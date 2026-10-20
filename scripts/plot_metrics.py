#!/usr/bin/env python3
"""Plot normalized protein-evaluation batch summaries. See README for CLI examples."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import statistics as stats
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "protein_eval_graph_mpl"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


BLUE = "#2563A6"
ORANGE = "#E87528"
RMSD_FIELDS = {"target", "binder", "placed", "irmsd"}


@dataclass(frozen=True)
class Metric:
    key: str
    label: str
    column: str
    unit: str
    higher: bool
    group: str

    @property
    def title(self) -> str:
        return f"{self.label} ({'↑' if self.higher else '↓'} better)"


METRICS = (
    Metric("target", "Target self-aligned backbone RMSD", "target_backbone_fit_rmsd_angstrom", "Å", False, "Folding"),
    Metric("binder", "Binder self-aligned backbone RMSD", "binder_self_aligned_backbone_rmsd_angstrom", "Å", False, "Folding"),
    Metric("lddt", "Global all-atom lDDT", "ost_lddt", "score", True, "Folding"),
    Metric("chirality", "Chirality accuracy", "chirality_all_accuracy_pct", "%", True, "Folding"),
    Metric("ilddt", "Interface lDDT", "ost_ilddt", "score", True, "Docking"),
    Metric("dockq", "DockQ", "ost_dockq_ave", "score", True, "Docking"),
    Metric("placed", "Target-aligned binder backbone RMSD", "target_aligned_binder_backbone_rmsd_angstrom", "Å", False, "Docking"),
    Metric("irmsd", "Interface RMSD", "ost_per_interface_json", "Å", False, "Docking"),
)
BY_KEY = {m.key: m for m in METRICS}
PREDICTION_FIELDS = {m.column for m in METRICS if m.key not in {"binder", "placed", "irmsd"}}


def number(value: object) -> float | None:
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def parse_source(spec: str) -> tuple[str, Path, str | None]:
    if "=" not in spec:
        raise ValueError("Use --source MODEL=CSV or MODEL=CSV::PATH_METHOD")
    model, tail = spec.split("=", 1)
    path, sep, path_method = tail.partition("::")
    if not model or not path or (sep and not path_method):
        raise ValueError(f"Invalid --source {spec!r}")
    return model, Path(path), path_method or None


def path_method(path: str) -> str | None:
    parts = Path(path).parts
    try:
        return parts[parts.index("predicted") + 1]
    except (ValueError, IndexError):
        return None


def seed_from_path(path: str, regex: str | None) -> str:
    if regex:
        match = re.search(regex, path)
        if not match or not match.groups():
            raise ValueError(f"Seed regex needs a capture group matching {path}")
        return match.group(1)
    for pattern in (r"(?:^|/)seed(\d+)(?:/|_)", r"(?:^|/)(\d{1,5})(?:/)", r"(?:^|[_/-])(\d{2,4})(?:[_./-]|$)"):
        match = re.search(pattern, path, re.I)
        if match:
            return match.group(1)
    return "single"


def read_confidence(path: Path | None) -> dict[str, float]:
    if path is None:
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or not {"prediction_path", "confidence_score"}.issubset(rows[0]):
        raise ValueError("Confidence CSV needs prediction_path and confidence_score columns")
    scores = {}
    for row in rows:
        value = number(row["confidence_score"])
        if value is None:
            raise ValueError(f"Invalid confidence score for {row['prediction_path']}")
        scores[row["prediction_path"]] = value
    return scores


def load_predictions(source_specs: list[str], seed_regex: str | None, confidence_csv: Path | None,
                     confidence_root: Path | None = None) -> list[dict]:
    confidence = read_confidence(confidence_csv)
    seen_models = set()
    predictions = []
    for spec in source_specs:
        model, path, method_filter = parse_source(spec)
        if model in seen_models:
            raise ValueError(f"Duplicate model label {model}")
        seen_models.add(model)
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        if not rows or not {"sample_id", "prediction_path", "status"}.issubset(rows[0]):
            raise ValueError(f"{path}: missing sample_id, prediction_path, or status")
        if method_filter is None:
            path_methods = {path_method(row["prediction_path"]) for row in rows}
            if len(path_methods) > 1 and not all(
                token and (re.fullmatch(r"seed\d+", token, re.I) or re.fullmatch(r".+_\d{2,5}", token))
                for token in path_methods
            ):
                raise ValueError(f"{path}: several path methods found ({sorted(str(x) for x in path_methods)}); add ::METHOD_TOKEN to --source")
        grouped = defaultdict(list)
        for row in rows:
            if method_filter and path_method(row["prediction_path"]) != method_filter:
                continue
            if row["status"].upper() == "SKIPPED":
                continue
            grouped[row["prediction_path"]].append(row)
        if not grouped:
            raise ValueError(f"{spec}: no included predictions")
        for pred_path, group in sorted(grouped.items()):
            first = group[0]
            if len({r["sample_id"] for r in group}) != 1:
                raise ValueError(f"Conflicting sample IDs for {pred_path}")
            item = {
                "model": model, "sample_id": first["sample_id"],
                "prediction_path": pred_path, "seed": seed_from_path(pred_path, seed_regex),
                "status": first["status"], "binder_groups": len(group),
            }
            for m in METRICS:
                if m.key in {"binder", "placed"}:
                    vals = [number(r.get(m.column)) for r in group]
                    finite = [v for v in vals if v is not None]
                    item[m.key] = stats.median(finite) if finite else None
                elif m.key == "irmsd":
                    raw = first.get(m.column)
                    interfaces = (json.loads(raw) if raw else None) or []
                    if not isinstance(interfaces, list):
                        raise ValueError(f"Invalid {m.column} for {pred_path}: expected an interface list")
                    vals = [number(i.get("irmsd")) for i in interfaces]
                    finite = [v for v in vals if v is not None]
                    item[m.key] = stats.median(finite) if finite else None
                else:
                    vals = {number(r.get(m.column)) for r in group}
                    if len(vals) != 1:
                        raise ValueError(f"Prediction-level {m.column} differs across binder rows: {pred_path}")
                    item[m.key] = vals.pop()
            for col in ("chirality_all_correct_count", "chirality_all_eligible_count", "chirality_all_skipped_count"):
                vals = {number(r.get(col)) for r in group}
                if len(vals) != 1:
                    raise ValueError(f"Prediction-level {col} differs across binder rows: {pred_path}")
                item[col] = vals.pop()
            if item["chirality"] is None:
                eligible = item["chirality_all_eligible_count"]
                skipped = item["chirality_all_skipped_count"] or 0
                correct = item["chirality_all_correct_count"]
                item["chirality"] = 100 * correct / (eligible - skipped) if eligible and correct is not None and eligible > skipped else None
            score = confidence.get(pred_path, number(first.get("confidence_score")))
            if score is None and confidence_root is not None:
                structure = confidence_root / pred_path
                candidate = structure.with_name(f"confidence_{structure.stem}.json")
                if candidate.is_file():
                    with candidate.open(encoding="utf-8") as handle:
                        score = number(json.load(handle).get("confidence_score"))
            item["confidence_score"] = score
            predictions.append(item)
    return predictions


def method_color(name: str, steering_names: set[str]) -> str:
    lower = name.lower()
    negative = any(token in lower for token in ("no_steer", "non_steer", "nonsteer", "unsteer", "baseline"))
    steering = name in steering_names or (not negative and ("steer" in lower or lower.endswith(("1x", "2x"))))
    return ORANGE if steering else BLUE


def selected_metric_keys(args: argparse.Namespace, predictions: list[dict]) -> list[str]:
    keys = args.metrics.split(",") if args.metrics else [m.key for m in METRICS]
    if len(keys) != len(set(keys)) or any(k not in BY_KEY for k in keys):
        raise ValueError(f"--metrics must be unique comma-separated keys from {', '.join(BY_KEY)}")
    if not keys:
        raise ValueError("No metrics selected")
    unavailable = [k for k in keys if not any(p[k] is not None for p in predictions)]
    if unavailable:
        raise ValueError(f"No data for metrics {unavailable}; use --metrics to omit unavailable panels")
    return keys


def require_values(points: list[dict], metric: str) -> list[dict]:
    available = [p for p in points if p[metric] is not None]
    if not available:
        raise ValueError(f"No values for {metric}")
    return available


def choose(candidates: list[dict], mode: str, metric: str) -> dict:
    if not candidates:
        raise ValueError("No candidates")
    if mode == "single":
        if len(candidates) != 1:
            raise ValueError(f"Expected one prediction; found {len(candidates)} for {candidates[0]['sample_id']} {candidates[0]['model']}. Choose top1 or oracle.")
        return candidates[0]
    if mode == "top1":
        if any(p["confidence_score"] is None for p in candidates):
            raise ValueError(f"Top-1 requires confidence_score for every seed of {candidates[0]['sample_id']} {candidates[0]['model']}")
        return min(candidates, key=lambda p: (-p["confidence_score"], p["seed"], p["prediction_path"]))
    if mode == "oracle":
        m = BY_KEY[metric]
        finite = require_values(candidates, metric)
        return min(finite, key=lambda p: (
            -p[metric] if m.higher else p[metric],
            -p["chirality"] if p["chirality"] is not None else 0,
            p["placed"] if p["placed"] is not None else float("inf"),
            p["seed"], p["prediction_path"],
        ))
    raise ValueError(f"Unknown selection {mode}")


def groups_by_input(predictions: list[dict]) -> dict[tuple[str, str], list[dict]]:
    groups = defaultdict(list)
    for p in predictions:
        groups[(p["model"], p["sample_id"])].append(p)
    return groups


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def save(fig, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
        "axes.spines.right": False, "axes.edgecolor": "#94A3B8",
        "axes.labelcolor": "#334155", "figure.facecolor": "white",
    })


def metric_axis(ax, metric: Metric) -> None:
    ax.set_title(metric.title, loc="left", fontsize=10, weight="bold")
    ax.set_ylabel(metric.unit)
    ax.grid(axis="y", color="#E2E8F0", linewidth=0.7, zorder=0)
    if metric.key in {"lddt", "ilddt", "dockq"}:
        ax.set_ylim(0, 1.05)
    elif metric.key == "chirality":
        ax.set_ylim(0, 105)
    else:
        ax.set_ylim(bottom=0)


def panels(n: int, title: str):
    cols = 4 if n > 6 else min(3, n)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(4.1 * cols, 4.0 * rows), squeeze=False)
    fig.suptitle(title, fontsize=15, weight="bold", y=1.01)
    for ax in axes.flat[n:]:
        ax.set_visible(False)
    return fig, list(axes.flat)


def describe_points(predictions: list[dict], keys: list[str], outdir: Path) -> None:
    rows = [{k: p.get(k) for k in ("model", "sample_id", "seed", "prediction_path", "status", "binder_groups", "confidence_score", "chirality_all_correct_count", "chirality_all_eligible_count", "chirality_all_skipped_count", *keys)} for p in predictions]
    write_csv(outdir / "normalized_predictions.csv", rows)


def compare(args, predictions, keys):
    points = [p for p in predictions if p["sample_id"] == args.input_id]
    if not points:
        raise ValueError(f"Input {args.input_id} not found")
    models = list(dict.fromkeys(p["model"] for p in predictions))
    groups = groups_by_input(points)
    fig, axes = panels(len(keys), f"{args.input_id}: folding and docking comparison")
    data = []
    for ax, key in zip(axes, keys):
        metric = BY_KEY[key]
        vals, labels, colors = [], [], []
        for model in models:
            candidates = groups.get((model, args.input_id), [])
            if not candidates:
                continue
            choice = choose(candidates, args.selection, key)
            if choice[key] is None:
                continue
            labels.append(model)
            vals.append(choice[key])
            colors.append(method_color(model, set(args.steering)))
            data.append({"model": model, "sample_id": args.input_id, "metric": key, "value": choice[key], "selection": args.selection, "seed": choice["seed"], "prediction_path": choice["prediction_path"]})
        if not vals:
            ax.set_visible(False)
            continue
        ax.bar(range(len(vals)), vals, color=colors, width=0.7, zorder=3)
        ax.set_xticks(range(len(vals)), labels, rotation=35, ha="right")
        metric_axis(ax, metric)
        if metric.key not in {"lddt", "ilddt", "dockq", "chirality"}:
            ax.set_ylim(0, max(vals) * 1.18 if max(vals) else 1)
    fig.tight_layout()
    save(fig, args.output_dir / f"{args.input_id}_comparison")
    write_csv(args.output_dir / "comparison_plot_data.csv", data)


def aggregate(chosen: list[dict], metric: str) -> tuple[float, float | None, int]:
    vals = [p[metric] for p in chosen if p[metric] is not None]
    if not vals:
        raise ValueError(f"No values for {metric}")
    if metric in RMSD_FIELDS:
        center = stats.median(vals)
    elif metric == "chirality":
        counts = [(p["chirality_all_correct_count"], p["chirality_all_eligible_count"], p["chirality_all_skipped_count"] or 0) for p in chosen]
        if all(c is not None and e is not None for c, e, _ in counts):
            den = sum(e - s for _, e, s in counts)
            center = 100 * sum(c for c, _, _ in counts) / den if den > 0 else stats.mean(vals)
        else:
            center = stats.mean(vals)
    else:
        center = stats.mean(vals)
    return center, stats.stdev(vals) if len(vals) > 1 else None, len(vals)


def selection_bars(args, predictions, keys):
    models = list(dict.fromkeys(p["model"] for p in predictions))
    grouped = groups_by_input(predictions)
    ids_by_model = {model: {p["sample_id"] for p in predictions if p["model"] == model} for model in models}
    ids = set.intersection(*ids_by_model.values()) if args.scope == "common" else set.union(*ids_by_model.values())
    if not ids:
        raise ValueError("No input IDs in the requested scope")
    fig, axes = panels(len(keys), f"Top-1 versus oracle · {args.scope} inputs")
    records, selections = [], []
    for ax, key in zip(axes, keys):
        metric = BY_KEY[key]
        bar_labels, centers, errors, colors = [], [], [], []
        for model in models:
            for mode in ("top1", "oracle"):
                chosen = []
                for sample_id in sorted(ids):
                    candidates = grouped.get((model, sample_id), [])
                    if not candidates:
                        continue
                    target_key = args.oracle_target if mode == "oracle" and args.oracle_target else key
                    item = choose(candidates, mode, target_key)
                    if item[key] is None:
                        continue
                    chosen.append(item)
                    selections.append({"metric": key, "model": model, "sample_id": sample_id, "selection": mode, "selection_metric": target_key, "seed": item["seed"], "prediction_path": item["prediction_path"], "value": item[key]})
                if not chosen:
                    continue
                center, sd, n = aggregate(chosen, key)
                label = f"{model}\n{mode}"
                bar_labels.append(label)
                centers.append(center)
                errors.append(sd or 0)
                colors.append(method_color(model, set(args.steering)))
                records.append({"metric": key, "model": model, "selection": mode, "scope": args.scope, "n_inputs": n, "value": center, "sample_sd": sd, "aggregation": "median" if key in RMSD_FIELDS else ("residue pooled" if key == "chirality" else "mean")})
        if not centers:
            ax.set_visible(False)
            continue
        x = list(range(len(centers)))
        ax.bar(x, centers, yerr=errors if args.error_bars else None, capsize=3 if args.error_bars else 0, color=colors, width=0.7, zorder=3)
        ax.set_xticks(x, bar_labels, rotation=35, ha="right")
        metric_axis(ax, metric)
        if key in RMSD_FIELDS:
            ax.set_ylim(0, max(c + (e if args.error_bars else 0) for c, e in zip(centers, errors)) * 1.18 or 1)
    fig.tight_layout()
    save(fig, args.output_dir / f"selection_bars_{args.scope}{'_sd' if args.error_bars else ''}")
    write_csv(args.output_dir / "selection_bar_summary.csv", records)
    write_csv(args.output_dir / "selection_plot_data.csv", selections)


def pearson(points: list[dict], key: str) -> float | None:
    pairs = [(p["chirality"], p[key]) for p in points if p["chirality"] is not None and p[key] is not None]
    if len(pairs) < 2:
        return None
    xs, ys = zip(*pairs)
    sx, sy = stats.pstdev(xs), stats.pstdev(ys)
    if sx == 0 or sy == 0:
        return None
    mx, my = stats.mean(xs), stats.mean(ys)
    return sum((x - mx) * (y - my) for x, y in pairs) / (len(pairs) * sx * sy)


def scatter(args, predictions, keys):
    models = list(dict.fromkeys(p["model"] for p in predictions))
    if len(models) > 2:
        raise ValueError("Scatter supports one or two models per figure")
    grouped = groups_by_input(predictions)
    fig, axes = panels(len(keys), f"Chirality accuracy versus structure · {args.selection}")
    fig._suptitle.set_y(1.075)
    plot_rows, pair_rows = [], []
    for ax, key in zip(axes, keys):
        selected = []
        if args.selection == "all":
            selected = [p for p in predictions if p[key] is not None and p["chirality"] is not None]
        else:
            for candidates in grouped.values():
                item = choose(candidates, args.selection, key)
                if item[key] is not None and item["chirality"] is not None:
                    selected.append(item)
        if not selected:
            ax.set_visible(False)
            continue
        for model in models:
            subset = [p for p in selected if p["model"] == model]
            ax.scatter([p["chirality"] for p in subset], [p[key] for p in subset], label=model,
                       color=method_color(model, set(args.steering)), alpha=0.75, s=35 if model == models[0] else 43,
                       marker="o" if model == models[0] else "^", zorder=3)
            r = pearson(subset, key)
            if r is not None:
                ax.text(0.98, 0.98 - models.index(model) * 0.07, f"{model}: r={r:.2f}", transform=ax.transAxes, ha="right", va="top", fontsize=8)
        if len(models) == 2:
            pair_keys = (lambda p: (p["sample_id"], p["seed"])) if args.selection == "all" else (lambda p: (p["sample_id"],))
            sides = [defaultdict(list) for _ in models]
            for p in selected:
                sides[models.index(p["model"])][pair_keys(p)].append(p)
            for pair_id in sorted(sides[0].keys() & sides[1].keys()):
                if len(sides[0][pair_id]) != 1 or len(sides[1][pair_id]) != 1:
                    raise ValueError(f"Ambiguous paired predictions for {pair_id}; choose one structure per input and seed")
                a, b = sides[0][pair_id][0], sides[1][pair_id][0]
                ax.plot([a["chirality"], b["chirality"]], [a[key], b[key]], linestyle=":", color="#64748B", alpha=0.47, zorder=2)
                pair_rows.append({"metric": key, "sample_id": a["sample_id"], "seed_1": a["seed"], "seed_2": b["seed"], "model_1": models[0], "model_2": models[1], "value_1": a[key], "value_2": b[key]})
        for p in selected:
            plot_rows.append({"metric": key, "model": p["model"], "sample_id": p["sample_id"], "seed": p["seed"], "selection": args.selection, "prediction_path": p["prediction_path"], "chirality_accuracy_pct": p["chirality"], "metric_value": p[key]})
        ax.set_xlim(-2, 102)
        ax.set_xlabel("Chirality accuracy (%)")
        ax.set_ylabel(BY_KEY[key].unit)
        ax.set_title(BY_KEY[key].title, loc="left", fontsize=10, weight="bold")
        ax.grid(color="#E2E8F0", linewidth=0.7, zorder=0)
        if key in {"lddt", "ilddt", "dockq"}:
            ax.set_ylim(-0.02, 1.02)
        elif key == "chirality":
            ax.set_ylim(-2, 102)
        else:
            ax.set_yscale("symlog", linthresh=2)
            ax.set_ylim(bottom=0)
    handles = [Line2D([0], [0], marker="o" if i == 0 else "^", linestyle="", color=method_color(model, set(args.steering)), label=model) for i, model in enumerate(models)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.035), ncol=len(models), frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    save(fig, args.output_dir / f"chirality_scatter_{args.selection}")
    write_csv(args.output_dir / "scatter_plot_data.csv", plot_rows)
    write_csv(args.output_dir / "scatter_pairs.csv", pair_rows)


def seed_range(args, predictions, keys):
    models = list(dict.fromkeys(p["model"] for p in predictions))
    if len(models) != 2:
        raise ValueError("Seed range requires exactly two sources")
    grouped = groups_by_input(predictions)
    ids = sorted({p["sample_id"] for p in predictions})
    rows = []
    for key in keys:
        fig, ax = plt.subplots(figsize=(12, max(5, len(ids) * 0.45 + 1.5)))
        for i, sample_id in enumerate(ids):
            for j, model in enumerate(models):
                points = [p for p in grouped.get((model, sample_id), []) if p[key] is not None]
                if not points:
                    continue
                if len({p["seed"] for p in points}) != len(points):
                    raise ValueError(f"Multiple predictions for a seed of {sample_id} {model}; select one structure per seed")
                y = i + (-0.16 if j == 0 else 0.16)
                values = [p[key] for p in points]
                color = method_color(model, set(args.steering))
                ax.plot([min(values), max(values)], [y, y], color=color, linewidth=2.5, zorder=2)
                ax.scatter(values, [y] * len(values), color=color, s=25, alpha=0.8, zorder=3)
                for p in points:
                    rows.append({"metric": key, "model": model, "sample_id": sample_id, "seed": p["seed"], "value": p[key], "min": min(values), "max": max(values), "prediction_path": p["prediction_path"]})
        ax.set_yticks(range(len(ids)), ids)
        ax.invert_yaxis()
        ax.set_xlabel(f"{BY_KEY[key].label} ({BY_KEY[key].unit})")
        ax.set_title(BY_KEY[key].title, loc="left", fontsize=12, weight="bold")
        ax.grid(axis="x", color="#E2E8F0")
        ax.legend(handles=[Line2D([0], [0], marker="o", color=method_color(model, set(args.steering)), label=model) for model in models], loc="upper right")
        fig.tight_layout()
        save(fig, args.output_dir / f"{key}_seed_minmax")
    write_csv(args.output_dir / "seed_range_plot_data.csv", rows)


def distribution(args, predictions, keys):
    models = list(dict.fromkeys(p["model"] for p in predictions))
    if len(models) != 2:
        raise ValueError("Distribution requires baseline and steering sources")
    baseline = [p for p in predictions if p["model"] == models[0] and p["sample_id"] == args.input_id]
    steering = [p for p in predictions if p["model"] == models[1] and p["sample_id"] == args.input_id]
    if len(baseline) < 2 or not steering:
        raise ValueError("Distribution needs at least two baseline seeds and one steering prediction for the input")
    fig, axes = panels(len(keys), f"{args.input_id}: {len(baseline)} unsteered seeds versus {len(steering)} steering predictions")
    rows = []
    for ax, key in zip(axes, keys):
        base = require_values(baseline, key)
        steer = require_values(steering, key)
        vals = [p[key] for p in base]
        violin = ax.violinplot(vals, positions=[0], widths=0.6, showextrema=False)
        for body in violin["bodies"]:
            body.set_facecolor(BLUE)
            body.set_alpha(0.25)
        ax.boxplot(vals, positions=[0], widths=0.17, showfliers=False, patch_artist=True,
                   boxprops={"facecolor": "white", "edgecolor": BLUE}, medianprops={"color": BLUE})
        ax.scatter([((i % 7) - 3) * 0.035 for i in range(len(vals))], vals, s=13, color=BLUE, alpha=0.52, zorder=3)
        ax.scatter([1 + (i - (len(steer) - 1) / 2) * 0.06 for i in range(len(steer))], [p[key] for p in steer], s=55, color=ORANGE, marker="D", zorder=4)
        for p in base + steer:
            rows.append({"metric": key, "model": p["model"], "sample_id": p["sample_id"], "seed": p["seed"], "value": p[key], "prediction_path": p["prediction_path"]})
        ax.set_xticks([0, 1], [f"{models[0]}\n(n={len(base)})", f"{models[1]}\n(n={len(steer)})"])
        ax.set_xlim(-0.5, 1.45)
        metric_axis(ax, BY_KEY[key])
        observed = vals + [p[key] for p in steer]
        low, high = min(observed), max(observed)
        span = max(high - low, 0.05 if key in {"lddt", "ilddt", "dockq"} else 0.1)
        cap = 103 if key == "chirality" else (1.02 if key in {"lddt", "ilddt", "dockq"} else float("inf"))
        ax.set_ylim(max(0, low - 0.14 * span), min(cap, high + 0.14 * span))
    fig.tight_layout()
    save(fig, args.output_dir / f"{args.input_id}_metric_distributions")
    write_csv(args.output_dir / "distribution_plot_data.csv", rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("compare", "selection-bars", "scatter", "seed-range", "distribution"):
        p = sub.add_parser(name)
        p.add_argument("--source", action="append", required=True, help="MODEL=CSV or MODEL=CSV::PATH_METHOD; repeat per model")
        p.add_argument("--output-dir", type=Path, required=True)
        p.add_argument("--metrics", help=f"Comma-separated keys; default all: {','.join(BY_KEY)}")
        p.add_argument("--steering", action="append", default=[], help="Model label to color orange; repeat as needed")
        p.add_argument("--seed-regex", help="Regex with one capture group for seed ID")
        p.add_argument("--confidence-csv", type=Path, help="CSV: prediction_path,confidence_score")
        p.add_argument("--confidence-root", type=Path, help="Evaluation directory containing adjacent Boltz confidence_<structure>.json files")
        if name in {"compare", "distribution"}:
            p.add_argument("--input-id", required=True)
        if name == "compare":
            p.add_argument("--selection", choices=("single", "top1", "oracle"), default="single")
        if name == "selection-bars":
            p.add_argument("--scope", choices=("common", "all"), default="common")
            p.add_argument("--oracle-target", choices=tuple(BY_KEY), help="Use one oracle metric across all panels; default metric-specific")
            p.add_argument("--error-bars", action="store_true", help="Show sample SD across inputs")
        if name == "scatter":
            p.add_argument("--selection", choices=("all", "top1", "oracle"), default="all")
    args = parser.parse_args()
    setup_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    predictions = load_predictions(args.source, args.seed_regex, args.confidence_csv, args.confidence_root)
    keys = selected_metric_keys(args, predictions)
    if args.command == "scatter" and not args.metrics:
        keys.remove("chirality")
    describe_points(predictions, keys, args.output_dir)
    {"compare": compare, "selection-bars": selection_bars, "scatter": scatter,
     "seed-range": seed_range, "distribution": distribution}[args.command](args, predictions, keys)
    manifest = {
        "command": args.command,
        "sources": args.source,
        "metric_keys": keys,
        "prediction_count": len(predictions),
        "predictions_by_model": {model: sum(p["model"] == model for p in predictions)
                                 for model in dict.fromkeys(p["model"] for p in predictions)},
        "missing_values_by_metric": {key: sum(p[key] is None for p in predictions) for key in keys},
        "settings": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
                     if k not in {"source", "output_dir"}},
    }
    (args.output_dir / "plot_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {args.command} plots for {len(predictions)} predictions to {args.output_dir}")


if __name__ == "__main__":
    main()

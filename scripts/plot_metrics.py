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
import textwrap
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "protein_eval_graph_mpl"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from rcg_style import (CONDITION_COLORS, DOUBLE_WIDTH, INK, PAIR_LINE, SINGLE_WIDTH,
                       SOURCE_URL, display_name, label_metric, letter_panels, model_color,
                       model_identity, model_marker, named_stem, ordered_models, panel_title, save_figure,
                       setup_style as rcg_setup_style)


BLUE = CONDITION_COLORS["baseline"]
ORANGE = CONDITION_COLORS["steering"]
RMSD_FIELDS = {"target", "binder", "placed", "irmsd"}
L_RMSD_CAPTION = (
    "L-RMSD denotes binder backbone (N, Cα, C, O) RMSD after jointly fitting all "
    "annotated target chains on their backbone atoms, without refitting the binder. "
)


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
    Metric("chirality", "D-Cα chirality correctness", "chirality_daa_accuracy_pct", "%", True, "Folding"),
    Metric("target", "Target self-aligned backbone RMSD", "target_backbone_fit_rmsd_angstrom", "Å", False, "Folding"),
    Metric("binder", "Binder self-aligned backbone RMSD", "binder_self_aligned_backbone_rmsd_angstrom", "Å", False, "Folding"),
    Metric("lddt", "Global heavy-atom lDDT", "ost_lddt", "score", True, "Folding"),
    Metric("placed", "L-RMSD", "target_aligned_binder_backbone_rmsd_angstrom", "Å", False, "Docking"),
    Metric("ilddt", "iLDDT", "target_binder_ilddt_macro", "score", True, "Docking"),
    Metric("dockq", "DockQ", "target_binder_dockq_mean", "score", True, "Docking"),
    Metric("irmsd", "Interface RMSD", "target_binder_dockq_per_interface_json", "Å", False, "Docking"),
)
BY_KEY = {m.key: m for m in METRICS}
DAA_COUNTS = ("chirality_daa_correct_count", "chirality_daa_eligible_count",
              "chirality_daa_skipped_count", "chirality_daa_assessable_count")


def metric_label_for_plot(args: argparse.Namespace, key: str) -> str:
    """Allow source-accurate lDDT wording without forking the plotter."""
    if key == "lddt" and args.lddt_label:
        return args.lddt_label
    if key == "ilddt" and args.ilddt_label:
        return args.ilddt_label
    return BY_KEY[key].label


def scatter_stat_name(model: str, models: list[str], steering_names: set[str]) -> str:
    """Use condition names in paired same-model panels; the legend has full names."""
    bases = {model_identity(name, steering_names)[0] for name in models}
    if len(models) == 2 and len(bases) == 1 and None not in bases:
        condition = model_identity(model, steering_names)[1]
        return {"baseline": "Baseline", "steering": "+ steering",
                "rcg": "+ RCG", "rcg_steering": "+ RCG + steering"}[condition]
    return display_name(model, steering_names)


def number(value: object) -> float | None:
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def parse_formats(raw: str) -> tuple[str, ...]:
    formats = tuple(dict.fromkeys(part.strip().lower() for part in raw.split(",")))
    if not formats or any(part not in {"svg", "png", "pdf"} for part in formats):
        raise argparse.ArgumentTypeError("--formats must be comma-separated svg,png,pdf values")
    return formats


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
                     confidence_root: Path | None = None, allow_single_method_seed_layout: bool = False,
                     stereochecks: str = "on") -> list[dict]:
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
            if len(path_methods) > 1 and not allow_single_method_seed_layout and not all(
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
            seed = seed_from_path(pred_path, seed_regex)
            if seed == "single" and allow_single_method_seed_layout:
                source_paths = {r.get("summary_path", "") for r in group}
                if len(source_paths) != 1 or not next(iter(source_paths)):
                    raise ValueError(f"Cannot recover seed from unique summary_path for {pred_path}")
                seed = seed_from_path(next(iter(source_paths)), seed_regex)
                if seed == "single":
                    raise ValueError(f"No seed in prediction_path or summary_path for {pred_path}")
            item = {
                "model": model, "sample_id": first["sample_id"],
                "prediction_path": pred_path, "seed": seed,
                "status": first["status"], "binder_groups": len(group),
                "warnings": first.get("warnings", ""),
                "unavailable_metrics_json": first.get("unavailable_metrics_json", ""),
            }
            for m in METRICS:
                if m.key == "chirality":
                    continue
                column = m.column + "_no_stereochecks" if (
                    stereochecks == "off" and m.key in {"lddt", "ilddt"}
                ) else m.column
                if m.key in {"binder", "placed"}:
                    vals = [number(r.get(column)) for r in group]
                    finite = [v for v in vals if v is not None]
                    item[m.key] = stats.mean(finite) if finite else None
                elif m.key == "irmsd":
                    raw = first.get(column)
                    interfaces = (json.loads(raw) if raw else None) or []
                    if not isinstance(interfaces, list):
                        raise ValueError(f"Invalid {m.column} for {pred_path}: expected an interface list")
                    vals = [number(i.get("irmsd")) for i in interfaces]
                    finite = [v for v in vals if v is not None]
                    item[m.key] = stats.mean(finite) if finite else None
                else:
                    vals = {number(r.get(column)) for r in group}
                    if len(vals) != 1:
                        raise ValueError(f"Prediction-level {column} differs across binder rows: {pred_path}")
                    item[m.key] = vals.pop()
            for col in (*DAA_COUNTS, "chirality_all_correct_count", "chirality_all_eligible_count",
                        "chirality_all_skipped_count", "chirality_all_assessable_count"):
                vals = {number(r.get(col)) for r in group}
                if len(vals) != 1:
                    raise ValueError(f"Prediction-level {col} differs across binder rows: {pred_path}")
                item[col] = vals.pop()
            if item["chirality_daa_assessable_count"] is None:
                eligible = item["chirality_daa_eligible_count"]
                skipped = item["chirality_daa_skipped_count"]
                if eligible is not None and skipped is not None:
                    item["chirality_daa_assessable_count"] = eligible - skipped
            item["chirality"] = item["daa_chirality"] = pooled_daa_accuracy([item])
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


def method_color(name: str, steering_names: set[str], all_models: list[str] | None = None) -> str:
    return model_color(name, all_models or [name], steering_names)


def selected_metric_keys(args: argparse.Namespace, predictions: list[dict]) -> list[str]:
    keys = args.metrics.split(",") if args.metrics else [m.key for m in METRICS]
    if len(keys) != len(set(keys)) or any(k not in BY_KEY for k in keys):
        raise ValueError(f"--metrics must be unique comma-separated keys from {', '.join(BY_KEY)}")
    if not keys:
        raise ValueError("No metrics selected")
    unavailable = [k for k in keys if not any(p[k] is not None for p in predictions)]
    if unavailable:
        raise ValueError(f"No data for metrics {unavailable}; use --metrics to omit unavailable panels")
    return ordered_panel_keys(keys, BY_KEY)


def ordered_panel_keys(keys: list[str], metric_by_key: dict[str, Metric]) -> list[str]:
    """Put the requested lead metric first in each folding/docking row."""
    folding = [key for key in keys if metric_by_key[key].group == "Folding"]
    docking = [key for key in keys if metric_by_key[key].group == "Docking"]
    folding = [key for key in folding if key in {"chirality", "daa_chirality"}] + [
        key for key in folding if key not in {"chirality", "daa_chirality"}]
    docking = [key for key in docking if key == "placed"] + [key for key in docking if key != "placed"]
    return folding + docking


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


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if not rows and fields is None:
        return
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{k: "NA" if v is None else v for k, v in row.items()} for row in rows])


def save(fig, stem: Path, formats: tuple[str, ...]) -> None:
    save_figure(fig, stem, formats)


def setup_style() -> None:
    rcg_setup_style()


def metric_axis(ax, metric: Metric) -> None:
    ax.set_title(panel_title(metric.label, metric.higher), loc="left", fontsize=9, weight="bold")
    ax.set_ylabel(label_metric(metric.label, metric.unit))
    if metric.key in {"lddt", "ilddt", "dockq"}:
        ax.set_ylim(0, 1.05)
    elif metric.key == "chirality":
        ax.set_ylim(0, 105)
    else:
        ax.set_ylim(bottom=0)


def selection_axis_label(metric: Metric) -> str:
    """Name the across-input aggregation used for a selection bar."""
    if metric.key in RMSD_FIELDS:
        return f"Median {label_metric(metric.label, metric.unit)}"
    if metric.key == "chirality":
        return "Pooled D-Cα chirality correctness (%)"
    return f"Mean {label_metric(metric.label, metric.unit)}"


def panels(keys: list[str], title: str, metric_by_key: dict[str, Metric] | None = None,
           discussion_title: bool = False):
    metric_by_key = metric_by_key or BY_KEY
    ordered = ordered_panel_keys(keys, metric_by_key)
    folding = [key for key in ordered if metric_by_key[key].group == "Folding"]
    docking = [key for key in ordered if metric_by_key[key].group == "Docking"]
    cols = max(len(folding), len(docking), 1)
    rows = 2 if len(keys) > 1 else 1
    fig, axes = plt.subplots(rows, cols, figsize=(DOUBLE_WIDTH if rows == 2 else SINGLE_WIDTH,
                                               5.7 if rows == 2 else SINGLE_WIDTH), squeeze=False)
    if discussion_title:
        fig.suptitle(title, fontsize=11, weight="bold", ha="left", x=0.04)
    mapped = {}
    for row, group in enumerate((folding, docking) if rows == 2 else (keys,)):
        for col, key in enumerate(group):
            mapped[key] = axes[row, col]
    for ax in axes.flat:
        if ax not in mapped.values():
            ax.set_visible(False)
    letter_panels([mapped[key] for key in (folding + docking) if key in mapped])
    return fig, mapped


def describe_points(predictions: list[dict], keys: list[str], outdir: Path) -> None:
    rows = [{k: p.get(k) for k in ("model", "sample_id", "seed", "prediction_path", "status", "binder_groups", "confidence_score", *DAA_COUNTS, "warnings", "unavailable_metrics_json", *keys)} for p in predictions]
    write_csv(outdir / "normalized_predictions.csv", rows)
    missing = [{"model": p["model"], "sample_id": p["sample_id"], "seed": p["seed"],
                "prediction_path": p["prediction_path"], "metric": key,
                "reason": p["unavailable_metrics_json"] or p["warnings"] or "metric value unavailable"}
               for p in predictions for key in keys if p[key] is None]
    write_csv(outdir / "missing_values.csv", missing,
              ["model", "sample_id", "seed", "prediction_path", "metric", "reason"])


def describe_skipped(source_specs: list[str], outdir: Path) -> None:
    skipped = []
    for spec in source_specs:
        model, path, method_filter = parse_source(spec)
        with path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if method_filter and path_method(row["prediction_path"]) != method_filter:
                    continue
                if row["status"].upper() == "SKIPPED":
                    skipped.append({"model": model, "sample_id": row["sample_id"],
                                    "prediction_path": row["prediction_path"],
                                    "reason": row.get("warnings") or row.get("unavailable_metrics_json") or "evaluation skipped"})
    unique = list({(r["model"], r["prediction_path"]): r for r in skipped}.values())
    write_csv(outdir / "skipped_predictions.csv", unique,
              ["model", "sample_id", "prediction_path", "reason"])


def compare(args, predictions, keys):
    points = [p for p in predictions if p["sample_id"] == args.input_id]
    if not points:
        raise ValueError(f"Input {args.input_id} not found")
    models = ordered_models(list(dict.fromkeys(p["model"] for p in predictions)), set(args.steering))
    groups = groups_by_input(points)
    fig, axes = panels(keys, f"{args.input_id}: folding and docking comparison", discussion_title=args.discussion_title)
    data = []
    for key in keys:
        ax = axes[key]
        metric = BY_KEY[key]
        vals, labels, colors = [], [], []
        for model in models:
            candidates = groups.get((model, args.input_id), [])
            if not candidates:
                data.append({"model": model, "sample_id": args.input_id, "metric": key, "value": None,
                             "selection": args.selection, "selection_metric": key, "reason": "input absent"})
                continue
            selection_metric = args.oracle_target if args.selection == "oracle" and args.oracle_target else key
            choice = choose(candidates, args.selection, selection_metric)
            if choice[key] is None:
                data.append({"model": model, "sample_id": args.input_id, "metric": key, "value": None,
                             "selection": args.selection, "selection_metric": selection_metric,
                             "seed": choice["seed"], "prediction_path": choice["prediction_path"], "reason": "metric unavailable"})
                continue
            labels.append(textwrap.fill(display_name(model, set(args.steering)), 13))
            vals.append(choice[key])
            colors.append(method_color(model, set(args.steering), models))
            data.append({"model": model, "sample_id": args.input_id, "metric": key, "value": choice[key], "selection": args.selection, "selection_metric": selection_metric, "seed": choice["seed"], "prediction_path": choice["prediction_path"]})
        if not vals:
            ax.set_visible(False)
            continue
        ax.bar(range(len(vals)), vals, color=colors, width=0.7, zorder=3)
        ax.set_xticks(range(len(vals)), labels)
        metric_axis(ax, metric)
        if metric.key not in {"lddt", "ilddt", "dockq", "chirality"}:
            ax.set_ylim(0, max(vals) * 1.18 if max(vals) else 1)
    fig.tight_layout()
    save(fig, named_stem(args.output_dir, f"{args.input_id}_comparison", args.figure_id, args.version), args.formats)
    write_csv(args.output_dir / "comparison_plot_data.csv", data)


def pooled_daa_accuracy(points: list[dict]) -> float | None:
    counts = [(p.get("chirality_daa_correct_count"), p.get("chirality_daa_assessable_count"))
              for p in points]
    valid = [(passed, assessable) for passed, assessable in counts
             if passed is not None and assessable and 0 <= passed <= assessable]
    total = sum(assessable for _, assessable in valid)
    return 100 * sum(passed for passed, _ in valid) / total if total else None


def aggregate(chosen: list[dict], metric: str) -> tuple[float, float | None, int]:
    vals = [p[metric] for p in chosen if p[metric] is not None]
    if not vals:
        raise ValueError(f"No values for {metric}")
    if metric in RMSD_FIELDS:
        center = stats.median(vals)
    elif metric == "chirality":
        center = pooled_daa_accuracy(chosen)
        if center is None:
            raise ValueError("D-Cα chirality requires assessable residue counts")
    else:
        center = stats.mean(vals)
    return center, stats.stdev(vals) if len(vals) > 1 else None, len(vals)


def selection_bars(args, predictions, keys):
    models = ordered_models(list(dict.fromkeys(p["model"] for p in predictions)), set(args.steering))
    grouped = groups_by_input(predictions)
    ids_by_model = {model: {p["sample_id"] for p in predictions if p["model"] == model} for model in models}
    ids = set.intersection(*ids_by_model.values()) if args.scope == "common" else set.union(*ids_by_model.values())
    if not ids:
        raise ValueError("No input IDs in the requested scope")
    fig, axes = panels(keys, f"Top-1 versus oracle · {args.scope} inputs", discussion_title=args.discussion_title)
    records, selections = [], []
    for key in keys:
        ax = axes[key]
        metric = BY_KEY[key]
        bar_labels, centers, errors, colors = [], [], [], []
        for model in models:
            for mode in ("top1", "oracle"):
                chosen = []
                for sample_id in sorted(ids):
                    candidates = grouped.get((model, sample_id), [])
                    if not candidates:
                        selections.append({"metric": key, "model": model, "sample_id": sample_id,
                                           "selection": mode, "value": None, "reason": "input absent"})
                        continue
                    target_key = args.oracle_target if mode == "oracle" and args.oracle_target else key
                    item = choose(candidates, mode, target_key)
                    if item[key] is None:
                        selections.append({"metric": key, "model": model, "sample_id": sample_id,
                                           "selection": mode, "selection_metric": target_key,
                                           "seed": item["seed"], "prediction_path": item["prediction_path"],
                                           "value": None, "reason": "metric unavailable"})
                        continue
                    chosen.append(item)
                    selections.append({"metric": key, "model": model, "sample_id": sample_id, "selection": mode, "selection_metric": target_key, "seed": item["seed"], "prediction_path": item["prediction_path"], "value": item[key]})
                if not chosen:
                    continue
                center, sd, n = aggregate(chosen, key)
                bar_labels.append("Top-1" if mode == "top1" else "Oracle")
                centers.append(center)
                errors.append(sd or 0)
                colors.append(method_color(model, set(args.steering), models))
                aggregation = {
                    **dict.fromkeys(RMSD_FIELDS, "median"), "chirality": "pooled assessable D sites"
                }.get(key, "system-equal mean")
                records.append({"metric": key, "model": model, "selection": mode, "scope": args.scope, "n_inputs": n, "value": center, "sample_sd": sd, "aggregation": aggregation})
        if not centers:
            ax.set_visible(False)
            continue
        x = list(range(len(centers)))
        ax.bar(x, centers, yerr=errors if args.error_bars else None, capsize=3 if args.error_bars else 0, color=colors, width=0.7, zorder=3)
        ax.set_xticks(x, bar_labels)
        metric_axis(ax, metric)
        ax.set_ylabel(selection_axis_label(metric))
        if key in RMSD_FIELDS:
            ax.set_ylim(0, max(c + (e if args.error_bars else 0) for c, e in zip(centers, errors)) * 1.18 or 1)
    fig.legend(handles=[Patch(facecolor=method_color(model, set(args.steering), models), label=display_name(model, set(args.steering))) for model in models],
               loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=len(models), frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, named_stem(args.output_dir, f"selection_bars_{args.scope}{'_sd' if args.error_bars else ''}", args.figure_id, args.version), args.formats)
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
    models = ordered_models(list(dict.fromkeys(p["model"] for p in predictions)), set(args.steering))
    if len(models) > 2:
        raise ValueError("Scatter supports one or two models per figure")
    grouped = groups_by_input(predictions)
    fig, axes = panels(keys, f"D-Cα chirality correctness versus structure · {args.selection}", discussion_title=args.discussion_title)
    # Reserve a separate band for correlations above each plotting area.
    fig.set_size_inches(DOUBLE_WIDTH if len(keys) > 1 else SINGLE_WIDTH,
                        6.4 if len(keys) > 1 else 4.1)
    plot_rows, pair_rows, unmatched_rows = [], [], []
    stat_artists = []
    for key in keys:
        ax = axes[key]
        selected = []
        if args.selection == "all":
            selected = [p for p in predictions if p[key] is not None and p["chirality"] is not None]
        else:
            for candidates in grouped.values():
                selection_metric = args.oracle_target if args.selection == "oracle" and args.oracle_target else key
                item = choose(candidates, args.selection, selection_metric)
                if item[key] is not None and item["chirality"] is not None:
                    selected.append(item)
        paired = {}
        if len(models) == 2:
            pair_key = (lambda p: (p["sample_id"], p["seed"])) if args.selection == "all" else (lambda p: (p["sample_id"],))
            sides = [defaultdict(list) for _ in models]
            for p in selected:
                sides[models.index(p["model"])][pair_key(p)].append(p)
            common = sides[0].keys() & sides[1].keys()
            for p in selected:
                if pair_key(p) not in common:
                    unmatched_rows.append({"metric": key, "model": p["model"], "sample_id": p["sample_id"], "seed": p["seed"], "reason": "no matched counterpart"})
            selected = [p for p in selected if pair_key(p) in common]
            for pair_id in sorted(common):
                if len(sides[0][pair_id]) != 1 or len(sides[1][pair_id]) != 1:
                    raise ValueError(f"Ambiguous paired predictions for {pair_id}; choose one structure per input and seed")
                paired[pair_id] = (sides[0][pair_id][0], sides[1][pair_id][0])
        if not selected:
            ax.set_visible(False)
            continue
        correlations = []
        for model in models:
            subset = [p for p in selected if p["model"] == model]
            ax.scatter([p["chirality"] for p in subset], [p[key] for p in subset], label=display_name(model, set(args.steering)),
                       color=method_color(model, set(args.steering), models), alpha=0.85, s=12.25,
                       marker=model_marker(model), zorder=3)
            r = pearson(subset, key)
            if r is not None:
                correlations.append(f"{scatter_stat_name(model, models, set(args.steering))}: r={r:.2f}")
            else:
                correlations.append(f"{scatter_stat_name(model, models, set(args.steering))}: r=NA")
        if len(models) == 2:
            for pair_id, (a, b) in paired.items():
                ax.plot([a["chirality"], b["chirality"]], [a[key], b[key]], linestyle=":", color=PAIR_LINE, linewidth=0.6, zorder=2)
                pair_rows.append({"metric": key, "sample_id": a["sample_id"], "seed_1": a["seed"], "seed_2": b["seed"], "model_1": models[0], "model_2": models[1], "value_1": a[key], "value_2": b[key]})
        for p in selected:
            plot_rows.append({"metric": key, "model": p["model"], "sample_id": p["sample_id"], "seed": p["seed"], "selection": args.selection, "selection_metric": args.oracle_target if args.selection == "oracle" and args.oracle_target else key, "prediction_path": p["prediction_path"], "d_ca_chirality_correctness_pct": p["chirality"], "metric_value": p[key]})
        ax.set_xlim(-2, 102)
        ax.set_xlabel("D-Cα chirality correctness (%)")
        label = metric_label_for_plot(args, key)
        ax.set_ylabel(label_metric(label, BY_KEY[key].unit))
        title = ax.set_title(panel_title(label, BY_KEY[key].higher),
                             loc="left", fontsize=9, weight="bold", pad=31)
        if correlations:
            stat = ax.text(0.5, 1.035, "\n".join(correlations),
                           transform=ax.transAxes, ha="center", va="bottom",
                           fontsize=8, clip_on=False)
            stat_artists.append((ax, stat, title))
        if key in {"lddt", "ilddt", "dockq"}:
            ax.set_ylim(-0.02, 1.02)
        elif key == "chirality":
            ax.set_ylim(-2, 102)
        else:
            ax.set_ylim(bottom=0)
    handles = [Line2D([0], [0], marker=model_marker(model), linestyle="", color=method_color(model, set(args.steering), models), label=display_name(model, set(args.steering))) for model in models]
    if not plot_rows:
        raise ValueError("No matched observations with both D-Cα correctness and selected metrics")
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=len(models), frameon=False)
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    stat_boxes = [stat.get_window_extent(renderer) for _, stat, _ in stat_artists]
    title_boxes = [title.get_window_extent(renderer) for _, _, title in stat_artists]
    for index, (ax, stat, title) in enumerate(stat_artists):
        bbox = stat_boxes[index]
        if (not fig.bbox.contains(*bbox.p0) or not fig.bbox.contains(*bbox.p1)
                or any(bbox.overlaps(other) for other in title_boxes)
                or any(bbox.overlaps(other) for other_index, other in enumerate(stat_boxes)
                       if other_index != index)
                or any(bbox.overlaps(other.bbox) for other in axes.values()
                       if other.get_visible())):
            raise RuntimeError("Scatter correlation text overlaps another label or plotting area; increase layout spacing")
    save(fig, named_stem(args.output_dir, f"chirality_scatter_{args.selection}", args.figure_id, args.version), args.formats)
    write_csv(args.output_dir / "scatter_plot_data.csv", plot_rows)
    write_csv(args.output_dir / "scatter_pairs.csv", pair_rows)
    write_csv(args.output_dir / "scatter_unmatched.csv", unmatched_rows,
              ["metric", "model", "sample_id", "seed", "reason"])


def seed_range(args, predictions, keys):
    models = ordered_models(list(dict.fromkeys(p["model"] for p in predictions)), set(args.steering))
    if len(models) != 2:
        raise ValueError("Seed range requires exactly two sources")
    grouped = groups_by_input(predictions)
    ids_by_model = {model: {p["sample_id"] for p in predictions if p["model"] == model} for model in models}
    ids = sorted(ids_by_model[models[0]] & ids_by_model[models[1]])
    if not ids:
        raise ValueError("Seed range requires inputs evaluated by both models")
    unmatched = [{"model": model, "sample_id": sample_id, "reason": "no matched input in other model"}
                 for model in models for sample_id in sorted(ids_by_model[model] - set(ids))]
    write_csv(args.output_dir / "seed_range_unmatched.csv", unmatched, ["model", "sample_id", "reason"])
    rows = []
    for key in keys:
        fig, ax = plt.subplots(figsize=(DOUBLE_WIDTH, max(SINGLE_WIDTH, len(ids) * 0.19 + 1.0)))
        for i, sample_id in enumerate(ids):
            for j, model in enumerate(models):
                points = [p for p in grouped.get((model, sample_id), []) if p[key] is not None]
                if not points:
                    continue
                if len({p["seed"] for p in points}) != len(points):
                    raise ValueError(f"Multiple predictions for a seed of {sample_id} {model}; select one structure per seed")
                y = i + (-0.16 if j == 0 else 0.16)
                values = [p[key] for p in points]
                color = method_color(model, set(args.steering), models)
                ax.plot([min(values), max(values)], [y, y], color=color, linewidth=1, zorder=2)
                ax.scatter(values, [y] * len(values), color=color, s=12.25, alpha=0.85, marker=model_marker(model), zorder=3)
                for p in points:
                    rows.append({"metric": key, "model": model, "sample_id": sample_id, "seed": p["seed"], "value": p[key], "min": min(values), "max": max(values), "prediction_path": p["prediction_path"]})
        ax.set_yticks(range(len(ids)), ids)
        ax.invert_yaxis()
        ax.set_xlabel(label_metric(BY_KEY[key].label, BY_KEY[key].unit))
        ax.set_title(panel_title(BY_KEY[key].label, BY_KEY[key].higher), loc="left", fontsize=9, weight="bold")
        ax.legend(handles=[Line2D([0], [0], marker=model_marker(model), color=method_color(model, set(args.steering), models), label=display_name(model, set(args.steering))) for model in models], loc="upper right", frameon=False)
        fig.tight_layout()
        save(fig, named_stem(args.output_dir, f"{key}_seed_minmax", args.figure_id, args.version), args.formats)
    write_csv(args.output_dir / "seed_range_plot_data.csv", rows)


def distribution(args, predictions, keys):
    models = list(dict.fromkeys(p["model"] for p in predictions))
    if len(models) != 2:
        raise ValueError("Distribution requires baseline and steering sources")
    if model_identity(models[0], set(args.steering))[1] != "baseline" or model_identity(models[1], set(args.steering))[1] != "steering":
        raise ValueError("Distribution sources must be baseline first and steering second; label models accordingly")
    baseline = [p for p in predictions if p["model"] == models[0] and p["sample_id"] == args.input_id]
    steering = [p for p in predictions if p["model"] == models[1] and p["sample_id"] == args.input_id]
    if len(baseline) < 2 or not steering:
        raise ValueError("Distribution needs at least two baseline seeds and one steering prediction for the input")
    fig, axes = panels(keys, f"{args.input_id}: {len(baseline)} unsteered seeds versus {len(steering)} steering predictions", discussion_title=args.discussion_title)
    base_color = method_color(models[0], set(args.steering), models)
    steering_color = method_color(models[1], set(args.steering), models)
    rows = []
    for key in keys:
        ax = axes[key]
        base = require_values(baseline, key)
        steer = require_values(steering, key)
        vals = [p[key] for p in base]
        violin = ax.violinplot(vals, positions=[0], widths=0.6, showextrema=False)
        for body in violin["bodies"]:
            body.set_facecolor(base_color)
            body.set_alpha(0.25)
        ax.boxplot(vals, positions=[0], widths=0.17, showfliers=False, patch_artist=True,
                   boxprops={"facecolor": "white", "edgecolor": base_color}, medianprops={"color": INK})
        ax.scatter([((i % 7) - 3) * 0.035 for i in range(len(vals))], vals, s=12.25, color=base_color, alpha=0.85, marker=model_marker(models[0]), zorder=3)
        ax.scatter([1 + (i - (len(steer) - 1) / 2) * 0.06 for i in range(len(steer))], [p[key] for p in steer], s=12.25, color=steering_color, alpha=0.85, marker=model_marker(models[1]), zorder=4)
        for p in base + steer:
            rows.append({"metric": key, "model": p["model"], "sample_id": p["sample_id"], "seed": p["seed"], "value": p[key], "prediction_path": p["prediction_path"]})
        ax.set_xticks([0, 1], [f"{display_name(models[0], set(args.steering))}\n(n={len(base)})", f"{display_name(models[1], set(args.steering))}\n(n={len(steer)})"])
        ax.set_xlim(-0.5, 1.45)
        metric_axis(ax, BY_KEY[key])
        observed = vals + [p[key] for p in steer]
        low, high = min(observed), max(observed)
        span = max(high - low, 0.05 if key in {"lddt", "ilddt", "dockq"} else 0.1)
        cap = 103 if key == "chirality" else (1.02 if key in {"lddt", "ilddt", "dockq"} else float("inf"))
        ax.set_ylim(max(0, low - 0.14 * span), min(cap, high + 0.14 * span))
    fig.tight_layout()
    save(fig, named_stem(args.output_dir, f"{args.input_id}_metric_distributions", args.figure_id, args.version), args.formats)
    write_csv(args.output_dir / "distribution_plot_data.csv", rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("compare", "selection-bars", "scatter", "seed-range", "distribution"):
        p = sub.add_parser(name)
        p.add_argument("--source", action="append", required=True, help="MODEL=CSV or MODEL=CSV::PATH_METHOD; repeat per model")
        p.add_argument("--output-dir", type=Path, required=True)
        p.add_argument("--formats", type=parse_formats, default=("svg",), help="Comma-separated output formats: svg,png,pdf (default: svg)")
        p.add_argument("--figure-id", help="Final manuscript figure number, e.g. Fig05 (adds figure/content/version filename)")
        p.add_argument("--version", default="v01", help="Final manuscript figure version (default: v01)")
        p.add_argument("--discussion-title", action="store_true", help="Add an 11 pt whole-figure title for standalone discussion figures")
        p.add_argument("--metrics", help=f"Comma-separated keys; default all: {','.join(BY_KEY)}")
        p.add_argument("--stereochecks", choices=("on", "off"), default="on",
                       help="Choose the saved OpenStructure lDDT and iLDDT versions")
        p.add_argument("--lddt-label", help="Source-accurate global lDDT label for scatter titles and axes")
        p.add_argument("--ilddt-label", help="Source-accurate interface lDDT label for scatter titles and axes")
        p.add_argument("--steering", action="append", default=[], help="Model label to color orange; repeat as needed")
        p.add_argument("--seed-regex", help="Regex with one capture group for seed ID")
        p.add_argument("--allow-single-method-seed-layout", action="store_true",
                       help="Treat a source CSV as one method despite seed directories or legacy path aliases; recover missing seed IDs from summary_path")
        p.add_argument("--confidence-csv", type=Path, help="CSV: prediction_path,confidence_score")
        p.add_argument("--confidence-root", type=Path, help="Evaluation directory containing adjacent Boltz confidence_<structure>.json files")
        if name in {"compare", "selection-bars", "scatter"}:
            p.add_argument("--oracle-target", choices=tuple(BY_KEY), help="Use one oracle selection metric across panels; default metric-specific oracle")
        if name in {"compare", "distribution"}:
            p.add_argument("--input-id", required=True)
        if name == "compare":
            p.add_argument("--selection", choices=("single", "top1", "oracle"), default="single")
        if name == "selection-bars":
            p.add_argument("--scope", choices=("common", "all"), default="common")
            p.add_argument("--error-bars", action="store_true", help="Show sample SD across inputs")
        if name == "scatter":
            p.add_argument("--selection", choices=("all", "top1", "oracle"), default="all")
    args = parser.parse_args()
    setup_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    predictions = load_predictions(args.source, args.seed_regex, args.confidence_csv,
                                   args.confidence_root, args.allow_single_method_seed_layout,
                                   args.stereochecks)
    keys = selected_metric_keys(args, predictions)
    if args.command == "scatter":
        keys = [key for key in keys if key != "chirality"]
        if not keys:
            parser.error("Scatter uses chirality on the X axis; select at least one other metric for a Y-axis panel")
    describe_points(predictions, keys, args.output_dir)
    describe_skipped(args.source, args.output_dir)
    {"compare": compare, "selection-bars": selection_bars, "scatter": scatter,
     "seed-range": seed_range, "distribution": distribution}[args.command](args, predictions, keys)
    manifest = {
        "command": args.command,
        "formatting_standard": SOURCE_URL,
        "sources": args.source,
        "metric_keys": keys,
        "prediction_count": len(predictions),
        "predictions_by_model": {model: sum(p["model"] == model for p in predictions)
                                 for model in dict.fromkeys(p["model"] for p in predictions)},
        "missing_values_by_metric": {key: sum(p[key] is None for p in predictions) for key in keys},
        "observation_unit": {"compare": "one selected prediction per system/model/metric",
                             "selection-bars": "one selected prediction per system/model/selection/metric",
                             "scatter": "system-seed pair" if getattr(args, "selection", None) == "all" else "one selected prediction per system/model/metric",
                             "seed-range": "system-seed pair", "distribution": "seed of one system"}[args.command],
        "settings": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
                     if k not in {"source", "output_dir"}},
    }
    (args.output_dir / "plot_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    model_ids = {model: {p["sample_id"] for p in predictions if p["model"] == model}
                 for model in manifest["predictions_by_model"]}
    matched = len(set.intersection(*model_ids.values())) if model_ids else 0
    selection = getattr(args, "selection", "top-1 versus metric-wise oracle" if args.command == "selection-bars" else "all evaluated seeds")
    if getattr(args, "oracle_target", None):
        oracle_note = f"Shared {args.oracle_target} oracle across panels"
    elif args.command == "selection-bars" or selection == "oracle":
        oracle_note = "Metric-wise oracle panels select independently; they are not a joint profile from one structure"
    else:
        oracle_note = "Each plotted multi-metric observation uses values from its stated prediction file"
    top1_note = " Top-1 means highest supplied confidence score." if args.command == "selection-bars" or selection == "top1" else ""
    caption = (
        f"# Caption draft — complete and verify before manuscript use\n\n"
        f"Figure {args.figure_id[3:] if args.figure_id else 'X'}. [Main finding or comparison].\n"
        f"Panels: {', '.join(metric_label_for_plot(args, key) for key in keys)}. "
        f"Observation unit: {manifest['observation_unit']}. "
        f"Input IDs present in every model before metric-specific exclusions: {matched}. "
        f"Selection: {selection}.{top1_note} {oracle_note}.\n\n"
        + (L_RMSD_CAPTION if "placed" in keys else "")
        + "[Verify source metric definitions, N/Cα/C/O alignment atoms, chain mappings, "
        "DockQ mode/interface scope, and lDDT implementation before stating them.] "
        f"OpenStructure lDDT and iLDDT use stereochemistry checks {args.stereochecks}. "
        "D-Cα correctness uses passed/assessable D sites at the evaluation cutoff. "
        "Missing metrics are NA in normalized_predictions.csv with reasons in missing_values.csv; "
        "skipped evaluations are in skipped_predictions.csv. "
        + ("Error bars show sample SD across selected system-level values. " if getattr(args, "error_bars", False) else "")
        + "[State metric-specific matched N, seeds per system, "
        "and any evaluation restrictions. For boxplots: median, IQR, 1.5×IQR whiskers, all raw points.]\n"
    )
    (args.output_dir / "caption_draft.md").write_text(caption, encoding="utf-8")
    print(f"Wrote {args.command} plots for {len(predictions)} predictions to {args.output_dir}")


if __name__ == "__main__":
    main()

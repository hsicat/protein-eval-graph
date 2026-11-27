#!/usr/bin/env python3
"""Compare metric-wise seed oracles across Single D, Partial D, and All D inputs."""

from __future__ import annotations

import argparse
import csv
import statistics as stats
from collections import defaultdict
from pathlib import Path

from plot_metrics import (
    BY_KEY as BASE_METRICS, L_RMSD_CAPTION, Metric, load_predictions, method_color,
    ordered_panel_keys, panels, parse_formats, pooled_daa_accuracy, save, setup_style, write_csv,
)
from matplotlib.patches import Patch
from rcg_style import INK, display_name, model_marker, named_stem, ordered_models, panel_title


GROUP_ORDER = ("Single D", "Partial D", "All D")
METRICS = (
    Metric("daa_chirality", "D-Cα chirality correctness", "chirality_daa_accuracy_pct", "%", True, "Folding"),
    BASE_METRICS["target"],
    BASE_METRICS["binder"],
    BASE_METRICS["lddt"],
    BASE_METRICS["placed"],
    BASE_METRICS["ilddt"],
    BASE_METRICS["dockq"],
)
BY_KEY = {metric.key: metric for metric in METRICS}


def read_groups(path: Path) -> dict[str, dict]:
    """Read either plot_daa.py chain_details.csv or the older research CSV."""
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No chain rows in {path}")
    sample_col = "sample_id" if "sample_id" in rows[0] else "sample"
    length_col = "chain_length" if "chain_length" in rows[0] else "length"
    required = {sample_col, length_col, "d_residues"}
    if not required.issubset(rows[0]):
        raise ValueError(f"{path}: requires {', '.join(sorted(required))}")
    by_sample = defaultdict(list)
    for row in rows:
        length, count = int(row[length_col]), int(row["d_residues"])
        if length <= 0 or not 1 <= count <= length:
            raise ValueError(f"Invalid D-residue count/chain length for {row[sample_col]}")
        label = "Single D" if count == 1 else "All D" if count == length else "Partial D"
        by_sample[row[sample_col]].append((label, length, count))
    groups = {}
    for sample, chains in by_sample.items():
        labels = {label for label, _, _ in chains}
        if len(labels) != 1:
            raise ValueError(f"{sample}: D-containing chains belong to different groups: {sorted(labels)}")
        length = sum(item[1] for item in chains)
        count = sum(item[2] for item in chains)
        groups[sample] = {"group": labels.pop(), "d_residues": count, "d_chain_count": len(chains),
                          "d_chain_residues": length, "d_fraction": count / length}
    return groups


def select_oracle(points: list[dict], key: str) -> dict:
    metric = BY_KEY[key]
    return min(points, key=lambda p: (
        -p[key] if metric.higher else p[key],
        -p["chirality"] if p["chirality"] is not None else 0,
        p["placed"] if p["placed"] is not None else float("inf"),
        p["seed"], p["prediction_path"],
    ))


def prepare(predictions: list[dict], groups: dict[str, dict], models: list[str],
            keys: list[str], expected_seeds: int | None) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    by_input = defaultdict(list)
    for point in predictions:
        by_input[(point["model"], point["sample_id"])].append(point)
    included, excluded = [], []
    for sample in sorted(groups):
        reasons = []
        for model in models:
            points = by_input.get((model, sample), [])
            if not points:
                reasons.append(f"{model}: no evaluated predictions")
            elif len({p["seed"] for p in points}) != len(points):
                reasons.append(f"{model}: duplicate seed IDs")
            elif expected_seeds is not None and len(points) != expected_seeds:
                reasons.append(f"{model}: {len(points)} evaluated seeds, expected {expected_seeds}")
            else:
                missing = [key for key in keys if any(p[key] is None for p in points)]
                if missing:
                    reasons.append(f"{model}: missing {','.join(missing)}")
        if reasons:
            excluded.append({"sample_id": sample, "group": groups[sample]["group"], "reason": "; ".join(reasons)})
        else:
            included.append(sample)
    if not included:
        raise ValueError("No common inputs with all requested metric values and seeds")
    seed_rows, input_rows = [], []
    for sample in included:
        for model in models:
            points = sorted(by_input[(model, sample)], key=lambda p: (p["seed"], p["prediction_path"]))
            for key in keys:
                metric = BY_KEY[key]
                for point in points:
                    seed_rows.append({"sample_id": sample, "group": groups[sample]["group"], "model": model,
                                      "seed": point["seed"], "metric": key, "value": point[key],
                                      "prediction_path": point["prediction_path"]})
                best = select_oracle(points, key)
                across_seeds = (pooled_daa_accuracy(points) if key == "daa_chirality"
                                else stats.mean(p[key] for p in points) if metric.higher
                                else stats.median(p[key] for p in points))
                input_rows.append({"sample_id": sample, **groups[sample], "model": model, "metric": key,
                                   "oracle_seed": best["seed"], "oracle_value": best[key],
                                   "oracle_prediction_path": best["prediction_path"], "n_seeds": len(points),
                                   "across_seed_value": across_seeds,
                                   "chirality_daa_correct_count": best.get("chirality_daa_correct_count") if key == "daa_chirality" else None,
                                   "chirality_daa_assessable_count": best.get("chirality_daa_assessable_count") if key == "daa_chirality" else None})
    group_rows = []
    for model in models:
        for key in keys:
            metric = BY_KEY[key]
            for group in GROUP_ORDER:
                selected = [row for row in input_rows if row["model"] == model and row["metric"] == key and row["group"] == group]
                values = [row["oracle_value"] for row in selected]
                if values:
                    statistic, center = (("mean", stats.mean(values)) if metric.higher
                                         else ("median", stats.median(values)))
                    if key == "daa_chirality":
                        statistic, center = "pooled assessable D sites", pooled_daa_accuracy(selected)
                    group_rows.append({"model": model, "metric": key, "group": group, "n_inputs": len(values),
                                       "statistic": statistic, "group_oracle_value": center})
    return excluded, seed_rows, input_rows, group_rows


def style_axis(ax, metric: Metric, values: list[float], groups: list[str], counts: dict[str, int]) -> None:
    ax.set_title(panel_title(metric.label, metric.higher), loc="left", fontsize=9, weight="bold")
    ax.set_ylabel(metric.label if metric.unit == "score" else f"{metric.label} ({metric.unit})")
    ax.set_xticks(range(len(groups)), [f"{g.replace(' ', chr(10))}\n(n={counts[g]})" for g in groups], fontsize=8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    if metric.unit == "%":
        ax.set_ylim(0, 105)
    elif metric.unit == "score":
        ax.set_ylim(0, 1.05)
    else:
        high = max(values) * 1.18 if values else 1
        ax.set_ylim(0, high)


def plot_bars(group_rows: list[dict], models: list[str], groups: list[str], keys: list[str],
              outdir: Path, formats: tuple[str, ...], discussion_title: bool = False,
              figure_id: str | None = None, version: str = "v01") -> None:
    fig, axes = panels(keys, "Metric-wise best seed per input, then D-content group summary", BY_KEY, discussion_title)
    lookup = {(r["model"], r["metric"], r["group"]): r for r in group_rows}
    counts = {group: lookup[(models[0], keys[0], group)]["n_inputs"] for group in groups}
    for key in keys:
        ax, metric = axes[key], BY_KEY[key]
        values = []
        for index, model in enumerate(models):
            heights = [lookup[(model, key, group)]["group_oracle_value"] for group in groups]
            values.extend(heights)
            positions = [i + (-0.18 if index == 0 else 0.18) for i in range(len(groups))]
            ax.bar(positions, heights, width=0.34, color=method_color(model, set(), models), zorder=3)
        style_axis(ax, metric, values, groups, counts)
    fig.legend(handles=[Patch(facecolor=method_color(model, set(), models), label=display_name(model)) for model in models],
               loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=len(models), frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, named_stem(outdir, "group_oracle_bars", figure_id, version), formats)


def plot_boxes(input_rows: list[dict], models: list[str], groups: list[str], keys: list[str],
               outdir: Path, formats: tuple[str, ...], discussion_title: bool = False,
               figure_id: str | None = None, version: str = "v01") -> None:
    fig, axes = panels(keys, "Metric-wise input oracles by D-content group · each dot is one input", BY_KEY, discussion_title)
    counts = {group: sum(r["group"] == group and r["model"] == models[0] and r["metric"] == keys[0] for r in input_rows) for group in groups}
    for key in keys:
        ax, metric = axes[key], BY_KEY[key]
        all_values = []
        for group_index, group in enumerate(groups):
            for model_index, model in enumerate(models):
                rows = sorted((r for r in input_rows if r["group"] == group and r["model"] == model and r["metric"] == key), key=lambda r: r["sample_id"])
                values = [r["oracle_value"] for r in rows]
                all_values.extend(values)
                position = group_index + (-0.18 if model_index == 0 else 0.18)
                color = method_color(model, set(), models)
                ax.boxplot(values, positions=[position], widths=0.27, patch_artist=True,
                           showfliers=False, manage_ticks=False,
                           boxprops={"facecolor": color, "alpha": 0.35, "edgecolor": color},
                           medianprops={"color": INK, "linewidth": 1.0},
                           whiskerprops={"color": color}, capprops={"color": color})
                jitter = [0 if len(values) == 1 else -0.085 + 0.17 * i / (len(values) - 1) for i in range(len(values))]
                ax.scatter([position + delta for delta in jitter], values, s=12.25, alpha=0.85,
                           color=color, marker=model_marker(model), edgecolors="white", linewidths=0.3, zorder=3)
        style_axis(ax, metric, all_values, groups, counts)
        ax.set_xlim(-0.55, len(groups) - 0.45)
    fig.legend(handles=[Patch(facecolor=method_color(model, set(), models), label=display_name(model)) for model in models],
               loc="lower center", bbox_to_anchor=(0.5, 0.01), ncol=len(models), frameon=False)
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    save(fig, named_stem(outdir, "group_oracle_boxplots", figure_id, version), formats)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chain-details", type=Path, required=True, help="chain_details.csv from plot_daa.py or legacy research composition")
    parser.add_argument("--source", action="append", required=True, help="MODEL=CSV or MODEL=CSV::PATH_METHOD; first source is non-steering, second is steering")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--metrics", help=f"Comma-separated subset of {','.join(BY_KEY)} (default all)")
    parser.add_argument("--stereochecks", choices=("on", "off"), default="on")
    parser.add_argument("--seed-regex", help="Regex with a capture group for seed ID")
    parser.add_argument("--expected-seeds", type=int, help="Require this many evaluated seeds per input and model")
    parser.add_argument("--formats", type=parse_formats, default=("svg",), help="Comma-separated svg,png,pdf (default: svg)")
    parser.add_argument("--figure-id", help="Final manuscript figure number, e.g. Fig05")
    parser.add_argument("--version", default="v01", help="Final manuscript version (default: v01)")
    parser.add_argument("--discussion-title", action="store_true", help="Add an 11 pt whole-figure title for standalone discussion figures")
    args = parser.parse_args()
    if len(args.source) != 2:
        parser.error("Exactly two --source arguments are required")
    if args.expected_seeds is not None and args.expected_seeds < 1:
        parser.error("--expected-seeds must be positive")
    keys = [part.strip() for part in args.metrics.split(",")] if args.metrics else list(BY_KEY)
    if not keys or len(keys) != len(set(keys)) or any(key not in BY_KEY for key in keys):
        parser.error(f"--metrics must contain unique keys from {','.join(BY_KEY)}")
    keys = ordered_panel_keys(keys, BY_KEY)
    setup_style()
    groups = read_groups(args.chain_details)
    predictions = load_predictions(args.source, args.seed_regex, None, stereochecks=args.stereochecks)
    models = ordered_models(list(dict.fromkeys(point["model"] for point in predictions)))
    excluded, seed_rows, input_rows, group_rows = prepare(predictions, groups, models, keys, args.expected_seeds)
    present_groups = [group for group in GROUP_ORDER if any(row["group"] == group for row in input_rows)]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "seed_metrics.csv", seed_rows)
    write_csv(args.output_dir / "input_metrics.csv", input_rows)
    write_csv(args.output_dir / "group_metrics.csv", group_rows)
    with (args.output_dir / "excluded_inputs.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "group", "reason"])
        writer.writeheader()
        writer.writerows(excluded)
    plot_bars(group_rows, models, present_groups, keys, args.output_dir, args.formats, args.discussion_title, args.figure_id, args.version)
    plot_boxes(input_rows, models, present_groups, keys, args.output_dir, args.formats, args.discussion_title, args.figure_id, args.version)
    counts = {group: sum(row["group"] == group and row["model"] == models[0] and row["metric"] == keys[0]
                         for row in input_rows) for group in present_groups}
    caption = (
        f"# Caption draft — complete and verify before manuscript use\n\n"
        f"Figure {args.figure_id[3:] if args.figure_id else 'X'}. [Main finding for D-content groups]. "
        f"Models: {display_name(models[0])} and {display_name(models[1])}. "
        f"Matched systems: {sum(counts.values())}; group sizes: "
        + ", ".join(f"{group} n={counts[group]}" for group in present_groups) + ". "
        f"Each dot is one system; predictions have {args.expected_seeds or 'variable'} evaluated seeds per system and model. "
        "The best seed is selected independently for each metric and system; these panels do not form a joint structure profile. "
        "Bars show system-level medians for RMSD, system-equal means for lDDT and DockQ, "
        "and residue-pooled D-Cα correctness. "
        "Boxes show median and interquartile range; whiskers reach the most extreme observations within 1.5×IQR; all system values are shown. "
        "D-Cα correctness uses passed/assessable D sites at the evaluation cutoff. "
        f"OpenStructure lDDT and iLDDT use stereochemistry checks {args.stereochecks}. "
        f"Excluded systems: {len(excluded)} (see excluded_inputs.csv). "
        + (L_RMSD_CAPTION if "placed" in keys else "")
        + "[Verify score provenance, chain mapping, lDDT implementation, and all other evaluation restrictions.]\n"
    )
    (args.output_dir / "caption_draft.md").write_text(caption, encoding="utf-8")
    print(f"Plotted {len(set(row['sample_id'] for row in input_rows))} common inputs; excluded {len(excluded)}. Figures: {', '.join(args.formats)}")


if __name__ == "__main__":
    main()

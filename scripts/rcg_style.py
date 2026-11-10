"""Shared RCG manuscript figure style (Notion rules, 2026-11-04)."""

from __future__ import annotations

import textwrap
import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager


SOURCE_URL = "RCG Paper Figure/Table Formatting (Notion, 2026-11-04)"
MM = 1 / 25.4
SINGLE_WIDTH = 90 * MM
DOUBLE_WIDTH = 180 * MM
INK = "#222222"
PAIR_LINE = "#B8B8B8"
REFERENCE = "#CCCCCC"
CONDITION_COLORS = {"baseline": "#007282", "steering": "#E69F00", "rcg": "#009E73", "rcg_steering": "#CC79A6"}
MODEL_COLORS = {
    "AF3": {"baseline": "#DACAB0", "steering": "#BCA987", "rcg": "#E6A45E", "rcg_steering": "#A96628"},
    "Boltz-1": {"baseline": "#C7E3EE", "steering": "#8AB2BE", "rcg": "#56B4E9", "rcg_steering": "#176A9C"},
    "Boltz-2": {"baseline": "#C6E5C6", "steering": "#8FB88F", "rcg": "#66C98B", "rcg_steering": "#237B4B"},
}
MODEL_MARKERS = {"AF3": "o", "Boltz-1": "s", "Boltz-2": "^"}
MODEL_ORDER = {"AF3": 0, "Boltz-1": 1, "Boltz-2": 2}
CONDITION_ORDER = {"baseline": 0, "steering": 1, "rcg": 2, "rcg_steering": 3}


def model_identity(label: str, steering_names: set[str] | None = None) -> tuple[str | None, str]:
    """Map common research aliases to the paper's base-model and condition names."""
    raw = label.lower().replace("_", "-").replace(" ", "")
    base = "AF3" if "af3" in raw or "alphafold3" in raw else (
        "Boltz-1" if "boltz-1" in raw or "boltz1" in raw else (
            "Boltz-2" if "boltz-2" in raw or "boltz2" in raw else None))
    rcg = "rcg" in raw
    negative = any(token in raw for token in ("no-steer", "non-steer", "nonsteer", "unsteer", "baseline"))
    steering = not negative and (label in (steering_names or set()) or "steer" in raw or raw.endswith(("1x", "2x", "+s", "-s")))
    condition = "rcg_steering" if rcg and steering else "rcg" if rcg else "steering" if steering else "baseline"
    return base, condition


def ordered_models(models: list[str], steering_names: set[str] | None = None) -> list[str]:
    return sorted(models, key=lambda name: (
        MODEL_ORDER.get(model_identity(name, steering_names)[0], 3),
        CONDITION_ORDER[model_identity(name, steering_names)[1]], name))


def model_color(label: str, all_models: list[str], steering_names: set[str] | None = None) -> str:
    base, condition = model_identity(label, steering_names)
    if base == "Boltz-2" and condition in {"rcg", "rcg_steering"}:
        raise ValueError("Boltz-2 + RCG is marked Not evaluated in the current RCG paper rules")
    bases = {model_identity(name, steering_names)[0] for name in all_models if model_identity(name, steering_names)[0]}
    if bases and any(model_identity(name, steering_names)[0] is None for name in all_models):
        raise ValueError("Every source label must identify its base model when an AF3 or Boltz source is present")
    if len(bases) > 1:
        if base is None:
            raise ValueError(f"Cannot assign a cross-model palette to {label!r}; use a model name containing AF3, Boltz-1, or Boltz-2")
        return MODEL_COLORS[base][condition]
    if not bases and len(all_models) > 1:
        conditions = [model_identity(name, steering_names)[1] for name in all_models]
        if len(set(conditions)) != len(conditions):
            raise ValueError("Model labels are ambiguous for the RCG palette; use recognizable AF3/Boltz names or --steering")
    return CONDITION_COLORS[condition]


def display_name(label: str, steering_names: set[str] | None = None) -> str:
    base, condition = model_identity(label, steering_names)
    if base is None:
        return label
    suffix = {"baseline": "", "steering": " + steering", "rcg": " + RCG", "rcg_steering": " + RCG + steering"}[condition]
    return base + suffix


def model_marker(label: str) -> str:
    return MODEL_MARKERS.get(model_identity(label)[0], "o")


def setup_style() -> None:
    try:
        font_manager.findfont("Arial", fallback_to_default=False)
    except ValueError as exc:
        raise RuntimeError("RCG figures require Arial at the final figure size; install a licensed Arial font on this machine") from exc
    plt.rcParams.update({
        "font.family": "Arial", "font.size": 9, "font.weight": "normal",
        "axes.titlesize": 9, "axes.titleweight": "bold", "axes.labelsize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
        "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 0.7,
        "lines.linewidth": 1.0, "lines.markersize": 3.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "figure.facecolor": "white", "savefig.facecolor": "white",
        "svg.fonttype": "none", "pdf.fonttype": 42,
    })


def label_metric(name: str, unit: str) -> str:
    return name if unit == "score" else f"{name} ({unit})"


def panel_title(name: str, higher: bool) -> str:
    return textwrap.fill(name, width=20, break_long_words=False) + f"\n({'↑' if higher else '↓'} better)"


def named_stem(outdir: Path, content: str, figure_id: str | None = None, version: str = "v01") -> Path:
    if figure_id is None:
        return outdir / content
    if not re.fullmatch(r"Fig\d{1,3}", figure_id):
        raise ValueError("--figure-id must look like Fig05")
    if not re.fullmatch(r"v\d{2}", version):
        raise ValueError("--version must look like v01")
    clean = re.sub(r"[^A-Za-z0-9_]+", "_", content).strip("_")
    return outdir / f"{figure_id}_{clean}_{version}"


def letter_panels(axes: list) -> None:
    if len(axes) < 2:
        return
    for i, ax in enumerate(axes):
        ax.text(-0.02, 1.02, chr(ord("a") + i), transform=ax.transAxes,
                fontsize=10, fontweight="bold", ha="right", va="bottom")


def save_figure(fig, stem: Path, formats: tuple[str, ...]) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        path = stem.with_suffix(f".{fmt}")
        fig.savefig(path, dpi=300, facecolor="white")
    plt.close(fig)

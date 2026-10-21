#!/usr/bin/env python3
"""Plot D-residue counts and fractions from Boltz YAML or AF3 JSON inputs."""

from __future__ import annotations

import argparse
import csv
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "protein_eval_graph_mpl"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml
from plot_metrics import parse_formats


D_CCD = frozenset("DAL DAR DAS DCY DGL DGN DHI DIL DLE DLY DPN DPR DSG DSN DTH DTR DTY DVA MED".split())
BLUE = "#2563A6"


def load_file(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle) if path.suffix.lower() == ".json" else yaml.safe_load(handle)


def sample_name(path: Path) -> str:
    stem = path.stem
    for suffix in ("_ccd_full", "_ccd", "_full"):
        if stem.endswith(suffix):
            return stem[: -len(suffix)]
    return stem


def find_inputs(folder: Path, pattern: str) -> list[Path]:
    if pattern == "auto":
        compact = sorted(folder.rglob("*_ccd.yaml")) + sorted(folder.rglob("*_ccd.yml"))
        if compact:
            return compact
        return sorted(p for p in folder.rglob("*") if p.is_file() and p.suffix.lower() in {".yaml", ".yml", ".json"} and not p.stem.endswith("_full"))
    return sorted(p for p in folder.rglob(pattern) if p.is_file())


def extract(folder: Path, pattern: str) -> tuple[list[dict], list[dict], list[str]]:
    files = find_inputs(folder, pattern)
    if not files:
        raise ValueError(f"No input files under {folder} matching {pattern}")
    samples, chains, without_d = [], [], []
    seen = set()
    for path in files:
        sample = sample_name(path)
        if sample in seen:
            raise ValueError(f"Duplicate sample ID {sample}; narrow --glob to one input version")
        seen.add(sample)
        data = load_file(path)
        if not isinstance(data, dict) or not isinstance(data.get("sequences"), list):
            raise ValueError(f"{path}: expected a sequences list")
        sample_chains = []
        for item in data["sequences"]:
            protein = item.get("protein") if isinstance(item, dict) else None
            if not protein:
                continue
            sequence = protein.get("sequence", "")
            if not isinstance(sequence, str) or not sequence:
                raise ValueError(f"{path}: protein sequence missing")
            mods = protein.get("modifications", [])
            d_positions = []
            for mod in mods:
                code = mod.get("ccd", mod.get("ptmType"))
                pos = mod.get("position", mod.get("ptmPosition"))
                if code not in D_CCD:
                    continue
                if not isinstance(pos, int) or not 1 <= pos <= len(sequence):
                    raise ValueError(f"{path}: invalid D-residue position {pos}")
                d_positions.append(pos)
            if len(d_positions) != len(set(d_positions)):
                raise ValueError(f"{path}: duplicate D-residue position")
            if not d_positions:
                continue
            ids = protein.get("id", [])
            ids = ids if isinstance(ids, list) else [ids]
            for chain_id in ids:
                if not chain_id:
                    raise ValueError(f"{path}: D-containing chain lacks an ID")
                row = {"sample_id": sample, "input_path": str(path), "chain_id": str(chain_id),
                       "chain_length": len(sequence), "d_residues": len(d_positions),
                       "d_positions": ";".join(map(str, sorted(d_positions)))}
                sample_chains.append(row)
                chains.append(row)
        if not sample_chains:
            without_d.append(sample)
            continue
        ids = [r["chain_id"] for r in sample_chains]
        if len(ids) != len(set(ids)):
            raise ValueError(f"{path}: duplicate chain ID")
        count = sum(r["d_residues"] for r in sample_chains)
        length = sum(r["chain_length"] for r in sample_chains)
        samples.append({"sample_id": sample, "input_path": str(path), "d_chain_count": len(sample_chains),
                        "d_residues": count, "d_chain_residues": length, "d_fraction": count / length,
                        "d_percent": 100 * count / length})
    if not samples:
        raise ValueError("No D-containing chains found")
    return samples, chains, without_d


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def plot(samples: list[dict], outdir: Path, formats: tuple[str, ...] = ("svg",)) -> None:
    ordered = sorted(samples, key=lambda r: (-r["d_fraction"], -r["d_residues"], r["sample_id"]))
    labels = [r["sample_id"] for r in ordered]
    ys = list(range(len(ordered)))
    fig, (left, right) = plt.subplots(1, 2, figsize=(14, max(5, len(ordered) * 0.38 + 1.5)), sharey=True)
    left.barh(ys, [r["d_residues"] for r in ordered], color=BLUE, height=0.72)
    right.barh(ys, [r["d_percent"] for r in ordered], color=BLUE, height=0.72)
    left.set_yticks(ys, labels)
    left.invert_yaxis()
    left.set_xlabel("Number of D residues")
    right.set_xlabel("D residues / residues in D-containing chains (%)")
    left.set_title("D-residue count", loc="left", weight="bold")
    right.set_title("D-residue fraction", loc="left", weight="bold")
    right.set_xlim(0, 105)
    for ax in (left, right):
        ax.grid(axis="x", color="#E2E8F0", zorder=0)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("D-amino-acid composition by input", weight="bold", y=1.015)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    for ext in formats:
        fig.savefig(outdir / f"daa_composition.{ext}", dpi=220, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--glob", default="auto", help="Recursive file pattern, or auto (prefer compact *_ccd.yaml files)")
    parser.add_argument("--formats", type=parse_formats, default=("svg",), help="Comma-separated output formats: svg,png,pdf (default: svg)")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    samples, chains, without_d = extract(args.inputs, args.glob)
    write_csv(args.output_dir / "sample_summary.csv", samples, list(samples[0]))
    write_csv(args.output_dir / "chain_details.csv", chains, list(chains[0]))
    write_csv(args.output_dir / "inputs_without_d_residues.csv", [{"sample_id": s} for s in without_d], ["sample_id"])
    plot(samples, args.output_dir, args.formats)
    print(f"Plotted {len(samples)} inputs; {len(without_d)} inputs had no D-containing chain")


if __name__ == "__main__":
    main()

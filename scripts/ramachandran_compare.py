#!/usr/bin/env python3
"""
ramachandran_compare.py

Compute backbone phi/psi angles for one native structure and any number of
predictions, then write CSV files and Ramachandran plots.

Designed to work with PDB or mmCIF files and to tolerate non-standard / D-amino
acid residue names, as long as each residue has backbone atoms N, CA, and C.

Dependencies:
    pip install biopython numpy pandas matplotlib

Example:
    python ramachandran_compare.py \
        --native 7YH8.cif B Native \
        --pred pred_boltz1.cif B Boltz1 \
        --pred pred_boltz2.cif B Boltz2 \
        --pred pred_af3.cif B AF3 \
        --outdir rama_7YH8 \
        --annotate
"""

from __future__ import annotations

import argparse
import math
import os
import re
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "protein_eval_graph_mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser, PDBParser


# ---------------------------------------------------------------------------
# Structure loading
# ---------------------------------------------------------------------------

def load_structure(path: Path):
    suffix = path.suffix.lower()

    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    elif suffix in {".pdb", ".ent"}:
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(
            f"Unsupported structure format: {path}\n"
            "Use .pdb, .ent, .cif, or .mmcif."
        )

    return parser.get_structure(path.stem, str(path))


def get_chain(structure, chain_id: str):
    """
    Return chain from model 0.

    Tries exact Bio.PDB chain ID first. If not found, prints available IDs.
    """
    model = structure[0]

    if chain_id in model:
        return model[chain_id]

    available = [chain.id for chain in model]
    raise KeyError(
        f"Chain '{chain_id}' not found. Available chain IDs: {available}"
    )


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------

def dihedral_deg(p0, p1, p2, p3) -> float:
    """
    Signed dihedral angle in degrees, range (-180, 180].
    """
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    p3 = np.asarray(p3, dtype=float)

    b0 = -(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    norm_b1 = np.linalg.norm(b1)
    if norm_b1 == 0:
        return np.nan

    b1 = b1 / norm_b1

    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1

    norm_v = np.linalg.norm(v)
    norm_w = np.linalg.norm(w)
    if norm_v == 0 or norm_w == 0:
        return np.nan

    v = v / norm_v
    w = w / norm_w

    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)

    return float(np.degrees(np.arctan2(y, x)))


def distance(a, b) -> float:
    return float(np.linalg.norm(np.asarray(a) - np.asarray(b)))


def circular_angle_difference_deg(a: float, b: float) -> float:
    """
    Smallest signed difference a-b for angles in degrees, in [-180, 180).
    """
    if pd.isna(a) or pd.isna(b):
        return np.nan
    return ((a - b + 180.0) % 360.0) - 180.0


# ---------------------------------------------------------------------------
# Residue extraction
# ---------------------------------------------------------------------------

def residue_display_name(residue) -> str:
    """
    Example: DAL1, DLY14, GLY5A
    """
    _, resseq, icode = residue.id
    icode = icode.strip()
    return f"{residue.resname.strip()}{resseq}{icode}"


def residue_key(residue) -> Tuple[int, str]:
    """
    Position key used for native-vs-prediction matching.
    """
    _, resseq, icode = residue.id
    return int(resseq), icode.strip()


def backbone_residues(chain):
    """
    Keep every residue that has N, CA, and C, regardless of residue name.

    This deliberately does NOT call is_aa(...), because D-amino acids and
    modified residues may be encoded with non-standard residue names.
    """
    residues = []
    for residue in chain.get_residues():
        if all(atom_name in residue for atom_name in ("N", "CA", "C")):
            residues.append(residue)
    return residues


def classify_region(phi: float, psi: float) -> str:
    """
    Deliberately broad, heuristic labels for quick inspection only.
    These are NOT a substitute for a formal favored/allowed-region definition.

    Useful here mainly for spotting alpha_L-like vs alpha_R-like basin switches.
    """
    if pd.isna(phi) or pd.isna(psi):
        return "undefined"

    # Broad alpha-like windows
    if -100 <= phi <= -30 and -90 <= psi <= 20:
        return "alpha_R_like"

    if 30 <= phi <= 100 and -20 <= psi <= 90:
        return "alpha_L_like"

    # Broad extended / beta-like area
    if -180 <= phi <= -60 and 60 <= psi <= 180:
        return "beta_like"

    return "other"


def compute_phi_psi(
    structure_path: Path,
    chain_id: str,
    label: str,
    peptide_bond_cutoff: float = 2.0,
) -> pd.DataFrame:
    """
    Compute phi/psi for all backbone-containing residues in a chain.

    A phi or psi angle is skipped across an apparent chain break. The chain-break
    test uses the peptide C-N distance with a default cutoff of 2.0 Å.
    """
    structure = load_structure(structure_path)
    chain = get_chain(structure, chain_id)
    residues = backbone_residues(chain)

    if len(residues) < 2:
        raise ValueError(
            f"{label}: chain {chain_id} contains fewer than 2 residues "
            "with N/CA/C backbone atoms."
        )

    rows = []

    for i, residue in enumerate(residues):
        phi = np.nan
        psi = np.nan

        current_n = residue["N"].coord
        current_ca = residue["CA"].coord
        current_c = residue["C"].coord

        # phi_i = C_(i-1) - N_i - CA_i - C_i
        if i > 0:
            prev = residues[i - 1]
            prev_c = prev["C"].coord

            if distance(prev_c, current_n) <= peptide_bond_cutoff:
                phi = dihedral_deg(prev_c, current_n, current_ca, current_c)

        # psi_i = N_i - CA_i - C_i - N_(i+1)
        if i < len(residues) - 1:
            nxt = residues[i + 1]
            next_n = nxt["N"].coord

            if distance(current_c, next_n) <= peptide_bond_cutoff:
                psi = dihedral_deg(current_n, current_ca, current_c, next_n)

        _, resseq, icode = residue.id

        rows.append(
            {
                "structure": label,
                "file": str(structure_path),
                "chain": chain_id,
                "residue_index_in_chain": i + 1,
                "resseq": int(resseq),
                "icode": icode.strip(),
                "resname": residue.resname.strip(),
                "residue": residue_display_name(residue),
                "phi_deg": phi,
                "psi_deg": psi,
                "region": classify_region(phi, psi),
            }
        )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def setup_ramachandran_axes(ax, title: str):
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks([-180, -120, -60, 0, 60, 120, 180])
    ax.set_yticks([-180, -120, -60, 0, 60, 120, 180])
    ax.axhline(0, linewidth=0.7, alpha=0.35)
    ax.axvline(0, linewidth=0.7, alpha=0.35)
    ax.grid(alpha=0.18)

    ax.set_xlabel(r"$\phi$ (degrees)")
    ax.set_ylabel(r"$\psi$ (degrees)")
    ax.set_title(title)

    # Approximate reference centers, intentionally labeled as "approx."
    ax.text(
        -60, -45, r"$\alpha_R$ approx.",
        ha="center", va="center", fontsize=9, alpha=0.45
    )
    ax.text(
        60, 45, r"$\alpha_L$ approx.",
        ha="center", va="center", fontsize=9, alpha=0.45
    )
    ax.text(
        -125, 130, r"$\beta$ approx.",
        ha="center", va="center", fontsize=9, alpha=0.45
    )

    ax.set_aspect("equal", adjustable="box")


def sanitize_filename(text: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text.strip())
    return text.strip("_") or "structure"


def plot_single(df: pd.DataFrame, out_path: Path, annotate: bool = False):
    valid = df.dropna(subset=["phi_deg", "psi_deg"])

    fig, ax = plt.subplots(figsize=(7, 7))
    setup_ramachandran_axes(
        ax,
        f"{df['structure'].iloc[0]} — chain {df['chain'].iloc[0]}"
    )

    ax.scatter(
        valid["phi_deg"],
        valid["psi_deg"],
        s=48,
        alpha=0.85,
    )

    if annotate:
        for _, row in valid.iterrows():
            ax.annotate(
                row["residue"],
                (row["phi_deg"], row["psi_deg"]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7,
                alpha=0.85,
            )

    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def plot_overlay(
    datasets: List[pd.DataFrame],
    out_path: Path,
    annotate_native: bool = False,
):
    fig, ax = plt.subplots(figsize=(8, 8))
    setup_ramachandran_axes(ax, "Native vs predictions — Ramachandran overlay")

    markers = ["o", "s", "^", "D", "P", "X", "v", "<", ">"]

    for i, df in enumerate(datasets):
        valid = df.dropna(subset=["phi_deg", "psi_deg"])
        marker = markers[i % len(markers)]

        label = str(df["structure"].iloc[0])
        lower = label.lower()
        is_steering = ("steer" in lower and not any(x in lower for x in ("no_steer", "nonsteer", "unsteer"))) or lower.endswith(("1x", "2x"))
        color = "#64748B" if i == 0 else ("#E87528" if is_steering else "#2563A6")
        ax.scatter(
            valid["phi_deg"],
            valid["psi_deg"],
            s=52 if i == 0 else 40,
            marker=marker,
            color=color,
            alpha=0.8,
            label=label,
        )

        # Usually only annotate native to avoid an unreadable plot.
        if i == 0 and annotate_native:
            for _, row in valid.iterrows():
                ax.annotate(
                    row["residue"],
                    (row["phi_deg"], row["psi_deg"]),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=7,
                    alpha=0.8,
                )

    ax.legend(frameon=True, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Native-vs-prediction comparison
# ---------------------------------------------------------------------------

def make_comparison_table(
    native_df: pd.DataFrame,
    pred_dfs: List[pd.DataFrame],
) -> pd.DataFrame:
    """
    Match residues primarily by (resseq, insertion code).

    Note:
      If a prediction renumbers the peptide differently from the native, this
      positional matching may fail. In that case use the per-structure CSVs or
      renumber chains consistently before comparison.
    """
    native = native_df.copy()
    native["key"] = list(zip(native["resseq"], native["icode"]))

    native_lookup = native.set_index("key")

    rows = []

    for pred in pred_dfs:
        pred_label = pred["structure"].iloc[0]

        for _, prow in pred.iterrows():
            key = (prow["resseq"], prow["icode"])

            if key not in native_lookup.index:
                continue

            nrow = native_lookup.loc[[key]].iloc[0]

            # Protect against duplicate residue IDs, though uncommon here.
            if isinstance(nrow, pd.DataFrame):
                nrow = nrow.iloc[0]

            dphi = circular_angle_difference_deg(
                prow["phi_deg"], nrow["phi_deg"]
            )
            dpsi = circular_angle_difference_deg(
                prow["psi_deg"], nrow["psi_deg"]
            )

            rows.append(
                {
                    "prediction": pred_label,
                    "resseq": prow["resseq"],
                    "icode": prow["icode"],
                    "native_resname": nrow["resname"],
                    "prediction_resname": prow["resname"],
                    "native_phi_deg": nrow["phi_deg"],
                    "native_psi_deg": nrow["psi_deg"],
                    "prediction_phi_deg": prow["phi_deg"],
                    "prediction_psi_deg": prow["psi_deg"],
                    "delta_phi_deg": dphi,
                    "delta_psi_deg": dpsi,
                    "abs_delta_phi_deg": abs(dphi) if not pd.isna(dphi) else np.nan,
                    "abs_delta_psi_deg": abs(dpsi) if not pd.isna(dpsi) else np.nan,
                    "native_region": nrow["region"],
                    "prediction_region": prow["region"],
                    "same_region": (
                        nrow["region"] == prow["region"]
                        if nrow["region"] != "undefined"
                        and prow["region"] != "undefined"
                        else np.nan
                    ),
                }
            )

    return pd.DataFrame(rows)


def print_summary(native_df, pred_dfs, comparison_df):
    print("\n=== Ramachandran summary ===")

    for df in [native_df] + pred_dfs:
        label = df["structure"].iloc[0]
        chain = df["chain"].iloc[0]
        counts = df["region"].value_counts().to_dict()
        valid_n = int(df[["phi_deg", "psi_deg"]].dropna().shape[0])

        print(
            f"{label:20s} chain={chain:>4s}  "
            f"valid(phi,psi)={valid_n:2d}  regions={counts}"
        )

    if not comparison_df.empty:
        print("\nMean absolute torsion differences vs native:")
        grouped = comparison_df.groupby("prediction")[
            ["abs_delta_phi_deg", "abs_delta_psi_deg"]
        ].mean()
        print(grouped.round(2).to_string())


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute and compare Ramachandran phi/psi angles for a native "
            "structure and multiple predictions."
        )
    )

    parser.add_argument(
        "--native",
        nargs=3,
        metavar=("FILE", "CHAIN", "LABEL"),
        required=True,
        help="Native structure specification, e.g. --native 7YH8.cif B Native",
    )

    parser.add_argument(
        "--pred",
        nargs=3,
        metavar=("FILE", "CHAIN", "LABEL"),
        action="append",
        required=True,
        default=[],
        help=(
            "Prediction specification. Repeat for multiple predictions, e.g. "
            "--pred pred1.cif B Boltz1 --pred pred2.cif B Boltz2"
        ),
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("ramachandran_output"),
        help="Output directory (default: ramachandran_output)",
    )

    parser.add_argument(
        "--annotate",
        action="store_true",
        help="Annotate residues on individual plots and native points on overlay.",
    )

    parser.add_argument(
        "--peptide-bond-cutoff",
        type=float,
        default=2.0,
        help=(
            "Maximum C(i)-N(i+1) distance in Å to treat residues as connected "
            "(default: 2.0)."
        ),
    )

    return parser.parse_args()


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    native_file, native_chain, native_label = args.native
    native_df = compute_phi_psi(
        Path(native_file),
        native_chain,
        native_label,
        peptide_bond_cutoff=args.peptide_bond_cutoff,
    )

    pred_dfs = []
    for pred_file, pred_chain, pred_label in args.pred:
        pred_df = compute_phi_psi(
            Path(pred_file),
            pred_chain,
            pred_label,
            peptide_bond_cutoff=args.peptide_bond_cutoff,
        )
        pred_dfs.append(pred_df)

    # Save per-structure CSV + plot
    all_dfs = [native_df] + pred_dfs
    for df in all_dfs:
        label = df["structure"].iloc[0]
        safe = sanitize_filename(label)

        csv_path = args.outdir / f"{safe}_phi_psi.csv"
        plot_path = args.outdir / f"{safe}_ramachandran.png"

        df.to_csv(csv_path, index=False)
        plot_single(df, plot_path, annotate=args.annotate)

    # Combined long-form CSV
    combined_df = pd.concat(all_dfs, ignore_index=True)
    combined_df.to_csv(
        args.outdir / "all_structures_phi_psi.csv",
        index=False,
    )

    # Overlay
    plot_overlay(
        all_dfs,
        args.outdir / "ramachandran_overlay.png",
        annotate_native=args.annotate,
    )

    # Native-vs-prediction residue table
    comparison_df = make_comparison_table(native_df, pred_dfs)
    comparison_df.to_csv(
        args.outdir / "native_vs_predictions_torsions.csv",
        index=False,
    )

    print_summary(native_df, pred_dfs, comparison_df)

    print(f"\nOutputs written to: {args.outdir.resolve()}")
    print("  - one *_phi_psi.csv per structure")
    print("  - one *_ramachandran.png/.svg per structure")
    print("  - all_structures_phi_psi.csv")
    print("  - native_vs_predictions_torsions.csv")
    print("  - ramachandran_overlay.png")
    print(
        "\nImportant: alpha_R_like / alpha_L_like / beta_like are broad "
        "heuristic basin labels for inspection, not formal validation regions."
    )


if __name__ == "__main__":
    main()

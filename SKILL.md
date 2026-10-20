---
name: protein-eval-graph
description: Generate D-amino-acid composition, folding/docking comparison, seed-selection, chirality scatter, seed-range, metric-distribution, and Ramachandran plots from protein evaluation data. Use for plotting existing protein prediction evaluations, not for calculating structure scores.
---

# Protein evaluation graphs

Use the bundled scripts to plot existing evaluation results. Read [references/graph-types.md](references/graph-types.md) for the selected graph's inputs, aggregation, selection rule, and example command. The repository [README.md](README.md) documents installation and all graph types.

## Workflow

1. Identify the requested graph and inspect the relevant inputs: `batch_summary.csv`, Boltz/AF3 input YAML or JSON, or native/predicted CIF/PDB structures. Confirm model labels, seed identifiers, sample IDs, and the applicable chain mapping. Do not infer biological chain roles from file names.
2. Run the matching `scripts/plot_metrics.py` subcommand, `scripts/plot_daa.py`, or `scripts/ramachandran_compare.py`. Use a separate output directory. For a mixed-method batch CSV, use each `--source NAME=CSV::PATH_METHOD` filter; for a one-method CSV, `--source NAME=CSV` suffices.
3. Check the written plot-data CSV and image files. Report exclusions, unavailable metrics, selection rule, aggregation, and any unmatched pairs. Never silently replace a requested metric or treat different seeds as independent biological inputs.

## Shared formatting and meaning

- Folding: **Target self-aligned backbone RMSD**, **Binder self-aligned backbone RMSD**, **Global all-atom lDDT**, **Chirality accuracy**.
- Docking: **Interface lDDT**, **DockQ**, **Target-aligned binder backbone RMSD**. **Interface RMSD** is the eighth comparison panel when available.
- Put `(↑ better)` or `(↓ better)` beside metric panel titles. Use blue for non-steering and orange for steering. Keep a metric's axis scale consistent across methods in a figure.
- RMSD uses Å; lDDT and DockQ use 0–1; chirality accuracy uses percent. RMSD values are lower-is-better. Other listed metrics are higher-is-better.
- Aggregate binder-group RMSDs to one median per prediction before comparing predictions. For Interface RMSD, use the median across interfaces. Across selected inputs, show medians for RMSDs, arithmetic means for lDDT and DockQ, and residue-pooled chirality accuracy when residue counts are available.
- Top-1 means highest supplied confidence score. A metric-specific oracle means the best observed seed for *that metric*, selected separately for each input and method; it is descriptive and not an achievable multi-metric prediction. If confidence scores are absent, request them or use a user-specified selection instead of inventing top-1.
- Pair two methods by input ID for one selected seed per input; for all-seed scatter plots, connect only equal input ID **and** equal seed ID. Do not draw arbitrary connections.
- For Ramachandran overlays, use explicitly selected native and predicted chain IDs. Sequence/residue numbering mismatches must be noted when interpreting per-residue differences.

The plotting scripts do not score structures. Use `protein-eval-skill` first when evaluation CSVs do not yet exist.

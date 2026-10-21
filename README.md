# protein-eval-graph

A reusable Codex skill and command-line plotting tools for protein prediction evaluations. It draws figures from existing `protein-eval-skill` `batch_summary.csv` files, D-amino-acid input files, and native/predicted structures. The scripts generate plots and plot-data CSVs; they do not evaluate structures.

## Install

Install the Codex skill from this standalone repository:

```bash
python3 ~/.codex/skills/.system/skill-installer/scripts/install-skill-from-github.py \
  --repo hsicat/protein-eval-graph --path . --name protein-eval-graph
```

For direct script use, clone the repository and install its Python dependencies:

```bash
git clone https://github.com/hsicat/protein-eval-graph.git
cd protein-eval-graph
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

Use `$protein-eval-graph` in Codex after installation. Example request:

> Use $protein-eval-graph to draw paired AF3 versus AF3 steering chirality scatter plots from these two `batch_summary.csv` files, using the metric-specific oracle seed per input. Save SVG and plot-data CSV files in `figures/`.

## Graph types

| Type | Description | Tool |
|---|---|---|
| DAA composition bars | For each input, compare D-residue count and D-residue fraction across its D-containing protein chains. | [plot_daa.py](scripts/plot_daa.py) |
| DAA-content group oracle bars and boxplots | Group inputs by the D proportion of their D-containing chains (Single D, Partial D, All D); compare metric-wise best seeds across two methods using group summary bars and per-input boxplots with dots. | [plot_daa_groups.py](scripts/plot_daa_groups.py) |
| Folding and docking comparison bars | For one input ID, compare model predictions across four folding and four docking metrics. | [plot_metrics.py](scripts/plot_metrics.py) `compare` |
| Seed-selection comparison bars | Compare confidence top-1 and oracle selections across models using median/mean bars, optionally with sample-SD error bars. | [plot_metrics.py](scripts/plot_metrics.py) `selection-bars` |
| Chirality versus structure scatter | Plot structure metrics against pooled D/L chirality accuracy using all seeds, top-1, or metric-specific oracle. Two methods can be paired with dotted lines. | [plot_metrics.py](scripts/plot_metrics.py) `scatter` |
| Seed min–max ranges | Show each seed score and the observed min–max range for each input and method, one graph per metric. | [plot_metrics.py](scripts/plot_metrics.py) `seed-range` |
| Metric distributions | Compare a many-seed unsteered violin/box/points with a few individual steering points for one input. | [plot_metrics.py](scripts/plot_metrics.py) `distribution` |
| Ramachandran overlay | Compare native and predicted chain backbone φ/ψ angles; export single plots, an overlay, and residue-level CSVs. | [ramachandran_compare.py](scripts/ramachandran_compare.py) |

All four tool scripts are in the repository's [scripts directory](scripts/). `plot_metrics.py` handles five graph types through the subcommands shown above.

The exact metric labels are **Target self-aligned backbone RMSD**, **Binder self-aligned backbone RMSD**, **Global all-atom lDDT**, **Chirality accuracy**, **Interface lDDT**, **DockQ**, and **Target-aligned binder backbone RMSD**. The older eight-panel comparison also contains **Interface RMSD**, so the tool includes it as the eighth metric when data are available. DAA-group plots specifically use **D-residue chirality accuracy**, while the other plots use pooled D/L chirality accuracy. Every metric panel states `(↑ better)` or `(↓ better)`. Multi-metric figures use two rows: folding metrics above docking metrics. Blue denotes non-steering and orange denotes steering; use `--steering MODEL` when the label alone does not identify steering. For DAA-group plots, pass the non-steering source first and steering source second. RMSDs use Å, chirality uses %, and lDDT/DockQ use 0–1 scores.

## Input conventions

`plot_metrics.py` reads evaluation CSVs with at least `sample_id`, `prediction_path`, and `status`, plus requested metric columns. It accepts the `batch_summary.csv` written by `protein-eval-skill`. Rows for multiple binder groups of one prediction are collapsed: binder RMSDs become their median; Interface RMSD is the median across entries in `ost_per_interface_json`. Skipped predictions are excluded. Prediction-level fields must agree across binder-group rows.

Pass one `--source` per model:

```text
--source MODEL=/path/to/batch_summary.csv
--source MODEL=/path/to/batch_summary.csv::METHOD_TOKEN
```

The `::METHOD_TOKEN` form filters a shared CSV by the component immediately after `predicted/` in `prediction_path`. **Use that filter whenever one CSV contains multiple models.** For separate one-model CSVs, omit it. For unusual file layouts, `--seed-regex 'pattern(with_capture_group)'` overrides seed extraction. For top-1 selection, scores must be present in a `confidence_score` CSV column, in `--confidence-csv` (`prediction_path,confidence_score`), or in Boltz `confidence_<structure>.json` files beside predicted structures under `--confidence-root`.

Default metric keys, in panel order: `target,binder,lddt,chirality,ilddt,dockq,placed,irmsd`. Use `--metrics target,binder,...` to choose a subset or omit unavailable iRMSD. Scatter plots omit chirality as a Y metric by default. **Every command saves SVG only by default**; use `--formats svg,png,pdf` to request additional image formats, or `--formats png` / `--formats pdf` for only that format. Data CSVs are still written. `plot_metrics.py` also writes `normalized_predictions.csv`.

## Examples

These commands assume the repository is the working directory. Replace the paths with your own evaluation files.

### 1. DAA composition

```bash
python scripts/plot_daa.py --inputs /data/inputs --output-dir figures/daa
```

The automatic file search prefers compact `*_ccd.yaml` Boltz inputs and otherwise reads JSON/YAML while excluding `_full` variants. Use `--glob` to narrow a mixed folder to one intended input version. The denominator includes full lengths of only the D-containing protein chains. Inputs without one are listed separately.

### 2. DAA-content group oracle bars and boxplots

```bash
python scripts/plot_daa_groups.py \
  --chain-details figures/daa/chain_details.csv \
  --source AF3=/data/af3_baseline/evaluation/results/batch_summary.csv \
  --source 'AF3+S=/data/af3_steering/evaluation/results/batch_summary.csv' \
  --expected-seeds 5 --output-dir figures/daa_groups
```

`chain_details.csv` comes from the DAA composition command above; the script also accepts the older research CSV with `sample`, `length`, and `d_residues`. Each D-containing chain must belong to the same group within an input: **Single D** means exactly one D residue; **Partial D** means more than one but fewer than the chain length; **All D** means every residue is D. Only inputs with evaluated predictions and complete requested metrics in both models are compared. `--expected-seeds 5` also requires five evaluated seeds on each side. Each metric selects its own best seed per input and model. Bar heights are medians across inputs for RMSDs and means for bounded scores and D-residue chirality; the boxplots show all input-level oracle values, one dot per input. Outputs are `group_oracle_bars.svg`, `group_oracle_boxplots.svg`, `seed_metrics.csv`, `input_metrics.csv`, `group_metrics.csv`, and `excluded_inputs.csv` by default. The first source is colored blue, the second orange.

### 3. One-input folding/docking comparison

```bash
python scripts/plot_metrics.py compare --input-id 6Q77 \
  --source 'Boltz1=/data/evaluation/results/batch_summary.csv::boltz1_no_steering' \
  --source 'Boltz1_steering=/data/evaluation/results/batch_summary.csv::boltz1_steering' \
  --source 'Boltz2=/data/evaluation/results/batch_summary.csv::boltz2_no_steering' \
  --source 'Boltz2_steering=/data/evaluation/results/batch_summary.csv::boltz2_steering' \
  --output-dir figures/6Q77
```

If a source contains multiple seeds for 6Q77, specify `--selection top1` (with confidence scores) or `--selection oracle` (best seed per metric).

### 4. Seed-selection bars

```bash
python scripts/plot_metrics.py selection-bars \
  --source 'Boltz1=/data/evaluation/results/batch_summary.csv::boltz1' \
  --source 'Boltz1x=/data/evaluation/results/batch_summary.csv::boltz1x' \
  --confidence-root /data/evaluation --scope common --error-bars \
  --output-dir figures/selection
```

Top-1 maximizes confidence score. Oracle selects the best seed separately for each metric; `--oracle-target chirality` uses one chirality-oracle seed across all panels. RMSD bar heights are medians across inputs; bounded scores are means; chirality is pooled over assessable residues. Each panel's Y axis names the corresponding aggregation. Error bars are sample SD across selected input-level scores, not confidence intervals.

### 5. Chirality scatter

```bash
python scripts/plot_metrics.py scatter --selection oracle \
  --source AF3=/data/af3_baseline/evaluation/results/batch_summary.csv \
  --source AF3_steering=/data/af3_steering/evaluation/results/batch_summary.csv \
  --output-dir figures/chirality
```

Use `--selection all` for every seed or `--selection top1` for confidence-selected seeds. With two models, selected points for the same input are connected. In all-seed mode, only equal input and seed IDs are connected. The X coordinate of an oracle point is the chirality accuracy of the prediction chosen for that Y metric.

### 6. Seed min–max ranges

```bash
python scripts/plot_metrics.py seed-range \
  --source AF3=/data/af3_baseline/evaluation/results/batch_summary.csv \
  --source AF3_steering=/data/af3_steering/evaluation/results/batch_summary.csv \
  --output-dir figures/seed_ranges
```

Each graph covers one metric. Input IDs are ordered alphabetically, with both methods on neighboring rows and one point per evaluated seed.

### 7. Metric distributions

```bash
python scripts/plot_metrics.py distribution --input-id 7YH8 \
  --source Boltz1=/data/boltz1_50seeds/batch_summary.csv \
  --source Boltz1x=/data/boltz1x_steering/batch_summary.csv \
  --output-dir figures/7YH8_distributions
```

Put the unsteered model first. The script records the actual numbers of baseline and steering predictions; the usual case is 50 versus 3–5.

### 8. Ramachandran plot

```bash
python scripts/ramachandran_compare.py \
  --native /data/GT/7YH8.cif B Native \
  --pred /data/predictions/7YH8.cif B Steering \
  --outdir figures/rama_7YH8 --annotate
```

Pass the confirmed chain ID for each structure. The script exports the native and prediction φ/ψ points, overlay, and residue comparison. A residue-level comparison assumes matching numbering and insertion codes; inspect the CSV if the structures were renumbered.

## Reproducibility

Keep the source CSVs, selected `--source` filters, confidence scores, script command, and exported plot-data CSV with each figure. Do not infer a model's seed ranking from evaluation metrics. The graph scripts do not change the source data. Their selection and aggregation rules are documented in [references/graph-types.md](references/graph-types.md).

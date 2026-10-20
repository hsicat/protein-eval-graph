# Plot selection and data rules

Read the matching section when a graph is requested. The [README](../README.md) has install and CLI examples.

## DAA composition

Use `plot_daa.py` on Boltz YAML or AF3 JSON inputs. Count a residue only when its modification CCD/PTM code is one of the 19 D-amino-acid codes in the script. Count every physical chain ID separately. For each input, sum D residues and full lengths of **D-containing protein chains only**. Skip inputs with none and report them in `inputs_without_d_residues.csv`. The `auto` file pattern prefers compact `*_ccd.yaml`; narrow `--glob` when a directory has several versions of one input.

## Folding/docking comparison

Use `plot_metrics.py compare --input-id ID` with one or more model sources. Eight panels are available: the four folding metrics, three requested docking metrics, and Interface RMSD as the eighth panel from the earlier comparison figures. By default each model must contribute exactly one prediction for that ID. If multiple seeds exist, explicitly choose `--selection top1` or `--selection oracle`; the latter selects a different seed per panel. If iRMSD is unavailable, omit it using `--metrics` and report that it was unavailable.

## Seed-selection bars

Use `selection-bars` for top-1 versus metric-specific oracle. Top-1 requires a `confidence_score` column, `--confidence-csv`, or adjacent Boltz confidence JSON files via `--confidence-root`. `--scope common` retains IDs with both methods; `--scope all` uses available IDs per method. RMSDs are medians across selected inputs; lDDT, iLDDT, and DockQ are arithmetic means; chirality is pooled over assessable residues when counts are present. `--error-bars` shows sample SD of selected input-level values, **not** uncertainty of the aggregate. `--oracle-target chirality` reproduces the older single-chirality-oracle selection on every panel.
Oracle ties prefer higher chirality accuracy, then lower target-aligned binder backbone RMSD, then the lowest seed ID.

## Chirality-versus-structure scatter

Use `scatter --selection all|top1|oracle`. The default six structure metrics plus Interface RMSD are plotted against pooled D/L chirality accuracy; chirality-versus-itself is omitted by default. For `oracle`, select the best seed separately for each Y metric, model, and input; its actual chirality accuracy supplies X. With two methods, connect same input IDs for selected points. For `all`, connect only points with the same input and seed IDs. Correlations printed on panels are descriptive; seeds from one input are repeated observations.
RMSD Y axes use a symmetric log scale, linear up to 2 Å, to show small errors and large outliers together.

## Seed min–max

Use `seed-range` with two model sources, normally non-steering and steering. One figure is written per metric. Each row is an alphabetically ordered input ID, with separate positions for each method. A point represents each seed; a horizontal line spans the smallest to largest observed value. Missing values are omitted from that metric and remain visible in `normalized_predictions.csv`.

## Metric distributions

Use `distribution --input-id ID` with the unsteered source first and steered source second. A violin, box, and individual blue points show the unsteered seeds; orange diamonds show each steered result. The typical comparison is 50 versus 3–5 seeds, but the script records the actual counts. Both groups must have values for every requested metric.

## Ramachandran overlay

Use `ramachandran_compare.py --native FILE CHAIN LABEL --pred FILE CHAIN LABEL`. Native and predicted chain IDs must be explicitly supplied from known chain mappings. The script computes φ/ψ from backbone atoms N, CA, C, omits absent or disconnected torsions, exports per-structure CSVs, PNG/SVG single plots and an overlay, and a residue comparison CSV. Per-residue comparison uses residue numbers and insertion codes; verify numbering before interpreting differences. Read labels as structure identities, not quality classes.

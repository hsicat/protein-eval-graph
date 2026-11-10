# Plot selection and data rules

Read the matching section when a graph is requested and apply the [RCG formatting and data rules](rcg-formatting.md). The [README](../README.md) has install and CLI examples.

## DAA composition

Use `plot_daa.py` on Boltz YAML or AF3 JSON inputs. Count a residue only when its modification CCD/PTM code is one of the 19 D-amino-acid codes in the script. Count every physical chain ID separately. For each input, sum D residues and full lengths of **D-containing protein chains only**. Skip inputs with none and report them in `inputs_without_d_residues.csv`. The `auto` file pattern prefers compact `*_ccd.yaml`; narrow `--glob` when a directory has several versions of one input.

## DAA-content group oracle bars and boxplots

Use `plot_daa_groups.py --chain-details PATH --source BASELINE=CSV --source STEERING=CSV --output-dir DIR`. Obtain `chain_details.csv` from `plot_daa.py` or supply the older research CSV with `sample`, `length`, and `d_residues`. Group each D-containing chain as **Single D** (one D residue), **Partial D** (more than one but fewer than chain length), or **All D** (every residue D). Require all D-containing chains of one input to have the same group. The comparison uses only inputs present and fully scored for every requested metric in both sources. Use `--expected-seeds N` to require an exact evaluated seed count per input and source; exclusions and reasons go into `excluded_inputs.csv`.

For each input and metric, choose its best seed independently: minimize RMSDs; maximize lDDT, DockQ, and **D-Cα chirality correctness**. The correctness denominator includes all targeted D sites, including unassessable sites. Group bars use the median of input oracle values for RMSDs and the arithmetic mean for bounded metrics, including system-equal D-Cα correctness. Boxplots show input oracle values (one dot per input), median and IQR, with 1.5-IQR whiskers. `seed_metrics.csv`, `input_metrics.csv`, and `group_metrics.csv` expose values, selected seeds, and aggregation. Colors follow the one-model or cross-model palette, according to source labels. The figure uses two metric rows, folding above docking.

## Folding/docking comparison

Use `plot_metrics.py compare --input-id ID` with one or more model sources. Eight panels are available: the four folding metrics, three requested docking metrics, and Interface RMSD as the eighth panel from the earlier comparison figures. Multi-metric figures put D-Cα chirality correctness first in the folding row and Target-aligned binder backbone RMSD first in the docking row, even when `--metrics` lists them later. By default each model must contribute exactly one prediction for that ID. If multiple seeds exist, explicitly choose `--selection top1` or `--selection oracle`; the latter selects a different seed per panel. If iRMSD is unavailable, omit it using `--metrics` and report that it was unavailable.

## Seed-selection bars

Use `selection-bars` for top-1 versus metric-specific oracle. Top-1 requires a `confidence_score` column, `--confidence-csv`, or adjacent Boltz confidence JSON files via `--confidence-root`. `--scope common` retains IDs with both methods; `--scope all` uses available IDs per method and is not a paired comparison. RMSDs are medians across selected inputs; lDDT, interface lDDT, DockQ, and D-Cα correctness are system-equal arithmetic means. Label each panel's Y axis with this aggregation. `--error-bars` shows sample SD of selected system-level values, **not** uncertainty of the aggregate. `--oracle-target chirality` selects the same D-Cα-oracle structure across panels; otherwise each metric panel has its own selected structure.
Oracle ties prefer higher D-Cα correctness, then lower target-aligned binder backbone RMSD, then the lowest seed ID.

## Chirality-versus-structure scatter

Use `scatter --selection all|top1|oracle`. The default six structure metrics plus optional Interface RMSD are plotted against **D-Cα chirality correctness**; chirality-versus-itself is omitted by default. For `oracle`, select the best seed separately for each Y metric, model, and input; its actual D-Cα correctness supplies X. With two methods, show only common input IDs for selected points; for `all`, only common input **and seed** IDs. Connect these with 0.6 pt light-gray dotted lines. `scatter_unmatched.csv` records unmatched observations. Correlations are descriptive; seeds from one input are repeated observations. Render Pearson r in the reserved band above each axis, below its panel title, and verify it does not cover data or overlap other labels at final width. Use `--lddt-label "Global lDDT" --ilddt-label "Interface lDDT"` for OpenStructure scores when the source has not been verified as the manuscript's heavy-atom lDDT implementation; document OpenStructure provenance in the caption. RMSD Y axes are linear from zero and cover all observations.

## Seed min–max

Use `seed-range` with two model sources, normally non-steering and steering. One figure is written per metric. Each row is an alphabetically ordered input ID, with separate positions for each method. A point represents each seed; a horizontal line spans the smallest to largest observed value. Missing values are omitted from that metric and remain visible in `normalized_predictions.csv`.

## Metric distributions

Use `distribution --input-id ID` with the baseline source first and steering source second. A violin, box, and individual points show baseline seeds; other points show each steered result. Colors follow the Notion palette, and marker shape identifies base model only. The typical comparison is 50 versus 3–5 seeds, but the script records actual counts. Both groups must have values for every requested metric; seeds remain repeated observations of one system.

## Ramachandran overlay

Use `ramachandran_compare.py --native FILE CHAIN LABEL --pred FILE CHAIN LABEL`. Native and predicted chain IDs must be explicitly supplied from known chain mappings. The script computes φ/ψ from backbone atoms N, CA, C, omits absent or disconnected torsions, exports per-structure CSVs, SVG single plots and an overlay, and a residue comparison CSV. Per-residue comparison uses residue numbers and insertion codes; verify numbering before interpreting differences. Read labels as structure identities, not quality classes.

All plotting scripts save SVG only by default, with editable Arial text and final-size layout. Pass `--formats svg,png,pdf` for additional PNG/PDF images; PNG renders at 300 dpi. Multi-metric figures place folding panels on the first row and docking panels on the second. Complete the manuscript caption and verify metric provenance before submission.

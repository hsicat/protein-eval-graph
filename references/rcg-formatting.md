# RCG paper figure rules

This is the bundled operational English guide for the user's Notion page “RCG Paper Figure/Table Formatting”, fetched 2026-11-04 and verified against the unchanged page on 2026-11-09. The [source snapshot](rcg-notion-source.md) preserves its text, omitting a private image reference and local script path. This local guide and snapshot are the formatting authority during plot generation; no Notion access is needed. Refresh them only when the user explicitly requests synchronization. These are manuscript working rules, subject to the target journal's final requirements. Apply them to main and supplementary figures.

## Group names, colors, and order

Choose colors by comparison purpose. **Within one base model**, use Baseline `#007282`, steering `#E69F00`, RCG `#009E73`, and RCG + steering `#CC79A6`. Preserve a group's color when plotting only two of the conditions. The PyMOL names are `marine`, `orange`, `aquamarine`, and `violet` respectively. Order four conditions Baseline → + steering → + RCG → + RCG + steering.

**When comparing base models in one figure**, use the exact cross-model palette:

| Model | Baseline | + steering | + RCG | + RCG + steering | Marker |
| --- | --- | --- | --- | --- | --- |
| AF3 | `#DACAB0` | `#BCA987` | `#E6A45E` | `#A96628` | circle `o` |
| Boltz-1 | `#C7E3EE` | `#8AB2BE` | `#56B4E9` | `#176A9C` | square `s` |
| Boltz-2 | `#C6E5C6` | `#8FB88F` | `#66C98B` | `#237B4B` | triangle `^` |

Ground truth is `#CCCCCC` in general plots. Boltz-2 RCG colors are **reserved**; never turn them into unevaluated results. Prefer separate panels for AF3, Boltz-1, and Boltz-2 (in that order) when comparing base models. If overlaid, retain the model marker. Markers identify base models only, never treatments or seeds. Include text or position as well as color so groups are distinguishable without color. Display full names in cross-model figures: `AF3`, `AF3 + steering`, `AF3 + RCG`, `AF3 + RCG + steering`, similarly for `Boltz-1` and `Boltz-2`. `Boltz 1x/2x` map to the corresponding `+ steering` names. When the base model is explicit in a panel title, `Baseline`, `+ steering`, `+ RCG`, `+ RCG + steering` are permitted. Never draw missing groups as zero; mark `Not evaluated` or omit.

## Final-size layout and typography

- Use **Arial**. At final size, a single-column figure is 90 mm wide; a double-column or multi-panel figure is 180 mm wide. Do not create at twice the width and shrink the finished figure without scaling text and lines accordingly.
- Default manuscript figures have left-aligned panel titles (9 pt bold) and lowercase panel letters (10 pt bold). Put the whole-figure title in the manuscript caption. An 11 pt bold, left-aligned whole-figure title is allowed for standalone discussion figures.
- Axis titles are 9 pt. Tick numbers, category labels, legends, in-panel statistics, structure labels, scale labels, and structure metrics are 8 pt. Colorbar title 9 pt and ticks 8 pt. A short bottom note may be 7.5 pt and at most two lines. Complete caption is 9 pt in the manuscript.
- Prefer horizontal category labels; wrap long labels instead of reducing font size. Keep equal panel heights within each row and equal panel widths within each column. Adjust figure height and margins before shrinking fonts. The user's two-row rule applies to multi-metric figures: folding metrics first, docking metrics second. A later user preference (2026-11-10) puts D-Cα chirality correctness first in the folding row and L-RMSD first in the docking row whenever those panels are present. Chirality scatter has no chirality-versus-itself panel.

## Standard metric names and provenance

The arrows `(↑ better)` and `(↓ better)` are display cues outside the metric names. Use these exact names and units on axes, wrapping long names as needed:

| Metric name | Unit | Direction | Meaning |
| --- | --- | --- | --- |
| D-Cα chirality correctness | % | ↑ | Correct target D-Cα sites divided by prespecified target D-Cα sites. |
| Systems with all D-Cα sites correct | % | ↑ | Systems with every target D-Cα site correct divided by included systems. |
| D-Cα margin pass rate | % | ↑ | Target D sites whose corrected directed angle has margin `q ≥ 15°`; distinct from general correctness. |
| Target self-aligned backbone RMSD | Å | ↓ | Target fitted independently to reference target. |
| Binder self-aligned backbone RMSD | Å | ↓ | Binder fitted independently to reference binder. |
| L-RMSD | Å | ↓ | Binder N/CA/C/O RMSD after jointly fitting all annotated target chains on their N/CA/C/O atoms; do not refit binder. |
| DAA-chain self-aligned backbone RMSD | Å | ↓ | Independently fit each relevant DAA chain to its reference. |
| DockQ | 0–1 | ↑ | Docking score; record software version, scoring mode, chain mapping and interface scope. |
| Global heavy-atom lDDT | 0–1 | ↑ | Local distance consistency on specified common heavy atoms. |
| Interface heavy-atom lDDT | 0–1 | ↑ | Local distance consistency on specified cross-chain heavy-atom pairs. |
| Selected peptide-bond constraint violation rate | % | ↓ | Violated selected bond constraints divided by assessable constraints. |
| Selected residue-geometry constraint violation rate | % | ↓ | Violated selected within-residue constraints divided by assessable constraints. |

For saved OpenStructure `ost_ilddt` plots, display the metric as **iLDDT** by default and define its target–binder interface scope and scoring provenance in the caption. This display name does not identify it with the separately defined Interface heavy-atom lDDT metric above.

The user's naming update (2026-11-24) replaces “Target-aligned binder backbone
RMSD” and “Binder pose RMSD” with **L-RMSD** in figures and tables, including
axes, panel titles, captions, and success-rate labels. Use **L-RMSD (Å)** and
**L-RMSD success**. This supersedes the old label in the preserved Notion source
snapshot; no Notion refresh is required. Keep the `placed` key and
`target_aligned_binder_backbone_rmsd_angstrom` source field unchanged. Define
the target fit and binder measurement in captions. Keep self-aligned RMSDs,
CA-only diagnostics, and OpenStructure per-interface `lrmsd` distinct.

The present graph scripts consume evaluation scores; they do **not** recompute them. Before publication, verify the score source implements the definition. The Notion page specifies N, Cα, C, O for backbone RMSD, common atoms, and proper rigid fits without mirroring. For the four-group analysis it specifies mean across included interfaces for system DockQ, median for interface-level target/binder RMSDs, and mean of independently fitted DAA-chain RMSDs. Document any different aggregation explicitly. Keep all-interface and DAA-related interface scores separate. Its current lDDT implementation is distance-only, pair-weighted on common heavy-atom pairs of different residues within about 15 Å reference distance, with 0.5/1/2/4 Å thresholds; interface lDDT uses cross-chain pairs. Do not claim these implementation details for an unverified CSV. These lDDT scores omit chemical symmetry and stereochemical penalties.

For chirality, the page defines corrected directed angle `q > 0°` as correct; `q = 0°` is not correct. The `q ≥ 15°` margin is a separate metric. Exclude sites without the specified chirality center. Keep target, assessable, correct and unassessable site counts. Missing target sites stay in the correctness denominator and must not be called confirmed inversions. The graph loader uses `chirality_daa_correct_count / chirality_daa_eligible_count`, not the evaluation CSV's assessable-site percentage or pooled D/L percentage. Check that input evaluation's DAA site set matches the intended target set.

## Statistical display

- State the observation unit: one system, one system–seed pair, one site, or one interface. For paired comparisons use the same systems; connect only corresponding observations. Multiple seeds from one system are repeated observations, not independent systems.
- State the structure selection rule: exported final structure, confidence top-1, a specified optimization state, or a named metric oracle. A metric-wise oracle plot selects independently per metric; it is **not** a profile from one structure. When metrics are claimed as a joint profile, they must come from the same selected structure file.
- Record missing results as `NA` with reasons. Compare matched available systems for paired analyses; never replace missing scores by zero. A boxplot shows median, 25th–75th percentiles, whiskers to the most extreme points within 1.5×IQR, and all raw points. Hiding duplicate outlier symbols is permitted only if raw points remain.
- Label error bars as SD, SEM, or 95% CI. Bootstrap system-level comparisons at the system level and retain pairing. Define change as new condition minus control. A percentage difference uses **percentage points**.
- Default display precision: correctness and rates 1 decimal, RMSD 2 decimals, DockQ/lDDT 3 decimals, counts integer. Compute at full precision; increase displayed precision consistently within a panel if a real effect would otherwise disappear.
- Across-system D-Cα correctness is the **system-equal mean**. If pooling sites instead, label it `Pooled-site D-Cα chirality correctness (%)`.

## Axes, marks, and structure figures

White `#FFFFFF` background; dark-gray `#222222` text/axes; 0.7 pt axes; 1.0 pt data outlines; 0.6 pt light-gray `#B8B8B8` paired lines; 0.7 pt dashed gray reference lines. Default points are about 3.5 pt and 0.85 opacity; distribution fills 0.25 opacity. Grid off unless a faint guide is necessary. Hide top/right statistical-plot spines and legend frames. Use 0–100% for correctness and 0–1 for DockQ/lDDT, allowing a small visual margin. Raw RMSDs start near zero and show all observations. No numeric-axis jitter. For a method A/B comparison use identical x/y ranges and a `y = x` reference line. Label zoom ranges and preserve full-range context.

For structure renderings: reference target `#C7C7C7`, reference binder `#555555`, predicted target `#969696`, predicted binder its treatment color. Use a black cross, ring or inset with legend for wrong-chirality sites. Hold view, zoom, crop and representation constant within a case. Align target for relative docking position; align binder separately for its own conformation and label that view `Binder self-aligned`. Put 8 pt metrics in consistent order, DockQ before RMSD. Ramachandran plots use φ/ψ axes from −180° to 180° and state the D-residue reference and site-filtering rule.

## Captions and delivery

Complete captions must identify the comparison, panels, observation unit, matched-system count, seeds, structure-selection rule, metric definitions and alignment, boxes/whiskers, error bars or CI, missing-data handling, and evaluation restrictions. Use the page's caption template. A short in-figure note may say `Dots: individual systems. Boxes: median and interquartile range.`

Keep the script, source data, chosen-structure list, editable figure, preview, and plot data together. Use a filename with figure number, content, version, e.g. `Fig05_AF3_RCG_DockQ_v01.svg`. Export statistics as SVG or PDF with editable text. PNG is a **requested preview** at 300 dpi; structure renders are at least 300 dpi, preferably 600 dpi. Use RGB and white background. Inspect final-width figures for readable text, clipping, complete points, and legend consistency. On this workstation the system Arial files carry restricted embedding flags: the scripts preserve editable Arial SVG text but do not embed the font; verify font availability/embedding requirements before external submission.

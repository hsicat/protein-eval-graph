"""Check aggregation and seed choices that materially affect scientific plots."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from plot_daa import extract  # noqa: E402
from plot_metrics import BY_KEY, choose, load_predictions, panels, parse_formats, selection_axis_label  # noqa: E402
from plot_daa_groups import prepare, read_groups  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


class DataRuleTests(unittest.TestCase):
    def test_metric_rows_and_requested_formats(self):
        fig, axes = panels(["dockq", "target", "placed", "lddt"], "Metric layout")
        try:
            self.assertEqual(axes["target"].get_subplotspec().rowspan.start, 0)
            self.assertEqual(axes["lddt"].get_subplotspec().rowspan.start, 0)
            self.assertEqual(axes["dockq"].get_subplotspec().rowspan.start, 1)
            self.assertEqual(axes["placed"].get_subplotspec().rowspan.start, 1)
            self.assertEqual(parse_formats("svg,png,pdf"), ("svg", "png", "pdf"))
            self.assertEqual(selection_axis_label(BY_KEY["target"]), "Median (Å)")
            self.assertEqual(selection_axis_label(BY_KEY["chirality"]), "Residue-pooled accuracy (%)")
        finally:
            plt.close(fig)

    def test_daa_groups_oracle_and_legacy_chain_csv(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chain_details.csv"
            path.write_text("sample,length,d_residues\nA,8,1\nB,8,8\nC,8,3\n", encoding="utf-8")
            groups = read_groups(path)
            self.assertEqual([groups[s]["group"] for s in "ABC"], ["Single D", "All D", "Partial D"])
            points = []
            for model in ("baseline", "steering"):
                for sample in "ABC":
                    for seed, rmsd, lddt in (("1", 4, .9), ("2", 2, .7)):
                        points.append({"model": model, "sample_id": sample, "seed": seed,
                                       "prediction_path": f"{model}/{sample}/{seed}", "chirality": 100,
                                       "target": rmsd, "lddt": lddt, "placed": 1})
            excluded, seeds, inputs, summary = prepare(points, groups, ["baseline", "steering"], ["target", "lddt"], 2)
            self.assertFalse(excluded)
            self.assertEqual(len(seeds), 24)
            self.assertEqual(next(r["oracle_seed"] for r in inputs if r["metric"] == "target"), "2")
            self.assertEqual(next(r["oracle_seed"] for r in inputs if r["metric"] == "lddt"), "1")
            self.assertEqual(next(r["group_oracle_value"] for r in summary if r["group"] == "All D" and r["metric"] == "target"), 2)

    def test_binder_groups_are_one_prediction_and_oracles_are_metric_specific(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "batch_summary.csv"
            fields = [
                "sample_id", "prediction_path", "status", "binder_self_aligned_backbone_rmsd_angstrom",
                "target_aligned_binder_backbone_rmsd_angstrom", "target_backbone_fit_rmsd_angstrom",
                "ost_lddt", "ost_ilddt", "ost_dockq_ave", "ost_per_interface_json",
                "chirality_all_accuracy_pct", "chirality_all_correct_count",
                "chirality_all_eligible_count", "chirality_all_skipped_count", "confidence_score",
            ]
            common = {"sample_id": "X001", "status": "SUCCESS", "target_backbone_fit_rmsd_angstrom": 1,
                      "ost_ilddt": 0.5, "ost_dockq_ave": 0.6,
                      "ost_per_interface_json": json.dumps([{"irmsd": 2}, {"irmsd": 4}]),
                      "chirality_all_eligible_count": 10, "chirality_all_skipped_count": 0}
            rows = [
                {**common, "prediction_path": "predicted/boltz1/00/X001/model.cif", "binder_self_aligned_backbone_rmsd_angstrom": 1,
                 "target_aligned_binder_backbone_rmsd_angstrom": 3, "ost_lddt": .8,
                 "chirality_all_accuracy_pct": 50, "chirality_all_correct_count": 5, "confidence_score": .9},
                {**common, "prediction_path": "predicted/boltz1/00/X001/model.cif", "binder_self_aligned_backbone_rmsd_angstrom": 5,
                 "target_aligned_binder_backbone_rmsd_angstrom": 7, "ost_lddt": .8,
                 "chirality_all_accuracy_pct": 50, "chirality_all_correct_count": 5, "confidence_score": .9},
                {**common, "prediction_path": "predicted/boltz1/01/X001/model.cif", "binder_self_aligned_backbone_rmsd_angstrom": 4,
                 "target_aligned_binder_backbone_rmsd_angstrom": 4, "ost_lddt": .9,
                 "chirality_all_accuracy_pct": 100, "chirality_all_correct_count": 10, "confidence_score": .7},
            ]
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(rows)
            points = load_predictions([f"Boltz1={path}::boltz1"], None, None)
            self.assertEqual(len(points), 2)
            self.assertEqual(points[0]["binder"], 3)
            self.assertEqual(points[0]["placed"], 5)
            self.assertEqual(points[0]["irmsd"], 3)
            self.assertEqual(choose(points, "top1", "lddt")["seed"], "00")
            self.assertEqual(choose(points, "oracle", "lddt")["seed"], "01")
            self.assertEqual(choose(points, "oracle", "binder")["seed"], "00")
            tied = [{**points[0], "chirality": 100}, {**points[1], "chirality": 100}]
            self.assertEqual(choose(tied, "oracle", "chirality")["seed"], "01")

    def test_daa_count_expands_physical_chains_and_excludes_other_chains(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "X001_ccd.yaml").write_text(
                "sequences:\n"
                "  - protein:\n"
                "      id: [A, B]\n"
                "      sequence: AAAAA\n"
                "      modifications: [{ccd: DAL, position: 2}]\n"
                "  - protein:\n"
                "      id: C\n"
                "      sequence: AAAAAAAAAA\n",
                encoding="utf-8",
            )
            samples, chains, skipped = extract(folder, "auto")
            self.assertEqual(skipped, [])
            self.assertEqual(len(chains), 2)
            self.assertEqual(samples[0]["d_residues"], 2)
            self.assertEqual(samples[0]["d_chain_residues"], 10)
            self.assertAlmostEqual(samples[0]["d_fraction"], .2)


if __name__ == "__main__":
    unittest.main()

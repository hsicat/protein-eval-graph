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
from plot_metrics import choose, load_predictions  # noqa: E402


class DataRuleTests(unittest.TestCase):
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

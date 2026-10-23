"""Regression checks for G6 candidate panel triage."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

import yaml


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "src/99_screen_g6_candidate_panel_004.py"
SPEC = importlib.util.spec_from_file_location("g6_panel_triage_004", SCRIPT)
assert SPEC and SPEC.loader
TRIAGE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(TRIAGE)
CONFIG = yaml.safe_load((PROJECT / TRIAGE.CONFIG).read_text(encoding="utf-8"))


class CandidatePanelTriageTest(unittest.TestCase):
    def test_below_25_percent_still_enters_blind_validation(self) -> None:
        status, next_step = TRIAGE.route("SIZE_PASS_PENDING_FULL_AUDIT", [], ["G1", "M5_Iron_Uptake"], CONFIG)
        self.assertEqual(status, "BLIND_VALIDATION_ELIGIBLE_LOW_COVERAGE")
        self.assertEqual(next_step, "RUN_EFFECT_BLIND_MEASUREMENT_VALIDATION")
        self.assertFalse(CONFIG["diagnosis_group_effects_allowed"])
        self.assertFalse(CONFIG["full_primary_eligibility_from_triage"])

    def test_unknown_panel_and_unknown_donor_count_stay_pending(self) -> None:
        self.assertEqual(TRIAGE.route("SIZE_PASS_PENDING_FULL_AUDIT", ["G1"], [], CONFIG)[0],
                         "PENDING_PANEL_EVIDENCE_FOR_BLIND_VALIDATION")
        self.assertEqual(TRIAGE.route("HOLD_DONOR_COUNTS", [], [], CONFIG)[0], "HOLD_DONOR_COUNTS")


if __name__ == "__main__":
    unittest.main()

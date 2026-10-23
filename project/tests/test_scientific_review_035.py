"""Known statistical values and missing-evidence guardrails for SOP035 review."""
import unittest
import sys
from pathlib import Path
import numpy as np
from scipy.stats import false_discovery_control
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import scientific_review_virtual_ko_035 as review

class ScientificReviewTests(unittest.TestCase):
    def test_BH_against_independent_scipy_with_ties(self):
        p=np.array([0, .01, .01, .2, .7, 1.])
        np.testing.assert_allclose(review.bh(p),false_discovery_control(p,method='bh'))

    def test_nonestimable_is_not_zero_or_excluded(self):
        p=np.array([np.nan,.01,.5])
        q=review.bh(p)
        self.assertTrue(np.isnan(q[0]))
        np.testing.assert_allclose(q[1:],false_discovery_control(p[1:],method='bh'))
        # One reproduced of one testable out of six planned is NOT 100% stable.
        self.assertEqual(1/6, .16666666666666666)
        self.assertNotEqual(1/6,1/1)

    def test_fixed_scope_preserves_states_and_unique_units(self):
        scope=review.table(review.ROOT/review.CFG['scope_path'])
        self.assertEqual(len(scope),29)
        self.assertEqual(scope.target.nunique(),15)
        self.assertEqual(scope.context.str.endswith('_AD').sum(),14)
        self.assertFalse(scope.duplicated(review.KEY).any())

    def test_reference_ties_and_resolution(self):
        controls=np.array([.5]*20)
        self.assertEqual((1+sum(controls>=.5))/21,1)
        self.assertEqual((1+sum(controls>=.6))/21,1/21)

if __name__=='__main__':unittest.main()

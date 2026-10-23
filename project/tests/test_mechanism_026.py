"""Focused safeguards for statistical and therapeutic interpretation boundaries."""
import importlib.util,json,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('audit026',ROOT/'src/133_report_audit_mechanism_026.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class Safeguards(unittest.TestCase):
    def test_BH_known_values(self):np.testing.assert_allclose(m.bh([.01,.04,.03,.2]),[.04,.0533333333333,.0533333333333,.2])
    def test_BH_permutation_invariant(self):
        p=np.array([.01,.5,.03,.07]);order=[3,0,2,1];np.testing.assert_allclose(m.bh(p)[order],m.bh(p[order]))
    def test_factorial_interaction_direction(self):self.assertEqual(m.contrast('(BD-B)-(D-C)',dict(BD=7.,B=2.,D=4.,C=1.)),2.)
    def test_unsafe_contrast_rejected(self):
        with self.assertRaises(ValueError):m.contrast('__import__("os")',{})
    def test_TQ_gate_stays_closed(self):
        d=json.loads((m.OUT/'TQ_D1_D2_audit.json').read_text());self.assertEqual(d['normalization_eligible'],0);self.assertEqual(d['primary_therapeutic_eligible'],0);self.assertFalse(d['reopen_017'])
    def test_original_sources_not_inferred_from_figures(self):
        d=json.loads((m.RAW/'FIGSHARE_30252505_metadata.json').read_text());self.assertEqual(len(d['files']),11);self.assertTrue(all(f['name'].lower().endswith('.tif') for f in d['files']))
if __name__=='__main__':unittest.main()

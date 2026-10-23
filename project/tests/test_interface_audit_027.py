"""Scientific guardrails against accidental qualification or pseudoreplication."""
import importlib.util,json,unittest
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/2026-10-03-disease-interface-perturbation-audit-027'
spec=importlib.util.spec_from_file_location('b',ROOT/'src/138_build_interface_qualification_027.py');b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
class Guardrails(unittest.TestCase):
    def test_unknown_exposure_cannot_pass(self):
        keys=['formula_identity','therapeutic_direction','matched_reference','exposure','viability_stress','context_qualification'];r={k:'PASS' for k in keys};self.assertTrue(b.p5_pass(r));r['exposure']='UNKNOWN';self.assertFalse(b.p5_pass(r))
    def test_genetic_not_compound(self):
        self.assertEqual(b.role('GSE290755','KD-NTC')[:2],('genetic_mechanism','PC-E'))
    def test_vehicle_not_rescue(self):
        self.assertEqual(b.role('PXD079947','DMSO-NoTreatment')[0],'vehicle_control');self.assertEqual(b.role('PXD079947','IronRSL3Fer-IronRSL3')[0],'external_rescue_compound')
    def test_g1_signs_preserved(self):
        genes=pd.read_csv(OUT/'D0_interface_genes.tsv',sep='\t');source=pd.read_csv(OUT/'D0_original_effects/G1.tsv',sep='\t');g=genes[genes.module=='G1_robust244'].set_index('gene');self.assertEqual(set(g.index),set(source.gene));self.assertEqual(len(g),244)
        for r in source.itertuples():self.assertEqual(g.loc[r.gene,'weight'],1 if r.meta_logFC>0 else -1)
    def test_no_automatic_therapy_direction(self):
        d=pd.read_csv(OUT/'D0_interface_registry.tsv',sep='\t');self.assertFalse(d.primary_therapeutic_eligible.any());self.assertFalse(d.normalization_eligible.any());self.assertFalse(d.activity_or_communication_direct_measurement.any())
    def test_spatial_units_not_143_animals(self):
        d=pd.read_csv(OUT/'candidate_sample_registry.tsv',sep='\t');s=d[d.accession=='GSE333352'];self.assertEqual(len(s),143);self.assertEqual(s.independent_unit.nunique(),5);self.assertLess(s.ROI.nunique(),143)
    def test_mixed_candidate_not_pooled(self):
        d=pd.read_csv(OUT/'candidate_sample_registry.tsv',sep='\t');s=d[d.accession=='GSE315101'];self.assertEqual(s.cohort.nunique(),4);self.assertEqual(len(s[s.cohort=='aged_cGAS_control_only']),3)
    def test_access_not_structure_validation(self):
        d=pd.read_csv(OUT/'candidate_file_access.tsv',sep='\t');self.assertTrue(d.accessible.all());self.assertFalse(d.body_validated.any())
    def test_provisional_samples_retained(self):
        d=pd.read_csv(OUT/'P4_sample_traceability.tsv',sep='\t');self.assertEqual((d.mapping_status=='UNRESOLVED').sum(),2)
    def test_p7_coverage_no_therapy_pass(self):
        d=pd.read_csv(OUT/'P7_shared_measurement_coverage.tsv',sep='\t');self.assertTrue((d.shared_genes<=d.interface_genes).all());self.assertFalse(d.primary_therapeutic_eligible.any());self.assertIn('G5_M5_Iron_Storage_ortholog_exploratory',set(d.interface))
if __name__=='__main__':unittest.main()

"""Input correctness and scientific rejection guards for 028."""
import importlib.util,json,unittest
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/2026-10-03-gse286094-mechanism-028'
spec=importlib.util.spec_from_file_location('audit',ROOT/'src/142_audit_gse286094_028.py');audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
class Checks(unittest.TestCase):
    def test_dimensions_reject(self):
        with self.assertRaises(ValueError):audit.validate_counts(csc_matrix([[1,2]]),2,2)
    def test_fractional_counts_reject(self):
        with self.assertRaises(ValueError):audit.validate_counts(csc_matrix([[1.5,2]]),1,2)
    def test_negative_counts_reject(self):
        with self.assertRaises(ValueError):audit.validate_counts(csc_matrix([[-1,2]]),1,2)
    def test_count_conservation(self):
        d=pd.read_csv(OUT/'all_unfiltered_barcode_QC.tsv.gz',sep='\t');s=pd.read_csv(OUT/'matrix_QC_diagnostics.tsv',sep='\t');observed=d.groupby('GSM').n_UMI.sum();self.assertEqual(len(d),int(s.barcodes.sum()))
        for r in s.itertuples():self.assertEqual(observed[r.GSM],r.UMI_sum)
    def test_namespaced_barcodes(self):
        d=pd.read_csv(OUT/'all_unfiltered_barcode_QC.tsv.gz',sep='\t',usecols=['cell_id']);self.assertFalse(d.cell_id.duplicated().any())
    def test_no_identity_imputation(self):
        d=pd.read_csv(OUT/'sample_identity_audit.tsv',sep='\t');self.assertTrue(d.declared_biological_replicate.all());self.assertFalse(d.independence_qualified.any());self.assertTrue((d.animal_id=='NOT_SUPPLIED').all())
    def test_no_pseudoreplicated_inference(self):
        d=json.loads((OUT/'analysis_status.json').read_text());self.assertEqual(d['independence_gate'],'NOT_PASSED');self.assertTrue(d['no_inferential_P_or_FDR']);self.assertFalse(d['thresholds_selected'])
if __name__=='__main__':unittest.main()

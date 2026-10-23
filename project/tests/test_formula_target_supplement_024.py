import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

SPEC=importlib.util.spec_from_file_location('targets024',Path(__file__).resolve().parents[1]/'src/113_supplement_formula_targets_024.py')
m=importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(m)
FSPEC=importlib.util.spec_from_file_location('final024',Path(__file__).resolve().parents[1]/'src/114_finalize_formula_targets_024.py')
f=importlib.util.module_from_spec(FSPEC); FSPEC.loader.exec_module(f)


class TargetRules(unittest.TestCase):
    def record(self):
        return dict(identity_exact=True,target_stable_id='UniProt:P35354',target_symbol='PTGS2',species='9606',single_protein=True,assay_mode='binding',source='ChEMBL',confidence=9,value='10',units='nM',relation='=',assay_id='A1',document='doi:example',activity_outcome='active',activity_type='IC50')

    def test_qualified(self): self.assertEqual(m.classify(self.record())[0],'E1')
    def test_similarity(self): self.assertEqual(m.classify({**self.record(),'similarity_only':True})[0],'S')
    def test_stereo(self):
        self.assertNotEqual(m.canonical('N[C@@H](C)C(=O)O'),m.canonical('N[C@H](C)C(=O)O'))
        self.assertEqual(m.classify({**self.record(),'identity_exact':False})[0],'E2')
    def test_nonhuman(self): self.assertEqual(m.classify({**self.record(),'species':'10090'})[0],'E2')
    def test_original_assay_origin_ambiguity(self): self.assertEqual(m.classify({**self.record(),'assay_description':'Inhibition of LYP (unknown origin)'})[0],'E2')
    def test_inactive(self): self.assertEqual(m.classify({**self.record(),'activity_outcome':'Inactive'})[0],'E2')
    def test_provider_inactive_comment_is_used(self):
        r={**self.record(),'activity_outcome':m.chembl_outcome('Not Active')}
        self.assertEqual(m.classify(r)[0],'E2'); self.assertFalse(f.positive(r))
    def test_source_no_measured_curve_not_quantitative_hit(self):
        r={**self.record(),'activity_outcome':m.chembl_outcome('Inhibition < 50% @ 10 uM and thus dose-reponse curve not measured')}
        self.assertEqual(m.classify(r)[0],'E2'); self.assertFalse(f.positive(r))
    def test_units(self): self.assertEqual(m.classify({**self.record(),'units':''})[0],'E2')
    def test_rate_constant_not_molar_endpoint(self): self.assertEqual(m.classify({**self.record(),'activity_type':'Kinact'})[0],'E2')
    def test_inactive_not_candidate_intersection(self): self.assertFalse(f.positive({**self.record(),'activity_outcome':'Inactive'}))
    def test_unknown_outcome_not_candidate_intersection(self): self.assertFalse(f.positive({**self.record(),'activity_outcome':''}))
    def test_lower_bound_not_positive_confirmation(self): self.assertFalse(f.positive({**self.record(),'relation':'>'}))
    def test_identity_hold_not_candidate_intersection(self): self.assertFalse(f.positive({**self.record(),'candidate_identity_verified':False}))
    def test_phenotype(self): self.assertEqual(m.classify({**self.record(),'target_type':'CELL-LINE','assay_mode':'functional_or_phenotypic'})[0],'E2')
    def test_confidence(self): self.assertEqual(m.classify({**self.record(),'confidence':8})[0],'E2')
    def test_association_not_prediction(self): self.assertEqual(m.classify({**self.record(),'association_only':True})[0],'H')
    def test_prediction_requires_method(self): self.assertEqual(m.classify({'identity_exact':True,'predicted':True})[0],'H')
    def test_cross_source_duplicate(self):
        a={**self.record(),'smiles':'C','original_measurement_id':'AID:100:SID:200','source':'ChEMBL'}
        self.assertEqual(m.evidence_key(a),m.evidence_key({**a,'source':'PubChem'}))
    def test_stable_ID_namespaces_normalized_before_dedup(self):
        a={**self.record(),'smiles':'C','target_hgnc_id':'HGNC:9605','original_measurement_id':'AID:100:SID:200','source':'ChEMBL'}
        b={**a,'source':'PubChem','target_stable_id':'NCBI_Gene:5743'}
        self.assertEqual(m.evidence_key(a),m.evidence_key(b))
    def test_different_assays_not_merge(self):
        a={**self.record(),'smiles':'C','activity_id':'100','compound_key':'C'}
        self.assertNotEqual(m.evidence_key(a),m.evidence_key({**a,'activity_id':'200'}))
    def test_different_accession_annotations_retained(self):
        a={**self.record(),'activity_id':'100','compound_key':'C','target_accession':'2O23_A'}
        self.assertNotEqual(m.evidence_key(a),m.evidence_key({**a,'target_accession':'2O23_B'}))
    def test_conflicting_activity_outcomes_retained(self):
        a={**self.record(),'activity_id':'100','compound_key':'C'}
        self.assertNotEqual(m.evidence_key(a),m.evidence_key({**a,'activity_outcome':'inactive'}))
    def test_parse_summary(self):
        d={'Table':{'Columns':{'Column':['AID','CID','Activity Outcome']},'Row':[{'Cell':['1','2','Inactive']}]}}
        self.assertEqual(m.pubchem_rows(d)[0]['Activity Outcome'],'Inactive')
    def test_unknown_summary_fails(self):
        with self.assertRaises(m.SourceError): m.pubchem_rows({})
    def test_incomplete_page_fails(self):
        c=object.__new__(m.Client)
        c.json=lambda *a: ({'activities':[{}],'page_meta':{'total_count':2,'next':None}},{'sha256':'x'})
        with self.assertRaises(m.SourceError): c.pages('https://example.org','activities')
    def test_cycle_fails(self):
        c=object.__new__(m.Client)
        c.json=lambda *a: ({'activities':[],'page_meta':{'total_count':0,'next':'https://example.org'}},{'sha256':'x'})
        with self.assertRaises(m.SourceError): c.pages('https://example.org','activities')
    def test_restart_cache_no_network(self):
        d=Path(__file__).resolve().parents[1]/'.cache/024-restart-test'
        d.mkdir(parents=True,exist_ok=True)
        try:
            c=m.Client('test',{'http':{'interval_seconds':1,'timeout_seconds':30,'backoff_seconds':[2,5,10]}},d)
            url='https://example.org/data'; key=m.digest(url.encode()); raw=b'{"ok":true}'
            (c.raw/(key+'.body')).write_bytes(raw)
            (c.raw/(key+'.json')).write_text(json.dumps({'status':200,'sha256':m.digest(raw)}))
            c.session.get=lambda *a,**kw: (_ for _ in ()).throw(AssertionError('network called'))
            self.assertEqual(c.get(url)[0],raw)
        finally:
            for p in (d/'test').glob('*'): p.unlink()
            (d/'test').rmdir(); d.rmdir()


if __name__=='__main__': unittest.main()

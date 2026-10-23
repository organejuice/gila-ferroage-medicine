import importlib.util
from pathlib import Path
import unittest

spec=importlib.util.spec_from_file_location('herb024',Path(__file__).resolve().parents[1]/'src/119_query_herb_candidates_024.py')
h=importlib.util.module_from_spec(spec); spec.loader.exec_module(h)


class TestHERBIdentity(unittest.TestCase):
    def summary(self,cid='6054',smiles='OCCc1ccccc1'):
        return [['PubChem id','Molecule smile'],[{'title':cid},smiles]]

    def test_exact(self):
        self.assertTrue(h.identity(self.summary(),{'smiles':'OCCc1ccccc1'},'6054')[0])

    def test_wrong_cid(self):
        self.assertFalse(h.identity(self.summary('1'),{'smiles':'OCCc1ccccc1'},'6054')[0])

    def test_wrong_structure(self):
        self.assertFalse(h.identity(self.summary(smiles='CCO'),{'smiles':'OCCc1ccccc1'},'6054')[0])

    def test_stereo(self):
        self.assertFalse(h.identity(self.summary(smiles='N[C@H](C)O'),{'smiles':'C[C@H](N)O'},'6054')[0])

    def test_missing_cid(self):
        self.assertFalse(h.identity(self.summary(''),{'smiles':'OCCc1ccccc1'},'6054')[0])

    def test_width(self):
        with self.assertRaises(ValueError): h.table([['A','B'],['a']])

    def test_duplicate_header(self):
        with self.assertRaises(ValueError): h.table([['A','A'],['a','b']])

    def test_h_remains_h(self):
        self.assertEqual(h.m.classify({'association_only':True,'identity_exact':True,'document':'1234'})[0],'H')

    def test_cell_list(self):
        self.assertEqual(h.cell([{'title':'A'}, {'title':'B'}]),'A;B')


if __name__=='__main__': unittest.main()

"""Audit prepared chemical identities and receptor cofactors; no docking claims."""
import csv
import hashlib
import json
import unittest
from pathlib import Path

import numpy as np
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
RUN = '2026-10-04-stage-b-tcmsp-docking-032'


class InputAudit(unittest.TestCase):
    def test_manifest_and_structures(self):
        directory = ROOT / 'results' / RUN
        with (directory / 'ligand_manifest.tsv').open(encoding='utf-8') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        self.assertEqual(len(rows), 23)
        self.assertEqual(sum(r['target'] == 'AKR1B1' for r in rows), 21)
        self.assertEqual(sum(r['target'] == 'LTA4H' for r in rows), 2)
        self.assertEqual(len({r['compound_key'] for r in rows if r['role'] == 'candidate'}), 16)
        for row in rows:
            with (ROOT / row['file']).open('rb') as handle:
                molecule = next(iter(Chem.ForwardSDMolSupplier(handle, removeHs=False)))
            self.assertIsNotNone(molecule)
            self.assertEqual(Chem.MolToInchiKey(Chem.RemoveHs(molecule)), row['inchikey'])
            self.assertEqual(Chem.GetFormalCharge(molecule), int(row['formal_charge']))
            self.assertEqual(molecule.GetNumHeavyAtoms(), int(row['heavy_atoms']))
            self.assertTrue(np.isfinite(molecule.GetConformer().GetPositions()).all())
        for relative, expected in json.loads((directory / 'input_hashes.json').read_text()).items():
            self.assertEqual(hashlib.sha256((ROOT / relative).read_bytes()).hexdigest(), expected)

    def test_receptor_and_reference(self):
        config = json.loads((ROOT / 'configs/stage_b_tcmsp_docking_032.json').read_text())
        for gene, target in config['targets'].items():
            lines = (ROOT / target['receptor']).read_text().splitlines()
            hetero = {line[17:20].strip() for line in lines if line.startswith('HETATM')}
            self.assertEqual(hetero, {target['retain_residue']})
            self.assertNotIn(target['reference_ligand'], hetero)
            self.assertTrue(all(line[21] == 'A' for line in lines if line.startswith(('ATOM', 'HETATM'))))
            ligand = ROOT / 'data/raw' / RUN / (target['pdb'] + '_' + target['reference_ligand'] + '_instance.sdf')
            with ligand.open('rb') as handle:
                molecule = next(iter(Chem.ForwardSDMolSupplier(handle, removeHs=True)))
            self.assertEqual(molecule.GetNumHeavyAtoms(), 24 if gene == 'AKR1B1' else 19)
            np.testing.assert_allclose(molecule.GetConformer().GetPositions().mean(axis=0), target['center'], atol=0.00051)


if __name__ == '__main__':
    unittest.main()

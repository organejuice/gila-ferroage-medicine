"""Prepare traceable inputs for the user-requested TCMSP GUI docking.

Actual docking is performed by TCMSP AI Agent, not by this script.
"""
import csv
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import requests
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

ROOT = Path(__file__).resolve().parents[1]
RUN = '2026-10-04-stage-b-tcmsp-docking-032'
RAW = ROOT / 'data/raw' / RUN
OUT = ROOT / 'data/processed' / RUN
RESULT = ROOT / 'results' / RUN


def fetch(url, path):
    if not path.exists():
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        path.write_bytes(response.content)
    return path


def embed(smiles, name):
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(name)
    molecule = Chem.AddHs(molecule)
    settings = AllChem.ETKDGv3()
    settings.randomSeed = 42
    if AllChem.EmbedMolecule(molecule, settings) != 0:
        raise RuntimeError(f'Embedding failed: {name}')
    if not AllChem.MMFFHasAllMoleculeParams(molecule):
        raise RuntimeError(f'MMFF parameters missing: {name}')
    status = AllChem.MMFFOptimizeMolecule(molecule, maxIters=2000)
    if status != 0:
        raise RuntimeError(f'MMFF not converged: {name}')
    molecule.SetProp('_Name', name)
    return molecule


def carboxylate(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    matches = molecule.GetSubstructMatches(pattern)
    if not matches:
        return None
    mutable = Chem.RWMol(molecule)
    for match in matches:
        atom = mutable.GetAtomWithIdx(match[2])
        atom.SetFormalCharge(-1)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(True)
    Chem.SanitizeMol(mutable)
    return Chem.MolToSmiles(mutable)


def prepare():
    for directory in (RAW, OUT, RESULT):
        directory.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'results/2026-10-03-stage-a-mechanism-030/mechanism_compound_target_matrix.tsv.gz'
    with gzip.open(source, 'rt', encoding='utf-8') as handle:
        candidates = [row for row in csv.DictReader(handle, delimiter='\t')
                      if row['gene'] in ('AKR1B1', 'LTA4H') and row['tier'] == 'E1'
                      and row['disease_set'] == 'G1_robust244']
    assert len(candidates) == 16
    manifest, targets = [], {}
    for gene, pdb, ligand, keep in [('AKR1B1', '2FZD', 'TOL', 'NAP'), ('LTA4H', '4L2L', '1V6', 'ZN')]:
        path = fetch(f'https://files.rcsb.org/download/{pdb}.pdb', RAW / f'{pdb}.pdb')
        ligand_file = fetch(f'https://models.rcsb.org/v1/{pdb.lower()}/ligand?auth_asym_id=A&auth_comp_id={ligand}&encoding=sdf', RAW / f'{pdb}_{ligand}_instance.sdf')
        with ligand_file.open('rb') as handle:
            reference = next(iter(Chem.ForwardSDMolSupplier(handle, removeHs=True)))
        if reference is None:
            raise ValueError('Invalid reference SDF')
        reference_smiles = Chem.MolToSmiles(reference)
        xyz = reference.GetConformer().GetPositions()
        center = xyz.mean(axis=0).round(3).tolist()
        # Retain dominant alternate location for each protein/cofactor atom.
        atom_records = {}
        for line in path.read_text().splitlines():
            if line[:6].strip() not in ('ATOM', 'HETATM') or line[21] != 'A':
                continue
            if line.startswith('HETATM') and line[17:20].strip() != keep:
                continue
            key = (line[22:27], line[17:20], line[12:16])
            occupancy = float(line[54:60])
            if key not in atom_records or occupancy > atom_records[key][0]:
                atom_records[key] = (occupancy, line[:16] + ' ' + line[17:])
        receptor = OUT / f'{gene}_{pdb}_receptor.pdb'
        receptor.write_text('\n'.join(record[1] for record in atom_records.values()) + '\nEND\n')
        assert ligand not in {line[17:20].strip() for _, line in atom_records.values()}
        assert keep in {line[17:20].strip() for _, line in atom_records.values()}
        target_dir = OUT / gene
        target_dir.mkdir(exist_ok=True)
        items = [row for row in candidates if row['gene'] == gene]
        for row in items:
            expected = Chem.MolToInchiKey(Chem.MolFromSmiles(row['smiles']))
            assert expected == row['inchikey'], (row['compound_key'], expected)
        items.append(dict(gene=gene, compound_key=f'CONTROL:{ligand}', cids='', herbs='', smiles=reference_smiles,
                          inchikey=Chem.MolToInchiKey(reference), names=f'{ligand}_cocrystal_positive_control'))
        for number, row in enumerate(items, 1):
            for state, smiles in [('catalogue', row['smiles']), ('carboxylate', carboxylate(row['smiles']))]:
                if smiles is None:
                    continue
                name = f'{gene}_{number:02d}_' + ('CONTROL_' + ligand if row['compound_key'].startswith('CONTROL') else 'CID' + row['cids']) + '_' + state
                molecule = embed(smiles, name)
                output = target_dir / f'{name}.sdf'
                with output.open('w', encoding='utf-8') as handle:
                    with Chem.SDWriter(handle) as writer:
                        writer.write(molecule)
                manifest.append(dict(target=gene, ligand_id=name, role='positive_control' if row['compound_key'].startswith('CONTROL') else 'candidate',
                                     compound_key=row['compound_key'], cids=row['cids'], names=row['names'], herbs=row['herbs'],
                                     state=state, smiles=smiles, inchikey=Chem.MolToInchiKey(Chem.RemoveHs(molecule)),
                                     formal_charge=Chem.GetFormalCharge(molecule), heavy_atoms=molecule.GetNumHeavyAtoms(), file=str(output.relative_to(ROOT))))
        targets[gene] = dict(pdb=pdb, reference_ligand=ligand, receptor=str(receptor.relative_to(ROOT)),
                             center=center, box=[24.0, 24.0, 24.0], retain_residue=keep,
                             rationale='Full resolved human WT inhibitor complex; centroid of crystallographic inhibitor; 24 A box covers ligand and pocket, not blind whole protein',
                             waters='removed; not modeled', alternate_locations='maximum occupancy per atom; ties first record')
    with (RESULT / 'ligand_manifest.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest[0].keys(), delimiter='\t'); writer.writeheader(); writer.writerows(manifest)
    config = dict(run=RUN, preparation_seed=42, conformer_method='RDKit ETKDGv3 then MMFF94 maxIters=2000, convergence required',
                  protonation='catalogue exact state plus carboxylate sensitivity; no claim of pH population or full tautomer enumeration',
                  targets=targets, docking_engine='TCMSP AI Agent GUI cloud AutoDock Vina; backend version/seed/exhaustiveness to inspect from results',
                  validation='symmetry-aware heavy atom RMSD to co-crystal in fixed receptor coordinates; no ligand fitting',
                  source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(), environment=dict(python=sys.version, rdkit=rdBase.rdkitVersion))
    (ROOT / 'configs/stage_b_tcmsp_docking_032.json').write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding='utf-8')
    files = list(RAW.glob('*')) + list(OUT.rglob('*.pdb')) + list(OUT.rglob('*.sdf'))
    (RESULT / 'input_hashes.json').write_text(json.dumps({str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}, indent=2), encoding='utf-8')
    print(json.dumps(dict(targets=targets, ligand_inputs=len(manifest)), indent=2))


if __name__ == '__main__':
    prepare()

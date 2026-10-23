"""Audit actual TCMSP output with Meeko export and RDKit fixed-coordinate RMSD."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from meeko import PDBQTMolecule, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[1]
RUN = '2026-10-04-stage-b-tcmsp-docking-032'
RESULT = ROOT / 'results' / RUN
INTERIM = ROOT / 'data/interim' / RUN


def run():
    with (RESULT / 'ligand_manifest.tsv').open(encoding='utf-8') as handle:
        manifest = {r['ligand_id']: r for r in csv.DictReader(handle, delimiter='\t')}
    config = json.loads((ROOT / 'configs/stage_b_tcmsp_docking_032.json').read_text(encoding='utf-8'))
    table, controls, receptor_audit = [], [], []
    uncharger = rdMolStandardize.Uncharger()
    for gene, target in config['targets'].items():
        folder = INTERIM / gene
        summary = json.loads((folder / 'docking-summary.json').read_text(encoding='utf-8'))
        expected_ids = {i for i, m in manifest.items() if m['target'] == gene}
        assert {r['ligand_id'] for r in summary['results']} == expected_ids
        assert summary['failed_tasks'] == 0
        np.testing.assert_allclose([summary['box']['center_' + axis] for axis in 'xyz'], target['center'])
        complex_files = list(folder.glob('complex_poses*.pdbqt'))
        kept_per_complex = []
        input_lines = (ROOT / target['receptor']).read_text(encoding='utf-8').splitlines()
        protein_lookup = {(l[12:16], l[17:27]): np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
                          for l in input_lines if l.startswith('ATOM')}
        coordinate_errors = []
        for complex_file in complex_files:
            lines = complex_file.read_text(encoding='utf-8').splitlines()
            kept_per_complex.append(sum(l.startswith(('ATOM', 'HETATM')) and l[17:20].strip() == target['retain_residue'] for l in lines))
            for line in lines:
                key = (line[12:16], line[17:27])
                if line.startswith(('ATOM', 'HETATM')) and key in protein_lookup:
                    coordinate_errors.append(float(np.linalg.norm(np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]) - protein_lookup[key])))
        receptor_audit.append(dict(target=gene, requested_cofactor=target['retain_residue'],
                                   complex_files=len(complex_files), cofactor_atoms_per_complex=kept_per_complex,
                                   cofactor_present_in_export=all(n > 0 for n in kept_per_complex),
                                   backend_prepared_receptor_file_delivered=False,
                                   matched_protein_atoms=len(coordinate_errors), maximum_protein_coordinate_difference_A=max(coordinate_errors),
                                   vina_params=summary['vina_params']))
        reference_path = ROOT / 'data/raw' / RUN / f"{target['pdb']}_{target['reference_ligand']}_instance.sdf"
        with reference_path.open('rb') as handle:
            reference = next(iter(Chem.ForwardSDMolSupplier(handle, removeHs=True)))
        reference = uncharger.uncharge(reference)
        receptor_lines = (ROOT / target['receptor']).read_text(encoding='utf-8').splitlines()
        protein_atoms = [(l[17:20].strip() + ':' + l[21:27].strip(), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
                         for l in receptor_lines if l.startswith('ATOM') and l[76:78].strip() != 'H']
        for record in summary['results']:
            ligand_id = record['ligand_id']
            metadata = manifest[ligand_id]
            files = list(folder.glob('*_' + ligand_id + '_out.pdbqt'))
            assert len(files) == 1
            text = files[0].read_text(encoding='utf-8')
            pdbqt = PDBQTMolecule(text, skip_typing=True)
            exported = RDKitMolCreate.from_pdbqt_mol(pdbqt)
            assert len(exported) == 1 and exported[0] is not None
            molecule = Chem.RemoveHs(exported[0])
            assert molecule.GetNumHeavyAtoms() == int(metadata['heavy_atoms'])
            assert molecule.GetNumConformers() == record['n_poses'] == 10
            energies = [float(l.split()[3]) for l in text.splitlines() if l.startswith('REMARK VINA RESULT:')]
            assert len(energies) == 10 and energies[0] == record['best_affinity_kcal_mol']
            output_key = Chem.MolToInchiKey(molecule)
            with (RESULT / (ligand_id + '_poses.sdf')).open('w', encoding='utf-8') as handle:
                with Chem.SDWriter(handle) as writer:
                    for cid in range(molecule.GetNumConformers()):
                        molecule.SetProp('vina_score_kcal_mol', str(energies[cid]))
                        writer.write(molecule, confId=cid)
            xyz = molecule.GetConformer(0).GetPositions()
            cofactor_xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in receptor_lines
                                     if l.startswith('HETATM') and l[17:20].strip() == target['retain_residue'] and l[76:78].strip() != 'H'])
            cofactor_min_distance = float(np.linalg.norm(xyz[:, None, :] - cofactor_xyz[None, :, :], axis=2).min())
            nearest = {}
            for residue, coordinate in protein_atoms:
                distance = float(np.linalg.norm(xyz - coordinate, axis=1).min())
                nearest[residue] = min(nearest.get(residue, 999), distance)
            contacts = ';'.join(f'{residue}={distance:.2f}' for residue, distance in sorted(nearest.items(), key=lambda p: p[1]) if distance <= 4.0)
            row = dict(target=gene, ligand_id=ligand_id, role=metadata['role'], cid=metadata['cids'],
                       names=metadata['names'], herbs=metadata['herbs'], input_state=metadata['state'],
                       score_kcal_mol=energies[0], n_poses=10, input_charge=int(metadata['formal_charge']),
                       output_charge=Chem.GetFormalCharge(molecule), output_inchikey=output_key,
                       exact_input_identity_preserved=output_key == metadata['inchikey'],
                       min_distance_to_requested_cofactor_A=round(cofactor_min_distance, 3), contacts_within_4A=contacts)
            table.append(row)
            if metadata['role'] == 'positive_control':
                normalized = uncharger.uncharge(molecule)
                rmsds = [rdMolAlign.CalcRMS(normalized, reference, prbId=cid) for cid in range(normalized.GetNumConformers())]
                self_rmsd = rdMolAlign.CalcRMS(reference, reference)
                assert self_rmsd < 1e-7
                shifted = Chem.Mol(reference)
                conformer = shifted.GetConformer()
                for atom in range(shifted.GetNumAtoms()):
                    conformer.SetAtomPosition(atom, conformer.GetAtomPosition(atom) + np.array([3., 0., 0.]))
                assert abs(rdMolAlign.CalcRMS(shifted, reference) - 3.0) < 1e-6
                controls.append(dict(target=gene, ligand_id=ligand_id, best_score_pose_rmsd_A=rmsds[0],
                                     minimum_rmsd_A=min(rmsds), minimum_rmsd_pose=1 + int(np.argmin(rmsds)), all_pose_rmsd_A=rmsds,
                                     method='Meeko bond-aware export; Uncharger for graph matching; RDKit CalcRMS heavy atoms, symmetry, no alignment',
                                     tests='reference self=0; translated reference=3 A'))
    with (RESULT / 'audited_docking_results.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=table[0].keys(), delimiter='\t')
        writer.writeheader(); writer.writerows(table)
    audit = dict(tasks=23, successful=23, poses=230, controls=controls, receptor_audit=receptor_audit,
                 changed_ligand_states=[r['ligand_id'] for r in table if not r['exact_input_identity_preserved']],
                 output_zip_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT / 'data/raw' / RUN).glob('*outputs.zip')})
    (RESULT / 'output_audit.json').write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding='utf-8')
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    run()


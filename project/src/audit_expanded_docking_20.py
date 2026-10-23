"""Audit reference pose reproduction without alignment and summarize all 20 targets."""
import csv, hashlib, io, json
from pathlib import Path
import numpy as np
from rdkit import Chem
from meeko import PDBQTMolecule, RDKitMolCreate
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-04-expanded-docking-20'
PRO=ROOT/'data/processed'/RUN;OUT=ROOT/'results'/RUN
C=json.loads((ROOT/'configs/expanded_docking_20.json').read_text(encoding='utf-8'))
rows=list(csv.DictReader((OUT/'job_manifest.csv').open(encoding='utf-8-sig')))
input_checks=[]
for gene in C['targets']:
    d=PRO/gene;md=json.loads((d/'receptor_manifest.json').read_text())
    selected=[l for l in (d/'receptor.pdb').read_text().splitlines() if l.startswith(('ATOM','HETATM')) and l[17:20].strip() in md['cofactors']]
    coords={tuple(float(l[a:b]) for a,b in [(30,38),(38,46),(46,54)]) for l in (d/'receptor.pdbqt').read_text().splitlines() if l.startswith(('ATOM','HETATM'))}
    for l in selected:
        if tuple(float(l[a:b]) for a,b in [(30,38),(38,46),(46,54)]) not in coords:raise ValueError(f'Lost cofactor atom {gene}')
    delta=None
    if (d/'reference.sdf').exists():
        with (d/'reference.sdf').open('rb') as handle:m=next(iter(Chem.ForwardSDMolSupplier(handle,removeHs=True)))
        delta=float(np.linalg.norm(m.GetConformer().GetPositions().mean(0)-md['center']))
        if delta>=0.01:raise ValueError(f'Wrong reference instance {gene}:{delta}')
    input_checks.append(dict(gene=gene,cofactor_heavy_atoms_preserved=len(selected),reference_centroid_delta_A=delta))
(OUT/'input_coordinate_audit.json').write_text(json.dumps(input_checks,indent=2),encoding='utf-8')
def write(name,rows):
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def export(path):
    pm=PDBQTMolecule(path.read_text(),skip_typing=True)
    sd,failed=RDKitMolCreate.write_sd_string(pm)
    if failed:raise ValueError(f'Failed molecule exports {failed}')
    path.with_suffix('.sdf').write_text(sd,encoding='utf-8')
    return [m for m in Chem.ForwardSDMolSupplier(io.BytesIO(sd.encode()),removeHs=True) if m]
aud=[]
for r in rows:
    path=OUT/r['gene']/f"{r['job']}.pdbqt"
    if not path.exists():continue
    mols=export(path);r=dict(r);r['poses']=len(mols)
    best=mols[0];xyz=best.GetConformer().GetPositions()
    with (PRO/r['gene']/f"{r['job']}.sdf").open('rb') as handle:
        input_mol=next(iter(Chem.ForwardSDMolSupplier(handle,removeHs=True)))
    r['input_identity_preserved']=Chem.MolToSmiles(best,isomericSmiles=True)==Chem.MolToSmiles(input_mol,isomericSmiles=True)
    if not r['input_identity_preserved']:raise ValueError(f'Output identity mismatch {r["gene"]}:{r["job"]}')
    receptor=PRO/r['gene']/'receptor.pdb'
    atoms=[l for l in receptor.read_text().splitlines() if l.startswith(('ATOM','HETATM'))]
    coords=np.array([[float(l[a:b]) for a,b in [(30,38),(38,46),(46,54)]] for l in atoms])
    distance=np.linalg.norm(xyz[:,None,:]-coords[None,:,:],axis=2)
    r['closest_receptor_heavy_A']=round(float(distance.min()),3)
    r['severe_pairs_below_1A']=int((distance<1).sum())
    contacts={}
    for i,l in enumerate(atoms):
        val=float(distance[:,i].min());key=f'{l[17:20].strip()}:{l[21]}:{l[22:26].strip()}'
        if val<=4:contacts[key]=min(contacts.get(key,999),val)
    r['contacts_within_4A']=';'.join(f'{k}={v:.2f}' for k,v in sorted(contacts.items(),key=lambda item:item[1]))
    md=json.loads((PRO/r['gene']/'receptor_manifest.json').read_text())
    co=np.array([[float(l[a:b]) for a,b in [(30,38),(38,46),(46,54)]] for l in atoms if l[17:20].strip() in md['cofactors']])
    r['closest_cofactor_A']=round(float(np.linalg.norm(xyz[:,None,:]-co[None,:,:],axis=2).min()),3) if co.size else ''
    r['reference_fixed_frame_symmetry_RMSD_A']=''
    r['reference_min_RMSD_all_modes_A']=''
    if r['role']=='reference':
        with (PRO/r['gene']/'reference.sdf').open('rb') as f:ref=next(iter(Chem.ForwardSDMolSupplier(f,removeHs=True)))
        matches=best.GetSubstructMatches(ref,uniquify=False,maxMatches=10000)
        if not matches:raise ValueError(f'Reference graph mismatch {r["gene"]}')
        original=ref.GetConformer().GetPositions()
        r['reference_fixed_frame_symmetry_RMSD_A']=round(min(float(np.sqrt(np.mean(np.sum((xyz[list(m)]-original)**2,axis=1)))) for m in matches),3)
        r['reference_min_RMSD_all_modes_A']=round(min(float(np.sqrt(np.mean(np.sum((pose.GetConformer().GetPositions()[list(m)]-original)**2,axis=1)))) for pose in mols for m in matches),3)
    aud.append(r)
if aud:write('pose_audit.csv',aud)

# Report reference-to-metal distances, including alternative instances.
distances=[]
for pdb,ref,met in [('3HS4','AZM','ZN'),('4NRE','C8E','FE2')]:
    lines=(ROOT/'data/raw'/RUN/f'{pdb}.pdb').read_text().splitlines()
    coords=lambda seq:np.array([[float(l[a:b]) for a,b in [(30,38),(38,46),(46,54)]] for l in seq])
    metal=coords([l for l in lines if l.startswith('HETATM') and l[17:20].strip()==met])
    for seq in sorted({l[22:26].strip() for l in lines if l.startswith('HETATM') and l[17:20].strip()==ref}):
        x=coords([l for l in lines if l.startswith('HETATM') and l[17:20].strip()==ref and l[22:26].strip()==seq])
        distances.append(dict(pdb=pdb,residue=ref,seq=seq,closest_metal_A=round(float(np.linalg.norm(x[:,None,:]-metal[None,:,:],axis=2).min()),3)))
write('reference_site_audit.csv',distances)
print(json.dumps(dict(audited_jobs=len(aud),reference=[{k:r[k] for k in ['gene','reference_fixed_frame_symmetry_RMSD_A']} for r in aud if r['role']=='reference'],sites=distances),ensure_ascii=False))

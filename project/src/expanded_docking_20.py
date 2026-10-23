"""Local exploratory docking with cached structures, explicit scope and provenance."""
import concurrent.futures, csv, hashlib, importlib.util, json, re, subprocess, sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

ROOT=Path(__file__).resolve().parents[1]
C=json.loads((ROOT/'configs/expanded_docking_20.json').read_text(encoding='utf-8'))
RAW=ROOT/'data/raw'/C['run_id']; PRO=ROOT/'data/processed'/C['run_id']; OUT=ROOT/'results'/C['run_id']
for p in (PRO,OUT): p.mkdir(parents=True,exist_ok=True)
spec=importlib.util.spec_from_file_location('prep032',ROOT/'src/stage_b_tcmsp_docking_032.py')
prep=importlib.util.module_from_spec(spec);spec.loader.exec_module(prep)
SOURCE=ROOT/'results/2026-10-04-g1-g5-docking-reuse/E1成分靶点去重清单.csv'
rows=list(csv.DictReader(SOURCE.open(encoding='utf-8-sig')))
assert len(rows)==80 and len({r['靶点'] for r in rows})==20

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write_csv(path, rows):
    if not rows:return
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def ligand(m,path):
    setups=MoleculePreparation().prepare(m)
    assert len(setups)==1
    text,ok,err=PDBQTWriterLegacy.write_string(setups[0])
    if not ok:raise ValueError(err)
    path.write_text(text,encoding='ascii')
    with path.with_suffix('.sdf').open('w',encoding='utf-8') as handle:
        w=Chem.SDWriter(handle);w.write(m);w.close()

def prepare():
    manifest=[];issues=[]
    only=sys.argv[sys.argv.index('--only')+1].split(',') if '--only' in sys.argv else None
    for gene,v in C['targets'].items():
        if only and gene not in only:continue
        pdb,chains,ref,refchain,refseq,keep,scope=v
        try:
            d=PRO/gene;d.mkdir(exist_ok=True)
            lines=(RAW/f'{pdb}.pdb').read_text().splitlines()
            coords=[list(map(float,(l[30:38],l[38:46],l[46:54]))) for l in lines if l.startswith('HETATM') and l[17:20].strip()==ref and l[21]==refchain and l[22:26].strip()==refseq and l[76:78].strip()!='H']
            if not coords:raise ValueError('No exact reference residue coordinates')
            center=np.array(coords).mean(0)
            selected={}
            for l in lines:
                if not l.startswith(('ATOM  ','HETATM')) or l[21] not in chains:continue
                if l.startswith('HETATM') and l[17:20].strip() not in keep:continue
                if l[76:78].strip()=='H':continue
                key=(l[21],l[22:27],l[17:20],l[12:16]);occ=float(l[54:60])
                if key not in selected or occ>selected[key][0]:selected[key]=(occ,l[:16]+' '+l[17:])
            rp=d/'receptor.pdb';rq=d/'receptor.pdbqt'
            rp.write_text('\n'.join(x[1] for x in selected.values())+'\nEND\n',encoding='ascii')
            run=subprocess.run([C['obabel'],'receptor.pdb','-O','receptor.pdbqt','-xr','-p','7.4'],cwd=d,capture_output=True,text=True)
            (d/'receptor_preparation.log').write_text(run.stdout+run.stderr,encoding='utf-8')
            if run.returncode or not rq.exists():raise ValueError('Receptor conversion failed')
            atoms=[l for l in rq.read_text().splitlines() if l.startswith(('ATOM','HETATM'))]
            if not atoms:raise ValueError('Empty receptor')
            # OpenBabel renames HET residues. Restore source labels by exact coordinates,
            # then verify cofactor heavy atoms individually instead of accepting renaming.
            original={tuple(round(float(l[a:b]),3) for a,b in [(30,38),(38,46),(46,54)]):l for _,l in selected.values()}
            converted=[]
            for l in rq.read_text().splitlines():
                if l.startswith(('ATOM','HETATM')):
                    xyz=tuple(round(float(l[a:b]),3) for a,b in [(30,38),(38,46),(46,54)])
                    if xyz in original:l=l[:12]+original[xyz][12:30]+l[30:]
                converted.append(l)
            rq.write_text('\n'.join(converted)+'\n',encoding='ascii')
            atoms=[l for l in converted if l.startswith(('ATOM','HETATM'))]
            retained={l[17:20].strip() for l in atoms}
            if not set(keep)<=retained:raise ValueError(f'Missing cofactors: {set(keep)-retained}')
            entry=json.loads((RAW/f'{pdb}_entry.json').read_text(encoding='utf-8'))
            citation=entry.get('rcsb_primary_citation',{})
            md=dict(gene=gene,pdb=pdb,chains=chains,reference=f'{ref}:{refchain}:{refseq}',center=center.tolist(),size=C['box_size'],cofactors=keep,scope=scope,citation_title=citation.get('title'),doi=citation.get('pdbx_database_id_DOI'),receptor_types=sorted({l.split()[-1] for l in atoms}),pdb_sha256=digest(RAW/f'{pdb}.pdb'),receptor_sha256=digest(rq))
            (d/'receptor_manifest.json').write_text(json.dumps(md,indent=2),encoding='utf-8')
            cfg='\n'.join([f'center_{a} = {val:.4f}' for a,val in zip('xyz',center)]+[f'size_{a} = {C["box_size"]}' for a in 'xyz']+[f'{k} = {C[k]}' for k in ['exhaustiveness','num_modes','energy_range']])
            (d/'vina_config.txt').write_text(cfg+'\n',encoding='ascii')
            # Reference is an instance from the chosen structure, with exact chemical graph.
            if gene not in ('ALOX12','GAPDH'):
                p=RAW/f'{pdb}_{ref}_{refchain}_{refseq}_exact.sdf'
                prep.fetch(f'https://models.rcsb.org/v1/{pdb.lower()}/ligand?auth_asym_id={refchain}&auth_comp_id={ref}&auth_seq_id={refseq}&encoding=sdf',p)
                if gene=='PKM':
                    # ModelServer drops explicit charges, despite N1+/O1- in source PDB.
                    # Restore the two crystal-record charges on exact SER instance.
                    with p.open('rb') as handle:
                        m=next(iter(Chem.ForwardSDMolSupplier(handle,removeHs=False,sanitize=False)))
                    if m.GetAtomWithIdx(0).GetAtomicNum()!=7 or m.GetAtomWithIdx(6).GetAtomicNum()!=8:
                        raise ValueError('Unexpected exact SER atom ordering')
                    m.GetAtomWithIdx(0).SetFormalCharge(1);m.GetAtomWithIdx(6).SetFormalCharge(-1)
                    Chem.SanitizeMol(m);mols=[m]
                else:
                    with p.open('rb') as handle:
                        mols=[m for m in Chem.ForwardSDMolSupplier(handle,removeHs=False) if m is not None]
                m=min(mols,key=lambda m:np.linalg.norm(m.GetConformer().GetPositions()[[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()>1]].mean(0)-center))
                delta=np.linalg.norm(m.GetConformer().GetPositions()[[a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum()>1]].mean(0)-center)
                if delta>0.01:raise ValueError(f'Exact reference centroid mismatch {delta}')
                ligand(Chem.AddHs(m,addCoords=True),d/'reference.pdbqt')
                manifest.append(dict(gene=gene,job='reference',role='reference',cid=ref,name=ref,scope=scope))
            for i,row in enumerate(r for r in rows if r['靶点']==gene):
                job=f'candidate_{i+1:02d}';m=prep.embed(row['SMILES'],row['原化学名称']);ligand(m,d/f'{job}.pdbqt')
                manifest.append(dict(gene=gene,job=job,role='candidate',cid=row['CID'],name=row['原化学名称'],scope=scope))
            print('prepared',gene,flush=True)
        except Exception as e:
            issues.append(dict(gene=gene,error=str(e)));print('PREPARATION_ERROR',gene,str(e),flush=True)
    if only:
        previous=list(csv.DictReader((OUT/'job_manifest.csv').open(encoding='utf-8-sig')))
        manifest=[r for r in previous if r['gene'] not in only]+manifest
    write_csv(OUT/'job_manifest.csv',manifest)
    (OUT/'preparation_issues.json').write_text(json.dumps(issues,indent=2),encoding='utf-8')
    return manifest

def dock(row):
    gene,job=row['gene'],row['job'];d=PRO/gene;o=OUT/gene;o.mkdir(exist_ok=True)
    pose=o/f'{job}.pdbqt';log=o/f'{job}.log'
    args=[C['vina'],'--receptor','receptor.pdbqt','--ligand',f'{job}.pdbqt','--config','vina_config.txt','--cpu',str(C['cpu_per_job']),'--seed',str(C['seed']),'--out',f'{job}_out.pdbqt']
    if not (pose.exists() and 'REMARK VINA RESULT:' in pose.read_text()):
        p=subprocess.run(args,cwd=d,capture_output=True,text=True);log.write_text(p.stdout+p.stderr,encoding='utf-8')
        if p.returncode:return dict(**row,status='failed',best_score='',error=(p.stdout+p.stderr)[-1000:])
        pose.write_bytes((d/f'{job}_out.pdbqt').read_bytes());(d/f'{job}_out.pdbqt').unlink()
    scores=re.findall(r'REMARK VINA RESULT:\s+(-?[\d.]+)',pose.read_text())
    print('docked',gene,job,scores[0] if scores else 'NO_SCORE',flush=True)
    return dict(**row,status='completed' if scores else 'failed',best_score=scores[0] if scores else '',error='')

if __name__=='__main__':
    manifest=prepare() if '--prepare' in sys.argv else list(csv.DictReader((OUT/'job_manifest.csv').open(encoding='utf-8-sig')))
    if '--prepare-only' not in sys.argv:
        result=[]
        with concurrent.futures.ThreadPoolExecutor(max_workers=C['workers']) as pool:
            for row in pool.map(dock,manifest):
                result.append(row);write_csv(OUT/'docking_results.csv',result)
        (OUT/'execution_provenance.json').write_text(json.dumps(dict(config=C,source_sha256=digest(SOURCE),code_sha256=digest(Path(__file__)),jobs=len(result),completed=sum(r['status']=='completed' for r in result)),indent=2),encoding='utf-8')

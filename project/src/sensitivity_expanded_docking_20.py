"""One independent-seed, increased-search run for each target's top candidate."""
import concurrent.futures,csv,io,json,re,subprocess,sys
from pathlib import Path
import numpy as np
from rdkit import Chem
from meeko import PDBQTMolecule,RDKitMolCreate
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-04-expanded-docking-20'
OUT=ROOT/'results'/RUN;PRO=ROOT/'data/processed'/RUN
C=json.loads((ROOT/'configs/expanded_docking_20.json').read_text(encoding='utf-8'))
rows=[]
for r in csv.DictReader((OUT/'job_manifest.csv').open(encoding='utf-8-sig')):
    p=OUT/r['gene']/f"{r['job']}.pdbqt"
    if p.exists():
        match=re.search(r'REMARK VINA RESULT:\s+(-?[\d.]+)',p.read_text())
        if match:r.update(status='completed',best_score=match[1]);rows.append(r)
jobs=[]
for gene in sorted(C['targets']):
    if '--only' in sys.argv and gene not in sys.argv[sys.argv.index('--only')+1].split(','):continue
    if '--exclude' in sys.argv and gene in sys.argv[sys.argv.index('--exclude')+1].split(','):continue
    x=[r for r in rows if r['gene']==gene and r['role']=='candidate' and r['status']=='completed']
    if x:jobs.append(min(x,key=lambda r:float(r['best_score'])))
def run(r):
    gene,job=r['gene'],r['job'];d=PRO/gene;o=OUT/gene
    path=o/f'{job}_sensitivity.pdbqt'
    if not path.exists():
        args=[C['vina'],'--receptor','receptor.pdbqt','--ligand',f'{job}.pdbqt','--config','vina_config.txt','--exhaustiveness',str(C['sensitivity_exhaustiveness']),'--seed',str(C['sensitivity_seed']),'--cpu',str(C['cpu_per_job']),'--out',f'{job}_sensitivity_out.pdbqt']
        p=subprocess.run(args,cwd=d,capture_output=True,text=True)
        (o/f'{job}_sensitivity.log').write_text(p.stdout+p.stderr,encoding='utf-8')
        if p.returncode:raise ValueError(p.stderr)
        path.write_bytes((d/f'{job}_sensitivity_out.pdbqt').read_bytes());(d/f'{job}_sensitivity_out.pdbqt').unlink()
    score=float(re.search(r'REMARK VINA RESULT:\s+(-?[\d.]+)',path.read_text())[1])
    mols=[]
    for p in [o/f'{job}.pdbqt',path]:
        sd,failed=RDKitMolCreate.write_sd_string(PDBQTMolecule(p.read_text(),skip_typing=True))
        if failed:raise ValueError(failed)
        if p==path:p.with_suffix('.sdf').write_text(sd,encoding='utf-8')
        mols.append(next(iter(Chem.ForwardSDMolSupplier(io.BytesIO(sd.encode()),removeHs=True))))
    a,b=mols;matches=b.GetSubstructMatches(a,uniquify=False,maxMatches=10000)
    if not matches:raise ValueError('Sensitivity graph mismatch')
    x=a.GetConformer().GetPositions();y=b.GetConformer().GetPositions()
    rmsd=min(float(np.sqrt(np.mean(np.sum((y[list(m)]-x)**2,axis=1)))) for m in matches)
    result=dict(gene=gene,job=job,cid=r['cid'],initial_score=float(r['best_score']),sensitivity_score=score,score_difference=round(score-float(r['best_score']),3),pose_RMSD_A=round(rmsd,3),seed=C['sensitivity_seed'],exhaustiveness=C['sensitivity_exhaustiveness'])
    print('sensitivity',gene,result['score_difference'],flush=True);return result
results=[]
workers=int(sys.argv[sys.argv.index('--workers')+1]) if '--workers' in sys.argv else C['workers']
output=sys.argv[sys.argv.index('--output')+1] if '--output' in sys.argv else 'sensitivity_top_candidates.csv'
with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
    for r in pool.map(run,jobs):
        results.append(r)
        with (OUT/output).open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(r));w.writeheader();w.writerows(results)

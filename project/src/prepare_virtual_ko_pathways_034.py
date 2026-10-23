"""Freeze audited KO inputs and reuse existing pathway snapshots, without rerunning KO."""
from pathlib import Path
import hashlib, json, subprocess
import pandas as pd
import yaml

ROOT=Path(__file__).resolve().parents[1]
CFG_PATH=ROOT/'configs/virtual_ko_pathways_034.yaml'
CFG=yaml.safe_load(CFG_PATH.read_text(encoding='utf-8'))
OUT=ROOT/CFG['output_dir']
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    jobs=[];assessments=[];files=[CFG_PATH]
    for analysis,rel in CFG['inputs'].items():
        directory=ROOT/rel
        for name in ['artifact_audit.json','completion_audit.json']:
            p=directory/name;a=json.loads(p.read_text(encoding='utf-8'))
            assert a['status'] in {'PASS','COMPLETED_WITH_NONESTIMABLE_TARGETS'}
            assert not a.get('issues',[]) and not a.get('errors',[]) and not a.get('missing_jobs',[])
            files.append(p)
        p=directory/'all_context_summary.tsv';s=pd.read_csv(p,sep='\t');files.append(p)
        assert not s.duplicated(['context','target','seed']).any()
        for r in s.to_dict('records'):
            r['analysis']=analysis
            assessments.append(r)
            if r['status']!='COMPLETED': continue
            p=directory/r['context']/f"seed_{r['seed']}"/f"{r['target']}_diffRegulation.tsv"
            dr=pd.read_csv(p,sep='\t')
            assert dr.gene.is_unique and r['target'] in set(dr.gene)
            assert len(dr)-1==r['n_tested']
            assert ((dr.gene!=r['target'])&(dr['p.adj']<CFG['fdr'])).sum()==r['n_significant']
            files.append(p)
            jobs.append(dict(job_id=f"{analysis}__{r['context']}__{r['target']}__{r['seed']}",
                analysis=analysis,context=r['context'],cell_type=r['cell_type'],state=r['state'],
                target=r['target'],seed=r['seed'],path=p.relative_to(ROOT).as_posix(),sha256=sha(p)))
    assert len(assessments)==198 and len(jobs)==84
    pd.DataFrame(assessments).to_csv(OUT/'assessment_registry.tsv',sep='\t',index=False)
    pd.DataFrame(jobs).to_csv(OUT/'job_registry.tsv',sep='\t',index=False)
    libraries=[]
    for name in CFG['libraries']:
        p=ROOT/CFG['library_dir']/f'{name}.gmt';files.append(p)
        rows=[line.rstrip('\n').split('\t') for line in p.read_text(encoding='utf-8').splitlines()]
        assert len({r[0] for r in rows})==len(rows)
        libraries.append(dict(library=name,terms=len(rows),genes=len(set(g for r in rows for g in r[2:])),sha256=sha(p),path=p.relative_to(ROOT).as_posix()))
    pd.DataFrame(libraries).to_csv(OUT/'pathway_library_registry.tsv',sep='\t',index=False)
    hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files}
    p=OUT/'frozen_input_hashes.json'
    if p.exists(): assert json.loads(p.read_text(encoding='utf-8'))==hashes, 'Frozen inputs changed; create a new run or document amendment'
    else: p.write_text(json.dumps(hashes,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    provenance={'run_id':CFG['run_id'],'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'jobs':84,'assessments':198,'reuse':['src/145_enrich_stage_a_030.R','existing three GMT snapshots','fgsea 1.32.4'],
        'resources':'local existing Conda environment; no install/download; no raw-count or KO recomputation',
        'parameters':CFG,'claim':'traceable exploratory analysis; DVC remote not configured; not strict cross-machine replication'}
    (OUT/'provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    trace=ROOT/'workspaces'/CFG['run_id']/'raw/traces'
    trace.mkdir(parents=True,exist_ok=True)
    with (trace/'events.jsonl').open('a',encoding='utf-8') as f:
        f.write(json.dumps({'event':'INPUTS_FROZEN','jobs':84,'assessments':198,'input_manifest_sha256':sha(p)})+'\n')
    print('Frozen 84 valid KO jobs, 198 assessments and three pathway libraries.')
if __name__=='__main__':main()

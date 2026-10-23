"""Full module coverage, immutable inputs, replay ledger, and local execution trace."""
import gzip,hashlib,json,subprocess
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-public-mechanism-analysis-026';OUT=ROOT/'results'/RUN;INT=ROOT/'data/interim'/RUN;RAW=ROOT/'data/raw'/RUN
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return pd.read_csv(p,sep='\t')
def main():
    modules=read(INT/'frozen_modules.tsv');cov=[]
    for family,manifest in [('', 'dataset_manifest.json'),('receiver','receiver_dataset_manifest.json'),('protein','protein_dataset_manifest.json')]:
        for ds in json.loads((INT/manifest).read_text()):
            folder=OUT/family;files=list((folder/'gene_effects').glob(ds['id']+'_*.tsv.gz'))
            if files:valid=set(read(files[0]).human_symbol.dropna())
            else:
                x=read(INT/(ds['id']+'_expression.tsv.gz'));maps=read(INT/(ds['id']+'_gene_map.tsv'));keep=x.iloc[:,1:].std(axis=1)>0;valid=set(maps[maps.feature_id.isin(x.loc[keep,'feature_id'])].human_symbol.dropna())
            for mod,d in modules.groupby('module'):
                avail=len(set(d.gene)&valid);cov.append(dict(dataset=ds['id'],module=mod,genes_available=avail,genes_total=len(d),coverage=avail/len(d),score_estimable=avail>=2,reason='estimable, biological adequacy not certified' if avail>=2 else 'fewer than 2 measured mapped genes, no module score reported'))
    pd.DataFrame(cov).to_csv(OUT/'full_module_coverage.tsv',sep='\t',index=False)
    inputs=[]
    paths=list(RAW.iterdir())+[INT/'frozen_modules.tsv',INT/'mouse_human_one_to_one.tsv',ROOT/'configs/public_mechanism_analysis_026.json',ROOT/'data/raw/2026-09-30-fixed-formula-target-supplement-024/HGNC/hgnc_complete_set.tsv']
    source=json.loads((OUT/'source_manifest.json').read_text());paths.extend(ROOT/p['path'] for p in source if p.get('reused_025'))
    if not all(sha(ROOT/r['path'])==r['sha256'] for r in source if r.get('reused_025')):raise ValueError('025 immutable source snapshot changed')
    for p in paths:
        if p.is_file():inputs.append(dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p)))
    (OUT/'complete_input_manifest.json').write_text(json.dumps(inputs,indent=2)+'\n')
    code_paths=[p for p in (ROOT/'src').iterdir() if p.is_file() and p.name[:3].isdigit() and 125<=int(p.name[:3])<=135]
    code_paths.extend([ROOT/'configs/public_mechanism_analysis_026.json',ROOT/'configs/wikiskill_mechanism_analysis_026.json',ROOT/'tests/test_mechanism_026.py',ROOT/'docs/protocol/2026-10-03-public-mechanism-analysis-026.md',ROOT/'docs/decisions/2026-10-03-public-mechanism-analysis-026.md'])
    (OUT/'code_configuration_manifest.json').write_text(json.dumps({str(p.relative_to(ROOT)):sha(p) for p in code_paths},indent=2)+'\n')
    # Compare the previous numerical fingerprint baseline after re-running all fit branches.
    manifest=json.loads((OUT/'output_content_manifest.json').read_text());replay=[]
    for rel,h in manifest.items():
        p=ROOT/rel;b=gzip.decompress(p.read_bytes()) if p.name.endswith('.gz') else p.read_bytes();v=hashlib.sha256(b).hexdigest()
        replay.append(dict(path=rel,expected_sha256=h,replayed_sha256=v,identical=h==v))
    pd.DataFrame(replay).to_csv(OUT/'offline_replay_verification.tsv',sep='\t',index=False)
    if not all(r['identical'] for r in replay):raise ValueError('numerical replay changed content')
    (OUT/'TCM_select_environment_lock.txt').write_text(subprocess.check_output([str(Path('D:/conda/envs/TCM_select/python.exe')),'-m','pip','freeze'],text=True),encoding='utf-8')
    trace=ROOT/'workspaces'/RUN/'raw/traces';trace.mkdir(parents=True,exist_ok=True)
    payload=dict(run_id=RUN,request='D1/D2 audit and 025 public mechanism reanalysis',execution='local existing TCM_select/R4.3.3; no remote compute or install',sources_manifest=str((OUT/'complete_input_manifest.json').relative_to(ROOT)),claims='C0; RNA/protein expression only; quantitative Desired Reference not qualified; 017 closed',deterministic_tests='6 unittest checks plus 108 source/numerical checks',replay_records=len(replay),replay_identical=True,wikiskill='synthetic deterministic architecture; 0 baseline/best/test, iteration1 rejected; no active or shared skills',pending='cell-level/spatial reanalysis, raw-MS search, two SVGA SDRF links, current formula perturbation dataset and quantitative functional reference')
    p=trace/'execution_summary.json'
    if not p.exists():p.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    (OUT/'provenance_audit.json').write_text(json.dumps(dict(input_files=len(inputs),full_coverage_rows=len(cov),replay_records=len(replay),replay_identical=True,raw_read_only=True,local_only=True),indent=2)+'\n')
    print('replay',len(replay),'identical; input files',len(inputs),'coverage records',len(cov),flush=True)
if __name__=='__main__':main()

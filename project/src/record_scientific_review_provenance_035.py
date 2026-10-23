"""Record complete read-only input fingerprints and final review provenance."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import json, subprocess, sys, importlib.metadata
import pandas as pd
import scientific_review_virtual_ko_035 as r

def main():
    assert json.loads((r.OUT/'final_review_status.json').read_text())['review_complete']
    reg=r.table(r.SRC/'job_registry.tsv');files=set()
    files.update(r.ROOT/p for p in reg.path)
    for jid in reg.job_id:
        files.update((r.SRC/'jobs'/jid).glob('*'))
    files.update(r.SRC.glob('*.tsv'));files.update(r.SRC.glob('*.json'))
    for d in (r.SRC/'networks').iterdir():
        if d.is_dir():files.update(d.glob('*.tsv'));files.add(d/'completion.json')
    for d in (r.SRC/'matched').iterdir():
        if d.is_dir():files.update(d.glob('*.tsv'));files.add(d/'completion.json')
    files.update(r.ROOT/p for p in json.loads((r.SRC/'frozen_input_hashes.json').read_text()))
    files.update([r.ROOT/r.CFG['scope_path'],r.ROOT/'docs/protocol/2026-10-05-virtual-ko-verification-035/核查SOP.md'])
    rows=[]
    for i,p in enumerate(sorted(files),1):
        if p.is_file():rows.append(dict(path=str(p.relative_to(r.ROOT)).replace('\\','/'),size=p.stat().st_size,sha256=r.sha(p),role='reviewed_source_asset_read_only'))
        if i%1000==0:print(f'SOURCE_FINGERPRINT {i}/{len(files)}',flush=True)
    r.save(pd.DataFrame(rows),'reviewed_source_artifact_manifest.tsv')
    fp=r.table(r.OUT/'frozen_input_review.tsv')
    for x in fp.itertuples(index=False):assert r.sha(r.ROOT/x.path)==x.current_sha256,x.path
    r.save(pd.DataFrame([dict(name=d.metadata['Name'],version=d.version) for d in importlib.metadata.distributions()]).sort_values('name'),'python_distribution_lock.tsv')
    envs={}
    for name in ['TCM_select','zhang2026-nicotine-copd-r']:
        records=[]
        for p in (Path('D:/conda/envs')/name/'conda-meta').glob('*.json'):
            x=json.loads(p.read_text(encoding='utf-8'));records.append({k:x.get(k) for k in ['name','version','build','subdir']})
        envs[name]=sorted(records,key=lambda x:x['name'])
    r.js('conda_environment_locks.json',envs)
    test=subprocess.run([sys.executable,'-X','utf8','-m','unittest','discover','-s','tests','-p','test_scientific_review_035.py'],cwd=r.ROOT,capture_output=True,text=True,encoding='utf-8')
    (r.OUT/'review_tests.txt').write_text(test.stdout+test.stderr,encoding='utf-8');assert test.returncode==0
    codes=['configs/virtual_ko_scientific_review_035.yaml','src/scientific_review_virtual_ko_035.py','src/review_virtual_ko_identity_035.R',
           'src/finalize_scientific_review_035.py','src/record_scientific_review_provenance_035.py','src/run_final_scientific_review_035.py',
           'tests/test_scientific_review_035.py','docs/decisions/2026-10-11-virtual-ko-final-scientific-review-035.md']
    manifest=pd.DataFrame([dict(path=p,size=(r.ROOT/p).stat().st_size,sha256=r.sha(r.ROOT/p),role='review_code_or_configuration') for p in codes])
    r.save(manifest,'review_code_manifest.tsv')
    r.js('provenance.json',dict(run_id=r.CFG['run_id'],recorded_at=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat(),
        git_head_at_analysis=json.loads((r.OUT/'environment.json').read_text())['git_head'],
        reviewed_source_assets=len(rows),source_assets_unchanged_after_review=True,reviewer='Codex_AI_not_human_signoff',
        local_resources='existing Windows Python/R environments; three local numerical audit processes; no remote computation',
        input_manifest='reviewed_source_artifact_manifest.tsv',code_manifest='review_code_manifest.tsv',
        dependency_installations=0,new_model_training=False,new_networks=False,DVC_remote_configured=False,
        strict_cross_machine_reproducibility=False,temporary_processes_closed=True,formal_source_or_results_deleted=False))
    output=[]
    for p in sorted(r.OUT.glob('*')):
        if p.is_file() and p.name!='output_manifest.tsv':output.append(dict(path=p.name,size=p.stat().st_size,sha256=r.sha(p)))
    r.save(pd.DataFrame(output),'output_manifest.tsv')
    print('PROVENANCE_COMPLETE',len(rows),flush=True)

if __name__=='__main__':main()

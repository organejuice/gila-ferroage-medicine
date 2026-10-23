"""Parallel matching and enrichment using the unchanged audited scientific adapters."""
import argparse
import concurrent.futures
import importlib.util
import json
from pathlib import Path
import subprocess
import threading
import time

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('parallel035',ROOT/'src/parallel_virtual_ko_verification_035.py')
runtime=importlib.util.module_from_spec(spec);spec.loader.exec_module(runtime)
base=runtime.base
LAUNCH_LOCK=threading.Lock()

def process_job(name,args,inherited_pid=0,completion=None):
    directory=base.OUT/'downstream_workers';directory.mkdir(exist_ok=True)
    marker=directory/f'{name}.json'
    if inherited_pid:
        base.write_json(marker,dict(stage=name,pid=inherited_pid,status='INHERITED_RUNNING'))
        while runtime.alive(inherited_pid):time.sleep(5)
        if completion and not completion.exists():raise RuntimeError(f'Inherited task {name} interrupted')
    else:
        with LAUNCH_LOCK:
            while runtime.available_memory()<2.5*1024**3:
                base.write_json(marker,dict(stage=name,status='WAITING_FOR_MEMORY'));time.sleep(10)
            log=base.OUT/'logs'/f'{name}.log'
            stream=log.open('ab')
            process=subprocess.Popen([base.R,*map(str,args)],cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT)
            base.write_json(marker,dict(stage=name,pid=process.pid,status='RUNNING',log=str(log.relative_to(ROOT))))
            time.sleep(30)
        try:code=process.wait()
        finally:stream.close()
        if code:
            base.write_json(marker,dict(stage=name,status='FAILED',exit_code=code))
            raise RuntimeError(f'{name} failed ({code}); inspect log')
    if completion:
        assert json.loads(completion.read_text())['status']=='PASS'
    base.write_json(marker,dict(stage=name,status='COMPLETE'))
    print('DOWNSTREAM_COMPLETE',name,flush=True)

def run(inherited,enrichment_workers=4,early_enrichment=False):
    assert 1 <= enrichment_workers <= 6
    base.freeze()
    config=str(base.CONFIG.relative_to(ROOT))
    for name in ['adapter_deterministic_tests.json','adapter_real_replay.json']:
        audit=json.loads((base.OUT/name).read_text())
        assert audit['status']=='PASS' and audit['adapter_sha256']==base.sha(ROOT/'src/verify_virtual_ko_035.R')
    tasks=base.pd.read_csv(base.OUT/'network_task_registry.tsv',sep='\t')
    assert len(tasks)==130
    # The previous coordinator validated all network signatures before matching.
    # Preserve all frozen inputs and additionally verify every closed task/output.
    for task in tasks.itertuples(index=False):
        marker=base.OUT/'networks'/task.task_id/'completion.json'
        assert json.loads(marker.read_text())['status']=='PASS'
        base.register(*base.task_rows(task))
    base.register(*base.ribosomes())
    manifest=base.OUT/'downstream_parallel_manifest.json'
    if manifest.exists():
        (base.OUT/('downstream_parallel_manifest_before_'+time.strftime('%Y%m%dT%H%M%S')+'.json')).write_bytes(manifest.read_bytes())
    base.write_json(manifest,dict(matching_workers=2,enrichment_workers=enrichment_workers,early_enrichment=early_enrichment,threads_per_process=1,minimum_available_memory_gib=2.5,inherited_tasks=inherited,scientific_parameters_changed=False,wrapper_sha256=base.sha(Path(__file__)),started=time.strftime('%Y-%m-%dT%H:%M:%S%z')))
    groups=list(base.SCOPE.groupby(['analysis','context']).groups)
    assert set(inherited)<=set(f'{a}__{c}' for a,c in groups)
    groups.sort(key=lambda pair:0 if '__'.join(pair) in inherited else 1)
    early_pool=None;early_futures=[]
    if early_enrichment:
        # Keep the registry immutable while early shards consume its current prefix.
        # No matched jobs are appended until all early shards join.
        early_hash=base.sha(base.OUT/'job_registry.tsv')
        early_pool=concurrent.futures.ThreadPoolExecutor(max_workers=enrichment_workers)
        early_futures=[early_pool.submit(process_job,f'enrich_early_{i}_of_{enrichment_workers}',
            ['src/enrich_virtual_ko_verification_035.R',config,i,enrichment_workers]) for i in range(1,enrichment_workers+1)]
    base.state('MATCHING_AND_EARLY_ENRICHMENT' if early_enrichment else 'PARALLEL_MATCHED_CONTROLS',status='RUNNING',matching_workers=2,enrichment_workers=enrichment_workers if early_enrichment else 0)
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures=[]
        for analysis,context in groups:
            key=f'{analysis}__{context}';completion=base.OUT/'matched'/key/'completion.json'
            if completion.exists():
                assert json.loads(completion.read_text())['status']=='PASS';continue
            futures.append(pool.submit(process_job,f'match_{analysis}_{context}',
                ['src/verify_virtual_ko_035.R',config,'match',analysis,context],inherited.get(key,0),completion))
        for future in concurrent.futures.as_completed(futures):future.result()
    if early_pool:
        try:
            for future in concurrent.futures.as_completed(early_futures):future.result()
            assert base.sha(base.OUT/'job_registry.tsv')==early_hash
        finally:early_pool.shutdown(wait=True)
    # Shared job registry is updated sequentially, then held fixed for all shards.
    for analysis,context in groups:
        directory=base.OUT/'matched'/f'{analysis}__{context}'
        statuses=base.pd.read_csv(directory/'KO_status.tsv',sep='\t');rows=[];metadata=[]
        assert (statuses.status=='COMPLETED').all()
        for r in statuses.itertuples(index=False):
            jid=f'match__{analysis}__{context}__{r.target}__{r.seed}'
            rows.append(dict(job_id=jid,analysis=analysis,context=context,cell_type=context.split('_')[0],state=context.split('_')[1],target=r.target,seed=r.seed,path=r.path,sha256=base.sha(ROOT/r.path)))
            metadata.append(dict(job_id=jid,stage='MATCHED',variant='',original_job_id='',replicate=0,task_id='baseline',n_removed=0))
        base.register(rows,metadata)
    registry_hash=base.sha(base.OUT/'job_registry.tsv')
    base.state('PARALLEL_PATHWAY_ENRICHMENT',status='RUNNING',workers=enrichment_workers,registry_sha256=registry_hash)
    with concurrent.futures.ThreadPoolExecutor(max_workers=enrichment_workers) as pool:
        futures=[pool.submit(process_job,f'enrich_parallel_{i}_of_{enrichment_workers}',['src/enrich_virtual_ko_verification_035.R',config,i,enrichment_workers]) for i in range(1,enrichment_workers+1)]
        for future in concurrent.futures.as_completed(futures):future.result()
    assert base.sha(base.OUT/'job_registry.tsv')==registry_hash
    base.state('FINAL_COMPUTATIONAL_AUDIT',status='RUNNING')
    n=base.summarize()
    for p,h in json.loads((base.OUT/'frozen_input_hashes.json').read_text()).items():assert base.sha(ROOT/p)==h
    base.state('COMPUTATIONS_AUDITED_PENDING_SCIENTIFIC_REVIEW',status='REQUIRES_FINAL_REVIEW',enrichment_jobs=n,direction_estimated=False)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--inherited-file',required=True)
    parser.add_argument('--enrichment-workers',type=int,default=4)
    parser.add_argument('--early-enrichment',action='store_true')
    args=parser.parse_args()
    try:run(json.loads(Path(args.inherited_file).read_text(encoding='utf-8')),args.enrichment_workers,args.early_enrichment)
    except Exception as exc:
        base.state('FAILED_REQUIRES_FIX',status='FAILED',error=str(exc));raise

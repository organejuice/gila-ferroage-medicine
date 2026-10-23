"""Run disjoint local single-core R jobs; retain logs and exact process exit states."""
from pathlib import Path
import json, os, subprocess, time
import yaml
ROOT=Path(__file__).resolve().parents[1]
CFG=yaml.safe_load((ROOT/'configs/virtual_ko_pathways_034.yaml').read_text(encoding='utf-8'))
OUT=ROOT/CFG['output_dir']
def main():
    assert os.environ.get('CONDA_PREFIX','').replace('\\','/').lower()=='d:/conda/envs/zhang2026-nicotine-copd-r', 'Run this launcher through conda run to inherit the correct R DLL paths.'
    handles=[];logs=[]
    for worker in range(1,CFG['execution_workers']+1):
        log=(OUT/f'worker_{worker}.log').open('ab');logs.append(log)
        p=subprocess.Popen(['D:/conda/envs/zhang2026-nicotine-copd-r/lib/R/bin/Rscript.exe','src/enrich_virtual_ko_pathways_034.R',str(worker),str(CFG['execution_workers'])],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        handles.append(p)
    state={'status':'RUNNING','pids':[p.pid for p in handles],'workers':len(handles)}
    (OUT/'run_state.json').write_text(json.dumps(state,indent=2)+'\n',encoding='utf-8');print(json.dumps(state),flush=True)
    try:
        previous=-1
        while any(p.poll() is None for p in handles):
            done=len(list((OUT/'jobs').glob('*/completion.json')))
            if done!=previous:print(f'Completed {done}/84 KO enrichment jobs',flush=True);previous=done
            time.sleep(10)
        codes=[p.wait() for p in handles]
        state={'status':'COMPUTATIONS_COMPLETE_PENDING_AUDIT' if not any(codes) else 'FAILED','exit_codes':codes,'expected_jobs':84,'completed_jobs':len(list((OUT/'jobs').glob('*/completion.json')))}
        (OUT/'run_state.json').write_text(json.dumps(state,indent=2)+'\n',encoding='utf-8');print(json.dumps(state),flush=True)
        if any(codes):raise SystemExit(1)
    finally:
        for log in logs:log.close()
if __name__=='__main__':main()

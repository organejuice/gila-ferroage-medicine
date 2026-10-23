"""Wait on the existing Windows R process handle, then audit and summarize once."""
from pathlib import Path
import argparse,ctypes,json,subprocess,sys,time
from ctypes import wintypes
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/2026-10-04-stratified-virtual-ko-033'
TRACE=ROOT/'workspaces/2026-10-04-stratified-virtual-ko-033-execution/raw/traces/lifecycle.jsonl'

def record(state,**fields):
    item={'time':time.strftime('%Y-%m-%dT%H:%M:%S'),'state':state,**fields}
    TRACE.parent.mkdir(parents=True,exist_ok=True)
    with TRACE.open('a',encoding='utf-8') as f:f.write(json.dumps(item,ensure_ascii=False)+'\n')
    (OUT/'lifecycle_state.json').write_text(json.dumps(item,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(item,ensure_ascii=False),flush=True)

def run(command,label):
    record(label)
    with (OUT/f'{label}.log').open('a',encoding='utf-8') as log:
        result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:
        record('FAILED',step=label,exit_code=result.returncode)
        raise SystemExit(result.returncode)

def main():
    global OUT,TRACE
    p=argparse.ArgumentParser();p.add_argument('--pid',type=int,required=True)
    p.add_argument('--config',default=None);p.add_argument('--stage-prefix',default='stratified_virtual_ko_033');a=p.parse_args()
    if a.config:
        import yaml
        cfg=yaml.safe_load((ROOT/a.config).read_text(encoding='utf-8'))
        OUT=ROOT/cfg['output_dir'];TRACE=ROOT/'workspaces'/f"{cfg['run_id']}-execution"/'raw/traces/lifecycle.jsonl'
    cfgarg=[a.config]if a.config else[]
    kernel=ctypes.WinDLL('kernel32',use_last_error=True)
    kernel.OpenProcess.argtypes=[wintypes.DWORD,wintypes.BOOL,wintypes.DWORD];kernel.OpenProcess.restype=wintypes.HANDLE
    kernel.WaitForSingleObject.argtypes=[wintypes.HANDLE,wintypes.DWORD];kernel.WaitForSingleObject.restype=wintypes.DWORD
    kernel.GetExitCodeProcess.argtypes=[wintypes.HANDLE,ctypes.POINTER(wintypes.DWORD)];kernel.GetExitCodeProcess.restype=wintypes.BOOL
    kernel.CloseHandle.argtypes=[wintypes.HANDLE]
    handle=kernel.OpenProcess(0x00100000|0x1000,False,a.pid)
    if not handle:
        record('PROCESS_HANDLE_UNAVAILABLE',pid=a.pid,error=ctypes.get_last_error());raise SystemExit(1)
    record('WAITING_EXISTING_COMPUTATION',pid=a.pid)
    try:
        while True:
            outcome=kernel.WaitForSingleObject(handle,10000)
            if outcome==0:break
            if outcome!=258:
                record('WAIT_FAILED',pid=a.pid,wait_result=outcome);raise SystemExit(1)
        code=wintypes.DWORD()
        if not kernel.GetExitCodeProcess(handle,ctypes.byref(code)):raise OSError(ctypes.get_last_error())
        record('COMPUTATION_TERMINATED',pid=a.pid,exit_code=code.value)
    finally:kernel.CloseHandle(handle)
    state=OUT/'run_state.json'
    if code.value or not state.exists() or json.loads(state.read_text())['status']!='COMPUTATIONS_COMPLETE_PENDING_AUDIT':
        record('COMPUTATION_NOT_COMPLETE');raise SystemExit(1)
    conda=['D:/conda/Scripts/conda.exe','run','--no-capture-output','-n','zhang2026-nicotine-copd-r','Rscript']
    # Verify original recorded input identities BEFORE refreshing implementation hashes.
    run(conda+['src/141_audit_virtual_ko_artifacts.R']+cfgarg,'artifact_audit')
    run(['D:/python/python.exe','src/137_virtual_ko_provenance.py']+cfgarg,'provenance_refresh')
    run(['D:/python/python.exe','src/133_summarize_virtual_ko.py']+cfgarg,'summary')
    audit=json.loads((OUT/'completion_audit.json').read_text(encoding='utf-8'))
    if audit['status']!='COMPLETED_WITH_NONESTIMABLE_TARGETS':
        record('COMPLETION_AUDIT_NOT_PASSED');raise SystemExit(1)
    for stage in [a.stage_prefix,a.stage_prefix+'_audit',a.stage_prefix+'_summary']:
        run(['C:/Users/administer/.local/bin/dvc.exe','commit','-f',stage],f'dvc_commit_{stage}')
    record('READY_FOR_SCIENTIFIC_REVIEW',expected_jobs=audit['expected_jobs'],recorded_jobs=audit['recorded_jobs'])

if __name__=='__main__':main()

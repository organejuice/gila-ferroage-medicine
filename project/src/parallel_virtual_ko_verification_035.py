"""Resource-only parallel continuation; scientific source/config hashes stay frozen."""
import argparse
import concurrent.futures
import ctypes
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('verification035', ROOT/'src/run_virtual_ko_verification_035.py')
base = importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
LOCK = threading.Lock()
MINIMUM_MEMORY_GIB = 4
ctypes.windll.kernel32.OpenProcess.restype = ctypes.c_void_p
ctypes.windll.kernel32.GetExitCodeProcess.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_ulong)]
ctypes.windll.kernel32.CloseHandle.argtypes = [ctypes.c_void_p]

class Memory(ctypes.Structure):
    _fields_ = [('length', ctypes.c_ulong), ('load', ctypes.c_ulong)] + [(n, ctypes.c_ulonglong) for n in ['total_phys','avail_phys','total_page','avail_page','total_virtual','avail_virtual','avail_extended']]

def available_memory():
    m = Memory(); m.length = ctypes.sizeof(m)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)):
        raise OSError('Cannot measure available memory')
    return m.avail_phys

def alive(pid):
    handle = ctypes.windll.kernel32.OpenProcess(0x1000, False, pid)
    if not handle: return False
    code = ctypes.c_ulong()
    ctypes.windll.kernel32.GetExitCodeProcess(handle, ctypes.byref(code))
    ctypes.windll.kernel32.CloseHandle(handle)
    return code.value == 259

def completed(task):
    p = base.OUT/'networks'/task/'completion.json'
    return p.exists() and json.loads(p.read_text())['status'] == 'PASS'

def worker(task, inherited_pid=0):
    marker = base.OUT/'parallel_workers'/f'{task}.json'
    marker.parent.mkdir(exist_ok=True)
    if inherited_pid:
        base.write_json(marker, dict(task=task,pid=inherited_pid,status='INHERITED_RUNNING'))
        while alive(inherited_pid): time.sleep(5)
        if not completed(task): raise RuntimeError(f'Inherited task {task} did not close')
    else:
        with LOCK:
            while available_memory() < MINIMUM_MEMORY_GIB * 1024**3:
                base.write_json(marker, dict(task=task,status='WAITING_FOR_MEMORY'))
                time.sleep(10)
            log = base.OUT/'logs'/f'{task}.log'
            f = log.open('ab')
            p = subprocess.Popen([base.R,'src/verify_virtual_ko_035.R',str(base.CONFIG.relative_to(ROOT)),'network',task],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
            base.write_json(marker, dict(task=task,pid=p.pid,status='RUNNING'))
            # Let initial allocations settle before another worker starts.
            time.sleep(30)
        try:
            code = p.wait()
        finally:
            f.close()
        if code or not completed(task):
            base.write_json(marker, dict(task=task,status='FAILED',exit_code=code))
            raise RuntimeError(f'{task} failed; inspect its log')
    base.write_json(marker, dict(task=task,status='COMPLETE'))
    return task

def execute(inherited_task, inherited_pid, workers, inherited_file=None):
    assert 2 <= workers <= 4
    inherited = json.loads(Path(inherited_file).read_text(encoding='utf-8')) if inherited_file else {}
    if inherited_pid:
        inherited[inherited_task] = inherited_pid
    assert all(isinstance(k,str) and isinstance(v,int) and v>0 for k,v in inherited.items())
    assert len(inherited) <= workers
    base.freeze()
    tasks = base.pd.read_csv(base.OUT/'network_task_registry.tsv',sep='\t').task_id.tolist()
    pending = [t for t in tasks if not completed(t)]
    assert set(inherited) <= set(tasks)
    pending.sort(key=lambda t: 0 if t in inherited else 1)
    base.state('PARALLEL_NETWORK_VERIFICATION',status='RUNNING',workers=workers,inherited_task=inherited_task)
    manifest=base.OUT/'parallel_execution_manifest.json'
    if manifest.exists():
        archive=base.OUT/('parallel_execution_manifest_before_'+time.strftime('%Y%m%dT%H%M%S')+'.json')
        archive.write_bytes(manifest.read_bytes())
    base.write_json(manifest,dict(workers=workers,threads_per_task=1,minimum_available_memory_gib=MINIMUM_MEMORY_GIB,inherited_tasks=inherited,scientific_parameters_changed=False,wrapper_sha256=base.sha(Path(__file__)),started=time.strftime('%Y-%m-%dT%H:%M:%S%z')))
    # Only task-local files are written by workers. Shared registries are written by
    # this coordinator after all network workers join.
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(worker,t,inherited.get(t,0)):t for t in pending}
        for f in concurrent.futures.as_completed(futures):
            t = f.result()
            print('PARALLEL_TASK_COMPLETE',t,flush=True)
    # Existing entry point validates all caches, registers all outputs, then runs
    # matching, enrichment and final audit with their original settings.
    base.run()

if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--inherited-task',required=True)
    parser.add_argument('--inherited-pid',type=int,required=True)
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--inherited-file')
    args=parser.parse_args()
    try: execute(args.inherited_task,args.inherited_pid,args.workers,args.inherited_file)
    except Exception as exc:
        base.state('FAILED_REQUIRES_FIX',status='FAILED',error=str(exc));raise

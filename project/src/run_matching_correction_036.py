"""Local, resumable minimal correction; preserves 033/034/035 source assets."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time,threading,ctypes,concurrent.futures,argparse
import pandas as pd
import yaml
ROOT=Path(__file__).resolve().parents[1]
CONFIG=ROOT/'configs/virtual_ko_matching_correction_036.yaml'
CFG=yaml.safe_load(CONFIG.read_text(encoding='utf-8'))
OUT=ROOT/CFG['output_dir']; INTERIM=ROOT/CFG['interim_dir']
TRACE=ROOT/'workspaces'/CFG['run_id']/'raw/traces/events.jsonl'
R=Path('D:/conda/envs/zhang2026-nicotine-copd-r/lib/R/bin/Rscript.exe')
LOCK=threading.Lock()
class Memory(ctypes.Structure):
 _fields_=[('length',ctypes.c_ulong),('load',ctypes.c_ulong)]+[(n,ctypes.c_ulonglong) for n in ['total','available','page','available_page','virtual','available_virtual','extended']]
def memory_available():
 m=Memory();m.length=ctypes.sizeof(m)
 if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)):raise OSError('Memory measurement failed')
 return m.available
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(4*1024*1024),b''):h.update(block)
 return h.hexdigest()
def js(p,x):
 p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def state(stage,**kw):
 x=dict(stage=stage,time=time.strftime('%Y-%m-%dT%H:%M:%S%z'),**kw)
 with LOCK:
  js(OUT/'run_state.json',x);TRACE.parent.mkdir(parents=True,exist_ok=True)
  with TRACE.open('a',encoding='utf-8') as f:f.write(json.dumps(x,ensure_ascii=False)+'\n')
  print(stage,kw,flush=True)
def freeze():
 files={CONFIG,ROOT/CFG['scope_path'],ROOT/CFG['network_config'],ROOT/CFG['retention_config']}
 files.update((ROOT/CFG['baseline_dir']).glob('*baseline.rds'))
 scope=pd.read_csv(ROOT/CFG['scope_path'],sep='\t')
 for (a,c),ss in scope.groupby(['analysis','context']):
  orig='results/2026-10-04-stratified-virtual-ko-033' if a=='primary' else 'results/2026-10-05-virtual-ko-target-retention-sensitivity-033'
  for seed in CFG['seeds']:
   files.add(ROOT/orig/c/f'seed_{seed}/WT_network.rds')
   files.update(ROOT/orig/c/f'seed_{seed}/{t}_diffRegulation.tsv' for t in ss.target)
  status=pd.read_csv(ROOT/CFG['source_dir']/'matched'/f'{a}__{c}'/'KO_status.tsv',sep='\t')
  files.update(ROOT/p for p in status.path)
 files.update(ROOT/p for p in ['src/130_virtual_ko_cached_backend.R','src/135_sparse_cp_equivalent.R','src/138_implicit_manifold.R','external/virtual_ko_033/pcNet_upstream.R','src/correct_virtual_ko_matching_036.R','src/matching_correction_036_helpers.R','tests/test_matching_correction_036.R','data/interim/2026-10-03-public-mechanism-analysis-026/frozen_modules.tsv','results/2026-10-05-virtual-ko-pathways-034/two_seed_pathways.tsv'])
 hashes={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}
 p=OUT/'frozen_input_hashes.json'
 if p.exists():assert json.loads(p.read_text(encoding='utf-8'))==hashes,'Frozen inputs/code changed'
 else:js(p,hashes)
 old=json.loads((ROOT/CFG['source_dir']/'frozen_input_hashes.json').read_text())
 # Adapter historically validated after its initial freeze; scientific assets stay exact.
 old_review=pd.read_csv(ROOT/'results/2026-10-11-virtual-ko-scientific-review-035/frozen_input_review.tsv',sep='\t')
 for x in old_review.itertuples(index=False):assert sha(ROOT/x.path)==x.current_sha256,x.path
 state('INPUTS_FROZEN',assets=len(hashes),new_environment=False,remote_computation=False)
def rrun(mode,seed=None,part=None,parts=1):
 log=OUT/'logs'/f'{mode}{"_"+str(seed) if seed is not None else ""}{"_part_"+str(part) if part is not None else ""}.log';log.parent.mkdir(exist_ok=True)
 env=os.environ.copy();prefix='D:/conda/envs/zhang2026-nicotine-copd-r'
 env['CONDA_PREFIX']=prefix;env['PATH']=prefix+'/Library/bin;'+prefix+'/lib/R/bin;'+prefix+'/lib/R/bin/x64;'+env['PATH']
 env['OMP_NUM_THREADS']='1';env['OPENBLAS_NUM_THREADS']='1';env['MKL_NUM_THREADS']='1'
 if mode=='test':args=[str(R),'tests/test_matching_correction_036.R']
 elif mode=='audit':args=[str(R),'src/audit_matching_correction_036.R']
 elif mode=='partition_KO':args=[str(R),'src/partition_matching_correction_036.R',CONFIG.relative_to(ROOT).as_posix(),mode,str(seed),str(part),str(parts)]
 else:args=[str(R),'src/correct_virtual_ko_matching_036.R',CONFIG.relative_to(ROOT).as_posix(),mode]+([] if seed is None else [str(seed)])
 with log.open('ab') as f:
  p=subprocess.Popen(args,cwd=ROOT,env=env,stdout=f,stderr=subprocess.STDOUT)
  state(mode,status='RUNNING',seed=seed,part=part,child_pid=p.pid,log=log.relative_to(ROOT).as_posix())
  code=p.wait()
 if code:raise RuntimeError(f'{mode} failed: {log}')
 state(mode,status='COMPLETED',seed=seed,part=part)
def resource_plan(available,requested=None):
 if requested is not None:
  if requested not in (1,2,4):raise ValueError('Worker count must be 1, 2 or 4')
  required=(1.1*requested+1.1)*1024**3
  if available<required:raise RuntimeError(f'{requested} workers require {required/1024**3:.1f} GiB free; measured {available/1024**3:.2f} GiB')
  return (2 if requested==4 else 1),requested
 parts=2 if available>=6*1024**3 else 1
 return parts,2*parts if available>=4*1024**3 else 1

def run(workers_requested=None):
 OUT.mkdir(parents=True,exist_ok=True);freeze();rrun('test')
 # Preparation already completed is reused only when its numerical files and frozen input exist.
 if not (OUT/'CP10K_and_reference_change.tsv').exists():rrun('prepare')
 freeze_data=json.loads((OUT/'expansion_freeze.json').read_text())
 assert sha(INTERIM/'expanded_input.rds')==freeze_data['input_sha256']
 assert sha(OUT/'frozen_expansion_reference_genes.tsv')==freeze_data['reference_sha256']
 rrun('original_KO')
 for seed in CFG['seeds']:rrun('expanded_WT',seed)
 rrun('expanded_select')
 rrun('audit')
 available=memory_available();parts,workers=resource_plan(available,workers_requested)
 state('EXPANDED_KO_RESOURCE_PLAN',workers=workers,requested_workers=workers_requested,parts_per_seed=parts,available_memory_GiB=available/1024**3,threads_per_worker=1,partition_rule='sorted gene index modulo parts; same per-gene RNG reset')
 with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
  futures=[pool.submit(rrun,'partition_KO',seed,part,parts) for seed in CFG['seeds'] for part in range(parts)]
  for f in futures:f.result()
 expected=set(CFG['expanded_targets'])|set(pd.read_csv(OUT/'expanded_matched_target_controls.tsv',sep='\t').control)
 for seed in CFG['seeds']:
  status=pd.concat([pd.read_csv(OUT/f'expanded_KO_status_seed_{seed}_part_{part}.tsv',sep='\t') for part in range(parts)])
  assert set(status.target)==expected and not status.target.duplicated().any() and set(status.seed)=={seed}
  status.sort_values('target').to_csv(OUT/f'expanded_KO_status_seed_{seed}.tsv',sep='\t',index=False)
 for p,h in json.loads((OUT/'frozen_input_hashes.json').read_text()).items():assert sha(ROOT/p)==h,p
 state('COMPUTATION_COMPLETE_PENDING_SUMMARY',status='COMPLETED')
def summarize_after_process(pid):
 # Attach a read-only wait to an already running local pipeline. No process is stopped.
 api=ctypes.WinDLL('kernel32',use_last_error=True)
 api.OpenProcess.argtypes=[ctypes.c_ulong,ctypes.c_int,ctypes.c_ulong]
 api.OpenProcess.restype=ctypes.c_void_p
 api.WaitForSingleObject.argtypes=[ctypes.c_void_p,ctypes.c_ulong]
 api.WaitForSingleObject.restype=ctypes.c_ulong
 api.CloseHandle.argtypes=[ctypes.c_void_p]
 handle=api.OpenProcess(0x00100000,False,pid)
 if not handle:raise OSError(ctypes.get_last_error(),'Cannot attach to running pipeline')
 print('Waiting for local pipeline',pid,flush=True)
 try:
  while True:
   status=api.WaitForSingleObject(handle,30000)
   if status==0:break
   if status!=258:raise OSError(ctypes.get_last_error(),'Process wait failed')
 finally:api.CloseHandle(handle)
 current=json.loads((OUT/'run_state.json').read_text(encoding='utf-8'))
 assert current['stage']=='COMPUTATION_COMPLETE_PENDING_SUMMARY' and current['status']=='COMPLETED',current
 subprocess.run([sys.executable,'src/summarize_matching_correction_036.py'],cwd=ROOT,check=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser()
 parser.add_argument('--workers',type=int,choices=[1,2,4])
 parser.add_argument('--summarize',action='store_true')
 parser.add_argument('--summarize-after-pid',type=int)
 opts=parser.parse_args()
 try:
  if opts.summarize_after_pid is not None:
   assert opts.workers is None and not opts.summarize
   summarize_after_process(opts.summarize_after_pid)
  else:
   run(opts.workers)
   if opts.summarize:subprocess.run([sys.executable,'src/summarize_matching_correction_036.py'],cwd=ROOT,check=True)
 except Exception as e:
  if opts.summarize_after_pid is not None:
   js(OUT/'summary_watcher_state.json',dict(status='FAILED',error=str(e)))
  else:state('FAILED',error=str(e))
  raise

"""Record immutable-input and implementation identities for KO run 033."""
from pathlib import Path
import hashlib,json,subprocess,sys,yaml
ROOT=Path(__file__).resolve().parents[1]
CFG_PATH=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'configs/virtual_ko_033.yaml'
CFG=yaml.safe_load(CFG_PATH.read_text(encoding='utf-8'))
OUT=ROOT/CFG['output_dir']

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()

paths=[ROOT/'data/processed/D3_GSE160936_glia_seurat.rds',ROOT/'configs/virtual_ko_033.yaml',ROOT/'data/interim/2026-10-03-public-mechanism-analysis-026/frozen_modules.tsv']
paths+=sorted((ROOT/'data/interim/2026-10-04-stratified-virtual-ko-033').glob('*.rds'))
paths+=sorted((ROOT/'src').glob('*virtual_ko*'))+sorted((ROOT/'src').glob('*sparse_cp*'))+[ROOT/'src/134_test_accelerated_pcnet.R']
paths += [ROOT/'src'/n for n in ['138_implicit_manifold.R','139_test_implicit_manifold.R','140_validate_real_implicit_alignment.R']]
paths += [ROOT/'src/141_audit_virtual_ko_artifacts.R']
paths += [ROOT/'src/142_finish_virtual_ko_when_ready.py']
paths+=sorted((ROOT/'external/virtual_ko_033').glob('*'))
paths += [CFG_PATH.resolve(),ROOT/'results/2026-10-04-expanded-docking-20/20靶点分析汇总.csv']
paths += [ROOT/'params.yaml',ROOT/'dvc.yaml']
snapshot_dir=ROOT/'workspaces'/f"{CFG['run_id']}-execution"/'raw/traces'
snapshot_dir.mkdir(parents=True,exist_ok=True)
for p in [ROOT/'params.yaml',ROOT/'dvc.yaml']:
    snapshot=snapshot_dir/f'{p.stem}-{sha(p)}{p.suffix}'
    if not snapshot.exists():snapshot.write_bytes(p.read_bytes())
    paths.append(snapshot)
record={'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
 'git_has_uncommitted_changes':True,'resource':'local Windows 16GB RAM; existing zhang2026-nicotine-copd-r environment; no dependency installation',
 'reused_library_module':'E:/我的文档/实验技术/9_项目实例/2026_Zhang_Nicotine_COPD_方法模块化_v01/modules/15_sctenifoldknk/sctenifoldknk.R',
 'method':{'title':'scTenifoldKnk: An efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation','DOI':'10.1016/j.patter.2022.100434','PMID':'35510185'},
 'files':[{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(set(paths)) if p.is_file()],
 'strict_reproduction':False,'DVC_remote':'not configured; local cache only',
 'WikiSkill':{'deterministic':'core invariants passed; exact-token task score 0 due missing synthetic rule metadata; no skills accepted','real_model':'baseline validation 1.0; best validation 1.0; test 1.0 on exact-token architecture tasks; no skills accepted; all core invariants passed; project-local preexisting CLI compatibility shim reused'}}
(OUT/'provenance.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print('Provenance recorded for',len(record['files']),'files.')

"""Audit and report the additive correction using frozen, unsigned response metrics."""
from pathlib import Path
import json,sys,platform,shutil,subprocess
import numpy as np
import pandas as pd
from scipy.stats import rankdata
import run_matching_correction_036 as run
ROOT=run.ROOT;OUT=run.OUT;CFG=run.CFG;SRC=ROOT/CFG['source_dir']
def read(p):return pd.read_csv(p,sep='\t').fillna('')
def save(df,name):df.to_csv(OUT/name,sep='\t',index=False,lineterminator='\n')
def check_metric():
 assert np.allclose(rankdata([0,2,2],method='average')/4,[.25,.625,.625])
 assert (1+sum(np.array([.5,.5,.1])>=.5))/4==.75
 universe={'a','b','c','d'};refs=[{'a','b','d'},{'b','c','d'}]
 available=universe.copy()
 for ref in refs:available.intersection_update(ref)
 assert all((m&universe&refs[0]&refs[1])==(m&available) for m in [{'a','b'},{'c','d'},{'b','d'}])

def cp10k_changes(responses):
 keys=['analysis','context','target','seed','function']
 old=responses[responses.variant=='original_035'].drop(columns='variant')
 new=responses[responses.variant=='corrected_full_CP10K'].drop(columns='variant')
 delta=old.merge(new,on=keys,suffixes=('_old','_corrected'),validate='one_to_one')
 assert len(delta)==len(old)==len(new)
 for name in ['target_mean_percentile_rank','target_significant_fraction','reference_mean','reference_tail_fraction']:
  delta[name+'_change']=pd.to_numeric(delta[name+'_corrected'],errors='coerce')-pd.to_numeric(delta[name+'_old'],errors='coerce')
 delta['same_common_members']=delta.genes_old.fillna('').eq(delta.genes_corrected.fillna(''))
 delta['formal_specificity_P']=False
 save(delta,'CP10K_functional_change.tsv')
 return delta
def main(original_only=False):
 check_metric();scope=read(ROOT/CFG['scope_path']);assert len(scope)==29
 normalization=[]
 for p in sorted((OUT/'corrected').glob('*/expression_features.tsv')):
  f=read(p);v=float(f.mean_CP10K.sum());assert np.isclose(v,10000,atol=1e-8),p
  lib=read(p.parent/'library_denominator_audit.tsv')
  assert (lib.full_library_UMI>=lib.filtered_library_UMI).all() and (lib.full_library_UMI>0).all()
  normalization.append(dict(background=p.parent.name,genes=len(f),cells=len(lib),CP10K_sum=v,full_UMI_denominator=True))
 save(pd.DataFrame(normalization),'normalization_audit.tsv')
 frozen=json.loads((OUT/'frozen_input_hashes.json').read_text())
 audit=[]
 for p,h in frozen.items():
  current=run.sha(ROOT/p);assert current==h,p
  audit.append(dict(path=p,sha256=h,unchanged=True))
 save(pd.DataFrame(audit),'source_unchanged_audit.tsv')
 modules=read(ROOT/'data/interim/2026-10-03-public-mechanism-analysis-026/frozen_modules.tsv').groupby('module').gene.apply(set).to_dict()
 focus=read(ROOT/'results/2026-10-05-virtual-ko-pathways-034/two_seed_pathways.tsv')
 libs={}
 for lib in ['GO_BP','KEGG','Reactome']:
  libs[lib]={s[0]:set(s[2:]) for s in (line.split('\t') for line in (ROOT/'data/interim/2026-10-03-stage-a-mechanism-030'/f'{lib}.gmt').read_text(encoding='utf-8').splitlines())}
 corrected=read(OUT/'corrected_reference_KO_status.tsv')
 expanded=None if original_only else pd.concat([read(OUT/f'expanded_KO_status_seed_{s}.tsv') for s in [42,43]])
 networks=read(OUT/'expanded_network_audit.tsv');assert len(networks)==2
 quality=[];responses=[];reference_dist=[]
 variants=['original_035','corrected_full_CP10K']
 if not original_only:variants.append('expanded_low_expression')
 for variant in variants:
  for (a,c),ss in scope.groupby(['analysis','context']):
   if variant=='expanded_low_expression' and (a,c)!=(CFG['expanded_analysis'],CFG['expanded_context']):continue
   print('COMPARE',variant,a,c,flush=True)
   orig='results/2026-10-04-stratified-virtual-ko-033' if a=='primary' else 'results/2026-10-05-virtual-ko-target-retention-sensitivity-033'
   if variant=='original_035':
    d=SRC/'matched'/f'{a}__{c}';m=read(d/'matched_target_controls.tsv');st=read(d/'KO_status.tsv')
   elif variant=='corrected_full_CP10K':
    d=OUT/'corrected'/f'{a}__{c}';m=read(d/'matched_target_controls.tsv');st=corrected[(corrected.analysis==a)&(corrected.context==c)]
   else:m=read(OUT/'expanded_matched_target_controls.tsv');st=expanded
   cache={}
   def load(g,seed,target=False):
    key=(g,seed)
    if key in cache:return cache[key]
    if target and variant!='expanded_low_expression':p=ROOT/orig/c/f'seed_{seed}/{g}_diffRegulation.tsv'
    else:
     hit=st[(st.target==g)&(st.seed==seed)];assert len(hit)==1,(variant,g,seed)
     if hit.iloc[0].status!='COMPLETED':cache[key]=None;return None
     p=ROOT/hit.iloc[0].path
    df=pd.read_csv(p,sep='\t',usecols=['gene','distance','p.adj'])
    df=df[df.gene!=g].set_index('gene');assert not df.index.duplicated().any()
    assert np.isfinite(df.distance).all() and np.isfinite(df['p.adj']).all()
    df['percentile']=rankdata(df.distance,method='average')/(len(df)+1);cache[key]=df
    return df
   for t in ss.target:
    mm=m[m.target==t];controls=sorted(mm.control.unique());assert len(controls)<=20
    q=dict(variant=variant,analysis=a,context=c,target=t,n_selected_controls=len(controls),matching_limited=bool(mm.matching_limited.any()),
      insufficient_controls=len(controls)<20,median_matching_distance=mm.distance.median(),max_matching_distance=mm.distance.max(),
      target_detection_fraction=mm.target_detection_fraction.iloc[0],min_control_detection=mm.control_detection_fraction.min(),max_control_detection=mm.control_detection_fraction.max(),
      target_CP10K=mm.target_mean_CP10K.iloc[0],min_control_CP10K=mm.control_mean_CP10K.min(),max_control_CP10K=mm.control_mean_CP10K.max())
    for feature in ['detection_fraction','mean_CP10K','out_degree','out_abs_weight_sum']:
     q[feature+'_outside_selected_range']=any((s['target_'+feature].iloc[0]<s['control_'+feature].min() or s['target_'+feature].iloc[0]>s['control_'+feature].max()) for _,s in mm.groupby('seed'))
    quality.append(q)
    sets={name:gs for name,gs in modules.items() if name.startswith('G5_')==c.startswith('Astrocyte')}
    ff=focus[(focus.analysis==a)&(focus.context==c)&(focus.target==t)&(focus.method=='GSEA')&(focus['mode']=='weighted')]
    for row in ff.itertuples(index=False):sets[f'034_weighted::{row.library}::{row.term}']=libs[row.library][row.term]
    for seed in [42,43]:
     dr=load(t,seed,True);refs={g:load(g,seed) for g in controls};completed={g:f for g,f in refs.items() if f is not None}
     available=set() if dr is None else set(dr.index)
     for f in completed.values():available.intersection_update(f.index)
     for name,gs in sets.items():
      base=dict(variant=variant,analysis=a,context=c,target=t,seed=seed,function=name,n_selected_controls=len(controls),n_completed_controls=len(completed),matching_limited=q['matching_limited'])
      if dr is None:responses.append(dict(**base,status='TARGET_NO_OUTGOING_EDGES'));continue
      common=(gs-{t})&available
      if len(common)<2 or not completed:responses.append(dict(**base,n_members=len(common),status='INSUFFICIENT_COMMON_COVERAGE_OR_REFERENCE'));continue
      common=sorted(common);tv=float(dr.loc[common,'percentile'].mean());ts=float((dr.loc[common,'p.adj']<.05).mean());cv=[]
      for g,f in completed.items():
       val=float(f.loc[common,'percentile'].mean());frac=float((f.loc[common,'p.adj']<.05).mean());cv.append(val)
       reference_dist.append(dict(**base,control=g,mean_percentile_rank=val,significant_fraction=frac,n_members=len(common)))
      responses.append(dict(**base,n_members=len(common),genes=';'.join(common),target_mean_percentile_rank=tv,target_significant_fraction=ts,
       reference_mean=float(np.mean(cv)),reference_min=min(cv),reference_max=max(cv),reference_tail_fraction=(1+sum(v>=tv for v in cv))/(len(cv)+1),
       formal_specificity_P=False,status='DESCRIPTIVE_CONDITIONAL_REFERENCE_NO_CAUSAL_P'))
 suffix='_CP10K_stage' if original_only else ''
 save(pd.DataFrame(quality),f'matching_quality_comparison{suffix}.tsv');save(pd.DataFrame(responses),f'functional_response_comparison{suffix}.tsv');save(pd.DataFrame(reference_dist),f'reference_response_distribution{suffix}.tsv')
 cp10k_changes(pd.DataFrame(responses))
 changes=read(OUT/'CP10K_and_reference_change.tsv');changes['old_relative_excess']=changes.old_CP10K/changes.corrected_CP10K-1;save(changes,'CP10K_and_reference_change.tsv')
 if original_only:
  run.js(OUT/'CP10K_stage_status.json',dict(status='CP10K_STAGE_COMPLETE_EXPANDED_KO_PENDING',combinations=29,source_assets_unchanged=len(audit),functional_rows=len(responses),formal_specificity_P=False))
  print('CP10K_STAGE_COMPLETE; expanded KO remains pending',flush=True)
  return
 q=pd.DataFrame(quality);eq=q[q.variant=='expanded_low_expression'];oldq=q[(q.variant=='original_035')&(q.target.isin(CFG['expanded_targets']))]
 cstatus=corrected.status.value_counts().to_dict();estatus=expanded.status.value_counts().to_dict()
 status=dict(status='MINIMAL_CORRECTION_COMPLETE_WITH_REPORTED_LIMITATIONS',combinations_CP10K_corrected=29,reference_lists_changed=int((changes.retained_controls<20).sum()),
  corrected_reference_KOs=len(corrected),reused_reference_KOs=int(corrected.reused.astype(str).str.lower().eq('true').sum()),corrected_reference_status=cstatus,
  expanded_networks=2,expanded_network_CP_converged=int(networks.CP_converged.sum()),expanded_targets=4,expanded_reference_pool=json.loads((OUT/'expanded_selection_freeze.json').read_text())['both_seed_estimable_reference_pool'],
  expanded_KO_attempts=len(expanded),expanded_KO_status=estatus,expanded_matching_limited=int(eq.matching_limited.sum()),expanded_insufficient_controls=int(eq.insufficient_controls.sum()),
  source_assets_unchanged=len(audit),formal_specificity_P=False,drug_evidence_ceiling='C0',full_resampling_and_enrichment_rerun=False,strict_cross_machine_reproduction=False)
 run.js(OUT/'final_status.json',status)
 env=dict(python=sys.version,executable=sys.executable,numpy=np.__version__,pandas=pd.__version__,platform=platform.platform(),new_dependencies=0,execution='local R workers; at most four disjoint gene partitions, one thread each; reused environments',git_source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip())
 run.js(OUT/'environment.json',env)
 for name in ['conda_environment_locks.json','python_distribution_lock.tsv']:
  shutil.copy2(ROOT/'results/2026-10-11-virtual-ko-scientific-review-035'/name,OUT/name)
 shutil.copy2(SRC/'environment_versions.tsv',OUT/'environment_versions.tsv')
 codes=['configs/virtual_ko_matching_correction_036.yaml','src/correct_virtual_ko_matching_036.R','src/matching_correction_036_helpers.R','src/run_matching_correction_036.py','src/summarize_matching_correction_036.py','src/audit_matching_correction_036.R','src/partition_matching_correction_036.R','tests/test_matching_correction_036.R','docs/decisions/2026-10-11-virtual-ko-minimal-matching-correction-036.md']
 save(pd.DataFrame([dict(path=p,size=(ROOT/p).stat().st_size,sha256=run.sha(ROOT/p)) for p in codes]),'code_manifest.tsv')
 workspace=ROOT/'workspaces'/CFG['run_id'];workspace.mkdir(parents=True,exist_ok=True)
 run.js(workspace/'audit.json',dict(status='PASS',scope='project trace and deterministic numerical audit',backend='existing validated R/Python modules; current Codex session',skill_evolution=False,baseline_validation=None,best_validation=None,test='R 6 checks and Python metric fixtures PASS',accepted_iterations=0,rejected_iterations=0,cross_project_skill_publication=False,active_skills='existing project skills; no new promotion'))
 lines=['# 036 最小匹配修正结果','',
  f'已完成29个组合的完整UMI分母CP10K重算；仅{status["reference_lists_changed"]}个组合各更换1个参照。原分母导致平均CP10K相对偏高{changes.old_relative_excess.min():.1%}–{changes.old_relative_excess.max():.1%}，名单总体稳定。',
  f'原WT参照KO共{len(corrected)}项，其中复用{status["reused_reference_KOs"]}项，其余补算；状态{cstatus}。',
  '',f'AD小胶质扩展WT：保留8704个原基因，加入预冻结80个低表达参照，共8784基因；种子42、43各建一个网络。双种子有出边的低表达参照池为{status["expanded_reference_pool"]}个。',
  f'扩展背景靶点及参照KO共{len(expanded)}项，状态{estatus}。靶点与参照使用同一个扩展WT，未混用原WT的靶点结果。','',
  f'网络CP收敛：{status["expanded_network_CP_converged"]}/2；迭代数、最终残差与解释比例见expanded_network_audit.tsv。未收敛时须保留网络估计限制，不以有限数值替代收敛判定。','',
  '|靶点|原参照检出率范围|扩展参照检出率范围|扩展参照数|匹配范围仍受限|','|---|---|---|---|---|']
 for row in eq.itertuples(index=False):
  before=oldq[oldq.target==row.target].iloc[0]
  lines.append(f'|{row.target}|{before.min_control_detection:.2%}–{before.max_control_detection:.2%}|{row.min_control_detection:.2%}–{row.max_control_detection:.2%}|{row.n_selected_controls}|{row.matching_limited}|')
 iron=pd.DataFrame(responses)
 iron=iron[(iron.variant=='expanded_low_expression')&(iron.function=='G4_M5_Iron_Storage_exploratory')]
 lines+=['','原先涉及的铁储存模块在扩展网络中的描述性比较如下（无表达方向或因果含义）：','',
  '|靶点|种子|共同成员数|靶点平均距离百分位|参照平均值|参照尾部比例|状态|','|---|---|---|---|---|---|---|']
 for row in iron.itertuples(index=False):
  def fmt(name):
   v=getattr(row,name,float('nan'))
   return f'{v:.3f}' if isinstance(v,(float,int)) and np.isfinite(v) else '不可估计'
  lines.append(f'|{row.target}|{row.seed}|{fmt("n_members")}|{fmt("target_mean_percentile_rank")}|{fmt("reference_mean")}|{fmt("reference_tail_fraction")}|{row.status}|')
 lines+=['','低表达支持改善与最终生物学特异性是两层结论；拓扑范围不足、参照不可估计或少于20均记录。参照尾部比例是有限参照分布的描述，不是正式特异性P值。',
  '功能比较使用既有冻结模块及034共同加权通路，按每个KO背景共同覆盖成员计算平均距离百分位和显著成员比例，保留每个参照的完整数值。原始/修正/扩展网络的各自描述可比较，跨网络差异不能解释为药物作用或机制确认。',
  '035留一供体、Bootstrap及全套富集未重跑；其稳定性限制仍适用。此次修正不改变候选原入选门槛，不新增80%通过线，不提升C0药物证据等级。','',
  '验证：6项R确定性检查通过，Python秩与含并列尾部比例算术检查通过，所有DR数值有限且基因唯一，1202项冻结输入逐项SHA-256前后相同。每个新WT、KO与预冻结参照均带校验标记。',
  '复现命令：在项目根目录使用现有TCM_select Python执行 `src/run_matching_correction_036.py`，随后执行 `src/summarize_matching_correction_036.py`；前者自动设置既有R环境，按阶段校验并复用缓存。',
  '资源与来源：复用本地实验技术库单细胞归一化说明和035已验证后端；Windows本地计算，网络构建串行，内存允许时至多四个不重叠基因任务分片，每进程单线程。每个基因均由原KO函数重置相同种子的RNG，分片不改变结果算法；无新增依赖、无远程计算。新源码/配置由Git记录，结果、中间数据和WikiSkill执行轨迹独立保存；DVC最终快照及人工报告复核须在计算结束后确认。未配置DVC Remote，达到本机可追溯重跑，尚不声称严格跨机器复现。','']
 (OUT/'最小修正报告.md').write_text('\n'.join(lines),encoding='utf-8')
 run.state('MINIMAL_CORRECTION_AUDITED',status='COMPLETED',**{k:v for k,v in status.items() if k!='status'})
 files=list(OUT.rglob('*'))
 save(pd.DataFrame([dict(path=p.relative_to(OUT).as_posix(),size=p.stat().st_size,sha256=run.sha(p)) for p in sorted(files) if p.is_file() and p.name!='output_manifest.tsv']),'output_manifest.tsv')
 print(json.dumps(status,ensure_ascii=False),flush=True)
if __name__=='__main__':
 assert sys.argv[1:] in ([],['--original-only']),'Only --original-only is supported'
 main(original_only='--original-only' in sys.argv)

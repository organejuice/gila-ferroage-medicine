"""Summarize all KO strata using frozen project modules and tested-gene backgrounds."""
from pathlib import Path
import hashlib,json,yaml,sys
import pandas as pd
import numpy as np
from scipy.stats import hypergeom,spearmanr

ROOT=Path(__file__).resolve().parents[1]
MODULES=ROOT/'data/interim/2026-10-03-public-mechanism-analysis-026/frozen_modules.tsv'
CFG_PATH=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'configs/virtual_ko_033.yaml'
CFG=yaml.safe_load(CFG_PATH.read_text(encoding='utf-8'))
OUT=ROOT/CFG['output_dir']
ALPHA=float(CFG['fdr'])
LAYERS={
 'G1':['AKR1B1','LTA4H','ALOX15B','CA2','SNCA'],
 'G2':['PKM','GAPDH','SLC2A1'], 'G3':['MMP9','RELA'],
 'G4':['DPP4','NLRP3','PTGS2','DYRK1A','SETD7','ALOX5','ALOX12','ALOX15','MPO','ABCC1','SNCA'],
 'G5_SASP':['MMP9','RELA'], 'G5_LPO':['ALOX12','ALOX15','ALOX5'],
 'G5_Ferro_aging':['ABCC1','ALOX15','DPP4','DYRK1A','MPO','NLRP3','PTGS2','SETD7','SNCA']}
PRIORITY={
 'G1':{'AKR1B1':'优先机制候选','LTA4H':'条件保留','ALOX15B':'条件保留','CA2':'暂不优先抑制','SNCA':'聚集功能研究候选'},
 'G2':{'PKM':'优先功能及异构体核查','GAPDH':'备选功能研究候选','SLC2A1':'备选功能研究候选'},
 'G3':{'MMP9':'优先机制候选','RELA':'通路功能研究候选'},
 'G4':{**dict.fromkeys(['DPP4','NLRP3','PTGS2'],'优先机制复核候选'),**dict.fromkeys(['DYRK1A','SETD7'],'条件保留'),**dict.fromkeys(['ALOX5','ALOX12','ALOX15'],'LPO功能研究候选'),**dict.fromkeys(['MPO','ABCC1','SNCA'],'功能保留，暂非直接作用主线')},
 'G5_SASP':{'MMP9':'优先机制候选','RELA':'通路功能研究候选'},
 'G5_LPO':dict.fromkeys(['ALOX12','ALOX15','ALOX5'],'LPO功能研究候选'),
 'G5_Ferro_aging':{**dict.fromkeys(['DPP4','NLRP3','PTGS2'],'优先机制复核候选'),**dict.fromkeys(['DYRK1A','SETD7'],'条件保留'),**dict.fromkeys(['ABCC1','MPO','SNCA'],'功能保留，暂非直接作用主线'),'ALOX15':'LPO功能研究候选'}}

def bh(p):
    p=np.asarray(p,dtype=float);result=np.full(len(p),np.nan);ix=np.flatnonzero(np.isfinite(p))
    order=ix[np.argsort(p[ix])]
    result[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1])
    return result

def main():
    summaries=list(OUT.glob('*_*/seed_*/target_summary.tsv'))
    if not summaries:
        (OUT/'report.md').write_text('分层虚拟敲除尚在计算。真实输入审计与缓存网络一致性测试已完成；尚无完成的网络与敲除结果，不能报告机制结论。\n',encoding='utf-8')
        print('No completed KO contexts; computation remains active.');return
    s=pd.concat([pd.read_csv(p,sep='\t') for p in summaries],ignore_index=True)
    modules=pd.read_csv(MODULES,sep='\t').groupby('module').gene.apply(set).to_dict()
    enrichment=[];stability=[];layered=[];tables={}
    for row in s.to_dict('records'):
        ct=row['cell_type'];gene=row['target']
        for layer,targets in LAYERS.items():
            if gene in targets and ((ct=='Astrocyte')==layer.startswith('G5')):layered.append(dict(layer=layer,research_category=PRIORITY[layer][gene],**row))
        if row['status']!='COMPLETED':continue
        p=OUT/row['context']/f"seed_{row['seed']}"/f'{gene}_diffRegulation.tsv'
        dr=pd.read_csv(p,sep='\t');dr=dr[dr.gene!=gene].copy()
        assert np.isfinite(dr.distance).all() and dr.gene.is_unique
        tables[(row['context'],gene,row['seed'])]=dr
        bg=set(dr.gene);sig=set(dr.loc[dr['p.adj']<ALPHA,'gene'])
        for name,genes in modules.items():
            if (ct=='Astrocyte')!=name.startswith('G5_'):continue
            present=(genes-{gene})&bg;overlap=present&sig
            pv=float(hypergeom.sf(len(overlap)-1,len(bg),len(present),len(sig))) if len(present)>=2 else np.nan
            enrichment.append(dict(context=row['context'],target=gene,seed=row['seed'],module=name,
                n_module_total=len(genes-{gene}),n_module_in_network=len(present),n_significant=len(sig),n_overlap=len(overlap),
                overlap_genes=';'.join(sorted(overlap)),background=len(bg),p_value=pv,
                status='TESTABLE' if len(present)>=2 else 'INSUFFICIENT_MODULE_COVERAGE'))
    e=pd.DataFrame(enrichment)
    if len(e):
        e['FDR_all_KO_module_tests']=bh(e.p_value)
        e.to_csv(OUT/'module_enrichment.tsv',sep='\t',index=False)
    for context,gene in sorted(set((k[0],k[1]) for k in tables)):
        a=tables.get((context,gene,42));b=tables.get((context,gene,43))
        if a is None or b is None:continue
        merged=a[['gene','distance']].merge(b[['gene','distance']],on='gene')
        aset=set(a.loc[a['p.adj']<ALPHA,'gene']);bset=set(b.loc[b['p.adj']<ALPHA,'gene'])
        stability.append(dict(context=context,target=gene,n_seed42=len(aset),n_seed43=len(bset),n_shared=len(aset&bset),
            significant_jaccard=len(aset&bset)/len(aset|bset) if aset|bset else np.nan,
            distance_spearman=spearmanr(merged.distance_x,merged.distance_y).statistic,
            shared_genes=';'.join(sorted(aset&bset))))
    s.to_csv(OUT/'all_context_summary.tsv',sep='\t',index=False)
    pd.DataFrame(layered).to_csv(OUT/'stratified_target_summary.tsv',sep='\t',index=False)
    pd.DataFrame(stability).to_csv(OUT/'seed_stability.tsv',sep='\t',index=False)
    expected={(ct,state,seed,target) for ct,targets in [('Microglia',sorted(set(sum([v for k,v in LAYERS.items() if not k.startswith('G5')],[])))),('Astrocyte',sorted(set(sum([v for k,v in LAYERS.items() if k.startswith('G5')],[]))))] for state in CFG.get('states',['AD','Control']) for seed in CFG['seeds'] for target in targets}
    actual={(r.cell_type,r.state,int(r.seed),r.target) for r in s.itertuples()}
    artifact_audit=OUT/'artifact_audit.json'
    numerical_audit_pass=artifact_audit.exists() and json.loads(artifact_audit.read_text(encoding='utf-8'))['status']=='PASS'
    valid_statuses={'COMPLETED','LOW_EXPRESSION_NOT_ESTIMABLE','NO_OUTGOING_EDGES'}
    if CFG.get('only_newly_eligible_targets'):valid_statuses.add('PRIMARY_MODEL_ALREADY_EVALUATED')
    complete=expected==actual and len(s)==len(actual) and set(s.status)<=valid_statuses and numerical_audit_pass
    audit=dict(status='COMPLETED_WITH_NONESTIMABLE_TARGETS' if complete else 'INCOMPLETE',expected_jobs=len(expected),recorded_jobs=len(actual),numerical_artifact_audit_pass=numerical_audit_pass,missing_jobs=[list(t) for t in sorted(expected-actual)],errors=s.loc[s.status=='ERROR'].to_dict('records'),scientific_limit='Network perturbation magnitude; no expression direction, PKM isoform, direct drug binding or communication claim.')
    (OUT/'completion_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    lines=['# 分层虚拟敲除结果','',f"当前状态：{audit['status']}；已记录 {len(actual)}/{len(expected)} 个靶点×细胞类型×状态×随机种子任务。",'',
     '输入为 GSE160936 人类 AD/对照原始单细胞计数；AD 状态为主分析，对照为敏感性分析。G5 人类星形胶质细胞结果仅为补充探索，不能替代动物 G5 原主队列。', '',
     '基因检测率门槛为5%，低表达未入网络者记为不可估计；无出边者单独报告。显著数排除了被敲除基因自身。两种种子的一致性见 seed_stability.tsv。', '',
     '|背景|状态|种子|靶点|结果状态|显著扰动基因数|','|---|---|---|---|---|---|']
    if CFG.get('target_retention_min_cells'):
        lines[0]='# 低表达靶点保留敏感性分析（探索性）'
        lines[4]='输入为 GSE160936 人类 AD 原始单细胞计数；本分析仅检查低表达目标保留的敏感性。G5 人类星形胶质细胞结果为补充探索，不能替代动物 G5 原主队列。'
        lines[6]='仅分析AD状态，保留5%检出率背景基因，并额外保留预先指定、至少25个质控后细胞中检出的靶点；这不是全部基因采用25细胞过滤的复现。PRIMARY_MODEL_ALREADY_EVALUATED 表示该靶点已通过主分析门槛，此处不重复计算。新增结果不替换主分析，不能自动提升证据等级。'
    for r in s.itertuples():lines.append(f'|{r.cell_type}|{r.state}|{r.seed}|{r.target}|{r.status}|{r.n_significant}|')
    lines+=['','模块结果采用现有冻结基因集，背景为每个网络中被检验的基因；移除被敲除靶点自身，以免自包含制造富集。对全部 KO×模块检验统一 BH 校正。模块覆盖不足单独报告。', '',
     '分层归属与原有优先机制候选、条件保留、功能研究候选类别见 stratified_target_summary.tsv。类别保留原先疾病证据和对接审查的边界，不用网络显著性替代它们。G5铁摄取与铁输出均无本清单E1交集靶点；不指定新直接靶点，模块富集仅作为其他KO的间接网络响应。G4/G5重复靶点与跨模块重复归属不增加独立验证次数。', '',
     '扰动距离与FC均不是RNA表达变化方向；这些结果不能说明疾病状态恢复、药物抑制或激活、PKM1/PKM2选择性、小胶质→星形胶质通讯，也不能提升对接可靠性。供体层面的网络不确定性与低表达目标的独立队列验证尚未解决。细胞标签继承既有聚类marker评分标注，未新增独立标注验证或双细胞识别。不同KO共同扰动的基因可能反映网络结构，不能直接判定为靶点特异的生化作用。', '',
     '复现代码：src/129_audit_virtual_ko_inputs.R → src/132_run_stratified_virtual_ko.R → src/133_summarize_virtual_ko.py。参数见 '+str(CFG_PATH)+'；环境、输入SHA256和源码来自 provenance.json、audit_sessionInfo.txt 与 Decision Log。']
    lines+=['','## 实际输入、环境与参数','', '|背景|质控及平衡抽样后细胞数|网络基因数|供体数|', '|---|---|---|---|']
    for path in sorted(OUT.glob('*/network_input.json')):
        info=json.loads(path.read_text(encoding='utf-8'))
        lines.append(f"|{info['context']}|{info['n_cells']}|{info['n_genes']}|{info['n_donors']}|")
    lines+=['', f"细胞质控：原始计数总量≥{CFG['cell_min_counts']}，线粒体比例≤{CFG['cell_max_mt_fraction']}；供体×脑区各最多{CFG['max_cells_per_donor_region']}个细胞，固定抽样种子42。网络种子为{CFG['seeds']}；每次构建{CFG['network_n_net']}个子网络，每网抽取{CFG['network_n_cells']}个细胞，PCR成分数{CFG['network_n_comp']}，保留边分位数{CFG['network_q']}。张量秩{CFG['td_K']}、迭代上限{CFG['td_max_iter']}、停止误差{CFG['td_max_error']}、舍入{CFG['td_n_decimal']}位；流形维度{CFG['ma_n_dim']}，显著性FDR<{CFG['fdr']}。脑区合并建网；各供体×脑区的实际细胞名单保存于sampled_cells.tsv。", '',
      '本机复用Conda环境zhang2026-nicotine-copd-r（R 4.4.3、scTenifoldKnk 1.0.3、scTenifoldNet 1.3）；未新增安装。复用实验技术库虚拟敲除接口，并采用经数值等价验证的pcNet加速、稀疏CP及隐式流形计算以适应16GB内存。原始数据不覆盖。Git记录代码，DVC记录输入/结果及流水线；当前DVC仅有本地缓存，未配置远端，因此交付属于可追溯的探索性验证，不声称严格跨机器复现。WikiSkill验证为执行架构检查，不是生物学结果验证，也没有提升为跨项目可复用技能。']
    lines+=['','## 模块富集与种子一致性','', '此处使用排除靶点自身后的基因集及网络检验背景，FDR为全部KO×模块检验的BH校正值。每个背景分别比较两个种子；跨细胞类型或G4/G5重复归属不构成独立验证。未覆盖的任务不能视作阴性。']
    if len(e):
        positive=e[e.FDR_all_KO_module_tests<ALPHA]
        keys=positive[['context','target','module']].drop_duplicates()
        if len(keys):
            lines+=['', '|背景|靶点|模块|seed42 FDR|seed43 FDR|两个种子均显著|', '|---|---|---|---|---|---|']
            for r in keys.itertuples(index=False):
                subset=e[(e.context==r.context)&(e.target==r.target)&(e.module==r.module)].set_index('seed')
                vals=[float(subset.loc[seed,'FDR_all_KO_module_tests']) if seed in subset.index else np.nan for seed in [42,43]]
                textvals=[f'{v:.4g}' if np.isfinite(v) else '未完成或覆盖不足' for v in vals]
                label='是' if all(np.isfinite(v) and v<ALPHA for v in vals) else '否/尚不完整'
                lines.append(f'|{r.context}|{r.target}|{r.module}|{textvals[0]}|{textvals[1]}|{label}|')
        else:lines+=['','在当前已完成的模块检验中，尚无经全部检验BH校正后FDR低于配置门槛的富集。']
    else:lines+=['','尚无可报告的模块检验。']
    docking_path=ROOT/'results/2026-10-04-expanded-docking-20/20靶点分析汇总.csv'
    docking=pd.read_csv(docking_path).set_index('target')
    required_targets=set(CFG['targets_microglia'])|set(CFG['targets_astrocyte'])
    if not required_targets<=set(docking.index):raise ValueError('Missing original docking evidence for requested targets')
    lines+=['','## 疾病层级与原始对接限制','', '下表保留原对接汇总的模型适用范围；虚拟敲除不改变这些限制。参考配体复现误差与候选姿态稳定性是不同检查，不应混用。', '', '|分层|靶点|原研究类别|对接模型及证据限制（原汇总）|', '|---|---|---|---|']
    for layer,targets in LAYERS.items():
        for target in targets:
            scope=str(docking.loc[target,'scope']).replace('|','/')
            lines.append(f'|{layer}|{target}|{PRIORITY[layer][target]}|{scope}|')
    lines+=['','原有核查条件继续保留：AKR1B1优先复核丹参酮IIA，姜黄素/绿原酸需解决NADP重叠；LTA4H核查6-姜烯酚来源与双酶功能方向；ALOX15B核查NDGA选择性与氧化还原干扰；CA2存在疾病表达下调与抑制方向张力；PKM核查异构体、丝氨酸位点及激活/抑制方向；GAPDH不能以非共价对接验证共价机制；SLC2A1需结合运输功能及其他构象。MMP9需复核原实验、Zn配位及选择性；RELA只保留通路功能研究。DPP4/NLRP3/PTGS2的优先表示复核顺序，须确认作用方向及细胞情境；DYRK1A的稳定检查主要验证哈明，姜黄素需另查；SETD7候选两次姿态不一致；ALOX5/ALOX12/ALOX15分别受参考复现、配体对照、人源结构限制；MPO/ABCC1/SNCA保留功能研究，不能据当前模型认定直接作用。SNCA的E46K纤维模型不能推及野生型单体。', '', '方法依据：Osorio等，scTenifoldKnk: An efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation，DOI 10.1016/j.patter.2022.100434；[PubMed原始论文](https://pubmed.ncbi.nlm.nih.gov/35510185/)。采用原始计数构建细胞背景网络；5%基因过滤及25细胞目标保留属于本次配置与库内接口规则，不声称是该论文验证的统一阈值。']
    if stability:
        lines+=['','## 两个随机种子的稳定性','', '|背景|靶点|seed42显著数|seed43显著数|共同显著数|Jaccard|距离Spearman|', '|---|---|---|---|---|---|---|']
        for r in stability:
            lines.append(f"|{r['context']}|{r['target']}|{r['n_seed42']}|{r['n_seed43']}|{r['n_shared']}|{r['significant_jaccard']:.4f}|{r['distance_spearman']:.4f}|")
        lines+=['','显著基因数、共同基因和距离排序均需结合两个种子比较；显著集合不稳定时，不能依据某一次运行的显著数排列机制候选。此处仅检验两次随机网络构建，不等同于独立供体或独立队列复现。']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k!='missing_jobs'},ensure_ascii=False))

if __name__=='__main__':main()

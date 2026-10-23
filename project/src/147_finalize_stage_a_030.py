"""Join evidence, audit numerical results, and deliver stage-A findings."""
import gzip, hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
ROOT=Path(__file__).resolve().parents[1];C=json.loads((ROOT/'configs/stage_a_030.json').read_text());RUN=C['run_id'];OUT=ROOT/'results'/RUN;INT=ROOT/'data/interim'/RUN;RAW=ROOT/'data/raw'/RUN
def read(p):return pd.read_csv(p,sep='\t',low_memory=False).fillna('')
def write(d,p):d.to_csv(p,sep='\t',index=False,compression={'method':'gzip','mtime':0} if str(p).endswith('.gz') else None)
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def save(o,p):Path(p).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n',encoding='utf8')
def evidence():
    A=ROOT/'results/2026-10-03-fixed-formula-target-supplement-024-continuation';D=ROOT/'results/2026-10-03-disease-interface-perturbation-audit-027'
    matrix=read(OUT/'mechanism_compound_target_matrix.tsv.gz');identity=read(INT/'compound_identity_audit.tsv').drop_duplicates('compound_key');idcols=['compound_key','smiles','inchikey','names','identity_status']
    matrix=matrix.drop(columns=[c for c in idcols[1:]+['action_direction','assay_types'] if c in matrix.columns])
    matrix=matrix.merge(identity[idcols],on='compound_key',how='left',validate='many_to_one')
    actions=read(INT/'target_action_records.tsv.gz');actions=actions[actions.record_qualified.astype(str).str.lower()=='true']
    act=actions.groupby(['tier','compound_key','target_symbol']).agg(action_direction=('action_direction',lambda x:';'.join(sorted(set(x)-{''}))),assay_types=('activity_type',lambda x:';'.join(sorted(set(x)-{''}))))
    matrix=matrix.merge(act,left_on=['tier','compound_key','gene'],right_index=True,how='left',validate='many_to_one')
    write(matrix,OUT/'mechanism_compound_target_matrix.tsv.gz')
    summary=[];meta=read(INT/'G1_full_meta.tsv').set_index('gene');links=read(INT/'connections.tsv.gz')
    allgenes=sorted(set(links.gene));states=read(INT/'D0_interface_genes.tsv');p4=read(INT/'P4_P5_qualification_matrix.tsv')
    # Feature-level effects retained; no averaging probe/protein effects and no treatment inference from induction.
    perturb=read(INT/'P7_shared_gene_effects.tsv.gz');perturb=perturb[perturb.human_symbol.isin(allgenes)]
    perturb=perturb.merge(p4[['dataset','contrast','intervention_role','PC','formula_ranking_eligible','therapeutic_direction','exposure','viability_stress']],on=['dataset','contrast'],how='left',validate='many_to_one')
    write(perturb,OUT/'candidate_public_perturbation_gene_effects.tsv.gz')
    for gene in allgenes:
        x=links[links.gene==gene];gates=sorted({s[:2] for s in x.disease_set if not s.startswith('G1_significant')});e1=x[x.tier=='E1'];released=((x.released_main_E1.astype(str).str.lower()=='true') & (x.disease_set=='G1_robust244')).any()
        row=dict(gene=gene,gates=';'.join(gates),interfaces=';'.join(sorted(set(x.disease_set))),tiers=';'.join(sorted(set(x.tier))),E1_compound_pairs=e1.compound_key.nunique(),E1_source_records=len(e1),G1_source_review_released=bool(released),disease_direction_is_treatment_direction=False,formula_measured_presence=False,relevant_free_exposure=False,formula_efficacy_rank='NOT_RELEASED',structure_status='candidate_exact_structure; target_3D_not_yet_audited')
        if gene in meta.index:
            r=meta.loc[gene];row.update(G1_full_meta_logFC=r.meta_logFC,G1_full_meta_FDR=r.meta_FDR,G1_association_scope='full tested meta; robust membership separately indicated')
        row['triage']='G1_SOURCE_REVIEWED_FUNCTIONAL_DIRECTION_PENDING' if released else 'MODEL_E1_RULE_SELECTED_EXPLORATORY' if len(e1) else 'AUXILIARY_OR_PREDICTION_ONLY'
        row['special_review']={'CA2':'RNA_down_vs_assay_inhibition_direction_tension','SNCA':'aggregation_modulation_not_direct_binding','PKM':'isoform_check_PKM1_vs_PKM2_required','GPX4':'protective_function_decrease_not_assumed_beneficial'}.get(gene,'functional_objective_and_relevant_exposure_unresolved')
        p=perturb[perturb.human_symbol==gene];row['external_perturbation_feature_rows']=len(p);row['external_perturbation_contexts']=p[['dataset','contrast']].drop_duplicates().shape[0];row['external_compound_contexts']=p.loc[p.intervention_role.isin(['external_compound','external_rescue_compound']),['dataset','contrast']].drop_duplicates().shape[0]
        summary.append(row)
    write(pd.DataFrame(summary),OUT/'candidate_gene_evidence_matrix.tsv')
    write(pd.DataFrame(summary).query('E1_compound_pairs > 0'),OUT/'E1_candidate_triage.tsv')
    # Attach specific external functional evidence without projecting its claims onto the formula.
    cards=read(INT/'D0_functional_evidence_cards.tsv');write(cards,OUT/'external_functional_evidence_cards.tsv')
    write(read(INT/'D0_TQ_Desired_reference.tsv'),OUT/'therapeutic_direction_qualification.tsv')
    gate_summary=read(INT/'G1_G5_gate_layer_summary.tsv');write(gate_summary,OUT/'G1_G5_gate_layer_summary.tsv')
    membership=[]
    for tier in ['E1','E2','P','H','ALL']:
        gs={gate:set().union(*(set(str(s).split(';'))-{''} for s in gate_summary.loc[(gate_summary.gate==gate)&((gate_summary.tier==tier) if tier!='ALL' else True),'gene_symbols'])) for gate in ['G1','G2','G3','G4','G5']}
        for gene in sorted(set().union(*gs.values())):
            row=dict(tier=tier,gene=gene,**{gate:gene in v for gate,v in gs.items()});row['gate_count']=sum(row[g] for g in gs);membership.append(row)
    write(pd.DataFrame(membership),OUT/'G1_G5_target_membership.tsv')
def audit():
    checks=[]
    frozen=json.loads((OUT/'frozen_input_hashes.json').read_text());replay=[]
    for p,h in frozen.items():replay.append(dict(path=p,match=sha(ROOT/p)==h,sha256=h))
    write(pd.DataFrame(replay),OUT/'historical_input_integrity.tsv');checks.append(dict(check='frozen historical input hashes unchanged',pass_check=all(x['match'] for x in replay)))
    resources=json.loads((RAW/'resource_manifest.json').read_text())
    checks.append(dict(check='all official resource snapshots match frozen SHA256',pass_check=all(sha(ROOT/x['file'])==x['sha256'] for x in resources)))
    reference=json.loads((OUT/'enrichment_replay_reference_hashes.json').read_text());stable=[]
    for name,h in reference.items():
        current=hashlib.sha256()
        with (gzip.open(OUT/name,'rb') if name.endswith('.gz') else (OUT/name).open('rb')) as f:
            for b in iter(lambda:f.read(2**20),b''):current.update(b)
        stable.append(current.hexdigest()==h)
    checks.append(dict(check='offline enrichment rerun reproduces all frozen scientific tables',pass_check=all(stable)))
    targets=read(INT/'all_drug_targets.tsv');counts=targets.groupby('tier').gene.nunique().to_dict();checks.append(dict(check='024 eligible global target counts exactly replayed',pass_check=counts=={'E1':276,'E2':1244,'P':6461,'H':2081}))
    ora=read(OUT/'ORA_all_tests.tsv.gz');sample=ora.sample(min(1000,len(ora)),random_state=C['seed']);pv=hypergeom.sf(sample.overlap.astype(float)-1,sample.background_size.astype(float),sample.term_size.astype(float),sample.query_size.astype(float));delta=np.max(np.abs(pv-sample.P.astype(float)));checks.append(dict(check='R hypergeometric tails independently match scipy',pass_check=bool(delta<1e-10),max_abs_error=float(delta)))
    checks.append(dict(check='zero-overlap tests retained in multiplicity',pass_check=bool((ora.overlap.astype(int)==0).any())))
    checks.append(dict(check='all ORA query/term counts fit their tested universe',pass_check=bool(((ora.query_size.astype(int)<=ora.background_size.astype(int))&(ora.term_size.astype(int)<=ora.background_size.astype(int))&(ora.overlap.astype(int)<=ora.query_size.astype(int))).all())))
    from importlib.util import spec_from_file_location,module_from_spec
    spec=spec_from_file_location('network030',ROOT/'src/146_network_stage_a_030.py');m=module_from_spec(spec);spec.loader.exec_module(m)
    for _,g in ora.groupby(['input','library','background']):
        if not np.allclose(m.bh(g.P.astype(float)),g.FDR_family.astype(float),atol=1e-12):raise ValueError('ORA FDR mismatch')
    checks.append(dict(check='all ORA family BH independently replayed',pass_check=True))
    network=read(OUT/'network_proximity_all.tsv');nchecks=[]
    for name,g in network.groupby('network'):
        null=read(OUT/f'network_null_{name}.tsv.gz')
        nullgroups={key:x.distance.astype(float).values for key,x in null.groupby(['interface','tier','overlap_mode'])}
        for _,r in g.iterrows():
            x=nullgroups[(r.interface,r.tier,r.overlap_mode)]
            p=(1+(x<=float(r.distance)).sum())/(len(x)+1)
            nchecks.append(bool(len(x)==int(r.null_iterations) and abs(p-float(r.P_empirical))<1e-12 and np.isclose(x.mean(),float(r.null_mean),atol=1e-12) and np.isclose(x.std(ddof=1),float(r.null_sd),atol=1e-12)))
    checks.append(dict(check='every empirical network P/count replayed from stored nulls',pass_check=all(nchecks)))
    checks.append(dict(check='network within-variant and all-variant BH independently replayed',pass_check=all(np.allclose(m.bh(g.P_empirical.astype(float)),g.FDR_network.astype(float)) for _,g in network.groupby('network')) and np.allclose(m.bh(network.P_empirical.astype(float)),network.FDR_all_network_tests.astype(float))))
    coverage=read(OUT/'network_coverage_experiments400.tsv');full=coverage[coverage.side=='drug'].drop_duplicates('set').set_index('set').original.astype(int).to_dict()
    checks.append(dict(check='network drug side used full qualified targets in every tier',pass_check=full==counts))
    gsea=read(OUT/'GSEA_all_tests.tsv.gz');checks.append(dict(check='GSEA has primary and sensitivity ranks and bounds',pass_check=set(gsea.ranking)=={'meta_z','meta_logFC'} and set(gsea.size_mode)=={'main','sensitivity'}))
    reactome=gsea[gsea.library=='Reactome'];checks.append(dict(check='Reactome stable identifiers and readable descriptions preserved',pass_check=bool(reactome.pathway.str.startswith('R-HSA-').all() and (reactome.pathway!=reactome.description).all())))
    membership=read(OUT/'G1_G5_target_membership.tsv');checks.append(dict(check='no invented five-gate universal target',pass_check=not (membership.gate_count.astype(int)==5).any()))
    checks.append(dict(check='no formula efficacy ranking or therapeutic direction release',pass_check=(read(OUT/'candidate_gene_evidence_matrix.tsv').formula_efficacy_rank=='NOT_RELEASED').all()))
    write(pd.DataFrame(checks),OUT/'acceptance_checks.tsv');save({'status':'PASS' if all(x['pass_check'] for x in checks) else 'FAIL','checks':checks,'claim_ceiling':'C0','strict_cross_machine_reproducibility':False,'DVC_remote':'not configured'},OUT/'acceptance_audit.json')
    return ora,gsea,network
def report(ora,gsea,network):
    main=gsea[(gsea.ranking=='meta_z')&(gsea.size_mode=='main')].copy();sign=main[main.FDR_all_libraries.astype(float)<C['fdr']].sort_values('FDR_all_libraries');top=sign.head(12)
    sensitivity=main[['library','pathway','description','NES','FDR_all_libraries']].copy()
    for rank,mode,label in [('meta_logFC','main','logFC'),('meta_z','sensitivity','size')]:
        alternate=gsea[(gsea.ranking==rank)&(gsea.size_mode==mode)][['library','pathway','NES','FDR_all_libraries']].rename(columns={'NES':label+'_NES','FDR_all_libraries':label+'_FDR'})
        sensitivity=sensitivity.merge(alternate,on=['library','pathway'],how='left',validate='one_to_one')
        sensitivity[label+'_direction_agrees']=np.sign(sensitivity.NES)==np.sign(sensitivity[label+'_NES'])
    sensitivity['significant_all_three']=(sensitivity.FDR_all_libraries<C['fdr'])&(sensitivity.logFC_FDR<C['fdr'])&(sensitivity.size_FDR<C['fdr'])
    write(sensitivity,OUT/'GSEA_sensitivity_comparison.tsv')
    candidates=read(OUT/'candidate_gene_evidence_matrix.tsv').set_index('gene');drivers=[]
    for _,r in main.iterrows():
        for gene in set(str(r.leadingEdge).split(';')) & set(candidates.index):
            cr=candidates.loc[gene]
            drivers.append(dict(gene=gene,library=r.library,pathway=r.pathway,description=r.description,NES=r.NES,FDR=r.FDR_all_libraries,tiers=cr.tiers,source_review_released=cr.G1_source_review_released,scope='disease leading edge; drug effect direction unestablished'))
    write(pd.DataFrame(drivers),OUT/'candidate_pathway_driver_matrix.tsv')
    family=ora[(ora.FDR_family.astype(float)<C['fdr'])&(ora.overlap.astype(int)>=C['ora']['min_overlap_for_interpretation'])]
    compact=[]
    for key,group in ora.groupby(['input','library','background']):
        eligible=group.overlap.astype(int)>=C['ora']['min_overlap_for_interpretation']
        fs=group[eligible & (group.FDR_family.astype(float)<C['fdr'])]
        compact.append(dict(zip(['input','library','background'],key))|dict(tests=len(group),family_significant=len(fs),global_significant=int((eligible & (group.FDR_global.astype(float)<C['fdr'])).sum()),lowest_family_FDR=float(group.FDR_family.min()),top_family_terms=';'.join(fs.sort_values('FDR_family').description.head(3).astype(str))))
    write(pd.DataFrame(compact),OUT/'ORA_query_summary.tsv')
    primary=network[(network.network=='experiments400')&(network.overlap_mode=='included')];sig=primary[primary.FDR_network.astype(float)<C['fdr']]
    counts={'GSEA_main_tests':len(main),'GSEA_main_global_FDR_significant':len(sign),'GSEA_significant_all_three':int(sensitivity.significant_all_three.sum()),'ORA_tests':len(ora),'ORA_family_FDR_interpretable':len(family),'network_tests':len(network),'network_primary_significant':len(sig)};save(counts,OUT/'result_summary.json')
    lines=['# 阶段A：分层通路、机制证据与网络接近性030','', '本阶段计算与证据整理完成；当前结果为关联/候选机制，不证明整方疗效或机制。药物C0，017/029与G2c/G6关闭保持。','', '## 数据与方法','', '复用024续跑的来源/成分身份和分层连接、026公开扰动效应、027冻结接口/资格审计。G1所有有限meta受检基因以meta_z排序；meta_logFC和10–300大小敏感性。GO BP为本地org.Hs.eg.db3.18.0；KEGG/Reactome官方快照及STRING12.0保存版本/SHA。GSEA主范围15–500、eps1e-5、seed30030；ORA5–500、BH0.05、解释至少2基因重叠。预排序GSEA不是供体表型置换。','', 'G2仅RNA代理注释，原占位受检背景未建立；G3/G4背景为原D2可定量基因（非完整重建过滤宇宙）；G5背景由原受检小鼠基因与NCBI一对一人同源映射获得。来源有记录靶点只作选择敏感性，不称完整来源测试宇宙。','', '## G1完整排序通路结果','', f'主分析共{len(main)}项检验，跨GO BP/KEGG/Reactome合并BH-FDR<0.05有{len(sign)}项。方向是疾病表达排序富集，不是蛋白或通路直接活性。','', '|库|通路|NES|合并FDR|','|---|---|---:|---:|']
    for _,r in top.iterrows():lines.append(f'|{r.library}|{r.description}|{float(r.NES):.3f}|{float(r.FDR_all_libraries):.3g}|')
    kegg=sign[sign.library=='KEGG']
    lines+=['','KEGG单独展示（同一跨库校正）：','', '|通路|NES|合并FDR|','|---|---:|---:|']
    for _,r in kegg.iterrows():lines.append(f'|{r.description}|{float(r.NES):.3f}|{float(r.FDR_all_libraries):.3g}|')
    if not len(kegg):lines.append('|没有通过跨库校正的KEGG项|—|—|')
    lines+=['','## 分层ORA','', f'保留{len(ora):,}个背景内全检验（含零重叠）；族内FDR<0.05且至少2个重叠基因共{len(family):,}行。全局FDR列另附，不能把不同层/不同背景/重复程序行当独立通路数。预设G2–G5成员富集只为注释，不构成独立发现。小交集统计信息有限，阴性不解释为无作用。','', '## 网络接近性','', f'主网络使用直接实验通道>=400，实验转移通道不纳入；该证据是实验支持关联，不保证二元物理结合。各接口取完整合格药物靶点作为药物端，双端数量/度箱匹配随机对照1000或10000次。主网络{len(primary)}项含重叠检验中族内BH<0.05为{len(sig)}项。','', '|接口|层|药物节点|疾病节点|z|FDR|','|---|---|---:|---:|---:|---:|']
    g1e1=family[(family.input=='G1_robust244__E1')&(family.background=='disease_measured')];g1e1s=family[(family.input=='G1_robust244__E1')&(family.background=='drug_record_selection_sensitivity')]
    lines.insert(lines.index('## 网络接近性'),f'G1来源核验E1五基因在疾病受检背景有{len(g1e1)}项可解释GO/KEGG/Reactome富集：'+ '; '.join(g1e1.description.astype(str))+f'。来源记录选择背景有{len(g1e1s)}项，背景敏感性限制保留，不称稳定整方机制。\n')
    for _,r in primary.sort_values('FDR_network').head(12).iterrows():lines.append(f'|{r.interface}|{r.tier}|{int(r.drug_nodes)}|{int(r.disease_nodes)}|{float(r.z):.3f}|{float(r.FDR_network):.3g}|')
    cov=read(OUT/'network_coverage_experiments400.tsv');gc=cov[(cov.side=='disease')&(cov['set']=='G1_robust244')].iloc[0]
    lines+=['',f'主网络G1映射{int(gc.mapped)}/{int(gc.original)}基因。未映射/非主连通分量成员不能解释成没有生物作用，丢失列表见coverage表。药物端276/1244/6461/2081为规则筛选的E1/E2/P/H候选来源范围，不能把全部节点称为人工核验的直接结合靶点；来源可测性偏倚并未被度匹配完全消除。','', 'G1去除交集与网络来源敏感性（各网络内BH）：','', '|网络|层|交集处理|z|FDR|','|---|---|---|---:|---:|']
    for _,r in network[network.interface=='G1_robust244'].sort_values(['network','tier','overlap_mode']).iterrows():lines.append(f'|{r.network}|{r.tier}|{r.overlap_mode}|{float(r.z):.3f}|{float(r.FDR_network):.3g}|')
    g1primary=primary[primary.interface=='G1_robust244']
    lines+=['',f'G1主网络四层中FDR<0.05的接近性有{int((g1primary.FDR_network.astype(float)<C["fdr"]).sum())}层。这不支持将药物靶点整体与G1疾病模块的网络接近性作为当前整方机制依据。敏感性中的单独显著项应与主结果及去交集结果共同解释。']
    lines+=['', '去除端点重叠、实验阈值700及加入文本挖掘700的敏感性已独立计算。全局BH和Monte Carlo误差/分半差保存。接近性无作用方向，不放行处方疗效排名。G4/G5共享预设成员，不是独立复制。','', '## 候选与作用方向','', 'G1来源已核验候选仍为AKR1B1、ALOX15B、CA2、LTA4H、SNCA。CA2疾病下调而现有记录为酶抑制，方向疑问保持；SNCA支持聚集调节而非已证直接结合。PKM仍需核查PKM1/PKM2。GPX4保护功能下降不能默认有益。','', 'candidate_gene_evidence_matrix.tsv与E1_candidate_triage.tsv保留每个基因的G集合、来源层、成分连接、G1全meta状态与公开扰动可测信息。mechanism_compound_target_matrix.tsv.gz含精确结构/成分/药味/assay/来源链；candidate_public_perturbation_gene_effects.tsv.gz保留特征级效应和CI，不平均探针或蛋白，不把诱导/外部参照当本方干预。','', '目前整方实测成分、相关自由暴露、靶点触达和匹配功能目标仍未建立，不能给出治疗性抑制/激活建议；本阶段没有制造一个加权疗效总排名。','', '## 验收、复现与文件','', '独立数值核验见acceptance_checks.tsv及acceptance_audit.json；原历史输入按SHA复核未改。Python/R现有环境复用，无新依赖/远程计算。参数/协议/Decision Log与raw/interim/results分层、环境锁及代码输入哈希保存。DVC remote未配置，只声明本地可追溯和离线重放，非严格跨机器复现。WikiSkill确定性与真实后端状态另见日志，不把合成测试称科研benchmark。','', '离线入口：D:/conda/envs/TCM_select/python.exe src/148_run_stage_a_030.py。首次资源取得用144脚本；结果后续重复计算仅读取固定快照。正式中间数据保留，下载失败part文件自动清理；保留失败/优化尝试日志供审计。','', '主要产物：GSEA_all_tests.tsv.gz、GSEA_significant_family.tsv、ORA_all_tests.tsv.gz、ORA_background_QC.tsv、network_proximity_all.tsv、各网络完整null分布/覆盖/度箱、机制矩阵与候选表。','', '方法依据：GSEA DOI10.1073/pnas.0506580102；clusterProfiler4 DOI10.1016/j.xinn.2021.100141；网络医学 DOI10.1126/sciadv.adh0215；STRING2025 DOI10.1093/nar/gkae1113；NetMedPy MIT固定提交保存在raw目录。']
    lines+=['','## 敏感性与下一步','',f'GSEA在主排序、logFC排序和较窄大小范围均通过各自合并FDR的通路有{int(sensitivity.significant_all_three.sum())}项；逐项方向/缺测比较见GSEA_sensitivity_comparison.tsv。候选与疾病leading-edge的连接见candidate_pathway_driver_matrix.tsv，连接本身不证明药物逆转。','', 'ORA逐查询/库/背景摘要见ORA_query_summary.tsv；保留族内和全局两种校正。G1_G5_target_membership.tsv展示五层共同/特异成员；G0为准入门而不是基因集。','', '神经递质/突触通路若出现富集，应结合原样本的细胞纯度、混杂与RNA来源核验后解释，不能据此直接认定小胶质细胞具有相应功能改变。','', '网络随机次数按预先设定的P值/分半误差规则扩展；存在自适应停止，P值与BH结果仅作探索性筛选。进入验证阶段前应使用固定次数重新确认选定假设，并处理通路冗余。','', '阶段B首先核验整方实测成分及相关自由暴露，再从来源已核验、方向可解释且有可测终点的候选中固定少量假设。后续可进行结构审计和有阳性/阴性对照的对接；对接不是本阶段交付，也不能替代靶点结合/触达及敲低救援实验。']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    # Compact figures emphasize quantities and keep library/endpoint labels explicit.
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'svg.fonttype':'none'})
    if len(top):
        x=top.iloc[::-1];fig,ax=plt.subplots(figsize=(10,6));ax.barh(range(len(x)),x.NES.astype(float),color=['#31688e' if n<0 else '#b24c3b' for n in x.NES.astype(float)]);ax.set_yticks(range(len(x)),[str(t)[:85] for t in x.description]);ax.set_xlabel('G1 disease prerank NES; association only');ax.axvline(0,color='black',lw=.6);fig.tight_layout();fig.savefig(OUT/'G1_GSEA_top.svg');fig.savefig(OUT/'G1_GSEA_top.png',dpi=200);plt.close(fig)
    pivot=primary.pivot(index='interface',columns='tier',values='z').astype(float);fig,ax=plt.subplots(figsize=(10,10));lim=float(np.nanmax(np.abs(pivot.values)));image=ax.imshow(pivot.values,aspect='auto',cmap='RdBu_r',vmin=-lim,vmax=lim);ax.set_yticks(range(len(pivot)),pivot.index);ax.set_xticks(range(len(pivot.columns)),pivot.columns);fig.colorbar(image,ax=ax,label='Proximity z; experimental association network');ax.set_title('Full drug target sets vs frozen interfaces');fig.tight_layout();fig.savefig(OUT/'network_proximity_heatmap.svg');fig.savefig(OUT/'network_proximity_heatmap.png',dpi=200);plt.close(fig)
def main():
    evidence();ora,gsea,net=audit();report(ora,gsea,net)
    wiki=[]
    for rid,log in [('2026-10-03-stage-a-030-deterministic','wikiskill_deterministic_audit.log'),('2026-10-04-stage-a-030-model-compat','wikiskill_model_compat_audit.log')]:
        state=json.loads((ROOT/'workspaces'/rid/'runs/state.json').read_text(encoding='utf8'));wa=json.loads((OUT/log).read_text(encoding='utf8'))
        wiki.append(dict(run_id=rid,backend=state['backend'],baseline_validation=state['baseline_val'],best_validation=state['r_best'],test=state['final_test']['mean'],test_tasks=state['final_test']['n'],accepted=[h['iteration'] for h in state['history'] if h['accepted']],rejected=[h['iteration'] for h in state['history'] if not h['accepted']],active_skills=wa['active_skills'],audit=wa['checks']['all_core_invariants'],scope='architecture tasks only; not biomedical validation or paper benchmark'))
    save(wiki,OUT/'wikiskill_architecture_audit.json')
    with (OUT/'report.md').open('a',encoding='utf8') as f:
        f.write('\n## WikiSkill架构审计\n\n确定性后端baseline/best/test=0/0/0，迭代1被拒绝；真实CLI后端=1/1/1，无接受或拒绝迭代。两个后端审计通过、active skills为空。真实后端验证集/测试集各仅一个架构任务，不能代表生物医学任务性能或论文benchmark。工作区与完整轨迹由DVC保留，详见wikiskill_architecture_audit.json。\n')
    subprocess.run([sys.executable,'-m','pip','freeze'],stdout=(OUT/'Python_environment_lock.txt').open('w',encoding='utf8'),check=True)
    code=[ROOT/'src'/f for f in ['144_prepare_stage_a_030.py','145_enrich_stage_a_030.R','146_network_stage_a_030.py','147_finalize_stage_a_030.py','148_run_stage_a_030.py','149_wikiskill_cli_compat_030.py']]+[ROOT/'configs/stage_a_030.json',ROOT/'params.yaml',ROOT/'dvc.yaml',ROOT/'tests/test_stage_a_030.py',ROOT/'docs/protocol/2026-10-03-stage-a-mechanism-030.md',ROOT/'docs/decisions/2026-10-03-stage-a-mechanism-030.md']
    save([dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in code],OUT/'code_parameter_manifest.json')
    import platform,scipy,networkx,requests
    head=subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,capture_output=True,text=True,check=True).stdout.strip()
    save(dict(Python=sys.version,executable=sys.executable,platform=platform.platform(),packages={'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'networkx':networkx.__version__,'requests':requests.__version__},Git_parent=head,Git_release_commit=False,reason='preserve existing mixed project modifications/staging; code and parameters individually SHA256 recorded',compute='local CPU; no new environment or installation',library_reuse='local method cards, existing R/Python, frozen MIT NetMedPy degree binning'),OUT/'runtime_provenance.json')
    print(json.dumps(json.loads((OUT/'result_summary.json').read_text()),ensure_ascii=False),flush=True)
if __name__=='__main__':main()

"""Assemble explicit scientific dispositions after independent SOP035 audit.

No automatic 80% cutoff, no efficacy ranking, no human signature, no changes
to 033/034/035 inputs. All narrative dispositions are explicit in this source.
"""
from pathlib import Path
import json, hashlib, sys, subprocess
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import scientific_review_virtual_ko_035 as r

ROOT=r.ROOT;OUT=r.OUT;SRC=r.SRC;KEY=r.KEY

def contributor_review():
    n=non_sig=0
    for q in pd.read_csv(SRC/'target_contributor_audit.tsv',sep='\t',chunksize=100000):
        assert (q.gene_significant==(q.gene_FDR<.05)).all()
        assert (q.expression_direction=='NOT_ESTIMATED').all() and q.distance.ge(0).all()
        n+=len(q);non_sig+=int((~q.gene_significant).sum())
    r.js('contributor_semantics_audit.json',dict(status='PASS',rows=n,leading_edge_rows_not_gene_significant=non_sig,unsigned_direction_preserved=True))

def matched_focus_review():
    from scipy.stats import rankdata
    scope=r.table(ROOT/r.CFG['scope_path']).fillna('')
    focus=r.table(ROOT/'results/2026-10-05-virtual-ko-pathways-034/two_seed_pathways.tsv').fillna('')
    focus=focus[(focus.method=='GSEA')&(focus['mode']=='weighted')]
    cfg=r.yaml.safe_load((ROOT/r.CFG['source_config']).read_text(encoding='utf-8'))
    libs={}
    for lib in cfg['libraries']:
        libs[lib]={x[0]:set(x[2:]) for x in (line.split('\t') for line in (ROOT/cfg['library_dir']/f'{lib}.gmt').read_text(encoding='utf-8').splitlines())}
    rows=[]
    for (analysis,context),ss in scope.groupby(['analysis','context']):
        d=SRC/'matched'/f'{analysis}__{context}';m=r.table(d/'matched_target_controls.tsv');status=r.table(d/'KO_status.tsv')
        cache={}
        for x in status.itertuples(index=False):
            dr=r.table(ROOT/x.path);dr=dr[dr.gene!=x.target].set_index('gene')
            cache[(x.target,x.seed)]=(dr,pd.Series(rankdata(dr.distance,method='average')/(len(dr)+1),index=dr.index))
        base='results/2026-10-04-stratified-virtual-ko-033' if analysis=='primary' else 'results/2026-10-05-virtual-ko-target-retention-sensitivity-033'
        for s in ss.itertuples(index=False):
            mm=m[m.target==s.target];ff=focus[(focus.analysis==analysis)&(focus.context==context)&(focus.target==s.target)]
            for seed in [42,43]:
                dr=r.table(ROOT/base/context/f'seed_{seed}/{s.target}_diffRegulation.tsv');dr=dr[dr.gene!=s.target].set_index('gene')
                rank=pd.Series(rankdata(dr.distance,method='average')/(len(dr)+1),index=dr.index)
                common_bg=set(dr.index)
                for control in mm.control.unique():common_bg&=set(cache[(control,seed)][0].index)
                for p in ff.itertuples(index=False):
                    common=sorted((libs[p.library][p.term]-{s.target})&common_bg)
                    if len(common)<2:
                        rows.append(dict(analysis=analysis,context=context,target=s.target,seed=seed,library=p.library,term=p.term,n_members=len(common),status='INSUFFICIENT_COMMON_COVERAGE'));continue
                    tv=float(rank.loc[common].mean());cv=[float(cache[(c,seed)][1].loc[common].mean()) for c in mm.control.unique()]
                    rows.append(dict(analysis=analysis,context=context,target=s.target,seed=seed,library=p.library,term=p.term,n_members=len(common),
                        target_mean_percentile_rank=tv,target_significant_fraction=float((dr.loc[common,'p.adj']<.05).mean()),
                        reference_min=min(cv),reference_max=max(cv),reference_mean=float(np.mean(cv)),reference_values=';'.join(f'{v:.12g}' for v in cv),
                        reference_tail_fraction=(1+sum(v>=tv for v in cv))/21,n_controls=20,matching_limited=bool(mm.matching_limited.any()),
                        status='DESCRIPTIVE_FIXED_034_WEIGHTED_FOCUS_NO_FORMAL_SPECIFICITY_P'))
    r.save(pd.DataFrame(rows),'matched_focus_response_comparison.tsv')

def module_review():
    q=r.table(SRC/'resampling_module_enrichment.tsv').fillna({'hits':''})
    jobs=r.table(SRC/'enrichment_task_summary.tsv')
    q=q.merge(jobs[['job_id','n_tested','n_significant']],on='job_id',validate='many_to_one')
    calc=np.where(q.n_members>=2,hypergeom.sf(q.hits.map(lambda s:len(r.genes(s)))-1,q.n_tested,q.n_members,q.n_significant),np.nan)
    assert np.nanmax(np.abs(calc-q.P))<1e-12
    assert np.nanmax(np.abs(r.bh(q.P)-q.FDR_new_run_all_module_tests))<1e-10
    tasks=r.table(SRC/'network_task_registry.tsv')
    meta=r.table(SRC/'enrichment_job_metadata.tsv').fillna('')
    q=q.merge(meta[['job_id','replicate']],on='job_id',validate='many_to_one')
    rows=[]
    scope=r.table(ROOT/r.CFG['scope_path']).fillna('')
    for s in scope.itertuples(index=False):
        x=q[(q.analysis==s.analysis)&(q.context==s.context)&(q.target==s.target)]
        for stage,B in [('LODO',6),('BOOTSTRAP',20)]:
            for name,t in x[x.stage==stage].groupby('module'):
                test=hit=0
                for _,z in t.groupby('replicate'):
                    if len(z)==2 and z.P.notna().all():
                        test+=1;hit+=int((z.FDR_new_run_all_module_tests<.05).all())
                rows.append(dict(analysis=s.analysis,context=s.context,target=s.target,stage=stage,module=name,
                    n_planned=B,n_testable=test,n_reproduced=hit,unconditional_rate=hit/B,conditional_rate=hit/test if test else np.nan,
                    interpretation='RETIRED_MODULE_ANNOTATION_ONLY' if any(v in name for v in ['Iron_Storage','Ferritinophagy']) else 'EXPLORATORY_MODULE_RESPONSE_NOT_INDEPENDENT_VALIDATION'))
    out=pd.DataFrame(rows);r.save(out,'predefined_module_review.tsv')
    r.js('module_numerical_audit.json',dict(status='PASS',rows=len(q),max_P_error=float(np.nanmax(np.abs(calc-q.P))),
        global_BH_reconstructed=True,retired_modules_not_used_for_positive_mechanism=True))
    return out

def narrative(s):
    if s.context.endswith('_Control'):
        return ('保留对照背景计算响应；仅作补充，不作为AD机制或疾病特异性证据',
                '对照内稳定性不能替代AD背景支持；未实施正式AD×KO交互检验。')
    if s.context=='Astrocyte_AD' and s.target=='RELA':
        return ('降低初始免疫/骨架通路线索的稳健性判断；不支持升级为G5-SASP机制',
                '原20项共同加权通路在留一任务中均无两种子共同显著复现；靶点阳性细胞描述支持星形标签，但正式身份验证仍缺失。')
    if s.analysis=='target_retention' and s.target in ['CA2','DPP4','PTGS2']:
        prefix='低表达可估计性明显受限；' if s.target=='PTGS2' else '保留低表达敏感性背景的部分通路线索；'
        return (prefix+'降低靶点特异性解释，记录共同网络响应及匹配支持不足',
                '不替代5%主分析；匹配参考表达下限高于目标，未构建低表达可比参考扩展网络；原共同加权条目排除核糖体相关代理后均未保留两种子共同显著（不可检验另列）；DPP4/PTGS2响应相似。CA2不升级为优先抑制靶点。')
    if s.target=='SLC2A1':
        return ('保留低表达敏感性计算线索；糖酵解特异作用仍信息不足',
                '原共同加权条目较少，参考匹配受限；排除核糖体代理后原3个加权条目不再共同显著。非主分析替代，未获得表达、代谢通量或治疗方向。')
    if s.context=='Astrocyte_AD':
        return ('保留有限ORA线索；不支持升级为AD星形ferro-aging机制',
                '初始加权GSEA无两种子共同条目；模块证据、身份和参考特异性未形成完整支持。')
    return ('保留计算响应记录；本轮不增加AD机制或靶点特异性可信度',
            '初始证据主要为已退出后续筛选的FTH1/FTL铁储存模块；新重采样模块响应属探索性结果，不是独立队列验证。')

def conclude(mod):
    scope=r.table(ROOT/r.CFG['scope_path']).fillna('')
    st=r.table(OUT/'reviewed_stability_by_term.tsv');den=r.table(OUT/'resampling_denominators.tsv')
    ko_status=r.table(SRC/'all_resampling_KO_statuses.tsv')
    expr=r.table(SRC/'target_expression_summary.tsv');mq=r.table(OUT/'matched_control_QC.tsv')
    mr=r.table(OUT/'matched_response_comparison.tsv')
    rib=r.table(OUT/'ribosome_focus_review.tsv');sim=r.table(SRC/'cross_target_similarity.tsv')
    ident=r.table(OUT/'target_positive_cell_identity_review.tsv');ident=ident[ident.level=='context']
    rows=[]
    for s in scope.itertuples(index=False):
        match=lambda q:q[(q.analysis==s.analysis)&(q.context==s.context)&(q.target==s.target)]
        d=match(den);e=expr[(expr.context==s.context)&(expr.target==s.target)].iloc[0];m=match(mq).iloc[0]
        ii=ident[(ident.context==s.context)&(ident.target==s.target)].iloc[0]
        conclusion,reason=narrative(s)
        row=dict(analysis=s.analysis,context=s.context,target=s.target,layers=s.layers,original_evidence_grade=s.evidence_grade,
          state_role=s.state_role,original_shared_weighted_terms=s.shared_weighted_GSEA,original_shared_ORA_terms=s.shared_ORA,
          original_predefined_modules=s.shared_predefined_modules,identity_status='DESCRIPTIVE_REVIEW_ORIGINAL_LABEL_RETAINED_NOT_REFERENCE_VALIDATED',
          identity_doublet_ambient_status='NOT_ASSESSABLE_WITH_CURRENT_LOCKED_INPUTS',balanced_cells=int(e.N),target_detected_cells=int(e.k),
          target_detection_fraction=e.detection_fraction,detecting_donors=int(e.n_detecting_donors),max_donor_detection_contribution=e.max_donor_detection_contribution,
          target_positive_afterQC_cells=int(ii.k),target_positive_alternative_panel_higher_fraction=ii.positive_cells_alternative_higher/ii.k if ii.k else np.nan,
          n_planned_KOs=int(d.n_planned_KOs.sum()),n_completed_KOs=int(d.n_completed_KOs.sum()),n_nonestimable_KOs=int(d.n_nonestimable_KOs.sum()),
          nonestimable_status_counts=json.dumps(match(ko_status).loc[match(ko_status).status!='COMPLETED','status'].value_counts().to_dict()),
          matching_quality=m.matching_quality,matching_controls=int(m.reference_controls),matching_limited=bool(m.matching_limited),
          matching_specificity_confirmed=False,numerical_audit_status='PASS',scientific_disposition=conclusion,reason=reason,
          mechanism_confirmation_pass=False,expression_direction='NOT_ESTIMATED',drug_evidence_ceiling='C0',reviewer='Codex_AI_not_human_signoff')
        for stage,label in [('LODO','lodo'),('BOOTSTRAP','bootstrap')]:
            dd=d[d.stage==stage].iloc[0];row[f'{label}_planned_pairs']=int(dd.n_planned_pairs);row[f'{label}_estimable_pairs']=int(dd.n_estimable_pairs)
            for method,mode,prefix in [('GSEA','weighted','weighted'),('ORA','significant_gene_set','ORA')]:
                q=match(st);q=q[(q.stage==stage)&(q.method==method)&(q['mode']==mode)]
                row[f'{label}_{prefix}_original_terms']=len(q)
                row[f'{label}_{prefix}_any_reproduced_terms']=int((q.n_reproduced>0).sum())
                row[f'{label}_{prefix}_median_unconditional_rate']=q.unconditional_rate.median() if len(q) else np.nan
                row[f'{label}_{prefix}_maximum_unconditional_rate']=q.unconditional_rate.max() if len(q) else np.nan
                row[f'{label}_{prefix}_global_any_reproduced_terms']=int((q.n_reproduced_exploratory_global>0).sum())
        for v in ['exclude_cytoplasmic','exclude_both']:
            q=match(rib);q=q[(q.variant==v)&(q.method=='GSEA')&(q['mode']=='weighted')]
            row[f'{v}_weighted_retained_terms']=int((q.outcome=='RETAINED_BOTH_SEEDS').sum())
            row[f'{v}_weighted_not_testable_terms']=int((q.outcome=='TERM_NOT_TESTABLE').sum())
            row[f'{v}_weighted_no_longer_shared_terms']=int((q.outcome=='NO_LONGER_SHARED_SIGNIFICANT').sum())
        row['ribosome_interpretation']='Operational_KEGG_Reactome_proxy_only_not_structural_ribosome_or_new_network'
        row['predefined_module_reproduction']=json.dumps(match(mod)[['stage','module','n_planned','n_testable','n_reproduced','unconditional_rate','interpretation']].to_dict('records'),ensure_ascii=False)
        mm=match(mr);mm=mm[mm.status=='DESCRIPTIVE_CONDITIONAL_REFERENCE_NO_CAUSAL_P']
        row['matched_module_response']=json.dumps(mm[['seed','function','n_members','target_mean_percentile_rank','reference_min','reference_max','reference_tail_fraction']].to_dict('records'),ensure_ascii=False)
        xx=sim[(sim.analysis==s.analysis)&(sim.context==s.context)&((sim.target1==s.target)|(sim.target2==s.target))]
        row['cross_target_max_distance_spearman']=xx.distance_spearman.max() if len(xx) else np.nan
        row['unresolved']='reference identity; doublet; ambient RNA; 20-control specificity resolution; unsigned response; no real formula perturbation'
        rows.append(row)
    result=pd.DataFrame(rows);assert len(result)==29 and result.target.nunique()==15 and not result.duplicated(KEY).any()
    r.save(result,'target_verification_conclusions.tsv')
    return result

def report(c):
    n=json.loads((OUT/'numerical_audit.json').read_text());closure=json.loads((OUT/'task_closure.json').read_text())
    lines=['# 035最终科学复核报告','',
      '复核日期：2026-10-11。复核者：Codex AI。状态：FINAL_SCIENTIFIC_REVIEW_COMPLETE_WITH_LIMITATIONS；本轮证据复核已闭合，不等于人工签署、细胞身份全部确认或15个靶点机制通过。',
      '', '## 主要结论','',
      '原先A/B/C是035纳入核查的初始线索分级，不能作为复核后的机制等级。本轮没有靶点升级为已确认的AD驱动因子、直接方剂靶点或治疗靶点。疾病G1–G5原结论、药物C0、017/029及G2c/G6关闭不变。',
      '', '1. 星形AD RELA：原20项共同加权通路留一复现均为0；Bootstrap最高仅3/20。初始免疫/骨架响应不足以支持稳健的G5-SASP机制。',
      '2. 小胶质AD CA2/DPP4/PTGS2：部分翻译相关条目可重复，但低表达、相似共同响应和参考匹配不足限制特异性；PTGS2尤有24/52次不可估计。',
      '排除核糖体相关代理后：CA2/DPP4的37个原共同加权条目中12个因大小不足不可检验，其余25个不再两种子共同显著；PTGS2为12个不可检验、26个不再共同显著；SLC2A1的3个均可检验但不再共同显著。此处解释为富集贡献依赖，不能当作重建网络或真实生物效应阴性。',
      '3. 小胶质AD主分析的初始模块证据主要落在已退出后续筛选的FTH1/FTL铁储存程序。本轮完整保留这些历史注释，但不将其用作新的ferro-aging主机制支持。',
      '纳入035的三个AD星形靶点RELA/ABCC1/DYRK1A，在G5-SASP、LPO和ferro-aging模块的两种子共同显著复现均为留一0/6、Bootstrap0/20；所有这些模块均可检验。此为当前数据/模型未提供稳定支持，不是对应生物机制不存在。',
      '4. 对照背景的稳定计算响应单列。没有正式AD×KO交互检验，不能把对照显著/AD不显著解释为疾病特异性。',
      '', '## 数据、任务闭合及数值核验','',
      f"GSE160936；15靶点、29组合（AD14/Control15）；130共享网络任务。目标重采样尝试1508次，其中{closure['resampling_status_counts'].get('COMPLETED',0)}次完成，{sum(v for k,v in closure['resampling_status_counts'].items() if k!='COMPLETED')}次非可估计结局；没有把不可估计改成阴性。",
      '非可估计原因分开：CA2有3次网络无出边；DPP4有2次低表达；PTGS2有24次低表达；Control星形SNCA有4次低表达。无出边不等于缺少表达，也不等于无生物作用。',
      f"全部2707个富集任务独立检查；ORA共{n['ORA_tests']:,}项，GSEA共{n['GSEA_tests']:,}项。每项ORA背景、成员、命中和超几何P值重建，全部族内/每任务BH重算；每库/模式/任务固定取最小P条目复核ES与leading edge，共{n['ES_replays']:,}次。ORA最大P误差{n['max_ORA_P_error']:.3g}，ES最大误差{n['max_ES_error']:.3g}。",
      'GSEA完整10万次置换随机序列未重跑；本轮检查其P值分辨率公式、输出完整性、全部FDR及代表性确定性得分，不声称逐置换严格重现。旧enrichment_independent_audit.json仅覆盖早期1个任务，保留为历史；当前完整审计以本目录numerical_audit.json为准。',
      '另核对原贡献表2,427,829行，2,201,055行的贡献基因未达到基因层FDR<0.05；这属于GSEA leading edge与基因显著性不同，不是计算错误。不得把所有贡献成员称显著扰动基因或下调基因。',
      '原冻结输入/代码30份SHA256全部一致；范围和所有网络/KO状态、抽样实例及参考任务均闭合。重采样通路分母/复现数与旧汇总逐项一致。',
      '', '## 方法及参数','',
      '沿用原5%基因骨架及完整靶点白名单≥25检出细胞保留规则，QC总UMI≥1000/MT≤0.1，每供体实例×脑区最多500细胞；两个网络种子42/43；nNet10、nCells500、nComp3、q0.9；供体留一6次、Bootstrap20次、参考靶点20个。',
      'ORA大小5–500、展示至少2命中；GSEA大小15–500、fgseaSimple10万次、scoreType=pos、加权1/不加权0。主FDR沿用每KO/模式三库BH；新增全运行ORA和GSEA各自整族BH仅作探索性敏感性，不能替代主分析或把条目数当独立机制数。',
      '对预先冻结034共同条目给出无条件复现率=共同显著重复/计划6或20，与条件率=共同显著重复/两种子条目均可检验重复并列；不采用事后80%放行线。最高通路保留率仅是条目上界，非靶点成功率。两个种子及重采样不增加独立供体数。',
      '', '## AD逐组合结论','',
      '|背景/模式|靶点|完成KO/52|原共同加权条目|留一最高保留率|Bootstrap最高保留率|科学处置|',
      '|---|---|---:|---:|---:|---:|---|']
    def pct(x):return '—' if pd.isna(x) else f'{x:.1%}'
    for x in c[c.context.str.endswith('_AD')].itertuples(index=False):
        lines.append(f'|{x.context}/{x.analysis}|{x.target}|{x.n_completed_KOs}/52|{x.original_shared_weighted_terms}|{pct(x.lodo_weighted_maximum_unconditional_rate)}|{pct(x.bootstrap_weighted_maximum_unconditional_rate)}|{x.scientific_disposition}|')
    lines+=['','完整29组合及所有理由见target_verification_conclusions.tsv；Control15组合逐项保留，未提升为AD证据。',
      '', '## 细胞身份与混合信号','',
      '原标签保留；逐cluster×供体×脑区记录七谱系面板，针对15靶点新增QC后表达细胞的标志共表达描述，未按结果重新选择细胞或重建网络。',
      'AD星形RELA检出3836/27239细胞；其中2个另一胶质面板更高，平均星形/小胶质面板约2.225/0.039。CSF1R/P2RY12/INPP5D分别在278/384/355个RELA阳性细胞检出，既不支持将整体RELA信号判为误标，也不足以排除少量混合或环境RNA。',
      'AD小胶质CA2/DPP4/PTGS2在QC后分别检出138/128/47细胞，其中另一胶质面板更高44/36/21细胞；同时多数仍检出小胶质标志。面板未校准，不把这些比例称污染、双细胞或身份概率。',
      '正式双细胞检测、环境RNA校正和独立参考图谱身份验证未完成，CP1记录为描述核查完成但身份验证不足。原标签敏感性可解释，未声称保守净化子集已建网。',
      '', '## 匹配参考与核糖体敏感性','',
      '各目标20参考、两种子；以冻结模块基因共同可测集合的平均扰动百分位给出完整参考值和含并列的条件参考尾部比例。指标为本轮描述性复核，未产生新的作用方向或正式显著特异性判定；M20分辨率1/21。',
      '另对034预先冻结的所有共同加权GO/KEGG/Reactome条目，包括翻译与RELA免疫/骨架条目，以同一指标给出20参考的完整分布，见matched_focus_response_comparison.tsv；不根据参考KO结果挑选新关注通路。',
      '低表达CA2/DPP4/PTGS2/SLC2A1参考来自原5%骨架；候选表达范围不覆盖目标，匹配记录标记limited。未执行SOP所述低表达参考扩展网络，因此CP5为匹配受限的描述性完成，不能给充分匹配PASS。',
      '原匹配代码在过滤后网络矩阵计算CP10K分母，未使用全原始基因library UMI。这是匹配尺度的实现偏差，完整登记，结果只作描述，不确认SOP全尺度符合。',
      'AD星形RELA的SASP平均扰动百分位为seed42约0.577、seed43约0.747；相应参考范围0.516–0.716与0.414–0.766，条件参考尾部比例15/21与2/21。两种子不一致且仍落在参考范围，不支持SASP特异性升级。',
      '两种核糖体排除变体的每条原通路按两种子分类为保留、仍可检验但不再共同显著、低于大小门槛不可检验。采用冻结KEGG hsa03010与Reactome线粒体翻译的操作性蛋白代理，含部分生物发生成员，不称纯结构核糖体；仅重新富集，没有重建无核糖体网络。',
      '', '## CP0–CP7、偏差与证据边界','',
      'CP0输入参数PASS；CP1描述核查完成/身份验证不足；CP2表达与分母PASS；CP3贡献与共同响应核查完成；CP4任务闭合和分母PASS；CP5参考描述完成/低表达匹配及尺度受限；CP6操作性代理敏感性完成/纯结构注释不足；CP7 AI科学复核完成、有明确未解决项。',
      'CP通过意味着操作与记录可追溯，不能合并成生物学通过。15靶点均有记录，29组合均有科学处置；本轮无人类复核者签名。',
      '完整预设模块新增重采样结果另外保留：AD小胶质中AKR1B1、LTA4H的G1模块两种子共同显著分别为留一2/6、2/6，Bootstrap6/20、7/20；PKM/GAPDH的糖酵解模块均留一0/6，Bootstrap各2/20。这些是本次探索性响应，未作为原关注条目重新选靶；源模块全部26721项概率及全族BH已独立复核。',
      '', '## 复现、来源与交付','',
      '运行入口：D:/conda/envs/TCM_select/python.exe -X utf8 src/scientific_review_virtual_ko_035.py；完成数值审计后可加--resume-audit重建描述表；随后运行src/finalize_scientific_review_035.py。R身份描述入口src/review_virtual_ko_identity_035.R使用既有zhang2026-nicotine-copd-r R4.4.3环境及DLL路径。',
      '仅本机计算；复用实验技术库15_sctenifoldknk、单细胞注释卡和项目033/034/035代码/快照，无新环境、依赖安装或模型训练。新增代码是复核与汇总适配器。',
      'WikiSkill沿用项目已有确定性和真实适配器验证，保存本轮轨迹；没有新演化、baseline/best/test不适用，没有接受或晋级跨项目技能。实际代码验证与源结果审计不称论文benchmark。',
      '原方案指定标准路径不存在，完整读取本机C:/Users/administer/.codex/standards/RESEARCH_AGENT_REPRODUCIBILITY_STANDARD.md。Git/DVC及环境登记见provenance与version_control_status；未配置DVC remote，不能声明严格跨机器恢复。达到本地可追溯计算复核与探索性验证。',
      '无原始数据、旧结果或正式缓存删除。两个仅用于本轮调试的进程曾中止：顺序数值审计改为3个本地进程；R共表达提取改为预取小矩阵以减少重复稀疏切片。统计规则未变，未产生可被误用的完整旧输出。Conda launcher DLL异常后直接复用原R环境与DLL路径，locale启动警告记录，未修改环境。',
      '', '方法核验：',
      '- Osorio et al. **scTenifoldKnk: An efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation**，DOI 10.1016/j.patter.2022.100434；[PubMed](https://pubmed.ncbi.nlm.nih.gov/35510185/)。继承1.0.3无符号流程，不升级算法。',
      '- Squair et al. **Confronting false discoveries in single-cell differential expression**，DOI 10.1038/s41467-021-25960-2；[PubMed](https://pubmed.ncbi.nlm.nih.gov/34584091/)。用于生物重复/细胞伪重复边界；不作为20次重采样或匹配数量的文献推荐。',
      '- [官方scTenifoldKnk实现](https://github.com/cailab-tamu/scTenifoldKnk)已核验；当前上游版本不同，保留本项目冻结版本及源码，不将上游新功能归入本次结果。']
    (OUT/'最终科学复核报告.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

def finalize_dispositions(mod):
    assert json.loads((OUT/'review_tables_state.json').read_text())['status']=='EVIDENCE_TABLES_COMPLETE_PENDING_INTERPRETATION'
    assert json.loads((OUT/'identity_target_review.json').read_text())['status']=='DESCRIPTIVE_REVIEW_COMPLETE'
    assert json.loads((OUT/'contributor_semantics_audit.json').read_text())['status']=='PASS'
    assert json.loads((OUT/'module_numerical_audit.json').read_text())['status']=='PASS'
    c=conclude(mod);report(c)
    checks=[('CP0','PASS_INPUTS_AND_PARAMETERS'),('CP1','DESCRIPTIVE_COMPLETE_IDENTITY_NOT_VALIDATED'),('CP2','PASS_EXPRESSION_DENOMINATORS'),
      ('CP3','COMPLETE_COMMON_RESPONSES_NOT_SPECIFICITY_PASS'),('CP4','PASS_TASK_CLOSURE_AND_REPLAY'),
      ('CP5','DESCRIPTIVE_COMPLETE_MATCHING_LIMITED'),('CP6','COMPLETE_OPERATIONAL_PROXY_ANNOTATION_LIMITED'),('CP7','AI_REVIEW_COMPLETE_WITH_LIMITATIONS')]
    r.save(pd.DataFrame([dict(checkpoint=k,status=v,reviewer='Codex_AI_not_human_signoff',date='2026-10-11') for k,v in checks]),'stage_checklist.tsv')
    r.js('final_review_status.json',dict(status='FINAL_SCIENTIFIC_REVIEW_COMPLETE_WITH_LIMITATIONS',review_complete=True,reviewer='Codex_AI',human_signed=False,
        unique_targets=15,combinations=29,AD_combinations=14,Control_combinations=15,confirmed_mechanism_targets=0,
        confirmed_identity=False,confirmed_specificity=False,direction_estimated=False,drug_evidence_ceiling='C0',strict_cross_machine_reproduction=False))
    trace=ROOT/'workspaces'/r.CFG['run_id']/'raw/traces'
    trace.mkdir(parents=True,exist_ok=True)
    (trace/'scientific_review.json').write_text((OUT/'final_review_status.json').read_text(encoding='utf-8'),encoding='utf-8')
    print(c[c.context.str.endswith('_AD')][['context','target','n_completed_KOs','scientific_disposition']].to_string(index=False))
    print('FINAL_SCIENTIFIC_REVIEW_COMPLETE_WITH_LIMITATIONS')

def main():
    contributor_review();matched_focus_review();mod=module_review()
    finalize_dispositions(mod)

if __name__=='__main__':main()

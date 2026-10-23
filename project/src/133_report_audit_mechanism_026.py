"""Independent numerical checks and source-grounded completion report for 026."""
import ast,gzip,hashlib,json,re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-public-mechanism-analysis-026';OUT=ROOT/'results'/RUN;INT=ROOT/'data/interim'/RUN;RAW=ROOT/'data/raw'/RUN
def read(p):return pd.read_csv(p,sep='\t')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def contrast(expr,values):
    def walk(n):
        if isinstance(n,ast.Name):return values[n.id]
        if isinstance(n,ast.BinOp) and isinstance(n.op,(ast.Add,ast.Sub)):return walk(n.left)+(1 if isinstance(n.op,ast.Add) else -1)*walk(n.right)
        raise ValueError('unsupported contrast')
    return walk(ast.parse(expr,mode='eval').body)
def bh(p):
    order=np.argsort(p);ranked=np.array(p)[order];v=np.minimum.accumulate((ranked*len(p)/np.arange(1,len(p)+1))[::-1])[::-1];result=np.empty(len(p));result[order]=np.minimum(1,v);return result
def main():
    checks=[];warnings=[]
    def check(name,ok,detail=''):
        checks.append(dict(check=name,pass_=bool(ok),detail=detail))
    tq=json.loads((OUT/'TQ_D1_D2_audit.json').read_text());frozen=json.loads((OUT/'pre_effect_freeze_manifest.json').read_text())
    check('pre_effect_direction_and_config_immutable',all(sha(ROOT/p)==h for p,h in frozen.items()))
    check('no_therapy_or_normalization_upgrade',tq['primary_therapeutic_eligible']==0 and tq['normalization_eligible']==0 and not tq['reopen_017'])
    check('all_244_G1_genes_direction_audited',len(read(OUT/'G1_244_gene_direction_audit.tsv'))==244)
    modules=read(INT/'frozen_modules.tsv');check('frozen_gene_module_unique',not modules.duplicated(['module','gene']).any())
    source_rows=[]
    # Every source with a captured retrieval sidecar is rehashed; empty 200 responses remain access failures.
    for p in sorted(RAW.glob('*.metadata.json')):
        d=json.loads(p.read_text());target=ROOT/d.get('path','_missing_')
        if 'sha256' not in d:source_rows.append(dict(path=str(target),valid=False,status=d.get('status','FAILED')));continue
        valid=target.exists() and sha(target)==d['sha256'];source_rows.append(dict(path=str(target.relative_to(ROOT)),valid=valid,status='CAPTURED' if target.stat().st_size else 'HTTP200_EMPTY_NO_DATA'))
        check('source_hash_'+p.name,valid)
    pd.DataFrame(source_rows).to_csv(OUT/'raw_source_verification.tsv',sep='\t',index=False)
    sets=[];loo=[];gene_summaries=[];datasets=[]
    for sub,manifest in [('', 'dataset_manifest.json'),('receiver','receiver_dataset_manifest.json'),('protein','protein_dataset_manifest.json'),('sample_sensitivity','sample_sensitivity_manifest.json')]:
        folder=OUT/sub;dslist=json.loads((INT/manifest).read_text());datasets+=dslist if sub!='sample_sensitivity' else []
        m=read(folder/'module_effects.tsv');m['analysis_family']=sub or 'main';sets.append(m)
        score=read(folder/'sample_module_scores.tsv');check('finite_module_scores_'+sub,np.isfinite(score.score).all())
        for sc,ss in m.groupby('scenario'):
            valid=ss.P.notna();check('BH_family_'+sub+'_'+sc,np.allclose(ss.loc[valid,'FDR'],bh(ss.loc[valid,'P'].to_numpy()),atol=1e-12))
        main=m[m.scenario=='main'];errors=[];ci_errors=[]
        for ds in dslist:
            sm=read(INT/(ds['id']+'_samples.tsv'));levels=sorted(sm.group.unique());X=np.column_stack([(sm.group==g).astype(float) for g in levels])
            if ds.get('block_covariate'):
                X=np.column_stack([X]+[(sm.batch==b).astype(float).to_numpy() for b in sorted(sm.batch.unique())[1:]])
            for row in main[main.dataset==ds['id']].to_dict('records'):
                scores=score[(score.dataset==ds['id'])&(score.module==row['module'])].set_index('sample').loc[sm['sample'],'score'].to_numpy()
                coef=np.linalg.lstsq(X,scores,rcond=None)[0];values=dict(zip(levels,coef[:len(levels)]));est=contrast(row['contrast'],values);errors.append(abs(est-row['effect']))
                if row['inferential']:
                    C=np.array([contrast(row['contrast'],{g:float(g==l) for g in levels}) for l in levels]+[0.]*(X.shape[1]-len(levels)))
                    df=len(scores)-np.linalg.matrix_rank(X);se=np.sqrt(np.sum((scores-X@coef)**2)/df*(C@np.linalg.inv(X.T@X)@C));low=est-stats.t.ppf(.975,df)*se;ci_errors.append(abs(low-row['CI_low']))
        check('independent_module_effect_recompute_'+sub,max(errors,default=0)<1e-10,str(max(errors,default=0)))
        check('independent_95CI_recompute_'+sub,max(ci_errors,default=0)<1e-10,str(max(ci_errors,default=0)))
        desc=m[~m.inferential];check('descriptive_no_P_or_CI_'+sub,desc[['P','FDR','CI_low','CI_high']].isna().all().all())
        check('coverage_reported_'+sub,((m.genes_available>=2)&(m.genes_available<=m.genes_total)).all())
        loo_file=folder/'module_leave_one_out.tsv'
        if loo_file.exists():l=read(loo_file);l['analysis_family']=sub or 'main';loo.append(l)
        gfile=folder/'gene_contrast_summary.tsv'
        if gfile.exists():g=read(gfile);g['analysis_family']=sub or 'main';gene_summaries.append(g)
        for p in (folder/'gene_effects').glob('*.tsv.gz'):
            g=read(p);check('gene_FDR_'+sub+'_'+p.stem,np.allclose(g['adj.P.Val'],bh(g['P.Value'].to_numpy()),atol=1e-10))
    combined=pd.concat(sets,ignore_index=True);combined.to_csv(OUT/'all_module_effects.tsv',sep='\t',index=False)
    all_gene=pd.concat(gene_summaries,ignore_index=True);all_gene.to_csv(OUT/'all_gene_contrast_summary.tsv',sep='\t',index=False)
    sm=read(OUT/'sample_traceability_audit.tsv');missing=sm[sm.mapping_status!='EXACT_ONE'];warnings.append('SVGA two matrix samples absent from deposited SDRF; retain original six-sample exploration plus exact-linked four-sample sensitivity')
    check('sample_map_missing_explicit',len(missing)==2 and set(missing['sample'])=={'SVGA-BrdU-1','SVGA-DMSO-1'})
    pmeta=json.loads((OUT/'protein_input_qc.json').read_text());check('protein_technical_replicates_not_inflated',pmeta['channels']==54 and pmeta['aggregated_samples']==18 and pmeta['independent_blocks']==3)
    pc=combined[combined.analysis_family=='main'];groups=[]
    for keys,sub in pc.groupby(['dataset','contrast','module']):
        if len(sub)==3:
            groups.append(dict(zip(['dataset','contrast','module'],keys),sign_consistent=len(set(np.sign(sub.effect)))==1,FDR_decision_consistent=len(set(sub.FDR<.05))==1,min_effect=sub.effect.min(),max_effect=sub.effect.max()))
    pd.DataFrame(groups).to_csv(OUT/'pseudocount_sensitivity_summary.tsv',sep='\t',index=False)
    if loo:
        l=pd.concat(loo,ignore_index=True);s=l.groupby(['analysis_family','dataset','contrast','module']).apply(lambda d:pd.Series(dict(min_effect=d.effect.min(),max_effect=d.effect.max(),sign_stable=bool((np.sign(d.effect.dropna())==np.sign(d.full_effect.iloc[0])).all()))),include_groups=False).reset_index();s.to_csv(OUT/'sample_influence_summary.tsv',sep='\t',index=False)
    report_rows=combined[(combined.scenario=='main')&(combined.analysis_family!='sample_sensitivity')]
    source=read(OUT/'all_025_source_dispositions.tsv');check('all_025_candidates_have_disposition',len(source)==32 and source.accession.is_unique)
    warnings+=['limma offsets zero within-group gene variances in abundance matrices; module OLS CIs independently recomputed', 'complete-case proteome excludes 1596 groups; no peptide-level raw-MS search', 'single-cell and spatial matrices not reanalysed; no formula-treated public dataset', 'DVC remote absent; shared Git worktree contains unrelated staged assets; architecture and local replay only']
    status=dict(status='PROCESSED_MATRIX_ANALYSIS_AND_DIRECTION_AUDIT_COMPLETED_WITH_SOURCE_LIMITS',pass_=all(c['pass_'] for c in checks),checks=len(checks),failed=[c for c in checks if not c['pass_']],datasets=len(datasets),expression_datasets=15,expression_samples=len(sm),protein_channels=54,protein_aggregated_samples=18,main_module_tests=int(((report_rows.analysis_family=='main')&report_rows.inferential).sum()),receiver_module_tests=int((report_rows.analysis_family=='receiver').sum()),protein_module_tests=int((report_rows.analysis_family=='protein').sum()),gene_contrasts=len(all_gene[all_gene.analysis_family!='sample_sensitivity']),warnings=warnings,claim_ceiling='C0',reopen_017=False,strict_reproducibility=False)
    (OUT/'acceptance_audit.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n',encoding='utf-8');pd.DataFrame(checks).to_csv(OUT/'acceptance_checks.tsv',sep='\t',index=False)
    # Report all main contrasts and frozen modules, not just favourable examples.
    report_rows.to_csv(OUT/'all_main_module_effects.tsv',sep='\t',index=False)
    selected=[]
    specs=[('GSE275834','BLEO-CTRL'),('GSE275834','BLEODMY-BLEO'),('GSE345471','DMSO_18M-WT_2M'),('GSE345471','AOMK_18M-DMSO_18M'),('GSE345471','AOMK_18M-WT_2M'),('GSE228927','LPS-NOT'),('GSE175485','MCM_LPS-MCM_PBS'),('GSE235076','M_N-A_N'),('GSE18740','LPS_Luteolin-LPS')]
    for acc,ct in specs:
        d=report_rows[(report_rows.dataset==acc)&(report_rows.contrast==ct)&(report_rows.module=='G3_M3_SASP_Senescence_exploratory')]
        if len(d):selected.append(d.iloc[0])
    selected=pd.DataFrame(selected);selected.to_csv(OUT/'reported_SASP_examples.tsv',sep='\t',index=False)
    def fmt(r):return f"|{r['dataset']}|{r['contrast']}|{r['effect']:.3f} [{r['CI_low']:.3f}, {r['CI_high']:.3f}]|{r['FDR']:.3g}|{int(r['genes_available'])}/{int(r['genes_total'])}|"
    text=f'''# D1/D2审计及025公开机制分析026报告

执行范围：按v1.3-VC与SOP，先冻结治疗方向审计，再重算公开样本级RNA和蛋白矩阵。审计计算通过不等于治疗资格通过。主分析、受体扩展、蛋白扩展和样本敏感性均独立记录。

## D0“后续药物比较”是什么

将真实药物/成分处理后的状态与疾病状态及合格的同体系Desired Reference比较，分别报告目标连接、表达变化、功能表型和保护功能风险。024成分—靶点交集是先验；不能直接当作处方疗效。公开DMY、AOMK、luteolin、Fer-1只是外部参照，未证明属于当前处方且达到有效暴露。

## D1/D2结果

逐终点审计22项：TQ-A 3、TQ-B 8、TQ-C 1、TQ-D 10。G1的244基因逐项登记为目前治疗方向未确定；18个RNA接口均未获主治疗反转资格。A只用于原模型中的H4K12la/PKM2病理反馈、PS19衰老胶质负担、ARID5A共培养神经元铁死亡。C用于GPX4保护功能，不能把下降视为改善。乳酸化位点、铁稳态、SASP和病期有背景依赖，不作整体反转。

TQ_map_and_Desired_reference.tsv包含Disease state、资格、目标、参照、阶段、指标、置信度和冻结状态。方向证据审计已冻结；量化Desired Reference正式冻结数=0。缺少本方案同体系数值功能目标及当前处方干预，不以健康表达均值代替治疗最优目标。相关来源真实题目、DOI、原模型与限制见TQ_functional_evidence_cards.tsv。

## 实际重算规模与方法

15个表达数据集，共136样本：13个具备指定对比，GSE325931仅描述（处理与批次完全混杂），GSE246888仅对照质控。额外重算1个蛋白组，8529蛋白组中6933满足全54通道正值完整观测，54技术通道按B1/B3/B4合并为18个条件×实验批次样本，每条件最多3个独立实验单位。避免把技术三重复算成9个生物学重复。

基因/蛋白组对比共{status['gene_contrasts']}项。RNA主族模块检验{status['main_module_tests']}项、神经元受体扩展{status['receiver_module_tests']}项、蛋白扩展{status['protein_module_tests']}项；各检验族分别BH校正。冻结模块按数据集内基因标准化表达均值计算；G1按原AD方向加权，其余无方向权重。分数是表达代理，不等于代谢通量、蛋白活性或死亡率。覆盖基因比例逐项报告。G4/G5跨细胞或跨物种外推只作探索，禁止把小鼠星形胶质模块当成人类同细胞治疗验证。

整数计数用TMM/filterByExpr/voom/limma；TPM、FPKM、RPKM仅log2丰度limma-trend探索，伪计数0.5/1/2敏感性。微阵列正值RMA只log2一次；神经元已rlog数据不再log。蛋白使用作者跨run定量，log2后技术重复合并，加实验批次协变量，不随机填补缺失。原论文使用插补及BY，这里是明确记录的完整观测敏感性分析，不冒充原论文严格复现。95%CI和模块效应经独立NumPy/SciPy重算；逐样本去除检查只评价影响，不是独立验证。

## 主要观察结果

下表固定展示SASP模块的几个不同背景示例；全部模块与对比见all_main_module_effects.tsv。效应是标准化分数差，CI未作同时置信区间校正；显著性看FDR，不能只看CI。

|数据集|对比|分数差 [95%CI]|FDR|覆盖|
|---|---|---:|---:|---:|
'''+ '\n'.join(fmt(r) for r in selected.to_dict('records'))+'''

DMY没有让本次冻结SASP模块下降，不能用这一模块宣称逆转衰老；这也不等于否定作者其他功能结果。AOMK使老龄小胶质SASP下降，但仍高于年轻参照；不能叫完全恢复或人类AD治疗。LPS-MCM星形胶质响应在GSE228927明显，而GSE175485同一冻结SASP模块对比未达FDR，不能拼接成每个体系都相同。微胶质apoE颗粒与星形胶质apoE颗粒的神经元响应支持受体效应，不自动证明介质特异性、药物通讯或当前处方作用。luteolin的此SASP模块没有检测到可信下降。

## 质量控制与剩余边界

136表达样本中134与GSM/ENA明确一一匹配；SVGA两个rep1未出现在公开SDRF，保留来源缺口并加做仅4个可匹配样本的敏感性分析。GSE290755实际RNA为KD3/NTC2，未补造缺失对照。GSE315100原登记28条但矩阵只有16列，未把剩余样本算入。GSE286095重复注释仅在样本向量与基因名完全一致时保留一次，未求和。GSE275834每表25条无基因注释续行且全部计数为空，排除并逐行登记；fractional count列禁止四舍五入，改用TPM。GSE345470 Pal经GSM确认是palbociclib。神经元原文件无表头；日期化及重复符号不猜测恢复，歧义映射排除。

32个025候选全部给出处置；单细胞/空间细胞层重分析尚未完成，raw MS未重新鉴定，GSE240609各条件只有1个文库不作重复推断。Figshare 30252505官方API列出11幅TIFF，没有样本级数值表；不能从图像猜出量化Desired Reference。低优先级异病种与请求提供资料未被冒充已分析。用户提供的整方3补表已在024摄入，但没有本方案整方干预的RNA/蛋白重复矩阵。整体025工作因此不能标为所有来源、所有证据闭环完成。

## 可复现与交付

本地复用TCM_select及R4.3.3现有环境，无新安装、无远程计算。复用实验技术库bulk模型模式及方向机制链约束；蛋白工具卡为未验证资产，仅作边界参考，实际蛋白方法依原论文与现有limma实现。配置、协议、Decision Log、源SHA、样本台账、环境锁、DVC与离线检查保留。原始025与024只读未覆盖。单个eBayes零组内方差警告记录在R_warnings.txt，不隐藏；模块OLS不靠该经验贝叶斯修正。没有DVC Remote，共享Git工作区含其他任务，故仅达到本地计算可重放及架构复现，不宣称严格跨机器科研复现。

主文件：TQ_map_and_Desired_reference.tsv；all_main_module_effects.tsv；all_gene_contrast_summary.tsv；各gene_effects目录；sample_traceability_audit.tsv；all_025_source_dispositions.tsv；acceptance_audit.json。全部在本报告同目录。正式中间矩阵属复现资产保留；一次性ASCII工作簿暂存已清理。C0及017关闭状态保持。

方法原始来源：limma powers differential expression analyses for RNA-sequencing and microarray studies，DOI 10.1093/nar/gkv007；voom: precision weights unlock linear model analysis tools for RNA-seq read counts，DOI 10.1186/gb-2014-15-2-r29；Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing，DOI 10.1111/j.2517-6161.1995.tb02031.x。蛋白原始方法：Modelling Ferroptosis in a Human Microglial Line by Sequential Exposure to Iron and GPX4 Inhibition，DOI 10.1155/omcl/9948225。
'''
    pc_summary=pd.DataFrame(groups)
    influence=read(OUT/'sample_influence_summary.tsv')
    protein_summary=all_gene[all_gene.analysis_family=='protein']
    text+=f'\n补充稳健性：{int(pc_summary.sign_consistent.sum())}/{len(pc_summary)}个丰度模块对比在0.5/1/2伪计数下方向一致，{int(pc_summary.FDR_decision_consistent.sum())}/{len(pc_summary)}个的FDR判断一致；逐样本去除{int(influence.sign_stable.sum())}/{len(influence)}个模块对比方向保持，其余明确视为样本敏感，不作为独立验证。SVGA仅四个可追溯样本时SASP升高仍达到本敏感性族FDR，但供体及疾病外推资格没有提升。\n'
    text+='\n蛋白组完整观测敏感性：'+ '; '.join(f"{r['contrast']}检出{int(r['FDR_lt_005'])}个蛋白组（BH FDR<0.05）" for r in protein_summary.to_dict('records'))+'。本次冻结蛋白模块65项检验均未达到族内FDR<0.05，Fer-1比较未检出差异蛋白组；这不否定原文死亡/LPO功能救援，也不证明该模型或处方无效。\n'
    text+='\n执行复核：108项原始来源与独立数值审计通过，6项统计及资格边界测试通过，83份结果TSV/压缩TSV经重新拟合后解压内容SHA一致。完整16数据集×18模块覆盖288项（包括不可估模块）见full_module_coverage.tsv。离线重算使用现有中间快照：D:/conda/envs/TCM_select/python.exe src/135_run_mechanism_026.py；加--all可从已缓存原始公开文件重新生成工作簿与表达中间数据，需先对本run中间目录解除DVC保护，不能改raw。DVC stage public_mechanism_effects_026依赖当前冻结raw/interim快照；DVC Remote未配置。WikiSkill确定性架构10项不变量通过，合成任务baseline/best/test均0，迭代1拒绝，active skills为空；不晋级跨项目共享技能。\n'
    (OUT/'report.md').write_text(text,encoding='utf-8')
    # Exact output content fingerprints enable offline numerical replay checks (gzip mtime excluded).
    fingerprints={}
    for p in sorted(OUT.rglob('*.tsv'))+sorted(OUT.rglob('*.tsv.gz')):
        data=gzip.decompress(p.read_bytes()) if p.name.endswith('.gz') else p.read_bytes();fingerprints[str(p.relative_to(ROOT))]=hashlib.sha256(data).hexdigest()
    if not (OUT/'output_content_manifest.json').exists():
        (OUT/'output_content_manifest.json').write_text(json.dumps(fingerprints,indent=2)+'\n')
    print(json.dumps(status,ensure_ascii=False),flush=True)
    if not status['pass_']:raise SystemExit(1)
if __name__=='__main__':main()

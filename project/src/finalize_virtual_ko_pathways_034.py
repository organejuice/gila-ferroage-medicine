"""Integrate complementary audits and add conservative, target-specific functional interpretation."""
from pathlib import Path
import hashlib,json,subprocess
import pandas as pd

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/2026-10-05-virtual-ko-pathways-034'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    a=json.loads((OUT/'acceptance_audit.json').read_text(encoding='utf-8'))
    b=json.loads((OUT/'GSEA_score_replay_audit.json').read_text(encoding='utf-8'))
    assert a['status']==b['status']=='PASS' and a['audited_jobs']==b['jobs']==84 and b['library_mode_groups']==504
    assert b['leading_edge_exact'] and b['max_ES_error']<1e-10
    marks=[json.loads(p.read_text(encoding='utf-8')) for p in (OUT/'jobs').glob('*/completion.json')]
    assert len(marks)==84 and len({m['signature'] for m in marks})==1
    state=json.loads((OUT/'run_state.json').read_text(encoding='utf-8'));assert state['exit_codes']==[0]*4 and state['completed_jobs']==84
    frozen=json.loads((OUT/'frozen_input_hashes.json').read_text(encoding='utf-8'))
    for p,h in frozen.items():assert sha(ROOT/p)==h,p
    s=pd.read_csv(OUT/'target_pathway_summary.tsv',sep='\t');st=pd.read_csv(OUT/'two_seed_pathways.tsv',sep='\t');r=pd.read_csv(OUT/'two_seed_weighting_robust_GSEA.tsv',sep='\t')
    assert len(s)==42 and len(list((OUT/'target_reports').glob('*.md')))==42
    # Main weighted and ORA findings precede large unweighted sensitivity lists.
    for row in s.to_dict('records'):
        x=st[(st.analysis==row['analysis'])&(st.context==row['context'])&(st.target==row['target'])]
        p=OUT/'target_reports'/f"{row['analysis']}__{row['context']}__{row['target']}.md"
        text=p.read_text(encoding='utf-8');before,rest=text.split('## 两种子均显著的通路',1);_,after=rest.split('## 每个种子的主要受扰动基因',1)
        lines=['## 两种子均显著的通路','','优先展示加权主分析和ORA；不加权结果作为敏感性单独展示。均为逐KO/模式跨库BH<0.05，非全运行统一FDR。','']
        for method,mode in [('GSEA','weighted'),('ORA','significant_gene_set'),('GSEA','unweighted_sensitivity')]:
            z=x[(x.method==method)&(x['mode']==mode)].sort_values('FDR42')
            lines+=[f'### {method} / {mode}（共{len(z)}项）','','|库|通路|FDR42|FDR43|共同贡献基因|','|---|---|---:|---:|---|']
            for rr in z.head(20 if mode!='unweighted_sensitivity' else 10).to_dict('records'):
                label=rr['shared_contributors'] if pd.notna(rr['shared_contributors']) else '—'
                lines.append(f"|{rr['library']}|{rr['description']}|{rr['FDR42']:.3g}|{rr['FDR43']:.3g}|{label}|")
            if not len(z):lines.append('|—|未检出两种子共同显著通路|—|—|—|')
            lines+=['','完整通路和贡献基因见two_seed_pathways.tsv；同名、嵌套或重叠条目不能计为独立机制。','']
        p.write_text(before+'\n'.join(lines)+'\n## 每个种子的主要受扰动基因'+after,encoding='utf-8')
    ad=r[r.context.str.endswith('_AD')]
    lines=['# 虚拟敲除后通路结果：科学解读','',
        '本轮将84个有效KO输出接入既有GO_BP、KEGG、Reactome快照，完成717,538项ORA和802,012项GSEA检验。所有198项原评估保留，包括低表达不可估计和敏感性中未重复主分析的记录。42份逐靶点/细胞背景报告均包含显著受扰动基因、扰动距离、网络FC、FDR、通路及贡献基因。','',
        '## AD背景重点结果','',
        '|靶点/背景|全库分析新增线索|解释边界|','|---|---|---|',
        '|G1 AKR1B1，小胶质主分析|没有两种子共同显著的加权GSEA或ORA通路；不加权GSEA有共同条目。|不能将只在不加权模式出现的结果称为强度稳健机制；维持原候选定位，不等于无生物学作用。|',
        '|G2 PKM，小胶质主分析|同样没有两种子共同显著的加权GSEA或ORA通路。|本轮没有补出稳定的强度加权糖酵解证据，仍需功能及异构体核查。|',
        '|G3/G5 RELA，星形主分析|两种子共同加权GSEA20项、ORA35项；20项加权通路也通过不加权检查，其中11项在四种组合的置换误差上限敏感性下仍通过。涉及免疫调控、细胞因子产生、细胞骨架/连接等功能。|统计条件下较一致的功能响应，不证明表达方向、SASP整体变化或药物直接结合。继承现有细胞标注限制。|',
        '|G4/G5 NLRP3|小胶质主分析及星形低表达敏感性均无共同加权GSEA条目。|仍存在不加权排序富集，不将其单独升级为稳定ferro-aging机制。|',
        '|G4 DPP4、PTGS2，小胶质低表达敏感性|共同加权GSEA分别37/38项，同样通过不加权；分别36/37项通过四组合置换误差上限检查，主要涉及翻译/核糖体等。|仅低表达敏感性背景；CA2也出现相似条目，提示共同网络响应，非DPP4/PTGS2特异铁死亡证据。|',
        '|G1 CA2、G2 SLC2A1，小胶质低表达敏感性|共同加权/不加权GSEA分别37/3项。|CA2原表达方向张力仍未解决；SLC2A1条目少，均不升级为已确认治疗靶点。|',
        '|MMP9及其他不可估计靶点|没有有效KO输出可供本轮全库分析。|低表达导致不可估计，不能写成敲除后无通路变化。|','',
        '以上“共同”使用逐KO/模式跨三库BH<0.05；如在全部KO中挑选任意显著结果，应同时查看全运行FDR_run_all_tests。共同条目仍包含数据库冗余，不能按数量排序靶点。置换区间为逐项计算误差，不是生物学置信区间或同时置信区间。','',
        '## 跨种子与加权方式均一致的AD通路','',
        '|分析|背景|靶点|共同条目数|四组合置换误差上限仍通过|','|---|---|---|---:|---:|']
    for keys,x in ad.groupby(['analysis','context','target']):lines.append('|'+ '|'.join(map(str,[*keys,len(x),int(x.all_four_MC_upper_FDR_pass.sum())]))+'|')
    lines+=['','所有条目及对应NES、两种子FDR、leading-edge genes和置换误差检查见two_seed_weighting_robust_GSEA.tsv；不是仅展示挑选的阳性结果。','',
        '## 如何读变化','',
        '基因distance是敲除前后网络表示的距离；网络FC是平方距离与背景平方距离之比，不是RNA表达倍数。ORA说明显著扰动基因是否集中于某功能；正NES说明该功能的成员集中于扰动排名高端。即使通路名称包含“positive regulation”，也不能写成通路被激活。KEGG疾病命名条目可能由共享的核糖体、免疫或代谢基因贡献，不证明存在该疾病。Leading edge不是已证实的直接调控或药物结合靶点。','',
        '本轮GO复用范围为生物过程（BP），未扩展到MF/CC。细胞注释、固定供体/脑区采样、两个随机种子、G5人类探索背景和对接来源/构象问题继承033限制。AD与对照分别计算，未做正式KO×疾病交互检验。保留原分层及优先机制/条件保留/功能候选分类。','',
        '## 审计与复现','',
        f"84/84结果通过受检背景、贡献基因和全检验范围核对，94份冻结输入哈希未变。全部ORA由独立scipy重算，最大P值误差{a['max_hypergeom_error']:.3g}；每个KO/库/模式均选取代表性GSEA条目重算，共{b['selected_terms_replayed']}项，ES最大误差{b['max_ES_error']:.3g}，leading edge完全一致。该ES重放是覆盖全部504组合的抽样核查，不宣称重算了全部GSEA条目。",'',
        '复用现有R/fgsea环境、GMT和项目富集方法；本机四个单核R进程，没有新增安装、下载或KO重算。初始自适应试算及Conda并发启动失败记录保留；正式输出经过单次环境激活、四个退出码0及结果审计确认。WikiSkill复用已有架构工作区：真实模型token验证1/1/1、无接受或晋级技能；确定性架构token分数0不代表生物学测试失败。本轮算法的人工构造集合和确定性GSEA测试通过。','',
        '代码/参数/Decision Log由Git管理，输入/正式结果由DVC本地缓存管理；未配置remote，交付为可追溯探索性验证，尚非严格跨机器复现。未删除原始计数、KO资产或验证证据。','',
        '详细入口：[完整统计报告](report.md)；[逐靶点汇总](target_pathway_summary.tsv)；[共同通路与贡献基因](two_seed_pathways.tsv)；[两种子与两种权重均一致的GSEA](two_seed_weighting_robust_GSEA.tsv)。']
    (OUT/'科学解读.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    a.update(GSEA_score_replay=b,production_signature_consistent=True,target_reports=42,status='PASS')
    (OUT/'acceptance_audit.json').write_text(json.dumps(a,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    state['status']='COMPLETED_AUDITED_AND_REPORTED'
    (OUT/'run_state.json').write_text(json.dumps(state,indent=2)+'\n',encoding='utf-8')
    sources=['configs/virtual_ko_pathways_034.yaml','src/prepare_virtual_ko_pathways_034.py','src/enrich_virtual_ko_pathways_034.R','src/run_virtual_ko_pathways_034.py','src/report_virtual_ko_pathways_034.py','src/audit_virtual_ko_GSEA_034.R','src/finalize_virtual_ko_pathways_034.py','tests/test_virtual_ko_pathways_034.R','docs/decisions/2026-10-05-virtual-ko-pathways-034.md']
    producer={'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'source_hashes':{p:sha(ROOT/p) for p in sources},'production_signature':marks[0]['signature'],'reused_libraries':'GO_BP/KEGG/Reactome stage-A frozen snapshots; all input hashes checked','reproduction':'DVC stage virtual_ko_pathways_034; existing Conda environment; no new installation','strict_cross_machine_reproduction':False}
    (OUT/'final_provenance.json').write_text(json.dumps(producer,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    outputs=[p for p in OUT.iterdir() if p.is_file() and p.name!='result_hashes.json']+list((OUT/'target_reports').glob('*.md'))
    (OUT/'result_hashes.json').write_text(json.dumps({p.relative_to(ROOT).as_posix():sha(p) for p in outputs},ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    trace=ROOT/'workspaces/2026-10-05-virtual-ko-pathways-034/raw/traces/events.jsonl'
    with trace.open('a',encoding='utf-8') as f:f.write(json.dumps({'event':'COMPLETED_AUDITED_AND_REPORTED','jobs':84,'target_reports':42,'input_hashes_checked':94})+'\n')
    print('PASS: all 84 KO outputs,42 contextual target reports,3 libraries,2 GSEA modes and complementary audits integrated.')
if __name__=='__main__':main()

"""Integrate audited KO results without promoting original evidence categories."""
from pathlib import Path
import hashlib, importlib.util, json, shutil
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = ROOT / 'results/2026-10-04-stratified-virtual-ko-033'
SENS = ROOT / 'results/2026-10-05-virtual-ko-target-retention-sensitivity-033'
OUT = ROOT / 'results/2026-10-05-virtual-ko-scientific-report-033'
spec = importlib.util.spec_from_file_location('ko_summary', ROOT/'src/133_summarize_virtual_ko.py')
model = importlib.util.module_from_spec(spec)
spec.loader.exec_module(model)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    evidence = []
    for directory, expected in [(PRIMARY,132),(SENS,66)]:
        for name in ['artifact_audit.json','completion_audit.json']:
            p=directory/name; a=json.loads(p.read_text(encoding='utf-8'))
            assert a['status'] in {'PASS','COMPLETED_WITH_NONESTIMABLE_TARGETS'}, (p,a)
            assert a['expected_jobs']==expected
            if name=='artifact_audit.json':
                assert a['audited_jobs']==expected and not a['issues']
            else:
                assert a['recorded_jobs']==expected and a['numerical_artifact_audit_pass'] and not a['errors'] and not a['missing_jobs']
            evidence.append(p)
    p=pd.read_csv(PRIMARY/'all_context_summary.tsv',sep='\t')
    q=pd.read_csv(SENS/'all_context_summary.tsv',sep='\t')
    assert len(p)==132 and len(q)==66
    assert not p.duplicated(['context','target','seed']).any()
    assert not q.duplicated(['context','target','seed']).any()
    assert p.status.value_counts().to_dict()=={'COMPLETED':66,'LOW_EXPRESSION_NOT_ESTIMABLE':66}
    assert q.status.value_counts().to_dict()=={'PRIMARY_MODEL_ALREADY_EVALUATED':34,'COMPLETED':18,'LOW_EXPRESSION_NOT_ESTIMABLE':14}
    mg=set(sum([v for k,v in model.LAYERS.items() if not k.startswith('G5')],[]))
    astro=set(sum([v for k,v in model.LAYERS.items() if k.startswith('G5')],[]))
    def expected(states):
        return {(ct,state,seed,t) for ct,targets in [('Microglia',mg),('Astrocyte',astro)] for state in states for seed in [42,43] for t in targets}
    def actual(table):
        return {(r.cell_type,r.state,int(r.seed),r.target) for r in table.itertuples()}
    assert actual(p)==expected(['AD','Control']) and actual(q)==expected(['AD'])
    docking_path=ROOT/'results/2026-10-04-expanded-docking-20/20靶点分析汇总.csv'
    docking=pd.read_csv(docking_path).set_index('target')
    assert set(p.loc[p.cell_type=='Microglia','target'])==set(docking.index) and len(docking)==20
    assert p.loc[p.cell_type=='Astrocyte','target'].nunique()==13
    st=pd.read_csv(PRIMARY/'seed_stability.tsv',sep='\t').set_index(['context','target'])
    sq=pd.read_csv(SENS/'seed_stability.tsv',sep='\t').set_index(['context','target'])
    rows=[]
    def pair(table,context,target,field):
        x=table[(table.context==context)&(table.target==target)].set_index('seed')
        assert set(x.index)=={42,43}
        return '/'.join(str(x.loc[seed,field]) for seed in [42,43])
    for layer,targets in model.LAYERS.items():
        ct='Astrocyte' if layer.startswith('G5') else 'Microglia'
        for target in targets:
            ctx=ct+'_AD'; a=p[(p.context==ctx)&(p.target==target)]
            b=q[(q.context==ctx)&(q.target==target)]
            if set(a.status)=={'COMPLETED'}: interpretation='主分析可估计；结合跨种子稳定性解释'
            elif set(b.status)=={'COMPLETED'}: interpretation='主分析低表达不可估计；仅敏感性分析可估计'
            else: interpretation='当前质控/抽样背景不可估计；不等于无生物学作用'
            primary_j=st.loc[(ctx,target),'significant_jaccard'] if (ctx,target) in st.index else None
            sens_j=sq.loc[(ctx,target),'significant_jaccard'] if (ctx,target) in sq.index else None
            rows.append(dict(layer=layer,target=target,cell_type=ct,research_category=model.PRIORITY[layer][target],
                primary_AD_status_42_43=pair(p,ctx,target,'status'),primary_AD_sig_42_43=pair(p,ctx,target,'n_significant'),
                primary_AD_jaccard=primary_j,sensitivity_AD_status_42_43=pair(q,ctx,target,'status'),
                sensitivity_AD_sig_42_43=pair(q,ctx,target,'n_significant'),sensitivity_AD_jaccard=sens_j,
                control_status_42_43=pair(p,ct+'_Control',target,'status'),control_sig_42_43=pair(p,ct+'_Control',target,'n_significant'),
                interpretation=interpretation,original_docking_scope=docking.loc[target,'scope']))
    table=pd.DataFrame(rows); OUT.mkdir(parents=True,exist_ok=True)
    assert len(table)==35 and not table.duplicated(['layer','target']).any()
    assert all(r.research_category==model.PRIORITY[r.layer][r.target] for r in table.itertuples())
    table.to_csv(OUT/'分层靶点_结果与证据边界.tsv',sep='\t',index=False)
    lines=['# 分层虚拟敲除最终报告','',
        '本轮完成20个小胶质细胞候选、13个星形胶质细胞候选的全部预定评估及数值审计。低表达不可估计与成功敲除分开报告，不能声称每个靶点都获得有效KO预测。结果不支持把原候选直接升级为已确认药物靶点、中心靶点或疾病驱动因子。','',
        '主分析：132项靶点×细胞类型×状态×种子评估，其中66项可计算敲除、66项低表达不可估计。敏感性分析：66项记录，其中18项新增可计算敲除、14项仍低表达不可估计、34项已在主分析评估而未重复计算。共保存84次有效网络敲除计算。两种分析分别进行模块检验和BH校正，不合并显著性或独立验证次数。','',
        '## 分层结果及研究定位','',
        '|层级|本轮结果|保留的研究定位|','|---|---|---|',
        '|G1|AKR1B1、LTA4H、ALOX15B、SNCA主分析可估计；CA2仅低表达敏感性可估计。小胶质AD主分析的显著集合跨种子重合度均低。|AKR1B1保留原优先机制复核入口；LTA4H/ALOX15B条件保留；CA2不升级为优先抑制靶点；SNCA保留聚集功能研究。|',
        '|G2|PKM/GAPDH主分析可估计，SLC2A1仅敏感性可估计；未得到AD背景糖酵解模块的两种子一致富集支持。|PKM继续做功能、异构体及作用方向核查；GAPDH/SLC2A1保留功能候选，不能确定直接结合。|',
        '|G3|MMP9在当前两类AD细胞背景中均不可估计；RELA可估计，但未得到AD背景SASP模块一致富集支持。|MMP9的原实验/对接复核优先级保留，本轮没有新增有效KO证据；RELA继续作为通路功能研究对象。|',
        '|G4|小胶质AD中NLRP3、DYRK1A、SETD7、ALOX5、ABCC1、SNCA主分析可估计；DPP4/PTGS2/ALOX12仅敏感性可估计；MPO/ALOX15仍不可估计。|DPP4/NLRP3/PTGS2的优先仍是复核顺序；不升级为ferro-aging核心驱动靶点。其余维持条件或功能类别。|',
        '|G5|星形AD中RELA、ABCC1、DYRK1A、SETD7、SNCA主分析可估计；ALOX12、ALOX5、MPO、NLRP3仅敏感性可估计；MMP9、ALOX15、DPP4、PTGS2仍不可估计。|分开处理SASP、LPO和综合程序。铁摄取/铁输出无本清单E1交集；不指定新靶点。人类星形结果是补充探索，不替代原小鼠G5队列。|','',
        '## 显著性与稳定性的实际含义','',
        '小胶质AD主分析12个可估计靶点的显著集合Jaccard为0–0.0182；低表达敏感性五靶点为0.0118–0.1111。星形AD主分析RELA为0.4201，其余四靶点为0–0.0237；新增四个星形敏感性靶点为0–0.0206。因此不能按单个种子的显著基因数对机制候选排序。两个种子仅检查随机网络构建，不等于独立供体或独立队列验证。','',
        '小胶质AD中，多种KO重复富集铁储存模块，但该冻结模块只有FTH1/FTL两个基因，seed42主要由FTH1、seed43由FTH1/FTL贡献。它支持报告共同网络响应，不能识别各候选特异的铁死亡或衰老机制。G1关联模块仅在DYRK1A的一个主分析种子及ALOX12的一个敏感性种子中显著，均缺乏另一种子一致支持。星形AD的SASP/LPO/综合程序未见校正后显著富集；ABCC1/SETD7的综合程序两种子富集发生在对照背景，不能推成AD治疗机制。','',
        '## 逐靶点结果','',
        '数字顺序为seed42/seed43，显著数排除KO靶点自身。—表示该背景不可估计或该敏感性任务未重复主分析。完整状态、对照结果及原对接限制见随附TSV。','',
        '|层级|靶点|原研究类别|AD主分析显著数|主分析Jaccard|AD低表达敏感性显著数|敏感性Jaccard|','|---|---|---|---|---|---|---|']
    for r in rows:
        def number(x): return '—' if x is None or pd.isna(x) else f'{x:.4f}'
        def counts(x): return x.replace('nan','—').replace('.0','')
        lines.append(f"|{r['layer']}|{r['target']}|{r['research_category']}|{counts(r['primary_AD_sig_42_43'])}|{number(r['primary_AD_jaccard'])}|{counts(r['sensitivity_AD_sig_42_43'])}|{number(r['sensitivity_AD_jaccard'])}|")
    lines+=['','## 方法、复现与边界','',
        '复用实验技术库scTenifoldKnk接口，以GSE160936原始计数按细胞类型及AD/对照分开建网。原始数据未覆盖。主要参数：细胞计数≥1000、线粒体比例≤0.1、供体×脑区上限500；主分析基因检出率≥5%；网络种子42/43、10个子网络、每网500细胞、PCR成分3、边分位数0.9、张量秩3、最多1000次迭代、误差1e-5、流形维度2、FDR<0.05。低表达敏感性仅在同一AD细胞背景的5%基因骨架中额外保留预定靶点且要求至少25个细胞检出；不等于所有基因采用25细胞过滤。阈值是本次配置/库内接口规则，不声称原论文验证了通用阈值。','',
        '全部计算在本机完成，复用zhang2026-nicotine-copd-r环境（R4.4.3、scTenifoldKnk1.0.3、scTenifoldNet1.3），没有新增安装。pcNet加速、稀疏CP及默认隐式流形实现保留数值验证记录。另测ncv80虽然真实数据重算一致，但模拟等价检查未通过，因此没有用于生产。主进程尾部解析错误经缓存恢复正常收尾；记录见Decision Log。','',
        '两个分析的数值审计检查原始输入SHA256、全部靶点/状态/种子覆盖、网络收敛、参数身份、低表达门槛、结果可读性、从保存流形重算差异调控及汇总一致性。Git/DVC用于追踪代码、输入与结果；DVC仅配置本地缓存。当前交付为可追溯探索性验证，未达到严格跨机器复现或原论文生物学复现。WikiSkill架构测试不等于生物学验证，没有提升跨项目复用技能。','',
        '细胞标签继承既有marker评分标注，未新增独立标签验证或双细胞识别；网络使用合并脑区、固定供体×脑区抽样，没有供体独立网络重复。扰动距离和FC不代表RNA上/下调，不能证明恢复疾病状态、成分抑制/激活、PKM1/PKM2选择性、直接药物结合或小胶质→星形通讯。对接的口袋、金属配位、参考复现、姿态稳定性和来源问题仍需按原记录解决。','',
        '方法论文：Osorio等，scTenifoldKnk: An efficient virtual knockout tool for gene function predictions via single-cell gene regulatory network perturbation，DOI 10.1016/j.patter.2022.100434，[PubMed](https://pubmed.ncbi.nlm.nih.gov/35510185/)。','',
        '结果明细：主分析与敏感性分析各目录中的all_context_summary.tsv、stratified_target_summary.tsv、seed_stability.tsv、module_enrichment.tsv、artifact_audit.json及provenance.json。参数、环境、源代码与输入身份详见两份provenance及配置。未删除原始数据、正式代码或验证记录；一次性核查采用内联命令，没有遗留一次性脚本。']
    (OUT/'分层虚拟敲除_最终报告.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    for d in [PRIMARY,SENS]:
        evidence += [d/n for n in ['all_context_summary.tsv','stratified_target_summary.tsv','seed_stability.tsv','module_enrichment.tsv','provenance.json']]
    evidence += [docking_path,ROOT/'src/143_build_virtual_ko_scientific_report.py',OUT/'分层靶点_结果与证据边界.tsv',OUT/'分层虚拟敲除_最终报告.md']
    # Preserve the existing audits and environment records as formal supporting assets.
    for directory in [PRIMARY,SENS]:
        support=OUT/'reproducibility'/directory.name
        support.mkdir(parents=True,exist_ok=True)
        for source in directory.iterdir():
            if source.is_file() and (source.name in {'provenance.json','input_audit.json','input_object_audit.txt','installed_knk_function.txt','target_input_coverage.tsv','audit_sessionInfo.txt','backend_tests.txt'} or 'tests.tsv' in source.name):
                dest=support/source.name
                shutil.copy2(source,dest)
                evidence.append(dest)
    support=OUT/'reproducibility';support.mkdir(exist_ok=True)
    for name in ['configs/virtual_ko_033.yaml','configs/virtual_ko_033_target_retention_sensitivity.yaml','docs/decisions/2026-10-04-stratified-virtual-ko-033.md']:
        source=ROOT/name;dest=support/source.name;shutil.copy2(source,dest);evidence.append(dest)
    manifest={'status':'PASS','primary_jobs':132,'sensitivity_records':66,'microglia_unique_targets':20,'astrocyte_unique_targets':13,
        'all_original_layers_preserved':set(table.layer)==set(model.LAYERS),'all_original_categories_preserved':True,
        'scope':'All prespecified assessments completed; nonestimable targets explicitly retained; no drug-binding, directional or communication claim.',
        'evidence':[{'path':x.relative_to(ROOT).as_posix(),'sha256':digest(x)} for x in evidence]}
    (OUT/'completion_requirements.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in manifest.items() if k!='evidence'},ensure_ascii=False))

if __name__=='__main__': main()

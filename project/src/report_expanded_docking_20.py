"""Join verified local docking, prior results and explicit noncomputed targets."""
import csv, hashlib, json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-04-expanded-docking-20';OUT=ROOT/'results'/RUN
C=json.loads((ROOT/'configs/expanded_docking_20.json').read_text(encoding='utf-8'))
def read(p,sep=','):return list(csv.DictReader(p.open(encoding='utf-8-sig'),delimiter=sep))
source=read(ROOT/'results/2026-10-04-g1-g5-docking-reuse/E1成分靶点去重清单.csv')
triage=read(ROOT/'results/2026-10-04-g1-g5-docking-reuse/E1重点候选核查表.csv')
audit=read(OUT/'pose_audit.csv');scores=read(OUT/'docking_results.csv')
old=read(ROOT/'results/2026-10-04-stage-b-tcmsp-docking-032/audited_docking_results.tsv','\t')
summary=[]
for gene in sorted({r['靶点'] for r in source}):
    pairs=[r for r in source if r['靶点']==gene];tier=next(r for r in triage if r['靶点']==gene)
    result=dict(target=gene,G_layers=tier['G层'],E1_pairs=len(pairs),calculation_source='',status='',candidate_jobs=0,best_candidate_cid='',best_candidate_score='',reference_score='',reference_RMSD_A='',scope='')
    if gene in C['reused']:
        x=[r for r in old if r['target']==gene and r['role']=='candidate'];best=min(x,key=lambda r:float(r['score_kcal_mol']))
        result.update(calculation_source='prior_TCMSP032',status='reused_exploratory',candidate_jobs=len(x),best_candidate_cid=best['cid'],best_candidate_score=best['score_kcal_mol'],scope='AKR1B1 best-score poses overlap requested NADP; tanshinoneIIA retained for review' if gene=='AKR1B1' else '6-shogaol formula source and dual enzyme function unresolved')
    elif gene in C['deferred']:
        result.update(calculation_source='none',status='deferred_no_qualified_model',scope=C['deferred'][gene])
    else:
        x=[r for r in scores if r['gene']==gene and r['role']=='candidate' and r['status']=='completed'];refs=[r for r in scores if r['gene']==gene and r['role']=='reference' and r['status']=='completed'];qa=[r for r in audit if r['gene']==gene and r['role']=='reference']
        best=min(x,key=lambda r:float(r['best_score'])) if x else None
        result.update(calculation_source='local_Vina1.2.7',status='completed_exploratory' if len(x)==len(pairs) else 'incomplete',candidate_jobs=len(x),best_candidate_cid=best['cid'] if best else '',best_candidate_score=best['best_score'] if best else '',reference_score=refs[0]['best_score'] if refs else '',reference_RMSD_A=qa[0]['reference_fixed_frame_symmetry_RMSD_A'] if qa else '',scope=C['targets'][gene][-1])
    summary.append(result)
with (OUT/'20靶点分析汇总.csv').open('w',encoding='utf-8-sig',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
count=sum(r['status']=='completed' for r in scores);candidates=sum(r['status']=='completed' and r['role']=='candidate' for r in scores)
text=f'''# 20个E1交集靶点对接分析

本轮覆盖20个去重靶点、80个唯一成分—靶点对；复用AKR1B1、LTA4H的既有平台结果，新增16个靶点的本地探索性计算。ALOX15与RELA未生成分数，原因单列。新任务完成{count}/{len(scores)}，其中候选成分任务{candidates}个。不能将20靶点范围核查表述为20靶点均已完成有效对接。

## 逐靶点结果

各行仅列该靶点内最低分的候选，供检索结果；不同受体或平台分数不用于疗效、靶点重要性或亲和力排名。

| 靶点 | G层 | 来源/状态 | 完成候选任务数 | 最低分CID | 分数(kcal/mol) | 参考配体RMSD(Å) |
|---|---|---|---:|---|---:|---:|
'''
for r in summary:
    text+=f"| {r['target']} | {r['G_layers']} | {r['calculation_source']}/{r['status']} | {r['candidate_jobs']} | {r['best_candidate_cid']} | {r['best_candidate_score']} | {r['reference_RMSD_A']} |\n"
text+='''
## 解释与限制

- **本轮没有据此确认高置信直接作用靶点。** 共晶对照按固定坐标与分子图对称重原子RMSD核查，完整连续值见pose_audit.csv。参考姿态复现较差的受体不能仅凭负分升级候选；接近共晶也不证明复方内直接结合。
- ALOX12铁中心盒无配体重对接对照；ALOX15B参考是底物模拟去污剂，不是选择性抑制剂。ALOX5使用稳定化构建体AKBA口袋，未覆盖其他作用位点。
- PKM使用PKM2-serine位点模型，不能代表PKM1或所有PKM作用位点。GAPDH的参考是共价抑制剂，本轮非共价计算只探索局部几何，不验证共价反应。
- SNCA使用E46K聚集纤维8X7L，保留整个所选晶体模型链组；不能等同野生型单体或人类AD中SNCA直接机制。ABCC1、SLC2A1只覆盖所选运输构象。
- 金属、辅因子保留，但通用Vina不专门模拟金属配位或共价化学；PTGS2保留原晶体钴血红素类似物，MPO保留原血红素及链组。距离和4Å接触仅作几何检查，不自动标为氢键。
- 既有AKR1B1最高分姜黄素/绿原酸姿态与请求保留的NADP重叠，不能作为优先证据；丹参酮IIA仍为复核候选。LTA4H的6-姜烯酚来源与双酶功能限制保留。
- AKR1B1已有19个候选任务包含15个成分结构及4个羧酸根状态，不能计为19个唯一成分—靶点对。全部唯一对仍为80个：本地59个、两既有靶点16个、两暂缓靶点5个。
- E1是已接纳实验记录连接，尚不等于实测整方成分、细胞内直接结合或有效自由暴露；模块成员身份与模块显著不能推导药理方向。
- G5按G3–SASP及G4–铁/LPO/ferro-aging程序标注，但保留星形胶质来源；重复模块成员不重复计为独立支持。铁储存及铁蛋白自噬不纳入后续筛选。

## 受体模型逐项范围

'''
text+='''
## 可继续复核的成分—靶点组合

以下为下一步结构/原assay核查顺序，不是高置信直接结合确认或疗效排名。

|组合|本轮候选分数(kcal/mol)|参考配体固定坐标RMSD(Å)|说明|
|---|---:|---:|---|
|β-谷甾醇（CID222284）—DPP4|−9.034|0.293|参考姿态复现较接近；仍需核对原assay与实际成分来源|
|姜黄素（CID969516）—DYRK1A|−8.887|0.314|候选结果可复核；不能由对接推断激酶选择性|
|槲皮素（CID5280343）—MMP9|−10.085|1.642|通用金属打分限制保留，尚未验证Zn配位化学|
|姜黄素（CID969516）—NLRP3|−7.946|0.763|仅所选NACHT/ADP构象口袋，原作用机制需核验|
|槲皮素（CID5280343）—PTGS2|−7.925|1.849|所选钴血红素结构的局部探索，PTGS1选择性未验证|

DYRK1A最低分候选哈明（Harmine，CID5280953）−9.272本身就是3ANR共晶参考配体，属于已知配体流程对照，不作为新发现靶点组合。SETD7的野黄芩素（Scutellarein，CID5281697）−7.684虽有较接近的参考配体复现，独立种子最高分姿态RMSD约6.63Å，应继续保留姿态不确定性，不能仅因表观遗传酶身份优先放行。

常用名称核对来自PubChem PUG REST缓存priority_candidate_pubchem_titles.json；精确输入仍以原SMILES/InChIKey为准。已知配体、非选择性或多靶点成分的结果不作为复方作用机制的独立证明。

'''
for r in summary:text+=f"- **{r['target']}**：{r['scope']}\n"
text+='''
## 方法与复现

本地复用技术库OpenBabel/Vina流程及032配体准备函数。Vina1.2.7、OpenBabel3.1.0、RDKit2025.09.6、Meeko0.7.1；无新增安装或远程对接。受体加氢-p7.4、配体ETKDGv3 seed42及MMFF94 maxIters2000并要求收敛。原结构状态未穷尽质子化/互变异构体。24Å盒中心来自指定共晶实例；ALOX12来自催化铁。exhaustiveness8、9模式、energy_range3、seed42020，每任务2CPU，4并发；为低成本探索初筛，不声称全面采样收敛。

方法来源：[AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization and multithreading](https://doi.org/10.1002/jcc.21334)，以及[AutoDock Vina 1.2.0: new docking methods, expanded force field, and Python bindings](https://doi.org/10.1021/acs.jcim.1c00203)。参数采用技术库初筛设置并通过参考重对接检验其局限，不因“常用”认定参数适配成功。每结构原始文章题目、DOI、坐标及受体SHA见data/processed同run目录中各靶点receptor_manifest.json。

完整任务见job_manifest.csv及docking_results.csv；所有姿态提供PDBQT及SDF，几何审计见pose_audit.csv，原位点核查见reference_site_audit.csv。配置见configs/expanded_docking_20.json；环境见Python_environment_lock.txt。

重跑：既有tcm_docking Python依次执行src/expanded_docking_20.py --prepare、src/audit_expanded_docking_20.py、src/report_expanded_docking_20.py。已完成任务采用输出缓存；参数变更需新run目录，不可沿用旧pose。正式输入、受体、姿态、日志保留；一次性重复输出已在转存后清理。raw不编辑。

交付等级为**可追溯探索性验证**，未达到跨机严格复现、普遍采样收敛或生物学验证；当前未配置DVC远端，跨机恢复未验证。WikiSkill只保存执行与审计，不晋升共享技能。
'''
if (OUT/'sensitivity_top_candidates.csv').exists():
    sens=read(OUT/'sensitivity_top_candidates.csv')
    text+='\n## 独立种子与搜索强度敏感性\n\n每个本地靶点首轮最低分候选使用seed42021、exhaustiveness16重复一次。该检查不覆盖所有候选，不证明全面收敛，也不补救参考配体复现差的模型。\n\n|靶点|CID|首轮分数|敏感性分数|差值|两次最高分姿态RMSD(Å)|\n|---|---|---:|---:|---:|---:|\n'
    for r in sens:text+=f"|{r['gene']}|{r['cid']}|{r['initial_score']}|{r['sensitivity_score']}|{r['score_difference']}|{r['pose_RMSD_A']}|\n"
text+='\n## 受体选择修正记录\n\nABCC1首选8VT4的Y01为脂质添加剂，已改用同篇研究8VUX的cGAMP底物位点；原试算位于ABCC1_lipid_pilot，不纳入当前结果。PKM参考配体下载起初混入聚合物SER，已改为游离SER A1532精确实例并按原PDB恢复N1+/O1-电荷后重算；原参考试算独立归档。SNCA参考下载起初选中另一实例，已限定IZ8 I202并重算参考；候选盒中心始终来自I202，候选结果不受这一参考下载修正影响。参考配体质心与指定实例逐项校验，防止错误实例进入最终审计。\n'
(OUT/'20靶点对接分析报告.md').write_text(text,encoding='utf-8')
hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for base in [ROOT/'data/processed'/RUN,OUT] for p in base.rglob('*') if p.is_file() and p.name!='delivery_audit.json'}
for p in [ROOT/'configs/expanded_docking_20.json',ROOT/'src/expanded_docking_20.py',ROOT/'src/audit_expanded_docking_20.py',ROOT/'src/report_expanded_docking_20.py',ROOT/'src/sensitivity_expanded_docking_20.py',ROOT/'results/2026-10-04-g1-g5-docking-reuse/E1成分靶点去重清单.csv']:
    hashes[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
record=dict(targets=len(summary),unique_pairs=len(source),covered_unique_pairs=75,local_completed_candidates=candidates,local_completed_jobs=count,reused_targets=C['reused'],deferred_targets=C['deferred'],all_20_successfully_docked=False,hashes=hashes)
(OUT/'delivery_audit.json').write_text(json.dumps(record,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in record.items() if k!='hashes'}))

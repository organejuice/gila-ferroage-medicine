"""Rebuild deliverables offline twice, test scientific guards, audit provenance."""
import gzip,hashlib,json,shutil,subprocess,sys
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-disease-interface-perturbation-audit-027';OUT=ROOT/'results'/RUN;RAW=ROOT/'data/raw'/RUN;PRIOR=ROOT/'results/2026-10-03-public-mechanism-analysis-026'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(name):return pd.read_csv(OUT/name,sep='\t')
def save(x,name):(OUT/name).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
def call(script):subprocess.run([sys.executable,'-X','utf8',str(ROOT/'src'/script)],cwd=ROOT,check=True)
def curated_hashes():
    prefixes=['D0_','P4_','P6_','P7_','candidate_','historical_input_hash_manifest']
    return {p.relative_to(OUT).as_posix():sha(p) for p in sorted(OUT.rglob('*')) if p.is_file() and any(p.relative_to(OUT).as_posix().startswith(x) for x in prefixes)}
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    for script in ['138_build_interface_qualification_027.py','139_audit_candidates_027.py']:call(script)
    baseline=curated_hashes()
    for script in ['138_build_interface_qualification_027.py','139_audit_candidates_027.py']:call(script)
    replay=curated_hashes();assert baseline==replay,'Offline replay changed curated content'
    pd.DataFrame([dict(file=k,sha256=v,second_sha256=replay[k],identical=v==replay[k]) for k,v in baseline.items()]).to_csv(OUT/'offline_replay_verification.tsv',sep='\t',index=False)
    test=subprocess.run([sys.executable,'-X','utf8','-m','unittest','discover','-s','tests','-p','test_interface_audit_027.py'],cwd=ROOT,capture_output=True,text=True,encoding='utf8');(OUT/'test_results.txt').write_text(test.stdout+test.stderr,encoding='utf8');assert test.returncode==0,test.stderr
    checks=[]
    def check(name,ok,detail=''):checks.append(dict(check=name,pass_check=bool(ok),detail=str(detail)))
    # Independently compare against the earlier accepted run's input/output manifests.
    for r in json.loads((PRIOR/'complete_input_manifest.json').read_text(encoding='utf8')):
        p=ROOT/r['path'];check('026_input_unchanged:'+r['path'],p.exists() and sha(p)==r['sha256'])
    for path,h in json.loads((PRIOR/'output_content_manifest.json').read_text(encoding='utf8')).items():
        p=ROOT/path;content=gzip.decompress(p.read_bytes()) if p.name.endswith('.gz') else p.read_bytes();check('026_output_unchanged:'+path,hashlib.sha256(content).hexdigest()==h,'026 uses decompressed-content SHA for gzip')
    for r in json.loads((OUT/'historical_input_hash_manifest.json').read_text(encoding='utf8')):
        check('curated_source_unchanged:'+r['path'],sha(ROOT/r['path'])==r['sha256'])
    for p in sorted(RAW.glob('*.retrieval.json')):
        r=json.loads(p.read_text(encoding='utf8'));check('new_snapshot_integrity:'+p.name,r['status']=='CAPTURED' and sha(ROOT/r['path'])==r['sha256'])
    reg=read('D0_interface_registry.tsv');qual=read('P4_P5_qualification_matrix.tsv');cov=read('P7_shared_measurement_coverage.tsv');dec=read('candidate_feasibility_decisions.tsv');samples=read('candidate_sample_registry.tsv')
    check('18_interfaces',len(reg)==18);check('244_G1',int(reg.loc[reg.gate=='G1','genes'].iloc[0])==244);check('all_three_candidates',set(dec.accession)=={'GSE286094','GSE315101','GSE333352'})
    check('file_header_sizes_no_mismatch',not (read('candidate_file_access.tsv').size_check=='MISMATCH').any(),'missing Content-Length remains unverified, not assumed matched')
    for gate in ['G1','G2','G3','G4','G5']:
        r=read('D0_effect_source_registry.tsv').set_index('gate').loc[gate];p=ROOT/r['source'];expected=gzip.decompress(p.read_bytes()) if p.name.endswith('.gz') else p.read_bytes();check('D0_original_text_preserved:'+gate,(OUT/'D0_original_effects'/f'{gate}.tsv').read_bytes()==expected)
    wk=json.loads((OUT/'wikiskill_architecture_audit.json').read_text(encoding='utf8'));check('WikiSkill_core_invariants',wk['checks']['all_core_invariants']);check('no_active_unvalidated_skill',wk['active_skills']==[])
    check('original_units_preserved',len(read('D0_effect_source_registry.tsv'))==5);check('all_files_listed',dec.filelist_all_matched.all());check('no_content_inferred_from_HEAD',not read('candidate_file_access.tsv').body_validated.any());check('no_qualified_formula_rank',not qual.formula_ranking_eligible.any());check('no_numeric_desired',not reg.normalization_eligible.any());check('all_guardrail_tests',test.returncode==0);check('offline_replay',baseline==replay)
    checksdf=pd.DataFrame(checks);checksdf.to_csv(OUT/'acceptance_checks.tsv',sep='\t',index=False)
    audit=dict(run_id=RUN,status='PASS' if checksdf.pass_check.all() else 'FAIL',checks=len(checks),passed=int(checksdf.pass_check.sum()),guardrail_tests=10,offline_replayed_files=len(baseline),interfaces=len(reg),qualification_rows=len(qual),context_rows=len(read('P6_context_comparability.tsv')),shared_coverage_rows=len(cov),candidate_samples_or_AOIs=len(samples),candidate_file_urls=len(read('candidate_file_access.tsv')),formal_therapy_ranking_conditions=int(qual.formula_ranking_eligible.sum()),strict_cross_machine_reproducibility=False,DVC_remote='UNCONFIGURED',claim_ceiling='C0',reopen_G6=False)
    save(audit,'acceptance_audit.json');assert audit['status']=='PASS',checksdf[~checksdf.pass_check].to_string()
    shutil.copyfile(PRIOR/'TCM_select_environment_lock.txt',OUT/'TCM_select_environment_lock.txt')
    configfiles=[ROOT/'src'/f'{i}_{name}.py' for i,name in [(136,'fetch_feasibility_027'),(137,'probe_files_027'),(138,'build_interface_qualification_027'),(139,'audit_candidates_027'),(140,'run_audit_027')]]
    configfiles += [ROOT/'tests/test_interface_audit_027.py',ROOT/'configs/disease_interface_perturbation_audit_027.json',ROOT/'configs/wikiskill_interface_audit_027.json',ROOT/'docs/protocol'/f'{RUN}.md',ROOT/'docs/decisions'/f'{RUN}.md']
    save([dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in configfiles],'code_configuration_manifest.json')
    save([dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(RAW.rglob('*')) if p.is_file()],'raw_snapshot_hash_manifest.json')
    report=f'''# D0正式接口、P4–P7资格与三候选审计027

本轮三项资料交付完成，验收PASS；“整理完成”不代表治疗资格通过。沿用v1.3-VC、017、026边界；没有重算疾病/药物效应，没有给处方疗效排名，没有重启G2c/G6。

D0正式包包含18个冻结接口、544条接口—基因成员记录、G1稳健244基因的原有方向与meta-SE、G2–G5原效应及原有不确定性、22项TQ/Desired Reference、原小鼠G5与一对一人类同源映射、153对模块重叠记录。原效应表保留已冻结来源原始文本。G2是H4K12la占位，026中的17基因RNA接口只是代理；两者不可等同。G3/G4的4/24/96小时比较分别保留；没有原CI的结果明确缺失。原模块成员与映射后的RNA接口成员可不同，原分数没有被重新命名为本轮重算分数。未建立TF/通路直接活性、完整M→A通讯或数值功能参照。原22终点TQ中A=3、B=8、C=1、D=10不变；18个RNA接口治疗方向均未获得主分析资格。

P4/P5共{len(qual)}行：36个主结果中的实际推断对比、1个SVGA四样本敏感性对比，以及3个仅描述对比和1个仅QC数据集；16个主数据集和1个敏感性数据集分开。遗传KD/KO为PC-E；DMY、AOMK、palbociclib、luteolin和ferrostatin救援保留为限定模型背景的外部化合物参照PC-B，不能自动算作当前16味处方干预。诱导、介质刺激、年龄与溶剂比较单列。原136个表达样本追溯保留，SVGA两份缺少SDRF匹配的样本仍未解决；四样本分析是敏感性，不冒充补齐。蛋白54技术通道按原实验块保留18个统计样本。名义剂量和时长的官方记录已导出；未知的相关自由暴露、应激/活力、治疗方向或匹配数值Desired Reference均未算通过。当前处方疗效排名条件=0。

P6共{len(read('P6_context_comparability.tsv'))}行，逐数据集×G1–G5登记物种、细胞、疾病模型、病期、区域、剂量暴露、时长、平台8维限制，不生成总分。同源映射不消除物种差异；细胞系与原代/iPSC/人类疾病细胞不自动等价。P7共{len(cov)}个对比×接口覆盖条目，附共同测量基因、原基因效应及CI/SE可用性、原模块各场景结果和已有样本影响/伪计数/删单基因敏感性。多探针/蛋白特征保留特征级，未人为合并；覆盖率不充当治疗资格，只有一个同源基因的模块也未被判充分。

|候选|文件与样本结果|可行性结论|
|---|---|---|
|GSE286094|11个声明生物重复文库：WT3、Atg7KO3、5xFAD-WT2、5xFAD-KO3；每份均有barcode/features/MTX记录；34条文件链接可访问，归档541,747,200字节|可准备新的模型范围分析；仍需核实动物ID、混样及性别，下载后核对矩阵/条码与注释。不能按细胞数作为动物重复。|
|GSE315101|27个文库，4个不同队列：AAV Ctrl4/DLK1-4；人iPSC少突细胞Ctrl3/sDLK1-3；TercKO5/WT5；老龄cGAS对照3。82条文件链接可访问，归档2,732,072,960字节|必须分队列；老龄cGAS队列无处理对照，不能单独估计cGAS作用。系列写“雄性”但AAV样本含雌性，依样本记录保留冲突；需核实动物/供体及混样。少突细胞不是星形胶质细胞。|
|GSE333352|143个GeoMx分段/AOI来自明确的5只小鼠；DCC及探针PKC记录齐全，145条链接可访问，DCC归档4,055,040字节。研究级设计28月龄雌鼠；各GSM未重复写年龄性别，未把研究描述改成逐样本已核实字段|仅支持老龄脑内区域/分段描述候选；独立单位是5只动物。缺少年轻/AD/药物组和已验证数值坐标；分段标记不是单细胞注释，不能建立完整M→A因果通讯。|

261/261条文件链接HEAD访问成功；181个样本/区域记录和官方归档文件清单完成审计。尚未下载大型矩阵或reads，未检查矩阵维度、barcode一一关系、数值空间坐标、细胞注释或实际QC。上述项目明确标为“未验证”；完成的是文件与样本可行性审计，不宣称已完成单细胞分析。两个单细胞队列均需确认生物独立单位，再冻样本层统计设计。

复用本地实验技术库单细胞QC卡和方向性机制链卡，以及项目06接口框架、127的SOFT解析、026冻结结果；没有新增QC阈值或统计算法，无环境重建、安装、大文件重复下载、远程计算。使用现有TCM_select Python3.14.3/pandas3.0.2；环境锁已保存。{audit['checks']}项来源/边界核查全部通过，10项科学防误判测试通过，{audit['offline_replayed_files']}个交付文件两次离线重放SHA完全一致；026原输入和已登记结果均验证未变。WikiSkill确定性架构审计核心不变量通过，合成验证/测试均0、拒绝1次迭代、active skills=0，不宣称科研基准或跨项目技能通过。

本轮Git/DVC独立登记，混合工作区旧修改保留；Decision Log、配置、环境锁、原始快照与执行轨迹可追溯。DVC Remote未配置，达到本地架构与离线重放验证，尚不满足严格跨机器复现。没有产生一次性大型缓存；正式源快照、代码和复现结果保留。

离线重放：`D:/conda/envs/TCM_select/python.exe -X utf8 src/140_run_audit_027.py`，无需网络。首次远端取证脚本136/137只取得小型元数据、文件清单及HEAD响应；现有快照按SHA复用，不覆盖原始数据。

主要文件：D0_formal_interface_package.json、D0_interface_registry.tsv、P4_P5_qualification_matrix.tsv、P6_context_comparability.tsv、P7_shared_measurement_coverage.tsv、P7_shared_gene_effects.tsv.gz、candidate_feasibility_decisions.tsv、candidate_group_independence_audit.tsv、candidate_file_access.tsv、acceptance_audit.json。全部位于本报告同目录。

建议下一步先处理GSE286094的动物独立单位与矩阵核验；GSE315101分队列后核验，GeoMx仅按5只动物的区域设计准备。后续须在新的模型范围内冻结分析，不改写历史疾病证据或把表达逆转默认当疗效。
'''
    (OUT/'report.md').write_text(report,encoding='utf8')
    trace=ROOT/'workspaces'/RUN/'raw/traces/curation_execution.json';trace.parent.mkdir(parents=True,exist_ok=True)
    if not trace.exists():trace.write_text(json.dumps(dict(run_id=RUN,backend='local_rule_based_curation_with_current_agent_review',scientific_effect_recalculation=False,data_scope='frozen_D0_P4_P7_metadata_feasibility',claim_ceiling='C0',skills_promoted=[],acceptance=audit),ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    save({p.relative_to(OUT).as_posix():sha(p) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='output_content_manifest.json'},'output_content_manifest.json')
    print(json.dumps(audit,ensure_ascii=False),flush=True)
if __name__=='__main__':main()

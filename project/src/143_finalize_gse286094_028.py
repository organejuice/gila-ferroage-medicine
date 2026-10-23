"""Offline input replay, documented gates and final audit report; no effect fitting."""
import hashlib,json,os,subprocess,sys
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-gse286094-mechanism-028';RAW=ROOT/'data/raw'/RUN;OUT=ROOT/'results'/RUN;OLD=ROOT/'results/2026-10-03-disease-interface-perturbation-audit-027'
CORE=['matrix_QC_diagnostics.tsv','sample_identity_audit.tsv','all_unfiltered_barcode_QC.tsv.gz','matrix_input_checks.tsv','analysis_status.json','feature_registry.tsv','input_QC_distributions.png']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(x,name):(OUT/name).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
def main():
    before={name:sha(OUT/name) for name in CORE}
    subprocess.run([sys.executable,'-X','utf8','src/142_audit_gse286094_028.py'],cwd=ROOT,check=True)
    after={name:sha(OUT/name) for name in CORE};assert before==after
    pd.DataFrame([dict(file=k,expected_sha256=v,replayed_sha256=after[k],identical=v==after[k]) for k,v in before.items()]).to_csv(OUT/'offline_replay_verification.tsv',sep='\t',index=False)
    checks=[]
    for name,h in json.loads((OLD/'output_content_manifest.json').read_text(encoding='utf8')).items():checks.append(dict(check='027_immutable:'+name,pass_check=sha(OLD/name)==h))
    for r in json.loads((OUT/'download_manifest.json').read_text(encoding='utf8')):checks.append(dict(check='028_raw_immutable:'+r['path'],pass_check=sha(ROOT/r['path'])==r['sha256']))
    tests=subprocess.run([sys.executable,'-X','utf8','-m','unittest','discover','-s','tests','-p','test_gse286094_audit_028.py'],cwd=ROOT,capture_output=True,text=True,encoding='utf8');(OUT/'test_results.txt').write_text(tests.stdout+tests.stderr,encoding='utf8');assert tests.returncode==0
    env=os.environ.copy();env['PYTHONUTF8']='1'
    wk=subprocess.run(['D:/conda/envs/wikiskill-reproduction/python.exe','E:/我的文档/实验技术/wikiskill-reproduction/scripts/audit.py','--project','.','--tasks','data/raw/wikiskill_tasks.json','--run-id',RUN+'-deterministic'],cwd=ROOT,env=env,capture_output=True,text=True,encoding='utf8',check=True);wiki=json.loads(wk.stdout);save(wiki,'wikiskill_architecture_audit.json')
    checks.extend([dict(check='offline_replay',pass_check=before==after),dict(check='negative_case_tests',pass_check=tests.returncode==0),dict(check='WikiSkill_core_invariants',pass_check=wiki['checks']['all_core_invariants'])]);assert all(r['pass_check'] for r in checks)
    pd.DataFrame(checks).to_csv(OUT/'provenance_acceptance_checks.tsv',sep='\t',index=False)
    st=json.loads((OUT/'analysis_status.json').read_text(encoding='utf8'));qc=pd.read_csv(OUT/'matrix_QC_diagnostics.tsv',sep='\t');matrixchecks=pd.read_csv(OUT/'matrix_input_checks.tsv',sep='\t')
    audit=dict(run_id=RUN,delivery_status='FROZEN_PROTOCOL_AND_INPUT_AUDIT_COMPLETE',audit_status='PASS',matrix_checks=len(matrixchecks),provenance_checks=len(checks),guardrail_tests=7,offline_replayed_files=len(CORE),model_mechanism_analysis_complete=False,scientific_gate='UNRESOLVED_ANIMAL_OR_POOL_IDENTITIES_AND_CELL_QC',raw_bytes=sum(p.stat().st_size for p in RAW.glob('GSM*.gz')),new_dependencies_installed=False,strict_cross_machine_reproducibility=False,reopen_G6=False,claim_ceiling='C0');save(audit,'acceptance_audit.json')
    (OUT/'TCM_select_environment_lock.txt').write_text(subprocess.check_output([sys.executable,'-m','pip','freeze'],encoding='utf8'),encoding='utf8')
    cmd=['D:/R/R-4.3.3/bin/Rscript.exe','-e',"print(sessionInfo());print(sapply(c('scDblFinder','SingleCellExperiment','edgeR','limma','Seurat'),requireNamespace,quietly=TRUE))"]
    r=subprocess.run(cmd,capture_output=True,encoding='utf8',errors='replace');(OUT/'R_environment_probe.txt').write_text(r.stdout+r.stderr,encoding='utf8')
    # Preserve error responses as retrieval evidence; they are not a paper or valid PDF/JSON.
    access=[dict(source='Europe PMC fullTextXML',file='paper_fulltext_europepmc_response.txt',result='HTTP500_NO_FULL_TEXT'),dict(source='PMC HTML',file='PMC11823820_page.html',result='HTTP403_NO_FULL_TEXT'),dict(source='PMC PDF',file='paper_pdf_response.bin',result='HTTP403_NOT_A_PDF'),dict(source='NCBI BioC',file='PMC11823820_bioc.json',result='HTTP200_ERROR_BODY_NOT_BIOC_JSON'),dict(source='RUP Crossref-listed PDF',file='publisher_fulltext_response_0.bin',result='HTTP403_NOT_A_PDF'),dict(source='Crossref DOI',file='paper_crossref_metadata.json',result='BIBLIOGRAPHIC_METADATA_ONLY'),dict(source='Europe PMC metadata',file='paper_search_metadata.json',result='BIBLIOGRAPHIC_METADATA_ONLY')]
    pd.DataFrame(access).to_csv(OUT/'paper_access_audit.tsv',sep='\t',index=False)
    blockers=[dict(item='animal_or_nonoverlapping_pool_identity',status='UNRESOLVED',needed='each GSM -> unique animal or independent nonoverlapping pool; whether libraries share animals; sex/batch if recorded'),dict(item='full_original_methods_and_sample_table',status='ACCESS_FAILED',needed='valid accessible original methods/supplement with exact sample mapping; HTTP200 error body is not source evidence'),dict(item='doublet_QC',status='NOT_EXECUTED',needed='compatible scDblFinder runtime; install/validate only once sample gate is ready'),dict(item='ambient_RNA',status='LIMITED_INPUT',needed='unfiltered droplet matrix for background-based correction, or explicit limitation plus qualified diagnostic strategy'),dict(item='cell_annotation_and_sample_level_effects',status='NOT_EXECUTED',needed='sample gate plus qualified cell QC, prespecified conservative annotation and sample-level model')]
    pd.DataFrame(blockers).to_csv(OUT/'mechanism_analysis_gate_and_required_inputs.tsv',sep='\t',index=False)
    methods=[dict(title='Loss of ATG7 in microglia impairs UPR, triggers ferroptosis, and weakens amyloid pathology control',DOI='10.1084/jem.20230173',URL='https://pubmed.ncbi.nlm.nih.gov/39945772/',purpose='original biological context; full methods not retrieved, animal mapping not inferred'),dict(title='Confronting false discoveries in single-cell differential expression',DOI='10.1038/s41467-021-25960-2',URL='https://pubmed.ncbi.nlm.nih.gov/34584091/',purpose='supports biological replicate aggregation; no pseudobulk inference executed'),dict(title='Doublet identification in single-cell sequencing data using scDblFinder',DOI='10.12688/f1000research.73600.2',URL='https://pubmed.ncbi.nlm.nih.gov/35814628/',purpose='planned required caller; GPL-3 upstream https://github.com/plger/scDblFinder; no installation or calls executed')]
    pd.DataFrame(methods).to_csv(OUT/'method_and_original_paper_sources.tsv',sep='\t',index=False)
    manifest=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(RAW.rglob('*')) if p.is_file()]
    save(manifest,'raw_source_hash_manifest.json')
    codes=[ROOT/'src'/f'{n}_{name}.py' for n,name in [(141,'capture_gse286094_028'),(142,'audit_gse286094_028'),(143,'finalize_gse286094_028')]]+[ROOT/'src/05_prepare_glial_pseudobulk.py',ROOT/'src/127_prepare_expression_026.py',ROOT/'configs/gse286094_mechanism_028.json',ROOT/'configs/wikiskill_gse286094_028.json',ROOT/'tests/test_gse286094_audit_028.py',ROOT/'docs/protocol'/f'{RUN}.md',ROOT/'docs/decisions'/f'{RUN}.md']
    save([dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in codes],'code_and_configuration_manifest.json')
    report=f'''# GSE286094冻结方案、输入与样本审计028

已完成冻结方案、33个官方processed 10x文件下载、11个文库矩阵审计和基础QC诊断；输入技术审计PASS。模型机制分析尚未执行：动物或独立混样身份门控未通过，细胞QC与注释也未完成。审计成功不等于生物分析资格通过。

研究问题：5月龄皮层CD45+细胞中ATG7条件基因型的表达变化是否依赖5xFAD背景；预定KO-WT、AD_KO-AD_WT、AD_WT-WT及二因素交互。没有方剂/药物干预，不重启017/G2c/G6。数据源GSE286094；原论文标题和DOI见method_and_original_paper_sources.tsv。原论文标题/摘要仅支持研究背景，不能替代样本独立性证据。

## 实际输入结果

33/33文件取得且gzip可完整解压，合计{audit['raw_bytes']:,}字节。11个文库均21,835行共同基因特征；列数与条码、行数与基因表相符，计数有限、非负、整数，条码在每份文库唯一，并用GSM前缀避免跨文库碰撞；UMI按列与总矩阵守恒。共96,376个已提交filtered条码，没有把它们称作本轮QC合格细胞。各样本gzip特征文件的压缩头可不同，共同基因目录比较使用解压后的文本；源压缩文件SHA单独保留。

|组别|文库数|提交条码数|
|---|---:|---:|
|WT|3|24,290|
|Atg7KO|3|25,544|
|5xFAD-WT|2|26,789|
|5xFAD-Atg7KO|3|19,753|

单文库条码范围3,791–16,276，mt比例中位数范围{qc.median_percent_mt.min():.2f}%–{qc.median_percent_mt.max():.2f}%。这些值只用于输入诊断，不能据此比较疾病丰度或宣布全部细胞质量良好。已保存UMI、基因数、mt比例分布图及96,376条逐条码记录；没有采用200基因/5%mt等固定阈值，没有过滤、双细胞替代算法或环境RNA校正。

## 样本及机制门控

GEO标题明确声明bio_rep，分组3/3/2/3已核对；但每GSM对应动物、是否混样、是否共享动物、性别和批次未给出。declared biological library与independently verified animal/pool分开，不能写“11只独立动物”。原始方法取证：PMC和出版社PDF返回403，Europe PMC全文接口500，BioC虽HTTP200却为错误说明正文；均记录失败，未当作有效论文/PDF。Crossref确认原论文DOI与出版社链接，但书目元数据不补齐实验身份。

因此四个预定条件对比和交互没有执行，未生成效应、P值/FDR、微胶质亚群标签或组成差异结论。CD45+分选集合含多种免疫细胞，不能直接当全体小胶质；分选/解离/捕获影响组成，条码数不等于脑组织绝对丰度。条件基因型包含Cx3cr1-CreERT2，未来须在原文核实Cre与Cx3cr1剂量等控制，不能把所有差异自动归于ATG7。

QC技能原文规定“Do not skip doublet calling or silently substitute another method.” 当前R4.3.3能读稀疏矩阵，并有Seurat、SingleCellExperiment、edgeR、limma，但scDblFinder未安装；本机R目录未发现可复用包。首先未通过样本门控，故本轮没有预先安装大量下游依赖或以其他方法替代，也没有声称QC通过。仅有filtered矩阵，无空滴背景，不能声称基于空滴的ambient校正已完成。注释、亚群及pseudobulk均待身份与QC门控后执行。

用户已明确回复没有额外动物/混样说明，本轮维持身份门控未通过，不反复要求补交。未来可恢复输入为各GSM对应唯一动物或独立不重叠混样说明、同一动物是否多文库、可用的原文方法/补表以及已记录的性别/批次。取得后在独立编号补充决策中确认推断单位，再完成scDblFinder、核心QC阈值冻结、保守注释及样本层设计。未知协变量不能补造；稀有类别缺失不能补零；细胞不是独立重复。

## 验证与复现

{len(matrixchecks)}项矩阵/文件核查全部通过，{len(checks)}项来源、历史027不变与门控审计通过，7项负例/计数/身份测试通过；7个核心文件（含逐条码QC和诊断图）离线重放SHA一致。raw只读，历史027结果不变。WikiSkill确定性合成架构核心不变量通过；baseline/best/test=0，迭代1拒绝、active skills=0；没有科研基准验证或跨项目技能晋级。当前Agent进行来源/门控审阅，未运行独立模型演化。

复用本地实验技术库单细胞QC、聚类注释和供体pseudobulk卡，复用项目05.read_lines和127.SOFT解析；没有复用旧D3人类阈值。原始矩阵逐文库稀疏读取，无全矩阵稠密化。本地TCM_select与R4.3.3，无新环境、依赖安装或远程计算；环境探测和锁文件已保存。Git/DVC、params及Decision Log登记本轮独立资产，混合工作区不整体提交；DVC Remote未配置，属于本地可追溯输入/离线重放验证，未达到严格跨机器复现。

下载使用原子临时part文件，完成后改名，重试的临时part清理；正式raw、QC表、源代码、图与取证失败响应作为复现/审计材料保留，不删除。未创建不必要的下游合并对象。

核心文件：analysis_status.json、matrix_QC_diagnostics.tsv、sample_identity_audit.tsv、all_unfiltered_barcode_QC.tsv.gz、input_QC_distributions.png、paper_access_audit.tsv、mechanism_analysis_gate_and_required_inputs.tsv、acceptance_audit.json。离线重放：D:/conda/envs/TCM_select/python.exe -X utf8 src/143_finalize_gse286094_028.py；仅输入审计，不触发下载或效应计算。
'''
    (OUT/'report.md').write_text(report,encoding='utf8')
    trace=ROOT/'workspaces'/RUN/'raw/traces/execution_audit.json';trace.parent.mkdir(parents=True,exist_ok=True)
    if not trace.exists():trace.write_text(json.dumps(dict(run_id=RUN,execution='local input audit',new_effects=False,acceptance=audit,scientific_gate=st,active_skills=[]),ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    save({p.relative_to(OUT).as_posix():sha(p) for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='output_content_manifest.json'},'output_content_manifest.json');print(json.dumps(audit,ensure_ascii=False),flush=True)
if __name__=='__main__':main()

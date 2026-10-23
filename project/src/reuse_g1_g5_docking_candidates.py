"""Deterministic joins of frozen stage-A evidence; no new biological inference.

Run with the existing tcm_docking Python environment. Uses standard library only.
"""
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CFG = ROOT / 'configs/g1_g5_docking_reuse.json'
config = json.loads(CFG.read_text(encoding='utf-8'))
SOURCE = ROOT / config['source_stage_a']
DOCK = ROOT / config['source_docking']
OUT = ROOT / 'results' / config['run_id']
GATES = ('G1', 'G2', 'G3', 'G4', 'G5')
inputs = []
source_hashes_before = {}

def read(path):
    inputs.append(path)
    source_hashes_before[path] = hashlib.sha256(path.read_bytes()).hexdigest()
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))

def tokens(values):
    return ';'.join(sorted({t for value in values for t in value.split(';') if t}))

def write(name, rows):
    if not rows:
        raise ValueError(f'Empty output: {name}')
    with (OUT / name).open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def require(condition, message):
    if not condition:
        raise ValueError(message)

members = read(SOURCE / 'G1_G5_target_membership.tsv')
genes = sorted({r['gene'] for r in members if any(r[g] == 'True' for g in GATES)})
evidence = read(SOURCE / 'candidate_gene_evidence_matrix.tsv')
by_gene = {r['gene']: r for r in evidence}
require(len(by_gene) == len(evidence), 'Duplicate evidence gene')
require(set(genes) <= set(by_gene), 'Membership gene lacks evidence')
connections = read(SOURCE / 'mechanism_compound_target_matrix.tsv.gz')
prior_summary = read(SOURCE / 'G1_G5_gate_layer_summary.tsv')
docking = read(DOCK / 'audited_docking_results.tsv')
directions = {r['gene']: r for r in read(SOURCE / 'G1_candidate_direction_review.tsv')}
require(len(docking) == 23, 'Unexpected docking task count')
e1 = {r['gene'] for r in members if r['tier'] == 'E1' and any(r[g] == 'True' for g in GATES)}
require(set(config['review_rules']) == e1, 'Review rules do not cover E1 candidates exactly')
active = [r for r in connections if r['gene'] in genes and any(
    r['disease_set'].startswith(g + '_') and (g != 'G1' or r['disease_set'] == 'G1_robust244')
    for g in GATES)]
pair_groups = defaultdict(list)
for row in active:
    if row['tier'] == 'E1':
        pair_groups[(row['gene'], row['compound_key'])].append(row)
require({g for g, _ in pair_groups} == e1, 'E1 membership/connection mismatch')
layer_summary = []
for old in prior_summary:
    actual = {r['gene'] for r in members if r['tier'] == old['tier'] and r[old['gate']] == 'True'}
    require(actual == set(old['gene_symbols'].split(';')), 'Frozen layer summary mismatch')
    layer_summary.append({'G层': old['gate'], '来源层': old['tier'], '靶点数': len(actual),
                          '靶点': ';'.join(sorted(actual)), '解释': 'G1为稳健244关联；G2–G5为预定义程序连接，非基因级独立验证'})

register = []
for gene in genes:
    old = by_gene[gene]
    membership = [r for r in members if r['gene'] == gene]
    gate_list = [g for g in GATES if any(r[g] == 'True' for r in membership)]
    tier_gate = ';'.join(f"{r['tier']}:{','.join(g for g in GATES if r[g] == 'True')}" for r in membership)
    pairs = [r for (g, _), rows in pair_groups.items() if g == gene for r in rows]
    docked = [r for r in docking if r['target'] == gene and r['role'] == 'candidate']
    rule = config['review_rules'].get(gene, ['辅助/预测层保留', old['structure_status'], '未取得该G1–G5范围的E1连接；不因网络或预测升级为直接靶点'])
    register.append({'靶点': gene, 'G层': ';'.join(gate_list), '逐证据层G归属': tier_gate,
        '证据层': tokens(r['tier'] for r in membership), '原始接口': old['interfaces'],
        'E1去重成分靶点对数': len({r['compound_key'] for r in pairs}),
        'E1候选CID': tokens(r['cids'] for r in pairs), 'E1原始assay_ID': tokens(r['assay_ids'] for r in pairs),
        'E1来源文献DOI或标识': tokens(r['documents'] for r in pairs),
        'G1稳健244成员': 'G1' in gate_list,
        'G1全受检meta_logFC': old['G1_full_meta_logFC'], 'G1全受检FDR': old['G1_full_meta_FDR'],
        '治疗方向': '未确定', 'G1原方向审阅': directions.get(gene, {}).get('direction_review', ''),
        '作用方向原字段': tokens(r['action_direction'] for r in pairs),
        '结构状态': rule[1], '已对接': bool(docked), '已对接候选任务数': len(docked),
        '已对接CID': tokens(r['cid'] for r in docked), '建议': rule[0], '核查事项': rule[2],
        '本方实测成分': '未建立', '相关自由暴露': '未建立', '药物证据上限': 'C0'})
compound_rows = []
for (gene, key), rows in sorted(pair_groups.items()):
    compound_rows.append({'靶点': gene, '精确结构键': key, 'CID': tokens(r['cids'] for r in rows),
        '原化学名称': tokens(r['names'] for r in rows), '目录药味关联_非实测': tokens(r['herbs'] for r in rows),
        '来源接口': tokens(r['disease_set'] for r in rows), '原始assay_ID': tokens(r['assay_ids'] for r in rows),
        '来源文献DOI或标识': tokens(r['documents'] for r in rows), 'assay类型': tokens(r['assay_types'] for r in rows),
        '来源数解释': '接口重复去重为一个结构键；不相加source_records作为独立证据数',
        '身份状态': tokens(r['identity_status'] for r in rows), 'SMILES': tokens(r['smiles'] for r in rows),
        'InChIKey': tokens(r['inchikey'] for r in rows), '作用方向原字段': tokens(r['action_direction'] for r in rows),
        '直接结合原状态': tokens(r['direct_binding'] for r in rows),
        '已有对接CID命中': any(r['target'] == gene and r['role'] == 'candidate' and r['cid'] in tokens(x['cids'] for x in rows).split(';') for r in docking)})
OUT.mkdir(parents=True, exist_ok=True)
write('G1-G5分层候选全表.csv', register)
write('E1重点候选核查表.csv', [r for r in register if r['靶点'] in e1])
write('E1成分靶点去重清单.csv', compound_rows)
write('各G证据层覆盖.csv', layer_summary)
excluded = [{'靶点': r['gene'], '原始接口': r['interfaces'], '理由': '仅属于其他/历史探索接口，未进入冻结G1–G5成员表'} for r in evidence if r['gene'] not in genes]
write('范围外历史候选.csv', excluded)
checks = {'membership_rows': len(members), 'candidate_genes': len(genes), 'E1_genes': len(e1),
          'E1_compound_target_pairs_deduplicated': len(compound_rows), 'E1_interface_rows': sum(map(len, pair_groups.values())),
          'excluded_historical_candidates': len(excluded), 'docked_targets': sorted({r['target'] for r in docking}),
          'stage_a_sources_unchanged': True, 'no_new_docking': True, 'no_new_model_training': True}
require(all(hashlib.sha256(p.read_bytes()).hexdigest() == digest for p, digest in source_hashes_before.items()),
        'Source changed during aggregation')
inputs.extend([CFG, Path(__file__)])
checks['input_sha256'] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
(OUT / 'reuse_audit.json').write_text(json.dumps(checks, ensure_ascii=False, indent=2), encoding='utf-8')
lines = ['# G1–G5分层候选复用结果', '',
    f'复用阶段A030和对接032：{len(genes)}个候选基因，{len(e1)}个E1候选靶点，{len(compound_rows)}个去重E1成分—靶点对。', '',
    'E1为阶段A已接纳的实验记录连接，不自动代表细胞内直接结合或本复方药效；E2为辅助实验、P为预测、H为待核验关联，均沿用原定义。', '',
    'G1为人类小胶质稳健244关联。G2–G5的交集是预定义模型程序成员连接，不代表每个成员已获得对应G层的独立基因级验证。尤其G4/G5模块重叠不可重复计为独立支持。', '',
    '## E1重点候选核查表', '', '|靶点|G层|去重成分对|已有对接|建议|核查事项|', '|---|---|---:|---|---|---|']
for r in register:
    if r['靶点'] in e1:
        lines.append(f"|{r['靶点']}|{r['G层']}|{r['E1去重成分靶点对数']}|{'是' if r['已对接'] else '否'}|{r['建议']}|{r['核查事项']}|")
lines += ['', '## 后续有限对接建议', '',
    '- 已有：AKR1B1优先复核；LTA4H条件保留，不重复计为新候选。',
    '- G2：先核PKM原assay异构体和作用方向，再决定是否能以PKM2进入有限对接；GAPDH、SLC2A1为备选。',
    '- G3/G5：MMP9优先补结构和原assay核查；RELA先确认是否直接作用及口袋适用性。',
    '- G4/G5：SETD7、DYRK1A、DPP4优先补结构和原assay核查。这里的优先仅是核查工作顺序，未放行对接、未确认药效；三者不需要全部计算。',
    '- ALOX家族、MPO、PTGS2及NLRP3保留条件核查；CA2、SNCA、ABCC1按表中功能与方向问题分别处理。',
    '- 新靶点的PDB、共晶对照、辅因子、口袋及成分实际来源尚未本轮新检索/验收，结构可用性仍待定。记录原文DOI只供追溯，不声称本轮重新核验了原始论文。', '',
    '## 复用、审计与复现', '',
    '- 全部计算在本机既有tcm_docking环境进行，仅Python标准库，无安装、下载或远程任务。',
    '- 使用实验技术库多维优先级与质量门槛原则，复用既有组件；不新建加权评分、算法或疗效排行榜。技术库卡片自身未本地验证，不当作方法性能证据。',
    '- 审计核查逐G逐层成员集合与旧汇总完全一致、E1成员与结构键清单一致、两已对接靶点映射、源文件SHA256。上游原始数据不修改；不重跑G1–G5或阶段A统计。',
    '- 216候选全表保留每个证据层的G归属，避免将低层的某G连接误写成该G的E1实验支持。25个范围外历史候选另表保留，不混入主清单。',
    '- 80个成分—靶点对按gene+compound_key去重；147条接口记录不是147个独立实验，assay/document ID保留而不相加来源数。',
    '- 治疗方向、实测成分、相关自由暴露未建立；药物证据上限C0、G2c/G6/017/029既有范围不改变。',
    '- 重跑入口：D:\\conda\\envs\\tcm_docking\\python.exe src/reuse_g1_g5_docking_candidates.py。规则冻结在configs/g1_g5_docking_reuse.json；输入和脚本哈希见reuse_audit.json。',
    '- 本次完成确定性复用与追溯，不代表新增结构验证、严格上游对接复现或生物学验证。DVC remote恢复未核验。正式脚本保留，无一次性缓存需清理。']
(OUT / '复用报告.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
print(json.dumps({k:v for k,v in checks.items() if k != 'input_sha256'}, ensure_ascii=False))

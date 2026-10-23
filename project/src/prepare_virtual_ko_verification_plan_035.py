"""Freeze a reviewable verification scope from existing KO evidence; no KO is run."""
from pathlib import Path
import csv
import hashlib
import json
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs/protocol/2026-10-05-virtual-ko-verification-035'
PATHWAYS = ROOT / 'results/2026-10-05-virtual-ko-pathways-034'
MODULES = {
    'primary': ROOT / 'results/2026-10-04-stratified-virtual-ko-033/module_enrichment.tsv',
    'target_retention': ROOT / 'results/2026-10-05-virtual-ko-target-retention-sensitivity-033/module_enrichment.tsv',
}


def read(path):
    with path.open(encoding='utf-8-sig', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summary = read(PATHWAYS / 'target_pathway_summary.tsv')
    robust = read(PATHWAYS / 'two_seed_weighting_robust_GSEA.tsv')
    mc_pass = defaultdict(int)
    for row in robust:
        if row['all_four_MC_upper_FDR_pass'].lower() == 'true':
            mc_pass[(row['analysis'], row['context'], row['target'])] += 1
    module_seeds = defaultdict(set)
    for analysis, path in MODULES.items():
        for row in read(path):
            fdr = row['FDR_all_KO_module_tests']
            if row['status'] == 'TESTABLE' and fdr and float(fdr) < 0.05:
                module_seeds[(analysis, row['context'], row['target'], row['module'])].add(row['seed'])
    common_modules = defaultdict(list)
    for (analysis, context, target, module), seeds in module_seeds.items():
        if seeds == {'42', '43'}:
            common_modules[(analysis, context, target)].append(module)
    records = []
    for row in summary:
        key = (row['analysis'], row['context'], row['target'])
        weighted = int(row['GSEA_weighted_shared'])
        ora = int(row['ORA_significant_gene_set_shared'])
        modules = sorted(common_modules[key])
        included = bool(weighted or ora or modules)
        grade = ('A_较一致通路线索' if weighted and mc_pass[key] else
                 'B_较弱通路线索' if weighted or ora else
                 'C_仅预设模块线索' if modules else 'D_未达到本轮纳入规则')
        records.append(dict(
            analysis=row['analysis'], context=row['context'], target=row['target'],
            layers=row['layers'], research_categories=row['research_categories'],
            included=included, evidence_grade=grade,
            state_role='AD主核查' if row['context'].endswith('_AD') else 'Control补充核查',
            shared_weighted_GSEA=weighted, shared_ORA=ora,
            shared_predefined_modules=';'.join(modules), MC_upper_pass_terms=mc_pass[key],
            shared_significant_genes=row['shared_genes'],
            limitation=('低表达靶点保留敏感性背景；非主分析替代' if row['analysis']=='target_retention' else '')
            + ('；Control线索不能移作AD机制' if row['context'].endswith('_Control') else '')
            + ('；两基因铁储存模块可由一个命中基因贡献，属探索性线索'
               if any('Iron_Storage' in m for m in modules) else ''),
            source_report=f"results/2026-10-05-virtual-ko-pathways-034/target_reports/{row['analysis']}__{row['context']}__{row['target']}.md",
        ))
    records.sort(key=lambda r: (r['state_role']!='AD主核查', r['evidence_grade'], r['context'], r['target']))
    selected = [r for r in records if r['included']]
    assert len({(r['analysis'],r['context'],r['target']) for r in records}) == len(records)
    assert all((ROOT/r['source_report']).is_file() for r in records)
    assert {r['target'] for r in selected} == {
        'AKR1B1','LTA4H','CA2','SNCA','PKM','GAPDH','SLC2A1','RELA',
        'DPP4','PTGS2','DYRK1A','SETD7','ABCC1','ALOX5','NLRP3',
    }
    assert len(selected)==29 and sum(r['state_role']=='AD主核查' for r in selected)==14
    write(OUT/'all_valid_result_dispositions.tsv', records)
    write(OUT/'verification_targets.tsv', selected)
    grouped=[]
    for target in sorted({r['target'] for r in selected}):
        rr=[r for r in selected if r['target']==target]
        grouped.append(dict(target=target, layers=';'.join(sorted(set(';'.join(r['layers'] for r in rr).split(';')))),
                            AD_contexts=';'.join(f"{r['analysis']}:{r['context']}" for r in rr if r['state_role']=='AD主核查'),
                            Control_contexts=';'.join(f"{r['analysis']}:{r['context']}" for r in rr if r['state_role']=='Control补充核查'),
                            n_selected_contexts=len(rr)))
    write(OUT/'unique_targets.tsv', grouped)
    sources = [PATHWAYS/'target_pathway_summary.tsv', PATHWAYS/'two_seed_weighting_robust_GSEA.tsv',
               PATHWAYS/'stratified_assessment_registry.tsv', *MODULES.values()]
    audit=dict(status='PASS', analysis_status='PLAN_ONLY_NO_NEW_KO_OR_CELL_AUDIT',
               valid_contexts_reviewed=len(records), selected_contexts=len(selected),
               unique_targets=len(grouped), AD_contexts=14, Control_contexts=15,
               source_hashes={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
               checks=['all valid results retained with disposition', 'unique target-context keys',
                       'all individual source reports exist', 'expected scope independently asserted',
                       'only two-seed common evidence qualifies', 'no unweighted-only qualification'])
    (OUT/'scope_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in audit.items() if k!='source_hashes'},ensure_ascii=False))


if __name__ == '__main__':
    main()

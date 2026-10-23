"""Replay frozen 024 evidence plus new source snapshots into a separate run."""
import importlib.util
import json
import shutil
from collections import Counter,defaultdict
from pathlib import Path

spec=importlib.util.spec_from_file_location('c024',Path(__file__).with_name('118_continue_formula_supplement_024.py'))
c=importlib.util.module_from_spec(spec); spec.loader.exec_module(c)
m=c.m


def main():
    c.OUT.mkdir(parents=True,exist_ok=True)
    files=['compound_identity_audit.tsv','evidence_ChEMBL.json','evidence_PubChem.json','ledger_ChEMBL.tsv','ledger_PubChem.tsv','prior_source_access.tsv','manual_assay_adjudication.tsv','source_sha256.tsv']
    manifest=[]
    for name in files:
        dest=c.OUT/name
        if not dest.exists(): shutil.copyfile(c.OLD/name,dest)
        actual=m.filehash(dest)
        if actual!=m.filehash(c.OLD/name): raise ValueError('baseline reusable input changed: '+name)
        manifest.append(dict(path=str(dest.relative_to(m.ROOT)),source=str((c.OLD/name).relative_to(m.ROOT)),sha256=actual))
    m.writetable(c.OUT/'reused_input_manifest.tsv',manifest)
    f_spec=importlib.util.spec_from_file_location('final024',Path(__file__).with_name('114_finalize_formula_targets_024.py'))
    f=importlib.util.module_from_spec(f_spec); f_spec.loader.exec_module(f)
    f.OUT=c.OUT
    # Existing raw snapshots are read through their original paths, without copies.
    f.RAW=m.ROOT/'data/raw/2026-09-30-fixed-formula-target-supplement-024'
    original_prior=f.prior_tables
    def joined_prior(queries):
        priors,ledger=original_prior(queries)
        hp=c.OUT/'evidence_HERB.json'; lp=c.OUT/'HERB_CID_query_ledger.tsv'
        extra=json.loads(hp.read_text(encoding='utf-8')) if hp.exists() else []
        prior_by_key=defaultdict(list)
        for r in extra: prior_by_key[r['compound_key']].append(r)
        statuses={r['cid']:r for r in m.readtable(lp)} if lp.exists() else {}
        for row in ledger:
            if row['source']!='HERB': continue
            q=next(q for q in queries if q['compound_key']==row['compound_key'])
            cids=[cid for cid in q['confirmed_cids'].split(';') if cid]
            if not cids: continue
            done=[statuses.get(cid,{}).get('status','RUNNING_NOT_YET_QUERIED') for cid in cids]
            if all(s.startswith('SUCCESS') for s in done): row['status']='SUCCESS_WITH_RECORDS' if prior_by_key[row['compound_key']] else 'SUCCESS_NO_RECORDS'
            elif any(s=='RUNNING_NOT_YET_QUERIED' for s in done): row['status']='RUNNING_NOT_YET_QUERIED'
            else: row['status']='QUERY_FAILED_OR_SOURCE_BLOCKED'
            row['detail']='public HERB CID search + exact CID/structure-verified details; associations remain H'
            row['pages_complete']=all(s.startswith('SUCCESS') for s in done)
        return priors+extra,ledger
    f.prior_tables=joined_prior
    f.main()
    state=json.loads((c.OUT/'machine_audit.json').read_text())
    state['run_id']=c.RUN
    state['parent_run_id']=c.OLD.name
    legacy=m.readtable(c.OUT/'legacy_structure_identity_followup.tsv') if (c.OUT/'legacy_structure_identity_followup.tsv').exists() else []
    state['legacy_identity_followup']=dict(Counter(r['status'] for r in legacy))
    state['new_legacy_candidate_releases']=sum(r.get('candidate_release')=='True' for r in legacy)
    paper_path=c.OUT/'publication_supplement_status.json'
    paper=json.loads(paper_path.read_text(encoding='utf-8')) if paper_path.exists() else None
    state['publication_supplements']=paper or {'status':'SOURCE_BLOCKED'}
    source_rows=m.readtable(c.OUT/'query_ledger.tsv')
    stage_running=any(r['status']=='RUNNING_NOT_YET_QUERIED' for r in source_rows) or len(legacy)<642
    state['all_available_source_runs_finished']=not stage_running
    state['status']='IN_PROGRESS' if stage_running else 'PARTIAL_SOURCE_BLOCKED'
    m.writejson(c.OUT/'machine_audit.json',state)
    groups=defaultdict(list)
    for r in m.readtable(c.OUT/'G1_244_layered_connections.tsv')+m.readtable(c.OUT/'exploratory_disease_connections.tsv'):
        if r['disease_set'].startswith('G1_significant'): continue
        groups[(r['disease_set'][:2],r['tier'])].append(r)
    summary=[dict(gate=g,tier=t,genes=len({r['gene'] for r in rows}),compound_target_pairs=len({(r['compound_key'],r['gene']) for r in rows}),source_records=len(rows),gene_symbols=';'.join(sorted({r['gene'] for r in rows}))) for (g,t),rows in sorted(groups.items())]
    m.writetable(c.OUT/'G1_G5_gate_layer_summary.tsv',summary)
    lines=['# 固定处方公共靶点补查024续跑报告','',f'Run ID：{c.RUN}。状态：{state["status"]}。药物C0；017关闭范围保持。','',
           '9月30日结果与原始响应保持；本目录独立合并新查询与历史成功记录。查询无记录只表示该公开检索范围无命中；失败/未查/身份暂停均不作为阴性。','',
           '## 来源覆盖','', '|来源|状态|结构数|','|---|---|---:|']
    lines.extend(f'|{r["source"]}|{r["status"]}|{r["structures"]}|' for r in state['source_coverage'])
    lines+=['','## 疾病连接','', '|门|层|基因数|成分靶点对|来源记录|','|---|---|---:|---:|---:|']
    lines.extend(f'|{r["gate"]}|{r["tier"]}|{r["genes"]}|{r["compound_target_pairs"]}|{r["source_records"]}|' for r in summary)
    lines+=['','E1人工放行范围仍限G1主集合；G2–G5保留物种、模型、RNA代理和同源映射限制。HERB原始数据库关联列H，PMID不等于靶点特异实验。BindingDB缺字段记录不升级。疾病表达方向不转换成治疗方向。','',
            '## 结构与论文待补','',f'旧结构追查状态：{json.dumps(state["legacy_identity_followup"],ensure_ascii=False)}。冻结候选新放行结构：{state["new_legacy_candidate_releases"]}。旧目录对应外部CID只用于纠错溯源，未经协议扩大候选。',
            '三个补充文件已由用户提供并归档：154条成分记录、37个PPI核心靶点及Western blot补图；处方等价性未建立，无逐成分—靶点对，独立H比较见publication_supplement_review.md。未升级E1或扩大候选。' if paper else '同名论文补表原链接重试仍受阻；Figshare DOI检索未命中。未取得文件，未声称同处方或精确靶点清单复现。','',
            '## 执行与复现','', '本机TCM_select复用；无安装或环境重建。方法复用113/114及实验技术库候选身份/分层约束。参数见configs/fixed_formula_supplement_024_continuation.json；决定见docs/decisions/2026-10-03-supplement-024-continuation.md。',
            '离线复现：TCM_select/python.exe -X utf8 src/121_finalize_continuation_024.py。网络查询须分别运行118/119/120；失败轮次与原始响应独立保存。该目录状态为IN_PROGRESS时不能当作最终验收。']
    (c.OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({'status':state['status'],'G1':state['G1_summary'],'legacy':state['legacy_identity_followup']},ensure_ascii=False),flush=True)


if __name__=='__main__': main()

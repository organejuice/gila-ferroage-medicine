"""Offline replay and source-aware audit for supplement 024 (no network)."""
from pathlib import Path
import importlib.util
import gzip
import json
import re
import ast
from collections import Counter, defaultdict
import yaml

spec=importlib.util.spec_from_file_location('s024',Path(__file__).with_name('113_supplement_formula_targets_024.py'))
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
ROOT=m.ROOT
cfg=json.loads(m.CONFIG.read_text(encoding='utf-8'))
OUT=ROOT/'results'/cfg['run_id']; RAW=ROOT/'data/raw'/cfg['run_id']


def positive(r):
    return (r.get('identity_exact') is True and r.get('candidate_identity_verified') is not False
            and r.get('species')=='9606' and bool(r.get('target_symbol')) and bool(r.get('target_stable_id'))
            and r.get('activity_outcome','').lower() in m.VALIDATED_RULES['positive_experimental_outcomes']
            and r.get('relation') not in m.VALIDATED_RULES['relations_not_positive_confirmation'])


def prior_tables(queries):
    hgnc=m.readtable(RAW/'HGNC/hgnc_complete_set.tsv')
    symbols={r['symbol']:r for r in hgnc if r['status']=='Approved'}
    entrez={r['entrez_id']:r for r in symbols.values() if r['entrez_id']}
    bycid=defaultdict(list)
    for q in queries:
        for cid in q['confirmed_cids'].split(';'):
            if cid: bycid[cid].append(q)
    access=m.readtable(OUT/'prior_source_access.tsv'); records=[]; unresolved=[]; hits=defaultdict(Counter)
    for a in access:
        if a['source']!='BATMAN' or 'browse_by_ingredients' not in a['url']: continue
        prediction='predicted_' in a['url']
        lines=gzip.decompress((ROOT/a['body_path']).read_bytes()).decode('utf-8-sig').splitlines()
        for line in lines[1:]:
            if not line.strip(): continue
            if prediction:
                first,_,rest=line.partition(' '); name,_,targets=rest.rpartition(' ')
            else:
                parts=line.split('\t')
                if len(parts)!=3: raise ValueError('BATMAN known schema changed')
                first,name,targets=parts
            if first not in bycid: continue
            for token in targets.split('|'):
                score=''; gene=''
                if prediction:
                    match=re.fullmatch(r'(\d+)\(([^()]*)\)',token)
                    if not match:
                        unresolved.append({'CID':first,'token':token,'reason':'unparsed BATMAN prediction ID/score'}); continue
                    gene,score=match.groups(); target=entrez.get(gene)
                else: target=symbols.get(token)
                if not target:
                    unresolved.append({'CID':first,'token':token,'reason':'not an approved HGNC ID/symbol'}); continue
                for q in bycid[first]:
                    r=m.base_evidence(q,'BATMAN'); r.update(identity_exact=True,actual_compound_id='CID:'+first,
                      target_symbol=target['symbol'],target_stable_id='NCBI_Gene:'+target['entrez_id'],species='9606',
                      predicted=prediction,association_only=not prediction,prediction_method='BATMAN-TCM2.0 similarity-based prediction; DOI10.1093/nar/gkad926' if prediction else '',
                      prediction_score=score,activity_id=('predicted' if prediction else 'known')+':'+first+':'+target['symbol'],
                      assay_id='',document='10.1093/nar/gkad926',response_sha256=a['sha256'],source_url=a['url'],activity_outcome='predicted' if prediction else 'database_known_association')
                    r['tier'],r['tier_reason']=m.classify(r); records.append(r); hits[q['compound_key']][r['tier']]+=1
    ledgers=[]
    for q in queries:
        has=bool(q['confirmed_cids'])
        ledgers.append(m.ledger(q,'BATMAN','SUCCESS_WITH_RECORDS' if hits[q['compound_key']] else 'SUCCESS_NO_RECORDS' if has else 'IDENTITY_HOLD', 'public known/predicted bulk parsed; no-record applies to this download only',pages_complete=has))
        ledgers.append(m.ledger(q,'HERB','MANUAL_QUERY_REQUIRED_NO_PUBLIC_BATCH_ENDPOINT_VERIFIED' if has else 'IDENTITY_HOLD', 'public page/JS captured; ingredient-target retrieval not completed; download control disabled'))
    m.writetable(OUT/'prior_identity_unresolved.tsv',unresolved,['CID','token','reason'])
    m.writetable(OUT/'ledger_prior_sources.tsv',ledgers)
    return records,ledgers


def disease_sets():
    robust=m.readtable(ROOT/cfg['inputs']['g1_244'])
    full=m.readtable(ROOT/cfg['inputs']['g1_full'])
    main={r['gene']:r for r in robust}
    sig={r['gene']:r for r in full if m.number(r.get('meta_FDR')) is not None and float(r['meta_FDR'])<.05}
    assert len(main)==244 and len(sig)==413, (len(main),len(sig))
    sets={'G1_robust244':main,'G1_significant413_exploratory':sig}; contexts={
        'G1_robust244':'human AD microglia; donor-aware RNA association',
        'G1_significant413_exploratory':'human AD microglia; full meta significant set, not robust244'}
    g2=yaml.safe_load((ROOT/'configs/g2_v3.2_execution.yaml').read_text(encoding='utf-8'))
    sets['G2_glycolysis_model_exploratory']={g:{} for g in g2['analysis']['glycolysis_target_set']['genes']}
    contexts['G2_glycolysis_model_exploratory']='5XFAD mouse microglia H4K12la target list; human G2c closed'
    g3=yaml.safe_load((ROOT/'configs/g3_v3.2_execution.yaml').read_text(encoding='utf-8'))
    for name,info in g3['modules'].items():
        key='G3_'+name+'_exploratory'; sets[key]={g:{} for g in info['genes']}; contexts[key]='A-beta iPSC microglia model; single source donor; RNA proxy'
    g5=m.readtable(ROOT/cfg['inputs']['g5_orthologs'])
    for r in g5:
        gene=r.get('human_one_to_one_symbol','')
        if not gene: continue
        key='G5_'+r['module']+'_ortholog_exploratory'; sets.setdefault(key,{})[gene]=r
        contexts[key]='APP/PS1-21 mouse astrocytes; frozen human one-to-one ortholog; not human astrocyte validation'
    tree=ast.parse((ROOT/'src/32_run_g4_v3_2.py').read_text(encoding='utf-8'))
    modules=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='MODULES' for t in n.targets))
    modules['M7_FerroAging']=[r['gene_symbol'].upper() for r in m.readtable(ROOT/'results/archive/2026-09-20-pre-v3-freeze/results/tables/D5_ferroaging_reference.tsv')]
    for name,genes in modules.items():
        key='G4_'+name+'_exploratory'; sets[key]={g:{} for g in genes if g!='DMT1'}
        contexts[key]='A-beta human iPSC microglia timecourse; original G4 frozen program; RNA proxy, not flux/death'
    return sets,contexts


def historical_disposition(evidence):
    baseline=m.readtable(ROOT/cfg['inputs']['baseline_023'])
    genes={r['gene'] for r in baseline}; old=m.readtable(ROOT/cfg['inputs']['priors']); result=[]
    hist=m.readtable(ROOT/cfg['inputs']['historical_assays'])
    similarities=[{**r,'tier':'S','historical_symbol_status':'unvalidated historical mapping; not a verified candidate target',
                   'disposition':'retain historical similar-compound hypothesis; exclude candidate experimental intersections'}
                  for r in hist if r['Evidence_Source']=='ChEMBL']
    m.writetable(OUT/'historical_similar_structure_hypotheses.tsv',similarities)
    bypair=defaultdict(list)
    for r in evidence:
        for cid in str(r.get('cids','')).split(';'):
            bypair[(cid,r.get('target_symbol',''))].append(r)
    for i,r in enumerate(old):
        if r['Target_Symbol'] not in genes: continue
        match=bypair[(r['PubChem_CID'],r['Target_Symbol'])]
        source=r.get('Final_Target_Source_DB') or r.get('Source','')
        similarity='similarity' in r.get('Target_DB_Notes','').lower()
        outcomes={x.get('activity_outcome','').lower() for x in match}
        inactive_only=bool(outcomes) and outcomes.issubset({'inactive','inconclusive','unspecified'})
        result.append({'historical_source_row':i+2,'cid':r['PubChem_CID'],'gene':r['Target_Symbol'],'source':source,
          'disposition':'WITHDRAW_EXPERIMENTAL_ATTRIBUTION_RETAIN_S' if similarity else 'WITHDRAW_POSITIVE_TARGET_ATTRIBUTION_RETAIN_TESTED_OUTCOME' if inactive_only else 'RETAIN_H_UNVERIFIED_ASSOCIATION',
          'reason':'historical similarity/name mapping does not establish exact candidate target; new records reviewed separately' if similarity else 'public source recheck shows only inactive/inconclusive outcomes for this CID/gene' if inactive_only else 'historical imported mapping lacks original positive target-specific assay',
          'new_same_CID_gene_records':len(match),'new_E1_records':sum(x['tier']=='E1' and positive(x) for x in match),
          'new_outcomes':';'.join(sorted({x.get('activity_outcome','') for x in match}))})
    bdb=[r for r in evidence if r['source']=='BindingDB' and '8768' in r.get('cids','').split(';')]
    result.append({'historical_source_row':'020 BindingDB CID8768','cid':'8768','gene':'PTGS2','source':'BindingDB',
      'disposition':'E2_OR_UNVERIFIED_NO_E1_RELEASE','reason':'local 19700 lacks verified unit/original assay; REST name-only results cannot establish stable target identity',
      'new_same_CID_gene_records':len(bdb),'new_E1_records':sum(r['tier']=='E1' and r.get('target_symbol')=='PTGS2' for r in bdb)})
    m.writetable(OUT/'historical_8_and_PTGS2_disposition.tsv',result)
    m.writetable(OUT/'historical_8_baseline_preserved.tsv',baseline)
    return result


def main():
    queries=m.readtable(OUT/'compound_identity_audit.tsv')
    candidates=m.readtable(ROOT/cfg['inputs']['candidates'])
    herbs_bycid={r['PubChem_CID']:r['Herbs'] for r in candidates}
    for q in queries:
        q['herbs']=';'.join(sorted({herb.strip() for cid in q['confirmed_cids'].split(';') for herb in herbs_bycid.get(cid,'').split(';') if herb.strip()}))
    m.writetable(OUT/'compound_identity_confirmed_attribution.tsv',queries)
    priors,ledgers=prior_tables(queries); evidence=[]
    for source in ('ChEMBL','BindingDB','PubChem'):
        p=OUT/f'evidence_{source}.json'
        if p.exists(): evidence.extend(json.loads(p.read_text(encoding='utf-8')))
        p=OUT/f'ledger_{source}.tsv'; have={r['compound_key']:r for r in m.readtable(p)} if p.exists() else {}
        ledgers.extend(have.get(q['compound_key']) or m.ledger(q,source,'RUNNING_NOT_YET_QUERIED','live fetch still pending') for q in queries)
    evidence.extend(priors)
    approved={r['symbol']:r for r in m.readtable(RAW/'HGNC/hgnc_complete_set.tsv') if r['status']=='Approved'}
    comments={}
    for path in (RAW/'ChEMBL').glob('*.json'):
        meta=json.loads(path.read_text(encoding='utf-8'))
        if '/activity.json' not in meta.get('url','') or meta.get('status')!=200: continue
        payload=json.loads((ROOT/meta['body_path']).read_text(encoding='utf-8'))
        comments.update({str(a['activity_id']):a.get('activity_comment') or '' for a in payload.get('activities',[])})
    unique={}
    for r in evidence:
        r['fetch_provisional_tier']=r.get('tier','')
        if r.get('source')=='ChEMBL':
            r['source_activity_comment']=comments.get(r['activity_id'],r.get('source_activity_comment',''))
            r['activity_outcome']=m.chembl_outcome(r['source_activity_comment'])
        if r.get('source')=='BindingDB' and r.get('actual_smiles') and not r.get('identity_exact'):
            r['similarity_only']=True
        r['herbs']=';'.join(sorted({herb.strip() for cid in r.get('cids','').split(';') for herb in herbs_bycid.get(cid,'').split(';') if herb.strip()})) if r.get('candidate_identity_verified') is not False else ''
        mapping=approved.get(r.get('target_symbol',''),{})
        r['target_hgnc_id']=mapping.get('hgnc_id',''); r['target_gene_id']=mapping.get('entrez_id','')
        r['tier'],r['tier_reason']=m.classify(r)
        key=m.evidence_key(r)
        if key not in unique: unique[key]=r
        else:
            unique[key]['duplicate_sources']=';'.join(sorted(set(unique[key].get('duplicate_sources',unique[key]['source']).split(';'))|{r['source']}))
    evidence=list(unique.values())
    sets,contexts=disease_sets(); reviews={}
    rp=OUT/'manual_assay_adjudication.tsv'
    if rp.exists(): reviews={r['assay_id']:r for r in m.readtable(rp)}
    for r in evidence:
        review=reviews.get(r.get('assay_id',''),{})
        if r['tier']=='E1' and review.get('decision')=='DOWNGRADE_E2':
            r['tier']='E2'; r['tier_reason']='manual original assay review: '+review['rationale']; r['manual_source_review']='DOWNGRADE_E2'
    cards={}
    for r in evidence:
        if r['tier']=='E1' and positive(r) and r['target_symbol'] in sets['G1_robust244']:
            aid=r['assay_id']; review=reviews.get(aid,{})
            r['manual_source_review']=review.get('decision','pending')
            r['review_note']=review.get('rationale','')
            cards[aid]={'assay_id':aid,'gene':r['target_symbol'],'document':r.get('document',''),
              'description':r.get('assay_description',''),'actual_compound_ids':';'.join(sorted({x['actual_compound_id'] for x in evidence if x.get('assay_id')==aid and x['tier']=='E1' and positive(x)})),
              'source_url':'https://www.ebi.ac.uk/chembl/api/data/assay/'+aid+'.json',
              'decision':review.get('decision','pending'),'rationale':review.get('rationale','')}
    m.writetable(OUT/'G1_E1_source_review_cards.tsv',list(cards.values()))
    connections=[]
    for r in evidence:
        tier=r['tier']
        if tier in ('E1','E2') and not positive(r): continue
        if tier not in ('E1','E2','P','H'): continue
        if not r.get('identity_exact') or r.get('candidate_identity_verified') is False: continue
        gene=r.get('target_symbol','')
        for key,genes in sets.items():
            if gene not in genes: continue
            disease=genes[gene]; released=tier=='E1' and r.get('manual_source_review')=='RETAIN' and key=='G1_robust244'
            connections.append({'gene':gene,'disease_set':key,'context':contexts[key],'disease_effect':disease.get('meta_logFC',disease.get('G1_meta_logFC','')),
              'disease_FDR':disease.get('meta_FDR',disease.get('G1_meta_FDR','')),'compound_key':r['compound_key'],'cids':r['cids'],'herbs':r['herbs'],
              'tier':tier,'source':r['source'],'target_stable_id':r.get('target_stable_id',''),'target_hgnc_id':r.get('target_hgnc_id',''),'target_gene_id':r.get('target_gene_id',''),'assay_id':r.get('assay_id',''),'activity_id':r.get('activity_id',''),
              'value':r.get('value',''),'units':r.get('units',''),'relation':r.get('relation',''),'document':r.get('document',''),
              'module_effect':disease.get('G5_module_logFC',''),'module_FDR':disease.get('G5_module_FDR_BH_7',''),'module_gate':disease.get('G5_module_gate',''),
              'effect_scope':'gene-level human G1 association' if key.startswith('G1_') else 'prespecified model program; module-level fields if available; no human gene-effect claim',
              'manual_review':r.get('manual_source_review','not_applicable'),'released_main_E1':released,'treatment_direction':'undetermined'})
    mainrows=[r for r in connections if r['disease_set']=='G1_robust244']
    m.writetable(OUT/'G1_244_layered_connections.tsv',mainrows)
    m.writetable(OUT/'exploratory_disease_connections.tsv',[r for r in connections if r['disease_set']!='G1_robust244'])
    bykey={q['compound_key']:q for q in queries}
    for r in ledgers:
        if r['source']=='BindingDB' and r['status']=='RETRY_EXHAUSTED' and 'JSON parse failure' in r['detail']:
            r['fetch_recorded_status']=r['status']; r['status']='PARSE_FAILED'
        q=bykey[r['compound_key']]
        r.update(smiles=q['smiles'],inchikey=q['inchikey'],confirmed_cids=q['confirmed_cids'],conflicting_cids=q['conflicting_cids'],herbs=q['herbs'])
        if r['source'] in ('PubChem','BATMAN','HERB'): r['query_key']=q['confirmed_cids']
    m.writetable(OUT/'query_ledger.tsv',ledgers)
    m.writetable(OUT/'pending_query_package.tsv',[r for r in ledgers if not str(r['status']).startswith('SUCCESS') and r['status']!='IDENTITY_HOLD'])
    m.writetable(OUT/'evidence_layered.tsv',evidence)
    historical=historical_disposition(evidence)
    frozen=m.readtable(OUT/'source_sha256.tsv'); checks=[]
    for r in frozen:
        if r['role']=='code': continue
        current=m.filehash(ROOT/r['path']); checks.append({'path':r['path'],'role':r['role'],'frozen_sha256':r['sha256'],'current_sha256':current,'match':current==r['sha256']})
    m.writetable(OUT/'frozen_input_verification.tsv',checks)
    coverage=[]
    for source in ('ChEMBL','BindingDB','PubChem','BATMAN','HERB'):
        rows=[r for r in ledgers if r['source']==source]; status=Counter(r['status'] for r in rows)
        for k,v in status.items(): coverage.append({'source':source,'status':k,'structures':v,'total_structures':len(queries)})
    m.writetable(OUT/'coverage_by_source.tsv',coverage)
    phase_coverage=Counter((r['source'],r['phase'],r['status']) for r in ledgers)
    m.writetable(OUT/'coverage_by_source_phase.tsv',[{'source':s,'phase':p,'status':t,'structures':n} for (s,p,t),n in sorted(phase_coverage.items())])
    summary=[]
    for tier in ('E1','E2','P','S','H'):
        rows=[r for r in mainrows if r['tier']==tier]
        summary.append({'tier':tier,'G1_genes':len({r['gene'] for r in rows}),'compound_target_pairs':len({(r['compound_key'],r['gene']) for r in rows}),'source_records':len(rows),
          'released_G1_genes':len({r['gene'] for r in rows if r['released_main_E1']}),'genes':';'.join(sorted({r['gene'] for r in rows}))})
    m.writetable(OUT/'G1_layer_summary.tsv',summary)
    overall=[]
    for tier in ('E1','E2','P','S','H'):
        rows=[r for r in evidence if r['tier']==tier]
        eligible=[r for r in rows if positive(r)] if tier in ('E1','E2') else [r for r in rows if r.get('identity_exact') and r.get('candidate_identity_verified') is not False and r.get('species')=='9606' and r.get('target_stable_id') and r.get('target_symbol')] if tier in ('P','H') else []
        overall.append({'tier':tier,'all_source_records':len(rows),'eligible_human_candidate_genes':len({r['target_symbol'] for r in eligible}),
          'eligible_compound_target_pairs':len({(r['compound_key'],r['target_symbol']) for r in eligible}),
          'eligible_source_records':len(eligible),'distinct_source_assay_IDs':len({(r['source'],r['assay_id']) for r in eligible if r.get('assay_id')}),
          'interpretation':'rule-selected records; original G1 assay reviews separately; assay IDs not independent experiments'})
    m.writetable(OUT/'all_target_layer_summary.tsv',overall)
    targets=[]
    for gene in sorted({r['gene'] for r in mainrows if r['released_main_E1']}):
        rows=[r for r in mainrows if r['gene']==gene and r['released_main_E1']]
        targets.append({'gene':gene,'G1_meta_logFC':sets['G1_robust244'][gene]['meta_logFC'],'G1_meta_FDR':sets['G1_robust244'][gene]['meta_FDR'],
          'candidate_CIDs':';'.join(sorted({cid for r in rows for cid in r['cids'].split(';') if cid},key=int)),
          'herbs_directory_level_only':';'.join(sorted({h for r in rows for h in r['herbs'].split(';') if h})),
          'compound_target_pairs':len({r['compound_key'] for r in rows}),'source_records':len(rows),'source_assay_IDs':len({r['assay_id'] for r in rows}),
          'experimental_context':'protein-specific aggregation modulation; direct binding not established' if gene=='SNCA' else 'target-specific biochemical enzyme inhibition',
          'interpretation':'public candidate compound assay connection; not measured formula presence or efficacy'})
    m.writetable(OUT/'G1_E1_released_target_summary.tsv',targets)
    success=sum(str(r['status']).startswith('SUCCESS') for r in ledgers)
    state={'status':'PARTIAL_SOURCE_BLOCKED' if success<len(ledgers) else 'COMPLETE','C0':True,'reopen_017':False,'candidate_CIDs':2276,
      'all_available_source_runs_finished':not any(r['status']=='RUNNING_NOT_YET_QUERIED' for r in ledgers),
      'query_structures':len(queries),'phase_structures':dict(Counter(q['phase'] for q in queries)),
      'CID_structure_conflict_queries':sum(bool(q['conflicting_cids']) for q in queries),'unconfirmed_queries':sum(not q['confirmed_cids'] for q in queries),
      'query_rows':len(ledgers),'query_rows_with_explicit_status':sum(bool(r['status']) for r in ledgers),
      'all_biological_inputs_hash_unchanged':all(r['match'] for r in checks),'source_coverage':coverage,
      'evidence_counts':dict(Counter(r['tier'] for r in evidence)),'G1_summary':summary,
      'all_target_layer_summary':overall,'historical_similarity_records_preserved':sum(r['Evidence_Source']=='ChEMBL' for r in m.readtable(ROOT/cfg['inputs']['historical_assays'])),
      'E1_G1_assays':len(cards),'E1_G1_assays_review_pending':sum(c['decision']=='pending' for c in cards.values()),
      'manual_original_assay_reviews':len(reviews),'manual_downgraded_assays':sum(r['decision']=='DOWNGRADE_E2' for r in reviews.values()),
      'cross_source_independent_experiment_count_not_established':True,'source_record_count_not_independent_experiments':True,
      'historical_disposition_rows':len(historical),'drug_efficacy_claim':False,'disease_direction_drug_direction_conversion':False}
    m.writejson(OUT/'machine_audit.json',state)
    runnote='本轮可执行来源补查已结束，受阻来源仍未完成。' if state['all_available_source_runs_finished'] else '抓取进行中，本报告为阶段结果。'
    lines=['# 固定处方公共靶点补查024报告','',f"状态：`{state['status']}`；药物C0，017关闭项保持。{runnote}",'',
      '## 输入与身份','',f"2276个候选CID均已取得公开结构。旧结构1299个（优先335，其余964）；保留差异并补取后共{len(queries)}个查询结构。{state['unconfirmed_queries']}个旧目录结构无确认CID；{state['CID_structure_conflict_queries']}个结构条目含冲突CID，不能视为已确认候选。",'',
      '## 来源覆盖','', '|来源|状态|结构数|','|---|---|---|']
    lines.extend(f"|{r['source']}|{r['status']}|{r['structures']}|" for r in coverage)
    lines += ['','查询状态包含失败/未完成/身份暂停；这些均不是无靶点。HERB公开页面已取得，但未核验成分靶点批量接口，保留逐项待补。BATMAN公开文件完成解析，文件范围无命中不代表无作用。','',
      '## 总体靶点覆盖','', '|层|规则合格人类候选靶点基因|成分靶点对|来源记录|','|---|---:|---:|---:|']
    lines.extend(f"|{r['tier']}|{r['eligible_human_candidate_genes']}|{r['eligible_compound_target_pairs']}|{r['eligible_source_records']}|" for r in overall)
    lines += ['','全体为规则筛选/分层结果，E1人工原始描述核验范围限新增G1交集；来源记录/assay编号不是独立实验数。旧534条相似结构假说单列historical_similar_structure_hypotheses.tsv，保留原字段及不可信靶点映射标识，不并入本轮候选实验交集。','',
      '## G1稳健244分层交集','', '|层|基因|成分靶点对|来源记录|人工放行基因|','|---|---:|---:|---:|---:|']
    lines.extend(f"|{r['tier']}|{r['G1_genes']}|{r['compound_target_pairs']}|{r['source_records']}|{r['released_G1_genes']}|" for r in summary)
    lines += ['', 'E1候选计数在人工核验前不放行；E2、P、H不是主实验验证。Inactive、不明确活性、身份未确认和下限删失记录保留原文，但不单凭其构造阳性实验交集。无效/速率终点单位不能混入浓度效力。','',
      '已放行主交集：'+ '、'.join(r['gene'] for r in targets)+'。核验依据是官方原始实验描述及结构/activity/target元数据；不宣称全文均已取得。SNCA为蛋白聚集调节，其余为靶点特异酶学抑制记录；都不能直接推出AD治疗方向。逐基因成分CID、药味目录关系、G1效应/FDR与实验背景见G1_E1_released_target_summary.tsv。全体E1是规则筛选记录，人工核验范围限新增G1交集。','',
      '## 历史与结论','', '原8个交集仍是历史映射先验；逐源处置见historical_8_and_PTGS2_disposition.tsv。CID8768—PTGS2的旧19700记录未因新抓取自动获得单位/原始实验放行。新增其他精确化合物对ALOX15B等的记录，不追认原相似化合物。', '',
      'G1全显著413及G2–G5分别另列探索接口，保留物种/细胞/模型限制。当前仅说明候选和疾病基因存在分层连接，不能判断整方有效、治疗方向、暴露可及性或状态恢复。缺数据或数据库无记录不代表药物无效。','',
      '## 复现与交付','', '参数configs/fixed_formula_target_supplement_024.json；脚本113/114/115；复用TCM_select现有环境及RDKit/requests/PyYAML，无依赖安装。原始响应按SHA核验，生物学输入核对结果见frozen_input_verification.tsv。代码修改与冻结输入区别见code_revision_log.json。来源记录数不等于独立实验数；跨库只有共同原始测量ID明确才去重。', '',
      '表：compound_identity_audit.tsv、query_ledger.tsv、evidence_layered.tsv、G1_244_layered_connections.tsv、exploratory_disease_connections.tsv、G1_E1_source_review_cards.tsv、pending_query_package.tsv及machine_audit.json。完整SOP见docs/protocol/2026-09-30-fixed-formula-target-supplement-SOP-024.md。最终验收、DVC与WikiSkill记录另附。']
    (OUT/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps({k:state[k] for k in ('status','query_structures','evidence_counts','G1_summary','E1_G1_assays_review_pending')},ensure_ascii=False),flush=True)


if __name__=='__main__': main()

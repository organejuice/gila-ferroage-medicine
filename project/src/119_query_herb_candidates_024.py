"""Sequential public HERB UI API queries with exact candidate identity gates.

The public UI's search_api and detail_api are used without disabled download
controls. Database cross references stay H; no literature/assay promotion.
"""
import importlib.util
import json
from collections import defaultdict
from pathlib import Path

spec=importlib.util.spec_from_file_location('continuation024',Path(__file__).with_name('118_continue_formula_supplement_024.py'))
c=importlib.util.module_from_spec(spec); spec.loader.exec_module(c)
m=c.m


def cell(value):
    if isinstance(value,dict): return str(value.get('title',''))
    if isinstance(value,list): return ';'.join(cell(v) for v in value)
    return str(value or '')


def table(data):
    if not isinstance(data,list) or not data: return []
    headers=[cell(h) for h in data[0]]
    if len(headers)!=len(set(headers)): raise ValueError('duplicate HERB table headers')
    rows=[]
    for row in data[1:]:
        if len(row)!=len(headers): raise ValueError('HERB table row width mismatch')
        rows.append(dict(zip(headers,[cell(v) for v in row])))
    return rows


def identity(summary,q,cid):
    rows=table(summary)
    if len(rows)!=1: return False,'summary not a unique ingredient'
    row=rows[0]
    if row.get('PubChem id','')!=str(cid): return False,'HERB PubChem id differs/missing'
    if m.canonical(row.get('Molecule smile',''))!=q['smiles']: return False,'exact stereo structure differs/missing'
    return True,'same CID and exact canonical isomeric SMILES'


def main():
    queries=m.readtable(c.OLD/'compound_identity_audit.tsv')
    bycid={cid:q for q in queries for cid in q['confirmed_cids'].split(';') if cid}
    hgnc=m.readtable(m.ROOT/'data/raw/2026-09-30-fixed-formula-target-supplement-024/HGNC/hgnc_complete_set.tsv')
    approved={r['symbol']:r for r in hgnc if r['status']=='Approved'}
    ledgerpath=c.OUT/'HERB_CID_query_ledger.tsv'
    old=m.readtable(ledgerpath) if ledgerpath.exists() else []
    done={r['cid']:r for r in old if r['status'].startswith('SUCCESS')}
    rows=list(done.values())
    evidencepath=c.OUT/'evidence_HERB.json'
    evidence=json.loads(evidencepath.read_text()) if evidencepath.exists() else []
    idpath=c.OUT/'HERB_candidate_identity.tsv'
    identities=m.readtable(idpath) if idpath.exists() else []
    failures=0
    for index,(cid,q) in enumerate(sorted(bycid.items(),key=lambda x:int(x[0])),1):
        if cid in done: continue
        if failures>=3:
            rows.append(dict(cid=cid,compound_key=q['compound_key'],status='SOURCE_UNAVAILABLE_NOT_QUERIED',detail='circuit open after three exhausted queries')); continue
        try:
            data,meta=c.herb_post(dict(func_name='search_api',label='Ingredient',keyword=cid))
            if not isinstance(data,dict) or 'res_data' not in data or data.get('error'): raise ValueError('HERB search schema/error')
            found=table(data['res_data'])
            matching=[]; start=len(evidence)
            for ingredient in found:
                iid=ingredient['Ingredient id']
                detail,dmeta=c.herb_post(dict(func_name='detail_api',label='Ingredient',key_id=iid,v=iid))
                if not isinstance(detail,dict) or 'summary' not in detail: raise ValueError('HERB detail schema')
                valid,reason=identity(detail['summary'],q,cid)
                identities.append(dict(cid=cid,compound_key=q['compound_key'],HERB_id=iid,exact=valid,reason=reason,response_sha256=dmeta['sha256']))
                if not valid: continue
                matching.append(iid)
                for section in ['ingredient_target','drug_paper_target']:
                    for target in table(detail.get(section,[])):
                        symbol=target.get('Target name',''); mapping=approved.get(symbol,{})
                        r=m.base_evidence(q,'HERB')
                        r.update(identity_exact=True,actual_smiles=q['smiles'],actual_compound_id=iid,
                            target_symbol=mapping.get('symbol',''),target_name=symbol,
                            target_stable_id='HGNC:'+mapping['hgnc_id'].split(':')[-1] if mapping else target.get('Target id',''),
                            target_gene_id=mapping.get('entrez_id',''),species='9606' if mapping else '',
                            association_only=True,activity_outcome='database_association_unverified',
                            activity_id=':'.join([iid,section,target.get('Target id',''),target.get('Paper id','')]),
                            assay_id='',document=target.get('PubMed id',''),HERB_target_id=target.get('Target id',''),
                            HERB_section=section,HERB_database_sources=target.get('Database sources',''),
                            HERB_paper_id=target.get('Paper id',''),source_url=dmeta['url'],response_sha256=dmeta['sha256'])
                        r['tier'],r['tier_reason']=m.classify(r); evidence.append(r)
            rows.append(dict(cid=cid,compound_key=q['compound_key'],status='SUCCESS_WITH_RECORDS' if len(evidence)>start else 'SUCCESS_NO_RECORDS',
                search_matches=len(found),exact_ingredient_ids=';'.join(matching),records=len(evidence)-start,
                detail='scope: public CID search; nonexact hits excluded; no target inference from absence',response_sha256=meta['sha256'],updated_at=m.now()))
            failures=0
        except (m.SourceError,ValueError,KeyError) as exc:
            failures+=1
            rows.append(dict(cid=cid,compound_key=q['compound_key'],status='QUERY_FAILED',detail=str(exc)[:500],updated_at=m.now()))
        # Checkpoint every candidate; recover successes only and retain original raw.
        m.writetable(ledgerpath,rows)
        m.writetable(idpath,identities,['cid','compound_key','HERB_id','exact','reason','response_sha256'])
        m.writejson(evidencepath,evidence)
        if index%20==0: print(f'HERB {index}/{len(bycid)} records={len(evidence)}',flush=True)
    m.writetable(ledgerpath,rows)
    m.writejson(evidencepath,evidence)
    print('HERB finished',len(rows),len(evidence),flush=True)


if __name__=='__main__': main()

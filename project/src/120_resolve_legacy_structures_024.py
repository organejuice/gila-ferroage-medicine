"""Trace unresolved legacy structures without enlarging frozen CID candidates."""
import importlib.util
import json
from pathlib import Path

spec=importlib.util.spec_from_file_location('c024',Path(__file__).with_name('118_continue_formula_supplement_024.py'))
c=importlib.util.module_from_spec(spec); spec.loader.exec_module(c)
m=c.m


def main():
    queries=m.readtable(c.OLD/'compound_identity_audit.tsv')
    frozen={cid for q in queries for cid in q['confirmed_cids'].split(';') if cid}
    pending=[q for q in queries if not q['confirmed_cids']]
    client=m.Client('PubChem_legacy_identity',c.CFG,c.RAW)
    path=c.OUT/'legacy_structure_identity_followup.tsv'
    old=m.readtable(path) if path.exists() else []
    done={r['compound_key']:r for r in old if r['status'].startswith('SUCCESS')}
    rows=list(done.values()); failures=0
    for index,q in enumerate(pending,1):
        if q['compound_key'] in done: continue
        if failures>=3:
            rows.append(dict(compound_key=q['compound_key'],status='SOURCE_UNAVAILABLE_NOT_QUERIED',candidate_release=False)); continue
        try:
            data,meta=client.json(c.CFG['sources']['PubChem']+'compound/fastidentity/smiles/cids/JSON',dict(smiles=q['smiles'],identity_type='same_stereo_isotope'),allow404=True)
            cids=[str(cid) for cid in (data or {}).get('IdentifierList',{}).get('CID',[])]
            verified=[]; standardization_mismatch=[]
            if cids:
                props,pmeta=client.json(c.CFG['sources']['PubChem']+'compound/cid/'+','.join(cids)+'/property/SMILES,InChIKey/JSON')
                for p in props.get('PropertyTable',{}).get('Properties',[]):
                    cid=str(p['CID'])
                    if m.canonical(p.get('SMILES') or p.get('IsomericSMILES',''))==q['smiles']: verified.append(cid)
                    else: standardization_mismatch.append(cid)
            release=set(verified)&frozen
            rows.append(dict(compound_key=q['compound_key'],legacy_cids=q['cids'],search_cids=';'.join(cids),exact_verified_cids=';'.join(verified),
                mismatched_standardized_cids=';'.join(standardization_mismatch),candidate_release=bool(release),frozen_candidate_matches=';'.join(sorted(release)),
                status='SUCCESS_WITH_EXACT_CID' if verified else 'SUCCESS_NO_EXACT_CID',response_sha256=meta['sha256'],
                disposition='outside frozen candidates; preserve provenance hold' if verified and not release else 'frozen candidate identity confirmed' if release else 'identity hold retained',updated_at=m.now()))
            failures=0
        except (m.SourceError,ValueError,KeyError) as exc:
            failures+=1; rows.append(dict(compound_key=q['compound_key'],status='QUERY_FAILED',detail=str(exc)[:500],candidate_release=False))
        m.writetable(path,rows)
        if index%20==0: print('Legacy identity',index,'/',len(pending),flush=True)
    m.writetable(path,rows)
    print('Legacy identity finished',len(rows),flush=True)


if __name__=='__main__': main()

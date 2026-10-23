"""Prospective continuation of SOP024; never overwrite the September snapshot."""
import importlib.util
import json
import sys
import shutil
import time
import requests
from pathlib import Path

spec = importlib.util.spec_from_file_location('s024', Path(__file__).with_name('113_supplement_formula_targets_024.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)
ROOT = m.ROOT
RUN = '2026-10-03-fixed-formula-target-supplement-024-continuation'
OLD = ROOT/'results/2026-09-30-fixed-formula-target-supplement-024'
OUT = ROOT/'results'/RUN
RAW = ROOT/'data/raw'/RUN
CFG = json.loads(m.CONFIG.read_text(encoding='utf-8'))
CONT = json.loads((ROOT/'configs/fixed_formula_supplement_024_continuation.json').read_text(encoding='utf-8'))
CFG['http'] = CONT['http']

HERB_SESSION=requests.Session()
HERB_LAST=0.0


def herb_post(payload):
    global HERB_LAST
    url=CONT['HERB']['url']
    request_text=json.dumps(payload,sort_keys=True,ensure_ascii=False)
    key=m.digest((url+request_text).encode())
    folder=RAW/'HERB'; folder.mkdir(parents=True,exist_ok=True)
    bodypath=folder/(key+'.body'); metapath=folder/(key+'.json')
    if metapath.exists():
        meta=json.loads(metapath.read_text())
        if bodypath.exists() and m.filehash(bodypath)==meta['sha256']:
            return json.loads(bodypath.read_bytes()),meta
        raise ValueError('cached HERB hash mismatch')
    attempts=[]
    for delay in [0]+CFG['http']['backoff_seconds']:
        time.sleep(delay)
        time.sleep(max(0,CFG['http']['interval_seconds']-(time.monotonic()-HERB_LAST)))
        HERB_LAST=time.monotonic()
        try:
            response=HERB_SESSION.post(url,json=payload,timeout=CFG['http']['timeout_seconds'])
            attempts.append({'status':response.status_code,'bytes':len(response.content),'retrieved_at':m.now()})
            if response.status_code==200:
                data=response.json()
                bodypath.write_bytes(response.content)
                meta={'source':'HERB','url':url,'method':'POST','request':payload,'status':200,'sha256':m.digest(response.content),'bytes':len(response.content),'body_path':str(bodypath.relative_to(ROOT)),'retrieved_at':m.now(),'attempts':attempts}
                m.writejson(metapath,meta)
                return data,meta
            if response.status_code not in [408,429,500,502,503,504]: break
            if response.headers.get('Retry-After','').isdigit(): time.sleep(int(response.headers['Retry-After']))
        except (requests.RequestException,ValueError) as exc:
            attempts.append({'error':str(exc),'retrieved_at':m.now()})
    m.writejson(folder/(key+'.failure.json'),{'url':url,'request':payload,'attempts':attempts})
    raise m.SourceError('HERB POST failed: '+str(attempts))


def herb_probe():
    records=[]
    for payload in [dict(func_name='browse_api',label='Ingredient',page=1,page_size=50),dict(func_name='search_api',label='Ingredient',keyword='quercetin')]:
        data,meta=herb_post(payload)
        print(json.dumps(data,ensure_ascii=False)[:6500],flush=True)
        records.append(meta)
    m.writejson(OUT/'HERB_API_validation.json',records)


def capture(source, url, params=None):
    client = m.Client(source, CFG, RAW)
    try:
        body, meta = client.get(url, params)
        return dict(meta, outcome='CAPTURED', preview=body[:160].decode('utf-8',errors='replace'))
    except m.SourceError as exc:
        return dict(source=source, url=url, outcome='SOURCE_BLOCKED', detail=str(exc), retrieved_at=m.now())


def probe():
    OUT.mkdir(parents=True, exist_ok=True)
    queries = {r['compound_key']:r for r in m.readtable(OLD/'compound_identity_audit.tsv')}
    failed = [r for r in m.readtable(OLD/'ledger_BindingDB.tsv') if not r['status'].startswith('SUCCESS')]
    m.writetable(OUT/'baseline_pending_BindingDB.tsv', failed)
    records=[]
    for row in failed[:3]:
        q=queries[row['compound_key']]
        rec=capture('BindingDB', CFG['sources']['BindingDB'], {'smiles':q['smiles'],'cutoff':'1.0','response':'application/json'})
        rec['compound_key']=q['compound_key']
        if rec['outcome']=='CAPTURED':
            try:
                data=json.loads((ROOT/rec['body_path']).read_bytes())
                rec['schema_valid']=bool(data.get('getLindsByUniprotResponse'))
            except (ValueError,AttributeError): rec['schema_valid']=False
        records.append(rec)
        m.writetable(OUT/'source_probe.tsv',records)
        print('BindingDB',q['compound_key'],rec['outcome'],rec.get('schema_valid'),flush=True)
    for source,url in [('HERB','http://herb.ac.cn/v2'),('HERB','http://herb.ac.cn/static/umi.js?version=762')]:
        rec=capture(source,url); records.append(rec)
        m.writetable(OUT/'source_probe.tsv',records)
        print(source,rec['outcome'],rec.get('bytes'),flush=True)
    for row in m.readtable(OLD/'publication_supplement_access.tsv'):
        rec=capture('Publication',row['url']); records.append(rec)
        m.writetable(OUT/'source_probe.tsv',records)
        print('Publication',rec['outcome'],flush=True)
    m.writejson(OUT/'probe_audit.json',{'run_id':RUN,'retrieved_at':m.now(),'records':records})


def bindingdb():
    OUT.mkdir(parents=True,exist_ok=True)
    for name in ['ledger_BindingDB.tsv','evidence_BindingDB.json','compound_identity_audit.tsv']:
        if not (OUT/name).exists(): shutil.copyfile(OLD/name,OUT/name)
    hgnc=ROOT/'data/raw/2026-09-30-fixed-formula-target-supplement-024/HGNC/hgnc_complete_set.tsv'
    _,uniprot=m.load_hgnc(hgnc)
    m.bindingdb_fetch(CFG,OUT,RAW,m.readtable(OUT/'compound_identity_audit.tsv'),uniprot)
    rows=m.readtable(OUT/'ledger_BindingDB.tsv')
    from collections import Counter
    print(json.dumps(Counter(r['status'] for r in rows)),flush=True)


def bindingdb_resume():
    """Re-open a circuit only after cooldown and a successful public control."""
    from collections import Counter
    queries=m.readtable(OLD/'compound_identity_audit.tsv')
    _,uniprot=m.load_hgnc(ROOT/'data/raw/2026-09-30-fixed-formula-target-supplement-024/HGNC/hgnc_complete_set.tsv')
    historypath=OUT/'BindingDB_resume_rounds.json'
    history=json.loads(historypath.read_text()) if historypath.exists() else []
    no_progress=0
    for roundid in range(len(history)+1,1000):
        rows=m.readtable(OUT/'ledger_BindingDB.tsv')
        before=sum(r['status'].startswith('SUCCESS') for r in rows)
        if before==len(queries): break
        time.sleep(CONT['BindingDB']['cooldown_seconds'])
        client=m.Client('BindingDB_health_round'+str(roundid),CFG,RAW)
        try:
            control,meta=client.json(CFG['sources']['BindingDB'],dict(smiles=CONT['BindingDB']['health_control_smiles'],cutoff='1.0',response='application/json'))
            if not isinstance(control,dict) or not control.get('getLindsByUniprotResponse'): raise m.SourceError('health schema failed')
        except m.SourceError as exc:
            history.append(dict(round=roundid,status='STOP_SOURCE_UNHEALTHY',reason=str(exc),time=m.now()))
            m.writejson(historypath,history); break
        roundraw=RAW/('round_'+str(roundid).zfill(3))
        m.bindingdb_fetch(CFG,OUT,roundraw,queries,uniprot)
        rows=m.readtable(OUT/'ledger_BindingDB.tsv')
        after=sum(r['status'].startswith('SUCCESS') for r in rows)
        history.append(dict(round=roundid,before=before,after=after,control_sha256=meta['sha256'],coverage=dict(Counter(r['status'] for r in rows)),time=m.now()))
        m.writejson(historypath,history)
        print('BindingDB round',roundid,'completed',after,'of',len(queries),flush=True)
        no_progress=no_progress+1 if after==before else 0
        if no_progress>=CONT['BindingDB']['max_no_progress_rounds']: break


if __name__=='__main__':
    if len(sys.argv)==1 or sys.argv[1]=='probe': probe()
    elif sys.argv[1]=='bindingdb': bindingdb()
    elif sys.argv[1]=='herb-probe': herb_probe()
    elif sys.argv[1]=='bindingdb-resume': bindingdb_resume()

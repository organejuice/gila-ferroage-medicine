"""Fetch official processed files atomically; no modifications of previous snapshots."""
import concurrent.futures,gzip,hashlib,json,time
from pathlib import Path
import requests
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-gse286094-mechanism-028';RAW=ROOT/'data/raw'/RUN;OUT=ROOT/'results'/RUN
OLD=ROOT/'data/raw/2026-10-03-disease-interface-perturbation-audit-027'
CFG=json.loads((ROOT/'configs/gse286094_mechanism_028.json').read_text(encoding='utf8'))
def capture(item):
    url,name=item;p=RAW/name;meta=RAW/(name+'.retrieval.json');p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists() and meta.exists():
        r=json.loads(meta.read_text());assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'];return r
    errors=[]
    for i in range(CFG['network']['attempts']):
        part=p.with_suffix(p.suffix+'.part')
        try:
            with requests.get(url,stream=True,timeout=CFG['network']['timeout_seconds']) as resp:
                resp.raise_for_status();h=hashlib.sha256();n=0
                with part.open('wb') as f:
                    for b in resp.iter_content(1048576):f.write(b);h.update(b);n+=len(b)
            if n==0:raise ValueError('empty body')
            if name.endswith('.gz'):
                with gzip.open(part,'rb') as f:
                    while f.read(1048576):pass
            part.replace(p);r=dict(url=url,path=p.relative_to(ROOT).as_posix(),bytes=n,sha256=h.hexdigest(),status='CAPTURED',attempts=i+1);meta.write_text(json.dumps(r,indent=2)+'\n',encoding='utf8');print(name,n,flush=True);return r
        except (requests.RequestException,OSError,ValueError) as e:
            errors.append(str(e));part.unlink(missing_ok=True);time.sleep(1)
    r=dict(url=url,path=p.relative_to(ROOT).as_posix(),status='UNRESOLVED',errors=errors);meta.write_text(json.dumps(r,indent=2)+'\n');return r
def main():
    RAW.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    head=json.loads((OLD/'supplementary_HEAD_snapshot.json').read_text());items=[(r['url'],r['url'].rsplit('/',1)[-1]) for r in head if r['accession']=='GSE286094' and 'GSM' in r['url'].rsplit('/',1)[-1]]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:records=list(pool.map(capture,items))
    (OUT/'download_manifest.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf8')
if __name__=='__main__':main()

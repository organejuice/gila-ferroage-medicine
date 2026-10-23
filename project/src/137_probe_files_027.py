"""Read observed supplementary URLs only; HEAD checks do not validate contents."""
import concurrent.futures,gzip,importlib.util,json,re,time
from pathlib import Path
import requests
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('capture',ROOT/'src/136_fetch_feasibility_027.py');capture=importlib.util.module_from_spec(spec);spec.loader.exec_module(capture)
RAW=capture.RAW
def probe(item):
    acc,url=item;errors=[]
    for attempt in range(2):
        try:
            r=requests.head(url,timeout=20,allow_redirects=True)
            return dict(accession=acc,url=url,status=r.status_code,accessible=r.status_code==200,bytes=r.headers.get('Content-Length',''),content_type=r.headers.get('Content-Type',''),last_modified=r.headers.get('Last-Modified',''),body_validated=False)
        except requests.RequestException as e:errors.append(str(e))
    return dict(accession=acc,url=url,status='UNRESOLVED',accessible=False,bytes='',errors=errors,body_validated=False)
def main():
    for acc in capture.CFG['candidate_accessions']:
        capture.fetch(f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/filelist.txt',acc+'_filelist.txt')
    dest=RAW/'supplementary_HEAD_snapshot.json'
    if dest.exists():return
    items=[]
    for p in sorted(RAW.glob('*family.soft.gz')):
        acc=p.name.split('_')[0];t=gzip.open(p,'rt',encoding='utf8').read()
        urls=sorted(set(re.findall(r'^!(?:Sample|Series)_supplementary_file[^=]* = (\S+)',t,re.M)))
        items.extend((acc,u.replace('ftp://','https://',1)) for u in urls if u.startswith('ftp://ftp.ncbi.nlm.nih.gov/'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:rows=list(pool.map(probe,items))
    dest.write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print('HEAD inventory',len(rows),'accessible',sum(x['accessible'] for x in rows),flush=True)
if __name__=='__main__':main()

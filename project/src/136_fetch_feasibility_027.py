"""Capture bounded official metadata/file listings; no large omics downloads."""
import hashlib,json,time
from pathlib import Path
import requests
ROOT=Path(__file__).resolve().parents[1];CFG=json.loads((ROOT/'configs/disease_interface_perturbation_audit_027.json').read_text());RUN=CFG['run_id'];RAW=ROOT/'data/raw'/RUN;OUT=ROOT/'results'/RUN
def fetch(url,name):
    RAW.mkdir(parents=True,exist_ok=True);p=RAW/name;mp=RAW/(name+'.retrieval.json')
    if p.exists() and mp.exists():
        r=json.loads(mp.read_text());assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'];return r
    attempts=[]
    for i in range(CFG['network']['attempts']):
        time.sleep(CFG['network']['interval_seconds'])
        try:
            with requests.get(url,stream=True,timeout=CFG['network']['timeout_seconds']) as response:
                response.raise_for_status();pieces=[];total=0
                for b in response.iter_content(65536):
                    total+=len(b)
                    if total>CFG['network']['max_metadata_bytes']:raise ValueError('metadata exceeds bounded download limit')
                    pieces.append(b)
                data=b''.join(pieces)
                if not data:raise ValueError('empty response')
                if name.endswith('.gz') and not data.startswith(b'\x1f\x8b'):raise ValueError('not gzip metadata')
                p.write_bytes(data);r=dict(url=url,path=str(p.relative_to(ROOT)),status='CAPTURED',bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),http_status=response.status_code,attempts=i+1)
                mp.write_text(json.dumps(r,indent=2)+'\n');return r
        except (requests.RequestException,ValueError) as e:attempts.append(str(e))
    r=dict(url=url,path=str(p.relative_to(ROOT)),status='ACCESS_UNRESOLVED',errors=attempts);mp.write_text(json.dumps(r,indent=2)+'\n');return r
def main():
    OUT.mkdir(parents=True,exist_ok=True);ledger=[]
    for acc in CFG['candidate_accessions']:
        base=f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}'
        for suffix,name in [('soft/'+acc+'_family.soft.gz',acc+'_family.soft.gz'),('suppl/',acc+'_suppl_listing.html')]:
            r=fetch(base+'/'+suffix,name);ledger.append(r);print(acc,name,r['status'],r.get('bytes',''),flush=True)
    (OUT/'access_manifest.json').write_text(json.dumps(ledger,indent=2)+'\n')
if __name__=='__main__':main()

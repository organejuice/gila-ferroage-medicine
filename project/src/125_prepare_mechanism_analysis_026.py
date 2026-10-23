"""Acquire processed public inputs for 026 without modifying 025 snapshots."""
import csv,gzip,hashlib,json,time
from pathlib import Path
import requests
ROOT=Path(__file__).resolve().parents[1]
RUN='2026-10-03-public-mechanism-analysis-026'
OLD=ROOT/'workspaces/2026-10-02-public-mechanism-bridge-search-025/raw/sources'
RAW=ROOT/'data/raw'/RUN
OUT=ROOT/'results'/RUN

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def writejson(p,d):
    p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def get(url,name):
    p=RAW/name;mp=RAW/(name+'.metadata.json')
    if p.exists() and mp.exists():
        meta=json.loads(mp.read_text());assert sha(p)==meta['sha256'];return meta
    attempts=[]
    for i in range(3):
        time.sleep(1)
        try:
            r=requests.get(url,timeout=(15,90));attempts.append(dict(status=r.status_code,bytes=len(r.content)))
            r.raise_for_status()
            if not r.content: raise ValueError('HTTP response has no data bytes')
            if r.content.lstrip().lower().startswith((b'<!doctype html',b'<html')) and not name.endswith('.html'): raise ValueError('HTML in data response')
            p.write_bytes(r.content)
            meta=dict(url=url,path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size,attempts=attempts)
            writejson(mp,meta);return meta
        except (requests.RequestException,ValueError) as e:
            error=str(e)
            if i<2:time.sleep(2**(i+1))
    meta=dict(url=url,path=str(p.relative_to(ROOT)),status='ACCESS_FAILED',error=error,attempts=attempts)
    writejson(mp,meta);return meta

def main():
    RAW.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    accs=['GSE228927','GSE175485','GSE345471','GSE275834','GSE286095','GSE18740','GSE275256','GSE345470','GSE235076','GSE246888','GSE240609']
    manifest=[]
    for acc in accs:
        meta=get(f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/soft/{acc}_family.soft.gz',acc+'_family.soft.gz')
        manifest.append(meta);print('metadata',acc,meta.get('status','CAPTURED'),flush=True)
    # Deposited processed tables. Matrix CEL/raw sequencing not downloaded.
    files={
      'GSE228927':['GSE228927_FPKM_anno.xls.gz'],
      'GSE175485':['GSE175485_Sup_file_2.xlsx','GSE175485_analysis.xlsx','GSE175485_sup_file_3.xlsx'],
      'GSE345471':['GSE345471_1_genes_fpkm_expression.xlsx'],
      'GSE275834':['GSE275834_BLEO_vs_CTRL.deseq2.annotation_HMC3_.xls.gz'],
      'GSE286095':['GSE286095_all.gene_counts.tsv.gz'],
      'GSE275256':['GSE275256_summary_star_genes_readcount.stranded.annotated.tsv.gz'],
      'GSE345470':['GSE345470_1_genes_fpkm_expression.xlsx'],
    }
    for acc,names in files.items():
        for name in names:
            meta=get(f'https://ftp.ncbi.nlm.nih.gov/geo/series/{acc[:-3]}nnn/{acc}/suppl/{name}',name)
            manifest.append(meta);print('data',name,meta.get('status','CAPTURED'),meta.get('bytes'),flush=True)
    # Series matrices retain the deposited normalized microarray signal and metadata.
    manifest.append(get('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18740/matrix/GSE18740_series_matrix.txt.gz','GSE18740_series_matrix.txt.gz'))
    for p in OLD.iterdir():
        if p.is_file(): manifest.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size,reused_025=True))
    writejson(OUT/'source_manifest.json',manifest)
    print('manifest',len(manifest),flush=True)

if __name__=='__main__':main()

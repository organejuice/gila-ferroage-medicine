"""Sparse input and core-QC diagnostics. No cell filtering, annotation or DE."""
import gzip,hashlib,importlib.util,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.io import mmread
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-gse286094-mechanism-028';RAW=ROOT/'data/raw'/RUN;OUT=ROOT/'results'/RUN;INT=ROOT/'data/interim'/RUN
OLD=ROOT/'data/raw/2026-10-03-disease-interface-perturbation-audit-027'
spec=importlib.util.spec_from_file_location('reader05',ROOT/'src/05_prepare_glial_pseudobulk.py');reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
spec2=importlib.util.spec_from_file_location('soft127',ROOT/'src/127_prepare_expression_026.py');softreader=importlib.util.module_from_spec(spec2);spec2.loader.exec_module(softreader)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(d,name):
    p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True);d.to_csv(p,sep='\t',index=False,compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)
def validate_counts(x,nfeatures,nbarcodes):
    if x.shape!=(nfeatures,nbarcodes):raise ValueError('matrix/feature/barcode dimensions differ')
    if not np.isfinite(x.data).all() or (x.data<0).any() or not np.equal(x.data,np.floor(x.data)).all():raise ValueError('nonfinite, negative or noninteger counts')
def main():
    OUT.mkdir(parents=True,exist_ok=True);INT.mkdir(parents=True,exist_ok=True);summaries=[];cells=[];checks=[];sampleaudit=[];feature_hashes=[]
    sources=softreader.soft(OLD/'GSE286094_family.soft.gz')
    inputmanifest=json.loads((OUT/'download_manifest.json').read_text());filelist=pd.read_csv(OLD/'GSE286094_filelist.txt',sep='\t')
    for r in inputmanifest:
        p=ROOT/r['path'];expected=filelist.loc[filelist.Name==p.name,'Size'];checks.append(dict(check='download_hash:'+p.name,pass_check=r['status']=='CAPTURED' and sha(p)==r['sha256']));checks.append(dict(check='filelist_size:'+p.name,pass_check=len(expected)==1 and p.stat().st_size==int(expected.iloc[0])))
    fig,axes=plt.subplots(11,3,figsize=(12,28),layout='constrained')
    for row,r in enumerate(sources):
        gsm=r['GSM'];base=list(RAW.glob(gsm+'*_matrix.mtx.gz'));assert len(base)==1
        p=base[0];prefix=p.name.replace('_matrix.mtx.gz','');feature=RAW/(prefix+'_features.tsv.gz');barcode=RAW/(prefix+'_barcodes.tsv.gz');flines=reader.read_lines(feature);features=[l.split('\t') for l in flines];genes=[x[1] for x in features];ids=[x[0] for x in features];bars=reader.read_lines(barcode)
        with gzip.open(p,'rb') as handle:x=mmread(handle).tocsc()
        validate_counts(x,len(genes),len(bars));assert len(set(bars))==len(bars);assert len(set(ids))==len(ids)
        total=np.asarray(x.sum(axis=0)).ravel().astype(np.int64);detected=np.asarray((x>0).sum(axis=0)).ravel();mt=np.asarray(x[np.array([g.lower().startswith('mt-') for g in genes]),:].sum(axis=0)).ravel();pct=np.divide(mt,total,out=np.zeros(len(total),dtype=float),where=total>0)*100
        group=prefix[len(gsm)+1:].rstrip('123');condition={'WT':'WT','KO':'KO','AD_WT':'AD_WT','AD_KO':'AD_KO'}[group]
        d=pd.DataFrame(dict(GSM=gsm,condition=condition,barcode=bars,cell_id=[gsm+':'+b for b in bars],n_UMI=total,n_genes=detected,percent_mt=pct,passes_QC='NOT_EVALUATED',doublet_status='NOT_RUN',annotation='UNASSIGNED'));cells.append(d)
        summary=dict(GSM=gsm,library=prefix[len(gsm)+1:],condition=condition,features=x.shape[0],barcodes=x.shape[1],stored_entries=x.nnz,UMI_sum=int(total.sum()),zero_UMI_barcodes=int((total==0).sum()),median_UMI=float(np.median(total)),median_genes=float(np.median(detected)),median_percent_mt=float(np.median(pct)),p95_percent_mt=float(np.quantile(pct,.95)),symbol_duplicates=len(genes)-len(set(genes)),feature_sha256=sha(feature),matrix_verified=True,cell_QC_qualified=False)
        summaries.append(summary);feature_hashes.append(hashlib.sha256('\n'.join(flines).encode('utf8')).hexdigest())
        chars=r.get('characteristics_ch1','');sampleaudit.append(dict(GSM=gsm,title=r['title'],library=summary['library'],condition=condition,declared_biological_replicate='bio_rep' in r['title'],animal_id='NOT_SUPPLIED',pooling_status='UNKNOWN',sex='NOT_SUPPLIED',batch='NOT_SUPPLIED',age='5-month-old per GEO',region='cortices per extract protocol; tissue field Brain',assay='CD45+ sorted 10x whole-cell 3prime v3',genotype_characteristics=chars,independence_qualified=False,source=f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}'))
        for col,(v,label) in enumerate([(np.log10(np.maximum(total,1)),'log10 UMIs'),(np.log10(np.maximum(detected,1)),'log10 detected genes'),(pct,'mitochondrial %')]):
            ax=axes[row,col];ax.hist(v,bins=50,color='#397c94');ax.set_title(summary['library']+' / '+label,fontsize=9);ax.set_ylabel('barcodes');ax.tick_params(labelsize=8)
        for name,ok in [('shape',x.shape==(len(features),len(bars))),('barcodes_unique',len(set(bars))==len(bars)),('feature_ids_unique',len(set(ids))==len(ids)),('nonnegative_integer_counts',np.equal(x.data,np.floor(x.data)).all() and (x.data>=0).all()),('UMI_conservation',int(total.sum())==int(x.sum()))]:checks.append(dict(check=gsm+':'+name,pass_check=bool(ok)))
        if row==0:write(pd.DataFrame(features,columns=['ensembl_id','symbol','feature_type']),'feature_registry.tsv')
        print(summary['library'],'features',x.shape[0],'barcodes',x.shape[1],'median mt',round(summary['median_percent_mt'],2),flush=True)
        del x
    checks.append(dict(check='common_feature_registry',pass_check=len(set(feature_hashes))==1));checks.append(dict(check='all_four_groups',pass_check={g:sum(r['condition']==g for r in summaries) for g in ['WT','KO','AD_WT','AD_KO']}=={'WT':3,'KO':3,'AD_WT':2,'AD_KO':3}))
    fig.suptitle('GSE286094 deposited filtered barcodes: diagnostic distributions only; no new filtering',fontsize=12);fig.savefig(OUT/'input_QC_distributions.png',dpi=140);plt.close(fig)
    write(pd.DataFrame(summaries),'matrix_QC_diagnostics.tsv');write(pd.DataFrame(sampleaudit),'sample_identity_audit.tsv');write(pd.concat(cells,ignore_index=True),'all_unfiltered_barcode_QC.tsv.gz');write(pd.DataFrame(checks),'matrix_input_checks.tsv')
    status=dict(run_id=RUN,file_integrity='PASS' if all(r['pass_check'] for r in checks) else 'FAIL',matrix_checks=len(checks),files=len(inputmanifest),libraries=len(summaries),deposited_barcodes=sum(r['barcodes'] for r in summaries),qualified_animals='UNRESOLVED',independence_gate='NOT_PASSED',cell_QC_gate='NOT_PASSED',mechanism_analysis='NOT_EXECUTED_PENDING_SAMPLE_IDENTITY_AND_QC',no_inferential_P_or_FDR=True,thresholds_selected=False,cell_annotation_available=False,ambient_correction='NOT_QUALIFIED_FILTERED_ONLY_NO_EMPTY_DROPLETS',doublet_calling='NOT_RUN_NO_SCDblFinder_LOCAL_PACKAGE',claim_ceiling='C0',reopen_G6=False)
    (OUT/'analysis_status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    assert status['file_integrity']=='PASS';print(json.dumps(status),flush=True)
if __name__=='__main__':main()

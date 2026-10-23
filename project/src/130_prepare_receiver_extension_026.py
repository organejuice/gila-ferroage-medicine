"""Prespecified deposited neuronal receiver extension, without pooling GEO series."""
import hashlib,importlib.util,json,re
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-public-mechanism-analysis-026'
INT=ROOT/'data/interim'/RUN;OUT=ROOT/'results'/RUN;RAW=ROOT/'data/raw'/RUN
spec=importlib.util.spec_from_file_location('fetch',ROOT/'src/125_prepare_mechanism_analysis_026.py');fetch=importlib.util.module_from_spec(spec);spec.loader.exec_module(fetch)
def main():
    meta=pd.read_csv(OUT/'official_sample_metadata.tsv',sep='\t').fillna('')
    ortho=pd.read_csv(INT/'mouse_human_one_to_one.tsv',sep='\t');mp=dict(zip(ortho.mouse_symbol,ortho.human_symbol))
    datasets=[];sources=[];samples=[]
    for acc in ['GSE235076','GSE246888']:
        sub=meta[meta.series==acc];matrices=[]
        for row in sub.to_dict('records'):
            url=row['supplementary_file_1'].replace('ftp://','https://');name=url.rsplit('/',1)[-1]
            sources.append(fetch.get(url,name));print(name,flush=True)
            d=pd.read_csv(RAW/name,sep='\t',index_col=0,header=None)
            if d.shape[1]!=1:raise ValueError('unexpected receiver columns')
            d.columns=[row['title']];matrices.append(d)
            samples.append(dict(dataset=acc,sample=row['title'],group=row['title'].rsplit('_',1)[0],GSM=row['GSM'],species='mouse',cell='primary_neuron',independent_unit='deposited_experimental_replicate',limitation='rlog only; single-series comparison, no cross-series batch pooling'))
        original=matrices[0].index.astype(str)
        if not all(list(m.index.astype(str))==list(original) for m in matrices):raise ValueError('receiver source row order differs')
        matrices=[m.reset_index(drop=True) for m in matrices]
        d=pd.concat(matrices,axis=1);assert not d.isna().any().any()
        # Deposited IDs are retained; use only exact mouse symbols from frozen orthologs.
        feature=[g+'__row'+str(i+1) for i,g in enumerate(original)]
        ambiguous=original.duplicated(keep=False)
        symbols=[mp.get(g,'') if not amb else '' for g,amb in zip(original,ambiguous)]
        if not any(symbols):
            info=pd.read_csv(ROOT/'data/raw/NCBI_gene_orthologs/Mus_musculus.gene_info.gz',sep='\t',dtype=str)
            ens={entry.split(':',1)[1]:r.Symbol for r in info.itertuples() for entry in str(r.dbXrefs).split('|') if entry.startswith('Ensembl:')}
            symbols=[mp.get(ens.get(g.split('.')[0],''),'') if not amb else '' for g,amb in zip(original,ambiguous)]
        d.insert(0,'feature_id',feature);d.to_csv(INT/(acc+'_expression.tsv.gz'),sep='\t',index=False)
        pd.DataFrame(dict(feature_id=feature,original_symbol=original,human_symbol=symbols,ambiguous_source_symbol=ambiguous)).to_csv(INT/(acc+'_gene_map.tsv'),sep='\t',index=False)
        ss=pd.DataFrame([s for s in samples if s['dataset']==acc]);ss.to_csv(INT/(acc+'_samples.tsv'),sep='\t',index=False)
        datasets.append(dict(id=acc,kind='rlog_normalized',species='mouse',cell='primary_neuron',formal=acc=='GSE235076',contrasts=['M_N-A_N'] if acc=='GSE235076' else [],limitation='deposited rlog abundance; neuronal response to 6ug/ml apoE particles, not proof of full communication chain; GSE246888 control-only cannot be pooled as matched control'))
    (INT/'receiver_dataset_manifest.json').write_text(json.dumps(datasets,indent=2)+'\n')
    pd.DataFrame(samples).to_csv(OUT/'receiver_analysis_samples.tsv',sep='\t',index=False)
    (OUT/'receiver_source_manifest.json').write_text(json.dumps(sources,indent=2)+'\n')
    (OUT/'receiver_pre_effect_freeze.json').write_text(json.dumps(dict(contrast='M_N-A_N',control_only='GSE246888 descriptive only',transformation='deposited rlog, no additional log transform',inference='limma-trend exploratory; existing frozen modules; receiver BH family separate; no between-series pooling',manifest_sha256=hashlib.sha256((INT/'receiver_dataset_manifest.json').read_bytes()).hexdigest()),indent=2)+'\n')
if __name__=='__main__':main()

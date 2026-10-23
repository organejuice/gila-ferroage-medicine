"""Trace samples to official metadata and assign a disposition to every 025 source."""
import importlib.util,json,re
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-public-mechanism-analysis-026';OUT=ROOT/'results'/RUN;INT=ROOT/'data/interim'/RUN
OLD=ROOT/'workspaces/2026-10-02-public-mechanism-bridge-search-025/raw/sources'
def main():
    samples=pd.read_csv(OUT/'analysis_samples.tsv',sep='\t').fillna('')
    meta=pd.read_csv(OUT/'official_sample_metadata.tsv',sep='\t').fillna('')
    old=pd.read_csv(ROOT/'results/2026-10-02-public-mechanism-bridge-search-025/sample_mapping.tsv',sep='\t').fillna('')
    spec=importlib.util.spec_from_file_location('prep',ROOT/'src/127_prepare_expression_026.py');prep=importlib.util.module_from_spec(spec);spec.loader.exec_module(prep)
    extra=pd.DataFrame([dict(series='GSE290755',**r) for r in prep.soft(ROOT/'data/raw'/RUN/'GSE290755_family.soft.gz')])
    meta=pd.concat([meta[meta.series!='GSE290755'],extra],ignore_index=True).fillna('')
    sdrf=pd.read_csv(OLD/'E-MTAB-15005.sdrf.txt',sep='\t').fillna('');rows=[]
    for s in samples.to_dict('records'):
        acc=s['dataset'];name=s['sample'];sub=meta[meta.series==acc]
        cand=[];rule=''
        if s['GSM']:cand=[s['GSM']];rule='existing exact 025 mapping'
        elif acc.startswith('E_MTAB'):
            key=re.sub(r'-(\d)$',r'\1',name.replace('-', '_',1)).replace('-', '_')
            # SDRF uses HMC3_BrdU3 while matrix uses HMC3-BrdU-3.
            key=name.split('-')[0]+'_'+name.split('-')[1]+name.split('-')[2]
            hit=sdrf[sdrf['Source Name']==key].drop_duplicates('Source Name')
            cand=hit['Comment[ENA_SAMPLE]'].tolist();rule='exact SDRF Source Name; ENA sample'
        elif acc=='GSE275834':cand=sub.loc[sub.title=='HMC3 '+name.replace('_tpm',''),'GSM'].tolist();rule='exact condition replicate title'
        elif acc in ['GSE345471','GSE345470']:cand=sub.loc[sub.title==name.replace('count.',''),'GSM'].tolist();rule='exact count/FPKM title'
        elif acc=='GSE175485':
            base=re.match(r'(ECM|MCM)_(PBS|LPS)_\d',name).group(0);cand=sub.loc[sub.title==base,'GSM'].tolist();rule='exact condition replicate title before sequencing suffix'
        elif acc=='GSE228927':
            group={'NOT':'Resting','LPS':'LPS','Hemin':'Hemin'}[name[:-1]];cand=sub.loc[sub.title==group+'-MCM, replicate '+name[-1],'GSM'].tolist();rule='official resting/LPS/hemin replicate title'
        elif acc=='GSE275256':
            base=name.replace('SLO_1-2','SLO 1/2').replace('SLO_1','SLO 1').replace('VEH','VEH control');base=base.rsplit('_',1);cand=sub.loc[sub.title=='HMC3 cells, '+base[0]+', replicate '+base[1],'GSM'].tolist();rule='official SLO dose/replicate title'
        elif acc=='GSE286095':
            lib=name.replace('sample.','').upper();cand=sub.loc[sub.description=='Library name: '+lib,'GSM'].tolist();rule='exact official Library name WT/KO'
        elif acc=='GSE18740':cand=sub.loc[sub.GSM==name,'GSM'].tolist();rule='matrix GSM column'
        elif acc=='GSE325931':cand=sub.loc[sub.description=='Library name: '+name,'GSM'].tolist();rule='exact official Library name'
        elif acc=='GSE290755':
            title='iTF_Microglia_RNA_'+('Arid5aKD' if 'arid5a_kd' in name else 'Nontargeting')+'_Rep'+name[-1]
            cand=sub.loc[sub.title==title,'GSM'].tolist();rule='full SOFT RNA condition replicate; RPF excluded'
        identifier=cand[0] if len(cand)==1 else ''
        rows.append(dict(**s,official_sample_id=identifier,mapping_rule=rule,mapping_status='EXACT_ONE' if len(cand)==1 else 'UNRESOLVED',candidate_count=len(cand)))
    if (OUT/'receiver_analysis_samples.tsv').exists():
        for s in pd.read_csv(OUT/'receiver_analysis_samples.tsv',sep='\t').to_dict('records'):rows.append(dict(**s,official_sample_id=s['GSM'],mapping_rule='deposited per-GSM file',mapping_status='EXACT_ONE',candidate_count=1))
    trace=pd.DataFrame(rows);trace.to_csv(OUT/'sample_traceability_audit.tsv',sep='\t',index=False)
    registry=pd.read_csv(ROOT/'results/2026-10-02-public-mechanism-bridge-search-025/candidate_registry.tsv',sep='\t').fillna('')
    analysed=set(samples.dataset)-{'GSE325931'};analysed.add('E-MTAB-15005');analysed.update(['GSE235076'])
    dispositions=[]
    for r in registry.to_dict('records'):
        acc=r['accession']
        if acc in analysed:state='PROCESSED_EXPRESSION_REANALYSED';reason='deposited sample-level matrix, replication and frozen contrasts; RNA only'
        elif acc=='GSE325931':state='DESCRIPTIVE_ONLY_BATCH_CONFOUNDED';reason='treatment and batch not identifiable; no inferential P/CI'
        elif acc=='GSE246888':state='CONTROL_ONLY_QC_COMPLETED';reason='3 control libraries; no within-series intervention; cross-series pooling prohibited'
        elif acc=='GSE240609':state='NOT_INFERENTIAL_SINGLE_LIBRARY_PER_CONDITION';reason='four conditions n=1, cannot estimate biological replication'
        elif acc=='YZDHD_2026':state='SUPPLEMENTS_INGESTED_024';reason='user-supplied WB figure, ingredients, PPI genes; not sample-level perturbation matrix'
        elif r['priority']=='EXISTING':state='EXISTING_REFERENCE_RETAINED';reason='existing disease/lactylation reference; no new formula intervention qualification'
        elif acc in ['GSE286094','GSE333352','GSE315101']:state='UNFINISHED_CELL_LEVEL_REANALYSIS';reason='single-cell/spatial requires separate cell/sample annotations, count QC and donor-aware workflow; no bulk substitution'
        elif acc=='PXD079947':state='PROCESSED_PROTEOME_REANALYSED';reason='54 Spectronaut quantity channels aggregated to 18 block-condition samples; no raw MS search or missing-value imputation; original peptide/protein FDR provenance recorded'
        elif acc=='FIGSHARE_30252505':state='FUNCTIONAL_FIGURE_SOURCE_ONLY';reason='API file inventory contains 11 TIFF figures and no sample-level numeric source table; original Srr study reviewed; do not infer quantitative reference from plotted images'
        elif acc=='TRAP1_VSMC':state='REQUEST_ONLY_CONTEXT_MISMATCH';reason='vascular smooth muscle, not present microglia/astrocyte endpoint; no public full matrix'
        else:state='CONTEXT_REFERENCE_NOT_REANALYSED';reason='lower-priority non-AD or nonmatching-system reference; not migrated to current therapeutic qualification'
        if acc=='GSE240609':r['title']='Microglial APOE3 Christchurch protects neurons from Tau pathology in a human iPSC-based model of Alzheimer\'s disease'
        dispositions.append(dict(accession=acc,priority=r['priority'],title=r['title'],doi=r['doi'],status=state,reason=reason,url=r['url']))
    pd.DataFrame(dispositions).to_csv(OUT/'all_025_source_dispositions.tsv',sep='\t',index=False)
    print('samples',len(trace),'unresolved',(trace.mapping_status!='EXACT_ONE').sum(),'sources',len(dispositions),flush=True)
if __name__=='__main__':main()

"""Conservative three-block proteome reanalysis; never count technical triplicates as N=9."""
import hashlib,json,re
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-public-mechanism-analysis-026';INT=ROOT/'data/interim'/RUN;OUT=ROOT/'results'/RUN;RAW=ROOT/'data/raw'/RUN
def main():
    # Record the approach before fitting any protein effects.
    decision=dict(source='PXD079947 raw_df.tsv; DOI 10.1155/omcl/9948225 methods 4.16.2',
      source_search='Spectronaut 19.7; deposited automatic cross-run quantity, peptide/protein 1% as reported, not independently re-searched',
      unit='B1/B3/B4 experimental blocks; aggregate A/B/C technical triplicates to one log2 mean; 3 per condition, never 9',
      missing='complete positive measurements in all 54 channels; no stochastic imputation; sensitivity restricted to complete-case proteins',
      normalization='deposited PG.Quantity, log2 once; no second cross-run normalization',model='limma-trend group + experimental batch; BH .05 exploratory, not identical to publisher BY/imputation analysis',
      contrasts=['IronDMSO-DMSO','IronRSL3-DMSO','IronRSL3Fer-IronRSL3','RSL3-DMSO','DMSO-NoTreatment'],
      modules='existing frozen G1-G5 mapped unambiguous gene symbols; mixed protein group genes excluded from modules',
      limitation='processed quantity analysis only; no raw MS search, missingness limits generalization, plate/sample provenance still restricted',normalization_eligible=False,reopen_017=False)
    (OUT/'protein_pre_effect_protocol.json').write_text(json.dumps(decision,indent=2)+'\n')
    d=pd.read_csv(RAW/'PXD079947_raw_df.tsv',sep='\t');cols=[c for c in d if c.endswith('.PG.Quantity')];assert len(cols)==54
    x=d[cols].apply(pd.to_numeric);keep=(x.notna() & (x>0)).all(axis=1);log=np.log2(x.loc[keep]);m=[];result={}
    for c in cols:
        batch,group,tech=re.search(r'(B\d)-02-HMC3-(\w+)-1500.*Elite-([ABC])_',c).groups()
        m.append(dict(dataset='PXD079947',input_column=c,batch=batch,group=group,technical_replicate=tech))
    md=pd.DataFrame(m);md.to_csv(OUT/'protein_54_channel_mapping.tsv',sep='\t',index=False)
    samples=[]
    for (batch,group),sub in md.groupby(['batch','group'],sort=True):
        assert set(sub.technical_replicate)==set('ABC');name=batch+'_'+group;result[name]=log[sub.input_column].mean(axis=1).to_numpy()
        samples.append(dict(dataset='PXD079947',sample=name,batch=batch,group=group,independent_unit='experimental batch, triplicates aggregated'))
    ds=d.loc[keep].reset_index(drop=True);ids=ds['PG.ProteinGroups'].astype(str);assert ids.is_unique
    approved=set(pd.read_csv(ROOT/'data/raw/HGNC/hgnc_complete_set.txt',sep='\t',dtype=str)['symbol']) if (ROOT/'data/raw/HGNC/hgnc_complete_set.txt').exists() else set()
    # Approved symbols already frozen in the 026 human mappings.
    if not approved:
        for p in INT.glob('*_gene_map.tsv'):approved.update(pd.read_csv(p,sep='\t').human_symbol.dropna())
    symbols=[str(g) if str(g) in approved else '' for g in ds['PG.Genes']]
    expr=pd.DataFrame(result);expr.insert(0,'feature_id',ids);expr.to_csv(INT/'PXD079947_expression.tsv.gz',sep='\t',index=False)
    pd.DataFrame(dict(feature_id=ids,original_gene=ds['PG.Genes'],human_symbol=symbols)).to_csv(INT/'PXD079947_gene_map.tsv',sep='\t',index=False)
    pd.DataFrame(samples).to_csv(INT/'PXD079947_samples.tsv',sep='\t',index=False)
    dset=dict(id='PXD079947',kind='protein_log2',species='human',cell='HMC3',formal=True,block_covariate=True,contrasts=decision['contrasts'],limitation=decision['limitation'])
    (INT/'protein_dataset_manifest.json').write_text(json.dumps([dset],indent=2)+'\n')
    (OUT/'protein_input_qc.json').write_text(json.dumps(dict(input_protein_groups=len(d),complete_case_protein_groups=int(keep.sum()),excluded_missing_or_nonpositive=int((~keep).sum()),channels=54,aggregated_samples=18,independent_blocks=3,unambiguous_mapped=int(sum(bool(s) for s in symbols)),protocol_sha256=hashlib.sha256((OUT/'protein_pre_effect_protocol.json').read_bytes()).hexdigest()),indent=2)+'\n')
    # Sample traceability sensitivity, fixed by missing SDRF entries rather than result.
    s=pd.read_csv(INT/'E_MTAB_15005_SVGA_samples.tsv',sep='\t');s=s[~s['sample'].str.endswith('-1')];s['dataset']='SVGA_verified4';s.to_csv(INT/'SVGA_verified4_samples.tsv',sep='\t',index=False)
    x=pd.read_csv(INT/'E_MTAB_15005_SVGA_expression.tsv.gz',sep='\t');x[['feature_id']+s['sample'].tolist()].to_csv(INT/'SVGA_verified4_expression.tsv.gz',sep='\t',index=False)
    (INT/'SVGA_verified4_gene_map.tsv').write_bytes((INT/'E_MTAB_15005_SVGA_gene_map.tsv').read_bytes())
    a=next(a for a in json.loads((INT/'dataset_manifest.json').read_text()) if a['id']=='E_MTAB_15005_SVGA');a['id']='SVGA_verified4';a['limitation']='sensitivity only: two matrix samples lack SDRF match and are excluded; 2 per group';(INT/'sample_sensitivity_manifest.json').write_text(json.dumps([a],indent=2)+'\n')
    print('protein complete groups',keep.sum(),'of',len(d),'channels 54 aggregated samples 18',flush=True)
if __name__=='__main__':main()

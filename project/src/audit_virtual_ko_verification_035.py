"""Independent ORA and representative positive-GSEA ES/leading-edge replay."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
import yaml
ROOT=Path(__file__).resolve().parents[1]
CFG=yaml.safe_load((ROOT/'configs/virtual_ko_verification_execution_035.yaml').read_text(encoding='utf-8'))
OUT=ROOT/CFG['output_dir']

def main():
    jobs=pd.read_csv(OUT/'job_registry.tsv',sep='\t');libs={}
    for lib in CFG['libraries']:
        libs[lib]={x[0]:set(x[2:]) for x in (line.split('\t') for line in (ROOT/CFG['library_dir']/f'{lib}.gmt').read_text(encoding='utf-8').splitlines())}
    results=[]
    for r in jobs.itertuples(index=False):
        d=OUT/'jobs'/r.job_id
        if not (d/'completion.json').exists():continue
        dr=pd.read_csv(ROOT/r.path,sep='\t');dr=dr[dr.gene!=r.target].sort_values(['distance','gene'],ascending=[False,True]);bg=set(dr.gene);sig=set(dr.loc[dr['p.adj']<.05,'gene'])
        o=pd.read_csv(d/'ORA_all.tsv.gz',sep='\t');g=pd.read_csv(d/'GSEA_all.tsv.gz',sep='\t');maxP=0;maxES=0;nReplay=0
        for lib,terms in libs.items():
            univ=bg&set().union(*terms.values());q=sig&univ
            tab=o[o.library==lib]
            expected={k for k,s in terms.items() if 5<=len(s&univ)<=500};assert set(tab.term)==expected
            for p in tab.itertuples(index=False):
                members=terms[p.term]&univ;hits=members&q
                value=hypergeom.sf(len(hits)-1,len(univ),len(members),len(q)) if q else 1
                error=abs(p.P-value);maxP=max(maxP,error);assert error<1e-12
                assert set(str(p.genes).split(';'))-{'','nan'}==hits
            for mode,power in CFG['gsea']['modes'].items():
                tab=g[(g.library==lib)&(g['mode']==mode)].sort_values('pval').head(3)
                for p in tab.itertuples(index=False):
                    indices=np.flatnonzero(dr.gene.isin(terms[p.pathway]).to_numpy());N=len(dr);Nh=len(indices)
                    weight=dr.distance.to_numpy()[indices]**power
                    weight=weight/weight.sum() if weight.sum()>0 else np.repeat(1/Nh,Nh)
                    tops=np.cumsum(weight)-(indices-np.arange(Nh))/(N-Nh);imax=int(np.argmax(tops));es=float(tops[imax]);maxES=max(maxES,abs(es-p.ES));assert abs(es-p.ES)<1e-10
                    edge=set(dr.gene.iloc[indices[:imax+1]]);assert edge==set(str(p.leadingEdge).split(';'))
                    nReplay+=1
        results.append(dict(job_id=r.job_id,ORA_tests=len(o),GSEA_ES_replays=nReplay,max_ORA_P_error=maxP,max_GSEA_ES_error=maxES))
    pd.DataFrame(results).to_csv(OUT/'independent_numerical_audit_jobs.tsv',sep='\t',index=False)
    state=dict(status='PASS_FOR_COMPLETED_JOBS' if results else 'NO_COMPLETED_ENRICHMENT_JOBS',registered_jobs=len(jobs),audited_jobs=len(results),coverage='partial until all registered and planned jobs complete',ORA_tests=sum(x['ORA_tests'] for x in results),GSEA_ES_replays=sum(x['GSEA_ES_replays'] for x in results),scientific_identity_validation=False)
    (OUT/'enrichment_independent_audit.json').write_text(json.dumps(state,ensure_ascii=False,indent=2),encoding='utf-8');print(json.dumps(state))

if __name__=='__main__':main()

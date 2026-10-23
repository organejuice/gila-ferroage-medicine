"""Final evidence review of completed SOP035; never alters source artifacts.

Reuses the recorded GMT, KO, enrichment, expression and donor tables. Audits
every enrichment output, replays representative ES, and derives review tables.
Run with the existing TCM_select Python environment; no network or R required.
"""
from pathlib import Path
import argparse, hashlib, json, sys, time, platform, subprocess
import numpy as np
import pandas as pd
import yaml
from scipy.stats import hypergeom, rankdata

ROOT = Path(__file__).resolve().parents[1]
CFG_PATH = ROOT / 'configs/virtual_ko_scientific_review_035.yaml'
CFG = yaml.safe_load(CFG_PATH.read_text(encoding='utf-8'))
SRC = ROOT / CFG['source_dir']; OUT = ROOT / CFG['output_dir']
KEY = ['analysis', 'context', 'target']

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()

def table(path, **kw):
    return pd.read_csv(path, sep='\t', **kw)

def save(df, name):
    df.to_csv(OUT/name, sep='\t', index=False, lineterminator='\n')

def js(name, value):
    (OUT/name).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')

def genes(s):
    return set(str(s).split(';'))-{'', 'nan', 'None'}

def bh(v):
    v = np.asarray(v, dtype=float); out = np.full(len(v), np.nan)
    ix = np.flatnonzero(np.isfinite(v)); order = ix[np.argsort(v[ix], kind='stable')]
    out[order] = np.minimum(1, np.minimum.accumulate((v[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1])
    return out

def fixtures():
    assert np.allclose(bh([.01,.04,.03]), [.03,.04,.04])
    assert np.isnan(bh([np.nan,.1])[0]) and bh([np.nan,.1])[1] == .1
    # Nonestimable replicates stay in unconditional denominators.
    assert 2/6 != 2/3
    assert (1+sum(np.array([.4,.5,.7]) >= .5))/4 == .75
    assert genes('FTH1;FTL;FTH1') == {'FTH1','FTL'}
    return {'BH_known_values':True,'missing_P_preserved':True,'nonestimable_denominator':True,
            'reference_ties_included':True,'gene_sets_deduplicated':True}

def init():
    OUT.mkdir(parents=True, exist_ok=True)
    js('deterministic_checks.json', {'status':'PASS','checks':fixtures()})
    cfg = yaml.safe_load((ROOT/CFG['source_config']).read_text(encoding='utf-8'))
    scope = table(ROOT/CFG['scope_path']).fillna('')
    assert len(scope)==29 and scope.target.nunique()==15 and not scope.duplicated(KEY).any()
    assert scope.context.str.endswith('_AD').sum()==14
    reg=table(SRC/'job_registry.tsv'); meta=table(SRC/'enrichment_job_metadata.tsv').fillna('')
    assert len(reg)==2707 and not reg.job_id.duplicated().any()
    assert set(reg.job_id)==set(meta.job_id)
    tasks=table(SRC/'network_task_registry.tsv').fillna(''); statuses=table(SRC/'all_resampling_KO_statuses.tsv').fillna('')
    assert len(tasks)==130 and not tasks.task_id.duplicated().any()
    assert len(statuses)==29*26*2 and not statuses.duplicated(['task_id','target','seed']).any()
    assert set(statuses.status)<= {'COMPLETED','LOW_EXPRESSION_NOT_ESTIMABLE','NO_OUTGOING_EDGES','INSUFFICIENT_CELLS'}
    for t in tasks.itertuples(index=False):
        d=SRC/'networks'/t.task_id; a=json.loads((d/'completion.json').read_text())
        assert a['status']=='PASS' and a['task_id']==t.task_id
        q=table(d/'KO_status.tsv').fillna(''); z=statuses[statuses.task_id==t.task_id]
        assert set(zip(q.target,q.seed,q.status))==set(zip(z.target,z.seed,z.status))
        inst=table(d/'cell_instances.tsv'); assert not inst.instance_cell_id.duplicated().any()
        counts=inst.groupby(['donor_instance_id','region']).size(); assert counts.max()<=500
        assert inst.groupby('donor_instance_id').source_donor_id.nunique().max()==1
        original_donors=json.loads(t.donors_json)
        assert inst.donor_instance_id.nunique()==len(original_donors)
        assert sorted(inst.groupby('donor_instance_id').source_donor_id.first())==sorted(original_donors)
        assert len(original_donors)==(5 if t.stage=='LODO' else 6)
    mf=[]
    frozen=json.loads((SRC/'frozen_input_hashes.json').read_text())
    for p,h in frozen.items():
        actual=sha(ROOT/p); role='scientific_input'
        unchanged=actual==h
        if p=='src/verify_virtual_ko_035.R':
            a=json.loads((SRC/'adapter_real_replay.json').read_text()); assert a['adapter_sha256']==actual and a['status']=='PASS'
            role='adapter_current_validated'
        if p.startswith('src/') or p==CFG['source_config']: role='code_or_execution_config'
        mf.append(dict(path=p,original_sha256=h,current_sha256=actual,unchanged=unchanged,role=role))
        if role=='scientific_input': assert unchanged, p
    save(pd.DataFrame(mf),'frozen_input_review.tsv')
    environment={'python':sys.version,'executable':sys.executable,'numpy':np.__version__,'pandas':pd.__version__,
                 'platform':platform.platform(),'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                 'source_R_versions':table(SRC/'environment_versions.tsv').to_dict('records')}
    js('environment.json',environment)
    (OUT/'frozen_review_params.yaml').write_text(CFG_PATH.read_text(encoding='utf-8'),encoding='utf-8')
    trace=ROOT/'workspaces'/CFG['run_id']/'raw/traces';trace.mkdir(parents=True,exist_ok=True)
    js('task_closure.json', {'status':'PASS','network_tasks':130,'target_resampling_attempts':len(statuses),
        'resampling_status_counts':statuses.status.value_counts().to_dict(),'enrichment_jobs':len(reg),
        'scope_combinations':29,'unique_targets':15,'direction_estimated':False})
    libs={}
    for lib in cfg['libraries']:
        libs[lib]={x[0]:(x[1],set(x[2:])) for x in (s.split('\t') for s in (ROOT/cfg['library_dir']/f'{lib}.gmt').read_text(encoding='utf-8').splitlines())}
    return cfg,scope,reg,meta,tasks,statuses,libs

AUDIT_LIBS=None
AUDIT_CONFIG=None
AUDIT_FOCUS=None

def audit_worker_init(cfg,libs,focus):
    global AUDIT_LIBS,AUDIT_CONFIG,AUDIT_FOCUS
    AUDIT_LIBS=libs;AUDIT_CONFIG=cfg;AUDIT_FOCUS=focus

def audit_one(job):
    from types import SimpleNamespace
    j=SimpleNamespace(**job);cfg=AUDIT_CONFIG;libs=AUDIT_LIBS
    d=SRC/'jobs'/j.job_id;a=json.loads((d/'completion.json').read_text())
    assert a['status']=='PASS' and a['source_sha256']==j.sha256==sha(ROOT/j.path)
    dr=table(ROOT/j.path);dr=dr[dr.gene!=j.target].sort_values(['distance','gene'],ascending=[False,True])
    assert not dr.gene.duplicated().any() and np.isfinite(dr.distance).all() and (dr.distance>=0).all()
    assert dr['p.adj'].between(0,1).all()
    bg=set(dr.gene);sig=set(dr.loc[dr['p.adj']<CFG['fdr'],'gene'])
    o=table(d/'ORA_all.tsv.gz');g=table(d/'GSEA_all.tsv.gz');qc=table(d/'QC.tsv')
    assert len(o)==a['ora_tests'] and len(g)==a['gsea_tests']
    assert not o.duplicated(['library','term']).any() and not g.duplicated(['library','mode','pathway']).any()
    assert o.P.between(0,1).all() and g.pval.between(0,1).all()
    calc=np.where(o.query_size.to_numpy()>0,hypergeom.sf(o.overlap.to_numpy()-1,o.background_size,o.term_size,o.query_size),1)
    op_error=float(np.max(np.abs(calc-o.P)));assert op_error<CFG['ora_tolerance']
    obh=float(np.max(np.abs(bh(o.P)-o.FDR_job_all_libraries)));assert obh<CFG['numeric_tolerance']
    gbh=0.; es_error=0.; n_es=0; focus_keys=AUDIT_FOCUS.get((j.analysis,j.context,j.target),set())
    for mode,q in g.groupby('mode'):
        gbh=max(gbh,float(np.max(np.abs(bh(q.pval)-q.FDR_job_all_libraries))))
    assert gbh<CFG['numeric_tolerance']
    assert set(g.leadingEdge.fillna('').str.split(';').explode())-{''}<=bg
    for lib,terms in libs.items():
        univ=bg&set().union(*(s for _,s in terms.values()));query=sig&univ
        oq=o[o.library==lib]
        assert set(oq.background_size)=={len(univ)} and set(oq.query_size)=={len(query)}
        expected_o={k:len(s&univ) for k,(_,s) in terms.items()};expected_o={k:n for k,n in expected_o.items() if cfg['ora']['min_size']<=n<=cfg['ora']['max_size']}
        assert set(oq.term)==set(expected_o)
        assert all(p.term_size==expected_o[p.term] for p in oq.itertuples(index=False))
        for p in oq.itertuples(index=False):
            hits=terms[p.term][1]&query
            assert p.overlap==len(hits) and genes(p.genes)==hits
        expected_g={k:len(s&bg) for k,(_,s) in terms.items()};expected_g={k:n for k,n in expected_g.items() if cfg['gsea']['min_size']<=n<=cfg['gsea']['max_size']}
        for mode,power in cfg['gsea']['modes'].items():
            q=g[(g.library==lib)&(g['mode']==mode)]
            assert set(q.pathway)==set(expected_g)
            assert all(p.size==expected_g[p.pathway] for p in q.itertuples(index=False))
            assert np.max(np.abs(q.pval-(q.nMoreExtreme+1)/(cfg['gsea']['nperm']+1)))<CFG['numeric_tolerance']
            assert np.max(np.abs(bh(q.pval)-q.FDR_library))<CFG['numeric_tolerance']
            for p in q.sort_values(['pval','pathway']).head(CFG['independent_ES_replays_per_library_mode']).itertuples(index=False):
                idx=np.flatnonzero(dr.gene.isin(terms[p.pathway][1]));N=len(dr);Nh=len(idx)
                weights=dr.distance.to_numpy()[idx]**power
                weights=weights/weights.sum() if weights.sum()>0 else np.full(Nh,1/Nh)
                tops=np.cumsum(weights)-(idx-np.arange(Nh))/(N-Nh);imax=int(np.argmax(tops));es=float(tops[imax])
                es_error=max(es_error,abs(es-p.ES));assert abs(es-p.ES)<CFG['numeric_tolerance']
                assert set(dr.gene.iloc[idx[:imax+1]])==genes(p.leadingEdge);n_es+=1
        assert np.max(np.abs(bh(oq.P)-oq.FDR_library))<CFG['numeric_tolerance']
    selected=[]
    for method,q in [('ORA',o),('GSEA',g)]:
        for pos,p in enumerate(q.itertuples(index=False)):
            mode='significant_gene_set' if method=='ORA' else p.mode
            term=p.term if method=='ORA' else p.pathway
            if (method,mode,p.library,term) in focus_keys:
                selected.append(dict(job_id=j.job_id,analysis=j.analysis,context=j.context,target=j.target,seed=j.seed,
                    method=method,mode=mode,library=p.library,term=term,p=float(p.P if method=='ORA' else p.pval),
                    main_FDR=p.FDR_job_all_libraries,overlap=p.overlap if method=='ORA' else np.nan,
                    NES=p.NES if method=='GSEA' else np.nan,contributors=p.genes if method=='ORA' else p.leadingEdge,
                    local_index=pos))
    sums=dict(job_id=j.job_id,analysis=j.analysis,context=j.context,target=j.target,seed=j.seed,n_tested=len(dr),n_significant=len(sig),
             ORA_rows=len(o),GSEA_rows=len(g),GSEA_MC_resolution_hits=int((g.nMoreExtreme==0).sum()))
    row=dict(job_id=j.job_id,ORA_tests=len(o),GSEA_tests=len(g),max_ORA_P_error=op_error,max_ORA_BH_error=obh,
             max_GSEA_BH_error=gbh,ES_replays=n_es,max_ES_error=es_error,status='PASS')
    return row,sums,selected,o.P.to_numpy(),g.pval.to_numpy()

def audit(cfg,reg,libs):
    from concurrent.futures import ProcessPoolExecutor
    focus=table(ROOT/'results/2026-10-05-virtual-ko-pathways-034/two_seed_pathways.tsv').fillna('')
    focus_sets={k:set(zip(q.method,q['mode'],q.library,q.term)) for k,q in focus.groupby(KEY)}
    rows=[];selected=[];sums=[];all_p={'ORA':[],'GSEA':[]};offset={'ORA':0,'GSEA':0}
    with ProcessPoolExecutor(max_workers=3,initializer=audit_worker_init,initargs=(cfg,libs,focus_sets)) as pool:
        for i,(row,summ,sel,op,gp) in enumerate(pool.map(audit_one,reg.to_dict('records'),chunksize=1),1):
            for r in sel:
                r['family_index']=offset[r['method']]+r.pop('local_index');selected.append(r)
            all_p['ORA'].append(op);all_p['GSEA'].append(gp);offset['ORA']+=len(op);offset['GSEA']+=len(gp)
            rows.append(row);sums.append(summ)
            if i%50==0:print(f'AUDIT {i}/{len(reg)}',flush=True)
    save(pd.DataFrame(rows),'independent_numerical_audit_jobs.tsv');save(pd.DataFrame(sums),'audited_job_summary.tsv')
    selected=pd.DataFrame(selected)
    for method in ['ORA','GSEA']:
        v=np.concatenate(all_p.pop(method));q=bh(v);ix=selected.method==method
        selected.loc[ix,'FDR_exploratory_all_run']=q[selected.loc[ix,'family_index'].astype(int)]
        print(f'GLOBAL_BH {method} {len(v)}',flush=True);del q,v
    save(selected.drop(columns='family_index'),'focus_term_evidence_all_jobs.tsv.gz')
    js('numerical_audit.json',dict(status='PASS',jobs=len(rows),ORA_tests=sum(x['ORA_tests'] for x in rows),
        GSEA_tests=sum(x['GSEA_tests'] for x in rows),ES_replays=sum(x['ES_replays'] for x in rows),
        max_ORA_P_error=max(x['max_ORA_P_error'] for x in rows),max_ES_error=max(x['max_ES_error'] for x in rows),
        full_ORA_background_and_hits_reconstructed=True,full_BH_reconstructed=True,
        GSEA_permutation_rng_replayed=False,identity_or_causal_pass=False,local_audit_processes=3))
    return selected

def identity():
    c=table(SRC/'cell_identity_audit.tsv'); panels=['Astrocyte','Microglia','Neuron','Oligodendrocyte','OPC','Endothelial','Pericyte']
    rows=[]
    for keys,q in c.groupby(['context','donor_id','region','seurat_clusters'],dropna=False):
        original=keys[0].split('_')[0];alternative='Microglia' if original=='Astrocyte' else 'Astrocyte'
        rows.append(dict(zip(['context','donor_id','region','cluster'],keys))|dict(N=len(q),
          original_panel_mean=q[original].mean(),alternative_panel_mean=q[alternative].mean(),
          n_alternative_higher=int((q[alternative]>q[original]).sum()),
          n_both_panels_two_markers=int(((q.Astrocyte_n_detected_markers>=2)&(q.Microglia_n_detected_markers>=2)).sum()),
          dominant_mean_panel=q[panels].mean().idxmax(),review_label=original,
          review_status='DESCRIPTIVE_ORIGINAL_LABEL_RETAINED_REFERENCE_AND_DOUBLET_UNRESOLVED'))
    save(pd.DataFrame(rows),'label_review_by_cluster_donor_region.tsv')
    save(table(SRC/'library_QC_availability.tsv'),'library_QC_availability.tsv')
    return table(SRC/'cell_identity_context_summary.tsv')

def matching(scope,cfg,libs):
    rows=[];quality=[];fm=table(ROOT/'data/interim/2026-10-03-public-mechanism-analysis-026/frozen_modules.tsv')
    modules=fm.groupby('module').gene.apply(set).to_dict()
    # Fixed pre-existing functional modules only; no new response-selected list.
    # Translation terms are separately displayed as existing focus evidence.
    for (a,context),ss in scope.groupby(['analysis','context']):
        d=SRC/'matched'/f'{a}__{context}';m=table(d/'matched_target_controls.tsv');status=table(d/'KO_status.tsv')
        assert json.loads((d/'completion.json').read_text())['status']=='PASS'
        assert (status.status=='COMPLETED').all()
        base='results/2026-10-04-stratified-virtual-ko-033' if a=='primary' else 'results/2026-10-05-virtual-ko-target-retention-sensitivity-033'
        cache={}
        for r in status.itertuples(index=False):
            dr=table(ROOT/r.path);dr=dr[dr.gene!=r.target];cache[(r.target,r.seed)]=dr.set_index('gene')
        for s in ss.itertuples(index=False):
            mm=m[m.target==s.target];assert mm.control.nunique()==20 and len(mm)==40
            limited=bool(mm.matching_limited.any())
            quality.append(dict(analysis=a,context=context,target=s.target,reference_controls=20,matching_limited=limited,
                min_reference_detection=mm.control_detection_fraction.min(),target_detection=mm.target_detection_fraction.iloc[0],
                median_matching_distance=mm.distance.median(),max_matching_distance=mm.distance.max(),
                target_expression_out_of_reference_range=bool((mm.target_detection_fraction<mm.control_detection_fraction.min()).any()),
                matching_quality='LIMITED_LOW_EXPRESSION_REFERENCE_SUPPORT' if limited else 'NEAREST_REFERENCE_DESCRIPTIVE_ONLY',
                empirical_resolution=1/21,formal_specificity_inference=False))
            for seed in CFG['network_seeds']:
                dr=table(ROOT/base/context/f'seed_{seed}/{s.target}_diffRegulation.tsv');dr=dr[dr.gene!=s.target].set_index('gene')
                for name,gs in modules.items():
                    if name.startswith('G5_') != context.startswith('Astrocyte'):continue
                    common=gs-{s.target}
                    common&=set(dr.index)
                    for control in mm.control.unique():common&=set(cache[(control,seed)].index)
                    if len(common)<2:
                        rows.append(dict(analysis=a,context=context,target=s.target,seed=seed,function=name,n_members=len(common),status='INSUFFICIENT_COMMON_COVERAGE'))
                        continue
                    def metric(df):
                        ranks=pd.Series(rankdata(df.distance,method='average')/(len(df)+1),index=df.index)
                        return float(ranks.loc[sorted(common)].mean()),float((df.loc[sorted(common),'p.adj']<CFG['fdr']).mean())
                    tv,ts=metric(dr);cv=[metric(cache[(control,seed)])[0] for control in mm.control.unique()]
                    rows.append(dict(analysis=a,context=context,target=s.target,seed=seed,function=name,n_members=len(common),
                        genes=';'.join(sorted(common)),target_mean_percentile_rank=tv,target_significant_fraction=ts,
                        reference_mean=np.mean(cv),reference_min=min(cv),reference_max=max(cv),
                        reference_tail_fraction=(1+sum(v>=tv for v in cv))/(len(cv)+1),n_controls=len(cv),
                        reference_values=';'.join(f'{v:.12g}' for v in cv),matching_limited=limited,
                        status='DESCRIPTIVE_CONDITIONAL_REFERENCE_NO_CAUSAL_P'))
    save(pd.DataFrame(quality),'matched_control_QC.tsv');save(pd.DataFrame(rows),'matched_response_comparison.tsv')
    return pd.DataFrame(quality)

def stability(scope,statuses,selected,meta):
    focus=table(ROOT/'results/2026-10-05-virtual-ko-pathways-034/two_seed_pathways.tsv').fillna('')
    used=['job_id','analysis','context','target','seed','method','mode','library','term','main_FDR','overlap','FDR_exploratory_all_run']
    relevant=meta[meta.stage.isin(['LODO','BOOTSTRAP'])][['job_id','stage','replicate']]
    sel=selected[used].merge(relevant,on='job_id',how='inner',validate='many_to_one')
    # Index the reviewed evidence once, rather than repeatedly scanning millions of rows.
    grouped={k:v for k,v in sel.groupby(KEY+['stage','method','mode','library','term'],sort=False,dropna=False)}
    rows=[];den=[]
    for s in scope.itertuples(index=False):
        k=(s.analysis,s.context,s.target);ss=statuses[(statuses.analysis==k[0])&(statuses.context==k[1])&(statuses.target==k[2])]
        f=focus[(focus.analysis==k[0])&(focus.context==k[1])&(focus.target==k[2])]
        for stage,B in [('LODO',6),('BOOTSTRAP',20)]:
            st=ss[ss.stage==stage]
            assert len(st)==B*2
            paired=sum((q.status=='COMPLETED').all() for _,q in st.groupby('replicate'))
            den.append(dict(analysis=k[0],context=k[1],target=k[2],stage=stage,n_planned_KOs=2*B,n_completed_KOs=int((st.status=='COMPLETED').sum()),
                n_nonestimable_KOs=int((st.status!='COMPLETED').sum()),n_planned_pairs=B,n_estimable_pairs=paired))
            for p in f.itertuples(index=False):
                q=grouped.get(k+(stage,p.method,p.mode,p.library,p.term),sel.iloc[:0])
                ntest=hit=0;global_hit=0
                for r,z in q.groupby('replicate'):
                    assert z.seed.nunique()==len(z)
                    if len(z)!=2:continue
                    ntest+=1
                    passing=(z.main_FDR<CFG['fdr'])
                    if p.method=='ORA':passing&=z.overlap>=2
                    hit+=int(passing.all());global_hit+=int((passing&(z.FDR_exploratory_all_run<CFG['fdr'])).all())
                rows.append(dict(analysis=k[0],context=k[1],target=k[2],stage=stage,method=p.method,mode=p.mode,library=p.library,term=p.term,
                    n_planned=B,n_estimable_pairs=paired,n_testable=ntest,n_reproduced=hit,unconditional_rate=hit/B,
                    conditional_rate=hit/ntest if ntest else np.nan,n_reproduced_exploratory_global=global_hit))
    df=pd.DataFrame(rows);save(df,'reviewed_stability_by_term.tsv');save(pd.DataFrame(den),'resampling_denominators.tsv')
    old=pd.concat([table(SRC/'leave_one_donor_out_summary.tsv'),table(SRC/'donor_bootstrap_summary.tsv')])
    on=KEY+['stage','method','mode','library','term'];comparison=old.merge(df,on=on,validate='one_to_one',suffixes=('_old','_new'))
    assert len(comparison)==len(old)==len(df)
    assert (comparison.n_reproduced_old==comparison.n_reproduced_new).all() and (comparison.n_testable_old==comparison.n_testable_new).all()
    js('stability_summary_replay.json',dict(status='PASS',rows=len(df),all_original_denominators_and_counts_equal=True))
    return df,pd.DataFrame(den)

def ribosome(scope,cfg,libs):
    focus=table(ROOT/'results/2026-10-05-virtual-ko-pathways-034/two_seed_pathways.tsv').fillna('');rows=[]
    mapping=table(SRC/'enrichment_job_metadata.tsv').fillna('')
    for s in scope.itertuples(index=False):
        f=focus[(focus.analysis==s.analysis)&(focus.context==s.context)&(focus.target==s.target)]
        for variant in ['exclude_cytoplasmic','exclude_both']:
            tables={}
            for seed in CFG['network_seeds']:
                jid=f'{s.analysis}__{s.context}__{s.target}__{seed}__{variant}';d=SRC/'jobs'/jid
                tables[seed]=(table(d/'GSEA_all.tsv.gz'),table(d/'ORA_all.tsv.gz'),int(mapping.loc[mapping.job_id==jid,'n_removed'].iloc[0]))
            for p in f.itertuples(index=False):
                ntest=hit=0;removed=[];edges=[]
                for seed in CFG['network_seeds']:
                    g,o,n_removed=tables[seed];removed.append(n_removed)
                    if p.method=='GSEA':
                        q=g[(g.library==p.library)&(g['mode']==p.mode)&(g.pathway==p.term)]
                        if len(q):hit+=int(q.FDR_job_all_libraries.iloc[0]<CFG['fdr']);edges.append(q.leadingEdge.iloc[0])
                    else:
                        q=o[(o.library==p.library)&(o.term==p.term)]
                        if len(q):hit+=int(q.FDR_job_all_libraries.iloc[0]<CFG['fdr'] and q.overlap.iloc[0]>=2);edges.append(q.genes.iloc[0])
                    ntest+=int(len(q)==1)
                rows.append(dict(analysis=s.analysis,context=s.context,target=s.target,variant=variant,method=p.method,mode=p.mode,
                    library=p.library,term=p.term,n_testable_seeds=ntest,n_significant_seeds=hit,
                    outcome='TERM_NOT_TESTABLE' if ntest<2 else 'RETAINED_BOTH_SEEDS' if hit==2 else 'NO_LONGER_SHARED_SIGNIFICANT',
                    removed_seed42=removed[0],removed_seed43=removed[1],contributors=' | '.join(map(str,edges)),
                    annotation_limit='frozen_KEGG_Reactome_operational_proxy_not_pure_structural_ribosome'))
    save(pd.DataFrame(rows),'ribosome_focus_review.tsv')
    return pd.DataFrame(rows)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--fixtures',action='store_true');parser.add_argument('--resume-audit',action='store_true');args=parser.parse_args()
    if args.fixtures: print(json.dumps(fixtures()));return
    cfg,scope,reg,meta,tasks,statuses,libs=init()
    print('INPUT_AND_TASK_CLOSURE_PASS',flush=True)
    if args.resume_audit:
        assert json.loads((OUT/'numerical_audit.json').read_text())['status']=='PASS'
        manifest=OUT/'reviewed_source_artifact_manifest.tsv'
        assert manifest.exists(),'A complete source fingerprint is required before audit reuse.'
        for x in table(manifest).itertuples(index=False):
            assert sha(ROOT/x.path)==x.sha256, f'Reviewed source changed: {x.path}; rerun full audit.'
        selected=table(OUT/'focus_term_evidence_all_jobs.tsv.gz')
    else:selected=audit(cfg,reg,libs)
    ids=identity();print('IDENTITY_DESCRIPTIVE_REVIEW_COMPLETE',flush=True)
    quality=matching(scope,cfg,libs);print('MATCHING_REFERENCE_REVIEW_COMPLETE',flush=True)
    stable,den=stability(scope,statuses,selected,meta);print('STABILITY_REPLAY_PASS',flush=True)
    ribo=ribosome(scope,cfg,libs);print('RIBOSOME_FOCUS_REVIEW_COMPLETE',flush=True)
    js('review_tables_state.json',dict(status='EVIDENCE_TABLES_COMPLETE_PENDING_INTERPRETATION',scope_combinations=29,unique_targets=15,
        numerical_audit='PASS',stability_replay='PASS',direction_estimated=False,identity_validated=False))

if __name__=='__main__':main()

"""Durable, single-network-worker execution of SOP035 with isolated outputs.

Use an activated existing R environment. Original 033/034 scripts are never run.
"""
from pathlib import Path
import argparse, csv, hashlib, itertools, json, os, subprocess, sys, time
import numpy as np
import pandas as pd
import yaml
from scipy.stats import spearmanr, hypergeom

ROOT=Path(__file__).resolve().parents[1]
CONFIG=ROOT/'configs/virtual_ko_verification_execution_035.yaml'
CFG=yaml.safe_load(CONFIG.read_text(encoding='utf-8'))
OUT=ROOT/CFG['output_dir']; INTERIM=ROOT/CFG['interim_dir']
R='D:/conda/envs/zhang2026-nicotine-copd-r/lib/R/bin/Rscript.exe'
SCOPE=pd.read_csv(ROOT/CFG['scope_path'],sep='\t')

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def write_json(path,x):
    path=Path(path);tmp=path.with_name(path.name+'.partial')
    tmp.write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8');tmp.replace(path)

def state(stage,**kwargs):
    write_json(OUT/'run_state.json',dict(stage=stage,updated_at=time.strftime('%Y-%m-%dT%H:%M:%S%z'),pid=os.getpid(),**kwargs))
    print(stage,kwargs,flush=True)

def rrun(name,*args):
    log=OUT/'logs'/f'{name}.log';log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('ab') as f:
        p=subprocess.Popen([R,*map(str,args)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
        state(name,status='RUNNING',child_pid=p.pid,log=str(log.relative_to(ROOT)))
        while p.poll() is None:
            time.sleep(5)
        code=p.wait()
    if code:raise RuntimeError(f'{name} failed ({code}); see {log}')
    state(name,status='COMPUTATION_FINISHED_PENDING_AUDIT')

def freeze():
    OUT.mkdir(parents=True,exist_ok=True);INTERIM.mkdir(parents=True,exist_ok=True)
    assert len(SCOPE)==29 and SCOPE.target.nunique()==15
    audit=json.loads((ROOT/'docs/protocol/2026-10-05-virtual-ko-verification-035/scope_audit.json').read_text(encoding='utf-8'))
    for p,h in audit['source_hashes'].items():assert sha(ROOT/p)==h,p
    files=[CONFIG,ROOT/CFG['scope_path'],ROOT/CFG['network_config'],ROOT/CFG['retention_config'],
           ROOT/'src/verify_virtual_ko_035.R',ROOT/'src/enrich_virtual_ko_verification_035.R',Path(__file__)]
    files+=list((ROOT/CFG['input_dir']).glob('*_raw_counts.rds'))
    files+=list((ROOT/CFG['library_dir']).glob('*.gmt'))
    files.append(ROOT/'data/interim/2026-10-03-public-mechanism-analysis-026/frozen_modules.tsv')
    files += [ROOT/x for x in ['src/130_virtual_ko_cached_backend.R','src/135_sparse_cp_equivalent.R','src/138_implicit_manifold.R','external/virtual_ko_033/pcNet_upstream.R']]
    for (a,context),group in SCOPE.groupby(['analysis','context']):
        original='results/2026-10-04-stratified-virtual-ko-033' if a=='primary' else 'results/2026-10-05-virtual-ko-target-retention-sensitivity-033'
        files += [ROOT/original/context/f'seed_{seed}/WT_network.rds' for seed in [42,43]]
    frozen={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in files}
    p=OUT/'frozen_input_hashes.json'
    if p.exists():assert json.loads(p.read_text(encoding='utf-8'))==frozen,'Run inputs/code changed: record a new run or documented invalidation.'
    else:write_json(p,frozen)
    (OUT/'frozen_params.yaml').write_text(yaml.safe_dump(CFG,allow_unicode=True,sort_keys=False),encoding='utf-8')
    trace=ROOT/'workspaces'/CFG['run_id']/'raw/traces/events.jsonl';trace.parent.mkdir(parents=True,exist_ok=True)
    with trace.open('a',encoding='utf-8') as f:f.write(json.dumps(dict(event='REAL_EXECUTION_STARTED_OR_RESUMED',bootstrap=20,matched=20,new_environment=False,skill_training=False))+'\n')

def baseline_registry():
    old=pd.read_csv(ROOT/'results/2026-10-05-virtual-ko-pathways-034/job_registry.tsv',sep='\t')
    merged=old.merge(SCOPE[['analysis','context','target']],on=['analysis','context','target'],validate='many_to_one')
    assert len(merged)==58
    return merged

def a3():
    jobs=baseline_registry();effects={};stats=[];contrib=[]
    pathways=pd.read_csv(ROOT/'results/2026-10-05-virtual-ko-pathways-034/two_seed_pathways.tsv',sep='\t').fillna('')
    selected=set(zip(SCOPE.analysis,SCOPE.context,SCOPE.target))
    for j in jobs.itertuples(index=False):
        assert sha(ROOT/j.path)==j.sha256
        d=pd.read_csv(ROOT/j.path,sep='\t');d=d[d.gene!=j.target].copy();effects[(j.analysis,j.context,j.target,j.seed)]=d.set_index('gene')
    for key,df in effects.items():
        analysis,context,target,seed=key
        pp=pathways[(pathways.analysis==analysis)&(pathways.context==context)&(pathways.target==target)]
        for p in pp.itertuples(index=False):
            text=getattr(p,f'contributors{seed}'); genes=set(str(text).split(';'))-{''}
            for gene in sorted(genes):
                assert gene in df.index
                contrib.append(dict(analysis=analysis,context=context,target=target,seed=seed,library=p.library,term=p.term,method=p.method,mode=p.mode,gene=gene,distance=df.loc[gene,'distance'],gene_FDR=df.loc[gene,'p.adj'],gene_significant=df.loc[gene,'p.adj']<.05,expression_direction='NOT_ESTIMATED'))
    for (analysis,context),g in SCOPE.groupby(['analysis','context']):
        for t1,t2 in itertools.combinations(sorted(g.target),2):
            for seed in (42,43):
                x=effects[(analysis,context,t1,seed)];y=effects[(analysis,context,t2,seed)];u=x.index.intersection(y.index)
                s1=set(x.loc[u].index[x.loc[u,'p.adj']<.05]);s2=set(y.loc[u].index[y.loc[u,'p.adj']<.05]);union=s1|s2
                stats.append(dict(analysis=analysis,context=context,target1=t1,target2=t2,seed=seed,common_background=len(u),n_significant1=len(s1),n_significant2=len(s2),n_shared=len(s1&s2),jaccard=len(s1&s2)/len(union) if union else np.nan,distance_spearman=spearmanr(x.loc[u,'distance'],y.loc[u,'distance']).statistic,shared_genes=';'.join(sorted(s1&s2))))
    pd.DataFrame(contrib).to_csv(OUT/'target_contributor_audit.tsv',sep='\t',index=False)
    pd.DataFrame(stats).to_csv(OUT/'cross_target_similarity.tsv',sep='\t',index=False)
    # Side-by-side continuous marker diagnostics, not an automatic relabeling.
    cells=pd.read_csv(OUT/'cell_identity_audit.tsv',sep='\t');q=[]
    for context,g in cells.groupby('context'):
        original=context.split('_')[0];alternative='Microglia' if original=='Astrocyte' else 'Astrocyte'
        q.append(dict(context=context,n_cells=len(g),original_panel_mean=g[original].mean(),alternative_panel_mean=g[alternative].mean(),fraction_alternative_panel_higher=(g[alternative]>g[original]).mean(),fraction_both_panels_at_least_two_markers=((g.Astrocyte_n_detected_markers>=2)&(g.Microglia_n_detected_markers>=2)).mean(),interpretation='DESCRIPTIVE_UNCALIBRATED_PANEL_NO_AUTOMATIC_DOUBLETS_OR_RELABEL'))
    pd.DataFrame(q).to_csv(OUT/'cell_identity_context_summary.tsv',sep='\t',index=False)
    write_json(OUT/'A_existing_data_audit.json',dict(status='COMPUTATIONS_COMPLETE_WITH_IDENTITY_REVIEW_LIMITATIONS',baseline_KOs=58,contributor_rows=len(contrib),similarity_rows=len(stats),doublet_ambient_full_QC='SEE_library_QC_availability.tsv',direction_estimated=False))

def ribosomes():
    sets={}
    for lib in CFG['libraries']:
        sets[lib]={x[0]:(x[1],set(x[2:])) for x in (line.rstrip('\n').split('\t') for line in (ROOT/CFG['library_dir']/f'{lib}.gmt').read_text(encoding='utf-8').splitlines())}
    rib=sets['KEGG']['hsa03010'][1]
    mito=set().union(*(v[1] for v in sets['Reactome'].values() if 'mitochondrial translation' in v[0].lower()))&rib
    # KEGG includes rRNAs/unresolved loci and some assembly components; preserve that limitation.
    protein={g for g in rib if not g.startswith(('RNA','RNR','LOC'))}
    mito=protein&mito;cyto=protein-mito
    assert len(mito)>20 and len(cyto)>40
    pd.DataFrame([dict(gene=g,compartment='mitochondrial' if g in mito else 'cytoplasmic_or_unresolved_KEGG_protein',source='frozen KEGG hsa03010; mitochondrial translation Reactome overlap',limitation='KEGG protein-gene proxy may include biogenesis components; not pure structural-complex annotation') for g in sorted(protein)]).to_csv(OUT/'ribosome_gene_sets.tsv',sep='\t',index=False)
    jobs=baseline_registry();new=[];mapping=[]
    for j in jobs.itertuples(index=False):
        dr=pd.read_csv(ROOT/j.path,sep='\t')
        for variant,excluded in [('exclude_cytoplasmic',cyto),('exclude_both',protein)]:
            d=dr[~dr.gene.isin(excluded)].copy();assert j.target not in excluded
            p=INTERIM/f'{j.job_id}__{variant}.tsv';d.to_csv(p,sep='\t',index=False)
            job=dict(j._asdict());job['job_id']=f'{j.job_id}__{variant}';job['path']=str(p.relative_to(ROOT)).replace('\\','/');job['sha256']=sha(p)
            new.append(job);mapping.append(dict(job_id=job['job_id'],stage='RIBOSOME',variant=variant,original_job_id=j.job_id,replicate=0,task_id='baseline',n_removed=int(dr.gene.isin(excluded).sum())))
    return new,mapping

def register(rows,mapping):
    p=OUT/'job_registry.tsv';mp=OUT/'enrichment_job_metadata.tsv'
    old=pd.read_csv(p,sep='\t') if p.exists() else pd.DataFrame()
    new=pd.DataFrame(rows)
    if len(new):
        new['enrichment_seed']=[CFG['seed']+int(hashlib.sha256(x.encode()).hexdigest()[:7],16) for x in new.job_id]
        if len(old):
            overlap=old.merge(new,on='job_id',suffixes=('_old','_new'));assert (overlap.sha256_old==overlap.sha256_new).all()
            new=new[~new.job_id.isin(old.job_id)]
        pd.concat([old,new],ignore_index=True).to_csv(p,sep='\t',index=False)
    oldm=pd.read_csv(mp,sep='\t') if mp.exists() else pd.DataFrame()
    m=pd.DataFrame(mapping)
    if len(m):pd.concat([oldm,m[~m.job_id.isin(oldm.job_id)] if len(oldm) else m],ignore_index=True).to_csv(mp,sep='\t',index=False)

def task_rows(task):
    d=OUT/'networks'/task.task_id
    statuses=pd.read_csv(d/'KO_status.tsv',sep='\t').fillna('');rows=[];mm=[]
    for r in statuses.itertuples(index=False):
        if r.status!='COMPLETED':continue
        rows.append(dict(job_id=f'{task.task_id}__{r.target}__{r.seed}',analysis=task.analysis,context=task.context,cell_type=task.context.split('_')[0],state=task.context.split('_')[1],target=r.target,seed=r.seed,path=r.path,sha256=sha(ROOT/r.path)))
        mm.append(dict(job_id=rows[-1]['job_id'],stage=task.stage,variant='',original_job_id=f'{task.analysis}__{task.context}__{r.target}__{r.seed}',replicate=task.replicate,task_id=task.task_id,n_removed=0))
    return rows,mm

def summarize():
    meta=pd.read_csv(OUT/'enrichment_job_metadata.tsv',sep='\t').fillna('');reg=pd.read_csv(OUT/'job_registry.tsv',sep='\t');summary=[];ribo=[];replay=[];lookup={};module_rows=[]
    focus=pd.read_csv(ROOT/'results/2026-10-05-virtual-ko-pathways-034/two_seed_pathways.tsv',sep='\t').fillna('')
    fm=pd.read_csv(ROOT/'data/interim/2026-10-03-public-mechanism-analysis-026/frozen_modules.tsv',sep='\t')
    modules=fm.groupby('module').gene.apply(set).to_dict()
    def bh(v):
        v=np.asarray(v,float);out=np.full(len(v),np.nan);ix=np.flatnonzero(np.isfinite(v));order=ix[np.argsort(v[ix],kind='stable')]
        out[order]=np.minimum(1,np.minimum.accumulate((v[order]*len(order)/np.arange(1,len(order)+1))[::-1])[::-1]);return out
    # Every per-job output is audited for hashes, finite P and gene/background membership.
    for j in reg.itertuples(index=False):
        d=OUT/'jobs'/j.job_id;a=json.loads((d/'completion.json').read_text())
        assert a['status']=='PASS' and a['source_sha256']==sha(ROOT/j.path)
        g=pd.read_csv(d/'GSEA_all.tsv.gz',sep='\t');o=pd.read_csv(d/'ORA_all.tsv.gz',sep='\t')
        assert g.pval.notna().all() and o.P.notna().all()
        dr=pd.read_csv(ROOT/j.path,sep='\t');dr=dr[dr.gene!=j.target];bg=set(dr.gene)
        assert set(g.leadingEdge.fillna('').str.split(';').explode())-{''}<=bg
        assert not g.duplicated(['library','mode','pathway']).any()
        m=meta[meta.job_id==j.job_id].iloc[0]
        fg=focus[(focus.analysis==j.analysis)&(focus.context==j.context)&(focus.target==j.target)]
        fk={(p.method,p.mode,p.library,p.term) for p in fg.itertuples(index=False)}
        tested=set();passed=set()
        for p in g.itertuples(index=False):
            key=('GSEA',p.mode,p.library,p.pathway)
            if key in fk:
                tested.add(key)
                if p.FDR_job_all_libraries<.05:passed.add(key)
        for p in o.itertuples(index=False):
            key=('ORA','significant_gene_set',p.library,p.term)
            if key in fk:
                tested.add(key)
                if p.FDR_job_all_libraries<.05 and p.overlap>=2:passed.add(key)
        lookup[j.job_id]=(tested,passed)
        sig=set(dr.loc[dr['p.adj']<.05,'gene'])
        for name,genes in modules.items():
            if name.startswith('G5_') != (j.cell_type=='Astrocyte'):continue
            members=(genes-{j.target})&bg;hits=members&sig
            module_rows.append(dict(job_id=j.job_id,stage=m.stage,analysis=j.analysis,context=j.context,target=j.target,seed=j.seed,module=name,n_members=len(members),hits=';'.join(sorted(hits)),P=hypergeom.sf(len(hits)-1,len(bg),len(members),len(sig)) if len(members)>=2 else np.nan))
        significant=g[(g['mode']=='weighted')&(g.FDR_job_all_libraries<.05)]
        summary.append(dict(job_id=j.job_id,analysis=j.analysis,context=j.context,target=j.target,seed=j.seed,stage=m.stage,replicate=m.replicate,weighted_significant_terms=len(significant),weighted_terms=';'.join(significant.pathway),n_tested=len(dr),n_significant=(dr['p.adj']<.05).sum()))
        if m.stage=='RIBOSOME':ribo.append(dict(job_id=j.job_id,analysis=j.analysis,context=j.context,target=j.target,seed=j.seed,variant=m.variant,n_removed=m.n_removed,n_tested=len(dr),weighted_terms=';'.join(significant.pathway)))
    pd.DataFrame(summary).to_csv(OUT/'enrichment_task_summary.tsv',sep='\t',index=False)
    pd.DataFrame(ribo).to_csv(OUT/'ribosome_sensitivity.tsv',sep='\t',index=False)
    # Denominators include low-expression/no-edge outcomes; original focus terms stay fixed.
    statuses=[]
    for p in (OUT/'networks').glob('*/KO_status.tsv'):statuses.append(pd.read_csv(p,sep='\t'))
    st=pd.concat(statuses,ignore_index=True);tasks=pd.read_csv(OUT/'network_task_registry.tsv',sep='\t');st=st.merge(tasks[['task_id','analysis','context','stage','replicate']],on='task_id')
    st.to_csv(OUT/'all_resampling_KO_statuses.tsv',sep='\t',index=False)
    accum=[]
    for r in SCOPE.itertuples(index=False):
        f=focus[(focus.analysis==r.analysis)&(focus.context==r.context)&(focus.target==r.target)]
        for stage,B in [('LODO',6),('BOOTSTRAP',CFG['bootstrap_replicates'])]:
            rr=st[(st.analysis==r.analysis)&(st.context==r.context)&(st.target==r.target)&(st.stage==stage)]
            for p in f.itertuples(index=False):
                tested=hit=0
                for b in range(1,B+1):
                    bb=rr[rr.replicate==b]
                    if len(bb)!=2 or not (bb.status=='COMPLETED').all():continue
                    passes=[];can=[]
                    for k in bb.itertuples(index=False):
                        jid=f'{k.task_id}__{k.target}__{k.seed}';key=(p.method,p.mode,p.library,p.term)
                        can.append(key in lookup[jid][0]);passes.append(key in lookup[jid][1])
                    tested+=all(can);hit+=all(passes)
                accum.append(dict(analysis=r.analysis,context=r.context,target=r.target,stage=stage,method=p.method,mode=p.mode,library=p.library,term=p.term,n_planned=B,n_testable=tested,n_reproduced=hit,unconditional_rate=hit/B,conditional_rate=hit/tested if tested else np.nan))
    aa=pd.DataFrame(accum)
    for stage,name in [('LODO','leave_one_donor_out_summary.tsv'),('BOOTSTRAP','donor_bootstrap_summary.tsv')]:aa[aa.stage==stage].to_csv(OUT/name,sep='\t',index=False)
    mr=pd.DataFrame(module_rows);mr['FDR_new_run_all_module_tests']=bh(mr.P);mr.to_csv(OUT/'resampling_module_enrichment.tsv',sep='\t',index=False)
    # Include module-only candidates in the final record rather than silently omitting them.
    conclusions=[]
    for r in SCOPE.itertuples(index=False):
        rr=st[(st.analysis==r.analysis)&(st.context==r.context)&(st.target==r.target)]
        bb=aa[(aa.analysis==r.analysis)&(aa.context==r.context)&(aa.target==r.target)]
        conclusions.append(dict(analysis=r.analysis,context=r.context,target=r.target,layers=r.layers,original_grade=r.evidence_grade,n_planned_KOs=len(rr),n_completed_KOs=(rr.status=='COMPLETED').sum(),n_nonestimable_KOs=(rr.status!='COMPLETED').sum(),maximum_shared_weighted_LODO_rate=bb[(bb.stage=='LODO')&(bb['mode']=='weighted')].unconditional_rate.max(),maximum_shared_weighted_bootstrap_rate=bb[(bb.stage=='BOOTSTRAP')&(bb['mode']=='weighted')].unconditional_rate.max(),original_modules=r.shared_predefined_modules,conclusion='COMPUTATIONAL_EVIDENCE_TABLES_READY_HUMAN_IDENTITY_AND_SPECIFICITY_REVIEW_REQUIRED',expression_direction='NOT_ESTIMATED'))
    pd.DataFrame(conclusions).to_csv(OUT/'target_verification_conclusions.tsv',sep='\t',index=False)
    return len(reg)

def run():
    assert os.environ.get('CONDA_PREFIX','').replace('\\','/').lower()=='d:/conda/envs/zhang2026-nicotine-copd-r'
    freeze()
    rrun('adapter_tests','src/verify_virtual_ko_035.R',str(CONFIG.relative_to(ROOT)),'test')
    for p in ['adapter_deterministic_tests.json','adapter_real_replay.json']:
        a=json.loads((OUT/p).read_text());assert a['status']=='PASS' and a['adapter_sha256']==sha(ROOT/'src/verify_virtual_ko_035.R')
    if not (OUT/'preparation.json').exists():rrun('prepare_existing_data','src/verify_virtual_ko_035.R',str(CONFIG.relative_to(ROOT)),'prepare')
    a3();state('A_EXISTING_DATA_DONE',status='COMPLETE_WITH_QC_INFORMATION_LIMITATIONS')
    rows,mm=ribosomes();register(rows,mm)
    # Prioritize complete donor task registration/execution; new enrichment is sequential and resumable.
    tasks=pd.read_csv(OUT/'network_task_registry.tsv',sep='\t')
    for t in tasks.itertuples(index=False):
        marker=OUT/'networks'/t.task_id/'completion.json'
        # Always call the adapter to verify cached signatures, rather than trusting file presence.
        rrun(t.task_id,'src/verify_virtual_ko_035.R',str(CONFIG.relative_to(ROOT)),'network',t.task_id)
        assert json.loads(marker.read_text())['status']=='PASS';register(*task_rows(t))
    for (analysis,context),group in SCOPE.groupby(['analysis','context']):
        rrun(f'match_{analysis}_{context}','src/verify_virtual_ko_035.R',str(CONFIG.relative_to(ROOT)),'match',analysis,context)
        d=OUT/'matched'/f'{analysis}__{context}';st=pd.read_csv(d/'KO_status.tsv',sep='\t');rr=[];mm=[]
        for r in st.itertuples(index=False):
            jid=f'match__{analysis}__{context}__{r.target}__{r.seed}'
            rr.append(dict(job_id=jid,analysis=analysis,context=context,cell_type=context.split('_')[0],state=context.split('_')[1],target=r.target,seed=r.seed,path=r.path,sha256=sha(ROOT/r.path)))
            mm.append(dict(job_id=jid,stage='MATCHED',variant='',original_job_id='',replicate=0,task_id='baseline',n_removed=0))
        register(rr,mm)
    rrun('enrich_all_registered_jobs','src/enrich_virtual_ko_verification_035.R',str(CONFIG.relative_to(ROOT)),1,1)
    n=summarize()
    for p,h in json.loads((OUT/'frozen_input_hashes.json').read_text()).items():assert sha(ROOT/p)==h
    state('COMPUTATIONS_AUDITED_PENDING_SCIENTIFIC_REVIEW',status='REQUIRES_FINAL_REVIEW',enrichment_jobs=n,direction_estimated=False)

if __name__=='__main__':
    try:run()
    except Exception as exc:
        OUT.mkdir(parents=True,exist_ok=True);state('FAILED_REQUIRES_FIX',status='FAILED',error=str(exc));raise

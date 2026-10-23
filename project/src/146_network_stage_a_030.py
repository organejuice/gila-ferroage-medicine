"""Degree-matched proximity using frozen NetMedPy bins and sparse shortest paths."""
import ast, hashlib, json, random, time
from pathlib import Path
import networkx as nx
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import dijkstra
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1];C=json.loads((ROOT/'configs/stage_a_030.json').read_text());RUN=C['run_id']
RAW=ROOT/'data/raw'/RUN;INT=ROOT/'data/interim'/RUN;OUT=ROOT/'results'/RUN
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
def write(d,p):d.to_csv(p,sep='\t',index=False,compression={'method':'gzip','mtime':0} if str(p).endswith('.gz') else None)
def bh(p):
    a=np.asarray(p,float);out=np.full(len(a),np.nan);idx=np.where(np.isfinite(a))[0];idx=idx[np.argsort(a[idx],kind='stable')]
    if len(idx):out[idx]=np.minimum(1,np.minimum.accumulate((a[idx]*len(idx)/np.arange(1,len(idx)+1))[::-1])[::-1])
    return out
def load_binner():
    source=(RAW/'NetMedPy/netmedpy/netmedpy.py').read_text(encoding='utf8');tree=ast.parse(source)
    funcs=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=='_get_degree_binning'];assert len(funcs)==1
    namespace={'nx':nx};exec(compile(ast.Module(body=funcs,type_ignores=[]),'frozen_NetMedPy_degree_binning','exec'),namespace)
    return namespace['_get_degree_binning']
def degree_groups(graph,bin_min):
    _,_,bins=load_binner()(graph,bin_min)
    bucket={node:i for i,values in enumerate(bins) for node in values}
    assert len(bucket)==graph.number_of_nodes()
    return [np.array(values,dtype=np.int32) for values in bins],bucket
def sample_set(template,bins,bucket,rng):
    counts={}
    for n in template:counts[bucket[int(n)]]=counts.get(bucket[int(n)],0)+1
    x=np.concatenate([rng.choice(bins[b],size=k,replace=False) for b,k in sorted(counts.items())])
    assert len(set(x))==len(template)
    return x

def sample_plan(template,bucket):
    values,counts=np.unique([bucket[int(n)] for n in template],return_counts=True)
    return list(zip(values,counts))

def draw_plan(plan,bins,rng):
    # Bins are disjoint; each draw is without replacement within its bin.
    return np.concatenate([rng.choice(bins[b],size=k,replace=False) for b,k in plan])
def graph_distance(adj,seeds):
    return dijkstra(adj,directed=False,indices=np.asarray(seeds,dtype=np.int32),unweighted=True,min_only=True)
def prepare_edges():
    dest=INT/'STRING_approved_edges.tsv.gz'
    if dest.exists():return read(dest)
    approved=set(read(ROOT/'data/raw/2026-09-30-fixed-formula-target-supplement-024/HGNC/hgnc_complete_set.tsv').query("status == 'Approved'").symbol)
    info=read(RAW/'STRING/9606.protein.info.v12.0.txt.gz');pid=info.columns[0]
    mapping=dict(zip(info[pid],info.preferred_name));approved_map={p:g for p,g in mapping.items() if g in approved}
    info['HGNC_approved']=info.preferred_name.isin(approved);write(info,OUT/'STRING_node_identifier_audit.tsv.gz')
    pieces=[]
    for part in pd.read_csv(RAW/'STRING/9606.protein.links.full.v12.0.txt.gz',sep=' ',chunksize=300000,usecols=['protein1','protein2','experiments','textmining','combined_score']):
        part=part[(part.experiments>=C['network']['threshold'])|(part.textmining>=C['network']['sensitivity_threshold'])].copy()
        part['a']=part.protein1.map(approved_map);part['b']=part.protein2.map(approved_map);part=part.dropna(subset=['a','b']);part=part[part.a!=part.b]
        part[['a','b']]=np.sort(part[['a','b']].values,axis=1);pieces.append(part[['a','b','experiments','textmining','combined_score']])
    edges=pd.concat(pieces).groupby(['a','b'],as_index=False)[['experiments','textmining','combined_score']].max();write(edges,dest)
    return edges
def run_network(name,edges,sets,drugsets,full_modules):
    print('network_start',name,len(edges),flush=True);g=nx.from_pandas_edgelist(edges,'a','b');lcc=max(nx.connected_components(g),key=len);g=g.subgraph(lcc).copy()
    nodes=sorted(g);idx={g:i for i,g in enumerate(nodes)};g=nx.relabel_nodes(g,idx,copy=True);adj=nx.to_scipy_sparse_array(g,nodelist=range(len(nodes)),format='csr',dtype=np.float64)
    adj=csr_matrix((adj.data,adj.indices.astype(np.int32),adj.indptr.astype(np.int32)),shape=adj.shape)
    bins,bucket=degree_groups(g,C['network']['degree_bin_min']);rng=np.random.default_rng(C['seed']);summary=[];nullrows=[];cover=[];binrows=[]
    for b,values in enumerate(bins):binrows.append(dict(network=name,bin=b,nodes=len(values),degree_min=min(dict(g.degree(values)).values()),degree_max=max(dict(g.degree(values)).values())))
    for interface,genes in sets.items():
        if not full_modules and interface not in ['G1_robust244','G2_GATE_UNION','G3_GATE_UNION','G4_GATE_UNION','G5_GATE_UNION']:continue
        S=np.array([idx[n] for n in sorted(genes) if n in idx],dtype=np.int32)
        cover.append(dict(network=name,side='disease',set=interface,original=len(genes),mapped=len(S),lost=';'.join(sorted(set(genes)-set(idx)))))
        if not len(S):continue
        observed=graph_distance(adj,S);pairs=[]
        for tier,genesT in drugsets.items():
            T=np.array([idx[n] for n in sorted(genesT) if n in idx],dtype=np.int32)
            cover.append(dict(network=name,side='drug',set=tier,original=len(genesT),mapped=len(T),lost=';'.join(sorted(set(genesT)-set(idx)))))
            for overlap_mode in ['included','removed']:
                Tuse=T if overlap_mode=='included' else np.setdiff1d(T,S).astype(np.int32)
                if len(Tuse):pairs.append(dict(tier=tier,T=Tuse,mode=overlap_mode,obs=float(observed[Tuse].mean()),overlap=len(np.intersect1d(T,S)),null=[]))
        splan=sample_plan(S,bucket)
        for pair in pairs:pair['plan']=sample_plan(pair['T'],bucket)
        maxiter=C['network']['null_iterations'];iteration=0
        while iteration<maxiter:
            Sr=draw_plan(splan,bins,rng);dist=graph_distance(adj,Sr)
            for pair in pairs:
                Tr=draw_plan(pair['plan'],bins,rng);pair['null'].append(float(dist[Tr].mean()))
            iteration+=1
            if iteration==C['network']['null_iterations']:
                for pair in pairs:
                    v=np.array(pair['null']);p=(1+(v<=pair['obs']).sum())/(len(v)+1);p1=(v[:len(v)//2]<=pair['obs']).mean();p2=(v[len(v)//2:]<=pair['obs']).mean()
                    pair['mc_split_gap_initial']=abs(p1-p2)
                    if .02<=p<=.08 or abs(p1-p2)>.02:maxiter=C['network']['extended_iterations']
        for pair in pairs:
            v=np.array(pair['null']);p=(1+(v<=pair['obs']).sum())/(len(v)+1);sd=v.std(ddof=1)
            summary.append(dict(network=name,interface=interface,tier=pair['tier'],overlap_mode=pair['mode'],drug_nodes=len(pair['T']),disease_nodes=len(S),original_overlap=pair['overlap'],distance=pair['obs'],null_mean=v.mean(),null_sd=sd,z=(pair['obs']-v.mean())/sd if sd>0 else np.nan,P_empirical=p,null_iterations=len(v),MC_SE=np.sqrt(p*(1-p)/(len(v)+1)),MC_split_gap_initial=pair['mc_split_gap_initial'],interpretation='NETWORK_ASSOCIATION_NOT_EFFICACY',direct_binding_claim=False))
            nullrows.extend(dict(network=name,interface=interface,tier=pair['tier'],overlap_mode=pair['mode'],iteration=i+1,distance=d) for i,d in enumerate(v))
        print('network_interface',name,interface,'iterations',maxiter,flush=True)
    write(pd.DataFrame(nullrows),OUT/f'network_null_{name}.tsv.gz');write(pd.DataFrame(cover),OUT/f'network_coverage_{name}.tsv');write(pd.DataFrame(binrows),OUT/f'network_degree_bins_{name}.tsv')
    write(pd.DataFrame([dict(network=name,nodes=len(nodes),edges=g.number_of_edges(),original_nodes=len(set(edges.a)|set(edges.b)),retained_component='largest connected',species='human STRING12.0',mapping='exact HGNC approved preferred name')]),OUT/f'network_QC_{name}.tsv')
    result=pd.DataFrame(summary);result['FDR_network']=bh(result.P_empirical);write(result,OUT/f'network_proximity_{name}.tsv');return result
def main():
    import sys
    only=sys.argv[sys.argv.index('--variant-only')+1] if '--variant-only' in sys.argv else None
    assemble='--assemble-existing' in sys.argv
    seal=OUT/'network_cache_manifest.json'
    def digest(p):
        h=hashlib.sha256()
        with p.open('rb') as f:
            for b in iter(lambda:f.read(2**20),b''):h.update(b)
        return h.hexdigest()
    dependencies=[INT/'all_drug_targets.tsv',INT/'D0_interface_genes.tsv',INT/'G1_full_meta.tsv',INT/'STRING_approved_edges.tsv.gz',RAW/'NetMedPy/netmedpy/netmedpy.py']
    tree=ast.parse(Path(__file__).read_text(encoding='utf8'))
    core=ast.dump(ast.Module(body=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name!='main'],type_ignores=[]))
    signature=dict(parameters=C['network'],seed=C['seed'],fdr=C['fdr'],core_sha256=hashlib.sha256(core.encode()).hexdigest(),full_source_sha256=digest(Path(__file__)),inputs={str(p.relative_to(ROOT)):digest(p) for p in dependencies if p.exists()})
    outputs=list(OUT.glob('network_*.tsv'))+list(OUT.glob('network_null_*.tsv.gz'))
    if '--seal-existing-verified' in sys.argv:
        audit=json.loads((OUT/'acceptance_audit.json').read_text())
        if audit['status']!='PASS':raise ValueError('cannot seal unaudited outputs')
        seal.write_text(json.dumps(dict(signature=signature,outputs={str(p.relative_to(ROOT)):digest(p) for p in outputs}),indent=2)+'\n')
        print('audited_network_cache_sealed');return
    if not only and not assemble and any(OUT.glob('network_proximity_*.tsv')):
        if not seal.exists():raise ValueError('existing network outputs have no validated cache manifest; explicitly recompute or audit and seal')
        frozen=json.loads(seal.read_text())
        if frozen['signature']!=signature or any(digest(ROOT/p)!=h for p,h in frozen['outputs'].items()):raise ValueError('network cache input/parameter/code/output hash mismatch; recompute required')
    edges=prepare_edges();members=read(INT/'D0_interface_genes.tsv');sets={k:set(x.gene) for k,x in members.groupby('module')}
    meta=read(INT/'G1_full_meta.tsv');sets['G1_significant413_exploratory']=set(meta.loc[meta.meta_FDR<C['fdr'],'gene'])
    for gate in ['G2','G3','G4','G5']:sets[gate+'_GATE_UNION']=set().union(*(v for k,v in sets.items() if k.startswith(gate+'_')))
    drug=read(INT/'all_drug_targets.tsv');drugsets={k:set(x.gene) for k,x in drug.groupby('tier')}
    assert {k:len(v) for k,v in drugsets.items()}=={'E1':276,'E2':1244,'P':6461,'H':2081},'invalid full drug target input'
    variants=[('experiments400',edges[edges.experiments>=C['network']['threshold']],True),('experiments700',edges[edges.experiments>=C['network']['sensitivity_threshold']],False),('experiments400_text700',edges,False)]
    if only:
        variants=[v for v in variants if v[0]==only]
        if not variants:raise ValueError('unknown network variant')
    results=[]
    for name,x,allmodules in variants:
        path=OUT/f'network_proximity_{name}.tsv'
        if assemble and not path.exists():raise ValueError('all three completed network variants are required for assembly')
        results.append(read(path) if path.exists() and not only else run_network(name,x,sets,drugsets,allmodules))
    if only:print('network_variant_complete',only,flush=True);return
    full=pd.concat(results,ignore_index=True);full['FDR_all_network_tests']=bh(full.P_empirical);write(full,OUT/'network_proximity_all.tsv')
    print('network_complete',len(full),flush=True)
    outputs=list(OUT.glob('network_*.tsv'))+list(OUT.glob('network_null_*.tsv.gz'))
    seal.write_text(json.dumps(dict(signature=signature,outputs={str(p.relative_to(ROOT)):digest(p) for p in outputs}),indent=2)+'\n')
if __name__=='__main__':main()

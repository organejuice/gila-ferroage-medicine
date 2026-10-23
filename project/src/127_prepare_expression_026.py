"""Build expression/sample manifests and frozen disease modules for 026."""
import ast,csv,gzip,importlib.util,io,json,re
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-public-mechanism-analysis-026'
RAW=ROOT/'data/raw'/RUN;OLD=ROOT/'workspaces/2026-10-02-public-mechanism-bridge-search-025/raw/sources'
INT=ROOT/'data/interim'/RUN;OUT=ROOT/'results'/RUN
CFG=json.loads((ROOT/'configs/public_mechanism_analysis_026.json').read_text())

def read(p,**kw): return pd.read_csv(p,sep='\t',**kw)
def write(d,p): p.parent.mkdir(parents=True,exist_ok=True);d.to_csv(p,sep='\t',index=False)
def soft(p):
    rows=[];r=None
    with gzip.open(p,'rt',encoding='utf-8') as f:
        for l in f:
            if l.startswith('^SAMPLE = '):
                if r:rows.append(r)
                r={'GSM':l.strip().split(' = ')[1]}
            elif l.startswith('^'):
                if r:rows.append(r);r=None
            elif r and l.startswith('!Sample_') and ' = ' in l:
                k,v=l.strip().split(' = ',1);r.setdefault(k.replace('!Sample_',''),[]).append(v)
    if r:rows.append(r)
    return [{k:' | '.join(v) if isinstance(v,list) else v for k,v in r.items()} for r in rows]

def orthology(hgnc):
    dest=INT/'mouse_human_one_to_one.tsv'
    if dest.exists(): return dict(zip(read(dest,dtype=str)['mouse_symbol'],read(dest,dtype=str)['human_symbol']))
    folder=ROOT/'data/raw/NCBI_gene_orthologs'
    pairs=set()
    for chunk in pd.read_csv(folder/'gene_orthologs.gz',sep='\t',dtype=str,chunksize=500000):
        chunk.columns=[x.lstrip('#').lower().replace(' ','_') for x in chunk.columns]
        sub=chunk[(chunk.tax_id.isin(['10090','9606'])) & (chunk.other_tax_id.isin(['10090','9606'])) & (chunk.tax_id!=chunk.other_tax_id) & (chunk.relationship=='Ortholog')]
        for r in sub.itertuples(index=False):pairs.add((r.geneid,r.other_geneid) if r.tax_id=='10090' else (r.other_geneid,r.geneid))
    mh=defaultdict(set);hm=defaultdict(set)
    for m,h in pairs:mh[m].add(h);hm[h].add(m)
    mouse=read(folder/'Mus_musculus.gene_info.gz',dtype=str,usecols=['GeneID','Symbol'])
    mice=dict(zip(mouse.GeneID,mouse.Symbol));human=dict(zip(hgnc.entrez_id,hgnc.symbol))
    rows=[dict(mouse_gene_id=m,human_gene_id=h,mouse_symbol=mice[m],human_symbol=human[h]) for m,h in sorted(pairs) if len(mh[m])==1 and len(hm[h])==1 and m in mice and h in human]
    frame=pd.DataFrame(rows);frame=frame[~frame.mouse_symbol.duplicated(keep=False)]
    write(frame,dest);return dict(zip(frame.mouse_symbol,frame.human_symbol))

def main():
    INT.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    h=read(ROOT/'data/raw/2026-09-30-fixed-formula-target-supplement-024/HGNC/hgnc_complete_set.tsv',dtype=str).fillna('');h=h[h.status=='Approved']
    ens={k:v for k,v in zip(h.ensembl_gene_id,h.symbol) if k};approved=set(h.symbol);mouse=orthology(h)
    datasets=[];samples=[];mappings=[];qc=[]
    def add(name,d,cols,groups,kind,species,cell,genes=None,ids=None,gsm=None,limit='',formal=True):
        assert len(cols)==len(groups) and len(set(cols))==len(cols)
        x=d[cols].apply(pd.to_numeric,errors='raise');assert np.isfinite(x.to_numpy()).all()
        if kind!='microarray_log2': assert (x.to_numpy()>=0).all()
        if kind=='counts': assert np.max(np.abs(x.to_numpy()-np.round(x.to_numpy())))<1e-7,'noninteger counts '+name
        ids=ids if ids is not None else d.iloc[:,0].astype(str)
        ids=[str(g).split('.')[0] if str(g).startswith(('ENSG','ENSMUSG')) else str(g) for g in ids]
        genes=list(genes) if genes is not None else ids
        if len(set(ids))!=len(ids):
            groups_by_id=defaultdict(list)
            for j,gid in enumerate(ids):groups_by_id[gid].append(j)
            duplicate_rows=[]
            for gid,indices in groups_by_id.items():
                if len(indices)<2:continue
                if len({str(genes[j]) for j in indices})!=1 or not np.all(x.iloc[indices].to_numpy()==x.iloc[indices[0]].to_numpy()):raise ValueError('conflicting duplicate feature '+name+' '+gid)
                duplicate_rows.append(dict(dataset=name,feature_id=gid,annotation_rows=len(indices),reason='identical sample vector and gene symbol; count measurement retained once, not summed'))
            write(pd.DataFrame(duplicate_rows),OUT/(name+'_annotation_duplicates.tsv'))
            indices=[v[0] for v in groups_by_id.values()];x=x.iloc[indices].reset_index(drop=True);ids=[ids[j] for j in indices];genes=[genes[j] for j in indices]
        symbols=[]
        for gid,g in zip(ids,genes):
            g=str(g).split('|')[0]
            sym=(g if g in approved else ens.get(gid,'')) if species=='human' else mouse.get(g,'')
            symbols.append(sym);mappings.append(dict(dataset=name,feature_id=gid,original_gene=g,human_symbol=sym,status='exact_approved' if species=='human' and sym else 'reciprocal_one_to_one' if sym else 'unmapped_or_ambiguous'))
        x.insert(0,'feature_id',ids);write(x,INT/(name+'_expression.tsv.gz'))
        mapping=pd.DataFrame(dict(feature_id=ids,human_symbol=symbols));write(mapping,INT/(name+'_gene_map.tsv'))
        sample_rows=[dict(dataset=name,sample=c,group=g,GSM=(gsm or {}).get(c,''),species=species,cell=cell,independent_unit='deposited_experimental_replicate; donor/generalization not implied',limitation=limit) for c,g in zip(cols,groups)]
        write(pd.DataFrame(sample_rows),INT/(name+'_samples.tsv'));samples.extend(sample_rows)
        datasets.append(dict(id=name,kind=kind,species=species,cell=cell,formal=formal,limitation=limit,contrasts=CFG['contrasts'].get(name,[])))
        qc.append(dict(dataset=name,features=len(x),samples=len(cols),mapped_features=sum(bool(s) for s in symbols),mapped_human_genes=len(set(symbols)-{''}),duplicate_human_symbols=sum(pd.Series(symbols)[pd.Series(symbols)!=''].duplicated()),minimum=float(x[cols].min().min()),maximum=float(x[cols].max().max()),kind=kind))
    sm=read(ROOT/'results/2026-10-02-public-mechanism-bridge-search-025/sample_mapping.tsv').fillna('')
    def old_meta(acc):return sm[sm.series==acc].to_dict('records')
    def wb(name,sheet=1,**kw):return read(INT/'workbooks'/f'{name}_sheet{sheet}.tsv',escapechar='\\',**kw)
    d=read(OLD/'GSE315100_matrix.txt.gz');cols=[f'bal_{i:03}_count' for i in range(1,17)];rows=old_meta('GSE315100');gsm={f"{r['title'].replace(' | ','').strip()}_count":r['sample'] for r in rows}
    add('GSE315100',d,cols,['Ctrl']*4+['BIBR']*4+['CtrlEt']*4+['BIBREt']*4,'counts','human','iPSC_microglia',ids=d.gene_id,gsm=gsm,limit='single male WTC11 line; only 16 of 28 deposited entries represented')
    a=wb('GSE275834_BLEODMY_vs_BLEO.xls');b=wb('GSE275834_BLEO_vs_CTRL.deseq2.annotation_HMC3_.xls')
    invalid=[]
    for label,frame in [('BLEODMY_vs_BLEO',a),('BLEO_vs_CTRL',b)]:
        keep=frame['Gene id'].astype(str).str.fullmatch(r'ENSG\d+(?:\.\d+)?')
        count_cols=[col for col in frame if col.endswith('_count')]
        if not frame.loc[~keep,count_cols].isna().all().all():raise ValueError('non-gene workbook rows contain counts and require review')
        invalid.extend(dict(table=label,source_excel_row=int(i)+2,first_cell=str(frame.loc[i,'Gene id'])[:200],reason='non-Ensembl annotation continuation row; all sample counts missing') for i in frame.index[~keep])
        if label=='BLEODMY_vs_BLEO':a=frame.loc[keep].reset_index(drop=True)
        else:b=frame.loc[keep].reset_index(drop=True)
    write(pd.DataFrame(invalid),OUT/'workbook_annotation_continuations.tsv')
    assert list(a['Gene id'])==list(b['Gene id'])
    for col in [f'BLEO0{i}_count' for i in range(1,4)]:assert np.array_equal(a[col].to_numpy(),b[col].to_numpy()),'shared BLEO count mismatch'
    cols=[f'{g}0{i}_tpm' for g in ['CTRL','BLEO','BLEODMY'] for i in range(1,4)]
    d=pd.concat([b[['Gene id','Gene name']+cols[:6]],a[cols[6:]]],axis=1)
    add('GSE275834',d,cols,['CTRL']*3+['BLEO']*3+['BLEODMY']*3,'abundance','human','HMC3',genes=d['Gene name'],ids=d['Gene id'],limit='published count columns are fractional, count provenance unresolved: TPM-trend exploration only; external DMY reference, not formula efficacy')
    d=read(OLD/'GSE290755_RNA_all_rpkm.txt.gz');cols=list(d.columns[1:]);add('GSE290755',d,cols,['KD']*3+['NTC']*2,'abundance','human','iPSC_microglia',genes=d.gene_name,limit='NTC third replicate absent; RPKM-only; no same-system healthy AD reference')
    for cell in ['HMC3','SVGA']:
        d=pd.read_csv(OLD/('HMC3_Microglia_TPM_values.csv' if cell=='HMC3' else 'SVGA_Astrocytes_TPM_values.csv'));cols=list(d.columns[1:])
        add('E_MTAB_15005_'+cell,d,cols,['BrdU']*3+['DMSO']*3,'abundance','human',cell,limit='separate cultures, not communication experiment; one cell line')
    d=wb('GSE228927_FPKM_anno.xls');cols=[g+str(i) for g in ['Hemin','LPS','NOT'] for i in range(1,4)]
    add('GSE228927',d,cols,['Hemin']*3+['LPS']*3+['NOT']*3,'abundance','mouse','primary_astrocytes',genes=d.id,limit='MCM receiver response; residual LPS/hemin and stroke context remain')
    d=wb('GSE175485_analysis.xlsx',2,header=1);cols=list(d.columns[27:39]);assert len(cols)==12
    add('GSE175485',d,cols,['ECM_PBS']*3+['ECM_LPS']*3+['MCM_PBS']*3+['MCM_LPS']*3,'counts','mouse','primary_astrocytes',genes=d.GeneSymbol,ids=d.GeneSymbol,limit='conditioned media from endothelial/microglia; sender-response and carryover require separate functional evidence')
    for acc in ['GSE345471','GSE345470']:
        d=wb(acc+'_1_genes_fpkm_expression.xlsx');cols=[c for c in d if c.startswith('count.')]
        groups=[re.sub(r'_[123]$','',c.replace('count.','')) for c in cols]
        add(acc,d,cols,groups,'counts','mouse','microglia' if acc=='GSE345471' else 'MG6_MG7',genes=d.gene_name,ids=d.gene_id,limit='aging, not AD; vehicle annotation conflicts resolved by title' if acc=='GSE345471' else 'Pal is palbociclib per GSM, not palmitate; cell line genotype labels incomplete')
    d=read(RAW/'GSE286095_all.gene_counts.tsv.gz');cols=[c for c in d if c.startswith('sample.')]
    add('GSE286095',d,cols,['KO' if '.ko' in c else 'WT' for c in cols],'counts','mouse','microglia',genes=d.external_gene_name,ids=d.ensembl_gene_id,limit='female 20mo mice, Atg7 conditional genotype; aging cohort not amyloid cohort')
    d=read(RAW/'GSE275256_summary_star_genes_readcount.stranded.annotated.tsv.gz');d=d[~d.iloc[:,0].astype(str).str.startswith('N_')];cols=list(d.columns[1:10]);groups=['SLO_half']*3+['SLO_full']*3+['VEH']*3
    add('GSE275256',d,cols,groups,'counts','human','HMC3',genes=d.Name,ids=d.iloc[:,0],limit='chemical senescence benchmark, no treatment rescue')
    s=gzip.open(RAW/'GSE18740_series_matrix.txt.gz','rt').read();t=s.split('!series_matrix_table_begin\n')[1].split('!series_matrix_table_end')[0];d=read(io.StringIO(t))
    sf=gzip.open(RAW/'GSE18740_family.soft.gz','rt',encoding='utf-8').read();platform=read(io.StringIO(sf.split('!platform_table_begin\n')[1].split('!platform_table_end')[0]),dtype=str).fillna('')
    probes={r['ID']:r['Gene Symbol'] if ' /// ' not in r['Gene Symbol'] else '' for r in platform.to_dict('records')};cols=list(d.columns[1:]);assert (d[cols].to_numpy()>0).all()
    # Deposited RMA signals are positive linear scale; explicitly log2 once.
    d[cols]=np.log2(d[cols]);add('GSE18740',d,cols,['Control']*3+['Luteolin']*3+['LPS']*3+['LPS_Luteolin']*3,'microarray_log2','mouse','BV2',genes=[probes.get(x,'') for x in d.ID_REF],ids=d.ID_REF,limit='50uM luteolin, 50ng/ml LPS; deposited linear RMA signal log2 transformed; no formula identity/exposure claim')
    # Confounded tri-culture file stays descriptive with exact deposited batch metadata.
    d=read(OLD/'GSE325931_matrix.txt.gz');cols=list(d.columns[1:]);batch={r['description'].split('a')[-1]:r for r in old_meta('GSE325931')}
    groups=['iron_RSL3' if c in ['a07','a08','a11','a12'] else 'tri_vehicle' if c in ['a27','a34','a35','a01','a02','a03'] else 'mono_vehicle' if c in ['a24','a25','a26'] else 'co_vehicle' for c in cols]
    add('GSE325931',d,cols,groups,'abundance','human','iPSC_microglia_triculture',genes=d.gene,formal=False,limit='treatment completely confounded with batch; co a31 missing; descriptive only')
    meta=[]
    for p in sorted(RAW.glob('*_family.soft.gz'))+list(OLD.glob('*_family.soft.gz')):
        meta.extend([dict(series=p.name.split('_')[0],**r) for r in soft(p)])
    write(pd.DataFrame(meta),OUT/'official_sample_metadata.tsv');write(pd.DataFrame(samples),OUT/'analysis_samples.tsv');write(pd.DataFrame(qc),OUT/'matrix_qc.tsv');write(pd.DataFrame(mappings),OUT/'gene_mapping_audit.tsv.gz')
    # Frozen disease sets are loaded through the existing 024 helper (no main execution).
    spec=importlib.util.spec_from_file_location('f024',ROOT/'src/114_finalize_formula_targets_024.py');f=importlib.util.module_from_spec(spec);spec.loader.exec_module(f)
    sets,contexts=f.disease_sets();modules=[]
    g1=read(ROOT/f.cfg['inputs']['g1_244']);sign=dict(zip(g1.gene,np.sign(g1.meta_logFC)))
    for name,genes in sets.items():
        if name=='G1_significant413_exploratory':continue
        for g in genes:modules.append(dict(module=name,gene=g,weight=sign.get(g,1) if name=='G1_robust244' else 1,context=contexts[name]))
    write(pd.DataFrame(modules).drop_duplicates(['module','gene']),INT/'frozen_modules.tsv')
    (INT/'dataset_manifest.json').write_text(json.dumps(datasets,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print('prepared',len(datasets),'datasets',len(samples),'samples',flush=True)

if __name__=='__main__':main()

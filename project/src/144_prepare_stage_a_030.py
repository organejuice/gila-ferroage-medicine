"""Capture public resource snapshots and prepare frozen stage-A inputs."""
import gzip, hashlib, importlib.util, json, shutil, subprocess, zipfile
from pathlib import Path
import pandas as pd
import requests
ROOT=Path(__file__).resolve().parents[1]
C=json.loads((ROOT/'configs/stage_a_030.json').read_text())
RUN=C['run_id']; RAW=ROOT/'data/raw'/RUN; INT=ROOT/'data/interim'/RUN; OUT=ROOT/'results'/RUN
A=ROOT/'results/2026-10-03-fixed-formula-target-supplement-024-continuation'
B=ROOT/'results/2026-10-03-public-mechanism-analysis-026'
D=ROOT/'results/2026-10-03-disease-interface-perturbation-audit-027'
INPUTS={}
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(2**20),b''):h.update(b)
    return h.hexdigest()
def read(p,**kw):
    p=Path(p); INPUTS[str(p.relative_to(ROOT))]=sha(p)
    return pd.read_csv(p,sep='\t',low_memory=False,**kw)
def write(d,p):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    d.to_csv(p,sep='\t',index=False,compression={'method':'gzip','mtime':0} if str(p).endswith('.gz') else None)
def save(v,p):Path(p).write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
def fetch(name,url):
    p=RAW/name; p.parent.mkdir(parents=True,exist_ok=True)
    if not p.exists():
        temp=p.with_suffix(p.suffix+'.part')
        try:
            with requests.get(url,stream=True,timeout=(30,120)) as r:
                r.raise_for_status()
                with temp.open('wb') as f:
                    for chunk in r.iter_content(2**20):f.write(chunk)
            temp.rename(p)
        except Exception:
            if temp.exists():temp.unlink()
            raise
    return dict(file=str(p.relative_to(ROOT)),url=url,sha256=sha(p),bytes=p.stat().st_size)
def resources():
    urls={
      'STRING/9606.protein.links.full.v12.0.txt.gz':'https://stringdb-downloads.org/download/protein.links.full.v12.0/9606.protein.links.full.v12.0.txt.gz',
      'STRING/9606.protein.info.v12.0.txt.gz':'https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz',
      'Reactome/ReactomePathways.gmt.zip':'https://reactome.org/download/current/ReactomePathways.gmt.zip',
      'Reactome/version.txt':'https://reactome.org/ContentService/data/database/version',
      'KEGG/link_hsa_pathway.tsv':'https://rest.kegg.jp/link/hsa/pathway',
      'KEGG/list_pathway_hsa.tsv':'https://rest.kegg.jp/list/pathway/hsa',
      'KEGG/info.txt':'https://rest.kegg.jp/info/kegg'}
    records=[]
    for name,url in urls.items():
        print('resource',name,flush=True);records.append(fetch(name,url))
    commit_path=RAW/'NetMedPy_commit.json'
    if not commit_path.exists():save(requests.get('https://api.github.com/repos/menicgiulia/NetMedPy/commits/main',timeout=30).json(),commit_path)
    commit=json.loads(commit_path.read_text())['sha']
    for name in ['LICENSE.txt','netmedpy/netmedpy.py','README.md']:
        records.append(fetch('NetMedPy/'+name,'https://raw.githubusercontent.com/menicgiulia/NetMedPy/'+commit+'/'+name))
    save(records,RAW/'resource_manifest.json')
    with zipfile.ZipFile(RAW/'Reactome/ReactomePathways.gmt.zip') as z:
        names=[n for n in z.namelist() if n.endswith('.gmt')];assert len(names)==1
        (INT/'Reactome_source.gmt').write_bytes(z.read(names[0]))
    print('resources_complete',flush=True)
def inputs():
    with zipfile.ZipFile(RAW/'Reactome/ReactomePathways.gmt.zip') as z:
        names=[n for n in z.namelist() if n.endswith('.gmt')];assert len(names)==1
        (INT/'Reactome_source.gmt').write_bytes(z.read(names[0]))
    connections=pd.concat([read(A/'G1_244_layered_connections.tsv'),read(A/'exploratory_disease_connections.tsv')],ignore_index=True).fillna('')
    write(connections,INT/'connections.tsv.gz')
    for name in ['D0_interface_registry.tsv','D0_interface_genes.tsv','D0_TQ_Desired_reference.tsv','D0_functional_evidence_cards.tsv','P4_P5_qualification_matrix.tsv','P6_context_comparability.tsv','P7_shared_measurement_coverage.tsv']:
        write(read(D/name),INT/name)
    cfg=json.loads((ROOT/'configs/fixed_formula_target_supplement_024.json').read_text())
    meta=read(ROOT/cfg['inputs']['g1_full']); assert not meta.gene.duplicated().any()
    write(meta,INT/'G1_full_meta.tsv')
    write(read(ROOT/cfg['inputs']['g1_244']),INT/'G1_244.tsv')
    write(read(A/'G1_E1_released_target_summary.tsv'),INT/'G1_E1_released_target_summary.tsv')
    write(read(A/'G1_G5_gate_layer_summary.tsv'),INT/'G1_G5_gate_layer_summary.tsv')
    write(read(A/'compound_identity_audit.tsv'),INT/'compound_identity_audit.tsv')
    write(read(D/'P7_shared_gene_effects.tsv.gz'),INT/'P7_shared_gene_effects.tsv.gz')
    # All positive qualified targets, rather than the disease intersections, form drug sets for proximity.
    layers={t:set() for t in ['E1','E2','P','H']};records=[]
    p=A/'evidence_layered.tsv';INPUTS[str(p.relative_to(ROOT))]=sha(p)
    cols=['tier','target_symbol','candidate_identity_verified','action_direction','activity_type','assay_description','activity_outcome','compound_key']
    spec=importlib.util.spec_from_file_location('final024',ROOT/'src/114_finalize_formula_targets_024.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
    # positive() also needs value/relation/outcome fields: read once with all columns in bounded chunks.
    layers={t:set() for t in layers};records=[]
    for part in pd.read_csv(p,sep='\t',chunksize=30000,low_memory=False,dtype=str):
        for row in part.fillna('').to_dict('records'):
            tier=row.get('tier');gene=row.get('target_symbol')
            row['identity_exact']=row.get('identity_exact','').lower()=='true'
            row['candidate_identity_verified']=row.get('candidate_identity_verified','').lower()!='false'
            eligible=tier in layers and bool(gene) and row['identity_exact'] and row['candidate_identity_verified'] and (tier in ['P','H'] or m.positive(row))
            if tier in layers and gene:records.append({k:row.get(k,'') for k in cols}|dict(record_qualified=bool(eligible)))
            if eligible:layers[tier].add(gene)
    assert {t:len(v) for t,v in layers.items()}=={'E1':276,'E2':1244,'P':6461,'H':2081},'global target eligibility drift'
    write(pd.DataFrame([dict(tier=t,gene=g) for t in layers for g in sorted(layers[t])]),INT/'all_drug_targets.tsv')
    write(pd.DataFrame(records),INT/'target_action_records.tsv.gz')
    # Mechanism matrix preserves exact compound keys, source assay and record links, without asserting measured exposure.
    keys=['disease_set','tier','gene','compound_key']; rows=[]
    for key,x in connections.groupby(keys,sort=True):
        row=dict(zip(keys,key)); row.update(context=x.context.iloc[0],cids=';'.join(sorted(set(';'.join(x.cids.astype(str)).split(';'))- {''})),herbs=';'.join(sorted(set(';'.join(x.herbs.astype(str)).split(';'))-{''})),source_records=len(x),sources=';'.join(sorted(set(x.source))),assay_ids=';'.join(sorted(set(x.assay_id.astype(str)))),documents=';'.join(sorted(set(x.document.astype(str))-{''})),human_G1_source_released=bool((x.released_main_E1.astype(str).str.lower()=='true').any()),disease_effect=x.disease_effect.iloc[0],disease_FDR=x.disease_FDR.iloc[0],treatment_direction='UNDETERMINED',measured_formula_presence='NOT_ESTABLISHED',relevant_free_exposure='NOT_ESTABLISHED',direct_binding='NOT_ESTABLISHED_FROM_LAYER_ALONE',structure_available='exact_candidate_key; linked to 024 compound_identity_audit.tsv',functional_reference='D0_functional_evidence_cards.tsv; context-specific; not formula proof',source_connection_file='G1_244_layered_connections.tsv' if key[0]=='G1_robust244' else 'exploratory_disease_connections.tsv')
        rows.append(row)
    write(pd.DataFrame(rows),OUT/'mechanism_compound_target_matrix.tsv.gz')
    # Annotation of nominated candidates; do not convert disease RNA direction to an intervention objective.
    nominated=read(A/'G1_E1_released_target_summary.tsv').copy()
    nominated['direction_review']=nominated.gene.map({'CA2':'RNA down; assay inhibition; direction tension, therapeutic objective unresolved','SNCA':'aggregation modulation; direct binding not established; no stable-pocket assumption','AKR1B1':'RNA up; enzyme inhibition record; functional therapeutic objective unresolved','ALOX15B':'RNA up; enzyme inhibition record; functional therapeutic objective unresolved','LTA4H':'RNA up; enzyme inhibition record; functional therapeutic objective unresolved'})
    write(nominated,OUT/'G1_candidate_direction_review.tsv')
    # Use source-measured backgrounds if reconstructible; otherwise no invented all-human universe.
    bg=[]
    for g in meta.gene.astype(str):bg.append(dict(gate='G1',gene=g,background_scope='G1 meta-tested genes'))
    # G3/G4 source is original D2 gene TPM; G5 uses frozen original mouse filter + one-to-one mapping.
    d2=read(ROOT/'results/archive/2026-09-20-pre-v3-freeze/results/tables/D2_gene_tpm.tsv')
    print('D2 columns',list(d2.columns[:4]),flush=True)
    gene_col=next((c for c in ['gene','gene_symbol','symbol'] if c in d2.columns),d2.columns[0])
    for gate in ['G3','G4']:
        for g in sorted(set(d2[gene_col].dropna().astype(str))):bg.append(dict(gate=gate,gene=g,background_scope='original D2 quantified genes; not reconstructed expression-filter universe'))
    filt=read(ROOT/'results/2026-09-22-disease-side-v3.4-g5/tables/GSE261698_filter_genes.tsv')
    mapping=read(ROOT/'results/2026-09-22-disease-side-v3.4-g5/registries/GSE261698_gene_mapping.tsv')
    write(filt,INT/'G5_filter.tsv');write(mapping,INT/'G5_mapping.tsv')
    # Reuse 054's reciprocal NCBI GeneID orthology rule; expand to the actual filtered G5 background.
    ort=ROOT/'data/raw/NCBI_gene_orthologs/gene_orthologs.gz';INPUTS[str(ort.relative_to(ROOT))]=sha(ort)
    pairs=set()
    for part in pd.read_csv(ort,sep='\t',dtype=str,chunksize=500000):
        part.columns=[c.lstrip('#').lower().replace(' ','_') for c in part.columns]
        x=part[((part.tax_id=='10090')&(part.other_tax_id=='9606'))|((part.tax_id=='9606')&(part.other_tax_id=='10090'))]
        for r in x.itertuples(index=False):
            if r.relationship=='Ortholog':pairs.add((r.geneid,r.other_geneid) if r.tax_id=='10090' else (r.other_geneid,r.geneid))
    from collections import defaultdict
    mh=defaultdict(set);hm=defaultdict(set)
    for m_id,h_id in pairs:mh[m_id].add(h_id);hm[h_id].add(m_id)
    human=read(ROOT/'data/raw/NCBI_gene_orthologs/Homo_sapiens.gene_info.gz',dtype=str)
    hs=dict(zip(human.GeneID,human.Symbol)); mapping['mouse_geneid']=mapping.raw_id.str.extract(r'_(\d+)$',expand=False)
    keep=set(filt.loc[filt.filtered_in.astype(str).str.lower()=='true','gene'])
    mapping['human_symbol']=mapping.mouse_geneid.map(lambda m:hs.get(next(iter(mh[m])), '') if len(mh[m])==1 and len(hm[next(iter(mh[m]))])==1 else '')
    mapping['filtered_in']=mapping.gene_symbol.isin(keep)
    write(mapping,INT/'G5_background_orthology_audit.tsv')
    for g in sorted(set(mapping.loc[mapping.filtered_in & (mapping.human_symbol!=''),'human_symbol'])):bg.append(dict(gate='G5',gene=g,background_scope='GSE261698 filtered mouse genes; reciprocal unique NCBI human orthologs'))
    # G2 RNA annotation is explicitly distinct from occupancy; use D2 quantified RNA as a sensitivity only.
    for g in sorted(set(d2[gene_col].dropna().astype(str))):bg.append(dict(gate='G2',gene=g,background_scope='D2 RNA-proxy sensitivity ONLY; mouse occupancy background unavailable'))
    write(pd.DataFrame(bg),INT/'disease_backgrounds.tsv')
    save(INPUTS,OUT/'frozen_input_hashes.json')
    print('inputs_complete',len(connections),len(rows),{k:len(v) for k,v in layers.items()},flush=True)
def main():
    for p in [RAW,INT,OUT]:p.mkdir(parents=True,exist_ok=True)
    import sys
    if '--inputs-only' not in sys.argv:resources()
    inputs()
if __name__=='__main__':main()

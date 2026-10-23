"""Offline curation of frozen evidence, not a new effect fit or efficacy ranking."""
import gzip,hashlib,importlib.util,itertools,json,re,shutil
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-disease-interface-perturbation-audit-027';OUT=ROOT/'results'/RUN;RAW=ROOT/'data/raw'/RUN
PRIOR=ROOT/'results/2026-10-03-public-mechanism-analysis-026';INT=ROOT/'data/interim/2026-10-03-public-mechanism-analysis-026'
REG=ROOT/'results/2026-09-24-disease-side-v3.5-g6-input-audit/registries'
spec=importlib.util.spec_from_file_location('p127',ROOT/'src/127_prepare_expression_026.py');p127=importlib.util.module_from_spec(spec);spec.loader.exec_module(p127)
INPUTS=set()
def read(p):
    p=Path(p);INPUTS.add(p);return pd.read_csv(p,sep='\t',low_memory=False)
def jread(p):
    p=Path(p);INPUTS.add(p);return json.loads(p.read_text(encoding='utf8'))
def write(d,name):
    p=OUT/name;p.parent.mkdir(parents=True,exist_ok=True)
    d.to_csv(p,sep='\t',index=False,compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)
def serial(d):return json.loads(d.to_json(orient='records'))
def jsave(value,name): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def role(dataset,contrast):
    if dataset in ['GSE290755','GSE286095']:return 'genetic_mechanism','PC-E','not a real compound contrast'
    if dataset in ['GSE228927','GSE175485','GSE235076']:return 'receiver_stimulus','NOT_APPLICABLE','receiver response; no therapeutic intervention'
    if dataset=='PXD079947':
        if contrast=='IronRSL3Fer-IronRSL3':return 'external_rescue_compound','PC-B','ferrostatin rescue in HMC3; outside human AD'
        if contrast=='DMSO-NoTreatment':return 'vehicle_control','NOT_APPLICABLE','vehicle effect is not drug efficacy'
        return 'ferroptosis_induction','NOT_APPLICABLE','iron/RSL3 induction, not candidate therapy'
    if dataset=='GSE275834' and 'DMY' in contrast:return 'external_compound','PC-B','DMY in chemically senescent HMC3; outside human AD'
    if dataset=='GSE345470':return 'external_compound','PC-B','palbociclib in murine microglial cell lines'
    if dataset=='GSE345471' and 'AOMK' in contrast:return 'external_compound','PC-B','AOMK in aged mouse microglia; not AD'
    if dataset=='GSE18740' and 'Luteolin' in contrast:return 'external_compound','PC-B','luteolin in BV2; LPS model, not AD'
    if dataset=='GSE345471':return 'aging_state_comparison','NOT_APPLICABLE','age contrast is not treatment'
    return 'senescence_or_inflammation_induction','NOT_APPLICABLE','induction/control; not therapeutic normalization'
def p5_pass(row):
    return all(row.get(k)=='PASS' for k in ['formula_identity','therapeutic_direction','matched_reference','exposure','viability_stress','context_qualification'])
def main():
    OUT.mkdir(parents=True,exist_ok=True)
    modules=read(INT/'frozen_modules.tsv');tq=read(PRIOR/'TQ_map_and_Desired_reference.tsv');g1=read(REG/'G1_microglia_signed_sender_244.tsv.gz')
    orth=read(REG/'G5_receiver_one_to_one_human_ortholog_registry.tsv.gz');original=read(REG/'G5_frozen_receiver_module_gene_registry.tsv.gz')
    originals={
        'G1':(REG/'G1_microglia_signed_sender_244.tsv.gz','human AD microglial meta-logFC; 5 cohorts; CI not stored in original'),
        'G2':(ROOT/'results/2026-09-20-disease-side-v3.2-g2/tables/g2b_target_set_h4k12la_effects.tsv','H4K12la ChIP log2 occupancy AD minus WT; mouse 5XFAD; M and whole-brain T kept separate; CI95 stored'),
        'G3':(ROOT/'results/2026-09-20-disease-side-v3.2-g3-model-scoped/tables/D2_model_temporal_evidence.tsv','human iPSC microglial Aβ vs h0 at 4/24/96h; RNA module mean difference; CI unavailable'),
        'G4':(ROOT/'results/2026-09-20-disease-side-v3.2-g4/tables/D2_G4_RNA_sample_level_module_effects.tsv','human iPSC microglial Aβ vs h0 at 4/24/96h; globally standardized RNA module score; CI unavailable'),
        'G5':(ROOT/'results/2026-09-22-disease-side-v3.4-g5/tables/GSE261698_primary_module_effects.tsv','mouse APP/PS1 astrocyte TRAP RNA; 24 vs 24 animals; adjusted mean voom-logCPM module difference; CI stored')}
    effectpack={};sources=[]
    for gate,(p,unit) in originals.items():
        d=read(p);dest=OUT/'D0_original_effects'/f'{gate}.tsv';dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(gzip.decompress(p.read_bytes()) if p.name.endswith('.gz') else p.read_bytes());effectpack[gate]={'source':p.relative_to(ROOT).as_posix(),'unit_and_scope':unit,'rows':serial(d)}
        sources.append(dict(gate=gate,source=p.relative_to(ROOT).as_posix(),rows=len(d),sha256=sha(p),unit_and_scope=unit))
    write(pd.DataFrame(sources),'D0_effect_source_registry.tsv');write(orth,'D0_G5_one_to_one_mapping.tsv');write(original,'D0_G5_original_mouse_registry.tsv');write(tq,'D0_TQ_Desired_reference.tsv');write(read(PRIOR/'TQ_functional_evidence_cards.tsv'),'D0_functional_evidence_cards.tsv')
    registry=[]
    for module,d in modules.groupby('module',sort=True):
        gate=module[:2];ref=tq[tq.endpoint==module].iloc[0].to_dict()
        selector=module[3:].replace('_ortholog_exploratory','').replace('_exploratory','') if gate in ['G3','G4','G5'] else 'all 244 signed genes' if gate=='G1' else 'target_set_id=M1_glycolysis_predeclared_v3_2; context=M separately from T'
        registry.append(dict(interface=module,gate=gate,genes=len(d),context=d.context.iloc[0],interface_quantity='RNA signed association' if gate=='G1' else 'RNA proxy; original occupancy is a distinct modality' if gate=='G2' else 'RNA module proxy',original_effects='D0_original_effects/'+gate+'.tsv',original_effect_selector=selector,original_membership_policy='effect refers to original measured module; 026 mapped human gene universe may differ, do not relabel original score as a refit',gene_direction_available=gate=='G1',individual_gene_direction_reason='G1 meta-logFC sign only; no therapeutic direction' if gate=='G1' else 'module-level disease direction must not be assigned to individual genes',activity_or_communication_direct_measurement=False,primary_therapeutic_eligible=False,normalization_eligible=False,desired_numeric_reference='UNAVAILABLE',TQ_class=ref['TQ_class'],therapeutic_objective=ref['Therapeutic_objective']))
    write(pd.DataFrame(registry),'D0_interface_registry.tsv');write(modules,'D0_interface_genes.tsv')
    sets={m:set(d.gene) for m,d in modules.groupby('module')};overlap=[]
    for a,b in itertools.combinations(sorted(sets),2):overlap.append(dict(interface_a=a,interface_b=b,shared=len(sets[a]&sets[b]),union=len(sets[a]|sets[b]),jaccard=len(sets[a]&sets[b])/len(sets[a]|sets[b]),independence_assumed=False))
    write(pd.DataFrame(overlap),'D0_overlap_audit.tsv')
    jsave(dict(version=RUN,status='FROZEN_FOR_CONTEXT_RESTRICTED_COMPARISON',interfaces=registry,genes=serial(modules),original_effects=effectpack,G5_mouse_gene_registry=serial(original),G5_orthology=serial(orth),TQ_and_Desired_reference=serial(tq),unsupported_quantities=['TF_activity','pathway_activity','direct_M_to_A_communication','numeric_functional_Desired_reference'],reopen_017=False,reopen_G6=False,claim_ceiling='C0'), 'D0_formal_interface_package.json')
    datasets=[]
    for name in ['dataset_manifest.json','receiver_dataset_manifest.json','protein_dataset_manifest.json','sample_sensitivity_manifest.json']:datasets+=jread(INT/name)
    ds={d['id']:d for d in datasets};summary=read(PRIOR/'all_gene_contrast_summary.tsv');samples=read(PRIOR/'analysis_samples.tsv');receiver=read(PRIOR/'receiver_analysis_samples.tsv');traces=read(PRIOR/'sample_traceability_audit.tsv')
    write(traces,'P4_sample_traceability.tsv');write(pd.concat([samples,receiver],ignore_index=True),'P4_sample_registry.tsv')
    # Preserve original sample protocols rather than infer missing exposure parameters.
    protocols=[]
    for p in sorted(list((ROOT/'data/raw/2026-10-03-public-mechanism-analysis-026').glob('*family.soft.gz'))+list((ROOT/'workspaces/2026-10-02-public-mechanism-bridge-search-025/raw/sources').glob('*family.soft.gz'))):
        INPUTS.add(p)
        for r in p127.soft(p):
            protocols.append(dict(dataset=p.name.split('_')[0],GSM=r['GSM'],title=r.get('title',''),characteristics=r.get('characteristics_ch1',''),treatment_protocol=r.get('treatment_protocol_ch1',''),source=p.relative_to(ROOT).as_posix()))
    write(pd.DataFrame(protocols),'P4_official_protocol_records.tsv')
    write(read(PRIOR/'official_sample_metadata.tsv'),'P4_inherited_official_sample_metadata.tsv')
    sdrf=ROOT/'workspaces/2026-10-02-public-mechanism-bridge-search-025/raw/sources/E-MTAB-15005.sdrf.txt'
    write(read(sdrf),'P4_BrdU_source_SDRF.tsv')
    jsave(jread(PRIOR/'protein_pre_effect_protocol.json'),'P4_protein_independence_protocol.json')
    qualification=[]
    for r in summary.to_dict('records'):
        d=ds[r['dataset']];purpose,pc,reason=role(r['dataset'],r['contrast']);row=dict(**r,intervention_role=purpose,PC=pc,PC_reference='current human AD M/A disease scope; model-restricted PC-B references only',PC_reason=reason,species=d['species'],cell=d['cell'],assay=d['kind'],source_limitation=d['limitation'],statistical_use='SENSITIVITY_ONLY' if r['analysis_family']=='sample_sensitivity' else 'PROVISIONAL_SAMPLE_MAPPING' if r['dataset']=='E_MTAB_15005_SVGA' else 'CONTEXT_RESTRICTED_INFERENCE',formula_identity='UNRESOLVED_NOT_ASSIGNED_TO_FORMULA',therapeutic_direction='NOT_QUALIFIED',matched_reference='NO_QUALIFIED_NUMERIC_DESIRED_REFERENCE',exposure='UNKNOWN_FREE_RELEVANT_EXPOSURE',viability_stress='NOT_QUALIFIED_IN_CURRENT_ANALYSIS',context_qualification='LIMITED_OR_UNMATCHED',compound_condition_qualified=False,formula_ranking_eligible=False,normalization_eligible=False)
        assert not p5_pass(row);qualification.append(row)
    # A descriptive batch-confounded dataset and a control-only dataset must also appear.
    for name in ['GSE325931','GSE246888']:
        d=ds[name];descriptive=read(PRIOR/'all_module_effects.tsv');contrasts=sorted(descriptive.loc[descriptive.dataset==name,'contrast'].unique()) if name=='GSE325931' else ['NO_INFERENTIAL_CONTRAST']
        for contrast in contrasts:
            qualification.append(dict(dataset=name,contrast=contrast,analysis_family='descriptive_or_QC',intervention_role='induction_or_culture_comparison' if name=='GSE325931' else 'control_only',PC='NOT_APPLICABLE',species=d['species'],cell=d['cell'],assay=d['kind'],source_limitation=d['limitation'],statistical_use='DESCRIPTIVE_ONLY' if name=='GSE325931' else 'QC_ONLY',formula_identity='UNRESOLVED_NOT_ASSIGNED_TO_FORMULA',therapeutic_direction='NOT_QUALIFIED',matched_reference='NO_QUALIFIED_NUMERIC_DESIRED_REFERENCE',exposure='UNKNOWN_FREE_RELEVANT_EXPOSURE',viability_stress='NOT_QUALIFIED_IN_CURRENT_ANALYSIS',context_qualification='LIMITED_OR_UNMATCHED',compound_condition_qualified=False,formula_ranking_eligible=False,normalization_eligible=False))
    for row in qualification:
        name=row['dataset'];contrast=row['contrast'];purpose=row['intervention_role']
        row['compound_name_observed']='ferrostatin rescue agent' if name=='PXD079947' and purpose=='external_rescue_compound' else 'DMY' if name=='GSE275834' and 'DMY' in contrast else 'AOMK' if name=='GSE345471' and 'AOMK' in contrast else 'palbociclib' if name=='GSE345470' else 'luteolin' if name=='GSE18740' and 'Luteolin' in contrast else 'NOT_A_THERAPEUTIC_COMPOUND_CONTRAST'
        row['compound_CID_identity_audit']='NOT_PERFORMED_NO_FORMULA_ASSIGNMENT';row['potency_vs_relevant_exposure']='UNQUALIFIED';row['nominal_regimen_source']='P4_official_protocol_records.tsv; P4_inherited_official_sample_metadata.tsv; P4_BrdU_source_SDRF.tsv; protein source has unqualified free exposure'
        row['stress_viability_interpretation']='No omnibus matched viability/stress criterion frozen for these expression contrasts; not inferred from RNA effect magnitude'
    write(pd.DataFrame(qualification),'P4_P5_qualification_matrix.tsv')
    contexts=[]
    for d in datasets:
        for gate in ['G1','G2','G3','G4','G5']:
            reference_species='mouse' if gate in ['G2','G5'] else 'human';targetcell='astrocyte' if gate=='G5' else 'microglia'
            cell=d['cell'];has_target=('astro' in cell or cell=='SVGA') if targetcell=='astrocyte' else ('microglia' in cell or cell in ['HMC3','BV2','MG6_MG7'])
            contexts.append(dict(dataset=d['id'],gate=gate,species='matched' if d['species']==reference_species else 'unmatched',species_reason=f"actual {d['species']} vs {reference_species}; orthology is not species equivalence",cell='partial' if has_target else 'unmatched',cell_reason=f'actual {cell} vs {targetcell}; line/state differences retained',disease_model='unmatched',model_reason='no dataset reproduces the frozen human AD / 5XFAD occupancy / Aβ time-course / APPPS1 astrocyte comparison',stage_pathology='unknown',stage_reason='no qualified matched pathology stage',brain_region='unknown',region_reason='no qualified region mapping for the comparison',dose_exposure='unknown',dose_reason='nominal protocols retained separately; relevant free exposure not qualified',duration='unknown',duration_reason='recorded durations do not establish correspondence to frozen disease history',assay='limited' if gate=='G2' or d['kind']=='protein_log2' else 'compatible',assay_reason='RNA/protein cannot directly measure site-specific occupancy; protein is a distinct modality' if gate=='G2' or d['kind']=='protein_log2' else 'gene expression quantity can be compared descriptively; score scales differ',overall_score='NOT_COMPUTED',therapeutic_transport_qualified=False))
    write(pd.DataFrame(contexts),'P6_context_comparability.tsv')
    h=read(ROOT/'data/raw/2026-09-30-fixed-formula-target-supplement-024/HGNC/hgnc_complete_set.tsv');approved=set(h.loc[h.status=='Approved','symbol']);coverage=[];shared=[]
    union=set(modules.gene)
    for folder in [PRIOR/'gene_effects',PRIOR/'receiver/gene_effects',PRIOR/'protein/gene_effects',PRIOR/'sample_sensitivity/gene_effects']:
        for p in sorted(folder.glob('*.tsv.gz')):
            d=read(p);mapped=d[d.human_symbol.isin(approved)];genes=set(mapped.human_symbol);dataset=d.dataset.iloc[0];contrast=d.contrast.iloc[0]
            for module,gs in sets.items():
                common=gs&genes;coverage.append(dict(dataset=dataset,contrast=contrast,interface=module,measured_approved_genes=len(genes),shared_genes=len(common),interface_genes=len(gs),coverage=len(common)/len(gs),shared_gene_list=';'.join(sorted(common)),effect_quantity='protein log2 difference' if dataset=='PXD079947' else 'source-specific RNA effect; see 026 normalization',disease_source_quantity='H4K12la occupancy; RNA is proxy only' if module.startswith('G2') else 'disease RNA association; model or human scopes retained',primary_therapeutic_eligible=False,coverage_adequacy='NOT_PREQUALIFIED_NO_NEW_THRESHOLD',uncertainty_source=p.relative_to(ROOT).as_posix()))
            common=mapped[mapped.human_symbol.isin(union)].copy();common['source']=p.relative_to(ROOT).as_posix();common['feature_to_gene_aggregation']='NOT_COLLAPSED_MULTIPLE_FEATURES_RETAINED';shared.append(common)
    for name in ['GSE325931','GSE246888']:
        matrix=read(INT/(name+'_expression.tsv.gz'));mapping=read(INT/(name+'_gene_map.tsv'));varying=set(matrix.loc[matrix.iloc[:,1:].std(axis=1)>0,'feature_id']);genes=set(mapping.loc[mapping.feature_id.isin(varying)&mapping.human_symbol.isin(approved),'human_symbol'])
        contrasts=sorted({r['contrast'] for r in qualification if r['dataset']==name})
        for contrast in contrasts:
            for module,gs in sets.items():
                common=gs&genes;coverage.append(dict(dataset=name,contrast=contrast,interface=module,measured_approved_genes=len(genes),shared_genes=len(common),interface_genes=len(gs),coverage=len(common)/len(gs),shared_gene_list=';'.join(sorted(common)),effect_quantity='DESCRIPTIVE_OR_QC_ONLY_NO_GENE_INFERENCE',disease_source_quantity='H4K12la occupancy; RNA proxy only' if module.startswith('G2') else 'disease RNA association; scope retained',primary_therapeutic_eligible=False,coverage_adequacy='NOT_PREQUALIFIED_NO_NEW_THRESHOLD',uncertainty_source='NO_GENE_EFFECT_OR_CI_AVAILABLE'))
    write(pd.DataFrame(coverage),'P7_shared_measurement_coverage.tsv');write(pd.concat(shared,ignore_index=True),'P7_shared_gene_effects.tsv.gz')
    write(read(PRIOR/'all_module_effects.tsv'),'P7_source_module_effects_and_scenarios.tsv')
    for name in ['module_leave_one_out.tsv','sample_influence_summary.tsv','pseudocount_sensitivity_summary.tsv','full_module_coverage.tsv']:
        write(read(PRIOR/name),'P7_'+name)
    inputs=[dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(INPUTS)]
    jsave(inputs,'historical_input_hash_manifest.json')
    print('D0',len(registry),'P4/P5',len(qualification),'P6',len(contexts),'P7',len(coverage),flush=True)
if __name__=='__main__':main()

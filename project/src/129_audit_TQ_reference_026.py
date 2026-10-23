"""Freeze D1/D2 scope from original functional evidence, not drug scores."""
import hashlib,json
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1];RUN='2026-10-03-public-mechanism-analysis-026';OUT=ROOT/'results'/RUN;INT=ROOT/'data/interim'/RUN
CARDS=[
 ('Pan2022','Positive feedback regulation of microglial glucose metabolism by histone H4 lysine 12 lactylation in Alzheimer\'s disease','10.1016/j.cmet.2022.02.013','https://pubmed.ncbi.nlm.nih.gov/35303422/','5XFAD mouse microglia; pharmacologic PKM2 interruption and microglial Pkm2 genetic ablation','H4K12la/glycolysis/PKM2 feedback contributes to dysfunction in the reported model','publication functional evidence; not all glycolytic genes or human AD intervention qualification'),
 ('Srr2026','Microglial serine racemase knockout alleviates Alzheimer-like neuropathology and behavioral deficit via lactylation-mediated anti-inflammation','10.1038/s42003-026-09772-y','https://www.nature.com/articles/s42003-026-09772-y','5XFAD Lyz2-Cre Srr deletion; H3K18la','reported increased H3K18la with anti-inflammatory and behavioral effects','H3K18la differs from H4K12la; Lyz2 specificity and sex/context restrictions'),
 ('Bussian2018','Clearance of senescent glial cells prevents tau-dependent pathology and cognitive decline','10.1038/s41586-018-0543-y','https://pmc.ncbi.nlm.nih.gov/articles/PMC6206507/','PS19 tauopathy mice; p16-positive cell clearance','senescent glial clearance improved specified tau/neurodegeneration outcomes','cell burden/clearance endpoint; does not justify reversing every SASP RNA or human AD signature'),
 ('Schafer2026','ARID5A RNA-binding coordinates microglial defense and ferroptosis in iPSC-derived models','10.1038/s41467-026-76131-0','https://www.nature.com/articles/s41467-026-76131-0','human iPSC microglia CRISPRi, neuronal co-culture and iron/RSL3 functional assays','ARID5A depletion affected lysosomal/iron/cytokine pathways and reduced co-cultured neuronal ferroptosis','single-model functional evidence; gene expression matrix 3 KD/2 NTC not full functional source data'),
 ('Yang2014','Regulation of ferroptotic cancer cell death by GPX4','10.1016/j.cell.2013.12.010','https://pmc.ncbi.nlm.nih.gov/articles/PMC4076414/','functional GPX4 inhibition/overexpression in ferroptosis models','GPX4 function opposes lipid peroxide-associated ferroptotic death','general functional protective evidence outside current microglia model; no RNA-to-enzyme equivalence'),
 ('BrainSenescence2025','Systematic profiling reveals distinct senescence signatures and regulators across human brain cell types','10.1038/s41467-025-66012-3','https://www.nature.com/articles/s41467-025-66012-3','BrdU treatment in separate human brain cell lines; additional functional experiments','senescence responses differ by cell type','BrdU induction not therapy; separate cell lines not communication or human AD rescue')
]

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    cards=pd.DataFrame(CARDS,columns=['evidence_id','title','doi','url','actual_context','finding_used','limits'])
    cards['verified_date']='2026-10-03';cards.to_csv(OUT/'TQ_functional_evidence_cards.tsv',sep='\t',index=False)
    modules=pd.read_csv(INT/'frozen_modules.tsv',sep='\t');rows=[]
    for mod,sub in modules.groupby('module',sort=True):
        mixed=mod.startswith('G2_') or (mod.startswith(('G3_','G4_')) and ('inflammation' in mod or 'SASP' in mod or 'Iron' in mod or 'Ferritinophagy' in mod or 'LPO' in mod))
        tq='TQ-B' if mixed else 'TQ-D'
        reason='mixed/context-dependent protective and maladaptive functions; RNA score alone cannot specify optimum' if mixed else 'current interface is disease/model RNA association or proxy; no endpoint-specific therapeutic direction established'
        refs='Pan2022;Srr2026' if mod.startswith('G2_') else 'Bussian2018;BrainSenescence2025' if mod.startswith('G3_') else 'Schafer2026;Yang2014' if mod.startswith('G4_') else 'frozen_G1_G5_registry'
        rows.append(dict(endpoint=mod,measured_quantity='signed RNA expression reference' if mod=='G1_robust244' else 'RNA module proxy',
            Disease_state=sub.context.iloc[0],TQ_class=tq,TQ_reason=reason,Therapeutic_objective='UNDETERMINED; no automatic global reversal',
            Desired_reference='NOT_FROZEN',stage_scope=sub.context.iloc[0],Primary_metric='descriptive frozen module/signature score only',
            Confidence='limited_model_or_association',Freeze_status='AUDIT_FROZEN_DIRECTION_UNRESOLVED',functional_evidence_ids=refs,
            primary_therapeutic_eligible=False,normalization_eligible=False,missing='endpoint-specific functional objective and qualified same-system disease/drug/desired reference'))
    extra=[
      ('H4K12la_PKM2_pathological_feedback','TQ-A','Pan2022','5XFAD microglia','reduce the pathological feedback in the reported model','matched WT and functional phenotype reference','site-specific H4K12la plus glycolytic/PKM2 activity and functional outcome'),
      ('senescent_glial_cell_burden','TQ-A','Bussian2018','PS19 tauopathy mice, p16-positive glial clearance','reduce pathology-associated senescent burden in this intervention context','matched nonpathological functional outcome reference','senescent cell burden plus tau/neurodegeneration outcome'),
      ('neuronal_ferroptotic_death_in_ARID5A_coculture','TQ-A','Schafer2026','reported human iPSC microglia-neuron co-culture with iron/RSL3','reduce neuronal ferroptotic death while preserve viable function','matched vehicle and ferroptosis-rescue control from same assay','direct neuronal death/rescue readout; not a transcript module'),
      ('GPX4_enzymatic_antiperoxidation','TQ-C','Yang2014','general ferroptosis cell models outside current human AD system','preserve protective GPX4 function; reduction is a potential risk','functional GPX4 activity and viable untreated reference','enzyme activity and lipid peroxide/death readout')]
    for endpoint,tq,ev,scope,obj,ref,metric in extra:
        rows.append(dict(endpoint=endpoint,measured_quantity='direct functional endpoint; not present in current RNA interfaces',Disease_state=scope,TQ_class=tq,
            TQ_reason='original functional intervention supports the specified endpoint only in its original context',Therapeutic_objective=obj,
            Desired_reference=ref,stage_scope=scope,Primary_metric=metric,Confidence='publication_functional_context_restricted',
            Freeze_status='QUALITATIVE_MODEL_OBJECTIVE_FROZEN_NUMERIC_REFERENCE_UNAVAILABLE',functional_evidence_ids=ev,
            primary_therapeutic_eligible=False,normalization_eligible=False,missing='same-system numeric functional data and prespecified quantitative criterion; current formula intervention absent'))
    frame=pd.DataFrame(rows);frame.to_csv(OUT/'TQ_map_and_Desired_reference.tsv',sep='\t',index=False)
    g1=modules[modules.module=='G1_robust244'][['gene']].copy();g1['TQ_class']='TQ-D';g1['Desired_direction']='UNDETERMINED';g1['reason']='human AD RNA association is not gene-level therapeutic qualification'
    g1.to_csv(OUT/'G1_244_gene_direction_audit.tsv',sep='\t',index=False)
    state=dict(status='EVIDENCE_AUDIT_COMPLETED_WITH_EXPLICIT_UNRESOLVED_REFERENCES',endpoints=len(frame),TQ_counts=frame.TQ_class.value_counts().to_dict(),
        G1_genes_audited=len(g1),current_RNA_interfaces_primary_eligible=0,primary_therapeutic_eligible=0,normalization_eligible=0,
        numerical_Desired_references_fully_frozen=0,C0=True,reopen_017=False,ranking_examined_before_freeze=False)
    (OUT/'TQ_D1_D2_audit.json').write_text(json.dumps(state,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [OUT/'TQ_map_and_Desired_reference.tsv',OUT/'TQ_functional_evidence_cards.tsv',ROOT/'configs/public_mechanism_analysis_026.json',ROOT/'docs/protocol/2026-10-03-public-mechanism-analysis-026.md']}
    mp=OUT/'pre_effect_freeze_manifest.json'
    if mp.exists() and json.loads(mp.read_text())!=manifest:raise ValueError('TQ or configuration changed after initial freeze')
    if not mp.exists():mp.write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(state,ensure_ascii=False),flush=True)

if __name__=='__main__':main()

from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "results" / "tables"
SUMMARY = ROOT / "data" / "interim" / "D3_discovery_summary.json"

summary = json.loads(SUMMARY.read_text(encoding="utf-8"))
rows = []
for cell_type, filename in [("Microglia", "DiseaseSignature_microglia.tsv"), ("Astrocyte", "DiseaseSignature_astrocyte.tsv")]:
    df = pd.read_csv(TABLE / filename, sep="\t")
    fdr_sig = df[df["FDR"] < 0.05]
    strong = fdr_sig[df.loc[fdr_sig.index, "logFC"].abs() > 1]
    top = df.iloc[0]
    rows.append({
        "dataset": "D3_GSE160936",
        "cell_type": cell_type,
        "n_pseudobulk_groups": int(top["n_pseudobulk_groups"]),
        "n_AD_groups": int(top["n_AD_groups"]),
        "n_Control_groups": int(top["n_Control_groups"]),
        "n_ranked_genes": int(len(df)),
        "n_FDR_lt_0.05": int(len(fdr_sig)),
        "n_abs_logFC_gt_1_and_FDR_lt_0.05": int(len(strong)),
        "top_nominal_gene": str(top["gene"]),
        "top_nominal_logFC": float(top["logFC"]),
        "top_nominal_P": float(top["P.Value"]),
        "top_nominal_FDR": float(top["FDR"]),
        "status": "exploratory_ranked_signature_only" if len(fdr_sig) == 0 else "primary_signature_supported",
        "interpretation": "No gene passed BH-FDR 0.05; retain signed ranking for hypothesis generation." if len(fdr_sig) == 0 else "At least one gene passed BH-FDR 0.05; review direction and sensitivity before promotion."
    })
acceptance = pd.DataFrame(rows)
acceptance.to_csv(TABLE / "D3_signature_acceptance.tsv", sep="\t", index=False, encoding="utf-8-sig")

m7 = pd.read_csv(TABLE / "D3_module_effects.tsv", sep="\t")
completion = pd.DataFrame([
    ["D3_matrix_ingestion", "complete", 24, "All 24 Cell Ranger filtered matrices passed archive structure and dimension audit."],
    ["D3_metadata_audit", "complete", 24, "24 samples, 12 donors, balanced AD/Control and EC/SSC; age, sex and RIN available."],
    ["D3_cell_QC", "complete", int(summary["cells_after_qc"]), f"{summary['cells_after_qc']}/{summary['cells_before_qc']} cells retained ({summary['qc_retention']:.4f}); sample-wise empirical thresholds; no doublet caller available."],
    ["D3_celltype_annotation", "provisional_complete", int(sum(summary["cell_type_counts"].values())), "Marker-module annotation completed; Microglia and Astrocyte retained as primary types; labels require sensitivity review."],
    ["D3_pseudobulk", "complete", int(summary["pseudobulk_groups"]), "Donor x region x cell type; 24 groups per primary cell type."],
    ["D3_gene_DE", "exploratory_only", int(acceptance["n_FDR_lt_0.05"].sum()), "Limma-voom disease contrast AD-Control with region, age, sex, RIN and donor block; zero genes passed FDR 0.05 in both primary cell types."],
    ["D3_M7_module", "exploratory_only", int((m7["FDR"] < 0.05).sum()), "All 95 M7 genes mapped, but module score was not significant in Microglia or Astrocyte."],
    ["I1_signed_disease_interface", "partial_exploratory", 2, "Signed ranked gene tables exist for Microglia and Astrocyte; not promoted to accepted disease signature."],
    ["I2_regulatory_activity", "pending", 0, "Requires D1 peak-to-gene annotation and regulatory activity construction."],
    ["I3_seven_module_vector", "partial", 1, "Only frozen M7 reference was tested here; M1-M6 are not yet frozen or evaluated."],
    ["C1_communication", "pending", 0, "Conditional on reproducible disease-state association; not run."],
    ["D6_independent_validation", "pending", 0, "Metadata audit exists; expression validation is not yet complete in this run."],
], columns=["component", "status", "n_units_or_positive_tests", "interpretation"])
completion.to_csv(TABLE / "DiseaseSideCompletion.tsv", sep="\t", index=False, encoding="utf-8-sig")

mechanism = pd.DataFrame([
    ["Overall", "D3 disease-side mechanism", 1, "Exploratory disease-state association", "D3 is cross-sectional and no gene-level or M7 module test passed BH-FDR 0.05; no causal perturbation or rescue."],
    ["M1", "Glycolysis", 0, "Not assessed", "Gene set not frozen in the current registry."],
    ["M2", "H4K12la regulatory interface", 0, "Not assessed", "D1 peak-to-gene annotation is pending; RNA cannot directly measure H4K12la."],
    ["M3", "SASP/senescence", 0, "Not assessed", "No frozen M3 gene-set composite was evaluated in D3."],
    ["M4", "Cell-cell communication", 0, "Not assessed", "C1 gate not passed and ligand-receptor analysis was not run."],
    ["M5", "Iron handling", 0, "Not assessed", "No stable empirical composite was frozen."],
    ["M6", "ACSL4/lipid peroxidation", 0, "Not assessed", "No frozen M6 composite and no direct lipid-flux measurement."],
    ["M7", "Ferro-aging reference", 1, "Hypothesis-compatible but unsupported", "95/95 reference genes mapped; Microglia logFC=-0.0898, FDR=0.8265; Astrocyte logFC=-0.0234, FDR=0.8468."],
], columns=["module", "mechanism_axis", "evidence_level_0_to_5", "conclusion", "reason"])
mechanism.to_csv(TABLE / "MechanismEvidenceAssessment.tsv", sep="\t", index=False, encoding="utf-8-sig")
print("Wrote D3_signature_acceptance.tsv, DiseaseSideCompletion.tsv, MechanismEvidenceAssessment.tsv")
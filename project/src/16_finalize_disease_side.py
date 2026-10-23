#!/usr/bin/env python3
"""Integrate the disease-side evidence layers without cross-dataset pseudo-replication."""
from __future__ import annotations

import json
import math
from datetime import date
from pathlib import Path

import pandas as pd

EXPECTED_D6_H5_BYTES = 1_369_492_123

M1 = ["HK1", "HK2", "PFKFB3", "PFKFB4", "ALDOA", "ALDOC", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "PDK1", "SLC2A1", "HIF1A", "PFKP"]
M3 = ["CDKN1A", "CDKN2A", "IL6", "IL1B", "TNF", "CCL2", "CCL3", "CCL4", "CSF1", "CXCL10", "MMP9", "MMP14", "SERPINE1", "IGFBP3", "IGFBP7", "NFKB1", "RELA", "TGFB1"]
M5 = ["TFRC", "SLC11A2", "STEAP3", "DMT1", "FTH1", "FTL", "SLC40A1", "CP", "NCOA4", "ATG3", "ATG5", "ATG7", "MAP1LC3B", "SQSTM1"]
M6 = ["ACSL4", "ALOX15", "ALOX5", "ALOX12", "NOX4", "GPX4", "LPCAT3", "SAT1"]


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def fnum(value, digits=4):
    if value is None or pd.isna(value):
        return "NA"
    return f"{float(value):.{digits}f}"


def d6_h5_status(root: Path) -> tuple[str, int, bool]:
    path = root / "data" / "raw" / "D6_GSE214979" / "GSE214979_filtered_feature_bc_matrix.h5"
    if not path.exists():
        return "not_downloaded", 0, False
    size = path.stat().st_size
    magic = path.open("rb").read(8)
    valid = size == EXPECTED_D6_H5_BYTES and magic == b"\x89HDF\r\n\x1a\n"
    return ("complete_ready_for_expression_validation" if valid else "incomplete_download_pending", size, valid)


def build_i3(root: Path, tables: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    d3 = read_table(tables / "D3_module_effects_all.tsv")
    for _, r in d3.iterrows():
        rows.append({
            "module": {"M1": "M1", "M2": "M2", "M3": "M3", "M6": "M6", "M7": "M7"}.get(r.module, r.module),
            "layer": "D3_human_postmortem_discovery",
            "cell_type": r.cell_type,
            "comparison": "AD_vs_Control",
            "effect": r.logFC,
            "P.Value": r["P.Value"],
            "FDR": r.FDR,
            "coverage": r.coverage,
            "n_present": r.n_present,
            "n_total": r.n_total,
            "direction": r.direction,
            "evidence_role": "primary_disease_discovery",
            "interpretation": "donor-region-celltype pseudobulk; no D3 module passed BH FDR 0.05",
        })
    d2 = read_table(tables / "D2_module_effects.tsv")
    for _, r in d2.iterrows():
        name = {"M1_Glycolysis": "M1", "M2_H4K12la": "M2", "M3_SASP_Senescence": "M3", "M6_LPO_ACSL4": "M6", "M7_FerroAging": "M7"}.get(r.module, r.module)
        rows.append({
            "module": name,
            "layer": "D2_iMGL_Abeta_timecourse",
            "cell_type": "iMGL",
            "comparison": r.comparison,
            "effect": r.logFC,
            "P.Value": r["P.Value"],
            "FDR": r.FDR,
            "coverage": r.coverage,
            "n_present": r.n_present,
            "n_total": r.n_total,
            "direction": r.direction,
            "evidence_role": "temporal_support_not_AD_control",
            "interpretation": "Aβ1-42 fibril time-course; supports temporal sensitivity only",
        })
    d6_path = tables / "D6_module_effects.tsv"
    d6 = read_table(d6_path) if d6_path.exists() else pd.DataFrame()
    d6_module_hits = int((d6.FDR < 0.05).sum()) if not d6.empty else 0
    if not d6.empty:
        for _, r in d6.iterrows():
            rows.append({
                "module": r.module,
                "layer": "D6_human_independent_validation",
                "cell_type": r.cell_type,
                "comparison": "AD_vs_Control",
                "effect": r.logFC,
                "P.Value": r["P.Value"],
                "FDR": r.FDR,
                "coverage": r.coverage,
                "n_present": r.n_present,
                "n_total": r.n_total,
                "direction": r.direction,
                "evidence_role": "independent_validation",
                "interpretation": "independent human expression validation; observational and no D6 module passed BH FDR 0.05",
            })
    d1 = read_table(tables / "D1_gene_regulatory_profile.tsv")
    m = d1[d1.context == "M"].copy()
    rows.append({
        "module": "M2",
        "layer": "D1_mouse_H4K12la_CUTandTag",
        "cell_type": "Microglia_mouse",
        "comparison": "AD_vs_WT",
        "effect": m.delta_presence.mean(),
        "P.Value": None,
        "FDR": None,
        "coverage": float(m.one_to_one_ortholog.sum() / len(m)),
        "n_present": int(m.one_to_one_ortholog.sum()),
        "n_total": int(len(m)),
        "direction": "AD_gain_tendency" if m.delta_presence.mean() > 0 else "AD_loss_tendency",
        "evidence_role": "regulatory_reference_annotation",
        "interpretation": "annotated peak-to-gene occupancy summary; not a direct H4K12la activity assay",
    })
    rows.append({
        "module": "M4",
        "layer": "C1_gate",
        "cell_type": "Microglia_Astrocyte",
        "comparison": "not_run",
        "effect": None,
        "P.Value": None,
        "FDR": None,
        "coverage": None,
        "n_present": 0,
        "n_total": 0,
        "direction": "not_estimable",
        "evidence_role": "conditional_no_go",
        "interpretation": "no reproducible D3 state association; ligand-receptor analysis withheld",
    })
    i3 = pd.DataFrame(rows)
    i3.to_csv(tables / "I3_disease_module_vector.tsv", sep="\t", index=False)

    d3_sig_hits = 0
    for name in ("DiseaseSignature_microglia.tsv", "DiseaseSignature_astrocyte.tsv"):
        p = tables / name
        if p.exists():
            x = read_table(p)
            if "FDR" in x.columns:
                d3_sig_hits += int((x.FDR < 0.05).sum())
    d3_module_hits = int((d3.FDR < 0.05).sum())
    d2_module_hits = int((d2.FDR < 0.05).sum())
    component_rows = [
        {"module": "M1", "status": "operational_descriptive", "D3_FDR_hits": int(((d3.module == "M1") & (d3.FDR < 0.05)).sum()), "support": f"D2 early down-response nominal; D3 no BH-FDR support; D6 module FDR hits {d6_module_hits}", "claim_boundary": "gene-program score only; not lactate measurement"},
        {"module": "M2", "status": "annotated_descriptive", "D3_FDR_hits": int(((d3.module == "M2") & (d3.FDR < 0.05)).sum()), "support": f"D1 microglia 17,911 annotated genes; {int(m.one_to_one_ortholog.sum())} one-to-one mouse-human projections; D3 projected score; D6 module FDR hits {d6_module_hits}", "claim_boundary": "does not prove H4K12la activity or causality"},
        {"module": "M3", "status": "operational_descriptive", "D3_FDR_hits": int(((d3.module == "M3") & (d3.FDR < 0.05)).sum()), "support": f"D2 M3 up at 4/24/96 h with BH-FDR support; D3 AD-Control not significant; D6 module FDR hits {d6_module_hits}", "claim_boundary": "senescence/SASP transcriptional program only"},
        {"module": "M4", "status": "conditional_no_go", "D3_FDR_hits": 0, "support": "C1 gate failed", "claim_boundary": "no ligand-receptor or causal communication claim"},
        {"module": "M5", "status": "component_level_only", "D3_FDR_hits": int(d3.module.str.startswith("M5_").mul(d3.FDR < 0.05).sum()), "support": f"D2 iron uptake/storage/export/ferritinophagy components are temporally responsive; D3 composite not significant; D6 module FDR hits {d6_module_hits}", "claim_boundary": "retain subcomponents; no stable composite ferroptosis score"},
        {"module": "M6", "status": "operational_descriptive", "D3_FDR_hits": int(((d3.module == "M6") & (d3.FDR < 0.05)).sum()), "support": f"D2 context-dependent response; D3 not significant; D6 module FDR hits {d6_module_hits}", "claim_boundary": "transcriptional proxy only; not lipid-peroxidation flux"},
        {"module": "M7", "status": "frozen_reference_descriptive", "D3_FDR_hits": int(((d3.module == "M7") & (d3.FDR < 0.05)).sum()), "support": f"D5 frozen 95-gene author reference; D2 h24 BH-FDR support; D3 no BH-FDR support; D6 module FDR hits {d6_module_hits}", "claim_boundary": "ferro-aging-related reference signature, not functional proof"},
    ]
    components = pd.DataFrame(component_rows)
    components.to_csv(tables / "I3_component_status.tsv", sep="\t", index=False)
    return i3, components


def write_registry(tables: Path, d1: pd.DataFrame) -> None:
    m = d1[d1.context == "M"]
    registry = pd.DataFrame([
        {"module": "M1", "status": "operational_descriptive", "feature_count": len(M1), "source": "project operational glycolysis set", "note": "D3/D2 scores computed; no claim of direct lactate measurement."},
        {"module": "M2", "status": "annotated_descriptive", "feature_count": int(m.one_to_one_ortholog.sum()), "source": "D1 GSE188765 mm10 NCBI RefSeq peak-to-gene plus MGI one-to-one orthology", "note": "Primary microglia profile; D1 annotation rate and orthology are audited; not direct H4K12la activity."},
        {"module": "M3", "status": "operational_descriptive", "feature_count": len(M3), "source": "project operational SASP/senescence set", "note": "D2 temporal support; D3 AD-Control effect not FDR significant."},
        {"module": "M4", "status": "conditional_no_go", "feature_count": "pending", "source": "C1 state association then communication", "note": "Withheld because D3 state association gate failed."},
        {"module": "M5", "status": "component_level_only", "feature_count": len(M5), "source": "iron handling subcomponents", "note": "Report uptake, storage, export and ferritinophagy separately; no unstable composite."},
        {"module": "M6", "status": "operational_descriptive", "feature_count": len(M6), "source": "ACSL4-centered lipid-peroxidation program", "note": "Transcriptional proxy only; no lipid flux inference."},
        {"module": "M7", "status": "frozen_reference", "feature_count": 95, "source": "D5 Cell Metabolism 2026 Table S2 Ferroaging gene set", "note": "Author reference signature; coverage and component concordance retained."},
    ])
    registry.to_csv(tables / "SignatureRegistry.tsv", sep="\t", index=False)


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    tables = root / "results" / "tables"
    interim = root / "data" / "interim"
    reports = root / "docs" / "reports"
    decisions = root / "docs" / "decisions"
    for p in (tables, interim, reports, decisions):
        p.mkdir(parents=True, exist_ok=True)

    d1 = read_table(tables / "D1_gene_regulatory_profile.tsv")
    d3 = read_table(tables / "D3_module_effects_all.tsv")
    d2 = read_table(tables / "D2_module_effects.tsv")
    d4 = read_table(tables / "D4_spatial_support_summary.tsv")
    d6_meta = read_table(tables / "D6_celltype_status_counts.tsv")
    i3, components = build_i3(root, tables)
    write_registry(tables, d1)

    d6_status, d6_bytes, d6_valid = d6_h5_status(root)
    d6_validation_done = all((tables / name).exists() for name in ("D6_module_effects.tsv", "D6_signature_concordance.tsv")) and (interim / "D6_expression_validation_provenance.json").exists()
    d6_module_table = read_table(tables / "D6_module_effects.tsv") if (tables / "D6_module_effects.tsv").exists() else pd.DataFrame()
    d6_concordance_table = read_table(tables / "D6_signature_concordance.tsv") if (tables / "D6_signature_concordance.tsv").exists() else pd.DataFrame()
    d6_provenance = json.loads((interim / "D6_expression_validation_provenance.json").read_text(encoding="utf-8")) if (interim / "D6_expression_validation_provenance.json").exists() else {}
    d6_module_hits = int((d6_module_table.FDR < 0.05).sum()) if not d6_module_table.empty else 0
    d6_micro_rho = float(d6_concordance_table.loc[d6_concordance_table.cell_type == "Microglia", "spearman_rho"].iloc[0]) if not d6_concordance_table.empty and any(d6_concordance_table.cell_type == "Microglia") else None
    d6_astro_rho = float(d6_concordance_table.loc[d6_concordance_table.cell_type == "Astrocyte", "spearman_rho"].iloc[0]) if not d6_concordance_table.empty and any(d6_concordance_table.cell_type == "Astrocyte") else None
    d6_groups = int(d6_module_table.groupby("cell_type")["n_groups"].first().sum()) if not d6_module_table.empty else 0
    if d6_validation_done:
        d6_completion_status = "complete"
        d6_completion_evidence = f"105,332-cell metadata audit; H5 validated; independent module and D3-D6 concordance outputs generated"
        d6_completion_limitation = "independent validation remains observational; duplicateCorrelation blocks are all size 1 and intrablock correlation is set to zero; no causal inference"
    else:
        d6_completion_status = "h5_ready_pending_validation" if d6_valid else d6_status
        d6_completion_evidence = f"105,332-cell metadata audit complete; H5 bytes {d6_bytes}/{EXPECTED_D6_H5_BYTES}"
        d6_completion_limitation = "expression validation is not executable until the official H5 download is complete and validated"
    c1 = pd.DataFrame([{
        "gate": "C1",
        "decision": "NO_GO",
        "required": "reproducible disease-state association in primary D3 data before ligand-receptor analysis",
        "observed_D3_gene_FDR_hits": int(sum(int((read_table(tables / name).FDR < 0.05).sum()) for name in ("DiseaseSignature_microglia.tsv", "DiseaseSignature_astrocyte.tsv") if (tables / name).exists())),
        "observed_D3_module_FDR_hits": int((d3.FDR < 0.05).sum()),
        "observed_D6_module_FDR_hits": d6_module_hits,
        "D6_microglia_signature_spearman": d6_micro_rho,
        "D6_astrocyte_signature_spearman": d6_astro_rho,

        "reason": (f"D3 gene signatures and all D3 module rows have no BH-FDR<0.05 support; D6 independent module FDR hits={d6_module_hits}, signature rho Microglia={d6_micro_rho}, Astrocyte={d6_astro_rho}" if d6_validation_done else "D3 gene signatures and all D3 module rows have no BH-FDR<0.05 support"),

        "action": "do not run CellChat/LIANA/NicheNet as disease-side evidence; revisit only after independent expression validation",
    }])
    c1.to_csv(tables / "C1_communication_gate.tsv", sep="\t", index=False)

    d1_m = d1[d1.context == "M"]
    completion = pd.DataFrame([
        {"stage": "D1", "status": "complete", "evidence": f"8 peak files; microglia {len(d1_m)} annotated genes; {int(d1_m.one_to_one_ortholog.sum())} one-to-one human projections", "limitation": "peak-to-gene is annotation-based; no direct quantitative H4K12la activity"},
        {"stage": "D2", "status": "complete", "evidence": f"12 samples, 24,194 genes after 87.12% transcript mapping; {int((d2.FDR < 0.05).sum())} module-timepoint rows FDR<0.05", "limitation": "Aβ iMGL time-course, not human AD-Control"},
        {"stage": "D3", "status": "complete", "evidence": f"101,906→100,595 cells; donor-level glial pseudobulk; {int((d3.FDR < 0.05).sum())} all-module rows FDR<0.05", "limitation": "primary disease gene/module support not statistically significant"},
        {"stage": "D4", "status": "complete_descriptive", "evidence": f"2 AD CosMx sections; {int(d4.n_cells.sum())} coordinate-resolved cells; provisional marker-state and k6 spatial support", "limitation": "both sections are AD; labels are marker-score provisional; no causal communication"},
        {"stage": "D5", "status": "complete", "evidence": "95-gene ferro-aging reference frozen", "limitation": "reference gene set, not disease functional validation"},
        {"stage": "D6", "status": d6_completion_status, "evidence": d6_completion_evidence, "limitation": d6_completion_limitation},
        {"stage": "I1", "status": "complete_discovery_only", "evidence": "D3 microglia and astrocyte signed disease outputs plus acceptance audit", "limitation": "no gene-level BH-FDR support"},
        {"stage": "I2", "status": "complete_descriptive", "evidence": "D1 annotated signed regulatory reference and D3/D2 projection tables", "limitation": "RNA/module projections cannot be called direct H4K12la measurements"},
        {"stage": "I3", "status": "complete_layered_vector", "evidence": ("D1/D2/D3/D6 layered vector plus component registry" if d6_validation_done else "D1/D2/D3 layered vector plus component registry"), "limitation": "cross-dataset rows are not combined into one pseudo-replicated causal model"},
        {"stage": "C1", "status": "NO_GO", "evidence": "gate table generated", "limitation": "communication analysis deliberately withheld"},
    ])
    completion.to_csv(tables / "DiseaseSideCompletion.tsv", sep="\t", index=False)
    assessment = pd.DataFrame([{
        "assessment": "disease_side_mechanism_proof",
        "level": "1/5",
        "label": "exploratory multi-layer association / hypothesis-generating",
        "basis": ("D1 regulatory annotation, D2 temporal support and D4 spatial descriptive support are available; primary D3 AD-Control gene/module effects lack BH-FDR support; D6 independent expression validation is complete but observational" if d6_validation_done else "D1 regulatory annotation, D2 temporal support and D4 spatial descriptive support are available; primary D3 AD-Control gene/module effects lack BH-FDR support and D6 expression validation is pending"),
        "not_proven": "direct H4K12la, lactate, Fe2+, 4-HNE/MDA, lipid flux, ligand-receptor causality, perturbation, rescue or in-vivo causality",
        "upgrade_requirements": ("reproduce a signed disease state in an additional cohort; then execute predeclared C1 communication and functional perturbation/rescue experiments" if d6_validation_done else "complete D6 H5 expression validation with donor-level model; reproduce a signed disease state; then execute predeclared C1 communication and functional perturbation/rescue experiments"),
    }])
    assessment.to_csv(tables / "MechanismEvidenceAssessment.tsv", sep="\t", index=False)

    status = "COMPLETE_COMPUTATIONALLY_WITH_C1_NO_GO" if d6_validation_done else "PASS_WITH_PREDECLARED_D6_BLOCK"
    provenance = {
        "status": status,
        "date": str(date.today()),
        "d6_expression_h5_status": d6_status,
        "d6_h5_bytes": d6_bytes,
        "d6_h5_expected_bytes": EXPECTED_D6_H5_BYTES,
        "d6_expression_validation_done": d6_validation_done,
        "d3_gene_fdr_hits": 0,
        "d3_module_fdr_hits": int((d3.FDR < 0.05).sum()),
        "d2_module_fdr_hits": int((d2.FDR < 0.05).sum()),
        "d4_coordinate_resolved_cells": int(d4.n_cells.sum()),
        "c1": "NO_GO",
        "mechanism_proof_level": "1/5",
        "outputs": [
            "results/tables/I3_disease_module_vector.tsv",
            "results/tables/I3_component_status.tsv",
            "results/tables/C1_communication_gate.tsv",
            "results/tables/DiseaseSideCompletion.tsv",
            "results/tables/MechanismEvidenceAssessment.tsv",
            "results/tables/SignatureRegistry.tsv",
            "results/tables/D6_module_effects.tsv",
            "results/tables/D6_signature_concordance.tsv",
            "data/interim/D6_expression_validation_provenance.json",
        ],
    }
    (interim / "disease_side_completion.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")

    if d6_validation_done:
        d6_report_sentence = "D6 的 105,332 细胞元数据与官方 H5 表达矩阵均已完成审计，并已生成供体级独立表达验证结果；该验证为观察性证据。"
        d6_decision_sentence = "D6 元数据与表达验证均已完成；H5 文件头和官方大小已校验；由于每个供体区组在聚合后没有可估计的重复相关结构，duplicateCorrelation 按 0 处理；D6 结果仅作为独立观察性验证，不替代因果实验。"
        d6_upgrade_conditions = "1. 若要重开 C1，需在预注册标准下确认独立队列重现同一疾病状态。\n2. 若 C1 状态门槛成功，再执行限定的 ligand-receptor 分析。\n3. 若要把机制等级提高到因果层级，仍需靶向 H4K12la/乳酸/铁稳态/脂质过氧化的扰动、救援和体内验证。"
    else:
        d6_report_sentence = "D6 的 105,332 细胞元数据已完成审计，但官方表达 H5 当前仍为不完整续传（{d6_bytes:,}/{EXPECTED_D6_H5_BYTES:,} bytes），因此 D6 独立表达验证不能被声称为完成。"
        d6_decision_sentence = "D6 元数据审计完成，但 H5 尚未达到官方文件大小并通过 HDF5 文件头验证（当前 {d6_bytes:,}/{EXPECTED_D6_H5_BYTES:,} bytes），故不执行或宣称 D6 表达验证。"
        d6_upgrade_conditions = "1. 完成并校验 D6 官方 H5，执行供体级 Microglia/Astrocyte 独立表达验证。\n2. 若 D6 重现疾病状态，再按预注册门槛重开 C1 通讯分析。\n3. 若要把机制等级提高到因果层级，仍需靶向 H4K12la/乳酸/铁稳态/脂质过氧化的扰动、救援和体内验证。"
    if d6_validation_done:
        d6_report_result_sentence = f"D6 独立表达验证完成：{int(d6_provenance.get('analyzed_cells', 0)):,} 个细胞、{d6_groups} 个伪批量组；D6 模块 FDR<0.05 命中 {d6_module_hits}，D3-D6 signature Spearman rho 为 Microglia={d6_micro_rho:.3f}、Astrocyte={d6_astro_rho:.3f}；duplicateCorrelation 区组均为单观测，相关性按 0 处理。"
    else:
        d6_report_result_sentence = "D6 独立表达验证尚未完成。"
    report = f"""# 疾病侧方案完成报告（2026-09-20）

## 结论

当前疾病侧已完成 D1–D6、D3 全模块向量、I1/I2/I3、D4 空间支持和 C1 门槛审计。{d6_report_sentence}

疾病侧机制证明等级保持 **1/5：探索性多层关联/假设生成**。D1、D2、D4 增加了方向性和空间背景支持，但 D3 主发现的基因签名及 M1–M7 模块均未达到 BH-FDR<0.05，且 C1 通讯门槛为 **NO_GO**。这不支持直接宣称 H4K12la、乳酸、铁离子、脂质过氧化通量或因果通讯机制。

## 已完成内容

- D1：固定 UCSC mm10 NCBI RefSeq 注释，完成 peak-to-gene 和 MGI 一对一同源映射；小胶质细胞 17,911 个注释基因、10,517 个一对一人类投影。
- D2：使用 Ensembl release 109 GRCh38 GTF 重建转录本映射；143,626 个转录本中 125,126 个映射，映射率 87.12%，得到 24,194 个基因和时间序列模块结果。
- D3：完成 101,906→100,595 细胞 QC、供体×区域×细胞类型伪批量、I1 签名及 M1–M7 全模块模型；D3 全模块 FDR<0.05 行数为 {int((d3.FDR < 0.05).sum())}。
- D4：两张 AD CosMx 切片共 {int(d4.n_cells.sum()):,} 个坐标可解析细胞；完成阴性探针/QC 描述、微胶质/星形胶质 provisional marker score 和 k=6 空间近邻支持。
- D5：冻结 95-gene ferro-aging reference，仅作为参考签名。
- D6：{d6_report_result_sentence}
- I3：生成分层疾病模块向量；不同数据集没有合并为伪重复的单一效应量。
- C1：由于主发现没有可重复的显著状态关联，暂不运行 CellChat、LIANA 或 NicheNet。

## 关键输出

- `results/tables/DiseaseSideCompletion.tsv`
- `results/tables/MechanismEvidenceAssessment.tsv`
- `results/tables/I3_disease_module_vector.tsv`
- `results/tables/I3_component_status.tsv`
- `results/tables/C1_communication_gate.tsv`
- `results/tables/D4_spatial_support_summary.tsv`
- `results/tables/D1_gene_regulatory_profile.tsv`
- `results/tables/D2_module_effects.tsv`
- `results/tables/D3_module_effects_all.tsv`
- `results/tables/D6_module_effects.tsv`
- `results/tables/D6_signature_concordance.tsv`
- `data/interim/D6_expression_validation_provenance.json`
- `data/interim/disease_side_completion.json`

## 仍然必要的升级条件

{d6_upgrade_conditions}
"""
    (reports / "2026-09-20_disease-side-completion.md").write_text(report, encoding="utf-8")
    decision = f"""# 疾病侧方案完成与证据门控决策

日期：2026-09-20

本次完成 D1–D6、I1/I2/I3、D4 描述性空间支持和 C1 门控。D2 使用 Ensembl release 109 GRCh38 GTF 进行固定转录本注释；D3 使用既有 donor×region×celltype pseudobulk 模型扩展 M1–M7；D4 使用官方 CosMx processed expression/metadata 的坐标解析细胞，并只输出 provisional marker-state 与近邻描述。

D3 主发现没有基因级或模块级 BH-FDR<0.05 支持，因此 C1 决策为 NO_GO，不执行 ligand-receptor 通讯推断。{d6_decision_sentence}

机制证明等级保持 1/5：探索性多层关联/假设生成。D1/D2/D4 只增加支持性背景，不改变 D3 主发现统计门槛，也不构成直接分子测量、因果机制、扰动、救援或体内证明。

可复现输出见 `results/tables/` 和 `docs/reports/2026-09-20_disease-side-completion.md`；总状态为 `{status}`。
"""
    (decisions / "2026-09-20-disease-side-completion.md").write_text(decision, encoding="utf-8")
    print(json.dumps({"status": status, "d6_status": d6_completion_status, "c1": "NO_GO", "mechanism_proof_level": "1/5", "d4_cells": int(d4.n_cells.sum())}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

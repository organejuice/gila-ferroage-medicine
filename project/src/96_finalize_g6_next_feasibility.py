#!/usr/bin/env python
"""Compile the locked, diagnosis-blind G6 feasibility outputs into an auditable report."""
from __future__ import annotations

import csv
import hashlib
import sys
from collections import Counter
from pathlib import Path

import yaml


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fmt_pct(value: float) -> str:
    return f"{100 * value:.2f}%"


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 96_finalize_g6_next_feasibility.py <project_root>")
    root = Path(sys.argv[1]).resolve()
    cfg_path = root / "configs/g6_next_blind_feasibility_001.yaml"
    config = yaml.safe_load(cfg_path.read_text(encoding="utf-8"))
    if config["diagnosis_group_effects"] is not False:
        raise RuntimeError("The locked feasibility config must prohibit diagnosis effects")

    res = root / "results/2026-09-27-g6-next-feasibility"
    cell = res / "cellbin"
    numi = res / "numi"
    candidate = res / "candidates"
    out_dir = root / "docs/reports/g6-next-feasibility-summary"
    out_report = out_dir / "execution_report.md"
    out_manifest = out_dir / "source_sha256.tsv"
    out_report.parent.mkdir(parents=True, exist_ok=True)

    integrity = read_tsv(numi / "fit_integrity.tsv")
    wp1_rows = read_tsv(numi / "mechanical_gate_by_donor.tsv")
    wp1_summary = read_tsv(numi / "donor_estimator_comparison.tsv")
    if len(integrity) != 4 or any(r["baseline_reproduced"] != "TRUE" or r["nUMI_matches_mode_definition"] != "TRUE" for r in integrity):
        raise RuntimeError("WP1 base-fit or nUMI-integrity audit failed")
    if len(wp1_rows) != 8:
        raise RuntimeError("WP1 must contain both donors, nUMI modes, and estimators")
    wp1_gate = {(r["GSM"], r["nUMI_mode"], r["estimator"]): r["gate_pass"] for r in wp1_rows}
    wp1_status = (numi / "report.md").read_text(encoding="utf-8")
    if "Status: STOP_NO_FULL_COHORT_REFIT" not in wp1_status:
        raise RuntimeError("WP1 report status differs from the locked mechanical result")
    max_full_weight_diff = max(float(r["full_weight_max_abs_diff_vs_locked"]) for r in integrity)

    cell_qc = read_tsv(cell / "cellbin_nucleus_counts_qc.tsv")
    post_qc = read_tsv(cell / "cellbin_nucleus_qc_after_150_genes.tsv")
    identity = read_tsv(cell / "identity_marker_gate.tsv")
    g1_gate = read_tsv(cell / "G1_coverage_gate_by_donor.tsv")
    g5_gate = read_tsv(cell / "G5_coverage_gate_by_module.tsv")
    marker_avail = read_tsv(cell / "marker_availability_summary.tsv")
    cluster_rows = read_tsv(cell / "exploratory_cluster_marker_annotation.tsv")
    annotation_report = (cell / "exploratory_annotation_report.md").read_text(encoding="utf-8")
    neighbor_report = (cell / "neighborhood_report.md").read_text(encoding="utf-8")
    if len(cell_qc) != 12 or len({r["donor_id"] for r in cell_qc}) != 12:
        raise RuntimeError("CellBin donor QC does not contain 12 unique donors")
    for key in ("cell_count_sum_matches_metadata", "gene_count_sum_matches_catalog", "cell_ID_gene_ID_crosswalk_unique", "polygon_coordinate_match"):
        if any(r[key] != "TRUE" for r in cell_qc):
            raise RuntimeError(f"CellBin integrity field failed: {key}")
    if len(post_qc) != 12:
        raise RuntimeError("CellBin post-QC table must contain 12 donors")
    if "Status: IDENTITY_NOT_CALIBRATED" not in annotation_report:
        raise RuntimeError("WP2 identity status does not match the locked gate")
    if "NEIGHBORHOOD_NOT_ELIGIBLE_MASK_UNAVAILABLE" not in neighbor_report:
        raise RuntimeError("WP2 neighborhood status does not match the predeclared mask rule")
    total_gt150 = sum(int(r["n_cells_ge150_genes"]) for r in cell_qc)
    valid_ge150 = sum(int(r["n_cells_ge150_and_polygon_valid"]) for r in post_qc)
    total_umi_fraction = [float(r["fraction_tissue_UMI_in_CellBin"]) for r in cell_qc]
    marker_donors: dict[str, int] = {}
    for marker_set in ("MG_training", "MG_holdout_disjoint_from_G1_G5", "ASC_training", "ASC_holdout_disjoint_from_G1_G5"):
        z = [r for r in marker_avail if r["marker_set"] == marker_set]
        marker_donors[marker_set] = sum(int(r["genes_meeting_1pct"]) >= 2 for r in z)
    annotated = read_tsv(cell / "annotated_nuclei_metadata.tsv")
    if len(annotated) != valid_ge150:
        raise RuntimeError("Annotated nucleus count does not match post-QC valid-polygon count")
    label_counts = Counter(r["label"] for r in annotated)
    g1_eligible_donors = sum(int(r["n_G1_genes_at_1pct"]) >= 61 for r in g1_gate)
    g5_best = min(int(r["donors_meeting_module_coverage"]) for r in g5_gate)
    g5_worst_module = next(r["module"] for r in g5_gate if int(r["donors_meeting_module_coverage"]) == g5_best)
    id_map = {r["lineage"]: r for r in identity}

    candidates = read_tsv(candidate / "candidate_design_summary.tsv")
    candidate_status = {r["candidate"]: r["status"] for r in candidates}
    if any(status == "PRIMARY_CANDIDATE" for status in candidate_status.values()):
        raise RuntimeError("Unexpected primary candidate in metadata-only screen")

    report = [
        "# G6 下一轮盲态可行性试验执行报告",
        "",
        "**方案编号：** `G6-FP-2026-09-27-001`  ",
        "**总体状态：** `G6_NOT_EXECUTABLE`  ",
        "**效应盲态：** 是。WP1/WP2 未读取诊断标签；WP3 只用诊断 metadata 审查技术设计。未计算任何 AD–对照表达、状态或邻域效应。",
        "",
        "## 结论",
        "",
        "本轮完成了模型 nUMI 口径比较、GSE269906 CellBin 核级可行性试验和替代队列 metadata 筛选。三条路线都没有获得进入新编号确认性 G6 方案的资格：WP1 停止全队列重拟合；WP2 未恢复可验证的 MG 身份且邻域掩膜缺失；WP3 没有合格主队列。因此维持 `G6_NOT_EXECUTABLE`，不解盲、不放宽门槛、不计算疾病组效应。",
        "",
        "## WP1：RCTD nUMI 口径",
        "",
        f"旧全库口径基线在两名供体的条码/分类顺序均一致，最大 full-weight 差为 `{max_full_weight_diff:.6g}`；两种 nUMI 定义与各自输入列和一致。拟合特征子集的重新拟合已完成。",
        "",
        "| 供体 | 口径 | 估计量 | MG≥0.5 单元 | G1≥1% 基因 | 留出 MG rho | 全部试验门槛通过 |",
        "|---|---|---|---:|---:|---:|---|",
    ]
    for r in wp1_summary:
        gate = wp1_gate[(r["GSM"], r["nUMI_mode"], r["estimator"])]
        report.append(
            f"| {r['GSM']} | {r['nUMI_mode']} | {r['estimator']} | {r['n_MG_ge_0_5']} | "
            f"{r['n_G1_genes_detected_at_least_1pct']} | {float(r['MG_holdout_marker_rho']):.3f} | {gate} |"
        )
    report.extend([
        "",
        "子集 nUMI 相比全库 nUMI 只相差 GSM8330063 的 172 UMI、GSM8330064 的 284 UMI；两种口径的身份支持和 G1 覆盖几乎不变。双细胞映射下每名供体仅 1 个 MG≥0.5 单元，全类型归一权重为 0–1 个；G1 仅 0–23 个基因达到 1%，均低于 ≥40 单元和 ≥61 个基因门槛。故状态为 `STOP_NO_FULL_COHORT_REFIT`。参考供体独立性仍未解决，此结果只用于技术诊断。",
        "",
        "## WP2：GSE269906 CellBin",
        "",
        f"12/12 供体均完成直接核级计数构建和输入 SHA-256 核对；`/cellBin/cellExp` 与转置存储 `/cellBin/geneExp` 的逐核×基因记录完全一致。共有 `{total_gt150:,}` 个核达到 ≥150 检出基因；在排除无效 polygon 后，进入盲态聚类/注释的核为 `{valid_ge150:,}`。CellBin 计数占 tissue bin1 总 UMI 的供体范围为 `{fmt_pct(min(total_umi_fraction))}`–`{fmt_pct(max(total_umi_fraction))}`；tissue-bin1 与 CellBin 矩阵的 UMI 差额单独记录，不将其解释为未分配 DNB，也不用于核身份判定。",
        "",
        f"聚类得到 `{len(cluster_rows)}` 个簇；QC 核标签数为 MG `{label_counts.get('MG', 0):,}`、ASC `{label_counts.get('ASC', 0):,}`、Ex `{label_counts.get('Ex', 0):,}`、Oli `{label_counts.get('Oli', 0):,}`。MG 身份验证方向为 `{id_map['MG']['donors_with_direction']}/{id_map['MG']['donors_total']}` 供体，ASC 为 `{id_map['ASC']['donors_with_direction']}/{id_map['ASC']['donors_total']}`，未达到 10/12。",
        f"全体 QC 合格核的标记可测性诊断显示：CD74/CD14 两个 MG 训练标记均在 ≥1% 核中检出的供体为 `{marker_donors['MG_training']}/12`；MG 独立留出标记中至少 2 个达到该比例的供体为 `{marker_donors['MG_holdout_disjoint_from_G1_G5']}/12`。对应 ASC 数字为 `{marker_donors['ASC_training']}/12` 和 `{marker_donors['ASC_holdout_disjoint_from_G1_G5']}/12`。该诊断不重标记细胞，也不改身份门槛；它支持 MG 标记稀疏/聚类未恢复 MG 的解释，但不能证明组织中没有小胶质细胞。",
        f"G1 只有 `{g1_eligible_donors}/12` 供体达到 61/244 覆盖；由于没有 MG 标签，零覆盖不能解释为 G1 基因在生物样本中缺失。最弱 G5 模块 `{g5_worst_module}` 仅 `{g5_best}/12` 供体通过。空间几何审计在 10/25/50 μm 半径均未找到合格 MG–ASC 对；同时无可审计组织空洞掩膜/配对组织学图像，邻域状态为 `NEIGHBORHOOD_NOT_ELIGIBLE_MASK_UNAVAILABLE`。",
        "",
        "本地作者 Scanpy 代码/锁定环境不可核验，因此核聚类采用单独标记的 Seurat/Harmony 探索性重建，不声称精确复现；冻结 G1/G5 与训练/留出标记未进入聚类特征。",
        "",
        "## WP3：候选队列 metadata",
        "",
        f"- GSE233208：AD `{next(r['n_AD'] for r in candidates if r['candidate']=='GSE233208')}`、对照 `{next(r['n_Control'] for r in candidates if r['candidate']=='GSE233208')}`；CaptureArea 与诊断完全混杂，设计矩阵不满秩，且参考独立性缺少精确 crosswalk；`{candidate_status['GSE233208']}`。",
        f"- GSE220442：3 AD + 3 对照，但 GEO GSM 到作者 patientID 映射为 0/6；`{candidate_status['GSE220442']}`。",
        f"- Zenodo 16938034 CosMx 与 Visium HD：逐供体诊断、脑区/技术设计和冻结面板覆盖仍不完整；状态分别为 `{candidate_status['Zenodo_16938034_CosMx']}`、`{candidate_status['Zenodo_16938034_VisiumHD']}`。",
        "",
        "没有候选达到 `PRIMARY_CANDIDATE`。本包只读 metadata，没有打开 Zenodo 表达对象。",
        "",
        "## 执行记录与复现",
        "",
        "WP3 metadata-only 筛选在 WP1/WP2 长时计算期间先行运行；它不读取表达矩阵，不向盲态分析传递诊断结果，因此不改变前两包。CellBin 首轮无效矩阵因 R 列名赋值错误未保留；修正后 12/12 基因名、核名、UMI 与基因数完整性检查均通过。GSM8330064 首次启动因未设置隔离 R 库而在载入 spacexr 前终止；按 SOP 设置 `R_LIBS_USER` 后完整重跑。原始 GEF 未修改，系统 TEMP 中的 ASCII 副本已由构建器清理。",
        "",
        f"代码、参数、结果表、日志和核级/模型中间对象纳入项目 Git/DVC 跟踪。报告汇总 stage 为 `g6_next_feasibility_summary`；计算 runner 按 SOP 命令本地执行，正式 interim/results 目录由 DVC 指针跟踪。DVC `{config['execution_result']['DVC_version']}` 使用本地 cache，`DVC_SITE_CACHE_DIR` 指向 `{config['execution_result']['DVC_SITE_CACHE_DIR']}`；本项目目前没有配置 remote。Git 基线为 `{config['execution_result']['git_base_revision']}`；开始执行时工作区已有变更，本轮未创建提交。只使用本机现有 R 4.3.3 / spacexr 2.2.1、Seurat 5.5.0、Harmony 2.0.5、irlba 与现有 Python 环境，未安装依赖、未下载大型表达矩阵、未使用远程计算。",
        "",
        "组件结果和运行签名见本目录下的 `results/2026-09-27-g6-next-feasibility/`；本报告的源文件/参数和组件审计表 SHA-256 列于配套 manifest。",
    ])
    out_report.write_text("\n".join(report) + "\n", encoding="utf-8")

    files = [
        "src/89_fit_g6_numi_sensitivity.R", "src/90_audit_g6_numi_sensitivity.R",
        "src/91_build_g6_cellbin_nucleus_counts.R", "src/92_annotate_g6_cellbin_nuclei.R",
        "src/93_audit_g6_cellbin_neighborhoods.py", "src/94_audit_g6_candidate_designs.py",
        "src/95_audit_g6_cellbin_marker_availability.R", "src/96_finalize_g6_next_feasibility.py",
        "configs/g6_next_blind_feasibility_001.yaml", "docs/protocol/2026-09-27-g6-next-feasibility-protocol.md",
        "docs/protocol/2026-09-27-g6-next-feasibility-sop.md", "docs/protocol/g6_next_feasibility_001.sha256",
        "docs/decisions/2026-09-27-g6-next-feasibility-execution.md",
        "data/interim/G6_next_feasibility_001.dvc", "results/2026-09-27-g6-next-feasibility.dvc",
        "results/2026-09-27-g6-next-feasibility/numi/fit_integrity.tsv",
        "results/2026-09-27-g6-next-feasibility/numi/mechanical_gate_by_donor.tsv",
        "results/2026-09-27-g6-next-feasibility/cellbin/cellbin_nucleus_counts_qc.tsv",
        "results/2026-09-27-g6-next-feasibility/cellbin/identity_marker_gate.tsv",
        "results/2026-09-27-g6-next-feasibility/cellbin/marker_availability_summary.tsv",
        "results/2026-09-27-g6-next-feasibility/cellbin/G1_coverage_gate_by_donor.tsv",
        "results/2026-09-27-g6-next-feasibility/cellbin/G5_coverage_gate_by_module.tsv",
        "results/2026-09-27-g6-next-feasibility/cellbin/neighborhood_geometric_gate.tsv",
        "results/2026-09-27-g6-next-feasibility/candidates/candidate_design_summary.tsv",
    ]
    manifest = []
    for relative in files:
        path = root / relative
        if not path.is_file():
            raise FileNotFoundError(f"Missing report dependency: {relative}")
        manifest.append({"path": relative, "bytes": path.stat().st_size, "sha256": sha256(path)})
    with out_manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Wrote report: {out_report.relative_to(root)}")
    print(f"Wrote SHA-256 manifest: {out_manifest.relative_to(root)}; files={len(manifest)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

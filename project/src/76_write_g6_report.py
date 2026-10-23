#!/usr/bin/env python3
"""Render the G6 audit and final analysis tables as a UTF-8 Chinese report."""
from __future__ import annotations

import csv
import datetime as dt
import json
import math
import sys
from pathlib import Path


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig")) if path.exists() else None


def read_csv(path: Path):
    if not path.exists():
        return None
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def number(value, digits=3):
    try:
        x = float(value)
        return f"{x:.{digits}g}" if math.isfinite(x) else "NA"
    except (TypeError, ValueError):
        return "NA"


def pvalue(value):
    try:
        x = float(value)
        if not math.isfinite(x):
            return "NA"
        return "<0.001" if x < 0.001 else f"{x:.3g}"
    except (TypeError, ValueError):
        return "NA"


def localize_blocker(message, preflight, g1_gate):
    if message.startswith("Frozen G1 coverage failed") and g1_gate:
        row = g1_gate[0]
        return (f"G1 发送端覆盖未通过：{row['G1_supported_genes']}/{row['G1_signature_genes']} 个冻结基因达到覆盖要求；"
                f"至少需要 {row['G1_minimum_supported_genes']} 个")
    if "lack at least one non-rejected bin with MG and ASC weights >=0.5" in message:
        failed = preflight.get("donors_failing_primary_weight_support", []) if preflight else []
        no_mg = sum(int(row.get("n_MG_weight_ge_0_5", 0)) == 0 for row in failed)
        no_asc = sum(int(row.get("n_ASC_weight_ge_0_5", 0)) == 0 for row in failed)
        return (f"主邻域权重支持未通过：{len(failed)}/12 位供体缺少至少一种权重 ≥0.5 的谱系 bin；"
                f"其中 MG 缺失 {no_mg} 位，ASC 缺失 {no_asc} 位")
    if message.startswith("Disjoint holdout-marker"):
        return "不重叠留出谱系标记门未通过"
    if message.startswith("No passing G5 module"):
        return "没有通过 G5 且达到 G6 覆盖要求的可估模块"
    if message.startswith("GSE233208"):
        return "GSE233208 与 GSE174367 的候选供体重叠经保守排除后，独立 AD 供体为 0 名"
    if message.startswith("GSE220442"):
        return "GSE220442 的 GSM—患者映射仍未解决，不能作为供体级验证"
    return message


def md_table(rows, columns, headers):
    if not rows:
        return "（无可报告记录）"
    lines = ["| " + " | ".join(headers) + " |", "|" + "|".join(["---"] * len(headers)) + "|"]
    lines.extend("| " + " | ".join(str(row.get(col, "NA")) for col in columns) + " |" for row in rows)
    return "\n".join(lines)


def main(root: Path):
    audit = root / "results/2026-09-24-disease-side-v3.5-g6-input-audit"
    out = root / "results/2026-09-24-disease-side-v3.5-g6-primary"
    (out / "tables").mkdir(parents=True, exist_ok=True)
    preflight = read_json(out / "G6_effect_blind_preflight.json")
    decision = read_json(out / "G6_decision.json")
    batch_status = read_json(audit / "logs/G6_RCTD_batch_resume_20260925_status.json")
    overlap = read_json(audit / "GSE269906_GSE174367_overlap_audit.json")
    input_audit = read_json(audit / "G6_input_audit.json")
    validation = input_audit["independent_validation"]
    section = read_tsv(audit / "tables/GSE269906_tissue_bin50_section_qc.tsv")
    mask = read_tsv(audit / "tables/GSE269906_bin50_CellBin_polygon_mask_qc.tsv")
    marker = read_csv(out / "tables/G6_RCTD_holdout_marker_gate.csv") or []
    donor_qc = read_csv(out / "tables/G6_RCTD_donor_fit_qc.csv") or []
    raw_weight_qc = read_csv(out / "tables/G6_RCTD_weight_sum_qc_all_donors.csv") or []
    g1_gate = read_csv(out / "tables/G6_G1_coverage_gate.csv") or []
    modules = read_csv(out / "tables/G6_G5_module_coverage_gate.csv") or []
    state = read_csv(out / "tables/G6_state_group_effects.csv") or []
    abundance = read_csv(out / "tables/G6_cell_type_abundance_group_effects.csv") or []
    neighborhood = read_csv(out / "tables/G6_neighborhood_group_effects.csv") or []
    reference_sensitivity = read_csv(out / "tables/G6_reference_donor_sensitivity_status.csv") or []
    status = decision["decision"] if decision else (
        "G6_ANALYSIS_INCOMPLETE" if preflight and preflight.get("preflight_pass") else "G6_NOT_EXECUTABLE"
    )
    total_bins = sum(int(row["n_bin50_ge400_and_centroid_supported"]) for row in section)
    total_polygons = sum(int(row["n_valid_polygons"]) for row in mask)

    lines = [
        "# G6：BA10 空间状态与微胶质细胞邻域分析",
        "",
        f"运行日期：{dt.date.today().isoformat()}；状态：**{status}**。",
        "",
        "## 问题与冻结决策",
        "",
        "检验 GSE269906 BA10 中，冻结的 G1 微胶质细胞 AD 方向状态是否改变，以及冻结的 G5 星形胶质细胞接收状态是否在空间上邻近高 G1 微胶质状态。统计单位为供体。25 µm Stereo-seq bin 是混合表达单元，RCTD 权重与状态均为推断值，不能解释为分割出的单细胞或因果效应。",
        "",
        "数据角色与放行规则见[数据集角色修订 001](../../docs/protocol/disease_side_v3.5_G6_dataset_roles_amendment_001.md)、[数据集选择决策](../../docs/decisions/2026-09-23-g6-spatial-dataset-selection.md)、[执行锁 001](../../docs/decisions/2026-09-24-g6-v3.5-execution-lock.md)、[执行附注 002](../../docs/decisions/2026-09-25-g6-execution-addendum-002.md)、[执行附注 003](../../docs/decisions/2026-09-25-g6-execution-addendum-003.md)、[执行附注 004](../../docs/decisions/2026-09-25-g6-execution-addendum-004.md)、[执行附注 005](../../docs/decisions/2026-09-25-g6-execution-addendum-005.md)和[执行附注 006](../../docs/decisions/2026-09-25-g6-execution-addendum-006.md)。父版冻结协议正文未被覆盖。",
        "空间零分布的固定接收点分母实现修正见[执行附注 007](../../docs/decisions/2026-09-25-g6-execution-addendum-007.md)。",
        "分开报告 RCTD 细胞组成与表达状态的预设敏感性见[执行附注 008](../../docs/decisions/2026-09-25-g6-execution-addendum-008.md)；效应盲预检入口及执行边界见[执行附注 009](../../docs/decisions/2026-09-25-g6-execution-addendum-009.md)；拒绝 bin 原始权重与分析权重的校验口径修订见[执行附注 010](../../docs/decisions/2026-09-25-g6-execution-addendum-010.md)；完整记录权重阈值支持并保持诊断标签盲态见[执行附注 011](../../docs/decisions/2026-09-25-g6-execution-addendum-011.md)。",
        "",
        "## 数据与盲态质控",
        "",
        f"主队列 GSE269906：{len(section)} 名 BA10 供体（作者标签 AD 6、Normal Aging 对照 6）。根据 tissue GEF 重建全基因 50×50 DNB、25 µm bin；保留检测基因数 ≥400、与有效 CellBin 多边形相交且至少有一个分割单元中心的 bin。共 {total_bins:,} 个合格 bin；12 名供体均达到每供体 5,000 bin 下限。",
        "",
        f"有效 CellBin 多边形共 {total_polygons:,} 个；cellTypeID 全部为默认 0。细胞组成因此由 GSE174367 的 18 供体人皮层 snRNA-seq 参考和 RCTD 双细胞模式推断。模型排除所有冻结 G1、G5 基因及预设不重叠谱系标记。RCTD 版本、固定提交、R 版本与 sessionInfo 见软件清单。",
        "RCTD 按冻结规则使用空间与参考共有、且剔除 G1/G5/留出标记后的特征矩阵，同时保留每个 bin 的原始全基因 UMI 总数作为归一化分母。spacexr 因 `nUMI` 与受限特征矩阵列和不同而发出提示；该差异由预设特征排除造成，分析保留全库 UMI 分母并在报告中披露。",
        "",
        f"GSE269906 与 GSE174367 有 {overlap['demographic_diagnosis_candidate_pairs']} 组仅按性别、诊断和五年年龄段匹配的候选组合；未找到作者交叉表，且没有候选组同时匹配可用 Braak 和斑块分期。个体是否重复无法确认。执行附注 002 保留锁定参考，但不称其为已证实独立参考。",
        "",
        "盲态谱系标记相关性门槛：",
        "",
    ]
    if raw_weight_qc:
        raw_medians = [float(row["raw_weight_sum_median"]) for row in raw_weight_qc]
        raw_outside = sum(int(row["n_nonrejected_weight_sum_outside_098_102"]) for row in raw_weight_qc)
        usable_bins = sum(int(row["n_bins"]) - int(row["n_rejected"]) for row in raw_weight_qc)
        lines.insert(18, f"原始 RCTD 权重和的供体中位数范围为 {min(raw_medians):.3f}–{max(raw_medians):.3f}；{usable_bins:,} 个非拒绝 bin 中有 {raw_outside:,} 个原始权重和不在 0.98–1.02。该描述性 QC 与使用前逐行归一化相分离；完整逐供体表见 `tables/G6_RCTD_weight_sum_qc_all_donors.csv`。")
    rctd_stderr_log = audit / "logs/GSM8330070_RCTD.stderr.log"
    if batch_status and batch_status.get("outcome") == "complete" and rctd_stderr_log.exists():
        rctd_log_text = rctd_stderr_log.read_text(encoding="utf-8", errors="replace")
        if "There were 13 warnings" in rctd_log_text:
            lines.insert(19, "GSM8330070 拟合日志另有 13 条警告的汇总提示，但没有保留逐条警告正文；批次状态为 exit code 0，13/13 分块及最终权重/QC 文件均已生成，因此警告的具体内容仍是未能核实的限制。")
    if marker:
        marker_rows = [{
            "lineage": "小胶质" if row["lineage"] == "MG" else "星形胶质",
            "donors_positive_rho": row["donors_positive_rho"],
            "minimum_disjoint_markers": row["minimum_disjoint_markers"],
            "median_rho": number(row["median_rho"]),
            "pass": "通过" if row["pass"].lower() == "true" else "未通过",
        } for row in marker]
        lines.append(md_table(marker_rows, ["lineage", "donors_positive_rho", "minimum_disjoint_markers", "median_rho", "pass"],
                               ["谱系", "正相关供体/12", "最少不重叠标记数", "中位 Spearman ρ", "门槛"]))
    else:
        lines.append("尚无完整 12 供体标记审计结果。")
    lines.extend(["", "G1 发送端覆盖门槛：每个基因须在 MG 权重 ≥0.5 的非拒绝 bin 中达到 ≥1% 检出，并在至少 10/12 供体满足；244 个冻结基因至少需 61 个通过：", ""])
    if g1_gate:
        row = g1_gate[0]
        g1_row = [{
            "G1_signature_genes": row["G1_signature_genes"],
            "G1_supported_genes": row["G1_supported_genes"],
            "G1_minimum_supported_genes": row["G1_minimum_supported_genes"],
            "G1_coverage_pass": "通过" if row["G1_coverage_pass"].lower() == "true" else "未通过",
        }]
        lines.append(md_table(g1_row,
                              ["G1_signature_genes", "G1_supported_genes", "G1_minimum_supported_genes", "G1_coverage_pass"],
                              ["冻结基因总数", "达到覆盖基因数", "最低要求", "门槛"]))
    else:
        lines.append("尚无 G1 覆盖门槛结果。")
    lines.extend(["", "主分析的 0.5 RCTD 权重阈值支持（全部供体均需有 MG 与 ASC 支持 bin）：", ""])
    if donor_qc:
        support_rows = [{
            "GSM": row["GSM"], "donor_id": row["donor_id"],
            "n_MG_weight_ge_0_5": row["n_MG_weight_ge_0_5"],
            "n_ASC_weight_ge_0_5": row["n_ASC_weight_ge_0_5"],
            "n_rejected": row["n_rejected"], "n_usable": row["n_usable"],
            "pass": "通过" if int(row["n_MG_weight_ge_0_5"]) > 0 and int(row["n_ASC_weight_ge_0_5"]) > 0 else "未通过",
        } for row in donor_qc]
        lines.append(md_table(support_rows,
                              ["GSM", "donor_id", "n_MG_weight_ge_0_5", "n_ASC_weight_ge_0_5", "n_rejected", "n_usable", "pass"],
                              ["GSM", "供体", "MG ≥0.5 bin", "ASC ≥0.5 bin", "拒绝 bin", "可用 bin", "门槛"]))
    else:
        lines.append("尚无完整 12 供体权重支持审计结果。")
    if preflight:
        blockers = preflight.get("preflight_blockers", [])
        if isinstance(blockers, str):
            blockers = [blockers]
        if blockers:
            lines.extend(["", "效应盲态预检未通过：" + "；".join(localize_blocker(x, preflight, g1_gate) for x in blockers) + "。"])
        else:
            lines.extend(["", "效应盲态预检：通过。"])
    lines.extend(["", "G5 模块覆盖门槛：至少 25% 的一对一基因在 ASC 权重 ≥0.5 的非拒绝 bin 中达到 ≥1% 检出，且该基因在 ≥10 个供体满足。只有 G5 v3.4 通过且 G6 覆盖可估的模块能支持 G6 主门；表格分别列出 G5 状态、G5 门和 G6 覆盖结果：", ""])
    if modules:
        module_rows = [{**row, "G6_state_estimable": "是" if row["G6_state_estimable"].lower() == "true" else "否"}
                       for row in modules]
        lines.append(md_table(module_rows,
                              ["module", "n_frozen_one_to_one_genes", "n_genes_supported_in_ge10_donors", "minimum_genes_required", "G5_module_gate", "G6_state_estimable"],
                              ["G5 模块", "一对一基因数", "≥10 供体达覆盖基因数", "最低基因数", "G5 门", "G6 可估"]))
    else:
        lines.append("尚未完成全部供体的状态覆盖审计。")

    lines.extend([
        "", "## 方法与统计", "",
        "盲态下按冻结规则将每个状态基因的 log2(CPM+1) 跨所有合格 bin 标准化。G1 按逐基因 AD 方向取带符号均值；G5 模块按冻结模块方向取均值。供体状态分数按对应连续 RCTD 权重加权。主模型为供体级 `状态 ~ 诊断 + 年龄中点`，并报告 Welch 未调整差及 95% CI；G1 加 7 个 G5 状态构成 BH-FDR 家族。另以每名供体非拒绝 bin 的平均归一化 MG/ASC RCTD 权重单独描述组成，按同一供体模型估计，并在两项组成指标内单独作 BH 校正；组成敏感性不进入 G6 判定门。",
        "",
        "高状态定义为：供体内、对应细胞类型权重 ≥0.5 的 bin 中分数最高四分位。邻域端点为高 G5 星形胶质 bin 中，25 µm 内至少有一个高 G1 微胶质 bin 的比例，减去 999 次受掩膜约束的整张接收状态图平移零分布均值。平移 250–1,000 µm，掩膜重叠至少 80%；连线的每个中间 25 µm 网格均须在组织掩膜内。每个模块的供体级 AD 减对照差以年龄调整线性模型估计；FDR 在每个半径和权重阈值的 7 个模块间校正。50/75 µm 半径和权重 0.4/0.7 为敏感性分析。",
        "", "## G1 与 G5 状态结果", "",
    ])
    if state:
        state_rows = [{**row,
                       "AD_minus_Control": number(row["AD_minus_Control"]),
                       "age_adjusted_effect": number(row["age_adjusted_effect"]),
                       "age_adjusted_CI_low": number(row["age_adjusted_CI_low"]),
                       "age_adjusted_CI_high": number(row["age_adjusted_CI_high"]),
                       "age_adjusted_p": pvalue(row["age_adjusted_p"]),
                       "FDR_BH_8": pvalue(row["FDR_BH_8"]),
                       "direction_concordant": "是" if row["direction_concordant"].lower() == "true" else "否",
                       "state_estimable": "是" if row["state_estimable"].lower() == "true" else "否"}
                      for row in state]
        lines.append(md_table(state_rows,
                              ["state", "n_AD", "n_Control", "AD_minus_Control", "age_adjusted_effect", "age_adjusted_CI_low", "age_adjusted_CI_high", "age_adjusted_p", "FDR_BH_8", "direction_concordant", "state_estimable"],
                              ["状态", "AD n", "对照 n", "未调整差", "年龄调整差", "95% CI 下限", "95% CI 上限", "年龄调整 p", "BH-FDR", "方向一致", "可估"]))
    else:
        lines.append("未计算诊断组间状态效应（盲态放行门槛未通过或分析尚未完成）。")

    lines.extend(["", "## 细胞组成敏感性（RCTD 权重）", ""])
    if abundance:
        abundance_rows = [{**row,
                           "AD_minus_Control": number(row["AD_minus_Control"]),
                           "age_adjusted_effect": number(row["age_adjusted_effect"]),
                           "age_adjusted_CI_low": number(row["age_adjusted_CI_low"]),
                           "age_adjusted_CI_high": number(row["age_adjusted_CI_high"]),
                           "age_adjusted_p": pvalue(row["age_adjusted_p"]),
                           "FDR_BH_2": pvalue(row["FDR_BH_2"])} for row in abundance]
        lines.append(md_table(abundance_rows,
                              ["cell_type", "n_AD", "n_Control", "AD_minus_Control", "age_adjusted_effect", "age_adjusted_CI_low", "age_adjusted_CI_high", "age_adjusted_p", "FDR_BH_2"],
                              ["推断组成", "AD n", "对照 n", "未调整差", "年龄调整差", "95% CI 下限", "95% CI 上限", "年龄调整 p", "两项 BH-FDR"]))
    else:
        lines.append("尚无完整的供体级 MG/ASC 平均 RCTD 权重对照结果。")

    lines.extend(["", "## 微胶质—星形胶质邻域结果", ""])
    primary = [row for row in neighborhood if float(row["radius_um"]) == 25 and float(row["cell_weight_threshold"]) == 0.5]
    if primary:
        spatial_rows = [{**row,
                         "AD_minus_Control": number(row["AD_minus_Control"]),
                         "age_adjusted_effect": number(row["age_adjusted_effect"]),
                         "age_adjusted_CI_low": number(row["age_adjusted_CI_low"]),
                         "age_adjusted_CI_high": number(row["age_adjusted_CI_high"]),
                         "age_adjusted_p": pvalue(row["age_adjusted_p"]),
                         "FDR_BH_7_within_radius": pvalue(row["FDR_BH_7_within_radius"]),
                         "spatial_gate_pass": "通过" if row["spatial_gate_pass"].lower() == "true" else "未通过"}
                        for row in primary]
        lines.append(md_table(spatial_rows,
                              ["module", "n_AD", "n_Control", "AD_minus_Control", "age_adjusted_effect", "age_adjusted_CI_low", "age_adjusted_CI_high", "age_adjusted_p", "FDR_BH_7_within_radius", "spatial_gate_pass"],
                              ["G5 模块", "AD n", "对照 n", "未调整差", "年龄调整空间超额差", "95% CI 下限", "95% CI 上限", "年龄调整 p", "BH-FDR", "主门"]))
    else:
        lines.append("未计算诊断组间邻域效应（盲态放行门槛未通过或分析尚未完成）。")

    final_decision = decision["decision"] if decision else "主分析尚未生成最终判定。"
    final_note = "；".join(localize_blocker(x, preflight, g1_gate) for x in decision.get("blockers", [])) if decision else ""
    lines.extend([
        "", ("本次效应盲态放行门未通过；`tables/` 保存完整盲态 QC 和覆盖审计。未生成诊断组间状态、组成或邻域效应，Not-ADNC、阈值/半径、逐供体及参考供体敏感性也未运行，`figures/` 不含组间效应图。" if decision and not decision.get("disease_effects_calculated", False) else "供体级状态与空间指标、覆盖审计、Not-ADNC 对照敏感性、阈值/半径敏感性和逐供体剔除结果见本目录 `tables/`；点图见 `figures/`。"),
        "", "## 验证与最终判定", "",
        f"GSE233208 的作者元数据更正和样本映射已审计；按与 GSE174367 的候选供体重叠执行 fail-closed 排除后，只剩 {validation['eligible_independent_AD_donors_after_fail_closed_exclusion']} 名独立 AD 与 {validation['eligible_independent_Control_donors_after_fail_closed_exclusion']} 名对照，不能验证 AD 方向。GSE220442 的 GSM—患者映射未解决，不纳入验证；AD-only 的 GSE263791 只保留为描述背景。",
        "", f"冻结协议判定：**{final_decision}**。{final_note}{'。' if final_note else ''}", "",
        "由于没有合格独立 AD 验证队列，即便主队列结果符合方向和显著性条件，也不能达到冻结协议的 G6_PASS。空间共现不代表因果。另有 6 对 6 小样本、潜在参考供体重叠未解、RCTD 原始平台验证未覆盖 Stereo-seq、无组织学图像或皮层层级/ROI、PMI/RIN 信息不可用于校正等限制。参考供体逐一剔除敏感性是否运行以对应状态表为准。",
        "", "## 数据来源与复现", "",
        "- 主队列：[GSE269906](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269906)；[原论文](https://doi.org/10.1038/s41467-024-54715-y)。",
        "- 验证候选：[GSE233208](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233208)；[作者处理说明](https://github.com/swaruplabUCI/DSAD_Spatial_Miyoshi_Morabito_2024)。",
        "- 跨脑区候选：[GSE220442](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220442)；[原论文](https://doi.org/10.1186/s40478-022-01494-6)。",
        "- 解卷积方法：[RCTD 方法论文](https://doi.org/10.1038/s41587-021-00830-w)；[固定 spacexr 提交](https://github.com/dmcable/spacexr/tree/9f5dc33c8060f946c6072a138b70e189636e1435)。",
        "- 本地实验技术库审阅了空间质量控制、邻域统计和细胞组成解卷积三份方法卡；均为未在本机验证的学习卡，无可直接复用代码。本分析复用了其方法建议并按冻结锁实现。",
        "- 本地 Windows 执行；R 4.3.3 和隔离 R 包库，RCTD 2.2.1 固定提交。原始数据与正式中间层已加入本地 DVC 缓存；未配置 DVC 远端。",
        "- 全量盲态解卷积：`powershell -ExecutionPolicy Bypass -File src/74_fit_g6_all_primary_rctd.ps1`；分析和报告：`dvc repro g6_primary_analysis`。",
        "- 报告生成器：`src/76_write_g6_report.py`。",
    ])
    if reference_sensitivity:
        ref_row = reference_sensitivity[0]
        if str(ref_row.get("performed", "")).lower() != "true":
            lines.extend(["", f"参考供体逐一剔除 RCTD 敏感性未运行：{ref_row.get('reason', '原因未记录')}。该项需要 18 组参考重建 × 12 个主队列供体的额外拟合，因此本地结果范围不包含该敏感性。"])
    report = out / "G6_report.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"Wrote G6 report: {report}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 76_write_g6_report.py <project_root>")
    main(Path(sys.argv[1]))

#!/usr/bin/env python
"""Re-audit GSE269906 against G6 plan 005 without linking diagnosis to expression."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import sys
from collections import Counter
from pathlib import Path

import yaml
from openpyxl import load_workbook


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def require_no_diagnosis_columns(paths: list[Path]) -> None:
    forbidden = {"diagnosis", "source_diagnosis", "group", "case_control", "disease_status"}
    for path in paths:
        rows = read_tsv(path)
        columns = set(rows[0]) if rows else set()
        if forbidden.intersection(column.lower() for column in columns):
            raise RuntimeError(f"Diagnosis/group column found in blind output: {path}")


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 100_reassess_g6_gse269906_plan005.py <project_root>")
    root = Path(sys.argv[1]).resolve()
    config_path = root / "configs/g6_gse269906_plan005_reassessment_001.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if config["scope"]["diagnosis_group_effects"] is not False:
        raise RuntimeError("Configuration must prohibit group effects")
    if config["scope"]["diagnosis_linked_cell_or_spatial_analysis"] is not False:
        raise RuntimeError("Diagnosis must remain disconnected from cell/spatial outputs")
    if config["protocol_scope"]["can_test_plan005_hippocampal_estimand"] is not False:
        raise RuntimeError("GSE269906 cannot be represented as a hippocampal test")

    out_dir = root / config["outputs"]["directory"]
    out_dir.mkdir(parents=True, exist_ok=True)

    # Source/design metadata are summarized separately from expression/spatial outputs.
    donor_manifest = root / config["dataset"]["donor_manifest"]
    donors = read_tsv(donor_manifest)
    donor_ids = {row["donor_id"] for row in donors}
    gsm_ids = {row["GSM"] for row in donors}
    diagnosis_counts = Counter(row["source_diagnosis"] for row in donors)
    regions = sorted({row["region"] for row in donors})
    if len(donors) != 12 or len(donor_ids) != 12 or len(gsm_ids) != 12:
        raise RuntimeError("Expected 12 unique donor-to-GSM rows")
    if diagnosis_counts != Counter({"AD": 6, "Control": 6}):
        raise RuntimeError(f"Unexpected source diagnosis counts: {diagnosis_counts}")
    if regions != ["BA10 prefrontal cortex"]:
        raise RuntimeError(f"Unexpected region metadata: {regions}")

    # Recheck source metadata separately from all cell/spatial measurements.
    soft_path, supplement_path = [root / value for value in config["dataset"]["source_metadata_files"]]
    with gzip.open(soft_path, "rt", encoding="utf-8", errors="replace") as handle:
        soft_text = handle.read()
    sample_starts = list(re.finditer(r"(?m)^\^SAMPLE = (GSM\d+)\s*$", soft_text))
    soft_records: dict[str, str] = {}
    technical_terms = re.compile(r"batch|chip|capture[_ ]?area|library[_ ]?(id|batch)|sequenc(e|ing)[_ ]?(run|batch)|flow[_ ]?cell|lane|slide|section|run[_ ]?id", re.I)
    protocol_descriptor_terms = re.compile(r"library[_ ]?(source|selection|strategy)", re.I)
    technical_fields: set[str] = set()
    protocol_descriptors: set[str] = set()
    for index, match in enumerate(sample_starts):
        end = sample_starts[index + 1].start() if index + 1 < len(sample_starts) else len(soft_text)
        gsm = match.group(1)
        block = soft_text[match.start():end]
        soft_records[gsm] = block
        for line in block.splitlines():
            if line.startswith("!Sample_"):
                field = line.split("=", 1)[0].strip()
                if technical_terms.search(field):
                    technical_fields.add(field)
                if protocol_descriptor_terms.search(field):
                    protocol_descriptors.add(field)
    if set(soft_records) != gsm_ids:
        raise RuntimeError("GEO SOFT sample records do not match the donor manifest")
    if any("BA10" not in re.search(r"(?m)^!Sample_title = (.*)$", block).group(1) for block in soft_records.values()):
        raise RuntimeError("GEO SOFT sample titles do not consistently identify BA10")

    workbook = load_workbook(supplement_path, read_only=True, data_only=True)
    sheet = workbook[workbook.sheetnames[0]]
    headers = [str(value).strip() if value is not None else "" for value in next(sheet.iter_rows(min_row=2, max_row=2, values_only=True))]
    supplement_technical_fields = [header for header in headers if technical_terms.search(header)]
    supplement_rows = list(sheet.iter_rows(min_row=3, values_only=True))
    supplement_rows = [row for row in supplement_rows if row and row[0] is not None and len(row) > 1 and row[1] in {"AD", "Control"}]
    workbook.close()
    supplement_ids = {str(row[0]).strip() for row in supplement_rows}
    if supplement_ids != donor_ids:
        raise RuntimeError("Supplementary clinical sample IDs do not match the 12-donor manifest")
    if len(technical_fields) == 0 and len(supplement_technical_fields) == 0:
        technical_status = "NOT_REPORTED_DESIGN_RANK_NOT_ASSESSED"
    else:
        technical_status = "TECHNICAL_FIELDS_PRESENT_REVIEW_REQUIRED"
    source_design_rows = [
        {"audit_item": "GEO_SOFT_sample_mapping", "status": "PASS", "result": f"{len(soft_records)} sample records; titles consistently BA10; GSMs match donor manifest"},
        {"audit_item": "Supplementary_Data1_sample_mapping", "status": "PASS", "result": f"{len(supplement_ids)} unique sample IDs match donor manifest; clinical fields only"},
        {"audit_item": "library_protocol_descriptors", "status": "AVAILABLE_NOT_BATCH_IDS", "result": f"GEO provides protocol descriptors={sorted(protocol_descriptors)}; these are not per-sample library/batch identifiers"},
        {"audit_item": "technical_batch_crosswalk", "status": "HOLD_METADATA", "result": f"GEO batch fields={sorted(technical_fields)}; supplement batch columns={supplement_technical_fields}; no auditable chip/run/library-ID/batch crosswalk"},
        {"audit_item": "design_matrix_rank", "status": "NOT_ASSESSED", "result": f"{technical_status}; diagnosis-versus-technical-factor identifiability cannot be assessed"},
    ]
    write_tsv(root / config["outputs"]["source_design_audit"], source_design_rows, ["audit_item", "status", "result"])

    # Verify every paired raw GEF against the hash manifest created by the prior blind run.
    input_manifest = read_tsv(root / config["dataset"]["raw_hash_manifest"])
    raw_rows: list[dict[str, object]] = []
    for row in input_manifest:
        if row["modality"] == "CellBin":
            raw_path = root / config["dataset"]["raw_cellbin_dir"] / row["file"]
        elif row["modality"] == "tissue_bin1":
            raw_path = root / config["dataset"]["raw_tissue_dir"] / row["file"]
        else:
            raise RuntimeError(f"Unregistered raw modality: {row['modality']}")
        if not raw_path.is_file():
            actual_bytes, actual_hash, state = 0, "MISSING", "FAIL_MISSING"
        else:
            actual_bytes = raw_path.stat().st_size
            actual_hash = sha256(raw_path)
            state = "PASS" if actual_bytes == int(row["bytes"]) and actual_hash == row["sha256"] else "FAIL_MISMATCH"
        raw_rows.append({
            "GSM": row["GSM"], "donor_id": row["donor_id"], "modality": row["modality"],
            "file": row["file"], "expected_bytes": int(row["bytes"]), "actual_bytes": actual_bytes,
            "expected_sha256": row["sha256"], "actual_sha256": actual_hash, "integrity": state,
        })
    if len(raw_rows) != 24 or any(row["integrity"] != "PASS" for row in raw_rows):
        raise RuntimeError("Raw GEF integrity verification did not pass for all 24 files")
    raw_hash_path = root / config["outputs"]["raw_hash_verification"]
    write_tsv(raw_hash_path, raw_rows, [
        "GSM", "donor_id", "modality", "file", "expected_bytes", "actual_bytes",
        "expected_sha256", "actual_sha256", "integrity",
    ])

    # Verify frozen signature registries against the exact hashes recorded in plan 005.
    frozen_hashes: dict[str, dict[str, str]] = {}
    for asset in ("G1", "G5"):
        asset_cfg = config["frozen_assets"][asset]
        asset_path = root / asset_cfg["registry"]
        actual_hash = sha256(asset_path)
        expected_hash = asset_cfg["expected_sha256"]
        if actual_hash != expected_hash:
            raise RuntimeError(f"Frozen {asset} registry hash mismatch")
        frozen_hashes[asset] = {"path": asset_cfg["registry"], "expected": expected_hash, "actual": actual_hash, "status": "PASS"}

    parent_freeze_path = root / config["run"]["parent_freeze_manifest"]
    parent_freeze_rows: list[dict[str, str]] = []
    for line in parent_freeze_path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        expected_hash, relative = line.split(maxsplit=1)
        relative = relative.strip()
        source_path = root / relative
        actual_hash = sha256(source_path)
        if actual_hash != expected_hash:
            raise RuntimeError(f"Parent plan freeze checksum mismatch: {relative}")
        parent_freeze_rows.append({"path": relative, "sha256": actual_hash, "status": "PASS"})
    if len(parent_freeze_rows) != 6:
        raise RuntimeError("Expected six frozen plan 005 checksum entries")

    input_tables = [root / config["existing_diagnostic_inputs"][key] for key in (
        "identity", "g1_by_donor", "g5_by_module", "neighborhood", "cellbin_qc",
        "cellbin_qc_after_150", "rctd_reference_status"
    )]
    for path in input_tables:
        if not path.is_file():
            raise FileNotFoundError(path)
    require_no_diagnosis_columns([path for path in input_tables if path.parent.name == "cellbin"])

    identity = {row["lineage"]: row for row in read_tsv(root / config["existing_diagnostic_inputs"]["identity"])}
    g1_rows = read_tsv(root / config["existing_diagnostic_inputs"]["g1_by_donor"])
    g5_rows = read_tsv(root / config["existing_diagnostic_inputs"]["g5_by_module"])
    geometry = read_tsv(root / config["existing_diagnostic_inputs"]["neighborhood"])
    cell_qc = read_tsv(root / config["existing_diagnostic_inputs"]["cellbin_qc"])
    cell_qc_after_150 = read_tsv(root / config["existing_diagnostic_inputs"]["cellbin_qc_after_150"])
    ref_rows = read_tsv(root / config["existing_diagnostic_inputs"]["rctd_reference_status"])
    if len(identity) != 2 or len(g1_rows) != 12 or len(cell_qc) != 12 or len(cell_qc_after_150) != 12 or len(geometry) != 3:
        raise RuntimeError("Unexpected row counts in prior blind outputs")
    if {row["donor_id"] for row in cell_qc} != donor_ids:
        raise RuntimeError("CellBin donor IDs do not match the audited source manifest")
    if {row["donor_id"] for row in cell_qc_after_150} != donor_ids:
        raise RuntimeError("Post-QC CellBin donor IDs do not match the source manifest")
    if {row["donors_total"] for row in identity.values()} != {"12"}:
        raise RuntimeError("Identity output does not cover all 12 donors")
    if {row["donor_id"] for row in g1_rows} != donor_ids:
        raise RuntimeError("G1 output donor IDs do not match the source manifest")
    if sum(int(row["n_G1_genes_at_1pct"]) for row in g1_rows) != 0:
        raise RuntimeError("Existing MG-targeted G1 table is not the expected zero-support result")
    if any(row["tissue_hole_mask_available"].lower() != "false" for row in geometry):
        raise RuntimeError("Existing geometry output unexpectedly reports an auditable tissue mask")
    if any(row["neighborhood_gate_status"] != "NEIGHBORHOOD_NOT_ELIGIBLE_MASK_UNAVAILABLE" for row in geometry):
        raise RuntimeError("Unexpected neighborhood gate status")
    if any("REFERENCE_IDENTITY_UNRESOLVED" not in row["reference_status"] for row in ref_rows):
        raise RuntimeError("RCTD reference status changed; new review required")

    mg_identity = identity["MG"]["donors_with_direction"]
    asc_identity = identity["ASC"]["donors_with_direction"]
    g1_mg_donors = sum(int(row["n_G1_genes_at_1pct"]) > 0 for row in g1_rows)
    g5_summary = {row["module"]: int(row["donors_meeting_module_coverage"]) for row in g5_rows}
    nuclei_ge150_genes = sum(int(row["n_cells_ge150_genes"]) for row in cell_qc)
    valid_nuclei = sum(int(row["n_cells_ge150_and_polygon_valid"]) for row in cell_qc_after_150)
    valid_polygons = valid_nuclei
    total_cells = sum(int(row["n_cells_total"]) for row in cell_qc)
    euclidean_pair_support = {row["radius_um"]: int(row["donors_with_ge20_Euclidean_pairs"]) for row in geometry}
    no_mask = all(row["tissue_hole_mask_available"].lower() == "false" for row in geometry)
    reference_status = sorted({row["reference_status"] for row in ref_rows})

    gates = [
        {"gate": "source_donor_mapping", "status": "PASS", "evidence": config["dataset"]["donor_manifest"],
         "result": f"{len(donor_ids)} unique donors; AD={diagnosis_counts['AD']}; Control={diagnosis_counts['Control']}; one GSM per donor"},
        {"gate": "technical_batch_crosswalk_and_design_rank", "status": "HOLD_METADATA", "evidence": config["outputs"]["source_design_audit"],
         "result": "Available GEO SOFT and clinical supplement do not provide an auditable chip/capture-area/run/library-ID/processing-batch crosswalk; design rank is not assessed"},
        {"gate": "plan005_common_hippocampal_region", "status": "FAIL_SCOPE", "evidence": config["dataset"]["donor_manifest"],
         "result": "Dataset region is BA10 prefrontal cortex; plan 005 requires a common hippocampal subregion"},
        {"gate": "raw_cellbin_and_tissue_integrity", "status": "PASS", "evidence": config["dataset"]["raw_hash_manifest"],
         "result": f"{sum(int(row['actual_bytes']) for row in raw_rows):,} bytes across {len(raw_rows)} files; exact byte and SHA-256 match"},
        {"gate": "frozen_G1_G5_registry_integrity", "status": "PASS", "evidence": "configs/g6_staged_research_plan_005.yaml",
         "result": f"G1={frozen_hashes['G1']['actual']}; G5={frozen_hashes['G5']['actual']}"},
        {"gate": "parent_plan005_freeze_integrity", "status": "PASS", "evidence": config["run"]["parent_freeze_manifest"],
         "result": f"{len(parent_freeze_rows)}/{len(parent_freeze_rows)} frozen plan files match their SHA-256 entries"},
        {"gate": "independent_MG_identity", "status": "NOT_CALIBRATED", "evidence": config["existing_diagnostic_inputs"]["identity"],
         "result": f"holdout-marker directional evidence in {mg_identity}/12 donors; no independently validated MG cell labels"},
        {"gate": "independent_ASC_identity", "status": "NOT_CALIBRATED", "evidence": config["existing_diagnostic_inputs"]["identity"],
         "result": f"holdout-marker directional evidence in {asc_identity}/12 donors; insufficient to establish calibrated ASC identity"},
        {"gate": "G1_MG_state_score_and_partial_score_fidelity", "status": "NOT_ESTIMABLE", "evidence": config["existing_diagnostic_inputs"]["g1_by_donor"],
         "result": f"no calibrated MG cell population; {g1_mg_donors}/12 donors have any G1 gene supported in the prior MG-targeted table; score fidelity not testable"},
        {"gate": "G5_ASC_state_score_and_partial_score_fidelity", "status": "NOT_ESTIMABLE", "evidence": config["existing_diagnostic_inputs"]["g5_by_module"],
         "result": "ASC identity was not calibrated; old 25% coverage summaries are not a substitute for 004/005 partial-score fidelity; descriptive legacy donor counts=" + json.dumps(g5_summary, sort_keys=True)},
        {"gate": "physical_neighborhood_and_tissue_mask", "status": "FAIL", "evidence": config["existing_diagnostic_inputs"]["neighborhood"],
         "result": f"auditable tissue-hole mask unavailable={no_mask}; donors with >=20 Euclidean pairs by radius={json.dumps(euclidean_pair_support, sort_keys=True)}"},
        {"gate": "RCTD_reference_independence", "status": "UNRESOLVED", "evidence": config["existing_diagnostic_inputs"]["rctd_reference_status"],
         "result": ";".join(reference_status)},
        {"gate": "independent_same_hippocampal_region_validation", "status": "NOT_FOUND", "evidence": "configs/g6_staged_research_plan_005.yaml",
         "result": "no eligible independent same-hippocampal-region validation cohort in the current registry"},
        {"gate": "AD_control_effect_analysis", "status": "NOT_RUN_PROHIBITED", "evidence": "G6-PLAN-2026-09-27-005 and amendment 006",
         "result": "No diagnosis labels were joined to cell/spatial outputs; no group effect, p-value, FDR, or group plot was generated"},
    ]
    write_tsv(root / config["outputs"]["gates"], gates, ["gate", "status", "evidence", "result"])

    audit = {
        "run_id": config["run"]["id"],
        "parent_plan": config["run"]["parent_plan"],
        "run_status": "G6_NOT_EXECUTABLE",
        "scope_status": "OUT_OF_REGION_TECHNICAL_FEASIBILITY_ONLY",
        "dataset": {"accession": "GSE269906", "region": regions[0], "unique_donors": len(donor_ids),
                    "source_diagnosis_counts": dict(diagnosis_counts)},
        "technical_design_audit": {"status": technical_status, "technical_fields_in_GEO_SOFT": sorted(technical_fields),
                                   "protocol_descriptors_in_GEO_SOFT": sorted(protocol_descriptors),
                                   "technical_columns_in_supplement": supplement_technical_fields,
                                   "design_matrix_rank": "NOT_ASSESSED"},
        "raw_input_integrity": {"files": len(raw_rows), "bytes": sum(int(row["actual_bytes"]) for row in raw_rows), "all_pass": True},
        "frozen_assets": frozen_hashes,
        "parent_plan_freeze": parent_freeze_rows,
        "blind_results": {"nuclei_total": total_cells, "nuclei_ge_150_genes": nuclei_ge150_genes,
                          "nuclei_ge_150_genes_and_valid_polygon": valid_nuclei,
                          "valid_polygons": valid_polygons, "MG_identity_direction_donors": int(mg_identity),
                          "ASC_identity_direction_donors": int(asc_identity), "G1_supported_donors": g1_mg_donors,
                          "G5_legacy_coverage_counts": g5_summary, "tissue_hole_mask_available": not no_mask,
                          "euclidean_pair_support": euclidean_pair_support, "RCTD_reference_status": reference_status},
        "effects_calculated": False,
        "g6_pass_allowed": False,
        "gate_rows": gates,
    }
    audit_path = root / config["outputs"]["audit_json"]
    audit_path.write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    report = f"""# GSE269906 按最新版 G6 方案 005 的盲态重评

**运行编号：** `{config['run']['id']}`  
**最终状态：** `G6_NOT_EXECUTABLE`  
**分析范围：** BA10 越区技术可行性；未计算 AD–对照效应。

## 结果

用户指定的 GSE269906 已按补充方案 006 重新纳入本轮资格评估。来源清单有 {len(donor_ids)} 个唯一供体、{diagnosis_counts['AD']} 名 AD 和 {diagnosis_counts['Control']} 名作者标注对照，区域均为 BA10 前额叶。原始 CellBin 与 tissue GEF 共 {len(raw_rows)} 个文件、{sum(int(row['actual_bytes']) for row in raw_rows):,} bytes；脚本对每个文件重算 SHA-256，结果全部匹配。冻结 G1/G5 registry 哈希亦与方案 005 一致，方案 005 的 6 个冻结文件也全部通过 checksum 核验。

**该数据不能作为方案 005 的海马主队列。** 方案 005 要求共同海马亚区；GSE269906 测的是 BA10 前额叶。本轮因此是越区的测量可行性重评，不能检验海马估计量，也不能作为独立海马验证。

此前的诊断盲态 CellBin 结果复核为：{nuclei_ge150_genes:,} 个核达到至少 150 个检出基因，其中 {valid_nuclei:,} 个同时有有效 polygon 并进入 QC 分析；独立留出标记的身份方向为 MG {mg_identity}/12 供体、ASC {asc_identity}/12 供体。两类身份都未能校准，MG 尤其没有可独立验证的细胞标签，因此 G1 的目标细胞分数以及 G5 在已验证 ASC 中的部分分数/保真度都不可估。G5 的旧覆盖汇总仅作描述性参考，不能以旧 25% 规则替代新方案要求的实测与分数保真检查。

邻域审计没有可审计的组织孔洞掩膜/组织学配准图像；10、25、50 μm 的记录均不能据此放行任何物理邻域。RCTD 参考供体独立性仍未解决，相关模型输出只作为诊断性试验材料。

## 门槛判定

| 门槛 | 状态 | 结果 |
|---|---|---|
| 供体和原始输入核验 | 通过 | 12 个唯一供体；24/24 原始 GEF 哈希匹配 |
| 母方案冻结完整性 | 通过 | 方案 005 冻结清单中的 6 个文件 SHA-256 全部匹配 |
| 技术批次与设计秩 | 暂缓 | 当前 GEO SOFT/补充表未提供芯片、建库或测序批次 crosswalk；设计秩未评估 |
| 方案 005 脑区范围 | 不通过 | BA10，不是共同海马亚区 |
| MG 身份与 G1 状态 | 未校准 / 不可估 | 独立身份方向 0/12；目标细胞分数及保真度未能估计 |
| ASC 身份与 G5 状态 | 未校准 / 不可估 | 独立身份方向 4/12；目标细胞状态与部分分数保真度不可评估 |
| 空间几何 | 不通过 | 缺少可审计 tissue-hole mask；无半径可放行 |
| 独立同海马验证 | 未找到 | 当前登记中没有合格候选 |
| AD–对照效应 | 未运行 | 受当前方案禁止；没有 P 值、FDR 或组间图 |

## 解释

004 修订取消的是来源级面板覆盖率的硬淘汰，并没有免除身份可验证、目标细胞实际检测、部分评分保真和空间几何要求。本次结论不是 G6 组间阴性，也不是“组织里没有小胶质细胞”；它表示在 BA10 这套现有测量/注释/几何输入下，G6 邻域效应不能被可靠估计。没有任何门槛被降低，也没有改变 G1/G5。

## 执行与复现

- 本轮使用本机已下载的 GSE269906 CellBin/tissue 文件；不联网、不重复下载、不重拟合、不重新聚类。
- 复用的盲态结果：`{config['dataset']['existing_blind_result_dir']}`；重审脚本：`{config['runner']['script']}`；参数：`{config['run']['config']}`。
- 来源级样本映射核对通过；现有 GEO SOFT 与 Supplementary Data1 未报告可用的 chip/run/library/batch crosswalk，设计矩阵秩未评估，记 `HOLD_METADATA`。
- 所有原始 GEF 保持只读。新生成的原始哈希核验、门槛表、JSON 和本报告位于 `{config['outputs']['directory']}`。
- 结果目录由配套 DVC 指针追踪；DVC 状态单独在 Decision Log 中登记，SHA-256 校验不替代 DVC 状态。
"""
    report_path = root / config["outputs"]["report"]
    report_path.write_text(report, encoding="utf-8")

    # Hash the immutable inputs and generated outputs; manifest itself is excluded.
    artifact_paths = [
        config_path,
        root / config["run"]["amendment"],
        root / config["run"]["decision_log"],
        root / config["runner"]["script"],
        parent_freeze_path,
        *[root / row["path"] for row in parent_freeze_rows],
        donor_manifest,
        *[root / value for value in config["dataset"]["source_metadata_files"]],
        root / config["dataset"]["raw_hash_manifest"],
        *[root / asset_cfg["registry"] for asset_cfg in config["frozen_assets"].values()],
        *input_tables,
        raw_hash_path,
        root / config["outputs"]["gates"],
        root / config["outputs"]["source_design_audit"],
        audit_path,
        report_path,
    ]
    manifest_rows = []
    for path in artifact_paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        manifest_rows.append({"path": str(path.relative_to(root)).replace("\\", "/"), "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_tsv(root / config["outputs"]["artifact_manifest"], manifest_rows, ["path", "bytes", "sha256"])
    print(f"Status: {audit['run_status']}")
    print(f"Verified raw GEF files: {len(raw_rows)}/{len(raw_rows)}")
    print(f"Wrote report: {report_path.relative_to(root)}")
    print(f"Wrote gate register: {(root / config['outputs']['gates']).relative_to(root)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

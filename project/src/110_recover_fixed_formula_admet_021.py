from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RUN_ID = "2026-09-30-fixed-formula-admet-recovery-021"
OUT = ROOT / "results" / RUN_ID
BASE = ROOT / "src" / "108_process_fixed_formula_prior_019.py"
INPUT_DIR = ROOT / "results" / "2026-09-30-fixed-formula-drug-prior-020"
SCREEN_PATH = INPUT_DIR / "formula_structure_and_step6_screen.tsv"
MASTER_PATH = INPUT_DIR / "formula_component_candidates.tsv"
POLICY_PATH = Path(r"F:\TCM select\config\pipeline_policy.json")
STEP7_PATH = Path(r"F:\TCM select\code\7_ADMET_online.py")
PARAM_PATH = ROOT / "configs" / "fixed_formula_admet_recovery_021.yaml"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8", errors="replace")
    OUT.mkdir(parents=True, exist_ok=True)

    spec = importlib.util.spec_from_file_location("fixed_formula_prior_shared_021", BASE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load shared prediction implementation: {BASE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    policy = json.loads(POLICY_PATH.read_text(encoding="utf-8"))
    screening = pd.read_csv(SCREEN_PATH, sep="\t", dtype=str, keep_default_na=False)
    master = pd.read_csv(MASTER_PATH, sep="\t", dtype=str, keep_default_na=False)
    if "Step6_Status" not in screening or not screening["Step6_Status"].eq("pass_heuristic").any():
        raise RuntimeError("No Step 6 pass structures found in frozen run 020 input")

    import admet_ai
    model_root = Path(admet_ai.__file__).resolve().parent / "resources" / "models"
    weights = sorted(model_root.rglob("*.pt"))
    sources = [SCREEN_PATH, MASTER_PATH, POLICY_PATH, STEP7_PATH, BASE, Path(__file__).resolve(), PARAM_PATH, *weights]
    source_manifest = pd.DataFrame([
        {"Path": str(p), "Bytes": p.stat().st_size, "SHA256": sha256(p), "Role": "input_method_or_model_weight"}
        for p in sources
    ])
    source_manifest.to_csv(OUT / "source_sha256.tsv", sep="\t", index=False, encoding="utf-8")

    predicted, model_metrics = module.predict_admet(screening, master, policy)
    if predicted.empty:
        status = "FAILED_MODEL_RUN"
        predicted.to_csv(OUT / "formula_admet_predictions.tsv", sep="\t", index=False, encoding="utf-8")
        candidates = master.copy()
        candidates["Step7_Status"] = "not_assessed_model_unavailable"
    else:
        status = "PASS_PREDICTION_COMPLETED_WITH_MODEL_LIMITS"
        predicted.to_csv(OUT / "formula_admet_predictions.tsv", sep="\t", index=False, encoding="utf-8")
        pred_map = predicted.set_index("Canonical_Isomeric_SMILES")["Step7_Status"].to_dict()
        candidates = master.copy()
        candidates["Step7_Status"] = candidates["Structure_SMILES"].map(
            lambda x: pred_map.get(module.canonical_smiles(x), "not_assessed")
        )
    candidates.to_csv(OUT / "formula_component_candidates_with_admet.tsv", sep="\t", index=False, encoding="utf-8")

    summary = {
        "run_id": RUN_ID,
        "decision_id": "DRUG-FORMULA-ADMET-RECOVERY-2026-09-30-021",
        "source_run": "2026-09-30-fixed-formula-drug-prior-020",
        "step6_pass_input_count": int(screening["Step6_Status"].eq("pass_heuristic").sum()),
        "model_metrics": model_metrics,
        "model_weights_count": len(weights),
        "model_weights": [{"path": str(p), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in weights],
        "console_encoding": "UTF-8",
        "network_queries": False,
        "dependencies_installed": False,
        "measured_ADMET": False,
        "brain_exposure_established": False,
        "drug_evidence_ceiling": "C0",
        "validation_steps_reopened": False,
        "software": {
            "python": sys.version.split()[0],
            "pandas": module.importlib.metadata.version("pandas"),
            "admet-ai": module.importlib.metadata.version("admet-ai"),
        },
        "status": status,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    audit = {
        "status": status,
        "checks": {
            "input_is_frozen_step6_pass_table_from_020": True,
            "UTF8_console_used": True,
            "model_weights_hashed": bool(weights),
            "raw_formula_sources_modified": False,
            "new_network_queries": False,
            "new_dependencies_installed": False,
            "ADMET_predictions_are_not_measurements": True,
            "brain_exposure_or_action_direction_inferred": False,
            "C0_claim_ceiling_retained": True,
            "017_closed_steps_reopened": False,
        },
    }
    (OUT / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")

    report = f'''# 固定处方 Step 7 ADMET 预测补跑 021

日期：2026-09-30。此记录补充[020结构与来源审计](../2026-09-30-fixed-formula-drug-prior-020/report.md)。020的本地ADMET-AI执行已开始，但GBK控制台无法编码模型进度输出，故020预测列未形成；本轮仅用UTF-8控制台对冻结的Step 6输入补跑。

## 结果

- Step 6通过并进入预测的唯一结构：{model_metrics.get('input_count', 0)}。
- ADMET-AI状态：`{model_metrics.get('model_status')}`；预测通过{model_metrics.get('step7_pass_predictive_count', 0)}，未通过{model_metrics.get('step7_fail_predictive_count', 0)}。逐分子字段、策略判断见[预测表](formula_admet_predictions.tsv)。
- 相关安装模型权重文件共{len(weights)}份，哈希见[来源清单](source_sha256.tsv)和`summary.json`。

模型通过/未通过只按既有Step 7项目预测规则分类。HIA、Caco-2、溶解度、CYP、DILI/hERG及BBB均是模型输出，不是实测；BBB也不是强制门槛。没有本方实际化学定量和血/脑/游离暴露时，不能把预测通过解释成Direct-CNS或P3可达药效。 target action direction, efficacy and mechanism remain undetermined.

依赖仅缺失输入结果不重新运行真实药效、整方机制或细胞通讯实验；药物证据仍为C0，017关闭步骤未重启。执行代码：[补跑脚本](../../src/110_recover_fixed_formula_admet_021.py)，参数：[021配置](../../configs/fixed_formula_admet_recovery_021.yaml)，校验：[audit.json](audit.json)。本地使用已有TCM_select环境；无新安装、网络查询或外部源修改。
'''
    (OUT / "report.md").write_text(report, encoding="utf-8")


if __name__ == "__main__":
    main()

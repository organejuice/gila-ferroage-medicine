#!/usr/bin/env python3
"""Export G6 signatures from frozen G1/G5 outputs without using G6 effects."""

from __future__ import annotations

import gzip
import hashlib
import json
from datetime import date
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/2026-09-24-disease-side-v3.5-g6-input-audit"
REG = OUT / "registries"
TABLES = OUT / "tables"
REG.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

G1 = ROOT / "results/archive/2026-09-20-disease-side-v3.1-g2-stop-historical/gates/G1_gene_signal.tsv"
G5_CONFIG = ROOT / "configs/g5_v3.4_execution.yaml"
G5_RESULT = ROOT / "results/2026-09-22-disease-side-v3.4-g5/tables/GSE261698_primary_module_effects.tsv"
M7_REFERENCE = ROOT / "results/archive/2026-09-20-pre-v3-freeze/results/tables/D5_ferroaging_reference.tsv"
MOUSE_MAP = ROOT / "results/2026-09-22-disease-side-v3.4-g5/registries/GSE261698_gene_mapping.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_gzip_tsv(frame: pd.DataFrame, path: Path) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        frame.to_csv(handle, sep="\t", index=False, na_rep="")


# G1's v3.1 audit defines strict human hits as meta-FDR < 0.05 and k == 5.
g1 = pd.read_csv(G1, sep="\t")
for required in ("gene", "k", "meta_logFC", "meta_FDR"):
    if required not in g1.columns:
        raise RuntimeError(f"Frozen G1 output is missing {required}.")
g1["k"] = pd.to_numeric(g1["k"], errors="coerce")
g1["meta_FDR"] = pd.to_numeric(g1["meta_FDR"], errors="coerce")
g1["meta_logFC"] = pd.to_numeric(g1["meta_logFC"], errors="coerce")
sender = g1.loc[(g1["k"] == 5) & (g1["meta_FDR"] < 0.05)].copy()
if len(sender) != 244 or sender["gene"].duplicated().any():
    raise RuntimeError(f"Expected 244 unique strict G1 senders, got {len(sender)}.")
sender["direction"] = sender["meta_logFC"].map(lambda value: "positive_AD" if value > 0 else "negative_AD" if value < 0 else "zero")
sender["source_rule"] = "G1 v3.1 strict hit: meta_FDR<0.05 and k==5; direction=meta_logFC sign"
write_gzip_tsv(sender, REG / "G1_microglia_signed_sender_244.tsv.gz")

# Resolve the seven modules from the immutable v3.4 execution config. M7 is
# expanded from its already-frozen author reference list.
config = yaml.safe_load(G5_CONFIG.read_text(encoding="utf-8"))
module_defs = config["analysis"]["modules"]
g5_results = pd.read_csv(G5_RESULT, sep="\t")
result_by_module = g5_results.loc[g5_results["analysis"] == "primary"].set_index("module").to_dict(orient="index")
mouse_map = pd.read_csv(MOUSE_MAP, sep="\t", dtype=str).fillna("")
mouse_by_symbol: dict[str, list[str]] = {}
for row in mouse_map.itertuples(index=False):
    symbol = getattr(row, "gene_symbol")
    raw_id = getattr(row, "raw_id")
    if symbol and raw_id:
        mouse_by_symbol.setdefault(symbol.upper(), []).append(raw_id)

module_rows: list[dict] = []
for module, definition in module_defs.items():
    if isinstance(definition, str) and definition.startswith("reference:"):
        reference_module = definition.split(":", 1)[1]
        if reference_module != "m7_reference":
            raise RuntimeError(f"Unexpected unresolved G5 module reference: {definition}")
        ref = pd.read_csv(M7_REFERENCE, sep="\t", dtype=str).fillna("")
        genes = sorted(set(ref.loc[ref["module"] == "M7", "gene_symbol"].str.upper()) - {""})
        gene_source = M7_REFERENCE.relative_to(ROOT).as_posix()
    else:
        genes = [str(gene).upper() for gene in definition]
        gene_source = G5_CONFIG.relative_to(ROOT).as_posix()
    if module not in result_by_module:
        raise RuntimeError(f"No frozen primary G5 result for module {module}.")
    effect = result_by_module[module]
    status = str(effect["status"])
    logfc = pd.to_numeric(pd.Series([effect["logFC"]]), errors="coerce").iloc[0]
    fdr = pd.to_numeric(pd.Series([effect["FDR_BH_7"]]), errors="coerce").iloc[0]
    if status == "non_estimable" or pd.isna(logfc):
        gate = "NON_ESTIMABLE"
        direction = "not_estimable"
    else:
        gate = "PASS" if pd.notna(fdr) and fdr < 0.05 else "NO_PASS"
        direction = "positive_model" if logfc > 0 else "negative_model" if logfc < 0 else "zero_model"
    for gene in genes:
        module_rows.append({
            "module": module,
            "gene_symbol_in_frozen_G5_list": gene,
            "gene_list_source": gene_source,
            "G5_primary_status": status,
            "G5_module_gate": gate,
            "G5_module_logFC": logfc,
            "G5_module_FDR_BH_7": fdr,
            "G5_module_direction_hypothesis": direction,
            "symbol_seen_in_GSE261698_mouse_model_mapping": gene in mouse_by_symbol,
            "GSE261698_mouse_raw_feature_ids_for_symbol": ";".join(sorted(set(mouse_by_symbol.get(gene, [])))),
            "human_one_to_one_ortholog_status": "PENDING_VERSIONED_ORTHOLOGY_AUDIT",
            "note": "Module direction is a gene-set-level model hypothesis; it is not an individual-gene direction or human G6 result.",
        })
modules = pd.DataFrame(module_rows)
if modules["module"].nunique() != 7:
    raise RuntimeError("The frozen G5 receiver family must contain exactly seven modules.")
if modules.loc[modules["module"] == "M5_Iron_Storage", "G5_module_gate"].nunique() != 1 or modules.loc[modules["module"] == "M5_Iron_Storage", "G5_module_gate"].iloc[0] != "NON_ESTIMABLE":
    raise RuntimeError("M5_Iron_Storage must remain non-estimable and cannot be promoted.")
write_gzip_tsv(modules, REG / "G5_frozen_receiver_module_gene_registry.tsv.gz")

module_summary = modules.drop_duplicates("module")[[
    "module", "G5_primary_status", "G5_module_gate", "G5_module_logFC",
    "G5_module_FDR_BH_7", "G5_module_direction_hypothesis"
]].to_dict(orient="records")
provenance = {
    "run_id": "2026-09-24-disease-side-v3.5-g6-input-audit",
    "generated_date": date.today().isoformat(),
    "G6_effects_used": False,
    "G1_sender": {
        "source": G1.relative_to(ROOT).as_posix(),
        "source_sha256": sha256(G1),
        "rule": "meta_FDR < 0.05 and k == 5, inherited from src/27_audit_disease_side_v3_1.py",
        "n_genes": len(sender),
        "directions": {str(key): int(value) for key, value in sender["direction"].value_counts().items()},
        "output": "registries/G1_microglia_signed_sender_244.tsv.gz",
    },
    "G5_receiver": {
        "config": G5_CONFIG.relative_to(ROOT).as_posix(),
        "config_sha256": sha256(G5_CONFIG),
        "primary_result": G5_RESULT.relative_to(ROOT).as_posix(),
        "primary_result_sha256": sha256(G5_RESULT),
        "M7_reference": M7_REFERENCE.relative_to(ROOT).as_posix(),
        "M7_reference_sha256": sha256(M7_REFERENCE),
        "gene_mapping": MOUSE_MAP.relative_to(ROOT).as_posix(),
        "gene_mapping_sha256": sha256(MOUSE_MAP),
        "module_count": int(modules["module"].nunique()),
        "gene_rows": len(modules),
        "modules": module_summary,
        "one_to_one_human_ortholog_mapping": "not emitted; a versioned one-to-one orthology source is required",
        "output": "registries/G5_frozen_receiver_module_gene_registry.tsv.gz",
    },
}
(REG / "G6_frozen_signature_provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps({"G1_signed_sender_genes": len(sender), "G5_modules": modules["module"].nunique(), "G5_gene_rows": len(modules), "G5_one_to_one_human_orthologs": "pending", "G6_effects_used": False}, ensure_ascii=False))

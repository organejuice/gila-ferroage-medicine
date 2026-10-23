#!/usr/bin/env python3
"""Map frozen G5 receiver genes via a dated NCBI mouse-human ortholog snapshot."""

from __future__ import annotations

import gzip
import hashlib
import json
import re
from collections import defaultdict
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/NCBI_gene_orthologs"
REG = ROOT / "results/2026-09-24-disease-side-v3.5-g6-input-audit/registries"
TABLES = ROOT / "results/2026-09-24-disease-side-v3.5-g6-input-audit/tables"
ORTHOLOG_FILE = RAW / "gene_orthologs.gz"
MOUSE_INFO_FILE = RAW / "Mus_musculus.gene_info.gz"
HUMAN_INFO_FILE = RAW / "Homo_sapiens.gene_info.gz"
G5_FILE = REG / "G5_frozen_receiver_module_gene_registry.tsv.gz"
G5_MOUSE_MAP = ROOT / "results/2026-09-22-disease-side-v3.4-g5/registries/GSE261698_gene_mapping.tsv"
PROVENANCE_FILE = REG / "G6_frozen_signature_provenance.json"
OUT = REG / "G5_receiver_one_to_one_human_ortholog_registry.tsv.gz"

for path in (ORTHOLOG_FILE, MOUSE_INFO_FILE, HUMAN_INFO_FILE, G5_FILE, G5_MOUSE_MAP):
    if not path.exists():
        raise FileNotFoundError(path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalize_col(name: str) -> str:
    name = name.lstrip("#").strip().lower()
    return re.sub(r"[^a-z0-9]+", "_", name).strip("_")


def read_ncbi_gene_info(path: Path, tax_id: str) -> pd.DataFrame:
    info = pd.read_csv(path, sep="\t", dtype=str, compression="gzip", low_memory=False)
    info.columns = [normalize_col(col) for col in info.columns]
    required = {"tax_id", "geneid", "symbol", "synonyms"}
    if not required.issubset(info.columns):
        raise RuntimeError(f"Unexpected NCBI gene_info columns in {path.name}: {info.columns.tolist()}")
    info = info.loc[info["tax_id"] == tax_id, ["tax_id", "geneid", "symbol", "synonyms"]].copy()
    info = info.drop_duplicates(["geneid"])
    return info


mouse_info = read_ncbi_gene_info(MOUSE_INFO_FILE, "10090")
human_info = read_ncbi_gene_info(HUMAN_INFO_FILE, "9606")
mouse_symbols = dict(zip(mouse_info["geneid"], mouse_info["symbol"]))
human_symbols = dict(zip(human_info["geneid"], human_info["symbol"]))

# NCBI records ortholog pairs as TaxID/GeneID and Other TaxID/GeneID. Keep
# only direct mouse-human records, then independently require reciprocal
# uniqueness so no ambiguous orthology group is labeled one-to-one.
pairs: set[tuple[str, str]] = set()
for chunk in pd.read_csv(ORTHOLOG_FILE, sep="\t", dtype=str, compression="gzip", chunksize=500_000, low_memory=False):
    original_columns = list(chunk.columns)
    chunk.columns = [normalize_col(col) for col in chunk.columns]
    required = {"tax_id", "geneid", "relationship", "other_tax_id", "other_geneid"}
    if not required.issubset(chunk.columns):
        raise RuntimeError(f"Unexpected NCBI gene_orthologs columns: {original_columns}")
    direct = chunk[
        ((chunk["tax_id"] == "10090") & (chunk["other_tax_id"] == "9606"))
        | ((chunk["tax_id"] == "9606") & (chunk["other_tax_id"] == "10090"))
    ]
    for row in direct.itertuples(index=False):
        if row.relationship != "Ortholog":
            continue
        if row.tax_id == "10090":
            mouse_id, human_id = row.geneid, row.other_geneid
        else:
            mouse_id, human_id = row.other_geneid, row.geneid
        pairs.add((str(mouse_id), str(human_id)))

mouse_to_human: dict[str, set[str]] = defaultdict(set)
human_to_mouse: dict[str, set[str]] = defaultdict(set)
for mouse_id, human_id in pairs:
    mouse_to_human[mouse_id].add(human_id)
    human_to_mouse[human_id].add(mouse_id)
accepted_pairs = {
    (mouse_id, next(iter(humans)))
    for mouse_id, humans in mouse_to_human.items()
    if len(humans) == 1 and len(human_to_mouse[next(iter(humans))]) == 1
}

# Resolve each frozen G5 symbol first against the actual GSE261698 feature map,
# then against NCBI's official mouse symbol/synonym fields for genes absent
# from that animal matrix. Ambiguity is retained as a non-mapping status.
g5 = pd.read_csv(G5_FILE, sep="\t", compression="gzip", dtype=str).fillna("")
model_map = pd.read_csv(G5_MOUSE_MAP, sep="\t", dtype=str).fillna("")
model_map["normalized_symbol"] = model_map["gene_symbol"].str.upper()
model_map["mouse_geneid"] = model_map["raw_id"].str.extract(r"_(\d+)$", expand=False).fillna("")
feature_ids_by_symbol: dict[str, set[str]] = defaultdict(set)
for row in model_map.itertuples(index=False):
    if row.normalized_symbol and row.mouse_geneid:
        feature_ids_by_symbol[row.normalized_symbol].add(row.mouse_geneid)

mouse_aliases: dict[str, set[str]] = defaultdict(set)
for row in mouse_info.itertuples(index=False):
    gene_id = str(row.geneid)
    symbol = str(row.symbol)
    if symbol and symbol != "-":
        mouse_aliases[symbol.upper()].add(gene_id)
    synonyms = str(row.synonyms)
    if synonyms and synonyms != "-":
        for alias in synonyms.split("|"):
            if alias and alias != "-":
                mouse_aliases[alias.upper()].add(gene_id)

rows: list[dict] = []
for gene in g5.itertuples(index=False):
    symbol = gene.gene_symbol_in_frozen_G5_list.upper()
    feature_ids = feature_ids_by_symbol.get(symbol, set())
    if len(feature_ids) == 1:
        candidate_mouse_ids = set(feature_ids)
        mouse_resolution = "GSE261698 feature Entrez suffix via exact normalized symbol"
    elif len(feature_ids) > 1:
        candidate_mouse_ids = set(feature_ids)
        mouse_resolution = "ambiguous GSE261698 feature mapping; no single mouse GeneID selected"
    else:
        candidate_mouse_ids = set(mouse_aliases.get(symbol, set()))
        mouse_resolution = "NCBI mouse official symbol/synonym" if candidate_mouse_ids else "no NCBI mouse symbol/synonym match"

    status = "no_mouse_geneid_match"
    mapped_mouse_id = ""
    human_candidates: set[str] = set()
    accepted_human_id = ""
    if len(candidate_mouse_ids) > 1:
        status = "ambiguous_mouse_geneid_match"
    elif len(candidate_mouse_ids) == 1:
        mapped_mouse_id = next(iter(candidate_mouse_ids))
        human_candidates = mouse_to_human.get(mapped_mouse_id, set())
        if len(human_candidates) > 1:
            status = "not_one_to_one_multiple_human_orthologs"
        elif len(human_candidates) == 0:
            status = "no_NCBI_mouse_human_ortholog_record"
        else:
            candidate_human_id = next(iter(human_candidates))
            if (mapped_mouse_id, candidate_human_id) not in accepted_pairs:
                status = "not_one_to_one_reciprocal"
            elif candidate_human_id not in human_symbols:
                status = "human_gene_info_missing"
            else:
                status = "accepted_NCBI_reciprocal_one_to_one"
                accepted_human_id = candidate_human_id

    rows.append({
        "module": gene.module,
        "frozen_G5_gene_symbol": symbol,
        "G5_gene_list_source": gene.gene_list_source,
        "G5_primary_status": gene.G5_primary_status,
        "G5_module_gate": gene.G5_module_gate,
        "G5_module_logFC": gene.G5_module_logFC,
        "G5_module_FDR_BH_7": gene.G5_module_FDR_BH_7,
        "G5_module_direction_hypothesis": gene.G5_module_direction_hypothesis,
        "mouse_symbol_resolution": mouse_resolution,
        "candidate_mouse_GeneIDs": ";".join(sorted(candidate_mouse_ids)),
        "mouse_official_symbols": ";".join(sorted({mouse_symbols.get(item, "") for item in candidate_mouse_ids if mouse_symbols.get(item, "")})),
        "human_ortholog_GeneIDs": ";".join(sorted(human_candidates)),
        "one_to_one_mapping_status": status,
        "human_one_to_one_GeneID": accepted_human_id,
        "human_one_to_one_symbol": human_symbols.get(accepted_human_id, ""),
        "GSE261698_feature_symbol_seen": gene.symbol_seen_in_GSE261698_mouse_model_mapping == "True",
        "note": "G5 direction is module-level animal-model evidence; it is not a gene-level direction or a G6 result.",
    })

mapped = pd.DataFrame(rows)
mapped["duplicate_human_gene_within_module"] = False
valid = mapped["one_to_one_mapping_status"] == "accepted_NCBI_reciprocal_one_to_one"
mapped.loc[valid, "duplicate_human_gene_within_module"] = mapped.loc[valid].duplicated(
    ["module", "human_one_to_one_GeneID"], keep=False
)
with gzip.open(OUT, "wt", encoding="utf-8", newline="") as handle:
    mapped.to_csv(handle, sep="\t", index=False, na_rep="")

updated_g5 = g5.merge(
    mapped[[
        "module", "frozen_G5_gene_symbol", "one_to_one_mapping_status",
        "human_one_to_one_GeneID", "human_one_to_one_symbol"
    ]],
    left_on=["module", "gene_symbol_in_frozen_G5_list"],
    right_on=["module", "frozen_G5_gene_symbol"],
    how="left",
    validate="one_to_one",
)
updated_g5["human_one_to_one_ortholog_status"] = updated_g5["one_to_one_mapping_status"]
updated_g5["human_one_to_one_ortholog_GeneID"] = updated_g5["human_one_to_one_GeneID"]
updated_g5["human_one_to_one_ortholog_symbol"] = updated_g5["human_one_to_one_symbol"]
updated_g5 = updated_g5.drop(columns=["frozen_G5_gene_symbol", "one_to_one_mapping_status", "human_one_to_one_GeneID", "human_one_to_one_symbol"])
with gzip.open(G5_FILE, "wt", encoding="utf-8", newline="") as handle:
    updated_g5.to_csv(handle, sep="\t", index=False, na_rep="")

module_summary = []
for module, group in mapped.groupby("module", sort=False):
    accepted = group.loc[group["one_to_one_mapping_status"] == "accepted_NCBI_reciprocal_one_to_one"]
    module_summary.append({
        "module": module,
        "G5_module_gate": group["G5_module_gate"].iloc[0],
        "G5_module_direction_hypothesis": group["G5_module_direction_hypothesis"].iloc[0],
        "n_frozen_input_genes": len(group),
        "n_NCBI_reciprocal_one_to_one_rows": len(accepted),
        "n_unique_human_GeneIDs_after_deduplication": accepted["human_one_to_one_GeneID"].nunique(),
        "coverage_one_to_one": accepted["human_one_to_one_GeneID"].nunique() / len(group),
        "n_missing_or_ambiguous": int((group["one_to_one_mapping_status"] != "accepted_NCBI_reciprocal_one_to_one").sum()),
    })
pd.DataFrame(module_summary).to_csv(TABLES / "G5_receiver_one_to_one_human_module_coverage.tsv", sep="\t", index=False)

provenance = {
    "mapping_source": "NCBI Gene FTP current snapshot",
    "retrieval_date": date.today().isoformat(),
    "NCBI_gene_orthologs": {
        "path": ORTHOLOG_FILE.relative_to(ROOT).as_posix(),
        "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz",
        "bytes": ORTHOLOG_FILE.stat().st_size,
        "sha256": sha256(ORTHOLOG_FILE),
        "source_last_modified_utc": "2026-09-23 (HTTP Last-Modified)",
        "method": "tax_id 10090 vs 9606; relationship Ortholog; reciprocal unique mouse GeneID and human GeneID",
        "n_mouse_human_pair_rows": len(pairs),
        "n_reciprocal_one_to_one_pairs": len(accepted_pairs),
    },
    "mouse_gene_info": {
        "path": MOUSE_INFO_FILE.relative_to(ROOT).as_posix(),
        "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
        "bytes": MOUSE_INFO_FILE.stat().st_size,
        "sha256": sha256(MOUSE_INFO_FILE),
        "source_last_modified_utc": "2026-09-23 (HTTP Last-Modified)",
    },
    "human_gene_info": {
        "path": HUMAN_INFO_FILE.relative_to(ROOT).as_posix(),
        "url": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
        "bytes": HUMAN_INFO_FILE.stat().st_size,
        "sha256": sha256(HUMAN_INFO_FILE),
        "source_last_modified_utc": "2026-09-23 (HTTP Last-Modified)",
    },
    "G5_registry": G5_FILE.relative_to(ROOT).as_posix(),
    "G5_registry_sha256": sha256(G5_FILE),
    "output": OUT.relative_to(ROOT).as_posix(),
    "output_sha256": sha256(OUT),
    "G6_effects_used": False,
    "modules": module_summary,
}
(REG / "G6_NCBI_one_to_one_orthology_provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

source_manifest = pd.DataFrame([
    {"accession": "NCBI_Gene", "source_url": provenance["NCBI_gene_orthologs"]["url"], "download_date": date.today().isoformat(), "source_last_modified_utc": "2026-09-23", "file": provenance["NCBI_gene_orthologs"]["path"], "bytes": provenance["NCBI_gene_orthologs"]["bytes"], "sha256": provenance["NCBI_gene_orthologs"]["sha256"], "processing_level": "NCBI mouse-human ortholog pairs; locally filtered reciprocal one-to-one", "license_access_conditions": "Public NCBI Gene FTP; NCBI terms apply"},
    {"accession": "NCBI_Gene", "source_url": provenance["mouse_gene_info"]["url"], "download_date": date.today().isoformat(), "source_last_modified_utc": "2026-09-23", "file": provenance["mouse_gene_info"]["path"], "bytes": provenance["mouse_gene_info"]["bytes"], "sha256": provenance["mouse_gene_info"]["sha256"], "processing_level": "NCBI official mouse symbol and synonym annotation", "license_access_conditions": "Public NCBI Gene FTP; NCBI terms apply"},
    {"accession": "NCBI_Gene", "source_url": provenance["human_gene_info"]["url"], "download_date": date.today().isoformat(), "source_last_modified_utc": "2026-09-23", "file": provenance["human_gene_info"]["path"], "bytes": provenance["human_gene_info"]["bytes"], "sha256": provenance["human_gene_info"]["sha256"], "processing_level": "NCBI official human symbol annotation", "license_access_conditions": "Public NCBI Gene FTP; NCBI terms apply"},
])
source_manifest.to_csv(REG / "G6_NCBI_orthology_source_manifest.tsv", sep="\t", index=False)

signature_provenance = json.loads(PROVENANCE_FILE.read_text(encoding="utf-8"))
signature_provenance["G5_receiver"]["one_to_one_human_ortholog_mapping"] = {
    "status": "audited_with_NCBI_Gene_FTP_snapshot",
    "source": "registries/G5_receiver_one_to_one_human_ortholog_registry.tsv.gz",
    "sha256": provenance["output_sha256"],
    "mapping_provenance": "registries/G6_NCBI_one_to_one_orthology_provenance.json",
}
PROVENANCE_FILE.write_text(json.dumps(signature_provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
print(json.dumps({
    "NCBI_mouse_human_pair_rows": len(pairs),
    "NCBI_reciprocal_one_to_one_pairs": len(accepted_pairs),
    "G5_input_gene_rows": len(mapped),
    "G5_accepted_one_to_one_rows": int(valid.sum()),
    "G5_unique_human_gene_ids": int(mapped.loc[valid, "human_one_to_one_GeneID"].nunique()),
    "G6_effects_used": False,
}, ensure_ascii=False))

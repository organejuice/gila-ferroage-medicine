from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = ROOT / "configs" / "fixed_formula_g1_overlap_scope_023.json"
SCRIPT_PATH = Path(__file__).resolve()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    config = json.loads(CONFIG_PATH.read_text(encoding="utf-8"))
    run_id = config["run_id"]
    output_dir = ROOT / "results" / run_id
    output_dir.mkdir(parents=True, exist_ok=True)
    inputs = {key: ROOT / config[key] for key in (
        "candidate_table", "broad_prior_table", "strict_assay_table",
        "g1_robust_table", "g1_full_meta_table", "previous_summary"
    )}
    inputs["config"] = CONFIG_PATH
    inputs["script"] = SCRIPT_PATH
    hashes_before = {key: sha256(path) for key, path in inputs.items()}

    candidates = read_tsv(inputs["candidate_table"])
    priors = read_tsv(inputs["broad_prior_table"])
    strict_records = read_tsv(inputs["strict_assay_table"])
    g1_rows = read_tsv(inputs["g1_robust_table"])
    full_g1_meta = read_tsv(inputs["g1_full_meta_table"])
    old_summary = json.loads(inputs["previous_summary"].read_text(encoding="utf-8"))

    candidate_by_cid = {row["PubChem_CID"].strip(): row for row in candidates if row["PubChem_CID"].strip()}
    candidate_cids = set(candidate_by_cid)
    g1_by_gene = {row["gene"].strip(): row for row in g1_rows if row["gene"].strip()}
    full_meta_by_gene = {row["gene"].strip(): row for row in full_g1_meta if row["gene"].strip()}
    matching_full_meta = sum(
        gene in full_meta_by_gene
        and abs(float(g1["meta_logFC"]) - float(full_meta_by_gene[gene]["meta_logFC"])) < 1e-10
        and abs(float(g1["meta_FDR"]) - float(full_meta_by_gene[gene]["meta_FDR"])) < 1e-10
        for gene, g1 in g1_by_gene.items()
    )
    prior_rows = [row for row in priors if row["Target_Species"].strip() == "Homo sapiens"]
    prior_cids = {row["PubChem_CID"].strip() for row in prior_rows}
    prior_symbols = {row["Target_Symbol"].strip() for row in prior_rows if row["Target_Symbol"].strip()}
    strict_cids = {row["PubChem_CID"].strip() for row in strict_records}
    strict_symbols = {row["Target_Symbol"].strip() for row in strict_records if row["Target_Symbol"].strip()}
    broad_overlap = sorted(prior_symbols & set(g1_by_gene))
    strict_overlap = sorted(strict_symbols & set(g1_by_gene))

    prior_by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in prior_rows:
        if row["Target_Symbol"].strip() in broad_overlap:
            prior_by_gene[row["Target_Symbol"].strip()].append(row)
    detailed_rows = []
    for gene in broad_overlap:
        matching = prior_by_gene[gene]
        cids = sorted({row["PubChem_CID"].strip() for row in matching}, key=int)
        source_tags = sorted({tag.strip() for row in matching for tag in row["Target_DB_Name"].split(";") if tag.strip()})
        g1 = g1_by_gene[gene]
        detailed_rows.append({
            "gene": gene,
            "G1_meta_logFC": g1["meta_logFC"],
            "G1_meta_FDR": g1["meta_FDR"],
            "G1_direction": g1["direction"],
            "candidate_CIDs": ";".join(cids),
            "candidate_names": ";".join(candidate_by_cid[cid]["Candidate_Name"] for cid in cids),
            "Step6_statuses": ";".join(candidate_by_cid[cid]["Step6_Status"] for cid in cids),
            "Step7_statuses": ";".join(candidate_by_cid[cid]["Step7_Status"] for cid in cids),
            "imported_source_types": ";".join(source_tags),
            "prior_table_rows": str(len(matching)),
            "evidence_tier": "imported_mapping_prior_not_qualified_direct_assay",
        })
    fields = list(detailed_rows[0]) if detailed_rows else [
        "gene", "G1_meta_logFC", "G1_meta_FDR", "G1_direction", "candidate_CIDs",
        "candidate_names", "Step6_statuses", "Step7_statuses", "imported_source_types",
        "prior_table_rows", "evidence_tier"
    ]
    write_tsv(output_dir / "broad_prior_G1_overlap.tsv", detailed_rows, fields)

    hub_rows = []
    for symbol in config["published_abstract_hub_symbols"] + sorted(strict_symbols):
        meta = full_meta_by_gene.get(symbol, {})
        hub_rows.append({
            "gene": symbol,
            "source_role": "published_abstract_hub" if symbol in config["published_abstract_hub_symbols"] else "022_strict_assay_target",
            "G1_meta_k": meta.get("k", ""),
            "G1_meta_logFC": meta.get("meta_logFC", ""),
            "G1_meta_FDR": meta.get("meta_FDR", ""),
            "in_G1_244": str(symbol in g1_by_gene).lower(),
        })
    write_tsv(output_dir / "published_hubs_and_strict_target_G1_status.tsv", hub_rows,
              ["gene", "source_role", "G1_meta_k", "G1_meta_logFC", "G1_meta_FDR", "in_G1_244"])

    counts = {
        "g1_genes": len(g1_by_gene),
        "g1_full_meta_genes": len(full_meta_by_gene),
        "g1_244_full_meta_matches": matching_full_meta,
        "candidate_cids": len(candidate_cids),
        "prior_rows": len(prior_rows),
        "prior_cids": len(prior_cids),
        "prior_target_symbols": len(prior_symbols),
        "broad_overlap_genes": len(broad_overlap),
        "published_hubs_in_g1_244": len(set(config["published_abstract_hub_symbols"]) & set(g1_by_gene)),
        "strict_assay_rows": len(strict_records),
        "strict_overlap_genes": len(strict_overlap),
    }
    checks = {key: counts[key] == expected for key, expected in config["expected"].items()}
    checks.update({
        "g1_gene_ids_unique": len(g1_by_gene) == len(g1_rows),
        "prior_cids_in_formula_candidates": prior_cids <= candidate_cids,
        "strict_cids_in_formula_candidates": strict_cids <= candidate_cids,
        "strict_target_present_in_prior_symbols": strict_symbols <= prior_symbols,
        "prior_overlap_rows_unique": len({row["gene"] for row in detailed_rows}) == len(detailed_rows),
        "previous_022_strict_overlap_reproduced": old_summary["candidate_compound_G1_overlap"]["exact_target_symbol_overlap_count"] == len(strict_overlap),
        "inputs_unchanged": hashes_before == {key: sha256(path) for key, path in inputs.items()},
    })
    summary = {
        "run_id": run_id,
        "scope": "descriptive comparison of two evidence tiers against the frozen G1 244 set",
        "counts": counts,
        "broad_prior_G1_overlap_genes": broad_overlap,
        "strict_assay_G1_overlap_genes": strict_overlap,
        "prior_CIDs_with_any_target_record": sorted(prior_cids, key=int),
        "prior_CIDs_in_G1_overlap": sorted({cid for row in detailed_rows for cid in row["candidate_CIDs"].split(";")}, key=int),
        "strict_target_symbols": sorted(strict_symbols),
        "published_abstract": config["published_abstract"],
        "status": "PASS_TWO_TIER_SCOPE_AUDIT" if all(checks.values()) else "FAIL_AUDIT",
        "interpretation": "8 imported target-mapping hypotheses overlap G1; the one exact-structure BindingDB record does not. Neither count establishes formula efficacy or target action direction.",
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    source_rows = [{"role": key, "path": str(path.relative_to(ROOT)), "bytes": str(path.stat().st_size), "sha256": hashes_before[key]} for key, path in inputs.items()]
    write_tsv(output_dir / "source_sha256.tsv", source_rows, ["role", "path", "bytes", "sha256"])
    audit = {"run_id": run_id, "status": summary["status"], "checks": checks, "all_pass": all(checks.values())}
    (output_dir / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    if not audit["all_pass"]:
        raise SystemExit("Audit failed; inspect audit.json")
    print(json.dumps({"status": summary["status"], **counts, "broad_overlap": broad_overlap}, ensure_ascii=False))


if __name__ == "__main__":
    main()

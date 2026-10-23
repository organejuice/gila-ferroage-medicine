"""Effect-blind audit of the publicly obtainable Avey section metadata.

Run from the project root: D:/python/python.exe src/102_audit_g6_avey_access.py
No Synapse credentials or expression matrices are read.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OLDER = ROOT / "data/raw/G6_source_metadata_2026-09-29"
NEWER = ROOT / "data/raw/G6_avey_access_probe_2026-09-29"
DESIGN = ROOT / "results/2026-09-29-g6-avey-dharshini-metadata-audit/avey_donor_design.tsv"
OUT = ROOT / "results/2026-09-29-g6-avey-access"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    required = [
        OLDER / "avey_supplementary_file9.csv",
        NEWER / "avey_author_BayesSpace_script.R",
        DESIGN,
    ]
    missing = [str(p) for p in required if not p.is_file()]
    if missing:
        raise FileNotFoundError("Missing required inputs: " + ", ".join(missing))

    with DESIGN.open(encoding="utf-8-sig", newline="") as stream:
        donor_rows = list(csv.DictReader(stream, delimiter="\t"))
    donor_ids = {int(row["donor_id"]) for row in donor_rows}
    if len(donor_ids) != 21:
        raise ValueError(f"Expected 21 unique donors, found {len(donor_ids)}")

    with required[0].open(encoding="utf-8-sig", newline="") as stream:
        spot_rows = list(csv.DictReader(stream))
    spot_sections = {(int(row["ID"]), row["section"].strip()) for row in spot_rows}
    code = required[1].read_text(encoding="utf-8")
    code_sections = {
        (int(donor), section + "1")
        for donor, section in re.findall(
            r'colnames\(Sample(\d+)([ABCD])\)\s*<-\s*paste0', code
        )
    }
    if not (spot_sections | code_sections) <= {
        (donor, section)
        for donor in donor_ids
        for section in ("A1", "B1", "C1", "D1")
    }:
        raise ValueError("Unexpected donor or section label in public source")

    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "partial_donor_section_crosswalk.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(
            [
                "donor_id",
                "section_label",
                "supplementary_file9_spot_evidence",
                "bayesspace_script_evidence",
                "evidence_scope",
            ]
        )
        for donor, section in sorted(spot_sections | code_sections):
            writer.writerow(
                [
                    donor,
                    section,
                    int((donor, section) in spot_sections),
                    int((donor, section) in code_sections),
                    "partial_public_only",
                ]
            )

    input_paths = sorted([*OLDER.glob("avey_*"), *NEWER.glob("avey_*")])
    with (OUT / "input_sha256.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["path", "bytes", "sha256"])
        for path in input_paths:
            writer.writerow([path.relative_to(ROOT).as_posix(), path.stat().st_size, sha256(path)])

    audit = {
        "audit_id": "G6-AVEY-ACCESS-2026-09-29-011",
        "effect_analysis_performed": False,
        "source_donors": len(donor_ids),
        "source_sections_reported_paper": 78,
        "supplementary_file9_spot_rows": len(spot_rows),
        "supplementary_file9_unique_donor_sections": len(spot_sections),
        "bayesspace_script_unique_donor_sections": len(code_sections),
        "public_partial_union_donor_sections": len(spot_sections | code_sections),
        "public_partial_union_donors": len({d for d, _ in spot_sections | code_sections}),
        "public_partial_union_fraction_of_paper_sections": len(spot_sections | code_sections) / 78,
        "full_donor_section_file_mapping_available": False,
        "library_batch_donor_mapping_available": True,
        "sequencing_batch_donor_mapping_available": False,
        "cdna_batch_and_slide_mapping_available": False,
        "synapse_access": "ANONYMOUS_READ_403_OWNER_APPROVAL_REQUIRED",
        "current_g6_status": "G6_NOT_EXECUTABLE",
    }
    (OUT / "audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

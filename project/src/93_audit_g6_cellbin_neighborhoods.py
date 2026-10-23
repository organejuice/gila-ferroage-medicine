#!/usr/bin/env python
"""Blind geometry-only CellBin neighborhood feasibility audit for G6-FP-2026-09-27-001."""
from __future__ import annotations

import csv
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 93_audit_g6_cellbin_neighborhoods.py <project_root>")
    root = Path(sys.argv[1])
    result_dir = root / "results/2026-09-27-g6-next-feasibility/cellbin"
    input_path = result_dir / "annotated_nuclei_metadata.tsv"
    if not input_path.exists():
        raise FileNotFoundError(f"Run the diagnosis-blind annotation step first: {input_path}")
    by_donor: dict[str, list[dict[str, str]]] = {}
    with input_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if any(key in row for key in ("diagnosis", "pathology", "braak", "cerad", "thal")):
                raise RuntimeError("Diagnosis/pathology field found in WP2 output")
            by_donor.setdefault(row["donor_id"], []).append(row)
    if len(by_donor) != 12:
        raise RuntimeError(f"Expected 12 donor maps, got {len(by_donor)}")

    radii = (10.0, 25.0, 50.0)
    rows: list[dict[str, object]] = []
    for donor, records in sorted(by_donor.items()):
        mg = [r for r in records if r["label"] == "MG"]
        asc = [r for r in records if r["label"] == "ASC"]
        unknown = sum(r["label"] == "unknown" for r in records)
        coords_mg = np.asarray([[float(r["x_um"]), float(r["y_um"])] for r in mg], dtype=float).reshape((-1, 2))
        coords_asc = np.asarray([[float(r["x_um"]), float(r["y_um"])] for r in asc], dtype=float).reshape((-1, 2))
        distances = np.empty(0, dtype=float)
        nearest = np.full(len(coords_mg), np.nan)
        if len(coords_mg) and len(coords_asc):
            tree_mg = cKDTree(coords_mg)
            tree_asc = cKDTree(coords_asc)
            nearest = tree_asc.query(coords_mg, k=1)[0]
            pair_lists = tree_mg.query_ball_tree(tree_asc, r=max(radii))
            distances = np.fromiter(
                (float(np.linalg.norm(coords_mg[i] - coords_asc[j]))
                 for i, js in enumerate(pair_lists) for j in js), dtype=float,
                count=sum(len(js) for js in pair_lists),
            )
        for radius in radii:
            rows.append({
                "donor_id": donor,
                "radius_um": radius,
                "n_MG_nuclei": len(mg),
                "n_ASC_nuclei": len(asc),
                "n_unknown_nuclei": unknown,
                "unique_Euclidean_MG_ASC_pairs": int(np.sum(distances <= radius)),
                "median_nearest_ASC_distance_from_MG_um": float(np.nanmedian(nearest)) if np.isfinite(nearest).any() else "",
                "fraction_MG_without_ASC_within_radius": float(np.mean(nearest > radius)) if len(nearest) else "",
                "bounding_box_width_um": float(np.ptp(np.vstack((coords_mg, coords_asc))[:, 0])) if len(mg) + len(asc) else "",
                "bounding_box_height_um": float(np.ptp(np.vstack((coords_mg, coords_asc))[:, 1])) if len(mg) + len(asc) else "",
                "tissue_hole_mask_available": False,
                "histology_validation_available": False,
            })

    eligible_counts = {
        radius: sum(1 for row in rows if row["radius_um"] == radius and row["unique_Euclidean_MG_ASC_pairs"] >= 20)
        for radius in radii
    }
    geometry_radius = next((r for r in radii if eligible_counts[r] >= 10), None)
    status = "NEIGHBORHOOD_NOT_ELIGIBLE_MASK_UNAVAILABLE"
    out_tsv = result_dir / "geometric_neighborhood_pair_counts.tsv"
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    gate_rows = [{
        "radius_um": r,
        "donors_with_ge20_Euclidean_pairs": eligible_counts[r],
        "required_donors": 10,
        "geometric_count_rule_pass": eligible_counts[r] >= 10,
        "tissue_hole_mask_available": False,
        "neighborhood_gate_status": status,
        "candidate_radius_by_count_rule_only": geometry_radius if geometry_radius is not None else "NONE",
    } for r in radii]
    with (result_dir / "neighborhood_geometric_gate.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(gate_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(gate_rows)
    report = [
        "# G6 CellBin geometric neighborhood feasibility (diagnosis blind)",
        "",
        f"Status: {status}",
        "Pair counts use the predeclared Euclidean 10/25/50 μm distances between exploratory MG and ASC nucleus centroids, within donor only.",
        "No paired histology image or independently auditable tissue-hole mask was available in the locked input audit. Therefore the Euclidean counts are descriptive geometry only; no radius is declared eligible because the no-crossing-tissue-holes condition cannot be verified.",
        f"Smallest radius meeting the count threshold before the mask condition: {geometry_radius if geometry_radius is not None else 'none'}.",
        "No disease labels or group effects were used.",
    ]
    (result_dir / "neighborhood_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

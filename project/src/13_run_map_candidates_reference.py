"""Run the official MAP checkpoint for each locally identified compound.

The fixed 0.1 uM input is a reference concentration only.  These outputs are
kept as VC model outputs and are not promoted to exposure-matched efficacy or
formula-level claims.
"""

from __future__ import annotations

import gc
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from omegaconf import OmegaConf


ROOT = Path(__file__).resolve().parents[1]
MAP_ROOT = ROOT / "external" / "MAP" / "MAP"
sys.path.insert(0, str(MAP_ROOT))

from demo import load_control_bulk  # noqa: E402
from model.model import MAPmodel  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    np.random.seed(42)
    raw = ROOT / "data" / "raw" / "MAP_official"
    minimal = ROOT / "data" / "interim" / "MAP_minimal_dataset" / "minimal_dataset"
    output = ROOT / "results" / "2026-09-20-drug-side-execution" / "map_reference_candidates"
    output.mkdir(parents=True, exist_ok=True)

    os.environ["MAP_MOLSTM_CKPT"] = str(raw / "molecule_model.pth")
    os.environ["MAP_MOLSTM_VOCAB"] = str(raw / "bart_vocab.txt")
    os.environ["TRANSFORMERS_CACHE"] = str(ROOT / "data" / "interim" / "hf_cache")

    cfg = OmegaConf.load(ROOT / "configs" / "MAP_se600m_local.yaml")
    se_ckpt = ROOT / "data" / "interim" / "MAP_se_from_epoch_3.safetensors"
    full_ckpt = raw / "epoch_3.pt"
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    compounds = pd.read_csv(ROOT / "results" / "tables" / "drug_compound_master.tsv", sep="\t")
    compounds = compounds[compounds["Canonical_SMILES"].fillna("").astype(str).str.len() > 0].copy()
    compounds = compounds.drop_duplicates(subset=["PubChem_CID"]).reset_index(drop=True)

    model = MAPmodel(
        se_ckpt=str(se_ckpt),
        se_cfg=cfg,
        smile_encoder="PrimeKG_ver1",
        hvg_info=None,
    )
    checkpoint = torch.load(full_ckpt, map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint.get("model_state_dict", checkpoint), strict=True)
    del checkpoint
    gc.collect()
    model = model.half().to(device)
    model.eval()

    gene_ids, expressions = load_control_bulk(
        "CVCL_0023",
        str(minimal / "preprocessed"),
        str(minimal / "preprocessed_se_inputs_memmap"),
        set_size=1,
    )
    ctrl_gene = gene_ids.unsqueeze(0).to(device)
    ctrl_expr = expressions.unsqueeze(0).to(device).half()

    rows = []
    first_pred_hvg = None
    for row in compounds.itertuples(index=False):
        cid = int(row.PubChem_CID)
        smiles = str(row.Canonical_SMILES)
        with torch.no_grad():
            drug_vec = model.pert_model.encode_drug([smiles]).cpu().float().numpy()[0]
        with torch.no_grad():
            pred_embs, pred_hvgs = model(
                ctrl_gene,
                ctrl_expr,
                [smiles],
                torch.tensor([0.1], dtype=torch.float16, device=device),
            )
        pred_emb = pred_embs.cpu().float().numpy()[0, 0]
        pred_hvg = pred_hvgs.cpu().float().numpy()[0, 0]
        if first_pred_hvg is None:
            first_pred_hvg = pred_hvg.copy()
        pred_hvg_maxdiff_vs_first = float(np.max(np.abs(pred_hvg - first_pred_hvg)))
        np.save(output / f"CID_{cid}_pred_embedding.npy", pred_emb)
        np.save(output / f"CID_{cid}_pred_hvg.npy", pred_hvg)
        rows.append(
            {
                "Herb": row.Herb,
                "Ingredient": row.Ingredient,
                "PubChem_CID": cid,
                "Canonical_SMILES": smiles,
                "cell_line": "CVCL_0023",
                "reference_concentration_uM": 0.1,
                "set_size": 1,
                "status": "completed_real_model_reference_run",
                "pred_embedding_l2": float(np.linalg.norm(pred_emb)),
                "pred_hvg_l2": float(np.linalg.norm(pred_hvg)),
                "pred_hvg_mean": float(pred_hvg.mean()),
                "pred_hvg_sd": float(pred_hvg.std(ddof=1)),
                "drug_vector_l2": float(np.linalg.norm(drug_vec)),
                "pred_hvg_maxdiff_vs_first": pred_hvg_maxdiff_vs_first,
                "interpretation_ceiling": "VC-only; no exposure match, no disease-state validation, no candidate ranking",
            }
        )

    table = pd.DataFrame(rows)
    table.to_csv(output / "drug_side_map_reference_outputs.tsv", sep="\t", index=False)
    meta = {
        "status": "completed_real_model_reference_runs",
        "n_unique_pubchem_cids": len(rows),
        "device": str(device),
        "cell_line": "CVCL_0023",
        "reference_concentration_uM": 0.1,
        "set_size": 1,
        "control_sampling_seed": 42,
        "full_checkpoint_sha256": sha256(full_ckpt),
        "se_checkpoint_sha256": sha256(se_ckpt),
        "interpretation_ceiling": "VC-only; no exposure match, no disease-state validation, no candidate ranking",
    }
    (output / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(meta, ensure_ascii=False, indent=2))
    print(table[["PubChem_CID", "Ingredient", "status", "pred_hvg_mean", "pred_hvg_sd"]].to_string(index=False))


if __name__ == "__main__":
    main()

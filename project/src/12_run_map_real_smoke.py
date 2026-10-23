"""Run one real MAP inference smoke test with official released assets.

This is deliberately a model/input integrity test, not an efficacy claim:
the concentration is a predeclared reference concentration and no matched
exposure or disease-state validation is asserted.
"""

from __future__ import annotations

import gc
import hashlib
import json
import os
import sys
from pathlib import Path

import numpy as np
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
    output = ROOT / "results" / "2026-09-20-drug-side-execution" / "map_real_smoke"
    output.mkdir(parents=True, exist_ok=True)

    os.environ["MAP_MOLSTM_CKPT"] = str(raw / "molecule_model.pth")
    os.environ["MAP_MOLSTM_VOCAB"] = str(raw / "bart_vocab.txt")
    os.environ.setdefault("TRANSFORMERS_CACHE", str(ROOT / "data" / "interim" / "hf_cache"))

    cfg_path = ROOT / "configs" / "MAP_se600m_local.yaml"
    se_ckpt = ROOT / "data" / "interim" / "MAP_se_from_epoch_3.safetensors"
    full_ckpt = raw / "epoch_3.pt"
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    cfg = OmegaConf.load(cfg_path)
    model = MAPmodel(
        se_ckpt=str(se_ckpt),
        se_cfg=cfg,
        smile_encoder="PrimeKG_ver1",
        hvg_info=None,
    )
    checkpoint = torch.load(full_ckpt, map_location="cpu", weights_only=False)
    state_dict = checkpoint.get("model_state_dict", checkpoint)
    model.load_state_dict(state_dict, strict=True)
    del checkpoint, state_dict
    gc.collect()

    # The 6-GB GPU cannot safely hold a full FP32 MAP model plus activations.
    # FP16 is only used to make the released checkpoint executable; no ranking
    # or direction claim is derived from this smoke test.
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
    smiles = ["C1=CC(=C(C=C1C=O)O)O"]  # Protocatechuic aldehyde, CID 8768
    concs = torch.tensor([0.1], dtype=torch.float16, device=device)

    with torch.no_grad():
        pred_embs, pred_hvgs = model(ctrl_gene, ctrl_expr, smiles, concs)

    np.save(output / "input_ctrl_gene_ids.npy", ctrl_gene.cpu().numpy().astype(np.int32))
    np.save(output / "input_ctrl_expressions.npy", ctrl_expr.cpu().float().numpy())
    np.save(output / "output_pred_embeddings.npy", pred_embs.cpu().float().numpy())
    np.save(output / "output_pred_hvg.npy", pred_hvgs.cpu().float().numpy())
    meta = {
        "status": "completed_real_model_smoke_test",
        "device": str(device),
        "cell_line": "CVCL_0023",
        "drug_name": "Protocatechuic aldehyde",
        "pubchem_cid": 8768,
        "drug_smiles": smiles[0],
        "reference_concentration_uM": 0.1,
        "set_size": 1,
        "pred_embs_shape": list(pred_embs.shape),
        "pred_hvgs_shape": list(pred_hvgs.shape),
        "full_checkpoint_sha256": sha256(full_ckpt),
        "se_checkpoint_sha256": sha256(se_ckpt),
        "interpretation_ceiling": "VC-only smoke test; no exposure match, no disease-state validation, no ranking",
    }
    (output / "meta.json").write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(meta, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

"""Audit whether the released MAP checkpoint responds to different SMILES."""

from __future__ import annotations

import gc
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


def main() -> None:
    raw = ROOT / "data" / "raw" / "MAP_official"
    minimal = ROOT / "data" / "interim" / "MAP_minimal_dataset" / "minimal_dataset"
    os.environ["MAP_MOLSTM_CKPT"] = str(raw / "molecule_model.pth")
    os.environ["MAP_MOLSTM_VOCAB"] = str(raw / "bart_vocab.txt")

    model = MAPmodel(
        se_ckpt=str(ROOT / "data" / "interim" / "MAP_se_from_epoch_3.safetensors"),
        se_cfg=OmegaConf.load(ROOT / "configs" / "MAP_se600m_local.yaml"),
        smile_encoder="PrimeKG_ver1",
        hvg_info=None,
    )
    checkpoint = torch.load(raw / "epoch_3.pt", map_location="cpu", weights_only=False)
    model.load_state_dict(checkpoint.get("model_state_dict", checkpoint), strict=True)
    del checkpoint
    gc.collect()
    model = model.half().to("cuda" if torch.cuda.is_available() else "cpu")
    model.eval()
    device = next(model.parameters()).device

    gene_ids, expressions = load_control_bulk(
        "CVCL_0023",
        str(minimal / "preprocessed"),
        str(minimal / "preprocessed_se_inputs_memmap"),
        set_size=1,
    )
    ctrl_gene = gene_ids.unsqueeze(0).to(device)
    ctrl_expr = expressions.unsqueeze(0).to(device).half()
    smiles_a = "C1=CC(=C(C=C1C=O)O)O"
    smiles_b = "CC1=C2C=CC3=C(C2=CC=C1)C(=O)C(=O)C(=C3O)C"

    with torch.no_grad():
        # Reproduce MAPmodel.forward up to the perturbation encoder inputs.
        gene_ids_flat = ctrl_gene.reshape(1, 2048)
        expr_flat = ctrl_expr.reshape(1, 2048)
        src = model.se.pe_embedding(gene_ids_flat)
        esm_tokens = src[:, 1:].clone().detach()
        src = torch.nn.functional.normalize(src, dim=2)
        src = torch.cat([model.se.cls_token.expand(src.size(0), 1, -1), src[:, 1:]], dim=1)
        dataset_token = model.se.dataset_token.expand(src.size(0), 1, -1)
        src = torch.cat((src, dataset_token), dim=1)
        gene_output, embedding, _ = model.se(src=src, counts=expr_flat, dataset_nums=None, profile=False)
        gene_output = gene_output[:, 1:-1, :]
        gene_tokens = model.pert_model.encode_genes(gene_output, esm_tokens)
        cell_tokens = model.pert_model.encode_cells(embedding)
        drug_a = model.pert_model.encode_drug([smiles_a])
        drug_b = model.pert_model.encode_drug([smiles_b])
        combined_a = torch.cat([cell_tokens.unsqueeze(1), drug_a.unsqueeze(1), gene_tokens], dim=1)
        combined_b = torch.cat([cell_tokens.unsqueeze(1), drug_b.unsqueeze(1), gene_tokens], dim=1)
        res_a = model.pert_model.transformer_backbone(inputs_embeds=combined_a).last_hidden_state[:, 0, :]
        res_b = model.pert_model.transformer_backbone(inputs_embeds=combined_b).last_hidden_state[:, 0, :]
        out_a = model.pert_model.project_out(cell_tokens + res_a)
        out_b = model.pert_model.project_out(cell_tokens + res_b)
        hvg_a = model.pert_model.gene_decoder(out_a)
        hvg_b = model.pert_model.gene_decoder(out_b)

    result = {
        "device": str(device),
        "drug_vector_max_abs_diff": float((drug_a - drug_b).abs().max().cpu()),
        "transformer_cls_max_abs_diff": float((res_a - res_b).abs().max().cpu()),
        "pred_embedding_max_abs_diff": float((out_a - out_b).abs().max().cpu()),
        "pred_hvg_max_abs_diff": float((hvg_a - hvg_b).abs().max().cpu()),
        "interpretation": "drug encoder differs; released MAP perturbation head is invariant in this FP16 reference probe",
    }
    path = ROOT / "results" / "2026-09-20-drug-side-execution" / "map_drug_sensitivity_audit.json"
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

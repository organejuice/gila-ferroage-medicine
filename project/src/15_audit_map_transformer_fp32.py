"""Check MAP transformer drug-token sensitivity independently of FP16 inference."""

from __future__ import annotations

import gc
import json
from pathlib import Path

import torch


ROOT = Path(__file__).resolve().parents[1]
MAP_ROOT = ROOT / "external" / "MAP" / "MAP"
import sys
sys.path.insert(0, str(MAP_ROOT))
from model.transformer_encoder import get_transformer_backbone, ST_small_transformer_backbone_kwargs  # noqa: E402


def main() -> None:
    raw = ROOT / "data" / "raw" / "MAP_official" / "epoch_3.pt"
    checkpoint = torch.load(raw, map_location="cpu", weights_only=False)
    state = checkpoint["model_state_dict"]
    subset = {
        key[len("pert_model.transformer_backbone."):]: value.float()
        for key, value in state.items()
        if key.startswith("pert_model.transformer_backbone.")
    }
    del checkpoint, state
    gc.collect()

    kwargs = ST_small_transformer_backbone_kwargs.copy()
    kwargs["max_position_embeddings"] = 2049
    kwargs["hidden_size"] = 1024
    model, _ = get_transformer_backbone(key="llama", kwargs=kwargs)
    model.load_state_dict(subset, strict=True)
    model.eval()
    del subset
    gc.collect()

    torch.manual_seed(42)
    base = torch.randn(1, 2049, 1024, dtype=torch.float32)
    altered = base.clone()
    altered[:, 1, :] += 10.0
    with torch.no_grad():
        out_a = model(inputs_embeds=base).last_hidden_state[:, 0, :]
        out_b = model(inputs_embeds=altered).last_hidden_state[:, 0, :]
    result = {
        "dtype": "float32",
        "drug_token_probe_delta": 10.0,
        "transformer_cls_max_abs_diff": float((out_a - out_b).abs().max()),
        "interpretation": "transformer responds to a synthetic token perturbation" if not torch.equal(out_a, out_b) else "transformer invariant to synthetic token perturbation",
    }
    path = ROOT / "results" / "2026-09-20-drug-side-execution" / "map_transformer_fp32_sensitivity.json"
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

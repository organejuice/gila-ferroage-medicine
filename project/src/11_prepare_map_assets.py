"""Prepare auditable MAP assets derived from the official released checkpoint."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import torch
from safetensors.torch import save_file


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw" / "MAP_official"
INTERIM_DIR = PROJECT_ROOT / "data" / "interim"
RESULT_DIR = PROJECT_ROOT / "results" / "2026-09-20-drug-side-execution"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    epoch_path = RAW_DIR / "epoch_3.pt"
    out_path = INTERIM_DIR / "MAP_se_from_epoch_3.safetensors"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    checkpoint = torch.load(epoch_path, map_location="cpu", weights_only=False)
    state = checkpoint["model_state_dict"]
    se_state = {
        key[3:]: value.detach().cpu().contiguous().clone()
        for key, value in state.items()
        if key.startswith("se.") and key != "se.pe_embedding.weight"
    }
    save_file(se_state, str(out_path), metadata={"source": "epoch_3.pt", "prefix": "se."})

    audit = {
        "source_checkpoint": str(epoch_path),
        "source_sha256": sha256(epoch_path),
        "derived_se_checkpoint": str(out_path),
        "derived_sha256": sha256(out_path),
        "source_state_keys": len(state),
        "derived_se_keys": len(se_state),
        "excluded_keys": ["se.pe_embedding.weight"],
        "checkpoint_epoch": checkpoint.get("epoch"),
        "checkpoint_global_step": checkpoint.get("global_step"),
    }
    INTERIM_DIR.joinpath("MAP_checkpoint_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    RESULT_DIR.joinpath("MAP_checkpoint_audit.json").write_text(
        json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(audit, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Run the prospective 4+4 metadata-only G6 candidate screen."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: 98_screen_g6_candidate_size_003.py <project_root>")
    core_path = Path(__file__).with_name("97_screen_g6_candidate_size_002.py")
    spec = importlib.util.spec_from_file_location("g6_candidate_size_screen", core_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load candidate screen: {core_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run(Path(sys.argv[1]), "configs/g6_candidate_size_screen_003.yaml")


if __name__ == "__main__":
    raise SystemExit(main())

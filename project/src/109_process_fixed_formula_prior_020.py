from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_SCRIPT = ROOT / "src" / "108_process_fixed_formula_prior_019.py"
RUN_ID = "2026-09-30-fixed-formula-drug-prior-020"

spec = importlib.util.spec_from_file_location("fixed_formula_prior_shared", BASE_SCRIPT)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not load shared implementation: {BASE_SCRIPT}")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
module.OUT = ROOT / "results" / RUN_ID
module.PARAM_PATH = ROOT / "configs" / "fixed_formula_drug_prior_020.yaml"
module.RUN_SCRIPT = Path(__file__).resolve()
module.RUN_ID = RUN_ID
module.DECISION_ID = "DRUG-FORMULA-PRIOR-CORRECTION-2026-09-30-020"
module.main()

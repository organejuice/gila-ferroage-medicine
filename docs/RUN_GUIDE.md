# Running and checking the frozen study source

## What can be checked immediately

From the repository root, run `python -B tools/check_source.py`. This standard-library check validates the publication SHA256 list, the original ESM manifest with its two declared privacy substitutions, and Python syntax. No biological analysis is run.

The original tests are retained for transparency. Several combine toy-data checks with checks of study-specific result files. Running all tests on this code-only archive is therefore expected to encounter missing-data or dependency errors. A syntax/integrity pass is not a numerical reproduction pass. No claim is made that all original tests were rerun for this public export.

## Environment records

The study reused separate existing local environments. `project/environment.yml` and `project/requirements.txt` are historical initial specifications; they are not sufficient to recreate every branch. `project/configs/TCM_select_024_pip_freeze.txt` is a pinned snapshot of one chemical-analysis environment, not a universal lock for the entire study. Inspect imports and the relevant source/configuration before installing anything; reproduce in an isolated environment rather than modifying a working environment.

The recorded manuscript environment includes scTenifoldKnk 1.0.3, scTenifoldNet 1.3, R 4.4.3 for the virtual-knockout branch and fgsea 1.32.4 for enrichment. The principal-component regression acceleration uses the separately pinned upstream `pcNet` source in `project/external/virtual_ko_033/`; the combined backend is not an unchanged official-package run. Some historical R launchers select other recorded local R installations. The structural branch records AutoDock Vina 1.2.7, OpenBabel 3.1.0, RDKit 2025.09.6 and Meeko 0.7.1. Other branches have their own historical versions. Do not silently substitute package or model versions and describe the result as the frozen study run.

Python branches use combinations of pandas, NumPy, SciPy, PyYAML, RDKit, requests, networkx, ADMET-AI, PyTorch and branch-specific libraries. R branches use combinations of Matrix, edgeR, limma, Seurat/SeuratObject, scTenifoldKnk/scTenifoldNet, fgsea, RSpectra, yaml/jsonlite and annotation packages. Editorial scripts additionally require document/PDF/plotting software. External model code, weights, binaries and database snapshots are not distributed here.

## Historical analysis entry points

Run scientific steps from **`project/`**, after obtaining and validating their registered inputs and adapting executable paths. The table identifies existing source entry points; it is not an instruction to run every historical branch sequentially.

| Analysis | Source entry points | Required context |
| --- | --- | --- |
| Human G1 associations | `src/22_prepare_disease_side_v3_1_inputs.py`, `src/23_prepare_d7_d8_pseudobulk.R`, `src/25_run_g1_v3_1.R`, `src/26_finalize_disease_side_v3_1.py` | Audited source counts, donor/region/cell annotations, prepared D3 counts and frozen cohort specifications; earlier `20_*` and later frozen definitions must not be mixed without checking provenance |
| G2/G3/G4 model evidence | `src/28_run_g2_v3_1.py`, `src/29_run_g2_v3_2.py`, `src/31_reassess_g3_v3_2_model_scoped.py`, `src/32_run_g4_v3_2.py` | Registered peak-to-gene/model inputs and frozen module definitions; human convergence remains unresolved |
| G5 astrocyte primary model | `src/47_run_g5_v3_4.R`, `src/48_finalize_g5_v3_4.py` | GSE261698 expression, sample metadata and the frozen G5 protocol; original fractional values must be retained |
| Layered formula targets | `src/113_supplement_formula_targets_024.py`, `src/114_finalize_formula_targets_024.py`, `src/121_finalize_continuation_024.py` | Exact-structure and assay/source registries, database exports and reviewed classifications |
| Stage-A enrichment/network | `src/148_run_stage_a_030.py` | Frozen stage-A inputs, gene maps, GMT snapshots and STRING network. `--recompute-network` changes the use of cached null outputs; inspect source before use |
| Candidate joins and docking | `src/reuse_g1_g5_docking_candidates.py`, `src/expanded_docking_20.py`, `src/audit_expanded_docking_20.py`, `src/sensitivity_expanded_docking_20.py` | Frozen E1 evidence, receptor/ligand structures, target-site/cofactor configuration and external executable paths |
| Primary virtual knockout 033 | `src/132_run_stratified_virtual_ko.R`, `src/133_summarize_virtual_ko.py`, `src/141_audit_virtual_ko_artifacts.R` | Frozen GSE160936 counts/metadata, four contexts, two seeds and registered target-retention sensitivity |
| Pathway analysis 034 | `src/run_virtual_ko_pathways_034.py` and associated preparation/audit/report scripts | Original KO outputs and frozen GMT files |
| Verification/review 035 | `src/run_virtual_ko_verification_035.py`, `src/run_final_scientific_review_035.py` | Registered 29-context scope, source networks, resampling/matched-reference outputs and audit inputs |
| Matching correction 036 | `src/run_matching_correction_036.py` | Full raw-count library denominators for the identical frozen cells, historical WT/KO inputs and frozen verification scope; Windows resource checks are built in |

There is no one-command supported `dvc repro` for a fresh clone. `project/dvc.yaml` and `.lock` preserve local dependency records, but the archive does not provide their data objects or a public DVC remote. In particular, neither `dvc pull` nor a cross-machine restore has been verified. Preserve the historical lock as provenance and inspect stages before an intentional rerun.

## Focused original tests

The retained tests include compound identity/evidence classification, candidate-panel gates, count conservation and rejection checks, BH correction, toy-network distance/sampling, matched-reference correction and deterministic fgsea examples. Examples after installing the corresponding branch dependencies and restoring required inputs:

```text
cd project
python -B -m unittest discover -s tests -p test_stage_a_030.py -v
Rscript tests/test_matching_correction_036.R
Rscript tests/test_virtual_ko_pathways_034.R
```

These commands are documented entry points, not newly certified cross-platform results. Some test modules read configuration at import time. Do not treat a missing frozen input as a zero result or remove failed scientific gates to obtain a passing run.

Legacy manuscript builders (`build_fig_yzdhd_submission.py`, etc.) require private author-only metadata and the complete editorial workspace. Their historical filename and email assertions have no role in the biological computation. This source archive is not a turnkey manuscript-generation package.

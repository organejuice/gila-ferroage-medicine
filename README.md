# gila-ferroage-medicine

Frozen analysis-source snapshot supporting the manuscript:

**Glial lactylation, senescence and iron homeostasis in Alzheimer's disease: an integrative transcriptomic and network pharmacology analysis of Yizhi Dihuang Decoction**

Hengjun Liu, Shi Jian, TingTing Zhang, Shangyi Cai, Yunze Shi, Zhengtong Xu, Yude Wang, Yu Zhang and Rui Qie. Correspondence: Rui Qie, qierui@hljucm.edu.cn.

This repository provides the custom analysis scripts, versioned configurations, recorded environment specifications, tests and provenance exported from Online Resource 3 of the Metabolic Brain Disease submission. It is a source snapshot of the completed study, rather than a tested, self-contained installation. Public raw count matrices, large frozen intermediates, full analysis outputs, model weights and a restorable DVC data remote are not included. Full end-to-end execution on another machine has not been validated.

## Start here

1. [Run guide](docs/RUN_GUIDE.md): historical analysis entry points, dependencies, tests and expected boundaries.
2. [Data and path requirements](docs/DATA_AND_PATHS.md): source accessions, absent inputs and platform-specific paths.
3. [Third-party notices](THIRD_PARTY_NOTICES.md): inherited GPL components, pinned solver and dependency licensing.
4. [Publication changes](PUBLICATION_CHANGES.json): source hashes, the two privacy-only substitutions and added release documentation.

To check the published source files with Python 3.11 or later and its standard library:

```text
git clone https://github.com/organejuice/gila-ferroage-medicine.git
cd gila-ferroage-medicine
python -B tools/check_source.py
```

This check validates file hashes and parses Python source. It does not install packages, access the network, execute scientific analyses, or validate the manuscript's biological conclusions.

## Repository map

| Path | Contents |
| --- | --- |
| `project/src/` | Frozen Python, R and PowerShell code for G1–G5, pharmacological evidence, docking, virtual knockout review and historical feasibility/editorial steps |
| `project/configs/` | Frozen parameters, historical pipeline snapshots and an environment-specific pip freeze |
| `project/tests/` | Original toy-data, numerical, evidence-gating and frozen-output regression tests |
| `project/environment.yml`, `project/requirements.txt` | Historical starting environment specifications, not complete unified lock files |
| `project/params.yaml`, `project/dvc.yaml`, `project/dvc.lock` | Original parameter and local data-pipeline records; their referenced data are not bundled |
| `project/external/virtual_ko_033/` | Pinned principal-component regression solver, original upstream metadata and GPL license text, restored from the study's local frozen source |
| `project/README.md` | Unchanged historical project status; its local links may refer to material outside this source archive |
| `claim_evidence.tsv`, `original_manuscript_provenance.json` | Claim-to-local-evidence map and original input provenance |
| `code_manifest.json` | Original 329-file ESM source manifest; two documented privacy substitutions have different published hashes |
| `SHA256SUMS.txt` | Current publication-file integrity manifest |

## Scientific scope

The manuscript integrates donor-level human AD microglial associations, model-scoped lactylation/senescence/iron-related programs, chemical and compound–target evidence, structural assessments, and virtual knockout sensitivity review. The frozen human association analysis identified 413 FDR-significant genes and a 244-gene complete-region subset. The study's primary network-proximity analysis was not significant. The final virtual knockout review considered 15 targets in 29 contexts.

These code and model outputs do not establish Yizhi Dihuang Decoction exposure, target engagement, efficacy, a confirmed human H4K12la mechanism, target specificity or causal microglia–astrocyte communication. Docking scores are not measured affinities. The historical G2c human-translation and G6 spatial-validation routes were closed for insufficient qualified evidence; their retained audit scripts do not represent completed spatial disease-effect analyses. Post-snapshot project extensions are outside this release.

Historical `FIG` filenames identify earlier editorial provenance and have been retained. They do not identify the current target journal. No scientific method, parameter, result or conclusion was changed when preparing this repository.

## Reuse and availability

Original scientific files are retained byte for byte, except that two legacy editorial wrapper scripts replace non-corresponding-author email assertions with `example.invalid` placeholders. The original source archive and per-file hashes remain identified in `PUBLICATION_CHANGES.json`. A previously omitted pinned third-party solver and its metadata/license are included without modification and match the archived DVC checksum.

This repository is publicly accessible. It does not apply a blanket permissive license to all files. The identified GPL-derived components remain under their stated GPL version-2-or-later terms; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and the included license text. Other files have no additional repository-wide license grant in this release. Third-party databases, packages, structures and datasets retain their own terms.

For a stable citation, record the exact Git commit or release tag used. No DOI or accepted-publication status is asserted by this source snapshot.

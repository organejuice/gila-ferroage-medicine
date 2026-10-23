# Data acquisition and path adaptation

## Primary public sources

| Source | Study role |
| --- | --- |
| [GSE160936](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160936) | Human discovery cohort; EC/SSC microglial association analysis and glial virtual-knockout source |
| [GSE148822](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148822) | Human primary validation; OC/OTC regions |
| [GSE174367](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174367) | Human secondary validation; PFC |
| [GSE188765](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188765) | Model H4K12la peak/annotation context, not validated human convergence |
| [GSE227221](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227221) | iMGL amyloid-fibril time course used for model-scoped transcriptional programs |
| [GSE261698](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261698) | Primary APP/PS1-21 astrocyte analysis |
| [GSE291239](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291239), [GSE74614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE74614), [GSE290101](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290101) | Registered G5 supporting/sensitivity contexts |
| [GSE214979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214979) | Exploratory human cohort only |

Other candidate accessions retained in configurations belong to historical feasibility or scope audits. Their presence is not evidence that they contributed disease-effect estimates. Source-study metadata and the manuscript's donor/region inclusion rules must be checked before analysis. An accession identifies a study, not a validated automatic sample mapping.

Additional inputs include the frozen ferro-aging reference (original DOI [10.1016/j.cmet.2026.02.010](https://doi.org/10.1016/j.cmet.2026.02.010); correction notice [10.1016/j.cmet.2026.08.011](https://doi.org/10.1016/j.cmet.2026.08.011)), BATMAN-TCM/TCMBank ingredient exports, PubChem-linked identities, reviewed ChEMBL/BindingDB assay records, STRING 12.0 channel snapshots, gene mappings, GO/KEGG/Reactome GMT files and registered PDB receptor/ligand structures. These sources retain their own access conditions and licenses; this repository does not redistribute their full datasets.

## Inputs absent from the repository

Raw counts, full single-cell objects, donor-level pseudobulk matrices, selected cell identities, frozen module tables, gene-set libraries, complete assay review tables, docking input structures and output models, WT networks, KO/regulation matrices, resampling outputs and complete result registries are not bundled. Some derived tables are supplied separately in the manuscript's Online Resource 2. `claim_evidence.tsv` and `original_manuscript_provenance.json` identify local evidence files and hashes; they are source maps, not data downloads. The latter preserves original local path strings solely as provenance.

The pinned `pcNet` file was absent from the original supplementary ZIP but has been recovered here from the existing study project. Its MD5 is `974ad35b48b4de423f51ecbab21730d3` (21,110 bytes), matching the archived 033 DVC dependency record. Its upstream commit is `ce4822c22a66ffcf1ed7ff3f9bb42f0de96dbbea`. See `project/external/virtual_ko_033/` for source, original authorship, license text and local provenance. It is one restored code dependency, not a replacement for the missing data.

## Reconstructing a run responsibly

1. Select the manuscript analysis branch and its exact frozen configuration. Do not infer chronology from numeric filenames alone or combine exploratory/retired branches with primary analyses.
2. Retrieve the corresponding public source files and permitted database exports. Record source URL/accession, original filename, retrieval date, size and checksum. Check the recorded hashes whenever the exact frozen source is available; current server exports can differ.
3. Rebuild the expected `project/data/raw/`, `project/data/interim/`, `project/data/processed/` and `project/results/` layout. Keep raw data read only. Restore reviewed protocol/sample mappings as required by the branch; several protocol documents are outside this code-only archive.
4. Adapt executable paths in wrappers and DVC commands. Historical examples include `D:/conda/...`, `D:/R/...`, `D:/python/...` and `F:/TCM select/...`. Some scripts use Windows `ctypes` resource checks. A path substitution alone does not establish portability or numerical equivalence.
5. Use compatible isolated environments and inspect toy/numerical tests before a full rerun. Preserve the original configurations, seeds, non-estimable outcomes, reference-matching constraints and inference boundaries. Record intentional environment or code changes in a separate run.
6. Compare outputs against the manuscript/supplementary evidence only after all data and numerical checks pass. A successful source hash check does not reproduce the biological analysis.

The repository supplies no credentials, private datasets, public DVC data endpoint, or complete cross-machine restoration recipe. Contact the corresponding author for clarification of frozen inputs that cannot be reconstructed from the referenced source records.

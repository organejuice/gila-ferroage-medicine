# Third-party notices and licensing scope

This file records identified reused components in the frozen study source. It does not replace their original notices or apply a blanket new license to the entire repository.

## scTenifold-derived source

The following files already identify their derivation and **GNU General Public License version 2 or later (GPL-2.0-or-later)** terms in their retained source headers:

- `project/src/130_virtual_ko_cached_backend.R`: adapted from scTenifoldKnk 1.0.3, [cailab-tamu/scTenifoldKnk](https://github.com/cailab-tamu/scTenifoldKnk). This adapter reuses a WT graph across deletions and calls installed scTenifold components.
- `project/src/135_sparse_cp_equivalent.R`: GPL >=2 adaptation of scTenifoldNet 1.3 `cpDecomposition`/`tensorDecomposition`, implemented through sparse slice-wise contractions.
- `project/src/138_implicit_manifold.R`: GPL >=2 adaptation of scTenifoldNet 1.3 `manifoldAlignment`, represented through an implicit operator.

Original source attribution and modification descriptions in these headers are preserved. The included [GNU GPL version 2 text](project/external/virtual_ko_033/COPYING) accompanies these identified components; their existing version-2-or-later grant is unchanged.

## Pinned scTenifoldNet principal-component regression solver

`project/external/virtual_ko_033/pcNet_upstream.R` is the unchanged local study copy of the upstream `R/pcNet.R` from [cailab-tamu/scTenifoldNet](https://github.com/cailab-tamu/scTenifoldNet) at commit **`ce4822c22a66ffcf1ed7ff3f9bb42f0de96dbbea`**. The local original `upstream_commit.txt`, `DESCRIPTION`, `COPYING` and `README.md` are included without modification. The corresponding upstream source reference is:

<https://github.com/cailab-tamu/scTenifoldNet/blob/ce4822c22a66ffcf1ed7ff3f9bb42f0de96dbbea/R/pcNet.R>

The retained DESCRIPTION records scTenifoldNet version 1.4.1, GPL (>=2), and authors Daniel Osorio, Yan Zhong, Guanxun Li, Jianhua Huang and James Cai. The study used this pinned solver with installed scTenifoldNet 1.3/scTenifoldKnk 1.0.3 components. These component versions must not be conflated. The restored local source checksum matches the 033 archived DVC dependency record (MD5 `974ad35b48b4de423f51ecbab21730d3`, 21,110 bytes). No upstream upgrade was performed for this publication export. The fixed upstream URL was not independently re-downloaded in this export; identity was verified against the existing frozen local source and lock.

The inherited README refers to historical numerical-equivalence tests in the full study workspace. Those outputs are not bundled and those analyses were not rerun for this export.

## Other code, software and data

No additional repository-wide license is assigned to the remaining author-written files in this release. Public visibility does not itself create a permissive license grant. This statement does not restrict rights already granted by an applicable third-party license.

The code calls external packages and tools including scTenifoldKnk/scTenifoldNet, fgsea, Seurat, edgeR/limma, RDKit, ADMET-AI, AutoDock Vina, OpenBabel, Meeko, and R/Python libraries. These dependencies are not relicensed by this repository. Their own source, attribution and licensing requirements apply. Full database exports, raw datasets, PDB structures and model weights are not included; consult each original provider before obtaining or redistributing them.

`project/configs/TCM_select_024_pip_freeze.txt` also records an editable TCMxPlore source dependency with its historical commit. This is an environment provenance record and is not a redistribution or new license grant for that project.

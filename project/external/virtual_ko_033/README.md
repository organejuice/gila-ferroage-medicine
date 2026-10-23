# Upstream algorithm reused for KO 033

`pcNet_upstream.R` is copied unchanged from cailab-tamu/scTenifoldNet at the commit in `upstream_commit.txt`. License GPL >=2, upstream authors as recorded in DESCRIPTION. Original source https://github.com/cailab-tamu/scTenifoldNet . Only its principal component regression solver is used; installed scTenifoldNet 1.3 supplies network sampling and manifold alignment, scTenifoldKnk 1.0.3 supplies direction enforcement and differential regulation.

Full backend equality and zero-edge-difference tests are recorded in results/2026-10-04-stratified-virtual-ko-033. The local sparse CP adapter derives from the official 1.3 dense implementation and matches its rounded outputs in deterministic numerical tests. This is a pinned composite implementation, not an unmodified official package run.

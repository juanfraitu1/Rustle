# ThesisVault methods audit — 2026-09-28

**Question:** does `ThesisVault/` (20 paper notes, 11 topic notes) hold a method we have not applied that would help attain O1/O2/O3?

**Answer:** no. Every method bearing on the three objectives is already in the code or already refuted. Method: grep of `src/`, `bench/`, `docs/` plus a read of the matching rows in `docs/NEGATIVE_RESULTS_REGISTER.md`. No experiment was rerun.

## Method-by-method status

| Vault method | Source | Status in Rustle |
|---|---|---|
| PSV discovery, read-correlation clustering, attraction/repulsion structure | Vollger 2019 (SDA) | Applied (`vg_family`, `psv_graph`). Soft-EM relaxation tried: changed 0/3,081 decisions. Depth-based collapse gate retired (`collapse_gate`, default off). |
| Best-hit plus mismatch margin as per-read baseline | Dishuck 2025 | Built as `--eichler-margin` and benchmarked; discards 83.7% of reads on Y ampliconic genes. |
| RNA-editing confound (A→G / T→C) | Clair3-RNA | Editing filter in `copy_assign`. |
| IsoCon significance test | Sahlin 2018 | Run at the assign gate's α: byte-identical output, so the δ = 0.005 threshold was inert. |
| Phase / linkage consistency across PSVs | longcallR | Covered by `origin_consistency_check` and the PSV co-observation analysis. K≥3 recombination obstruction is machine-checked. |
| WSSD-style read depth for collapsed or absent copies | Bailey 2002 | A read-depth proxy outscored every structural statistic. WGS trio k-mer copy-number check done (`docs/archive/2026-09/O3_WGS_TRIO_CN_2026-09-25.md`). |
| Single-linkage protein-identity families | Makova 2024, Pal 2026 | Superseded by the MCL and identity-lattice definition. |
| SQANTI3 categories, polyA / 5′ end support, SIRV spike-ins | LRGASP | SQANTI3 benchmarked; `--polish-tes` shipped. Only the human testis library is spiked. |
| StringTie2 / FLAIR / isoseq comparators | Kovaka 2019, Tang 2020 | Bakeoff done. |
| Phylogenetic paralog groups, intra-allelic-variation threshold | Guitart 2024, Dishuck 2025 | Not adopted. The threshold is arbitrary (Canzar objects), and register row 778 found TBC1D3's clusters are positional with no sequence signal. |
| Facility location / max-flow assignment | Canzar 2016 | Rejected on purpose: no abstain option. |

## Not in the code (none moves an objective)

1. **%LRC** (LRGASP): fraction of a transcript model covered by aligned reads. Cheap to compute. It would only tighten assembly-level claims, and the chain-rule search is closed. Optional reporting metric.
2. **BUSCO completeness**: reference-free completeness proxy. Not needed under the "no genome-only discovery" scope rule.
3. **Coverage-power calculation** (longcallR: power above 80% only over about 50× at SOR = 2): could justify the O2 abstain rate. A framing point, not a lever.

## Caveats

- Absence of a grep hit was read as "not applied" only after checking the register; a method under a different name could have been missed.
- Nothing here was prereg'd or measured; if any of the three items above is taken up, write a `docs/PREREG_*.md` first.

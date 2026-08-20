# NUMBERS — every quotable figure with its substrate

**Look a number up here before quoting it.** This file exists because provenance kept living in a
section header instead of next to the number, and that produced two wrong labels on the same document
in two days. A number without its substrate is not a result.

## Before you quote anything

1. **Name the species and the data.** ⚠ The two headline O1 rates are on **different species**:
   false-merge is HUMAN, false-omission is GORILLA. **Never pool them, never average them, never put
   them in one sentence without both substrates.**
2. **Name what the number IS.** A specificity is not a precision. A bound is not an estimate. A
   description of a known set is not a rate.
3. **Name the catalog, if it is a catalog property.** The shipped **494**-family catalog is *not*
   reproducible by the current binary, which emits **627** (`o1_catalog_provenance.md`).
4. **An arm is not the population.** A rule measured on the 164-pair FP/TP arms costs 0–6% there and
   3.67–12.80% genome-wide. Measure edge rules on the **whole edge set**.
5. **Report a count and a rate together, or neither.** On 2026-08-20 a certified-false edge *rate* rose
   0.2500 → 0.3333 while the *count* fell 7 → 1.
6. **Record `-p` and `-N` with any copy count.** A default `-p 0.8` silently discarded 8 of MAPKBP1's
   9 copies.

---

## O1 — family definition

| number | what it IS | substrate | notes |
|---|---|---|---|
| **2/150 = 1.33%** [0.0037, 0.0473] false-merge | **specificity, LOWER bound** — *not* a precision | **HUMAN** CHM13 v2.0 / `A119b.t2t.bam`; 150 gene-tight **single-locus windows** from 1,630 eligible, seed 101 | No positive stratum ⟹ no prevalence ⟹ **no precision**. Power demonstrated: 108/150 windows carry a merge opportunity over 2,466 node pairs. ⭐**Re-measured 2026-08-20 under the new defaults: unchanged, same two windows** (`o1_false_merge_remeasured.md`) |
| **28 → 3** spurious E_r edges (−89%); overlapping-span **26 → 1** | absolute counts on a negative panel | same panel as above | The orientation guard's benefit measured on an **independent** panel, not the GGO arm it was derived from |
| **9/162 = 5.6%** [0.0295, 0.1022] false-omission | rate, ARM 3 unbiased; upper bound 24.1% | **GORILLA** mGorGor1 + **matched fibroblast** IsoSeq; whole-genome excision of one copy from 162 two-copy families | ⚠ **Different species and design from the false-merge row. Never pool.** Not re-measured under the new defaults |
| **0/728** identity-clause failures | count | clause census | ⚠ **Cite the substrate LIST, never an ordinal** ("sixth independent substrate" is unaudited and used twice): 171/171 NPIP pairs · 194/194 records over 4,287 within-family pairs · 0/728 clause census |
| **22/40 = 0.5500** [0.3983, 0.6929] reach | recovery at τ=0.50; ~0.55 defensible genome-wide | **HUMAN** CHM13 **chr1**; 40-family GFF denominator fixed before the run | chr1 is representative under a **matched** clustering rule (40.0% vs 36.3% genome-wide, **Fisher p = 0.6090**). ⚠ The "chr1 is 80% clustered so it flatters the method" objection is **REFUTED** — the 80% and 36.3% came from *different rules*. Register entry; do not re-raise |
| **identical partition, 7/7** DNA vs RNA | agreement count | same loci, with RNA evidence | Small n. The RELATION transcends substrate; the NODE SET does not |
| **41/494 = 8.30%** definitional exposure ceiling | ceiling; **30 is a floor** | **GORILLA**, the **494**-family catalog | ⚠ Not re-measured on the 627 catalog. Needs the case census re-classified, not a recount |

### Structure — GORILLA, the **627**-family catalog, from the pipeline's own certificates (2026-08-20)

| | value |
|---|---:|
| 2-copy share — no split possible | **0.7018** |
| n ≥ 3 — the hierarchy ceiling | 0.2982 |
| complete graphs — γ provably inert | 0.0893 |
| **γ inert overall** (2-copy + complete) | **0.7911** |
| real reach (density < 1) | 0.2089 |
| λ = 1 with n ≥ 3 | 0.1786 |
| `cut_certified` (λ ≥ 2) | **75/627 = 0.1196** |

⭐ These agree with the offline re-derivation over the 494 catalog to **±0.019**, across a different
catalog, transcript set, code path and method — so the structural claims are substrate-robust.

## O2 — copy assignment

| number | what it IS | notes |
|---|---|---|
| **21.75%** near-ties | the **target population** | ⚠ **NOT MAPQ-0**, which is **0.0004** inside the multi-copy loci — ~500× smaller |
| **TPR 0.5066 / FPR 0.0280, AUC 0.7995** | held-out abstention | vs **MAPQ at chance, AUC 0.4944** — minimap2 is confidently wrong and its confidence is uninformative |
| **98.4%** restatement of the primary flag | 30/5,378 novel | ⚠ **Never claim "assigns better than minimap2"** — net headroom ~**0.1%**. Defend O2 on **abstention** |

## O3 — reference-absent copies

| number | what it IS | notes |
|---|---|---|
| **M ≤ 6.4** | bound, **unique**-sequence stratum | the route that works |
| **π = 1/35 = 0.0286**, 0/26 at cov ≥ 0.8 | **formally unbounded ⟹ vacuous** | paralogous stratum, unmapped-read route — **and O3's target class lives here** |
| **TPR 0.2703 / FPR 0.0200** | held-out S2 detector | set by **divergence not abundance**: 0.4500 above 0.01 vs **0.0588** below, and **45.78%** of positives lie below 0.01 |
| **ABSORBED 64.2% / ORPHANED 33.3%**, depth **1.75×** | excision fates | the signature is **UNMAPPED READS, not clipping** |
| ~~8/9/9 vs 5/6/8~~ | ⚠⚠ **DO NOT QUOTE** | does not reproduce at either `-p`; deficits shrink 3,3,1 → ~1,2,0 (`o3_haplotype_cnv_result.md`) |

## Operational

| | |
|---|---|
| genome-wide GGO catalog | **~2h20m, ~25 GB peak**, 627 families. State `D` + swap is **normal**, not failure |
| test baseline | **795 passed / 0 failed / 11 ignored**, `cargo test --release --all-targets` |

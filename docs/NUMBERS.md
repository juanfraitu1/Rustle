# NUMBERS — every quotable figure with its substrate

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the
> **superseded** 2026-07-17 catalog, which **no invocation of the current binary reproduces** — it was
> built with `refine` on, a default removed on 2026-08-20. The current default emits **627 families /
> 2,019 copies**. Re-measure before quoting: see [`NUMBERS.md`](NUMBERS.md) and
> [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

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
| **2/150 = 1.33%** [0.0037, 0.0473] false-merge | **specificity, LOWER bound** — *not* a precision | **HUMAN** CHM13 v2.0 / `A119b.t2t.bam`; 150 gene-tight **single-locus windows** from 1,630 eligible, seed 101 | No positive stratum ⟹ no prevalence ⟹ **no precision**. Power demonstrated: 108/150 windows carry a merge opportunity over 2,466 node pairs. ⭐**Re-measured 2026-08-20 under the new defaults: unchanged, same two windows** (`o1_investigations.md#false-positive-hardening-rules-that-survived-falsification`) |
| **28 → 3** spurious E_r edges (−89%); overlapping-span **26 → 1** | absolute counts on a negative panel | same panel as above | The orientation guard's benefit measured on an **independent** panel, not the GGO arm it was derived from |
| **9/162 = 5.6%** [0.0295, 0.1022] false-omission | rate, ARM 3 unbiased; upper bound 24.1% | **GORILLA** mGorGor1 + **matched fibroblast** IsoSeq; whole-genome excision of one copy from 162 two-copy families | ⚠ **Different species and design from the false-merge row. Never pool.** Not re-measured under the new defaults |
| **0/728** identity-clause failures | count | clause census | ⚠ **Cite the substrate LIST, never an ordinal** ("sixth independent substrate" is unaudited and used twice): 171/171 NPIP pairs · 194/194 records over 4,287 within-family pairs · 0/728 clause census |
| **22/40 = 0.5500** [0.3983, 0.6929] reach | recovery at τ=0.50; ~0.55 defensible genome-wide | **HUMAN** CHM13 **chr1**; 40-family GFF denominator fixed before the run | chr1 is representative under a **matched** clustering rule (40.0% vs 36.3% genome-wide, **Fisher p = 0.6090**). ⚠ The "chr1 is 80% clustered so it flatters the method" objection is **REFUTED** — the 80% and 36.3% came from *different rules*. Register entry; do not re-raise |
| **identical partition, 7/7** DNA vs RNA | agreement count | same loci, with RNA evidence | Small n. The RELATION transcends substrate; the NODE SET does not |
| **41/494 = 8.30%** definitional exposure ceiling | ceiling; **30 is a floor** | **GORILLA**, the **494**-family catalog | ⚠ Not re-measured on the 627 catalog. Needs the case census re-classified, not a recount |
| **275 loci / 120 of 627 families = 19.14%** | O1's **defensible miss** — loci with exon homology above `E_r`'s own floor AND primary read support | **GORILLA**, 627 catalog; `RUSTLE_ER_EDGE_DUMP` on a determinism-gated re-run (copies byte-identical to 2026-08-20) | Decomposes **0.8067 node construction / 0.1919 `E_r` / 0.0014 γ**. ⚠ **T8** — offline re-derivation, not pipeline-confirmed. Unit is **merged loci**, from 714 intervals |
| **107/17,924 = 0.60%** reps with `degree > 0` that never ship | γ's total cost, rep unit | same run | 96 reachable, 2,942,231 bp = **0.083% of genome**. ⭐ Replaces the retracted "γ costs 1/3,928" |
| ~~**45.28%** zero-primary but >=3 real alignments~~ ⛔**RETRACTED 2026-08-22** | **FAILED ITS MATCHED CONTROL** | — | Candidate 109/882 = 0.1236 vs size/compartment-matched control 104/881 = 0.1180, ⛔**REFUTED 08-23 on a VOLUME-INSENSITIVE test.** Per SECONDARY read, `de` at the locus vs `de` at that read's OWN primary: median per-locus fraction preferring the locus is **0.0000** in BOTH arms; pooled cand **233/108,688 = 0.0021** vs ctrl **151/37,559 = 0.0040**; **MWU p = 0.2685**. Both arms sit far below the genome-wide **0.0196**. The reads are SPILLOVER. **C5 stays OFF; O1 ⊥ O2 stands.** ⚠the earlier p = 0.7701 gate control was itself defective (a `de<=0.05` filter dropped 66% of secondaries) — quote THIS test, not that one. See `o1_ledger.md` §4i |

### ⚠⚠ RETRACTED 2026-08-21 — the SD "node gap" (audit: 17 agents, 0 of 4 angles survived)

| do NOT quote | why |
|---|---|
| **0.8984 "never a node"** as a gap | it **IS the base rate** for catalog-absent SD sequence (matched comparator 0.9086, n = 73,324; cleanest 0.9092, n = 69,396). The set is if anything node-*enriched* |
| **"γ costs 1/3,928"** | tie-break-degenerate (0 or 1 by tie policy); p flips across four nulls (0.0588 / 0.0274 / 0.2338 / 0.6517). Use **107/17,924** instead |
| **n = 3,928** as a sample size | 78.77% self-overlap, 49.97% strictly nested ⟹ **effective n = 1,556**; block-unit rate 0.8464 not 0.8984 |
| **~282** or **~1,374** expressed loci absent | scaled estimates on a denominator conditioned on the prediction. Use **275 loci / 120 families** |
| any **depletion factor** | **sign flips with the unit** — interval OR 1.32 (enriched), block OR 0.71 (depleted). Not an identified effect |
| **"398 nodes"** | 398 is **intervals**; the rep count is **248** |

⚠ **The truth set cannot measure O1 recall at all**: 21.89% of those intervals have *zero* bases
homologous to the anchor copy's exons, and 53.39% fall below the **0.50 coverage floor `E_r` itself
requires**. An SD partner is a duplicated **SEGMENT**, not a duplicated **GENE**.

### ⛔⛔ RETRACTED 2026-08-23 — the O1 loss decomposition (sign reversed)

| do NOT quote | why |
|---|---|
| **"γ costs 3.0× the coverage clause"** (246 = 0.2162 vs 82 = 0.0721) | scored on an edge set that is **NOT `E_r`** — `bench/o1_gamma_adjudicate.py:50` uses the min-length coverage form and **never reads the PAF strand field**, admitting minus-strand records the definition rejects at `denovo_pipeline.rs:4473` |
| **"the 114 questionable γ over-splits"** | **the set is EMPTY as a γ target**: 0/114 share an `E_r` component (Wilson [0.0000, 0.0326]); 113–114/114 have only a MINUS-STRAND qualifying record. γ never saw them |
| **"530 copies / 26.25% / 56 families"** at γ→0 | **466 copies / 38 families** (+14/4) = **480/2,019 = 0.2377** across **42/627 = 0.0670** |
| **"the partition does not transport, ARI 0.7064"** | **blob-driven** — drop the largest component and flat **ARI = 0.9707** (coarse 0.9578). ⭐ the **87.06% edge transport is unaffected** |

⭐ **CORRECT decomposition** — shipped `-k11 -w5` PAF, shipped coverage form, copy-pair unit, denom 1,135:
**orientation guard 167 = 0.1471 · coverage clause 147 = 0.1295 · identity 0 · γ 11 = 0.0097.**
⟹ the EDGE RULE costs **~28×** what γ costs. The orientation guard is the largest single loss channel — it buys precision (spurious edges **28 → 3**, antisense families **7.09% → 0.64%**, HUMAN panel) and that trade is now quantified on both sides. See `o1_ledger.md` §4m.

### ⚠ Anything quoted per **/1415** or **/494** is on the SUPERSEDED catalog

The 494-family catalog had **1,415** copies; the 627-family catalog has **2,019**. Re-measure before
quoting. Worked example — the short-copy population that the scale-free coverage mechanism exploits:

| | 494 catalog | **627 catalog** |
|---|---:|---:|
| copies with rep ≤ 2,000 bp | 352/1415 = **0.2488** | 432/2019 = **0.2140** |

The mechanism claim survives (still ~1 copy in 5), but **0.2488 is stale — quote 0.2140.**

### ⭐ Orientation guard, measured END-TO-END on the real binary (2026-08-20)

| | guard OFF (494) | **guard ON (627)** |
|---|---:|---:|
| families containing an overlapping same-family pair | 35/494 = 0.0709 (35/35 opposite-strand) | **4/627 = 0.0064** |
| spurious E_r edges, HUMAN 150-window negative panel | 28 | **3** |

Two independent substrates, both through the shipped binary. **No longer T8.**

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
| ~~8/9/9 vs 5/6/8~~ | ⚠⚠ **DO NOT QUOTE** | does not reproduce at either `-p`; deficits shrink 3,3,1 → ~1,2,0 (`o3_missing_copy_evidence.md`) |

## Operational

| | |
|---|---|
| genome-wide GGO catalog | **~2h20m, ~25 GB peak**, 627 families. State `D` + swap is **normal**, not failure |
| test baseline | **795 passed / 0 failed / 11 ignored**, `cargo test --release --all-targets` |

## Missing-copy detection (2026-08-26)

| number | what it is | the trap it avoids |
|---|---|---|
| **NPIP 3/31 vs 1/31 vs 12/31** | loci in a family: different animal / same animal 41% depth / same animal full depth | **depth beats individual identity** — the same animal at low depth loses to a different animal |
| **ABSORBED 64.2% / ORPHANED 33.3%** | fate of an excised copy's reads | ⚠ the 34.53% figure is a **READ** fraction, this is a **COPY** fraction |
| **113/162 = 69.8%** | excised copies with a homologous landing site | predicted from sequence alone; matches the 64.2% measured from reads, unfitted |
| **TPR 0.4248 / FPR 0.0000, AUC 0.8034** | depth caller, single genome, 40× sim | ⚠ threshold 1.5× is **NOT held out** |
| **0.7333 vs S2's 0.0588** | depth caller below 0.01 divergence | **12.5× on the stratum S2 cannot see**; Wilson95 [0.4805, 0.8910] |
| **15× is the knee** | coverage requirement | ⚠ **refutes my own ~25× estimate** — TPR is flat in coverage; depth buys PRECISION, not sensitivity |
| **456.91× / 1,389 Gbp** | WGS available under SAMN04003007 | ⚠ excludes the 62 RS II runs — those are 2015 **Y flow-sorted** DNA, a different experiment |


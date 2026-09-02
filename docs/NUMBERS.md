# NUMBERS — every quotable figure with its substrate

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the
> **superseded** 2026-07-17 catalog, which **no invocation of the current binary reproduces** — it was
> built with `refine` on, a default removed on 2026-08-20. The current default emits **627 families /
> 2,019 copies**. Re-measure before quoting: see [`NUMBERS.md`](NUMBERS.md) and
> [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

> ⚠⚠ **SECOND PROVENANCE BREAK — 2026-08-30. THE O1 NODE FLOOR CHANGED FROM 3 TO 2**
> (`NODE_MIN_READS`, ledger §6ac). **Every catalog figure recorded before 2026-08-30 — 627 families /
> 2,019 copies, 83 families / 484 copies, and every NPIP rate — was computed at a floor of 3 and the
> current binary does NOT reproduce them.** Set `RUSTLE_GATE_MIN_READS=3` to reproduce one.
> ⚠ `GATE_MIN_READS` is still 3 and still governs `copy_assign`'s tie-invariance bar — a different
> question (O1 ⊥ O2). Do not quote one for the other.

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
| **2/150 = 1.33%** [0.0037, 0.0473] false-merge — ⚠⚠**RE-MEASURED 2026-09-02 (§6bt.1). This figure is REPRODUCIBLE — it needs `RUSTLE_GATE_MIN_READS=3`, the node floor in force when it was measured.** At the shipped default (floor 2, §6ac) the same panel and binary give **3/150 = 2.00%** [0.0068, 0.0571] on a **DISJOINT** window set (W063/W106 → W033/W034/W065, no window in both). 2 vs 3 events is not a significant difference; the controlled swap is the finding, and **6/6 floor-2 edges carry a 2-read endpoint**. Quote the floor with the rate.** | **specificity, LOWER bound** — *not* a precision | **HUMAN** CHM13 v2.0 / `A119b.t2t.bam`; 150 gene-tight **single-locus windows** from 1,630 eligible, seed 101 | No positive stratum ⟹ no prevalence ⟹ **no precision**. Power demonstrated: 108/150 windows carry a merge opportunity over 2,466 node pairs. ⭐**Re-measured 2026-08-20 under the new defaults: unchanged, same two windows** (`o1_investigations.md#false-positive-hardening-rules-that-survived-falsification`) |
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

---

## ⚠⚠ STRATIFY RECALL BY SUBFAMILY — the coarse denominator mis-scores the rule (2026-08-28, ledger §6l)

**Any O1 recall figure quoted against a whole-family denominator penalises the rule for correctly
declining to link genes that are not recent duplicates.** On the 72-gene annotated DNA panel, the
"missing record" rate decomposes as:

| stratum | pairs with an alignment record | **missing** |
|---|---|---|
| **within-SUBFAMILY** | 188/201 = 0.9353 | **13/201 = 0.0647** |
| cross-subfamily, same coarse family | 198/605 = 0.3273 | 407/605 = 0.6727 |
| cross-coarse-family | 195/1,750 = 0.1114 | — |

⭐⭐ **The headline 420/806 = 0.5211 "no nucleotide record" is almost entirely a CROSS-SUBFAMILY
phenomenon. Within a true subfamily only 6.47% of pairs lack a record.** ⟹ **quote 6.47%, not 52%, as
the within-subfamily miss rate, and always name which stratum a recall number belongs to.**

⚠ This is the same failure mode already recorded for the dispersed-family stratum — *"the 8 were
mis-stratified domain superfamilies under a chr1-scoped rule"*, which moved the rate from 0/8 to
**1/24 = 4.2%**. A wrong stratum in the denominator has now cost two headline numbers.

⚠⚠ **The same edge set scores pair recall 0.2767 or 0.8209 depending only on the truth granularity**
(coarse product-family vs subfamily; component recall 0.6667 vs 0.9194, purity 1.0000 vs 0.9286).
**Granularity is a MODELLING CHOICE the data does not determine — every recall/purity figure must state
which one produced it.**

⭐ **Validate subfamily claims with GENOMIC POSITION, not product names.** NCBI product text is
similarity-derived ("golgin subfamily A member 6-**LIKE** protein 9"), and the protein channel alone
reproduces the DNA headline (ARI 0.9042, purity 1.0000), so a label-shuffle null tests against *no
structure*, not against the circular alternative. Position was never an input to minimap2: against
position-only tandem-array blocks the partition scores **ARI 0.6333 — exactly the annotated labels' own
0.6333** — beating protein 0.4073, identity 0.3930 and length 0.3770.

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

## Where O1 loses copies, and the seeded closure (2026-08-27)

| number | what it is | the trap it avoids |
|---|---|---|
| **30/31 vs 12/31** | oracle nodes vs real nodes, same edge rule | ⟹ the definition owns **3%** of the loss, node construction **58%** |
| **0.1840 vs 0.0580** | edge formation, single-exon vs ≥15-exon reps | completeness is PENALISED; survives a repeat attack at 3.67× |
| **50.5% / 49.0%** | random 30 kb genomic pairs that align / clear identity ≥0.60 | ⟹ identity never binds; the COVERAGE clause is what separates a repeat from a paralogue. ⭐**MECHANISM (§6m)**: minimap2 cannot emit below its scoring floor `B/(A+B)` = **0.667** default / **0.800** asm20, both ABOVE the 0.60 clause — **0 of 984,574 records across 8 PAFs fall below 0.60**, and the minimum tracks the PRESET (default 0.6313–0.7621, asm20 0.8291/0.8295) while the biology varies. The clause is inert as a discriminator on the E_r path. ⚠`tier2_rescue` uses `nm/bl` instead, where the same constant DOES fire (1.2–6.9%, min 0.1478) |
| **+0.0018 FPR for +0.4476 TPR** | lowering the coverage floor, on 87,990 real pairs | ⚠ **strong edge evidence that FAILS end-to-end** — NPIP stays at 12/31 |
| **25/31, converged round 3** | seeded closure, gorilla NPIP, ONE annotated seed | non-member hits constant at 1 every round ⟹ no repeat-chaining |
| **23/23 = 1.000** | expressed NPIP members found by that closure | vs 13/23 for the pipeline's own discovery |
| **65/65 converge, median recall 1.000** | human Soto panel, 65 families | ⚠ HUMAN — never pool with the gorilla figures |
| **SD-like 0.885 vs gene-like 0.895** | closure recall by family type | the SD-drift concern was well-founded and **does not manifest** |
| **r = −0.043** | closure recall vs family size, 65 families | large families are not harder; 14/16/17-member families all reach 1.000 |
| **8/248 = 0.0323 vs 6/256 = 0.0234** | canonical motifs among footprint gaps vs chance | Poisson p ≈ 0.24 ⟹ the old gate admitted footprints **at random**; fixed |



---

## 2026-08-30 — O1 recall, and two retractions

**Substrate for every figure below:** gorilla **fibroblast** `npip3.bam` (SAMN04003007 / KB3781, the
assembly animal), 1,009,396 primaries + 1,489,757 secondaries, 3 contigs. The **31 "NPIP loci" are a
minimap2 projection of HUMAN NPIP onto the gorilla genome** (`hsa2ggo.paf`), NOT gorilla annotation —
every rate below inherits that. `-p`/`-N` are the catalog defaults.

| configuration | families | copies | NPIP loci | pure families |
|---|---|---|---|---|
| floor 3, no flagfree (pre-08-30 default) | 83 | 484 | 12/31 | 3 |
| **floor 2** (shipped 08-30) | 121 | 678 | **14/31** | 3 |
| floor 2 + `FLAGFREE_SITES` (measured, NOT default) | 501 | 3,905 | **26/31** | **6** |

**⛔ RETRACTION 1 — "the remaining loci are not expressed."** That read a POST-ALIGNMENT PRIMARY COUNT as
transcription. Missing loci carry a **median of 404 secondary alignments vs 495 at recovered loci**; one
has 2 primaries and **1,119 secondaries**. A low primary count is exactly what read redistribution
predicts, so it can never evidence silence.

**⛔ RETRACTION 2 — the "23/31 expressed ceiling."** Derived the same way, and the flag-free arm reaches
**26/31**, past it. **There is currently no established expression ceiling for this panel.**

**Precision figures (size-robust — read these, not raw edge backing):**

| | pairs per NEW edge | density 2-5 / 6-15 / 16-39 / 40+ |
|---|---|---|
| floor 2 vs floor 3 | 1.5 | 0.81 / 0.44 / 0.34 / 0.38 |
| flagfree vs floor 2 | 2.0 | **0.96 / 0.60 / 0.50 / 0.53** |
| `COLLAPSE_EXONIC` (**refuted**) | **6.5** | falling |

⚠ **Raw direct-edge backing is CONFOUNDED WITH FAMILY SIZE** (pairs ~n², edges ~n) and produced a wrong
"do not ship" verdict for floor 2. Quote amplification and size-matched density instead
(`bench/merge_precision_arms.py`).

**Refuted 08-30, with the number that did it:** `COLLAPSE_EXONIC` (6.5 pairs/new edge, 62% of new pairs in
one family) · shared-read `E_c` edges as a definition tier (**fusion fraction ~72% INVARIANT from 2 to 100
shared reads** ⟹ no threshold exists).


## 2026-08-30 (later) — overlapping multimappers as site evidence

**Rule:** don't seed candidate sites from the primary flag (a coin toss between near-identical copies);
let every alignment propose a site, and require that **an unspliced site have a primary** — a spliced site
needs no such guard, because agreeing on an exact intron chain is already corroboration.

Same substrate as above (`npip3.bam`; the 31 NPIP loci are a HUMAN projection):

| | NPIP loci | largest family | single-copy controls (n=3) |
|---|---|---|---|
| primary-seeded (shipped default) | 14/31 | 54 | 3/3 PASS |
| every placement, no guard | 26/31 | 312 | ATP5F1A **FAIL** (115-copy family) |
| every placement + guard | **26/31** | 50 | 3/3 PASS |

⚠ `RUSTLE_FLAGFREE_SITES` is **opt-in, not the default**; the guard is intrinsic to it. Quote 26/31 only
with the n=3 caveat on the control and the note that the NPIP panel is a human projection.


## 2026-09-01 — reads that span two catalogued copies

**Substrate:** gorilla fibroblast `npip3.bam` (SAMN04003007/KB3781), one region query over
`NC_073244.2:69,492,851-69,524,989`. No pipeline run. The pair is GWFAM118:0 (4 exons, 275 reads) and
GWFAM111:3 (1 exon, 2,096 bp, 11 reads, `stub=true`), separated by 16,341 bp.

| measure | value | denominator |
|---|---|---|
| bridging N gaps on a CANONICAL motif | **329 / 329 = 100%** | all bridging gaps in spanning reads |
| distinct junctions canonical (all GT..AG) | **19 / 19** | distinct (donor, acceptor) pairs |
| spanning reads with aligned bases inside the second copy | **95 / 95 = 100%** | reads spanning both |
| spanning reads sharing the SAME last acceptor (69,522,690) | **95 / 95 = 100%** | reads spanning both |
| that acceptor's distance upstream of the second copy's start | **203 bp** | — |

⟹**GWFAM111:3 is the TERMINAL EXON of the GWFAM118:0 gene, not a copy.** 93 of the 99
`readthrough_span` rows are this one artifact.

### ⭐ THE DISCRIMINATOR — how to diagnose ANY read spanning two catalogued copies

Reusable, and it is the junction structure that decides — not the coverage, not the identity:

| observation | meaning | whose problem |
|---|---|---|
| **all spanning reads share ONE canonical junction**, and align inside both | one gene; the second "copy" is an **exon** of it | **O1** — a stub wrongly promoted to a copy |
| breakpoints **scatter** across reads | **mis-chaining** | O1 — the mis-chain filter's target |
| **no junction at all**, PSV alleles uniform along the read | **tandem-array mis-alignment**; the molecule came from ONE copy | O2 — assign or abstain |
| junction present **and PSV alleles SWITCH partway** | genuine **readthrough / fusion** of two real copies | biology — keep both, disclose |

⚠ measured on **n = 1 locus**. The distribution across the 71 engulfed copies is UNKNOWN.


## 2026-09-01 — the engulfment arms (§6as)

**Substrate:** `npip_cat/npip3.bam`, catalog `arm_f2` (node floor 2). Each arm ~17–28 min, 17.2 GB peak,
run strictly serially. A0 (OFF, new binary) is **byte-identical** to `arm_f2` on all four artifacts
(`cat.copies.tsv` md5 `9849dcb45b63e48e7b9b4d4358113a10`).

| arm | copies | families | NPIP | SharedAcrossFam | strictly engulfed |
|---|---|---|---|---|---|
| baseline `arm_f2` | 678 | 121 | 14/31 | 98 | 60 |
| `RUSTLE_READ_STRAND=1` | 634 | 125 | **15/31** | 34 | 20 |
| `RUSTLE_COLLAPSE_UNSTRANDED=1` | 497 | **87** | 14/31 | **5** | **1** |

⛔**NEITHER SHIPPED** — both fired pre-registered kill criteria. The second loses **34 families (25 losing
every copy span)** and **897 net E_r edges**; the largest family fragments. **A target hit at the cost of a
third of the catalog is not a fix.**

**Root cause, measured:** 58 of 60 engulfments are `('+','-')`, **all 1,888 single-exon nodes carry the
placeholder `'+'`** (`strand.unwrap_or('+')`), so `strand_conflict` is the only clause blocking the
collapse that would absorb the stub.

⚠**do NOT quote the split-only certificate as seed-invariance evidence** — it is a theorem of the code and
cannot fail (register, §6as). ⚠the family-purity/housekeeping control was **absent from the arm scorer**
and **fires** against the `READ_STRAND` arm; re-run every arm with it before revisiting.


## 2026-09-01 — non-canonical junctions and `RUSTLE_JUNCTION_MAJORITY`

**The advisor's objection is CORRECT.** Three independent lines of evidence; all figures on
`npip_cat/npip3.bam` (fibroblast KB3781, 3 contigs) unless stated.

**1. Depth — the motif test does real work at low support, and only there.**

| min reads | canonical | non-canonical | nc share |
|---|---|---|---|
| 1 | 47,640 | 18,844 | 28.34% |
| 10 | 16,586 | 207 | 1.23% |
| 500 | 3,462 | 7 | **0.20%** |

Non-canonical are **8.03% of DISTINCT junctions but 0.51% of OBSERVATIONS** — shallow, not massive.

**2. But the deep tail is real.** 18 non-canonical junctions at ≥100 reads → **7 distinct sites** (±30 bp);
**3 have NO canonical junction within ±30 bp**. ⚠the other 4 ARE jitter off a canonical site (one hotspot
supplies 12 of the 18 coordinates) — **a motif filter would be partly right for the WRONG reason.**

**3. Cross-substrate recurrence** — `GGO_ds.bam`, testis, **OR6737: different ANIMAL and TISSUE**:

| site | motif | fibroblast | testis | share of covering reads |
|---|---|---|---|---|
| `NC_073244.2:59965384-59967477` | CT..AT | 2,407 | 220 | **87.6%** |
| `NC_073241.2:52758208-52758462` | GT..GG | 462 | 65 | **92.9%** |
| `NC_073244.2:60228246-60228362` | CT..CC | 128 | 1 | 1.4% |

⭐**two are the DOMINANT form at their locus in a second animal.** `CT..AT` is the NPIPB12 motif class —
the rule already discards a 109-read NPIPB12 model over one such junction.

**`RUSTLE_JUNCTION_MAJORITY` arm** (OFF byte-identical, `copies.tsv` 9849dcb4):

| | OFF | MAJORITY |
|---|---|---|
| families / copies | 121 / 678 | 117 / 700 |
| NPIP loci | 14/31 | **14/31 (no gain)** |
| strictly engulfed | 60 | 63 |
| housekeeping control | 3/3 | **3/3** |
| the 2 lost sites | 0 nodes | **1 node each** |

⛔**STILL OFF BY DEFAULT. The blocker is named: the chr16 harm (copies 66→34, families 20→11) was measured
ON chr16, and this arm ran on NC_073241/242/244 — THE FAILURE IS NOT REACHABLE THERE, so a clean pass is
not a green light. A chr16 arm is required before any flip.**

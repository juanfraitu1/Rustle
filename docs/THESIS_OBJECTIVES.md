# Thesis objectives — the current three

The scope agreed with the advisor. **This file is the live status; nothing here is superseded.**
The 2026-08-07 five-objective VG/EM decomposition this file used to carry is archived in
[`o1_investigations.md`](o1_investigations.md#superseded-five-objective-vg-em-decomposition).


| # | objective |
|---|---|
| **O1** | Define a multi-copy gene family topologically at the RNA level (quasi-clique in E_r; MCC = χ(H)) |
| **O2** | Decide, for a read at a multi-copy locus, whether the evidence warrants assigning it to a COPY AT ALL — and abstain when it does not. Contested set = **alignment-score near-ties**, not MAPQ-0. No 1/k. |
| **O3** | Detect + flag expressed transcript paths not explained by represented reference copies, **STRATIFIED by whether the orphaned reads have anywhere to go**. Detect-and-flag with a measured FPR per stratum; **completeness is never claimed**. (this was O4) |

**The allele-specific-junction objective is DROPPED** — cut for time, and because it does not connect
to the others. ⚠ It was the only objective the 2026-06-25 audit rated ATTAINED; the ASJ result itself
is de-scoped, not retracted.

⚠ **"O3" denotes three different things across this repo's history** — the EM below, ASJ in the
06-25 audit, and reference-absent copies now. Resolve the scheme before quoting any objective number.
The numbered sections below are kept for provenance only.

### ⭐ O3 RESTATED (2026-08-19) — stratify the target, bound each stratum

O3 as originally posed — *"find copies missing from the genome"* — is **not achievable**, and the
reason is measured, not conjectural. In the whole-genome excision control (one copy of 162 two-copy
families deleted, matched IsoSeq) a deleted copy has **two fates**:

| fate | rate | what the reads do |
|---|---:|---|
| **ORPHANED** | 33.3% | median **92.7% of reads unmapped** — detectable |
| **ABSORBED** | 64.2% | reads land on the best paralogue at **1.75× depth**, concordance 0.967 — invisible to any unmapped-read method |

**Expression is not the constraint** (99.34% of copies clear the floor); **where the reads go** is.
⚠ And O3's original target class — a collapsed paralogue — sits in the ABSORBED stratum.

This is the same move that rescued O2: restate the target population and the claim, then bound each
stratum honestly.

**Target.** Copies absent from the assembly, **stratified by whether the orphaned reads have anywhere
to go.** **Claim.** Detect-and-flag with a measured FPR, stated **per stratum**:

| stratum × route | status | bound |
|---|---|---|
| unique sequence, **unmapped-read** route | **works** | **M ≤ 6.4** missing expressed copies |
| paralogous sequence, **unmapped-read** route | ⚠⚠ **vacuous** | π = 1/35 = 0.0286, 0/26 at cov ≥ 0.8, **formally unbounded** |
| paralogous sequence, **depth (S2)** route | **partial** | held-out **TPR 0.2703 / FPR 0.0200**; sensitivity set by **divergence, not abundance** (0.4500 above 0.01 divergence vs **0.0588** below) — ⚠ and **45.78% of positives lie below 0.01** |

**Say this, not more:** *the instrument flags candidates with a measured false-positive rate in a
named stratum, and has explicitly no power for unmapped-read detection in the collapsible stratum.*
**Never claim completeness.** The signature is **UNMAPPED READS, not clipping** (34.53% pooled,
MAPQ-60 before deletion) — no published collapse detector uses clipping.

⚠ **The one real candidate does not close:** STON1+GTF2A1L, ~116.7 kb absent from mGorGor1, 125
near-full-length unmapped reads, gapless chromosome, 0 GFF lines, present in chimp and orangutan —
but **single-copy, n = 1, p = 0.055, UNCONFIRMED**. It supports the instrument; it is not a result.

⭐ **The niche is empty**, which is the thesis value: nobody has found a reference-absent copy from
transcriptome data. The field standard is S1 re-assemble / S2 depth+PSV / S3 peptides. A *bounded
negative* is therefore itself publishable.

⚠ **The advisor's "the reference is an average" objection fails on PROVENANCE, not on measurement:**
mGorGor1 is a haplotype-resolved assembly of **one animal** and the fibroblast IsoSeq is that
animal's **own cell line**. Keep the two arguments separate — the objection dies on how the substrate
was built, and the copy-number-polymorphism rate is **still unmeasured** (the 2026-08-19 mat-vs-pat
run was uninformative: control floor 0.1512 vs signal 0.0278). ⚠⚠ And **do not quote the 8/9/9-vs-5/6/8
haplotype deficit** — it does not reproduce; see [`o3_missing_copy_evidence.md`](o3_missing_copy_evidence.md).

**Current O3 avenues and claim boundary (2026-08-16):** see
[`docs/o3_missing_copy_evidence.md` §8](o3_missing_copy_evidence.md#8-possible-o3-avenues--decision-record-2026-08-16).
Liftoff or a second genome is not required for the main experiment; natural RNA-only findings remain
candidates unless independently validated with donor DNA.

**Current O1 purity rules and expanded known-family graphs (2026-08-16):** see
[`docs/o1_investigations.md#false-positive-hardening-rules-that-survived-falsification`](o1_investigations.md#false-positive-hardening-rules-that-survived-falsification) and the 19-graph
[`expanded audit`](../bench/o1_expanded_family_audit/README.md). Soto SD membership is discovery
evidence, not automatic gene-family membership; primary and audit graphs are emitted separately.

# ⭐ SCOREBOARD — what is WON, what is DEAD, what is OPEN

> **Read this before proposing anything.** Every row is measured unless marked *estimated*. Numbers carry
> their unit and, where it matters, their comparator. `NUMBERS.md` holds provenance; the 631-entry
> `NEGATIVE_RESULTS_REGISTER.md` holds the long tail; `o1_ledger.md` §3–§4l holds O1's route-by-route detail.
> Last updated **2026-08-23**.

## ⚠⚠ THE META-LEDGER — the failure modes that keep recurring

These have cost more than any single wrong idea. **Four headlines died in one week to #1 alone.**

| # | failure mode | the scar |
|---|---|---|
| **1** | **A denominator conditioned on the prediction.** | Killed 7 findings, then 4 headlines in one week: the SD "node gap" (0.8984 *was* the base rate 0.9086); a degree-0 rate of 0.9947 against a 0.9945 base rate; a 93.06% self-overlap that inverts on its control; the primary-flag finding. **A rate without its comparator is a stop condition, not a to-do.** |
| **2** | **A CONTROL can be broken too.** | The gate control that "refuted" the primary-flag finding filtered secondaries at `de ≤ 0.05`, silently discarding **66%** of the evidence and preferentially the diverged reads that mattered. **Check what a control's filters remove.** |
| **3** | A null matched on **COUNT** proves nothing — match the **SIZE distribution** (and, for haplotype work, **COMPOSITION**). | The 2026-08-19 hapcnv run died on span- vs composition-matched controls (floor 0.1512 vs signal 0.0278). |
| **4** | Judging a change to **what a NODE or EDGE IS** on node-level metrics. | Failed **3× end-to-end**. Judge on the emitted PARTITION. |
| **5** | **T8** — offline re-derivation is a hypothesis generator, **never a test**. | 4 prior offline-proxy errors on the E_r/PSV line alone. |
| **6** | **argmax / best-hit counting is not counting** (T19); always record **`-p` and `-N`** with any copy count. | MAPKBP1 gives 1/1 at `-p 0.8` and 9/8 at `-p 0.1`. |
| **7** | Never pool **human** and **gorilla**; state the **unit** every time. | The two headline O1 rates are on different species. |
| **8** | **Substrate provenance.** | The shipped catalog was built from `GGO_ds.bam` = **OR6737 testis**, a *different animal* from the assembly (KB3781). Same-animal RNA exists: `fibroblasts/GCA_029281585.2_flnc_mm.bam`. |

---

## O1 — define a multi-copy gene family topologically at the RNA level

### ✅ WON
| result | number | substrate |
|---|---|---|
| Shipped catalog | **627 families / 2,019 copies** | GORILLA mGorGor1 |
| **P1 seed-invariance is a THEOREM** | — | definitional |
| **The RELATION transports across ANIMAL *and* TISSUE** ⭐ *(08-23, not T8)* | **87.06% of `E_r` edges recovered** (Jaccard 0.7397). Partition **ARI 0.9707** once the single largest component is dropped (0.7064 raw — **blob-driven**) | testis(OR6737) vs depth-matched fibroblast(KB3781), 10,143 shared loci |
| λ (edge connectivity) shipped as a per-family certificate | `lambda` / `cut_certified` | — |
| Orientation guard shipped as default | spurious edges **28 → 3**; antisense families **7.09% → 0.64%** | HUMAN negative panel + GORILLA |
| ONE path — `--refine` removed, hard error on O1 | — | — |
| False-**merge** | **2/150 = 1.33%** [0.37, 4.73] | ⚠ **HUMAN**, a *specificity LOWER bound*, NOT a precision |
| False-**omission** | **9/162 = 5.6%** [2.95, 10.22] | ⚠ **GORILLA**, excision ARM 3. Different species/design from the row above — **never pool** |
| Reach | **22/40 = 0.5500** [0.398, 0.693], ~0.55 defensible genome-wide | HUMAN chr1; chr1 representative, Fisher p = 0.6090 |
| Identity clause **never binds** | **0/728** | clause census |
| DNA and RNA give the **identical partition** on the same loci | **7/7** | small n |
| Per-rep **exon arrays** shipped ⭐ *(08-22)* | **45.06%** of "a node is here" was an **INTRON** | rep spans are 90.83% intron by bp |

### ⛔ DEAD — do not retry
| route | the number that killed it |
|---|---|
| **Fixing the coverage hole** | HERC2 (a real family) splits at `c_long ≥ 0.034` while the FIRST false positive dies only at 0.05 ⟹ **TP loss starts before FP rejection; no constraint-satisfying threshold exists** |
| Reads as a third source for it | tiling breadth significant in the **WRONG direction** (AUC 0.1259, p = 0.0005); reads are **redundant** with the rep alignment (median gain **+0.000**) |
| Catalog-anchored repeat gate (R5) | **breaks P1** — same pair accepted/rejected by universe (seed-invariance 94/147 moves) |
| Genome-anchored `min_shared_gmult` **as an admission rule** | rejects **12.80%** of shipped edges; parity with γ needs M≈3 and discards 48% of the edge set. **VETO only** |
| Graph theory — density / Fiedler / min-cut | **AUC 0.36 / 0.35 / 0.35** — all *below* chance |
| **Shared-exon / isoform pooling** *(08-22)* | **8.29%** `E_r` recovery, **25.04%** family retention against λ=1 in 88.04%; genuine new payload **9** node pairs at MIN_COUNT≥2, **0** at ≥3 |
| `RUSTLE_LOCUS_EXON_UNION` | **−20 recall points** — concatenation inflates its own denominator |
| `--single-copy-baseline` as a rep source | wrong object (`E_c` not `E_r`); **83.3%** of the catalog's own copies structurally unreachable |
| Facility location / bipartite matching for loci | **BANNED** — and bipartite assignment *is* the primary flag (99.10%) |
| FamilyGraph · k-mer Jaccard · `vg` for ties | dead / linear |


- **The shared-exon denominator repair** (`RUSTLE_SHARED_EXON`) — refuted §4w. Every endpoint moves the
  wrong way; the completeness deficit "improves" 3.16× → 1.32× only because the STUB rate collapses too.
  ⚠ a ratio of two rates is not an endpoint when both may fall.
- **The repeat-hub gate ported from `family_define`** (`RUSTLE_ER_REPEAT_GATE`) — refuted §4x, actively
  harmful: NPIP loci 12/31 → 7/31, pure families 3 → 1. **Genome multiplicity cannot distinguish a repeat
  from a high-copy family, and high-copy families are O1's subject.** Its apparent signal was tautological
  (mult ≥ 20 ⟹ edge rate **1.0000**).
- **Lowering the E_r coverage floor** — refuted end-to-end §5c, *despite* the strongest edge-level
  evidence O1 has produced (FPR +0.0018 for TPR +0.4476 against 87,990 real genomic pairs). NPIP recovery
  is UNCHANGED at 12/31 at every floor while pure families fall 3 → 1 and the largest balloons 39 → 104.

### 🔄 OPEN
- **The named definitional hole**: the min-length coverage denominator is **scale-free** (a ~1 kb repeat is ≥0.50 of ANY node < 2 kb; 24.88% of gorilla copies are ≤ 2 kb). Ceiling **41/494 = 8.30%** on the **superseded 494-family catalog**, ⚠ *not re-measured on the current 627*.
- ⛔ **NODE CONSTRUCTION — the stub defect is REAL but NOT FIXABLE by representative choice** *(08-25, ledger §4s)*. `pick_locus_rep` keeps ONE chain; **46%** of reps covering a known family member are single-exon stubs and **53%** of those loci have a discarded gate-passing spliced chain.
  - ⭐ **`E_r` has NO UPGRADE PATH: 0/215** — for stub-incident edges whose stub has a unique containing spliced model at its own locus, the container inherits the edge **0/215 = 0.0000** (validity checks: **1,416/1,416** containers are themselves OFF reps, and **116/215** have edges of their own). *"Edge is a property of the locus"* predicts 215/215; *"of the exact representative sequence"* predicts ~0.01. **Any representative change rewrites the edge set from scratch.**
  - ⛔ **`RUSTLE_SPLICED_REP` was measured end-to-end and REGRESSED**: chr7 F1 **0.570 → 0.411**, chr16 **0.910 → 0.761**, loses NPIPB12. The counterfactual on the strand arms is a null (annotated-intron precision 0.9700 vs 0.9714).
  - ⚠ **60.51%** of the edge set falls below the 0.50 coverage floor under a uniform **1.49×** rep-length inflation ⟹ **any intervention that lengthens a representative pays `RUSTLE_LOCUS_EXON_UNION`'s −20 recall points**, selected or constructed.
  - ✅ **What ships instead: a FLAG** — `frac_pure_intron` per copy against the gorilla-native GFF. Within-family paired effect **1.37×** (22/118 vs 16/118); burden 7.4–12.8% of copies. ⚠ quote 1.37×, **not** the unclustered 13.55×; ⚠ the GFF holds **0 NOTCH2NL\*, 1 SRGAP2, 1 LIMS1** — structurally blind at the flagships.
- ⭐⭐ **THE ORIENTATION GUARD IS THE LARGEST LOSS CHANNEL** *(08-23, ledger §4m)*. On the shipped `-k11 -w5` PAF with the shipped coverage form (copy-pair unit, denom 1,135): **orientation guard 167 = 0.1471 · coverage clause 147 = 0.1295 · identity 0 · γ 11 = 0.0097** ⟹ the edge rule costs **~28× γ**.
  - ⛔ **RETRACTED:** "γ costs 3.0× the coverage clause" — scored on a **non-`E_r`** edge set. ⛔ The **"114 questionable γ over-splits" set is EMPTY**: 0/114 share an `E_r` component; 113–114/114 have only a minus-strand record, so γ never saw them.
  - ⛔ **A HIERARCHY IS A NO-GO** *(§4m)*: the level that recovers the γ-cut pairs **is** the level that re-forms the blob; new coarse pairs are duplication-supported at **0.0016** vs **0.1252** — **78× less**. What ships is a *disclosure*, not a hierarchy.
- ⭐⭐⭐ **THE GUARD WAS FILTERING ON AN UNMEASURED FIELD — RESOLVED** *(08-24/25, ledger §4n–§4r)*. `strand.unwrap_or('+')` gave **all 5,928 single-exon reps `'+'` and zero `'-'`** (spliced split 0.4867/0.5133), so a third of the rep set carried a **placeholder**. **98.55%** of guard-blocked pairs involve such a rep (base 0.3943, **2.50×**), leaving only **58** genuine antisense candidates genome-wide. Read orientation measures the strand at **0.9650** (0.9934 at margin 0.90, 1.0000 at unanimity) against a constant-`'+'` floor of **0.4867**.
  - Shipped **opt-in** (`RUSTLE_READ_STRAND`, `RUSTLE_READ_STRAND_MARGIN`); the OFF gate is **byte-identical** to the shipped catalog. All pre-registered criteria PASS (0/58 antisense admitted; antisense-family rate 0.0048 both arms; **HUMAN 150-window false-merge 2/150 = 0.0133 UNMOVED**).
  - ⛔ **ABSTENTION DOES NOT RESCUE IT** *(§4q)*: the dominant dissolution cause is the **same 32 families in both arms**, so it is not vote-related. **Line CLOSED after four iterations, stop pre-registered.**
  - ⭐⭐ **THE REAL EFFECT IS NODE CONSTRUCTION** *(§4r)*: measuring the strand **re-picks 2,306 representatives and turns 1,928 = 0.8361 of them from a single-exon STUB into a SPLICED model** — the largest named node-construction defect — at a cost of 32 dissolved families. ⚠ That a spliced rep is *better* is **inferred, not measured**.
  - ⭐ **WHAT IT SETTLES:** `E_r` edges are computed on the **representative's sequence**, so improving a node breaks the edges formed with its worse version. **That is a property of ONE-REPRESENTATIVE-PER-LOCUS**, not of the strand fix, and any future change to representative selection pays the same cost.
  - ⚠⚠ **A matching lesson worth more than the result**: reading a property off *"the rep with the largest overlap"* returns a **neighbouring** rep when the exact span is gone — the case under test. That produced a phantom class of 18 families, retracted **twice**. **Match loci EXACTLY.**
- ⚠ The blob is **466 copies / 38 families** (+14/4 in a second) = **480/2,019 = 0.2377** across **42/627 = 0.0670**. γ is **non-monotone under edge addition**.

---

- **Node construction, not the edge rule.** Measured at NPIP (§5c/§5d): the edge rule has **1 locus** of
  headroom, node construction **4**, read fragmentation **≤ 7**, and **7 loci are not transcribed at all**.
  ⚠ §5c first claimed 18 behind node construction; that used reads overlapping by ≥ 1 bp and is
  **withdrawn** — requiring reads ≥ 50% inside gives 23/31 expressed, not 28/31, and 419 reads, not 5,544.
- **The completeness deficit itself has NO known repair** (§4v, survives a repeat attack at 3.67×). The
  denominator route is measured and dead; the failure surface is a hard cliff at exactly 50% shared
  sequence, divergence-independent (§5b).

## O2 — assign a read to a copy, or abstain

### ✅ WON
| result | number |
|---|---|
| **The abstention defence** — first non-circular validation (excision, truth by design) | held-out **TPR 0.5066 / FPR 0.0280, AUC 0.7995** |
| ⭐⭐ **MAPQ is AT CHANCE on the same reads** | **AUC 0.4944**, median 60 vs 60 ⟹ *minimap2 is confidently wrong and its confidence is uninformative* |
| Target population **restated** | contested set = **alignment-score near-ties (21.75%)**, not MAPQ-0 (**0.0004** inside multi-copy loci) — ~500× larger |
| PSV engine ships with a **derived** bound, not a tuned threshold | gate `n_decisive ≥ 1`; **τ = ln((1−p)/p)**, τ=6.9 ⟹ stated **P(misassign) ≤ 1e-3** |
| Given the copy set, **assignment DECOMPOSES** ⟹ per-read argmax is optimal | — |
| **O1 → O2 is a CONDITIONING, not a correlation** ⭐ | `copy_assign.rs` has **0 hits** for λ/min-cut/density/γ; it consumes only **cardinality**, and n→2 moves **14/104,147 = 0.013%** of verdicts |

### ⛔ DEAD — do not retry
| route | the number |
|---|---|
| **Reassignment as the headline** | a secondary fits better by `de` in **1.96%** of reads-with-secondary; net headroom **~0.1%**. **NEVER claim "assigns better than minimap2"** |
| The EM | changes **ZERO** evidenced decisions (**0/3,081**) |
| Swapping votes → likelihood | **identical in 16/16** configs (monotone-equivalent under flat error). *The lever is the GATE* |
| A better aligner (winnowmap) | **identical** to minimap2; DSFAM42 stays 95% MAPQ-0 under both — the limit is **identifiability** |
| λ / density / min-cut / `de`-weighted capacity as assignability predictors | all fail; capacity **equals `min_de` on 82.14%** |
| Structural-anchor reassignment truth | **cannot scale** — filter L1 rejects **29/29** anchors; **0/53** scored reads are near-ties |
| Secondary-read corroboration / C5 admission | median per-locus preference **0.0000 in BOTH arms**, MWU **p = 0.2685**; both ~5–9× *below* the genome-wide 0.0196 ⟹ **spillover**. **C5 stays OFF** |

### 🔄 OPEN
- `n_decisive = 0` (exonically identical copies) is a **genuine, aligner-invariant floor** — abstention there is correct, not a gap.
- The abstention defence rests on **one** excision panel; a second independent validation would strengthen it.

---

## O3 — detect + flag reference-absent copies (**completeness never claimed**)

### ✅ WON
| result | number |
|---|---|
| A positive control exists at all | whole-genome excision, **162** two-copy families, truth **by design** |
| **The two fates, measured** | **ABSORBED 64.2%** (concordance 0.967, depth **1.75×**) · **ORPHANED 33.3%** (median **92.7%** of reads unmapped) |
| ⟹ the signature is **UNMAPPED READS, not clipping** | no collapse detector in the literature uses clipping |
| S2 detector | held-out **TPR 0.2703 / FPR 0.0200** |
| Sensitivity is set by **DIVERGENCE, not abundance** | **0.4500** at ≥0.01 diverged vs **0.0588** below |
| DNA-side instrument validated | FP floor **0/817**; fires on **11.42%** of random intervals; recovers a published expansion exactly |
| Literature position established | **nobody** has found a reference-absent copy from transcriptome data |

- **A depth caller for the ABSORBED stratum** (§5a, simulation). TPR 0.4248 / FPR 0.0000, AUC 0.8034 —
  and **TPR 0.7333 below 0.01 divergence where the S2 detector scores 0.0588**, so the two are
  COMPLEMENTARY rather than competing. 15× is the coverage knee. ⚠ FPR 0.0000 is measured under uniform
  simulated coverage and **must not be quoted as a real-data figure**; the 1.5× threshold is not held out.

### ⛔ DEAD — do not retry
| route | the number |
|---|---|
| Unmapped reads in the **paralogous** stratum | **π = 1/35 = 0.0286**, **0/26** at cov ≥ 0.8 — **VACUOUS, and O3's target class lives there** (a collapsed copy is absorbed, not orphaned) |
| **The haplotype route** *(08-23)* | parameter stability **Jaccard 0.5263** (47% turnover on one flag); genes **DEPLETED** vs a composition-matched control (0.96% vs 2.12%, CMH p = 0.5749); a mirror arm with the *opposite* truth value is **symmetric** (56 vs 51, p = 0.6992) |
| …and its "larger compartment" premise | absent fraction is **0.3295% of genic span** ⟹ 10.9–42.1 Mb vs the existing ~37.0 Mb — **~1.1× at best, 3.4× SMALLER in the gene frame** |
| **O3 as a residual of O1+O2** | degenerates: the excision panel is **two-copy**, so masking leaves \|C\|=1 — no assignment, hence no abstention. The PSV residual at \|C\|=1 **is the shipped S2 detector renamed** |
| The **depth** residual | refuted — 16/104 destinations had **zero baseline reads**; 1.75× needs a "before" a real case does not have |
| DNA-side collapse detection at the current compartment | **0 collapses at 816/817**, but **underpowered BY CONSTRUCTION** (1.1224% of genome ⟹ 0.47–0.94 predicted) |

### 🔄 OPEN
- **One unconfirmed candidate**: STON1+GTF2A1L (~116.7 kb absent, 125 near-full-length unmapped reads) — **n = 1, p = 0.055**.
- ⚠ **45.78% of real positives sit BELOW the divergence line** where S2 is near-blind — a property of the problem, not the detector.
- An **n ≥ 3 excision panel** would unblock the residual test. ⚠ Feasibility: 14 families at ≥20 reads, 56 at ≥5, 157 at ≥3 — but only **6/43 = 13.95%** of scenarios put the deleted copy below de 0.01 against **45.78%** of real positives, so it would sample the **wrong divergence regime**. Measure that distribution first; it is free.
- A **larger compartment** remains the named lever for the DNA side.

---

## ⚠ DROPPED — allele-specific junctions (ASJ)
Cut 2026-08-07 for time and because it does not connect to the others. ⚠ It was the **only** objective the
2026-06-25 audit rated ATTAINED — say what was given up. The result is **de-scoped, not retracted**.

---

**Thesis claim:** assemble the isoforms of *individual paralogous copies* of multi-copy gene families from
long reads — recovering, quantifying, and resolving copies that primary-alignment assembly (StringTie) cannot.
The variation graph (`--vg`) shares evidence across copies; an EM assigns ambiguous reads to copies.

Status key: ✅ done · 🔄 partial / in progress · ⛔ blocked.

---

## O1 — Formally define a multi-copy gene family in the VG  ✅
**Foundational** — makes O2–O5 measurable (you cannot score copy recovery without defining a copy).

A multi-copy family is a **connected component of the read-multi-mapping graph** (a *splice-graph
generalization*: nodes = exon-classes carrying per-copy sequence, paths = copy × isoform) satisfying:
**M** ≥2 copies · **H** copies share a backbone (reads multi-map) · **X** ≥2 copies independently expressed
(edit-distance anchored) · **I** graded identifiability (full / partial / none), reported not assumed.

- Spec: `docs/superpowers/specs/2026-06-01-multicopy-family-definition-design.md` (jargon scrubbed for the advisor).
- Operational predicate `classify_family` shipped (`src/rustle/vg.rs`), emitted as `--vg` GTF attributes
  (`family_verdict`, `family_identifiability`, `family_n_copies/n_expressed`, `family_locus_rel`).
  Commits `d286400`..`6d91bf1`. Expression counted per **identifiability class** (`3609b60`), runs
  independent of EM eligibility (`6d91bf1`).
- Validated on **14 known families** (`bench/paralog_secondary_scan/validate_known_families.py`): DAZ / NBPF /
  RBMY / amylase → `family`; TSPY (97% ties) / MAGEA → `family_nonidentifiable`; SORD+LOC → spillover;
  β-defensin / protocadherin / PRAME correctly excluded (not expressed here).
- **Fresh O1 purity challenge (2026-08-16):** current Rustle was rerun on predeclared regional
  extracts from the original whole-genome-aligned GGO/HSA BAMs, without old node or family ids as
  inputs. It re-emits 124/133 audited loci and 72/75 independently named targets. Although 14/16
  unrelated conflicting-gene loci are real/reproducible emissions, only 1/16 rejoins the target;
  all nine re-emitted non-NBPF loci from the adversarial repeat bridge remain outside fresh NBPF.
  GOLGA2 is now separated as a documented broad-family/recent-subfamily outgroup rather than
  mis-scored as an unrelated false positive; an RNA identity-0.80 view removes it but damages MAGEA
  and NBPF, so this is a typed hierarchy rather than a new global threshold.
  The cost is explicit: 69/75 named targets land in the modal family (three not emitted, three
  split). The HSA run also discloses one node-pair decision delegated to O2 `reads_distinguish`,
  so this particular node set is not sequence-only. Direct Rustle tables, rule certificates,
  logs, and actual fresh `E_r` GFAs:
  `bench/o1_fresh_emission_validation/`; interpretation: `docs/o1_investigations.md#false-positive-hardening-rules-that-survived-falsification`.
- **Deferred implementation:** emit the current RNA homology family as the broad family plus an
  opt-in, nested DNA-supported recent-copy subfamily (`RECENT_COPY`, `BROAD_ONLY`, or
  `DNA_UNRESOLVED`). The annotation-free algorithm, flag-off byte-identity requirement, output
  schema, GOLGA discriminator, and cross-family safety tests are specified in
  `docs/o1_investigations.md#block-aware-duplication-provenance-graph`. Production Rustle does not yet emit these fields.
- **Provenance-model avenue:** represent loci as ordered paths through homologous duplication blocks,
  with separate RNA-homology, DNA-duplication, read-conflict, and optionally rooted ancestry edges.
  This can express the mosaic GOLGA2 + ITSN2-UTR origin of the chr15 GOLGA expansion without calling
  ITSN2 a GOLGA family member. A single genome yields an unrooted network; directional “ancestral”
  claims require outgroup sequence/synteny, without using the outgroup as the assembly reference.
  Formal model, five-family local pairwise-witness prototype, empirical concordance results, and
  proof-of-concept criteria: `docs/o1_investigations.md#block-aware-duplication-provenance-graph`. Durable typed tables and GFA
  projections: `bench/o1_provenance_witness_prototype/`. This is evidence that the representation
  separates coherent cores from repeat bridges; stable multi-locus block-class construction is
  still deferred and the current graphs remain `UNROOTED`. A deferred single-outgroup extension now
  specifies direct minimap2 block/flank alignment to both phased gorilla haplotypes, two-sided
  synteny rooting, explicit abstention states, and flag-off invariance. It uses no annotation
  projection and cannot change human family membership. Local ape assembly paths and WSL mount
  instructions: `docs/REFERENCE.md`. A GOLGA proof of concept now finds recurrent human
  intervals from GOLGA2 into 8 audited family loci and from ITSN2 into 6, while both proposed source
  loci have unique two-sided synteny in both gorilla haplotypes. These are retained as
  `ROOT_CANDIDATE_SINGLE_OUTGROUP` with `direction_status=UNROOTED`, because stable multi-locus
  block classes are not implemented. Evidence and rerun script:
  `bench/o1_outgroup_rooting_poc/`.
- **Open:** the §9 "report the range of consistent splits" quantification is deferred (shipped behavior:
  label non-identifiable + abstain from splitting); scope gate is `any_spliced` (lenient vs per-copy wording).

## O2 — Recover copies that primary-only assembly misses  🔄 capability shown · ⛔ no validated real-data recovery
Copies whose annotated transcripts are covered overwhelmingly by **secondary** alignments (their primary went
to a sibling paralog) — StringTie drops them; the VG *can instantiate a model* for them. No end-to-end run yet
demonstrates a dropped copy actually recovered on real data.

- Genome-wide scan (`bench/paralog_secondary_scan/`, GGO IsoSeq vs RefSeq): of 846 expressed multi-copy
  candidates, an edit-distance copy-anchor **classifies** ~93 expressed_real_copy (NBPF, LRPAP1-like, …) vs
  ~89 pure_spillover loci that naive secondary-use would *fabricate* (the honest counterexamples;
  SORD/LETM1/CES1 pseudogenes). This is *classification*, not assembly — no run shows `--vg` emitting a
  expressed_real_copy locus that primary-only output drops.
- **DAZ3 is a RETRACTED false positive (corrected 2026-06-02), not a recovery showcase.** Of its ~158–164
  +strand reads, the overwhelming majority are *secondary* alignments of DAZ1's reads (median NM ~6 at DAZ1
  vs ~88 at DAZ3); only ~2–6 reads genuinely prefer DAZ3. The earlier "5 isoforms at cov ~113" was phantom
  mass echoed from DAZ1. The committed anchored-prior + joint-strand EM collapses it to **cov 4.04 (2 tx,
  `low_confidence`, identifiability `none`)** — DAZ3 is real-but-very-lowly-expressed and the tool correctly
  **abstains**. (The cov-113 result reproduces only with the fix disabled:
  `RUSTLE_VG_JOINT_STRAND_EM=0 RUSTLE_VG_ANCHOR_PRIOR=0` → 5 tx, top cov 112.77.) DAZ3's value is now a
  **correctness** example — the EM refuses to fabricate a copy from a sibling's reads — guarded by
  `daz3_isoforms_max` in the oracle.
- **Open:** (1) the 93/89 rate is unreproducible from committed state (stalled verification, inputs
  uncommitted); (2) the decisive recovery test — mask a known copy out of the reference, confirm `--vg`
  recovers + assembles it — is not done; (3) no validated real-data O2 recovery exists.

## O3 — Assign ambiguous reads to the right copy (the EM)  ✅ synthetic
Fingerprint-EM: latent = copy of origin, parameters = per-copy abundance, per-read likelihood from
match/mismatch at distinguishing positions → the raw **ΔNM** anchor (T=2 decisive; ΔNM=0 non-identifiable).
Ties are **apportioned by the prior, never fabricated** — the answer to "a single SNP makes this trivial."

- Synthetic oracle: read-to-copy accuracy **100% / 100% decisive** (Obj 4).
- Derivation figures: `bench/paralog_secondary_scan/em_derivation*.png` (worked nucleotide example →
  E/M steps → convergence on DAZ-like numbers → the tie/decisive boundary).
- **Open:** real-data per-read validation on GGO is thin.

## O4 — Assemble each copy's distinct isoforms / structural variants  ✅ synthetic (capability) · 🔄 real-data
Per-copy isoform recovery — each copy gets its own source→sink paths through the family splice graph.

- Synthetic oracle: per-copy isoform **Sn/Pr 100/100** (Obj 3) — but this fixture is spatially
  separable, so the non-VG baseline ALSO scores 100/100; it certifies non-regression, not a
  StringTie-beating capability.
- **Discriminating benchmark built (2026-06-03):** `bench/tandem_attribution/o4o5_copy_benchmark.sh`
  (spec `docs/superpowers/specs/2026-06-03-o4o5-copy-resolution-benchmark.md`). Key finding:
  **gffcompare chain Sn/Pr cannot discriminate** — when copies have full-length reads the baseline
  always ties or beats VG. The honest discriminator is **copy ATTRIBUTION** (which copy each
  ambiguous read came from): on a merged-but-separable fixture VG attributes multimappers to their
  true source copy at **100% (default settings, id 0.97, 5/5 seeds)** while the baseline produces
  **no copy metric at all** (undefined, not worse). VG **abstains 100% at identical copies** (the
  DAZ limit — fabricates nothing). Measured calibration gap surfaced: at the **1-SNP boundary VG is
  overconfident** (dec_acc 0.44 at dec_frac 0.75) — fix-and-remeasure item.
- **Open:** real-data breadth (per-copy isoform recovery + copy attribution on curated GGO families);
  the synthetic truth is read-name prefixes.

## O5 — Share evidence across copies via the graph  🔄 partial · ⛔ one blocker
Cross-family borrowing (coverage + junction propagation) and the *intended* **structural-linkage channel** —
copy-specific *junctions* as distinguishing positions, adding to ΔNM to resolve reads SNPs alone cannot
(figure B's DAZ1−/DAZ3+ illustration). **Status caveat:** the junction channel (`lambda_j`) is currently
**inert/redundant** in the EM (ablation is a no-op — see `project_color_cgroup_parity`), and DAZ3 is
near-silent (~2–6 genuine reads, see O2), so this is the design concept, not a demonstrated resolution.

- Coverage/junction borrowing implemented in the VG pipeline.
- **The "intra-bundle splitter" blocker was a mirage (corrected 2026-06-01).** GOLGA6L7 is **antisense-silent**
  (0 own-strand reads; its 53+17 reads are +strand, belonging to an overlapping +strand lncRNA) — it correctly
  emits nothing; the "miss" was a gffcompare overlap artifact, not a tool bug. Checked on 3 genuine same-strand
  expressed paralog pairs (gap 15 bp–3 kb): they fail by **bundle rejection / partial assembly in separate
  bundles, never merge-collapse** — so an intra-bundle splitter has **no validated target**.
- **🔄 Real open gap:** paralog bundles (low-cov, family-secondary-stripped) get **rejected/partially
  assembled** rather than each copy emitted — a depth-aware-bundling / family-read-handling problem, not a
  splitter. This is the genuine O5 work item.

---

## Cross-cutting honest gaps
- Most scores are **synthetic**; real-data per-copy validation (O3/O4 on GGO) is the weakest link.
- O2's **93/89 rate** is classification-only and unreproducible from committed state; no validated real-data
  recovery exists, and the former DAZ3 showcase is a retracted false positive (see O2).
- O5's **GOLGA6L7 splitter** is unsolved.
- Default de-novo (non-`--vg`) headline held at **95.6 / 90.5** throughout — none of the VG work regressed it.

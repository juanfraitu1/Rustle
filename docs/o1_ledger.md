# O1 LEDGER — what has and has not worked

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the
> **superseded** 2026-07-17 catalog. The current default emits **627 families / 2,019 copies**. See
> [`NUMBERS.md`](NUMBERS.md) and [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

One page. Every route attempted against O1's definition, its verdict, and the number that decided it.
Detail lives in the linked documents; the [register](NEGATIVE_RESULTS_REGISTER.md) holds the full
639-row history. **Read this before proposing anything.**

## 1. What SHIPPED and holds

| | evidence |
|---|---|
| **`E_r` = identity ≥ 0.60 ∧ coverage ≥ 0.50 of the shorter, one record, exon-sum reps** | false-merge **2/150 = 1.33%** (HUMAN, 150 windows), **unchanged** when re-measured under the new defaults |
| **P1 seed-invariance is a THEOREM** | and it was **defended at cost** — R5 and R2 were rejected *for* violating it |
| **γ-quasi-clique at 0.20** | not decorative: at γ→0 one component holds **530 copies = 26.25% of the catalog / 56 families** |
| **λ / `cut_certified` certificates** | populated **627/627** on the current catalog; 75 are `cut_certified` |
| **transcript-orientation guard, DEFAULT since 2026-08-20** | spurious edges on the human negative panel **28 → 3**; antisense-overlap families **7.09% → 0.64%** on the real binary |
| **one path** — `--refine` rejected on the O1 catalog | a flag that changes what the catalog *is* hid a 494-vs-627 discrepancy for six weeks |

## 2. What is TRUE and must be disclosed

* **The named hole.** ~30 of 105 classified bad-family cases are definitional, one mechanism: the
  min-length coverage denominator is **scale-free**. Exposure ceiling **41/494 = 8.30%**, and 30 is a
  **floor**.
* **The discriminating clause does not order the classes.** True GFPT1×GFPT2 = **0.5353** scores below
  false ATP1A1×ATP4A = **0.5689** ⟹ no threshold on it separates them.
* **Identity never binds** — 0/728, 0/674, 245/245, 171/171. Coverage does 100% of the work.
* **Scope, not generality.** The definition covers **expressed loci with well-delineated transcript
  models**. It does not reach unexpressed copies, fragmentary models, or dispersed/diverged families
  (**0/8** at every rung).

## 3. What FAILED — by family, with the number that killed it

### 3a. The threshold / denominator / substrate space — CLOSED

| route | verdict |
|---|---|
| raise or lower `c` | **impossibility**: HERC2 splits at c_long ≥ 0.034 while the first FP dies only at 0.05 |
| coverage-of-longer on exon-sums (R2) | makes `E_r` a function of **which representative you extract**; 19/494 families fall below γ |
| coverage-of-longer on genomic spans | **no operating point** — c=0.10 costs 58.7% of TPs. Spans are mostly **intron**; introns don't align |
| genomic span as the default substrate | at `c=0.50`: **0/14 FP but 30/150 TP** — 80% recall loss. **Substrate and threshold are coupled** |
| absolute aligned bp (R3) | kills 100% of short-rep pairs at B ≥ 2000 |
| read-supported core as denominator | pipeline F1 **0.704 → 0.632** |

**All four cells of substrate × denominator fail, for *opposite* reasons: the exon-sum is too short, the genomic span is mostly intron.**

### 3b. Statistics of a different kind

| route | verdict |
|---|---|
| junction-crossing / exon structure | **12.80%** of edges genome-wide, **100× monotone bias by exon count** |
| full-length FLNC read, **edge** | FP 14/14 but TP **130/144 = 0.9028** — fires on ~90% of true pairs |
| full-length FLNC read, **node** | **inert**, 1,412/1,415 pass (T20) |
| soft clipping | **algebraic**: `clip_frac = 1 − cov`. It *is* the coverage clause |
| hard clips / supplementary linkage | 37 reads span ≥2 copies ⇒ 5 links **at one read each** |
| block internal-consistency (divergence contrast) | **AUC 0.4336**; the ancestry baseline is contaminated by repeats |
| block count / density | **0.7548 → 0.6410** per kb — half was span length |
| catalog-counted repeat promiscuity (R5) | **breaks P1** — MRPS17 scores 50 whole, 1 from a seed |

⭐ **The pattern:** the FP class **is** the short-node class, so any statistic computed **from the two
nodes** re-encodes length and collapses when normalised.

### 3c. Graph theory — every standard concept, all below chance

| concept | result |
|---|---|
| density / **Fiedler (spectral)** / min-cut | **AUC 0.36 / 0.35 / 0.35** — size dominates |
| edge betweenness | 0.683 → **0.531** residualised — "degree in disguise" |
| k-core | Simpson's paradox, p 0.0086 → **0.8125** |
| connected components | 238- and 196-member blobs at density 0.08 |
| cliques | only 81% of structural families are cliques |
| edge weighting by identity/coverage | same partition — the hard floors left **no dynamic range** |
| VG / minimizer-multiplicity | AUC 0.686, and **catalog-dependent** |
| per-family VG + flow decomposition | **Eichler's group built it and abandoned it** |

⭐ **Why below chance, not merely weak:** the false merges are **dense hubs**, not sparse anomalies —
GWFAM210's AluY hub joins 127–190 reps across 72–107 families. Community detection **optimises for
exactly that structure**, so density methods *reinforce* the error.

## 4. What WORKS but only as a CERTIFICATE

Every high-sensitivity signal is **internal** to the pair and carries the false-merge problem; every
signal that cleanly rejects false merges is **external** and has low sensitivity.

| signal | reference | FP admitted | TP reach |
|---|---|---:|---:|
| exon-sum coverage ≥ 0.50 | internal | 14/14 | **150/150** |
| genomic-span coverage ≥ 0.50 | genome | **0/14** | 30/150 |
| SD containment | SD catalog | **0/14** | 24–53/150 |
| genome-anchored repeat veto | genome | rejects 10/12 | **0/135** cost |
| transcript orientation | strand | rejects 6/14 | 4/9,032 edges |

⟹ **The architecture is one high-sensitivity RELATION plus low-sensitivity, high-precision
CERTIFICATES.** Every attempt to *replace* the relation with a certificate failed — they are not the
same kind of object.

## 5. OPEN

| | |
|---|---|
| **weight `E_r` by EXTERNAL evidence and re-partition** | the one graph-theoretic move left: external evidence never passed through `E_r`'s floors, so unlike identity/coverage it retains dynamic range. The code's own note licenses the revisit |
| re-measure the **8.30%** ceiling on the 627 catalog | needs the case census re-classified, not a recount |
| **114** SD-supported pairs γ splits *outside* the hairball | the residual over-split, down from 246 once forced cuts are removed |
| genomic span at **its own** re-fitted floor | a project (held-out fit, genome-wide measurement), not a flag flip |

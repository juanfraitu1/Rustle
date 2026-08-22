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

*Full treatment, including what each certificate licenses and how to reproduce it:*
[`o1_corroboration_layer.md`](o1_corroboration_layer.md).


Every high-sensitivity signal is **internal** to the pair and carries the false-merge problem; every
signal that cleanly rejects false merges is **external** and has low sensitivity.

| signal | reference | FP admitted | TP reach |
|---|---|---:|---:|
| exon-sum coverage ≥ 0.50 | internal | 14/14 | **150/150** |
| genomic-span coverage ≥ 0.50 | genome | **0/14** | 30/150 = 0.2000 |
| SD containment ≥ 0.90 | SD catalog | **0/14** | 24/150 = 0.1600 |
| **flank homology** ⭐ *portable* | **genome only** | **0/14** | 30/150 = 0.2000 |
| **UNION of the three** | — | **0/14** | **41/150 = 0.2733** |
| genome-anchored repeat veto | genome | rejects 10/12 | **0/135** cost |
| transcript orientation | strand | rejects 6/14 | 4/9,032 edges |

⭐ **Flank homology is the portable certificate.** Take 5 kb outside each locus, **exclude the gene
bodies**, and ask whether the extra-genic sequence aligns (≥ 300 bp at ≥ 0.90). A segmental duplication
copies a neighbourhood, so a real pair shares flanks; a mobile element inserted into two unrelated
genes does not. It reaches SD containment's specificity **using only the genome** — no SD catalog, no
species-specific third input — which is exactly the portability objection that stopped SD containment
becoming a membership condition. ⚠ A *negative* is weak (a duplication with boundaries inside the gene,
or an old one with diverged flanks, looks the same); the **positive** is the informative direction.

⚠ The three certificates **overlap substantially** — 30 + 24 + 30 = 84 if disjoint, 41 in fact — so
they measure related things (all are genomic co-duplication evidence). The union still admits **0/14**,
Wilson95 **[0.0000, 0.2153]** — bounded, not proven zero.

⟹ **The architecture is one high-sensitivity RELATION plus low-sensitivity, high-precision
CERTIFICATES.** Every attempt to *replace* the relation with a certificate failed — they are not the
same kind of object.

## 4a. ⚠ The search for more FP rules has hit diminishing returns — structurally

Four rules already reject all or nearly all of the 14 false merges; the binding constraint has always
been **genome-wide TP cost** (2.05% / 3.67% / 12.80% / 80%). And the space of *references* a new rule
could use is now largely enumerated:

| reference | status |
|---|---|
| the two nodes themselves | **exhausted** — length in disguise |
| the graph | **exhausted** — every concept scores below chance |
| the genome, via shared-anchor multiplicity | `gmult` — works, as a veto |
| the genome, via the neighbourhood | flank homology — works, portable |
| the genome, via whole-locus rank | **subsumed by `gmult`**; the tier cannot see repeats at rep scale |
| a duplication catalog | SD containment — works, not portable |
| strand | orientation — **shipped** |
| the annotation, as repeat content | softmask — works, library-based |
| the annotation, as gene content | ⚠ **circular** — entailed by the truth predicate |

**What remains genuinely untried is a different KIND of reference, not another statistic over the
genome:** another individual, another species' independently-built catalog, or experimental evidence.
Each is a project rather than a rule.

## 4b. ⭐ RECALL — what is missing so ALL members are found (2026-08-21)

Measured against **independent duplication evidence** (SEDEF `final.bed`) rather than only against
excision. Three tiers, in order of size — and the effort so far has gone to the smallest.

### Tier 1 — the copy must become a NODE *(largest, and outside `E_r` entirely)*

| | |
|---|---|
| excision destinations that are not even nodes | **54/56 = 96.4%** |
| families with a member unservable by reads | **23.08%** |

⚠⚠ **CORRECTED 2026-08-21 — tier 1 is NOT wholly an expression limit.** This section previously said
an unexpressed copy is structurally invisible and *"no edge rule, threshold or partition can recover a
locus that never became a node"*. Only part of that holds. Measured directly: segmental duplications
predict **3,928** copy locations where the catalog has **nothing**, and sampling 400 of them for read
support (`-F 2308`):

| | of 400 | |
|---|---:|---:|
| **no reads** — genuinely invisible to RNA | 170 | **0.4250** |
| 1–2 reads | 90 | 0.2250 |
| **≥ 3 reads** — clears the pipeline's own `MIN_READS`, a node *could* have been built | **140** | **0.3500** |

⟹ scaled, **~1,374 expressed, duplication-supported locations are absent from the catalog** — against
2,019 copies present, i.e. **68% as many candidate copies missing as present.** About **a third** of
tier 1 is a **node-construction gap**, not an expression limit; only ~42.5% is truly unrecoverable.

⚠ **Upper bound.** An SD partner interval need not contain a gene, and clearing `MIN_READS` is not the
node builder's only criterion — the catalog has **0/1415 single-exon nodes**, so an unspliced
read-pile correctly does not become one. Sampling was systematic (stride over a sorted key list), not
random.

### ⭐ Why this is the thesis's binding constraint

**O2 assigns a read to a copy. A copy O1 misses can never receive a read** — the possibility is not on
the table — whereas a copy O1 wrongly adds, O2 can decline. The asymmetry favours **recall**, and
essentially all effort has gone to **precision**: seven closed repair routes on the coverage clause,
which costs **7.2%**, while tier 1 leaves ~1,374 expressed candidates unlaid.

### Tiers 2 and 3 — given that both loci ARE nodes

| | pairs | share |
|---|---:|---:|
| SD-supported pairs with both loci already nodes | 1,484 | |
| ├ co-duplicated **different** genes (no alignment; correctly apart) | 346 | |
| **recall denominator** — pairs that should be one family | **1,138** | |
| **FOUND** | 810 | **0.7118** |
| **lost to the COVERAGE clause** | 82 | 0.0721 |
| **lost to γ** — the edge existed, the partition cut it | **246** | **0.2162** |

⭐ **γ costs 3.0× what the coverage clause costs.** The named definitional hole is real, but it is the
**smaller** of the two recoverable losses, and essentially all effort has gone to it.

Of the 246, **132 lie inside the 530-copy hairball** — forced cuts, where γ had to partition a blob no
definition could emit whole. The genuinely questionable over-splits are the remaining **114**.

⚠ **The 71.18% is an over-estimate of general recall.** SD-supported pairs are recent, high-identity
duplications — the easy end. And SD co-membership is not family membership, so the 1,138 denominator
may still contain co-duplicated distinct genes that happen to align.

### ⟹ What is missing, in order of size

1. **Expression / node construction** — unrecoverable within `E_r`; the substrate's limit.
2. **γ's partition** — 21.6%, of which 114 pairs are questionable rather than forced. **Untouched.**
3. **The coverage clause** — 7.2%. Seven repair routes closed; the hole is named and bounded.

## 4c. The 114 questionable γ splits, adjudicated (2026-08-21)

`bench/o1_gamma_adjudicate.py`. Each locus mapped to its best-overlapping annotated gene, then
compared. ⚠ The annotation is legitimate here — it produced neither the γ partition nor the SD calls.

| | pairs | |
|---|---:|---:|
| **UNINFORMATIVE** — different `LOC` ids | 75/114 | 0.6579 |
| UNRESOLVED — no overlapping gene | 26/114 | 0.2281 |
| **γ CORRECT** — different named genes | 11/114 | 0.0965 |
| **OVER-SPLIT** — same named gene | **2/114** | 0.0175 |

**Among the 13 adjudicable pairs: over-split 2/13 = 0.1538, Wilson95 [0.0432, 0.4229].**
The single confirmed over-split is **DHRSX**, one gene's copies separated into `GWFAM578` / `GWFAM579`.

### ⚠⚠ The methodological trap, recorded because it nearly went the other way

The first run counted **"different `LOC` ids" as evidence γ was CORRECT**, giving over-split
**2/88 = 2.27%** — a clean exoneration. It is wrong. **RefSeq assigns a distinct `LOC` id per LOCUS**,
so two copies of one *unnamed* family carry different ids **by construction**. That class cannot
distinguish "two copies of one unnamed family" from "two different unnamed genes"; it is
**uninformative**, and counting it as correct inflated the exoneration **~7×**.

(An earlier version made the opposite error — excluding `LOC` entirely — which discarded the 2
over-splits' comparability along with everything else. A `LOC` id counts for **equality** and never
for the **stem** test.)

### What this establishes, and what it does not

* ✅ **γ is not obviously broken.** Where the annotation can adjudicate, it is right 11 of 13 times.
* ❌ **It is not exonerated either.** **88.6% of the questionable set is un-adjudicable by annotation**
  — precisely because γ acts in unannotated, repeat-rich, recently-duplicated territory. The same
  blocker as the hairball: **no truth where the hard cases are.**
* ⟹ The γ question cannot be settled with the annotation. It needs truth of a different kind for
  unnamed loci — cross-species orthology, or manual curation of a sample.

## 4d. Can γ be SIMPLIFIED? (2026-08-21)

γ looks expensive: inert on **79.11%** of families, splitting only **16** components, **NP-hard** so the
code certifies *a valid witness* rather than the optimum, and documented as having **two divergent
implementations**. Two questions, both now answered.

### Can an edge veto replace it? NO — measured

γ's one real job is preventing the 530-copy hairball. If that blob were held together by repeat-bridge
edges, a pair-local veto would dissolve it and the partition step could go. Dropping hairball edges by
genome-anchored multiplicity:

| drop edges with `gmult` ≥ | edges kept | largest component |
|---:|---:|---:|
| none (γ→0) | 4,019 | 508 |
| 100 | 2,866 | 371 |
| 20 | 2,513 | 289 |
| **5** | **1,492 = 37%** | **206** |

**Dropping 63% of the edges still leaves a 206-node blob**, against γ's 56 families. The hairball's
connectivity is **distributed**, not carried by a cuttable set of bridges. ⟹ **γ's job is a PARTITION
job and no edge filter substitutes for it.**

### Is the two-implementation complexity real? NO — already gone on the O1 path

`ONE_METHOD.md`'s consistency table (2026-07-21, marked HISTORICAL) says *"Still two:
`family_definition::refine_component` (CNM) vs `family_split::gamma_quasi_clique_partition` (Louvain)…
deferred to deliverable C."* **That overstates the current state.** Traced 2026-08-21:

* `gw_family_catalog` calls **`refine_families_exon_sum`**, which routes to **Louvain**;
* the CNM `refine_component` is reached only through `family_definition::refine_families`, whose sole
  caller is **`driver.rs:386`** — the legacy `family_define` binary.

⟹ **Only one γ implementation is reachable from the O1 catalog.** The divergence is real but confined
to a legacy parity fixture, not to the shipped definition.

### What is irreducible

γ-quasi-clique partitioning is **NP-hard**, so the emitted partition is *a certified valid witness*,
not the optimum. That is inherent to the object, not an implementation choice — and it is the honest
caveat to state, rather than the two-implementations one.

| | |
|---|---|
| **weight `E_r` by EXTERNAL evidence and re-partition** | ⚠**PRICED 2026-08-20, headroom is MODERATE.** On the 530-copy hairball, γ's cuts already land preferentially on repeat-bridge edges — median `gmult` **20 for cut vs 5 for kept**, ≥50 rate **42.4% vs 22.3%**, **AUC 0.6690**. So γ captures ~⅔ of the ordering external evidence provides, **purely from density**. Residual headroom: **22.3% of KEPT edges carry gmult ≥ 50** and **57.6% of CUT edges are below it**. ⚠**Blocker: the hairball has no truth labels**, so a re-partition cannot be scored as better — it can only be shown to differ. Build only alongside truth for those 56 families |
| re-measure the **8.30%** ceiling on the 627 catalog | needs the case census re-classified, not a recount |
| **114** SD-supported pairs γ splits *outside* the hairball | the residual over-split, down from 246 once forced cuts are removed |
| genomic span at **its own** re-fitted floor | a project (held-out fit, genome-wide measurement), not a flag flip |

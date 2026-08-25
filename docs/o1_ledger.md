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
2,019 copies present.

### ⚠⚠ Both readings of that 1,374 were wrong. Corrected the same day.

**(i) The threshold was the wrong one — ~282, not ~1,374.** "≥ 3 reads in the interval" is not the
pipeline's criterion; it requires ≥ 3 reads on a **consistent intron chain**. Re-measured against the
builder's own rule (`bench/o1_why_no_node.py`, 250 sampled):

| | of 250 | |
|---|---:|---:|
| below `MIN_READS` overall | 165 | 0.6600 |
| spliced but **no chain reaching 3 reads** — correctly not a node | 67 | 0.2680 |
| **≥ 3-read spliced chain and still no node** | **18** | **0.0720** |

⟹ **~282**, a **5× reduction**. Two thirds of the "expressed" intervals do not clear `MIN_READS` at
all once measured properly.

**(ii) And those 282 cannot be blamed on the node builder.** `copies.tsv` lists **family members
only**, so "no node here" and "a node that joined no family" are **conflated**. The well-supported
examples settle it by inspection — `NC_073224.2:145475666-145481637` carries **133 primaries, a
106-read chain and 22 introns**; a locus like that unquestionably becomes a node. It is a **singleton**,
not a missing node.

⟹ ~~**The gap is in `E_r`/γ — family assignment — not in node construction.**~~
⚠⚠⚠ **RETRACTED THE SAME DAY — see §4e.** This was inference from three inspected examples, and the
measurement reverses it. `--single-copy-baseline` is also the WRONG instrument (it returns before the
`o1_homology` dispatch and builds from the read-conflict graph `E_c`; 83.3% of the shipped catalog's
own copies are structurally unreachable by it). The rep set was obtained instead from
`RUSTLE_ER_EDGE_DUMP` on the shipped run, and on the honest gated subset **node construction dominates
at 0.8067**, not `E_r`/γ.

⚠ **Upper bound.** An SD partner interval need not contain a gene, and clearing `MIN_READS` is not the
node builder's only criterion — the catalog has **0/1415 single-exon nodes**, so an unspliced
read-pile correctly does not become one. Sampling was systematic (stride over a sorted key list), not
random.

### ⭐ Why this is the thesis's binding constraint

**O2 assigns a read to a copy. A copy O1 misses can never receive a read** — the possibility is not on
the table — whereas a copy O1 wrongly adds, O2 can decline. The asymmetry favours **recall**, and
essentially all effort has gone to **precision**.

~~After both corrections the recall picture is: **~282 …** and the evidence points at **family
assignment** rather than node construction.~~ ⚠⚠⚠ **RETRACTED — §4e.** The ~282 is a scaled estimate
whose denominator is conditioned on the prediction, and both the direction and the size change once a
base rate is computed. The defensible replacement is **275 loci across 120/627 families**, decomposing
**0.8067 node construction / 0.1919 `E_r` / 0.0014 γ**.

### Tiers 2 and 3 — given that both loci ARE nodes

| | pairs | share |
|---|---:|---:|
| SD-supported pairs with both loci already nodes | 1,484 | |
| ├ co-duplicated **different** genes (no alignment; correctly apart) | 346 | |
| **recall denominator** — pairs that should be one family | **1,138** | |
| **FOUND** | 810 | **0.7118** |
| **lost to the COVERAGE clause** | 82 | 0.0721 |
| **lost to γ** — the edge existed, the partition cut it | **246** | **0.2162** |

⛔⛔ ~~**γ costs 3.0× what the coverage clause costs.**~~ **RETRACTED 2026-08-23, SIGN REVERSED — see §4m.**
The edge set these 1,138 pairs were scored against was NOT `E_r`: `bench/o1_gamma_adjudicate.py:50` uses
the **min-length** coverage form `(qe-qs)/min(ql,tl)` instead of the shipped axis-following form, and it
**never reads `f[4]`, the PAF strand field**, so it admitted minus-strand records the definition rejects
at `denovo_pipeline.rs:4473`. On the shipped instrument the ordering inverts.

 The named definitional hole is real, but it is the
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

---

## 4e. ⭐⭐⭐ The SD "node gap" AUDITED — it was the base rate (2026-08-21)

`RUSTLE_ER_EDGE_DUMP` on a re-run of the shipped invocation (2h19m, exit 0). **Determinism gate
passed**: 94,257 skeletons → 17,924 reps → 16,483 γ-blocks → 627 families → 2,019 copies, the copy set
**byte-identical** to the 2026-08-20 catalog. So `dump/ggo.nodes.tsv` (17,924 rows + `degree`) and
`dump/ggo.edges.tsv` (4,778 edges) X-ray **the shipped object**, not a lookalike. Params self-recorded
in `dump/ggo.params.tsv`: `homology_genomic_span=false`, `core_substrate=exon-sum`,
`alignment_orientation=forward-only(+)`, `min_coverage=0.50`, all `RUSTLE_*` overrides unset.

The classification that motivated §4b was then audited by 17 agents across four angles, each
adversarially verified by three independent lenses. **0 of 4 angles survived.** The arithmetic is
correct — five reimplementations plus `bedtools` reproduce **3,529 / 398 / 1 / 0** over the 3,928
intervals, and the suspected fixed-window bug never fires (max required backscan 4). **The counting was
right; the inference was not.**

### The dissolution

| comparator | never-a-node | |
|---|---:|---|
| observed (3,928 SD-absent intervals) | **0.8984** | the claimed gap |
| all non-copy SD sides (n = 73,324) | **0.9086** | matched compartment |
| SD sides with no catalog copy at either mate (n = 69,396) | **0.9092** | cleanest |

**0.8984 IS the base rate for catalog-absent SD sequence.** The set is, if anything, node-*enriched*.
The apparent deficit came from an implicit **uniform-genome** comparator, which is 13.75% SD while the
truth set is 100% SD; the two biases cancel almost exactly (uniform-null depletion −0.0755 in AUC,
anchor-status enrichment +0.0765).

### RETRACTED — do not quote

| number | why |
|---|---|
| "the gap is node construction" (a=0.8984) | it is the base rate; **no gap to attribute** |
| **"γ costs 1/3,928"** | tie-break-degenerate (0 or 1 by tie policy); p-value flips across four nulls (0.0588 / 0.0274 / 0.2338 / 0.6517) |
| **n = 3,928** as a sample size | 78.77% self-overlap, 49.97% strictly nested ⟹ **effective n = 1,556**; block-unit rate **0.8464** |
| "398 nodes" | 398 is **intervals**; the rep count is **248** |
| any depletion factor | **sign flips with the unit** — interval OR 1.32 (enriched) vs block OR 0.71 (depleted) |
| class (c) = 0/3,928 | circular — the BED was built by deleting SD intervals overlapping a copy |

### SURVIVES

- **γ is exonerated, on a clean number:** **107/17,924 = 0.60%** of reps have `degree > 0` and never
  ship (96 reachable, 2,942,231 bp = **0.083% of genome**). Third independent exoneration of γ.
- **The truth set cannot measure O1 recall.** **21.89%** of intervals have *zero* bases homologous to
  the anchor copy's exons, and **53.39%** fall below the **0.50 coverage floor `E_r` itself requires**
  (median coverage 0.4208 vs control 0.8971; positive control 40/1,180 = 3.39%, Fisher p = 4.5e-62).
  For a majority, catalog-absence is `E_r` behaving exactly as defined. **An SD partner is a duplicated
  SEGMENT, not a duplicated GENE** — the same mechanism as the false-omission work.

### ⭐⭐⭐ The real finding: the binding constraint is the PRIMARY FLAG

| | |
|---|---:|
| intervals with **zero primary reads** (`-F 2308`) yet ≥3 genuinely aligning reads | **37.14%** |
| same, restricted to exon-homologous intervals (gate defined **without** reads) | **829/1,831 = 45.28%** |
| sampled zero-primary intervals carrying a **≥3-read exact intron chain in secondary space** | **136/235 = 57.87%** |
| …of those, all-canonical (would clear `PASS1_MIN_READS`, `GATE_MIN_READS` *and* the junction filter) | **29.36%** |

Spot-checked secondaries: `de` **0.011–0.020** (98–99% identity) at **MAPQ 0**. Builder default **C5
(admit candidates from all-secondary regions) is OFF** (`src/rustle/types.rs:1044`).

**This is structurally forced.** The 3,928 are *defined* as the SD mate that **loses minimap2's
primary/secondary tie-break**. Applying `-F 2308` removes, **from one arm only**, exactly the evidence
the set's own definition selected against.

⚠ **The house rule "always `-F 2308`" is right for per-read CIGAR statistics and WRONG for asking
"is this locus transcribed" at a locus that lost a tie-break.**

⟹ ~~**O1 AND O2 ARE NOT SEPARABLE AT THE COPIES O1 MISSES.**~~
⚠⚠⚠ **RETRACTED 2026-08-22 (same day it was written) — IT FAILED ITS OWN MATCHED CONTROL. See §4i.**

 Whether a locus becomes a node depends
on an **assignment** decision made upstream by the aligner. This is a second — and far larger — scoped
exception to O1 ⊥ O2 than the co-located junction-less pair of 2026-08-13.

⚠ **45.28% is an UPPER BOUND**; primary-only 0.3465 is the LOWER bound. A ≥90%-identical duplicate of
an expressed gene attracts secondaries whether or not it is itself transcribed. Separating them
requires **PSVs, not read counts**.

⚠ C5 also inherits **`-p`/`-N` dependence**: this BAM is `-N 50 -p 0.1 --secondary=yes`; at minimap2's
default `-p 0.8` paralogues scoring under 0.8× the self-hit are discarded outright — the same
parameter-sensitivity that invalidated the haplotype-CNV copy counts. **Record `-p` and `-N` with any
copy count derived this way.**

### The defensible number

Gated on exon homology above `E_r`'s own floor **and** primary read support:
**714 intervals → 275 merged loci → 120/627 = 19.14% of shipped families.**

| | of 714 | |
|---|---:|---:|
| no node at all → **node construction** | 576 | **0.8067** |
| node with degree 0 → **`E_r` sensitivity** | 137 | **0.1919** |
| node with degree > 0 → **γ** | 1 | **0.0014** |

⟹ on the honest subset **`E_r`'s share is 19%, not 10%**, and node construction still dominates —
but on a **5.5× smaller denominator** at a **2.6× coarser unit** than §4b claimed. Anything larger runs
through the unvalidated secondary-visible tier and is bounded by the primary-flag problem, which makes
it an **O1 × O2 joint** quantity, not an O1 recall number.

---

## 4f. Design: permissive admission + PSV corroboration (2026-08-21)

The question raised: *admit loosely, tolerate some weird merges, then corroborate with the RNA.* The
architecture is right; the corroborating statistic must be **allelic, not quantitative**.

**Why read presence cannot corroborate.** Already refuted 2026-08-16: reads are **redundant with the
rep alignment** (median gain **+0.000**, > 0.10 in **0/14**), and tiling breadth discriminates
significantly in the **wrong direction** (AUC 0.1259, p = 0.0005) with a *passing* positive control.
Here it fails harder: **a phantom node's rep is built from the very reads that would corroborate it**,
so "≥ 3 spliced reads" is true *by construction* — it is the admission criterion.

**Why `E_r` cannot filter phantoms either.** A node admitted from spillover has a rep ≈ the paralog's
sequence, so it clears identity ≥ 0.60 and coverage ≥ 0.50 **trivially** — it *is* a copy at the
sequence level. `E_r` compares rep to rep and cannot see that one rep was assembled from the other's
reads. Every phantom joins the family it spilled from, attacking false-merge (1.33% [0.37, 4.73])
invisibly, because a phantom looks like an *excellent* member.

**What escapes.** The rep is a **consensus**, and a consensus averages away the discriminating signal:
spillover reads carry the **source copy's alleles at positions where the two references differ**. So:

> **A phantom retains ZERO reads after PSV assignment. A real copy retains its own.**

The machinery already ships — `psv_linkage.rs::assign_read_to_copy`, gate `n_decisive ≥ 1` (default
since 2026-06-24) — and carries a *derived* bound rather than a tuned threshold:
**τ = ln((1−p)/p)**, τ = 6.9 ⟹ **stated P(misassign) ≤ 1e-3**. On real gorilla data τ = 2 and τ = 6.9
give an identical operating point, so the principled choice is free. ⚠ **The lever is the GATE, not the
scoring** — votes ≡ LLR under flat error, identical in 16/16 configs.

**Ordering — corroboration cannot run at admission.** PSVs are defined **between copies**, so a node in
isolation has none: the merge must happen first, provisionally. Nor can it run after the partition —
removing a node would perturb a γ-block and λ certificate computed *with it present*, so the
certificate would no longer certify the object built. The only valid slot:

> **nodes → `E_r` edges → PSV corroboration → drop refuted → γ → λ**

**Three outcomes, not two.** `n_decisive = 0` — exonically identical copies with no PSVs — is a genuine
floor, **aligner-invariant** (winnowmap tested: identical; DSFAM42 stays 95% MAPQ-0 under both). Those
can be neither corroborated nor refuted, and it is exactly where phantoms are most likely. So
**corroborated / refuted / ABSTAINED**, with the abstain class flagged rather than silently merged —
the same shape as O3's detect-and-flag, and consistent with defending O2 on abstention.

### Next step

1. ✅ **Emit per-rep exon arrays** — `catalog_input::exon_blocks_str`, `nodes.tsv` gains `exons` +
   `exon_bp`. 801 tests pass; validated against 2,019 real copies (758 minus-strand) with **zero**
   convention violations. Independently reproduces the motivating figure: catalog copies are
   **0.9078 intron by bp** (audit: 0.9083 over all reps).
2. **▶ IN FLIGHT** — re-run with the new binary, then redefine (a)/(b1)/(b2) on **exons instead of
   spans**. Rep spans are **90.83% intron by bp**, so today's class (a) cannot distinguish *"the
   interval lies inside a rep's intron"* (correct rejection) from *"real transcript sequence no node
   covers"* (a real miss). That ambiguity is what invalidated the audit's expression angle.
3. **Then the only non-circular test.** Everything above is **T8** (offline re-derivation on both
   sides). Run the real binary restricted to the flagged loci in **three arms — off / admit-only /
   admit + PSV-corroborate** — scoring all three outcome classes with the abstain rate reported
   separately, never folded into either side. Two arms would not do: without admit-only, a PSV gate is
   untestable, because you cannot show a filter works without measuring what it filtered.
   Pre-registered rule: if the 137 degree-0 intervals stay degree-0 under **exon-level** coverage
   ≥ 0.50, `E_r` has a real sensitivity hole worth a theorem; if they gain edges, the residual was a
   span artefact.

⚠ Before that run: fix the duplicate rep triple `NC_073229.2:140449427-140453573` (idx 3036 and 13514
differ in strand, exon count and degree) — inert for this classification, but any `copies.tsv` ↔
`nodes.tsv` join on the triple silently picks one of two rows.

---

## 4g. ⭐⭐ The exon arrays land: 45% of "a node is here" was an INTRON (2026-08-22)

Re-ran the catalog with the exon-array build. **Determinism gate passed on every check**: 94,257 →
17,924 reps, 16,483 blocks → 627 families, copies **identical** to the shipped catalog, and the
cross-path invariant `exon_bp == exon_sum_len` holds on **17,924/17,924 reps, 0 mismatches** — so the
malformed-chain guard never fires on real data and the arrays are complete.

Reclassifying the 3,928 SD intervals on **exons** instead of spans (`bench/o1_classify_exonic.py`,
exact max-end-prefix sweep replacing the fixed 800-row window):

| | of 3,928 | |
|---|---:|---:|
| (a) never a node — ⚠ **the BASE RATE, not a gap** (§4e) | 3,529 | 0.8984 |
| **(b-EXONIC)** overlaps a spliced rep's **exon** — transcript really is there | **189** | 0.0481 |
| **(b-INTRONIC)** inside a spliced rep's **intron** only — **CORRECT rejection** | **155** | 0.0395 |
| (s) stub rep only — **INDETERMINATE** (exon array == span) | 55 | 0.0140 |

⭐ **Of the 344 intervals with a spliced rep, 45.06% are INTRONIC.** The span-based classification
counted all 344 as "a node is here"; only **189** have transcript sequence in the interval. That is the
ambiguity the exon arrays were built to remove, and it removes **nearly half** of the apparent
node-level signal. Median exonic bp in-interval where exonic: **571**.

⭐ **The stub worry, measured: small HERE.** Single-exon reps are 33.07% of the rep set, but only
**55/3,928 = 1.40%** of these intervals are stub-only. The stub problem is real at the rep level
(§ stubs memory) and does **not** materially bind this particular analysis.

γ appears once more at **1/189**, consistent with §4e's 107/17,924 = 0.60%. **Fourth exoneration.**

### ⚠ RETRACTED WITHIN THE HOUR — "E_r is the binding constraint, 188/189 = 0.9947"

Generated and killed in the same session. Spliced reps that are **not** catalog copies have degree 0 at
rate **10,460/10,518 = 0.9945** — the observed 0.9947 **is the base rate to four decimals**. The
interval set was built by excluding catalog copies, and non-catalog ⟹ degree 0 almost surely, so the
figure is entailed by construction. **The same defect as §4e's dissolved headline, one level down.**
(Sound-join check: 0/1,478 spliced catalog copies have degree 0, as required.)

⟹ **This analysis cannot rank `E_r` against node construction.** What it establishes is narrower and
solid: **45% of the "node present" class is a correct rejection**, and the residual candidate set is
**189 intervals** carrying real exonic sequence — inspectable, with coordinates, not an estimate.

---

## 4h. ⛔ SHARED-EXON / ISOFORM POOLING — CLOSED, NO-GO (2026-08-22)

The proposal: stop discarding the non-representative isoform chains, so a locus's full exonic content
participates. Designed and killed **before** spending a run (17 agents, 4 angles, 3 adversarial lenses each).
Everything below is **T8** — offline re-derivation on existing artifacts; it bounds the payload from above,
which is enough to stop, and is not enough to have proceeded.

### Dead on BOTH readings — measured on GORILLA, not transferred from human

| | measured | unit |
|---|---:|---|
| shipped `E_r` edges the rule recovers | **396/4,778 = 0.0829** | edge |
| shipped families retaining ≥1 shared-exon edge | **157/627 = 0.2504** | family |
| families with λ = 1 (a single cut destroys them) | **552/627 = 0.8804** | family |

⟹ **dead as a REPLACEMENT.** And after deleting exon matches on the *same genomic interval*:

| `MIN_COUNT` | genuine NEW node pairs |
|---:|---:|
| ≥1 | 376 (194 with fully disjoint spans) |
| **≥2** | **9** |
| ≥3 | **0** |

⟹ **dead as an ADDITION.** Pre-registered launch floor was ≥100 disjoint-exon new pairs at the target
`MIN_COUNT`; measured 9. Nine pairs cannot move a 627-family partition behind a ≥2-distinct-loci gate and
a γ splitter.

### ⚠⚠ RETRACTED: "MIN_COUNT ≥ 2 is the untested axis"

Quoted from the source docstring on 2026-08-21 and **wrong** — it was tested 2026-08-03 on HUMAN
chr1+chr15 **through the real binary**, against a matched control, on the EMITTED PARTITION:

| `MIN_COUNT` | copies | families |
|---:|---:|---:|
| 1 | 545 | 121 |
| 2 | 323 (−40.73%) | 62 (−48.76%) |
| 4 | 243 | 43 |
| 8 | 197 | 43 |

Monotone loss. The docstring has been replaced with both measurements so the lead is not inherited again.

### ⭐ TWO REAL BUGS — and they are why the rule looked alive

**1. The distinct-locus guard was missing.** The only predicate was `ra != rb`, so two representatives over
the SAME genomic span linked to each other — their exons being the same physical DNA, aligning to
themselves at ~100%.

| of the 2,795 node pairs the rule reported that `E_r` does not | | in-`E_r` control |
|---|---:|---:|
| overlapping rep spans | **2,601 = 0.9306** | 50/396 = 0.1263 |
| every matched exon on the same physical interval | 2,419 = 0.8655 | 49/396 = 0.1237 |
| overlapping **and opposite-strand** (a locus matching its own antisense model) | 2,348 = 0.8401 | 3/396 = 0.0076 |

The control **inverts**, so this is a property of what the rule ADDS. **FIXED** (guard mirroring
`distinct_locus_reps_grouped`, with a positive control in the test). 802 tests pass.

**2. `MIN_COUNT` is not counting independent exons.** Under `ISOFORMS=1`, `support` accumulates exon
INDEX pairs over an UN-DEDUPLICATED pool, so one genomic exon pair contributes `n_iso_a × n_iso_b`. It
measures **isoform multiplicity**. ⟹ every `MIN_COUNT` figure ever quoted as evidence about *independent
exons* is void (the partition-cost table above is unaffected, and is what actually kills the axis).
Dedup fix specified, **not landed** — the path is a recorded dead end.

### ⭐⭐ CORRECTION — a GORILLA-NATIVE RefSeq GFF exists

`/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_genomic.gff` — NCBI RefSeq `GCF_029281585.2-RS_2024_02`,
genome-build `NHGRI_mGorGor1-v2.0_pri`, taxon 9595, **41,193 gene+pseudogene records**, contigs matching the
substrate. The standing claim *"there is no gorilla RefSeq GFF on this machine"* was **false**, over-generalised
from the (true) fact that the HUMAN `HSA_genomic.gff` drops 29.1% of loci. It had been constraining truth-set
designs. Call results **"gorilla-native loci, human-curated relation"** — never "gorilla ground truth".

### ⟹ Where the finding actually routes

**93.06% self-overlap means this belongs to the NODE BUILDER, not the edge rule** — the same destination
`o1_error_case_census.md` reached for the coverage hole (47/105 cases node-construction). The real defect the
investigation surfaced is upstream and already documented: rep selection keeps ONE chain, **46% of
representatives covering a known family member are single-exon stubs**, and **53% of those loci have a
gate-passing ≥3-read spliced chain that was discarded** (NOTCH2NLA / SRGAP2C / LIMS1 at 92 / 65 / 124 reads).

⚠ Note also, for any future additive tier: **γ-quasi-clique is NON-MONOTONE under edge addition** — 4 shipped
families / 63 copies drop below γ=0.20 from a single +1-node/+1-edge attachment, and 626/627 families have at
least one such partner. A union tier can LOSE copies; never pre-register a one-sided readout.

### Sharp edge found while testing

`exon_seqs_of` derives exon LENGTHS from genomic coordinates and slices `seq` by them, so a rep whose span
and sequence length disagree contributes **ZERO exons, silently**. This is why the guard test's positive
control failed on first run. ⚠ It also means `pooled_isoform_exons_recover_what_a_stub_representative_lacks`
passes its rep-only assertion because its `stub` (span 1,200, seq 600) emits no exons — **not** for the
"carries only filler" reason its comment states. That test is weaker than it reads.

---

## 4i. ⛔⛔ §4e's PRIMARY-FLAG FINDING IS REFUTED BY ITS OWN CONTROL (2026-08-22)

§4e concluded that the binding constraint on O1 recall is minimap2's primary flag, and therefore that
**O1 and O2 are not separable at the copies O1 misses**. That was committed, pushed, and written to
memory. **It is wrong.** The statistic was never compared against a matched control; when one is
computed the effect vanishes.

**Candidate** = SD mate of a catalog copy, exon-homologous, no rep. **Control** = SD side with NO catalog
copy on *either* mate, no rep, size- and compartment-matched. Unit = merged genomic locus.

| gate stage | candidate (n=882) | control (n=881) | |
|---|---:|---:|---|
| zero-primary | 497 = **0.5635** | 679 = **0.7707** | ⚠ **control HIGHER**, p ≈ 0 |
| + ≥3 good secondaries | 192 = 0.2177 | 211 = **0.2395** | control higher, p = 0.2817 |
| + **≥3-read EXACT intron chain** | **109 = 0.1236** | **104 = 0.1180** | **Fisher p = 0.7701** |

⟹ *"reads land here as secondaries carrying a shared intron chain"* is **the base rate of duplicated
genic sequence** — not evidence that a copy is transcribed there. **O1 ⊥ O2 STANDS. No second scoped
exception is warranted.** Reinforcing it: **21,770/23,807 = 0.9144** of the secondary reads at those loci
carry the **anchor's** alleles, i.e. the class is overwhelmingly spillover.

### ⚠⚠ CORRECTION 2026-08-23 — §4i's OWN CONTROL WAS DEFECTIVE. STATUS: OPEN, NOT REFUTED.

The gate script filtered secondaries at `de <= 0.05`, which discards **two thirds of all overlapping
secondaries** (97,980/290,942 = 0.3368 kept) and discards them NON-RANDOMLY: it keeps only alignments
>= 95% identical, which is backwards for detecting a DIVERGED second copy. Sweeping the threshold:

| `de <=` | CAND | CTRL | Fisher p | cand secondaries kept |
|---:|---:|---:|---:|---:|
| **0.05** (as run) | 109/882 = 0.1236 | 104/881 = 0.1180 | **0.7701** | 0.3368 |
| 0.10 | 247/882 = 0.2800 | 164/881 = 0.1862 | **< 0.0001** | 0.8090 |
| **0.20+ (unfiltered)** | **270/882 = 0.3061** | **181/881 = 0.2054** | **< 0.0001** | 1.0000 |

⟹ **the tie was a filter artefact**, so "p = 0.7701" and "it IS the base rate" are withdrawn as the
BASIS of the retraction.

⚠⚠ **But the finding is NOT restored, because the control is not matched on the thing that generates
reads.** A CANDIDATE is the SD mate of a **catalog copy** — an expressed locus. A CONTROL has no catalog
copy at EITHER mate. So candidates have a nearby expressed paralogue to spill from and controls do not,
which predicts more secondaries at candidates **with no transcription at the candidate required** —
consistent with 21,770/23,807 = 0.9144 of those secondaries carrying the ANCHOR's alleles.

⟹ the gate alone cannot decide it: it measures read VOLUME, and volume is confounded with the label.

### ⭐⭐⭐ CLOSED 2026-08-23 — the volume-insensitive test says SPILLOVER. §4e IS REFUTED.

`bench/o1_competitor_resolved.py`, pre-registered before the result. For every SECONDARY read at a
gate-passing locus L, compare minimap2's own `de` at L against `de` at **that read's own primary locus P**
— the alternative it actually contested. This asks whether the reads that arrive PREFER L, and is
therefore **insensitive to how many arrive**, which is what the gate could not separate. Nothing is
re-derived: both `de` values are minimap2's, so this is not the offline-PSV route that has failed 4x.

| | CAND (270 loci) | CTRL (181 loci) |
|---|---:|---:|
| median per-locus fraction of reads fitting **L** better | **0.0000** | **0.0000** |
| pooled READ unit | 233/108,688 = **0.0021** | 151/37,559 = **0.0040** |

**Mann-Whitney over per-locus fractions: U = 23,271, p = 0.2685 — TIE.**

Three comparators, all failed:
1. **the matched control arm** — tie (p = 0.2685), and the candidate rate is if anything the LOWER of
   the two (0.0021 vs 0.0040), i.e. pointing AWAY from the hypothesis;
2. **the genome-wide rate** at which a secondary wins by `de` at all — **0.0196**. Both arms sit ~5-9x
   BELOW it, which is what high-identity SD mates should look like if the reads are spillover;
3. **a constant predictor** — "never prefer L" scores 0.0021 wrong on candidates. The median locus has
   ZERO reads preferring it.

⟹ **The secondaries at these loci fit the locus they were assigned to, essentially always. They are
spillover** — corroborating the 21,770/23,807 = 0.9144 anchor-allele figure by an independent statistic.

⭐ **VERDICT (pre-registered): C5 stays OFF permanently. No second O1 ⊥ O2 exception. §4e's primary-flag
finding is REFUTED** — not by the defective p = 0.7701 gate control, but by a test the volume confound
cannot reach. **O1 ⊥ O2 STANDS.**

⚠ Note the sequence honestly: the finding was asserted (08-22), retracted on a broken control (08-22),
the retraction's basis was withdrawn (08-23), and the claim was then refuted on a sound one (08-23). The
conclusion is the same as the first retraction; **the reasoning that got there the first time was not
sound, and only the pre-registered test settled it.**

⚠ **This is the FOURTH headline in one week to die to the same defect — a rate quoted without its
comparator.** The other three: §4e's own 0.8984 "never a node" (= the base rate 0.9086); the 0.9947
degree-0 rate (= a 0.9945 base rate); the shared-exon 93.06% self-overlap (inverts on its control).
**Treat "I have a rate but not its comparator" as a stop condition, not a to-do.**

**Blast radius that would have been at stake**, recorded so the cost is known when this is proposed
again: 109 merged loci, **47/627 = 0.0750 of families**, 201/2,019 copies; under worst-case pendant
attachment **4/47 receiving families drop below γ = 0.20** and **13/75 = 0.1733 of all λ ≥ 2 certificates
are lost**. **C5 stays OFF** (`src/rustle/types.rs:1044`).

⚠ Register the standing recommendation *"keep a site if ≥3 reads share an intron chain by ANY
alignment"* as **REFUTED — zero specificity**.

### ⚠⚠ AND §4f's premise is gone with it

The permissive-admission + PSV-corroboration design was motivated by that finding. Its ordering argument
and its phantom analysis remain correct as *reasoning*, but there is now **no measured reason to admit
anything**. Do not build it.

## 4j. O1 → O2 IS A CONDITIONING, NOT A CORRELATION (2026-08-22)

Tested whether an O1 graph certificate bounds O2's assignability — the theorem-shaped claim. **It is
false by construction, and I checked rather than assumed.**

- `grep -E "lambda|stoer|min_cut|density|gamma|quasi_clique|cut_certified" src/rustle/vg_family/copy_assign.rs`
  → **0 hits.** The ONLY quantity O2 consumes from O1 is the family **cardinality**, at
  `copy_assign.rs:453`: `let thr = p.alpha / (n.saturating_sub(1).max(1) as f64);`
- It could not matter much if it did: on 104,147 read-records **90.95% sit outside any movable zone**
  (min_p exactly 0 in 39.82%, < 1e-12 in 26.73%, ≥ 1.0 in 29.40%), and collapsing **every** family to
  n = 2 — the maximal relaxation of that denominator — moves **14/104,147 = 0.013%** of verdicts.
- Every retrofit failed: density / Fiedler / min-cut **AUC 0.36 / 0.35 / 0.35** (§3c); `de`-weighted
  capacity **equals `min_de` on 515/627 = 0.8214** and is ≥ it by construction on the rest; `E_r` degree
  is **at chance after depth+length matching** (0.4607 / 0.4899 / 0.4839 vs 0.5000).

### ⚠ ENTAILED — do not sell as findings

| claim | why it is entailed |
|---|---|
| **λ = 1 on 552/627 = 0.8804** | **440/627 = 0.7018 are 2-copy, and λ = 1 there is a THEOREM.** The informative residue is **112/187 = 0.599 of families with ≥3 copies**. Quote *that*. |
| "λ certifies a graph blind to 42.57% of within-family pairs" | = 1 − density; γ = 0.20 *permits* 80% of pairs edgeless by definition. The parameter restated. |
| "family order beats α by 20×" | matched perturbation of n and of α gives **0/104,147** disagreements. |
| "n_decisive = 0 ⟹ abstain is an identifiability floor" | 100% of n_decisive = 0 rows are tied — it *is* the `if`. |
| δ-bound "copies below divergence δ are unassignable" | pigeonhole on the definition of divergence; the headline moves 3.07% → 8.92% → 13.08% across δ = 0.0005 / 0.005 / 0.01. |

### ⟹ The defensible statement

**O1 → O2 is a CONDITIONING, not a correlation.** Given the copy set, assignment decomposes and per-read
argmax is optimal; the combinatorial difficulty lives entirely in **choosing the set** — which is O1,
with O3 discharging the precondition that the set is complete. That is also why O2 is defended on
**abstention**, its only non-circular validation: excision-held-out **TPR 0.5066 / FPR 0.0280,
AUC 0.7995** against **MAPQ at AUC 0.4944** (median 60 vs 60) on the same reads.

⚠ Do NOT present "O1 and O2 are views of one graph object" — already REFUTED (register §9.1:994),
reframed as **two graphs, one arrow**.

## 4k. O3 AS A RESIDUAL OF O1+O2 — DEGENERATES ON ITS OWN POSITIVE CONTROL

- The **abstention** residual is undefined at **162/162** excision positives: the panel is TWO-COPY, so
  masking leaves |C| = 1 — no assignment, hence no abstention.
- The **PSV** residual at |C| = 1 **IS the shipped S2 within-pile mixture detector** — a rename.
- The **depth** residual was already refuted: 16/104 destinations had zero baseline reads; the 1.75×
  needs a "before" that a real case does not have.

**O3's honest status is unchanged.** Bar to beat: **TPR 0.2703 / FPR 0.0200** on the excision panel.
Sensitivity is set by **divergence, not abundance** (0.4500 at ≥0.01 diverged vs 0.0588 below) and
**45.78% of real positives sit below that line** — a property of the problem, not the detector.
Two levers, in order: an **n ≥ 3 excision panel** (⚠ feasibility: 14 families at ≥20 reads, 56 at ≥5,
157 at ≥3 — and only **6/43 = 0.1395** of scenarios put the deleted copy below de 0.01 against
**45.78%** of real positives, so the panel would sample the WRONG divergence regime; measure that
distribution offline first, it is free), and a **larger compartment**.

---

## 4l. ⭐⭐⭐ CROSS-SUBSTRATE REPLICATION — the RELATION transports, the PARTITION less so (2026-08-23)

**The first test of whether O1's object is GENOMIC or SAMPLE-CONDITIONAL**, and the first result in this
line that is **NOT T8** — it ran through the real binary twice.

Every thesis number rests on a catalog built from `GGO_ds.bam` = **OR6737 TESTIS**, a *different animal*
from the assembly (KB3781). A depth-matched KB3781 **FIBROBLAST** catalog was built with the same binary,
reference and parameters, so exactly two things differ: **animal and tissue**.

⚠ `fibro_ds.bam` was NOT depth-matched to the shipped substrate despite the note saying so — 3,962,281
primaries against 1,628,629, i.e. **2.43x**. Subsampled by READ-NAME hash (`-s 42.4110`, so every
alignment of a kept read survives together and the primary/secondary structure `E_r` depends on is
preserved) to **1,627,629 primaries — 0.06% off target.**

| | testis (OR6737) | fibroblast (KB3781) |
|---|---:|---:|
| primaries | 1,628,629 | 1,627,629 |
| skeletons | 94,257 | **55,366** |
| reps | 17,924 | **13,196** |
| families | 627 | **356** |

⚠ **THE COUNTS ARE NOT THE READOUT.** Testis transcribes far more of the genome, so fewer fibroblast loci
is expected biology. Judging the definition on per-sample yields would repeat the error made 3x here
(judging what a NODE is on node-level counts). The question is agreement **on the loci both see**.

### Agreement on shared loci — reciprocal-50% matched, n = 10,143 reps

(0.5659 of testis, 0.7686 of fibroblast. ⚠ The matcher's lookback window is **inert**: 60 / 400 / 2000
give byte-identical results, so it is not truncating — the defect that nearly invalidated the SD work.)

| | measured | comparator |
|---|---:|---|
| **`E_r` edges recovered** — testis edges on shared loci also found by fibroblast | **1,171/1,345 = 0.8706** | edge density ~2.7e-5, so chance ≈ 0 |
| edge **Jaccard** | 1,171/1,583 = **0.7397** | as above |
| partition **ARI** on 724 shared catalog copies | **0.7064** | size-preserving label permutation: mean **-0.0000**, 95th pct 0.0038, **p ≤ 0.005** (200 perms) |
| **co-membership pair Jaccard** (the harsher, uncorrected framing) | 1,192/2,174 = **0.5483** | — |

⭐⭐ **THE RELATION TRANSPORTS ACROSS ANIMAL *AND* TISSUE: 87.06% of `E_r` edges are recovered.** That is
the strongest substrate-independence evidence O1 has, and it is consistent with the standing conclusion
*"the RELATION transcends; the LEVEL can be made to; the NODE SET cannot."*

### ⚠ The partition transports LESS well, and the instability is DIRECTIONAL

Of 261,726 copy pairs on shared loci: together in both **1,192**, apart in both 259,552,
**testis-together/fibroblast-apart (fibroblast SPLITS) 290 = 0.2953 of disagreements**, and
**testis-apart/fibroblast-together (fibroblast MERGES) 692 = 0.7047**.

**Fibroblast merges 2.4x more often than it splits** — and its graph on the shared loci is *denser*
(1,409 edges vs 1,345) despite having 26% fewer reps overall, which points the mechanism at γ: a denser
graph yields larger quasi-clique blocks. ⚠ Report ARI 0.7064 **and** the pair Jaccard 0.5483; ARI is
chance-corrected but 259,552/261,726 = 99.2% of pairs are trivially apart in both.

### SELECTION CHECK — the shared stratum is NOT the easier one

Testis reps that MATCHED have **mean degree 0.429** and 10.46% with degree > 0; the UNMATCHED have
**0.669** and 13.69%. So the shared loci are, if anything, *less* connected — the agreement above is not
bought by restricting to a well-connected core.

### What this does and does not establish

- ⭐ It **removes the cross-individual objection** as a blanket criticism: the relation survives a change
  of animal *and* tissue at 0.8706.
- ⚠ **Animal and tissue are CONFOUNDED** by this design; it cannot say which drives the residual.
- ⚠ One fibroblast run; no determinism gate exists for it (no prior fibroblast catalog to compare).
- ⚠ The 627 vs 356 family difference is **expected biology, untested** — do not quote it as instability.

---

## 4m. ⛔⛔⛔ γ IS NOT THE PROBLEM — THE ORIENTATION GUARD IS (2026-08-23)

Commissioned to design the deferred hierarchical emission against the "114 questionable γ over-splits".
**The target set is empty, and the premise that produced it was measuring a different object.**

### The finding: 0/114

**0/114 of the recovery target is recoverable at ANY level, because 0/114 share an `E_r` connected
component** — and coarsening can only merge *within* a component. Wilson 95% **[0.0000, 0.0326]**.
⭐ The comparator shows the lookup is not empty by construction: **129–132/132** of the "forced"
inside-hairball pairs *are* same-component.

**Mechanism: 113–114/114 have a qualifying MINUS-STRAND record and no qualifying forward record.**
`ggo.rule.tsv` records `alignment_orientation = forward-only (+)`; the guard is
`denovo_pipeline.rs:4473-4474`. **γ never saw these pairs.** They were removed one stage earlier, by the
orientation guard shipped as a default on 2026-08-20.

### ⛔ RETRACTED — and the sign reverses

Re-scored on the **shipped** `-k 11 -w 5` PAF with the **shipped** coverage form (the only instrument the
definition actually uses), copy-pair unit, denominator 1,135:

| loss channel | pairs | share |
|---|---:|---:|
| **orientation guard** | **167** | **0.1471** |
| **coverage clause** | **147** | **0.1295** |
| identity floor | 0 | 0.0000 |
| **γ — the partition** | **11** | **0.0097** |

⟹ **the EDGE RULE costs ~28× what γ costs**, not the reverse. Guard:coverage is **1.14×, near parity** —
not the 2.80× one arm claimed. ⚠ Three re-derivations gave three answers because two used the offline
`gw_ava_*.paf` (different presets) instead of the shipped PAF. **Quote only the row above.**

| also retracted | correction |
|---|---|
| "the honest recovery target is the remaining 114" | **the set is empty as a γ target** |
| "530 copies = 26.25% of the catalog / 56 families" at γ→0 | **466 copies / 38 families** (+14/4 in a second component) = **480/2,019 = 0.2377** across **42/627 = 0.0670** of families |
| §4l's "the PARTITION does not transport (ARI 0.7064)" | **blob-driven.** Drop the single largest component (~172/724 = 23.8% of shared copies) and flat **ARI = 0.9707**, coarse 0.9578 — **the flat partition transports BETTER.** The 87.06% edge-transport figure is unaffected. |

### ⛔ NO-GO on a hierarchy, on its own arithmetic

The level that recovers the genuine γ-cut pairs **is** the level that re-forms the blob — they are the
same object, and there is no threshold in between. Matched-permissiveness pricing of the coarse level:
new co-membership pairs are duplication-supported at **167/104,125 = 0.0016** against **810/6,469 =
0.1252** for the relation the catalog already asserts and **507/1,926,577 = 0.0003** for a random
cross-component pair. **78× less supported than production**; to surface 167 supported pairs it asserts
**103,958 unsupported** ones.

A *recent* (refinement) level adds **0** new assertions by theorem, but is informative on only
**69/627 = 11.00%** and **splits 36/523 = 6.88% of SEDEF-certified co-contained pairs at id ≥ 0.90** —
measurably wrong on the only stratum with independent truth. Its "SEGMENTAL dissolved 0/136"
corroboration is **circular** (selected at SD id ≥ 0.90 while DISSOLVED *is* max within-family identity
< 0.90).

### ✅ MODIFIED-GO — a disclosure, and it must NOT be called a hierarchy

`copies.tsv` / `copies.fa` **byte-identical**; **627 families / 2,019 copies** remains the headline.
Three columns appended to `families.tsv` after `cut_certified` (the same additive contract as
`lambda`): **`er_component_id`** (min `node_key` of the family's `E_r` component; 627/627 families lie in
exactly one), **`er_component_n_families`** (42/627 = 0.0670 carry > 1), **`gamma_forced`** (true iff γ
actually split this component; 42/627, 2 components genome-wide).

⚠ **Membership is still answered by `family_id`. There is no query.** These columns are a *disclosure* —
"this family is one γ-block of a 466-copy component" — and must not be consumed by any predicate, gate or
filter. Name it `er_component_id`, **never `superfamily_id`**: nobody may write "a 466-copy KZNF family".

### ⟹ Where the effort should go

**The orientation guard is now the largest single measured loss channel (0.1471).** It is not wrong — it
cut spurious edges **28 → 3** on a HUMAN negative panel and antisense families **7.09% → 0.64%** — but
its recall cost had never been measured, and it is at parity with the coverage clause. **That trade is
now quantified on both sides and is the open question.** γ, at 0.0097, is not worth further work.

---

## 4n. ⭐⭐⭐ THE ORIENTATION GUARD IS FILTERING ON AN UNMEASURED FIELD (2026-08-24)

§4m established the guard as O1's largest loss channel (167/1,135 = **0.1471** of SD-supported pairs,
against coverage 0.1295 and γ 0.0097). This is the mechanism, and it is a **defect, not a trade-off**.

### A third of the rep set has a PLACEHOLDER strand, not a measured one

`denovo_assemble.rs:1010`: `let strand = strand.unwrap_or('+');` — a rep's strand comes from the gate's
canonical-junction classification, so a **single-exon rep has no junctions and no determinable strand**,
and the code stamps `'+'`.

| rep class | n | `'+'` | `'-'` |
|---|---:|---:|---:|
| **single-exon** | 5,928 | **5,928 = 1.0000** | **0 = 0.0000** |
| spliced (≥2 exons) | 11,996 | 5,839 = 0.4867 | 6,157 = 0.5133 |

Not one single-exon rep is `'-'`. Spliced reps split ~50/50, which is what biology looks like.
⚠ **33.07% of the rep set carries an unmeasured field.** If single-exon reps resemble spliced ones,
roughly **3,000 of the 5,928 currently carry the WRONG strand** *(estimated from the spliced ratio)*.

### Why that costs edges

A rep's `seq` is stored **in its `strand` orientation**, so a rep wrongly marked `'+'` is stored
**reverse-complemented**; its alignment to a correctly-oriented paralogue comes out MINUS-strand, and the
guard drops it at `denovo_pipeline.rs:4473`. Measured on the shipped `-k11 -w5` PAF under the shipped rule
(`bench/o1_guard_cost.py`), rep-pair unit:

| | measured | |
|---|---:|---|
| guard-blocked pairs (a qualifying `'-'` record, no qualifying `'+'`) | **4,009** | |
| …involving ≥1 **single-exon** rep | **3,951 = 0.9855** | ⭐ |
| **…both spliced — genuine antisense candidates** | **58 = 0.0145** | genome-wide |
| **COMPARATOR** — kept (`'+'`) pairs involving ≥1 single-exon rep | 1,884/4,778 = **0.3943** | **2.50× enrichment** |

⟹ **the guard is overwhelmingly discarding artefacts of its own placeholder, not antisense biology.**

### ⭐ The fix is to MEASURE the strand, not to relax the guard

Relaxing it (e.g. to `n_exon ≥ 2` both sides) would admit 3,951 rep pairs wholesale and forfeit the
precision the guard was shipped for (spurious edges **28 → 3**, antisense families **7.09% → 0.64%**,
HUMAN panel). Instead, determine the strand the reads already carry — then those reps are stored in the
right orientation, align `'+'`, and **pass the guard without it being weakened.**

**Validated on a LABELLED set** (`bench/o1_strand_recovery.py`): spliced reps have junction-determined
strand, i.e. ground truth. On 400 sampled at random from 11,996:

| signal | agreement with junction-determined strand |
|---|---:|
| `ts:A:` majority (BAM is `-ax splice:hq -uf`) | **386/400 = 0.9650** |
| FLAG 0x10 majority | **386/400 = 0.9650** |
| **COMPARATOR — constant `'+'` predictor** | **0.4867** |

⟹ read orientation recovers the strand at **0.9650** vs a 0.4867 floor. Applied to single-exon reps this
replaces an ~50%-wrong placeholder with a ~3.5%-wrong measurement — roughly a **14× reduction in
mis-stranded reps** *(estimated)*.

### ⚠ How this MUST be judged

This changes **what a node IS** (the orientation its sequence is stored in). Judging such a change on
node- or edge-level counts has **failed 3× end-to-end**. Acceptance must be on:
1. the emitted **PARTITION** (families/copies), not edge counts;
2. the **HUMAN 150-window negative panel** — false-merge must not regress from 2/150 = 1.33%;
3. the **antisense-family rate** — must not regress from 0.64% back toward 7.09%;
4. the 58 genuine antisense pairs must **still be rejected**.

⚠ Un-run. Requires a source change plus a full catalog rebuild (~2h20m), and it is a **definition-level
change**, so the negative panel is mandatory, not optional.

---

## 4o. ⭐⭐⭐ READ-STRAND MEASUREMENT — RUN, AND EVERY PRE-REGISTERED CRITERION PASSES (2026-08-24)

§4n found that `strand.unwrap_or('+')` gives a third of the rep set an UNMEASURED strand. Both arms are
now run through the real binary. **This is NOT T8.**

### The OFF gate — the opt-in is genuinely inert

`copies.tsv` md5 **35eb5c51c141cda25fc7c1866d310f1a**, identical to the shipped catalog; 94,257 skeletons
→ 17,924 reps → 627 families; single-exon reps still 5,928 `'+'` / 0 `'-'`. `params.tsv` differs by
**exactly one row** (`env.RUSTLE_READ_STRAND <unset>`) plus the output path. 2:11:27, 24.4 GB peak.

### The ON arm — the placeholder is gone and the distribution looks like biology

| | OFF | ON |
|---|---:|---:|
| single-exon reps `'+'` / `'-'` | **5,928 / 0** | **2,060 / 1,926** |
| skeletons → reps | 94,257 → **17,924** | 94,257 → **15,626** |
| families / copies | 627 / 2,019 | **630 / 1,958** |
| wall / peak RSS | 2:11:27 / 24.4 GB | 2:44:23 / 24.4 GB |

The measured split **0.5168 / 0.4832** closely tracks the spliced (junction-determined) distribution
**0.4867 / 0.5133** — a measured field that looks like biology, not noise.

### The TWO-SIDED ledger (`bench/o1_guard_delta.py`) — the reason a one-sided count would mislead

| | measured |
|---|---:|
| guard-blocked pairs OFF → ON | **4,009 → 551 (−86%)** |
| **GAINED** (934 of them were guard-blocked in OFF) | **972** |
| **LOST** (only 214 newly guard-blocked; the rest are pairs whose reps MERGED) | **864** |
| **NET** | **+108** |

⚠ **A one-sided ledger would have reported "934 recovered".** The net is +108. Report both arms.

### ⭐ ACCEPTANCE — all pre-registered before any result, all PASS

| criterion | OFF | ON | |
|---|---|---|---|
| the 58 frozen genuine-antisense pairs admitted (keyed by `node_key`) | — | **0/58** | **PASS** |
| families with an overlapping OPPOSITE-STRAND pair | 3/627 = 0.0048 | **3/630 = 0.0048** | **PASS** (pre-guard 0.0709) |
| **HUMAN 150-window negative panel, false-merge** | **2/150 = 0.0133** | **2/150 = 0.0133** | **PASS — unmoved** |
| partition shape | 627 fam, size-2 440, max 42 | 630 fam, size-2 441, max 40 | stable |

⚠ **Negative-panel E_r edges rose 3 → 5.** Both new classes are ALREADY-NAMED pathologies, not new ones,
and neither reaches `copies.tsv` (which is why family-level false-merge is unmoved):
- **W038** `chr14:64934142-64989490(-) × chr14:64928524-64936834(-)`, identity **exactly 1.000000**,
  spans OVERLAP — the **self-identity** pathology, one locus emitted as two (same class as the known W063).
- **W106 ×2**, targeting a **206 bp** rep at `chr5:141029012-141029218`, identity 0.834/0.831,
  coverage 0.966–1.000 — the **scale-free coverage denominator**: a 206 bp rep is ≥0.50 covered by almost
  anything. This is O1's named definitional hole, surfaced (not created) by the strand fix.

### ⟹ WHAT THIS ACTUALLY IS: a NODE-CONSTRUCTION fix, not an edge-recovery win

The partition barely moves (627 → 630) while the **rep set moves a lot (17,924 → 15,626, −12.8%)**:
co-located reps that a placeholder strand had artificially kept apart now collapse correctly
(`family_detect.rs:670` gates span-aware collapse on `a.strand != b.strand`). So the dominant effect is
**fixing over-fragmentation of loci** — which is exactly where the ledger says O1's largest unaddressed
problem lives — and NOT the edge recovery the guard analysis implied.

⚠ **STILL OPEN before flipping the default**: **61 copies left the catalog** (2,019 → 1,958). At 0.9650
strand accuracy some merges are wrong, and the human panel cannot see gorilla-side merges. Understand
those 61 before making this the default. **Ships opt-in (`RUSTLE_READ_STRAND`); the default is unchanged.**

---

## 4p. ⚖️ WHERE THE 61 COPIES WENT — ROUGHLY HALF CORRECTION, HALF LOSS. **DO NOT FLIP THE DEFAULT.**

§4o passed every pre-registered criterion but lost 61 copies net. This resolves them
(`bench/o1_strand_copy_delta.py`).

### The net hides the churn

**224 lost, 163 gained** — not 61 of anything. Classifying the 224 (⚠ the "is this locus still
represented" question is ASYMMETRIC: a merged rep is LARGER than its parts, so a reciprocal-overlap rule
would manufacture losses):

| class | n | | |
|---|---:|---:|---|
| **B DEMOTED** — an ON rep covers it, but it is in no family | 112 | 0.5000 | the concerning class |
| **A ABSORBED** — an ON CATALOG COPY covers it (merged, still in a family) | 61 | 0.2723 | benign; the change working |
| C "VANISHED" — no ON rep covers ≥50% of it | 51 | 0.2277 | ⚠ see below |

⚠ **My own ≥50% rule manufactured most of class C.** Re-examined under ANY overlap, **38/51 = 0.7451
still have an ON rep** (the locus survives under a smaller/different representative, because
`pick_locus_rep` re-picks after a merge) and 19/51 overlap an ON catalog copy. **Only 13/51 are truly
absent from the rep set.**

**The strand fix is demonstrably the cause**: single-exon copies are **144/224 = 0.6429** of the losses
against a base rate of **541/2,019 = 0.2680** — **2.40× enriched**.

### ⚖️ The verdict: the dissolved families split almost evenly

112 demotions touch 58 OFF families; **35 dissolved entirely** (median OFF size 2). Reading each
dissolved family's members at their ON rep strand:

| | n | |
|---|---:|---|
| members now **OPPOSITE-strand** ⟹ the family was an **ANTISENSE ARTEFACT** | **16/35 = 0.4571** | **dissolving is CORRECT** |
| members still **SAME-strand** ⟹ the edge was lost for another reason | **18/35 = 0.5143** | **a genuine cost** |
| indeterminate | 1/35 | |

Examples correctly dissolved: **GWFAM46** (`116689952-116693606` OFF`+`→ON`-` vs `130657657-130685127`
ON`+`), **GWFAM65**, **GWFAM109**.

### ⟹ DECISION: **RUSTLE_READ_STRAND STAYS OPT-IN. The default is UNCHANGED.**

The change fixes a real defect — a third of the rep set carried an unmeasured strand — and it corrects
**16** genuine antisense artefacts while cutting guard-blocked pairs 86%. But it costs **18** real
families, and at the family level the catalog is roughly break-even (627 → 630) with large churn in both
directions. **That is not a basis for changing the definition's default.**

### ⭐ The actionable refinement: ABSTAIN instead of guessing

The 18 lost families are the price of a **0.9650** vote propagating through locus collapse
(`family_detect.rs:670` collapses on strand, so one wrong call merges two distinct loci or splits one).
The fix is not a better signal but a **decision rule**: when the read vote is not decisive, **keep the
`'+'` placeholder rather than flipping**. Abstention converts a wrong call into the status quo ante,
which is exactly the shape O2 is already defended on (assign-or-abstain, `n_decisive >= 1`) — and it
makes the two objectives consistent in their treatment of insufficient evidence.

⚠ Un-measured: what majority margin buys what accuracy, and how many of the 16 correct dissolutions
survive a stricter rule. That trade is the next measurement, and it is cheap — the vote tallies are
already computable from the BAM without a catalog rebuild.

---

## 4q. ⛔ ABSTENTION DOES NOT RESCUE IT — THE 18 LOSSES ARE MARGIN-INVARIANT. **LINE CLOSED.** (2026-08-25)

§4p left the read-strand fix break-even (16 correctly-dissolved antisense artefacts against 18 genuine
losses) and proposed abstention: don't flip on a bare majority, keep the `'+'` placeholder when the vote
is not decisive. The offline curve supported it — accuracy 0.9600 → **0.9934** at margin 0.90, a 6.1×
drop in wrong calls, while still flipping 0.93 of the reps a bare majority would.

**Run at margin 0.90 through the real binary. It does not work.**

| arm | families | copies | dissolved | **CORRECT** (now antisense) | **LOST** (still same-strand) | indet |
|---|---:|---:|---:|---:|---:|---:|
| OFF | 627 | 2,019 | — | — | — | — |
| ON, bare majority | 630 | 1,958 | 35 | 16 | **18** | 1 |
| **ON, margin 0.90** | **624** | 1,948 | 37 | 18 | **18** | 1 |

The margin **is** binding — single-exon `'-'` calls fall 1,926 → 1,765 as abstainers return to the
placeholder (OFF 5,928/0 · ON 2,060/1,926 · ON90 2,293/1,765). So the rule fired and did what it was
built to do.

⭐ **But the genuine losses are EXACTLY 18 in both arms.** They are **invariant to the strand-vote
margin**, so they are **not caused by low-confidence calls** — the mechanism §4p hypothesised is wrong.
Abstention changed *which* families dissolve (correct 16 → 18) without protecting a single lost one, and
left the catalog with **fewer** families than either the bare-majority arm or the baseline (624 vs 630 vs
627).

### ⟹ CLOSED. Four iterations is enough, and the stop was pre-registered

Recorded before the run: *"if both classes fall together, the honest conclusion is that this defect is
real but its correction doesn't pay for itself at this read depth — written down, not iterated past a
fifth time."* The outcome is worse than that: the loss class did not move **at all**.

**`RUSTLE_READ_STRAND` and `RUSTLE_READ_STRAND_MARGIN` both stay OPT-IN, default OFF.** The OFF gate is
byte-identical to the shipped catalog (§4o), so nothing is at risk.

### What SURVIVES from this line — the defect, not the fix

- ⭐ **A third of the rep set carried an UNMEASURED strand**: `strand.unwrap_or('+')` gave all **5,928**
  single-exon reps `'+'` and **zero** `'-'`, against a spliced split of 0.4867/0.5133. That is a real,
  named defect in the node builder and it is now documented in the source.
- ⭐ **Read orientation measures it at 0.9650** (386/400 on a junction-determined labelled set, against a
  constant-`'+'` floor of 0.4867), rising to **0.9934** at margin 0.90 and **1.0000** at unanimity.
- ⭐ **98.55%** of the orientation guard's blocked pairs involve such a rep (vs a 0.3943 base rate among
  kept pairs, 2.50×), leaving only **58** genuine antisense candidates genome-wide — all of which
  **remain rejected** in every arm.
- ⚠ **Correcting the strand is roughly break-even on the emitted catalog**, and the residual 18-family
  loss has an unknown cause that is NOT vote confidence. **That is the open question if anyone returns
  to this** — and it should be attacked directly, not through the strand rule.

---

## 4r. ⭐⭐⭐ WHY THOSE FAMILIES DISSOLVE — THE REPRESENTATIVE IS RE-PICKED, STUB → SPLICED (2026-08-25)

§4q left one open question: the dissolutions have *"an unknown cause that is NOT vote confidence."*
Answered. **It is not antisense, not vote confidence, and not the coverage floor — the family's
representative is replaced by a better one, and the `E_r` edge was formed with the old one.**

### ⛔ FIRST — §4p AND §4q's CLASSIFICATION IS RETRACTED (my error, twice)

Both sections classified a dissolved family by reading its members' strand from *"the ON rep with the
largest overlap"*. That matcher returns a **neighbouring rep** when the exact span is gone — which is
precisely the case under test. Re-derived with **exact (chrom, start, end) matching**:

| dissolved family, cause | ON (bare majority) | ON90 (margin 0.90) |
|---|---:|---:|
| **member span is NO LONGER A NODE — rep re-picked** | **32** | **32** (the same 32) |
| OPPOSITE strand in ON — a true antisense artefact | 3 | 5 |
| same-strand with **both** spans intact | **0** | **0** |
| total | 35 | 37 |

⟹ **"16 correct / 18 lost" (§4p) and "the 18 losses are margin-invariant" (§4q) are BOTH WITHDRAWN.**
There is no same-strand-both-intact class at all: it was an artefact of the overlap lookup.
⭐ **§4q's CONCLUSION STANDS AND IS STRENGTHENED** — the dominant cause is **the identical 32 families in
both arms**, and it is not vote-related, so abstention could never have addressed it.

### The mechanism, measured

**The locus survives; its representative moves.** Of the 41 changed member spans, **39/41 = 0.9512 still
have an overlapping ON rep** (median overlap of the new rep with the old span **1.0000** — it *contains*
it; median boundary shift **52,957 bp**); only **2** are truly absent from the rep set.

**And the new representative is better.** At those 39 loci the median `n_exon` goes **1 → 5**:

| | measured | comparator |
|---|---:|---|
| old rep was a single-exon **STUB**, new one is **SPLICED** | 27/39 = **0.6923** | **1,928/2,306 = 0.8361** genome-wide among ALL re-picked reps ⟹ **0.83×, NOT enriched** |
| the new rep IS a catalog copy in ON (family relabelled) | 11/39 = 0.2821 | |
| the new rep is in NO ON family (a real membership loss) | 28/39 = 0.7179 | |

⚠ The comparator matters and reverses the reading: stub→spliced is **what the strand fix does
everywhere** (83.61% of 2,306 re-picked reps), not something peculiar to the dissolved families. The 32
dissolutions are **collateral**, not the phenomenon.

⚠ An earlier comparator here was **vacuous by construction** — "families that survived" was *defined* as
all copies present at their exact OFF spans, so it can contain no re-picked rep (n = 0). Not evidence.

### ⭐⭐ THE REAL HEADLINE: the strand fix repairs 1,928 STUB REPRESENTATIVES

The ledger's largest named node-construction defect is that `pick_locus_rep` keeps ONE chain, and **46%
of representatives covering a known family member are single-exon stubs** while **53% of those loci have
a gate-passing ≥3-read spliced chain that was discarded**. Measuring the strand lets a stub group with
the spliced transcripts at its own locus (they were kept apart by a placeholder `'+'` against a real
`'-'`), and `pick_locus_rep` then selects the spliced model.

**That is 1,928 loci genome-wide whose representative goes from a stub to a spliced model** — and it is a
far larger effect than the ±3 families the catalog headline moves by. ⚠ **UNVALIDATED as an improvement**:
that a spliced rep is *better* is an inference from the stub census, not measured here. The catalog-level
consequence is 32 families dissolving.

### ⟹ What this settles, and the honest position

The read-strand line stays **CLOSED and OPT-IN** (§4q). But the reason is now precise: correcting a
placeholder strand **re-picks 2,306 representatives**, and `E_r` edges are computed on the
representative's sequence — so improving a node breaks the edges that were formed with its worse version.
**This is a property of the one-representative-per-locus design, not of the strand fix**, and any future
change to representative selection will pay the same cost. That is the finding worth carrying forward.


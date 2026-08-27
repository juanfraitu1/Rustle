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
[`o1_ledger.md`](o1_ledger.md).


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
`o1_investigations.md#census-of-incorrectly-called-families` reached for the coverage hole (47/105 cases node-construction). The real defect the
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

> ⛔ **RETRACTED IN PART — see §4r.** The dissolution classification below reads each member's
> strand from *the ON rep with the largest overlap*, which returns a NEIGHBOURING rep when the exact
> span is gone — exactly the case under test. **"16 correct / 18 lost" is withdrawn.** Exact matching
> gives **32 rep-re-picked / 3 opposite-strand / 0 same-strand-both-intact**. The 224-lost/163-gained
> churn figures and the single-exon enrichment (2.40×) are unaffected.

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

> ⛔ **RETRACTED IN PART — see §4r.** "The 18 genuine losses" inherits §4p's overlap-matching
> defect and is withdrawn; there is no same-strand-both-intact class. ⭐ **The CONCLUSION stands and is
> STRONGER**: the dominant cause is the *same 32 rep-re-picked families* in both arms, so it is not
> vote-related and abstention could never have addressed it.

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


---

<!-- merged from o1_ledger.md on 2026-08-25 -->

### The corroboration layer — external certificates for `E_r`

> ⚠ **CATALOG PROVENANCE.** Figures quoted per **/494** (families) or **/1415** (copies) describe the
> **superseded** 2026-07-17 catalog. The current default emits **627 families / 2,019 copies**. See
> [`NUMBERS.md`](NUMBERS.md) and [`o1_catalog_provenance.md`](o1_catalog_provenance.md).

**Status 2026-08-20. All candidates are T8 — offline, not through the shipped binary — except the
transcript-orientation guard, which is the default.** Companion to
[`o1_investigations.md#false-positive-hardening-rules-that-survived-falsification`](o1_investigations.md#false-positive-hardening-rules-that-survived-falsification) (the veto side) and
[`o1_ledger.md`](o1_ledger.md) (the one-page verdict list).

### 1. Why a layer, and not a better rule

Two measurements force the architecture.

**Every statistic computed from the two nodes is length in disguise.** The FP class *is* the
short-node class, so an internal statistic half-works and then collapses when length-normalised —
which is what happened to coverage on all four substrate × denominator cells, to block count
(0.7548 → 0.6410 per kb), to junction crossing (100× exon-count bias) and to read tiling.

**And the discriminating clause does not order the classes at all.** True GFPT1×GFPT2 scores
**0.5353**, false ATP1A1×ATP4A scores **0.5689** — a true pair below a false one, so no threshold on
it separates them.

⟹ **The object is one high-sensitivity RELATION plus low-sensitivity, high-precision CERTIFICATES that
mark which edges independent evidence corroborates.** Every attempt to *replace* the relation with a
certificate failed, because they are not the same kind of object: a relation must reach everything, a
certificate must be right.

### 2. The certificates

Unit = **pair**, GGO, FP arm 14 / TP arm 150. ⚠ The FP arm is **not held out** — this line of work was
derived on it — so the TP column is the load-bearing one and every 0/14 carries **Wilson95
[0.0000, 0.2153]**: bounded, not proven zero.

| certificate | references | FP admitted | TP covered | portable? |
|---|---|---:|---:|---|
| exon-sum coverage ≥ 0.50 *(the relation, for contrast)* | nothing external | 14/14 | **150/150** | — |
| genomic-span coverage ≥ 0.50 | the genome | **0/14** | 30/150 = 0.2000 | ✅ genome only |
| SD containment ≥ 0.90 | an SD catalog | **0/14** | 24/150 = 0.1600 | ❌ species-specific |
| **flank homology** | the genome | **0/14** | 30/150 = 0.2000 | ✅ **genome only** |
| **UNION of the three** | — | **0/14** | **41/150 = 0.2733** | |

#### 2a. Flank homology — the portable one

Take **5 kb outside each locus, exclude the gene bodies entirely**, and ask whether the extra-genic
sequence aligns (**≥ 300 bp at ≥ 0.90 identity**). A segmental duplication copies a *neighbourhood*, so
a real pair shares flanks; a mobile element inserted into two unrelated genes does not.

| min bp | min id | FP admitted | TP kept |
|---:|---:|---:|---:|
| 300 | 0.90 | **0/14** | 30/150 |
| 1000 | 0.80 | 1/14 | 36/150 |
| 2000 | 0.80 | **0/14** | 27/150 |

It reaches SD containment's specificity **using only the genome**, which removes the portability
objection that stopped SD containment from being folded in. ⚠ **A negative is weak** — a duplication
whose boundaries fall inside the gene, or an old one with diverged flanks, looks identical. **The
positive is the informative direction.**

#### 2b. SD containment — what it licenses

**Sayable:** for **136/627 = 21.7%** of families every member locus lies entirely inside one segmental
duplication at ≥ 0.90 — *the family and the duplicated segment are the same event.*
**Not sayable:** which copy is ancestral. With one genome every duplication edge is symmetric, so this
is `DNA_DUPLICATION`, never `DERIVED_FROM`. ⚠ Nor is later **gene conversion** excluded — containment
dates the **segment**, not the family.

#### 2c. The vetoes, for completeness

| veto | references | FP rejected | cost |
|---|---|---:|---|
| **transcript orientation** — SHIPPED DEFAULT | strand | 6/14 | 4/9,032 edges; **2.05%** genome-wide |
| genome-anchored repeat multiplicity ≥ 50 | the genome | 10/12 | 0/135 arms; **3.67%** genome-wide |

⚠ The repeat veto is a **veto, never an admission criterion**: TP median `gmult` is 2, so a criterion
strict enough to exclude repeats excludes most real paralogues.

### 3. Limits, stated plainly

* **They overlap heavily.** 30 + 24 + 30 = 84 if disjoint; **41** in fact. All three are genomic
  co-duplication evidence, not three independent lines — do not present them as convergent proof.
* **Coverage is a quarter.** 27.3% of true pairs are corroborated; the rest are **uncorroborated, not
  rejected.** The layer adds confidence where it fires and says nothing where it does not.
* **T8.** None of the certificates has run through the shipped binary. The orientation *veto* has.
* **0/14 is bounded, not zero**, and the arm is not held out.

### 4. Reproduce

```bash
python3 bench/o1_flank_extension.py          # flank homology
python3 bench/o1_sd_anchor.py 0.90           # SD containment, identity floor as argv
python3 bench/o1_sd_certificate.py           # duplication-mechanism classes on a catalog
python3 bench/o1_substrate_denominator.py    # genomic-span coverage
```

---

## 4s. ⛔⛔ REPRESENTATIVE SELECTION IS CLOSED — `E_r` HAS NO UPGRADE PATH (2026-08-25)

§4r ended by naming one-representative-per-locus as the root cause and pointed at the stub defect.
**Attacked directly. NO-GO on changing the representative; the defect ships as a FLAG.**

### ⭐ The carrying number: 0/215 (verified independently by the controller)

For every OFF `E_r` edge whose **shorter** endpoint is a single-exon rep that has a **unique containing
spliced model at its own locus** (n = **215** edges, from **1,632** such loci), does that spliced model
form an `E_r` edge with the stub's partner?

| | |
|---|---:|
| **the container inherits the stub's partner edge** | **0/215 = 0.0000** [Wilson 95% upper 0.0170] |
| *"the edge is a property of the LOCUS"* predicts | 215/215 |
| *"the edge is a property of the exact REPRESENTATIVE SEQUENCE"* predicts | ~0.01 (base rate 4,778/C(17,923,2) = **2.97e-05**) |

**Two must-pass validity checks, both run before believing it:**
- **1,416/1,416 = 1.0000** of the containers are themselves OFF representatives, so they sat in the same
  all-vs-all and *could* have formed the edge. (Had this been ~0, the result would be a lookup artefact.)
- The containers are not inert: **116/215** have ≥1 `E_r` edge of their own; genome-wide 151/1,416 =
  0.1066 have degree > 0.

⟹ **`E_r` has NO UPGRADE PATH. Any representative change rewrites the edge set from scratch.** The
pre-registered kill was *"if the spliced model inherits < 20% of the stub's edges, representative
selection is closed"*. Observed **0.0000**.

### ⛔ And my brief was wrong: "is a spliced rep better?" HAS been measured — it regressed

I told the investigation this had *"never been measured"*. Wrong in both directions:

| level | measurement | result |
|---|---|---|
| **partition** (the one that counts) | `RUSTLE_SPLICED_REP`, **two real-binary e2e runs** (register row 357) | chr7 family-detection **F1 0.570 → 0.411**; chr16 **0.910 → 0.761**; **loses NPIPB12** |
| node (association, not intervention) | register row 990: in-band fraction at spliced-rep loci 69% (n=45) vs stub-rep loci 28% (n=32) | 2.46× — but this stratifies ONE catalog by an outcome-correlated covariate |
| counterfactual on the arms on disk | per-copy **annotated-intron precision** OFF **0.9700** (13,415/13,830) vs ON **0.9714** (13,243/13,633) | **a null** |

⟹ the premise is true at node level, where it does not matter, and **false at partition level, where it
does**.

### ⚠ The generalisation: the edge set sits within 1.49× of the coverage floor

Projecting a uniform **1.49×** inflation of the shorter rep's exonic length (the *observed median*
stub→container ratio; q1 1.14, q3 2.04, n = 1,632) onto the whole OFF edge set:
**2,891/4,778 = 0.6051 of edges fall below the 0.50 coverage floor.** The stub class is not special —
110/215 = 0.5116 against that 0.6051 comparator.

⭐ **This restates `RUSTLE_LOCUS_EXON_UNION`'s −20 recall points as a general property of the coverage
clause: ANY intervention that lengthens a representative pays it**, whether the longer sequence is
selected or constructed. That prices `RUSTLE_COTHREAD_REP` too — it constructs a max-weight path, so it
escapes the selection entailment, but not this one. ⚠ Its documented Soto gains (correctly-sized loci
15 → 26; SRGAP2 0.05 → 0.99) are about **locus SIZE**, a different objective needing its own
pre-registration — do not run it as a stub fix.

### ✅ MODIFIED-GO — ship the defect as a FLAG, not a fix

Emit **`frac_pure_intron`** per catalog copy: the fraction of a copy's asserted exonic bp lying inside an
annotated intron with **no** annotated exon overlap, against `GGO_genomic.gff` (277,703 distinct introns).
Genuinely new and not a closed route in disguise — it touches no node, no edge, no denominator and no
partition; **emission ≠ definition**, so P1 is untouched.

| | measured | how to quote it |
|---|---:|---|
| single-exon copies placing > 50% of asserted bases in pure intron | 149/431 = **0.3457** vs 37/1,450 = 0.0255 multi-exon | 13.55× — **unclustered upper reading only** |
| clustered by family | 55/185 = 0.2973 vs 31/551 = 0.0563 | 5.3× |
| **within-family paired**, 118 families holding both classes | **22/118 = 0.1864 vs 16/118 = 0.1356** | ⭐ **1.37×** — *this* is the effect size (sign test 20 vs 7, 91 ties) |

Catalog burden **149/2,019 = 0.0738** (floor) to **259/2,019 = 0.1283** (ceiling). Corroborates §4g's
45.06% intron rate by an independent route.

⚠ **Mandatory caveat, every time**: differential exclusion — **110/541 = 20.3%** of single-exon copies
fall outside annotated space vs **28/1,478 = 1.9%** of multi-exon; and `GGO_genomic.gff` holds
**0 NOTCH2NL\*, 1 SRGAP2, 1 LIMS1**, so it is **structurally blind at the named flagship families**.

### ⟹ The thesis statement about node construction

**The stub defect is real, it is measurable, and it is not fixable by choosing a different
representative** — because `E_r` edges are properties of the exact representative sequence (0/215), and
because 60.51% of the edge set is within 1.49× of the coverage floor. O1 should **disclose** it
(`frac_pure_intron`) rather than claim to have solved it.

---

## 4t. ⛔⛔ THE SCOPED ORIENTATION GUARD — SOUND IN PRINCIPLE, CATASTROPHIC IN EFFECT (2026-08-25)

§4s closed representative selection. This closes the remaining route into the stub defect: **relaxing the
orientation guard where it has no information.** Implemented, then **killed offline** — no catalog run
was spent.

### Why it looked right

Stratifying §4s's 0/215 by whether the "container" contains the stub's **exons** (not merely its span):

| stratum | n | dominant failure |
|---|---:|---|
| **genuine exonic upgrade** | 124 | **123/124 = 0.9919 the ORIENTATION GUARD** |
| stub lies in the container's **introns** | 91 | 84/91 = 0.9231 **no alignment** — different objects |

and **123/129 = 0.9535** of those guard blocks have a **single-exon (placeholder-strand) partner**.
⟹ for real stub→spliced upgrades the blocker is **the guard acting on an unmeasured field, on the
PARTNER side** — not representative choice, and not coverage (**1/215 = 0.0047**).

⛔ **This CORRECTS §4s**, which named the coverage floor as the mechanism from a *projection*
(60.51% below floor at 1.49× inflation). Measured on the actual pairs, coverage is **0.5%**.
The projection was hypothetical; the guard is what fires.

### Why it is fatal

Exempting minus-strand records where a rep's strand was never measured **and** the spans are disjoint
(both clauses needed — the unmeasured clause alone admits 3,951 pairs of which **1,727 = 0.4371** have
overlapping spans against **0.0109** in the shipped set, **40×**) leaves **2,224** pairs. Then:

| | measured |
|---|---:|
| SEDEF support vs the shipped edge set | **0.2350 vs 0.3025 = 0.78×** |
| pairs joining two reps that are **already catalog copies** | 1,913/2,224 = 0.8602 |
| …of those, joining **two DIFFERENT families** | **1,912/1,913 = 0.9995** |
| distinct family pairs fused | **112** |
| catalog touched | **77/627 = 0.1228** |
| ⟹ collapses into | **17 blocks, the largest swallowing 43 families** |
| pairs that would form a genuinely NEW family | 77/2,224 = **0.0346** |

**The trade is 1,912 family merges to gain 77 new families.** That is the hairball pathology, arrived at
from a new direction.

### ⚠⚠ THE INSTRUMENT LESSON: the negative panel is BLIND to this change

The pre-registered cheap test was the HUMAN 150-window panel. It ran and returned **false-merge 2/150 =
0.0133 and 3 E_r edges — both unchanged.** ⛔ **That is NOT a pass.** The panel offers **0 qualifying
pairs**: its windows are gene-tight and single-locus, so all 25 of its guard-blocked pairs are
**co-located**, which the disjoint-span clause excludes by construction.

⭐ **An acceptance instrument that cannot express the change under test returns a null for the wrong
reason.** Check the candidate count BEFORE reading the verdict — the same discipline as the mutation
test that proved the read-strand ingest was wired.

### ⟹ What this settles about the orientation guard

**The guard is doing correct work for a wrong-sounding reason.** §4n established that it filters on an
unmeasured field — **98.55%** of what it blocks involves a placeholder strand. That remains true, and it
is *still* not a reason to relax it: those pairs are overwhelmingly **cross-family bridges**, and
admitting them on the principled ground that the field is unmeasured would destroy the partition.

⟹ **The stub defect has no remaining route.** Representative selection is closed (§4s, 0/215); the
guard cannot be scoped (here); exon-union, shared-exon and the hierarchy were already NO-GO. **O1 should
DISCLOSE the defect (`frac_pure_intron`, §4s) and not claim to have solved it.** `RUSTLE_ER_GUARD_SCOPED`
ships **off**, retained only as the record of this experiment.

---

## 4u. ⛔ WGS CANNOT GIVE A COPY-NUMBER BASELINE FOR O1 — AND THE ONE IT COULD GIVE NEEDS NO WGS (2026-08-25)

Asked whether the downloaded WGS can inform the O1 definition, or at least supply an expected copy
number per family. **Killed at the design stage, before any counting.**

### ⭐ The physics: family paralogues are too diverged for exact k-mer matching

Measured on the shipped edge dump (`o1_strand/off/dump/ggo.edges.tsv`), within-family `E_r` edges,
n = 3,715 — **verified by the controller, not re-derived**:

| within-family identity | ⟹ P(a 21-mer is identical) |
|---|---:|
| **median 0.8234** | **0.0169** |
| q1 0.7635 | 0.0035 |
| q3 0.8922 | 0.0911 |

**1,524/3,715 = 0.4102 of within-family edges join copies sharing under 1% of their 21-mers.**

⟹ a k-mer multiplicity readout returns **CN ≈ 1 for a copy regardless of its family size**, because a
copy's k-mers are essentially private to it. ⚠ This is not a tunable: lowering k to recover sharing
destroys specificity, and the identity distribution is a property of the biology, not of the method.
A supporting probe reports **216/627 = 0.3445 of families share ZERO genomic 21-mers across their
copies** (48.01% share ≤ 10) — consistent with the arithmetic above. ⚠ that figure is UNVERIFIED (its
verifier died on a session limit); the identity table is the controller-verified part.

### ⭐⭐ And the baseline the question actually wants is FREE — the WGS is not needed

**The WGS animal IS the assembly animal** (both SAMN04003007 / KB3781). So the decomposition is

> **CN_WGS = CN_assembly + collapse_deficit**

- **CN_assembly** — how many times a copy's sequence occurs in `GGO.fasta` — needs **no WGS at all**, and
  the machinery is already in-tree (`min_shared_gmult`: genome-anchored, seed-invariance **0/147**,
  AUC 0.9429). ⚠ Any such count must record **`-p` and `-N`** (MAPKBP1 gives 1/1 at `-p 0.8` and 9/8 at
  `-p 0.1`).
- **collapse_deficit** — the WGS-only residual — **is O3's ABSORBED stratum** (64.2% of excised copies,
  depth 1.75×), and the project's own O1/O3 boundary puts it outside O1: *O1 finds copies missing from
  the ANNOTATION; O3 finds copies missing from the GENOME.* The o3 work pre-registers **< 1 collapse**
  in this compartment, so the residual is expected to be empty.

⟹ **everything the WGS adds beyond the assembly is, by this project's own boundary, not O1's.**

### ⛔ NOT in the definition — three independent reasons

1. **It measures the wrong thing** (above): CN ≈ 1 for diverged paralogues.
2. **It would break substrate portability.** The cross-substrate replication that gives O1 its strongest
   generalisation claim — **87.06% of `E_r` edges recovered** across a different ANIMAL *and* TISSUE —
   ran on a fibroblast sample with **no matching WGS**. A definition consuming DNA copy number could not
   be evaluated on its own best evidence.
3. **It repeats a closed entry.** "Joint DNA/RNA" was already RETRACTED as a definition and survives as a
   PROPERTY; the measured gain there was **CONTIGUITY, not jointness** (a length-matched DNA-only arm
   recovered 91/91).

⭐ The defensible use is the third of the three the design considered: **an external validation set,
never consumed by the pipeline** — the same standing as `frac_pure_intron` (§4s) and the λ certificate,
under the settled principle **emission ≠ definition**.

### ⚠ And the comparison would be confounded even if it worked

O1 sees only **EXPRESSED** copies, from **testis of a different animal**; DNA sees **all** copies. So
**O1 < DNA is EXPECTED and is not evidence of an O1 defect.** A probe reports that the stratum where the
comparison would be clean (core ≥ 200 shared k-mers AND no copy in an SD) contains **2 families, none
with ≥ 3 copies** ⚠ (unverified — same session-limit loss). Even granting the method, there is nothing
to measure on.

### ⟹ What to do instead, if a baseline is wanted

Count each catalog copy's occurrences **in the assembly** (`GGO.fasta`), recording `-p` and `-N`, and
publish it as an **annotation** beside each family — not as a term in the definition, and not from the
WGS. Cost: one genome pass, no k-mer database, no 42.5 GB of counting.


## §4v — the completeness deficit is REAL, not a repeat artefact (2026-08-26)

**Substrate:** the 3-contig full-depth fibroblast NPIP catalog (2,847 reps / 83 families). ⚠ NOT
comparable to the genome-wide testis catalog; whole-contig subset, so off-contig partners cannot form
edges.

**The finding under test (§4v.1).** Edge formation falls as the transcript model gets more complete:
single-exon **227/1,234 = 0.1840**, 2–4 exons 0.2706, 5–14 exons 0.1560, **≥15 exons 16/276 = 0.0580**.
Restricted to models both complete and well supported (≥15 exons AND ≥100 reads): **7/153 = 0.0458**.
Mechanism: coverage ≥ 0.50 of the SHORTER sequence on a scale-free denominator — a stub needs half of
~2 kb, a 21-exon model half of ~30 kb across diverged UTRs and alternative exons. **Completeness raises
the bar the model itself must clear.** Sharpest single case: `NC_073242.2:29391569-29415846`, 21 exons,
**4,187 reads**, `degree = 0`.

**The attack (§4v.2): are stub edges just repeats?** If stubs pair via Alu/L1-class sequence the deficit
would be an artefact of junk, not of the coverage clause. Measured genome multiplicity of every rep —
`minimap2 -x asm20 -c --eqx -N 200 -p 0.1` against the three contigs, a hit being ≥90% identity and
≥50% coverage OF THE REP, distinct non-overlapping target loci counted. **`-p` and `-N` recorded.**

⚠⚠ **THE CROSS-CLASS COMPARISON IS INVALID AND IS RETRACTED HERE BEFORE USE.** An exon-sum is a
CONCATENATION of exons; it cannot align contiguously to genomic DNA across an intron, so the ≥50%
single-record coverage filter systematically fails for spliced reps — hence their median multiplicity of
**0**, which is a measurement artefact, not biology. **Only the within-stub arms are comparable**, both
being single-exon where exon-sum == genomic segment.

**Result (valid arm).** Among single-exon reps:

| arm | n | median mult | q90 | frac ≥ 20 loci |
|---|---:|---:|---:|---:|
| **with an edge** | 227 | 1 | 6 | **0.0881** |
| **no edge** | 1,007 | 1 | 1 | **0.0000** |

⭐ **The majority of stub edges are NOT repeat-driven** (median multiplicity 1 — family-sized, not
repeat-sized).

⛔⛔ **RETRACTED 2026-08-26 — the second half of this section was WRONG.** It read the 20/227-vs-0/1,007
contrast as "high multiplicity is diagnostic of a stub edge", called those edges removable junk, and
proposed the veto as a free precision gain. **All three claims are refuted by §4x.** Multiplicity and
degree measure THE SAME THING: a rep with N near-identical genomic copies has N potential partners.
Stratified, the stub edge rate is **mult 1 → 0.1443 · 2–4 → 0.4574 · ≥20 → 1.0000 (20/20)** — a perfect
rate, because those are the highest-copy-number loci, which is precisely what O1 exists to find. The
"junk" was NPIP.

**Verdict.** Applying the veto takes the stub rate **0.1840 → 207/1,234 = 0.1677**, so the deficit
**SURVIVES**: **3.67×** against well-supported spliced models (was 4.02×), **2.89×** against all ≥15-exon
reps (was 3.17×). ⟹ **the completeness deficit is a property of the EDGE RULE, not of repeat contamination**,
and ~9% of stub edges are separately removable junk.

**Consequence for O1's claim.** The definition's detection power is **anti-correlated with evidence
quality**. This is not "rigid" — a rigid definition would miss hard cases while keeping easy real ones.
This one keeps the DEGENERATE cases and drops the well-resolved ones. Every pure-NPIP family in this
catalog (`GWFAM30` 6/6, `GWFAM32` 3/3, `GWFAM31` 2/2) is built from single-exon stubs.
⟹ **claim reproducibility and falsifiability, NOT adequacy.** See `docs/THESIS_OBJECTIVES.md`.

## §4w — the shared-exon edge rule does NOT fix the completeness deficit (2026-08-26) — REFUTED

**Motivation.** §4v established that edge formation is anti-correlated with model completeness, mechanism
= the scale-free coverage denominator (half of the SHORTER sequence). The ledger's own prescription was
"change the DENOMINATOR, not the threshold". `shared_exon_edges` is the one in-tree implementation of
that and its precision had never been measured: exons are matched **any-to-any** with an **absolute
aligned-bp floor** (`>= 100 bp`), so no sequence lengthens and no denominator moves. It also does not
touch representative selection, so it cannot regress the way `RUSTLE_SPLICED_REP` did.

**Arms.** Same substrate, same binary, same 2,847 reps (node construction untouched — confirmed identical
in all three). OFF = incumbent (coverage >= 0.50 of shorter, id >= 0.60). SE98 = `RUSTLE_SHARED_EXON=1`
at its default id >= 0.98. SE60 = the same at **id >= 0.60**, matched to the incumbent, which is the arm
that isolates the denominator. ⚠ SE98 alone is CONFOUNDED (it moves identity AND denominator together);
it is reported only for shape.

**Endpoint: reaches a family, as a fraction of all reps in that exon class.**

| class | reps | OFF | SE98 | SE60 |
|---|---:|---:|---:|---:|
| single-exon | 1,234 | 0.1718 | 0.0429 | 0.0900 |
| 2–4 exons | 510 | 0.2627 | 0.0235 | 0.1431 |
| 5–14 exons | 827 | 0.1487 | 0.0278 | 0.0798 |
| **≥15 exons** | 276 | **0.0543** | **0.0326** | **0.0399** |
| **deficit** | | **3.16×** | 1.32× | 2.26× |

⛔⛔ **THE RATIO IMPROVES WHILE THE TERM WE CARE ABOUT GETS WORSE.** The deficit narrows 3.16× → 2.26× →
1.32×, and quoting that alone would read as a fix. It is not: the **≥15-exon rate FALLS in every arm**
(0.0543 → 0.0399 → 0.0326; absolute copies **15 → 11 → 9**). The deficit closes because the *stub* rate
collapses (0.1718 → 0.0900 → 0.0429), not because complete models do better.
⟹ **a ratio of two rates is not a valid endpoint when both are free to fall.** Same family as the
already-registered "ratio-to-truth ⟹ use the IN-BAND FRACTION, never the median".

**Every pre-registered endpoint moves the wrong way** (OFF → SE60): families **83 → 59**, copies
**484 → 261**, NPIP loci in a family **12/31 → 10/31**, 100%-NPIP families **3 → 1**, largest family
39 → 22. There is no arm in which complete models gain.

**Verdict: REFUTED.** Replacing the fractional denominator with an absolute aligned-bp floor does not
rescue complete transcript models; it is uniformly stricter and costs recall and NPIP purity together.
⟹ **the designed denominator fix is now MEASURED and dead**, closing the route §4s/[[project_o1_single_exon_stubs]]
left open ("precision UNMEASURED"). The completeness deficit of §4v stands with **no known repair**.

⚠ Candidate count checked before reading the verdict: SE60 linked **966 locus pairs** from 11,385 exons
(SE98: 164) against OFF's 1,847 E_r edges — the instrument had ample candidates, so this is a real
negative, not a blind panel.


## §4x — the repeat-hub veto is REFUTED for the catalog path (2026-08-26)

**What was built.** `family_define` has carried a repeat-hub gate (`min_shared_mult >= 20`) ON BY DEFAULT
since the annotation era; the de novo catalog path never had one. Ported it as
`RUSTLE_ER_REPEAT_GATE=<M>`: per-rep genome multiplicity (distinct non-overlapping reference loci at
>= 90% identity, >= 50% coverage of the rep, `-p 0.1 -N 200` recorded), vetoing an edge when BOTH
endpoints are hubs. **Genome-anchored, so P1 holds by construction** — multiplicity is a property of the
rep and the reference alone, unlike R5's catalog-anchored count which broke seed-invariance at 94/147.

**Discipline checks PASS.** Flag unset ⟹ `copies.tsv` md5 `2c002d7c…` and `families.tsv` `03e8497e…`
**byte-identical** to baseline, and `repeat_hub_gate` is emitted into the params certificate (the M2
defect, shipped twice before, does not recur).
⚠ **The first ON run was a silent no-op** — `intron_fasta` is populated only under
`--homology-genomic-span`, so the veto had no reference and returned the baseline exactly. **"Identical
to baseline" would have read as "safe and inert" had the veto not logged its own skip.** Same failure
class as the human 150-window panel scoring 2/150 on zero qualifying pairs. **Make a no-op announce
itself.**

**Result at M=20 (veto fired: 19/2,847 reps are hubs; edges 1846 → 1826).**

| | OFF | GATE20 |
|---|---:|---:|
| families | 83 | 82 |
| copies | 484 | 473 |
| **NPIP loci in a family** | **12/31** | **7/31** |
| **100%-NPIP families** | **3** | **1** |

⛔⛔ **REFUTED — ACTIVELY HARMFUL ON THE TARGET FAMILY.** All 11 lost copies are single-exon and
**7/11 sit at NPIP loci**. The gate destroys 5 of 12 recovered NPIP loci and 2 of 3 pure NPIP families to
remove 20 edges.

⭐⭐⭐ **WHY, AND IT GENERALISES: GENOME MULTIPLICITY CANNOT DISTINGUISH A REPEAT FROM A HIGH-COPY GENE
FAMILY — AND HIGH-COPY GENE FAMILIES ARE O1'S ENTIRE SUBJECT.** NPIP has 31 loci, so an NPIP exon-sum
trips any "occurs in ≥ 20 places" test **by being what it is**. The gate's predictor is confounded with
the target it is supposed to protect. ⟹ **the veto's apparent signal in §4v was TAUTOLOGICAL**: mult ≥ 20
⟹ edge rate 1.0000, because both quantities ask "does this sequence have paralogues?".

⚠ Why it is nevertheless correct in `family_define`: that path is annotation-driven, so a named gene is
never at risk of being reclassified as a repeat by its own copy number. **Do not port gates across paths
on the strength of the threshold agreeing.**

**Disposition:** flag RETAINED, default OFF, **DO NOT ENABLE** — it is the record of the experiment, like
`RUSTLE_ER_GUARD_SCOPED`. ⟹ **§4v's "free precision gain" does not exist; there is no queued win.**

## §4y — reference bias is NOT the binding constraint; DEPTH is (2026-08-26)

**The claim under test** (advisor, 2026-08-26): *all copies are found only when the fibroblast reads are
aligned to the T2T assembly of that same animal; align any other individual and copies are missed.*

**Direct test on NPIP** — the same 31 homology-defined loci scored against three catalogs:

| substrate | animal | depth / scope | NPIP loci in a family |
|---|---|---|---:|
| testis OR6737 | **different** | genome-wide | **3/31** |
| fibroblast KB3781 | **same as assembly** | 41% subsample, genome-wide | **1/31** |
| fibroblast KB3781 | **same as assembly** | full depth, 3 contigs | **12/31** |

⛔ **THE PREMISE IS FALSE TWICE OVER.** (1) Same animal, full depth, aligned to its own assembly recovers
**12/31 — not all**. (2) The same animal at reduced depth recovers **1/31, WORSE than a different
animal's testis at 3/31** ⟹ individual identity is not the ordering variable.

⭐ **The clean comparison is rows 2 vs 3: same animal, same tissue, same pipeline, same 31-locus
denominator — a 12× difference from READ DEPTH ALONE.** Row 3 is if anything handicapped (a 3-contig
catalog cannot form off-contig edges; all 31 NPIP loci are on those contigs).

⚠ Rows 1 and 3 differ in tissue AND scope, so only the depth comparison is clean. Corroborating, from
§4l: **87.06% of E_r edges replicate** across a different animal *and* tissue. Reference bias is real and
second-order; O1's measured reach is **0.55** even on favourable data.

### Where the reads of a missing copy go — measured, not assumed

Whole-genome excision, 162 families:

| fate | share | detail |
|---|---:|---|
| **ABSORBED** — redistributed onto surviving paralogues | **104/162 = 64.2%** | leaves a **depth ghost of 1.75×** |
| **ORPHANED** — become unmapped | **54/162 = 33.3%** | median **92.73%** of that copy's reads unmapped |

⚠ **UNIT TRAP (T12):** the **34.53%** figure that circulates is a **READ** fraction; the 33.3% is a
**COPY** fraction. They nearly coincide and are different units — always say which.

### Can the redistribution be reverse-engineered into a "there was probably another copy" call?

**Partly, and the limit is DIVERGENCE, not abundance.** The held-out S2 detector runs at
**TPR 0.2703 / FPR 0.0200**, but its power is **0.4500 above 0.01 divergence vs 0.0588 below**, and
**45.78% of true positives lie below 0.01** ⟹ recently-duplicated near-identical copies are close to
invisible — exactly the class that gets collapsed. The unmapped-read route is separately **VACUOUS** for
this stratum (π = 1/35, 0/26 at cov ≥ 0.8): an ABSORBED copy is not orphaned, so there is no pile to find.

⟹ **the unexploited signal is the 1.75× depth ghost**, and nothing currently keys on it. See §4z for why
a depth caller cannot be built from the WGS on hand.

## §4z — a WGS depth caller for the ABSORBED stratum: sound design, BLOCKED ON DEPTH (2026-08-26)

**The idea.** §4y leaves the 1.75× depth ghost as the only unexploited signal for the **64.2% ABSORBED**
stratum — the copies that both existing routes miss (the unmapped route is vacuous there; the S2 detector
runs at 0.0588 TPR below 0.01 divergence, and 45.78% of positives are below 0.01). If a collapsed copy's
reads pile onto its surviving paralogue, that locus should carry **two copies' worth of DNA depth**, and a
caller could flag it.

⚠ **§4u does NOT refute this.** §4u killed **k-mer multiplicity → CN** (paralogues at median identity
0.8234 share ~1.7% of their 21-mers, so a copy's k-mers are private and CN reads as 1). **Depth is a
different mechanism** — it needs reads to MAP to the locus, which is exactly what "absorbed" means. The
two routes fail and succeed independently.

### ⛔→✅ Blocker 1 WITHDRAWN SAME DAY: the 0.65× was a TRUNCATED DOWNLOAD, not the data

⚠⚠ **RETRACTED 2026-08-26, hours after writing.** The measurement below was taken on the local file and
reported as the RUN's depth. It is not: `gzip -t` returns **`unexpected end of file`**, and ENA reports
`SRR26152581` at **34.62 Gbp = 11.39×**. The local copy holds **1.98 Gbp = 5.7% of the run.**
⟹ **"the WGS is too shallow" was a statement about a broken download.** **`gzip -t` EVERY archive
BEFORE deriving a number from it.**

⭐⭐ **The real picture: 108 runs exist under SAMN04003007, 102 of them WGS, totalling 1,389 Gbp =
456.91× of a 3.04 Gb assembly** — **18× the requirement**, not 1/38th of it.

| platform / model | runs | Gbp | coverage |
|---|---:|---:|---:|
| OXFORD_NANOPORE / PromethION | 6 | 561.2 | 184.61× |
| PACBIO_SMRT / Sequel II (HiFi) | 22 | 375.6 | 123.56× |
| ILLUMINA / NovaSeq 6000 | 8 | 279.7 | 92.00× |
| ILLUMINA / HiSeq 2500 | 3 | 92.1 | 30.30× |
| PACBIO_SMRT / RS II (CLR) | 62 | 44.2 | 14.52× |
| ILLUMINA / NextSeq 500 | 1 | 36.2 | 11.92× |

**Two Sequel II HiFi runs reach the target**: `SRR26152597` (13.70×) + `SRR26152596` (13.39×) =
**27.09×**. HiFi is the right substrate here — the whole problem is resolving collapsed paralogues, where
long reads map unambiguously and short reads do not.

⚠ The ENA metadata caveat still holds: `sample_title` LIES for every run under this BioSample (it labels
fibroblast IsoSeq as "Y flow-sorted DNA"). The STRUCTURED fields used above — `library_strategy`,
`instrument_model`, `base_count` — are the trustworthy ones.

**The original (WRONG) measurement, kept as the record:** `SRR26152581_subreads.fastq.gz` locally holds
163,171 reads / 1,978,946,468 bp, mean 12,128 bp ⟹ 0.65×.

Power for a 1-copy vs 2-copy call, stated as arithmetic rather than opinion. Reads overlapping a locus of
length *L* per 1× is `(L + readlen)/readlen`; separation of *n* vs 2*n* Poisson counts is
`z = n/sqrt(3n)`:

| locus | reads/locus at 1× | coverage for n≥10 (z=1.83) | for n≥30 (z=3.16) |
|---|---:|---:|---:|
| 1.5 kb | 1.10 | 9.1× | **27.3×** |
| 3 kb | 1.20 | 8.3× | **25.0×** |
| 10 kb | 1.66 | 6.0× | 18.0× |
| 30 kb | 2.99 | 3.3× | 10.0× |

⟹ **z = 1.83 at n = 10 is NOT separable**; a usable call needs n ≈ 30, i.e. **~25× at the 1.5–3 kb loci
that dominate** (24.88% of gorilla copies are ≤ 2 kb). At 0.65× a 3 kb locus expects **0.78 reads**.
**Shortfall ≈ 38×.**

### ⛔ Blocker 2: the compartment is pre-registered EMPTY

Per §4u, `CN_WGS = CN_assembly + collapse_deficit`, the WGS animal IS the assembly animal
(SAMN04003007/KB3781), and the residual — which is exactly this absorbed stratum — is pre-registered at
**< 1 collapse**. So even at 25× the expected yield **in this genome** is ~0. A null result would confirm
the pre-registration, not validate the caller.

### ⚠ And RNA cannot substitute

The 1.75× ghost was measured on RNA (`o3_excise/panel_primary.bam` is `splice:hq` FLNC) in a **PAIRED**
design — same locus, reference with and without the sibling, expression held constant. **A caller gets no
counterfactual**: given one BAM it must compare observed depth to an EXPECTATION, and on RNA that
expectation is expression, which is unknown and varies over orders of magnitude. The pairing is what made
the measurement work and is precisely what a caller cannot have.

**Verdict (REVISED): the DEPTH blocker is GONE; only Blocker 2 remains, and it is scoped.** 27.09× is two
downloads away. Blocker 2 says the NATIVE compartment is pre-registered empty (< 1 collapse), so the
caller must be **built and validated on the excision ablation** — 162 arms, 104 absorbed, known ground
truth — where the positives are constructed and plentiful. Running it on the native genome afterwards is
a TEST OF THE PRE-REGISTRATION, and a null there is a confirmation, not a failure.

**Cost to proceed:** ~82 Gbp of `fastq.gz` to fetch (order 25–30 GB compressed; 579 GB free on
`/mnt/linuxdisk`), then a genome-wide `minimap2 -x map-hifi` of ~76 Gbp on 4 threads — a multi-hour to
overnight run, and the largest single compute commitment in this project so far. ⚠ ONE AT A TIME.

## §5a — a DEPTH CALLER FOR THE ABSORBED STRATUM WORKS, and it is complementary to S2 (2026-08-26)

**First positive result on the missing-copy problem in this session.** §4y left the 1.75× depth ghost as
the only unexploited signal for the **64.2% ABSORBED** stratum. Built and measured it in simulation.

### Design — the reason it is not circular

> **Reads are simulated from the COMPLETE `GGO.fasta` and aligned to the MASKED `GGO.masked.fasta`.**

The redistribution therefore happens inside the aligner and is OBSERVED, not modelled. Simulating from
the masked genome would assume the answer. Library matched to the real one: mean **17,266 bp**, identity
**89%**, both measured from `SRR26152581` against this same reference (⚠ those reads are **CLR, not
HiFi** — a stored note saying otherwise is wrong; two presets agree at identity 0.8376–0.8910).

**Ground truth:** 162 excised copies / 6.85 Mb, recovered by diffing masked vs original (⚠ **NOT** from
N-runs: the assembly's own gaps are 2.00 Mb of the 8.85 Mb of N, and a line-seam bug in a first attempt
reported **645 Mb** of "excision" — 18% of the genome — which looked plausible. **Verify an interval is
actually all-N before trusting it.**)
**Positives:** 113 landing sites (mean identity 0.9482). **Negatives:** 226 **size-matched** controls.

⭐ **UNPLANNED VALIDATION:** sequence alone predicts **113/162 = 69.8%** have somewhere to be absorbed;
the real-read excision measured **104/162 = 64.2% ABSORBED**, and the 49 with no landing site match the
54 observed ORPHANED. Prediction and observed fate agree without being fitted to each other.

### The caller, and what a real one can actually use

⚠ **The paired ratio is NOT available to a real caller.** Masked-vs-complete depth gives a beautiful
signal (landing sites median **1.2744**, controls **1.0000** at every quartile, MWU **p = 1.35e-47**) but
it needs the complete genome, which is the thing we do not have. **The caller must compare a locus to
BACKGROUND estimated from the same genome.** All numbers below do that.

| coverage | background | TPR@1.5× | FPR | AUC |
|---:|---:|---:|---:|---:|
| 5× | 5.41 | 0.4867 | 0.1416 | 0.7011 |
| 10× | 11.00 | 0.4159 | 0.0575 | 0.7337 |
| **15×** | 16.53 | **0.3982** | **0.0177** | 0.7526 |
| 25× | 27.36 | 0.4071 | 0.0177 | 0.7758 |
| 40× | 44.32 | 0.4248 | 0.0000 | 0.8034 |

⭐⭐ **TPR IS FLAT IN COVERAGE; DEPTH BUYS PRECISION, NOT SENSITIVITY.** This **refutes my own
pre-registered ~25× requirement**, which came from a per-READ Poisson argument (n ≈ 30 for z = 3.16).
That framing was wrong: depth is measured **per base over a ~7 kb locus**, which aggregates far more
information than a read count. ⟹ **15× is the knee**, and **ONE WGS run (SRR26152597 = 13.70×) suffices**
— the second download was cancelled on this result, saving ~29.6 GB and ~12 hours.

### ⭐⭐⭐ The finding that matters: it is strongest where S2 is blind

| divergence to landing site | n | depth TPR@1.5× | S2 detector |
|---|---:|---:|---:|
| **< 0.01 (near-identical)** | 15 | **0.7333** [Wilson95 0.4805–0.8910] | **0.0588** |
| 0.01–0.05 | 49 | 0.3265 | — |
| 0.05–0.10 | 36 | 0.3889 | — |
| ≥ 0.10 | 13 | 0.5385 | 0.4500 |

**12.5× better on exactly the stratum S2 cannot see**, and S2's 0.0588 lies outside the Wilson interval,
so the comparison survives n = 15. **The mechanism explains the direction:** a near-identical copy's
reads redistribute MOST completely, so the depth ghost is strongest precisely where sequence divergence
offers nothing. ⟹ **the two detectors are COMPLEMENTARY, not competing** — S2 above 0.01 divergence
(0.4500), depth below it (0.7333). Overall single-genome: **TPR 0.4248 / FPR 0.0000, AUC 0.8034** vs S2's
**0.2703 / 0.0200**.

### ⚠ What is NOT established

- ⚠⚠ **FPR 0.0000 IS OPTIMISTIC AND MUST NOT BE QUOTED AS A REAL-DATA FIGURE.** The simulated controls sit
  at ratio **exactly 1.0000 at every quartile** — that is a property of uniformly-sampled simulated reads,
  not evidence about real libraries. Real WGS depth varies with GC content, mappability and library prep,
  all of which inflate FPR by an amount that is **UNMEASURED**. The real-data arm was cancelled
  2026-08-26 (see `wgs/README_WGS.md`), so this gap is open by choice, not by oversight. **TPR and AUC are
  far less exposed** — they rest on the depth ghost, whose magnitude the real-read excision independently
  measured at 1.75×.
- **The 1.5× threshold is NOT held out** — it was chosen as the point where FPR reaches 0 in this data.
- **FPR 0.0177 is not free genome-wide**: on a 20,000-locus scan that is ~350 false calls. ⟹ **apply the
  caller to a CANDIDATE SET (known multi-copy families), not as a genome-wide sweep.**
- The native gorilla compartment is still pre-registered at **< 1 collapse** (§4u), so running this on the
  real WGS is a **test of that pre-registration**, where a null CONFIRMS rather than fails. **UNRUN** —
  the download was cancelled at 7.53/28.6 GB once the simulation had answered the design question.
  Resume instructions, checksum and the CLR/Y-flow-sorted warnings are in `wgs/README_WGS.md`.

## §5b — simulation characterises O1's RECALL exactly, and CANNOT touch its precision (2026-08-26)

**Question:** can the §5a simulation approach repair O1's caveats — chiefly the §4v completeness deficit,
whose one designed repair §4w refuted?

**Design.** The deficit hypothesis is about the EDGE RULE, not about reads, so it is testable directly on
constructed sequence pairs — no read simulation, no catalog run. 210 pairs, each a TRUE paralogue pair by
construction: both members carry the same conserved CORE (diverged 0.02–0.20), each padded with its own
private flank cut from a different real gorilla representative. **Sweeping the core fraction IS sweeping
"how complete is this model relative to what it shares".** Scored with the verbatim shipped invocation
(`-c -X --no-long-join -k 11 -w 5`) and the verbatim rule (identity ≥ 0.60, coverage ≥ 0.50 of the
shorter, forward only).

### ⭐⭐⭐ The failure surface is a clean step function, and divergence plays NO part

| core fraction | pairs | edge rate | median best coverage |
|---:|---:|---:|---:|
| 0.90 | 30 | 1.0000 | 0.8998 |
| 0.75 | 30 | 1.0000 | 0.7497 |
| 0.60 | 30 | 1.0000 | 0.5999 |
| **0.50** | 30 | **0.3667** | 0.4999 |
| 0.40 | 30 | **0.0000** | 0.3998 |
| 0.30 | 30 | 0.0000 | 0.2999 |
| 0.20 | 30 | 0.0000 | 0.1999 |

**Coverage tracks the shared fraction to three decimals.** ⟹ the clause is a direct readout of *how much
of the model is shared*, and the floor is a **hard cliff at exactly 50% shared** — deterministic, not
probabilistic. **Divergence 0.02–0.20 is irrelevant**, confirming with exact truth the standing result
that IDENTITY NEVER BINDS.
⚠ **The pooled 48.1% loss is NOT a real-world rate** — it is an artefact of sweeping core fraction
uniformly. The transferable claim is the **cliff position**, not the percentage.

⭐ Recall available by relaxing the floor, on exact truth: **0.40 → 0.6095 · 0.30 → 0.7810 ·
0.20 → 0.9286** (shipped 0.50 → 0.4810).

### ⛔ But the precision half is UNANSWERABLE by simulation — three constructions, three answers

A recall number alone is the registered trap "a denominator conditioned on the prediction". Scoring the
same rules against a negative class gave a **different winner every time the negatives were rebuilt**:

| negative class | what happened |
|---|---|
| 210 pairs sharing one fixed **1 kb** element | any bp-floor above 1 kb scores FPR 0.0000 — **the threshold was above the element by construction** |
| 96 pairs sharing an **identical 300–6000 bp** element | **nothing separates** (best Youden 0.1071) — but negatives share IDENTICAL sequence while positives share DIVERGED cores, so the negatives align better than the positives |
| 210 **real** cross-contig rep pairs, never same family | **7/210 clear identity ≥ 0.60**, best coverage max 0.1163 ⟹ **FPR 0.0000 for every rule is VACUOUS** |

⚠ The third is the same blind instrument as the human 150-window panel scoring 2/150 on zero qualifying
pairs. **CANDIDATE COUNT CHECKED BEFORE READING THE VERDICT** — which is the only reason it was caught.

### ⭐⭐ The limit, and it is structural

**Simulation supplies exact POSITIVES; it cannot supply the NEGATIVE class, because "a pair that looks
homologous but is not one family" is precisely what O1 lacks ground truth about in the first place.**
That is why §5a worked and this does not: O3's controls were REAL genomic loci, whereas O1's negatives
must be invented, and the invention decides the answer.

⟹ **ROUTE, not yet run:** pair the simulation's exact positives with a REAL negative panel that has real
truth — the human 150-window false-merge panel (measured 1.33%). Recall from simulation, precision from
the panel. Neither half can be got from the other.
⛔ **Do NOT lower the coverage floor on the strength of the recall column alone.** The 0.9286 at 0.20 is
real and the precision cost is UNMEASURED.

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

⚠ **SCOPE CORRECTED 2026-08-27 — these arms ran the REP-ONLY variant.** `pool_isoforms` defaults to
false, so `RUSTLE_SHARED_EXON=1` alone never exercised `shared_exon_edges_pooled`, in which each locus
contributes the exons of EVERY isoform collapsed into it. That variant targets the strongest stated
motivation for the whole route (46% of reps covering a known member are stubs while 53% of those loci
have a discarded gate-passing spliced chain), and the original §4w wording implied a coverage it did not
have. **Now run at matched identity 0.60** (`RUSTLE_SHARED_EXON_ISOFORMS=1`): 109,135 exons over 2,847
loci against 11,385 rep-only, **1,494 locus pairs linked against 966** — and **worse on every endpoint**:

| | families | NPIP loci | pure NPIP |
|---|---:|---:|---:|
| shipped exon-sum | 83 | **12/31** | **3** |
| SE60 rep-only | 59 | 10/31 | 1 |
| **SE60 + isoforms** | **134** | **10/31** | **0** |

⟹ pooling every isoform's exons **fragments harder** (134 families) and **destroys NPIP purity entirely
(3 → 0)**. **The verdict stands and now covers both variants.**
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

## §5c — the coverage floor CAN be lowered on edge evidence, and it FAILS end-to-end (2026-08-26)

**A fourth instance of the standing rule: a change to what an EDGE IS cannot be judged on edge metrics.**

### ⭐⭐ First, a real negative class — found by accident

Attempting to build a background of *verified-distinct* loci for an end-to-end benchmark, **417 of 420
random 30 kb gorilla regions were discarded for mutual homology.** That is not a broken experiment:
**50.5% of random genomic region pairs produce an alignment and 49.0% clear identity ≥ 0.60.**
⟹ **repeat-driven cross-homology between random genomic regions is UNIVERSAL**, and the coverage clause
is precisely the mechanism separating *"shares an Alu"* from *"is a paralogue"*. It also re-explains, from
a third direction, why **identity never binds** — identity is satisfied by every repeat in the genome.

Those 87,990 rejected pairs are the **real negative class §5b could not construct**: real sequence, real
repeat content, ≥ 50 kb clear of any known family.

| coverage floor | FPR (87,990 real genomic pairs) | TPR (§5b constructed positives) |
|---:|---:|---:|
| **0.50 (shipped)** | 0.026287 | 0.4810 |
| 0.40 | 0.026787 | 0.6095 |
| 0.30 | 0.027446 | 0.7810 |
| **0.20** | **0.028094** | **0.9286** |

**Between 0.20 and 0.50 there is a DEAD ZONE**: FPR moves +0.0018 while TPR moves +0.4476, a 250:1 trade.
The genomic pairs that pass are near-complete matches that pass at ANY floor, so lowering it admits
almost nothing new. On edge evidence alone this is the strongest repair candidate O1 has produced.

### ⛔ And it does not survive the partition

| coverage | families | copies | largest | NPIP loci | pure NPIP fams |
|---:|---:|---:|---:|---:|---:|
| **0.50 shipped** | 83 | 484 | 39 | **12/31** | **3** |
| 0.30 | 97 | 691 | 94 | **12/31** | **1** |
| 0.20 | 140 | 1020 | 104 | **12/31** | **1** |

⛔ **NPIP recovery is UNCHANGED at 12/31 in every arm**, while pure NPIP families fall **3 → 1**, the
largest family balloons **39 → 104**, and **87–94 baseline families fuse** (28 and 25 new families each
absorbing more than one). Cost with no benefit. *(Comparator: the scoped guard that was killed fused 112
pairs with the largest swallowing 43 — this is smaller, and still all cost.)*

### ⭐⭐⭐ Why — and it redirects the whole line

Node counts are **identical (2,847) in all three arms**: node construction is untouched, so only edges
moved. At NPIP:

| | loci |
|---|---:|
| NPIP loci total | 31 |
| **with a node at all** | **13** |
| reaching a family (every arm) | 12 |
| ⟹ **headroom for ANY edge-rule change** | **1** |
| ⟹ **headroom behind NODE CONSTRUCTION** | **18** |

**The edge rule has ONE locus of headroom at NPIP; node construction has EIGHTEEN.** §4v, §4w, §5b and
this section all optimised the edge rule — the wrong half of the pipeline for this family. A perfect edge
rule cannot connect a locus that never became a node.

⟹ **VERDICT: do NOT lower the floor.** The edge-level evidence is genuine and the dead zone is real, but
the binding constraint is upstream. ⭐ **The productive target is NODE CONSTRUCTION** — the 18 NPIP loci
that have reads (28/31 carry ≥ 3) yet produce no node, which is the same stub/one-rep-per-locus territory
as §4s/§4r, now with a number attached to what it is worth.

## §5d — ⚠ CORRECTION to §5c: node construction has 4 loci of headroom, not 18 (2026-08-26)

§5c redirected the line at node construction on the strength of "18 NPIP loci have reads but no node".
**Diagnosing those 18 individually shows the premise was inflated by a loose read criterion**, the same
error class that has cost this project several retractions.

### ⛔ The read counts were counting the neighbours

`samtools view <region>` returns every read overlapping by **≥ 1 bp**, so reads belonging to adjacent
genes that merely clip a locus boundary were counted as that locus's evidence. Requiring a read to lie
**≥ 50% inside** the locus:

| | ≥1 bp overlap | ≥50% inside |
|---|---:|---:|
| NPIP loci with ≥ 3 reads | **28/31** | **23/31** |
| total reads across the 31 loci | **5,544** | **419** |

**Only 7.6% of the reads survive.** The worst case is `NC_073242.2:29415572` — **4,784 reads by overlap,
19 inside**; the other 4,765 belong to the adjacent gene whose 21-exon / 4,187-read node sits 24 kb
upstream and touches this locus by 223 bp. ⚠ **"28/31 expressed" and "5,544 reads" are WITHDRAWN.**

### The corrected diagnosis of the 18 node-less loci

Per locus, reads required ≥ 50% inside, chains pooled over the junction-incidence component exactly as
`locus_support` does:

| cause | loci | is it a METHOD failure? |
|---|---:|---|
| **no reads inside at all** | **7/18** | ⛔ no — not transcribed in fibroblast |
| **pooled chain support < 3** | **7/18** | ⚠ partly — reads exist but fragment across chains that share no junction |
| **passes the read gate, still no node** | **4/18** | ✅ yes — genuine downstream failure |

### ⭐ Revised headroom at NPIP, replacing §5c's table

| layer | loci |
|---|---:|
| edge rule | **1** |
| node construction, downstream of the read gate | **4** |
| read fragmentation (chains that never reach 3 pooled) | **≤ 7** |
| **not transcribed — outside the method entirely** | **7** |

⟹ **§5c's "18 loci behind node construction" is WITHDRAWN; the correct figure is 4, or ≤ 11 if the
fragmentation cases are counted as addressable.** The redirection still points the right way — node
construction outranks the edge rule 4:1 rather than 18:1 — but the dominant limit at NPIP is **evidence,
not method**: 12/31 recovered against a ceiling of 23/31 expressed, and 7 loci carry no fibroblast
transcription at all.

⚠ This is consistent with the standing result that **O1 sees only EXPRESSED copies**, and with the
register entry that already refuted a "node construction is the binding constraint" claim on exactly this
kind of loose criterion (`~1,374 candidates` → ~282 after re-measurement, a 5× reduction).

## §5e — ORACLE ABLATION: the definition costs ~3%, node construction costs ~58% (2026-08-26)

**Design.** Substitute each pipeline stage with its perfect version and see what survives. This is the
first measurement that separates *"the definition is wrong"* from *"the definition never got the input"*.

| rung | nodes | edges | NPIP loci grouped |
|---|---|---|---:|
| 0 | oracle | oracle | 31/31 — ceiling by construction |
| **1** | **oracle (true locus sequence)** | **REAL shipped rule** | **30/31** |
| **1b** | **oracle EXON-SUMS** (real substrate) | **REAL shipped rule** | **5/5 available** |
| 3 | real | real | **12/31** — shipped |

**Rung 1**: the 31 true NPIP locus sequences taken straight from the genome (25.7 kb mean), run through
the verbatim shipped invocation and rule. **463/465 pairs align; 199 edges form; 30 of 31 loci land in
ONE component, 1 singleton.**

**Rung 1b** repeats it on the substrate the pipeline actually uses. A real gorilla NPIP transcript was
built from annotation — `NPIPB11`'s longest mRNA, **26 exons, 12,878 bp spliced** — and projected onto
each locus to yield genuine exon-sums. It projects cleanly to only **5/31** loci (a limit of
single-transcript projection across 0.82-identity paralogues, not of the rule), but **all 5 group, median
coverage 1.0000**.

⚠ Rung 1 uses GENOMIC spans, which retain introns and so favour a family that arose by segmental
duplication; rung 1b is the like-for-like substrate but is **n = 5**. Neither alone is conclusive; they
agree, and they agree with the shipped catalog's own behaviour (12 of the 13 loci that HAVE a node reach
a family = 0.92).

### ⭐⭐⭐ The decomposition

| where the loss happens | loci | share |
|---|---:|---:|
| **the definition (E_r + γ)** | **1** | **3%** |
| not transcribed in fibroblast | 7 | 23% |
| reads present but chains never reach pooled 3 | 7 | 23% |
| passes the read gate, still no node | 4 | 13% |
| **⟹ everything upstream of the definition** | **18** | **58%** |

⟹ **GIVEN A NODE, THE DEFINITION GROUPS IT — 30/31 on perfect input, 12/13 in the shipped run. The
open problem in O1 is NOT the definition; it is that 18 of 31 loci never become nodes.**

This retro-explains four refutations in a row (§4w, §4x, §5b, §5c): all four optimised the edge rule,
which owns **3%** of the loss. §5c's lowered coverage floor could not move NPIP recovery off 12/31 for
exactly this reason — a perfect edge rule cannot connect a locus that has no node.

⟹ **Where the remaining work is, in order of size:** the **7 fragmentation** cases (reads exist, chains
never agree — the shattered-locus problem `locus_support` pooling was built for and evidently does not
reach), then the **4** that clear the read gate and still produce nothing. The 7 untranscribed loci are
outside the method by construction — **O1 sees only expressed copies.**

## §5f — the canonical-junction rule is nearly INERT, and NPIP recall is at the EVIDENCE ceiling (2026-08-26)

**Proposal tested (user, 2026-08-26): do not enforce canonical junctions at all.** §5e had isolated 4
loci as the only genuinely addressable failures, and all four carried non-canonical junctions
(4/39, 8/38, 5/16, 2/11), which made the rule the obvious suspect — it is all-or-nothing, so one bad
junction kills a 26-exon model.

**Run end-to-end**, with `junction_majority` / `junction_nc_max_bp` first added to the params certificate:

| arm | reps | families | NPIP nodes | in a family | of the 4 rescued |
|---|---:|---:|---:|---:|---:|
| shipped (all-canonical) | 2,847 | 83 | 13/31 | 12/31 | — |
| `JUNCTION_MAJORITY=1` (< 10 kb) | 2,846 | 83 | 13/31 | 12/31 | **0** |
| + `NC_MAX_BP=1e9` (any length) | 2,845 | 83 | 13/31 | 12/31 | **0** |

⛔ **The rule is nearly INERT: it moves 2 reps out of 2,847 and changes NOTHING else** — identical family
count, identical NPIP recovery, zero of the four rescued. Skeleton count is **23,672 in every arm**,
because the canonical test runs in `build_spliced_seq`, downstream of a stage the four loci never reach.
⟹ **do not change it.** No benefit, and the code documents it as doing DOUBLE DUTY — it is also the
principal mis-chain filter, since a spuriously chained junction is usually non-canonical.

### Four hypotheses for the 4 loci, all refuted

| # | hypothesis | test | verdict |
|---|---|---|---|
| A | non-canonical junctions block them | both arms above | ⛔ 0/4 rescued |
| B | reads agree but junction coordinates jitter | re-chain with ±5 bp snapping | ⛔ chains@exact **==** chains@±5 bp |
| C | chains are compatible windows of one structure (read spanning introns 1–5 vs 3–8) | pairwise compatibility over the SHARED coordinate range | ⛔ genuinely **conflicting**: 42/91, 13/15, 11/15, 13/21 |
| D | the reads were mis-assigned here from sibling copies | MAPQ and `de` | ⛔ not supported — **1/45 MAPQ 0**, median MAPQ **60**, median `de` 0.0012–0.0064 (good alignments) |

### ⭐⭐ What is actually true

These loci carry **6–19 reads, and no two agree on a splice structure** — 14 reads yielding 14 distinct
chains, 6 yielding 6. `PASS1_MIN_READS = 2` is the binding constraint, and it is a *reasonable minimum*:
with no structure observed twice, there is nothing to assemble a model from.

⟹ **NPIP recall of 12/31 is at or near the EVIDENCE ceiling, not an algorithmic one.** Final accounting:

| | loci | fixable? |
|---|---:|---|
| not transcribed in fibroblast | 7 | ⛔ outside the method |
| 1–6 reads inside | 7 | ⛔ too few reads |
| 6–19 reads, no two agreeing on a structure | 4 | ⛔ nothing repeats |
| node exists, joins no family | 1 | the definition's 3% |
| **recovered** | **12** | |

⚠ **One architectural route remains and is NOT recommended**: a splice-GRAPH assembler would emit some
model from 14 conflicting reads where exact-chain matching emits none. That is what StringTie did, and it
was deliberately removed (`RETIREMENT_AND_MIGRATION.md`). It would trade "no model" for "a model of
unknown correctness" at exactly the loci with the least evidence.

⟹ **O1's node-construction line is CLOSED for NPIP.** §5e's "4 addressable loci" is revised to **~0**:
the 4 are evidence-limited too, just less starkly than the 7.

## §5g — the DNA substrate on RNA-found loci is WORSE, and it explains §5e (2026-08-27)

**Proposal tested (user): run the definition on DNA, add RNA afterwards.** ⚠ This is NOT the retracted
DNA∪RNA union — it is the substrate decision the joint-definition investigation itself *recommended*:
*"assembly supplies every base and therefore the WHOLE EDGE RELATION; RNA supplies WHICH LOCI EXIST."*
It is half-implemented as `--homology-genomic-span`, and had never been run end-to-end.

| arm | families | copies | NPIP nodes | **in a family** | **pure NPIP** | runtime |
|---|---:|---:|---:|---:|---:|---:|
| exon-sum (shipped) | 83 | 484 | 13/31 | **12/31** | **3** | 15:08 |
| **genomic span (DNA)** | 91 | 445 | 13/31 | **10/31** | **1** | **2:25:54** |

⛔ **Worse on every endpoint that matters, at 9.6× the wall clock** (peak RSS 16.4 GB; the substrate carries
~12× the sequence — 25.7 kb spans against ~2 kb exon-sums). 18 baseline families fuse.

### ⭐⭐ Why this contradicts §5e's 30/31 — and the answer is the whole point

§5e ran the same rule on genomic spans and grouped **30/31**. This arm gets 10/31. **The rule is identical;
the NODE SET is not.**

| | nodes | result |
|---|---|---|
| §5e rung 1 | the **31 TRUE locus spans**, complete and correctly bounded | 30/31 |
| §5g arm | the genomic spans of the **13 RNA-detected reps**, bounded wherever the reads fell | 10/31 |

⟹ **Taking the genomic span of a WRONG interval does not recover the gene** — it brings more sequence
with the same bad boundaries. The 24 kb-upstream node stays 24 kb upstream. **RNA still decides WHERE, so
the DNA substrate inherits every boundary error and cannot reach the 18 loci that have no node.**

⭐ **The register's recommendation is therefore REAL BUT BOUNDED**: DNA does supply a better relation
(§5e, 30/31 on correct nodes), and it is worthless while locus DISCOVERY stays with RNA. The binding
constraint is unchanged from §5e — 18/31 loci never become nodes.

### The version that WOULD work, and what it costs

A genuinely DNA-first pipeline must discover loci by **self-homology** (segmental-duplication detection),
not from RNA and not from known families — ⚠ `--from-genome` seeded from truth windows is **CIRCULAR** and
the register bans quoting it. Such a route would find all 31 NPIP loci; the homology search in this
session did exactly that from human orthologs.

⚠⚠ **BUT IT CHANGES WHAT O1 CLAIMS.** O1's objective is to *define a multi-copy gene family topologically
at the RNA level*. A pipeline that discovers loci by DNA self-homology and then annotates expression is
**segmental-duplication detection plus an expression filter** — a well-trodden method (SEDEF and
relatives), on which this project already holds cross-ape catalogs. The recall gain is real (31 vs 13
loci at NPIP) and the novelty claim moves. **That is a thesis-scope decision, not a technical one.**

⟹ **Do not flip `--homology-genomic-span`.** Its measured effect on the shipped node set is negative, and
the substrate objection the register already recorded — *"a span is mostly INTRON and paralogue introns do
not align"* — is joined by a second: **the span cannot fix a boundary the reads got wrong.**

## §5h — ⭐ REDEFINE THE NODE, NOT THE EDGE: footprint nodes recover +6 loci (2026-08-27)

**Question (user): redefine at the RNA level, or define at both levels?** The measurements answer it.
**"Both levels" is not needed** — §5g showed the DNA substrate inherits RNA's boundary errors, and the
DNA∪RNA union is long retracted. **The productive change is at the RNA level, and it is to the NODE, not
the EDGE**: §5e established the edge rule recovers 30/31 given correct nodes, so the edge rule is not what
is failing.

**What a node currently requires:** ≥3 pooled reads agreeing on ONE consistent intron chain. §5f showed
that at the failing loci **no two reads agree on a splice structure** (14 reads → 14 chains), which is
unreachable by any pooling or junction-tolerance rule.

**FOOTPRINT NODE (prototype):** the union of bases covered by ≥2 reads, exonic only (the reads' `N`
operations exclude introns), ≥300 bp, **with no chain-consistency requirement and no assembly**. Still
purely RNA — only read-covered bases.

| ladder, identical edge rule throughout | NPIP loci grouped |
|---|---:|
| oracle true locus spans (§5e) | 30/31 |
| **footprint nodes** | **18/31** |
| shipped assembled chains | 12/31 |
| genomic span of the shipped nodes (§5g) | 10/31 |

Footprints exist at **24/31** loci against 13/31 for assembled chains, median **6,175 bp**, and 18 group
into 2 components.

**PRECISION, measured — not the recall column alone.** The same construction applied at 60 non-NPIP
expressed loci as a control:

| edge class | n |
|---|---:|
| NPIP ↔ NPIP (wanted) | 30 |
| **NPIP ↔ control (false)** | **1** |
| control ↔ control | 1 |

⟹ **contamination 1/31 = 0.0323** for **+6 loci = +50% relative recall**.

### ⚠ What this is NOT yet

- **OFFLINE (T8).** A prototype over PAF, not the shipped binary. This project's own rule: offline
  re-derivation is a hypothesis generator, **never a test**.
- **Node-level.** Three prior definitional changes passed node metrics and failed end-to-end (§4x, §5b,
  §5c). **The verdict must come from the PARTITION** — families, copies, NPIP purity, cross-family fusion.
- The 60-locus control is small, and its 1 control↔control edge suggests it contains real paralogue pairs,
  so contamination is an **upper** bound on that panel and an **underestimate** of what a genome-wide
  control would contain.
- Footprints are ~3× longer than exon-sums, so the scale-free coverage denominator (§4v) will act on them
  differently at genome scale.

⟹ **Recommended next step, costed:** implement the footprint node behind a flag (default OFF, OFF arm
byte-identical, flag recorded in the params certificate), run the NPIP catalog end-to-end, and judge on
the partition. This is the first change in this session whose evidence points the right way on **both**
axes — every other candidate (§4w, §4x, §5c, §5f, §5g) failed on one or the other.

## §5i — footprint nodes IMPLEMENTED and measured end-to-end: +1 locus, ships OPT-IN (2026-08-27)

§5h's prototype is now a real flag (`RUSTLE_FOOTPRINT_NODES`, default OFF), run through the binary.
**The OFF arm reproduces the baseline BYTE-IDENTICALLY** (`copies.tsv` and `families.tsv`) and the knob is
in the params certificate.

| | OFF (shipped) | FOOTPRINT ON |
|---|---:|---:|
| families | 83 | 87 |
| copies | 484 | 491 |
| largest family | 39 | **38** |
| NPIP loci with a node | 13/31 | **15/31** |
| **NPIP loci in a family** | **12/31** | **13/31** |
| **pure NPIP families** | **3** | **3** |
| cross-family fusion | — | 8 families absorb >1 baseline (largest 5) |
| NPIP loci lost | — | **0** |

⭐ **The first change this session positive end-to-end on BOTH axes**: +1 locus recovered, 0 lost, purity
held 3/3, largest family did not grow. Fusion is modest against the comparator that killed the scoped
guard (112 pairs, largest swallowing 43).

### ⚠⚠ The prototype over-predicted 6×, and two of three bugs were invisible to it

§5h predicted **18/31**; the binary gives **13/31**. The prototype scored footprint SEQUENCES against each
other in isolation, while the real pipeline also runs locus collapse, representative selection and the
γ partition, where a footprint competes with existing reps and can be absorbed. **T8 vindicated
concretely.**

1. **The block set.** The pass runs INSIDE pass1, so `existing` still held the singleton chains the gate
   discards — and at exactly the target loci every read forms its own chain (§5f), so those doomed
   singletons reserved the region. Fixed to block only against `n_reads >= min_reads`.
2. ⚠ **The grouping threshold — the real one.** `max_gap = 100 kb` applied GENOME-WIDE chains every
   covered region on a contig together: on `NC_073244.2`, **73 islands, max span 4.6 Mb**. That is why the
   first two ON arms added 14 nodes and moved nothing. The prototype applied the same constant inside a
   30 kb window, where it could not misbehave.

**Re-chosen from data, not from the outcome:** 2 kb → 3,245 islands / median 2.4 kb · **5 kb → 1,443 /
7.1 kb** · 10 kb → 721 / 20 kb · 20 kb → 353 / 47 kb. 5 kb matches a typical gene span and independently
matches §5h's 6,175 bp footprint median, at the right scale for ~950 reps per contig. ⚠ **It was still
chosen on the substrate it is scored on** — a held-out confirmation is owed before any default change.

⟹ **SHIPS OPT-IN, default OFF.** The gain is real but small (+1 locus), it costs 8 family fusions, and its
one free parameter was tuned here. Same disposition as the read-strand fix (§4o/§4p): passed its criteria,
stayed opt-in because the ledger is two-sided.

## §5j — NPIP is ONE family fragmented into 5–6, not three real subfamilies (2026-08-27)

**Question (user): is NPIP three families, and should each family get its own graph?** Both halves have
clean answers, and the first is a defect rather than structure.

### The 3 "pure NPIP families" are FRAGMENTS of one family

| evidence | result |
|---|---|
| recovered NPIP loci, and how many families they occupy | 12 loci across **5 families** |
| does the split follow the real NPIPA / NPIPB subfamily boundary? | ⛔ **NO** — `GWFAM30` holds 5 NPIPB + 1 **NPIPA**, `GWFAM31` holds one of each |
| loci landing in MORE THAN ONE family | **3 of 12** — `21081589` (30+31), `35184501` (29+30), `21077139` (29+30) |
| the 17 NPIP reps graphed among THEMSELVES, shipped rule | **6 components**, sizes 6/4/3/2/1/1 |
| the 31 TRUE locus spans, same rule (§5e) | **1 component of 30** |

⟹ **The fragmentation is not subfamily structure** — it cuts across the NPIPA/NPIPB split, which is the
deep, real division in this family. And it is **not an artefact of the global graph**: graphing the NPIP
reps in isolation still yields 6 components. ⭐ **With correct nodes the same rule yields ONE component of
30 (§5e), so the fragmentation is caused by the NODES, not by the rule or by the partition** — the third
independent arrival at §5e's conclusion.

⚠ **3 of 12 loci are represented TWICE, in different families.** That is simultaneous over-splitting and
double-counting: two reps at one genomic locus that the collapse did not merge, then assigned to
different blocks. It also means a "copies" count over-states distinct loci at NPIP by ~25%.

### ⛔ One graph per family is CIRCULAR

The families are the OUTPUT of the graph; conditioning the graph on them assumes the answer. This project
has already registered that failure mode — *"a denominator conditioned on the prediction"* killed 7
claims, and *"prediction ⊆ its own truth is tautological"* killed 3 more. The pipeline is correctly
structured: ONE global all-vs-all → γ-partition → families.

⭐ **What IS legitimate, and already shipped:** per-family analysis AFTER the partition, as a certificate
rather than a definition — `families.tsv` carries `n_edges`, `density`, `lambda` and `cut_certified` per
family. Diagnosing a block is fine; letting the diagnosis choose the block is not.

⟹ **The actionable item is the double representation**, not the graph architecture: 3/12 NPIP loci carry
two reps that never collapsed. That is `family_detect`'s locus-collapse stage, upstream of `E_r` — the
same layer §5e/§5i already identified.

## §5k — ⛔ A FAMILY/SUBFAMILY HIERARCHY CANNOT RECONCILE THE FRAGMENTATION (2026-08-27)

**Proposal (user): define a SUBFAMILY level to reconcile §5j's fragmentation** — family = a coarse
grouping, subfamily = the current γ-blocks. NPIP would then read as *1 family, 6 subfamilies* rather than
5–6 unrelated families.

⚠ **First, γ can only SPLIT a connected component, never merge one.** So a hierarchy can only help if the
fragments are γ-splits of one component. They are not: §5j measured the 17 NPIP reps as **6 genuinely
DISCONNECTED components** at the shipped floor. Reconciling them therefore requires a MORE PERMISSIVE
edge floor, not a coarser partition of the same edges.

### The permissive tier does connect NPIP — by absorbing it into a hairball

17 NPIP reps in isolation: 6 components at floor 0.50, **1 component at 0.10**. But run genome-wide on
all 2,847 reps, the component NPIP lands in is:

| floor | NPIP components | their sizes | **non-NPIP reps in them** |
|---:|---:|---|---:|
| **0.50 (shipped)** | 5 | [384, 10, 2, 1, 1] | **381** |
| 0.30 | 3 | [856, 10, 2] | 851 |
| 0.20 | 2 | [1093, 10] | 1086 |
| **0.10** | **1** | **[1359]** | **1342** |
| 0.05 | 1 | [1440] | 1423 |

⛔ **THERE IS NO FLOOR AT WHICH NPIP IS ONE COMPONENT *AND* THAT COMPONENT IS NPIP-SPECIFIC.** The
single-component state always IS the blob — at 0.10 the component holds 1,359 reps of which **1,342 are
not NPIP**, and genome-wide that one component absorbs 1,359 of the 1,474 reps that are in any component
at all. The very permissiveness that joins NPIP's fragments joins everything else first.

⟹ **the hierarchy has no meaningful coarse level to offer.** Strict ⟹ fragmented; permissive ⟹ one blob.
Same shape as the coverage-repair impossibility (*"TP loss starts BEFORE FP rejection ⟹ no
constraint-satisfying threshold exists"*).

### ⭐⭐ And it corrects how §5j should be read

**At the SHIPPED floor, NPIP's largest component already contains 381 NON-NPIP reps.** So γ is not
over-splitting a clean NPIP component — **γ is carving NPIP out of a 384-rep hairball**, which is a much
harder problem and exactly the job γ exists to do (§4m: without γ one component holds 466 copies over 38
families). The 5–6 fragments are γ's best effort on a badly contaminated component, not a gratuitous split.

⟹ **Do not add a subfamily level.** It would relabel a node-construction defect as biological structure
while offering no level at which the coarse grouping means anything. ⭐ **The fragments are also NOT
subfamilies on the evidence** (§5j): they cut across the real NPIPA/NPIPB boundary, and 3 of 12 loci sit
in two of them at once — a locus cannot belong to two subfamilies.

⟹ **The finding stands where §5e/§5i/§5j put it: fix the NODES.** With correct nodes the same rule and the
same γ give ONE component of 30/31 (§5e) — no hierarchy required.


## §5l — graph-to-graph similarity: the implemented half is refuted, the other half targets the wrong axis (2026-08-27)

**Proposal (user): compare loci as GRAPHS rather than as one representative sequence each.** It splits
into two distinct ideas, and they have different fates.

**(a) Use every isoform, not one representative — IMPLEMENTED, now REFUTED.** This is
`shared_exon_edges_pooled`, and §4w had never actually run it (see the scope correction there). Run at
matched identity it is worse than both the shipped rule and the rep-only variant: 134 families,
**10/31 NPIP loci, 0 pure NPIP families**. ⟹ the isoform information is real (9.6× more exons) and using
it this way costs more than it buys.

**(b) Compare splice-graph STRUCTURE — unimplemented, and aimed at the wrong constraint.** Two loci would
be linked when their junctions CORRESPOND, not merely when their sequence aligns. ⚠ Paralogues sit at
different coordinates, so their junctions are not directly comparable: the comparison requires an
alignment first, then a check that junctions correspond in the alignment frame. That is strictly STRONGER
than the current sequence test, so it can only REMOVE edges ⟹ it buys **precision**.
⛔ **Precision is not the binding constraint**: §5c measured FPR moving **+0.0018** across the entire
coverage range while TPR moved **+0.4476**, and §5e attributes only **3%** of the loss to the edge rule.
A structural test would sharpen an axis that is already not costing anything.

⟹ **Best arm remains the footprint node (§5i): 13/31 loci, 3 pure, +1 over shipped.** Every route that
enriches the EDGE (shared exon, pooled isoforms, DNA substrate, lower floor, graph structure) has now been
measured, and the ones that moved the needle at all moved the NODE.

## §5m — ⭐⭐⭐ ONE SEED FINDS EVERY EXPRESSED MEMBER (2026-08-27)

**Question (user): starting from minimal annotation or seeds, can we find all family members — or at
least the expressed ones?** Measured on NPIP. **Yes, and the expressed answer is complete.**

**Design, and why it is not circular.** The seed is the **single annotated gorilla NPIP copy**,
`NPIPB11` at `NC_073242.2:31337556-31462161` (124.6 kb) — the only one the gorilla GFF names. The truth
set is the 31 loci located from **human** orthologs. Different sources, so recovery is a real measurement.
`minimap2 -x asm20 -c --eqx -N 200 -p 0.1` — ⚠ `-p` and `-N` recorded, copy counts are sensitive to both.

| | true NPIP loci found |
|---|---:|
| round 1 — the single seed | **16/31** |
| round 2 — iterate, using round 1's hits as seeds | **25/31** (+9) |
| 19 human seeds (this session's original search) | 31/31 |
| **the RNA pipeline's own node construction** | **13/31** |

### ⭐⭐⭐ Against the EXPRESSED subset, the seeded search is complete

| | of the 23 expressed loci |
|---|---:|
| **seeded search (1 seed + 1 iteration)** | **23/23 = 1.000** |
| RNA pipeline's own discovery | **13/23 = 0.565** |
| currently reach a family | 12 |

⭐ **RNA discovery is a strict SUBSET of seeded discovery**: seeding finds **12 loci RNA misses** and
misses **0** that RNA finds. **A single annotated copy plus one iteration recovers EVERY expressed member
of the family** — the pipeline's own discovery reaches 56.5% of them.

⟹ **combined with §5e** (with correct nodes the shipped rule groups 30/31), seeded discovery + the
EXISTING definition projects to ~23 recovered against the current **12**, from ONE annotation.
⟹ this sits inside the O1 reframing already accepted on 2026-08-12 (*"genome + MINIMAL ANNOTATION +
reads"*), so it is **not** the scope change that DNA-first discovery would be (§5g).

### ⭐⭐ CONVERGENCE AND TRANSITIVITY — the closure is well-behaved (2026-08-27)

⚠ §5m stopped at round 2 **arbitrarily**, not on convergence, so 25/31 was "what one iteration gave".
Run to a fixed point it is a **convergence result**, and precision does NOT degrade:

| round | loci hit | true NPIP | new NPIP | non-NPIP | precision |
|---:|---:|---:|---:|---:|---:|
| 1 | 17 | 16 | 16 | 1 | 0.941 |
| 2 | 26 | 25 | 9 | 1 | 0.962 |
| **3** | 26 | **25** | **0** | **1** | **0.962** |

**Converged at round 3.** ⭐ **The non-NPIP count is 1 at EVERY round** — transitive closure does NOT
chain through shared repeats here, which was the obvious failure mode and the reason precision was
tracked per round rather than only at the end.

⭐⭐⭐ **THE SIX LOCI IT NEVER REACHES ALL HAVE ZERO READS** (`28300719`, `29228842`, `21285416`,
`15670981`, `28701609`, `5715342`). Every one is unexpressed ⟹ **one seed + closure recovers
23/23 = 1.000 of the EXPRESSED members**, and everything it misses is invisible to RNA by construction.

⟹ **"SHOULD THERE BE MORE THAN ONE SEED?" — NO, not for an RNA-level objective.** One annotated copy
saturates the expressed set; additional seeds can only add loci with no RNA evidence. 19 human orthologs
reach 31/31, but the extra 6 over the closure's 25 are exactly the zero-read ones.
⟹ **"IS IT TRANSITIVE?" — YES, AND IT TERMINATES:** 3 rounds to a fixed point, precision flat at 0.962.


### ⚠ What this is not

- **n = 1 family.** A demonstration, not a rate. The 2026-08-12 minimal-annotation work measured
  Δ+0.1000 at **P = 0.0635 (not significant)** across 22 families, surviving only at m≥3 (5/22 vs 0/22),
  and left *"coordinate-only baseline never computed"* open.
- **NPIP is a recent hominoid expansion at high identity**, which is the easy case for homology seeding.
  An older family whose copies have diverged would be harder — the same age-dependence §5g raised.
- The seed is a **124.6 kb whole-gene span**; a shorter seed (CDS only) is untested and likely weaker.
- "23/23 expressed" is relative to the human-derived 31-locus truth set, and to §5d's corrected
  expression criterion (≥3 reads lying ≥50% INSIDE the locus, not ≥1 bp overlapping).

⟹ **NEXT, COSTED (~20 min):** feed the 25 seeded loci to the catalog as search windows and measure how
many reach a family. §5e predicts ~24; the current pipeline gets 12. ⚠ Use the SEEDED windows, never the
truth windows — `--from-genome` seeded from truth is circular and the register bans quoting it.

## §5n — seeded windows through the real binary: 26/31, and the fragmentation dissolves (2026-08-27)

**§5m's confirmatory run.** The 25 loci found by seeding from `NPIPB11` alone (plus 10 kb flanks, 26
windows) were fed to the catalog via `--from-genome`. ⚠ **SEEDED windows, never the truth windows** —
`--from-genome` seeded from truth is circular and the register bans quoting it.

⚠ **Scope: `--from-genome` is mutually exclusive with `--bam`, so this arm groups on DNA and builds NO RNA
nodes.** It measures seeded DISCOVERY + grouping, not the RNA-level objective.

| route — identical edge rule throughout | NPIP loci in a family |
|---|---:|
| §5e oracle, TRUE locus spans | 30/31 |
| **seeded windows (this)** | **26/31** |
| shipped, RNA discovery + RNA nodes | **12/31** |
| §5g, DNA spans of RNA-FOUND nodes | 10/31 |

### ⭐⭐ Three prior findings confirmed at once

1. **§5g's diagnosis is right.** Same DNA substrate, different node set: spans bounded by *where reads
   fell* give **10/31**; spans bounded by *homology to a real gene* give **26/31**. The substrate was
   never the problem — **the boundaries were.**
2. **§5m's projection lands.** Predicted ~23–24 from seeded discovery plus the existing grouping;
   measured **26**.
3. **§5k's fragmentation dissolves.** The catalog emits **2 families** — `GWFAM0` with 31 copies and
   `GWFAM1` with 3 — instead of the shipped 5–6 fragments. NPIP comes back as ONE family, without any
   change to the edge rule, γ, or a subfamily layer.

### Precision, and an honest asymmetry in the truth set

**5 of GWFAM0's 31 copies lie outside the 31-locus truth set (purity 0.839).** Their nearest true NPIP
locus is 63 kb – 1.87 Mb away, so they are not boundary spill. ⚠ **The truth set was built from HUMAN
orthologs, so a GORILLA-SPECIFIC NPIP copy is absent from it by construction and is scored here as an
error.** These 5 are therefore an upper bound on false positives and a plausible lower bound on
gorilla-specific copies — exactly O3's target class. They should not be quoted as 5 false positives.

⚠ **Cost:** 29 min and **22.6 GB peak** for 26 windows — the DNA substrate's ~12× sequence volume again
(§5g). Genome-wide this route would not fit on this machine.

⟹ **What is still unmeasured, and it is the actual objective:** seeded windows with **RNA nodes** built at
them. `--from-genome` cannot express that (no `--bam`), so it needs the footprint pass restricted to given
windows — a small feature, and the natural next step. **The DNA number above must not be quoted as an
RNA-level result.**

## §5o — ⚠⚠ A FOURTH CONSTRUCTION ERROR, SHIPPED BEHIND A FALSE COMMENT; and seeded windows REJECTED (2026-08-27)

Investigated by an 8-agent workflow with an adversarial verifier that re-measured every load-bearing
number from real binary runs and run logs rather than trusting the diagnoses.

### ⛔⛔ THE BUG: footprints were NOT exempt from the canonical-motif test

`denovo_assemble.rs` documented footprints as *"marked `footprint: true` so the canonical-motif test is
skipped … sequence taken as-is."* **That comment was FALSE.** Footprints routed to
`build_spliced_seq_with(majority = true, nc_max = u64::MAX)`, and `nc_max` lifts only the per-junction
SIZE veto — the majority path still ran `if !introns.is_empty() && plus == 0 && minus == 0 { return None }`,
demanding at least one canonical junction. A footprint's 50 bp-binned **coverage gaps** have no reason to
carry `GT-AG` / `GC-AG` / `AT-AC`.

**Measured:** across 22 candidate footprints, **8 of 248 gaps were canonical = 0.0323**, against a uniform
expectation of **6/256 = 0.0234** — 5.81 expected vs 8 observed, **Poisson p ≈ 0.24, INDISTINGUISHABLE
FROM CHANCE**. Every surviving multi-gap footprint passed on **exactly one accidental dinucleotide** out of
4–40 gaps. Genome-wide the guard killed **225 of 288 = 78%** of footprints.
⟹ **the feature's entire yield was a lottery ticket**, and §5i's +1 locus was drawn from it.

**FIXED**: a dedicated `build_footprint_seq` — no motif test, strand from the read vote as for an unspliced
model — with a regression test that pins BOTH halves (the footprint path accepts a motif-free gap; the
spliced path still rejects it even at `nc_max = MAX`).

**Effect of the fix**, scored on the objective (a locus counts when a family copy lies ≥50% INSIDE it):

| arm | families | nodes | **loci in a family** | pure |
|---|---:|---:|---:|---:|
| shipped | 83 | 13/31 | **12/31** | 3 |
| footprint, bug present (§5i) | 87 | 15/31 | **13/31** | 3 |
| **footprint, gate fixed** | 100 | **16/31** | **13/31** | 3 |

Post-gate transcripts +63 → **+210**; footprint survival roughly tripled. ⟹ **+1 NODE, 0 additional loci
on the objective** — the bottom of the verifier's pre-registered range (0 to +2, point estimate +1).
**Keep the fix** (the code now does what it claims, and node admission is no longer a coin flip), but it
is **NOT a gain on the objective**. ⚠ families 87 → 100 and copies 491 → 534: the extra blocks are outside
NPIP and their correctness is **UNMEASURED**.

### ⛔ SEEDED WINDOWS → RNA NODES: DO NOT IMPLEMENT

| arm | loci in a family | pure NPIP |
|---|---:|---:|
| shipped | 12/31 | 3 |
| **unseeded footprint (needs no seeds, no truth set, no geometry)** | **13/31** | **3** |
| seeded 26 windows | **12/31** | **2** |

**Seeding buys 0 loci over shipped and costs 1 locus and 1 pure family against the free feature.**
Two further findings from the verifier:
- ⚠ **The 22/31 ceiling I was aiming at is an ORACLE** — obtained by feeding the 31 TRUE locus spans as
  windows. Not reachable from any seed set buildable without the truth annotation, and tuning window
  extent against `ggo_loci.json` is this project's own banned move.
- ⚠ **The chain-capable ceiling is 11/31, BELOW the shipped 13/31** ⟹ a third of the shipped NPIP nodes
  are unspliced stubs, not chain agreement.

### ⚠⚠ TWO HEADLINE NUMBERS WERE COINCIDENCES

Reproduced from run logs alone, no modelling: *"reps unchanged at 2,847"* is **+6 / −6**, and *"NPIP
unchanged at 13/31"* is **+1 / −1**. Of the 6 reps the footprints destroyed, **4 were 100% inside a true
NPIP locus**, including a 3,637 bp model displaced by a 137 kb footprint only 20.46% inside — which then
merged 4 shipped reps, 2 of them with zero NPIP overlap. **An unchanged total is not evidence of
inertness; diff the sets.**

⟹ **§5m and §5n STAND** (one seed finds 23/23 expressed; seeded windows give 26/31 on DNA). What is
refuted is only that seeded boundaries help **RNA node construction**.

## §5p — the one-seed closure GENERALISES: 65 human families, 65/65 converge (2026-08-27)

⚠⚠ **HUMAN (CHM13 v2.0, Soto's 80-family panel). NEVER pool with the gorilla NPIP numbers.**

§5m's closure result was **n = 1 family**. Tested on Soto's human panel: 65 families with ≥ 3 members,
one seed each (the largest member), iterated to a fixed point. All families in ONE genome pass per round,
query names carrying family lineage. Same command as NPIP — `minimap2 -x asm20 -c --eqx -N 200 -p 0.1`,
hits ≥ 500 aligned bp. ⚠ `-p`/`-N` recorded.

**⚠ Coordinate validity checked first, not assumed**: `soto_replication/` holds CHM13 **v1.0**, so v2.0
compatibility was verified by requiring members of a family to be homologous TO EACH OTHER — 11, 81 and
30 cross-member alignments on three test families. They are.

### ⭐⭐ Stratified by the SD-vs-gene-family caveat (raised by the user, and measured)

| stratum | n | median recall | mean recall | converged | non-members grew |
|---|---:|---:|---:|---:|---:|
| gene-family-like (≤25% accession names) | 42 | **1.000** | 0.895 | **42/42** | 10/42 |
| mixed | 11 | **1.000** | 0.845 | **11/11** | 4/11 |
| **SD-like (≥50% accessions)** | 12 | **1.000** | **0.885** | **12/12** | 4/12 |

⭐ **65/65 families CONVERGE — no closure exploded.** Median recall is **1.000 in every stratum**.
⭐⭐ **THE SD CAVEAT DOES NOT MANIFEST: SD-like 0.885 vs gene-like 0.895.** The concern was well-founded a
priori — SDs are repeat-rich and closure should chain there — but measured, it does not. Non-member hits
move a median of **1 → 2** across all rounds, and only **18/65 = 0.277** of families show ANY growth.

⭐ **Recall is NOT size-dependent**: Pearson **r = −0.043** over 65 families; the families of 14, 16 and 17
members each reach **1.000**. So a large family is not inherently harder to close.

### ⚠ What may NOT be concluded

⛔ **Do NOT read this as "NPIP is harder than a typical Soto family"** even though 0.806 < 1.000 median.
**The truth sets are constructed differently**: NPIP's 31 loci come from HOMOLOGY to human orthologs — a
permissive, method-adjacent set that may include marginal loci — whereas Soto's members are a CURATED
panel. A permissive denominator lowers recall by construction. The two recalls are not commensurable.
⚠ Species differ (human vs gorilla), assemblies differ, and the seeds differ (largest member vs the one
annotated copy). **Only the WITHIN-panel comparisons above are safe.**

⟹ **The transferable claim, and it is now n = 66 rather than n = 1:** *a homology search seeded from one
family member and iterated to a fixed point converges, does not chain through repeats, and recovers the
median family completely — and this holds for segmental-duplication-like groupings as well as for gene
families.*

## §5q — EDGE WEIGHTS: the partition discards a 0.9491-AUC signal, and using it changes nothing (2026-08-27)

**Question (user): would other graph-theory concepts improve precision and recall?** Most are already
refuted in the register — connected components (238-member blobs at density 0.08), pure cliques
(over-splits tubulins/LILRs), min-cut/Fiedler/density as over-merge discriminators (AUC 0.35–0.52; **size
dominates and connectivity INVERTS once size is controlled**), λ as a membership clause (entailed at n=2,
and ~57% of families are pairs), and a pair-local edge veto replacing the partition (**the hairball's
connectivity is DISTRIBUTED — dropping 63% of its edges still leaves a 206-node component**).

⭐ **One thing was genuinely unexploited.** `denovo_pipeline.rs` flattened every edge to weight **1.0**
before the partition — while Louvain underneath is a *weighted*-modularity algorithm. Identity and
coverage were computed, used once for the threshold, then discarded.

**The discarded value carries strong signal** (unit = E_r edge, 3-contig NPIP catalog):

| edge class | n | median identity | median coverage |
|---|---:|---:|---:|
| NPIP ↔ NPIP (wanted) | 22 | **0.9878** | 0.9995 |
| NPIP ↔ non-NPIP (suspect) | 42 | **0.8281** | 0.9525 |
| non-NPIP ↔ non-NPIP | 1,782 | 0.7757 | 0.6642 |

**AUC(identity) = 0.9491**, AUC(coverage) = 0.6418.

### ⛔ And it does nothing end-to-end — the FOURTH such failure this session

`RUSTLE_ER_WEIGHTED_PARTITION=identity`, OFF arm byte-identical:

| arm | families | copies | NPIP loci | NPIP families | pure |
|---|---:|---:|---:|---:|---:|
| shipped (unweighted) | 83 | 484 | **12/31** | **5** | **3** |
| **weighted by identity** | 88 | 484 | **12/31** | **5** | **3** |

Identical on every objective endpoint. It splits 5 more blocks elsewhere at an unchanged copy count, and
does not touch NPIP's recovery or its 5-way fragmentation.

### ⭐⭐ THE MECHANISM, and it bounds the whole idea

`family_split.rs:511` — `induced_density` computes `edges.iter().filter(...).count()`. **It DISCARDS the
weight.** So the γ TEST (*whether* a block is split) is unweighted; only Louvain's `split_once` (*how* it
is split) sees weights. Since **γ is inert on 79% of families**, weighting reaches only a fraction of a
minority.

⚠ **And making `induced_density` weighted would make things WORSE, not better.** Identity weights are
< 1, so every weighted density is LOWER than its unweighted counterpart, more blocks fall under γ = 0.20,
and the partition splits MORE — while NPIP is already over-split 5 ways (§5j). Recovering the current
behaviour would need γ re-tuned to the weight distribution, which is fitting a threshold on the arms it
is scored against — the move this project forbids.

⟹ **VERDICT: edge weights are not the lever.** Retained as `RUSTLE_ER_WEIGHTED_PARTITION` (default OFF,
byte-identical, in the params certificate) as the record of the experiment.
⟹ **Answer to the original question: no graph-theoretic concept now looks promising.** The partition is
not where O1 loses copies — §5e put **~3%** of the loss on the edge rule and everything downstream of it,
and **~58%** on node construction. A better partition operates on the 3%.

## §5r — ⭐⭐⭐ TIER-2 ADMISSION: expression and assemblability are different questions (2026-08-27)

**Proposal (user): separate coverage so it says which loci are REALLY expressed, and let SIMILARITY admit
the rest.** This is the first proposal this session aimed at the layer that actually loses copies.

### Why NPIP is not fully recovered — the funnel, each row measured

| step | loss | remaining |
|---|---:|---:|
| 31 NPIP loci | | 31 |
| no expression (< 3 reads lying ≥ 50% INSIDE) | −8 | 23 |
| **expressed but NO NODE built** | **−10** | 13 |
| node built, joins no family | −1 | **12** |

⭐ **The dominant loss is loci that ARE expressed and cannot be assembled.** Reads 3–19 per locus, but
their reads never agree on a splice structure (worst case 14 reads → 14 distinct chains), so no chain
reaches `PASS1_MIN_READS = 2` and nothing pools to `GATE_MIN_READS = 3`.
⟹ the pipeline **conflates "is this locus expressed?" with "can I assemble a transcript model here?"**,
and failing the second currently kills the first.

### The measurement

For each of the 10 expressed-but-nodeless loci, its best match among the 13 loci that DO have a node:

| substrate | clears the SHIPPED `E_r` rule |
|---|---:|
| oracle (true locus spans) | **10/10** — identity 0.7890–0.9980, coverage 0.999–1.000 |
| **pipeline-knowable** (span of the READ CLUSTER vs a node-bearing locus's span) | **10/10** — identity 0.7890–0.9979, coverage 0.8802–1.000 |

⭐ **It works on sequence the pipeline can actually obtain** — the read cluster's own extent, not the
truth interval. ⟹ projected **12/31 → 22/31**, against the 23-locus expressed ceiling.

### ⭐⭐ NEGATIVE CONTROL — and this is what the previous four candidates failed

The same rule applied to 200 non-NPIP read clusters (≥ 3 reads):

| | rate |
|---|---:|
| true positives admitted | **10/10 = 1.0000** |
| **non-NPIP controls admitted** | **4/200 = 0.0200** |

153/200 controls produced a forward alignment to tier-1 and only 4 cleared the rule, so **the coverage
clause is what rejects them**, not the absence of homology.
⚠ **0.0200 is an UPPER bound**: "non-NPIP" is defined by the 31-locus homology-derived truth set, so a
genuine NPIP-related locus outside those 31 counts here as a false positive.

### ⚠ Why this is not the same pattern as §4x / §5b / §5c / §5q

Those four had strong edge-level evidence and failed end-to-end **because they changed the EDGE rule,
which §5e puts at ~3% of the loss** — a better edge cannot reach a locus that has no node. **This targets
the ~58%**: it admits loci that never became nodes. It is the first candidate whose mechanism matches the
measured failure.

⚠ **STILL A PROJECTION, NOT A RESULT.** No end-to-end run yet. Open questions the run must answer: how the
read-cluster span behaves when reads scatter (the `max_gap` failure of §5i), what the genome-wide
admission rate is against 2,847 reps rather than 200 controls, and whether admitted copies survive locus
collapse and the γ partition.

⟹ **The output must carry the EVIDENCE TIER** — tier 1 assembled its own model, tier 2 was admitted by
similarity plus weak expression. They are different claims and a single count would hide that.

## §5s — ⭐⭐⭐ TOP-UP: coverage is NOT the limit — LOCUS COLLAPSE is (2026-08-27)

**Proposal (user): top up every locus with simulated reads, so the definition can be tested with
expression removed as a variable.** Run, and it overturns the standing attribution.

### The construction, and its verification

40 simulated reads per locus from a 2%-error template, **added to** (not substituting for) the real
fibroblast BAM: 1,240 reads, merged total 1,010,636 primaries.
⚠ **Stratified by necessity**: only **11/31** loci accept a transcript projected from the one annotated
gorilla copy even at `-N 50 -p 0.05` — the copies are too diverged (median within-family identity 0.82)
for a sibling's splice structure to transfer. The other 20 were simulated from the genomic span.

**The delivery was verified, not assumed**: **93.71%** of simulated reads landed at their SOURCE locus,
5.56% were absorbed by a paralogue, 0.73% fell outside. **30/31 loci received ≥ 20 reads; 0 received
none.** So coverage really was delivered.

### ⛔ And recovery barely moves

| arm | families | copies | NPIP nodes | NPIP in a family |
|---|---:|---:|---:|---:|
| real reads only | 83 | 484 | 13/31 | **12/31** |
| **topped up (40 reads each)** | 76 | 445 | **17/31** | **15/31** |

Per arm (never pooled — the unspliced arm builds genomic-substrate nodes, worse per §5g):
**spliced 5/11 in a family · unspliced 10/20.**
⚠ **Two-sided: +8 loci gained, 5 LOST**, net +3. Families fall 83 → 76 and copies 484 → 445.

⟹ **With ~40 clean reads at essentially every locus, recovery is 15/31, not 31/31. COVERAGE IS NOT THE
LIMITING FACTOR.**

### ⭐⭐⭐ Where it actually fails — and it is not where §5e implied

Of the **14** loci still nodeless, **13 PASS the read gate** — 38–94 reads, pooled support up to **43**,
best single chain up to **23**. No drop reason is logged for any of them (the readthrough and mis-chain
filters name their coordinates; 160 and 14 drops respectively, none at these loci), and the transcripts
survive to the 18,088 that reach collapse.

⚠ **Checked that this is not my measurement rule again: 13/14 have NO overlapping rep AT ALL**, not a rep
failing a ≥50% overlap test. The nearest rep sits **0.7–141.6 kb away** and does not overlap.

⟹ **The transcripts are built and then ABSORBED by locus collapse.** `collapse_loci_span_aware` merges
the locus with a neighbour, `pick_locus_rep` keeps ONE transcript, and the winner's span is its own — so
an absorbed locus leaves **no trace at its own coordinates**.

### What this changes

⚠ **§5e attributed ~58% of the loss to "node construction". This localises it further: the transcripts
EXIST. The loss is in LOCUS COLLAPSE and representative selection, downstream of the gate and upstream of
the definition.** Same mechanism as §5j's double representation (3/12 loci carrying two reps in different
families) and §5o's absorbed footprints, now seen from the other side.

⚠ §4s closed *"which transcript becomes the representative"* (0/215 edge inheritance). **This is a
different question — whether a LOCUS survives collapse at all — and it is not covered by that closure.**

## §5t — ⚠ CORRECTION: the ceiling is NOT expression, and the chain disagreement is an ALIGNER artefact (2026-08-27)

**Two corrections, both prompted by the user asking why the topped-up run should be bounded by expression
at all. It should not be, and I said otherwise.**

### ⛔ RETRACTED: "31/31 is unreachable because 8 loci have no expression"

That applies to REAL-READ runs only. §5s topped **every** locus to ~40 reads, so in that run no locus is
unexpressed. Measured on the 8 loci that had no real expression:

| after top-up | count |
|---|---:|
| became a NODE | **6/8** |
| reached a FAMILY | **4/8** |
| still no node **despite pooled support 40** | **2/8** |

⟹ **expression was not their binding constraint** — six became nodes the moment they had reads.
⟹ **the "23-locus expressed ceiling" is a property of the real-read substrate, NOT of the method.** In the
topped-up run every locus has reads and recovery is still 15/31, so **the ceiling is set by LOCUS COLLAPSE
and by the fragmentation below — 31/31 is not excluded.**

### ⭐⭐⭐ The chain disagreement §5f attributed to biology is produced by the ALIGNER

Reads simulated from ONE identical template, 40 per locus, after `minimap2 -ax splice:hq`:

| arm | loci | median reads | **median distinct chains** | median best chain |
|---|---:|---:|---:|---:|
| **spliced** | 11 | 49 | **28** | 10 |
| unspliced | 20 | 47 | **3** | 1 |

⭐ **40 reads from a SINGLE template yield a median of 28 distinct intron chains.** There is no isoform
diversity in that input — the template is one sequence. The aligner is placing junctions inconsistently
across near-identical paralogues. Unspliced reads, having no junctions to disagree about, cluster cleanly
at a median of 3.

**How far apart are the placements?** Donor positions within 500 bp of a neighbour (n = 171):
median gap **64 bp**, q25 8, q75 187 — **only 0.187 of gaps are ≤ 5 bp**, 0.368 ≤ 20 bp, 0.456 ≤ 50 bp.

⚠⚠ **THIS QUALIFIES §5f.** That section tested ±5 bp snapping, found chains@exact == chains@±5bp, and
concluded the reads "genuinely disagree on splice structure — not coordinate jitter". **The ±5 bp window
was far too small: 81.3% of neighbouring donor placements differ by MORE than 5 bp.** The disagreement is
real at 5 bp and is nonetheless an ALIGNMENT artefact, visible here because the ground truth is a single
template. §5f's conclusion that the loci are evidence-limited stands for the real reads; its mechanism —
"genuinely different splice structures" — does **not**.

⚠ Scope: measured on SIMULATED reads. Whether real reads fragment by the same mechanism is **NOT
MEASURED**, though §5f's real-read pattern (14 reads → 14 chains) is the same shape.

⟹ **Two repairs are now on the table and they are independent**: locus collapse (§5s) and junction
placement tolerance at a realistic scale (tens to low hundreds of bp, not 5). Neither is the definition.

## §5u — ⛔⛔ RETRACTION of §5t's aligner claim: the simulated templates were not transcripts (2026-08-27)

**Prompted by the user asking whether the simulated reads were simply not varied enough, or the splicing
not convincing. The splicing was not convincing. It was not splicing at all.**

### The defect

§5t's spliced arm simulated reads from "transcripts" produced by projecting the one annotated gorilla NPIP
copy onto each locus with `minimap2 -x splice`. Inspected:

| | |
|---|---:|
| blocks under 20 bp | **6,429 / 6,919 = 0.9292** |
| `NPIP0` | 2,088 blocks, **median 2 bp** |
| `NPIP1` | 2,043 blocks, **median 1 bp** |
| projected junctions that are canonical | 81/94 = **0.8617** (real splice sites are ~0.99) |

⟹ **these are SHATTERED ALIGNMENT BLOCKS, not gene models.** Splice-aware alignment across
0.82-identity paralogues fragmented the transcript into thousands of micro-blocks, and the simulated
reads were drawn from that. A read carrying 2,000 one-base "exons" is not an mRNA, and no aligner
behaviour can be inferred from how it maps.

### ⛔ What is RETRACTED

- **"40 reads from ONE template yield a median of 28 distinct chains — the ALIGNER fragments them."**
  The input was garbage; the fragmentation is the simulation's, not the aligner's.
- **The donor-gap distribution** (median 64 bp, q25 8, q75 187, 0.187 of gaps ≤ 5 bp). Measuring scatter
  from an invalid template measures nothing.
- **§5t's qualification of §5f.** §5f tested ±5 bp snapping on REAL reads and concluded they genuinely
  disagree on splice structure. **That conclusion STANDS, unqualified.** The retraction is mine, not its.

### ✅ What SURVIVES, and why

The two arms had different template quality, and only one was invalid:

| arm | template | loci | nodes | in a family |
|---|---|---:|---:|---:|
| spliced | projected — **92.9% of blocks < 20 bp, INVALID** | 11 | 5 | 5 |
| **unspliced** | the locus **genomic span** — a real contiguous sequence, **VALID** | 20 | 12 | **10** |

⟹ **on the valid arm alone, 10/20 recovered — so 10 loci still fail with ~40 reads from a legitimate
template.** §5s's core finding is unaffected: **coverage is not the limiting factor and locus collapse
absorbs loci that pass the read gate.** Of the 14 nodeless loci in §5s, **8 are from the unspliced (valid)
arm**, including cases with pooled support up to 43.

⟹ **§5t's OTHER finding also survives** — that 6/8 previously-unexpressed loci became nodes once topped
up, so the "23 expressed loci" ceiling is a property of the real-read substrate rather than of the method.
That count does not depend on template quality.

### ⚠ The lesson

**A simulation is only as good as its generative model, and I never inspected mine.** The 28-chain result
was striking enough to feel like a discovery, and it was an artefact of an input I had not looked at. The
check that caught it — *are these blocks the size of exons?* — costs one command and should have run
before the first read was simulated. ⚠ **The user's second objection (no isoform diversity) also stands
and is unfixed**: one template per locus means the simulation cannot reproduce real isoform structure,
which makes it a LOWER bound on difficulty, not an upper one.

---

## §6a — AMY: the family definition stated on DNA alone (2026-08-27)

> ⛔⛔ **SUPERSEDED IN PART BY §6c (2026-08-28).** An adversarial audit broke this section's two
> headlines: the `-p`/seed-invariance result (⛔ `-p` is INERT at this tier — `-X` implies
> `--dual=no`) and the human comparator (⛔ CHM13 has 9–11 amylase units, not 5 — the LOC-name trap
> named in this very section). The identity/coverage reading is also withdrawn: the 0.60 identity
> floor sits below minimap2's emission floor and cannot fire. **Read §6c before quoting anything
> here.** What stands: gorilla = 3 AMY copies, the complete triangle, and MGAM/MGAM2 excluded.


**Why this section exists.** Every §5 experiment defined the family on the RNA representation — exon-sum
representatives, one rep per locus, micro-exon guards, read-derived strand. Each of those is a place the
definition can fail for reasons that have nothing to do with the definition. This section asks the prior
question: **stated on DNA, with genomic intervals as nodes, does the definition work?** Nothing here uses
reads, exon sums, rep picking, or a strand vote.

### The substrate

Gorilla annotation (`GGO_genomic.gff`, GCF_029281585.2, NHGRI_mGorGor1-v2.0_pri) carries **3 coding
amylase genes plus 1 pseudogene** in ~112 kb of `NC_073224.2`. ⚠ **AMY1 is present but under a `LOC`
name** — `LOC101133335` has product `alpha-amylase 1B`; searching on `gene=AMY` alone finds only AMY2A
and AMY2B and would have reported the family as size 2. Human CHM13 carries **5** (AMY2B, AMY2A, AMY1A,
AMY1B, AMY1C) in ~392 kb of `NC_060925.1`.

⚠ **Annotation gene spans are unusable as DNA nodes here.** `AMY2A`'s gene span is 50,111 bp and
*overlaps* `AMY2B`'s, and the pseudogene sits inside `AMY2B` — a readthrough transcript variant inflates
the span. The **CDS envelope** is clean and, notably, **uniform**:

| | gene | CDS envelope | span | CDS blocks |
|---|---|---|---|---|
| GGO | LOC101133335 (AMY1B) | 136,222,106–136,230,308 | 8,202 | 10 |
| GGO | AMY2A | 136,262,941–136,271,157 | 8,216 | 10 |
| GGO | AMY2B | 136,309,116–136,317,270 | 8,154 | 10 |
| HSA | AMY2B / AMY2A / AMY1A / AMY1B / AMY1C | — | 7,898–8,227 | 10 each |

⭐ **Every AMY unit, both species, is 10 CDS blocks and ~8.2 kb.** This matters for one specific reason:
O1's single named definitional hole is that the min-length coverage denominator is **scale-free**
(§4/def_failures — a ~1 kb repeat clears 0.50 of any node under 2 kb). With every node at 8.2 kb, **that
hole cannot bite on this family**, so anything that fails here fails for another reason.

### Stage 1 — discovery (seed → genome)

One seed, the AMY1B genomic interval, against the whole genome, `-x asm20 -p 0 -N 500`, keeping hits at
identity ≥ 0.60 and ≥ 0.50 of the **query** covered (the chromosome is not a node, so the denominator is
the query):

| copy | identity | cov_q | annotation |
|---|---|---|---|
| NC_073224.2:136,222,099–136,230,351 | 1.0000 | 1.000 | AMY1B (self) |
| NC_073224.2:136,262,934–136,271,200 | **0.9278** | 1.000 | AMY2A |
| NC_073224.2:136,309,109–136,317,310 | **0.9303** | 1.000 | AMY2B |

**3/3 annotated copies, at full query coverage.** Genome-wide the seed produces **only these 3 records**;
relaxing the coverage floor 0.50 → 0.20 → 0.05 adds **nothing**, and there are no off-cluster hits.

### Stage 2 — definition (copies all-vs-all)

The discovered intervals, all-vs-all under the shipped rule (`-c -X --no-long-join -k 11 -w 5`; identity
≥ 0.60 **and** coverage ≥ 0.50 **of the shorter**, forward-only), with MGAM and MGAM2 — which carry the
α-amylase domain but are a different family — as the **negative control**:

| pair | strand | identity | cov of shorter | edge |
|---|---|---|---|---|
| cpA–cpB | + | 0.9282 | **1.0000** | YES |
| cpA–cpC | + | 0.9308 | **1.0000** | YES |
| cpB–cpC | + | 0.9360 | **1.0000** | YES |
| MGAM–cpA | − | 0.7698 | 0.0342 | no |
| MGAM2–cpB | + | 0.6814 | 0.0635 | no |
| MGAM–MGAM2 | − | 0.7613 | 0.0127 | no |

⭐ **The AMY copies form a complete triangle (γ = 1.0); MGAM and MGAM2 are separate singleton
components.** ⭐⭐ **The controls pass the identity floor and are excluded by coverage alone** — 0.68–0.78
identity against a 0.60 floor, killed at 0.013–0.064 coverage. This is the **same signature** recorded on
RNA ("IDENTITY NEVER BINDS", 0.749–0.803 vs a 0.60 floor, 0/728): it is a property of the rule, not of
the RNA representation.

### ⭐⭐ `-p` decides whether the definition is even well-posed

`-p` is minimap2's secondary-to-primary score ratio. Running all three copies as seeds against the
cluster (real minimap2 runs at each `-p`, not an offline re-derivation):

| `-p` (with `-N 500`) | seed cpA | seed cpB | seed cpC | seed-invariant? |
|---|---|---|---|---|
| 0.0 | 3 | 3 | 3 | ✅ |
| 0.2 | 3 | 3 | 3 | ✅ |
| 0.5 | 1 | 2 | 2 | ❌ |
| **0.8 (minimap2 default)** | 1 | 1 | 1 | family invisible |

⭐⭐ **P1 seed-invariance holds at `-p ≤ 0.2` and breaks at `-p ≥ 0.5`.** At the default `-p 0.8` the
family collapses to 3 singletons. This reproduces the MAPKBP1 lesson (hapcnv: 1/1 at `-p 0.8`, 9/8 at
`-p 0.1`) on an independent family, and it converts `-p` from a tuned knob into a **derived** one:
choose the largest `-p` at which the partition is seed-invariant. That is a threshold fixed by a
property of the output, not by fitting to truth.

### What this establishes, and what it does not

⭐ **On DNA the definition is exact on this family**: 3/3 recovered, a perfect clique, a domain-sharing
negative control correctly excluded, and no unannotated copies. The RNA-side failures catalogued in §5
(2,088-block projections, stub representatives, placeholder strand, one-rep-per-locus) are therefore
**representation failures, not definition failures** — the same rule, given clean DNA nodes, is exact.

⚠ **Gorilla 3 vs human 5** is consistent with the human-specific AMY1 expansion and is *not* an O3
detection: no reference-absent copy was found, at any coverage floor down to 0.05.
⚠ **Scope**: one family, one genome, n=1 — this is an existence proof that the DNA statement is clean,
not a measurement of its precision or recall.
⚠ **Sensitivity caveat**: `asm20` addresses ~20% divergence, so copies between the 0.60 identity floor
and ~0.80 could in principle be missed; the `-k 11 -w 5` genome pass is the control for that and is
recorded separately.

### §6a addendum — the sensitivity control, and two traps it contained

The `-k 11 -w 5` genome pass intended as the divergence control **was killed at 24,943,544 KB RSS**
(24.9 GB of a 25 GB machine) after 5 minutes. ⚠ **OPERATIONAL: `k=11` against the 3.5 Gb genome is an
OOM, not a slow job** — the index reports *average occurrences 736.4*, and the shipped pipeline only ever
uses `-k 11 -w 5` on the small representative FASTA, never on a genome. Do not repeat it.

The control was instead run **independently of minimap2**, against the SEDEF segdup catalog
(`GGO_segdup_putative.bedpe`), asking: does any genome-wide duplication block overlap an AMY gene body?

⚠ **TRAP 1 — a single garbage record made every gene look 100% duplicated.** The catalog contains
`NC_073224.2:73,990,166-145,002,185 <-> NC_011120.1:0-16,411`: a **71,012,019 bp** span paired with the
**16,411 bp mitochondrial genome**. It spuriously overlaps every gene on the arm, and reported *100.0% of
the gene* for all three AMY copies. **A 71 Mb "duplication block" is not a duplication** — check the
block length before reading an overlap fraction.

⚠ **TRAP 2 — the one real off-cluster candidate was a repeat, and it looked convincing.** After removing
the mtDNA record, exactly one informative row survived: AMY2A paired over **7,840 bp = 95.4% of the gene**
with `NC_073227.2:75,086,075-75,096,529`. Directly tested, that interval is **99.4% soft-masked repeat,
carries no annotated gene, and produces ZERO alignment to AMY2A** at `-k 11 -w 5 -p 0 -N 50`. It is a
dispersed repeat whose coordinates happen to overlap AMY2A, not an AMY copy.

⭐ **Result: two independent methods agree.** The minimap2 seed search finds 3 copies and nothing
off-cluster down to 0.05 coverage; the SEDEF catalog's only off-cluster candidate is refuted on direct
test. **Gorilla has 3 AMY copies; human CHM13 has 5.** No reference-absent copy — this family presents
no O3 case.

---

## §6b — the DNA definition across four families: precision 1.0000, and orientation is the whole recall gap (2026-08-27)

> ⛔⛔ **SUPERSEDED IN PART BY §6c (2026-08-28).** ⛔ The title's "precision 1.0000" is not
> supportable (denominator = the prediction); ⛔ the scale-free narrowing was a ZERO-CANDIDATE
> result and the hole is real (36/150 on repeat-rich small nodes); ⛔ the recall gain is
> misattributed — it measures a node-extraction bug (reference vs gene orientation), not the
> forward-only guard; ⛔ "+33 recall points" is really +23.08 POINTS; ⛔ the answer key is wrong in
> 3 of 8 assertions. **Read §6c before quoting anything here.**


§6a was n=1. This section runs the same DNA-level definition over **30 annotated CDS-envelope nodes**
drawn from AMY, MAGEA, RFPL, GOLGA, HERC, SRGAP and NPIP, scored against an **answer key written before
the run** (`amy/PREREGISTER.md`).

⚠ **The prefix trap was avoided by design.** `GOLGA*` and `HERC*` are not families — GOLGA6L7/6L10 are
the multi-copy pair while GOLGA1-5/7 are distinct genes, and HERC1/2 (giant) are unrelated to HERC3-6
(small). ZNF already sprang this trap once (231 "members" that are a domain superfamily). The key
therefore separates **recent near-identical duplicates** (expected to form families) from **ancient
domain-sharing paralogs** (expected to stay apart, as MGAM/MGAM2 correctly did in §6a).

### The pre-registered failure mode did NOT occur

Node size spans **951 bp (MAGEA1) to 302,838 bp (SRGAP1) — a 319× range**, so I pre-registered that the
scale-free coverage denominator would fire: ~476 bp of the 951 bp MAGEA1 landing anywhere inside the
302 kb SRGAP1 gives coverage-of-the-shorter ≥ 0.50 and a spurious edge.

⭐ **Zero such edges appeared.** 226 pairs had *some* alignment; 8 became edges; **0 were cross-family.**
⟹ **The scale-free hole is not triggered by size disparity alone — it needs a shared REPEAT.** Unrelated
gene bodies do not contain ~476 bp of contiguous ≥60% homology. This narrows the hole recorded in
[[project_o1_definitional_failures]]: it is a *repeat* problem, not a *length-ratio* problem, and a
length-ratio-based guard would therefore be aimed at the wrong variable.

### Which clause fails, over every within-family pair

| outcome | count |
|---|---|
| edge formed | 8 |
| fails COVERAGE | 1 |
| fails IDENTITY | **0** |
| **NO FORWARD ALIGNMENT** | **8** |

⭐⭐ **Identity fails zero times** — the fourth independent confirmation of "identity never binds"
(0/728 in false_omission; 0.749–0.803 vs the 0.60 floor in def_failures; 0.68–0.78 for MGAM/MGAM2 in
§6a). ⭐⭐⭐ **But 8 of the 10 misses have no forward alignment at all, and every one of them is an
INVERTED DUPLICATION**: MAGEA1/4/9 are (+) while MAGEA10/12 are (−); RFPL1/3 are (+) while RFPL2 is (−);
GOLGA6L7 is (−) while GOLGA6L10 is (+). On the reverse strand **7/8 are recoverable**, including
**RFPL2–RFPL3 at coverage 0.999, identity 0.959** — a near-perfect paralogue pair killed purely by
orientation. The 8th (GOLGA6L7–GOLGA6L10) fails coverage at 0.294 even on the reverse strand, so it is a
genuine coverage miss, not an orientation one.

| arm | edges | cross-family pairs | members placed |
|---|---|---|---|
| forward-only | 8 | 0 | 9/13 = **0.6923** |
| both strands | 16 | 0 | 12/13 = **0.9231** |

per family, both-strands: AMY **1.0000**, MAGEA **1.0000** (was 0.6000), RFPL **1.0000** (was 0.6667),
GOLGA6L 0.5000 (unchanged — the real coverage miss).

### ⚠ This validates a decision already in the code; it is NOT a bug report

`RefineParams::forward_only_active()` is `require_forward_alignment && substrate ==
Substrate::TranscriptOriented`, documented as *"meaningful only between transcript-oriented reps, so it
is inert on a reference-oriented substrate no matter what `require_forward_alignment` says"*, and
`denovo_pipeline.rs:4985` sets `dna_params.substrate = Substrate::ReferenceOriented` for `--from-genome`.
⟹ **The shipped DNA path already disables the guard.** What is new here is the *measurement*: the design
choice is worth **+33 recall points (0.6923 → 0.9231) at zero cross-family cost** on a pre-registered
key. ⚠ **This says nothing about the RNA default** — on transcript-oriented reps an antisense overlap is
a different biological object, and §4o's frozen-antisense panel (0/58 admitted, antisense-family rate
0.0048) is the instrument that governs that decision, not this one.

### §6b precision control — a size-matched negative panel, and why the both-strands arm passed a HARDER test

The §6b precision figure rested on only 17 singleton nodes as possible false partners, which cannot
support the words "precision 1.0000". Replaced with a proper control: **600 unrelated gorilla genes drawn
size-matched** to the family nodes (median 16,663 bp vs the families' 17,140 bp; range 763–347,015 vs
945–302,838) — matched on the SIZE DISTRIBUTION, not on count, because a count-matched null proves
nothing here.

⚠ **Scoring rule fixed in advance**: a false merge is counted **only** as a family↔negative edge.
Negative↔negative merges are NOT scored, because two randomly drawn genes may be genuine paralogues and
calling such a merge an error would assume the answer.

30 family nodes as query against the 600 negatives as target (18.8 s, peak RSS 1.68 GB — the earlier
630×630 self-alignment was abandoned at 27 min for doing ~11× the work to answer the same question):

| arm | eligible pairs (any alignment) | false merges | rate | Wilson 95% |
|---|---|---|---|---|
| forward-only | 205 | **0** | 0.000000 | [0, **0.0184**] |
| both strands | **389** | **0** | 0.000000 | [0, **0.0098**] |

**389 of 18,000 possible pairs had any alignment at all (0.0216)** — the candidate count is real, so the
zero is not the vacuous kind.

⭐⭐ **The load-bearing observation: allowing both strands nearly DOUBLES the eligible pair count
(205 → 389).** The both-strands arm was therefore exposed to ~1.9× as many opportunities to false-merge,
produced **zero**, and ends with a **tighter** upper bound (0.0098 vs 0.0184). So the +33 recall points
of §6b are not bought by loosening: the harder-tested arm is also the better-bounded one.

⚠ **Do NOT compare this to the RNA false-merge rate** (human 150-window panel, 2/150 = 1.33%
[0.37, 4.73]). Different substrate (DNA intervals vs RNA exon-sum reps), different species panel, and
different construction. The numbers must not be pooled or ranked against each other.
⚠ **n is still 4 families / 13 members** on the recall side; only the precision side is now well-powered.

---

## §6c — ADVERSARIAL AUDIT OF §6a/§6b: seven of eight claims broken or narrowed (2026-08-28)

A 25-agent audit (2 adversarial lenses per claim, then adjudication) attacked every load-bearing claim
in §6a and §6b. **Seven of eight were REFUTED or OVERSTATED.** The arithmetic reproduced everywhere —
every count in §6a/§6b is correct — but the *inferences* drawn from the counts mostly do not survive.
This section is the correction of record; §6a and §6b must be read through it.

### ⛔⛔⛔ 1. THE IDENTITY CLAUSE IS INERT — AND THIS REACHES BEYOND DNA

⛔ **RETRACT "identity never binds" as a finding about the definition, on every substrate.** The E_r
identity floor is 0.60. minimap2 cannot emit an alignment below its own scoring floor `B/(A+B)`:
**default `A=2,B=4` ⟹ 0.667; asm20 `A=1,B=4` ⟹ 0.800.**

⭐ **The observed minima track the preset, with the biology held constant**: the three default-scored
PAFs bottom out at **0.6500 / 0.6313 / 0.6577**, while `amy_vs_genome.asm20.paf` bottoms out at
**0.8118**. Over 71,950 non-self records, **0 fall below 0.60**. A log-linear fit to the [0.66, 0.76)
bins predicts 34.2 records in [0.60, 0.65) against **0** observed (Poisson P = 1.4e-15) ⟹ the
distribution is **TRUNCATED, not exhausted**.

⭐ **Verified independently here** (synthetic control, seed 7): a 6 kb sequence against copies at 5–45%
divergence emits an alignment **only for the 5% copy (identity 0.947)**; 15/25/30/35/40/45% give **no
alignment at all**, under default *and* asm20 scoring, with **0 records below 0.60 identity**.

⟹ ⭐⭐ **E_r is not a two-clause rule. It is COVERAGE (+ orientation); the identity clause has never been
able to fire under default scoring.** Deleting the 0.60 floor leaves the 30-node partition
byte-identical. ⚠⚠ **This reinterprets the RNA results too** — "IDENTITY NEVER FAILS 0/728" and
"0.749–0.803 vs a 0.60 floor" are measurements of minimap2's emission floor, **not** findings about
paralogue biology. Four "independent confirmations" were four measurements of the same constant.

⛔ Also misattributed in §6b: **orientation blocks 9/9 and coverage does ZERO independent
discriminative work in the forward arm.** The pair booked as "fails COVERAGE" (RFPL1–RFPL2) is itself an
inverted duplication — its only forward record is a **292 bp fragment at coverage 0.0386**, recovering at
**0.7682 / 0.9537 on reverse**. Coverage's only genuine bite is **1/16**, visible solely in the
both-strands arm.

### ⛔⛔ 2. §6b's SCALE-FREE NARROWING WAS A ZERO-CANDIDATE RESULT — AND THE HOLE IS REAL

⛔ **RETRACT "the scale-free hole needs a shared REPEAT, not a length ratio."** It rested on an **empty
cell**: **0 of 54 small (<2 kb) × giant (>100 kb) pairs produced any alignment record**, and no record
among 89,447 joins a ≤2 kb node to a ≥100 kb node. **The quoted 319× ratio appears in NO pair the rule
ever evaluated** (max ratio actually scored: **49.08×**; and the smallest node is MAGEA12 at 945 bp, not
MAGEA1 at 951). This is precisely the project's own *check the candidate count before reading a verdict*
trap — committed one section after restating it.

⛔ **The stated cause was also false.** Unrelated gene bodies DO carry long homology: **133/209 = 63.6%**
of aligned cross-family pairs span ≥476 bp on the shorter node, with contiguous match runs to **1,161 bp**
(span 5,013 bp, HERC1×HERC2, identity 0.9469). They fail because the shorter node's **median length among
aligned pairs is 18,977 bp**, where a 5 kb block is coverage 0.03.

⭐⭐ **When the missing cell was supplied, the hole fired immediately**: 150 real repeat-rich 1,000 bp
gorilla windows from a chromosome carrying no family node give **36/150 = 0.2400 spurious E_r edges**
(coverage to 1.0000, identity 0.714–0.965) versus **0/150 repeat-poor windows** (Fisher **p = 2.7e-12**).
⭐ **What actually protects the 30-node panel is ANNOTATION, not the rule**: node length and repeat
content are confounded (Spearman ρ = +0.735), and **0/120 annotated genes <2 kb carry a ≥476 bp masked
run, versus 261/286 = 0.9126 of genes ≥20 kb**. ⟹ the hole is real, ~24% on repeat-rich small nodes, and
annotated CDS envelopes are structurally shielded from it.

### ⛔⛔ 3. THE RECALL GAIN IS MISATTRIBUTED — IT IS MY EXTRACTION BUG, NOT THE GUARD

⛔ **`fam_nodes.fa` was extracted in REFERENCE orientation even though `fam_nodes.bed` carries the gene
strand (16 of 30 nodes are '−').** All 9 within-family pairs failing the forward arm are
opposite-annotated-strand; all 8 passing are same-strand — **a perfect 17/17 confound.** Re-extracting
the nodes in **GENE orientation** yields **16 edges / 12-13 placed / no-forward 8 → 0 WITH THE GUARD
STILL ON** (confirmed two ways, including an alignment-free algebraic relabelling of the original PAF).

⟹ ⭐⭐ **The correct fix is to store DNA nodes 5'→3' in gene orientation, NOT to drop the forward-only
guard.** §6b's framing ("the guard is an RNA rule with no justification on DNA") is withdrawn as the
*explanation*; the +23 points are real but they measure my node construction, not the guard.
⛔ **Arithmetic correction: the gain is +23.08 percentage POINTS** (9/13 = 0.6923 → 12/13 = 0.9231).
"+33 recall points" conflated the **relative** gain (+33.3%) with points — an inflation of 1.44×.
⚠ **Power**: the +3 members come from 2 of 4 families and 3 inverted genes, so **n_eff ≈ 2**; a cluster
bootstrap over families gives 95% CI **[0.0, 37.5] points, P(Δ=0) = 0.059**. ⚠ The metric is monotone
non-decreasing under edge addition, so the SIGN was guaranteed a priori.

### ⛔⛔ 4. `-p` AND `-N` ARE INERT AT THE E_r TIER — §6a's HEADLINE CONTRADICTS THE TREE

⛔ **RETRACT "P1 seed-invariance makes `-p` a derived parameter."** Under `-c -X --no-long-join -k 11
-w 5`, **`-X` implies `--dual=no` and `-p`/`-N` do nothing**: the shipped `fam_allvall.paf` is
**byte-identical (89,447/89,447 records, md5 4b5920241b41e7e4a0022aee33387c07)** to a rerun at the
minimap2 **default `-p 0.8`**. ⚠⚠ **`denovo_pipeline.rs:4194` already documents this** ("`-N`/`-p` are
INERT at this tier; `-X` is the operative difference, because it implies `--dual=no`") — the §6a headline
contradicted an in-tree finding. §6b's 8 forward / 16 both-strand / 0 cross-family result was **already a
default-`-p` result.**

⚠ The `-p` sweep did measure something real — the **stage-1 discovery pass** (no `-X`, target contains the
seed's own locus, so a perfect self-hit sets the reference) — but there the true breakpoint is
**0.471627 = 3898/8265, not 0.5**, the {0, 0.2, 0.5, 0.8} grid never probed 0.3–0.47, and count-invariance
is **NON-MONOTONE** (holding on [0, 0.4716) and again on [0.5248, 1.0]) ⟹ "largest `-p` preserving
invariance" is not well-defined.

### ⛔⛔ 5. HUMAN AMY: THE LOC-NAME TRAP, COMMITTED ONE PARAGRAPH AFTER NAMING IT

⛔ **"Human CHM13 has 5 AMY genes" is FALSE.** `NC_060925.1:103.3–103.9 Mb` carries **12 CDS-bearing
loci — 11 protein-coding α-amylase plus the AMYP1 pseudogene**; applying §6a's OWN unit criterion (10 CDS
blocks, 7.0–9.5 kb envelope) gives **9 units**. **"5" is exactly what `gene=AMY[0-9]` returns** — the
identical LOC-name trap §6a documents for gorilla *one paragraph earlier*. ✅ **GGO = 3 stands** (same
script). ⟹ the gorilla-vs-human contrast is **3 vs 9–11**, not 3 vs 5.

⛔ **"Two independent methods agree" is withdrawn.** `GGO_segdup_putative.bedpe` carries
**identity=nan in 399,405/399,405 rows**, score 0 throughout, and **FAILS its positive control 0/3**: it
contains **no row pairing any two of the three real AMY duplications**. It is not a working instrument
here. ⚠ Both arms also read `GGO.fasta`, so they were never independent, and **neither can address
"reference-absent" by construction.** The off-cluster candidate refutation stands, but on the 600-gene
panel **39/114 = 0.3421 of length-matched genes have ≥1 such off-locus repeat partner** — that outcome is
the *negative control's* expected behaviour, not evidence.

### ⛔⛔ 6. THE ANSWER KEY ITSELF IS WRONG IN 3 OF 8 ASSERTIONS — SAME ROOT CAUSE EVERY TIME

`PREREGISTER.md` was built by a **prefix grep over NAMED genes**, which silently drops LOC-named
paralogues. Both auditors independently found the same three errors:

- ⛔ **GOLGA6L7 and GOLGA6L10 are NOT a family — they are different subfamilies.** 6L7 forms a clique with
  LOC115931294 / LOC134757625 / LOC101137218 (all "golgin subfamily A member 6-like protein 7",
  id 0.9668–0.9675, **cov 1.0000**); 6L10 with LOC115930840 / LOC101138066 ("6-like protein 9") and
  LOC129523543 (id 0.9555–0.9922). Between the two clusters: **0/28 edges, max cov 0.3341.**
  ⟹ ⭐ **§6b's ONE declared recall miss (GOLGA6L 0.5000) is a KEY ERROR, not a limitation of the rule.**
- ⛔ **RFPL4A is not a singleton** — it heads an **8-locus tandem array** on NC_073244.2 (~9,142 bp
  periodicity, 8 distinct GeneIDs, all "ret finger protein-like 4A") at **id 0.9738–0.9825, cov 1.0000** —
  *tighter than the declared RFPL family*.
- ⛔ **NPIPB11 is not a singleton, and "only NPIP member annotated" is factually false** — 5 NPIPB loci
  exist; NPIPB11 edges LOC101141990 (id 0.9682, cov 0.8007) and LOC101134557 (id 0.9518, cov 0.8344).

✅ **Survives**: AMY = {LOC101133335, AMY2A, AMY2B}; MAGEA = {1,4,9,10,12} as one family (id 0.7743–0.8723,
cov 0.6339–1.0000, MAGEA10 the weakest link but clearing the floor); RFPL{1,2,3} with RFPL4A excluded
(563 bp 3'-domain block, cov 0.4086); and the ancient sets correctly non-merging — **SRGAP1/2/3 max cov
0.0059, HERC1–6 max 0.0414, GOLGA1–5/7 max 0.0725.** ⭐ The "fake pass" worry about SRGAP is **unfounded**:
SRGAP1/2/3 have no annotated paralogue in the GFF.

### ⚠ 7. PRECISION: THE ZERO IS REAL, THE WORD "1.0000" IS NOT — AND 0.50 IS UNTESTED

The zero is genuine and not luck: within-family base rate among the 226 aligned pairs is **17/226 =
0.0752 (13.3× lift)**, and P(a random 16-of-226 draw is all within-family) = **1.3e-23**.
⛔ But **"precision 1.0000" has the prediction as its denominator** (n=8, n=16 ⟹ Wilson-95 *lower* bounds
only **0.6756 / 0.8064**). ✅ **Correct statement: 0 cross-family edges among 209 eligible cross-family
pairs, Wilson 95% upper 0.0180** (0/189, upper 0.0199 forward-only), plus **0/389 family×negative,
upper 0.0098.**
⚠⚠ ⭐ **The decision gap is EMPTY: highest cross-family coverage 0.4086 vs lowest true-edge coverage
0.6339** ⟹ **any coverage floor in [0.41, 0.64] returns the identical partition, so the shipped 0.50 is
UNTESTED on this panel.** ⚠ Also **every predicted component is already a clique**, so transitive closure
added 0 pairs and single-linkage chaining — the usual false-merge source — was never exercised.
⚠ The negative panel is **blind at the nodes supplying 2 of the 3 gained members**: all 5 MAGEA nodes have
**zero** cross-family alignment records.

---

## §6d — the node set rebuilt from SEQUENCE, not names (2026-08-28)

> ⛔⛔ **LARGELY SUPERSEDED BY §6f (2026-08-28).** A second audit refuted this section's three
> headlines as TRUE BY CONSTRUCTION: the "guard costs nothing" result is a 0-of-0 (0/202 pairs could
> even supply a reverse-strand candidate) and its recommendation is backwards — the shipped
> `--from-genome` path already force-disables the guard; the "7 clean components" are forced by the
> stage-1 gate (58/58 within-seed pairs entailed); and the NPIP "SD blocks" are in fact 6 UNANNOTATED
> GENE COPIES (85–99% of NPIPB11's CDS union). The hyphen mechanism came from a truncated index.
> **Read §6f before quoting anything here.**

§6c showed every name-based error in §6a/§6b came from one cause: a **prefix grep over gene symbols**,
which silently drops LOC-named paralogues. This section rebuilds the node set name-blind and re-runs.

**Construction.** One seed per family, **CDS envelope extracted in GENE orientation** (minus-strand genes
revcomped — the fix for §6c's 17/17 confound). Seeds vs whole genome, `-x asm20 -p 0 -N 500`. Every hit
at identity ≥ 0.60 and ≥ 0.50 **query** coverage becomes a copy interval, **extracted in the orientation
of the hit**. Gene names are attached **afterwards, as a readout only** — they never define membership.

### ⭐⭐ BOTH name-based routes fail, and they fail the SAME way

| product keyword | genes | LOC-named | invisible to a symbol grep |
|---|---|---|---|
| ret finger protein-like | 9 | 9 | **100%** |
| golgin subfamily A member 6 | 14 | 14 | **100%** |
| alpha-amylase | 1 | 1 | 100% |
| MAGE family member A | 5 | 0 | 0% |

⚠⚠ **The mirror trap: a PRODUCT-keyword census is no better.** NCBI phrases the same family differently
for symbol-bearing vs LOC-named genes, and **a single hyphen splits it**:

| symbol gene | its product | LOC gene | its product |
|---|---|---|---|
| RFPL1 | "ret finger protein **like** 1" | LOC101151087 | "ret finger protein**-like** 4A" |
| GOLGA6L7 | "golgin **A6 family like 7**" | LOC101137218 | "golgin **subfamily A member 6-like protein** 7" |
| NPIPB11 | "nuclear pore complex **interacting**" | LOC101141990 | "nuclear pore complex**-interacting**" |

⟹ symbol grep finds one half of each family, product keyword finds the other, and **both split along the
same symbol/LOC boundary** — which is why every §6b key error took the form *"X is a singleton"*.
⭐ Some members are invisible to **both**: `LOC101137871` (MAGEA component) and `LOC101149946 /
LOC101141018 / LOC101140740` (GOLGA6L component) carry **no product description at all**.
⟹ **Only sequence is a reliable membership test.** MAGEA is the one family the old key got right, and it
is the one family that is 0% LOC-named.

### Result: 7 components, one per seed, zero cross-seed contamination

| seed | symbol grep gave | sequence discovery gives |
|---|---|---|
| GOLGA6L | 2 | **7** |
| NPIP | 1 | **8** (⚠ but see below) |
| AMY | 3 | 3 |
| MAGEA | 5 | 3 (⚠ preset-limited, see below) |
| RFPL | 3 | 3 |
| SRGAP / HERC | 1 | **1 each** ✅ correctly singleton |

The GOLGA6L component contains exactly the three LOC "6-like protein 7" genes §6c named
(LOC115931294 / LOC134757625 / LOC101137218) plus three product-less LOC genes.

### ⭐⭐⭐ The forward-only guard costs NOTHING once nodes are co-oriented

```
202 aligned pairs -> 58 E_r edges (both strands) | 58 forward-only   (difference: 0)
```

⭐⭐ **This confirms §6c's diagnosis and closes the question.** §6b's "+23.08 points from dropping the
guard" measured **my reference-orientation extraction**, not the guard. Extract nodes in gene/hit
orientation and the guard is **free**: it blocks 0 within-family edges while remaining available to
reject genuine antisense. ⟹ **the fix is node construction, and the guard stays on.**

### ⚠⚠ TWO TRAPS THIS RUN EXPOSED

⛔ **"NPIP: 1 → 8 copies" would have been WRONG.** The NPIPB11 seed is a **104,913 bp gene span**, and
**7 of the 8 discovered intervals contain NO NPIP-family gene at all** (they are annotated titin-like,
sortilin-related-receptor-like, multidrug-resistance-protein-1); only the seed's own interval carries
NPIPB11. ⟹ **the 105 kb seed found 8 homologous ~100 kb SEGMENTAL DUPLICATION BLOCKS, not 8 gene
copies.** ⭐⭐ **NODE GRANULARITY DETERMINES WHAT "COPY" MEANS**: an 8 kb CDS envelope (AMY) returns gene
copies; a 105 kb gene span returns SD blocks. State the node unit with every copy count.

⚠⚠ **THE DISCOVERY PRESET DETERMINES THE COPY COUNT — record it like `-p`/`-N`.** From the same MAGEA1
seed on the same 4 Mb window: **`asm20` returns 3 copies, `-k 11 -w 5` returns ~11–12** at
**identity 0.76–0.85 and query coverage up to 0.965**. That band is *inside* asm20's nominal ~20%
divergence range, but its `k=19 / w=10` seeding does not find it. ⟹ **extend the standing rule: ALWAYS
RECORD THE PRESET (`-x`/`-k`/`-w`) ALONGSIDE `-p` AND `-N`.** ⚠ **§6a's "gorilla has exactly 3 AMY
copies, no unannotated copy" used `asm20` and is therefore under-powered by exactly this gap** — a
sensitive per-contig re-run is the test, recorded separately.

---

## §6e — the sensitivity test: AMY survives, MAGEA is 26 copies, and presets are NON-NESTED (2026-08-28)

> ⚠ **PARTLY CORRECTED BY §6f.** The AMY, RFPL and MAGEA-multi-copy results STAND. But the preset
> gap is larger than stated (asm20 3 vs -k11 -w5 **18** on the window, not ~11–12) and its MECHANISM is
> wrong: seeding contributes **0** of the gap — the gate is SCORING (`-s` with the mismatch penalty).
> The naming breakdown here is computed from the same truncated index §6f retracts.

§6d left one live threat: every genome-wide count used `asm20`, which was shown to miss MAGEA copies at
identity 0.76–0.85. Test: re-run all four small seeds at **`-k 11 -w 5 -p 0 -N 500`, one contig at a
time** (per-contig keeps memory ~1.6 GB; `k=11` on the whole genome OOMs at 24.9 GB). 644,905 records.

| seed | records | copies `-k11 -w5` | copies `asm20` |
|---|---|---|---|
| **AMY** | 20,418 | **3** | **3** |
| RFPL | 587,948 | 3 | 3 |
| GOLGA6L | 36,510 | 6 | 7 |
| **MAGEA** | 29 | **26** | 3 |

### ✅ §6a's AMY claim SURVIVES the sensitivity attack

The sensitive scan returns **the same 3 intervals** — LOC101133335 / AMY2A / AMY2B, all annotated,
**nothing unannotated, nothing off-cluster** — despite 20,418 raw records. ⟹ the 0.60–0.80 identity gap
§6c flagged is **closed for AMY**: **gorilla has 3 AMY copies under both `asm20` and `-k11 -w5`,
genome-wide, `-p 0 -N 500`.**
⭐ RFPL is the strongest demonstration of the coverage clause doing real work: **587,948 raw records
collapse to 3 copies.** The repeat noise is enormous and the clause removes essentially all of it.

### ⭐⭐⭐ MAGEA HAS 26 COPIES, AND A SYMBOL GREP SEES 26.9% OF THEM

26 intervals, all on `NC_073247.2` (X-linked, as MAGEA should be), identity 0.7275–1.0000,
coverage 0.513–1.000, **0/26 above 50% soft-masked (max 9.7%) — none is a repeat artefact.** About half
(13/26) are near-full-length at coverage ≥ 0.95; the rest are partial copies.

| how the copy is annotated | n |
|---|---|
| symbol-named `MAGEA*` | 7 |
| LOC, product **"melanoma-associated antigen"** | 9 |
| LOC, **no product at all** (pseudogene) | 6 |
| **NO annotation whatsoever** | **3** |
| other product text | 1 |

⚠⚠⚠ **A THIRD NAMING CONVENTION — and it defeated the fix §6d proposed.** The same family, in the same
GFF, is written `MAGEA9` → *"MAGE family member A9"* for symbol genes and `LOC129530018` →
*"melanoma-associated antigen 9"* for LOC genes: **different words entirely, not a hyphen**. §6d's
product-keyword census searched `"MAGE family member A"` and returned exactly **5** — it missed **all 9**
melanoma-associated-antigen genes. ⟹ **the census I built to escape the symbol trap fell into the same
trap one level down.** A symbol grep sees **7/26 = 26.9%**; the product grep sees 5/26.
⟹ ⭐⭐ **No name-based membership test is trustworthy. Only sequence.**

⭐ **O3 material**: **3 copies carry no annotation at all**, and 6 more are pseudogenes with no product —
present in the reference, absent from the naming.

### ⚠⚠ PRESETS ARE NON-NESTED — NEITHER DOMINATES

⛔ **The sensitive scan MISSED a GOLGA6L copy that `asm20` found** (`NC_073240.2:24,974,895-24,979,577`,
LOC101153208 "golgin subfamily A member 6-like protein 24"), while finding 6 of the other 7. ⟹ **`-k11
-w5` is not a strict superset of `asm20`**; each preset finds copies the other misses.
⟹ **a copy count is a property of (seed, preset, `-p`, `-N`, node unit) — report all five, and prefer the
UNION over any single preset.** This extends the standing rule beyond `-p`/`-N`.
⚠ It also means **§6e's own "AMY = 3" is a two-preset agreement, not a proof** — the strongest available
statement, but a third preset could still differ.

---

## §6f — SECOND AUDIT: §6d is largely an artefact of its own construction (2026-08-28)

A second 15-agent audit attacked §6d/§6e. **Four claims REFUTED, one OVERSTATED.** Every count again
reproduced exactly; the inferences again did not. The single deepest finding:

### ⛔⛔⛔ THE TWO-STAGE DESIGN IS CIRCULAR — STAGE 2 MADE NO DECISION AT ALL

Stage 1 admits an interval iff it covers **≥ 0.50 of the seed**. Stage 2 then asks whether pairs of those
intervals cover ≥ 0.50 of each other. ⭐⭐ **All 26/26 nodes are 100.00% covered by their own seed's
alignment blocks (the node IS the merged target span), so 58/58 within-seed pairs are FORCED to clear
0.50** — verified by interval-overlap lower bound, min stage-1 query coverage 0.5084. **All 58 possible
within-seed pairs aligned and all 58 became edges; 0/267 cross-seed.** ⟹ **every family is a complete
clique because fragmentation was impossible, and the constant predictor "edge iff same seed" reproduces
both arms exactly.** ⛔ **RETRACT §6d's "7 components, one per seed, zero cross-seed contamination" as
evidence for anything about E_r** — the cross-seed zero is real but low-information (0/144 aligned
cross-seed pairs, Wilson-95 upper 0.0260).

⟹ **A discover-then-define pipeline cannot validate its own edge rule unless stage 2 applies a criterion
stage 1 did not already enforce.** This is the general lesson and it should be checked against the
shipped RNA path.

### ⛔⛔⛔ "THE GUARD COSTS NOTHING" IS A 0-OF-0 — AND IT IS BACKWARDS

⛔ **RETRACT §6d's "58 vs 58, the forward-only guard costs nothing."** The guard's **candidate count is
ZERO**: **0/202 aligned pairs supply a reverse-strand record clearing identity 0.60 AND coverage 0.50**;
max reverse coverage in the panel is **0.1532** (0.0592 among the 58 edges) — **3.26×** below the floor.
`discover.py` extracts every interval in the orientation of its hit to ONE per-family seed (audited:
0/26 strand-mixed merges, 0 loci discovered under two seeds), so same-seed nodes are **co-oriented by
fiat**. It is the same *check the candidate count* trap as §6c, in a new costume.

⭐⭐ **The measurement this panel DOES support points the OTHER WAY.** Relabelling the same PAF to
reference orientation (alignment-free; coverage and identity are revcomp-invariant) drops the forward arm
**58 → 26**, and the 32 pairs the guard rejects are **32/32 WITHIN-SEED, 0/32 cross-seed** — *every*
rejection it can make here **destroys a true edge** (null: E[lost] = 29.00, **P(lost = 0) = 1.9e-06**).

⛔⛔ **§6d's recommendation "the fix is node construction, and the guard stays ON" is WRONG IN
DIRECTION.** The shipped `--from-genome` path already stores reference-oriented reps and **force-disables**
the guard, locked in on 2026-08-19 by `genome_mode_grouping_keeps_an_inverted_duplication` and
`refine_params_default_is_orientation_agnostic` — whose message reads *"if you are here to flip this
default, that is the bug."* ⚠ **I proposed flipping exactly that default.** The tree was already right.

⚠ **And the proposed fix is not available to a de novo caller**: co-orientation needs one seed per family
— *the assignment the graph is meant to compute* — and deriving orientation from the graph instead
(2-colouring each component by alignment-strand parity) requires precisely the reverse-strand records the
guard forbids.
⚠ **Provenance**: no minimap2 CMD log survives for `disc_allvall.paf`, and it retains **530 full-length
diagonal self hits**, which is NOT `-X` behaviour — the quoted preset is unverified.

### ⛔⛔ NPIP: THE INFERENCE WAS BACKWARDS — THEY *ARE* GENE COPIES

⛔ **RETRACT "7 of 8 NPIP intervals contain no NPIP gene, therefore they are SD blocks not gene copies."**
Both halves fail. **(1) A name-negative is zero-information here**: only **5** NPIP-family gene records
exist genome-wide, so P(a random 104,913 bp window contains one) = **0.000198**; expected count among the
7 blocks is **0.0014** and observed is 0. **The constant predictor "contains no annotated NPIP gene"
scores 8/8** against the claim's 7/8 — the readout measured NCBI's annotation sparsity, not block content.
**(2) The blocks ARE gene copies**: CIGAR-projecting NPIPB11's 27-block, 12,609 bp CDS union through each
discovery record recovers **0.8520–0.9935 of the CDS union in 6 of the 7 non-seed intervals** (26–27 of 27
blocks). ⟹ ⭐⭐ **these are 6 unannotated NPIP gene copies — an O3 result — not segmental duplication
blocks.** ⚠ The "node granularity determines what a copy means" moral was drawn from a false premise;
it may still be true, but this panel is not evidence for it.

### ⛔⛔ THE NAMING STORY: RIGHT PHENOMENON, WRONG MECHANISM, BROKEN INSTRUMENT

⛔ **RETRACT the hyphen mechanism and the 7/26 figure.** `gene_index.pkl` was **truncated**: 22,706
entries holding **one mRNA-level `product=` per gene**, with **no CDS products, no pseudogenes (7,027),
no lncRNAs (7,330)**. ⛔ **The hyphen claim is false on its own exemplars** — read from the GFF's **CDS**
`product=`, the symbol-bearing genes carry the hyphenated form themselves: RFPL1 = *"LOW QUALITY PROTEIN:
ret finger protein-like 1"*, RFPL2/3 = *"ret finger protein-like 2/3 isoform X1"*, GOLGA6L7 = *"LOW
QUALITY PROTEIN: golgin subfamily A member 6-like protein…"*. ⭐ With a complete census, **symbol-only
intervals = 0/26 = 0.0000** (symbol 11/26, keyword 16/26), against the **7/26** I reported off the broken
index. ✅ **What survives**: RFPL is 9/9 and GOLGA6 14/14 LOC-named, so a *symbol* grep genuinely misses
most members, and 3 MAGEA copies have no annotation at all. ⛔ **But "both name routes fail identically,
split by a hyphen" is withdrawn** — it was an artefact of my own truncated parse.

### ⚠ PRESET SENSITIVITY: TRUE AND UNDERSTATED, BUT THE MECHANISM IS WRONG

✅ The headline holds and is **stronger** than stated: on the 4 Mb MAGEA window, `-x asm20` → **3** copies
and `-k11 -w5` → **18** (not "~11–12"; that figure could not be reproduced at any `-p` and was most
likely a **unit swap** — annotated loci touched, counted as copies).
⛔ **But "asm20's `k=19/w=10` seeding misses them" is FALSE: seeding contributes 0 of the 15-copy gap.**
`-x asm20 -k11 -w5` still gives **3**, while `-k19 -w10` with default scoring gives **16** and
`-x asm20 -A 2` gives **16**. ⟹ **the gate is SCORING (the `-s` minimal-chain-score threshold interacting
with the mismatch penalty), not seeding.** The practical rule from §6e stands — record the preset — but
the stated cause must be corrected.

---

## §6g — an honest precision/recall table, and the 0.50 coverage floor TESTED (2026-08-28)

§6f established that §6b/§6d could not measure the rule: their node sets were built by discovery from a
seed, so stage 2's edges were entailed by stage 1's gate. **This panel is built differently and is
non-circular by construction**:

- **Nodes come from ANNOTATION, not discovery** — the CDS envelope of every annotated member, extracted
  in **gene orientation**. Stage 2 can therefore disagree with stage 1, and it does.
- **Truth comes from normalised `product=` text over the FULL GFF** (2,374,329 feature lines; CDS
  products, pseudogenes and lncRNAs included — the parse §6f showed was missing before).
- **Negatives are SIZE-MATCHED** unrelated genes (576 of them; median 2,421 bp vs the families' 2,547).

⚠ **The rule could fail here, and did**: it misses 10 of 72 members and produces 3 candidate false
merges. That is what makes the numbers meaningful.

### Families, by product text (symbol-blind)

| family | members | LOC-named |
|---|---|---|
| MAGEA | 35 | 40% |
| GOLGA6 | 16 | 88% |
| RFPL | 13 | 69% |
| NPIP | 5 | 80% |
| AMY | 3 | 33% |

### Precision / recall at the SHIPPED floor (identity ≥ 0.60, coverage ≥ 0.50 of the shorter)

| family | members | annotated pairs | edges | pair recall | Wilson 95% | component recall |
|---|---|---|---|---|---|---|
| AMY | 3 | 3 | 3 | **1.0000** | [0.4385, 1.0000] | **3/3** |
| RFPL | 13 | 78 | 32 | 0.4103 | [0.3078, 0.5211] | 8/13 |
| NPIP | 5 | 10 | 4 | 0.4000 | [0.1682, 0.6873] | 3/5 |
| GOLGA6 | 16 | 120 | 32 | 0.2667 | [0.1957, 0.3521] | 7/16 |
| MAGEA | 35 | 595 | 152 | 0.2555 | [0.2221, 0.2920] | 27/35 |
| **ALL** | **72** | **806** | **223** | **0.2767** | **[0.2469, 0.3086]** | **48/72 = 0.6667** |

**Precision**: 223 edges, **0 cross-family**; cross-family false-merge rate **0/195 eligible** of 1,750
possible, Wilson 95% **[0, 0.0193]**. Family × size-matched negative: **3/184 eligible** = 0.0163,
Wilson 95% [0.0056, 0.0468]. **Component purity 1.0000 — 0 components mix two annotated families.**

⚠ **The 3 "false merges" are probably TRUE positives.** All three are the same partner,
`LOC115931109` at `NC_073242.2:29,435,860` — chr16, in the SULT1A / CLN3 / PDXDC1 neighbourhood, i.e. the
classic **NPIP segmental-duplication region** — hitting three different NPIP genes at **identity
0.9411–0.9674, coverage 0.7284–0.8862**. Annotated *"serine/arginine repetitive matrix protein 2-like"*.
A 95%-identity relationship over 3/4 of a 20 kb span is a paralogue, not an error. ⟹ **report as
CANDIDATE false merges; the true rate is plausibly 0/184.**

### ⭐⭐⭐ THE 0.50 COVERAGE FLOOR IS TOO CONSERVATIVE — 0.30 STRICTLY DOMINATES IT

§6c found the decision gap EMPTY, so 0.50 was untested. **On this panel the band is populated on the
positive side and empty on the negative side**, which is exactly what makes the floor testable:

- **positives with coverage in [0.30, 0.70]: 142**
- **negatives with coverage in [0.30, 0.70]: 0 / 379** (Wilson 95% upper on the in-band negative rate
  **0.0151**) — the negatives are **bimodal**: repeat noise below 0.30, and the 3 suspect NPIP pairs
  at 0.73–0.89.

⟹ **Across the whole [0.30, 0.70] band the floor trades recall against NOTHING.**

| floor | pair recall | precision | component recall | purity | mixed comps |
|---|---|---|---|---|---|
| 0.20 | 0.4355 | 0.9832 | 0.8611 | 1.0000 | 0 |
| **0.30** | **0.3821** | **0.9904** | **0.8611** | **1.0000** | **0** |
| 0.40 | 0.3151 | 0.9883 | 0.7639 | 1.0000 | 0 |
| **0.50 (shipped)** | 0.2767 | 0.9867 | 0.6667 | 1.0000 | 0 |
| 0.70 | 0.2060 | 0.9822 | 0.6528 | 1.0000 | 0 |

⭐⭐ **0.50 → 0.30 gains +85 true edges and +19.4 points of component recall (0.6667 → 0.8611) at
UNCHANGED purity (1.0000, 0 mixed components) and NO additional false positives** — precision actually
*rises* to 0.9904. **GOLGA6 goes 7/16 → 16/16 and NPIP 3/5 → 5/5 (both complete); RFPL 8/13 → 11/13.**
[0.20, 0.30] is a plateau, so **0.30 is the conservative end of the optimum.**

⭐ **The coverage clause itself is validated: AUC = 0.9808** (n_pos 386, n_neg 379). It is the right
quantity; only the threshold is miscalibrated.

⚠ **Scope, stated plainly.** (i) n = **5 families / 72 members / one genome** — this is a calibration
signal, not a universal constant. (ii) Truth is annotation-derived and annotation is unreliable — 3 MAGEA
copies genome-wide carry no annotation at all (§6e), so the recall denominator is itself incomplete.
(iii) **MAGEA stays 27/35 at every floor** — 8 members are unreachable by coverage at any threshold, so
the residual miss is not a threshold problem. (iv) ⚠⚠ **This does NOT transfer to the RNA path
unexamined**: §4s measured that **60.51% of the RNA edge set falls below the 0.50 floor under a 1.49×
rep-length inflation**, so the same change there has a different and unmeasured consequence. Test the DNA
arm first.

---

## §6h — ⛔ RETRACT §6g's coverage-floor recommendation: the register already refuted it (2026-08-28)

⛔⛔ **§6g's "lower the E_r coverage floor 0.50 → 0.30" is WITHDRAWN before it reached the pipeline.**
`docs/NEGATIVE_RESULTS_REGISTER.md` §4.1 (**REFUTED / RETRACTED — High redo risk**) already contains the
identical experiment:

> *"Relax the E_r coverage floor from 0.50 to 0.30 to recover lost true pairs — cov≥0.30: **402 T / 132 F
> (P 0.753)** vs shipped **340 T / 62 F (P 0.846)** — **62 true bought for 70 false**. And 127/180 of the
> target pairs have coverage EXACTLY 0.00 — no floor can admit a missing record."*

⚠⚠ **The register flags this row "High redo risk" — i.e. as a thing likely to be re-proposed — and I
re-proposed it.** The standing rule *CONSULT THE REGISTER BEFORE PROPOSING AN APPROACH* exists for
exactly this, and I consulted it only after recommending the change.

### Why §6g's panel disagreed: its negative set was the weak link

§6g's decisive observation was *"142 positives in [0.30, 0.70] against **0 of 379 negatives**"*, which
made the floor look free to lower. That zero is a property of **how I drew the negatives** — 576
size-matched *randomly chosen annotated genes* — not of the rule. Re-scoring against a genuinely
adversarial negative set (every genome-wide hit of the four seeds at identity ≥ 0.60, **644,905
records**, which contains every repeat bridge in the genome):

| coverage band | positives | negatives |
|---|---|---|
| [0.30, 0.40) | 2 | **1** |
| [0.40, 0.50) | 0 | **1** |
| [0.50, 0.60) | 4 | **3** |
| [0.60, 0.70) | 1 | 0 |

⭐ **5 negatives fall in [0.30, 0.70]** where the size-matched panel had **0/379**, and moving the floor
0.50 → 0.30 admits **2 more positives against 2 more negatives — a 1:1 trade**, consistent in sign with
the register's better-powered 62-for-70. ⟹ **the shipped 0.50 stands.**

⚠ **Neither negative set is clean, and this bounds the whole exercise.** The size-matched panel is too
EASY (0 in the decision band); the genome-wide set is contaminated the OTHER way — it labels a hit
"negative" unless it lands on a same-family *annotated* gene, but §6e showed **3 MAGEA copies carry no
annotation at all**, so the 19 negatives at [0.8, 0.9) and 5 at [0.9, 1.0) are very likely **unannotated
true copies**. ⟹ **a defensible floor experiment needs a negative set that is adversarial AND
annotation-independent; neither of mine is.** That, not the floor value, is the open problem.

### What §6g still supports

✅ The **P/R table itself stands** — it is the first non-circular measurement in §6 (nodes from
annotation, truth from full-GFF product text, negatives size-matched, and the rule visibly fails):
**pair recall 223/806 = 0.2767 [0.2469, 0.3086], component recall 48/72 = 0.6667, component purity
1.0000, 0 cross-family edges (0/195 eligible, Wilson 95% upper 0.0193).**
✅ **Coverage is the right quantity — AUC 0.9808** (n_pos 386, n_neg 379).
✅ **MAGEA stays 27/35 at EVERY floor**, so its residual miss is not a threshold problem — the register's
companion line explains it: *pairs with coverage exactly 0.00 cannot be admitted by any floor.*
⛔ Only the **recommendation** to move the floor is withdrawn.

⚠ **LESSON (third instance this session of the same class).** §6b's precision rested on 17 possible
partners; §6d's guard result was a 0-of-0; §6g's floor result rested on a negative set with an empty
decision band. **All three were zero/near-zero-candidate results read as strong ones.** Before quoting a
rate: state the candidate count, and state what the control's construction EXCLUDES.

---

## §6i — ⛔ RETRACT §6c's "the identity clause is inert": I measured the WRONG ESTIMATOR (2026-08-28)

⛔⛔⛔ **§6c's headline — "E_r's 0.60 identity floor sits below minimap2's emission floor and is
structurally incapable of firing", and its corollary "E_r is really a one-clause rule" — is WITHDRAWN.**
It was computed with **gap-compressed identity (`1 − de:f`)**. **The shipped rule does not use that
estimator.** `denovo_pipeline.rs:3091` reads:

```rust
if nm / bl >= params.sensitive_identity && cov >= params.min_coverage {
```

i.e. **matches / block-length**, in which every gap BASE counts as a mismatch, whereas `de:f` charges a
gap of any length once. minimap2's scoring floor `B/(A+B)` bounds the **gap-compressed** quantity, not
`nm/bl` — so the emission-floor argument simply does not apply to the clause the code evaluates.

Measured across every PAF in this investigation:

| PAF | records | **< 0.60 by `nm/bl` (SHIPPED)** | < 0.60 by `1 − de:f` | min `nm/bl` |
|---|---|---|---|---|
| `pr/ff.paf` | 14,739 | **735** | 0 | 0.2373 |
| `pr/fn.paf` | 268 | **11** | 0 | 0.4655 |
| `seq/disc_allvall.paf` | 234,086 | **16,139** | 0 | 0.1478 |
| `fam_allvall.paf` | 89,447 | **2,950** | 0 | 0.2014 |
| `seq/sensitive_all.paf` | 644,905 | **7,519** | 0 | 0.2071 |

⟹ ⭐⭐ **Under the shipped estimator the identity clause fires on 1.2–6.9% of records and reaches as low
as 0.1478. It is NOT inert, it is NOT dominated by the aligner, and E_r remains a genuine TWO-CLAUSE
rule.** ✅ **What the clause actually does is now clear and is a real function: it rejects INDEL-HEAVY
alignments** — high gap-compressed identity but low base-level identity — which is exactly the
repeat-bridge / spurious-chain shape the coverage clause alone would not catch.

⚠⚠ **Consequences, stated explicitly:**
- ⛔ **"Four independent confirmations were four measurements of the same constant" is withdrawn** for any
  result computed with the shipped `nm/bl` rule. ⚠ The RNA figures (`0/728` in false_omission,
  0.749–0.803 in def_failures) must each be re-checked for **which estimator produced them** before being
  either defended or re-retracted — I have not done that, and must not assume either way.
- ⛔ The proposal to **verify inertness end-to-end by expecting a byte-identical catalog is dropped** —
  the prediction is now expected to FAIL, and running it would have burned the run to learn what one
  `grep` of line 3091 shows.
- ⭐ **Cause: `amy/dna_define.py::ident()` prefers the `de:f` tag and only falls back to `nmatch/blen`.**
  Every offline number in §6a–§6h used it. ⚠ **The estimator must MATCH the shipped rule, or the offline
  panel is measuring a different definition than the pipeline implements.** Re-derived offline results
  are hypothesis generators (T8) — this is the failure mode that makes them so.
- ⚠ **Not everything moves**: for the small within-family tallies the audit checked, the 17-pair result
  was identical under both estimators (8 / 8 / 1 / 0), so §6g's P/R table is only mildly affected —
  735/14,739 = 5.0% of its records shift. **It should be recomputed under `nm/bl` before being quoted.**

---

## §6j — the P/R table, corrected: recall depends on TRUTH GRANULARITY, and the protein tier is real (2026-08-28)

Three corrections to §6g, in order of size.

### 1. ✅ The shipped-estimator recompute changes NOTHING

§6i retracted §6c because the offline helper used `1 − de:f` while the shipped rule uses `nm/bl`.
Rescoring §6g's whole panel with the shipped predicate copied verbatim from `denovo_pipeline.rs:3091`:
**223 edges, pair recall 0.2767, component recall 48/72, purity 1.0000, 0/195 cross-family, 3/184
family × negative — every number identical.** The 735 sub-0.60 `nm/bl` records in `ff.paf` either fail
the coverage clause anyway or are not the best-coverage record for their pair. ⟹ **§6c's CLAIM was
wrong; §6g's MEASUREMENTS survive it intact.**

### 2. ⛔ The preset-union hypothesis is REFUTED — and it un-does §6e's "non-nested"

Pre-registered: *presets are non-nested (§6e), so union recall > max(arm) at precision drop ≤ 0.01.*

| arm | edges | pair recall | component recall |
|---|---|---|---|
| `-k11 -w5` (sensitive) | 223 | 0.2767 | 0.6667 |
| `asm20` (core) | 86 | 0.1067 | 0.3889 |
| **UNION** | **223** | **0.2767** | **0.6667** |

⭐ **asm20-only edges: 0. `asm20` ⊂ `-k11 -w5` STRICTLY at the edge level — the union buys nothing.**
⚠ **This also undermines §6e's "presets are NON-NESTED"**: that comparison ran `-k11 -w5` **per contig**
against `asm20` **genome-wide**, so the one GOLGA6L copy `asm20` "uniquely" found is most likely a
chunking artefact, not a preset property. ⟹ **§6e's non-nestedness claim is downgraded to UNVERIFIED**;
the union rule ("prefer the union of presets") is withdrawn.

### 3. ⭐⭐⭐ THE RECALL NUMBER IS A PROPERTY OF THE TRUTH SET'S GRANULARITY

§6g's "MAGEA 27/35" scored the rule for failing to merge genes my truth set had lumped. `"MAGEA"` is in
fact **eight distinct MAGE subfamilies** — MAGE-B 12, MAGE-A 12, MAGE-D 4, MAGE-C 3, MAGE-E 2, MAGE-F 1,
MAGE-H 1 — and MAGE-D/E/F/H are the divergent non-CT MAGEs. **The rule SHOULD refuse those merges.**
GOLGA6 likewise splits into 6L / 6C / 6A (the first audit had already found 6L7 and 6L10 are different
subfamilies). Rescoring at subfamily granularity, same nodes, same rule:

| truth granularity | families | genes | pair recall | Wilson 95% | component recall | purity | cross edges |
|---|---|---|---|---|---|---|---|
| **coarse** (product family) | 5 | 72 | 0.2767 | [0.2469, 0.3086] | 48/72 = 0.6667 | **1.0000** | **0** |
| **subfamily** | 12 | 62 | **0.8209** | [0.7620, 0.8677] | **57/62 = 0.9194** | 0.9286 | 40 |

⭐⭐ **The SAME edge set scores pair recall 0.2767 or 0.8209 depending only on what you call a family.**
Coarse truth buys purity 1.0000 at recall 0.2767; fine truth buys recall 0.8209 at purity 0.9286. ⟹
**neither number is "O1's recall" — the granularity is a MODELLING CHOICE the data does not determine,
and every recall/purity figure in this project must state which one it used.** The 40 cross-subfamily
edges are not obviously errors: MAGE-A↔MAGE-B is a real relationship at the family level.

### ⭐⭐ The protein tier is a genuine, additive win

`protein_tail` (mmseqs, **DEFAULT OFF**) is documented to reach *"coding paralogs past the ~0.65
nt-identity floor where both nucleotide tiers find no edge"* — exactly the **420/806 = 52.11%** of
within-family pairs that have **no nucleotide alignment at all** under either preset (vs 163/806 = 20.22%
that have a record but fall below the coverage floor). Scored on this panel (annotated CDS translated
in-frame — an **UPPER BOUND** on the shipped tier, which must guess the ORF by 6-frame translation):

| granularity | nucleotide | **+ protein** | new edges |
|---|---|---|---|
| coarse, component recall | 0.6667 | **0.7361** | 19 — **19 within-family, 0 cross** |
| subfamily, component recall | 0.9194 | **0.9516** | 11 — 8 within, 3 cross |
| subfamily, pair recall | 0.8209 | **0.8607** | |

⭐ **NPIP goes 3/5 → 5/5 (complete) and RFPL4A-family stays 8/8; coarse purity stays 1.0000.**
⚠ At subfamily granularity it adds 3 cross-subfamily edges (purity 0.9286 → 0.9231), so it is additive
but **not free** — and mmseqs **is installed here**, so the tier is runnable.
⟹ **Worth an end-to-end pipeline A/B with `protein_tail` ON.** ⚠ Untested through the binary; this panel
is annotation-derived, and the shipped tier's 6-frame ORF guess will do worse than the in-frame CDS used
here.

---

## §6k — TBC1D3: subfamily as a LEVEL of one tree, selected by the widest merge gap (2026-08-28)

> ⚠ **HEDGE FROM §6l**: the merge-height *divergence ladder* below must NOT be claimed as
> information beyond identity — partial correlation of height with an independent protein divergence
> proxy, controlling for DNA identity, is **−0.0465 (p = 0.8919)**. The ladder is a presentation of
> `1 − identity`, not a second measurement. The LEVEL SELECTION (widest merge gap) is unaffected.

The user's scoping: subfamilies are meaningless at n=2 and meaningful for clusters like TBC1D3. The
register agrees arithmetically — a split needs n≥3 and **348/494 = 70.45% of GGO families are 2-copy**,
which is why hierarchy-as-catalog-refinement was refuted (inert on 80.57%). This section applies it to
**one clustered family**, which is the domain where it can work.

### Locating the family: gorilla TBC1D3 has 18 copies to human's 9

⚠ **Gorilla TBC1D3 is 0/9 by symbol** (the three `TBC1D3*` names are TBC1D30/31/32, different genes) and
**no gorilla TBC1D3 coordinates existed anywhere in the project**. Located here by homology from the 9
human CHM13 copies (`tbc1d3/hsa_tbc1d3.fa`).
**Preset recorded in full** (§6e's rule): `-x asm20 -p 0 -N 500`, per contig (bounded ~1.1 GB RSS),
identity ≥ 0.60 by **`nm/bl`** (the shipped predicate, §6i), query coverage ≥ 0.50, merge slop 5 kb.

| tier | n | identity | coverage | location |
|---|---|---|---|---|
| **full-length** | **11** | 0.9713–0.9785 | 1.000 | NC_073228.2 |
| partial (TBC1D3G-like) | 5 | 0.8422–0.9054 | 0.507–0.678 | NC_073228.2 |
| distal | 2 | 0.8792–0.9043 | ~1.000 | NC_073243.2 |

⭐ **The identity distribution is BIMODAL with an empty gap from 0.9054 to 0.9713** — a two-level
structure visible before any clustering is run.

### ⭐⭐⭐ The widest merge gap selects the level, and the level is right

Copies extracted **in hit orientation** (co-oriented), all-vs-all `-c -X --no-long-join -k11 -w5`,
average linkage on **d = 1 − identity × coverage** (the composite that won the parallel variant test;
plain 1 − identity did **worse**, contradicting `family_hierarchy.py`'s design note that coverage "would
only add noise" — on these genomic spans it does not).

| rank | gap | between | selects |
|---|---|---|---|
| **1** | **0.0427** | h = 0.0534 → 0.0961 | **k = 6** |
| 2 | 0.0418 | 0.1662 → 0.2080 | k = 2 |
| 3 | 0.0339 | 0.0961 → 0.1300 | k = 5 |

⭐⭐ **At the gap-selected k = 6, all six clusters are PURE** with respect to an independent
full-length/divergent tier label, and **the 11 full-length copies form EXACTLY one cluster (set equality
verified, |cluster| = |full| = 11).** The level is chosen by a property of the tree — the widest gap in
merge heights — not by a tuned threshold, which is the Canzar-shaped form.

**Divergence ladder**: the recent expansion coalesces entirely **below h = 0.0534**; then the widest gap
in the tree; the older copies join between **h = 0.0961 and h = 0.2080**. ⟹ merge heights order the
expansion, and the 11-copy recent burst is separated from the ancient copies by the largest single step
in the tree.

### ⚠ Why this works here and NOT on the pooled panel

On the 72-gene / 12-subfamily panel (§6j), **no variant achieved label purity at or below the true group
count**: the NPIP-validated `1 − identity` distance was pure only at **k = 44–72** — degenerate, since
k = 72 is all-singletons and pure by construction — and the best variant reached purity 1.0000 only at
**k = 20** for 12 true subfamilies (0.5806 at k = 12–13, **0.0000** at k = 2–5). The 2026-08-05 NPIP
result held *below* the true group count (pure at k = 2..8 for 2 groups); pooled, the floor sits at 3.7×
it.
⟹ ⭐⭐ **A single global tree level does NOT transfer across families — subfamilies in different families
formed at different depths, so one cut cannot serve all.** The unit is ONE CLUSTERED FAMILY, exactly as
the user scoped it. Report the level per family, selected by that family's widest merge gap.

⚠ **Scope**: n = 1 family. The purity label (full-length vs divergent) comes from the *human*-vs-gorilla
alignment while the tree comes from *gorilla*-vs-gorilla, so the tree never saw it — but both are
sequence-based, so this is a genuine external check, not a fully independent one. Copy counts remain
preset-dependent (§6e). ⚠ **NC_073228.2 is a 195 Mb contig; whether it is the human-chr17 ortholog is NOT
established here** and must not be asserted.

---

## §6l — the hierarchy audit: what the tree does and does not add (2026-08-28)

Nine agents built five hierarchy variants on the 72-gene panel and adversarially checked the winner.
Every headline number reproduced from raw `ff.paf` via the shipped `nm/bl` predicate. Four findings, two
of which matter more than the hierarchy itself.

### ⭐⭐⭐ 1. THE 52% "MISSING RECORDS" GAP IS MOSTLY CORRECT BEHAVIOUR

§6j reported **420/806 = 0.5211** of within-coarse-family pairs having no nucleotide record at all, and
read it as the recall ceiling. Decomposed by stratum:

| stratum | pairs with a record | missing |
|---|---|---|
| **within-SUBFAMILY** | 188/201 = 0.9353 | **13/201 = 0.0647** |
| cross-subfamily, within coarse family | 198/605 = 0.3273 | **407/605 = 0.6727** |
| cross-coarse-family | 195/1,750 = 0.1114 | — |

⭐⭐ **The 52% is almost entirely a CROSS-SUBFAMILY phenomenon. Within a true subfamily only 6.47% of
pairs lack a record.** ⟹ **most of what §6j called O1's recall ceiling is the rule correctly declining to
link genes that are not recent duplicates.** The real within-subfamily miss rate is 6.47%, not 52%.

### ⭐⭐ 2. IDENTITY ALONE IS NOT A DISTANCE — COVERAGE CARRIES THE SIGNAL

Medians of the best-coverage record: **within-subfamily identity 0.7915 vs cross-COARSE identity
0.7906 — indistinguishable — while coverage differs 23× (0.9811 vs 0.0421).** Of the 100 shortest
distances under `1 − identity`, only 63 are within-subfamily; **21 are cross-COARSE**, and average
linkage merges those first (e.g. `MAGEC3—LOC115933515`, MAGE vs GOLGA6L26, d = 0.0363 — shorter than
MAGE-A's own formation height).
⟹ this is why `1 − identity` scored ARI 0.5849 while `1 − identity×coverage` scored **0.8484**, and it
is a **fourth independent confirmation that the coverage clause is the load-bearing half of E_r**.
⭐ **Gating pairs by the SHIPPED edge rule before building the tree lifts ARI 0.5849 → 0.8374** at the
shipped constants (not tuned) — the edge rule is exactly the missing distance filter.

### ⛔ 3. MERGE HEIGHTS DO NOT TRACE DIVERGENCE BEYOND IDENTITY

Heights correlate strongly with an independent protein proxy (**Spearman −0.8909, p = 0.0002, n = 11**),
so they are not noise. **But the partial correlation controlling for DNA identity is −0.0465,
p = 0.8919**, while the reverse partial is +0.6823, p = 0.0207.
⟹ ⛔ **merge height is essentially `1 − identity` relabelled; "tracing divergence" is a RESTATEMENT, not
an additional measurement.** A divergence ladder may still be a useful *presentation*, but it must not be
claimed as information the identity column does not already carry. ⚠ **This hedge applies to §6k's
TBC1D3 divergence ladder too.**

### ⚠ 4. THE HIERARCHY'S LIFT OVER A TRIVIAL BASELINE IS SMALL, AND PURITY IS NEARLY VACUOUS

- **Lift over the best trivial baseline is +0.0351 ARI** (average linkage 0.8484 at k=20 vs
  **single linkage 0.8133** — i.e. one pairwise threshold on identity×coverage). That difference is
  **exactly two gene attachments** (MAGEB5, MAGEC1).
- ⚠ **"17/17 pure clusters" is not the honest denominator.** At k=20, 4 of those 17 are size-1 and cannot
  be impure; **RFPL4B is pure BY CONSTRUCTION** (the any-alignment graph has 6 components and RFPL4B is an
  isolated 2-node one, so no linkage rule could merge it) and MAGE-E likewise. ⟹ **the honest score is
  0 errors / 12 genuine opportunities, and 2/12 non-singleton subfamilies are recovered by construction.**
- ⚠ **best_k is ORACLE-selected** by maximising ARI against the labels being scored ⟹ 0.8484 is an
  UPPER BOUND, not a prediction. (Null control passes: 1,000–2,000 label permutations, null mean 0.0253,
  p95 0.0569, p ≤ 0.0005 — the signal is real; it is the LEVEL that fails.)
- ✅ **Missing-pair fill policy is INERT**: four policies (worst×1.0, ×1.5, 1.0, ×3.0) give identical
  best_k, ARI and purity, though only ×1.5 reproduces the merge topology byte-for-byte.

### ⭐⭐ 5. CIRCULARITY IS REAL AND LARGE — BUT A SEQUENCE-FREE CHECK EXISTS AND PASSES

⚠⚠ The labels are **similarity-derived by construction** ("golgin subfamily A member 6-**LIKE** protein
9"). Pushing the **protein** channel — the very evidence those names were written from — through the
identical pipeline reaches **ARI62 = 0.9042, purity62 = 1.0000, numerically identical to the DNA
headline**, and **192/201 = 0.9552** of the within-subfamily pairs are visible to it. ⟹ **the label-shuffle
p-value tests against NO STRUCTURE, not against the circular alternative**, and "labels.tsv is independent
of any alignment" is true only of *this* alignment.
⭐ **Not refuted, because a genuinely sequence-free criterion exists**: **genomic position**, which was
never an input to minimap2. Against position-only 2 Mb tandem-array blocks the k=20 composite scores
**ARI 0.6333 — exactly equal to the annotated labels' own 0.6333** — beating protein 0.4073, identity
0.3930 and length-only 0.3770. ⟹ **position is the non-circular validator; use it, not the product names,
whenever a subfamily claim must be defended.**

### No transferable height (confirms §6k)

Only **3 of 5** coarse families have ≥2 annotated subfamilies (AMY and NPIP have exactly one, so the
question is vacuous there). Of the 3: GOLGA6 admits a window (0.06846, 0.42226], RFPL (0.02594, 0.97024],
and **MAGE is DEGENERATE — its subfamilies form at or above the height at which the whole 72-gene tree
closes.** **Intersection over the informative families: EMPTY.** ⟹ **the level must be read PER FAMILY**,
as §6k does for TBC1D3.

---

## §6m — ⛔ RETRACT §6i, RESTORE §6c: I read the wrong function (2026-08-28)

⛔⛔ **§6i is withdrawn.** It retracted §6c's "the identity clause cannot fire" on the grounds that the
shipped rule uses `nm/bl`. **That line — `denovo_pipeline.rs:3091` — is inside `tier2_rescue`, the
tier-2 admission path, NOT the E_r edge rule.**

**The main E_r decision is `denovo_pipeline.rs:4665–4673`:**

```rust
let de = f[12..].iter().find_map(|x| x.strip_prefix("de:f:")...);
let ident = match de { Some(d) => 1.0 - d, None => nmatch / alnlen };
```

i.e. **`1 − de`, with `nmatch/blocklen` only as a fallback** — matching what the params certificate has
always emitted (`identity_metric = "1-de (fallback nmatch/blocklen when de:f: absent)"`,
`denovo_pipeline.rs:4407`) and what `docs/METHOD_PSEUDOCODE.md:177` has always documented.

⭐ **The fallback is DEAD CODE: 0 of 984,574 records across 8 PAFs are missing the `de:f:` tag
(0.000000).** minimap2 `-c` always emits it. ⟹ the main path always uses `1 − de`, which is exactly the
quantity §6c measured.

### ✅ §6c is RESTORED, and now rests on 8 PAFs

**0 of 984,574 records fall below the 0.60 floor**, and the observed minimum **tracks the scoring
preset** with the biology varying freely:

| preset | emission floor `B/(A+B)` | observed min `1 − de` |
|---|---|---|
| default (`A=2, B=4`) | 0.667 | 0.6313 / 0.6437 / 0.6500 / 0.6591 / 0.6862 / 0.7621 |
| `asm20` (`A=1, B=4`) | 0.800 | **0.8291 / 0.8295** |

⟹ ⭐⭐ **on the E_r path the identity clause is structurally incapable of firing: its floor sits below
minimap2's emission floor.** This independently confirms `NUMBERS.md`'s existing line — *"identity never
binds; the COVERAGE clause is what separates a repeat from a paralogue"* (49.0% of random 30 kb pairs
clear identity ≥ 0.60) — and now supplies the mechanism for it.

### ⚠⚠ A REAL DEFECT THIS EXPOSED: two paths, two estimators

`tier2_rescue` (line 3091) applies `nm/bl >= sensitive_identity` **directly**, with no `de:f:` branch,
against the same constant 0.60. The two estimators are not interchangeable — `de` is gap-compressed while
`nm/bl` charges every gap base — and under `nm/bl` the same constant **does** fire, on **1.2–6.9%** of
records, reaching **0.1478**.
⟹ **the tier-2 admission path applies a materially stricter identity test than the definition it is
admitting into, while appearing to use the same threshold.** Tier-2 is currently uncommitted/opt-in, so
nothing shipped is wrong — but this must be reconciled before tier-2 is enabled, either by routing it
through the same `1 − de` branch or by documenting the difference as deliberate.

⚠ **Process note.** §6i was written from ONE `grep` hit without checking which function contained it, and
it retracted a correct finding. The lesson from §6i stands in corrected form: **the offline panel must
compute the shipped predicate — and "the shipped predicate" means the one on the path under test, which
must be identified by reading the enclosing function, not by matching a line.**

---

## §6n — the precision/recall table with error DIAGNOSTICS (2026-08-28)

Computed with the **main E_r predicate** (`1 − de`, `denovo_pipeline.rs:4665`; §6m). ⭐ **The `nm/bl`
variant used in §6g–§6l gives the IDENTICAL 223-edge set — 0 pairs differ** — so those tables are
unaffected by the §6i/§6m confusion. Panel: 72 annotated gorilla CDS-envelope DNA nodes, gene-oriented,
576 size-matched negatives.

| truth granularity | pair recall | Wilson 95% | pair precision | component recall |
|---|---|---|---|---|
| **subfamily** (12 non-singleton) | **165/201 = 0.8209** | [0.7620, 0.8677] | 165/223 = 0.7399 | **57/62 = 0.9194** |
| **coarse family** (5) | 223/806 = 0.2767 | [0.2469, 0.3086] | **223/223 = 1.0000** | 48/72 = 0.6667 |

⚠⚠ **Precision and recall move in OPPOSITE directions with granularity because both are computed on the
SAME 223 edges.** The 58 edges that count as false positives at subfamily level are true positives at
coarse level. **Neither row is "O1's precision/recall" — always name the granularity.**

### Where the 36 missed within-subfamily pairs go

| outcome | n | fraction of 201 |
|---|---|---|
| edge recovered | 165 | 0.8209 |
| **FN — record exists, fails COVERAGE** | 23 | 0.1144 |
| **FN — no alignment record at all** | 13 | 0.0647 |
| FN — fails IDENTITY | **0** | 0.0000 |

⭐ **Identity accounts for zero misses**, consistent with §6m (it cannot fire on this path). The
recoverable stratum is the 23 coverage failures; the 13 with no record are beyond any threshold.

### ⭐⭐ The "false positives" are NOT repeats and NOT errors

Of 223 edges, **58 cross a subfamily boundary and 0 cross a coarse-family boundary** — i.e. **zero
genuine errors on this panel.** All 58 join subfamilies of the same family, at **identity 0.9458–0.9922,
coverage 1.000, alignment blocks 2,785–7,868 bp** (e.g. `GOLGA6L10~GOLGA6L19` at 0.9922 / 1.000 /
7,834 bp; `RFPL1~RFPL3` at 0.9458 / 1.000 / 3,383 bp).
⭐ **They are not repeat-driven: 0/58 have a node above 50% soft-masked** (median max-masked 0.254).
⟹ **these are genuine near-identical paralogue pairs that the ANNOTATION splits into different
subfamilies** — the numbering ("6-like protein 9" vs "…19") is finer than the sequence supports. **The
label is the questionable object here, not the edge.**

### The 3 family × negative merges: repeat-rich, but the blocks are too long for repeats

| family node | partner | identity | cov | block | masked (node / partner) |
|---|---|---|---|---|---|
| NPIPB11 | LOC115931109 | 0.9411 | 0.886 | **23,932 bp** | 0.531 / 0.651 |
| LOC101134557 | LOC115931109 | 0.9500 | 0.763 | 19,949 bp | 0.488 / 0.651 |
| LOC101141990 | LOC115931109 | 0.9674 | 0.728 | 14,878 bp | 0.422 / 0.651 |

⚠ **Both readings are partly right, and the block length decides.** The partner is **65.1% soft-masked**,
so "repeat-rich" is true — but a dispersed element cannot produce a **14.9–23.9 kb** collinear block at
94–97% identity (Alu ≈ 300 bp, L1 ≈ 6 kb). ⟹ **this is a segmental duplication, i.e. a real paralogous
relationship that the annotation names "serine/arginine repetitive matrix protein 2-like"**, sitting in
the chr16 SULT1A/CLN3/PDXDC1 NPIP SD region. **Classify as annotation disagreement, not rule error** —
and note the three family nodes are 1.9 Mb away, 74.5 Mb away, and on a different contig, so this is not
a local artefact.

### ⭐ Truly NEW unannotated copies: 3, and only in MAGEA

Genome-wide sensitive scan, per seed, copies passing E_r stage-1 (identity ≥ 0.60, cov_q ≥ 0.50):

| seed | copies | overlap an annotated gene | **unannotated** |
|---|---|---|---|
| AMY | 3 | 3 | **0** |
| GOLGA6L | 6 | 6 | **0** |
| RFPL | 3 | 3 | **0** |
| **MAGEA** | **26** | 23 | **3** ⭐ |

⟹ **the only candidate reference-present-but-unannotated copies in this panel are 3 MAGEA copies.** AMY,
GOLGA6L and RFPL are fully accounted for by the annotation — which is itself the useful negative result:
**the method is not manufacturing copies.**

---

## §6o — "no alignment record" is mostly HARMLESS; the coverage floor is what actually breaks families (2026-08-28)

§6n reported the 201 within-subfamily pairs as 165 edges / 23 coverage failures / 13 no-record. That
reads as though the 13 are the hopeless stratum. **It is the other way round.**

**What "no record at all" means mechanically**: minimap2 emitted no line for the pair — not a weak line,
none. Its minimizer seeding found too few shared exact k-mers to build a chain, so **the edge rule never
evaluated the pair**. That is a *seeding* failure, upstream of both clauses, and it is why no threshold
can address it.

### ⭐⭐ But a missed PAIR is not a lost MEMBER — the family routes around it

| missed stratum | pairs | **members still in the same component** |
|---|---|---|
| no alignment record at all | 13 | **12/13 = 0.9231** |
| record exists, fails COVERAGE | 23 | **12/23 = 0.5217** |

⟹ ⭐⭐⭐ **the pairs the aligner never saw are almost entirely benign, while the pairs it DID see and the
coverage floor rejected are what actually separate members.** Of the 12 genuinely separated pairs,
**11 align but fall below coverage** and only **1** has no record. This inverts the natural reading of
§6n's table and redirects any future repair effort at the coverage stratum, which is at least *visible*.

### ⚠ The damage is about REDUNDANCY, not subfamily size

⛔ **My "small subfamilies break" hypothesis is refuted: only 1/12 = 0.08 of separated pairs sit in a
subfamily of ≤3 members.**

| subfamily | members | pairs | missed | separated | placed |
|---|---|---|---|---|---|
| MAGE-B | 12 | 66 | **23** | **0** | **12/12** |
| MAGE-A | 12 | 66 | 0 | 0 | 12/12 |
| RFPL4A | 8 | 28 | 0 | 0 | 8/8 |
| **NPIP** | 5 | 10 | 6 | **6** | **3/5** ⛔ |
| **MAGE-D** | 4 | 6 | 5 | **5** | **2/4** ⛔ |
| **MAGE-E** | 2 | 1 | 1 | **1** | **1/2** ⛔ |

⭐ **MAGE-B misses 23 of its 66 pairs (34.8%) and loses NOBODY**; NPIP misses 6 of 10 (60%) and loses 2.
⟹ **what matters is whether the surviving edges keep the subfamily CONNECTED, not how many pairs are
missed.** A pair-level recall of 0.8209 therefore understates member-level recovery (component recall
0.9194), and the two must not be used interchangeably.

### Are these members impossible to find? Almost never — n = 1 of 201

- **12 of 13** no-record pairs are recovered **transitively**, through other members.
- The **protein tier rescues 2 more** — both NPIP (`LOC101141990 ~ LOC109028586` and `~ LOC115934567`),
  which is the mechanism behind §6j's NPIP 3/5 → 5/5.
- ⭐ **Exactly one pair resists everything tried: `MAGEE1 ~ MAGEE2`** — no nucleotide record at any
  preset **and no protein hit**. It is the only genuinely undetectable pair in the panel, and it is
  fatal only because MAGE-E has **2 members**, so there is no third copy to route through.

⟹ **the honest statement is: 1 of 201 within-subfamily pairs (0.50%) is beyond every method tried here.
"No record" is not a synonym for "impossible" — it is a synonym for "not seen by THIS seeding", and in a
family with redundancy it costs nothing.**

---

## §6p — O2 preconditions on the fibroblast substrate: threading is possible, the population is 40 reads (2026-08-28)

Question: with O1 in better shape, does O2 work — can PSVs on the same molecule still be threaded through
the variation graph? **Threading is implemented and the molecules can carry it. The problem is that
almost nothing needs it on this substrate.**

⚠ **Stale reference corrected**: the register cites `psv_linkage.rs:292` for "`--vg` emits ZERO
copy-assignment by default". **That file no longer exists**; the PSV machinery is in
`src/rustle/vg_family/copy_assign.rs` (1,070 lines) and `src/bin/copy_assign.rs` (2,392).

⭐ **Threading IS implemented**: `read_psv_obs(read, psv_positions)` returns the read's allele at **every**
PSV column, `logl` combines them, and `n_decisive` counts how many discriminate, gating on
`resolvable = n_decisive >= 1`. Multi-PSV molecules are the designed input, not a missing feature.

### ⛔ TBC1D3 has the structure but no expression — the standing scope is EMPTY for O2

Among the 11 full-length gorilla TBC1D3 copies (§6k): median pairwise identity **0.9784** over
**12,790 bp** ⟹ **~276 differences per pair, one every ~46 bp**. A 500 bp read would span ~10.
**Ideal threading substrate.**
⛔ **But the fibroblast FLNC BAM carries 9 primary reads across all 18 copies, with 12 copies at zero.**
TBC1D3 is testis-biased. ⟹ ⚠⚠ **the standing scope (TBC1D3 + fibroblast) is viable for O1, which works
on genomic sequence, and EMPTY for O2, which needs molecules.** This must be stated when the scope is
quoted.

### Where reads exist, there is no ambiguity to resolve

FLNC primaries at the annotated family nodes (fibroblast, `-F 2308`): MAGE **11,156** over 9 expressed
nodes, NPIP **2,623** over 5/5, RFPL 315, GOLGA6 62, AMY 21.
⭐ Note the two deepest are **MAGE-D (4/4 expressed)** and **NPIP (5/5)** — precisely the two subfamilies
O1 separated in §6o. But score margins kill the O2 opportunity there:

| locus | reads | MAPQ 0 | **near-tie s2/s1 ≥ 0.95** | ≥ 0.80 |
|---|---|---|---|---|
| MAGED1 | 6,019 | 0 | **0 (0.0000)** | 0 |
| LOC101137992 | 3,121 | 0 | **0 (0.0000)** | 0 |
| LOC109028586 | 2,423 | 0 | **0 (0.0000)** | 0 |
| **NPIPB11** | 114 | 4 | **40 = 0.3509** | 113 = 0.9912 |

⭐⭐ **11,563 reads at the three deep loci have ZERO score competition — every one is MAPQ 60 with no
second chain within 20%.** There is no assignment problem to solve at these loci; minimap2 has already
placed every molecule unambiguously. **NPIPB11 is the only locus with a genuine O2 target**, and its
**0.3509 near-tie rate matches the reframed O2 population** (alignment-score near-ties, ~21.75%
genome-wide) rather than MAPQ-0.

### ⭐ At NPIPB11 threading would work — on 40 molecules

The 40 near-tie reads are **mean 2,551 bp (688–5,681), spanning a mean of 2.4 introns**, i.e. long
spliced molecules. Their competitors are other NPIP copies on the same contig (1,989 secondary records in
the window). NPIPB11's near-tie partners sit at identity **0.9518–0.9682**, i.e. 3–5% divergent ⟹ a
2,551 bp exonic read crosses on the order of **85–125 discriminating positions**. ⟹ **the molecules
amply carry multi-PSV evidence; nothing about read length or divergence blocks threading here.**

### ⚠⚠ The structural finding: what needs O2 and what is expressed are anti-correlated

**The family with the PSV structure to thread (TBC1D3 — 11 copies at 0.9784, one difference per 46 bp)
has 9 reads. The families with thousands of reads (MAGE-D, NPIP's deep nodes) have ZERO near-ties,
because the expressed copies are the divergent ones minimap2 can already separate.** Their intersection
on this substrate is **one locus and 40 molecules.**
⟹ **O2 cannot be evaluated at scale on fibroblast FLNC.** Any O2 claim needs either a tissue that
expresses the near-identical copies (testis, for TBC1D3 — a different animal, so out of the current
scope) or a substrate where similar copies are transcribed. ⚠ This is a POWER statement, not a refutation
of threading: n = 40 supports a demonstration, never a rate.

---

## §6q — O1↔O2 synergy: NOT a coupled definition, but O1's edge identity as an O2 PLANNER (2026-08-28)

### ⛔ First, what is already closed: coupling the two definitions

Every attempt on record returns **zero**:

| attempt | result |
|---|---|
| Joint DNA/RNA object as **definitional** (union or intersection of E_dna, E_rna) | **0 blocks in O1, 0 of 7,691 reads in O2, 0 in O3** |
| Intron-chain comparison to "fix the O1⊥O2 leak" | branch fires **0 times on every substrate on record** |
| Iterative joint estimator (EM, longcallR-style) | **0 disagreements among 3,081** co-committed reads |
| Combining orthogonal levers in one classifier | in-sample 0.915 → grouped LOFO-CV **0.80 ± 0.05 vs 0.840** for the best single feature |

✅ **Verified in code**: `copy_assign.rs` contains **zero** references to λ, min-cut, γ or `E_r`. O1 → O2 is a
pure **conditioning** with no leakage, and that orthogonality is load-bearing — it is what stops the
definition from presupposing its own members (`a VG presupposes its members`, register §4.1).
⟹ **do not couple the definitions.**

### ⭐⭐ What DOES work: O1 already computes the quantity that predicts where O2 has work

`METHOD_PSEUDOCODE.md` records that the graph γ runs on is **unweighted** — "identity and coverage decide
whether an edge exists" and are then discarded. **That discarded identity is an excellent predictor of
O2's difficulty.** Per node, take the maximum E_r identity to any same-family partner, and compare with
the node's observed alignment-score near-tie rate (`s2/s1 ≥ 0.95`) among its fibroblast FLNC reads:

| node class | nodes | reads | **read-weighted near-tie rate** |
|---|---|---|---|
| max within-family identity **≥ 0.95** | 14 | 1,919 | **0.7577** |
| max within-family identity **< 0.95** | 6 | **11,606** | **0.0000** |

**Spearman ρ = 0.7173** (node-level, n = 20, threshold-free); Pearson r = 0.5718.
⭐⭐ **The negative direction is the strong one: 11,606 reads at low-max-identity nodes produced NOT ONE
near-tie.** ⭐ The positive extreme is equally sharp — the five nodes at rate **1.0000** all have a partner
at identity ≥ 0.9779, including the `LOC101144987 ~ LOC101146541` pair at **0.9997** (both nodes: every
read a near-tie) and `LOC101148313 ~ LOC115935025` at **1.0000**.

⚠ **Why this is not circular**: O1's identity is measured **rep vs rep** over CDS envelopes; O2's near-tie
is measured **read vs genome**, over shorter, spliced, exon-only molecules competing genome-wide rather
than within the family. Related through divergence, but different objects and different measurements —
so this is a prediction, not an entailment.

⭐ **The subfamily median is the WRONG statistic; the per-node MAX is the right one.** MAGE-D has median
within-subfamily identity **0.8000** — which reads as easy — yet **max 0.9997**, and exactly the two nodes
in that near-identical pair have near-tie rate **1.0000 (687/687 and 655/655)** while the other two MAGE-D
nodes, with 9,140 reads between them, sit at **0.0000**. A family-level summary hides this completely.

### What the synergy buys, without touching either definition

- **Skip work that provably has none**: **11,606 of 13,525 reads (85.8%)** sit at nodes where no near-tie
  exists. Running copy_assign there cannot change an assignment.
- **Predict where abstention will be needed** before running O2, from a number O1 has already computed.
- **Ship it the λ way** — as a per-node REPORTED certificate that decides no membership and no assignment,
  which is the only form of O1/O2 contact this project's record supports.

⚠ **Limits.** n = 20 nodes on one substrate; **the 0.95 cut was chosen after seeing the data**, so the
two-group contrast is fitted and only the threshold-free ρ = 0.7173 is honest without a held-out check.
Two real exceptions with depth sit above the cut at rate 0.0000: **RFPL1 (279 reads, max id 0.9528)** and
**LOC115934567 (50 reads, 0.9545)** — so ≥ 0.95 is a *prior*, not a guarantee.

---

## §6r — HELD-OUT test of the §6q planner: the signal holds, the clean zero does NOT (2026-08-29)

§6q fitted a 0.95 cut on a 20-node panel and got a read-weighted near-tie rate of **0.7577 vs 0.0000**.
I flagged the threshold as fitted. Tested on a genuinely held-out set — **the shipped fibroblast catalog
(`o1_replicate`): 356 families, 1,070 copies, 2,075 E_r edges, built by the pipeline, not by me** — with
**nothing refitted**:

| group (rule fixed in advance) | copies | reads | near-ties | rate | Wilson 95% |
|---|---|---|---|---|---|
| max within-family identity **≥ 0.95** → predict YES | 217 | 200,558 | 110,218 | **0.5496** | [0.5474, 0.5517] |
| max within-family identity **< 0.95** → predict NO | 853 | **1,661,251** | 144,043 | **0.0867** | [0.0863, 0.0871] |

✅ **The signal replicates**: 6.3× separation, **Spearman ρ = 0.5804** threshold-free over 1,070 copies,
**permutation p = 0.0005** (2,000 shuffles). O1's discarded edge identity genuinely predicts where O2
faces alignment-score near-ties, on data that had no part in choosing the rule.

⛔ **But the strong form is REFUTED.** §6q's "below 0.95, not one near-tie in 11,606 reads" does **not**
hold: the held-out low group carries **144,043 near-ties at rate 0.0867**. That zero was a small-panel
artefact — the training group had **6 nodes**. ⚠ The effect is also weaker throughout (ρ 0.7173 → 0.5804;
0.7577/0.0000 → 0.5496/0.0867).

⚠⚠ **Consequence for the planner: it may RANK, it must not SKIP.** 89.23% of reads sit at copies the rule
predicts to be free of near-ties, but **8.67% of them are not** — skipping that stratum would silently
drop **144,043** contested reads. ⟹ **§6q's "skip work that provably has none" is withdrawn.** The
defensible use is a *prior* on where abstention will concentrate, never a gate.

### What this licenses, and what it does not

✅ Emit the weight. `denovo_pipeline.rs` flattens every edge to 1.0 before the partition (§5q), and the
per-copy output has no identity column at all — so a quantity with **AUC 0.9491** for wanted-vs-suspect
edges (§5q) and **ρ = 0.5804** against O2 difficulty (here) is computed and then dropped from the primary
artifact. It should travel with the copy, in the `families.tsv` idiom: *"APPENDED last, so existing
header-keyed readers keep working unchanged… a REPORT on the emitted family, never an input to it."*
⛔ **Do NOT re-weight the partition.** §5q settled this twice: `induced_density` (`family_split.rs:511`)
discards weights so the γ TEST is unweighted, γ is inert on 79% of families, and since identity weights
are < 1 every weighted density is LOWER — more blocks fall under γ = 0.20 and the partition splits MORE,
with NPIP already over-split 5 ways. Recovering current behaviour would need γ re-tuned to the weight
distribution, i.e. fitting a threshold on the arms it is scored against.
⟹ **the edge weight leaves O1 as a REPORTED per-copy certificate and enters O2 as a prior. It never
re-enters the partition.**

---

## §6s — the three objectives end to end: the chain works, the substrate has no natural instance (2026-08-29)

**User's framing**: O1 says a locus has 2 tandem copies; O2 threads PSVs and assigns reads to them; O3
finds the individual actually had 3 — *"if we had a sample from an individual where the sample matches
the genome we would [know]"*. That condition **holds here**: the fibroblast sample and the assembly are
the same animal (SAMN04003007 / KB3781).

### ⭐⭐ The O3 leg is VALIDATED — by excision, at FPR 0

`o3_excise/` removes a copy from the reference, so the sample genuinely carries one more copy than the
reference shows, and asks where its reads go. Depth at the surviving paralogue, masked genome ÷ complete
genome:

| stratum | n | median ratio | q1 | q3 | >1.10 | max |
|---|---|---|---|---|---|---|
| **copy excised (positive)** | 113 | **1.2744** | 1.0186 | 1.8233 | **0.6018** | **2.45×** |
| control (nothing excised) | 225 | **1.0000** | 1.0000 | 1.0000 | **0.0000** | 1.0015× |

⭐ **AUC = 0.8924**, and at the ratio admitting **zero** controls (> 1.0015): **TPR = 0.7965 at FPR =
0.0000.** ⟹ **when a copy is genuinely absent from the reference, its reads relocate onto the surviving
paralogue and the excess is detectable.** The architecture the user describes is sound and demonstrated.

### ⛔ But O2's ambiguity is NOT a proxy for the missing copy

Lining all three objectives up on the 356-family fibroblast catalog:

| family stratum | n | median depth skew (max/median reads) | median max identity |
|---|---|---|---|
| O2 near-tie rate ≥ 0.20 | 58 | **1.49** | **0.9949** |
| O2 near-tie rate < 0.20 | 298 | **1.56** | **0.7720** |

⭐⭐ **Depth skew is the SAME in both groups (1.49 vs 1.56) while identity differs enormously (0.9949 vs
0.7720).** ⟹ **O2's near-tie rate tracks copy SIMILARITY, not a missing copy** — confirming §6r — and the
hypothesised chain "an unmodeled copy causes both O2 ambiguity and depth excess" **does not appear as a
correlation.** The two signals are orthogonal on this catalog.

### ⭐⭐⭐ Why RNA depth cannot carry the O3 leg — the numbers that decide it

**The absorbed-copy signal is ~1.27× (median, above).** Between-copy expression variance in the same
catalog reaches **109.96× (GWFAM225) and 55.75× (GWFAM163, 13 copies)**. ⟹ **a 1.27× effect cannot be
recovered from a between-copy RNA depth comparison; it is buried two orders of magnitude down.**
⭐ **The excision panel works precisely because it compares the SAME locus to ITSELF** (masked reference
vs complete), which cancels expression exactly. ⟹ **the O3 leg must be a WITHIN-LOCUS comparison (one
locus, two references) or DNA depth — never a between-copy RNA comparison.** This is a design conclusion,
not a tuning one.

### ⚠ And on this animal the natural yield is pre-registered at < 1

Because the WGS animal **is** the assembly animal, **CN_WGS = CN_assembly + collapse_deficit**, and the
WGS-only residual **is** O3's absorbed stratum — **pre-registered at < 1 collapse** (§4u). ⟹ **the "there
were really 3" case is expected to be essentially absent in this genome**, which is why the chain has
never had a natural instance to run on and why the excision panel is the right instrument.
⟹ **Report the architecture as demonstrated-by-construction (TPR 0.7965 / FPR 0.0000 on 113 excisions),
never as a discovery rate on this substrate.**

### Fixed here

⭐ **`tier2_rescue` now uses the E_r estimator.** It applied `nm/bl >= sensitive_identity` while the
definition applies `1 − de` with `nm/bl` only as a fallback (§6m) — against the same 0.60 constant, so
tier-2 was silently STRICTER than the definition it admits into (the raw ratio drops 1.2–6.9% of records,
down to identity 0.1478). Its comment claimed "the SAME rule", which is how it hid. Now assembles the
`de`-preferring identity exactly as `nucleotide_edges_scored` does.
⭐ **`tier2_rescue` no longer re-types the E_r tier flags** — it builds the command from `ER_TIER_FLAGS`
and `ER_SENSITIVE_SEED` (excluding `-X`, which does not apply when aligning candidates against reps). The
anti-drift guard had caught this; it was a genuine second literal copy.

---

## §6t — ⚠ CORRECT §6s: the excision ratio is an ORACLE statistic, not a detection rate (2026-08-29)

§6s reported the O3 leg as "validated by excision: **TPR 0.7965 at FPR 0.0000, AUC 0.8924**". **That is
not a detection rate and must not be quoted as one.**

`o3_excise/align.sh` states the design: *"PAIRED alignment: the SAME simulated reads against the COMPLETE
genome and the MASKED one. Depth at a landing site with the paralogue present (complete) vs absent
(masked) isolates the redistribution."*

⟹ ⛔ **the masked/complete ratio requires BOTH references — i.e. it requires already knowing where the
copy was removed.** A real O3 caller possesses only the actual assembly, which *is* the masked arm. The
statistic cannot be computed in the situation it is meant to detect. This is the same class of error the
session caught three times already: a quantity conditioned on the answer.
⚠ Two further scope limits: the reads are **simulated at 40×** (`sim_40x.fastq.gz`), not real fibroblast;
and the run is `-ax map-pb --eqx -N 50 -p 0.1` — record the preset with any figure from it.

### ✅ What the panel DOES establish — and it is still the load-bearing result

**The redistribution mechanism is real and quantified.** Same reads, both references, so library
composition cannot explain the difference:

| stratum | n | median masked/complete | >1.10 | max |
|---|---|---|---|---|
| copy excised | 113 | **1.2744** | 0.6018 | 2.45× |
| control | 225 | **1.0000** | **0.0000** | 1.0015× |

⭐ **When a copy is absent from the reference its reads are ABSORBED by the surviving paralogue**, and the
controls move not at all (max 1.0015×). ✅ This is a **mechanism confirmation**, and it agrees with the
independent figures already in `docs/o3_missing_copy_evidence.md`: FixItFelix's ~1.5× for collapsed
regions and this project's own §5 measurement of 1.75×. ⚠ It is whole-genome, so the doc's
"mini-reference manufactured the signal" critique does **not** apply to this panel — but the doc's other
finding stands and is the decisive one: **the shipped detector fires 0/915 on the matched whole-genome
arm**, because a read that does not fit locus *u* on 3.6 Gb finds a better primary home elsewhere instead
of piling onto *u*.

### ⭐⭐ This STRENGTHENS the design conclusion rather than weakening it

Because the paired-reference statistic is unavailable in practice, the O3 leg **cannot** be a depth
comparison at all — not between copies (expression variance reaches **109.96×** against a **1.27×**
signal, §6s) and not between references (requires the answer). ⟹ **it must read structure WITHIN the one
locus it can see.**

✅ **The shipped code already does exactly this**, and no change is needed: `absent_copy.rs` admits a
collapsed copy from **co-varying allele clusters inside the host locus** (gate 1: `n_clusters >= 3`;
gate 2: `min_p_distinct`; gate 3: strand-symmetry, rejecting A→G editing artefacts; gate 4: alleles
placeable; gate 5: the synthetic copy must NOT remap to its host at ≥ 98%), and anything RNA cannot
settle returns **`DnaNeeds`** rather than a guess — the module's own header says the separation *"needs
DNA"*. ⟹ **the architecture the user described is implemented correctly; what the excision panel supplies
is the mechanism behind it, not a rate for it.**

---

## §6u — ⭐⭐⭐ DIRECT PROOF of the advisor's phenomenon, inside ONE animal (2026-08-29)

**Question**: is there proof that a reference misses copies an individual carries? **Yes — and it does not
need a second individual.** KB3781 is haplotype-resolved, so its maternal and paternal assemblies are two
haploid genomes of the same animal. Same probe, same parameters, both haplotypes
(`o3_hapcnv/{mat,pat}.chr.fa`, `-x asm20 -p 0 -N 500`, per contig; identity ≥ 0.60 by `1 − de`,
cov_q ≥ 0.50, 5 kb merge).

⚠ **Design note — why this avoids the 2026-08-19 failure.** That run was uninformative because its
span-matched controls were not composition-matched (null 0.1512 vs signal 0.0278). Here there are **no
random-interval controls**: the same probe is counted on two genomes, so preset sensitivity cancels in
the DIFFERENCE, and single-copy genes calibrate it.

### ⛔ First, the confound the calibration caught

Four probes gave MAT *n* / PAT **0**: MAGEA1 3→0, MAGED1 1→0, **MAGEE1 1→0, MAGEH1 1→0**. The last two
were my single-copy CONTROLS — and a control showing the same effect as the test is how you know the
effect is not the one you are testing for. **All four are X-linked and have ZERO alignment records in the
paternal haplotype** (autosomal MAGEF1 gives 1 and 1). ⟹ **KB3781 is male; the paternal haplotype carries
Y, so the X's content is absent by SEX, not by copy number.** Discard all X-linked probes.

### ⭐⭐⭐ The autosomal result: a tandem array of 3 on one haplotype, 2 on the other

| probe | MAT | PAT |
|---|---|---|
| AMY (LOC101133335) | 3 | 3 |
| NPIPB11 | 8 | 8 |
| RFPL1 | 3 | 3 |
| **GOLGA6L7** | **6** | **7** |
| calibration: APOBR / CLN3 / NUPR1 / MAGEF1 | 1 | 1 |
| calibration: PDXDC1 / SULT1A2 (not single-copy, but concordant) | 5 | 5 |

**6/6 autosomal controls concordant; 1 of 4 autosomal families differs, by exactly one unit.** The copies
correspond one-to-one across haplotypes — matching sizes and identities to 3 decimal places
(4,682/4,686 bp both at id 0.924, cov 0.748; 3,182 bp both at id 0.8944, cov 0.508) — **except in one
cluster**:

| MAT (2 units) | PAT (3 units) |
|---|---|
| 100,405,873 · 6,035 bp · id 0.9672 | 104,789,835 · 6,033 bp · id 0.9670 |
| 100,446,711 · 6,036 bp · id 0.9675 | 104,830,654 · 6,037 bp · id 0.9677 |
| — | **104,871,462 · 6,035 bp · id 0.9673 · cov 1.000** |

⭐ **The extra unit is not a fragment**: full length, coverage 1.000, and identity indistinguishable from
its two siblings. ⭐⭐ **The spacing is a perfectly regular tandem array — 40,819 and 40,808 bp between
consecutive units on PAT** (MAT's two are 40,838 apart). ⟹ **the paternal haplotype carries a third
tandem repeat unit that the maternal haplotype does not.**

### ⟹ What this proves, and why "sample matches genome" does not rescue it

**This is exactly the scenario: a locus where one haploid genome shows 2 copies and the individual
actually has 3.** And it needs no second animal — it is within KB3781.
⭐⭐ **The primary assembly is a MOSAIC of 16 PAT / 9 MAT chromosomes** (independently re-derived from
exact chromosome lengths, [[project_o3_haplotype_cnv_uninformative]]). It therefore carries **one
haplotype per chromosome** — so on any chromosome drawn from one parent, the other parent's
haplotype-specific copies are **structurally invisible, for the very animal the assembly was built
from.** ⟹ **"the sample matches the genome" does NOT mean "the genome holds all the sample's copies":**
a diploid individual has two haplotypes and a reference has one per chromosome. Note the probe's own
locus returns identity **1.0000** in PAT, confirming GOLGA6L7 is a paternal-haplotype locus in the
primary.

⚠ **Scope.** One family, one array, n = 1 difference among 4 autosomal families tested. GOLGA6 is the most
SD-rich family in the panel, so an assembly-resolution difference between haplotypes is the alternative
explanation — argued against by the unit's full length, full coverage, sibling-matched identity and the
regular 40.8 kb periodicity, but not excluded. ⚠ Record the preset with the count (above); §6e showed copy
counts are preset-dependent, which is why only the mat-vs-pat DIFFERENCE at identical settings is quoted.

---

## §6v — scaling the haplotype scan: a rate, and why it must be read as an UPPER BOUND (2026-08-29)

§6u established the phenomenon at one locus. Next step: how often? One gene-scale probe per **autosomal
multi-copy family** from the shipped fibroblast catalog, both haplotypes, identical parameters
(`-x asm20 -p 0 -N 500`, per contig; identity ≥ 0.60 by `1 − de`, cov_q ≥ 0.50, 5 kb merge).

⚠ **Probe selection matters and is stated**: of 356 catalog families, **307 are autosomal multi-copy**;
**206 were dropped because their longest copy exceeds 20 kb** — §6d showed a large genomic span returns
**segmental-duplication blocks rather than gene copies**, so counting them would answer a different
question. **101 probes remain, median 8,187 bp.**

| outcome | n |
|---|---|
| concordant (MAT == PAT) | **93** |
| **different** | **8** |

**Rate = 8/101 = 0.0792, Wilson 95% [0.0407, 0.1486].**

### ⚠⚠ But the DIRECTION test — which I built as a control — fires

| direction | n |
|---|---|
| MAT > PAT | **7** |
| PAT > MAT | **1** |

Sign test **p = 0.0703**. Real heterozygous CNV should be roughly balanced; this is not. The cause is
visible in the totals:

| quantity | value |
|---|---|
| total copies summed over all 101 families | **MAT 245 / PAT 231** → ratio **1.0606** |
| assembly size | mat 3.507 Gb / pat 3.352 Gb → ratio **1.0464** |

⭐⭐ **The copy-count ratio tracks the assembly-size ratio (1.0606 vs 1.0464).** ⟹ **most of the MAT excess
is plausibly assembly bulk/redundancy, not biology**, and the largest differences (−5, −3, −3) sit in the
≥10-copy families where assembly resolution bites hardest (GWFAM327 24/19, GWFAM300 14/11, GWFAM321 6/3).
⟹ **8/101 = 7.92% is an UPPER BOUND on the haplotype CNV rate, not an estimate.** Do not quote it as a
rate. A defensible figure needs assembly-quality-matched haplotypes or per-locus validation.

### ⭐ What survives, and why §6u's case is the credible one

**GOLGA6L7 is PAT > MAT — the MINORITY direction, running AGAINST the assembly bias**, as is the only
other PAT>MAT family (GWFAM259, 21/22). ⟹ **a difference in the under-represented direction is far harder
to explain as assembly redundancy**, which is exactly why §6u's case — full-length 6,035 bp unit, coverage
1.000, identity 0.9673 against siblings at 0.9670/0.9677, regular 40.8 kb periodicity — remains the
load-bearing instance. **The phenomenon is established; its frequency is not.**

⚠ **Reusable design note.** The direction test cost nothing and caught a confound the headline rate would
have hidden. **Any two-assembly comparison must report the direction split and the assembly-size ratio
beside the rate** — a rate whose sign distribution tracks a size asymmetry is measuring the assemblies,
not the biology.

---

## §6w — validating the 8: half are fragment artefacts, and the rate falls (2026-08-29)

§6v's **8/101 = 7.92%** was flagged as an upper bound because the direction split (7 MAT>PAT vs 1) tracked
an assembly-size asymmetry. Each of the 8 was then profiled per copy. **The decisive test: restrict to
FULL-LENGTH copies (coverage ≥ 0.90).** A difference carried by fragments is assembly noise; one carried
by full-length units in a regular array is the §6u signature.

| verdict | n | families |
|---|---|---|
| **ARTEFACT** — difference vanishes on full-length copies | **4** | GWFAM168 8/7→**2/2**, GWFAM259 21/22→**20/20**, GWFAM29 4/3→**1/1**, GWFAM321 **6/3→1/1** |
| unresolved — full-length but dispersed, not tandem | 3 | GWFAM140 6/5, GWFAM300 7/8, GWFAM304 7/6 |
| **REAL** — full-length units in a regular tandem array | **1** | **GWFAM327 24/19** |

⭐ **GWFAM321's headline −3 was entirely fragments** (6/3 → 1/1), as was GWFAM259's +1 (21/22 → 20/20).
⟹ **half the differences were the fragment confound.**

| statistic | rate | Wilson 95% |
|---|---|---|
| §6v, before validation | 8/101 = **0.0792** | [0.0407, 0.1486] |
| after removing fragment artefacts | **4/101 = 0.0396** | [0.0155, 0.0977] |
| tandem-array-confirmed only | 1/101 = **0.0099** | [0.0017, 0.0536] |

### The two credible instances, and why they differ in strength

**GWFAM327 — 24 (MAT) vs 19 (PAT) full-length units**, **zero fragments on either haplotype**, uniform
unit size **1,779–1,823 bp**, identical identity range (0.9103–1.0000 both), and near-identical
periodicity: **median tandem gap 36,440 bp (MAT, n=13) vs 36,386 bp (PAT, n=8)**. ⚠ **But it runs WITH
the assembly-size bias** (MAT > PAT), and a ~1.8 kb unit at 36.4 kb spacing is a **higher-order repeat**,
exactly the structure where assembly length differences are notorious. ⟹ **large and clean, but bias-aligned:
cannot be separated from an assembly-length difference on this evidence.**

⭐ **GOLGA6L7 (§6u) remains the load-bearing case** precisely because it is **+1 in the MINORITY
direction** — against the assembly bias — with a full-length 6,035 bp unit at coverage 1.000 and identity
0.9673 against siblings at 0.9670/0.9677, on a regular 40.8 kb period. **A difference in the
under-represented direction cannot be explained by the redundancy that produces the majority direction.**

### ⟹ Where this leaves the claim

✅ **The phenomenon is established** (§6u), and it needs no second individual.
⚠ **The frequency is not.** The defensible statement is **≤ 4/101 = 3.96% [1.55, 9.77] of autosomal
multi-copy families differ between one animal's haplotypes**, of which **exactly one is confirmed as a
tandem-array difference in the bias-resistant direction**. ⛔ **Do not quote 7.92%** — half of it was
fragments.
⚠ The 3 dispersed survivors (GWFAM140/300/304, all ±1) are genuinely unresolved: full-length, but not in
an array, so an assembly placement difference and a real CNV look the same. Settling them needs read-level
depth at those loci, not more assembly comparison.

## §6x — tier-2 admission: the mechanism finally fires, and the trade is real (2026-08-29)

**§5r's projection was never tested**, because the tier-2 candidate generator inherited a defect.

### ⚠⚠ THE FIFTH CONSTRUCTION ERROR — AND A REPEAT OF ONE ALREADY FIXED

`footprint_skeletons` has two branches. The **windowed** branch blocks a region only when existing reps
cover **≥ 50%** of it (fixed earlier). The **genome-wide** branch — the one `tier2_rescue` uses — still
blocked on **ANY overlap**. At the target loci reps sit **0.7–3 kb away**, so a candidate whose 5 kb
grouping merely touched one was discarded whole.
⟹ **tier-2 had never fired where it was designed to**: every rep at an NPIP locus was tier-1, and its
earlier +3 came indirectly, from tier-2 reps added elsewhere perturbing collapse.
⭐ **LESSON: when a predicate is wrong in one branch, check its siblings before moving on.** The correct
form was 40 lines away in the same function. A test asserting both branches agree on a shared case is
cheaper than a 17-minute run.

### After the fix

| arm | families | copies | **NPIP loci** | **pure** | T2 copies | NPIP reached via T2 |
|---|---:|---:|---:|---:|---:|---:|
| shipped | 83 | 484 | **12/31** | **3** | 0 | 0 |
| tier-2, any-overlap block | 98 | 557 | 15/31 | 2 | 47 | 4 |
| **tier-2, coverage block** | 117 | 631 | **17/31** | **0** | **140** | **17** |

Admissions rise **50 → 238** (the block was suppressing 79% of candidates). ⭐ **The mechanism now works
as designed** — 17 NPIP loci reached via tier-2, against 4 before.

### The trade, stated both ways

⭐ **Recall up and fragmentation DOWN**: 12 → **17/31** loci, and NPIP now spans **3 families instead of
5** (§5j's fragmentation partly resolved).
⛔ **Purity 3 → 0.** Of the 33 impurities in those 3 families: **19/33 lie within 50 kb of a true NPIP
locus (median 6.0 kb, q1 2.8 kb)** — adjacent sequence, plausibly same-locus flanks or gorilla-specific
copies **absent from a human-ortholog-derived truth set by construction**; but **11/33 are >1 Mb away or
on another contig**, which is genuine contamination.

⟹ **NOT shippable as a default.** 140 tier-2 copies genome-wide (22% of the shipped catalog size) are
unvalidated, and 11 distant impurities are real false merges. **Ships OFF, with the `T2~` tid prefix so
every admitted copy is identifiable** — the evidence tier must travel with the copy.
⟹ **31/31 remains unreachable**: 8 loci have NO reads, and §5s showed locus collapse absorbs loci even at
40× coverage. Tier-2 routes around that defect; it does not repair it. **The repair belongs in
`collapse_loci_span_aware`.**

## §6y — The remaining lever, closed; §5s RETRACTED (08-30)

**Shipped win.** `RUSTLE_COLLAPSE_EXONIC` (pre-existing at `family_detect.rs:662`, never measured, and
missing from the params certificate — the recurring M2 defect, now added). Measures collapse containment
on EXONS not the genomic span. Real fibroblast BAM: reps 2,847 → 3,303, families 83 → 94, copies
484 → 550, **NPIP purity HELD at 3/3**, NPIP loci 12/31 unchanged. OFF arm reproduces baseline
byte-identically once the `max_family_identity` column (commit 0d8d606) is excluded — md5 03e1ddec both.

**⛔ RETRACTION 1 — §5s's mechanism.** §5s concluded "the transcripts are built and then absorbed by
locus collapse". **Refuted by ablation.** `RUSTLE_LOCUS_JUNCTION_ONLY=1` removes BOTH sequence-based
clauses (b) span and (c) POA-core; it recovers **exactly the same single locus (28)** as the exonic fix
— 16/31 either way — while families explode 76 → 153 and copies 445 → 699. ⟹ **locus collapse owns
1 of 15 residual loci, not 13 of 14.** The exonic fix therefore captures 100% of what that mechanism can
yield at 94 families rather than 153: the efficient frontier of an exhausted lever.

**⛔ RETRACTION 2 — the top-up substrate is INVALID at the loci it was built to test.** The injected
reads (named `NPIP<n>_`, 1,240 = 31×40, 91.3% at source, 29/31 with ≥20 — §5s's 93.7%/30-of-31 broadly
reproduce) are **18,467 bp median at the MISSING loci vs 6,013 bp at recovered ones**, while REAL reads
at those same loci are **1.2–2.9 kb, indistinguishable from recovered loci**. Long reads span many
junctions ⟹ 27/23/20/28 distinct exact intron chains from ~24 spliced reads ⟹ no chain reaches
`GATE_MIN_READS`, so no node can form BY CONSTRUCTION. ⟹ **"15/31 with 40 reads everywhere" never
measured what it claimed.** ⚠ a top-up must match the real read-length distribution or it tests the
simulator.

**⭐ THE REAL RESIDUAL IS NODE CONSTRUCTION.** Of the 15 missing loci, **13 have NO overlapping node at
all** (any-overlap, not the 50% attribution rule) and 1 more overlaps only 0.01 of a 24.3 kb node — this
despite real support: **L3 = 62 real reads with 16 sharing ONE exact intron chain** (vs a floor of 3),
**L5 = 4,784 real reads with 4,187 on a single chain.** Not the edge rule, not γ, not collapse, not
depth. ⟹ the next question is WHICH construction step drops them; the cause is NOT yet established.

⚠ purity is 0 on BOTH topped-up arms vs 3 on the real BAM — the simulated reads contaminate purity, so
**purity is only interpretable on the real substrate.**

## §6z — Stage attribution: WHICH construction step drops the missing NPIP loci (08-30)

Instrumented the SHIPPED path. `RUSTLE_DEBUG_LOCUS` existed only in `detect_conflict_catalog_genome_wide_xchrom`;
`--homology-primary` runs `detect_homology_catalog_genome_wide`, which had NO trace — a first run produced 0
dbg lines. Ported the hook there (comma-separated loci) plus `gate_reject_reason`, which recomputes the
SHIPPED gate predicate in the shipped order INCLUDING POOLED locus support. Pure logging; build clean.

**The funnel over all 31 projected loci (real BAM, shipped config):**

| stage | loci | median reads |
|---|---|---|
| 1. PASS-1 — no exact chain reaches 2 reads | **11** | **4** |
| 2. GATE rejects every skeleton | 4 | 16 |
| 3. LOCUS COLLAPSE absorbs it | 2 | 86 |
| 4. reaches a REP | 14 | 35 |

14 reps → shipped `copies.tsv` reports 12, so **2 more are lost at the edge/family stage**. 11+4+2+2 = 19 = 31−12. ✔

**⭐ THE DOMINANT CAUSE IS EXPRESSION, NOT A DEFINITIONAL DEFECT.** The 11 Pass-1 losses carry a **median
of 4 reads**; only 2 loci die at Pass-1 with ≥15 reads (37 and 18). This is consistent with the independently
derived **expressed ceiling of 23/31** — these copies are not transcribed in fibroblast at detectable depth.
⟹ **do not report the 19-locus gap as method error.**

**⭐⭐ THE ONE SHARP LEVER: all 26 `read_floor` rejections sit at pooled = 2 against `GATE_MIN_READS` = 3 —
every one is short by EXACTLY ONE READ.** (Gate verdicts over 157 traced skeletons: 44 KEPT, 26 read_floor,
7 noncanonical_or_mixed_strand.) A floor of 2 would admit all 26. ⚠UNTESTED for precision: it must be
scored genome-wide on false merges + NPIP purity before any claim, and the register warns that lowering a
read floor multiplies spurious skeletons.

⚠ collapse costs **2 loci** here (median 86 reads) — consistent with §6y, where the exonic fix recovered 1.
⚠ this supersedes four offline hypotheses that each died to their comparator (collapse-absorption,
all-canonical, projection offset, chain fragmentation): none could model POOLED support. **T8 confirmed —
only the instrumented binary settled it.**

**⚠ PROVENANCE CORRECTION (08-30).** §6y and §6z were first produced in `/mnt/linuxdisk/.../mut_repo`, a
NON-git tree ~3,000 lines behind canonical (`denovo_pipeline.rs` 7,765 vs 10,855 lines, no tier2). The
instrumentation was ported here and the trace RE-RUN on the canonical binary: the funnel reproduces
**EXACTLY** — 11 / 4 / 2 / 14 loci and 44 KEPT / 26 read_floor / 7 noncanonical — so §6z stands on the
shipped pipeline. ⚠canonical's gate has a FOOTPRINT branch and a READ-STRAND argument the stale tree
lacked; `gate_reject_reason` mirrors the shipped clause order including both. **Check `git rev-parse` before
trusting a measurement's tree.**

## §6aa — FLAG-DEFAULT AUDIT: nothing qualifies to flip (08-30)

**Scope.** 124 `RUSTLE_*` flags exist; only 13 are named in any doc. Classified every boolean behaviour
flag on the O1 path by MEASURED verdict.

| flag | measured | default |
|---|---|---|
| `COLLAPSE_EXONIC` | see the retraction below | **stays OFF** |
| `LOCUS_EXON_UNION` | −20 recall points | off |
| `SPLICED_REP` | e2e regression chr7 F1 .570→.411, chr16 .910→.761 | off |
| `TIER2_ADMIT` | 17/31 loci but purity 3→0 | off |
| `ER_WEIGHTED_PARTITION` | refuted; all edges weight 1.0 ⟹ runs unweighted | off |
| `READ_STRAND` | break-even, 224 lost / 163 gained | off (opt-in) |
| `COTHREAD_REP`, repeat gate | priced negative / P1-breaking | off |
| `ER_SUM_COVERAGE`, `SHARED_EXON`, `FLAGFREE_SITES`, `TSS_SNAP` | **UNMEASURED** | off |
| `*_DUMP` `*_AUDIT` `*_TRACE` `LOCUS_JUNCTION_ONLY` | diagnostics | never |

⟹ **every measured flag is negative or break-even; the rest are unmeasured. The shipped defaults are
already correct and NOTHING should be flipped.**

**⛔ RETRACTION — §6y's "purity HELD 3/3" was too narrow to license `COLLAPSE_EXONIC`.** Co-membership
pairs, OFF vs ON: **2,778 → 4,481**. The flag **CREATES 1,991 new co-membership pairs and destroys only
288** — it merges ~7× more than it splits, and the largest family goes **39 → 60**. The rising family count
(83 → 94) hid this because 66 copies were added at the same time. New merges are **333 cross-chromosome**
and 1,658 same-chrom at a **median separation of 14.9 Mb** (max 111 Mb). Dispersed families do span such
distances, so this is not proof of error — but it IS 1,991 unvalidated merges, and a 3-family NPIP purity
check cannot speak to them. ⚠**"purity held" on a 3-family subset is not a precision measurement.**

⚠**AND THE AVAILABLE PANEL CANNOT ADJUDICATE IT.** `bench/family_def_airtight_panel.py` computes edges from
ANNOTATION-derived cDNA homology and never invokes the pipeline (no subprocess, no `RUSTLE_` env) ⟹ it is
**BLIND to a LOCUS-COLLAPSE flag**; running both arms would return "identical", the 2/150 trap exactly.
⚠also: exonic containment is **NOT monotonically stricter** — exon overlap shrinks but so does the
denominator, so it can CREATE merges span-containment refuses. That is what the 1,991 are.
⟹ **to ship it, price those 1,991 merges against an instrument that actually runs the pipeline.**

## §6ab — `COLLAPSE_EXONIC` REFUTED as a default: its merges are hub fusion (08-30)

**New instrument: `bench/merge_precision_arms.py`** — prices the merges a flag CREATES against the
pipeline's OWN output, filling the gap §6aa exposed (the airtight panel recomputes edges from annotation
and never runs the pipeline, so it is blind to any collapse/representative flag).

**The statistic.** A family is an E_r component, so co-membership is transitive BY CONSTRUCTION: two copies
being co-family is not evidence they resemble each other. A few new edges can fuse large components and
manufacture a quadratic number of pairs. So: what fraction of co-membership pairs carries a DIRECT E_r
edge — **read against the baseline arm's own rate.**

| pair set | direct E_r edge |
|---|---|
| OFF baseline (2,778 pairs) | **49.53%** ← the comparator |
| ON, kept from OFF (2,490) | 52.29% |
| **ON, newly created (1,991)** | **14.92%** |

⭐**+154 edges produce +1,991 co-membership pairs, and 85% of them have NO direct homology backing** — a
3.3× drop against the arm's own baseline. **1,244 of the 1,991 (62%) land in ONE family, `GWFAM8`, which
grows to 60 copies.** ⟹**this is the hairball failure mode, not recall**: exonic containment fuses
components through hubs. ⛔**`RUSTLE_COLLAPSE_EXONIC` STAYS OFF — refuted as a default, not merely
unvalidated.**

⚠the 297 backed new pairs are legitimate (median identity 0.7639, min 0.6804) — the flag is not
uniformly wrong, it is that 85% of what it adds is unsupported. ⚠**quote 14.92% only WITH the 49.53%.**

## §6ac — O1 NODE FLOOR: 3 → 2, SHIPPED (08-30)

**The inconsistency.** Pass-1 forms a skeleton at **2** reads; the gate then demanded **3**, discarding
chains the previous stage had accepted. §6z found EVERY read-floor rejection at the NPIP panel sitting at
pooled support = 2 — all 26 short by exactly one read.

**Arms on the real fibroblast BAM.** Control first: floor 3 with `RUSTLE_GATE_MIN_READS` unset is
**BYTE-IDENTICAL** to the shipped baseline (md5 2c002d7c both), so the treatment arm is interpretable.

| | floor 3 | floor 2 |
|---|---|---|
| NPIP loci | 12/31 | **14/31** |
| family purity | 3 | **3 (HELD)** |
| families / copies | 83 / 484 | 121 / 678 |
| reps | 2,847 | 3,598 |

**⚠⚠ THE FIRST VERDICT WAS WRONG — `merge_precision_arms.py` HAD A SIZE CONFOUND.** Raw direct-edge
backing (22.44% new vs 49.53% baseline) printed "do not ship". **Pairs grow as n² while edges grow
linearly, so backing FALLS MECHANICALLY as families grow**, sound merges or not; newly-created pairs sit in
the LARGER families, so the comparison was not like-for-like. This is the register's own "match the SIZE
distribution, not the edge count" trap, walked into while holding the instrument built to avoid it.

**Size-robust statistics (now what the script reports):**

| | pairs per NEW EDGE | 2–5 | 6–15 |
|---|---|---|---|
| floor 3 (baseline) | — | 0.89 | 0.45 |
| **floor 2** | **1.5** | **0.81** | **0.44** |
| COLLAPSE_EXONIC (refuted) | **6.5** | 0.89 | 0.44 |

⟹ floor 2 adds **1,446 GENUINELY NEW EDGES at 1.5 pairs each** — real homology, not amplification — and its
families are built like the shipped ones at matched size. ⭐**SHIPPED as `NODE_MIN_READS = 2`.**

⚠**`GATE_MIN_READS` STAYS 3** — it governs `copy_assign`'s tie-invariance certificate, a different
question (O1 ⊥ O2). The knob feeds only `GateParams::default()`.
⚠⚠**EVERY CATALOG NUMBER RECORDED BEFORE 08-30 WAS COMPUTED AT 3** (the 627-family catalog included); set
`RUSTLE_GATE_MIN_READS=3` to reproduce one.
⚠**A TEST ALREADY ENCODED THE OLD DEFAULT** (`gate_rejects_below_min_reads`) and failed on the flip —
**grep the tests before changing a default.** It now pins the CLAUSE at an explicit floor, and
`o1_node_floor_is_two_and_is_overridable` pins the VALUE. 815 pass / 0 fail.

## §6ad — THE GAP WAS NEVER EXPRESSION: flag-free site construction recovers 26/31 (08-30)

**⛔⛔ RETRACTION — "the missing loci are not expressed" (§6z, and everything built on it).** That read a
POST-ALIGNMENT PRIMARY COUNT as if it were transcription. Measured on the real BAM (which retains
**1,489,757 secondary alignments, MORE than its 1,009,396 primaries**):

| | median primaries | median SECONDARIES |
|---|---|---|
| missing loci (17) | 6 | **404** |
| recovered loci (14) | 33 | 495 |

**Reads align at the missing loci at essentially the recovered rate** — locus 0 has 2 primaries against
**1,119 secondaries**. They lose the primary-flag COIN TOSS among near-identical paralogs. Low primary
count is exactly what redistribution predicts, so it can never evidence silence. ⚠**the "23/31 expressed
ceiling" was derived the same way and is RETRACTED with it** — the flag-free arm reaches 26/31, past it.

**The fix was already in the tree, unmeasured: `RUSTLE_FLAGFREE_SITES`** (the advisor's own objection —
"which placement minimap2 flags primary is close to arbitrary, yet Pass-1 seeds candidate sites from
primaries ALONE"). Every placement proposes a site; a site survives only if the EXISTING gates find reads
agreeing on a chain THERE. Deliberately an ABSTENTION from read→locus assignment: **O1 ("which copies
exist") must not require first answering O2 ("where did this molecule come from")**.

**Arm vs the just-shipped floor-2 default:**

| | floor 2 | + FLAGFREE |
|---|---|---|
| NPIP loci | 14/31 | **26/31** |
| PURE families | 3 | **6** |
| families / copies | 121 / 678 | 501 / 3,905 |
| pairs per NEW EDGE | 1.5 | **2.0** |
| density 2-5 / 6-15 / 16-39 / 40+ | 0.81 / 0.44 / 0.34 / 0.38 | **0.96 / 0.60 / 0.50 / 0.53** |

⭐**Families are DENSER than baseline in EVERY size band** — the opposite of hub fusion (COLLAPSE_EXONIC
ran at 6.5 pairs/edge with falling density). ⭐**Phantom-copy test PASSES**: 3,905 copies over **3,361
DISTINCT genomic footprints = 1.16 copies/footprint vs the baseline's 1.17**, within-family overlapping
pairs 0.005% — flag-free seeding is not duplicating placements, it is finding distinct loci.

⚠**NOT yet flipped**: a 5.8× copy increase genome-wide is a large behavioural change, the NPIP panel is a
HUMAN PROJECTION, and runtime rises (2.5M records vs 1M). Every available precision proxy improved, but
they are proxies. ⟹**recommend ON; awaiting an explicit call.**

## §6ae — SHARED-READ EDGES (E_c) AS AN ADDITIVE TIER: no depth threshold exists (08-30)

**The proposal** (user): emit a definition edge when two loci share reads — the starved-locus logic that
motivated flag-free. Measured on the flag-free arm's 6,818 reps, E_c = "both reps carry an AS-tied
placement of the SAME read" (`read_conflict.rs`'s own definition), 1,558,151 rep-overlapping placements.

**It is NOT circular, and it is NOT redundant.** Circularity check: tied-≥2 among pairs sharing ANY read is
**54.5%**, nowhere near the ~100% that would mean flag-free co-creation ENTAILS the edge — the tie
requirement does independent work. And E_c overlaps E_r by only **4,084 / 35,258 = 11.6%**: read-conflict
and sequence-homology capture largely DIFFERENT relations.

**⛔ But it destroys the partition at EVERY depth.** Of 31,174 additive edges, **22,485 FUSE two families**
(501 families exist). Sweeping the shared-read floor:

| min shared reads | additive | FUSE families | fusion fraction |
|---|---|---|---|
| 2 | 31,174 | 22,485 | **72.1%** |
| 10 | 9,435 | 6,272 | 66.5% |
| 100 | 1,714 | 1,228 | **71.6%** |

⭐⭐**THE FUSION FRACTION IS INVARIANT TO DEPTH (~72% at 2 and at 100).** Depth does not discriminate; it
only scales both classes down together. Structurally identical to §-coverage-repair's "TP loss starts
BEFORE FP rejection ⟹ no constraint-satisfying threshold exists". ⟹**E_c cannot be an additive definition
tier at any depth.**

**Mechanism** — the same one §5b measured: **50.5% of random genomic region pairs align and 49.0% clear
identity ≥ 0.60**; repeat-driven cross-homology is universal. E_r survives it because of the COVERAGE
clause, which separates "shares an Alu" from "is a paralogue". A shared tied read IS the Alu case, and
read DEPTH is not a substitute for span coverage. ⟹**the shipped division stands: E_r DEFINES membership,
E_c SCOPES co-resolution** (and per the register, promoting E_c to "the principled family definition" was
the root cause of the advisor's "inconsistent approaches" complaint).

⚠**INSTRUMENT NOTE**: the first run returned **0 E_c edges from 2,499,153 records** — a `split("\t", 12)`
that lumps every tag past field 12 into one string, so the `AS:i:` scan saw none and every record was
skipped. It read as a clean null ("shared reads add nothing") and was pure artefact; only the same run
reporting **61,941 E_r edges over the identical reps** made it impossible. **A null needs a positive
control in the same run.**

## §6af — FLAGFREE_SITES CLOSED: recall real, but it fails the negative control (08-30)

§6ad measured flag-free at **NPIP 14/31 → 26/31, pure families 3 → 6**, amplification 2.0, density HIGHER
in every size band, phantom-copy test passing. Every one of those is a POSITIVE-panel or INTERNAL measure.
The missing evidence was a NEGATIVE control, and it reverses the recommendation.

**Clean negatives — canonical single-copy housekeeping genes** (`bench/negative_control/`, the panel's own
FAIL criterion: gene overlaps a copy in a family of ≥ 2 distinct loci):

| gene | floor 3 | floor 2 | FLAGFREE |
|---|---|---|---|
| ALDOA | PASS | PASS | PASS |
| **ATP5F1A** | PASS | PASS | **FAIL — GWFAM17, n=115** |
| RPL13A | PASS | PASS | PASS |

⚠**CANDIDATE COUNT: only 3 of the 30 housekeeping genes lie on this catalog's 3 contigs** — quote n=3, and
do not read 2/3 PASS as a pass. But **a canonical single-copy gene inside a 115-copy family is an
unambiguous false positive**, not a marginal call.

**Large-n corroboration** (2,207 uniquely-named protein-coding genes on the 3 contigs, in a ≥2-copy
family): floor 3 **12.78%** · floor 2 **16.63%** · flagfree **33.08%**. ⚠this negative class is
CONTAMINATED — a uniquely-named gene can be a real family member whose paralogues are LOC-named (RFPL 9/9,
GOLGA6 14/14), which is why the floor is already 12.78% at baseline. Read the DELTA, not the level: flag-free
adds **+20.3 points over floor 3**, i.e. 448 more genes acquire a family.

**⭐⭐ THE METHODOLOGICAL POINT: SIZE-MATCHED DENSITY DID NOT CATCH THIS.** Flag-free's families were denser
than baseline in every band (0.96/0.60/0.50/0.53) — and GWFAM17 (115 copies, swallowing a housekeeping
gene) is *part of* that. **A clique of shared repeats is DENSE.** Density is an internal-consistency
measure and cannot distinguish a tight true family from a tight repeat cluster; only an external negative
control can. ⟹**never close a recall win on internal structure alone.**

⛔**CLOSED: `RUSTLE_FLAGFREE_SITES` STAYS OPT-IN, NOT DEFAULT.** The 26/31 recall is real and the
primary-flag critique it answers is correct — the defect is that seeding from EVERY placement has no
analogue of E_r's coverage clause, so repeat-driven placements seed sites. ⟹**the reopenable question is a
site-level admission bar** (≥2 independent tied placements, or a coverage floor on the seeded site), NOT
the flag itself.

## §6ag — Overlapping multimappers as site evidence: the cheap win (08-30)

**The rule, in one sentence.** Don't decide where a copy exists from the primary flag — between
near-identical copies that flag is a coin toss — so let **every** alignment propose a site, and keep a
site when **more than one read's placements land on it**.

That is the whole idea. The advisor's objection all along was that seeding candidate sites from primary
alignments alone makes "which copies exist" (O1) depend on an answer to "where did this molecule come
from" (O2). Abstaining from that assignment is simpler than the machinery it replaces, not more complex.

**One guard is needed, and only one.** The two Pass-1 paths differ in KIND:

* **Spliced** reads must agree on an **exact intron chain** — already a real quorum, because a shared
  repeat does not manufacture a shared splice structure. Safe under flag-free, no change needed.
* **Unspliced** reads cluster on **span overlap alone** — no corroboration at all, so every repeat
  position accumulates secondary placements into one giant node.

⟹ **an unspliced site needs a primary; a spliced site does not.** That is the entire guard
(`RUSTLE_FLAGFREE_UNSPLICED=1` restores the unguarded behaviour for A/B). It preserves the intronless
pseudogene/retrocopy class that a blanket "spliced only" rule would have discarded.

**What it costs and what it buys** (fibroblast `npip3.bam`; ⚠ the 31 NPIP loci are a minimap2 projection
of HUMAN NPIP, not gorilla annotation):

| | NPIP loci | largest family | 1-exon share | single-copy controls |
|---|---|---|---|---|
| primary-seeded (shipped default) | 14/31 | 54 | 1.9% | 3/3 PASS |
| every placement, NO guard | 26/31 | 312 | 95.2% | **ATP5F1A FAIL, n=115** |
| **every placement + the guard** | **26/31** | **50** | 68.0% | **3/3 PASS** |

⭐**All of the recall (+12 loci), none of the pile-up**: 496,290 unspliced secondary-only placements
barred, largest family back to 50 vs the primary-seeded 54. The guard costs nothing on NPIP because the
sites carrying those loci are spliced, and spliced sites were never the ones barred.

⚠**Honest residuals.** The largest family is still **68% single-exon** (vs 1.9% primary-seeded), so
pile-up is reduced, not eliminated. On the large-n class (2,207 uniquely-named protein-coding genes) the
guarded arm sits **+6.7 points above the shipped default** (23.29% vs 16.63%) — the unguarded arm was
+16.5 — and part of that rise is legitimate, since flag-free is meant to find more real families and the
class is contaminated by LOC-named paralogues. ⚠**the clean control is n=3**: only ALDOA, ATP5F1A and
RPL13A lie on this catalog's contigs. ⟹**flag-free stays OPT-IN; the guard is intrinsic to it. A ship
decision needs the full-contig run where all 30 housekeeping genes are in scope.**

## §6ah — The genome-anchored repeat veto is a NULL on top of the site bar (08-30)

Ran `RUSTLE_ER_REPEAT_GATE=50` on top of flag-free + the asymmetric site bar (§6ag), to target the
residual unspliced pile-up the bar leaves (largest family still 68% single-exon).

| | families | copies | NPIP | largest | 1-exon | housekeeping |
|---|---|---|---|---|---|---|
| floor 2 (shipped default) | 121 | 678 | 14/31 | 54 | 1.9% | 3/3 PASS |
| + flagfree + bar | 291 | 1,584 | 26/31 | 50 | 68.0% | 3/3 PASS |
| + flagfree + bar + repeat veto | 285 | 1,557 | **26/31** | **50** | **68.0%** | **3/3 PASS** |

The veto flags **92 of 4,507 reps** as genome-anchored hubs and drops **188 of 6,554 edges (2.9%)**;
every endpoint is unchanged. ⟹**the bar had already removed what the veto would have caught. The residual
68% single-exon share is NOT repeat-hub driven.**

**⛔ THIS ALSO UNDERCUTS THE ML EXPERIMENT IT WAS RUN TO ENABLE.** The plan was: label ← genome-anchored
multiplicity (from the REFERENCE, independent of the RNA features), features ← RNA, task ← predict the DNA
verdict from RNA evidence. That label is the only non-circular one available (the catalogue is circular;
cross-substrate replication is confounded by expression — the trap behind the retracted "not expressed").
Measured, **it identifies 92 reps whose removal changes no endpoint** ⟹ a model trained on it would be
redundant with a rule already in the tree. ⚠**the binding constraint on ML here is LABELS, not models.**

⭐**What would break it**: an assembly-derived per-copy truth set of the kind that gave GOLGA6L7 = 6 on
`_mat` / 7 on `_pat` — DNA-derived, per-copy, independent of the RNA evidence. Bounded work for a few dozen
families, and useful to the thesis whether or not a model is ever trained on it.


## §6ai — k+1 MEC for reference-absent copies: NO-GO, four independent kills (08-30)

Four-lens design panel + a synthesizer that RAN the objective itself on the real dumped PSV matrix
(79,176 rows / 101 copies / 36,411 columns / 12 families), then adversarial verification. Verdict
**NO-GO**, on four grounds each independently sufficient.

**1. THE OBJECTIVE IS DEGENERATE UNDER MISSING DATA.** "Cost = cells disagreeing with their group's
consensus" is ill-posed when a group observes no cell in a column: the optimiser minimises by SHRINKING
THE SCORED CELL SET, and best-of-R selects exactly that. Measured: GWFAM115 at k=6 returns **MEC = 0
with 15.3% of observed cells scored**, against a genuine k=5 fit of MEC = 1,076 at 100% scored — an
apparent Δ of +1,076 that is **100% artifact**. GWFAM111 k=5 → MEC = 0 at 22.5% scored; GWFAM83 → 0 at
both k=5 and k=6. Reseeding empty clusters does not fix it (the attractor is PARTIAL, not empty). The
only repair is a hard min-group-size floor — and it is **load-bearing**: at n_min=5, GWFAM83's k=6
solution is [348,15,10,5,5,5], three groups sitting exactly ON the floor. ⟹ the answer becomes a
function of a tuned constant, which is the "no arbitrary thresholds" standard the advisor set.

**2. THE CLEAN NEGATIVE DROPS HARDEST, AND THERE IS NO ELBOW.** With the floor and 100% of cells
scored: **GWFAM83 — the ONLY `reference_resolved` family, `collapsed_copies = 0`, i.e. the one family
with NO missing copy — has the LARGEST k→k+1 relative drop (0.576)**, against 0.419 (GWFAM111) and
0.322 (GWFAM115), both collapsed. And k+1→k+2 is comparable or larger in 2/3 ⟹ smooth monotone decay,
no k is selectable.

**3. THE k = n_copies CONTROL DOES NOT PASS.** ARI of the MEC partition against the shipped assignment
on evidenced reads: GWFAM115 **1.0000**, GWFAM111 **0.9841**, GWFAM83 **0.7621**, GWFAM104 **0.2718** —
2 of 4 below the pre-registered 0.80 ⟹ stop before k+1.

**4. THE SUBSTRATE HOLDS ≈0.2 EXPECTED TRUE POSITIVES.** `o3_missing_copy_evidence.md` §4.2: M ≤ 8.5
expressed reference-absent copies genome-wide, point estimate 1.80, CI [0.09, 8.54], over 915 catalog
copies ⟹ this 101-copy batch expects **0.20**. And the truth set cannot substitute: `o3_excise` deletes
one copy of **162 TWO-copy families**, so every positive is the k=1→2 regime the incumbent already
handles **40/40** — **the panel contains zero instances of the object under test**, and on it depth
alone already gives AUC 0.8034 (p=0.0005) with no MEC at all.

**⛔⛔ TWO OF MY OWN CLAIMS CORRECTED BY THE PANEL.**
* **"Read sparsity caps any phasing-style method here" — WRONG.** Measured **26,177/26,183 informative
  reads span ≥2 columns**, pairwise read-read overlap **28–1,044 columns**. Linkage is not the blocker.
* **"The EM (0/3,081) is the soft version of MEC, so k=n is a control" — WRONG.**
  `em_copy_assign::m_step` re-estimates only `pi` while `logl` is an **immutable input built from the
  reference `CopyProfiles`** ⟹ the EM structurally CANNOT form a novel consensus, so register:576 does
  not bound MEC in either direction. k=n is a control on different grounds — and it fails.

**⭐ AND THE CONCEPT ALREADY SHIPS.** `denovo_pipeline.rs:1223-1235` calls `split_locus_copies` with NO
reference input and emits every non-host cluster as a reference-absent candidate — **289 such candidates
already exist across these 12 families, 0 admitted**. MEC's gain over it is MECHANICAL (error- and
span-tolerance vs exact-vector equality), and mechanical gain is not detection gain: the incumbent's
multi-copy failure is a **SPECIFICITY** failure (identical statistics with and without an absent copy),
which MEC does not address. Also: the PSV column space is **reference-conditioned**
(`copy_assign_pipeline.rs:441-442` inserts a column only where two COPY sequences differ), so a k+1
consensus is necessarily a mosaic over reference alleles — an object `copy_assign.rs:1449` already models
as `copy_conversions`, and register:520 already refutes identifiability at K≥3 via cross-copy
recombination (our k = 5–19).

⟹**MEC is closed for O3 detection.** What survives is a narrower, gated POWER question (can a
read-inferred group be told from the family's own reference copies at all, on a within-family matched
pair) — **never a detector, and no TPR/FPR/detection rate may be quoted from it.**

## §6aj — The advisor's two PSV objections, measured (08-31)

Four investigative lenses + a synthesizer + three adversarial refuters, all recomputing from the real
12-family dump (36,410 columns / 79,175 read-rows / 101 copies). **All three refuters fired**; what
follows separates what survived independent recomputation from what did not.

### Objection 1 — "you are not doing variant calling; give me a VCF"

**✅ SURVIVES — the frame is partly a category error.** `discover_psvs`
(`copy_assign_pipeline.rs:385-441`) star-aligns every copy to copy 0 and inserts a column only where
`is_acgt(ca) && is_acgt(cb) && ca != cb`. **A PSV is a FIXED DIFFERENCE BETWEEN ASSEMBLED COPY
SEQUENCES**, not a call from a pileup. There is no DP/AF/QUAL to compute *for the column*.

**✅ SURVIVES — and he is factually wrong twice.**
1. **A read pileup gate ALREADY runs by default**: `read_supported_columns` (`:652`, wired `:1470-74`,
   disable with `RUSTLE_PSV_READFILTER=0`) keeps a column only if reads show ≥2 alleles at ≥2 reads each.
   Counterfactual through the real binary: **36,414 candidates → 36,410 shipped = 4 dropped (0.011%)**.
   Real, default-on, near-inert here.
2. **Sequencing error IS modelled** — per-base Phred QVs are read from the BAM (`fill_psv_obs :751`) and
   enter the likelihood and the per-read certificate, phred→error clamped to [1e-4, 0.25].

**✅ SURVIVES — read support, confirmed by an independent recompute.** **0 / 36,410 columns are
read-unobserved** (min depth 11, median 677) — and not by filter artifact: the predicate at `:669`
explicitly passes columns with cov<4. **36,170 / 36,410 = 0.9934** have every copy-claimed allele seen in
≥2 reads. ⚠quote the third-allele rate **stratified**: 105,554/12,439,427 = **0.008485** at 2-allele
columns, NEVER the pooled 0.005543 (the 4-allele stratum is tautologically 0 over 24.7% of that
denominator).

**✅ SURVIVES — the RNA-editing confound is REFUTED with a mirror control.** Copy alleles are stored in
transcription orientation (verified 70/70 direct, 50/50 complement vs GGO.fasta), so A-to-I would appear
uniformly as A→G. Measured: **A→G 15,450/3,234,432 = 0.004777 vs G→A 20,607/1,953,346 = 0.010550, ratio
0.453 — A→G is DEPLETED, not enriched.** Reference-level A/G columns 3,102 vs the C/T strand-complement
control 2,959 (ratio 1.048, n.s.); Ti/Tv 0.657.

**⛔ THE PROPOSED DELIVERABLE (build a proper VCF) IS REFUTED.** `genome_pos` **IS NOT A LOCUS**:
`copy_assign_pipeline.rs:1495-1500` sets it via `col_canon[col].get_or_insert(g)` — the FIRST copy
carrying the column — an internal canonicalization device so switch breakpoints share one frame. The other
copies sit **33 kb to 81 Mb away** (GWFAM55 max 81,042,216 bp). ⟹ every ALT a VCF emitted would be a base
from sequence **megabases from the POS the row asserts**. A VCF would make the misrepresentation WORSE.
⚠the existing `bench/igv_tracks.py::write_psv_vcf` is independently broken: `chrom_of` (`:92-96`) scans
regions by coordinate ignoring contig ⟹ **6,039/36,410 = 0.1659 wrong CHROM**, and REF is the
transcription-strand base ⟹ **19,034/36,410 = 0.5228 contradict the genome** and fail
`bcftools norm --check-ref`.

### Objection 2 — "similar copies cannot carry many PSVs; usually there is no difference"

**✅ SURVIVES — the strong form is REFUTED by direct count.** **0 / 490 within-family copy pairs have zero
distinguishing columns.** Minimum 4; minimum **11** over the 469 pairs with ≥200 mutually-ungapped columns.
Not a majority — not one pair.

**✅ SURVIVES — but his prediction is a TAUTOLOGY, not an empirical claim.** In the panel frame,
PSV-per-kb **= 1000 × (1 − identity) to machine precision** (max deviation 0.0000000000 over 490 pairs;
Spearman exactly **−1.0000** pooled and within all 12 families). Density and divergence are the same
number, so the prediction cannot fail. ⚠**presenting a correlation here as confirmation would be
presenting arithmetic as evidence.**

**⛔ EVERY IDENTITY-REGIME CLAIM IS REFUTED — THE AXIS IS NOT AN IDENTITY.** "Panel-frame identity" has the
PSV panel as its own denominator — **a denominator conditioned on the prediction**, this project's named
trap. Against the pipeline's OWN identity metric (`-k11 -w5`, `1-de`) it correlates at **rho = +0.109,
n = 438, p = 0.022, r² = 0.012**. ⟹ **retract "our families are ~75%-identical paralogues" and "only
1.3% of pairs sit in the too-similar regime"** — both are artifacts of that axis.

### ⭐⭐ A REAL BUG SURFACED BY THE THIRD REFUTER

**1,587 read names are genotyped in ≥2 family matrices** (79,175 rows over 77,372 distinct names).
**517 carry `status='assigned'` in two or more families simultaneously**, giving **518 double-counted
assignment rows** (53,715 assigned rows over 53,197 distinct assigned names). Of 519 assigned-copy pairs,
**309 overlap genomically but 210 are DISJOINT** — median separation **16,341 bp**, max **37,727,971 bp**.
**A single molecule cannot originate from two disjoint loci, so at least one member of each of those 210
is wrong.** Not previously known; not handled anywhere in the pipeline.

## §6ak — Cross-family double assignment: fixed, and INDEPENDENTLY verified (08-31)

`RUSTLE_XFAM_RECONCILE` (off = default / report / abstain), committed a8079cf. The workflow's verify
phase crashed the terminal twice (3 xhigh agents each re-running `cargo test --release` + `copy_assign`;
one `copy_assign` alone is **9.6 GB RSS**, so three concurrently exceed the 25 GB box — an authoring
error, heavy work must never sit inside a `parallel()` block). Verified SERIALLY instead, by hand.

**⭐ MY ORIGINAL FRAMING WAS THE WRONG CUT.** I posed the bug as "210 assignments at DISJOINT loci". The
right axis is whether two assigned rows resolve to ONE alignment record or TWO. Measured stratification
of all 519 contested pairs: **309 shared_locus** (one placement, two family labels — an O1 naming defect
in an O2 costume) · **99 readthrough_span** (one record spanning both copies — biology) · **111
cross_family_contradiction** (two independent records — genuinely impossible). ⚠**an interval-geometry
gate fires on all 519 and would abstain 415 correct single-placement assignments.**

**Independent verification (my own recount from the three arms, not the implementer's numbers):**

| arm | rows | assigned | reads in ≥2 fams | disjoint copy pairs |
|---|---|---|---|---|
| off | 79,175 | 53,715 | 517 | 210 |
| report | 79,175 | 53,715 | 517 | 210 |
| abstain | 79,175 | **53,494** | **407** | **98** |

* **OFF is byte-identical to the pre-fix dump** — `assignments.tsv` fe563640, `quant.tsv` 219f9007,
  `families.tsv` 033b315a, all SAME (my md5s, not the implementer's).
* **Row count is INVARIANT at 79,175 across all three arms** ⟹ the fix changes STATUS, it does not hide
  rows. 221 demotions, all `assigned → ambiguous`.
* **111/111 genuine contradictions demoted, 0 remain.** The 112th demotion is a readthrough pair whose
  molecule is in ≥3 families and also carries a contradiction — the rule working, not a leak.
* **The 98 residual disjoint-by-geometry pairs are 98/98 `readthrough_span` with `same_record = true`** —
  ONE record spanning two copies, median separation 16,341 bp, max 76,199 bp. Correctly kept.
  ⟹**geometry alone cannot separate a two-record contradiction from a one-record readthrough**, which is
  exactly why the record axis is the right one.

**Discipline checks, all pass:** `xfam_reconcile` IS in the params certificate (`copy_assign.rs:2456` —
the M2 defect avoided); an unknown mode returns `Err` rather than silently running as off; the tests pin
the CLAUSE (`abstain_demotes_only_the_two_record_contradiction`) and the VALUE
(`abstain_row_counts_are_exactly_two_demotions`) separately. Suite **824 passed / 0 failed / 11 ignored**
(baseline 815; delta = exactly the 9 new tests, no existing test edited).

⭐**A SECOND BUG WAS FOUND AND FIXED EN ROUTE (f7ed833): `bench/mechanism/byte_identity_gate.sh` pointed
at a `BIN` path that did not exist**, so `check` ran no binary and md5'd stale files — the gate had been
passing VACUOUSLY. Same class as the blind airtight panel and the missing trace hook: an instrument
reporting success while structurally unable to fail.

⚠**default stays OFF** pending a genome-wide arm; `quant.tsv` and `families.tsv` are byte-identical
across all three modes, so abundance is unmoved by the demotion.

## §6al — VACUOUS INSTRUMENTS: a defect class, one confirmed, scope unassessed (08-31)

Three instruments were found this week that **reported success while structurally unable to fail**. That
is a class, not three coincidences, and it is the most dangerous defect shape in this project because a
vacuous PASS is indistinguishable from a real one in the record.

**CONFIRMED, all three:**
1. `bench/mechanism/byte_identity_gate.sh` — `BIN` pointed at a path that did not exist, so `check` ran
   NO binary and md5'd stale files. **The gate had been passing vacuously.** Fixed (f7ed833): default
   corrected to the real `CARGO_TARGET_DIR`, hard abort on a missing binary, second corpus added.
2. `RUSTLE_DEBUG_LOCUS` existed only in `detect_conflict_catalog_genome_wide_xchrom`, but
   `--homology-primary` runs `detect_homology_catalog_genome_wide`. A full trace run produced **0 dbg
   lines** — the instrument was absent from the path under test. Fixed (fef1317).
3. `bench/family_def_airtight_panel.py` — computes edges from ANNOTATION-derived cDNA homology and never
   invokes the pipeline ⟹ **structurally blind to any collapse/representative flag**. Both arms return
   "identical", which reads as a pass. Not fixed; documented in §6aa.

⚠ A fourth, in my own analysis code: `split("\t", 12)` lumped every SAM tag past field 12 into one
string, so an `AS:i:` scan found none and **every record was skipped — 0 E_c edges from 2,499,153
records**. It read as a clean null ("shared reads add nothing"). Caught only because the SAME run
reported 61,941 E_r edges over identical reps, making the zero impossible. ⟹**A NULL NEEDS A POSITIVE
CONTROL IN THE SAME RUN.**

**SCOPE — FOUND, NOT ASSESSED. Do not quote these as defect counts.**
* **`target/release/` DOES NOT EXIST in the repo** (builds go to `CARGO_TARGET_DIR` on linuxdisk), yet
  **36 bench files invoke `target/release/<binary>`**, 8 of them named gate/check/verify/validate.
* ⛔**CORRECTED BY §6am — THIS CLAIM IS WRONG.** `/home/juanfra/winloci_scratch` is a **SYMLINK that
  RESOLVES**, so these scripts RUN rather than aborting, which is WORSE. Original text kept for the
  record: ~~349 files reference `/home/juanfra/winloci_scratch`, a path that no longer exists (data moved to
  `/mnt/linuxdisk/.../_from_wsl/winloci_scratch`); **37 of those are gates or validators**.~~
* **Only 25 of 64 bench shell scripts use `set -e`**, so a failed command does not abort — the script
  continues to its comparison and can still print a verdict. That is exactly the `byte_identity_gate`
  mechanism.
* Sampled 3 (`core_gate_gw.sh`, `gate_onoff_gw.sh`, `soto/salvage_validation.sh`): all three combine a
  stale BAM path, the nonexistent `target/release` binary, and no `set -e`.

⚠⚠**WHAT IS NOT ESTABLISHED: how many of the 37 return a FALSE PASS versus dying loudly.** Failing loudly
is safe; passing vacuously invalidates whatever it certified. **One** is confirmed (`byte_identity_gate`).
Extrapolating from a sample of 3 would repeat the "210 disjoint" over-reach corrected in §6ak. The audit
is scoped but UNRUN.

⟹**STANDING RULE: an instrument must be shown capable of FAILING before its PASS is evidence.** Check the
candidate count (§4t's human panel offered 0 qualifying pairs), check the binary exists, check a null has
a positive control in the same run.

## §6am — The vacuous-gate audit: 10 confirmed false PASSes, and two of my premises corrected (08-31)

40 verdict-emitting bench instruments classified STATICALLY (nothing executed: cat/grep/ls/`test -e`/
`bash -n` only, on a box where one pipeline binary is 9.6 GB and parallel builds crashed the terminal
twice). **10 FALSE_PASS · 19 FAILS_LOUD · 11 OK · 0 OVERTURNED** — every false-pass candidate was upheld
under independent re-derivation by a second agent instructed to default to overturning it.

### ⛔ TWO OF MY BRIEFING PREMISES WERE WRONG

1. **`/home/juanfra/winloci_scratch` IS NOT GONE — it is a SYMLINK** (Aug 13) to
   `/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch`. ⟹**retract "349 files on a stale path"**;
   the path RESOLVES. **This makes things worse, not better**: those scripts RUN rather than aborting, so
   the failure that would have saved them never fires. Three of the ten false PASSes
   (`validate_exon_sum`, `soto_family_validate`, `rna_reframe_validate`) exist *because* of this.
2. **`target/` does not exist AT ALL**, not merely `target/release/`. **33 executable bench scripts
   invoke a binary beneath it; only 5 were in the audited 40 ⟹ 28 UNAUDITED scripts carry the canonical
   `byte_identity_gate` shape.**
⚠count drift: measured **26/65** shells with `set -e` (register says 25/64) and **34** files invoking
`target/release` (register says 36). Same conclusion — **re-derive before quoting either as exact**.

### The ten, worst first

* **`soto/gate_sensitivity_sweep.sh` — HIGHEST RISK, AND DESTRUCTIVE.** `set -u` only. `GWCAT` resolves
  via `soto/_detect_unit.sh:6` to the nonexistent repo-local binary. Line 16 **TRUNCATES the committed
  result** `$CACHE/gate_sensitivity.tsv` (342 B, Jul 23 17:05) and `run_config:37` `rm -f`s the cached
  per-chrom outputs — **both BEFORE any binary is attempted** — then scores whatever stale file survives
  and appends a verdict row. Quoted in `bench/soto/gate_robustness.md`. **It deletes the evidence and
  writes a fabricated replacement in the same run.**
* **`soto/salvage_validation.sh`** — the `byte_identity_gate` bug reproduced **three times in one file**,
  including a byte-identity gate of its own: `BIN` missing ⟹ `timeout 600 "$BIN" >/dev/null 2>&1` returns
  127, nothing aborts, and it md5-compares leftover outputs and prints PASS. Cited `o1_ledger.md:5289`.
* **`family_def_airtight_panel.py`** — already known (§6aa); **still runs to completion today** because
  its inputs are live. Never invokes the pipeline ⟹ blind to every `RUSTLE_*` flag it was used to judge;
  its negative criteria pass **by default on an empty evidence set**. Quoted in the register and ledger.
* **`gate_onoff_gw.sh`** — a silent exec failure leaves a PRIOR run's GTFs in place, so the keyed diff
  reports **lost=0 gained=0**: the strongest possible claim ("flipping the gate ON loses no real
  transcript") produced by **zero execution**.
* **`family_def_newbam_validate.py`** — ⭐**the two arms are the SAME INODE**: `GGO.bam` is a symlink to
  `GGO_mm.bam`, and `:27/:28` set `OLD=GGO.bam`, `NEW=GGO_mm.bam`. The comparison is against itself and
  cannot fail. ⚠**BUT the recorded table (`FAMILY_DEF.md:526-530`, secondaries 4,818→313,883) PREDATES the
  symlink and IS a real measurement — re-running today would silently null a TRUE result. DO NOT
  re-freeze that table from this script.**
* **`soto_family_validate.py`** — **MAD over an all-zero famCN vector is exactly 0.0**, so a map-back
  returning nothing scores as "100% famCN-CONSISTENT". Every span is cut from the same genome used as the
  minimap2 target, so a zero-hit copy is *proof of failure* and is scored as a pass. `cn_ok = mad < 1.0`
  is unfalsifiable. Quoted as **"famCN-CONSISTENT 67/82 = 81.7%"**.
* **`o4_diploid_validate.py`** — both `.mmi` inputs verified absent ⟹ the COPY branch of `classify()` is
  unreachable ⟹ **"COPIES confirmed: 0" is true by construction**; the verdict and JSON are still written.
* **`validate_exon_sum.py`** · **`rna_reframe_validate.py`** · **`psv_sizing/psv_phase_validate.py`** —
  unchecked subprocesses / `max(1,…)` denominators / silent `if exists else []` comparators, each turning
  an empty computation into a printed rate.

### ⭐ One recorded number partially rescued by an internal positive control

`VALIDATION_AND_STATUS.md` §F2's **81.7%** is *probably real*: the same run reports **DNA shared-exon
25/82**, which is NONZERO and requires real PAF hits, so the aligner demonstrably worked at record time.
No such internal control exists for `gate_robustness.md`'s sweep rows or `o1_ledger.md:5289`'s
salvage/gate_onoff citations — **their provenance is UNRESOLVED.**

⚠**SCOPE: 40 of ~657 bench scripts = ~6%**, and the selection was mine, not random. **Absence from this
report is not evidence about a script.** Static reading establishes *prospective capability to fail*; it
does NOT establish that any recorded number was produced vacuously. That needs a run, per script.

## §6an — The ten false-PASS gates repaired; every one retains residual holes (09-01)

10 scripts / **12 files, 721 insertions, 89 deletions**; all 12 parse (`bash -n` / `py_compile`, my own
check). **Zero data or doc files touched** — `bench/soto/gate_sensitivity.tsv` is byte-identical at
md5 `18156526`, so the committed evidence quoted in `gate_robustness.md` is intact. Nothing was executed.

Every repair follows `byte_identity_gate.sh` (f7ed833): a real overridable default path, a **hard abort
with a nonzero exit BEFORE any work or truncation**, and a stated reason. Verified independently: all 10
**can now fail** before printing a verdict, and all 10 **preserve behaviour** on a fully-provisioned run.

**Notable repairs.** `gate_sensitivity_sweep.sh` now aborts above the truncation, **copies the previous
$SUM to `.prev.<timestamp>` before rewriting**, stamps binary+commit+REASON into the new file, and adds a
freshness backstop requiring all NBED+1 units to produce a `copies.tsv` newer than a pre-run stamp —
because `recompute_perchrom.sh` exits 0 even when every unit dies. `family_def_newbam_validate.py` now
**aborts when its two arms resolve to the same inode** (the symlink defect). ⚠the agent deliberately did
NOT change `recompute_perchrom.sh`'s exit contract (other callers depend on it) and did NOT touch
`soto_cache_score.py`'s 4-leg `hit()` — flagged as science, not a guard.

**⚠⚠ RESIDUAL HOLES — THE REPAIRS CLOSED THE HEADLINE MECHANISM, NOT EVERY ROUTE.** Each was found by an
independent verifier reading the repaired file; none is fixed:

* **`family_def_airtight_panel.py`** — the positive control is **panel-wide**, so **3 of the 7
  counterexamples still pass on ABSENCE of evidence**: all 9 pairs of the name-coincidence cases (NDUFS,
  NDUFA, COX) have ZERO Hd records. The under-resolution abort also fires only AFTER every verdict is
  printed. ⟹**its PASS is still not citable as a flag verdict.**
* **`salvage_validation.sh`** — Part A **prints the three md5s but never COMPARES them**: if OFF run1 ≠
  OFF run2, or ON == OFF (salvage inert), nothing aborts and it proceeds to "DONE".
* **`soto_family_validate.py`** — the repair closed famCN=0 but **NOT famCN=1**: every query span is cut
  from the minimap2 TARGET, so a self-hit is guaranteed and `min(famcn) > 0` is **satisfied by
  construction**. `SOTO_MIN_HIT_FRAC=0` also disables the evidence floor by default.
* **`gate_onoff_gw.sh`** — **the gate is never required to have FIRED**: `gated_merges` is recorded but
  never checked, so if 0 merges are gated genome-wide the ON and OFF arms are the same computation and
  "lost=0 gained=0" is again vacuous. Also truncates $SUM on every run while `continue`-ing past latched
  contigs, so a resumed run scores only what it ran this time.
* **`o4_diploid_validate.py`** — the positive control and empty-block guard are **floors of ONE, not
  proportions**; the reuse stamp binds the query but **never the targets**, so PAFs built against the
  wrong targets are accepted; and a **failed determinism check never aborts**.
* **`validate_exon_sum.py`** — per-family aligner vacuity survives (the guard is genome-wide only), so a
  family whose minimap2 exits 0 with empty stdout still scores purity 1/n and **strengthens** the verdict.
* **`rna_reframe_validate.py`** — `except KeyError: continue` on the `de` tag is the largest hole: a
  tag-stripped BAM silently empties the computation.
* **`psv_phase_validate.py`** — guards stop at the inputs and never reach the scoring loop; every locus
  can be discarded by the size filter and a verdict is still printed.

⟹**a repaired gate is capable of failing, which is strictly better than before, but "PASS" from any of
these still needs its evidence count read.** ⚠**28 dead-binary scripts remain UNAUDITED** (§6am).

## §6ao — Provenance of the numbers the false-PASS gates produced (09-01)

§6am left two artifacts' provenance UNRESOLVED. Both now resolve, and the concern is much narrower than
I stated.

* **`salvage_validation.sh` and `gate_onoff_gw.sh` have NO quoted numbers anywhere.** The audit cited
  "o1_ledger.md:5289", but that line is **my own §6al sampled-3 sentence** — not a place their output is
  used as evidence. Grep confirms the only other "salvage" hits in the register are a different subject
  (mis-chain read salvage; a Gotoh traceback). ⟹**nothing rests on those two scripts.**
* **`gate_robustness.md`'s sweep numbers STAND, on an internal positive control the audit missed.** The
  four configs report **four DISTINCT copy counts — 467 / 376 / 483 / 513 (a 137-copy spread)** — each
  moving as its gate predicts. **A dead binary scores the SAME stale files for every config**, so a
  spread that size requires a genuine per-config rebuild. Same rescue shape as §F2's 81.7% (whose control
  was a nonzero DNA shared-exon 25/82). A provenance note is now IN `gate_robustness.md` itself, beside
  the numbers, not only here.
  ⚠**the recall column is the weak half**: `soto_cache_score.py`'s `hit()` unions FOUR legs, three of
  them static files independent of the rebuild, so the **81.2–82.3% band is insensitive by construction**
  and is not evidence of stability. Quote the copy swing; discount the recall band.

⛔**ALSO CORRECTED IN PLACE: §6al's "349 files reference a path that no longer exists" is WRONG** and is
now struck through there with a pointer to §6am. `/home/juanfra/winloci_scratch` is a **symlink that
resolves**. Left visible rather than deleted, because the wrong version was committed and quoted.

⟹**REMAINING TRULY OPEN: the 28 dead-binary scripts (§6am), and the residual holes in all ten repaired
gates (§6an).** Neither is a provenance question about a recorded number.

## §6ap — Reads spanning two near-identical copies: keep current for O2, and a LIVE O1 defect found (09-01)

Four lenses + a synthesizer + three adversarial refuters, all re-deriving from the real dumps.
**All three refuters fired.** What follows separates what survived independent recomputation.

### ✅ O2: KEEP CURRENT — and it is now PROVEN, not assumed

A record spanning two disjoint copies is assigned in BOTH families, and that is **not** double-counting.
`fill_psv_obs` writes an observation only inside an M/=/X block and `allele_at` maps genomic→read offset
**injectively**, so on disjoint copies the two families interrogate **disjoint genomic positions and
therefore disjoint READ BASES**. Measured on the dominant pair: **294 + 1,512 = 1,806 observed columns at
disjoint intervals = 1,806 DISTINCT read bases of ~3,255 aligned.** ⟹**each copy's claim is independently
earned on bases the other never saw.** One molecule, two segments, two independent evidence sets.

**⭐ ABUNDANCE PROVABLY CANNOT MOVE**: `quant.tsv` md5 **219f9007fe2cfa7e39e14ca00da62207 across ALL FOUR
arms** (pre-reconcile baseline, off, report, abstain). No reconcile decision — existing or proposed —
moves a single abundance figure; only the `status` column can change. ⟹**this was never an abundance bug.**

**⛔ AND IT NEEDS NO NEW MACHINERY.** The population is **99 rows / 77,372 distinct molecules = 0.128%**,
and **93 of the 99 are ONE copy pair** (GWFAM118:0 × GWFAM111:3, separation 16,341 bp) — essentially one
stub. A 3-flag / 3-arm proposal was **refuted on necessity**: the honest deliverable is keep-current, add
diagnostic columns, write the paragraph. ⚠also `same_record` is **logically incapable of failing on 94%
of the stratum**, so its agreement is not evidence.

### ⛔⛔ THE REAL FINDING: THE DAZ FAILURE MODE IS ALREADY LIVE IN THE SHIPPED CATALOG

The claim that "the shipped stance is already correct" is **FALSE**. Running the in-tree
`catalog_overlaps` logic over `arm_f2/cat.copies.tsv` (678 copies, 121 families):
* **98 `SharedAcrossFamilies` overlaps** — the enum's own doc says this "is ALWAYS a defect: a readthrough
  transcript in one family spans loci that are separate copies of another. Observed at GSTM."
* **71 copies are STRICTLY ENGULFED by a copy of a DIFFERENT family, touching 34/121 families.**
* In `dump/e.nodes.tsv`, **441/3,598 nodes engulf ≥1 other node (910 engulfments), and 115 of those carry
  ≥1 E_r edge.**
⟹**this is not a hypothetical guarded against; it is shipping.** ⚠the proposed "extension" would have made
node creation strictly worse.

**The mechanism, located precisely:** `is_chimeric_bridge` (`denovo_assemble.rs:1550`) — "two spans either
intersect or they do not", threshold-free — is applied at the GATE but **NOT at collapse**.
`collapse_parent` phase 1 (`family_detect.rs:613-624`) unions any two transcripts sharing an exact
`(chrom, donor, acceptor)` **with no chimera test**, so a spanning object that clears the gate on its own
reads still merges the two loci — the GSTM/DAZ mechanism, one stage later.

### ⚠ THE USER'S LITERAL CASE IS INVISIBLE TO THE DETECTOR

Near-identical **and** spannable is **0/3,141 E_r edges at identity ≥0.98**, but at ≥0.90 there are **27
pairs within 100 kb and 4 within 10 kb** — and **two are WITHIN a family in this batch**: GWFAM111 c0~c1
(identity **0.9400**, gap **8,166 bp**) and GWFAM118 c4~c5 (**0.9603**, **41,885 bp**).
⚠**`copy_assign.rs:659` skips same-family pairs, so the xfam detector structurally CANNOT see them**, and
`best_overlap_copy` picks ONE frame per record ⟹ **the second segment's bases are silently discarded.**
That is the real gap the question points at — not the 99.

### ⚠ WHAT WE CANNOT TELL FROM WHAT IS DUMPED

**Whether the ~13 kb gap is a REAL splice junction or an aligner-invented N gap.** No per-read CIGAR,
intron chain, junction coordinate or read strand exists in any dump. **That is exactly the test separating
a genuine readthrough from a chimera or mis-chain, and it is UNRUN** — so "biology" in the stratum comment
is an ASSUMPTION. ⭐**The settling statistic is breakpoint concordance: a real readthrough reuses ONE exact
donor/acceptor across the 93 records; a mis-chain scatters them.**
⚠also unresolved: whether GWFAM111 copy3 is a distinct copy at all, or the alternative terminal exon of
the GWFAM118 copy0 gene 21 Mb away (single-exon, 2,096 bp, 11 reads, `stub=true`).

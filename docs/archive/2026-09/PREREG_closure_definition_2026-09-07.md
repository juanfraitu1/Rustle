# PREREG — the family as a CLOSURE (fixed point of "shares a core with half the family"), 2026-09-07

**Written before the prototype.** md5 in `mcl_ann/adj/closure/PREREG.md5`.

## The definition under test
> **A multi-copy gene family is a set of genomic intervals that all carry copies of the same segment, where
> that segment is what at least half of them share.**

Operationally a closure, iterated to a fixed point:
1. **Seed** — any set of loci that align to each other. Here: **connected groups** of the annotation's
   all-vs-all links (identity >= 0.70, coverage of the longer >= 0.30, >= 300 bp). **No MCL.**
2. **Core** — for each current member, the positions linked by duplication pairs to **at least half of the
   other current members**. Member iff core >= half its span, or >= half the set's median core (shipped rule).
3. **Extend** — align the members' core hulls to the genome; every interval carrying the core that is not
   already a member is admitted if it passes step 2 against the current members.
4. Repeat 2-3 until the member set stops changing.

This replaces two named objects with one: MCL (inflation 2.8, prune 1e-9) disappears, and the annotation
becomes a seed rather than the node set. Preparatory measurements that motivate it, both from 2026-09-07:
connected groups + the core rule reproduce MCL + the core rule on NPIP (59 candidates pruned to 35 records
over 25 loci, precision 35/35); and leave-one-out core projection recovered a deleted member in **25/25**
folds, admitting nothing else.

## Substrates and truth
Gorilla three contigs (`allgenes.asm20.paf`, `GGO_sedef_final.bed`, `npip3_contigs.fa`), truth = the 26 LCR16a
loci; and human chr16+chr18 (`npip_hsa/hsa.paf`, `HSA_sedef_pairs.bed`, `chm13v2.0.fa`), truth = the 26
description-complete NPIP records. Human and gorilla never pooled. `lcr16a.bed` is a core projection of the
shipped catalog (register 727): it is **valid for membership** questions and **invalid for boundary** ones,
so no boundary claim is scored on it.

## Predictions
| # | prediction |
|---|---|
| **P1** | the closure **converges in <= 5 iterations** on both substrates |
| **P2** | at the fixed point, precision (members that are truth) **>= 0.95** on each substrate |
| **P3** | sensitivity **>= 25/26** (gorilla) and **>= 24/26** (human) |
| **P4** | **seed-independence**: started from a random half of the seed set, and from a single locus's direct neighbours, the closure reaches the **same** member set as from the full connected group |
| **P5** | the iteration is **not monotone** - at least one record leaves and later rejoins (observed once already on NPIP), so a convergence proof cannot rest on a decreasing set |
| **P6** | **exonic-core refinement** (require the shared segment to carry read-supported exonic bases) raises chrY DAZ specificity above the amplicon result of 4/15 |

## Interpretation fixed in advance
- P1+P2+P4 holding ⟹ the closure is the definition to write in the thesis and MCL is dropped to a footnote.
- P4 failing ⟹ the fixed point is seed-dependent; the definition is then **not** well posed as stated and must
  name its seed, which is a materially weaker claim and will be reported as such.
- P6 failing ⟹ the amplicon limit stands and is stated as a limit, not patched.

---

## AMENDMENT 1 — the MINIMAL closed set (2026-09-07, written after P4 failed on chrY, before the runs)

Register 735: the closure is **not unique**. On chrY the AZFc amplicon (15 units) and {DAZ1-4} are both closed;
the group seed reaches the first, a DAZ seed the second. The candidate repair, now under test:

> **The family of a locus is the SMALLEST closed set containing it.**

Operationally: start the same closure from a **single locus** and iterate. One sentence, still no MCL, and it
makes the definition single-valued if the least fixed point reachable from a locus is the same wherever inside
the family you start.

### Predictions
| # | prediction |
|---|---|
| **Q1** | from a single DAZ copy the closure reaches **exactly the four DAZ copies**, not the amplicon |
| **Q2** | from a single gorilla NPIP locus it reaches the **same 25 truth loci** as the group seed |
| **Q3** | **well-definedness**: every single-locus start inside one family reaches the same set (tested on all four DAZ copies and on five gorilla NPIP loci) |
| **Q4** | each minimal set is a **subset** of the set the group seed reaches |

### Interpretation fixed in advance
- Q1+Q3 holding ⟹ "the smallest closed set containing the locus" is the definition to write, and it disposes
  of register 735.
- Q3 failing ⟹ the family depends on WHICH member you start from, which is worse than depending on a seed set;
  the repair fails and the definition must name its seed explicitly. That will be reported, not patched again.
- ⚠ A single-locus seed makes the majority condition vacuous on the first round (there are no "other members"),
  so round 1 admits whatever aligns; the test is what the SECOND round prunes it to. If a single locus grows to
  the amplicon and stays there, Q1 fails and the amplicon is simply the smallest closed set on chrY.

---

## AMENDMENT 2 — a member is a LOCUS, not an alignment fragment (2026-09-07, before the runs)

Register 739 diagnosed: at the 4-copy DAZ level two members are alignment fragments inside DAZ1 and DAZ4
(23.2 kb and 20.8 kb) that hold **9 % and 5 % of those loci's exons**. The extend step admitted the piece that
happened to align rather than the locus it lies in, so members of one level are not comparable and the
exon-based level rule preferred the 2-copy level.

**Change:** an admitted interval is grown to the locus containing it — the overlapping annotated span where one
exists (the annotation proposes extent, which the ablation showed is robust to ±5 kb), otherwise the interval
is kept as-is and marked **unannotated**. The duplicon stays what the loci SHARE; the member becomes the locus.

### Predictions
| # | prediction |
|---|---|
| **R1** | at the 4-copy DAZ level all four members are whole loci |
| **R2** | exon coverage of the shared segment at that level is **> 0.5 for every member** |
| **R3** | the level rule ("largest level whose shared segment still carries the transcribed exons") then picks the **4-copy** level over the 2-copy one |
| **R4** | gorilla NPIP is **unchanged**: 25 loci from a single-locus seed, precision 1.000 |

R3 failing ⟹ level selection is not an exon question and the chain is reported without a chosen level.

# Empirical rules for near-identical copies — don't merge what RNA can still separate

Derived from the Soto-family recovery on A119b: 83 human segdup families, 295 well-covered members
(≥20 reads), 9,436 pairwise member alignments, cross-referenced with which members our de-novo detection
resolved as distinct copies vs collapsed.

## Finding 1 — sequence %identity does NOT predict separability (the headline)

Resolved and collapsed copies have the **same** identity distribution:

| well-covered members | median nearest-sibling id | reaches |
|---|---|---|
| **resolved** as a distinct copy | 99.7 % | 100.0 % |
| **collapsed** (not resolved) | 99.7 % | 100.0 % |

Collapse rate is **flat (~30 %) across every identity band from 95 % to 100 %** — no threshold separates the
two populations. Concrete proof: **AMY1A / AMY1B / AMY1C are 99.89–99.95 % identical yet ALL resolved as
separate copies.**

> **Rule 1: never use a sequence-identity threshold (SD98-style) to decide whether two loci are one copy or
> two.** It would merge separable copies (AMY) and fail to flag the truly-unresolvable ones. Identity is the
> wrong axis.

## Finding 2 — resolution is set by the read-visible EXONIC signal, not genomic similarity

Genomic identity includes introns and flanks, which RNA never sees. Two copies can be ~100 % identical
genome-wide and still separate on a handful of **exonic PSVs** or a **copy-specific junction** (the AMY case,
and the DAZ2 junction-only case); or collapse when their **expressed exon sequence is identical** regardless
of how divergent their introns are.

> **Rule 2: compute the keep-separate decision on the EXON-SUM sequence + junctions only — never on the full
> genomic span.** The distinguishing evidence must be (i) a read-supported PSV inside a shared exon, or
> (ii) a copy-specific splice junction. One is enough.

## Finding 3 — the genuine-collapse floor is ~15–17 %, and it's the K=0 wall, not a tuning miss

Of well-covered members: **71 % resolved (own copy), ~11 % have a copy nearby (offset/adjacent —
recoverable by better locus stitching), ~15–17 % genuinely collapse** (no separate locus even within 50 kb).
That residual is the copies whose *expressed* sequence is identical — the irreducible K=0 floor.

> **Rule 3: the genuine-collapse residual is where DNA parCN is required — do NOT loosen thresholds to
> "recover" it.** Loosening breaks the 100 % precision we measured (every called copy a real member). Flag it
> as "K=0 unresolvable from RNA" and route to DNA, rather than fabricating a merge or a split.

## Finding 4 — a near-identical copy needs its OWN assignable reads to be seeded

Collapse persists even at ≥20 primary reads, so it isn't raw depth. Part of the ~11 % "nearby copy" tail is
copies under-seeded because minimap2 placed their reads as *primaries at a sibling* (a secondary-sink). The
multimapper ASSIGNMENT (PSV/junction) must feed the detection step, not just the aligner's arbitrary primary.

> **Rule 4: seed de-novo loci from assignable reads, not from aligner primaries alone** — otherwise a copy
> whose reads minimap2 parked on its twin never gets its own locus.

## The operational "don't-merge" test (all four, distilled)

Two co-expressed loci are **two copies** (keep separate) iff, on their **shared exons**, the reads support
**≥1 PSV** (a column where ≥2 alleles each clear the read-support floor) **OR** the loci carry **≥1
copy-specific junction** — *independent of % sequence identity*. If neither holds (exon-sum + junctions
identical), they are **K=0 unresolvable from RNA**: report the copy *number* (χ_H / genome projection) and
route the split to DNA parCN. Never merge on identity; never split without an exonic witness.

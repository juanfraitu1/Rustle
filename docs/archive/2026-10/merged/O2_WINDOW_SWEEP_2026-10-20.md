# O2's window past Eichler's 10 AS: a development sweep on gorilla KB3781 (2026-10-20)

**User:** "we could also go past the 10 AS difference from Eichler, if that is too close to the fully tied impossible
cases"; "Sounds good, please proceed".

**Status.** Development families of `PREREG_o2_personal_copies_2026-10-20` (140 families), descriptive and exploratory. No
verdict; the held-out half was not used. Code: `bench/o3_graph/o2window.py` (+ `test_o2window.py`). Work dir:
`/mnt/linuxdisk/tmp/o3_graph/personal/window/`.

## What was swept

- **The window.** A read is tied in window W when its alignments within W of its best `ms` lie on 2 or more copies of one
  family, none off them, minimap2's primary among them. `ms` is minimap2's intron-free score; on these reads it tracks the
  edit count, one substitution = 5 points. W = 10, 20, 30, 50: **89,777** development reads are tied in some window.
- **Stitched alignments.** An alignment laying >= 50 bp on two copies' exons (SM392's chimeras) is no placement for O2.
- **O2's pick:** the tied copy with the unique fewest edits over the common window (Amendment 1's rule).
- **O2's set:** the tied copies that no other tied copy beats at every substitution column where the two differ
  (`pareto_set`). It drops a copy only when the read contradicts it and never supports it against that rival.
- **The answer key.** Every read was aligned to his maternal and paternal assemblies (6 batches, about 5 min each).
  - The registered edit key (NM + clipped bases, strict minimum) **over-penalizes clipped read ends**. In the windows up to
    20, 344 of minimap2's 754 "errors" are decided only by minimap2's alignment clipping more bases (5-19 bp in 241 of
    them). Example: SM375.0 matches the read with 0 mismatches and clips 4 end bases; SM375.1 aligns all of it with 2
    mismatches. The edit key gives the read to SM375.1.
  - The sweep therefore also uses an **intron-free key**: the copy with the highest `ms` in his genome, at least one
    substitution (5 points) above any other copy or locus (`truth_ms`).

## Results under the intron-free key (each band = the reads first tied in that window)

| band (`ms` points from the best) | scored | minimap2 wrong | O2 pick fixes / breaks | set holds the truth | set of one copy (right) | mean pool → set |
|---|---|---|---|---|---|---|
| < 10 | 519 | 225 | 3 / 60 | 517 | 218 (216) | 2.01 → 1.59 |
| 10-19 | 1,257 | 191 | 1 / 22 | 1,255 | 980 (980) | 2.03 → 1.24 |
| 20-29 | 7,700 | 97 | 0 / 11 | 7,698 | 7,155 (7,153) | 2.01 → 1.08 |
| 30-49 | 80,151 | 18 | 1 / 55 | 80,150 | 79,520 (79,519) | 2.00 → 1.01 |

**By origin, bands < 30:**

| band | reads whose copy matches the reference: scored / minimap2 wrong / set of one (right) | reads carrying the other parent's copy: scored / minimap2 wrong / set holds truth / set of one (right) |
|---|---|---|
| < 10 | 199 / 1 / 188 (187) | 320 / 224 / 319 / 30 (29) |
| 10-19 | 1,033 / 3 / 976 (976) | 224 / 188 / 222 / 4 (4) |
| 20-29 | 7,583 / 0 / 7,146 (7,146) | 117 / 97 / 115 / 9 (7) |

minimap2's errors in these three bands are almost all **SM251**: 212, 188 and 94. Its reads come from his maternal SM251.0,
which is about 41 edits from the reference's own copy.

## Reading

- **Widening adds reads minimap2 already gets right.** Past 20-30 points minimap2 is almost never wrong (18 of 80,151 in
  the 30-49 band). Nearly every read it gets wrong is inside 30.
- **The set does what the user asked.** It shrinks the pool without dropping the true copy:
  - it holds the truth on 99.6-100% of reads in every band;
  - it narrows to one copy for about 95% of reads whose copy matches the reference, and that copy is right;
  - for reads from his divergent maternal SM251.0 it keeps both copies (sets of one: 43, 40 of them right). An allele-blind
    rule cannot name that copy, and the set does not guess.
- **O2's single pick (fewest edits) breaks more than it fixes in every band.** The set is the usable output.
- **Caveat: the reference-matching reads are close to self-scoring.** For them the intron-free key on his own haplotype is
  the reference's own score ordering. minimap2's near-perfect score there, and the set holding the truth, are largely by
  construction. The independent evidence is the other-parent reads: the set holds the truth on 656 of 661 of them, but
  narrows to one copy on only 43.

## Next

Pre-register "intron-free window 30 + the set rule" and test it on fresh reads:
- the mirror reference (his other haplotype as the reference), or HG002 with its own assembly;
- simulations with a known origin.

Score the set by:
- truth held;
- set size against the pool;
- singletons right.

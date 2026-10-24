# Graph O2 with a structure-free fallback, and cross-family-guarded O3, on a sixth handpicked family set (2026-10-16)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User goal:** "improve the simulation and use handpicked families (maybe soto or NPIP or other mainstays) to improve the graph
based O2 and O3". Round 4 (`docs/PREREG_o3_graph_round4_2026-10-15.md`, register 1302) claimed O3 for families without
sister-family hosting, and did not claim O2. Both failures had one root: the graph method gave up on reads it could not
explain within a family.
- **O2 failed A1 (0.969):** in CROSS it abstained on paralog-structure reads it does not thread, which the current O2 places
  correctly.
- **O3:** on development set 2 it pooled another family's reads (NPIPA.5 hosted on NPIPB) into a family's unexplained set
  and called them.

## Methods under test

**O2: graph O2 with a fallback** (`o2graph.py run WORLD unionf 2`, after `unionp 2` and `seqsupport`), `fallback_call`:
- padded union O2 decides first (round 4);
- a read it cannot explain in this family is `other_family` when another family's union explains it (within delta);
- otherwise the read is decided on its spliced alignments to each catalogued copy's genomic locus (`locus`: any canonical
  junction, the same margin and delta).

**O3: K** (`detect.py`): as round 4's H, but its pool is the reads the fallback O2 leaves `unexplained`. So another family's
reads and reads the family's own loci explain stay out.

Beside them: the current O2 and CUR; padded union O2 and H (round 4); locus O2; S23, S3, S123m; Kn (K before the genome
check).

## Worlds

As round 4, seed 61: INTACT, OUTSIDE, REMOVED (special copy hard-masked, own index), BRANCH (novel 150-bp exon), and CROSS
(`--projection liftover`).

## Families (held out, picked before any run)

The candidates were taken in a fixed order: the TRIM and LILR subfamilies first, so that sister families could occur, then
SPDYE, STAG3L, MT1, GSTM, CEACAM, FAM230, HERC2P, CYP2D, BTN3A, SPANXN, GUSBP, ANKRD20A, FRG2, GOLGA6 and PMS2. Too few
qualified at first, so the last five were added in order, still before any run. The set is the first 8 with >= 3 copies.

| family | copies / isoforms |
|---|---|
| TRIM64 (64, 64B, 64C) | 3 / 3 |
| TRIM51 (TRIM51 and 8 pseudogenes, chr11 and chr2) | 9 / 12 |
| SPDYE | 10 / 17 |
| STAG3L | 3 / 12 |
| FAM230 | 9 / 20 |
| HERC2P | 3 / 4 |
| GUSBP | 3 / 3 |
| ANKRD20A (8 chromosomes) | 10 / 12 |

`/mnt/linuxdisk/tmp/o3_graph/families5.json` (the earlier sets identical to `families4.json`).

## Known before running (development: the five earlier family sets)

O2, the fallback graph O2 against round 4's bars (re-labelled C1-C4 below):

| set | C1: lowest world ratio of correct calls to current O2's | C2: lowest precision | C3: REMOVED and BRANCH fewer; pooled non-catalogued ratio | C4: vs locus O2 | correct: fallback / padded union / current |
|---|---|---|---|---|---|
| dev | 1.025 | 1.0 | yes; 0.52 | x1.049 | 11,389 / 11,315 / 10,939 |
| dev-2 | 1.083 | 1.0 | yes; 0.41 | x1.126 | 12,702 / 12,701 / 11,530 |
| dev-3 | 1.065 | 1.0 | yes; 0.45 | x1.055 | 7,115 / 7,106 / 6,678 |
| dev-4 | 0.997 | 1.0 | yes; 0.54 | x1.056 | 9,220 / 9,185 / 9,230 |
| dev-5 (round 4's set) | 1.002 | 1.0 | yes; 0.29 | x1.059 | 8,928 / 8,814 / 8,906 |

- Every paralog-structure read the fallback assigns is right: 957, 612, 270, 222 and 471 of them; there are no wrong-family
  calls.
- The fallback places more non-catalogued reads than padded union (dev-4: 479 vs 375). C3's pooled ratio rises towards 0.54.

O3, K on dev-2, the only development set with sister-family hosting (550-790 reads per world hosted on another family):

| world | H false-call families | K false-call families | K detected |
|---|---|---|---|
| INTACT | 1 (NPIPB) | 0 | - |
| CROSS | 1 | 0 | - |
| OUTSIDE | 1 | 0 | - |
| REMOVED | 1 | 0 | 2 / 2 detectable |
| BRANCH | 0 | 0 | 8 / 8 |

K on the other development sets was still running at commit; it is reported under Results.

## Bars (held out: the 8 families above, seed 61; wrong-family calls count as wrong)

| id | bar |
|---|---|
| C1 | fallback graph O2 correct >= 0.98 x the current O2's, in every world |
| C2 | fallback graph O2 precision >= 0.995 over reads of catalogued copies, in every world |
| C3 | in REMOVED and in BRANCH it places fewer reads of non-catalogued copies than the current O2, and <= 0.6 x as many over OUTSIDE, REMOVED and BRANCH together |
| C4 | fallback graph O2 correct >= 1.03 x locus O2's, pooled over the five worlds |
| D1 | K detects >= 6 of the 8 BRANCH absent copies |
| D2 | K has a false call in <= 1 of the 24 INTACT, OUTSIDE and CROSS families |
| D3 | K detects every detectable REMOVED copy CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |
| D4 | K detects at least as many absent copies (REMOVED + BRANCH) as S23, with no more false-call families |
| D5 | in every family-world with >= 100 reads hosted from another family, K has no false call; not measurable without one |

Detectable (D3): the removed copy's reads have a median primary `de` on the masked genome > delta (`summary.py o3`).

## Readings (fixed now)

1. **C1-C4 all pass:** claimed, "on simulated reads of handpicked real families, graph O2 that threads every isoform structure
   through every copy, and falls back to each copy's genomic locus for reads no structure explains:
   - keeps >= 98% of the current O2's correct assignments in every world, including when isoforms are annotated on another
     paralog;
   - has >= 99.5% precision;
   - places <= 60% as many reads of non-catalogued copies;
   - keeps >= 3% more correct assignments than a locus-only O2."

   Any failing bar: not claimed.
2. **D1-D4 all pass:** claimed, round 4's O3 sentence for K. D3 not measurable: claimed without the removed-copy clause. D3
   measurable and failing, or any other bar failing: not claimed. D5 is reported. If it is measurable and fails, the claim
   keeps round 4's scope (families without sister-family hosting).
3. Each claim stands or falls as a whole. A failed bar is reported by name.

Scope either way: simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

## Results (held out, 2026-10-16)

Run: code and prereg `1c49489a`, the frozen copy_assign binary, seed 61. Outputs in `/mnt/linuxdisk/tmp/o3_graph/real/heldout5/`
(`r5_*.json`, `o3_score.json`). Summaries from `bench/o3_graph/summary.py` (`o2bars-unionf`, `o3`).

**O2** (correct / assigned; wrong-family calls counted):

| world | current | **fallback graph O2** | locus O2 | non-catalogued placed: current / fallback |
|---|---|---|---|---|
| INTACT | 2,357 | 2,546 / 2,546 (x1.080) | 2,456 | - |
| CROSS (liftover) | 2,793 | 3,005 / 3,006 (x1.076) | 2,922 | - |
| OUTSIDE | 2,066 | 2,260 / 2,260 (x1.094) | 2,237 | 199 / 449 |
| REMOVED | 2,066 | 2,260 / 2,260 (x1.094) | 2,237 | 233 / 449 |
| BRANCH | 2,357 | 2,546 / 2,546 (x1.080) | 2,456 | 104 / 0 |

- On this set the uncatalogued (OUTSIDE) and removed (REMOVED) special copies are all within delta of a catalogued paralog: the
  removed copies' median `de` is 0.0013-0.0077, so none is detectable.
- Their reads are explained by that paralog. The fallback, like any distance-based O2, places 449 of them on it (wrong by
  construction). The current O2 abstains on more of them: 199 and 233.

**O3** (8 families per world):

| | CUR | S23 | S3 | H | **K** | Kn (before the genome check) |
|---|---|---|---|---|---|---|
| BRANCH detected | 0 | 8 | 0 | 8 | **8** | 8 |
| INTACT + OUTSIDE + CROSS false-call families (of 24) | 2 | 0 | 0 | 0 | **0** | 5 |
| REMOVED + BRANCH false-call families | 2 | 0 | 0 | 0 | **0** | 4 |
| REMOVED detected (0 detectable) | - | - | - | - | - | - |

**Sister-family hosting (D5).** TRIM64 hosts 143-230 reads from another family (TRIM51) in every world, and TRIM51 hosts 60 in
REMOVED. K has no false call in any of them; Kn has TRIM64 or TRIM51 false calls in all five worlds, which the guard and the
genome check remove. H, without the guard, has none either, so on this set the genome check alone sufficed. On development
set 2 it did not: H had NPIPB false calls, and K none.

| bar | value | |
|---|---|---|
| C1 | lowest world ratio 1.076 | pass |
| C2 | lowest precision 0.9997 (CROSS 3,005 / 3,006) | pass |
| C3 | REMOVED 449 vs 233 (more); pooled 898 / 536 = 1.68 | **fail** |
| C4 | 12,617 vs 12,308: x1.025 | **fail** |
| D1 | 8 / 8 | pass |
| D2 | 0 of 24 | pass |
| D3 | 0 detectable removed copies | not measurable |
| D4 | K 8 = S23 8, false-call families 0 = 0 | pass |
| D5 | measurable (TRIM64, TRIM51); K no false call | pass |

### Verdict (2026-10-16)

**O3: claimed without the removed-copy clause (reading 2: D1, D2, D4 pass; D3 not measurable).** "Graph O3 that clusters what
no family's graph, nor the family's own loci, explains detects novel-exon absent copies (8 of 8; the current method 0). It
makes no false call over intact, uncatalogued-paralog and paralog-annotation families (0 of 24; the current method 2),
including families that host another family's reads (D5)."

Together with round 4, which claimed the removed-copy clause for H (6 of 6) on families without sister-family hosting:
- the graph-based O3 holds on two fresh sets;
- K detects everything H detects on all five development sets, and has no false call on any of them.

Attribution, as in round 4: the novel exons come from the branch calls, and the false-call drop from the genome check plus,
for sister families, the cross-family guard.

**O2: not claimed (C3, C4 fail).** The fallback graph O2 assigns 7.6-9.4% more reads correctly than the current O2 in every
world, at precision >= 0.9997, including the paralog-annotated isoforms (C1, C2). But:
- where an uncatalogued copy is within delta of a catalogued one, it places that copy's reads on the paralog, as any
  distance-based O2 would, while the current O2 abstains on more of them (C3);
- it keeps only 2.5% more correct calls than locus O2 here (C4).

Across six fresh family sets graph O2 has traded one failure for another:
- round 4: abstaining on structures it does not thread;
- round 5: assigning reads of near-identical uncatalogued copies.

No version has met all of its pre-registered bars.

## Review corrections (round-5 review, 2026-10-16; post hoc)

A fresh reviewer re-derived every number above with independent code, and all of them reproduce. Timing is clean: commit
00:52:34.052, first held-out file 7 ms later. The family selection reproduces exactly (first 8 with >= 3 copies, in GROUPS
order, from `families5_all.json`). No code bug changes a number. The corrections are about what the held-out set showed.

**D5 passed without exercising the guard (I1).**
- None of the TRIM51 reads hosted on TRIM64 (148-230 per world) aligns to TRIM64's loci, so no sister-family PSV cluster forms
  there. Development set 2 had 30-43 per world.
- Hn's and Kn's sister-family calls are identical. All of them are branch calls, and the genome check drops all of them
  (median `de` 0.0011-0.0041, no novel bases).
- So "which the guard and the genome check remove" was wrong: the guard removed nothing.
- "Measurable (TRIM64, TRIM51)" was wrong too: TRIM51 hosts 60 foreign reads, under the 100 needed; only TRIM64 qualifies.
- D5 is passed but says nothing about the guard. The sister-family evidence for K is development set 2 only.

**The guard covers PSV clustering only (I2).** K's branch calls are built from all of a family's hosted reads, as S2's are, not
from K's pool. So another family's reads with >= 50 genome-unaligned bases, or `de` > delta, could still be called; only the
genome check stands in the way. The prereg's "another family's reads ... stay out" holds for K's PSV clustering only. Not
changed after the fact: it would be a method change for a new prereg.

**Wording, attribution and history (I3, I4).**
- The claim sentence should be the pre-registered one (below).
- On this set K, H and S23 are equal in every cell: 8 / 8 BRANCH, 0 / 24 false-call families, 0 in REMOVED + BRANCH. S3,
  the current method plus the genome check, is also 0 / 24.
- All 8 detections are branch calls that S23 also makes, and the removed-copy clause had no data. So none of K's own parts
  (the fallback pool, the guard, the mosaic rule, pooling) was tested on held-out data.
- On development data K adds over S23 only in REMOVED (GAGE, GOLGA8, PSG, CGB, KIR, PRAMEF), exactly the clause this set
  could not measure.
- This is round 5, and the detector changed after each round. The guard fixes false calls that round 4's mosaic extension
  introduced; S23 has none on dev-2.
- "Holds on two fresh sets" mixed H (round 4) and K (round 5).

**Where C3's failure comes from (I5).** Placements of the 481 non-catalogued reads per world:

| world | current O2 | padded union | locus O2 | fallback |
|---|---|---|---|---|
| OUTSIDE | 199 | 303 | 437 | 449 |
| REMOVED | 233 | 303 | 437 | 449 |

- The fallback adds 146 per world: its own `locus` calls on reads padded union left `unexplained` (TRIM51.7 60, HERC2P 59,
  STAG3L 27). It cannot tell a catalogued copy's unthreaded structure from a near-identical uncatalogued copy.
- Its own calls on held-out: 48 correct against 292 non-catalogued reads placed. On the development sets: 233 against 259.
- Counting non-catalogued placements, precision is 0.834 in OUTSIDE and REMOVED, against the current O2's 0.912 and 0.899.
  The ">= 0.9997" was C2's catalogued-only figure.
- On this set padded union alone would also fail C3 (REMOVED 303 > 233; pooled 1.13) and C4 (x1.021). It would pass C1
  (lowest 1.062, CROSS).
- So C3 fails mostly because of this family set, and the fallback makes it worse.

**The O2 history, by round (I6).** Five fresh held-out sets, each failing on a different bar:
- **Round 1:** path O2 passed its own bars (O2a-O2d). The review then found the gain annotation-decided, and Amendment 1's
  CROSS failed (precision 0.42 on paralog-structure reads).
- **Round 2:** union O2, U4 0.981.
- **Round 3:** padded union, no gain (x0.995).
- **Round 4:** padded union, A1 0.969 in CROSS.
- **Round 5:** the fallback, C3 and C4.

"Traded one failure for another" shows only on development data: the fallback raises CROSS (dev-5: 2,330 vs 2,252) and
raises non-catalogued placements (dev-4: 479 vs 375).

**Smaller.**
- The "dev" row pools six worlds (CROSS and CROSSL).
- Not fixed, recorded:
  - `summary.py` labels the fallback's output "padded union" and its bars A1-A4;
  - the guard does not stop a `locus` call when the read's own family's union also fails (none occurred);
  - there is no table-level test of `run_fallback` (the reviewer re-derived every row);
  - the selection step is not a committed command (it reproduces exactly from the rule above);
  - a future prereg should require >= 3 detectable removed copies before running, since detectability is computable
    without reads.

### Verdict, revised after the review (2026-10-16)

**O3: claimed without the removed-copy clause (reading 2), in the pre-registered words.** "On simulated reads of handpicked real
families, graph O3 that clusters a family's O2-unexplained reads together, with the mosaic rule, detects novel-exon absent
copies with at most one false call over intact, uncatalogued-paralog and paralog-annotation families." Observed: 8 of 8 and
0 of 24.

What this set tested: the branch calls (novel exons) and the genome check (no false call), exactly as S23 and S3 do. It did not
test K's own parts. The sister-family evidence for the guard is development set 2 only (K 0 false calls, H 4).

**O2: not claimed (C3, C4 fail).**
- The fallback assigns 7.6-9.4% more reads correctly per world, at precision >= 0.9997 over reads of catalogued copies.
- But its own `locus` calls are mostly reads of near-identical uncatalogued copies. Counting those, its precision is 0.834,
  against the current O2's 0.90-0.91 in OUTSIDE and REMOVED.
- Across five fresh held-out sets no graph O2 version has met all of its pre-registered bars.

## Post hoc: can any operating point of graph O2 dominate the current O2? (2026-10-16; all sets seen, no claim)

Worlds with an uncatalogued copy (OUTSIDE, REMOVED); a placement of a non-catalogued copy's read counts as wrong. Correct /
wrong:

| set | world | current O2 | fallback graph O2 |
|---|---|---|---|
| round 2's | REMOVED / OUTSIDE | 1,034 / 162; 1,034 / 119 | 1,102 / 119 (dominates both) |
| round 3's | REMOVED / OUTSIDE | 1,624 / 373; 1,621 / 241 | 1,619 / 235 (5 fewer correct, far fewer wrong) |
| round 4's | REMOVED / OUTSIDE | 1,437 / 121; 1,434 / 60 | 1,440 / 60 (dominates both) |
| round 5's | REMOVED / OUTSIDE | 2,066 / 233; 2,066 / 199 | 2,260 / 449 (trade-off) |

Could a per-read error test separate a near-identical uncatalogued copy's reads from sequencing error? Reads assigned to their
own copy sit at a median 1.6-1.7 differences per kb on the copy's locus, with a heavy tail (95th percentile 6-14). The
uncatalogued copies' reads sit at 4-5 per kb. No cutoff from 0.5 to 9.5 per kb makes graph O2 dominate where it does not
already.

So graph O2 is at least as good as the current O2 on most sets and worlds, and better in precision. It loses precision exactly
where an uncatalogued copy is within sequencing-error range of a catalogued one (round 5's set), and nothing tried here tells
those reads apart.

# Graph O2 as assign-or-abstain, and family-pooled mosaic-aware O3, on a fifth handpicked family set (2026-10-15)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User goal (2026-10-15):** "improve the simulation and use handpicked families (maybe soto or NPIP or other mainstays) to
improve the graph based O2 and O3".

## Where rounds 1-3 left it (registers 1298-1301, after their reviews)

No O2 or O3 claim has survived as pre-registered.
- **O2.** Graph O2's gains in correct assignments were decided by the annotation: per-copy isoform lists (round 1), one
  family's non-canonical intron (FCGR, round 2), and annotated end positions (round 2's development families). With
  padded paths, nothing was gained on round 3's fresh set (x0.995).
- **O2, what held every time:** graph O2 keeps the current O2's correct calls within 2% at precision >= 0.999. It also
  abstains on reads no catalogued copy explains within delta, placing far fewer reads of non-catalogued copies (round 2:
  223 vs 608; round 3: 375 vs 891). An annotation-free locus O2 abstains the same way, but it loses 5-6% of the correct
  calls.
- **O3.** Branch calls and the genome check were right every time: 8 / 8 novel-exon absent copies on three fresh sets, no
  false call after the genome check. Removed copies failed the bar twice.
  - Round 2: PSG's reads were split over four hosts.
  - Round 3: KIR was missed because the mosaic rule was coded more narrowly than documented. It skipped linked clusters
    with few copy-private alleles; fixed in `b0a3e59f`.

This round tests the two things that held or were fixed, on a fresh family set, as one claim each.

## Methods under test

**O2: padded union O2** (`o2graph.py run WORLD unionp 2`), unchanged from round 3. Each copy's paths are its own curated
isoforms plus every paralog isoform spliced-projected onto its locus, with both ends padded by 50 bp of the copy's genome.
A read is assigned to the closest copy at margin 2 and abstains (`unexplained`) beyond delta. Beside it: the current O2
(copy_assign `--no-as-tied-only`) and locus O2 (the same decision on each copy's genomic locus).

**O3: H** (`detect.py`), as round 3, with two fixes from the round-3 review (`b0a3e59f`):
- the mosaic rule, as documented: every cluster no copy explains is called when its nearest single copy is > delta away;
- one call per pooled member set.

H clusters the reads padded union O2 cannot explain, pooled over all of the family's hosts. The mosaic rule and the
branch calls apply, and the genome check keeps a call only when the genome does not explain its reads.

## Worlds

As round 3, seed 51: INTACT, OUTSIDE, REMOVED (special copy hard-masked, own index), BRANCH (novel 150-bp exon), and CROSS
(`--projection liftover`). Worlds share INTACT's reads of catalogued copies (round-2 review M4), so the O2 bars are set per
world.

## Families (held out, picked before any run)

Candidates in a fixed order: TRIM, LILR, DEFB10, CSAG, MAGEB, CTAG, SERF1, CBWD, MAGEC, CYP2A, NBPF, GTF2I. Too few qualified,
so more candidates were added in a fixed order before any run: PRAMEF, VCX, SPDYE, FAM86, GPR89, STAG3L, POM121, CCDC144. The
set is the first 8 with >= 3 copies, so that a removed copy leaves >= 2 catalogued hosts (round 3's CYP4F had none).

| family | copies / isoforms |
|---|---|
| TRIM (TRIM49, 49B, 49C) | 3 / 5 |
| DEFB10 (DEFB103A/B) | 3 / 3 |
| CSAG | 4 / 13 |
| CYP2A (2A6, 2A7, 2A13) | 3 / 9 |
| NBPF (7, 8, 9, 11, 12) | 5 / 15 |
| GTF2I (pseudogenes GTF2IP2, P3, P8) | 3 / 3 |
| PRAMEF | 9 / 11 |
| VCX | 4 / 9 |

`/mnt/linuxdisk/tmp/o3_graph/families4.json` (earlier sets identical to `families3.json`; `families4_all.json` holds every
candidate).

## Known before running (development: the four earlier family sets)

Development sets: dev (the original development families), dev-2 (round 1's held-out set), dev-3 (round 2's), dev-4
(round 3's), all seen now.

O2, the bars on development (padded union O2 against the current O2 and locus O2):

| set | A1: lowest world ratio of correct calls | A2: lowest world precision | A3: non-catalogued placed, padded union vs current (OUTSIDE, REMOVED, BRANCH; pooled ratio) | A4: vs locus O2 |
|---|---|---|---|---|
| dev | 1.017 | 1.0 | 118 vs 118, 118 vs 228, 8 vs 314; 0.37 | x1.041 |
| dev-2 | 1.083 | 1.0 | 87 vs 67, 87 vs 119, 0 vs 240; 0.41 | x1.125 |
| dev-3 | 1.060 | 1.0 | 93 vs 119, 93 vs 162, 37 vs 327; 0.37 | x1.053 |
| dev-4 | 0.983 (CROSS) | 1.0 | 183 vs 241, 183 vs 373, 9 vs 277; 0.42 | x1.052 |

In OUTSIDE the uncatalogued paralog is in the genome and mostly within delta of a catalogued copy, so no abstention can
tell its reads apart (by O3's own definition they are explained). Padded union O2 places as many as the current O2 on dev
(118 vs 118) and more on dev-2 (87 vs 67). A3's "fewer" clause therefore covers REMOVED and BRANCH; OUTSIDE is reported.

O3, H with the post-review fixes (`b0a3e59f`), re-run on dev-3 and dev-4:

| set | BRANCH detected | REMOVED detected / detectable | false-call families, INTACT + OUTSIDE + CROSS | Hn false-call families (before the genome check) |
|---|---|---|---|---|
| dev-3 | 8 / 8 | 3 / 3 (FCGR, LRRC37A, PSG) | 0 / 24 | 8 world-families |
| dev-4 | 8 / 8 | 3 / 4 (CGB, KIR, UGT2B; CYP4F unreachable) | 0 / 24 | 6 world-families |

The same re-run on dev and dev-2 was still running when this was committed. Its numbers are added under Results; they do
not change any bar.

## Bars (held out: the 8 families above, seed 51; wrong-family calls count as wrong)

| id | bar |
|---|---|
| A1 | padded union O2 correct >= 0.98 x the current O2's, in every world |
| A2 | padded union O2 precision >= 0.995 over reads of catalogued copies, in every world |
| A3 | in REMOVED and in BRANCH padded union O2 places fewer reads of non-catalogued copies than the current O2, and <= 0.6 x as many over OUTSIDE, REMOVED and BRANCH together |
| A4 | padded union O2 correct >= 1.03 x locus O2's, pooled over the five worlds |
| B1 | H detects >= 6 of the 8 BRANCH absent copies |
| B2 | H has a false call in <= 1 of the 24 INTACT, OUTSIDE and CROSS families |
| B3 | H detects every detectable REMOVED copy CUR detects, and >= 0.75 of the detectable ones; needs >= 3 detectable |
| B4 | H detects at least as many absent copies (REMOVED + BRANCH) as S23, with no more false-call families |

Detectable (B3): the removed copy's reads have a median primary `de` on the masked genome > delta.

## Readings (fixed now)

1. **A1-A4 all pass:** claimed, "on simulated reads of handpicked real families, O2 on a family graph that threads every isoform
   structure through every copy keeps >= 98% of the current O2's correct assignments at >= 99.5% precision. It abstains on
   reads no catalogued copy explains, so it places fewer reads of non-catalogued copies (<= 60% as many). It keeps >= 3% more
   correct assignments than an annotation-free locus O2 with the same abstention." Any failing bar: not claimed.
2. **B1-B4 all pass:** claimed, "graph O3 that clusters a family's O2-unexplained reads together, with the mosaic rule, detects
   novel-exon absent copies and >= 75% of detectable removed copies, with at most one false call over intact,
   uncatalogued-paralog and paralog-annotation families". B3 not measurable (< 3 detectable): claimed without the
   removed-copy clause. B3 measurable and failing, or any other bar failing: not claimed.
3. Each claim stands or falls as a whole. A failed bar is reported by name.

Reported beside the bars:
- path O2;
- G, S23, S123m and Hn (H without the genome check);
- O2 per family;
- the liftover CROSS structures whose intron chain union does not thread (round-2 review I2).

Scope either way: simulated reads (curated CAT isoforms, HiFi errors) on the real CHM13 genome. Real data's other species,
fragments and allele background are not modelled.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

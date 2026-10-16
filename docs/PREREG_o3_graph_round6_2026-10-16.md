# O2 + an unannotated-copy flag, with heterozygous alleles, on a seventh handpicked family set (2026-10-16)

Pre-registration, written before any world of the new held-out set was built. Code: branch `o3-graph` (`bench/o3_graph/`),
the commit that adds this file.

**User question (2026-10-16):** "can we check if O2 is doing that because it is detecting unannotated copies? In that case that
could be useful". The answer was to pre-register the flag on a fresh family set, with alleles so that the copy-versus-allele
cost is measured; the user agreed ("yes please").

## Background (post hoc, round 5's set and earlier, no alleles)

O2's wrong placements in OUTSIDE and REMOVED worlds are reads of the uncatalogued copy, placed on its nearest catalogued
paralog. As a group they often form a pure PSV cluster whose copy-private alleles are linked. That is an unannotated copy's
fingerprint, usually at 0.1-0.5% divergence, below delta (0.96%).

Across six family sets, the flag defined below:
- covered 24 of 48 removed copies, 13 of them below delta;
- held 45% of O2's misplaced reads;
- falsely flagged 4 of 48 intact families.

There are no alleles in those worlds. In real data a heterozygous site gives the same signature (linked differences shared by
about half of a copy's reads). This round measures that cost.

## Simulation: heterozygous alleles (`realworld.py --het RATE`)

The donor is diploid. Each copy (catalogued or not) has heterozygous sites at RATE of its exonic bases, and its haplotype 2
carries the alternative bases. Each read comes from haplotype 1 or 2 with probability 1/2 (truth column `hap`). **RATE = 0.001**,
human heterozygosity. With `--het 0` a world is byte-identical to the earlier simulator's.

## Methods under test

- **The flag** (`o2graph.flag_reads`): a PSV cluster of the current method's per-host clustering (`copy_assign
  --psv-clusters`) whose copy-private alleles are linked (`unexplained`), with >= 20 reads, whatever its divergence.
- **O2 + flag** (`o2graph.py run WORLD flag unionp`): padded union O2 (round 3), with every read of a flagged cluster reported
  as an unannotated-copy candidate instead of being assigned.

Padded union is the base because, without alleles, padded union + flag had at least the current O2's correct calls (within
0.4%) and fewer wrong ones in all 18 development worlds.

## Worlds

INTACT, REMOVED (special copy hard-masked, own index) and OUTSIDE, all with `--het 0.001`, seed 71.

## Families (held out, picked before any run)

The candidates were taken in a fixed order: FRG2, GOLGA6, PMS2, CCL3/4, DEFA, BCRP, GGT, CFHR, PRB, ANKRD30, FRG1, GSTT2, CTAGE,
RIMBP3. FRG2, GOLGA6 and PMS2 were round-5 candidates past its first 8 and never built into a world. The set is the first 8
with >= 3 copies:

| family | copies / isoforms |
|---|---|
| FRG2 | 5 / 9 |
| GOLGA6 (6B, 6C, 6D) | 3 / 5 |
| PMS2 (PMS2 + 9 pseudogenes) | 10 / 20 |
| CCL34 (CCL3, CCL3L1, CCL3L3) | 3 / 9 |
| DEFA (DEFA1, DEFA1B, DEFA3) | 3 / 4 |
| BCRP | 5 / 5 |
| GGT | 5 / 10 |
| FRG1 | 10 / 10 |

`/mnt/linuxdisk/tmp/o3_graph/families6.json` (the earlier sets identical to `families5.json`).

## Known before running (development: round 5's families with `--het 0.001`, seed 61; `real/dev_het/`)

| | INTACT | REMOVED | OUTSIDE |
|---|---|---|---|
| flagged clusters: copy / allele / other | 0 / 13 / 8 | 6 / 10 / 7 | 6 / 10 / 7 |
| uncatalogued copies flagged | - | 2 of 8 | 2 of 8 |
| families with a false flag | 7 of 8 | 7 of 8 | 7 of 8 |
| correct calls: current O2 / padded union / padded union + flag | 2,356 / 2,542 / 1,886 | 2,065 / 2,256 / 1,727 | 2,065 / 2,256 / 1,727 |
| wrong calls (incl. non-catalogued): current / padded union / + flag | 9 / 0 / 0 | 226 / 297 / 165 | 192 / 297 / 165 |

With alleles the flag fires on one haplotype of most catalogued copies. It costs 26% of padded union's correct calls in
INTACT, and the uncatalogued copies' reads split by haplotype, so fewer copies are covered. The same families without alleles
had 2 false-flag families and 4 of 8 copies flagged.

## Bars (held out: the 8 families above, seed 71, `--het 0.001`)

| id | bar |
|---|---|
| H1 | the flag covers >= 50% of the uncatalogued copies (REMOVED and OUTSIDE together, 16): >= 80% of a copy's hosted reads in flagged clusters that are >= 90% its reads |
| H2 | the flag falsely flags (a cluster not >= 90% one uncatalogued copy's reads) <= 1 of the 8 INTACT families |
| H3 | O2 + flag keeps >= 0.98 x padded union O2's correct calls in INTACT |
| H4 | O2 + flag places fewer wrong reads (non-catalogued + wrong copy + wrong family) than the current O2 in REMOVED and in OUTSIDE |

Truth classes come from `summary.py flags`: `copy`, `allele` (>= 90% one catalogued copy's reads, >= 90% of those from one
haplotype) and `other`.

## Readings (fixed now)

1. **H1-H4 all pass:** claimed, "with realistic heterozygosity, flagging PSV clusters of linked copy-private alleles below
   delta turns O2's misplaced reads of unannotated copies into candidates. It covers >= half of them, falsely flags <= 1 of 8
   intact families, and keeps >= 98% of O2's correct calls."
2. **H2 or H3 fails:** not claimed. "Without allele information (DNA of the same individual, or phasing), the flag cannot
   tell an unannotated copy below delta from a heterozygous allele." The development data predicts this reading.
3. **H1 fails alone:** not claimed; the flag misses unannotated copies even where alleles are not the problem.

Scope: simulated reads (curated CAT isoforms, HiFi errors, heterozygous sites at 0.1%) on the real CHM13 genome. No allelic
expression imbalance, no somatic or editing variation.

## Results

Added after the runs. Nothing above this heading changes once a statistic has been seen.

## Results (held out, 2026-10-16)

Run: code and prereg `ef78d6ed`, the frozen copy_assign binary, seed 71, `--het 0.001` (60 heterozygous sites). Outputs in
`/mnt/linuxdisk/tmp/o3_graph/real/heldout6/` (`flags.json`, `o2_*.json`); summaries from `summary.py flags` and `o2graph.py
score`.

| | INTACT | REMOVED | OUTSIDE |
|---|---|---|---|
| flagged clusters: copy / allele / other | 0 / 11 / 2 | 5 / 7 / 2 | 5 / 7 / 2 |
| uncatalogued copies flagged | - | 2 of 8 | 2 of 8 |
| families with a false flag | 4 of 8 (FRG2, GGT, GOLGA6, PMS2) | 4 | 4 |
| correct / wrong: current O2 | 1,983 / 0 | 1,683 / 157 | 1,683 / 157 |
| correct / wrong: padded union O2 | 2,104 / 0 | 1,790 / 118 | 1,790 / 118 |
| correct / wrong: padded union + flag | 1,701 / 0 | 1,505 / 117 | 1,505 / 117 |

"Wrong" counts reads of non-catalogued copies placed, wrong-copy and wrong-family calls.

| bar | value | |
|---|---|---|
| H1 | 4 of 16 uncatalogued copies covered (0.25) | **fail** |
| H2 | 4 of 8 INTACT families falsely flagged | **fail** |
| H3 | 1,701 / 2,104 = 0.81 of padded union's correct calls kept | **fail** |
| H4 | 117 < 157 in REMOVED and in OUTSIDE | pass |

### Verdict (2026-10-16)

**Not claimed (reading 2: H2 and H3 fail).** Without allele information (DNA of the same individual, or phasing), the flag
cannot tell an unannotated copy below delta from a heterozygous allele.
- At human heterozygosity (0.1%), 11 of the 13 clusters it flags in intact families are one haplotype of a catalogued copy.
- It costs 19% of O2's correct calls.
- The unannotated copies' own reads split by haplotype, so it covers only 4 of 16.

So the answer to the question that started this round: O2's misplaced reads do carry an unannotated-copy signature, but at
the divergence where O2 misplaces them, alleles carry the same one.

**Reported beside the bars (no claim).** On this fresh set, padded union O2 without the flag beats the current O2 on both
counts in every world:
- INTACT 2,104 vs 1,983 correct (no wrong);
- REMOVED and OUTSIDE 1,790 vs 1,683 correct, and 118 vs 157 wrong.

Elsewhere, from the round-6 development tables (corrected after a first draft overstated it):
- it dominates in every world on the original development set and on round 2's set;
- on round 1's set it dominates except in OUTSIDE (87 vs 81 wrong);
- on round 3's and round 4's sets it is within 0.4% of the current O2's correct calls with fewer or equal wrong ones;
- on round 5's set, with or without alleles, it places more wrong reads (its uncatalogued copies are near-identical to
  catalogued ones). That is a candidate O2 claim for its own
pre-registration on a new family set: "padded union O2 dominates the current O2". It was not this round's test, so it is not
claimed here.

## Review corrections (round-6 review, 2026-10-16; post hoc)

A fresh reviewer re-derived every number above with independent code; all reproduce. Timing is clean (commit 09:32:19.491,
`heldout6/` 18 ms later; nothing above "## Results" changed). No Critical finding. The corrections are about why the bars
failed and what the flag did.

**Why H1 fails (I1).** It was not mainly haplotype splitting:
- Five of the six copies not covered have no linked copy-private alleles at all, so with MIN_ALT = 2 they can never be
  flagged, with or without alleles:
  - BCRP6: 6 alts, none unexplained;
  - CCL3L1: one unexplained alt;
  - DEFA1B: one private site, and DEFA got no het site;
  - FRG1.2 and PMS2.8: their haplotype-1 clusters are explained.
- Only GOLGA6B is a near miss: 47 / 60 = 78%, split over two hosts.
- GGT2 did split by haplotype (32 + 27) and was still covered.

So H1 also fails for reading 3's reason (copies with no PSV signature).

**The flag contributed almost nothing to O2 here (I2).**
- Padded union misplaces 118 reads: FRG1.2 60, CCL3L1 57, GOLGA6B 1. Only 1 of them is in a flagged cluster, against 45% on
  the six earlier sets.
- Of the 286 reads labelled `flagged` in REMOVED and OUTSIDE, 285 were correct calls of catalogued copies.
- So H4 passes because of padded union, not the flag.

**Dominance scope (I3).**
- The side observation holds for INTACT, REMOVED and OUTSIDE only. In CROSS, padded union keeps 0.983 (round 3's set) and
  0.969 (round 4's: that round's A1 failure) of the current O2's correct calls; round 4's REMOVED is 0.9958.
- This set has no CROSS world, and its REMOVED and OUTSIDE are identical for O2 (2,644 reads, same hosting).
- A dominance prereg needs a CROSS world.

**Realism (I4).** The verdict's direction is robust, but "19%", "4 of 8" and "11 of 13" depend on these design choices:
- Haplotype 1 is the reference exactly, and every alt allele sits on haplotype 2. A real donor's homozygous differences from
  the reference would make whole copies flaggable, so the true cost is probably higher. Random phase would lower the flag
  count; the reviewer's estimate is H2 still about 4 of 8 and H3 about 0.87, both failing.
- Het sites are a deterministic count (round(0.001 x exonic length)). Copies under 1.5 kb of exon get at most one site and can
  never be flagged, and the CCL3 copies sit just under. A Poisson draw would give about 6 of 8 families with a flaggable copy.
- Expression is uniform: about 30 reads per haplotype, above the 20-read minimum.

**Smaller.**
- **M1.** With het 0 every product is byte-identical except the truth file, which gains a `hap` column (all 1).
- **M2.** `classify_cluster`'s allele rule counts reads that cover no het site, and labels a >= 90% haplotype-1 cluster as
  an allele. By the truth, all 13 INTACT flags here are allele clusters (two were labelled `other`).
- **M3.** "13 below delta" in the background counts by PSV est_div; by O3's genome-`de` measure it is 10.
- **M4.** Pre-run statements:
  - "fewer wrong in all 18 development worlds": 5 tie;
  - "within 0.4%": the worst is 0.42%;
  - "fires on most catalogued copies": 13 of 44 held-out copies, which is most of those with >= 2 het sites.
- **M5.** `apply_flag` labels only reads O2 had assigned. Members of flagged `copy` clusters that O2 abstained on are not
  listed as candidates.
- **M6.** REMOVED and OUTSIDE are identical here, so "4 of 16" is 2 of 8 distinct copies counted twice.

### Verdict, revised after the review (2026-10-16)

**Not claimed (reading 2: H2 and H3 fail; H1 also fails for reading 3's reason).**

The answer to the question that started this round:
- On the six earlier sets, about 45% of O2's misplaced reads carried an unannotated-copy signature.
- On this fresh set, O2's misplacements come from copies with no PSV signature at all, while the flag fires on heterozygous
  alleles: all 13 intact flags are allele clusters.
- So the flag is not usable without allele information (DNA of the same individual, or phasing), and it is not what explains
  O2's misplacements.

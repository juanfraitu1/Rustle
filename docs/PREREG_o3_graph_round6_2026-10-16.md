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

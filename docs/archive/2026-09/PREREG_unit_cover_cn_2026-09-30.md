# Pre-registration: finer units for the family cover — duplicon x copy-number class, sliver-free usage (2026-09-30, KEY=unitcovercn)

Written before any statistic below was computed. Follows KEY=unitcover (HOLDS, held-out mean Jaccard 0.473 vs location 0.277, but exact
sets 10.7%, 30.3% of single-family genes falsely multi, NPIP shared genes 1-3 of 6 families). Same inputs, same split, same multi /
clean gene sets, same statistics. Script: `bench/soto_m2/soto_m2_unit_cover.py` extended with the arms below (arm A must reproduce
KEY=unitcover's numbers exactly, else stop).

## 1. The two refinements (fixed now)

- **Sliver-free usage (E):** a gene uses duplicon D only if D is the dominant duplicon (most bases; ties: ID order) of at least one of
  its exons (exons merged). No threshold: a sliver at an exon edge never dominates that exon. (KEY=unitcover used any overlap >= 1 bp.)
- **Copy-number classes (K):** the units are (duplicon, class). For a duplicon D, the clean genes using D that have a famCN are sorted by
  famCN and cut into classes wherever two consecutive values differ by 2 or more (single linkage at Soto's own rule: a pair joins when
  |famCN difference| < 2, i.e. MAD < 1). Owner of a class = the family with the most exonic bases on D among the class's clean genes
  (ties: lower family number). Clean genes without famCN do not vote; if D has no clean gene with famCN, D is one class owned as in
  KEY=unitcover. **A gene with famCN c** uses the classes of D holding a member within |c - x| < 2; none, no unit from D. **A gene
  without famCN** uses every class of D (Soto's pair rule never cuts a pair with only one famCN). Clean genes are predicted
  leave-one-out (removed from D's classes before their own prediction).

## 2. Arms

- **A:** KEY=unitcover (1 bp usage, duplicon units). **B:** E only. **C:** K only (1 bp usage). **D:** E + K (the proposed structure).
- Each with **S1C famCN** (Soto's copy numbers) and with **our famCN** (recomputed from the same 268 SGDP tracks; less circular).
- Per arm: mean Jaccard on multi genes per half, location baseline (unchanged), permutation p (owners shuffled among owned units,
  1,000 permutations, seed 20260930); exact set, recall, precision; clean genes leave-one-out: exact own family, two or more, empty.

## 3. Decision rule (primary: arm D with S1C famCN, held-out half)

- **REFINES:** permutation p < 0.01, held-out mean Jaccard above arm A's (0.473), and the clean-gene false-multi share below arm A's
  (0.303).
- **PARTIAL:** p < 0.01 and exactly one of the two improvements.
- **NO GAIN:** otherwise.
Arm D with our famCN, and arms B and C, are reported beside it without a verdict.

## 4. Circularity (stated now)

With S1C famCN the classes are built from the copy numbers Soto used to cut their families, so arm D with S1C famCN asks whether the
structure can **express** Soto's table, not whether it discovers it independently. Our famCN is an independent measurement of the same raw
data (r = 0.977 with S1C; 211 of 1,793 genes differ by 2 or more) and is the less circular check.

## 5. Seen before (disclosed)

KEY=unitcover's full output (section 7 of its pre-registration and `docs/SOTO_UNIT_COVER_2026-09-30.{md,tsv}`), including the NPIP-side
lines. Nothing about arms B-D has been computed.

## 6. Result

(Filled in after the run, below this line, without editing anything above.)

Run on 2026-09-30 after this file was committed (`af100ada`, sha1 `f1cfdab2f50b4d841ee5d36b47fd2da6c73d5d45`). Deviation (form only): the
arms live in a new script, `bench/soto_m2/soto_m2_unit_cover_cn.py`, so that `soto_m2_unit_cover.py` stays as KEY=unitcover ran it; arm A
reproduces KEY=unitcover exactly (asserted: 0.439 / 0.473, same exact / recall / precision / clean shares). 9 s, light lock. Output
`docs/archive/2026-09/SOTO_UNIT_COVER_CN_2026-09-30.md`, per gene `docs/SOTO_UNIT_COVER_CN_2026-09-30.tsv`.

**VERDICT (section 3): REFINES.** Arm D with S1C famCN, held-out: mean Jaccard 0.660 (arm A 0.473), clean-gene false-multi 0.157 (arm A
0.303), permutation p 0.001.

| arm | famCN | Jaccard dev / held-out | exact | recall | precision | clean: own family / 2+ / empty |
|---|---|---|---|---|---|---|
| A (KEY=unitcover) | - | 0.439 / 0.473 | 0.107 | 0.555 | 0.718 | 0.546 / 0.303 / 0.020 |
| B sliver-free | - | 0.440 / 0.464 | 0.121 | 0.496 | 0.800 | 0.607 / 0.201 / 0.023 |
| C CN classes | S1C | 0.548 / 0.624 | 0.215 | 0.917 | 0.617 | 0.644 / 0.224 / 0.045 |
| C CN classes | ours | 0.492 / 0.601 | 0.161 | 0.817 | 0.608 | 0.626 / 0.225 / 0.053 |
| **D both** | **S1C** | **0.693 / 0.660** | **0.349** | **0.869** | **0.734** | **0.697 / 0.157 / 0.052** |
| D both | ours | 0.540 / 0.654 | 0.195 | 0.768 | 0.739 | 0.668 / 0.165 / 0.061 |

All permutation p = 0.001 (the floor); location baseline 0.221 / 0.277.

The two refinements do different jobs: sliver-free usage mostly removes false memberships (precision 0.718 -> 0.800, clean false-multi
0.303 -> 0.201), copy-number classes mostly add the missing ones (recall 0.555 -> 0.917); together they beat both. NPIP: with S1C famCN all
five shared genes get exactly Soto's six (or five) families; with our famCN 3-6 of them.

**Reading.** The structure — a partition of (duplicon, copy-number class) units, genes as paths, a gene without a copy number using every
class — expresses most of Soto's cover, including the NPIP shared genes exactly. With S1C famCN this is expressiveness, not independent
discovery (section 4); with our independently recomputed famCN the held-out Jaccard is nearly the same (0.654), while exact sets drop
(0.195), following the copy-number disagreements. Still incomplete: 65% of multi genes are not exact, 16% of single-family genes would be
called multi, 5% get no family.

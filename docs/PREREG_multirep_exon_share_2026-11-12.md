# PREREG row 1654: multiple representatives per locus under the exon-sharing edge rule

## Why
User question 2026-11-12: is there a better way to call families than one representative per locus? Rows 1575/1576
(`locus_reps.py --disjoint` + attach) added true pairs and lost none but ran under the old coverage edge rule
(`--min-shared-exon-frac 0.60`) and missed the precision bar by 0.005 on A119b. Exon sharing (row 1593, default since
2026-11-09) removes mostly false pairs. The two were never combined. Both samples are seen (A119b dev, HG002 used as the
held-out of 1576 and 1593): a pass here is development evidence; transfer = gorilla OR6737, a later row.

## Arms (binary `target_sddefault/release/mcl_families`, the one of the row 1593 port check; no new code)
Flags of the shipped default: `--min-exonic-bp 1 --emit-units --min-shared-exon-frac 0.0 --min-exon-share 0.5 --inflation 2.0`.
- Baselines (existing, `port1593/`): A119b damped `Ab`, plain `Aa`; HG002 damped `Hb`, plain `Ha`.
- F: the same flags on the existing multi-rep inputs `r1571/F.families.gtf` (A119b) and `r1576/F.families.gtf` (HG002),
  then `piece_merge.py`, then `piece_attach.py` onto the matching baseline's unchanged families:
  F2 = `--keep-all` (registered arm), Fa = without (descriptive).
Tracks: damped = the shipped default where a SEDEF BED exists (registered); plain = descriptive.
Scoring: pairscore (simiiformes default, PAIRSCORE_SAME_SYMBOL=ignore; HG002 `--species hg002`), universe = the
sample's unchanged `families.gtf`. pairdiag on baseline and F2 (damped) re-counts the miss classes under the new default.

## Bars: F2 vs baseline, damped track, each sample
1. Precision >= baseline (no tolerance).  2. Recall >= baseline + 0.02.  3. Baseline hits lost <= 5.
4. Copies <= 1.05 x baseline.
All four on both samples -> proposed as an opt-in families post-step (default = the user's call). Anything else = FAIL
as registered, no retuning on these samples.

## Prediction
HG002 passes all four. A119b: bars 2-4 pass, bar 1 borderline FAIL (baseline precision .918; attached pairs entered at
about .6 under the old rule).

## Result (2026-11-12, job 4467136, 1,217 s, `/anvme/workspace/v138da13-rustle/r1654_multirep_share/`)
**FAIL as registered (bar 1 on A119b); HG002 passes all four.** Prediction right on both samples.

| sample, track | arm | precision | recall | hits | false | hits lost | copies |
|---|---|---|---|---|---|---|---|
| A119b damped | baseline | 0.918 | 0.782 | 1,707 | 153 | - | 6,495 |
| A119b damped | F2 | 0.896 FAIL | 0.820 | 1,789 | 208 | 0 | 6,782 (1.044 x) |
| HG002 damped | baseline | 0.703 | 0.355 | 147 | 62 | - | 3,418 |
| HG002 damped | F2 | 0.718 | 0.498 | 206 | 81 | 0 | 3,540 (1.036 x) |
| A119b plain | baseline / F2 | 0.868 / 0.829 | 0.783 / 0.825 | 1,709 / 1,800 | 261 / 372 | 0 | 7,256 / 7,593 |
| HG002 plain | baseline / F2 | 0.545 / 0.583 | 0.396 / 0.541 | 164 / 224 | 137 / 160 | 0 | 3,927 / 4,060 |

Fa (attach without --keep-all), damped: A119b 0.903 / 0.800, HG002 0.730 / 0.490.
Added pairs (F2, damped): A119b 82 true / 55 false (0.60); HG002 59 true / 19 false (0.76). The exon-sharing rule did not
lift the added pairs' precision on A119b above what the coverage rule gave.
Miss classes under the exon-sharing default (pairdiag, damped): A119b baseline 476 missed = absorbed 71, split 205,
unlinked 200; F2 394 = absorbed 44, split 150, unlinked 200. HG002 baseline 267 = absorbed 42, split 91, unlinked 134;
F2 208 = absorbed 25, split 62, unlinked 121. Multiple representatives act on absorbed and split pairs and do not touch
the unlinked class on A119b.

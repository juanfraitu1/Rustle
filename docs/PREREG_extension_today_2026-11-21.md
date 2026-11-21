# Row 1883: transfrag extension at today's defaults, twenty datasets (pre-registered 2026-11-21)

Row 1811 (2026-11-19) switched extension off: +3.6% / +2.5% correct chains on A119b / OR6737, under 0.4% on five non-primate
datasets, `x.gtf` precision lower on all. Since then four defaults changed (junction ratio, spliced share, every read counted,
filter r1813 + truncation rule), and row 1882 found that with extension `final.gtf` now gains on BOTH axes on three ape testis
libraries and loses precision on mouse cortex. This row measures it everywhere. `RUSTLE_EXTEND_TRANSFRAGS=1`; nothing is built.

## Seen before registration
Row 1882's four datasets (Sumatran, Bornean, chimpanzee testis: both axes up; mouse cortex: 256 chains added, 0 annotated).

## Runs
The row 1880 job (`assemble` + `oneread`, same driver snapshot `3e35a0ed` and binaries) with `RUSTLE_EXTEND_TRANSFRAGS=1` on the
sixteen datasets row 1882 did not run; `off` = the row 1880 runs. Workspace `r1883_extension_today`.

## Measures (`final.gtf`, multi-exon chains)
chains / correct / precision off and on; NEW and LOST chains with their correct counts; the row 1880 comparison (list against
each rival) with extension on.

## Reading, fixed now
- Extension PAYS on a dataset when `final.gtf` has at least as many correct chains and at least the precision, one strictly
  better. HARM = more than 0.1% of the correct chains lost, or more than 0.001 of precision.
- **GENERAL** = pays on at least 6 of 7 non-primate and at least 11 of 13 primate datasets, no HARM anywhere.
- **PRIMATE-SIDE, HARMLESS ELSEWHERE** = pays on at least 11 of 13 primate datasets and no HARM on any non-primate one.
- **CONDITIONAL** = pays on at least 8 datasets and HARMs on at least one: then the table is printed with two library
  statistics computed without a reference from the `off` run, fixed now: (a) the share of `x.gtf` chains with 2 reads;
  (b) the median number of introns of `x.gtf` chains. A statistic is a LEAD for a gate only if one threshold puts every PAYS
  dataset on one side and every HARM dataset on the other; a gate would be its own row on libraries not used here.
- Otherwise: extension stays the opt-in it is, with the table.
- In every case the default is the user's call (the user switched it off in row 1811).

## Limits, stated in advance
- Real-library truths are incomplete: an extended chain that is real and unannotated counts as wrong.
- On the simulated sets extension has little to act on if full-length reads cover the truth (row 1811: under 0.4%).
- The four datasets of row 1882 were seen; they count in the tallies and are marked.

## Result (2026-11-22): CONDITIONAL by the letter, no gate; extension pays on the ten non-human ape samples and nowhere else

Jobs 4494713-4494728 (assembly stage exit 0 on all sixteen) plus row 1882's four; `r1883_extension_today/ext.NAME.txt` (`ext.py`).
`final.gtf` = chains / correct / precision; new and lost chains with their correct counts in brackets.
| dataset | group | off | with extension | correct | precision | new | lost | verdict | 2-read share | median introns |
|---|---|---|---|---|---|---|---|---|---|---|
| fly simulated | N | 6940 / 6753 / 0.9731 | 6955 / 6743 / 0.9695 | -10 (-0.15%) | -0.0035 | 31 (1) | 16 (11) | HARM | 0.087 | 3 |
| mouse simulated | N | 25864 / 25077 / 0.9696 | 26007 / 24626 / 0.9469 | -451 (-1.80%) | -0.0227 | 642 (10) | 499 (461) | HARM | 0.065 | 7 |
| zebrafish simulated | N | 13300 / 12129 / 0.9120 | 13501 / 11929 / 0.8836 | -200 (-1.65%) | -0.0284 | 448 (12) | 247 (212) | HARM | 0.095 | 5 |
| mouse ES | N | 41641 / 17073 / 0.4100 | 41817 / 17047 / 0.4077 | -26 (-0.15%) | -0.0023 | 767 (28) | 591 (54) | HARM | 0.180 | 8 |
| mouse heart | N | 21160 / 12922 / 0.6107 | 21239 / 12943 / 0.6094 | +21 (+0.16%) | -0.0013 | 309 (49) | 230 (28) | HARM | 0.217 | 8 |
| mouse adrenal | N | 23018 / 14162 / 0.6153 | 23129 / 14228 / 0.6152 | +66 (+0.47%) | -0.0001 | 463 (92) | 352 (26) | neither | 0.217 | 7 |
| mouse cortex (row 1882) | N | 29767 / 15490 / 0.5204 | 30023 / 15490 / 0.5159 | +0 (+0.00%) | -0.0044 | 659 (49) | 403 (49) | HARM | 0.227 | 9 |
| WTC11 r2 human iPS | P | 31336 / 20016 / 0.6388 | 31548 / 20113 / 0.6375 | +97 (+0.48%) | -0.0012 | 630 (150) | 418 (53) | HARM | 0.275 | 7 |
| A119b human testis | P | 147948 / 34882 / 0.2358 | 153648 / 35920 / 0.2338 | +1038 (+2.98%) | -0.0020 | 15795 (1510) | 10095 (472) | HARM | 0.230 | 6 |
| KB3781 gorilla fibroblast | P | 53153 / 23407 / 0.4404 | 53364 / 23710 / 0.4443 | +303 (+1.29%) | +0.0039 | 1879 (450) | 1668 (147) | PAYS | 0.189 | 8 |
| HG002 human | P | 94368 / 29150 / 0.3089 | 97885 / 29857 / 0.3050 | +707 (+2.43%) | -0.0039 | 7462 (917) | 3945 (210) | HARM | 0.183 | 5 |
| OR6737 gorilla testis | P | 57830 / 24315 / 0.4205 | 58306 / 24774 / 0.4249 | +459 (+1.89%) | +0.0044 | 2937 (622) | 2461 (163) | PAYS | 0.263 | 8 |
| Sumatran fibroblast | P | 51147 / 22378 / 0.4375 | 51349 / 22835 / 0.4447 | +457 (+2.04%) | +0.0072 | 1942 (542) | 1740 (85) | PAYS | 0.196 | 7 |
| Sumatran testis (row 1882) | P | 54748 / 15552 / 0.2841 | 56157 / 16327 / 0.2907 | +775 (+4.98%) | +0.0067 | 5254 (862) | 3845 (87) | PAYS | 0.327 | 4 |
| bonobo testis | P | 48548 / 17383 / 0.3581 | 49329 / 17734 / 0.3595 | +351 (+2.02%) | +0.0014 | 2728 (448) | 1947 (97) | PAYS | 0.327 | 9 |
| bonobo fibroblast | P | 69404 / 18951 / 0.2731 | 70011 / 20062 / 0.2866 | +1111 (+5.86%) | +0.0135 | 5472 (1228) | 4865 (117) | PAYS | 0.205 | 6 |
| Bornean fibroblast | P | 55156 / 21370 / 0.3874 | 55479 / 22068 / 0.3978 | +698 (+3.27%) | +0.0103 | 3132 (811) | 2809 (113) | PAYS | 0.208 | 7 |
| Bornean testis (row 1882) | P | 83473 / 19525 / 0.2339 | 85484 / 20478 / 0.2396 | +953 (+4.88%) | +0.0056 | 9255 (1195) | 7244 (242) | PAYS | 0.338 | 6 |
| chimpanzee cell line | P | 61218 / 26268 / 0.4291 | 61484 / 26776 / 0.4355 | +508 (+1.93%) | +0.0064 | 2747 (724) | 2481 (216) | PAYS | 0.203 | 7 |
| chimpanzee testis (row 1882) | P | 46969 / 20309 / 0.4324 | 47756 / 21106 / 0.4420 | +797 (+3.92%) | +0.0096 | 3455 (919) | 2668 (122) | PAYS | 0.341 | 8 |

- **Tally.** PAYS on 10 of 20, all ten of them the non-human ape samples (five testis, five cell lines; +1.3% to +5.9%
  correct chains, +0.001 to +0.014 precision). HARM on 9: the three simulated sets, mouse ES, heart, cortex, and all three human
  samples (WTC11, A119b, HG002). Mouse adrenal neither (+66 correct, -0.0001).
- **By the registered reading: CONDITIONAL** (pays on >= 8, harms on >= 1); not GENERAL, not harmless outside primates.
  Neither library statistic is a lead for a gate: the 2-read share of the PAYS datasets is 0.19-0.34 and of the HARM datasets
  0.07-0.28 (KB3781 pays at 0.189, WTC11 harms at 0.275); the median intron count overlaps entirely.
- **On complete truth extension destroys true chains.** Mouse simulated: 499 chains lost, 461 of them true, for 642 new ones of
  which 10 are true (-1.8% correct, -0.023 precision). Zebrafish 212 true lost for 12 gained, fly 11 for 1.
- **The split is not primate against non-primate, it is the annotation** (my reading). The three human samples, scored against
  the curated RefSeq set on CHM13, behave like mouse: more correct chains (+0.5% to +3.0%) at lower precision. The ten samples
  that pay are scored against RefSeq annotations of the non-human apes, which are in large part computed models. A chain made
  longer than any read matches a computed model more often than it matches a curated or a simulated one. This row cannot prove
  that reading; what it shows is that the gain exists only against those annotations, and that where the truth is known the
  extension is wrong 98% of the time it replaces a chain.
- Consequence for row 1880: StringTie's one lead (Sumatran testis) and its top-quarter leads on short-read testis libraries are
  on exactly this kind of truth. StringTie builds chains longer than its reads; the same caution applies to that lead.
- The row 1880 comparison with extension on was scored for row 1882's four datasets only; the other sixteen one-read runs were
  still finishing when this was written and add nothing to the verdict.

**Decision: none needed.** Extension stays off by default (row 1811, user). Do not re-open it on the strength of ape RefSeq
gains; a test that could re-open it needs complete truth or a curated annotation.

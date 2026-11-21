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

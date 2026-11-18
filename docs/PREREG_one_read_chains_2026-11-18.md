# Row 1805: can a one-read chain be accepted on evidence other than its read count? (pre-registered 2026-11-18)

User (2026-11-18): one read should be the bare minimum, without the flood of wrong chains a floor of one brings (row 1793:
A119b +3,509 correct, +130,800 wrong). Also the user's ruling: datasets already used are reused, not "spent"; a rule found here
is checked on one more unseen dataset. Registered before the script below was run. No assembly is rerun.

**Data.** The `floor1` and `base` arms of row 1793 (`PREFIX.gtf`): fly simulated, mouse simulated, mouse ES (+ SIRV, shown),
OR6737, A119b. **One-read chain** = a multi-exon transcript of `floor1` with `reads "1"` whose intron chain is not in `base`.
Truth and unit as in row 1793 (exact intron chains).

**Tests, reference-free** (evidence = the transcripts of `base`, i.e. chains with at least two reads):
- **A, introns proven:** every intron of the chain is an intron of some `base` transcript (same contig, strand).
- **B, ends proven:** the chain's first base and last base each lie within 50 bp of the same-side end of some `base` transcript on
  that strand.
- **A+B.**

**Bar, per test.** On a dataset the test is *usable* when the chains it accepts have precision at least that of `base`
`PREFIX.gtf` on that dataset (so adding them does not lower precision) and it keeps at least 30% of the correct one-read chains.
CEILING REACHED = usable on all five datasets; otherwise the datasets and the side that failed are named.
Reported besides: counts for "neither test"; SIRV; the share of the simulated one-read truth chains (mouse 7,949, fly 1,734)
that `floor1` emits at all, since the other rules also remove them.

**Limits.** A measurement of what is separable, not a rule: nothing is implemented and no default changes. On real data "wrong"
includes unannotated true chains. A passing test still needs an implementation and the check on an unseen dataset.

## Result (2026-11-18): CEILING NOT REACHED; neither test is usable on any dataset

Tables: workspace `r1805_one_read/{dme,mmu,es,sirv,or,a1}.txt` (`oneread.py`). Accepted one-read chains / correct / precision;
kept = share of the correct one-read chains.

| dataset (`base` precision) | all one-read | A introns proven | B ends proven | A+B | neither |
|---|---|---|---|---|---|
| fly sim (0.833) | 770 / 412 / 0.535 | 31 / 10 / 0.323, kept 0.02 | 360 / 116 / 0.322, kept 0.28 | 23 / 5 / 0.217 | 402 / 291 / 0.724, kept 0.71 |
| mouse sim (0.821) | 3,881 / 1,881 / 0.485 | 179 / 100 / 0.559, kept 0.05 | 1,488 / 266 / 0.179, kept 0.14 | 101 / 43 / 0.426 | 2,315 / 1,558 / 0.673, kept 0.83 |
| mouse ES (0.320) | 22,193 / 1,037 / 0.047 | 1,294 / 48 / 0.037 | 7,855 / 250 / 0.032 | 838 / 32 / 0.038 | 13,882 / 771 / 0.056 |
| OR6737 (0.347) | 57,638 / 3,317 / 0.058 | 2,664 / 362 / 0.136 | 13,305 / 1,058 / 0.080 | 1,459 / 260 / 0.178 | 43,128 / 2,157 / 0.050 |
| A119b (0.170) | 130,241 / 3,395 / 0.026 | 7,322 / 362 / 0.049 | 26,547 / 1,069 / 0.040 | 3,699 / 238 / 0.064 | 100,071 / 2,202 / 0.022 |

SIRV: 3 one-read chains, none correct.

**Reading.** Evidence borrowed from chains with two or more reads does not separate one-read chains: A accepts 2-11% of the correct
ones, and no accepted set reaches the precision of `base` anywhere. The correct one-read chains are mostly (65-83%) the ones that
pass *neither* test: a lone transcript of a locus with no other assembled chain, where nothing can vouch for it. On simulated
reads that class is 0.67-0.72 precise; on real libraries 0.02-0.06, so what makes a real lone one-read chain wrong is absent
from the simulation (and part of it is unannotated truth). One-read chains next to an expressed transcript (B) are the worse
class on the simulated sets (0.18-0.32): variants of a supported isoform.
Share of the simulated one-read truth chains that `floor1` emits: mouse 1,881 + 46 already in `base` of 7,949; fly 412 + 10 of
1,734. So three quarters of them are removed by the other rules or never assembled, whatever the floor.

# Row 1880: one ranked Rustle list against every other assembler, primates and non-primates (pre-registered 2026-11-21)

User goal (2026-11-21): rules that raise sensitivity and precision and outperform all other assemblers, not fitted to primates.
State at today's defaults (public comparison page, row 1849): `final.gtf` is more precise than StringTie and FLAIR on 12 of 12
ape samples, has more correct chains than StringTie on 9; StringTie has more on three testis libraries (chimpanzee, both
orangutans), IsoQuant and bambu are more precise on some samples at about half the correct chains. Row 1838 showed where
StringTie's extra chimpanzee chains come from: 77% have StringTie coverage under 2, i.e. about one read. Rustle keeps one-read
chains out of `final.gtf` and ranks them in a side file (`PREFIX.oneread.ranked.gtf`, row 1806; default product of `all`).

**The rule tested here:** Rustle's output read as ONE ranked list.
1. the chains of `final.gtf`, by `reads`, most first;
2. then the one-read chains of `PREFIX.oneread.ranked.gtf` that `final.gtf` lacks, by `oneread_score`, best first.
A cut of the list at N chains takes whole groups of equal reads; inside the group that N splits, the correct count is taken
in proportion (the expected value of a random order). Nothing is fitted in this row: both rankings exist and are frozen
(`reads`; `tools/data/oneread_r1806.json`).

Workspace `r1880_ranked_vs_rivals`. Registered before any run of this row.

## Seen before registration (development for this row)
Post hoc on existing files (`r1876_whats_left/tiercurve.py`, `strictcurve.py`): chimpanzee testis and WTC11 r2 with the one-read
tier against StringTie's point (reached at 8,594 and 270 added chains, precision 0.388 and 0.635 against StringTie's 0.350 and
0.543); the read-cut curve of `final.gtf` on chimpanzee testis, WTC11 r2, A119b, OR6737 against IsoQuant's and bambu's points
(at 3 reads Rustle is ahead of IsoQuant on both axes on all four; bambu needs 4 reads on A119b). These four are development;
the other sixteen datasets have not been looked at with the list.

## Data and runs
Twenty datasets. Per dataset a fresh `assemble` + `oneread` (no partner library) at today's defaults: driver =
`git show HEAD:tools/rustle_pipeline.sh` at this commit, binaries `~/rustle_target/release`, `RUSTLE_PTC=0`.
| group | datasets | truth |
|---|---|---|
| complete truth | fly, mouse, zebrafish simulated | rows 1779 / 1800 truths |
| non-primate real | mouse ES (mouse contigs), heart, adrenal, cortex | GENCODE vM27 |
| primate | the twelve ape samples of the comparison page; WTC11 r2 | RefSeq (the page's files; RefSeq on CHM13 for human) |
Rivals, de novo, as they are on disk (commands of `apetools/job.sbatch`): StringTie 3.0.3 `-L`; FLAIR 3.0.1 collapse; IsoQuant
4.1.0; bambu 3.12.1 at NDR 1. Where a rival has no output yet (StringTie and IsoQuant on zebrafish, heart, adrenal, cortex) it is
run with the same command. FLAIR and bambu are scored where they exist (FLAIR: apes, fly, mouse simulated, mouse ES; bambu: apes).

## Measure
Multi-exon intron chains (`bench/famcmp/curve.py load`). For a rival with N chains and C correct: Rustle's list cut at N, its
correct count R(N). If the list is shorter than N, R = the whole list's correct count.
- **AHEAD** of the rival on a dataset when R(N) >= C: at the rival's own output size Rustle reports at least as many correct
  chains, so it is at least as sensitive and at least as precise.
- For StringTie also along its own ranking (`cov`, top 25%, 50%, 75% of its chains): AHEAD at each of the three.

## Reading, fixed now
- Per rival: **OUTPERFORMED ACROSS THE BOARD** = AHEAD on every dataset where the rival has an output, in both groups. One miss
  = not across the board; the table says where, and by how many chains.
- **The list transfers** when StringTie and IsoQuant (the two rivals with outputs in both groups) are both outperformed on at
  least 6 of the 7 non-primate datasets and at least 11 of the 13 primate ones.
- Reported, no bar: the list cut at `final.gtf` plus the best 10% of the one-read tier (a fixed "sensitive" cut), and at 3 reads
  (a fixed "strict" cut): chains, correct, precision per dataset, beside each rival's point.
- Expectation written now: FLAIR is NOT outperformed on A119b and HG002 (its correct count there is beyond what Rustle's one-read
  tier adds); StringTie, IsoQuant and bambu are, with the simulated mouse set the closest call for IsoQuant.
- **Decision.** A table for the user. If the list transfers, the candidate is a product, not a filter: one ranked file (or a
  rank attribute) written by `all`, with the two fixed cuts named. Building it is the user's call.

## Limits, stated in advance
- A rival's point is its default; only StringTie is followed along its own ranking. IsoQuant has no knob; bambu's and FLAIR's
  other cuts are not scored here.
- Real-library truths are incomplete: a one-read chain that is real and unannotated counts as wrong, for every tool alike.
- The one-read tier costs a second assembly at a one-read floor (A119b + 1,628 s on 36 threads, row 1806).
- AHEAD at equal size is a statement about two lists, not about a default output: `final.gtf` itself is unchanged.

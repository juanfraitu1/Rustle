# Row 1842: near misses at splice sites, acceptor against donor (pre-registered 2026-11-20)

User (2026-11-20): look at the rules that judge near misses at the intron acceptor sites. Measurement only, GTF level; no
assembler code, no default. Workspace `r1842_acceptor_near` (`near.py`). Committed before the script was written or run.

## What exists (code reading, 2026-11-20)
No rule treats the acceptor apart from the donor. Every near-miss rule uses one window of 10 bp on both sides:
| rule | state | what it does |
|---|---|---|
| `--polish-junction-snap reads` | default (row 1800) | junction K -> neighbour J (donor and acceptor within 10 bp, same locus) when fewer than 2 of K's reads align better at K than at J |
| `--polish-junction-ratio R` | opt-in (rows 1819, 1835; 20 proposed) | K -> J within 10 bp when J has R times K's support |
| `--polish-fuzzy-junction`, `--read-junction-snap` | opt-in, failed (rows 866, 1527, 1618) | merges within a window, no evidence test |
| chain score features (wobble class, non-GT-AG count) | in the `final.gtf` filter | no side |
| support profile flags | failed as a filter (row 1830) | shifts of 1-3 bp, no side |
Known: 18-27% of the wrong chains of `final.gtf` have an intron within 10 bp of a true one (row 1818); Rustle keeps 4-7 times
StringTie's near-junction wrong chains on real libraries (row 1838); weaker shifted-site chains were measured on apes only,
without side or distance (rows 1725-1729: frame-keeping shifts are more often annotated); splice-site sequence strength
failed as a score feature (row 1789). Not measured anywhere: the split by side and by distance, on non-primates.

## Measured (`x.final.gtf` of the row 1839 `base` arms, twelve datasets, truths of row 1839; `x.gtf` beside it)
A PAIR = two emitted chains of one contig and strand with the same number of introns, equal except in ONE intron whose
donor and acceptor each differ by at most 10 bp. The WEAKER chain is the one with fewer `reads` (ties: both left out).
Cell of a pair = side (acceptor only, donor only, both) x distance class (1-2, 3, 4-5, 6 or 9, 7-8 or 10 bp; for `both` the
larger of the two). Per cell: weaker chains, the share of them that is correct, their share with reads ratio < 0.05 to the
stronger chain, and whether the stronger chain is correct. For acceptor-only pairs also: is the weaker chain's acceptor an AG.
Beside it, table E: the wrong chains that differ from a TRUTH chain in one intron by <= 10 bp, by the same cells (the errors,
whether or not the true partner was emitted).

## Reading, fixed now
- **The sides differ** when, on at least 9 of 12 datasets, the weaker chains of acceptor-only pairs are at least twice or at
  most half as precise as those of donor-only pairs (all distances pooled). Otherwise a rule has no reason to treat them apart.
- **A cell is a rule candidate** when it holds at least 50 weaker chains and they are at most a quarter as precise as the
  dataset's `final.gtf`, on at least 6 of 7 non-primate and 4 of 5 primate datasets. A candidate goes to the user as a
  proposal for a registered rule row (drop or re-place the weaker chain of that cell); nothing is built here.
- **A cell is protected** (a rule must spare it) when its weaker chains are at least as precise as the dataset's `final.gtf`
  on at least 6 of 7 and 4 of 5.

## Expectation, written before the script
Acceptor-only pairs outnumber donor-only ones and the 3 bp acceptor cell (NAGNAG) is the most precise one, likely protected on
the simulated sets. Distances 1-2 and 4-5 at the acceptor are the artifact cells; I am unsure they reach "a quarter" on real
libraries, since rows 1725-1727 found the weaker shifted chains only mildly worse than the read count predicts on human.

## Limits
- Pairs differ in exactly one intron; chains with two shifted sites are not seen. The default snap has already acted, so the
  pairs are what it left.
- Real-library truths are incomplete: an unannotated true alternative acceptor counts as wrong.
- Simulated reads: junction errors come from the simulator and the aligner. Twelve reused datasets, said.

## Result (2026-11-20): the sides do not differ; no cell is a rule candidate; the distance matters, the side does not

Job 4490045, twelve runs exit 0; `r1842_acceptor_near/near.NAME.txt`, cross-dataset view `agg.py`. `x.final.gtf` below;
`x.gtf` gives the same picture (in the files).

Weaker chain of a pair, correct / all. Order: fly, mouse, zebrafish simulated | mouse ES, heart, adrenal, cortex | OR6737,
A119b, chimpanzee, bonobo, KB3781.
| | fly | mouse sim | zebrafish | ES | heart | adrenal | cortex | OR6737 | A119b | chimp | bonobo | KB3781 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `final.gtf` precision | 0.964 | 0.953 | 0.891 | 0.347 | 0.577 | 0.578 | 0.487 | 0.378 | 0.219 | 0.415 | 0.337 | 0.373 |
| acceptor, 1-2 bp | 0/4 | 9/38 | 0/20 | 6/508 | 1/150 | 2/125 | 2/220 | 1/62 | 11/147 | 0/36 | 0/40 | 3/91 |
| acceptor, 3 bp | 22/28 | 319/373 | 16/26 | 441/2,837 | 246/832 | 269/915 | 318/1,425 | 1,114/2,230 | 1,398/3,971 | 677/1,184 | 717/1,481 | 1,336/2,702 |
| acceptor, 4-5 bp | 4/4 | 69/75 | 5/8 | 78/923 | 32/221 | 46/284 | 45/415 | 150/754 | 246/1,482 | 94/353 | 88/402 | 176/877 |
| acceptor, 6 or 9 bp | 8/9 | 44/46 | 5/5 | 37/433 | 22/103 | 21/111 | 29/163 | 80/346 | 109/578 | 36/135 | 40/181 | 83/361 |
| acceptor, all | 38/49 | 455/546 | 28/66 | 574/4,866 | 309/1,343 | 349/1,476 | 404/2,297 | 1,372/3,571 | 1,814/6,588 | 832/1,806 | 868/2,215 | 1,635/4,197 |
| donor, 1-2 bp | 0/8 | 6/34 | 1/18 | 8/254 | 4/52 | 3/51 | 2/57 | 4/90 | 16/212 | 3/53 | 3/57 | 6/80 |
| donor, 3 bp | 6/12 | 17/29 | 1/10 | 23/248 | 13/71 | 15/64 | 13/107 | 27/109 | 41/206 | 18/38 | 15/53 | 32/123 |
| donor, 4-5 bp | 6/25 | 75/125 | 9/45 | 74/791 | 31/166 | 45/190 | 45/254 | 149/594 | 252/1,375 | 92/317 | 74/375 | 164/714 |
| donor, all | 23/56 | 175/281 | 19/88 | 205/1,991 | 87/452 | 120/483 | 129/718 | 377/1,314 | 560/2,913 | 243/707 | 218/846 | 482/1,649 |
| both sides, all | 3/17 | 19/157 | 2/42 | 8/195 | 2/44 | 4/32 | 3/40 | 9/70 | 18/249 | 7/25 | 7/33 | 13/88 |

- **The sides do not differ** (registered reading): acceptor-only weaker chains are 0.9-2.0 times as precise as donor-only
  ones, a factor 2 on none of twelve (closest fly 0.78 / 0.41 and zebrafish 0.42 / 0.22). A rule has no reason to treat the
  acceptor apart. What differs is the number: on real libraries acceptor-only pairs are 2-3 times the donor-only ones, and
  more than half of them are the 3 bp cell.
- **No cell is a rule candidate.** The 1-2 bp cells are almost never correct on real libraries (acceptor 26 of 1,379 over the
  nine real libraries, donor 49 of 906) but hold under 50 chains on the three simulated sets and on chimpanzee and bonobo,
  and A119b is over the quarter bar (0.075 against 0.055): acceptor 4 of 7 and 2 of 5, donor 4 of 7 and 4 of 5.
- **No cell is protected across both groups.** The 3 bp acceptor cell (NAGNAG) is at least as precise as the assembly on 5 of
  5 primate datasets (0.35-0.57) and on 0 of 7 non-primate ones: 0.62-0.86 on the simulated sets (below their 0.89-0.96) and
  0.16-0.30 on the mouse libraries (about half the assembly's precision).
- Expectation: acceptor pairs outnumber donor pairs, yes; the 3 bp acceptor cell is the most precise acceptor cell on all nine real
  libraries (on the simulated sets 4-5 and 6 / 9 bp are above it); "likely protected on the simulated
  sets" was wrong. The 1-2 bp cells are far under a quarter wherever they have the numbers; the 4-5 bp cells are not artifact
  cells (0.08-0.92, near the pair average).
- The weaker chain's acceptor is an AG in 99.5-100% of the acceptor-only pairs: the assembly already holds no non-AG variant,
  so the dinucleotide cannot judge them.
- Table E (errors against the truth): a wrong chain one near intron away from a truth chain is 26%, 30%, 16% of the wrong
  chains on the simulated sets, 13-14% on the mouse libraries, 2-6% on the apes; the true partner is emitted in 81-94% of
  them on real libraries (the default snap could have moved them and did not: at least 2 reads prove the weaker site).
  Simulated only: the largest error cell on mouse simulated is an intron moved at BOTH sides by 4-5 bp (127 of 396, partner
  emitted in 108), 9 on fly, 35 on zebrafish, 60 on mouse ES.

**Reading (mine).** The near misses that are left are judged by distance, not by side. Three things could be made rules, all
small: (1) a weaker chain 1-2 bp from a stronger one, either side: 0.02-0.05 correct on real libraries, 2.2% of the wrong
chains of mouse ES, 0.2% of A119b's; (2) an intron moved at both sides, 0.04-0.15 correct, a few hundred chains at most;
(3) the 3 bp acceptor cell must be spared on apes and simulated data and is the bulk of the pairs. None is acceptor-specific
and none is primate-overfit; the default snap's 10 bp window treats 1-2 bp like 9 bp, and a distance-aware version is the only
change this table supports. Each would be its own registered row; the user decides.

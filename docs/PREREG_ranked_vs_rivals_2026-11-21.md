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

## Result (2026-11-21): the list transfers; IsoQuant and bambu outperformed across the board, StringTie on 19 of 20, FLAIR on 14 of 15

Jobs 4494627-4494646 (Rustle, every stage exit 0) and 4494647-4494650 (StringTie and IsoQuant on zebrafish, heart, adrenal,
cortex); driver `3e35a0ed`, binaries of 2026-11-20; `r1880_ranked_vs_rivals/score.NAME.txt` (`rank.py`), tables `tables.md`.
`final.gtf` equals row 1849's default numbers on the twelve datasets both rows ran.

Each rival cell: the rival's correct chains of its chains -> Rustle's correct chains at that size (difference).
`final.gtf` = chains / correct / precision. N = non-primate, P = primate.
| dataset | group | `final.gtf` | StringTie | IsoQuant | FLAIR | bambu |
|---|---|---|---|---|---|---|
| fly simulated | N | 6940 / 6753 / 0.9731 | 6,549 of 6,801 -> 6,622 (+73) | 6,592 of 6,780 -> 6,603 (+11) | 6,614 of 9,183 -> 8,101 (+1487) | - |
| mouse simulated | N | 25864 / 25077 / 0.9696 | 23,705 of 24,730 -> 24,026 (+321) | 24,098 of 24,924 -> 24,206 (+108) | 27,048 of 50,555 -> 29,095 (+2047) | - |
| zebrafish simulated | N | 13300 / 12129 / 0.9120 | 11,826 of 13,270 -> 12,104 (+278) | 11,503 of 12,571 -> 11,518 (+15) | - | - |
| mouse ES | N | 41641 / 17073 / 0.4100 | 16,401 of 44,591 -> 17,404 (+1003) | 13,443 of 28,508 -> 15,039 (+1596) | 16,686 of 172,968 -> 18,872 (+2186) | - |
| mouse heart | N | 21160 / 12922 / 0.6107 | 12,577 of 22,478 -> 13,121 (+544) | 11,755 of 18,760 -> 12,170 (+415) | - | - |
| mouse adrenal | N | 23018 / 14162 / 0.6153 | 13,741 of 24,761 -> 14,388 (+647) | 12,880 of 20,415 -> 13,369 (+489) | - | - |
| mouse cortex | N | 29767 / 15490 / 0.5204 | 14,935 of 34,611 -> 16,207 (+1272) | 13,416 of 26,060 -> 14,672 (+1256) | - | - |
| WTC11 r2 | P | 31336 / 20016 / 0.6388 | 20,064 of 36,963 -> 21,134 (+1070) | 16,609 of 26,191 -> 18,019 (+1410) | - | - |
| A119b human testis | P | 147948 / 34882 / 0.2358 | 31,781 of 211,865 -> 36,480 (+4699) | 18,903 of 84,539 -> 28,766 (+9863) | 39,230 of 626,983 -> 41,662 (+2432) | 12,793 of 41,745 -> 20,187 (+7394) |
| KB3781 gorilla fibroblast | P | 53153 / 23407 / 0.4404 | 22,087 of 58,164 -> 23,910 (+1823) | 15,785 of 36,929 -> 20,031 (+4246) | 24,734 of 197,662 -> 26,771 (+2037) | 14,232 of 26,430 -> 16,546 (+2314) |
| HG002 human | P | 94368 / 29150 / 0.3089 | 27,456 of 115,579 -> 30,146 (+2690) | 16,087 of 67,703 -> 26,859 (+10772) | 38,489 of 728,951 -> 32,318 (-6171) | 11,073 of 25,126 -> 16,167 (+5094) |
| OR6737 gorilla testis | P | 57830 / 24315 / 0.4205 | 23,248 of 67,491 -> 25,408 (+2160) | 16,677 of 37,089 -> 19,830 (+3153) | 24,293 of 138,125 -> 29,369 (+5076) | 15,381 of 31,947 -> 18,236 (+2855) |
| Sumatran fibroblast | P | 51147 / 22378 / 0.4375 | 21,510 of 55,400 -> 22,819 (+1309) | 15,667 of 36,794 -> 19,450 (+3783) | 23,993 of 183,790 -> 26,261 (+2268) | 14,021 of 27,270 -> 16,514 (+2493) |
| Sumatran testis | P | 54748 / 15552 / 0.2841 | 16,981 of 72,495 -> 16,295 (-686) | 12,866 of 41,631 -> 13,531 (+665) | 16,513 of 111,097 -> 18,297 (+1784) | 10,451 of 37,259 -> 12,858 (+2407) |
| bonobo testis | P | 48548 / 17383 / 0.3581 | 17,288 of 61,623 -> 18,643 (+1355) | 13,109 of 35,490 -> 15,045 (+1936) | 15,934 of 69,488 -> 19,240 (+3306) | 12,429 of 29,640 -> 13,771 (+1342) |
| bonobo fibroblast | P | 69404 / 18951 / 0.2731 | 18,780 of 96,488 -> 20,451 (+1671) | 13,139 of 60,515 -> 17,993 (+4854) | 21,171 of 337,026 -> 22,866 (+1695) | 9,128 of 25,929 -> 11,879 (+2751) |
| Bornean fibroblast | P | 55156 / 21370 / 0.3874 | 20,787 of 66,155 -> 22,164 (+1377) | 14,803 of 42,147 -> 19,297 (+4494) | 23,246 of 231,157 -> 25,303 (+2057) | 12,203 of 25,866 -> 14,985 (+2782) |
| Bornean testis | P | 83473 / 19525 / 0.2339 | 20,561 of 123,272 -> 21,227 (+666) | 14,094 of 57,528 -> 16,810 (+2716) | 19,881 of 181,913 -> 23,052 (+3171) | 10,653 of 36,021 -> 13,396 (+2743) |
| chimpanzee cell line | P | 61218 / 26268 / 0.4291 | 25,433 of 69,534 -> 27,076 (+1643) | 16,579 of 41,003 -> 21,872 (+5293) | 30,658 of 235,531 -> 31,510 (+852) | 13,650 of 25,649 -> 16,407 (+2757) |
| chimpanzee testis | P | 46969 / 20309 / 0.4324 | 21,532 of 61,509 -> 22,310 (+778) | 15,761 of 35,718 -> 17,487 (+1726) | 19,730 of 72,154 -> 23,501 (+3771) | 13,835 of 30,781 -> 16,205 (+2370) |

| rival | datasets with an output | AHEAD, non-primate | AHEAD, primate | verdict |
|---|---|---|---|---|
| StringTie | 20 | 7 of 7 | 12 of 13 | not across the board: Sumatran testis, 686 chains behind |
| IsoQuant | 20 | 7 of 7 | 13 of 13 | **OUTPERFORMED ACROSS THE BOARD** |
| bambu (NDR 1) | 12 | - | 12 of 12 | **OUTPERFORMED ACROSS THE BOARD** (primates only: no non-primate output) |
| FLAIR | 15 | 3 of 3 | 11 of 12 | not across the board: HG002, 6,171 chains behind |

- **The list transfers** (registered bar: StringTie and IsoQuant both outperformed on >= 6 of 7 non-primate and >= 11 of 13
  primate datasets): 7 of 7 and 12 of 13.
- **The one StringTie miss** is the Sumatran orangutan testis library, and it is a miss along StringTie's whole ranking (-673,
  -460, -62 at 25 / 50 / 75%). It is one of the two short-read orangutan testis libraries (the other, Bornean, is ahead by 666
  at StringTie's size and behind by 235 at its top quarter). Chimpanzee testis and the Bornean fibroblast line are behind at
  StringTie's top quarter only (-99, -35) and ahead from there on.
- **FLAIR** is outperformed on A119b against the expectation written above (+2,432: the whole list has more correct chains than
  FLAIR's 626,983 chains hold) and not on HG002, as expected: FLAIR reports 38,489 correct chains there in 728,951; Rustle's whole
  list holds 32,318 in 318,884. Where FLAIR has more chains than Rustle's list (9 of 15 datasets) "AHEAD" means more correct
  chains in a shorter list.
- The closest calls for IsoQuant are the simulated sets, as expected: fly +11, zebrafish +15, mouse +108.

Along StringTie's own ranking (difference in correct chains at its top 25%, 50%, 75%), and the two fixed cuts
(chains / correct / precision):
| dataset | StringTie top 25% | top 50% | top 75% | sensitive cut | strict cut |
|---|---|---|---|---|---|
| fly simulated | +14 | +50 | +65 | 7,132 / 6,933 / 0.9721 | 6,336 / 6,185 / 0.9762 |
| mouse simulated | +25 | +96 | +212 | 26,373 / 25,493 / 0.9666 | 24,191 / 23,527 / 0.9726 |
| zebrafish simulated | +65 | +141 | +228 | 13,720 / 12,455 / 0.9078 | 12,009 / 11,047 / 0.9199 |
| mouse ES | +314 | +617 | +747 | 45,683 / 17,516 / 0.3834 | 34,554 / 16,131 / 0.4668 |
| mouse heart | +111 | +230 | +343 | 23,598 / 13,322 / 0.5645 | 16,747 / 11,539 / 0.6890 |
| mouse adrenal | +145 | +182 | +424 | 25,716 / 14,566 / 0.5664 | 18,248 / 12,709 / 0.6965 |
| mouse cortex | +349 | +664 | +882 | 33,423 / 16,054 / 0.4803 | 23,244 / 14,051 / 0.6045 |
| WTC11 r2 | +32 | +307 | +714 | 36,611 / 21,068 / 0.5755 | 23,110 / 16,823 / 0.7280 |
| A119b human testis | +3140 | +4797 | +5460 | 189,978 / 36,035 / 0.1897 | 111,463 / 31,887 / 0.2861 |
| KB3781 gorilla fibroblast | +147 | +712 | +1395 | 60,000 / 24,074 / 0.4012 | 43,433 / 21,647 / 0.4984 |
| HG002 human | +943 | +2296 | +2468 | 116,819 / 30,176 / 0.2583 | 76,462 / 27,801 / 0.3636 |
| OR6737 gorilla testis | +267 | +1075 | +1609 | 70,201 / 25,649 / 0.3654 | 42,779 / 21,342 / 0.4989 |
| Sumatran fibroblast | +15 | +286 | +798 | 58,611 / 23,096 / 0.3941 | 41,410 / 20,597 / 0.4974 |
| Sumatran testis | -673 | -460 | -62 | 71,899 / 16,272 / 0.2263 | 36,220 / 12,698 / 0.3506 |
| bonobo testis | +208 | +580 | +995 | 61,357 / 18,628 / 0.3036 | 33,081 / 14,614 / 0.4418 |
| bonobo fibroblast | +331 | +1297 | +1742 | 86,501 / 20,058 / 0.2319 | 54,540 / 17,349 / 0.3181 |
| Bornean fibroblast | -35 | +886 | +1349 | 65,338 / 22,107 / 0.3383 | 43,823 / 19,642 / 0.4482 |
| Bornean testis | -235 | +326 | +842 | 111,942 / 20,828 / 0.1861 | 55,008 / 16,546 / 0.3008 |
| chimpanzee cell line | +249 | +964 | +1351 | 71,517 / 27,216 / 0.3806 | 49,284 / 23,991 / 0.4868 |
| chimpanzee testis | -99 | +313 | +706 | 60,223 / 22,158 / 0.3679 | 31,255 / 16,367 / 0.5237 |


**The two fixed cuts (reported, no bar).**
- SENSITIVE (`final.gtf` + the best 10% of the one-read tier) has at least StringTie's correct chains AND at least its precision
  on 19 of 20 datasets (not Sumatran testis). This is one fixed cut, not the rival's size.
- STRICT (3 reads) is more precise than StringTie on 20 of 20, but as a fixed cut it does not dominate the precise tools: it
  has both more correct chains and higher precision than IsoQuant on 13 of 20 (on the seven misses it is more precise with
  fewer correct chains, or the reverse on Sumatran testis) and than bambu on 5 of 12 (on the other seven bambu is more precise
  at about 60% of the correct chains). Against those two only the list cut at their size is ahead everywhere.

**What this is and is not (my reading).** No transcript changed: `final.gtf` is the same file. The finding is that Rustle's
two existing rankings, joined, give one list that at any other tool's output size holds more annotated chains than that tool,
on 65 of 67 comparisons, equally in mouse, fly, zebrafish and the apes. The join matters at the long end: without the one-read
tier Rustle cannot reach StringTie's size at all on any real library. The one-read ranking was fitted on seven datasets (row 1806: fly, mouse and
zebrafish simulated, mouse ES, heart, OR6737, A119b) and is used here on thirteen it never saw.
Limits as registered; in addition the ranking inside `final.gtf` is the raw read count, so the list's top can trail a
coverage ranking by a few dozen chains on short-read libraries (the three top-quarter misses).

**Candidate, the user's call:** write the list as a product of `all` (one ranked GTF or a `rank` attribute plus the tier of each
transcript), with SENSITIVE as a named cut. Not proposed: STRICT at a fixed read count.

## Post hoc (2026-11-22, not registered): the Sumatran testis miss depends on the annotation; row 1884 closes the FLAIR miss

- Rescored against the CAT annotation of each species (`genomes/TAG/cat.gff3` converted with gffread; `rank.py` with
  `RANK_TRUTH`, files `r1880_ranked_vs_rivals/cat.score.NAME.txt`). Difference in correct chains, Rustle minus StringTie:
  | testis sample | at StringTie's size, RefSeq -> CAT | top 25% | top 50% | top 75% |
  |---|---|---|---|---|
  | Sumatran | -686 -> **+237** | -673 -> -493 | -460 -> -163 | -62 -> +267 |
  | Bornean | +666 -> +550 | -235 -> -69 | +326 -> +288 | +842 -> +651 |
  | chimpanzee | +778 -> +669 | -99 -> +40 | +313 -> +320 | +706 -> +573 |
  Against CAT the list is ahead of StringTie at StringTie's size on all three; the Sumatran sample stays behind inside
  StringTie's top half. Every other rival stays AHEAD on the three samples under CAT. One annotation against another: this
  shows the single lost comparison is not robust to the choice of truth, not that it is won.
- Row 1884 (`docs/PREREG_list_tail_2026-11-22.md`) appends a third tier and closes the FLAIR miss on HG002 (-6,171 -> +6,718).
  Tally with it, against RefSeq: IsoQuant 20 of 20, bambu 12 of 12, FLAIR 15 of 15, StringTie 19 of 20.

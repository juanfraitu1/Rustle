# Row 1808: one-read chains seen in a second library and ranked high, against the two-read chains (pre-registered 2026-11-18)

User (2026-11-18), after row 1807: take the bar "as good as the chains the floor already lets in" and check the candidate rule
on libraries not used in row 1807. Registered before any download was mapped or any run started. No default changes.

**Rule under test (fixed, from row 1807).** Promote a one-read chain of `PREFIX.oneread.ranked.gtf` when (T1) its exact intron
chain is carried by at least one primary alignment in another library of the same species and (T3) it is in the best 20% of the
ranking (`oneread_rank`, model `tools/data/oneread_r1806.json`, unchanged).

**Libraries** (none used in row 1807 or in the ranking's fit; partners may be used ones):
| target | partner |
|---|---|
| Sumatran orangutan testis (`reads/PAB/PAB_mm.bam`, RefSeq) | its fibroblast line (`PAB_AG06213_fib_mm.bam`) |
| Bornean orangutan testis (`reads/PPY/PPY_mm.bam`, RefSeq) | its fibroblast line (`PPY_AG05252_fib_mm.bam`) |
| mouse adrenal gland, ENCODE ENCSR852HRM (ENCFF429RCV + ENCFF387GTW pooled, 1,942,036 reads; GENCODE vM27) | mouse heart, cortex, ES cells |

Shipped mapping for the new library; `assemble` at today's defaults; stage `oneread`.

**Bar, per target.** The promoted chains are at least as precise as the default's chains with exactly two reads, and their correct
chains number at least 2% of the default's correct chains. RULE CONFIRMED = all three targets.
Reported: the default by read count (2, 3, 4-5, >= 6); T1 and T3 alone; gained chains only (nothing is removed by promotion).

**Limits.** Needs a second library of the species: never a single-sample default. The two orangutan pairs are one individual
each (testis and fibroblast of different animals would be a harder test). Annotation as truth; two-read chains are judged by the
same incomplete annotation, so the comparison is like for like. The strict tier's rules were tuned with orangutan testis data
(row 1766 default test); the rule under test does not use them.

## Result (2026-11-18): RULE CONFIRMED, PASS on all three libraries

Jobs 4482760 (Sumatran orangutan), 4482761 (Bornean orangutan), 4482763 (mouse adrenal; mapping 4482762), every stage exit 0.
Tables: workspace `r1808_one_read_promote_check/{pab,ppy,adrenal}.txt`. Chains / correct / precision.

| | Sumatran orangutan testis | Bornean orangutan testis | mouse adrenal |
|---|---|---|---|
| default | 75,238 / 16,928 / 0.225 | 112,924 / 21,346 / 0.189 | 26,160 / 14,468 / 0.553 |
| default, 2 reads | 28,140 / 0.133 | 43,518 / 0.091 | 6,206 / 0.258 |
| default, 3 reads | 12,438 / 0.172 | 18,347 / 0.128 | 3,089 / 0.337 |
| default, 4-5 reads | 11,502 / 0.212 | 17,080 / 0.171 | 3,355 / 0.439 |
| default, >= 6 reads | 23,158 / 0.372 | 33,979 / 0.356 | 13,510 / 0.766 |
| **promoted (T1+T3)** | **1,773 / 492 / 0.277** | **1,821 / 487 / 0.267** | **1,446 / 444 / 0.307** |
| share of the default's correct chains | 2.9% | 2.3% | 3.1% |
| T1 alone (second library) | 6,467 / 1,603 / 0.248 | 6,304 / 1,203 / 0.191 | 5,991 / 1,290 / 0.215 |
| T3 alone (best 20%) | 14,708 / 1,202 / 0.082 | 19,711 / 1,018 / 0.052 | 2,877 / 556 / 0.193 |

The promoted chains are 2.1, 2.9 and 1.2 times as precise as the two-read chains, and sit between the 3-read and 4-5-read chains
on mouse adrenal and above the 4-5-read chains on both orangutans. T1 alone also clears the two-read bar on the orangutans
(yield 5.6-9.5%) but not on mouse adrenal (0.215 against 0.258). Effect on `PREFIX.gtf` + promoted: Sumatran 0.225 -> 0.226,
Bornean 0.189 -> 0.190, mouse adrenal 0.553 -> 0.540.

**Implemented, opt-in.** Stage `oneread` with `RUSTLE_ONEREAD_PARTNERS=A.bam:B.bam` (other libraries of the same species) writes
`PREFIX.oneread.promoted.gtf` (`oneread_promoted "1"`; `bench/famcmp/oneread_tier.py promote`); with `RUSTLE_ONEREAD=1` it is
appended to `PREFIX.final.gtf` after the strict tier. `PREFIX.gtf` is never changed.

# Row 1807: which one-read chains can be promoted into the main output? (pre-registered 2026-11-18)

User (2026-11-18): find a way to promote one-read chains. Rows 1805 / 1806: neither borrowed evidence nor the score brings them to
the precision of the default output. Registered before the script ran. Measurement on existing `oneread` outputs; no default changes.

**Data** (stage `oneread` outputs of row 1806): mouse cortex (outside the score's fit), chimpanzee testis (outside), mouse heart
(inside the fit, so its T3 figures are in-sample). One-read chain = a chain of `PREFIX.oneread.ranked.gtf`.

**Tests.**
- **T1, a second library:** the exact intron chain is carried by at least one primary alignment in another library of the same
  species: cortex <- heart + ES cells; heart <- cortex + ES cells; chimpanzee testis <- chimpanzee LCL (`PTR_AG18354_lcl_mm.bam`).
  Two molecules from two libraries instead of two from one.
- **T2, the strict tier** (`PREFIX.oneread.gtf`, row 1766 rule S3b).
- **T3, the best 20% of the ranking** (`oneread_rank`).
- T1+T2, T1+T3, T2+T3, all three.

**Bar.** A test is PROMOTABLE on a library when the chains it accepts are at least as precise as the default `PREFIX.gtf` of that
library and number at least 1% of the default's correct chains. A WAY TO PROMOTE = promotable on all three libraries.
Reported: accepted, correct, precision for every test; for T1 the count by partner library.

**Limits.** T1 needs a second library of the species, so it is not a single-sample rule. A recurrent alignment artefact repeats
across libraries. Truth is an annotation: chains seen in two libraries and still unannotated count against the test.

## Result (2026-11-18): NO WAY TO PROMOTE as registered; promotable on chimpanzee only

Job 4482741 (388 s). Tables: workspace `r1807_one_read_promote/{cortex,ptr,heart}.txt`. Accepted / correct / precision.

| test | mouse cortex (default 0.467) | chimp testis (0.386) | mouse heart (0.557) |
|---|---|---|---|
| T1 second library | 6,149 / 1,284 / 0.209 | 7,926 / 2,325 / 0.293 | 6,174 / 1,329 / 0.215 |
| T2 strict tier | 872 / 181 / 0.208 | 1,365 / 491 / 0.360 | 640 / 163 / 0.255 |
| T3 best 20% of ranking | 4,155 / 726 / 0.175 | 10,576 / 1,925 / 0.182 | 2,973 / 674 / 0.227 |
| T1+T2 | 449 / 142 / 0.316 | 454 / 285 / **0.628** | 420 / 136 / 0.324 |
| T1+T3 | 1,699 / 554 / 0.326 | 2,340 / 944 / **0.403** | 1,468 / 518 / 0.353 |
| T2+T3 | 460 / 123 / 0.267 | 925 / 384 / **0.415** | 297 / 89 / 0.300 |
| all three | 265 / 100 / 0.377 | 325 / 222 / **0.683** | 199 / 72 / 0.362 |

Bold = promotable. T1 by partner: cortex <- heart 2,407 (743 correct), ES 5,335 (1,111); heart <- cortex 3,362 (940), ES 5,114
(1,122); chimpanzee <- LCL 7,926 (2,325). No test reaches the default's precision on either mouse library.

**Post hoc, not the registered bar: the chains the floor already admits.** Precision of the default chains by read count:

| reads | mouse cortex | chimp testis | mouse heart |
|---|---|---|---|
| 2 | 8,507 chains, 0.181 | 21,332, 0.227 | 5,507, 0.275 |
| 3 | 4,421, 0.273 | 8,912, 0.311 | 2,841, 0.352 |
| 4-5 | 4,487, 0.349 | 8,377, 0.379 | 3,121, 0.453 |
| >= 6 | 16,212, 0.702 | 18,133, 0.612 | 12,136, 0.760 |

The default's overall precision is carried by its well-supported chains; its two-read chains match the annotation 0.18-0.28.
Against those, T1+T3 (second library and best 20% of the ranking) is higher on all three libraries (0.326 / 0.403 / 0.353, about
the level of the default's 3- to 5-read chains) with 554 / 944 / 518 correct chains (3.5% / 4.3% / 3.9% of the default's);
T1 alone is higher on cortex and chimpanzee, lower on heart. Row 1766 used this two-read comparison as its bar.

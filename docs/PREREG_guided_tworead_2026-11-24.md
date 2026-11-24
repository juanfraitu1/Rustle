# PREREG: guided mode, chains of guide introns at a two-read floor (row 1939, 2026-11-24)

## Current state
Scored 2026-11-24: FAIL (4 of 10 judging samples pass; 9 needed). Outcome at the end.

User (2026-11-24), after row 1934: "lets try arm B with two reads". Row 1934's arm B (one-read chains made only of guide
introns) added 13-40% hidden transcripts but failed the precision bar on 10 of 10 samples. Written before any hidden share
of the two-read class was computed. Seen so far, on the development sample OR6737 only: of the 3,963 chains of arm B's
class 3,388 have one read and 575 have two or more.

## Rule
Arm B2 = arm A of row 1934 (today's guided mode since the port) + every chain of the floor-1 assembly
(`PREFIX.oneread.d/f1.gtf`) that `PREFIX.gtf` lacks, whose introns are ALL guide introns, that is neither a contiguous run
nor an intron-retaining version of a guide chain, and that has `reads` >= 2 in the floor-1 assembly.
Reported beside, no bar: the same at >= 3 reads.

## Data and decision (as row 1934; no new assembly, the runs of `r1934_guided_final_oneread/runs/` are rescored)
Development: OR6737, A119b. Judging: the other ten samples. Truth = hidden intron chains of the guide split.
- B1: hidden share of the added chains >= arm A's precision outside the guide.
- B2: hidden added >= 2% of arm A's hidden found.
- PASS = both bars on >= 9 of the 10 judging samples. New and lost chains separately (nothing is lost by construction).

## Predictions (mine)
The added chains are far fewer than at one read (hundreds per sample) and more precise. B1 holds on 7 or more of the 10,
I cannot say whether on 9. B2 holds on most but is close to its floor on the samples with a large arm A. Overall: uncertain,
leaning FAIL by one or two samples.

## Limits
One guide split per species (the hidden 31%); the floor-1 assembly costs a second assembly run (within 2 x `assemble`, row
1934). A PASS makes a candidate; adding it to the guided product is the user's call.

## Outcome (2026-11-24; workspace `r1939_guided_tworead/{score2.py,all.txt}`)
OR6737 and A119b are development; the other ten judge.

| sample | arm A: hidden / precision | B2 (>= 2 reads): added (hidden, share) | >= 3 reads, no bar |
|---|---|---|---|
| A119b | 7722 / 0.1588 | 2773 (415, 0.1497) FAIL | 1926 (332, 0.1724) PASS |
| OR6737 | 5421 / 0.3923 | 575 (181, 0.3148) FAIL | 369 (117, 0.3171) FAIL |
| HG002_kinnex.t2t | 7207 / 0.2011 | 4175 (843, 0.2019) PASS | 3206 (734, 0.2289) PASS |
| GGO_KB3781_fib_mm | 5822 / 0.3575 | 853 (289, 0.3388) FAIL | 617 (230, 0.3728) PASS |
| PTR_mm | 3894 / 0.4604 | 242 (125, 0.5165) PASS | 111 (56, 0.5045) FAIL |
| PTR_AG18354_lcl_mm | 6362 / 0.3671 | 1218 (606, 0.4975) PASS | 850 (440, 0.5176) PASS |
| PPA_mm | 3441 / 0.3263 | 212 (79, 0.3726) PASS | 103 (48, 0.4660) FAIL |
| PPA_PR00251_fib_mm | 4291 / 0.2286 | 807 (182, 0.2255) FAIL | 483 (112, 0.2319) PASS |
| PPY_mm | 3947 / 0.1966 | 251 (48, 0.1912) FAIL | 134 (24, 0.1791) FAIL |
| PPY_AG05252_fib_mm | 5132 / 0.3364 | 794 (243, 0.3060) FAIL | 555 (182, 0.3279) FAIL |
| PAB_mm | 3047 / 0.2202 | 180 (43, 0.2389) FAIL | 68 (16, 0.2353) FAIL |
| PAB_AG06213_fib_mm | 5451 / 0.3600 | 790 (272, 0.3443) FAIL | 544 (195, 0.3585) FAIL |

- **FAIL: both bars on 4 of the 10 judging samples** (HG002, chimp testis, chimp cell line, bonobo testis). B1 (precision
  of the added chains) holds on 5 of 10, B2 (yield) on 8 of 10.
- What the two-read floor does: the added chains are now about as precise as arm A's own chains outside the guide (share
  minus arm A's precision between -0.030 and +0.130, 7 of 10 within 0.025), so the overall precision moves by at most
  0.009; the yield is +1.2% to +11.7% hidden transcripts (one read: +13% to +40%).
- At three reads (no bar): 5 of 10 would pass; the two testis samples that pass at two reads fail the yield bar at three.
- Predictions: "B1 on 7 or more" WRONG (5); "B2 on most, close to the floor" held (8, two below); overall FAIL as leaned.
- Reading: the class has no floor at which it is both clearly more precise than the guided assembly and large enough to
  matter on every sample. At two reads it is precision-neutral and small. Not added to the guided product.

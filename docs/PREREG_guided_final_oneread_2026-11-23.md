# PREREG r1934: guided mode on final.gtf, and one-read chains made only of guide introns (2026-11-23)

User (2026-11-23): "make the guided mode a bit better ... we are not really taking full advantage of guided mode"; approved
registering the two arms below. Written before any run of this row. Workspace `r1934_guided_final_oneread`, scorer
`bench/famcmp/guided_arms.py`, job `run.sbatch` there.

## Current state
Run 2026-11-23/24 (jobs 4498896, 4498934-4498944). Arm A PASSES on 10 of 10 judging samples (candidate; the default is the user\'s call). Arm B FAILS on 10 of 10 (precision bar B1). Outcome at the end.

## What is known before (step 0, read-only, runs of `artifact_refresh_2026-11-21`, OR6737 and A119b only; workspace `guided_step0_2026-11-23`)
- Guided mode (rule G1d, row 1622) filters `PREFIX.gtf` and adds guide chains with a full-chain read, so it cannot find a hidden
  transcript the de novo assembly lacks, and it does not see the chain score filter (`PREFIX.final.gtf`).
- The same verdicts on `final.gtf`: A119b 7,705 hidden at 0.158 (today 7,776 at 0.137); OR6737 5,422 at 0.393 (5,429 at 0.371).
- Non-guide chains whose introns are all guide introns in a new combination are the best class: 3,928 hidden of 10,957 on A119b
  (0.36), 2,467 of 4,666 on OR6737 (0.53); chains with a new intron beside guide introns 0.07 / 0.14.
- Gated one-read files of other runs (not the pool of this row): that class at 25 of 76 (A119b) and 123 of 278 (OR6737).
Both samples are therefore development for this row. KB3781 was the test sample of the G1d floors (row 1622) but has not been
looked at for these arms.

## Arms (frozen; no parameter is chosen from the results)
- **Arm A:** G1d's verdicts applied to the chains of `PREFIX.final.gtf` instead of `PREFIX.gtf`. Floors unchanged (10 reads with
  a new intron beside guide introns, 3 reads otherwise; fragments and intron-retaining versions of guide chains dropped). Guide
  chains with a full-chain primary alignment as today.
- **Arm B:** arm A + every chain of the floor-1 assembly (`PREFIX.oneread.d/f1.gtf`, the `oneread` stage's own run at today's
  defaults) that `PREFIX.gtf` lacks, whose introns are ALL guide introns, and that is neither a contiguous run nor an
  intron-retaining version of a guide chain. No read floor, no score: the guide confirms every junction, a read the combination.

## Design
Twelve primate samples, `assemble --seed-pool primary --guide` then `oneread`, driver and binaries of today's HEAD, defaults.
Guide = RefSeq with 31% of the transcripts hidden by hash (`GUIDE_HIDE=tx`, the guides of row 1763 and the pages); truth = the
hidden chains. Development: OR6737 (run first, alone, to check the mechanics), A119b. **Judging: the other ten** (HG002,
KB3781, chimpanzee, bonobo, Bornean and Sumatran orangutan, testis and cell line each). New and lost chains against today's
`PREFIX.guided.gtf` of the same run are counted separately.

## Bars (per sample; a rule PASSES when it passes on at least 9 of the 10 judging samples)
- **A1** hidden chains lost against today <= 2% of today's hidden found.
- **A2** the lost chains are mostly wrong: hidden share among them <= half of today's precision outside the guide.
- **B1** the added chains do not dilute: hidden share among them >= arm A's precision outside the guide on that sample.
- **B2** the addition is worth a rule: hidden chains added >= 2% of arm A's hidden found.
- **Cost (both):** wall of `oneread` <= 2 x wall of `assemble` on the sample (reported; a miss makes B "PASS, too slow").
Reported without a bar: hidden found and precision outside the guide of today / A / B beside StringTie -G, FLAIR guided, IsoQuant
guided and bambu guided on the same guide (their GTFs of row 1763); the floor-1 pool by class; the two development samples.

## My prediction (stated before the runs)
A passes everywhere (it only removes chains the chain score already rejects). B1 passes on most samples (0.3 to 0.5 against
0.15 to 0.4); B2 is the doubtful one: I expect +1% to +4% hidden chains, so B may fail on yield.

## What a PASS means
A candidate rule for the opt-in guided mode, the user's call; nothing changes in the de novo products. A Rust port
(`src/guide.rs` reading `final.gtf`; the one-read class behind a switch) follows only on the user's word. A FAIL closes the arm.
Limits stated now: RefSeq is the only truth, and on the non-human apes it is computed models; a new combination of known
introns that RefSeq lacks counts as wrong although some are real; one-read chains have no second read by definition.

## Outcome (2026-11-24; workspace `r1934_guided_final_oneread/runs/STEM/score.txt`)
OR6737 and A119b are the development samples, the other ten judge.

| sample | today: hidden / precision | arm A | A lost chains (hidden) | arm B | B added chains (hidden, share) |
|---|---|---|---|---|---|
| A119b | 7794 / 0.1377 | 7722 / 0.1588 | 7986 (72) | 8793 / 0.1467 | 11300 (1071, 0.0948) |
| OR6737 | 5428 / 0.3709 | 5421 / 0.3923 | 816 (7) | 6330 / 0.3560 | 3963 (909, 0.2294) |
| HG002_kinnex.t2t | 7238 / 0.1842 | 7207 / 0.2011 | 3468 (31) | 8418 / 0.1891 | 8679 (1211, 0.1395) |
| GGO_KB3781_fib_mm | 5834 / 0.3437 | 5822 / 0.3575 | 690 (12) | 6570 / 0.3419 | 2930 (748, 0.2553) |
| PTR_mm | 3895 / 0.4369 | 3894 / 0.4604 | 458 (1) | 5450 / 0.4253 | 4356 (1556, 0.3572) |
| PTR_AG18354_lcl_mm | 6372 / 0.3553 | 6362 / 0.3671 | 603 (10) | 7816 / 0.3651 | 4081 (1454, 0.3563) |
| PPA_mm | 3447 / 0.3157 | 3441 / 0.3263 | 371 (6) | 4221 / 0.3147 | 2865 (780, 0.2723) |
| PPA_PR00251_fib_mm | 4305 / 0.2034 | 4291 / 0.2286 | 2397 (14) | 4935 / 0.2190 | 3769 (644, 0.1709) |
| PPY_mm | 3962 / 0.1770 | 3947 / 0.1966 | 2310 (15) | 4687 / 0.1920 | 4334 (740, 0.1707) |
| PPY_AG05252_fib_mm | 5143 / 0.3095 | 5132 / 0.3364 | 1363 (11) | 5902 / 0.3109 | 3729 (770, 0.2065) |
| PAB_mm | 3062 / 0.1938 | 3047 / 0.2202 | 1962 (15) | 3736 / 0.2153 | 3514 (689, 0.1961) |
| PAB_AG06213_fib_mm | 5465 / 0.3433 | 5451 / 0.3600 | 776 (14) | 6226 / 0.3389 | 3227 (775, 0.2402) |

- **Arm A (rule G1d on `final.gtf`): PASS on 10 of 10 judging samples** (and on both development samples). It drops 371 to
  3,468 chains per judging sample, of which 1 to 31 are hidden transcripts (share 0.2% to 1.8%); precision outside the guide
  rises by 0.010 to 0.027 on every sample; no chain is added. Bars A1 and A2 hold everywhere.
- **Arm B (A + one-read chains made only of guide introns): FAIL on 10 of 10**, always on B1: the added chains are less often
  hidden transcripts than arm A's chains outside the guide (closest: PTR_AG18354 0.356 against 0.367). B2 holds everywhere:
  +12% to +40% hidden transcripts found (PTR testis 3,894 -> 5,450). Against TODAY's guided mode arm B still finds more hidden
  transcripts on all twelve samples, at a precision above today's on 6 of the 10 judging samples and below on 4.
- Cost: `oneread` within 2 x `assemble` on every sample.
- Predictions: A passes everywhere, held. B1 "passes on most" WRONG (0 of 10). B2 "doubtful, +1% to +4%" WRONG the other way
  (+12% to +40%).
- Not done, the user's call: port arm A to `src/guide.rs`; arm B as an opt-in "sensitive" guided output (it is a trade, not a
  gain by the registered bar); a stricter class inside B (e.g. two reads) would be a new row.

# PREREG r1934: guided mode on final.gtf, and one-read chains made only of guide introns (2026-11-23)

User (2026-11-23): "make the guided mode a bit better ... we are not really taking full advantage of guided mode"; approved
registering the two arms below. Written before any run of this row. Workspace `r1934_guided_final_oneread`, scorer
`bench/famcmp/guided_arms.py`, job `run.sbatch` there.

## Current state
Registered, jobs not yet run. Outcome section is appended below when they end.

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

# PREREG: rules that let one-read chains into the pipeline without much loss of precision (row 1990, 2026-11-27)

User goal: one or more empirical rules to fold the one-read minimum into the pipeline without losing much precision.

## Where things stand (not to be rediscovered)
- Default since this morning (row 1987): one-read LOCI whose every intron is in another read. Real genes (introns all
  annotated more often than two-read chains) but mostly fragments: exact chain annotated 8-17%, so `PREFIX.gtf` loses
  0.02-0.03 of exact-chain precision.
- The strict one-read tier (`PREFIX.oneread.gtf`, row 1766 S1: lone read, no longer read contains it, every junction in
  another read, ORF >= 900 nt, no early stop) is as precise as the sample's two-read chains (0.19 against 0.16), not as the assembly.
- Rows 1805 / 1806 / 1770 / 1771: junction and end proofs, a 23-feature score, BAM tags and read-end features do not
  bring one-read chains to the assembly's precision.
- Caveat that sets the bar: on real libraries the annotation holds few low-count chains at all. Exact chain annotated,
  WTC11 r1: chains with 2 reads 0.30, 3 reads 0.39, 4-5 0.47, > 10 0.84. A one-read class at 0.30 is as good as what the
  floor already lets in.

## What is new here: count the reads a chain EXPLAINS, not the reads that equal it
The floor asks for two reads with the identical chain. IsoQuant and StringTie count the reads COMPATIBLE with a model.
For a one-read chain C (`bench/o2_truth/one_read_rules.py`), over the primary reads sharing an intron with it:
- `sub` = reads whose chain is a contiguous proper part of C (cut-short copies of the same molecule kind);
- `sub_only` = those of them that no assembled transcript explains (only C does);
- `longer` = reads whose chain contains C (then C is itself the fragment); `conflict` = reads sharing an intron, neither;
- whether C is a part of / extends / conflicts with an assembled transcript, or lies where the assembly has none;
- longest ORF and early stop (the strict tier's test), introns.
Units: every spliced chain of the one-read-floor assembly that the main assembly lacks and exactly one primary read has.
Development: WTC11 r1 and mouse heart. Judged once, rule fixed: WTC11 r2, mouse adrenal, cortex, ES; gorillas as a check.

## Claim and bars
- **R1:** a rule written on the development libraries selects, on at least 3 of the 4 judging libraries, >= 300 chains
  whose exact chain is annotated at least as often as that library's TWO-READ chains (the floor's own bar, user 2026-11-18).
- **R2 ("without losing much"):** adding the selected chains lowers the exact-chain precision of `PREFIX.gtf` (without the
  row 1987 loci) by at most 0.005 on at least 3 of the 4, while adding correct chains.
- **R3:** the same rule applied to the row 1987 supported loci (keep only the ones it selects) brings their cost under 0.005 too.
Prediction: R1 holds for a rule built on `sub_only` >= 2 with no `longer` and no `conflict`; R2 holds only because the
class is small; R3 holds.

## Development result (WTC11 r1, mouse heart; job 4515685) and the rule written from it, before the judging libraries are opened
One-read chains the assembly lacks: 59,303 and 23,350; exact chain annotated 0.112 and 0.123 (two-read chains of the
assembly 0.300 and 0.289; whole assembly 0.562 and 0.578).
- Counting explained reads lifts less than predicted: `sub_only` >= 2 gives 0.18-0.25 (1,378 / 427 chains); with the frame
  test 0.27 / 0.24. Under the floor's bar on WTC11. My prediction for this feature was too optimistic.
- The strong signal is CONFLICT: reads that share an intron with the chain but splice differently elsewhere.
  More than 5 such reads: 0.12 / 0.08 (38,822 / 12,900 chains, two thirds of all). None: 0.10 / 0.24, and with a complete
  reading frame 0.33 / 0.49 (781 / 890 chains); with >= 3 introns as well 0.42 / 0.53 (375 / 644).
- A longer read containing the chain, or an assembled transcript containing it: fragments (exact 0.03-0.11, part 0.48-0.78).
**Rule B (judged for R1 and R2):** a one-read chain enters when (1) no read's chain and no assembled transcript's chain
contains it, (2) no read conflicts with it (every read sharing an intron has the same chain or a contiguous part of it),
(3) its longest ORF is >= 900 nt and does not stop more than 50 nt before the last junction.
**Rule A (reported beside it):** B and >= 3 introns.
On development: B adds 259 / 439 correct chains for -0.004 / -0.003 of exact-chain precision; A 157 / 341 for -0.001 / -0.001.
R3 is judged as: the chains at empty places with every intron in another read (the row 1987 class) that also pass B.

## Judged once (rule fixed; `judge.py` in the workspace)
Exact intron chains. "cost" = change of the exact-chain precision of `PREFIX.gtf` when the class is added.
| library | assembly precision | two-read chains | rule B: chains, correct, exact, cost | rule A: chains, correct, exact, cost | row 1987 class: chains, exact, cost |
|---|---|---|---|---|---|
| WTC11 r1 (development) | 0.5624 | 0.300 | 781, 259, 0.332, -0.0040 | 375, 157, 0.419, -0.0012 | 1,929, 0.121, -0.018 |
| mouse heart (development) | 0.5781 | 0.289 | 890, 439, 0.493, -0.0032 | 644, 341, 0.530, -0.0013 | 1,211, 0.173, -0.021 |
| WTC11 r2 | 0.5932 | 0.345 | 918, 380, 0.414, -0.0047 | 558, 249, 0.446, -0.0024 | 2,258, 0.124, -0.029 |
| mouse adrenal | 0.5749 | 0.273 | 1,005, 446, 0.444, -0.0051 | 721, 344, 0.477, -0.0027 | 1,707, 0.140, -0.028 |
| mouse cortex | 0.4963 | 0.206 | 788, 308, 0.391, -0.0026 | 467, 218, 0.467, -0.0004 | 1,464, 0.161, -0.015 |
| mouse ES | 0.3858 | 0.121 | 859, 211, 0.246, -0.0026 | 397, 142, 0.358, -0.0002 | 1,315, 0.100, -0.008 |
| gorilla OR6737 (check) | 0.3539 | 0.165 | 1,148, 181, 0.158, -0.0032 | 383, 92, 0.240, -0.0006 | 4,319, 0.078, -0.016 |
| gorilla KB3781 (check) | 0.3902 | 0.157 | 793, 146, 0.184, -0.0027 | 330, 98, 0.297, -0.0005 | 1,679, 0.094, -0.008 |
- **R1 HOLDS, 4 of 4** (B is above the two-read chains on every judging library; on the gorillas 0.158 against 0.165 and 0.184 against 0.157).
- **R2 HOLDS, 3 of 4** (adrenal misses the 0.005 bound by 0.0001). Rule A stays within 0.003 on all eight and is above the
  two-read chains on all eight.
- **R3 is EMPTY:** 1 to 8 chains of the row 1987 class pass B. The two are different populations: B = a complete
  transcript nothing else at the place contradicts; row 1987 = a real gene seen through reads that disagree (fragments).
  B cannot clean the row 1987 class; that class costs 0.008-0.029 whatever is done.
Counting explained reads (the idea this row started from) did not carry the rule; conflict and the reading frame did.

## Decision (user, 2026-11-27) and the build
Rule A into `PREFIX.gtf` by default; the row 1987 supported loci leave `PREFIX.gtf` for the side file `PREFIX.oneloci.gtf`.
`stage_oneloci` (after the main assembly of `all`): `tools/one_read_loci.py clean` (tag `one_read_chain "clean"`, ids `R1_`;
`RUSTLE_ONE_READ_CHAINS=0` = off), then the supported loci on top into the side file (`RUSTLE_ONE_READ_LOCI=0` = off).
Check (job 4515730): the built pipeline adds exactly the judged chains.
| data | chains added | `PREFIX.gtf` correct | precision | side file loci |
|---|---|---|---|---|
| mouse heart | 644 | 13,086 -> 13,427 | 0.5781 -> 0.5767 | 1,046 |
| WTC11 r2 | 558 | 20,095 -> 20,344 | 0.5932 -> 0.5908 | 2,029 |
| fly simulated (whole `all`) | 401 | 6,783 -> 7,175 | 0.9393 -> 0.9414 | 15 |
No warning, `PREFIX.final.gtf` made, the one-read tier reused the floor-1 run, a second `all` redid nothing; with both
switches off the assembly is byte-identical to the one before today's defaults; cargo test 1,591 passed.
None of StringTie, FLAIR, IsoQuant, bambu uses a reading-frame test to admit transcripts (their source; FLAIR only labels
productivity afterwards): the frame condition is ours and keeps one-read NON-CODING chains out. Without it the no-conflict
class is 0.10 / 0.24 exact on the development libraries.
Not measured: the added chains' effect on families and later stages; `PREFIX.final.gtf` does not get them.

## Addendum before running (2026-11-27, user: is anything added today detrimental for gene families?)
What families read: `PREFIX.gtf`. Added to it this week: own-evidence seeding (row 1984, already tested on planted
families) and today's clean one-read chains (not tested there). The side file is read by no stage.
Known risk before looking: rule A has no test of the lone read's placement. At a silent copy a read of the expressed copy can
sit as an arbitrary primary with nothing to contradict it; that would be a clean chain on a copy that is not expressed.
**Planted families, fresh seeds 1992 / 1993** (scenarios E all expressed, S half silent, thin = 3 reads per copy), three arms:
`before` = `--seed-pool good`, chains off; `own` = default pool, chains off (yesterday); `now` = today's default.
- **F1:** `now` against `own`: isoform precision of the mode and of the main assembly not lower by more than 0.01, exact
  isoforms not fewer, silent copies with a transcript not more, the same loci labelled expressed, three-read recall not
  lower by more than 0.02.
**Real, gorilla KB3781:** families stage with and without the chains.
- **F2:** families and copies change by under 1%; the clean chains that lie in a family copy are exact at least as often
  as the others minus 0.05.
Prediction: F1 fails on silent copies (the risk above) and holds elsewhere; F2 holds. If silent copies gain transcripts, the
repair to test is the rival test of row 1987 (no other alignment within 10 ms) added to rule A.

## Result of the family addendum (jobs 4515976, 4516082)
Planted families, seeds 1992 / 1993; `own` (yesterday) -> `now` (today's default); `before` in brackets.
| measure | 1992 | 1993 |
|---|---|---|
| clean chains added, E / S / thin | 0 / 0 / 2 | 0 / 1 / 4 |
| E: isoform precision, exact isoforms | 0.742 -> 0.742, 196 -> 196 (0.622, 196) | 0.767 -> 0.767, 198 -> 198 (0.647, 198) |
| S: isoform precision, mode / main | 0.528 -> 0.528 / 0.595 -> 0.595 (0.472 / 0.477) | 0.490 -> 0.488 / 0.560 -> 0.557 (0.449 / 0.451) |
| S: silent copies with a transcript, mode / main | 44 -> 44 / 31 -> 31 (44 / 44) | 47 -> 47 / 34 -> 34 (47 / 46) |
| loci labelled expressed, E / S | 83 -> 83 / 44 -> 44 | 83 -> 83 / 45 -> 45 |
| thin: isoform recall, mode / main | 0.645 -> 0.649 / 0.557 -> 0.566 (0.702 / 0.623) | 0.684 -> 0.689 / 0.553 -> 0.566 (0.719 / 0.627) |
- **F1 HOLDS on every part.** The silent-copy risk I named did not occur (prediction wrong); one wrong transcript in one S set.
Gorilla KB3781, families stage without -> with the 330 chains: families with >= 2 copies 260 -> 260, copies 724 -> 724,
the 260 families identical as sets of loci. Three clean chains overlap a family copy (none exact); 327 elsewhere (0.300).
- **F2 first part HOLDS** (no change at all); the second part cannot be judged on three chains.
(The first gorilla submission failed in the families stage because my job exported the human duplication file; rerun without it.)
The change that does move family numbers is yesterday's seeding (row 1984), as already recorded: precision up, silent
copies in the main assembly down, three-read recall of the mode down 0.035-0.057.

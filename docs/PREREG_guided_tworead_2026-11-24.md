# PREREG: guided mode, chains of guide introns at a two-read floor (row 1939, 2026-11-24)

## Current state
Registered, not scored.

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

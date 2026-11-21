# Row 1866: one-read candidates with the pooling guard bypassed, promoted by a partner library (pre-registered 2026-11-21)

User (2026-11-21): do the next step of row 1864. There: at a 1-read floor a single-read chain passes the assembly gate only
on pooled locus support, and the anti-bridging guard marks full-length single reads of long genes as chimeric, so the
`oneread` stage (whose candidates come from a floor-1 run) never sees them: on WTC11 2,401 correct single-read chains
beside the 4,270 it does see. This row feeds the stage with the guard bypassed and asks whether the existing promotion rule
(row 1808: a ranked one-read chain in the best 20% that a read of a second library also carries; as precise as the
default's two-read chains) then promotes more correct chains at the same quality. Workspace `r1866_oneread_guard`.
Committed before the switch was written and before any run. The guard in the main path is not touched.

## The switch, fixed
`RUSTLE_ONEREAD_GATE=1` (opt-in; unset = today's command, byte-identical): the floor-1 run inside stage `oneread` gets
`RUSTLE_GATE_MIN_READS=1`, so every single-read spliced chain of pass 1 reaches the candidate set. Ranking model
(`tools/data/oneread_r1806.json`), the best-20% cut and the partner test are unchanged. `PREFIX.gtf` is untouched.

## Runs
Per library: `assemble` at today's defaults, then `oneread --partners ...`, once without and once with the switch (arms A, B).
| | library | partner(s) | role |
|---|---|---|---|
| D | WTC11 replicate 2 | WTC11 replicate 1 | development (rows 1862 / 1864 were read on it) |
| C1 | WTC11 replicate 3 | WTC11 replicate 1 | check; replicate 3 has never been assembled |
| C2 | mouse adrenal (ENCSR852HRM) | heart, cortex, ES (row 1808's partners) | check; was row 1808's own check library |
| C3 | chimpanzee testis | chimpanzee LCL (`PTR_AG18354_lcl_mm.bam`) | check, primate |
Replicate 1 is used as a partner only (its reads are looked up, it is not assembled or scored).

## Measured, per library and arm
`PREFIX.oneread.promoted.gtf`: chains, correct (exact chain in the truth), precision; the ranked candidate set's size; and the
comparator of row 1808: precision of the default's two-read chains (`final.gtf`, `reads` = 2).

## Reading, fixed now (on C1, C2, C3; D is shown)
- **GAIN:** arm B promotes at least 10% more correct chains than arm A.
- **QUALITY:** arm B's promoted chains are at least as precise as the library's two-read chains, and not more than 0.02
  below arm A's promoted chains.
- **Candidate** = GAIN and QUALITY on all three check libraries. A candidate makes the switch the proposed default of the
  `oneread` stage; the user decides. GAIN without QUALITY anywhere = the guard was doing useful work for this stage.

## Expectation
GAIN on the WTC11 libraries (the hidden candidates are 56% of the visible correct ones there). I am unsure about QUALITY: the
best-20% cut is relative to the candidate set, so more candidates move the cut; the hidden chains are long (20 or more
introns are the most precise on WTC11) and the partner test should be strict for long chains. On mouse adrenal I expect a
smaller gain (mouse ES: hidden chains 76% of visible, at 0.044).

## Limits
- WTC11 replicates share the cell line: a partner read confirms a chain more easily than a different tissue would.
- The ranking model was fitted (row 1806) on candidate sets that lacked these chains.
- The floor-1 run uses collapsed read counts (row 1850's note).

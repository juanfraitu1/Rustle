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

## Result (2026-11-21): GAIN on every library, QUALITY on none: not a candidate by the bars

Jobs 4491366-4491369, all steps exit 0; code `d1d6d873`; `r1866_oneread_guard/NAME/{A,B}`.
| library | `final.gtf` | two-read chains | A (today): candidates; promoted correct / all | B (switch): candidates; promoted correct / all | GAIN | QUALITY |
|---|---|---|---|---|---|---|
| WTC11 r2 (development) | 20,016 / 0.639 | 0.388 | 32,744; 574 / 1,731 = 0.332 | 52,756; 865 / 2,906 = 0.298 | +50.7% | no |
| WTC11 r3 (check) | 25,423 / 0.584 | 0.309 | 34,492; 423 / 1,393 = 0.304 | 63,941; 771 / 2,771 = 0.278 | +82.3% | no (0.278 < 0.309; 0.026 under A) |
| mouse adrenal (check) | 14,162 / 0.615 | 0.305 | 14,910; 294 / 1,351 = 0.218 | 26,989; 484 / 2,109 = 0.229 | +64.6% | no (0.229 < 0.305) |
| chimpanzee testis (check) | 20,309 / 0.432 | 0.251 | 53,298; 576 / 2,325 = 0.248 | 132,547; 1,145 / 4,900 = 0.234 | +98.8% | no (0.234 < 0.251) |
Correct chains the promotion adds to `final.gtf`, A -> B: WTC11 r3 +1.66% -> +3.03%; mouse adrenal +2.08% -> +3.42%;
chimpanzee +2.84% -> +5.64%. The chains only B promotes are 0.25-0.29 correct (493 / 1,813; 306 / 1,221; 825 / 3,412).

- **GAIN: yes on all three check libraries** (65-99% more correct promoted chains). The guard was hiding about as many good
  candidates as the stage could see.
- **QUALITY: no on all three.** B is within 0.03 of A everywhere (-0.026, +0.011, -0.014), but below the library's two-read
  chains on all three, and on WTC11 r3 more than 0.02 under A.
- **Row 1808's own claim no longer holds under today's defaults:** arm A, the promotion as it ships, is also below the
  two-read chains on all three (0.304 / 0.309, 0.218 / 0.305, 0.248 / 0.251). The comparator moved: with every read counted
  (row 1850) and the new rules (rows 1845, 1848) the default's two-read chains are a cleaner set than in row 1808
  (mouse adrenal then 0.258, now 0.305).
- Expectation: GAIN as expected; the doubt about QUALITY was justified, for a reason I had not named (the comparator).
- Cost: the switch adds 5-39% to the stage's wall time.

**Decision:** not a candidate as registered. The switch stays opt-in (`RUSTLE_ONEREAD_GATE=1`, default byte-identical): it
roughly doubles what the partner promotion finds at 0.01-0.03 lower precision. For the user: the promoted tier as a whole
is now a recall option at about 0.22-0.30 precision, no longer "as precise as two-read chains"; its documentation should
say so. WTC11 replicate 3 is now assembled and has been used for this check.

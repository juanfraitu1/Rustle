# PREREG: can a minimap2 parameter change remove or decide the near-tied reads (row 1940, 2026-11-24)

## Current state
Registered, not run.

User (2026-11-24): the advisor is inquisitive; "a quick sweep with minimap2 to prove that no parameter sweep can avoid the
near tie multimappers". A sweep cannot prove this for every setting; it can show it for a stated grid, and one part of
the claim does not depend on the aligner at all (see "What no setting can change"). Written before any cell was run.

## Setup
- HG002 Kinnex on CHM13, the 42,770 selected reads and the repaired truth table of row 1785 (reused).
- Baseline B0 (the command behind every earlier row): `minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes`,
  version 2.30. Every cell changes ONE thing in B0; `-N 50 -p 0.1` stay, so rivals are always reported.
- Cells (37): seeding `-k 11 13 17 19 21`, `-w 3 10`, `-f 0.00002 0.002 0.02`; chaining `-n 2 6`, `-m 20 80`,
  `-g 2000 10000`, `-r 200,20000 1000,50000`, `-G 100k 500k`; scoring `-B 2 6 8`, `-O 2,24 12,32`, `-C 0 9`;
  alignment `-z 200,100 800,400`, `-s 40 200`, `--end-bonus 10`, `-J 0`; presets `-ax splice` and `-ub` in place of
  `-uf`; annotated junctions `--junc-bed` (CHM13 annotation). A cell minimap2 refuses is reported as refused.
- Per cell: (a) the number of reads whose best rival at another locus is within 10 alignment-score points of the
  primary, and with exactly the same score (the second count does not depend on the score scale; cells that change the
  scoring are marked); (b) the primary's pick on the decidable reads: near ties (600, baseline 553 / 47) and clear-lead
  reads (14,596, baseline 18 wrong), pick = the candidate locus of the table the primary overlaps on its strand.

## Decision
- Split by read name (sha256 even = tuning, odd = judging), as row 1936.
- DECIDES THE NEAR TIES: the cell with the fewest wrong near ties on the tuning half (ties: fewer clear-lead wrong, then
  nearest B0 in the list) has, on the judging half, wrong <= two thirds of B0's (rounded down) AND clear-lead wrong
  <= B0 + 3. Otherwise NO CELL DECIDES THEM.
- AVOIDS THE NEAR TIES: a cell with unchanged scoring has at most half of B0's reads with an exactly tied rival, while
  its clear-lead wrong stays <= B0 + 3 on all decidable reads. Otherwise NO CELL AVOIDS THEM.

## What no setting can change (stated, not swept)
Row 1784 / 1791: 4,249 of the 7,770 near ties match several places of HG002's OWN two haplotype assemblies equally
(2,950 of them several copies of one haplotype). For those reads the tie is in the individual's genome, not in the
alignment to CHM13; an aligner setting can only hide it (fewer secondaries reported, another primary picked).

## Predictions (mine)
NO CELL DECIDES and NO CELL AVOIDS. Seeding and chaining cells leave both counts within about 10% of B0; scoring cells
move the within-10 count mechanically (the scale changes) and the exact-tie count little; a few cells (large k, `-f`
small, `-s 200`) lose alignments and get worse on the clear-lead reads. Annotated junctions change the primary of a few
reads in either direction.

## Limits
One sample, one individual, 42,770 reads on family genes; one parameter at a time, no combinations; the truth is the
catalogue-free test with its known limits.

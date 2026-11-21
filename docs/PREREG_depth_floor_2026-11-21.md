# Row 1852: a read floor that scales with library depth (pre-registered 2026-11-21)

User (2026-11-21): do it, after row 1851 (a chain of N reads is worth far less in a deep library: 2 reads are 13% annotated
on mouse ES against 31% on heart; ES cut to heart's reported count has heart's precision). Measurement only, as a post-cut on
today's `final.gtf`; no assembler code. A build follows only on the user's word. Workspace `r1852_depth_floor` (`floor.py`).
Committed before the depths were counted and before any cut was made.

## The rule measured
Depth M = millions of mapped alignment records in the BAM (`samtools idxstats`, all contigs). Floor F = max(2, ceil(c x M));
a multi-exon transcript of `final.gtf` is dropped when its `reads` (every record counted, the default since row 1850) is
below F. Arms c = 0.25, 0.5, 1 reads per million, chosen a priori (1 per million is the round figure; row 1851's cuts on ES,
reads >= 5 to 10 at 8.7 M reads, fall between 0.5 and 1). Single-exon transcripts are not cut and not scored.

## Data
`r1849_dup_count_rerun/NAME/dc/x.final.gtf` (today's default) on the twelve datasets, truths of row 1849. Reused, said.

## Reading, fixed now
- An arm **pays** on a dataset when the chains it removes are at most half as precise as the dataset's `final.gtf` AND it
  keeps at least 90% of the correct chains. On a dataset where F = 2 the arm changes nothing and counts as neutral (left out
  of the count, reported).
- **Candidate** = an arm pays on at least 6 of 7 non-primate and 4 of 5 primate datasets, neutral ones aside, with at least
  three datasets changed. A candidate is put to the user as a rule to build (in the polish, before the filter).
- Reported per arm and dataset: F, chains removed, their precision, correct kept, precision after; and the spread of
  precision across the four mouse libraries before and after (does the floor make libraries comparable).

## Expectation
The simulated sets decide it. Mouse simulated has the depth of mouse ES by construction and a complete truth in which 2-read
transcripts are real, so a floor of 5-9 reads there removes true chains at the assembly's own precision: I expect it to FAIL
on mouse simulated (and fly, if deep) and to pay on the deep real libraries (ES, A119b, the gorillas). If so, the floor is a
statement about real-library noise that simulated reads do not carry, and the verdict by the bar is FAIL or undecided with
the table.

## Limits
- A post-cut on `final.gtf` is not the same as a floor inside the assembly (a removed chain's reads are not re-assigned).
- `idxstats` counts supplementary and secondary records too; depth is approximate.
- Real-library truths are incomplete: row 1851 showed most low-count "wrong" chains are canonical and seen by other tools.

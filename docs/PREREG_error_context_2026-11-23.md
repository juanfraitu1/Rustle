# PREREG: sequencing-error context against minimap2's choice, real reads with own-haplotype truth (row 1923, 2026-11-23)

User, after row 1922: "we also have the haplotypes for the gorilla fibroblasts". Row 1922 left one thing untested:
minimap2 charges a sequencing error like a true difference between copies, and real HiFi errors are mostly short indels
in homopolymers. KB3781 (gorilla fibroblast) is the assembled individual: at the read's true copy every mismatch is a
sequencing error or the other haplotype, so it is the right data to ask whether discounting such errors changes the
placement, and whether the change is right.

## Data (existing tables, nothing new aligned)
- KB3781 on its own reference: `r1785_truth_repair/kb/reads.tsv` (repaired truth, row 1785: 23,710 decidable reads;
  primary wrong on 1 of 20,611 with a clear lead and on 5 of 45 near ties), alignments `r1775_eichler_baseline/mm.sel.bam`.
- HG002 on CHM13 (sample differs from the reference; reported beside it, same rules): `r1785_truth_repair/hg/reads.tsv`,
  alignments `r1776_hg002_room/mm.sel.bam`.
Known before any run: the room on KB3781 is 6 reads.

## Measurement
1. Error profile: on decidable reads whose primary is on the origin (KB3781), events per kb of aligned read:
   substitutions, indels inside a homopolymer run (>= 3 of one base in the reference, the indel made of that base),
   other indels. This is what minimap2 pays for errors alone.
2. Context score per alignment, two fixed forms: A = `ms` + (6 + length) for every homopolymer indel (they cost
   nothing); B = `ms` + 7 for every one-base indel anywhere.
3. Rule: the locus with the highest context score among the read's reported alignments; when it is one locus and not
   the primary's, the read moves. Scored against the origin: fixed / broke, by stratum (clear lead, near tie).

## Claims
C1 (either form, KB3781): fixes at least twice what it breaks AND fixed - broke >= 3 (half of the room).
C2 (same form, HG002): fixes at least twice what it breaks.
C1 + C2 = worth a registered design inside the `o2` labels. Otherwise the idea is closed with the counts.

## Predictions (mine)
Homopolymer indels are the majority of error events. The rule moves few reads and fails C1: where the primary is wrong
the cause is not an error in the read (rows 1778, 1783: the copies' reference differences do not separate them).

## Result (job 4496786, workspace `r1923_error_context`, `ctx.py`)
Error profile of real HiFi reads at their true copy (events per kb of aligned read):

| | substitutions | homopolymer indels (one-base) | other indels | all |
|---|---|---|---|---|
| KB3781 (23,704 reads, 56.6 Mb) | 0.452 | 0.467 (0.432) | 0.200 | 1.12 |
| HG002 (18,232 reads, 38.3 Mb) | 0.571 | 0.607 (0.571) | 0.184 | 1.36 |

The rule (A = homopolymer indels free, B = one-base indels free), decidable reads:

| | reads | primary right | A: moved / fixed / broke | B: moved / fixed / broke |
|---|---|---|---|---|
| KB3781 clear lead | 20,611 | 20,610 | 0 / 0 / 0 | 0 / 0 / 0 |
| KB3781 near ties | 45 | 40 | 0 / 0 / 0 | 0 / 0 / 0 |
| HG002 clear lead | 14,596 | 14,578 | 13 / 1 / 12 | 13 / 1 / 12 |
| HG002 near ties | 600 | 553 | 31 / 3 / 27 | 28 / 1 / 26 |

- **C1 FAILS: on the gorilla fibroblast the rule moves no read at all**, in 20,656 decidable reads with a rival locus.
  None of the 6 wrong primaries is caused by an error the context could discount.
- **C2 FAILS and the rule is harmful where the sample differs from the reference:** on HG002 it moves 44 reads, fixes
  4 and breaks 39. Short indels in homopolymers are real differences between copies as often as they are errors, so
  making them free removes evidence.
- Predictions: "homopolymer indels are the majority of error events" was wrong (42-45%; substitutions are as common);
  "few moves, C1 fails" held.
- **Closed.** With row 1922: minimap2's choice is not wrong because it mis-weighs errors. A real HiFi read carries about
  one error per kb, a tenth of what separates most copies that can be separated at all; where copies differ by less,
  the read usually does not cover the difference. Do not offer error-aware re-scoring of the alignments again.

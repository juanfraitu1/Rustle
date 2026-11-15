# PREREG row 1710 (2026-11-14): allele vs copy on short reads, step 0 on real data (measurement, then close or hand over)

The open line of `LEDGER_O3.md` since 2026-11-01. Its source (archive `LEDGER_O3.md`, simulated world, MRGPRX): the
allele cutoff delta = 0.00958 is calibrated on whole copies, so two heterozygous sites within about 200 bp of a short
read exceed it and a candidate made of such reads is called a copy when it is an allele. It was seen in simulation
with reads of 190-506 bp. Never measured: whether real flagged candidates that are alleles are short.

## Measurement (tables that exist; nothing is rerun)
Flagged candidate contigs of `candidate_copies` default pool on KB3781 with their haplotype class from `pool_score.py`:
`o2pool/scoreA.contigs.tsv` (row 1513 panel) and `r1703_secpool/scoreA.contigs.tsv` (row 1703 panel), 80 families.
Classes `allele_of_known_copy` (the false flag this confound would produce) against `absent_from_ref`. Length = the
table's `len` column (contig length).
The confound needs a short sequence: k heterozygous sites exceed delta below k / 0.00958 bp (209 bp for two sites, 313
for three; with the error floor a little more). SHORT = len < 500 bp, a generous bound.

## Reading rule (written before the lengths of the allele class were looked at)
- Fewer than 8 allele-class contigs over both panels: UNDERPOWERED, reported, the line stays open.
- Fewer than 20% of the allele-class contigs SHORT: the short-read confound is not what makes real allele false flags
  at Iso-Seq read lengths. The open line is CLOSED as a simulation-only effect; no rule is built.
- 20% or more SHORT: the confound is real on real reads; a length-scaled cutoff is the registered next row, and
  whether to run it is the user's decision.
Reported without a bar: length quantiles per class, the `d` of the allele-class candidates where the candidates table
has it.
Seen before writing: the 12 + 7 `absent_from_ref` rows of row 1703 (all 1.4 kb or longer). Not seen: any allele-class row.

## Outcome (2026-11-14): UNDERPOWERED as registered (7 allele-class contigs, floor 8)
| class | contigs | under 500 bp | min / median / max length |
|---|---|---|---|
| allele_of_known_copy | 7 (1 + 6) | 0 | 1,542 / 1,570 / 9,921 |
| absent_from_ref | 15 | 0 | 1,145 / 2,093 / 6,272 |
| not_in_genome | 41 | 0 | 645 / 2,261 / 6,698 |
| between | 8 | 0 | 1,203 / 1,498 / 6,445 |

By the rule the line is not closed by this row: one contig short of the floor. What the tables do show, outside the
rule: none of the 71 flagged contigs of any class is under 500 bp (shortest 645), so the precondition of the confound
(a sequence short enough for two or three sites to exceed delta) does not occur among flagged candidates on these 80
families. The 7 allele-class candidates sit at d 0.017 to 0.147 from their nearest reference locus over 1.5 kb or
more: they are not two-site effects. Recommendation: close the line as a simulation-only effect; that is the user's
decision, since the registered floor was missed.

User ruling 2026-11-14: "Close as simulation-only". The item leaves the O3 open list; the UNDERPOWERED verdict of step 0 stands as recorded.

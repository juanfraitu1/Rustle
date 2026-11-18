# PREREG r1773: which loci do a family's multimapping reads reach? (2026-11-17)

Descriptive, the reverse direction of r1772 (`docs/PREREG_sibling_reach_2026-11-17.md`), written after r1772's result
and before any count of this row. Workspace `r1773_sibling_reverse`, script `bench/famcmp/sibling_reach.py` with
`SIB_REVERSE=1` (r1772's tables are unchanged by the switch).

## Question
r1772: nearly every sibling locus carries secondary records of the other members' reads. For those records to FIND
the siblings (advisor), the loci they reach must mostly be siblings. What share are?

## Definitions
- Truth families, strand and exon rules, samples and models as r1772 (OR6737, KB3781, PTR). Families with >= 10
  source molecules (primaries on the family's members).
- One pass over every secondary record of those molecules, genome-wide. Each record goes to one target, first match:
  1. member: exons of a member of the molecule's family on its strand (not the member holding the molecule's primary);
  2. older paralog: a gene that is a Compara paralog of a member at the eutheria level but not in the family;
  3. pseudogene: an annotated gene whose biotype contains "pseudogene" (family not known);
  4. other coding / other non-coding annotated gene (exons, strand);
  5. unannotated: no annotated exon on the strand; records in adjacent 1 kb bins of one chromosome and strand = one locus.
  A record on the exons of the member holding its own primary is dropped.
- A target is reached at (K, band) when >= K molecules of the family have a record there with AS >= band x the
  molecule's primary AS. Bands any / 0.5 / 0.98. K rules: 2; 5; "5%" = at least 5% of the family's source molecules
  and at least 2.
- Member share = reached member targets / all reached targets, pooled over families. "Relative share" also counts
  older paralogs and pseudogenes as correct (pseudogenes cannot be attributed to a family here, so it is an upper bound).

## Reading rule (per band and K rule, both gorillas; PTR reported, not read)
- member share >= 0.80: the reads single out the siblings at that setting;
- < 0.50: they do not;
- between: partial. The setting's recall is r1772's reach in the same band.

## Predictions
Band any, K = 2: member share below 0.30 (fails). Band 0.98: above 0.80 (holds), at the recall r1772 already showed
(0.15-0.43). No setting has both member share >= 0.80 and r1772 reach >= 0.80.

## Limits stated beforehand
- Unannotated targets may be real unannotated copies; the annotation decides the class.
- "5%" and K = 5 were not used in r1772, so recall at those settings comes from this run's member counts
  (members reached / members of judged families that hold no primary of their own >= 2), reported beside the share.

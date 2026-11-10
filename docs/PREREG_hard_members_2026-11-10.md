# PREREG r1601: the family members StringTie and FLAIR both miss — are ours real copies, and do family layers add more? (2026-11-10)

User question (2026-11-10): for known multi-copy families, do we find the copies CAT says are expressed where StringTie and
FLAIR (which drop or demote multimappers) do not, and does the family-wise mode add more? Written before any re-scoring.

## Hard set H (fixed)
`~/famcmp/A119b_any.genes.tsv` (r1521 panel, Amendment 2, T definition: >= 2 molecules with a primary or a tied secondary on the
member's exons): the members with `stringtie == 0` and `flair == 0`. n = 150. Seen already (counted 2026-11-10 before this file):
the r1521 run's GTF (`fusedgw/A/x.gtf`) overlaps 68 of them; 82 are missed by all three. Nothing finer was looked at.

## Levels, per member per GTF (models on the member's strand)
- L1 overlap: a model with an exon overlapping the member's CAT exon blocks (the r1521 hit).
- L2 copy-specific: an L1 model whose exons overlap no other member of the same family in the genes table (same strand).
- L3 structure: an L2 model that, for a spliced member, carries >= 1 CAT intron of the member exactly (both ends, from the
  CAT/Liftoff slim GFF transcripts of that gene); for a single-exon member, has >= 50% of its exonic bases inside the
  member's exon blocks.
- L4 exact: an L2 model whose intron chain equals one CAT transcript's chain of the member (spliced members only).

## Arms
- `fusedgw/A/x.gtf` (r1521 base) and `fusedgw/A/x.cut.gtf` (r1523 family-wide cuts on it).
- `audit/e2e2/a1.gtf` (current default assembly) and `audit/e2e2/a1.final.gtf` (current layered output).
- StringTie and FLAIR as controls (L1 = 0 on H by construction).

## Bars (registered)
- Q1 "ours recovers the hard members as real copies": PASS if, on `a1.gtf`, L3 >= 0.7 x L1 and L3 >= 40.
- Q2 "family-wise output adds copies": per pair (x.gtf -> x.cut.gtf; a1.gtf -> a1.final.gtf), PASS if L3 rises by >= 5 and
  no member falls out of L3.
- Descriptive, no bar: the read class (T-only / MAPQ0 >= 0.5 / ordinary) of members at each level; the misses by family name.

## Outcome
(empty until the run)

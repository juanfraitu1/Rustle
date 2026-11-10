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

## Outcome (2026-11-10, `~/famcmp/r1601.hard.tsv`)
| arm | L1+ | L2+ | L3+ | L4 |
|---|---|---|---|---|
| x.gtf (r1521 base) | 68 | 39 | 35 | 26 |
| x.cut.gtf (r1523 cuts) | 68 | 66 | 35 | 26 |
| a1.gtf (current default) | 68 | 39 | 35 | 26 |
| a1.final.gtf (current layered) | 38 | 37 | 36 | 26 |
| StringTie / FLAIR | 0 | 0 | 0 | 0 |

- **Q1 FAIL as registered**: L3 35 < 0.7 x 68 and < 40. Of the 35, 32 are multimapper members (MAPQ0 18, T-only 14): TSPY 14,
  GAGE 4, CT45 3, SPDYE 3, PRR 2, CFC 2, SPANXA/SSX/LGALS/GGTLC/PGA/VCY/TCAF 1 each; 26 match a CAT intron chain exactly.
- The 33 at L1-L2 are 27 TAF11L members (26 T-only): one model spans several copies of the chr5 tandem array (tandem
  readthrough, register 1516/1517); plus 2 OR, GAGE12E, NOTCH2NLR, FAM, IFNA.
- **Q2 FAIL** on both pairs: cuts +0 L3 (they make 27 more models copy-specific but carry no CAT intron); final.gtf +1 (H2AB3),
  none lost; final drops the 30 overlap-only models.
- Missed by all three: 82, of them 63 ordinary (OR 7, TAF11L 6, DEFB 6, DEFA 4, SIGLEC 4, ...), 19 multimapper.

# What is a corrected primary worth? The oracle ceiling (2026-11-19, PSU, register 1824)

User's question: tied multimappers pick their primary alignment arbitrarily, so the primary used for seeding may be wrong;
seeding with every secondary misleads (rows 1767, 1820). For the reads that can be told apart, does correcting the primary gain
anything? Development only: all four HG002 quarters pooled (spent; reuse per the user's ruling). Workspace
`goal_multicopy_dev/corr` (linked from `r1824_primary_correction`); script `bench/famwise_eval/correct_primary.py`.

## Design
The ORACLE: HG002's own-haplotype truth (`neutral_truth.py`, repaired rule) decides which reads are corrected. For every read
the truth decides (one origin copy) whose primary record is off its origin while another record is on it, the best-ms record
on the origin becomes the primary; nothing is removed. Assemble with `--seed-pool primary` (Rustle) and StringTie `-L` on the
primaries, before and after. Any real correction rule can gain at most this.

## Result
- Decidable molecules 73,616; primary off its origin **241 (0.33%)**, 161 of them near-tied (two best places within 10 ms);
  240 movable (a record on the origin exists).
| | resolved units (of 3,322) | precision (copy claims) | family members found (of 409 with own reads) |
| Rustle primaries | 2,708 | 0.873 | 307 |
| corrected, all 240 | 2,714 (+6) | 0.874 | 307 |
| corrected, 160 near ties | 2,714 (+6) | 0.873 | 307 |
| StringTie -> corrected | 2,266 -> 2,274 (+8) | 0.852 -> 0.852 | 343 -> 343 |
- Rustle, all 240 corrected, against the uncorrected: 13 new chains (6 EQUAL, 4 NARROW, 2 wrong copy, 1 unjudged), 13 lost (5 EQUAL,
  1 NARROW, **7 wrong-copy chains that 397 misplaced molecules had built**).

## Reading
Correcting primaries, even perfectly, adds about 0.2% resolved units and no family member; its one real benefit is removing a
few wrong-copy chains built by misplaced reads. The primary is wrong too rarely (0.33% of decidable molecules) for the seeding
to depend on it. Practical rules recover at most 5-9 of the primary's errors (rows 1778-1785) and copy_phasing none (row 1823),
so a real rule would get a fraction of this ceiling. Not worth building for assembly; the reads that matter are the ones that
cannot be told apart (copy sets, `MULTICOPY_FAMILIES_2026-11-19.md`).

Caveats: HG002 is a cell line on another individual's reference (CHM13); the truth decides only reads with one near-perfect home
on HG002's haplotypes (ties, `outside` and unjudged reads are untouched); no testis sample has such a truth.

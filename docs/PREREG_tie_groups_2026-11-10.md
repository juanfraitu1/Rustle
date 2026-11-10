# PREREG r1603: tie groups — a family-level product for the copies no read can tell apart (2026-11-10)

Follows r1601/r1602: 72 of the 150 A119b truth-family members that StringTie and FLAIR both miss have no read that prefers
them over another placement. A per-copy model there is unsupported; the honest product is the group. `assign` cannot give it:
its catalogue has no copy in the TAF11L array (the assembly has one 113 kb model there). So the groups are built from the
reads' own tie sets. Written and committed before the script is run. Panel first (fast-feedback rule); A119b is a spent
sample, so a PASS here is development evidence only and the held-out run (KB3781) is a later row.

## Product (`bench/famcmp/tie_groups.py`, prototype)
- **Molecules**: every read with a non-supplementary record overlapping the span of any of the 1,136 panel members
  (`~/famcmp/A119b_any.genes.tsv`). All their records are taken from the whole BAM (one `samtools view -N` pass).
- **Tied record**: AS equal to the molecule's best AS over its records; the molecule is tied when >= 2 of its records are.
- **Locus**: per (chrom, strand), tied records sorted by start; a record joins the open locus when it overlaps the locus's
  first record by >= 50% of the shorter of the two, else it opens a new one.
- **Group**: connected components of loci, two loci linked when >= 2 molecules have a tied record in both.
- **Per locus model**: the most common intron chain of its records (unspliced counts as a chain); start and end = medians
  over the records with that chain.
- Outputs: `PREFIX.gtf` (one model per locus, attributes `group`, `locus_molecules`, `group_molecules`, `group_loci`) and
  `PREFIX.groups.tsv` (one row per group; each molecule counted once per group, never split across loci).

## Bars (registered; M = the 72 members with `signal == indistinguishable` in `~/famcmp/r1602.signal.tsv`)
A member is **reached** when a locus model on its strand has an exon overlapping its CAT exon blocks.
- B1 reach: >= 58 of 72 reached.
- B2 purity: among groups with a locus reaching any panel member, >= 0.90 have no foreign locus. A locus is foreign when
  it overlaps exons of a CAT gene (same strand, `chm13v2.0_CAT_Liftoff.genes.tsv`) but of no gene in the panel family that
  most of the group's reached members belong to. Loci overlapping no CAT gene are reported, not judged.
- B3 unit structures: >= 20 of the 31 indistinguishable TAF11L members are reached by a model of span <= 5 kb.
- B4 cost: the job (names + one BAM pass + grouping) finishes within 1 h wall on one node.
- PASS = all four. Any miss = FAIL as registered, with the numbers.
- Descriptive: group size distribution, loci with no CAT gene, and reach on the other 78 hard members.

## Not in this row
A Rust port, a pipeline stage, per-copy abundance. Decided after the outcome.

## Outcome (2026-11-10, job 4456161, `/anvme/workspace/v138da13-rustle/r1603/`): FAIL as registered
| bar | result | verdict |
|---|---|---|
| B1 reach >= 58 of 72 | 39 | FAIL |
| B2 purity >= 0.90 | 273 of 318 = 0.858 | FAIL |
| B3 TAF11L <= 5 kb models >= 20 of 31 | 11 | FAIL |
| B4 wall <= 1 h | 1,957 s (names 234 s, BAM pass + build 1,721 s) | pass |

- 1,072,265 molecules -> 76,500 tied records, 2,658 loci, 1,615 groups; 1,210 groups have one locus (their partner tied
  with < 2 molecules, or lies on the other strand of the seed rule). 247 loci in judged groups overlap no CAT gene.
- Unreached: 20 of 31 TAF11L, TSPY10 x2, and 11 singles (CENPVL2, DEFA1B, DEFB103A, H2AB1, LGALS9C, NOTCH2NLR, OR11H12,
  PPIAL4E, PRAMEF7, PRR20E, SPDYE13).
- Post hoc (not a registered result): over the TAF11L array the group models are unit-length (1.9-4.8 kb, 25 of 27) but
  carry 1-2 molecules each. The array has 1,818 secondary records; exact AS equality keeps few of them. The panel's own T
  definition is AS >= 0.98 x best; the exact-tie rule is the likely cause of the reach miss. A near-tie arm is a new row.

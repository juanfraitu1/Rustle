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

## Outcome
(empty until the run)

# PREREG — isoform-level phantoms in the shipped GTF: the same isoform emitted at several copies on abstaining reads (2026-09-09, before the measurement)

**Question (user):** the GTF's transcripts sit where the aligner placed their reads, not where O2 did. An isoform
whose reads coin-toss over near-identical copies is therefore emitted once per copy the aligner scattered
them to. How often? And how often are two copies' transcripts the SAME isoform with independent evidence at
each (genuinely shared) — the object the isoform-level certificate (step 2) will have to model?

**Objects.** Human MCL0, `bakeoff/human/ours_final2.gtf` (684 family transcripts: 27 `assigned`, 137
`undecidable`, 520 `unadjudicated`), copy spans `copies16.tsv`, copy-to-copy alignments
`scratchpad/human_gspans.paf` (minimap2 asm20 `--eqx`, genomic unit spans). An isoform = an intron chain.
Two transcripts at copies A and B are **the same isoform** when A's chain, lifted base-by-base through the A→B
alignment, matches B's chain at every intron boundary within ± 5 bp (indel jitter; exact-after-lift also
reported). Isoform groups = connected components of "same isoform" over the 26 copies.
Member reads of a transcript = primaries (-F 2308) in the region whose intron chain equals the transcript's
(the collapse rule). Evidence at copy B for a group = ≥ 1 member read at B that is either a UNIQUE mapper
(not admitted by the AS-tied gate, i.e. absent from `ours_final2.assignments.tsv`) or O2-`assigned` to B.
Denominator = family transcripts with ≥ 2 introns (single-intron and unspliced chains are excluded: they match
too promiscuously).

| # | prediction | refuted by |
|---|---|---|
| P1 | **10–40 %** of the multi-intron family isoforms are emitted at ≥ 2 copies | < 5 % or > 60 % |
| P2 | of the isoform groups emitted at ≥ 2 copies, **≥ 50 % have ≥ 1 copy with NO evidence** (all its member reads abstain) — a phantom placement | < 30 % |
| P3 | **≥ 5 groups are genuinely shared**: evidence at ≥ 2 copies (the object the isoform certificate must represent, not collapse) | < 2 |
| P4 | the phantom copies concentrate in the near-identical groups (6/7/8, 10/11, 12/19, 23/24/25, 2/8): ≥ 80 % of phantom transcripts sit in a copy with a sibling ≥ 0.985 identity in the same group | < 60 % |
| P5 | the 27 `assigned` transcripts are never phantoms (their copy has a certified read by construction) | any |
Script `bench/isoform_copy_lift.py`. Human only. No shipped behaviour changes; the number is the baseline
for step 3's emission rule.

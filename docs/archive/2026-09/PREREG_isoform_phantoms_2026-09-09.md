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

## Outcome (2026-09-09) — `bench/isoform_copy_lift.py`, groups in `bakeoff/human/hard/isoform_groups.tsv`
542 multi-intron family isoforms; 156 same-isoform transcript pairs across copies (154 exact after lift), 452 groups.
| # | verdict |
|---|---|
| P1 | ✓ **143/542 = 26.4 %** of multi-intron isoforms are emitted at ≥ 2 copies (53 groups) |
| P2 | ⛔ **13/53 = 25 %** of those groups have a copy without evidence (predicted ≥ 50): most multi-copy isoforms are backed at EVERY copy by unique mappers — paralogs share exon structure, and their unique mappers say so |
| P3 | ✓ **43 genuinely shared groups** (evidence at ≥ 2 copies), up to 7 copies for one chain (9/10/13/14/18/20/21) |
| P4 | ✓ 14 phantom copies, **13 with a sibling ≥ 0.985** in the group: copies 7 (5), 6 (4), 23 (4), 24 (1) — all `undecidable` transcripts |
| P5 | ✓ 0 `assigned` transcripts are phantoms |
Extra: **39 single-copy isoforms sit on abstaining reads only** (an address the aligner's tie-break chose); member
reads of multi-intron transcripts: 4,029 unique / 2,806 abstaining / 109 assigned. The phantom problem is
therefore small and sharply localized — **53 transcripts (14 + 39) of 542 have no evidence-backed address**, all in
the near-identical groups — and the shared-isoform question has a measured answer: paralogs do share isoforms
(43 groups), and at the twin copies the reads cannot say which twin(s), which is what step 2 pools for.

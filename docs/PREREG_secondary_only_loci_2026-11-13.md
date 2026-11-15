# PREREG row 1671: loci supported only by secondary alignments, real or shadow (2026-11-13)

User request: loci found only when secondary alignments are considered are probably shadows of other loci, but some
may be unfairly excluded because no tied multimapper got the primary there. Check which are real and, if real ones
exist, derive rules that certify them. Assembly layer (which loci are emitted); no read is assigned to a copy.

**Truth.** The annotation alone cannot decide this: a shadow lands on a paralog and paralogs are annotated (row 1060:
97% of secondary-seeded loci are annotated). The annotation gives the class of the locus (coding / pseudogene / other
gene / none). Whether the molecules come from the locus is read from the animal's own haplotypes (KB3781 = Jim, mat and
pat assemblies; method of `bench/o2_truth/haplotype_truth.py`).

**Step 0 (this row, descriptive, KB3781 default assembly `r1665_monofeed/k_base/k.gtf`).**
`bench/famcmp/sec_only_loci.py` (selftest). Population = assembly loci with no primary record on an exon and >= 1
molecule with a secondary record there inside the default pool (AS >= 0.98 x the primary's). Per locus up to 30 molecules
(hash-sampled) and the longest transcript of the locus L and of each assembly locus P holding those molecules' primaries
are aligned to both haplotypes. Molecule class: L (matches L's haplotype loci, not its P's), P, both (the haplotypes
cannot tell them apart), other. Locus class = the class of >= half of >= 3 judged molecules, else mixed / thin.
Reading: L = real and unfairly excluded; both = real candidate (identical copies, the aligner's pick is arbitrary);
P = shadow.

**Features recorded before labels are seen** (candidates for a certificate, none preferred): share of molecules whose
record here ties the primary's AS exactly; share with NM here <= NM at the primary; median AS ratio; molecule count;
share of molecules whose primaries sit on one locus; spliced or not; annotation class.

**Stop rule.** Loci of class L or both < 10 -> report the table, nothing to certify, row closes.
Otherwise a second registration freezes a rule on the odd-numbered chromosomes (by BAM header order) with bars on fixed
counts (certified L/both vs certified P), tests it on the even ones, then on HG002 (truth in `r1644_hg002_unspliced`);
OR6737 and A119b by annotation only. A rule would add an attribute; it would not add or drop models by default.

Not for `advisor_share` until the user says so. Repeats (RepeatMasker) are not classed in step 0.

## Result (2026-11-13, jobs 4469111 and 4469120, `/anvme/workspace/v138da13-rustle/r1671_seconly/{kb,kb0}`)

KB3781: 14,860 assembly loci (merged by span), 194 with no primary record on an exon; all have a pool molecule.

| truth | L | both | P | other | mixed | thin |
|---|---|---|---|---|---|---|
| as registered (NM slack 2) | 2 | 24 | 139 | 7 | 1 | 21 |
| strict, slack 0 (post hoc, `SEC_ONLY_SLACK=0`) | 0 | 8 | 155 | 10 | 0 | 21 |

By annotation (registered run): coding 152 (L 1, both 19, P 114), pseudogene 29, other gene 2, none 11. Annotation does
not separate the classes.
Stop rule as registered: L + both = 26 -> GO. Strict re-read: 0 loci where the molecules match the locus better than
their primaries' locus, 8 where the haplotypes cannot tell the two apart -> under the bar of 10. The two L loci of the
registered run (5 and 8 molecules) are an effect of the slack.
Features: no recorded feature separates `both` from `P` (NM here <= NM at the primary: median 0.00 in both; exact tie
0.43 vs 0.67; primaries on one locus 1.00 in both). No rule derived; the second registration is not written.
Reading: on this sample the default assembly emits ~155 shadow loci of 194 that rest on secondary alignments alone, and
no locus was found where the aligner's primary choice hid a real origin. Limits: one sample, own-genome reference
(a cross-individual reference could differ); loci merged by span, so nested ones are not counted; loci never
assembled (secondaries outside the 0.98 pool) are not in the population.

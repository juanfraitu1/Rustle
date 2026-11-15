# PREREG row 1672: flag for loci supported only by secondary alignments (2026-11-13)

User request ("register it") after row 1671. Assembly layer; no read is assigned to a copy; nothing is dropped.

**Rule (frozen, no parameter):** an assembly locus (merged same-strand overlapping transcripts) with no primary record
on any exon is flagged "supported only by secondary alignments". Found on KB3781 (194 loci: strict truth shadow 155,
indistinguishable 8, real 0, elsewhere 10, thin 21).

**Held-out test: HG002** (not used for row 1671), default assembly `r1669_monofeed_hg/base/hg.gtf`, BAM
`hg002.pat.bam`, haplotype indexes `r1644_hg002_unspliced/{mat,pat}.mmi`, annotation `~/data/hg002/genes.gff3`,
`bench/famcmp/sec_only_loci.py` with `SEC_ONLY_SLACK=0` (the strict truth, fixed here before the run).

**Bars** (flagged loci with a truth class, i.e. not thin):
1. class L (molecules match the locus, not their primaries' locus) <= 5 loci and <= 0.03 of judged;
2. class P (shadow) >= 0.75 of judged;
3. validity: flagged loci >= 30 and judged >= 20, else the test is void.

PASS = the flag is wired as an output of `assemble` (a table of flagged loci; models and defaults unchanged), with a
driver test. FAIL = no wiring, row closes. The wording of the flag stays "supported only by secondary alignments";
it does not say "artifact", since the indistinguishable class is real candidates.
Not for `advisor_share` until the user says so.

## Result (2026-11-13, job 4469136, `r1672_shadowflag/hg/loci.tsv`) — FAIL

HG002: 1,518 of 19,893 default loci have no primary record on an exon; 970 judged, 548 thin.

| annotation | L | both | P | other | mixed | thin | total |
|---|---|---|---|---|---|---|---|
| coding | 0 | 2 | 47 | 32 | 0 | 1 | 82 |
| pseudogene | 2 | 2 | 85 | 17 | 2 | 1 | 109 |
| other gene | 1 | 2 | 3 | 20 | 1 | 25 | 52 |
| none | 103 | 82 | 98 | 311 | 160 | 521 | 1,275 |
| total | 106 | 88 | 233 | 380 | 163 | 548 | 1,518 |

1. class L 106 loci, 0.109 of judged (bar <= 5 and <= 0.03): FAIL.
2. class P 233, 0.240 of judged (bar >= 0.75): FAIL.
3. validity 1,518 flagged, 970 judged: valid.

No wiring; the row closes. "No primary record" is not a shadow mark on HG002: 106 flagged loci carry molecules
that match the locus and not their primaries' locus under the strict truth. Those 106 are single-exon
(median spliced models 0), 103 unannotated, median 5 molecules, exact AS tie 0.75. The annotated flagged loci
(coding + pseudogene, 191) do behave as on KB3781: 132 shadow, 2 real.
Not checked: why HG002 flags 7.6% of loci against 1.3% on KB3781 (the BAM is aligned to the paternal haplotype
only, so maternal-origin molecules have no true locus to be primary on; this may produce both the "other" class
and part of L); the 106 were not inspected individually. No rule was fitted after seeing this.

# Row 1686 — step 0 ceiling: do segmental duplications hold the gene pairs fused by readthrough? (2026-11-13)

Goal (user, 2026-11-13): explore other ways to include segmental duplications. Assembly layer; a count only; nothing
wired; no read is assigned to a copy. Last SD use not yet measured (rows 1682-1685 closed the others).

**Background.** Rows 1516/1517: 73 fused loci on A119b join paralogs or same-stem genes (69 spliced, 4 intronless
blobs); de novo cut rules reached 5-6 of 69 because nothing told a tandem paralog pair from two exons of one gene.
A tandem SD pair with one mate over each gene would be that signal, and its boundary the cut point.

**Count (frozen).** Table `~/fused/tandem_fused.tsv`; gene spans from the A119b benchmark RefSeq GTF (CHM13); CHM13
SEDEF. A fused locus is "SD-held" when one SEDEF pair has a mate over each of its two listed genes (>= 300 bp each,
`edgesd.sd_id`).
**Go:** >= 35 of the 69 spliced fused loci SD-held (half the class). Fewer -> STOP: SDs cannot locate the cut for
most of the class, and no cut rule is built.

## Result (2026-11-13, `r1682_sdid/tandem_ceiling.out`) — STOP as registered (28 of 69 < 35)

Spliced fused loci: 28 SD-held, 41 not. Intronless blobs: 2 held, 2 not. All 73 loci had coordinates for both genes.
Held examples: POLR2J2/POLR2J3 (0.998), TBC1D3D/E (0.996), GAGE1/GAGE12D, RBMY1A1/B, HSPA1A/B, CDK11A/B, UGT2B15/17,
LILRA6/LILRB3, KLRC1/3 and 13 KRAB-ZNF pairs at 0.83-0.93.

The registered bar (half the class) is not met, so no cut rule is built under this row. Stated plainly for the user:
28 is five times the reach of rows 1516/1517 (5-6 of 69), so SDs do carry the "two units" signal for two fifths of
the class; the other 41 are neighbouring genes of one family that no SEDEF pair links (diverged paralogs: readthrough
between distinct genes). Whether a cut rule on the 28 is worth a registered row is the user's decision; it would need
a count of SD-held gene pairs inside correctly unfused loci as its false-cut control.

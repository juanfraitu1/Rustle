# PREREG r1524: recover the base-mode family misses (strand vote, pass-1 drop, catalog-family cut)

User's goal text (2026-11-04), verbatim:
> Mechanism: (1) single-exon models take their strand from the read-orientation vote (RUSTLE_READ_STRAND on) instead of the
> "+" placeholder; (2) find and fix the pass-1 stage that drops a two-read canonical intron chain with no competing model
> (PRAMEF8, chr1:12700000-12740000, two identical MAPQ 22 primaries, polish off, primaries only); (3) the fused-model cut
> from r1523 runs on catalog families (x.fam.copies.tsv), not truth genes, as a final.gtf layer.
> Panel: r1521, 420 primates-level families, A119b, judged under T, plus the 34 missed members with a >= 2-read local chain
> (~/famcmp/assemblable.tsv).
> Bars: of the 34, >= 20 gain a same-strand model; complete families >= 204 + 8; fused families <= 40; ordinary-member hit
> rate >= 0.819 (base value); genome-wide spliced chains identical to base (bar 5 of r1518, >= 0.99); wall: one genome-wide
> assembly run.
> PASS -> (1) and (3) default on in final.gtf layers, (2) fixed in base.
> FAIL -> each part judged alone; a part that passes its bars ships alone.

## Fixed before any number
- Base = fusedgw/A/x.gtf (r1518 arm A, layers off). Arm = one genome-wide assemble run with (1) on and (2) fixed, layers off,
  then (3) applied to its GTF; the cut layer is also scored on base alone so each part is judged alone.
- "Same-strand model" for the 34 = any model (any exon count) on the gene's strand overlapping its exons (r1521 Amendment 2).
- The 34 are listed in ~/famcmp/assemblable.tsv (ours = 0, maxP_chain >= 2); the 3 tied-only members are not in the 34.
- Bar 5: spliced intron chains of the arm GTF vs base, identical share >= 0.99 of base chains (bench/seed_unspliced_score.py bar 5).
- Panel scorer: bench/famcmp/family_panel.py, hit = any model (default), definition T.
- (2) is a code fix only if its cause is found; the fix must keep the 17 family tests and add one test from the PRAMEF8 reads.

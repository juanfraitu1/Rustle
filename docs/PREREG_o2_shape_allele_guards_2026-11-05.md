# PREREG O2 SHAPE GUARD (register 1540) AND ALLELE-SITE GUARD (register 1541), 2026-11-05, before any number

Why (user 2026-11-05 "register and try the shape guard and the allele versus paralog case"): the two failure classes the
r1514 ledger names and never measured. Both are measured as POST-FILTERS on the r1514 run (KB3781, 494 families, `--tie-ms 10
--primary-on-family`), re-run only to add `--dump-psv --dump-star` (assignments expected byte-identical: dump flags); the strict
truth is r1514's (`pof_heldout.strict3_dump.tsv` 151 decidable: O2 129 right / 22 wrong; `pof_panel.strict3_dump.tsv` 419:
393 / 26). A guard that is worth porting into the binary (where it would also refit the PSV model) is a later row.

## r1540 shape guard
Rule: an assigned molecule whose record on the family windows is SPLICED (an N op >= 30 bp in the primary record, from the BAM)
may not be assigned to a SINGLE-EXON copy (`n_exon` = 1 in `kb.fam.copies.tsv`); such a row becomes abstained. Nothing else moves.
Step 0 (counted first, reported): rows of the class among the decidable held-out and panel rows, split right / wrong.
Bars (held-out, then panel reported beside it): the guard removes >= 2 wrong rows AND removes more wrong than right rows.
Both hold -> port into the binary and re-register with r1514's bars; else the class is too small or the guard is wrong: closed.

## r1541 allele-site guard (the MCL183 allele-vs-paralog trap, r1413)
Mechanism: the mat allele of copy A matches copy B's reference at a PSV column, so the column is ALLELIC within this individual,
not a paralog difference. Blind rule (no haplotype assembly): a PSV column is ALLELIC for copy A when A's own UNIQUELY placed
reads (primary, MAPQ 60, not contested, overlapping the column; from the BAM) carry >= 2 alleles with >= 2 reads each and the
minor one at >= 20%. A column that is allelic for either candidate of a molecule's deciding pair does not count; a molecule whose
deciding columns fall below 2 (`n_decisive` as the binary counts them) becomes abstained; a molecule whose remaining columns
favour the other candidate is counted as MOVED (reported, scored as its new copy).
Step 0: among the 22 wrong and 129 right held-out rows, how many have >= 1 allelic deciding column; MCL183's rows listed.
Bars (held-out, panel beside): the guard removes or moves-to-right >= 2 wrong rows AND costs fewer right rows than it fixes.
Both hold -> port (the column filter belongs in the certificate) and re-register; else closed as a within-sample rule and the
haplotype-aware version (needs assemblies) stays unregistered.
Inputs: `<out>.psv_reads.tsv` (read x column base + assignment), `<out>.psv_cols.tsv` (column -> genome position),
`<out>.star_reads.tsv` (the deciding columns per molecule). Script: `bench/o2_contested/guards.py`.
Scope limits: post-filters cannot show the PSV-model refit effect r1514 bar 3 found; the truth is slack 0 on one individual.

# PREREG row 1565: parallel BAM pass for `--polish-mono-primed` (SnapBam::mono_primed)

## Why (r1564 stamps, read before this row)
After r1564 the assemble stage is 1,685 s; the serial per-contig polish loop is ~1,170 s of it, and in BOTH runs ~1,000 s
of the loop is the mono-primed evidence pass (main: TSS PROOF -> MONO PRIMED 1,015 s; final_run: TES EXTEND -> MONO
PRIMED 1,001 s). `SnapBam::mono_primed` queries ONE region spanning [min, max) of a contig's single-exon models, i.e.
nearly the whole chromosome, and decodes every record with one reader (4 inflate workers; ~3 cores observed).

## Change
Split the query range into 8 equal sub-ranges, each on its own reader, in scoped threads. A record is counted by the
sub-range holding its alignment start (sub-range 0 also takes starts before the range), so every record of the old
single query is seen exactly once; per-record contributions to (reads, primed) are independent sums, so the output is
the same by construction. Nothing else changes (the loop over contigs stays serial).

## Run
Fresh build of main with the change; same BAM, genome, flags and harness as r1564 (shipped --region-threads 8, assemble
budget 20M), singlenode, wall limit 40 min.

## Bars (all must hold)
1. a1.gtf, a1.families.gtf, a1.final.gtf, a1.bridges.tsv, a1.bridge_junctions.tsv, a1.quant.tsv byte-identical to r1564.
2. Stage wall <= 1,100 s (r1564: 1,685 s).
3. RSS <= 120 GB.

## Prediction
PASS; mono-primed pass ~1,000 -> ~200-300 s per run (uneven read density across equal spans, centromere gaps), stage
wall ~950 s.

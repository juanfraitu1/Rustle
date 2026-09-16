# PREREG — what junction-jitter tolerance should RUSTLE_JUNCTION_FUZZ_BP use? (2026-09-15, written before any gffcompare/SQANTI3 score from this feature is computed)

## Question
On real chr20 IsoSeq reads, how far do individual reads' own CIGAR-derived splice junction coordinates scatter around the nearest real annotated RefSeq intron boundary -- and what tolerance, chosen by ONE fixed rule stated here, should `merge_fuzzy_skeletons` (docs/superpowers/specs/2026-09-15-fuzzy-junction-merge-design.md) use?

## Method
`bench/measure_junction_jitter.py` against `/mnt/linuxdisk/home/juanfraitu/bakeoff/human_chr20/chr20.bam` and the same substrate's `chr20_ref.gtf`. Every spliced primary read's own junction is matched to the nearest annotated intron by donor-site distance (capture window 500bp); the signed offset at both donor and acceptor is pooled into one distribution.

## Result
reads scanned: 25341
junctions seen: 101697
junctions matched to an annotated intron within 500bp: 100479
pooled donor+acceptor |offset| samples: 200958
median |offset|: 0
90th percentile |offset|: 0
95th percentile |offset|: 672

## Decision (fixed now, before Task 5 runs)
`RUSTLE_JUNCTION_FUZZ_BP` = 0 (the 90th percentile value above, an integer number of base pairs). This rule (90th percentile of real, freshly-measured jitter) was chosen in the design spec BEFORE this script ran and is not changed based on Task 5's result.

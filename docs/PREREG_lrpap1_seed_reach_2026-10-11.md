# How many seeds does the DNA mode need to find every LRPAP1 copy? (gorilla mGorGor1), 2026-10-11

**Pre-registered before any seed is aligned.** User (2026-10-11): the advisor is not sure the DNA (guided) mode finds all copies; show, for the
gorilla LRPAP1 expansion, how many copies are found from one or more of the other copies used as seeds, i.e. how many seeds are needed.

**Truth.** The 11 gorilla LRPAP1 loci of the LRPAP1 artifact (`/mnt/linuxdisk/tmp/lrpap1/browser/lrpap1_browser.json`): 8 full-length copies
(c00-c07) and 3 copies of the 5′ three-quarters (p12, p14, p16); 10 of them are Yoo 2025's list, the 11th (p16) was found 2026-10-04.

**The DNA-mode finder (as shipped in `bench/guided_pipeline.py`, unchanged).** A seed copy contributes its transcript (its model's exons,
spliced; `minimap2 -x splice -N 100 -p 0.1` on the primary genome index `GGO.splice.mmi`) and its gene body (model span, first to last exon;
`-x asm20`, same flags). Passing hits: `transcript_hits` (identity >= 0.80, >= 50% of the seed aligned) and `gene_body_chains` (identity >= 0.80,
aligned >= 50% of the shorter of seed and extrapolated span). Seed i FINDS copy j when a passing hit of i overlaps j's span on its chromosome.
Hits outside the 11 are counted per seed (candidate burden), not judged.

**Reported.** The 11 x 11 finds matrix; for every seed set size k = 1..11 and every set of that size: copies found in ONE pass (union of the seeds'
finds) and with ITERATION (found copies become seeds until nothing new: reachability in the finds graph); the minimum number of seeds that finds
all 11 in one pass (smallest set cover) and with iteration (number of source components of the finds graph).

**Expectation (registered).** E1: with iteration, ONE seed suffices, whichever copy is chosen (the copies are 95.7-99% identical, far above 0.80).
E2: in one pass, every full-length seed finds all 11; a fragment seed finds every copy carrying its exons 1-5, i.e. all 11.
If either fails, the copies a seed misses and why (coverage or identity) are reported; no threshold is changed after looking.

# DNA vs RNA, one engine: reproducing Soto's families from the genome vs from transcripts

**Date:** 2026-07-27. Benchmark: Soto 2025, 83 families / 362 members on T2T-CHM13v2.0 (`bench/soto/80_fams.chr.bed`).

## The result

| substrate | front-end | member recovery vs Soto |
|-----------|-----------|-------------------------|
| **DNA** — genome only, no annotation, no reads | self-align windows → genomic-locus reps | **358 / 362 = 98.9%** |
| **RNA** — human A119b Iso-Seq | reads → skeletons → spliced reps | **306 / 362 = 84.5%** (strict own-copy); 313/362 = 86.5% counting the coverage-recoverable members |

Both numbers come from the **same grouping engine** — `denovo_pipeline::homology_blocks`
(E_r homology edges → γ-quasi-clique partition). The RNA path (`--homology-primary`) and the DNA path
(`--from-genome`) call the identical function; the *only* thing that differs is the front-end that produces
the reps it groups:

- **RNA rep** = a spliced transcript (`seq` = exon-sum, introns removed).
- **DNA rep** = a genomic locus (`seq` = the full locus incl. introns, empty intron chain).

Reproduce: `gw_family_catalog --from-genome bench/soto/80_fams.chr.bed --fasta <chm13v2.0.fa> --min-identity 0.90 --out dna_mode`
then `python3 bench/soto/score_from_genome.py dna_mode.copies.tsv` (see `bench/soto/run_from_genome.sh`).
The run: 362 windows → 577 duplicated-locus reps → 49 families.

## What this shows

Give the *same method* the genome and it reproduces ~99% of Soto; give it RNA and it reaches ~85%. The
family-grouping engine is not the bottleneck — it is correct on both substrates. **RNA's shortfall is the
substrate**: splicing removes the intronic/flanking sequence that distinguishes near-identical copies, so
RNA reps are confined to the most-conserved exons and the copies collapse. DNA reps carry that sequence
(measured elsewhere: NCF1 copies differ at ~2 exonic positions but ~75 genomic positions), so they separate.

That is the concrete form of "the advisor is asking a different question": Soto (and this DNA mode) classify
copies the **assembly already separated**; the RNA method must *deconvolve* copies from ambiguous reads — a
strictly harder problem, bounded by identifiability, not by the engine.

## What this does NOT show (read before quoting the 98.9%)

- **Not "RNA copy-detection = Soto."** Impossible by construction — Soto consumed an assembled genome +
  WGS copy-number; RNA has neither. No experiment can deliver that; this one does not claim it.
- **The genome hands the method pre-separated copies.** DNA mode reaching 99% is *expected*: it consumes the
  same substrate Soto did (an assembled genome where the copies are already distinct loci). The value is the
  *contrast* with RNA under the same engine, not the absolute DNA number.
- **Region-window circularity (disclosed).** DNA mode uses Soto's family *regions* as search windows (scope:
  the 83-family benchmark). It discovers the copies + families *within* those windows from sequence, not from
  the member list; the Soto BED is used only to score. Genome-wide discovery (no window prior) is a later
  extension.
- **Member recovery, not family-count reproduction.** DNA mode emitted **49** families (30 cross-chromosome)
  covering 358 members — it does not reproduce Soto's **83** family *partition*. Sequence homology alone
  over-merges near-identical paralogs across loci (the same tendency the de-novo DNA reconstruction showed,
  `dna_family_pr.tsv`: 81.5% with 13 over-merges). Soto splits those with **parCN** (paralog-specific
  copy-number from WGS read depth) — the ingredient neither this DNA mode nor RNA has. So "98.9%" is
  *member recovery*, and the family granularity is coarser than Soto's.

## Bottom line

One engine, two substrates: **DNA 98.9%, RNA 85%.** The gap is the question, not the method — and the one
axis where RNA is not redundant with the genome is transcript-level *resolution* of the copies (read PSVs /
junctions), the analog of Soto's parCN.

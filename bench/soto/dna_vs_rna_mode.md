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

## DNA and RNA recover DIFFERENT members (same-engine cross-tab)

Scoring each of the 362 members by whether the **DNA-mode** pipeline recovers it (≥2-copy family overlap)
and whether the **RNA** pipeline recovers it (verified attribution = FOUND / resolved-elsewhere / coverage-recoverable):

| | count |
|--|------:|
| both recover | 310 |
| **DNA-only** (in the genome, RNA can't see) | **48** |
| DNA-pipeline gap (in the genome; RNA caught, `--from-genome` missed) | 3 |
| neither yet | 1 |

**RNA is a subset of the genome — it cannot recover what DNA lacks.** So the genome is the superset (existence
~100%); the DNA-mode *pipeline* recovers 358/362; RNA recovers the expressed-and-resolvable 313/362. The "3" are
**not** RNA-exclusive: they are large multi-copy genes (TCAF2, NPIPB2, NPIPB15) present in the genome that the
`--from-genome` pass failed to group — a fixable implementation limit of its self-alignment preset — plus one
218 bp member (AC239809.1) below the discovery floor.

The 48 DNA-only members — in the genome but RNA-invisible — by *why RNA missed them*:

| RNA-miss reason | DNA-only members |
|-----------------|-----------------:|
| not-expressed (silent copies) | 22 |
| mis-chain | 8 |
| collapse-K0 (indistinguishable sibling) | 6 |
| genuine-miss | 6 |
| seeding-gap | 5 |
| thin-single-exon | 1 |

The DNA advantage is **exactly** the RNA-invisible set: silent copies (22 not-expressed — the genome carries
them; RNA can't see what isn't transcribed) plus K=0 collapses. The genome is a superset of RNA — the 3
"DNA-pipeline gap" members are in the genome (the pipeline missed them, not RNA reach), so RNA never recovers
anything DNA lacks.

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

# PREREG r1788: assembly recovery split by Ensembl transcript support level, human A119b (2026-11-18)

Written before any count of this row. Workspace `r1788_tsl_split`, script `bench/famcmp/tsl_split.py`. A reporting split:
it changes no rule and no default.

## Question
User 2026-11-18: is the Ensembl transcript support level (TSL) useful, and where. Proposed use: report recovery of
reference transcripts by TSL on a human sample, so a missed model with no supporting evidence weighs less than a missed
well-supported one.

## Data
- Reference: public CHM13 v2.0 annotation (CAT + Liftoff of GENCODE v35, `chm13.draft_v2.0.gene_annotation.gff3`), which
  carries `transcript_support_level` per lifted transcript. This is not the RefSeq truth of the assembly tables; counts
  are not comparable with them.
- Sample: A119b (human testis, CHM13). Spent for assembly rules; nothing is tuned here.
- Arms (existing GTFs, nothing rerun): Rustle default-day runs of `artifact_refresh_2026-11-16` (`h.pri/a1.gtf`, and its
  `a1.final.gtf`), StringTie, FLAIR, IsoQuant as in `after_h.sh`. The other tools are trusted baselines.

## Counting
- Unit: a multi-exon reference intron chain (chrom, strand, introns). Several transcripts with one chain count once, at
  their best (lowest) level. Classes: 1, 2, 3, 4, 5, NA (annotated "NA" or absent), novel (CAT models with no source
  transcript). Single-exon models are left out (TSL does not grade them).
- Recovered = an arm has a transcript with exactly that intron chain.
- Expressed-gene restriction: reference chains of genes in which at least one arm recovers at least one chain. The
  reading below is on this subset, since level 1 models are also the highly expressed ones.

## Reading rule (fixed)
- The split is INFORMATIVE if, for Rustle on the expressed-gene subset, recovery at level 1 is at least twice recovery
  at level 5 and recovery does not rise from level 1 to level 5 at any step by more than 2 points. Then TSL is proposed
  as a standing column of the human benchmark tables (user decides).
- NOT INFORMATIVE if level 1 recovery is under 1.5 times level 5: the split adds nothing.
- Between: reported, no proposal.
- Classes with fewer than 200 chains are not read.

## Predictions
1. Recovery falls from level 1 to level 5 for every arm; INFORMATIVE.
2. The order of the arms by recovered chains is the same at every level.
3. Of the reference chains an arm misses in expressed genes, more than half are level 3-5, NA or novel.
4. FLAIR recovers the largest share of level 4-5 chains (it keeps the most models).

## Limits known beforehand
TSL is GENCODE v35 evidence from mRNA / EST alignments, human only; it says nothing for apes, single-exon models or
copy assignment. Exact intron-chain identity only; ends are not scored.

## Amendment 1 (2026-11-18, before any HG002 count): second human sample
User: "should we do this for the other human sample too?" Same script, reference, counting and reading rule on HG002
(Kinnex, CHM13, cell line). Arms, existing de novo GTFs: Rustle `r1692_trunc5_hg/hg.gtf` (and `hg.final.gtf`), StringTie
`apetools/HG002_kinnex.t2t/st_dn.gtf`, FLAIR `fl_dn.isoforms.gtf`, IsoQuant `iq_dn/x/x.transcript_models.gtf`.
Purpose: does the A119b shape hold outside testis. Predictions: (5) level 1 is at least 4 times level 5 for Rustle;
(6) level 5 is again above level 4 in every arm; (7) the level 1 share is higher than on A119b for Rustle (a cell line
expresses fewer oddities). If (6) fails, the broken ladder is a testis effect and the "between" verdict is reopened.

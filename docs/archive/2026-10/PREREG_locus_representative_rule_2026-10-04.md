# PREREG — the locus representative must carry the locus's structure: most-junctions vs most-reads (written before any run, 2026-10-04)

## Why (from `docs/archive/2026-10/SPLICED_COPY_SUPPORT_2026-10-04.md`, register rows 1232-1233)

The de novo locus (a `gene_id` group of the assembled GTF) is represented downstream by ONE transcript — the one with the most reads
(`mcl_families --from-gtf`, ties: longer span, then last `transcript_id`; `src/bin/mcl_families.rs` "GtfLocus"). Its exons are the locus's
exons in `loci.gff3` (the family step's exonic filter and edge weight), the copy's exons and spliced sequence in the copy table O2 reads,
and the structure every recovery claim is judged on. In a 5'-truncated library the most-read transcript is a 3' fragment: at NPIPA9 the
locus holds 42 transcripts / 38 read-supported junctions and its representative has 1; under the strict "found" rule the NPIP page's
23 / 21 / 24 own nodes are 13 / 9 / 7 found copies while the loci themselves contain 20 / 19 / 16. The 09-27 representatives study
(rows 1121-1127) measured the rep's exon-bp coverage (98.7-100% of the locus's exon bp carried) but never its junction count; its
"rep fine for families / O2" does not hold at NPIP. This prereg tests the one obvious alternative.

## Rules

- **R_M (shipped):** representative = the transcript with the most `reads`; ties: longer span, then the last `transcript_id`.
- **R_J (candidate):** representative = the transcript with the most junctions (gaps >= 50 bp between consecutive exons — the same
  junction floor as the strict "found" rule); ties: most `reads`, then longer span, then the last `transcript_id`. Rationale: the
  assembler already admits only read-supported junctions (floor 2, strict / majority junctions), so the transcript with the most junctions
  is the locus's most complete read-supported chain; no new constant.
- Not run: the exon UNION of a locus's transcripts as its representative (retained-intron isoforms would fuse exons into one giant exon,
  register T5), and the longest-span transcript (a readthrough or a retained intron wins). Both are named so they are not re-proposed
  without a reason.
- Implementation: `mcl_families --representative most-reads|most-junctions` (default `most-reads`; unset = byte-identical to today),
  driver knob `RUSTLE_REPRESENTATIVE=most-junctions` on `tools/rustle_pipeline.sh families` (the pattern of `RUSTLE_MIN_COV_SHORTER`).
  Nothing else in the family step changes: the all-vs-all aligns whole locus bodies, so `loci.paf` is the same under both rules; only the
  exonic filter / edge weight, `loci.gff3`, the copy table and the units move.

## Substrates (de novo mode only; the guided mode's representatives are the annotation)

- **Dev (rules were developed here; reported, never decides):** human A119b chr16 (NPIP; the page), gorilla OR6737 NC_073244.2.
- **Held-out (decides):** human chr6 (untouched), chr2, chr8, chr10 (the reused verdict set of Figure 7), gorilla NC_073234.2 (untouched).
- Inputs: the per-contig de novo GTFs of the Figure 7 dev tables (`/mnt/linuxdisk/tmp/rustle_figures/fig7/current/<species>_<contig>.denovo.gtf`,
  the genome-wide assembly restricted to the contig), the driver's `families` stage run twice per contig (R_M / R_J) with the same binary
  and the stage's shipped defaults (`--min-exonic-bp 1 --min-shared-exon-frac 0.60 --emit-units`, `--bridge-regroup f1v2`,
  `--min-cov-shorter 0.70`), into a scratch dir (the figure caches are not touched); references and scorers exactly as Figure 7's
  (`figures/_o1_recovery.py`: `family_score` vs Compara and Soto for human, the Liftoff self-lift pairs for every species, `npip_u2` on chr16).

## Pre-registered readouts

- **H1 (families, held-out):** per contig and rule, the bipartite F / sensitivity / precision against the primary reference (human:
  Compara; gorilla: Liftoff pair sensitivity, where precision is not defined) — Soto and the protein referee reported beside.
- **H2 (copy table structure, held-out):** per contig and rule, over the copies of `<prefix>.fam.copies.tsv`: median junctions per copy
  (`n_exon - 1`), fraction of copies with >= 2 junctions, total exon bp.
- **H3 (found genes, strict rule, held-out and dev):** copies = the contig's annotated protein-coding genes (human CAT/Liftoff v2.0,
  gorilla RefSeq) that are spliced-expressed in the sample (>= 2 reads with >= 2 read-supported junctions, `bench/copy_support.py`'s
  rule); found = a same-strand locus of the arm whose representative carries >= 2 of the gene's supported junctions. Reported: found
  under R_J vs R_M, the genes that change, and on chr16 the 25 NPIP copies (the page's question) with the locus-level reading beside.
- **Decision (held-out only):** R_J becomes the default iff, on every one of the 5 held-out contigs, (a) bipartite F against the primary
  reference >= R_M − 0.005 and precision >= R_M − 0.01 where defined; (b) H3 found genes >= R_M; (c) H2 median junctions per copy >= R_M.
  A single violation keeps R_M as the default and R_J stays opt-in; the dev results never override a held-out violation. The outcome,
  the per-contig tables and the flipped (or not) default go in `docs/archive/2026-10/LOCUS_REPRESENTATIVE_RULE_2026-10-04.md` and the register.
- Reported beside, no rule: run time of the two arms; the number of loci whose representative changes; O2 is not re-run here (the copy
  table it reads changes with the rule; its effect is a separate measurement).

## Not changed

Locus formation, node admission (single-exon loci stay), the family edge rule and its thresholds, MCL, the bridge regroup, O2's certificate.

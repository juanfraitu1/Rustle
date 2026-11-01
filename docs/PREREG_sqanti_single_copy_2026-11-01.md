# Pre-registration: SQANTI3 on single-copy genes, our plain assembly vs StringTie, FLAIR and IsoCon (A119b chr20, 2026-11-01)

**Written before any number.** User: "can we at some point make an analysis based on sqanti of how wrong it was [IsoCon] with our
method, vs stringtie and flair for single copy transcripts (... the mode that assembles without need of forming gene families)".
Register row 1482. Run dir `~/sqanti_sc_A119b/`, code `bench/sqanti_single/sc_panel.py`, SQANTI3 6.1 from bioconda in the pixi
project `~/tools_home/sqanti3/`.

## What is known (not predicted)
- Register 1478 (OR6737): IsoCon modeled 27 of 30 single-copy genes vs ours 30, StringTie 29, FLAIR 30; junction precision 0.88 vs
  ours 0.83 with 1.6-1.8x our transcript count.
- SQANTI3 5.5.4 on a six-chromosome testis panel (other library; `docs/archive/2026-10/SQANTI3_CAT_2026-10-01.md`): ours FSM 47.4%,
  artifact 6.1%, PASS 93.9%; StringTie 45.7% / 7.1% / 93.3%; FLAIR 32.1% / 12.1% / 71.6%. IsoCon was never run through SQANTI3.

## Universe (A119b, chr20; every eligible gene, no sampling)
CAT/Liftoff protein-coding gene with an Ensembl id, no Compara paralog at any level (`famtruth.truth_pairs`), expressed (a de novo
locus of `a1.gtf` overlaps its exons, `famtruth.expressed_from_gtf`), overlapping no copy of `a1.fam.copies.tsv` (locus columns), and
50-5,000 BAM records on its span (`isocon_targets.py select`'s depth window). The gene SPAN (first exon start to last exon end) is the
unit: each arm's GTF is cut to the transcripts overlapping a universe span on the same strand or unstranded, so the four arms are
scored on identical loci.

## Arms
- **ours** = `a1.gtf` (the genome-wide plain assembly, `rustle_pipeline.sh assemble`, 2026-10-29; family copies are a tag, no
  family step is needed to produce it).
- **StringTie** 3.0.1 `-L`, no annotation, and **FLAIR** 3.0.1: the lab's own genome-wide A119b runs on the same BAM
  (`~/work/jfim/rustle_data/analysis/benchmark_collapse/A119b.{stringtie,flair.isoforms}.gtf.gz`).
- **IsoCon** 0.3.3 per gene: pool = every primary and secondary record on the span (`isocon_targets.py pool`), `IsoCon pipeline`
  defaults, final candidates mapped back with `minimap2 -ax splice:hq` (`isocon_targets.py gtf`), the 1478 protocol unchanged.

## Scoring
SQANTI3 6.1 `sqanti3_qc.py` against the CAT chr20 GTF (`bench/annotation/cat_ref_gtf.py`) and `chm13v2.0.fa`, then
`sqanti3_filter.py rules` with the default rules, per arm on its cut GTF. Reported per arm, over the universe: transcripts, genes
with >= 1 transcript, FSM / ISM / NIC / NNC / artifact (genic + antisense + intergenic + genic_intron + fusion) counts and shares,
rules PASS share, transcripts per modeled gene, and per gene the FSM count so a reader can see the 1478 genes' shape.

## Decision rule (fixed now)
An arm is **WORSE than ours** when its FSM + ISM share is below ours by > 5 points, or its artifact share is above ours by > 5
points; **BETTER** when FSM + ISM share is above ours by > 5 points and artifact share is not above ours by > 5 points; **TIED**
otherwise. Genes modeled is reported, not a bar (every tool models nearly every single-copy gene). Nothing is adopted or changed
on any outcome: this is a report that places IsoCon on the same SQANTI axes as the other three tools.

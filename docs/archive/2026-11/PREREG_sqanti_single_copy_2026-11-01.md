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

## RESULT (2026-11-01, run once; `~/sqanti_sc_A119b/{targets.tsv,table.tsv,sq3/}`; jobs 4403511 (IsoCon, 370 pools in 22 min), 4403550 (mapping + SQANTI), srun 4403625/4403636 (IsoCon arm alone))
370 universe genes (1,508 chr20 CAT genes -> 545 protein-coding -> 543 with an Ensembl id -> minus 37 with a paralog, 89 unexpressed,
8 overlapping a copy, 39 outside the depth window). Cut GTFs: ours 3,259 of 5,319 chr20 transcripts, StringTie 2,928 of 5,731, FLAIR
11,218 of 20,917, IsoCon 19,067 of 23,872 placed candidates.

| arm | transcripts | genes modeled / 370 | FSM | ISM | NIC | NNC | artifact | FSM+ISM % | artifact % | rules PASS % | tx / gene | genes with >= 1 FSM | distinct ref tx matched | distinct intron chains | mono-exon |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ours | 3,259 | 370 | 631 | 353 | 765 | 1,346 | 164 | 30.2 | 5.0 | 70.5 | 8.8 | 320 | 631 | 3,259 | 10 |
| StringTie | 2,927 | 369 | 533 | 147 | 683 | 1,024 | 540 | 23.2 | 18.4 | 58.5 | 7.9 | 297 | 533 | 2,927 | 208 |
| FLAIR | 11,218 | 370 | 668 | 254 | 2,139 | 5,284 | 2,873 | 8.2 | 25.6 | 47.8 | 30.3 | 296 | 641 | 11,210 | 2,388 |
| IsoCon | 19,067 | 363 | 5,065 | 7,983 | 2,027 | 1,892 | 2,100 | 68.4 | 11.0 | 74.2 | 52.5 | 300 | 645 | 9,153 | 5,484 |

**Verdict by the registered rule (ours = FSM+ISM 30.2, artifact 5.0):**
- **StringTie WORSE** (FSM+ISM 23.2 < 25.2; artifact 18.4 > 10.0).
- **FLAIR WORSE** (8.2; 25.6).
- **IsoCon WORSE** (artifact 11.0 > 10.0), although its FSM+ISM share, 68.4, is far above ours. The rule's two clauses disagree and
  the artifact clause decides, as written.

**Reading (after the verdict).** IsoCon's FSM+ISM share is a per-candidate share over 52 candidates per gene: 19,067 candidates carry
9,153 distinct intron chains (half are redundant), 5,484 are mono-exon, and 7,983 (42%) are ISM = fragments of a reference transcript.
On per-gene and per-reference measures the four tools are close and ours leads: genes with >= 1 FSM 320 (ours) / 300 (IsoCon) / 297 /
296; distinct reference transcripts matched 631 / 645 / 533 / 641. IsoCon is not "more wrong" than ours on single-copy genes in the
sense of the 1478 question: it recovers the reference structures about as well, at 6x the transcripts, with twice the artifact share
and 7 genes without any model (1478's shape, 27 / 30). Our own FSM+ISM share (30%, NNC 41%) is lower than on the 2026-10-01 SQANTI
panel (47%): different library (A119b testis), reference (CAT chr20) and universe; not comparable across panels, only across arms here.

**Procedure notes.** (1) `isocon_targets.py gtf` re-indexed the whole genome per gene (26 genes in 50 min); job 4403511 was cancelled
after IsoCon and the mapping re-run from a prebuilt `splice:hq` index in job 4403550 (IsoCon outputs untouched; the step now says so
in its usage). (2) IsoCon transcript ids contain `;` `,` `(` `)`, which gffread (inside `sqanti3_qc.py`) rejects; the writer now
replaces them with `_` and the IsoCon arm was re-run alone on the cleaned GTF (`isocon.cut.clean.gtf`), twice (the first srun
launched 4 tasks on the same directory; cleaned and re-run with `-n 1`, identical numbers). Nothing adopted or changed.

## gffcompare on the same panel, updated 2026-11-01 (descriptive; after registers 1485 / 1487)
Reference = the 370 panel genes' CAT transcripts (`~/sqanti_sc_A119b/gffc/ref370.gtf`), each GTF cut to the universe spans as above;
`ours` = the panel's `a1.gtf`; `ours_tssr` = chr20 re-assembled with today's default (`--polish-tss rescue`, r1485); `ours_med` = that
plus `--polish-tss-end median` (r1487, opt-in). Sensitivity / precision in %; the last three columns are the `=` matched multi-exon
transcripts, the share whose 5' end lies within 50 bp of the reference's, and the median 5' distance (`bench/assembly_ends/end_score.py`).
Table: `~/sqanti_sc_A119b/gffc/table.tsv`.

| arm | transcripts | base Sn / Pr | exon Sn / Pr | intron Sn / Pr | chain Sn / Pr | transcript Sn / Pr | locus Sn / Pr | missed exons | novel introns | `=` multi | 5' within 50 bp | median 5' |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ours | 3,259 | 71.4 / 54.7 | 57.4 / 63.4 | 70.0 / 70.3 | 22.9 / 19.2 | 22.7 / 19.1 | 87.9 / 71.2 | 11.8% | 4.4% | 623 | 43.0% | 63 bp |
| ours_tssr | 3,298 | 71.5 / 54.5 | 57.5 / 63.2 | 70.0 / 70.3 | 23.0 / 19.1 | 22.8 / 19.1 | 87.9 / 71.2 | 11.8% | 4.4% | 628 | 42.8% | 64 bp |
| ours_med | 3,297 | 69.6 / 58.7 | 58.7 / 64.8 | 70.0 / 70.3 | 23.0 / 19.1 | 22.8 / 19.1 | 87.9 / 70.9 | 12.1% | 4.4% | 628 | 69.8% | 26 bp |
| ours_all (r1485 + r1487 + r1488 + r1489) | 3,385 | 71.9 / 58.5 | (not run) | (not run) | 23.0 / 19.1 | 22.8 / 18.6 | 87.9 / 70.6 | (not run) | (not run) | 628 | 69.8% | 61.3% 3' |
| StringTie | 2,928 | 75.1 / 28.7 | 51.0 / 57.5 | 68.6 / 69.8 | 19.0 / 19.1 | 18.8 / 17.7 | 83.3 / 45.8 | 10.7% | 7.2% | 519 | 40.7% | 78 bp |
| FLAIR | 11,218 | 85.3 / 19.0 | 63.7 / 33.4 | 76.9 / 37.6 | 21.8 / 6.7 | 21.8 / 5.3 | 80.2 / 29.7 | 6.9% | 9.8% | 599 | 70.3% | 25 bp |
| IsoCon | 19,067 | 69.1 / 35.1 | 57.8 / 37.3 | 63.1 / 70.6 | 22.4 / 4.5 | 22.3 / 3.2 | 81.4 / 37.7 | 16.8% | 5.6% | 4,914 | 58.6% | 36.5 bp |

With the two end layers ours keeps the best transcript-, chain- and locus-level precision and now places 5' ends as well as FLAIR
(69.8% vs 70.3% within 50 bp), with the highest base precision (58.7). The remaining gap is base and exon sensitivity against FLAIR
(69.6 vs 85.3 base Sn; FLAIR emits 3.4x the transcripts). IsoCon's `=` count holds its redundant candidates (median ~4 per reference).
ours_med has one transcript fewer than ours_tssr in the cut: a moved start can leave the universe span.

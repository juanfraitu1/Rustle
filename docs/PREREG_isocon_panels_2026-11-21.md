# Row 1872: IsoCon gene panels for all six species (pre-registered 2026-11-21)

User (2026-11-21): add IsoCon panels for all species to the assembler comparison page; IsoCon has no guided mode and stays out
of the guided and with-secondary comparisons. Until now the panel existed for gorilla only (rows 1478, 1515, 1764). This row
runs the same procedure on the ten other samples. Descriptive (page content); no rule, no default. Workspace
`r1872_isocon_panels` (`isocon.sbatch`). Committed before any target was selected or any IsoCon run started.

## Panel, the same definition for every species (`bench/denovo_eval/isocon_targets.py select`, unchanged)
- **Multi-copy class:** every gene of the species' annotation whose symbol or description matches one of the four families of
  the gorilla panel (`famtruth.GGO_MAINSTAYS`: LRPAP1, MAPKBP1, JMJD7-PLA2G4B, SPTBN5). These are expansions described for
  gorilla (Yoo 2025); another species may carry one copy of each. Said on the page.
- **Single-copy class:** 30 protein-coding genes with no Compara paralog, with a Rustle locus, 50-5,000 alignment records on
  the gene, drawn with the fixed seed of the script.
- Annotation: RefSeq GFF of the species' T2T assembly (`genomes/TAG/TAG_genomic.gff`); human: CAT / Liftoff on CHM13
  (`chm13v2.0_CAT_Liftoff.slim.gff3.gz`), which has no descriptions, so the families match by symbol only there.
- One panel per species, selected on the first sample (testis; human A119b) and used for both, as for gorilla.

## Procedure per sample (as `apetools/isocon_kb.sbatch`)
Read pools per target from the BAM, `IsoCon pipeline` per pool, candidates placed with minimap2 (`isocon_targets.py gtf`), then
`toolcompare.py --targets` (gene copies modelled, junction precision against exon-block boundaries, spliced transcripts) for
Rustle (`final.gtf` of the 2026-11-21 refresh, primary pool), StringTie, FLAIR, IsoQuant, bambu (NDR 1) and IsoCon, and
`isocon_junc.py` (Rustle's junctions cut by support to IsoCon's count on copies both model, RefSeq introns as truth).

## Reported on the page, per sample
The three columns of the existing gorilla table, and the matched-count junction precision of Rustle against IsoCon. No bar.
Statements about IsoCon in the page's text are written from the twelve tables after they exist and are checked by the build
script (asserted, not typed).

## Limits, stated in advance
- The four families were chosen on gorilla; elsewhere the multi-copy class may be small or single-copy in fact.
- IsoCon sees only the reads of the panel's targets; the other tools see the whole library.
- Junction precision against annotation boundaries counts unannotated real junctions as wrong, for every tool alike.

## Result (2026-11-21)
Jobs 4491490-4491499, ten samples, every step exit 0 (IsoCon 128-471 s per sample; two pools failed in IsoCon: DAZL on
Sumatran fibroblast, DNAJB8 on chimpanzee LCL, which then have no IsoCon model). `r1872_isocon_panels/ID/{targets.tsv,
score.txt,junc.txt}`; gorilla rescored with today's Rustle in `artifact_refresh_2026-11-21/panel/`.
Genes of the four families in the panel: human 4 (no LRPAP1 by symbol), gorilla 49, chimpanzee 44, bonobo 33, Bornean
orangutan 4, Sumatran orangutan 5; 30 single-copy genes each.

Gene copies modelled / junction precision / spliced transcripts (Rustle; IsoCon), and junction precision with Rustle cut to
IsoCon's junction count on copies both model (Rustle / IsoCon):
| sample | Rustle | IsoCon | matched count |
|---|---|---|---|
| human A119b | 32 / 0.56 / 160 | 31 / 0.55 / 763 | 0.887 / 0.795 |
| human HG002 | 24 / 0.59 / 104 | 23 / 0.48 / 1,230 | 0.844 / 0.649 |
| gorilla OR6737 | 54 / 0.74 / 203 | 47 / 0.75 / 508 | 0.905 / 0.878 |
| gorilla KB3781 | 47 / 0.77 / 187 | 44 / 0.79 / 663 | 0.939 / 0.892 |
| chimpanzee testis | 39 / 0.77 / 135 | 27 / 0.88 / 204 | 0.928 / 0.901 |
| chimpanzee LCL | 37 / 0.80 / 129 | 30 / 0.74 / 840 | 0.881 / 0.794 |
| bonobo testis | 33 / 0.85 / 128 | 28 / 0.87 / 202 | 0.934 / 0.895 |
| bonobo fibroblast | 35 / 0.80 / 152 | 30 / 0.78 / 1,420 | 0.924 / 0.810 |
| Bornean testis | 34 / 0.81 / 150 | 31 / 0.80 / 411 | 0.904 / 0.835 |
| Bornean fibroblast | 30 / 0.87 / 101 | 28 / 0.81 / 958 | 0.936 / 0.851 |
| Sumatran testis | 34 / 0.81 / 117 | 32 / 0.83 / 260 | 0.900 / 0.848 |
| Sumatran fibroblast | 32 / 0.89 / 79 | 28 / 0.82 / 1,429 | 0.954 / 0.801 |
Across the twelve: IsoCon writes the most transcripts on 11, models fewer gene copies than Rustle on 12, has the higher raw
junction precision on 5; at IsoCon's junction count Rustle is the more precise on 12 of 12. The page section is built from
these tables by `artifact_refresh_2026-11-21/build.py` (statements computed, not typed).

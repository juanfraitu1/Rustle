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

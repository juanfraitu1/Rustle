# PREREG r1610: tie groups, frozen recipe, transfer to gorilla OR6737 (2026-11-10)

The recipe is frozen from the human A119b rows r1603-r1609 (user, 2026-11-10: freeze r1607 and register the gorilla row;
reusing the samples is accepted, the standard is transfer to another species). Nothing here is tuned on gorilla. Written
and committed before any gorilla step runs.

## Frozen recipe (r1607, stranded)
`tie_groups.py build PREFIX as98 0 1 RMSK_BED:0.5`: tied = AS >= 0.98 x the molecule's best; loci per (chrom, strand); link
= >= 2 shared tied molecules; a locus >= 50% interspersed repeat is set aside. Scorer `full` (full primates truth, same-family
pseudogenes) with the co-duplicated category from SEDEF.

## Gorilla inputs
- Reads `~/work/jfim/rustle_data/reads/GGO_mm.bam` (OR6737 testis); AS table `~/runs/gorilla_OR6737/or.molecules.tsv`.
- Genes: RefSeq `DENOVO_GGO_GFF` through `famtruth.load_genes('gorilla')` (a gene takes the human Compara identity of its
  symbol; LOC genes have none). Panel: `family_panel.py` with `PANEL_SPECIES=gorilla`, tools ours = `or.gtf` (the base
  run's assembly, not rebuilt), StringTie and FLAIR = `benchmark_collapse/GGO.*.gtf.gz`.
- RepeatMasker: UCSC GenArk `GCF_029281585.2.repeatMasker.out.gz`, converted to BED
  (`/anvme/workspace/v138da13-rustle/annot/dl_rmsk_ggo/GGO.repeatMasker.bed`, 5,462,308 rows).
- SEDEF: `/anvme/workspace/v138da13-rustle/sedef_dl/gorilla/final.bed`.

## Member set M (same definitions as human)
Panel members (T definition) with no StringTie model and no FLAIR model (r1601 rule), of which those with 0 decisive
molecules (r1602 rule, `hard_signal.py`). If M has fewer than 20 members the row STOPS as uninformative, with the count.

## Bars
- B1g grouped reach >= 0.80 of M (the human bar, 58 of 72).
- B2s purity, co-duplicated-only groups counted as pure, >= 0.90.
- B5r (replaces the fixed 60-locus limit, user-approved 2026-11-10 after r1607-r1609 showed the largest groups are
  single-chromosome arrays): the largest group lies on one chromosome, or >= 0.90 of its annotated loci are in one family.
- B4 the whole job (panel, records, build, score) <= 1 h wall.
- PASS = all four. No TAF11L bar (human-specific); the largest groups are named in the outcome.
- Descriptive: human figures beside gorilla (70/72 = 0.972, 0.982, strict 0.884), loci set aside as repeat, strict B2f.

## Outcome
(empty until the run)

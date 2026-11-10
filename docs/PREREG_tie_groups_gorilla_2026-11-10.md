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

## Outcome (2026-11-10, job 4456437, `/anvme/workspace/v138da13-rustle/r1610/`): STOP as registered, M = 0
- The gorilla panel is small: 75 judged families, 198 members listed, 176 expressed under T (human: 324 families, 1,136
  members). 12 members have no StringTie and no FLAIR model; 11 of them have >= 2 decisive molecules and 1 has one, so none
  is indistinguishable: **M = 0 < 20**. The bars are not judged.
- The panel has no member with MAPQ0 share >= 0.5 and one T-only member. Member hit rates: ours 0.858, StringTie 0.869,
  FLAIR 0.898.
- The recipe ran end to end in 341 s (panel 99 s, records 78 s): 70,008 molecules, 2,710 tied records, 151 loci, 80
  groups, 181 loci set aside as repeat. 23 judged groups: 17 pure, 6 impure only through co-duplicated loci (B2s 23/23,
  strict 0.739). Largest group 17 loci on 4 chromosomes, none of its loci in one truth family (B5r would fail).
- Cause of the small panel: a gorilla gene takes its family from the human gene of the same symbol, and the copies of
  multi-copy families in the RefSeq gorilla annotation are mostly `LOC` genes with no symbol, so they are outside the truth.
  A gorilla transfer needs a truth that names those copies (gorilla's own paralogue table, or sequence-based membership).

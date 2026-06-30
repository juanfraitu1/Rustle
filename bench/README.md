# bench/ docs index

The analysis record, consolidated into 10 thematic docs + this index (git history retains every original).
Each doc merges several former `*.md` verbatim under `## <name>` sections.

- **FAMILY_DEF.md** — the family-definition body: read-conflict/de-tie definition, formal proof,
  BAM/junction/VG signals, genome-wide PSV graph, criterion bake-offs, detection validation.
- **COPY_ASSIGNMENT_AND_GATE.md** — copy-split, the significance gate, PSV LLR-vs-votes, τ sweeps,
  phasing, identifiability bound, resolution boundary, quant, StringTie head-to-head, **the definitive
  O2 recompute on GGO_mm**.
- **THEORY.md** — the machine-checked theory (Lemmas/Theorems 1–7, MWCA facility-location capstone,
  F4 scope, identifiability limits).
- **VALIDATION_AND_STATUS.md** — the fully-simulated ground-truth benchmark, flagship case studies,
  paper-grounded reviews / airtight fixes / loose-ends audit, the defense-readiness audit, and the
  O1–O5 objective-status partition table.
- **FAMILY_LEVELS_AND_RELATED.md** — RNA/DNA/protein three-level cross-tab, the methods-to-find-
  multi-copy-families table, DNA/RNA overlay + Compara validation.
- **REFERENCE_ABSENT_AND_UNMAPPED.md** — reference-absent/collapsed copies, unmapped rescue,
  cross-chrom discovery, the reference-absent milestone, de-novo zero-annotation paralog copies.
- **DENOVO_PIPELINE.md** — de-novo family pipeline, two-pass, read-coherence, intron-chain discovery.
- **ASJ.md** — allele-specific junctions.
- **PERFORMANCE_AND_IO.md** — pipeline speedup, input formats/ties, winnowmap vs minimap2.
- **PANELS_AND_NOTES.md** — case panels + recent lever notes (conflict-edge unification, Bayesian
  posterior, intron-retention rescue, POA-core completion, loci/copies/isoforms, SEDEF build).

To recover any original: `git log --follow --diff-filter=D -- bench/<old-name>.md` then
`git show <rev>^:bench/<old-name>.md`.

## IGV visualization
- **`copy_assign --gtf`** — FLAIR-like assembly emit: writes `<out>.gtf`, an IGV-loadable transcriptome of all
  de-novo isoforms in the swept regions (intron-chain collapse → gene grouping, annotation-free). Multi-copy
  family copies appear as separate genes tagged `family_id`/`copy_index`/`multicopy "true"`; everything else is
  an ordinary isoform. Purely additive — the assignment outputs are byte-identical with or without it.
- **igv_tracks.py** — turns a `copy_assign` result into IGV tracks: `<out>.tagged.bam` (reads auto-coloured by
  assigned copy via the `YC` tag, grouped by `cp`, tied reads grey) + `<out>.copies.bed` (per-copy loci track).
  `python bench/igv_tracks.py --assignments OUT.assignments.tsv --bam reads.bam --regions regions.txt --out OUT`
- Together: load `<out>.gtf` (the transcripts, families flagged) + `<out>.tagged.bam` (reads coloured by copy) +
  the genome FASTA in IGV — the family copies as a transcript track above the copy-coloured reads.

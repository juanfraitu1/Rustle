# bench/ docs index

The analysis record, consolidated to **13 docs + this index** (down from 115; git history retains every
former file). Each canonical home has absorbed its former one-off notes and current-result anchors; recover
any removed file with `git log --follow --diff-filter=D -- bench/<old-name>.md` then
`git show <rev>^:bench/<old-name>.md`.

## Canonical topic homes

- **FAMILY_DEF.md** — **O1**, the RNA multi-copy family definition: read-conflict/de-tie vs E_r homology,
  γ-quasi-clique, formal proof, BAM/junction/VG signals, criterion bake-offs, and the shipped refine gates
  (reciprocal-identity divergence floor, multi-repeat-bridge gate, catalog artifact audit).
- **COPY_ASSIGNMENT_AND_GATE.md** — **O2**, copy-assignment (⚠ contested set is alignment-score
  near-ties, NOT MAPQ-0 — restated 2026-08-15, see `docs/copy_assignment_definition.md` §0) + the significance/de-tie
  gate (assign-or-abstain, never 1/k): copy-split, PSV LLR-vs-votes, τ sweeps, identifiability/resolution
  bounds, primary/secondary invariance, reference-free copy number, the SUN ladder, the recombinant-abstain
  gate, and the gene-conversion-vs-RT discriminator (biology caveated).
- **THEORY.md** — machine-checked theory: MCC = χ(H), NP-hardness, Strong Separation, the K-frontier,
  the facility-location capstone, identifiability limits, SUN, and EM = soft SDA PSV-clustering (consistency).
- **REFERENCE_ABSENT_AND_UNMAPPED.md** — **O4**, reference-absent/collapsed copies, unmapped rescue,
  cross-chrom discovery, de-novo zero-annotation paralog copies.
- **ASJ.md** — **O3**, allele-specific junctions (PSI / ΔPSI).
- **DENOVO_PIPELINE.md** — the de-novo family + copy-assign pipeline: two-pass, read-coherence,
  intron-chain discovery, the R4 readthrough filter, the DAZ2 locus_support fix, the containment coverage-floor.
- **FAMILY_LEVELS_AND_RELATED.md** — RNA/DNA/protein three-level cross-tab, the methods-to-find-multi-copy-
  families table, DNA/RNA overlay + Compara validation, and the **assembly-based parCN validation** (per-copy /
  diploid famCN on the mGorGor1 haplotypes).
- **VALIDATION_AND_STATUS.md** — the validation/reviews/objective-status home: the fully-simulated
  ground-truth benchmark **and the family-detection demo** (100%/100%), flagship case studies (incl. the TSPY
  honest-tie sim), the adversarial defense-readiness scorecard, cross-modal (Liftoff + SEDEF) confirmations,
  the GW false-positive audit + `--refine` default, the false-negative audit, **known-family
  sensitivity/precision** (incl. RFPL + EEF1A1/SRGAP2 controls), the **human cross-species** anti-overfit
  result, and the O1–O5 objective-status table.
- **PERFORMANCE_AND_IO.md** — pipeline speedup, SAM/BAM/CRAM I/O + ties, the alignment recipe, winnowmap vs minimap2.
- **PANELS_AND_NOTES.md** — case panels + lever notes, including the GSTM worked example.

## Definitions

- **DEFINITIONS_FORMAL.md** — the defense-grade formalization: the four-oracle homology lattice
  (E_a/E_b/E_r/E_p, χ(H)), the paralog / segmental-duplication / family / expansion / reference-absent
  predicates, **and the one-line glossary** of defense terms.
- **OBJECTIVES_FLOW.md** — the O1→O4 reads→VG→result walkthrough on real data.

## The benchmark

- **SOTO_BENCHMARK.md** — the headline real-data benchmark (state at 2026-07-15, commit 6fbc0e0): RNA copy
  recovery on the Soto human segmental-duplication catalog — the recovery panels, per-member
  sensitivity/precision (76.2% / 93%), the `--collapse-enumerate` measured effect + EEF1A1 control, the
  empirical near-identical keep-separate rules, and the "we do not over-call" specificity concordance (13/13).

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

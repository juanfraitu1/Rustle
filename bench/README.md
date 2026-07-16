# bench/ docs index

The analysis record, consolidated to **21 docs + this index** (down from 115; git history retains every
former file). Docs are grouped by role: the **canonical topic homes** (one per objective/theme, each having
absorbed its former one-off notes), the **definitions**, and the **current shipped-state anchors** (the
newest real-data results). Superseded experiments, probes, prototypes, and the retired-assembler-era docs
were removed — recover any with
`git log --follow --diff-filter=D -- bench/<old-name>.md` then `git show <rev>^:bench/<old-name>.md`.

## Canonical topic homes

- **FAMILY_DEF.md** — **O1**, the RNA multi-copy family definition: read-conflict/de-tie vs E_r homology,
  γ-quasi-clique, formal proof, BAM/junction/VG signals, criterion bake-offs, and the shipped refine gates
  (reciprocal-identity divergence floor, multi-repeat-bridge gate, catalog artifact audit).
- **COPY_ASSIGNMENT_AND_GATE.md** — **O2**, copy-assignment under MAPQ-0 ambiguity + the significance/de-tie
  gate (assign-or-abstain, never 1/k): copy-split, PSV LLR-vs-votes, τ sweeps, identifiability/resolution
  bounds, primary/secondary invariance, reference-free copy number, the genome-wide SUN ladder, the
  recombinant-abstain gate, and the gene-conversion-vs-RT discriminator (biology caveated).
- **THEORY.md** — machine-checked theory: MCC = χ(H), NP-hardness, Strong Separation, the K-frontier,
  the facility-location capstone, identifiability limits, SUN, and EM = soft SDA PSV-clustering (consistency).
- **REFERENCE_ABSENT_AND_UNMAPPED.md** — **O4**, reference-absent/collapsed copies, unmapped rescue,
  cross-chrom discovery, de-novo zero-annotation paralog copies.
- **ASJ.md** — **O3**, allele-specific junctions (PSI / ΔPSI).
- **DENOVO_PIPELINE.md** — the de-novo family + copy-assign pipeline: two-pass, read-coherence,
  intron-chain discovery, the R4 readthrough filter, the DAZ2 locus_support fix, the containment coverage-floor.
- **FAMILY_LEVELS_AND_RELATED.md** — RNA/DNA/protein three-level cross-tab, the methods-to-find-multi-copy-
  families table, DNA/RNA overlay + Compara validation.
- **VALIDATION_AND_STATUS.md** — the fully-simulated ground-truth benchmark, flagship case studies (incl.
  the TSPY honest-tie sim), the paper-grounded reviews, the adversarial defense-readiness scorecard, the
  cross-modal (Liftoff + SEDEF) confirmations, the GW false-positive audit + `--refine` default, the
  false-negative audit, and the O1–O5 objective-status table.
- **PERFORMANCE_AND_IO.md** — pipeline speedup, SAM/BAM/CRAM I/O + ties, the alignment recipe, winnowmap vs minimap2.
- **PANELS_AND_NOTES.md** — case panels + lever notes, including the GSTM worked example.

## Definitions

- **DEFINITIONS_FORMAL.md** — the defense-grade formalization: the four-oracle homology lattice
  (E_a/E_b/E_r/E_p, χ(H)), and the paralog / segmental-duplication / family / expansion / reference-absent
  predicates. Supersedes the four former per-level formal-definition files.
- **GLOSSARY.md** — one-line working definitions of the defense terms.
- **OBJECTIVES_FLOW.md** — the O1→O4 reads→VG→result walkthrough on real data.

## Current shipped-state anchors (state at 2026-07-15, commit 6fbc0e0)

- **PARCN_VALIDATION.md** — assembly-based per-copy / diploid famCN validation on the mGorGor1 haplotypes.
- **soto/SOTO_A119B_RECOVERY.md** — Soto human-family recovery in A119b (presence / de-novo / enumerate panels
  + the "we do not over-call" specificity concordance).
- **soto/SOTO_MEMBER_DETECTION.md** — Soto per-member sensitivity **and** precision (76.2% / 93%).
- **soto/NEAR_IDENTICAL_RULES.md** — the empirical keep-separate rules (%identity does not predict separability).
- **soto/COLLAPSE_ENUMERATE_MEASURE.md** — the `--collapse-enumerate` measured effect + EEF1A1 control.
- **SIM_DETECTION_DEMO.md** — 100% member detection / precision on a planted, non-circular ground-truth genome.
- **KNOWN_FAMILY_SENSITIVITY_PRECISION.md** — known-family sensitivity/precision (precision 1.00; the RFPL
  flagged-failure and EEF1A1/SRGAP2 negative controls).
- **HUMAN_CROSSSPECIES.md** — the identical binary/recipe on human testis Iso-Seq, tracking species-specific
  copy numbers (anti-overfit).

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

# Segmental-duplication extent / breakpoint calling (VG)

**Status (2026-06-03):** BUILT, opt-in `RUSTLE_VG_SEGDUP_EXTENT` (default OFF, analysis-only —
emits a trace, no GTF/assembly change). Validated on planted ground truth + real GOLGA8.

## Idea
A bare gene paralog shares homology only across the gene; a true SEGMENTAL DUPLICATION copied
the gene PLUS flanking sequence, so cross-copy homology extends into the intergenic flanks and
the duplication BREAKPOINTS are where that homology ends. IsoSeq reads are spliced mRNA and do
NOT cover intergenic flanks, so the signal is GENOME-based (the VG already loads the genome for
family discovery), anchored at the gene the family discovered.

## Method
- Pure `flank_homology_extent(a, b, window, min_identity)` (vg_hmm/segdup.rs): anchored sliding-
  window homology scan outward from the gene boundary; returns the breakpoint distance (where
  windowed identity drops below the floor). `call_segdup_extent` runs it on the upstream
  (reversed) + downstream flanks and classifies segdup vs bare paralog (flank total ≥ min_flank).
- `detect_and_report_segdups` (vg.rs): per family, for each copy pair, fetch the genome flanks
  (window default 8 kb) anchored at the copies' gene loci, take the max-extent pair, emit
  `[VG-SEGDUP] family=.. gene~Wbp up_flank=.. down_flank=.. duplicated_extent=.. => segdup|bare`.

## Validation
- 5 unit tests: identical flanks extend fully; homologous-prefix-then-random breaks near the
  boundary; immediately-divergent ≈ 0; segdup-vs-bare classification.
- Pipeline ground truth (gen_synthetic --segdup-flank): planted flank 600 → measured up/down
  **600/600 bp**; 400 → **399/401 bp** (window=30; breakpoints = gene boundary ± flank, EXACT).
- Negative: bare paralog (flank 0) → up/down ≈ window/2 noise, is_segdup=false.
- Real GOLGA8: **5 segmental duplications** (families 3/6/8/17 — flank homology hits the 8 kb
  fetch cap, i.e. multi-kb shared flanks) vs **4 bare paralogs** (families 0/2/11/13, gene-only).
  The main 7-copy GOLGA8 family reads as bare paralog (no flanking homology beyond the gene).

## Limitations / notes
- Genome-based (not read-based) — fundamental: mRNA doesn't transcribe the flanks.
- The anchored scan assumes the copies' gene boundaries correspond; large indels in the flanks
  would drift the offset (a banded re-anchor would harden it).
- A synthetic with IDENTICAL shared flanks ≥ ~700 bp trips the genome-based discovery quality
  filter (extra cross-copy k-mers → family dropped); a real segdup's flanks are at gene-level
  identity, not identical, so this is a synthetic-only artifact (real GOLGA8 segdups detect fine).

## UPDATE: hardened + surfaced to GTF (2026-06-03, commit f004e2f)
From an adversarial perf/hardening audit (4 auditors → filtered; ~10 premature/churn items rejected):
- HONESTY: flank_homology_extent returns 0 (not phantom window/2) when the anchor window itself
  fails — bare paralogs now read up=0/down=0 (was 25/25). The reported bp ARE the deliverable.
- DEFINITION: is_segdup now requires homology on BOTH flanks (up.min(down) ≥ min_each_flank,
  default 50; RUSTLE_VG_SEGDUP_MIN_EACH) — a true segdup has two breakpoints. Real GOLGA8
  reclassifies one one-sided family (up=0/down=520) segdup→bare: now **4 segdup / 5 bare**.
- GTF: the call is surfaced as family_segdup / family_segdup_extent / family_flank_up /
  family_flank_down attrs (opt-in; default-off byte-identical).
- DEFERRED (audit agreed): banded re-anchor for flank indels is a FEATURE, not a hardening
  fix — the offset-lock limitation stays documented; revisit if needed for indel-rich segdups.

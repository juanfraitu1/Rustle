# Local variation graph for tandem mega-bundles (intra-bundle copy resolution)

**Date:** 2026-06-02
**Status:** Design DRAFT — scoping option 3 (per advisor discussion). Not yet approved for implementation.
**Scope:** VG / `--vg`-gated only. Default de-novo (95.6/90.5) must stay byte-identical.
**Depends on:** [`2026-06-01-flow-capacity-apportionment-design.md`](2026-06-01-flow-capacity-apportionment-design.md) (EM apportionment supplies abundance/abstention). Supersedes the structure-borrowing premise of [`2026-06-02-family-bundle-o5-design.md`](2026-06-02-family-bundle-o5-design.md) (corrected: borrowing already exists; the real gap is upstream — see below).

## Problem (established empirically 2026-06-02)

VG copy-resolution works for **dispersed** paralogs and is **structurally inert for tandem arrays**. Root cause, traced end-to-end on RBMY1 (6-copy testis Y array, `NC_073248.2:19.56–19.74 Mb`):

- The 6 copies are 93–99.8% identical, ~23 kb apart. Long reads' **secondary** alignments (512 of them, RBMY is near-identical) scatter each read across **all 6 copies**, bridging the 9 kb intergenic gaps — even though primary coverage is 86.5% zero across the span. → the whole array collapses into **one bundle locus** (`19,558,214–19,739,408`, ×2 as a strand-split).
- VG family discovery is **inter-bundle**: it links ≥2 *separate*, non-overlapping, non-strand-mirror bundles by shared multimappers (`discover_family_groups`, vg.rs:264). With one locus there is nothing to link — supplementary scan = **1** cross-bundle link (488/512 secondaries land in the same bundle as their primary), sequence-similarity = the 2 same-span bundles are overlap-skipped.
- Result: **0 families discovered → EM / completion / borrow never engage.** Ordinary within-bundle assembly still carves the copies into genes (RSTL.1–5 by position); the read-starved outlier c6 (2 reads) never clears the assembly threshold → 0 transcript.

Contrast (the proof): **DAZ** (2 dispersed inverted loci → 2 bundles) gives **179** cross-bundle links → family discovered → VG engages. The difference is *dispersed vs tandem*, nothing else.

> **Honesty correction (must land regardless of this design):** the earlier RBMY "validation" (5/6 copies, `capacity_confidence 1.000`) is *ordinary assembly*, not VG copy-resolution — `cc 1.000` is the no-multimapper default; no RBMY family was ever discovered. The `rbmy_analysis`/`rbmy_validation` figures must be corrected. The DAZ result stands.

## Key insight — most of this is already built and tandem-aware

`build_family_graph` (family_graph.rs:411) takes a `FamilyGroup` whose `bundle_indices` are **one-bundle-per-copy**, and its own comments target "RBMY / TSPY scattered across a chromosome" with a sequence-similarity fallback (`merge_singletons_by_sequence`, family_graph.rs:~460) and `cluster_by_position` (family_graph.rs:183). The downstream is all present:

- `run_fingerprint_em(families, bundles, family_graphs, max_iter)` (vg.rs) — per-read apportionment across copies, with the anchored-prior + score-gap gate (flow-capacity spec).
- `forward_against_path_for_copy(fg, read, path, copy_id)` (scorer.rs:389) — per-copy forward score of a read against one paralog's path. **This is the read-threading + abstention primitive.**
- `build_bundle_completion_nodes` / `build_bundle_borrow_junctions` / `build_bundle_borrow_coverage` (vg.rs) — dark-exon / junction / coverage borrowing, already default-on.
- `FamilyVerdict` carries a `LocusRel::Tandem` classification (vg.rs:1719) — the vocabulary exists.

**The machinery is fed inter-bundle families and never receives a tandem one.** So the genuinely-new work is narrow: **decompose a single tandem mega-bundle into a synthetic per-copy `FamilyGroup`**, then hand it to the existing pipeline. This is the local variation graph: `build_family_graph` over the decomposed copies *is* the graph — shared ExonClasses across copies are the tandem cycle; copy-specific exons/SNVs are the bubbles.

## What a "true variation graph" buys here (and its ceiling)

Representing the mega-bundle as a graph (shared backbone + copy-specific bubbles + a tandem cycle) gives three things a linear bundle can't:
1. **Decompose, not collapse** — copies become distinct graph paths even though they share a bundle.
2. **Copy number from traversal depth** — depth through the tandem cycle estimates how many copies are present.
3. **Honest read threading** — a read covering a copy-specific bubble is assigned; a read covering only shared backbone is *routed to the backbone and abstains*.

**Ceiling (non-negotiable):** the graph **represents** ambiguity; it does not **manufacture** identity. RBMY copies (real copy-specific positions) → genuine recovery. DAZ-tied reads (99.97% identical, intronic-only diagnostics) → abstain, same conclusion as the linear analysis. This design must not imply otherwise, and must never fabricate a silent copy (the DAZ3 rule).

**Out of scope here:** *inversions* (the dispersed-inverted DAZ case — handled by strand logic, not tandem decomposition) and *TE-exon* bubble detection (a separate novel-isoform capability). This design targets **tandem expansions / mega-bundle polishing** only.

## Design — intra-bundle tandem decomposition → existing FamilyGraph/EM

Five components; only #1–#2 are new code, #3–#5 are wiring of existing functions.

1. **Tandem mega-bundle detector** `vg::detect_tandem_bundle(bundle, genome) -> Option<TandemDecomposition>`.
   A bundle is a tandem-array candidate when:
   - its reads/assembled exons form **≥2 genomically-separated position clusters** (reuse `cluster_by_position` on the bundle's own exons), AND
   - the clusters are **mutually sequence-similar** above a threshold (reuse the min-hash/Jaccard from `discover_sequence_similar_bundles`, applied *within* the bundle to the per-cluster reference sequence), AND
   - separation ≥ a min gap (avoid splitting a single multi-exon gene).
   Returns the per-cluster genomic spans (the copy loci) + a similarity matrix. Gate: `RUSTLE_VG_TANDEM` (default OFF).

2. **Synthetic per-copy FamilyGroup builder** `vg::decompose_tandem_to_family(bundle, decomp, &mut bundles) -> (FamilyGroup, Vec<usize>)`.
   For each detected copy cluster, emit a **synthetic sub-bundle** containing that cluster's reads (mirror `rescue::synthesize_bundles_*`, which already build `Bundle`s from read/exon sets). Append to `bundles`; collect their indices as `bundle_indices` of a new `FamilyGroup`. Populate `multimap_reads` from reads whose alignments span >1 cluster (the secondaries that caused the collapse — now they're *evidence*, not a collapse force).

3. **Local variation graph** = `build_family_graph(synthetic_family, bundles, genome, …)` — **unchanged call**. Position-clustering + the singleton sequence-merge produce ExonClasses spanning the copies (the tandem cycle); copy-specific exons stay `copy_specific=true` (the bubbles). `annotate_per_copy_exon_coverage` fills `per_copy_cov`.

4. **Read threading + abstention** = `run_fingerprint_em([synthetic_family], bundles, [Some(fg)], max_iter)` — **unchanged call**. Per-read apportionment writes `read.weight` per copy; `em_anchored`/`em_weight_gap` (flow-capacity) mark decisive vs abstained reads. A read whose best-vs-next per-copy `forward_against_path_for_copy` gap is below the score-gap gate stays at shared weight → low `capacity_confidence` downstream (honest abstention). No new scoring code.

5. **Assembly + copy number + labelling.** The synthetic sub-bundles assemble through the normal path (they're in `bundles`); completion/borrow now fire (a family exists), so a starved copy (c6) gets sibling structure at its own apportioned abundance. Copy-number = number of populated clusters; emit as `family_n_copies` (already in `FamilyVerdict`). Transcripts carry `vg_family_id`/`vg_copy_id`/`capacity_confidence`/`abstain` as for dispersed families; add `tandem_copy "true"` so the provenance is explicit.

## Algorithm

```
for each bundle B in --vg mode (RUSTLE_VG_TANDEM on), with genome:
  decomp = detect_tandem_bundle(B, genome)            # NEW (component 1)
  if decomp is None: continue                          # not a tandem array
  (fam, sub_bis) = decompose_tandem_to_family(B, decomp, &mut bundles)  # NEW (2)
  fg = build_family_graph(fam, bundles, genome, …)     # EXISTING (3)
  annotate_per_copy_exon_coverage(&mut fg, …)          # EXISTING
  run_fingerprint_em([fam], &mut bundles, [Some(fg)], …)  # EXISTING (4)
  # sub-bundles now assemble via the normal loop with completion/borrow active (5)
  # tag resulting transcripts: vg_family_id, vg_copy_id, capacity_confidence,
  #   tandem_copy=true; copies with no decisive reads -> low_confidence (abstain)
```

The original mega-bundle B is **replaced** by its synthetic sub-bundles for assembly (so we don't double-emit). De-novo path untouched (entire block gated by `vg_mode && RUSTLE_VG_TANDEM`).

## c6 — the concrete win

c6's 2 reads form a position cluster → a synthetic copy in the family → the family graph borrows the conserved structure from c1–c5 → c6 assembles full-length at its **own** apportioned abundance (~2), tagged `tandem_copy`/`capacity_confidence` low. RBMY copies are distinguishable (the outlier c6 at 93–96% is *more* so), so threading genuinely assigns its reads — this is a real recovery, not a DAZ-style tie.

## Honesty contract

- **No fabricated copies.** A cluster needs ≥1 own read to become a copy. Depth-based copy-number is reported as an *estimate/range* with confidence — never materialized as a silent transcript.
- **Abstain, don't guess.** Reads with no decisive per-copy variation stay shared-weight → low `capacity_confidence` → `low_confidence` tag (kept, not dropped).
- **Inferred vs observed.** Borrowed structure on a starved copy is labelled (`tandem_copy` + low `capacity_confidence`); expression is the copy's own apportioned reads.
- **DAZ unchanged.** DAZ is dispersed/inverted, not a tandem mega-bundle → `detect_tandem_bundle` returns None → no behaviour change.

## Gating

- `RUSTLE_VG_TANDEM` (default **OFF**), `RUSTLE_VG_TANDEM_MIN_SEP` (default e.g. 5 kb), `RUSTLE_VG_TANDEM_MIN_JACCARD` (default 0.20, reuse).
- Entirely inside `config.vg_mode`; new transcript fields written only in VG mode, emitted only when `Some`. De-novo byte-identical (hard gate).
- Requires `--genome-fasta` (needed for per-copy sequence + similarity). If absent, detector returns None.

## Test & validation plan

1. **Regression (first):** default de-novo on GGO_19, no `--vg` → 95.6/90.5, GTF byte-identical to HEAD.
2. **RBMY primary success:** `--vg --genome-fasta <chrY> RUSTLE_VG_TANDEM=1` → the array is discovered as a **1 family / 6-ish copies**; c6 emits a full-length transcript at abundance ~2 with `tandem_copy`/low `capacity_confidence`; c1–c5 gain `vg_family_id`/`vg_copy_id` attribution (vs the current bare assembly).
3. **DAZ unchanged:** `detect_tandem_bundle` returns None on DAZ; DAZ GTF identical to the non-tandem `--vg` run.
4. **No-fabrication:** a synthetic tandem with a truly silent copy (0 reads) yields no transcript for it; abstaining reads get `low_confidence`.
5. **False-split guard:** a single large multi-exon gene with internal repeats is **not** split into copies (separation + cross-cluster-similarity thresholds); add a targeted fixture.
6. **Genome-wide:** re-run `bench/paralog_secondary_scan` with the flag; assert no spurious tandem families on non-array loci and no precision regression.

## Risks

1. **False decomposition** (splitting one gene into "copies") — guard with min-separation + high cross-cluster sequence similarity + min per-cluster support; test #5.
2. **Performance** — O(n²) similarity within a bundle; cap cluster count and total exons (reuse the existing 2000-exon guard in `build_family_graph`).
3. **Copy-number over/under-call** — report as estimate+confidence; never fabricate; depth-based refinement is optional and clearly labelled.
4. **Double-emission** — the original mega-bundle must be replaced (not assembled alongside its sub-bundles); assert no overlapping duplicate transcripts in test #2.
5. **Identifiability over-claim** — abstention path (#4) is mandatory; do not assign reads lacking copy-specific signal.

## Effort & staging

- **Phase A (detector + decomposition, components 1–2):** the only new algorithms; unit-testable on synthetic tandem fixtures without the full pipeline. ~the bulk of the work.
- **Phase B (wiring 3–5):** call existing `build_family_graph` / `run_fingerprint_em` on the synthetic family; tag transcripts. Mostly plumbing.
- **Phase C (validation):** RBMY + DAZ + de-novo + false-split + genome-wide.

A separate TDD implementation plan follows once this design is approved. In parallel and independent of approval: **correct the RBMY figures** (option 2) so the validation narrative is honest immediately.

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

## Validation results (2026-06-03, prototype on RBMY1)

Implemented (branch `vg/flow-capacity-apportionment`, opt-in `RUSTLE_VG_TANDEM=1`). Two integration-driven corrections to the T6/T7 design were required and made:
1. **detect** clusters PRIMARY-only footprints — secondary/chimeric alignments of near-identical copies bridge the copies (RBMY footprints up to 56 kb), collapsing gap-clustering into one. Primaries (one best placement each) resolve the 6 copy loci cleanly.
2. **decompose** assigns ALL alignments by start-in-span (not a fixed read→cluster map), so a read's primary+secondary across copies becomes a multimapper the EM apportions — borrowed structure, honest abundance.

**Outcome on RBMY1** (`--vg --genome-fasta <chrY> RUSTLE_VG_TANDEM=1`):
- `[VG-TANDEM] bundle 0 (NC_073248.2:19,558,214-19,739,408) → 6 copy sub-bundles`; discovered as **1 family / 6 copies**, `family_locus_rel "tandem"`. The full FamilyGraph + fingerprint-EM engaged (4354 diag-sites, 110 reads apportioned, decisive=29/moderate=43/uncertain=38) — none of this ran before (0 families on the collapsed mega-bundle).
- **c6 (LOC129530242) recovered**: 0 transcripts → 1 transcript at 19,720,060-19,730,933, cov 1.7, **capacity_confidence 0.956**, `tandem_copy "true"`, full `family_id/copy_id/copy_independent_support` attribution. This is the headline result.
- `capacity_confidence` now *varies per copy* (0.08–0.96) — a real EM signal (vs the pre-feature default 1.000 = no apportionment).

**Gates passed:**
- **De-novo unaffected**: no `--vg` → identical 16 transcripts, zero tandem/family attrs (all changes are `config.vg_mode`-gated; the `tandem_copy` field defaults None and emits only when `Some(true)`).
- **DAZ unchanged**: `RUSTLE_VG_TANDEM=1` on DAZ fires **0** `VG-TANDEM` (detect returns None — DAZ is dispersed/inverted, 1 cluster per bundle); DAZ still discovers its 1 dispersed family via the normal path.

**Known limitations (honest, for follow-up — NOT blockers):**
- **c4 (LOC129530238) absorbed into c5.** Baseline: c4=1 tx, c5=6 tx. Tandem: c4=0 tx, c5=6→9 tx. c4's reads (99.8% identical to the dominant 45-read c5, with no decisive copy-specific sites) are apportioned by the EM to c5. This is *arguably correct* (c4 is non-identifiable from c5) but it loses a copy baseline kept, and inflates c5's isoform set (over-enumeration). The fix is in the **EM apportionment** (flow-capacity territory, not the tandem split): preserve a read's PRIMARY-placement weight as anchoring evidence for its own copy, or abstain-don't-merge (emit c4 at low cc rather than vanish). The tandem decomposition *exposes* this EM behavior; it does not cause it.
- **Mild over-enumeration** (20 vs 16 transcripts) from borrowed sibling structure — same root as c4 absorption; the apportioned-coverage gate should filter copy-absent isoforms more aggressively.

**Remaining validation (not yet run):** full GGO_19 de-novo byte-identical regression (the hard gate — structurally guaranteed but unmeasured here); false-split guard on a single multi-exon gene; genome-wide spillover scan. The feature stays opt-in (`RUSTLE_VG_TANDEM` default OFF) pending these + the c4/over-enum EM refinement.

## c4 fix + regression gates (2026-06-03)

**Regression gates established (the safety net for the EM-adjacent change):**
- DAZ3-suppression baseline: DAZ resolves as a family with DAZ3 at cov 4.0/1.4 (the validated flow-capacity residual; a DAZ-style phantom copy shows independent_support ~0.6).
- De-novo reference: GGO_19.bam → 1943 transcripts, 0 tandem/family attrs. (Byte-level diff is unusable on the large bam — rustle's parallel assembly emits transcripts in nondeterministic ORDER; the same binary re-run differs by ~40k lines but is sorted-identical. The correct gate is **sorted-identical content**.)
- EM regression suite: vg_hmm_em 6/6, vg_hmm_negative_control 3/3, vg_hmm_diagnostic 5/5 green. `vg_fingerprint_em_correctness` has 1 PRE-EXISTING failure (`fp_em_decisiveness_scales_with_diagnostic_site_coverage`, fixture drift on `test_data/synthetic_family` — never touched by the tandem commits; unrelated).

**The c4 problem (grounded):** RBMY1 c4 is genuinely identifiable (NM=35 / de=0.0022 from the dominant c5) and certified (independent_support **0.895**) — NOT a DAZ3-style phantom. Yet it had `this_cov` ~100 (raw) but **0 transcripts**: the fingerprint-EM apportioned its own primary reads to c5 (at the diagnostic sites c4's reads fit c5 marginally better, 0.0021 vs 0.0054), starving c4 of assemblable weight. The support metric cleanly separates the cases — c4 = 0.895, DAZ1 (real) = 0.989, DAZ phantom = 0.596 — so a support gate at 0.75 includes c4 and excludes phantoms.

**Fix:** after `compute_copy_independent_support` populates `vg_copy_support` (pipeline.rs:10894, pre-assembly, weights still post-EM), a **support-gated weight floor**: for each `RescueClass::TandemCopy` sub-bundle whose `independent_support ≥ RUSTLE_VG_TANDEM_SUPPORT_FLOOR` (0.75), raise its OWN primary reads' weight to ≥ `RUSTLE_VG_TANDEM_WEIGHT_FLOOR` (0.5) where it was under that. Only RAISES under-weighted primaries (already-emitting copies untouched). Scoped to TandemCopy sub-bundles → DAZ/dispersed and de-novo are on different code paths.

**Validated (all gates pass):**
- RBMY1 c4 **recovered**: 0 → 2 transcripts (cov 1.1/1.8) at **capacity_confidence 0.220** — honest abstain-in-place (low cc reflects the contested evidence), not a silent merge or a confident claim. All **6 copies** now emit; c6 retained.
- DAZ GTF **byte-identical** to baseline (DAZ3 stays cov 4.0/1.4; `[VG-TANDEM] weight floor` never fires — DAZ has no TandemCopy bundles). **No DAZ-regression possible** by construction.
- GGO_19 de-novo **content-identical** (1943 transcripts, sorted-identical).

Net RBMY1 result with `RUSTLE_VG_TANDEM=1`: the collapsed mega-bundle is resolved into all 6 tandem copies with per-copy attribution; the outlier c6 (was unassemblable) and the core c4 (was EM-absorbed) both recovered, each at a capacity_confidence honestly reflecting its evidence (c6 high 0.956, c4 low 0.220). Feature remains opt-in pending false-split/spillover guards + the pre-existing fixture-test fix.

## False-split / spillover scan + promotion decision (2026-06-03)

**Genome-region spillover scan** — `--vg` baseline vs `--vg RUSTLE_VG_TANDEM=1` on `bench/multi_copy_eval/ggo_Y.bam` (11,333 reads, all on NC_073248.2 / chrY, spanning many multi-copy regions: TSPY ~9 Mb, RBMY ~19 Mb, DAZ ~42 Mb, …) with the chrY genome:
- **316 transcripts in BOTH runs — identical count, zero regression.**
- **0 `[VG-TANDEM]` detections, 0 false splits, weight floor never fired** across all 28 discovered families. So tandem mode is **SAFE genome-wide**: it does not spuriously split single genes or perturb non-tandem loci. The false-split guard (min-gap + min-jaccard similarity filter + the `detect_dissimilar_clusters_returns_none` unit test) holds empirically.

**Sensitivity gap (the decisive finding for promotion).** In the *full-chrY* context, tandem did **not** fire on RBMY — yet RBMY is **not** correctly resolved by normal discovery either: all 6 RBMY genes come out tagged `copy_id "3"`, `capacity_confidence 1.000` (the no-apportionment default), swept into an unrelated cross-mapping family (`FAM_16`). Because the pre-pass **skips bundles already in a discovered family** (`discovered.contains(&bi)`), and RBMY's bundles were mis-assigned to FAM_16, the tandem pass never runs where it is needed. So the feature currently fires only when an array collapses to a bundle that escaped (incorrect) family assignment — i.e. on **isolated / region bams** (like the RBMY-only `tspy.bam`), not on whole-chromosome input.

**DECISION: do NOT promote to default-on. Keep `RUSTLE_VG_TANDEM` opt-in.**
Rationale: the promotion criterion "on a real multi-copy benchmark, adds recovered real copies without adding FPs" is **not met** — on whole-chromosome input the feature doesn't fire (sensitivity gap), and on isolated bams it works but that is not the production path. Promoting now would ship a default-on feature that is inert (or worse, interacts with mis-familying) on real input. It is **safe** (spillover-clean), so shipping it opt-in is correct: it is a validated capability for isolated-array analysis and a research result (per-copy confidence on a tandem testis array), not yet a whole-genome default.

**Path to promotion (future work):** make the tandem pre-pass fire on bundles whose family assignment is a *cross-mapping artifact* (e.g. low within-family junction/sequence concordance, or a bundle that itself spans ≥2 sequence-similar position clusters) rather than blanket-skipping any discovered-family bundle — so RBMY-in-full-genome (mis-familied into FAM_16) is still decomposed. Re-run this spillover scan + a real multi-copy benchmark after that fix; promote only if it then adds real copies without FPs.

## Sensitivity gap CLOSED (2026-06-03, commit 8554608)

The diagnostic (`RUSTLE_VG_TANDEM_TRACE`) settled the cause: in full-chrY input RBMY is a single mega-bundle (`bundle 179`, 19,558,214–19,739,408, 530 reads) that `detect_tandem_bundle` correctly fires on (**6 copy spans**) but the pre-pass **skipped** because `in_discovered_family=true`. All **4** chrY ampliconic arrays were detected-but-skipped (RBMY=6, TSPY≈34.7 Mb=5+5, 25.7 Mb=2) — none decomposed.

**Fix:** run detection on EVERY non-empty bundle (drop the blanket discovered-family skip). `detect_tandem_bundle` returns `None` for a genuine dispersed-family member (single position cluster), so this only adds REAL collapsed arrays; dispersed families (DAZ) are untouched. A mis-claimed bundle is neutralized (NOT removed from its family — that would shift the family's `multimap_reads` position indices), so its stale family slot emits nothing and the array's sub-bundles form their own tandem family.

**Validated (chrY `ggo_Y.bam`, all gates pass):**
- All 4 ampliconic arrays now decompose; RBMY resolves with **real per-copy capacity_confidence** (0.068–0.435; c4 region = 0.220), 41 tandem-copy transcripts (19 honestly flagged cc<0.5), vs the prior all-`copy_id 3`/`cc 1.000` default.
- **No false splits** — exactly the 4 known ampliconic arrays fired, nothing else.
- Change is **array-localized**: RBMY 54→29 transcripts (mega-bundle jumble → cleaner per-copy), TSPY count unchanged. Non-array regions differ only within `--vg`'s own run-to-run nondeterminism (two baseline `--vg` runs already differ; count stable 227=227) — no real side-effect.
- DAZ GTF **byte-identical** (dispersed, never fires); GGO_19 de-novo **sorted-identical**; isolated RBMY (`tspy.bam`) still recovers c4 (cc 0.220) + c6.

**Promotion status (updated):** the gap is closed — the feature now fires on whole-chromosome input and resolves all real tandem arrays, safely. It is a **strong promotion candidate**. Remaining before flipping default-on: a **precision/recall benchmark against annotation** to confirm the per-copy resolution (e.g. RBMY 54→29) is more *accurate*, not merely more granular (no ground-truth scoring was available in this environment). Until then it stays opt-in, but now genome-wide-functional rather than isolated-bam-only.

## Benchmark harness (2026-06-03)

`bench/tandem_benchmark.sh <reads.bam> <genome.fa> [annotation.gtf] [region.bed] [out_dir]`
runs `--vg` ± `RUSTLE_VG_TANDEM=1` and reports (a) structural metrics always
(transcript counts, arrays decomposed, per-copy cc distribution, per-array deltas)
and (b) transcript-level Sn/Pr vs a supplied annotation (gffcompare `-R -Q`,
optionally region-restricted to the arrays). This is the harness for the remaining
promotion gate. `bench/chrY_tandem_arrays.bed` lists the chrY arrays (RBMY/TSPY/
25.7 Mb). Validated: structural mode on `ggo_Y.bam` (4 arrays, RBMY 54→29);
accuracy plumbing confirmed (baseline-vs-self Sn/Pr 100.0). **To run the gate:**
supply an NCBI/RefSeq GTF covering NC_073248.2 as arg 3. **For per-copy
ATTRIBUTION accuracy** (did each read go to the right copy — the sharper thesis
metric), a synthetic with per-read source-copy truth is needed; score the EM's
max-weight placement (via `RUSTLE_VG_FP_ATTR_TSV`) against truth — a complementary
harness, not yet built.

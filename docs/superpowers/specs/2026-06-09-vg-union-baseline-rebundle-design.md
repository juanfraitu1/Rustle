# VG over-collapse recovery — clean re-bundle + primary-gated union

**Date:** 2026-06-09
**Status:** Design (not yet implemented)
**Context:** `bench/flow_recall_phase0/VG_REGRESSION_FINDING.md` (root cause + validated direction)

## Problem

VG mode (`--vg`) is **net-neutral vs baseline** genome-wide (baseline 24,374 FSM vs VG 24,373):
it recovers 108 annotated isoforms baseline misses (real copy/multimapper wins) but **drops 109**
that plain baseline `-L` assembles fine. The `VG⊇baseline` invariant does not hold.

**Root cause (confirmed):** at secondary-bearing family bundles, cross-mapped secondary reads
pollute the splice graph, and the assembly drops PRIMARY-supported isoforms that the primary-only
(baseline ≡ StringTie) path keeps. A mega-bundle variant: secondaries bridge a gene to its
neighbour into an unsplit read-through. The dropped isoform and the polluting secondaries are
structurally indistinguishable (no junction-support strip separates them). Documented at
`pipeline.rs:12127-12176`. My attribution confirms it: 94% of regressions are `context_only`
(recover when the family doesn't form on a slice); removing the decisive gate recovers only 1/17.

**Why prior fix failed:** the existing `RUSTLE_VG_UNION_BASELINE` produces baseline transcripts by
**cloning sub-bundles** (`b.clone()`, pipeline.rs:12198) that inherit the mega-bundle's stale
cached graph structures and junction-pair stats. The clone-and-patch is incomplete → mis-assembly
→ recovers only ~6/13, net-harmful on mega-bundles.

## Goal & validated direction

Recover the 109 over-collapse regressions so VG ⊇ baseline, taking VG from net −1 to **+108** over
baseline, **net-F1-positive**. Measured ceilings (`bench/flow_recall_phase0/`):
- Per-chrom (NC_073247.2): post-process union = 1043 FSM, recovers 17/17 regressions, loses no VG win.
- Genome-wide: VG 24,373 → 24,482 (+108 over baseline, +109 over VG).
- **Targeted** union (over-collapse zones) recall:FP = **1.27** (net-positive — first extraction
  lever >1.0); naive genome-wide union = 0.10 (net-negative — why UNION_BASELINE was abandoned).
  Targeting + novelty are load-bearing.

## Approach: clean re-bundle + primary-gated union (internal, opt-in)

A targeted, strictly-additive recovery pass inside the VG run, reusing the existing union stage
and replacing only the broken transcript-production step.

```
VG run (unchanged) ──► all_transcripts (VG output)
  ├─ identify secondary-bearing bundles  (existing: bundle has ≥1 secondary read)
  ├─ NEW: re-bundle their PRIMARY reads via the normal splitter into FRESH bundles,
  │       run the normal assembly closure → clean baseline transcripts (holdout)
  └─ end-stage union (existing, pipeline.rs:19337): add a holdout tx iff
        · intron chain ABSENT from all_transcripts  (novel — never displaces a VG tx)
        · multi-exon (single-exon skipped)
        · NEW: primary-support gate — longcov ≥ N
```

Three load-bearing properties: **targeting** (only secondary-bearing bundles feed the holdout),
**novelty** (only chains VG dropped), **primary-support gate** (only real primary-backed dropped
isoforms — the over-collapse signature — not baseline's thin alt-isoform collateral).

### Component 1 — clean re-bundle (the core fix)

Replace clone-and-patch with **build-fresh**. For each secondary-bearing bundle:
- take its **primary reads only**; split by the runoff gap (`split_spans_by_runoff`) — un-bridges
  the mega-bundles the secondaries papered over;
- for each read group, construct a **fresh `Bundle` from defaults** (NOT `b.clone()`) carrying
  only those reads, all cached state empty (`bundlenodes`/`read_bnodes`/`bnode_colors`=None,
  `synthetic`=false, `vg_family_id`=None) so the assembly rebuilds everything from the reads;
- run through the **same per-bundle assembly** the normal `-L` path uses; tag outputs
  `UnionBaseline`, route to the existing holdout vector.

**Principle:** a transcript's correctness must not depend on state carried from the polluted
parent bundle. Build-fresh guarantees this; clone-and-patch could not. Implements the
"re-bundle the primary reads (run the splitter)" fix the code comment names as correct.

**Open implementation detail (first plan task):** confirm the per-bundle assembly is cleanly
callable on an arbitrary `Bundle`. If not trivially reusable, factor it into a function first.

### Component 2 — primary-gated union (extend the existing stage)

The end-stage union (pipeline.rs:19337) already: strictly additive, dedup by exact intron chain,
skip single-exon. **Add one condition:** union a holdout tx only if `longcov ≥ N`
(`RUSTLE_VG_UNION_BASELINE_MIN_LONGCOV`, **default 2** as starting point, confirmed/adjusted in
validation, env-overridable). Keep the existing `[VG-UNION-BASELINE] unioned N …` log line.

## Env gating

- Reuse `RUSTLE_VG_UNION_BASELINE` (this is its corrected implementation; build-fresh + gate
  replace the broken clone path). Default-OFF → VG output **byte-identical** when unset.
- `RUSTLE_VG_UNION_BASELINE_MIN_LONGCOV` — the gate threshold (default 2; confirmed/adjusted in validation).

## Validation (mirrors the strand-bundling discipline)

1. Per-locus confirm: NC_073224.2 `XM_063708549.1` — full `--vg` + flag recovers the FSM chain;
   without flag, byte-identical to current default.
2. Genome-wide ±flag across 26 contigs: measure recovered regressions, FP cost, net recall:FP.
3. Promote to default only if net-positive genome-wide (separate decision, like the strand flip).

## Testing

- **Unit (TDD):**
  - re-bundle helper: synthetic secondary-bearing bundle (primary reads spanning two genes bridged
    by a secondary) → build-fresh re-splits into two clean bundles whose assembly yields both
    genes' chains (clone path would yield a read-through).
  - gate: a below-threshold holdout tx is NOT unioned.
  - novelty: a holdout chain already in VG is NOT duplicated.
- **Integration:** NC_073224.2 `XM_063708549.1` — flag recovers the FSM chain; default byte-identical.
- **Regression guard:** default-off run byte-identical (no behavioral change unless opted in).

## Out of scope

- Promoting to default (separate post-validation decision).
- The deeper in-graph fix (stop secondaries polluting the graph in the first place) — the
  prior in-flow attempts at that failed; the re-bundle/union is the validated path.
- Single-exon recovery (FP-prone; gated separately by `--read-chain-single`).

## Risk

- Re-running assembly on secondary-bearing bundles' primary reads adds compute (bounded — only
  secondary-bearing bundles, and only their primary reads). Per-locus `--vg` was 2.4s/720MB; the
  extra pass is a fraction of the run.
- Build-fresh assembly must reproduce baseline output for the sub-bundle — the whole point; the
  unit test (two-gene bridge) guards the mega-bundle case that broke clone-and-patch.
- Strictly additive + novelty dedup ⇒ cannot displace a VG transcript (the prior in-flow
  displacement leak is structurally excluded).

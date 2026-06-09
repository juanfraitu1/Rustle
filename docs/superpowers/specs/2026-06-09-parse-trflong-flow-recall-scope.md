# parse_trflong flow work — recall-direction scope (diagnostic-gated)

**Date:** 2026-06-09
**Status:** Design (Phase 0 not yet started)
**Author:** scoping session
**Supersedes framing of:** `2026-05-28-pathpat-flow-parity-scope.md` (that doc scoped the
*precision* direction — over-enumeration FPs; this scopes the *recall* direction —
non-generation of annotated chains). Read that doc first for the flow-internals map.

## Problem

Genome-wide, rustle misses **668 annotated isoforms that StringTie recovers** (gffcompare
FSM `=`); see `bench/copy_recovery_eval/results_genomewide/ST_ONLY_MISS_ANALYSIS.md`. The
container sub-bucket was traced with the cross-tool `path_extracted` infra: **~50% are pure
non-generation** — rustle's flow (`parse_trflong` / `path_extract.rs`) never enumerates the
chain that ST produces. The coverage *scale* divergence (1000× tlen-inflation) is a red
herring (the gates are cov-*ratios*); the `retained_intron` gate is a faithful ST port and
rarely the divergence point. The recall gap is in **flow path-enumeration**, the same
architectural lever the prior precision scope isolated — now confirmed genome-wide on
annotated transcripts.

rustle already has a partial `parse_trflong_st.rs` port (env-gated). The two pivotal slices
are stubbed: `seed_order_st` ("NOT YET IMPLEMENTED") and `update_abundance_st` ("TRUE STUB —
always returns None"). "The parse_trflong flow work" concretely = completing the implicated
stub(s).

## Goal & success criterion

**Net-F1-positive.** Recall gains must not regress precision overall. The prior precision
scope established that rustle's over- and under-extracted chains are a MIX of FP and TP that
is **indistinguishable by any available feature** (longcov, min_jct_mm, cov all overlap).
Therefore: generating ST's missing real chains will also generate spurious ones, and net-F1
is only achievable if a **discriminator exists** that separates them. Proving that
discriminator exists — or doesn't — is the central job of Phase 0, BEFORE any flow rewrite.

## Approach: diagnostic-gated targeted completion

Two stages. No flow code is written until Phase 0's gate passes.

### Key infrastructure decision

Phase 0 runs cross-tool comparison on **per-locus BAM slices** (one per st_only miss,
~150 reads, sub-second/tool), NOT genome-wide whole-chromosome runs. This avoids the WSL2
OOM that bit the genome-wide `--vg` runs (see `project_genomewide_oom_protocol`), runs the
668-case sweep in batched minutes, and reuses the validated `/tmp/xtool_*.py` harness.
Aggregate Sn/Pr oracle numbers come from the per-chrom GTFs already on disk
(`results_genomewide/perchrom/*/{rustle,st}.gtf`, `r_fsm/s_fsm`) — zero recompute.

## Phase 0 — diagnostic (low-risk, read-only)

Four measurements over existing per-chrom artifacts + locus-sliced cross-tool runs.

**(a) Generation census.** For each of the 668 st_only misses: slice both BAMs, run both
tools with `path_extracted` logging, classify in rustle as *generated* (chain in
`path_extracted`) vs *non-generated*. Output: genome-wide split into "flow never produced
it" vs "produced then filtered", broken down by the existing category
(`st_only_classified.jsonl`: partial_multi / altsplice / container / contained).

**(b) Recall-oracle ceiling.** Compute Sn/Pr if rustle adopted ST's exact extracted chain
set at these loci (add every ST-generated chain rustle misses). Mirror of the prior
precision oracle (−2.2 Sn / +4.2 Pr). Yields the upper bound on recall AND the precision
cost — how many *non-annotated* chains ST also extracts at these loci that ride along.

**(c) Separability test — THE GATE INPUT.** For chains ST generates but rustle doesn't,
build feature vectors (`longcov`, `entry_abund`, `seed_tf` provenance, pathpat reachability,
junction read-support, n-exons) and test whether any feature/threshold separates the
**real-missing** (annotated → st_only) from the **spurious-missing** (ST-extracted but
unannotated) chains. The prior work found no discriminator for the precision direction; we
need one for the recall direction or net-F1-positive is impossible.

**(d) Mechanism attribution.** For the non-generated set, instrument rustle's seed loop
(`path_extract.rs` ~5863 seed entry / ~6507 extraction) to log WHY each ST chain is not
produced:
  - seed never selected → `seed_order_st`
  - extension stopped short → `onpath_long` / back-fwd extension gates
  - seed starved of capacity → `update_abundance_st` stub
Identifies which stub Phase 1 completes.

## Go/No-Go gate (pre-registered)

Proceed to Phase 1 **only if BOTH**:
1. **Material recall ceiling** — recall oracle (b) recovers ≥ ~150–200 annotated isoforms
   genome-wide (enough to matter at scale), AND
2. **A discriminator exists** (c) — a threshold where recovered-real / added-spurious is
   favorable enough to hold precision flat-or-up (net-F1-positive).

**If the gate fails** (no discriminator / small ceiling / pure Sn↔Pr reshape): Phase 0 IS
the deliverable — a rigorous negative result closing the question, like the prior bucket-A
lowintron oracle (+0.1 Pr → falsified, saved high-risk surgery). This is a success, not a
failure.

## Phase 1+ — contingent implementation

Only if the gate passes. Deliberately narrow.

- Complete ONLY the stub the attribution (d) implicates (`seed_order_st`, the extension
  gates, or `update_abundance_st`) — not a wholesale `parse_trflong` rewrite.
- Env-gated (`RUSTLE_*`), default-off, default output **byte-identical**.
- Validated by **parity-diff to zero divergence** against ST's `path_extracted` on a frozen
  locus set BEFORE trusting any F1 number (prior lesson: per-phase F1 is unreliable —
  redist/firstcov each regressed alone).
- If the (c) discriminator is needed to stay F1-positive, it ships as a gate on the
  newly-generated chains.
- **Abort condition:** if the stub achieves path parity but does NOT yield the predicted
  recall, the divergence is deeper (graph construction) — stop and re-scope.

## Validation & risk

- **Determinism:** all harness scripts sort deterministically; no `Date`/`random` (matches
  the genome-wide protocol constraints).
- **OOM safety:** locus slices only, serial / small batches — no concurrent whole-chrom
  `--vg`. Arm a `free -m` watchdog per `project_genomewide_oom_protocol` if any larger run
  is needed.
- **Risk:** Phase 0 low (read-only). Phase 1 touches core flow (can disturb the ~1700
  existing matches) — mitigated by default-off gating, zero-divergence parity validation,
  and the abort condition.
- **Reuse:** `parity_decisions` infra (`STRINGTIE_PARITY_LOG` / `RUSTLE_PARITY_LOG`,
  `tools/parity_decisions/diff.py`), the `/tmp/xtool_*.py` harness, per-chrom artifacts.
  ST must be built `make clean release` (10s/run; `touch rlink.cpp` to force WSL recompile);
  `-o /dev/null` breaks ST — use a real output path.

## Cross-mode note

The premise is "fix baseline non-generation → propagates to VG." Phase 0 measures on
baseline `-L`; if the gate passes, a secondary VG check confirms propagation on the
paralog/family loci (RABL2/RBMY) where VG recall matters most.

## Out of scope

- Full `parse_trflong` byte-parity rewrite (the prior doc's Approach B) — only the implicated
  stub.
- The precision-direction over-enumeration gap (covered by the 2026-05-28 doc).
- The non-flow recall levers (strand-aware bundling `RUSTLE_STRAND_PURE_MINORITY`) — separate,
  parallel track.

# Scope: rustle == StringTie baseline parity (the alt-splice path-extraction re-port)

**Date:** 2026-06-14. **Status:** scope for decision — NOT yet a task-by-task plan.
**Goal (user principle):** baseline rustle (no flags) emits the *same transcripts* as StringTie; rustle's enhanced sensitivity moves behind a flag (`RUSTLE_PRECISE` / `--vg`). **Success metric:** rustle-vs-StringTie transcript divergence (`bench/divergent_loci.py`) → 0, accepting transient regression-vs-annotation per stage.

## 1. The honest magnitude (read this first)

This is **the central hard problem of the rustle port**, not a fresh tractable task. It has been investigated across many prior sessions and is documented as a **large, coupled, multi-session port** — not a flag-flip and not a filter tune. Every surgical filter/threshold lever has been falsified. And it is, by the user's own framing, a **lateral move on truth**: ~90% of rustle's "extra" isoforms are read-backed real splices, so converging to StringTie *deliberately discards real recall* (acceptable under the faithful-baseline principle, but that is the trade).

## 2. What the divergence actually is (from the 2026-06-14 understanding sweep)

- **It is over-enumeration of MINOR alt-isoforms, not dominant-isoform selection.** The dominant/backbone isoform is shared with StringTie at **82–93% of divergent genes**; rustle emits 1.04–1.13× StringTie's transcript count, the extras being low-coverage, ~70%-novel-junction minor isoforms.
- **The root cause is UPSTREAM of the selection filter.** rustle accepts junctions *strict-early* (its accepted-junction set is a ~42% SUBSET of StringTie's), builds a finer graph, and fragments reads into ~2× more thin transfrags/seeds → a thin dominant flux → weak depletion → a different surviving isoform set. StringTie is *permissive-early + heavy-downstream*. This single architectural inversion drives **both** the precision gap (rustle-only extras) and the recall gap (ST-only chains rustle apportions too little coverage to and then kills).
- **The selection filters are already ported** line-for-line (`predcluster_st.rs` / `isofrac_st` == StringTie `rlink.cpp:19149-19217`; `transcript_filter.rs` print_predcluster). They are **default-off** and, turned on standalone, **over-fire** — because the upstream flow over-enumerates, not because the predicates differ.

## 3. What is already done / in place

- **Graph-layer convergence (clean, shipped, byte-identical default):** ST donor-snap + acceptor-snap alt-junction snapping → chain divergence 186/80 → 162/76 (−28). Graph-node decisions aren't flow-coupled, so these were clean.
- **Pre-assembly wins:** strand-bundling (`RUSTLE_STRAND_PURE_MINORITY`, +48 FSM genome-wide), micro-node smart-trim.
- **Filters ported** (predcluster_st, isofrac_st, kill_included, kill_retained_intron) — present, default-off.
- **The flip architecture:** `precise_mode()` (`RUSTLE_PRECISE` = strict-early escape hatch byte-identical to commit 4705ab1; default = ST-faithful convergence target; ~18 `!precise_mode()` sites).
- **Bidirectional parity instrumentation:** rustle + the committed instrumented StringTie both emit `parity_decisions` JSONL (junction_accept, path_extracted, node_flux, transfrag_define, seed_reject, checktrf_result, pred_kill, …), diffed by `tools/parity_decisions/diff.py` on (step,start,end,strand). Convergence measured by `bench/divergent_loci.py` (RU-only vs ST-only intron chains).
- **Prior scope specs:** `docs/superpowers/specs/2026-06-11-coupled-flow-downstream-port-design.md`, `2026-06-12-coupled-faithful-port-scope.md`.

## 4. The unsolved core (the blocker)

**Flow apportionment.** rustle's thin-seed over-enumeration produces a weak dominant flux that under-depletes shared coverage, so minor isoforms retain flow and survive (rustle-only) while ST-only chains get too little coverage and die. The ST-faithful flow builders that would fix this (`back_to_source_fast_long_st`, `long_max_flow_st`, `RUSTLE_FLOW_ST`) are **"close-but-broken"** — they under-deplete / over-extract low-cov backbone, so routing through them **regresses both directions** (186→223 default; 294→371 with junction-layer on). This is the documented flow-decomposition ceiling, and it is **coupled**: partial mimicry composes negatively (ST-faithful junction acceptance *alone* regresses chain divergence 186→322), so junction acceptance + flow + filters must move *together*.

## 5. Phased roadmap (coupled — stages validated together, transient regression expected)

- **Phase A — junction-acceptance convergence.** Flip rustle strict-early → StringTie-permissive (the `RUSTLE_ST_JUNC` umbrella + graph snaps already converge accepted junctions 7318 → ~12535 of ST's 17487). Chain divergence regresses transiently (intended) — most of the new ST-only gap becomes *flow*, not junction-blocked.
- **Phase B — flow-apportionment re-port (THE core).** Fix the `_st` flow builders so the converged-junction graph apportions coverage like StringTie (heaviest-path order, depletion arithmetic, checktrf redistribute/independent-store). This is the unsolved, broad-blast-radius work (re-apportions all ~1741 shared-chain coverages). Without it, Phases A and C only regress.
- **Phase C — filter activation.** Turn on the already-ported predcluster_st / isofrac_st under the converged flow; they stop over-firing once coverage apportions correctly.
- Each phase tracked by `bench/divergent_loci.py` (drive RU-only + ST-only → 0) + the parity JSONL diff; `RUSTLE_PRECISE` byte-identity to 4705ab1 re-verified each step.

## 6. Recommended bounded first slice (before committing to the whole grind)

Do **not** open the full multi-session port blind. The tractable, measurable first step is to **establish the live convergence baseline + isolate the Phase-B flow gap on one locus**:
1. Run the parity capture (`bench/capture_parity.sh` on chr19) + `bench/divergent_loci.py` to get the *current* RU-only/ST-only counts at HEAD (the prior numbers predate this session's work).
2. Pick the single highest-divergence multi-exon locus, diff the `parity_decisions` `path_extracted` + `node_flux` events rustle-vs-StringTie, and pin *exactly* where the flow apportionment diverges (seed order vs depletion vs checktrf).
3. Attempt ONE targeted flow-apportionment fix at that locus; measure chain-divergence movement + RUSTLE_PRECISE byte-identity. If it converges cleanly → the Phase-B approach is viable, proceed; if it regresses (as prior `_st`-builder attempts did) → the core is confirmed intractable-by-slice and we reassess scope.

## 7. Decision

Committing to this is committing to the project's hardest, multi-session work, for a lateral-on-truth outcome (StringTie-reproducibility, not better transcripts). The bounded first slice (§6) de-risks that decision: it costs little and tells us whether Phase B is viable before the full grind. Recommend: **do the bounded first slice, then decide on the full port from the live numbers.**

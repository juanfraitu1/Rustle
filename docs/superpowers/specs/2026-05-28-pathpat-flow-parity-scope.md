# Scope: PATHPAT / flow path-selection parity with StringTie

Status: SCOPING (not started). Author session: 2026-05-28.
Branch context: builds on `parity/isofrac-chain-dedup` (chain-dedup committed; Sn 96.5 / Pr 90.7).

## 1. Goal & target

Close the dominant precision gap: **~180 over-enumeration false positives** Rustle emits that
StringTie does not.

- 120 `j`-class (novel-junction-combo) extras + ~60 `c/m/k/n/o` (containment/retained-intron-class)
  extras. These are the SAME phenomenon — alt-isoform over-enumeration — gffcompare just labels them
  differently by containment relationship to the reference.
- Spread across **77 loci**, concentrated in complex multi-isoform genes (RSTL.415: 8, RSTL.398: 6).
- Source split (current code): **82 from `flow`** (main path-extraction), **38 from `checktrf_rescue`**.

Target: eliminate these without losing matched chains. Upper bound if fully closed: Pr ~90.7 → ~95+,
F1 well above the current 93.5 ceiling. This is the only remaining lever with that magnitude
(all post-hoc filter levers verified exhausted — see project_coverage_metrics_deviation memory).

## 2. Mechanism (why Rustle over-enumerates)

Rustle's flow path-extraction generates alt-junction paths ST never traverses. Root cause is in the
path-selection logic, NOT in any downstream filter (retained-intron/isofrac filters already fire at
ST's scale and count — verified 2026-05-28: Rustle 1798 retained-intron kills vs ST 1601, same
per-base cov scale).

Two contributing engines:
1. **Flow seed extraction** (`parse_trflong`): `back_to_source` / `fwd_to_sink` build a path pattern
   (`pathpat`) by OR-ing in compatible transfrags' patterns. Rustle's `onpath_long` gate is more
   permissive than ST's, so chimeric/stub transfrags pollute `pathpat`, and the flow decomposes into
   extra alt-isoform paths. (82 extras.)
2. **checktrf rescue + abundance redistribution**: seeds deferred to checktrf get rescued and their
   abundance redistributed onto sub-seeds, which grow into alt-isoforms ST doesn't emit. (38 extras.)

## 3. Why prior / local fixes failed (READ THIS before touching gates)

Existing env-gated gates are **NO-OPS on the headline** in current code (measured 2026-05-28):

| gate | tx | Sn | Pr | j-extras |
|---|---|---|---|---|
| default | 1974 | 96.5 | 90.7 | 120 |
| RUSTLE_STRINGTIE_EXACT=1 (strict onpath_long) | 1973 | 96.5 | 90.8 | 120 |
| RUSTLE_STUB_GATE_PASTMAXP=1 | 1974 | 96.5 | 90.7 | 120 |
| RUSTLE_BACK_SKIP_PAST_SEED_MAX=1 | 1974 | 96.5 | 90.7 | 120 |

The `onpath_long` strictness divergence (ST `rlink.cpp:7434` `if(trnode.Last()>maxp) return false`
vs Rustle `path_extract.rs:3528-3530` `node_can_reach` relaxation) is real but **insufficient alone**:
the over-enumeration is a multi-seed cascade with redundant entry points. Fixing one seed's decision
leaves the alt-isoform produced by other seeds' redistribution. Conclusion from RSTL.182 trace
(memory project_rstl182_pathpat_or_cascade): "don't search for the one fix; coordinate across all
contributing seeds." => This must be a **coherent stage-by-stage port**, not a knob.

## 4. Surface area (functions to bring to ST-parity)

| Concern | StringTie (rlink.cpp) | Rustle |
|---|---|---|
| read→transfrag build (pattern, abundance) | `update_abundance` 4367 | scattered: transfrag_process.rs, path_extract.rs (`update_abundance_st` 995 is a stub) |
| path-compat gate | `onpath_long` 7416 | `onpath_long` path_extract.rs:3466 |
| seed extraction loop | `parse_trflong` 9681 | `extract_transcripts` seed loop path_extract.rs:6080+ |
| back/fwd path build + PATHPAT_OR | back/fwd inside parse_trflong | `back_to_source_fast_long` 4913, `fwd_to_sink_fast_long` 3943, OR sites 2043/3295/4867/5794 |
| max-flow | bfs flow ~8578 | `long_max_flow_seeded_with_used_pathpat_st` max_flow.rs:428 (already faithful — verified) |
| checktrf rescue + redistribution | checktrf handling in get_trf_long | path_extract.rs:9216+ (`original_abundances` 6517, `best_trf_match` port 2761, redistribution 1429/2983) |

Note: max-flow math is already line-for-line faithful (verified stage-3, 2026-05-28). The divergence
is in the INPUTS (transfrag set/abundance/pattern) and the path-building gates — i.e. rows 1, 2, 4, 6.

## 4b. PHASE 0 RESULTS (2026-05-28) — attribution table built; scope REVISED

Harness: `bench/pathpat_phase0.py` / `phase0c.py`. Captured path_extracted + pred_kill from
both tools (full chr19), attributed each of the 186 multi-exon extras. **Output reshapes the plan:**

| bucket | count | what it is |
|---|---:|---|
| **A: filter divergence** | **97 (52%)** | ST extracts the EXACT chain then kills it; Rustle keeps it. ST kill reason: retained_intron 31+, isofrac 16+, included_drop 4+ (46 coord-shifted, same 3 reasons). |
| B: path-enum | 89 (48%) | ST never extracted the chain. Of these: |
| — alt-splice diverge | 48 | Rustle picks a shifted donor/acceptor vs ST (a miss AND an extra) |
| — subset/contained | 32 | Rustle drops a junction ST keeps |
| — genuine extra-junction | **9** | true over-generation |

**Key correction**: only ~9 extras are genuine over-generation. The original "Rustle generates ~180
paths ST doesn't" framing was WRONG. Half (97) is FILTERING (both tools produce the path; ST kills,
Rustle doesn't), and most of bucket B is junction SELECTION (wrong donor/acceptor), not extra paths.

**Caveat on bucket A**: Rustle's retained_intron filter already fires MORE than ST globally (1798 vs
1601) and isofrac fires — so the divergence is PER-TRANSCRIPT (filter doesn't fire for THESE 97). Root
is the per-transcript filter INPUTS: coverage values (flow-allocation divergence), lowintron flags, or
the n1-pairing (which higher-cov neighbor). Phase 1 must determine which.

Full table: `/tmp/phase0_attribution.tsv` (regenerate via the bench scripts).

## 5. Phased plan (REVISED per Phase 0 — each phase env-gated, parity-diff validated before F1)

REVISED ORDER (was: transfrag-build → flow → checktrf). Now target the biggest+most-tractable bucket
first:
- **Phase 1 (NEW): filter parity on co-extracted paths (bucket A, 97 extras).** For the 97 chains ST
  extracts-and-kills but Rustle keeps, determine per-transcript why Rustle's retained_intron/isofrac/
  included_drop doesn't fire (coverage value? lowintron flag? n1-pairing?). Align the firing. Most
  tractable — no flow change, both tools already produce the path. Sub-investigate retained_intron
  (31) first (largest, and the included_pred/has_retained_intron logic is well-localized).
- **Phase 2: alt-splice junction selection (bucket B, 48).** Rustle picks shifted donor/acceptor. This
  is the documented alt-splice direction bias (project_missed_tx_breakdown: 78 cases, smart-trim got
  some). Junction-selection in the flow / graph layer. Recovers misses AND removes extras (2-for-1).
- **Phase 3: subset/contained (32) + genuine over-enum (9).** The residual flow path-selection work
  (the original onpath_long/PATHPAT_OR scope). Smallest payoff; do last.

### (superseded) original phase plan below

**Phase 0 — Diagnostic harness (1 session).** Extend the parity-decisions trace with a per-seed
path-decision record: seed idx, `maxpath`, `pathpat` membership, `onpath_long` pass/fail per candidate
transfrag, back/fwd edge choices. Add the same to ST. Freeze a test set: the 77 over-enum loci +
existing `bench/path_decision_diff/RSTL182_trace/`. Build a **per-extra → originating-decision
attribution table** (which seed + which divergent decision produced each of the 180 extras).
EXIT: every extra traced to a specific Rustle-vs-ST decision divergence.

**Phase 1 — transfrag-build parity (rows 1).** The transfrag set/abundance/pattern feed everything;
if they differ, the flow diverges regardless of downstream. Diff `transfrag_define` events ST vs
Rustle on test loci. Port divergences (`update_abundance_st` is currently a stub).
EXIT: transfrag sets match ST on test loci.

**Phase 2 — flow path-selection parity (rows 2,4).** Port ST's `onpath_long` + back/fwd + PATHPAT_OR
semantics as a UNIT. The bet: with the transfrag inputs correct (Phase 1) and the gate strict, the
matches that `RUSTLE_STRINGTIE_EXACT` lost in isolation are retained (they came from ST-faithful
behavior, not the relaxation). EXIT: the 82 `flow`-source extras eliminated, no match loss, on test loci.

**Phase 3 — checktrf rescue parity (row 6).** Port ST's checktrf semantics so rescue + redistribution
doesn't synthesize alt-isoforms ST doesn't. EXIT: the 38 `checktrf_rescue` extras eliminated.

**Phase 4 — integration + full-genome sweep.** Run full GGO_19, confirm net F1-positive, no Sn
regression. EXIT: Pr materially up, Sn held.

## 6. Validation methodology

The engine is the parity-decisions JSONL diff (both tools instrumented; ST binary built at
`tools/stringtie/stringtie`). For each phase: run both tools with `*_PARITY_LOG`, join on
`(step,chrom,start,end,strand)`, drive the divergence count to ZERO on the frozen test set BEFORE
measuring full-genome F1. Do NOT trust F1 alone per phase (each phase may regress in isolation — the
redistribution stage proved this). The parity-match is the per-phase gate; F1 is the final gate.

## 7. Risk & mitigation

- **High risk**: touches the core flow algorithm; can disturb the 1738 currently-matched chains.
- **Each phase may be F1-negative alone** (proven pattern: redist stage, firstcov). Mitigation:
  env-gate every phase; keep `parity/isofrac-chain-dedup` as the safe fallback; only flip a phase
  default-on after the full stack is net-positive.
- **The "lost 7 matches" trap**: tightening onpath_long alone loses matches because downstream
  behaviors depend on the looseness. Phase ordering (transfrag-build FIRST) is designed to remove that
  dependency before tightening.
- **Cascade redundancy**: a single seed fix is a no-op. Phase 2 must port the whole gate+build unit,
  validated by the attribution table from Phase 0 (confirm ALL contributing decisions flip).

## 8. Effort & exit criteria

Multi-session / multi-PR. Rough: Phase 0 ≈ 1 session, Phases 1–2 ≈ 2–4 (the hard core), Phase 3 ≈ 1–2,
Phase 4 ≈ 1. Abort criterion: if after Phase 2 the flow extras don't drop with parity-match achieved,
the divergence is deeper (graph construction) and the project should be re-scoped.

## 9. Existing assets

- Parity-decisions infra wired both sides (reference_parity_decisions_infra memory).
- Instrumented ST binary built (env-gated traces ST_TRACE_COV_NODES, parity log).
- `bench/path_decision_diff/RSTL182_trace/` slice + TSVs.
- `best_trf_match` already ported (path_extract.rs:2761 for checktrf, `best_trf_match_redist` from
  the redist stage).
- No-op gates (STRINGTIE_EXACT, STUB_GATE_PASTMAXP, BACK_SKIP_PAST_SEED_MAX) as scaffolding/reference.

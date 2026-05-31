# Design: StringTie-faithful source/sink (terminal) edge-wiring alignment

Status: DESIGN (approved 2026-05-30). Next: implementation plan (writing-plans).
Grounding: the foundation parity audit (findings §6j "Foundation parity audit" + "Graph edges now
MEASURED"; memory `project_color_cgroup_parity`); the over-enumeration floor
(`project_jfp_missr_characterization`); the coverage/flow divergence
(`project_coverage_metrics_deviation`, `project_pathpat_flow_parity_scope`).

## 1. Problem & insight

The foundation parity audit established that Rustle starts from StringTie's graph substrate — bundles,
bundlenodes, junction-spanning connectivity, colors, and raw junctions all match (every divergence is
the known mm_negative/support floor, an inert artifact, or RU being more precise). The DIRECT
graph-edge measurement (`bench/edge_diff.py`) pinpointed the FIRST real point of divergence:
**source/sink (terminal) edge wiring**.

- On the 3226 identical-node shared bundles, edge sets are identical on 3037 (94.2%). Junction-spanning
  connectivity is the same standpoint (exactly 1 RU-missing junction edge, and it is a CORRECT
  alt-donor omission — RU more precise).
- The divergence is terminal wiring: **RU emits 2708 EXTRA source/sink edges** vs ST (1308 on
  identical-node bundles, a near-perfect superset — only 1 missing).
- Same coverage-drop threshold (`ERROR_PERC * DROP = 0.05`) but a DIFFERENT RULE:
  - RU `add_coverage_source_sink_edges` (`graph_build.rs:68-272`): a blanket per-node parent/child
    coverage-ratio test on EVERY node, PLUS a phantom-zero-cov recursion pass (`:146-164`, `:185-202`)
    documented as compensating for phantom zero-cov nodes RU creates. Toggle exists:
    `RUSTLE_COVLINK_RECURSE_ZERO_OFF=1` disables the recursion; `RUSTLE_COVLINK_THRESHOLD_LOOSE`
    changes the threshold (off by default = 0.05, matching ST).
  - ST spreads terminal wiring across THREE narrower mechanisms: structural (source iff no real parent
    `parent.Count()==1 && !parent[0]`, sink iff no children — `rlink.cpp:2539,2581-2583`); `find_trims`
    (intra-exon coverage cliff via a `CHI_WIN=100`/`CHI_THR=50` sliding window — `rlink.cpp:1500-1601`);
    and a `create_graph` coverage-based source/sink check.

**Insight:** RU's single blanket rule + phantom recursion is more permissive → more flow start/end
points → feeds the known over-enumeration. This is the on-ramp to the flow stage and the first place
RU stops matching ST. BUT: the over-enumeration *junction* edges (161 RU-extra) sit on
node-MISMATCHED bundles, which the injection oracle cannot cleanly reach — so the realizable
terminal-wiring prize may be modest. The oracle measures this before any build.

## 2. Goal & non-goals

- **Goal:** measure — oracle-first, default-OFF — whether aligning RU's source/sink wiring to ST's
  reduces over-enumeration FPs and at what TP cost; decide ship/shelve from the oracle.
- **Decision rule (user-chosen):** oracle-bound first, then decide. Hard abort gate after Phase 0.
- **Non-goals:** node-decomposition parity on node-mismatched bundles; the full flow/path-extraction
  port; the coverage *substrate* divergence (documented elsewhere); `find_trims`/structural re-port
  unless Phase 0 shows it is the lever.

## 3. Architecture — three phases, hard abort gate after Phase 0

- **Phase 0a — Free probe** (no new code; existing toggles + committed `bench/edge_diff.py` +
  `bench/capture_parity.sh`). Run `RUSTLE_COVLINK_RECURSE_ZERO_OFF=1`; measure how many of the 2708
  extra terminal edges vanish and the F1/FP delta. Isolates the phantom-recursion contribution and the
  F1 direction.
- **Phase 0b — Injection oracle** (one env-gated analysis harness, default OFF). On identical-node
  bundles, override RU's terminal edges with ST's exact captured source/sink set; run flow+selection;
  measure causal FP-reduction / TP-cost. The prize ceiling for the injectable subset.
- **Abort gate.** Abort if net (FP-reduction − TP-loss) ≤ 0, or realizable net < ~5 and concentrated
  where injection cannot reach. Proceed only on a clear F1-positive net with FPs causally tied to
  removed terminal edges.
- **Phase 1 (only if gate clears) — align the rule** behind `RUSTLE_ST_TERMINAL` (default OFF),
  targeting whatever 0a/0b identified as the over-firing mechanism.

## 4. Components

### Phase 0a — Free probe (existing toggles)
- Capture: `RUSTLE_COVLINK_RECURSE_ZERO_OFF=1 RUSTLE_PARITY_LOG=/tmp/ru_edge_norecurse.jsonl
  RUSTLE_PARITY_FILTER_STEPS=graphnode_list,graph_edge RUSTLE_PARITY_FILTER_CHROM=NC_073243.2
  ./target/release/rustle -L GGO_19.bam -o /tmp/ru_norecurse.gtf`.
- `bench/edge_diff.py /tmp/ru_edge_norecurse.jsonl /tmp/st_edge.jsonl` → extra-edge reduction.
- gffcompare `/tmp/ru_norecurse.gtf` vs reference → F1/FP delta vs 95.6/90.5.
- Confirm `RUSTLE_COVLINK_THRESHOLD_LOOSE` is off (0.05 == ST), so divergence is rule/recursion not
  threshold.

### Phase 0b — Injection oracle (`RUSTLE_TERMINAL_ORACLE`, default OFF)
- New env-gated mode in the graph stage: after `add_coverage_source_sink_edges`, for each bundle whose
  node set is byte-identical to ST's (match by envelope + node-coordinate set), DISCARD RU's SRC/SNK
  edges and INSTALL ST's (parse the captured `graph_edge` SRC/SNK tokens from
  `RUSTLE_TERMINAL_ORACLE=<path to st_edge.jsonl>`, map node-coordinate tokens to RU node ids).
  Node-mismatched bundles keep RU's wiring (cannot inject cleanly). Must also reconcile the synthetic
  source/sink GraphTransfrag abundances so flow has consistent inputs (mirror the
  `add_coverage_source_sink_edges` abundance derivation for installed edges).
- Run flow + selection unchanged → final GTF → gffcompare.
- `bench/terminal_oracle_report.py`: edges injected, bundles touched, FP/FN/F1 delta vs baseline, and
  how many removed FPs are path-endpoint-attributable to a removed terminal edge (causal check).

### Phase 1 — align the rule (`RUSTLE_ST_TERMINAL`, default OFF; only if gate clears)
- Targeted at the identified over-firing mechanism — most likely: disable the phantom-zero recursion
  under the flag, and tighten the blanket interior-node check toward ST's structural rule
  (no-real-parent→source / no-child→sink) + `find_trims` intra-exon cliff. Composes with the existing
  `RUSTLE_COVLINK_*` toggles.

## 5. Data flow (oracle/flag ON)

reads → graph (nodes identical to ST) → **[terminal wiring: ST-aligned via injection (0b) or rule (1)]**
→ flow → selection → final GTF.

## 6. Validation (parity + F1)

- **Phase 0a:** diagnostic — edges-removed + F1 delta; pinpoints the culprit mechanism. No pass/fail.
- **Phase 0b (abort gate):** net = FP-reduction − TP-loss on the identical-node subset. Abort if ≤ 0,
  or realizable net < ~5 concentrated where injection can't reach. Proceed only on a clear F1-positive
  net with FPs causally tied to removed terminal edges.
- **Phase 1 (if reached):** `edge_diff.py` shows RU terminal edges converging toward ST (2708 extra →
  small); F1 measured. Ship only on a confirmed F1 gain; else keep `RUSTLE_ST_TERMINAL` opt-in + record
  the cost.
- **Default-unchanged guard** at every step (all flags default OFF; baseline 95.6/90.5 verified).

## 7. Safety & exit

- Free probe = existing default-OFF toggles → zero default risk. Injection oracle behind
  `RUSTLE_TERMINAL_ORACLE` (reads a captured file) → zero default risk. Phase 1 behind
  `RUSTLE_ST_TERMINAL` (default OFF). Default verified each step.
- **Abort at 0b (plausible):** record the terminal-wiring prize ceiling + which mechanism over-fires +
  where the over-enum FPs actually live (likely node-mismatched bundles = a deeper divergence). A
  successful bounding result that redirects to the real lever or confirms the floor.
- **Proceed → Phase 1:** measure convergence + F1; ship only on F1 gain, else opt-in.
- **Honest expectation (recorded up front):** prior flow experiments were F1-negative and the over-enum
  junction edges sit on node-mismatched bundles the injection can't reach — so the realizable
  terminal-wiring prize may be small. That outcome is still valuable: it would prove the over-enum root
  is the node-decomposition / flow-traversal on those bundles, not the terminal wiring per se,
  sharpening the next target.

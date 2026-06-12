# Coupled Faithful Port (junction acceptance → flow → filtering) — Scope

**Status:** SCOPE (multi-session, user-committed 2026-06-12)
**Branch:** vg/flow-capacity-apportionment · HEAD 1da762e
**Goal philosophy:** [[feedback_baseline_equals_stringtie]] — rustle reproduces StringTie EXACTLY as the default baseline (incl. ST's permissive junctions, FP/missed isoforms); rustle's strictness = opt-in levers. Escape hatch `RUSTLE_PRECISE` = today's stable rustle, must stay byte-identical to 4705ab1.

## Why coupled (the core finding)

rustle's accepted junctions are a STRICT SUBSET of ST's (chr19: rustle 7318 ⊂ ST 17487, rustle-only=0). ST accepts 10169 more, 85% single-read. **Aligning junction acceptance ALONE converges the junction SET but regresses chains 186/80 → 322/150** (tested, `RUSTLE_KEEP_MM_NEG=1`): rustle emits chains from the new junctions that ST's downstream flow/filtering kills but rustle's doesn't. The three stages are co-designed in ST and cannot be flag-separated (memory project_flow_input_parity §6u). So the faithful baseline = re-port all three together, accepting **transient** chain regression per layer; **success metric per layer = SET/decision convergence vs ST, NOT chain divergence** (chains recover only when all three layers land).

## Layers

### Layer 1 — Junction acceptance (accept ST's superset)
Make rustle accept ST's 17487 junctions. Two rejection mechanisms (chr19, the 10169 ST-only):
- **mm_negative (3470):** rustle sets mm=-1 → rejects (graph_build.rs:843 `reject_reason="mm_negative"`). Toggle EXISTS: `RUSTLE_KEEP_MM_NEG` (keeps mm_neg w/ nreads_good>0) → converges 10169→6699 ST-only. mm=-1 is set across killed_junctions.rs + pipeline.rs apply_bad_mm_neg_stage (2159) — KEEP_MM_NEG at the graph stage suffices to keep them.
- **single-read filter (6679):** single-read junctions (nreads_good=1) reach `junction_raw` but are dropped BEFORE `good_junction` (NOT `filter_weak_junctions`/min_junction_reads=1.0 — that's already ST-faithful). ⚠ EXACT FILTER LOCATION = the remaining Layer-1 sub-task (in rustle good_junc / anchor-support pipeline, bundle.rs count_good_junctions leftsupport/rightsupport vs ST's leftsupport=1 acceptance; ST runs the single-read through good_junction→junction_accept reason="ok"). Find it, relax it gated.
- **Gate:** unify under `RUSTLE_ST_JUNC` (= KEEP_MM_NEG + single-read relax), default off, `!precise_mode()`.
- **Metric:** accepted junction-set convergence (ST-only 10169 → 0; rustle-only stays 0). Expect chain regression (~186→350+); that is the intended transient baseline.

### Layer 2 — Flow on ST's junction set
Once rustle accepts ST's junctions, port ST's flow (parse_trflong/long_max_flow) faithfully on the new (denser) graph so the per-isoform apportionment matches ST. Builds on the prior flow investigation (the 2x-seed / depletion / nodeflux analysis, memory project_coupled_port_viability). Metric: path-population convergence vs ST path_extracted.

### Layer 3 — Filtering
Port ST's downstream pred filtering (readthr/isofrac/RI/included_drop) faithfully on the ST-faithful flow output, so the chains ST kills get killed. The filter PORT is already ~complete (memory: rustle downstream is a superset); the work is making the cov/longcov INPUTS faithful (cov is aligned r=0.95; bpcov/longcov diverges). Metric: chain divergence → 0.

## Gating / invariants (every layer)
- Behind `!precise_mode()` + a per-layer env flag (RUSTLE_ST_JUNC, RUSTLE_ST_FLOW, RUSTLE_ST_FILTER), independently toggleable for isolation.
- `RUSTLE_PRECISE=1` byte-identical to 4705ab1 (checked: `bench/mini3/check_precise.sh`).
- Suite green (`RUSTLE_PRECISE=1 cargo test --release --lib`).

## Metrics & harnesses
- Junction-set parity: junction_accept(accepted=True) sets, rustle vs `/tmp/st_parity.jsonl`.
- Graph node-parity: graphnode_list node-sets (the 2 snaps already shipped, 306/103→284/92).
- Chain divergence: `bench/gtf_chain_diff.py out.gtf /tmp/st_all.gtf` (the eventual all-layers metric).
- ST parity log: `STRINGTIE_PARITY_LOG=/tmp/st_parity.jsonl tools/stringtie/stringtie -L GGO_19.bam`.

## Already done (this branch)
- 2 graph-node faithful snaps (donor b038f8b + acceptor 1da762e): alt-junction snapping, chain 186/80→162/76, gated RUSTLE_ST_GRAPH_SNAP[_ACC]. (Graph-node decisions are NOT flow-coupled → ported clean; junction ACCEPTANCE is coupled → this scope.)

## Next concrete step
Layer 1 sub-task: pin the single-read raw→good_junction filter in rustle's good_junc/anchor pipeline; build `RUSTLE_ST_JUNC` (= KEEP_MM_NEG + single-read relax); measure ST-only junctions 10169→~0; document the transient chain regression; ship gated.

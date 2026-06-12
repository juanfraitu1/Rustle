# StringTie-Parity Milestone & Continuation Plan

**Status:** MILESTONE BANKED (2026-06-12) · Branch `vg/flow-capacity-apportionment` · HEAD **41f1d71**
**Assumption going forward:** parity (rustle-vs-StringTie chain divergence → 0) is *eventually attainable* via the documented continuation steps below. We bank here and pivot to the VG-improvement track.

## 1. Goal & philosophy

rustle (Rust port of StringTie) must reproduce StringTie EXACTLY as the **default baseline** — including ST's quirks, FP, read-through artifacts, and missed isoforms. rustle's better precision/recall is **opt-in** (gated flags), never default. Success metric: rustle-vs-StringTie intron-chain divergence → 0 on chr19 (`GGO_19.bam` vs `/tmp/st_all.gtf`).

Hard invariant (must always hold): `RUSTLE_PRECISE=1` is **byte-identical to commit 4705ab1** (the stable pre-parity rustle). Every faithful-port change is gated behind `!precise_mode()` + a per-layer env flag, so it is independently toggleable and roll-back-able. This is what lets us "prove a decision then roll it back to isolate individual causes."

## 2. What is DONE (committed, faithful, gated)

The coupled port has three layers: (1) junction acceptance, (2) flow on the denser junction set, (3) downstream filtering. **Layer 1 is shipped; Layer 2 is characterized as a ceiling (see §3).**

### Graph-node snaps (committed b038f8b, 1da762e)
- `compute_st_donor_snaps` (`RUSTLE_ST_GRAPH_SNAP`) + `compute_st_acceptor_snaps` (`RUSTLE_ST_GRAPH_SNAP_ACC`): port ST's alt-junction snapping (build_graphs rlink.cpp:14486-14708) using JunctionStat nm/leftsupport/rightsupport/guide_match. Node-parity 306/103→284/92; chr19 chain div 186/80→162/76.

### Layer 1 — junction acceptance (committed 41f1d71)
- **`RUSTLE_ST_JUNC` umbrella** (`!precise_mode`): the single Layer-1 lever. Turns on:
  - the good_junc **alt-donor/acceptor witness exemption** (`killed_junctions.rs` good_junc): skip the left-witness kill when `rightsupport > 5*junctionthr` (alt-donor signature), mirror for alt-acceptor — ST snaps+keeps these and never witness-evaluates them (verified by instrumenting ST's good_junc directly).
  - the two higherr floors relaxed to 1.0 (`RUSTLE_HIGHERR_UNRELIABLE_FLOOR`, `RUSTLE_HE_SMALL_SHIFT_FLOOR`).
  - `KEEP_MM_NEG` (accept mm<0 demotion-marker junctions on raw read support, `graph_build.rs`).
  - Individual env flags still override for isolation.
- **Effect** (`RUSTLE_ST_JUNC=1` + the two snaps): junction acceptance **7318 → 12161 of ST's 17475** (ST-only 10169 → 5316); target alt-donor 20539503 accepted; rustle-only stays 2 (no precision explosion). chr19 chains regress 186/80 → 294/163 = the **intended coupled-port transient** (Layer 1 converges the junction SET; chains recover only once Layers 2–3 land).
- Default (no flags) byte-identical 186/80; `RUSTLE_PRECISE=1` byte-identical to 4705ab1; suite 288/0.

### Dead levers (proven, do not re-attempt as filter slices)
- **strand_mismatch** (the 644/1091 junctions ST accepts, rustle rejects): rustle is one-strand-per-bundle; the filter is load-bearing (relaxing injects wrong-strand edges). ST routes junctions to per-strand graphs — there is no bundle-strand test to relax. The only faithful endpoint is a per-strand dual-graph rearchitecture. Real recoverable miss ~10 weak junctions. (workflow wf_11a5bf5e, unanimous adversarial verdict)

## 3. Layer 2 (flow port) — characterized as the documented broad CEILING

Decisive reframe (workflow wf_2a73b504, flow-trace at concrete loci): the 163 ST-only chains (with Layer-1 flags) are **100% flow-blocked but NOT flow-enumeration failures** — rustle's flow already extracts paths covering every intron, with longcov identical to ST. Three candidate flow slices were implemented + measured, **all falsified**:

| slice | gate | result |
|---|---|---|
| per-path noderate cov | `RUSTLE_ST_FLOW_NODERATE` | catastrophe (in-both 1658→625) — per-path is the flow-capacity rate, not cov noderate; rustle cov already faithful (per-node `n.noderate`) |
| back-extension port | `RUSTLE_ST_BACK_EXT` (endpath=i→maxpath + maxpath bump) | no-op (294→296) — over-decomposition is a multi-seed cascade |
| full canonical `_st` port | `RUSTLE_FLOW_ST` | regresses vs ST (294→371, 163→202) — `_st` builders under-deplete ("close-but-broken") |

The two real divergences are **non-flow**: (a) the 55 extracted-then-killed lose a cov-**ratio** contest = nodeflux apportionment (the documented seed-over-fragmentation broad ceiling); (b) the 103 "never-extracted" are actually extracted-as-**superset** carrying one extra ~100–240bp micro-intron ST never built (junction/transfrag-level).

## 4. Continuation plan (to eventually reach parity)

These are the documented broad-blast-radius ports. Each is gated `!precise_mode` + a per-layer flag, measured by SET convergence (not chains) until all land. Order by tractability:

1. **103 micro-intron supersets (junction-level, UNVERIFIED — verify first).** rustle's path carries a rustle-only ~100–240bp intron ST never built (ST_refs=0). ⚠ The junction-set diff says rustle-only-accepted = 2, which *contradicts* "103 rustle-only micro-introns" — so `ST_refs=0` may be transfrag-refs, not junction acceptance. **Root-cause the mechanism before building any fix** (the strand_mismatch / bucket-A lesson: promising-looking levers here repeatedly evaporate on inspection). If it is a clean ST min-intron / junction filter rustle is missing, porting it recovers ~103 chains without touching flow.
2. **Downstream pairwise RI / included_drop decision parity** (the §3a/§3b coupled pair, prior doc 2026-06-11-coupled-flow-downstream-port-design.md). Cross-tool RI instrumentation already wired (ST `pred_ri_eval` vs rustle `RUSTLE_RI_TRACE`). Diff `ri_result` per (victim,killer) on the readthr-surplus, align the diverging branch, then ship the readthr cov-scale exemption as its safe partner. Reaches chr19 ~253 in prior estimates.
3. **lowintron coverage-accumulation port** (oracle target chr19 258, −6). Bounded + oracle-validated (`/tmp/st_lowintron.oracle`, `RUSTLE_LOWINTRON_DUMP` vs ST `pred_intron_low`) but needs read→bpcov pipeline parity (count intronic-overhang/soft-clip bases like ST's add_read_to_cov; broad blast radius).
4. **nodeflux apportionment / `_st` under-deplete fix** (the deepest, the 55). The `_st` builders over-extract low-cov ~1.0 backbone sub-paths ST depletes; fixing depletion-among-competitors (rlink.cpp:8627-8665) is the real flux-gap lever but is the broad ceiling. Validate against `RUSTLE_TRACE_COV_NODES` vs `ST_TRACE_COV_NODES`.

## 5. How to validate (harnesses & invariants)

- **Chain divergence:** `rustle -L GGO_19.bam -o out.gtf` then `python3 bench/gtf_chain_diff.py out.gtf /tmp/st_all.gtf`.
- **Junction-set convergence:** run with `RUSTLE_PARITY_LOG=/tmp/ru.jsonl`, diff vs `/tmp/st_parity.jsonl` via `bench/junction_set_diff.py`. ST parity log: `STRINGTIE_PARITY_LOG=/tmp/st_parity.jsonl tools/stringtie/stringtie -L GGO_19.bam`.
- **path_extracted / flow convergence:** both tools emit `path_extracted` (with `introns`), `node_flux`, `pred_kill`. Worklist of the 163 misses: `/tmp/layer2_worklist.tsv`.
- **Escape hatch (must always pass before commit):** `RUSTLE_PRECISE=1 bash bench/mini3/check_precise.sh` prints "byte-matches 4705ab1"; suite `cargo test --release --lib` green (288/0).
- **Ops invariants:** ⚠ never `pkill -f rustle` (cwd contains "rustle"; use `pgrep -af target/release/rustle`); whole-genome `-L` OOMs ~18GB → per-chrom serial only; `*.bam`/`*.gtf` gitignored → `git add -f`; `RAYON_NUM_THREADS=1` for clean parity dumps; commit only on request; commit messages end `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>`; **never `git add -A`** — the tree carries a large unrelated VG-family track (family_manifest.rs, global_flow.rs, vg_family/*, etc.) that must NOT be committed; stage source files explicitly.

## 6. Memory pointers
`project_coupled_port_viability` (full detail), `feedback_baseline_equals_stringtie` (philosophy), `project_st_faithful_flip_progress`, `project_flow_input_parity`, `project_flow_parity_scope`, `project_junction_parity`.

# ST-Faithful Baseline (rustle-vs-StringTie → 0) — Design

**Status:** Approved (design choice: *flip-now*). In execution — Milestone 1 partial. Revert point: git `4705ab1` (`main`/current HEAD) — restore today's behavior by reverting to this commit, or in-code via `RUSTLE_PRECISE=1`.

**Progress (2026-06-10 → 06-11):** `precise_mode()` gate live (`bc58844`). Two generalizing default-flips shipped, each gated behind `RUSTLE_PRECISE`, escape-hatch byte-identical throughout, full suite 284/0 under `RUSTLE_PRECISE=1`:
1. **Witness-check gate** (`2385cec`): LR splice-pair witness OFF by default (reproduces ST's cross-gene read-through). chr19 divergence 291→273; mini3 3/16 → 2/7 (locus C converged).
2. **Single-exon-multinode flux depletion** (`cab62d4`): keep SE-multinode seeds (extract → deplete flux) instead of skipping them, so their abundance isn't stranded as residual and dumped onto flow fragments by the redistribution. chr19 273→269; mini3 2/7 → 1/7 (the inflated fragment RI-killed). No SE-FP increase.

**Cumulative: mini3 3/16 → 1/7, chr19 291 → 269.** ST+rustle redistribution instrumentation committed (ST submodule `72a71a1`: `trf_redistribute`/`trf_residual`; rustle `ea8adb9`). Unified root proven: rustle's flow under-depletes → residual → redistribution inflates fragments.

**Tried + reverted (all regress genome-wide):** ST seed-order (`4d92176`→`33b6951`, mini3-overfit), geometry-RI kill, endpoint-demotion, canonical-fold, redistribution-scoping, RI-relax. **Methodology:** validate every mini3 win on chr19 (~7s) before defaulting.

**Milestone 2 (2026-06-11): flow-decomposition co-port — 2 more gates shipped.** Discovery (`docs/superpowers/plans/2026-06-11-st-faithful-m2-flow-decomposition.md`) pinned the remaining-7 root via ST checktrf instrumentation: ST stores the 36063554-terminus chains via `independent_store` at their natural terminus; rustle's **checktrf JAB + fwd-extend** push them past the junction-acceptor node split (36063554 → 36063658 → 36065683), making them retained-intron victims of the dominant chain. Two gates shipped behind `precise_mode()`:
1. **JAB-extension gate** (`...`): chr19 269→268, mini3 → 1/6 (STRG.2.12 converged).
2. **fwd-extend gate**: chr19 268→**267** (−1 rustle-only FP).

3. **Surgical terminal-boundary gate**: the 3rd layer (`apply_terminal_boundary_evidence_to_longread_txs`) over-extends checktrf alt-TES chains to 36065683. Gating it broadly is neutral (un-extends real 3' UTRs); restricting to **checktrf_rescue-source** chains is net-positive — chr19 267→**266**, mini3 → 2/4.

4. **Endpoint-demotion** (`a235a0d`): the fwd soft-gate now also defers a flow seed to checktrf when ST's `_st` clone stops at a strictly-earlier 3' node (rustle over-extends past ST's stop). This was a **no-op before the M2 checktrf gates** (over-extended chains were shared) but **converges now** because checktrf stores the deferred seed at its natural terminus — the **coupled flow+checktrf co-port** validated. chr19 266→**264**.

**Cumulative session: chr19 291→264 (6 gates), mini3 3/16→2/4.** The remaining 4 mini3 ST-only are alt-splice **combination** differences (different intron counts) — `ST_ORDER` converges them (4→1) but overfits (chr19 regresses); canonical still regresses even with the M2 gates. Geometry-based truncation proven dead (post-hoc trim and fwd junction-acceptor-stop both blow up — junctions are too common). The lever is the alt-splice combination selection that generalizes (seed-order/depletion coupling) or more coupled-gate pairs like endpoint-demotion.

**Date:** 2026-06-10

---

## 1. Goal

Make **baseline rustle byte-match StringTie** — including StringTie's *mistakes* (read-through artifacts, FP chains, missed isoforms). rustle's better precision and its extra recall (the 32 real isoforms ST misses, the read-backed novel splices) remain available but become a **deliberate opt-in** (`RUSTLE_PRECISE` flag and/or VG mode), **not** the default.

**Success metric (flipped):** `bench/gtf_chain_diff.py rustle.gtf stringtie.gtf` → **0 rustle-only / 0 ST-only** multi-intron chains, measured **vs StringTie** (not vs the annotation). Target order: mini3 fixture → full GGO_19.bam → genome-wide per-chrom.

This inverts today's architecture: currently rustle's *default* is strict-early/precise and ST-faithful is a partial flag (`RUSTLE_ST_SHADOW` ~layer 2/6, `RUSTLE_PARSE_TRFLONG_ST_CANONICAL`). The target is **ST-faithful default, precision opt-in**.

## 2. Approach: flip-now

A single master gate **`RUSTLE_PRECISE`** is introduced:

- **`RUSTLE_PRECISE=1`** → today's exact strict-early/precise behavior (the in-code escape hatch / revert; must stay byte-identical to commit `4705ab1`).
- **default (unset)** → the **ST-faithful** path (the converging baseline).

Day 1 the default flips to ST-faithful. Because the ST port is incomplete, the **default will regress** (toward the partial-ST-faithful state — `RUSTLE_ST_SHADOW` alone is ~1494 rustle-only) and then **converge back to 0/0** as the coherent stage-groups are ported. This transient regression is **accepted** (per `feedback_full_stringtie_mimicry`: pursue exact-ST decisions per stage, accept transient regressions, do not tune knobs to compensate). The user explicitly chose flip-now with `4705ab1` as the safety net.

**Invariant throughout:** `RUSTLE_PRECISE=1` output stays byte-identical to commit `4705ab1` (the escape hatch never regresses). Verified each task via `gtf_chain_diff(RUSTLE_PRECISE=1 output, 4705ab1 output) == 0/0`.

## 3. Architecture

### 3.1 The `RUSTLE_PRECISE` gate
A helper `precise_mode() -> bool` (read `RUSTLE_PRECISE` once, cache) gates every divergence point. Each precise-vs-ST-faithful fork becomes:
```rust
if crate::stringtie_parity::precise_mode() { /* today's strict-early behavior */ }
else { /* ST-faithful behavior */ }
```
Existing partial-ST flags (`RUSTLE_ST_SHADOW`, `RUSTLE_PARSE_TRFLONG_ST_CANONICAL`, `RUSTLE_KEEPTRF_RET2_LEGACY`, the soft gate `RUSTLE_PARSE_TRFLONG_ST_GATE_OFF`) are folded so that, with `RUSTLE_PRECISE` unset, the ST-faithful side is active coherently; with `RUSTLE_PRECISE=1`, today's side is active.

### 3.2 Coherent stage-groups (the unit of work)
Stages are **coupled** (permissive junctions ↔ heavy filters ↔ coverage), so they are ported as **groups that converge together** — porting one stage alone diverges *more* (the `ST_SHADOW`→1494 lesson). Dependency order:

| Group | Stage | What "ST-faithful" means here | Coupled to |
|---|---|---|---|
| **G1** | Junction acceptance (`graph_build.rs`, `junction_graph*.rs`, `killed_junctions.rs`) | accept StringTie's permissive junction set (finish `ST_SHADOW` layers 3–6: support/consensus-less acceptance, the `mm`/`nreads_good`/`1.25*junctionthr` logic matched to `count_good_junctions` rlink.cpp:14292–14653) | G3 |
| **G2** | Coverage / flow (`max_flow.rs`, `path_extract.rs`, node coverage) | the coverage the permissive junctions + flow produce; already mostly faithful (CPath/usepath ported) — confirm it matches ST given G1's inputs | G1, G3 |
| **G3** | Prediction filtering (`transcript_filter.rs` print_predcluster analog) | StringTie's exact filter cascade — `included_drop`, `retained_intron`, `isofrac`, readthr, polymerase-runoff — matched to rlink.cpp print_predcluster, acting on G2 coverage | G1 |
| **G4** | Extraction (`parse_trflong_st.rs`) | already mostly ported (canonical `_st`); fold in the remaining "extension admits low-cov ~1.0 sub-paths" gap (§6q–s) | — |

**G1 and G3 converge together**: permissive junctions are only "safe" once ST's downstream filters are matched. Expect the default to regress when G1 lands and recover when G3 lands.

### 3.3 Where precision moves
Everything that makes today's rustle diverge to be *more correct* goes behind `RUSTLE_PRECISE` (or stays in VG): the strict-early junction rejection (strand_zero/mm_negative), the read-through avoidance (locus C), `RUSTLE_MIN_MULTI_INTRON_COV`, the precision-oriented kill rules, etc. None are *deleted* — they are *gated*.

## 4. Methodology: tracer-driven on mini3

Every iteration:
1. Run both tools on the **mini3 fixture** (`bench/mini3/mini3.bam`, 0.09s/run vs ~30s full).
2. `gtf_chain_diff(rustle, st)` → list of rustle-only + ST-only chains.
3. For each divergent chain, run the **first-divergence tracer** (junction_accept → transfrag_collapse → pred_filter_stage stages → pred_kill) to find the exact fork stage + decision (e.g. `AFTER_pairwise_overlap_filter` / `included_drop` / `retained_intron`).
4. Port that decision exactly (within its coherent group), gating today's behavior behind `RUSTLE_PRECISE`.
5. Re-measure rustle-vs-ST on mini3; assert `RUSTLE_PRECISE=1` still byte-matches `4705ab1`.
6. Converge mini3 to 0/0, then full-BAM, then genome-wide.

## 5. Fixture

Promote `/tmp/mini3.bam` to a committed fixture:
- `bench/mini3/mini3.bam` + `.bai` (284 reads, 119 KB, the 3 loci: A 17160875-17175957 / B 36012544-36069321 / C 68848317-68871262).
- `bench/mini3/expected_st.gtf` (StringTie's output on mini3 — the convergence target).
- `bench/mini3/check.sh` — runs rustle on mini3 + `gtf_chain_diff` vs expected_st, prints rustle-only/ST-only counts (current baseline before any work: **3 / 16**).

mini3 reproduces the full-BAM per-locus behavior exactly because the 3 loci are isolated bundles.

## 6. Scope boundaries

**In:** junction acceptance (G1), coverage confirmation (G2), prediction filtering (G3), extraction tail (G4) — i.e. the stages causing the 3-loci divergence and the genome-wide 187/104.

**Out (gated, not removed):** VG-mode features, the precision filters, the +recall — all become `RUSTLE_PRECISE`/VG opt-ins. No new precision work; this project only *relocates* precision behind a flag and *adds* ST-faithful behavior.

## 7. Risks & mitigations

- **Transient default regression** (default goes worse before converging). *Mitigation:* accepted by design; `RUSTLE_PRECISE=1` and commit `4705ab1` are escape hatches; the test suite uses `RUSTLE_PRECISE=1` for the precise-behavior assertions so they don't all break during the flip.
- **Coupling / non-monotonic progress** (a fix in G1 needs G3 to land before it helps). *Mitigation:* work in coherent stage-groups; measure rustle-vs-ST, accept intra-group regression, gate at group boundaries.
- **Long tail of scattered decisions.** *Mitigation:* tracer-driven on mini3 first (cheap), only escalate to full/genome-wide per group.
- **Genome-wide validation OOM** (whole-genome -L ~18 GB). *Mitigation:* per-chrom serial (`gw_run.sh` JOBS=1); ⚠ never `pkill -f rustle` (cwd path).
- **The escape hatch silently drifts.** *Mitigation:* a CI-style invariant test: `RUSTLE_PRECISE=1` on mini3 must byte-match the frozen `4705ab1` mini3 output, asserted every task.

## 8. Success criteria

1. mini3: rustle-vs-ST = 0/0.
2. full GGO_19.bam: rustle-vs-ST = 0/0 (or a documented, characterized residual).
3. genome-wide per-chrom: rustle-vs-ST → 0 (characterized residual acceptable).
4. `RUSTLE_PRECISE=1` byte-matches commit `4705ab1` throughout.
5. The precision/recall improvements remain reachable via `RUSTLE_PRECISE`/VG and are documented.

# Design: `RUSTLE_ST_SHADOW` — coherent end-to-end StringTie-faithful mode

Status: DESIGN (approved 2026-05-28). Next: implementation plan (writing-plans).
Context: see `docs/STRINGTIE_PARITY_FINDINGS.md` (the per-layer divergences this composes).

## 1. Problem & insight

Every StringTie-parity change tried this session regresses **in isolation** (junction acceptance,
read-split, lowintron, multicov, longcov, flow allocation) because Rustle's defaults are a *coherent
system* tuned to make sense of incomplete long-read data. Each Rustle relaxation/strictness is
load-bearing given the others. Aligning one layer toward StringTie breaks the compensation the
others rely on → F1 drops.

**Insight (the project premise):** make every layer StringTie-faithful *together*, in a separate
shadow mode, so the behaviors reinforce instead of fighting. The mid-stack F1 regressions stop
mattering because we never optimize F1 mid-stack — we validate **parity per layer** and only check
F1 once the whole stack converges.

## 2. Goal & non-goals

- **Goal:** bit-exact StringTie parity as a *proof* — when `RUSTLE_ST_SHADOW=1`, Rustle reproduces
  StringTie's output as closely as possible (ideally chain-identical), per-layer validated. The win
  is "Rustle can be StringTie when it chooses," and a precise account of every divergence.
- **Non-goal (now):** changing the default. Default (flag off) stays today's Rustle operating point
  (Sn 96.5 / Pr 90.7 — higher sensitivity, a deliberate point). Not a production `--stringtie-compat`
  mode yet; not a best-F1 hybrid.

## 3. Architecture — dispatch, not duplication

- New predicate `st_shadow()` in `stringtie_parity.rs`, **default OFF** (`RUSTLE_ST_SHADOW=1` to
  enable). Distinct from the always-on partial `stringtie_exact()` (which is the *default* Rustle
  behavior and stays as-is).
- Each layer's existing `stringtie_exact()` / `canonical` checks are OR'd with `st_shadow()`.
- The *new* layer gates (the ones that regress in isolation) are gated on `st_shadow()` ALONE
  (so they never affect the default). Existing experiment flags
  (`RUSTLE_KEEP_MM_NEG`, `RUSTLE_NO_KILLED_SPLIT`, `RUSTLE_LEFTOVER_REDIST`,
  `RUSTLE_ISOFRAC_CHAIN_DEDUP`, …) become *components* — `st_shadow()` implies them.
- No pipeline duplication: same orchestration; ST-faithful logic dispatched per layer via the
  existing `_st` modules (`parse_trflong_st`, `junction_graph_st`, `long_max_flow_st`).

## 4. Layers (build + validate bottom-up)

Order matters: each layer depends on the previous being ST-identical. Each has a parity event
already wired as its gate.

| # | layer | ST-faithful behavior under shadow | code site(s) | parity gate |
|---|---|---|---|---|
| 1 | junction acceptance | accept `mm_negative` / single-read junctions (no reject on `stat.mm<0`) | graph_build.rs:834 (`RUSTLE_KEEP_MM_NEG` logic) | `junction_accept`: Rustle-rej ∩ ST-acc → 0 |
| 2 | graph construction | ST node/edge set (junction_graph_st), no extra low-cov islands | junction_graph_st.rs, graph_build.rs | `bundlenode_list` / `graphnode_list` structural diff |
| 3 | read→transfrag | no kill-split/orphan; ST `update_abundance` segmentation | map_reads.rs (`RUSTLE_NO_KILLED_SPLIT` logic), transfrag_process.rs | `transfrag_pre_depl`: Rustle-only chains → 0 |
| 4 | flow | `long_max_flow_st` nodeflux + leftover-transfrag redistribution | max_flow.rs, path_extract.rs (`RUSTLE_LEFTOVER_REDIST`) | `path_extracted`: extracted-chain set matches |
| 5 | filters | ST isofrac multicov (chain-dedup) + retained-intron on ST coverage | transcript_filter.rs | `pred_kill`: kill set matches by reason |
| 6 | abundance | ST `Cov_Sum` denominator + seed-abundance longcov | transcript_filter.rs:481, path_extract.rs:9113 | GTF cov/TPM/longcov vs ST |

## 5. Per-layer validation (the discipline)

For each layer N, in order:
1. Enable layer N's shadow gate (layers 1..N all on).
2. Run both tools with the relevant `*_PARITY_LOG` step.
3. Drive layer N's parity-diff to **zero divergence on the frozen test set** (full GGO_19 chr19, or
   a fixed locus subset) before enabling layer N+1.
4. A layer is "done" only when its parity event matches ST.

**We do NOT measure F1 mid-stack** (it will look bad — that's expected and is the whole point).
F1 is checked only once all six layers converge. The per-layer parity-diff is the gate; the final
chain-identity-to-`GGO_19.gtf` + F1 is the proof.

Tooling: `bench/transfrag_parity_diff.py` (layer 3), the junction/graph/pred_kill diffs via the
parity-decisions JSONL, `RUSTLE_TRACE_READ` / `ST_TRACE_READ_START` for read-level checks. ST built
`make clean release` (~10s).

## 6. Build order & increments

Bottom-up. **First increment (first plan/session):** scaffold `st_shadow()` + wire Layer 1 (junction
acceptance) + drive `junction_accept` Rustle-rej∩ST-acc → 0. Each subsequent layer is its own
plan/session with its own parity gate. The project converges when the shadow GTF is chain-identical
to `GGO_19.gtf` (or the residual is fully accounted for).

## 7. Safety & exit

- `st_shadow()` default-off → zero risk to the default 96.5/90.7 operating point. Every shadow gate
  is additionally individually inspectable via its existing `RUSTLE_*` env.
- **Abort/re-scope criterion:** if after Layers 1–3 the `path_extracted` set still diverges with all
  three gates parity-matched, the divergence is deeper than these layers model (e.g. read alignment
  ingestion) and the layer decomposition needs revisiting.
- This is explicitly a multi-session effort; each session lands one converged layer.

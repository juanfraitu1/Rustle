# Design: `RUSTLE_ST_SHADOW` Layers 2.5 + 3 — graph-node parity + read→transfrag construction

Status: DESIGN (approved 2026-05-29). Next: implementation plan (writing-plans).
Parent spec: `docs/superpowers/specs/2026-05-28-stringtie-shadow-mode-design.md`.
Findings: `docs/STRINGTIE_PARITY_FINDINGS.md` (§5b, §6b, §6c).

## 1. Problem & insight

Layers 1+2 of the shadow mode made junction *acceptance* match StringTie (mm_negative → 0; the
`good_junc` `mm=-1` demotion removed under shadow). But that *increased* read→transfrag
over-segmentation: `transfrag_pre_depl` Rustle-only chains rose 3157 → **5756** under shadow L1+2 (vs
ST 4383 total). This is the expected mid-stack regression — keeping more junctions alive gives
Rustle's unchanged read→transfrag construction more material to build distinct chains from.

The read-level investigation (2026-05-29) found the 5756 Rustle-only chains are **49% truncations
(subsets) of an ST chain, 21% extensions, 80% single-read** — dominated by Rustle *over-splitting*
reads into partial chains, all using junctions ST already has (zero novel junctions in the truncation
set). Crucially, suppressing the split alone (two existing flags) made parity **worse** (5756 →
6623 / 7203), because **the graph node set itself differs**: even reads Rustle keeps whole produce
different intron coordinates than ST (e.g. a 12-intron Rustle chain vs ST's 11-intron chain over the
same span). Identical transfrag chains are impossible until the underlying graph node boundaries
match.

**Insight:** "Layer 3" decomposes into two coupled mechanisms that must land together:
**Layer 2.5 (graph node boundary parity)** is a hard prerequisite for **Layer 3 (read→transfrag
split + trim)**. Parity is measured only after both land — consistent with the full-StringTie-mimicry
strategy (accept transient regressions; align every coupled layer together).

## 2. Goal & non-goals

- **Goal:** under `RUSTLE_ST_SHADOW=1`, make Rustle's graph node set and read→transfrag construction
  StringTie-faithful so `transfrag_pre_depl` Rustle-only chains drive from 5756 → 0 (or a fully
  accounted residual). Default (flag OFF) stays at today's operating point (Sn 96.5 / Pr 90.7).
- **Non-goal:** changing the default; measuring F1 mid-stack (only at the end of the 6-layer
  project); the downstream flow/filter/abundance layers (4–6).

## 3. Architecture — dispatch, not duplication

- All behavior gated on the existing `crate::stringtie_parity::st_shadow()` (default OFF). With the
  flag off, every code path is byte-identical to today.
- **Layer 2.5 (graph nodes): dispatch-in-place** in `src/rustle/graph_build.rs`. The node-boundary /
  intron-split decision is expected to be a localized rule; branch on `st_shadow()` to match ST's
  node-splitting (where ST introduces a node boundary vs where Rustle does).
- **Layer 3 (read→transfrag): parallel ST-faithful path** in `src/rustle/map_reads.rs`, matching the
  established `_st` convention (`parse_trflong_st`, `junction_graph_st`, `long_max_flow_st`). Add an
  ST-faithful split/pattern routine dispatched under shadow; leave Rustle's default functions intact.
  `map_reads.rs` currently has zero `st_shadow()` references — greenfield wiring.

## 4. Components

### Layer 2.5 — graph node parity
- **Site:** `src/rustle/graph_build.rs` node/bundlenode construction (where introns become node
  boundaries). The exact divergence (extra/shifted split nodes) is characterized in implementation
  Task 2 (the oracle in Task 1 confirms node-parity is sufficient first).
- **ST-faithful behavior:** Rustle's node boundaries and intron-split points match ST's graph so the
  same read yields the same intron coordinates.
- **Gate:** `graphnode_list` / `bundlenode_list` structural diff between the tools → 0 (or accounted).

### Layer 3 — read→transfrag construction
- **Sites:** `src/rustle/map_reads.rs` `split_read_segments` (~1339), `add_or_update_transfrag`
  (~1666), the single-node-fragment drop (~652); reference `collect_read_nodes_exact` (~250).
- **ST-faithful behavior (port of `get_fragment_pattern` rlink.cpp:4688–4708 + `update_abundance`
  4367–4496):**
  1. Split a read's node path ONLY at consecutive non-contiguous nodes where no read junction
     exactly matches the node boundary (`juncs[i]->start==prev_end && juncs[i]->end==curr_start`).
     Keep all contiguous / junction-matched nodes in one pattern. (Replaces Rustle's killed/BADJUNC
     split condition under shadow.)
  2. Port the long-read end-node trim (rlink.cpp:4378–4496): trim spanning end nodes back to
     source/sink so partial reads canonicalize onto the full chain (the `findtrf_in_treepat`
     collapse).
  3. Do not create/drop single-node fragments under shadow — ST never creates the fragment in the
     first place (so the map_reads.rs:652 drop has nothing to drop).
- **Gate:** `transfrag_pre_depl` Rustle-only chains (`bench/transfrag_parity_diff.py`) → 0.

## 5. Data flow

reads → bundle → **graph nodes (Layer 2.5)** → **read node-paths / transfrags (Layer 3)** →
[deferred: flow (4), filter (5), abundance (6)]. Layer 2.5 output (node set) is Layer 3's input;
that ordering is why 2.5 must be correct before Layer 3 can converge.

## 6. Build order & validation (the discipline)

Parity-diffs are the gates (not unit tests — this is graph-level behavior, same as Layers 1–2). The
`st_shadow()` predicate already has unit-test coverage. Plan task order:

1. **Oracle bound (de-risk first).** Verify `graphnode_list`/`bundlenode_list` events are wired on
   both tools (parent spec says they are — confirm). Ingest ST's `graphnode_list` and force Rustle to
   adopt ST's node boundaries (oracle), then re-measure `transfrag_pre_depl` Rustle-only. **Exit
   check:** if it approaches 0, node-parity is the right prerequisite → proceed. If it does NOT, the
   divergence is deeper than nodes+split → replan before any porting (cheap kill).
2. **Layer 2.5 real implementation** → drive `graphnode_list` diff → 0 under shadow.
3. **Layer 3 real implementation** → drive `transfrag_pre_depl` Rustle-only → 0 under shadow.
4. **Validate combined** (both gates ≈ 0) + record in findings doc. F1 NOT measured.

We do NOT measure F1 mid-stack. Each step's parity-diff is its gate; the combined transfrag-chain
identity is the proof of Layers 2.5+3.

## 7. Safety & exit

- `st_shadow()` default-OFF → zero risk to the 96.5/90.7 default; verified by the unchanged-default
  regression check at each implementation step (as in Layers 1–2).
- **Abort/re-scope criterion:** if the Task-1 oracle shows ST node boundaries still don't bring
  Rustle-only near 0, nodes+split is not the full story (e.g. read-ingestion / boundary-snap upstream)
  — revisit the decomposition before investing in the port.
- Transient regression of the live shadow path (e.g. Rustle-only rising during Layer 2.5 before
  Layer 3 lands) is expected and acceptable per the full-mimicry strategy; only the default path is
  guarded.
- Explicitly a multi-session effort; 2.5 and 3 land as a coherent pair.

# Transfrag-provenance layer — SCOPE (not yet implemented)

Goal: trace each **ST_ONLY** transfrag (StringTie constructs it, rustle never
does — surfaced by `transfrag_construction_diff.py`) to the **exact rustle
lifecycle step** that failed to produce it. This is the entry point for the
deferred transfrag-construction rewrite that is the terminal root cause of
the cov gap (`project_cov_estimation_divergence_2026_05_15`).

## Transfrag lifecycle (both tools, verified)

```
reads
  → read→transfrag build
      rustle: map_reads_to_graph[_bundlenodes]  (pipeline.rs ~12431-12453)
              then tag_transfrags_origin_if_missing(.., "read_map")
      ST:     get_read_to_transfrag / CTransfrag build (rlink.cpp ~2913/5141)
  → process_transfrags                            (rustle transfrag_process.rs:1703 ;
                                                    ST rlink.cpp:5769)
      - eliminate_transfrags_under_thr  (rustle :1803; abund<min OR > MAX_TRF_NUMBER)
      - keeptrf grouping/absorb          (rustle :2038+; absorbed → weak=1, NOT deleted)
  → transfrag_define emit  (EXISTING, POST-process: rustle pipeline.rs ~12704;
                            ST rlink.cpp ~16352 — all transfrags, coord-aligned)
  → seeding / get_trf_long / flow
```

Key fact: the existing `transfrag_define` fires **after** process_transfrags.
keeptrf-absorbed transfrags are marked `weak=1` but still emitted → they show
up as MATCH/ABUND_DIVERGE, not ST_ONLY. Therefore an **ST_ONLY** transfrag
(absent from rustle's POST emit) was either:

| class | meaning | rustle site |
|---|---|---|
| `NEVER_CONSTRUCTED` | absent from rustle even *pre*-process | read→transfrag map: junction-correction window, `long_read_min_len`, exon→node mapping in `map_reads_to_graph` |
| `DROPPED_IN_PROCESS` | present pre, gone post | `eliminate_transfrags_under_thr` (abund<min_abundance, or MAX_TRF_NUMBER cap) |

(keeptrf absorb is NOT an ST_ONLY cause — emitted weak=1.)

The current scaffold cannot distinguish these. The provenance layer adds the
one missing observation point + reason tags to do so.

## Provenance layer = 3 additions per tool (all env-gated, parity-only)

1. **`transfrag_define_pre` emit** — identical payload shape to
   `transfrag_define`, emitted *immediately after read→transfrag build,
   before `process_transfrags`*.
   - rustle: pipeline.rs right after `tag_transfrags_origin_if_missing(&mut
     transfrags,"read_map")` (~12454).
   - ST: rlink.cpp right before the `process_transfrags(s,...)` call
     (~16338), same `for t in transfrag[s][b]` loop body as the existing
     post emit.
   - Lets the diff 3-way join PRE-rustle / POST-rustle / POST-ST →
     classify every ST_ONLY as NEVER_CONSTRUCTED vs DROPPED_IN_PROCESS.

2. **Enrich both `transfrag_define` + `_pre` payloads** with provenance
   fields (already available on the structs — no algorithm change):
   - `origin` — rustle `tf.origin_tag` (read_map / guide / guide_edge /
     source_connector / sink_connector / terminal_edge_inject_* / oracle);
     ST: a construction-site classifier (read / guide / addsource /
     addsink / future) set where each `CTransfrag` is built.
   - `read_count` — rustle `tf.read_count`; ST `abundance` (long-read
     proxy, ~1.0/read). Distinguishes single-long-read transfrags (the
     ST_ONLY class is abund≈1 long interior-start spans).
   - `real`, `longread` — already on both structs.

3. **`transfrag_drop` emit** wrapping `eliminate_transfrags_under_thr`
   (rustle :1803) and the ST equivalent: emit (chain, abund, reason ∈
   {under_thr, max_trf_cap}) for every transfrag it removes. Turns
   DROPPED_IN_PROCESS into a reasoned attribution.

## Diff-tool extension (`transfrag_construction_diff.py --pre`)

New 3-way mode: inputs PRE-rustle, POST-rustle, POST-ST logs. For each
ST_ONLY (POST-ST chain absent from POST-rustle):
- in PRE-rustle?  no  → **NEVER_CONSTRUCTED** (+ report ST origin/read_count;
  the read→transfrag mapping diverged — deepest class).
- in PRE-rustle?  yes → **DROPPED_IN_PROCESS**; join `transfrag_drop` by
  chain → reason (under_thr / max_trf_cap), with PRE abund vs ST abund.

Output: per ST_ONLY transfrag a one-line verdict
`chain | NEVER_CONSTRUCTED | st_origin=read_map st_reads=1` or
`chain | DROPPED_IN_PROCESS reason=under_thr pre_abund=1.0 st_abund=1.0`.

## Effort / risk

- Per tool: 1 new emit (reuse existing transfrag_define loop body), 1 payload
  enrichment, 1 drop-wrap. Tool: ~60 lines. Mirrors the existing
  transfrag_define scaffold exactly.
- **Risk: low.** Pure instrumentation, env-gated by the parity-log vars
  already in use; zero algorithm change → rustle default byte-identical
  1746/1948 (verify each step, as with every prior layer). ST changes
  parity-payload-only, local (gpertea upstream, not pushed).
- ~1 implementation session. Decisive output: splits the 9 STRG.15.1
  ST_ONLY transfrags into "read-map never built it" vs "process dropped it
  (reason)", which picks the fix locus (map_reads_to_graph junction/length
  handling vs eliminate-under-thr threshold) for the eventual — still
  multi-session, maximal-blast-radius — construction rewrite.

## Status

SCOPE ONLY. No code changed. The fix itself remains the documented deferred
rewrite; this layer only makes the divergence step-attributable so that
rewrite can be targeted and verified against 1746/1948.

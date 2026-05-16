# transfrag-construction scaffold — cross-tool diff of CONSTRUCTED transfrags

## Why this exists

The rustle<StringTie post-flow coverage gap that kills ~28 missed refs at
the isofrac/pairwise filters was traced (see
`project_cov_estimation_divergence_2026_05_15` memory) to a **transfrag
construction divergence upstream of max_flow**: rustle does not build the
long-range *interior-start* spanning transfrags StringTie builds, so its
capacity network lacks chord edges → fewer augmenting iters → ~0.5× flux →
low cov → isofrac kill. This scaffold makes that divergence directly
observable and diffable, as the entry point for any future
transfrag-construction work.

## Components (all env-gated; default byte-identical 1746/1948)

1. **`transfrag_define` parity event — both tools, construction-time.**
   - rustle: `pipeline.rs` ~12704, `for tf in transfrags.iter()`
     immediately after `process_transfrags(...)` (BEFORE seeding / flow
     depletion). Payload: `abund, longread, n_introns, trflong_seed, weak,
     usepath, n_nodes, guide, introns`.
   - StringTie: `rlink.cpp` ~16352, `for t in transfrag[s][b]` immediately
     after `process_transfrags(...)` (BEFORE `get_trf_long`). Payload:
     `abund, longread, n_introns, n_nodes, introns` (n_nodes added
     2026-05-15 for parity).
   - **Coordinate convention is ALIGNED across tools**: both emit introns
     as `donor_end+1 .. acceptor` where rustle's 0-based acceptor ==
     StringTie's 1-based−1. The emitted strings match directly — unlike raw
     node-coord dumps (seed/keeptr/capedge), which carry a 0/1-based
     off-by-one that fakes "ABSENT". **Always diff at transfrag_define /
     by intron-chain, never by raw node coords.**

2. **`transfrag_construction_diff.py`** — joins the two parity logs by
   (strand, canonical intron-chain); buckets into MATCH / ABUND_DIVERGE /
   ST_ONLY / RUSTLE_ONLY. ST_ONLY sorted by span = the missing
   long-range/interior-start transfrags.

## Running (canonical case: STRG.15.1, NC_073243.2:16598646-16622552)

```bash
RUSTLE_PARITY_LOG=/tmp/tc_r.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 \
  RUSTLE_PARITY_FILTER_RANGE=16598646-16622552 RUSTLE_FLOW_RESIDUAL_SE=1 \
  ./Rustle/target/dev-opt/rustle -L GGO_19.bam -o /tmp/tc_r.gtf
STRINGTIE_PARITY_LOG=/tmp/tc_s.jsonl STRINGTIE_PARITY_FILTER_CHROM=NC_073243.2 \
  STRINGTIE_PARITY_FILTER_RANGE=16598646-16622552 \
  ./stringtie/stringtie -L GGO_19.bam -o /tmp/tc_s.gtf
python3 Rustle/tools/trace_analysis/transfrag_construction_diff.py \
  --rustle /tmp/tc_r.jsonl --stringtie /tmp/tc_s.jsonl \
  --range 16598646-16622552 --min-introns 1 --top 20
```

## Validated result (STRG.15.1 locus, 2026-05-15)

`MATCH=48 ABUND_DIVERGE=7 ST_ONLY=9 RUSTLE_ONLY=5`. The 9 ST_ONLY are the
long-range interior-start spanning transfrags rustle never constructs, e.g.
`16601927-16612061 nI=18`, `16602409-16622552 nI=27`,
`16604715-16629514 nI=23`, `16599290-16628462 nI=30` — exactly the
capacity-chord transfrags whose absence starves STRG.15.1's flux. ABUND_DIVERGE
shows shared chains with construction-time abundance deltas (mostly ±1–2).

## How to use this for the fix (future, multi-session)

1. Pick an ST_ONLY transfrag (e.g. `16601927-16612061`). It is a long read
   (or merged read group) spanning interior→downstream that StringTie's
   `process_transfrags` keeps as a transfrag but rustle's does not.
2. The next instrumentation layer (not yet built): tag each
   `transfrag_define` with its construction PROVENANCE — which read(s) /
   merge step produced it — on both sides, so an ST_ONLY transfrag can be
   traced to the rustle `process_transfrags` step that dropped/merged it
   away. Candidate rustle sites: `process_transfrags` (pipeline.rs),
   transfrag collapse/merge, the long-read→transfrag construction in
   `transfrag_process.rs`.
3. Any construction change has MAXIMAL blast radius (transfrags feed the
   whole pipeline). Gate opt-in, diff with this tool first, then full
   chr19 regression vs 1746/1948 before considering default-on.

## Status

Scaffold only — NO algorithm change. rustle default byte-identical
1746/1948 F1=92.21%. StringTie change is parity-payload-only
(`n_nodes` field), getenv-gated by `STRINGTIE_PARITY_LOG`. The fix itself
(reconciling rustle transfrag construction with StringTie) is NOT attempted
— it is the documented terminal root cause and a multi-session rewrite.

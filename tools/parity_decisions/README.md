# Parity Decisions Scaffold

Cross-tool decision logger for diffing rustle vs StringTie at graph-build time.

## How it works

Both rustle and StringTie are instrumented to emit one JSONL line per decision
to a log file. The two logs can then be diffed by `(step, start, end, strand)`
key to surface divergences.

## Files

| File | Purpose |
|---|---|
| `Rustle/src/rustle/parity_decisions.rs` | Rust emitter (env-gated, BufWriter, locked) |
| `stringtie/parity_decisions.{h,cc}` | C/C++ emitter (line-buffered fprintf) |
| `Rustle/tools/parity_decisions/diff.py` | JSONL diff tool |

## Running

### Rustle

```bash
RUSTLE_PARITY_LOG=/tmp/p_rustle.jsonl \
RUSTLE_PARITY_FILTER_RANGE=22149137-22155355 \
    /scratch/jxi21/Assembler/Rustle/target/release/rustle \
    /scratch/jxi21/Assembler/GGO_19.bam -o /tmp/x.gtf
```

### StringTie

```bash
STRINGTIE_PARITY_LOG=/tmp/p_stringtie.jsonl \
STRINGTIE_PARITY_FILTER_RANGE=22149137-22155355 \
    /scratch/jxi21/Assembler/stringtie/stringtie -L \
    /scratch/jxi21/Assembler/GGO_19.bam -o /tmp/x_st.gtf
```

### Diff

```bash
python3 tools/parity_decisions/diff.py /tmp/p_rustle.jsonl /tmp/p_stringtie.jsonl \
    --step junction_accept --range 22149137-22155355
```

## Env vars (both sides have parallel knobs)

| rustle | stringtie | meaning |
|---|---|---|
| `RUSTLE_PARITY_LOG` | `STRINGTIE_PARITY_LOG` | output JSONL path (enables logging) |
| `RUSTLE_PARITY_FILTER_RANGE` | `STRINGTIE_PARITY_FILTER_RANGE` | `LO-HI` coord range (overlap test) |
| `RUSTLE_PARITY_FILTER_CHROM` | `STRINGTIE_PARITY_FILTER_CHROM` | chromosome name |
| `RUSTLE_PARITY_FILTER_STEPS` | `STRINGTIE_PARITY_FILTER_STEPS` | comma-sep step name whitelist |

## Schema

```jsonl
{"step":"junction_accept","tool":"rustle","start":22140526,"end":22149725,"strand":"-","payload":{"accepted":true,"bundle_strand":"-","jstrand":-1,"mm":2.0,"reason":"ok"}}
```

`step` and `tool` mandatory. `chrom` optional. `start`/`end`/`strand` are conventional;
`payload` is step-specific.

## EM scope note

"EM" has two distinct meanings in this codebase. **Do not conflate them:**

| term | what it is | active by default |
|---|---|---|
| **Native flow estimation** | Both tools' built-in flow optimization (min-cost flow in ST, Rustle's own flow EM). Always runs. | **Yes — always active in both tools** |
| **BundlenodeGraphInferer EM** | Experimental Rust re-implementation of ST's graph/EM pipeline, invoked via `RUSTLE_JUNCTION_GRAPH=1`. | **No — opt-in flag only** |

Parity comparisons run the **native** flow estimation on both sides. The BundlenodeGraphInferer is out of scope until specifically enabled. Pre/post EM snapshots are valuable for debugging when EM is active, but adding those parity events is deferred until the BundlenodeGraphInferer is promoted to default.

## Pipeline stages and parity coverage

```
bundlenode construction
    └─ [bundlenode_list] ← wired ✓
trim-point detection (chi-square) + read→node coverage accumulation
    └─ [graphnode_list]  ← wired ✓ (2026-05-26)
flow graph EM (native, always active)
    └─ [transfrag_pre_depl] ← wired ✓
depletion / isofrac filtering
path assembly (get_trf_long / checktrf)
    └─ [path_extracted] ← wired ✓
pred filtering
    └─ [pred_intron_low] ← wired ✓
GTF emit
    └─ [path_emit] ← wired ✓
```

### STRG.334.2 corrected diagnosis (2026-05-26, fully traced)

**All prior diagnoses were wrong.** Full stage-by-stage trace shows:

- `graphnode_list`: **identical** — both tools have the 60bp microexon (52980399-52980458, cov=96884) as a graphnode. Rustle has 1 extra 4bp low-cov island node (52966246-52966249) ST merges away; otherwise identical.
- `transfrag_pre_depl`: both tools have transfrags through the microexon (`52971296-52980398,52980459-52981544`) AND skip-path transfrags (`52971296-52981544`).
- `path_extracted`: both tools emit 13-exon paths through the microexon. Both tools emit some 13-exon paths with the correct last intron `53042103-53043292` (ST: 12 such paths; Rustle: 10, including 4 with nexons=13).
- **Actual divergence point: pred filtering — confirmed (2026-05-26 path_extracted diff).** All 4 Rustle paths with last intron `53042103-53043292` (nexons=12/13, lc=2–19) are killed. Kill mechanism traced via `pred_kill` events:
  - Primary kill: `stage=isofrac`, `isofraclong=0.01` (cov ~10–20 / usedcov ~2516 ≈ 0.4–0.8% < 1% threshold)
  - Secondary kills: `stage=pairwise`, `reason=retained_intron` or `included_drop` (minor-start-variant copies)
  - Survivor: dominant path `(52954348, 53043514)` nexons=13, cov=365, uses last intron `53042103-53043336`
- ST path_extracted shows 5 such paths (same coords), but ST's pred-filter keeps them. ST path_emit not instrumented in these logs, but STRG.334.2 is present in StringTie's final GTF.

**The STRG.334.2 miss is confirmed as an isofrac divergence**: Rustle's 1% isofrac gate kills the minority `53043292`-acceptor paths (cov < 1% of the dominant `53043336`-acceptor cluster). StringTie either applies a lower effective threshold or uses per-cluster (not global) isofrac normalization for these paths.

### Full GGO_19 chr19 graphnode_list diff (2026-05-26)

```
3351 common bundles (3405 Rustle / 3430 ST)
  2865 structurally identical (node coords + cov within 0.5 tolerance)
   486 cov-only differences (minor fp drift)
   116 structural n_nodes mismatches
     109: Rustle has MORE nodes (extra low-cov islands ST merges)
       7: ST has MORE nodes (ST's chi-square trim-point detection fires; Rustle doesn't split)
54 only-Rustle / 79 only-ST (bundle boundary divergences — same loci, different start/end)
```

The 7 cases where ST has more graphnodes were investigated (2026-05-26) and are **not MISS_R sources**.
All 7 splits fall at terminal exon 3'/5' endpoints, not at intron junctions. Since gffcompare's
intron-chain matching ignores outermost exon coordinates, Rustle's slightly different terminal
exon length is invisible to the benchmark: 33/34 reference transcripts at these loci are matched
by Rustle. The 1 MISS_R (STRG.267.3 at 38791833-38890087) is a path-enumeration miss unrelated
to the graphnode split. The largest ST-extra case: +3 nodes at 44054883-44094254 (-).

## Wired-up steps

| step | rustle site | stringtie site |
|---|---|---|
| `bundlenode_list` | `pipeline.rs` — top of `process_graph` closure (~line 12789) | `rlink.cpp` — top of `create_graph()` (~line 3230) |
| `graphnode_list` | `pipeline.rs` — after `map_reads_to_graph`, before `process_transfrags` (~line 13700) | `rlink.cpp` — in `build_graphs`, after `process_refguides`, before `process_transfrags` (~line 15401) |
| `junction_accept` | `graph_build.rs::filter_junctions_for_bundle` | `rlink.cpp:14441` post-checkfeat loop |
| `bundle_define` | `pipeline.rs` | `rlink.cpp:15557` |
| `transfrag_define` | `pipeline.rs` after `process_transfrags` | `rlink.cpp:16003` |
| `transfrag_seed` | `pipeline.rs` — same loop as `transfrag_define`, `tf.trflong_seed==true` | `rlink.cpp` — after each `trflong.Add()` in `process_transfrags` |
| `seed_reject` | `path_extract.rs` — `reject_seed!` macro at early-exit gates in `extract_transcripts` | (rustle-only for now) |
| `path_extracted` | `path_extract.rs` — before `out.push(Transcript{..})` in `extract_transcripts` | `rlink.cpp` — after each `parse_trflong` block in `get_trf_long` and `get_trf_long_mix` |
| `pred_filter_stage` | `transcript_filter.rs` | `rlink.cpp:18609,18613,...` |
| `pred_kill` (pairwise) | `transcript_filter.rs` — `kill!` macro in `pairwise_overlap_filter_with_summary` | `rlink.cpp` — short_first_last_exon (18359), retained_intron (18369), included_drop (18420) — **wired 2026-05-26** |
| `pred_kill` (isofrac) | `transcript_filter.rs:2550` — isofrac_with_summary loop | `rlink.cpp` — longunder branch at line 18620 — **wired 2026-05-26** |
| `path_emit` | `gtf.rs:209` | `rlink.cpp:19681` |

## Adding new decision points

1. **Rustle**: `crate::parity_decisions::emit("step_name", chrom, start, end, strand, &payload_json)`
2. **StringTie**: `pd_emit("step_name", chrom, start, end, strand, payload_json)` (after `#include "parity_decisions.h"`)
3. Format payload as JSON object body (no surrounding braces). Use `format!` in Rust, `snprintf` in C.

## Suggested next steps to wire

- `pred_filter_stage` — per-surviving-transcript events at each filter stage (Rustle already emits these; ST side not yet wired; needed to diff which paths survive the full pred-filter pipeline)
- `edge_create` — parent → child with abundance (junction edges only)

## How to use path_extracted vs path_emit

`path_extracted` fires when the path is assembled, before predcluster filtering.
`path_emit` fires when the path is written to GTF, after all filtering.

To classify a stringtie-only miss:
1. Check `path_extracted` — if present in stringtie but absent in rustle, the divergence is in assembly (path-building) not filtering.
2. Check `path_emit` — if absent in both but the ref produces it, the transcript was lost at filtering on both sides.
3. Check `pred_kill` on the rustle side with `stage:"pairwise"` — shows what pairwise_overlap_filter eliminated.

```bash
# Example: trace a specific intron chain
python3 tools/parity_decisions/diff.py /tmp/r.jsonl /tmp/st.jsonl \
    --step path_extracted --range 22140000-22155000
```

## Known divergences surfaced

### Junction-acceptor indexing off-by-one — FIXED (2026-04-30)

Rustle stores 0-based half-open exon `[s, e)`. StringTie uses 1-based
inclusive `[s, e]`. The donor numerically coincides (rustle 0-based
exclusive-end == StringTie 1-based inclusive-last-base). The acceptor
differs by 1 (rustle 0-based first-base-of-next-exon vs StringTie
1-based first-base-of-next-exon = rustle + 1).

Fix: in `graph_build.rs::filter_junctions_for_bundle`, emit
`j.acceptor + 1` to the parity log. Internal storage unchanged.

After fix: 22 common keys, 13 stringtie-only, 15 rustle-only, 22 common
with payload differences (different field names between sides; not
errors).

### Junction `mm`-sentinel kill — INVESTIGATED, fix-in-isolation regresses F1 (2026-05-01)

For junction `22140526-22150015 (-)`:
```
rustle:    accepted=False, mm=-1.0, reason=mm_negative
stringtie: accepted=True,  mm=1.0
```

Traced to `killed_junctions.rs:697` (`higherr_low_support_bad`) and a
parallel path at line 1194. Rustle preemptively sets `mm=-1` as a
"deferred" mark intending `apply_bad_mm_neg_stage` to re-evaluate, but
that stage skips `mm<0` entries → mark is permanent. StringTie's actual
pattern (rlink.cpp:14914): `if !good_junc(...) jd.mm=-1` — set ONLY on
rejection, not preemptively.

**Removing the preemptive mark recovers the target junction but
catastrophically regresses overall F1:**

| variant | Sn | Pr | F1 | matches |
|---|---|---|---|---|
| baseline (with preemptive mark) | 92.1 | 86.3 | **89.11** | 1693 |
| remove mark | 81.0 | 48.4 | 60.59 | 1489 |
| `RUSTLE_RESCUE_LOWSUPPORT_BAD=1` | 82.4 | 51.1 | 63.08 | 1516 |

The mark is doing necessary FP-suppression work because rustle's bpcov
"rescue" check at `killed_junctions.rs:1188` is too loose (fires at
virtually every splice site). Proper fix needs a stricter bpcov check
matching StringTie's `good_junc` step 6 (rlink.cpp:14823) — with
`mult = 1/ERROR_PERC^2 = 100` and asymmetric donor/acceptor ratios.

**Partial port (2026-05-01):** added bw=5 window-based bpcov-noise check
to the `higherr_low_support_bad` mark decision in
`killed_junctions.rs::good_junc_stats`. F1 preserved at 89.11.
Bypass via `RUSTLE_HIGHERR_NO_BPCOV_GUARD=1`. The same change on the
`cjunctions` parallel path (line 1196-1238) regresses F1 to 60.75 —
that path runs at a later pipeline stage with different `leftsupport`
values, requiring different thresholds. Committed for the safe path
only; real fix requires consolidating the two parallel paths.

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

## Wired-up steps (initial scaffold)

| step | rustle site | stringtie site |
|---|---|---|
| `junction_accept` | `graph_build.rs::filter_junctions_for_bundle` | `rlink.cpp` post-checkfeat loop |

## Adding new decision points

1. **Rustle**: `crate::parity_decisions::emit("step_name", chrom, start, end, strand, &payload_json)`
2. **StringTie**: `pd_emit("step_name", chrom, start, end, strand, payload_json)` (after `#include "parity_decisions.h"`)
3. Format payload as JSON object body (no surrounding braces). Use `format!` in Rust, `snprintf` in C.

## Suggested next steps to wire

- `node_create` — when graph builder creates a new node (start, end, hardstart/hardend flags)
- `edge_create` — parent → child with abundance
- `transfrag_seed` — when a long-read becomes a trflong_seed (abund, node list)
- `seed_reject` — when a seed is dropped pre-flow (reason)
- `path_emit` — when a final path is output (cov, longcov, exon list)

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

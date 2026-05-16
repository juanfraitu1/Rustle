# Read→transfrag divergence — focused scope (the fix locus)

The provenance layer proved all 9 STRG.15.1 ST_ONLY transfrags are
**NEVER_CONSTRUCTED** (origin `read_long`, not dropped by
`eliminate_transfrags_under_thr` / keeptrf). The divergence is in the
read→transfrag mapping. This doc localizes it precisely.

## The mechanism (design-level divergence)

**StringTie** (`get_read_to_transfrag` rlink.cpp:5141 → `get_read_pattern`
→ `CTransfrag(nodes,trpat,...)` / long `longtr` ~2928/2974): builds **one
full-length transfrag per long read's pattern**. It does NOT fragment a
long read's node chain at killed junctions or coverage valleys. (The
"split pairs" code paths are paired-end short-read group splitting, not
long-read chain fragmentation.)

**Rustle** (`map_reads_to_graph_bundlenodes` map_reads.rs:774) fragments a
long read into multiple sub-transfrags via TWO rustle-specific mechanisms
StringTie has no analogue for:

1. **`split_read_segments`** (map_reads.rs:1048) — splits a read's
   node chain at:
   - **killed junctions** (the "V99" behavior: segments right of a killed
     junction become `killed_junction_orphan`)
   - unitig source/sink boundaries
   plus single-node fragments from such reads are then *skipped*
   (map_reads.rs:920-922).
2. **`split_chimeric_transfrags`** (called map_reads.rs:947) — further
   splits a transfrag at coverage valleys.

Net: a long **interior-start spanning** read whose chain crosses a killed
junction or a coverage valley is fragmented; rustle never forms its
full-length transfrag → it is ST_ONLY / NEVER_CONSTRUCTED → its
capacity-chord edge is missing → STRG.15.1-class flux/cov shortfall →
isofrac kill of the ~28 cov-gated refs.

This is consistent with the option-2 path-build finding (97% pphash
divergence) and the cov chain — all roll up to: **rustle fragments long
reads at construction; StringTie keeps them whole.**

## Why the fragmentation exists (do not naively delete it)

`split_read_segments` killed-junction splitting is "V99" — it was added to
stop killed-junction reads from forming spurious chains; the single-node
skip (920-922) suppresses ~23 spurious source/sink transfrags per locus.
`split_chimeric_transfrags` removes genuine chimeric artifacts. These
protect precision. Removing them wholesale will regress (same risk class as
every targeted lever this arc). The fix must keep the protections for true
artifacts while NOT fragmenting reads whose full chain is graph-valid.

## Fix options (increasing risk)

A. **Targeted, opt-in probe (recommended first):** in
   `split_read_segments`, do NOT split at a killed junction when the read's
   *entire* node chain consists of accepted graph edges AND the read is
   long AND spans ≥N introns (i.e. the "killed" junction is actually a
   valid graph edge for this read — the kill was a different-context
   decision). Env-gated default-off. Measure: ST_ONLY count on STRG.15.1
   (via `transfrag_construction_diff.py --pre`), then full chr19 F1 vs
   1746/1948.
B. Gate `split_chimeric_transfrags` to skip long reads whose full chain is
   graph-valid (coverage valley inside a real long isoform is not a
   chimera).
C. Full reconciliation of rustle read→transfrag with StringTie
   `get_read_pattern` (one transfrag per long-read pattern). Maximal blast
   radius — the deferred rewrite.

## Recommended next experiment (bounded, ~1 session)

Implement option A as `RUSTLE_NO_SPLIT_VALID_LONG_CHAIN=1` (default off):
skip killed-junction split when `is_long && read.exons.len()>=K &&
every junction of the read maps to an accepted graph edge`. Then:
1. `transfrag_construction_diff.py --pre` on STRG.15.1 → expect the 9
   ST_ONLY to drop toward 0 (NEVER_CONSTRUCTED → constructed).
2. Full chr19 GGO: matching/query vs 1746/1948. Sweep K.
3. If F1-positive and STRG.15.1 recovers as `=` → candidate default-on,
   precision-checked. If F1-negative (expected risk) → keep opt-in,
   document, and the construction rewrite (option C) remains the only path.

Decisive instrumentation already in place: `transfrag_construction_diff.py
--pre` (NEVER_CONSTRUCTED count is the direct readout), plus the 1746/1948
regression gate. No code changed by this scope doc.

## Status

SCOPE ONLY. The two fragmentation sites (`split_read_segments`
map_reads.rs:1048; `split_chimeric_transfrags` map_reads.rs:947) are the
exact fix locus. Option A is the smallest testable step; the full
read→transfrag reconciliation stays the deferred maximal-blast-radius
rewrite. Practical cov-gated ceiling unchanged: **1746/1948 F1=92.21%**.

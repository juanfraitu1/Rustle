# Rustle Multi-Copy Meeting Prep — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a glass-box live-demo pipeline for tomorrow's advisor meeting that mechanically generates 5 visuals from real Rustle runs on (Tier 1) `synthetic_family` and (Tier 2) the chr19 GOLGA6L paralog triple in `GGO_19.bam`.

**Architecture:** Two new Rust dumps emit splice-graph edges and per-iteration flow state. Five Python renderers (graphviz + matplotlib + networkx) consume the existing + new TSVs to produce PNGs. A shell orchestrator runs everything end-to-end on both tiers. No source-code citations are used in the meeting; every claim is backed by a TSV the advisor can inspect.

**Tech stack:** Rust 1.x (existing Cargo project), Python 3 (matplotlib, graphviz, networkx, pysam), bash for orchestration.

**Spec:** `docs/superpowers/specs/2026-04-26-rustle-multicopy-meeting-prep-design.md`

**Time budget:** ~12-14 hours of work across today + early tomorrow. Drop priority order if compressed: V5 → V4 → Tier 3 closing slide. V1+V2+V3 are non-negotiable.

---

## File Structure

**New files:**
- `src/rustle/parity_graph_edges_dump.rs` — graph-edges dump emitter (mirrors `parity_junction_dump.rs` pattern)
- `src/rustle/parity_flow_iter_dump.rs` — per-iteration flow state dump
- `tools/demo/__init__.py` — package marker
- `tools/demo/common.py` — shared TSV reader, graphviz helpers, color palette
- `tools/demo/render_splice_graph.py` — V1: single splice graph PNG from nodes+edges TSVs
- `tools/demo/compare_graphs_isomorphism.py` — V2: 3-way isomorphism + side-by-side render
- `tools/demo/render_flow_iterations.py` — V3: per-iteration grid PNG
- `tools/demo/render_multimapper_pinning.py` — V4: read pinning before/after EM
- `tools/demo/render_phase3_vg_sketch.py` — V5: variation graph conceptual diagram
- `tools/demo/identify_locus_bundle.py` — find the bundle ID covering chr19 GOLGA6L triple
- `tools/demo/run_demo_synthetic.sh` — Tier 1 orchestrator
- `tools/demo/run_demo_chr19_triple.sh` — Tier 2 orchestrator
- `tools/demo/runbook.md` — meeting-time terminal recipe with timing notes
- `tests/regression/test_parity_graph_edges_dump.rs` — Rust integration test
- `tests/regression/test_parity_flow_iter_dump.rs` — Rust integration test
- `tools/demo/tests/test_render_splice_graph.py` — Python smoke test
- `tools/demo/tests/test_compare_graphs_isomorphism.py` — Python isomorphism test
- `tools/demo/tests/test_render_flow_iterations.py` — Python smoke test

**Files modified:**
- `src/rustle/lib.rs` — add `mod parity_graph_edges_dump;` + `mod parity_flow_iter_dump;`
- `src/rustle/path_extract.rs:5435` — call into graph_edges_dump alongside existing `RUSTLE_PARITY_GRAPH_TSV`
- `src/rustle/max_flow.rs:1243` (and parallel sites at 1333, 2064, 2138) — call into flow_iter_dump inside Edmonds-Karp loops

**Outputs (artifacts produced by the run scripts):**
- `tools/demo/out/synthetic/{nodes,edges,flow,reads,transcripts}.tsv`, `tools/demo/out/synthetic/*.png`
- `tools/demo/out/chr19_triple/{nodes,edges,flow,reads,transcripts}.tsv`, `tools/demo/out/chr19_triple/*.png`

---

## Phase 1 — New Rust dumps

These are the prerequisites for V1, V2, V3. Sequential within the phase: 1 then 2.

### Task 1: Graph-edges dump (`RUSTLE_PARITY_GRAPH_EDGES_TSV`)

**Files:**
- Create: `src/rustle/parity_graph_edges_dump.rs`
- Modify: `src/rustle/lib.rs` (add module declaration)
- Modify: `src/rustle/path_extract.rs` line 5435 (call site, after the existing `RUSTLE_PARITY_GRAPH_TSV` block)
- Test: `tests/regression/test_parity_graph_edges_dump.rs`

**Reference pattern:** The existing `RUSTLE_PARITY_GRAPH_TSV` emitter at `path_extract.rs:5393-5435` is the template. It uses `OnceLock<Mutex<Option<File>>>` for thread-safe file handling and a `OnceLock<Mutex<bool>>` for header-once tracking.

#### Step 1.1 — Read the existing `Graph` struct to confirm edge representation

Run: `grep -n -A 5 "pub struct Graph " /mnt/c/Users/jfris/Desktop/Rustle/src/rustle/graph.rs | head -20`

Expected: confirms `Graph` has `nodes: Vec<GraphNode>`, with each `GraphNode` having `children: SmallBitset` and `parents: SmallBitset` (graph.rs:126-127). Edges are implicit via the children bitsets.

Run: `grep -n -B 1 -A 8 "fn build_capacity\|let.*capacity.*Vec" /mnt/c/Users/jfris/Desktop/Rustle/src/rustle/max_flow.rs | head -30`

Expected: locates how edge capacities are computed in max_flow. Note the function/formula. (Capacity is likely `min(node.coverage, transfrag_pattern_count)` or similar — record the exact expression used at the time the dump fires.)

For the dump, we will emit per-edge: `from_node`, `to_node`, `from_start`, `from_end`, `to_start`, `to_end`, `edge_kind` (collinear if `to_start == from_end + 1`, else junction), `from_coverage`, `to_coverage`, `bottleneck_coverage` (= `min(from_cov, to_cov)`). We do NOT need the actual flow capacity matrix here — bottleneck coverage is the right pedagogical proxy and can be computed from existing GraphNode fields. The advisor will not see "max-flow capacity" as a separate concept; he sees coverage at each end of each edge.

- [ ] **Step 1.2 — Write failing integration test**

```rust
// tests/regression/test_parity_graph_edges_dump.rs
use std::process::Command;
use std::path::Path;

#[test]
fn graph_edges_dump_produces_expected_schema() {
    let bam = "test_data/synthetic_family/reads_sorted.bam";
    let out_gtf = "/tmp/test_graph_edges_dump.gtf";
    let edges_tsv = "/tmp/test_graph_edges_dump.tsv";
    let _ = std::fs::remove_file(edges_tsv);

    let status = Command::new(env!("CARGO_BIN_EXE_rustle"))
        .args(["-L", "-o", out_gtf, bam])
        .env("RUSTLE_PARITY_GRAPH_EDGES_TSV", edges_tsv)
        .status()
        .expect("rustle run failed");
    assert!(status.success(), "rustle exited non-zero");

    assert!(Path::new(edges_tsv).exists(), "edges TSV not created");
    let content = std::fs::read_to_string(edges_tsv).expect("read edges TSV");
    let mut lines = content.lines();

    let header = lines.next().expect("missing header");
    let expected_cols = "source\tchrom\tbdstart\tbdend\tstrand\tfrom_idx\tto_idx\tfrom_start\tfrom_end\tto_start\tto_end\tedge_kind\tfrom_cov\tto_cov\tbottleneck_cov";
    assert_eq!(header, expected_cols, "header mismatch");

    let body_count = lines.filter(|l| !l.trim().is_empty()).count();
    assert!(body_count > 0, "no edge rows written");
}
```

- [ ] **Step 1.3 — Verify test fails**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && cargo test --release --test test_parity_graph_edges_dump 2>&1 | tail -20`

Expected: FAIL — file not created (env var not yet observed by code).

- [ ] **Step 1.4 — Create `parity_graph_edges_dump.rs`**

```rust
// src/rustle/parity_graph_edges_dump.rs
use std::io::Write;
use std::sync::{Mutex, OnceLock};

use crate::graph::Graph;

static WR: OnceLock<Mutex<Option<std::fs::File>>> = OnceLock::new();
static HDR: OnceLock<Mutex<bool>> = OnceLock::new();

/// Dumps graph edges to RUSTLE_PARITY_GRAPH_EDGES_TSV when set.
/// Mirrors the path_extract.rs:5393 pattern. Emits per-edge rows usable
/// alongside RUSTLE_PARITY_GRAPH_TSV for splice-graph rendering.
pub fn emit(
    graph: &Graph,
    bundle_chrom: &str,
    bundle_strand: &str,
    bundle_id: &str,
) {
    let path = match std::env::var("RUSTLE_PARITY_GRAPH_EDGES_TSV") {
        Ok(p) if !p.is_empty() => p,
        _ => return,
    };

    let wr = WR.get_or_init(|| {
        Mutex::new(
            std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&path)
                .ok(),
        )
    });
    let hdr = HDR.get_or_init(|| Mutex::new(false));

    let (bs, be) = bundle_id
        .split_once(':')
        .and_then(|(_, rng)| rng.split_once('-'))
        .map(|(s, e)| (
            s.trim().parse::<u64>().unwrap_or(0),
            e.trim().parse::<u64>().unwrap_or(0),
        ))
        .unwrap_or((0, 0));

    let (Ok(mut f_opt), Ok(mut hdr_w)) = (wr.lock(), hdr.lock()) else { return };
    let Some(f) = f_opt.as_mut() else { return };

    if !*hdr_w {
        let _ = writeln!(
            f,
            "source\tchrom\tbdstart\tbdend\tstrand\tfrom_idx\tto_idx\tfrom_start\tfrom_end\tto_start\tto_end\tedge_kind\tfrom_cov\tto_cov\tbottleneck_cov"
        );
        *hdr_w = true;
    }

    for (from_idx, n) in graph.nodes.iter().enumerate() {
        if from_idx == graph.source_id || from_idx == graph.sink_id {
            continue;
        }
        for to_idx in n.children.iter() {
            if to_idx == graph.source_id || to_idx == graph.sink_id {
                continue;
            }
            let to = &graph.nodes[to_idx];
            let kind = if to.start == n.end + 1 {
                "collinear"
            } else {
                "junction"
            };
            let bottleneck = n.coverage.min(to.coverage);
            let _ = writeln!(
                f,
                "rustle\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}",
                bundle_chrom, bs, be, bundle_strand,
                from_idx, to_idx,
                n.start, n.end, to.start, to.end,
                kind,
                n.coverage, to.coverage, bottleneck,
            );
        }
    }
}
```

> **Note on `n.children.iter()`:** verify by reading `src/rustle/bitset.rs` that `SmallBitset` has an `iter()` returning set indices. If the API is named differently (e.g. `iter_set()`, `bits()`, `ones()`), adjust accordingly. If no iterator exists, fall back to `for to_idx in 0..graph.nodes.len() { if n.children.contains(to_idx) { ... } }`.

- [ ] **Step 1.5 — Register module in `lib.rs`**

```rust
// src/rustle/lib.rs (add alphabetically with other parity_* modules)
pub mod parity_graph_edges_dump;
```

- [ ] **Step 1.6 — Wire call site in `path_extract.rs`**

Locate the closing brace of the `RUSTLE_PARITY_GRAPH_TSV` block (around line 5435), and immediately after it (before the `RUSTLE_PARITY_TF_TSV` block at 5436), add:

```rust
        crate::parity_graph_edges_dump::emit(
            graph,
            bundle_chrom,
            bundle_strand,
            bundle_id,
        );
```

- [ ] **Step 1.7 — Verify test passes**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && cargo build --release 2>&1 | tail -10 && cargo test --release --test test_parity_graph_edges_dump 2>&1 | tail -20`

Expected: build green, test PASS.

If the build fails on the bitset iter — check the actual API and adjust per the note in Step 1.4.

- [ ] **Step 1.8 — Smoke check the output manually**

Run: `RUSTLE_PARITY_GRAPH_EDGES_TSV=/tmp/edges.tsv ./target/release/rustle -L -o /tmp/synth.gtf test_data/synthetic_family/reads_sorted.bam && head -10 /tmp/edges.tsv && wc -l /tmp/edges.tsv`

Expected: header row + ~6-12 edge rows for the 2-bundle synthetic dataset (3-4 nodes each with 1-3 outgoing edges).

- [ ] **Step 1.9 — Commit**

```bash
git add src/rustle/parity_graph_edges_dump.rs src/rustle/lib.rs src/rustle/path_extract.rs tests/regression/test_parity_graph_edges_dump.rs
git commit -m "feat(parity): add RUSTLE_PARITY_GRAPH_EDGES_TSV dump for splice-graph edges"
```

---

### Task 2: Flow per-iteration dump (`RUSTLE_FLOW_ITER_TSV`)

**Files:**
- Create: `src/rustle/parity_flow_iter_dump.rs`
- Modify: `src/rustle/lib.rs` (add module declaration)
- Modify: `src/rustle/max_flow.rs` lines 1243, 1333, 2064, 2138 (4 augmenting-path call sites)
- Test: `tests/regression/test_parity_flow_iter_dump.rs`

#### Step 2.1 — Locate the data we'll emit per iteration

Run: `sed -n '1230,1290p' /mnt/c/Users/jfris/Desktop/Rustle/src/rustle/max_flow.rs`

Expected: the loop reads `pred` (predecessor array from BFS), computes the augmenting path from sink back to source, finds the bottleneck capacity, and updates `flow_mat`. We emit *after* the path is determined and before the flow matrix is mutated, so we can record the path and bottleneck cleanly.

We need:
- `iter_idx` — counter for this bundle's max-flow run
- `bundle_id` — passed in from the caller
- `path` — node sequence (source → ... → sink)
- `bottleneck` — min remaining capacity along the path
- `total_flow_after` — sum of flow already pushed including this iteration

- [ ] **Step 2.2 — Write failing integration test**

```rust
// tests/regression/test_parity_flow_iter_dump.rs
use std::process::Command;
use std::path::Path;

#[test]
fn flow_iter_dump_emits_monotone_iterations() {
    let bam = "test_data/synthetic_family/reads_sorted.bam";
    let out_gtf = "/tmp/test_flow_iter_dump.gtf";
    let flow_tsv = "/tmp/test_flow_iter_dump.tsv";
    let _ = std::fs::remove_file(flow_tsv);

    let status = Command::new(env!("CARGO_BIN_EXE_rustle"))
        .args(["-L", "-o", out_gtf, bam])
        .env("RUSTLE_FLOW_ITER_TSV", flow_tsv)
        .status()
        .expect("rustle run failed");
    assert!(status.success());

    assert!(Path::new(flow_tsv).exists());
    let content = std::fs::read_to_string(flow_tsv).expect("read flow TSV");
    let mut lines = content.lines();
    let header = lines.next().unwrap();
    assert_eq!(
        header,
        "source\tchrom\tbdstart\tbdend\tstrand\tbundle_run\titer_idx\tpath_len\tpath_nodes\tbottleneck\ttotal_flow_after"
    );

    let row_count = lines.filter(|l| !l.trim().is_empty()).count();
    assert!(row_count >= 2, "expected ≥2 flow iterations, got {}", row_count);
}
```

- [ ] **Step 2.3 — Verify test fails**

Run: `cargo test --release --test test_parity_flow_iter_dump 2>&1 | tail -10`

Expected: FAIL — file not created.

- [ ] **Step 2.4 — Create `parity_flow_iter_dump.rs`**

```rust
// src/rustle/parity_flow_iter_dump.rs
use std::io::Write;
use std::sync::{Mutex, OnceLock};

static WR: OnceLock<Mutex<Option<std::fs::File>>> = OnceLock::new();
static HDR: OnceLock<Mutex<bool>> = OnceLock::new();
static RUN: OnceLock<Mutex<u64>> = OnceLock::new();

/// Per-iteration emit: path, bottleneck, total flow.
/// Called once per augmenting path found by Edmonds-Karp BFS, after
/// the path is determined but before flow_mat is mutated.
#[allow(clippy::too_many_arguments)]
pub fn emit(
    bundle_chrom: &str,
    bundle_strand: &str,
    bundle_id: &str,
    bundle_run: u64,
    iter_idx: usize,
    path: &[usize],
    bottleneck: f64,
    total_flow_after: f64,
) {
    let outpath = match std::env::var("RUSTLE_FLOW_ITER_TSV") {
        Ok(p) if !p.is_empty() => p,
        _ => return,
    };

    let wr = WR.get_or_init(|| {
        Mutex::new(
            std::fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(&outpath)
                .ok(),
        )
    });
    let hdr = HDR.get_or_init(|| Mutex::new(false));

    let (bs, be) = bundle_id
        .split_once(':')
        .and_then(|(_, rng)| rng.split_once('-'))
        .map(|(s, e)| (
            s.trim().parse::<u64>().unwrap_or(0),
            e.trim().parse::<u64>().unwrap_or(0),
        ))
        .unwrap_or((0, 0));

    let (Ok(mut f_opt), Ok(mut hdr_w)) = (wr.lock(), hdr.lock()) else { return };
    let Some(f) = f_opt.as_mut() else { return };

    if !*hdr_w {
        let _ = writeln!(
            f,
            "source\tchrom\tbdstart\tbdend\tstrand\tbundle_run\titer_idx\tpath_len\tpath_nodes\tbottleneck\ttotal_flow_after"
        );
        *hdr_w = true;
    }

    let path_str = path
        .iter()
        .map(|n| n.to_string())
        .collect::<Vec<_>>()
        .join(",");

    let _ = writeln!(
        f,
        "rustle\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}",
        bundle_chrom, bs, be, bundle_strand,
        bundle_run, iter_idx, path.len(), path_str,
        bottleneck, total_flow_after,
    );
}

/// Allocate a fresh per-bundle run id so multiple Edmonds-Karp invocations
/// against the same bundle can be told apart in the TSV.
pub fn next_run_id() -> u64 {
    let m = RUN.get_or_init(|| Mutex::new(0));
    let Ok(mut g) = m.lock() else { return 0 };
    *g += 1;
    *g
}
```

- [ ] **Step 2.5 — Register module in `lib.rs`**

```rust
// src/rustle/lib.rs
pub mod parity_flow_iter_dump;
```

- [ ] **Step 2.6 — Read max_flow.rs around 1243 to identify caller context**

Run: `sed -n '1200,1260p' /mnt/c/Users/jfris/Desktop/Rustle/src/rustle/max_flow.rs`

Expected: identify the function name containing line 1243 (likely `weight_bundle` / `solve_flow` / similar), and what variables are in scope at that point: `n` (node count), `capacity`, `flow_mat`, `link`, `pred`, `bundle_chrom`, `bundle_strand`, `bundle_id` (or similar — record exact names).

- [ ] **Step 2.7 — Wire emission inside the augmenting-path loop**

For EACH of the 4 sites (lines 1243, 1333, 2064, 2138), the structure is:

```rust
// before the while loop, hoist a run_id and an iter counter:
let __flow_run_id = crate::parity_flow_iter_dump::next_run_id();
let mut __flow_iter = 0usize;
let mut __flow_total = 0f64;

while bfs_augmenting_path(/* ... */) {
    // existing code that walks pred from sink to source and finds bottleneck
    // ...

    // After bottleneck `b` is determined and BEFORE flow_mat is updated, insert:
    {
        let mut __path: Vec<usize> = Vec::new();
        let mut __cur = sink_id; // or whatever the sink variable is named
        while __cur != usize::MAX {
            __path.push(__cur);
            __cur = pred[__cur];
        }
        __path.reverse();
        __flow_iter += 1;
        __flow_total += b;
        crate::parity_flow_iter_dump::emit(
            &__flow_bundle_chrom, &__flow_bundle_strand, &__flow_bundle_id,
            __flow_run_id, __flow_iter, &__path, b, __flow_total,
        );
    }

    // existing flow_mat update follows here
}
```

The bundle metadata (`__flow_bundle_chrom` etc.) needs to be hoisted from the caller. Read the function signature and trace where these are available.

> **Heuristic:** if the function only has access to indices and not the bundle's chrom/strand/id strings, pass those in as new arguments. Update all call sites of the modified function. Use a quick `grep -n "fn weight_bundle\|fn edmonds_karp" max_flow.rs` to find the function names, then `grep -rn "weight_bundle(\|edmonds_karp(" src/rustle/` to find callers.

If passing strings through is too disruptive, an acceptable fallback is to pass placeholder strings (`"unknown"`, `"."`) and accept that the per-iter TSV won't have the bundle coords. The downstream visualizer joins on `path_nodes` against the graph dump anyway. This degrades the TSV but keeps the patch small.

- [ ] **Step 2.8 — Verify test passes**

Run: `cargo build --release 2>&1 | tail -10 && cargo test --release --test test_parity_flow_iter_dump 2>&1 | tail -10`

Expected: build green, test PASS with ≥2 iterations recorded for synthetic.

- [ ] **Step 2.9 — Smoke-check the output**

Run: `RUSTLE_FLOW_ITER_TSV=/tmp/flow.tsv ./target/release/rustle -L -o /tmp/s.gtf test_data/synthetic_family/reads_sorted.bam && head -20 /tmp/flow.tsv && wc -l /tmp/flow.tsv`

Expected: ~3-6 iteration rows (one or two iterations per bundle, two bundles). Bottleneck values strictly positive; total_flow_after monotonically non-decreasing within a bundle_run.

- [ ] **Step 2.10 — Commit**

```bash
git add src/rustle/parity_flow_iter_dump.rs src/rustle/lib.rs src/rustle/max_flow.rs tests/regression/test_parity_flow_iter_dump.rs
git commit -m "feat(parity): add RUSTLE_FLOW_ITER_TSV dump for max-flow augmenting paths"
```

---

## Phase 2 — Python visualizers

Independent of each other (after Phase 1 done). V4 and V5 don't depend on Phase 1 dumps; can be parallelized via subagents.

### Task 3: Tools-demo dir setup + shared utilities

**Files:**
- Create: `tools/demo/__init__.py` (empty)
- Create: `tools/demo/common.py` — shared TSV loaders, color palette, graphviz wrappers
- Create: `tools/demo/tests/__init__.py` (empty)

#### Step 3.1 — Create the package structure

- [ ] **Step 3.1.1**

```bash
mkdir -p /mnt/c/Users/jfris/Desktop/Rustle/tools/demo/tests
touch /mnt/c/Users/jfris/Desktop/Rustle/tools/demo/__init__.py
touch /mnt/c/Users/jfris/Desktop/Rustle/tools/demo/tests/__init__.py
```

- [ ] **Step 3.1.2 — Create `common.py`**

```python
# tools/demo/common.py
"""Shared utilities for demo visualizers.

All renderers consume the same TSV schemas. This module centralizes the schemas,
palette, and graphviz wrapping so the renderers can stay focused.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Sequence

# Color palette — keep accessible (deuteranopia-safe).
PALETTE = {
    "node_default":     "#E8F1F8",
    "node_source_sink": "#F8E5C8",
    "edge_collinear":   "#1F4E79",
    "edge_junction":    "#C0504D",
    "edge_residual":    "#9D9D9D",
    "edge_augmenting":  "#E37222",
    "iso_pair":         ["#1F77B4", "#2CA02C", "#D62728", "#9467BD", "#FF7F0E"],
}


@dataclass(frozen=True)
class GraphNode:
    bdstart: int
    bdend: int
    strand: str
    chrom: str
    node_idx: int
    start: int
    end: int
    cov: float
    hardstart: bool
    hardend: bool
    nchildren: int
    nparents: int


@dataclass(frozen=True)
class GraphEdge:
    bdstart: int
    bdend: int
    strand: str
    chrom: str
    from_idx: int
    to_idx: int
    from_start: int
    from_end: int
    to_start: int
    to_end: int
    edge_kind: str   # "collinear" or "junction"
    from_cov: float
    to_cov: float
    bottleneck_cov: float


def _bundle_match(row: dict, locus: tuple[str, int, int] | None) -> bool:
    if locus is None:
        return True
    chrom, start, end = locus
    if row["chrom"] != chrom:
        return False
    bs, be = int(row["bdstart"]), int(row["bdend"])
    return be >= start and bs <= end


def read_nodes_tsv(path: Path, locus: tuple[str, int, int] | None = None) -> list[GraphNode]:
    out: list[GraphNode] = []
    with Path(path).open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if not _bundle_match(row, locus):
                continue
            out.append(GraphNode(
                bdstart=int(row["bdstart"]),
                bdend=int(row["bdend"]),
                strand=row["strand"],
                chrom=row["chrom"],
                node_idx=int(row["node_idx"]),
                start=int(row["start"]),
                end=int(row["end"]),
                cov=float(row["cov"]),
                hardstart=row["hardstart"] in ("1", "true", "True"),
                hardend=row["hardend"] in ("1", "true", "True"),
                nchildren=int(row["nchildren"]),
                nparents=int(row["nparents"]),
            ))
    return out


def read_edges_tsv(path: Path, locus: tuple[str, int, int] | None = None) -> list[GraphEdge]:
    out: list[GraphEdge] = []
    with Path(path).open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if not _bundle_match(row, locus):
                continue
            out.append(GraphEdge(
                bdstart=int(row["bdstart"]),
                bdend=int(row["bdend"]),
                strand=row["strand"],
                chrom=row["chrom"],
                from_idx=int(row["from_idx"]),
                to_idx=int(row["to_idx"]),
                from_start=int(row["from_start"]),
                from_end=int(row["from_end"]),
                to_start=int(row["to_start"]),
                to_end=int(row["to_end"]),
                edge_kind=row["edge_kind"],
                from_cov=float(row["from_cov"]),
                to_cov=float(row["to_cov"]),
                bottleneck_cov=float(row["bottleneck_cov"]),
            ))
    return out


def group_by_bundle(nodes: Sequence[GraphNode]) -> Iterator[tuple[tuple[str, int, int, str], list[GraphNode]]]:
    """Yield ((chrom, bdstart, bdend, strand), nodes_in_bundle)."""
    keyed: dict[tuple, list[GraphNode]] = {}
    for n in nodes:
        keyed.setdefault((n.chrom, n.bdstart, n.bdend, n.strand), []).append(n)
    for k, v in keyed.items():
        yield k, sorted(v, key=lambda n: n.node_idx)
```

- [ ] **Step 3.1.3 — Commit**

```bash
git add tools/demo/__init__.py tools/demo/common.py tools/demo/tests/__init__.py
git commit -m "scaffold(demo): tools/demo package + shared TSV schemas"
```

---

### Task 4: V1 — Single splice-graph renderer

**Files:**
- Create: `tools/demo/render_splice_graph.py`
- Test: `tools/demo/tests/test_render_splice_graph.py`

#### Step 4.1 — Write failing test

```python
# tools/demo/tests/test_render_splice_graph.py
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TOOL = REPO / "tools" / "demo" / "render_splice_graph.py"
NODES_FIXTURE = Path(__file__).parent / "fixtures" / "small_nodes.tsv"
EDGES_FIXTURE = Path(__file__).parent / "fixtures" / "small_edges.tsv"


def test_render_produces_png_for_small_graph(tmp_path):
    out_png = tmp_path / "graph.png"
    out_dot = tmp_path / "graph.dot"

    cp = subprocess.run([
        "python3", str(TOOL),
        "--nodes-tsv", str(NODES_FIXTURE),
        "--edges-tsv", str(EDGES_FIXTURE),
        "--output-png", str(out_png),
        "--output-dot", str(out_dot),
        "--locus", "chr_test:0-20000",
    ], capture_output=True, text=True)

    assert cp.returncode == 0, cp.stderr
    assert out_png.exists() and out_png.stat().st_size > 100
    assert out_dot.exists()

    dot = out_dot.read_text()
    assert "digraph" in dot
    assert "label=" in dot
    # Three nodes + 2 edges in fixture, expect at least 3 node entries.
    assert dot.count("[shape=") >= 3
```

- [ ] **Step 4.2 — Create test fixtures**

```python
# tools/demo/tests/fixtures/_create_fixtures.py — run once to seed
# (also commit the resulting TSVs)
import csv
from pathlib import Path

OUT = Path(__file__).parent

NODES = [
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","node_idx":0,"start":1000,"end":1200,"cov":12.5,"hardstart":1,"hardend":0,"nchildren":1,"nparents":0},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","node_idx":1,"start":2000,"end":2300,"cov":15.0,"hardstart":0,"hardend":0,"nchildren":2,"nparents":1},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","node_idx":2,"start":3000,"end":3150,"cov":11.0,"hardstart":0,"hardend":0,"nchildren":1,"nparents":1},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","node_idx":3,"start":4000,"end":4500,"cov":13.0,"hardstart":0,"hardend":1,"nchildren":0,"nparents":2},
]
EDGES = [
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","from_idx":0,"to_idx":1,"from_start":1000,"from_end":1200,"to_start":2000,"to_end":2300,"edge_kind":"junction","from_cov":12.5,"to_cov":15.0,"bottleneck_cov":12.5},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","from_idx":1,"to_idx":2,"from_start":2000,"from_end":2300,"to_start":3000,"to_end":3150,"edge_kind":"junction","from_cov":15.0,"to_cov":11.0,"bottleneck_cov":11.0},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","from_idx":1,"to_idx":3,"from_start":2000,"from_end":2300,"to_start":4000,"to_end":4500,"edge_kind":"junction","from_cov":15.0,"to_cov":13.0,"bottleneck_cov":13.0},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","from_idx":2,"to_idx":3,"from_start":3000,"from_end":3150,"to_start":4000,"to_end":4500,"edge_kind":"junction","from_cov":11.0,"to_cov":13.0,"bottleneck_cov":11.0},
]

def write(name, fields, rows):
    p = OUT / name
    with p.open("w") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

write("small_nodes.tsv", list(NODES[0].keys()), NODES)
write("small_edges.tsv", list(EDGES[0].keys()), EDGES)
```

Run once: `python3 tools/demo/tests/fixtures/_create_fixtures.py` — produces `small_nodes.tsv` and `small_edges.tsv`. Commit them.

- [ ] **Step 4.3 — Run test, verify it fails**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && python3 -m pytest tools/demo/tests/test_render_splice_graph.py -v`

Expected: FAIL — tool not yet implemented.

- [ ] **Step 4.4 — Implement `render_splice_graph.py`**

```python
#!/usr/bin/env python3
"""Render a single splice graph from RUSTLE_PARITY_GRAPH_TSV + RUSTLE_PARITY_GRAPH_EDGES_TSV.

Mechanical layout via graphviz dot. No hand-coded positioning. Output is PNG
plus the underlying .dot file (for reproducibility — paste into a graphviz
viewer if PNG is contested).

Usage:
    render_splice_graph.py --nodes-tsv ... --edges-tsv ... \\
        --output-png out.png --output-dot out.dot \\
        [--locus chrom:start-end] [--title "Title"]
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from demo.common import (
    PALETTE,
    GraphEdge,
    GraphNode,
    group_by_bundle,
    read_edges_tsv,
    read_nodes_tsv,
)


def _parse_locus(s: str | None) -> tuple[str, int, int] | None:
    if not s:
        return None
    chrom, rng = s.split(":")
    start, end = rng.split("-")
    return (chrom, int(start), int(end))


def _node_label(n: GraphNode) -> str:
    badge = ""
    if n.hardstart:
        badge += "▶"
    if n.hardend:
        badge += "◀"
    return f"#{n.node_idx}{badge}\\n[{n.start}-{n.end}]\\ncov={n.cov:.1f}"


def _edge_color(e: GraphEdge) -> str:
    return PALETTE["edge_junction"] if e.edge_kind == "junction" else PALETTE["edge_collinear"]


def _edge_label(e: GraphEdge) -> str:
    return f"{e.edge_kind}\\nbn={e.bottleneck_cov:.1f}"


def to_dot(nodes: list[GraphNode], edges: list[GraphEdge], title: str = "") -> str:
    out = ['digraph G {']
    out.append('  rankdir=LR;')
    out.append('  node [fontname="Helvetica" fontsize=10];')
    out.append('  edge [fontname="Helvetica" fontsize=9];')
    if title:
        out.append(f'  labelloc="t"; label="{title}"; fontsize=14;')

    for n in nodes:
        fill = PALETTE["node_source_sink"] if (n.hardstart or n.hardend) else PALETTE["node_default"]
        label = _node_label(n).replace('"', '\\"')
        out.append(
            f'  n{n.node_idx} [shape=box style="rounded,filled" fillcolor="{fill}" label="{label}"];'
        )
    for e in edges:
        out.append(
            f'  n{e.from_idx} -> n{e.to_idx} '
            f'[color="{_edge_color(e)}" label="{_edge_label(e)}" '
            f'penwidth={1.0 + min(4.0, e.bottleneck_cov / 5.0):.1f}];'
        )
    out.append('}')
    return "\n".join(out)


def render(dot_text: str, output_png: Path, output_dot: Path | None) -> None:
    if output_dot:
        output_dot.write_text(dot_text)
    p = subprocess.run(
        ["dot", "-Tpng", "-o", str(output_png)],
        input=dot_text, text=True, capture_output=True,
    )
    if p.returncode != 0:
        raise SystemExit(f"graphviz failed: {p.stderr}")


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodes-tsv", required=True, type=Path)
    ap.add_argument("--edges-tsv", required=True, type=Path)
    ap.add_argument("--output-png", required=True, type=Path)
    ap.add_argument("--output-dot", type=Path, default=None)
    ap.add_argument("--locus", default=None)
    ap.add_argument("--title", default="")
    args = ap.parse_args(argv)

    locus = _parse_locus(args.locus)
    nodes = read_nodes_tsv(args.nodes_tsv, locus)
    edges = read_edges_tsv(args.edges_tsv, locus)

    if not nodes:
        raise SystemExit(f"No nodes match locus {args.locus} in {args.nodes_tsv}")

    bundles = list(group_by_bundle(nodes))
    if len(bundles) > 1:
        # If multiple bundles match, render only the first and warn.
        print(f"WARN: locus matched {len(bundles)} bundles; rendering first only", file=sys.stderr)
    (key, bnodes) = bundles[0]
    chrom, bs, be, _ = key
    bedges = [e for e in edges if e.chrom == chrom and e.bdstart == bs and e.bdend == be]

    title = args.title or f"{chrom}:{bs}-{be} ({len(bnodes)} nodes, {len(bedges)} edges)"
    dot = to_dot(bnodes, bedges, title=title)
    render(dot, args.output_png, args.output_dot)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
```

- [ ] **Step 4.5 — Run test, verify pass**

Run: `cd /mnt/c/Users/jfris/Desktop/Rustle && python3 -m pytest tools/demo/tests/test_render_splice_graph.py -v`

Expected: PASS.

- [ ] **Step 4.6 — Visual smoke check**

Run:
```bash
python3 tools/demo/render_splice_graph.py \
  --nodes-tsv tools/demo/tests/fixtures/small_nodes.tsv \
  --edges-tsv tools/demo/tests/fixtures/small_edges.tsv \
  --output-png /tmp/v1_smoke.png --output-dot /tmp/v1_smoke.dot \
  --locus chr_test:0-20000
file /tmp/v1_smoke.png && wc -c /tmp/v1_smoke.png
```

Expected: PNG image, > 1 KB.

- [ ] **Step 4.7 — Commit**

```bash
git add tools/demo/render_splice_graph.py tools/demo/tests/test_render_splice_graph.py tools/demo/tests/fixtures/
git commit -m "feat(demo): V1 single splice-graph renderer (graphviz)"
```

---

### Task 5: V2 — Three-way isomorphism comparator

**Files:**
- Create: `tools/demo/compare_graphs_isomorphism.py`
- Test: `tools/demo/tests/test_compare_graphs_isomorphism.py`

#### Step 5.1 — Write failing test for isomorphism logic

```python
# tools/demo/tests/test_compare_graphs_isomorphism.py
import subprocess
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TOOL = REPO / "tools" / "demo" / "compare_graphs_isomorphism.py"
FIXTURES = Path(__file__).parent / "fixtures"


def test_isomorphic_three_copies_produces_bijection(tmp_path):
    # All 3 copies have the same shape; iso must be found.
    out_json = tmp_path / "iso.json"
    out_png = tmp_path / "panel.png"

    cp = subprocess.run([
        "python3", str(TOOL),
        "--nodes-tsv", str(FIXTURES / "iso_triple_nodes.tsv"),
        "--edges-tsv", str(FIXTURES / "iso_triple_edges.tsv"),
        "--locus", "copy1:0-10000",
        "--locus", "copy2:0-10000",
        "--locus", "copy3:0-10000",
        "--output-json", str(out_json),
        "--output-png", str(out_png),
    ], capture_output=True, text=True)

    assert cp.returncode == 0, cp.stderr
    data = json.loads(out_json.read_text())
    assert data["isomorphic"] is True
    # bijection: {pair_label: {a_node_idx: b_node_idx}}
    assert "bijection" in data
    assert all(len(v) >= 3 for v in data["bijection"].values())


def test_non_isomorphic_reports_partial(tmp_path):
    out_json = tmp_path / "iso.json"
    out_png = tmp_path / "panel.png"

    cp = subprocess.run([
        "python3", str(TOOL),
        "--nodes-tsv", str(FIXTURES / "non_iso_nodes.tsv"),
        "--edges-tsv", str(FIXTURES / "non_iso_edges.tsv"),
        "--locus", "copyA:0-10000",
        "--locus", "copyB:0-10000",
        "--output-json", str(out_json),
        "--output-png", str(out_png),
    ], capture_output=True, text=True)

    assert cp.returncode == 0, cp.stderr
    data = json.loads(out_json.read_text())
    assert data["isomorphic"] is False
    assert "diff" in data
```

- [ ] **Step 5.2 — Create iso/non-iso fixtures**

Use the same fixture-creator pattern as Task 4 to produce four files with three identical-shape graphs (copy1/copy2/copy3) and two different-shape graphs (copyA: 4 nodes, copyB: 5 nodes).

```python
# tools/demo/tests/fixtures/_create_iso_fixtures.py (run once, commit the TSVs)
import csv
from pathlib import Path

OUT = Path(__file__).parent
NODE_FIELDS = ["source","chrom","bdstart","bdend","strand","node_idx","start","end","cov","hardstart","hardend","nchildren","nparents"]
EDGE_FIELDS = ["source","chrom","bdstart","bdend","strand","from_idx","to_idx","from_start","from_end","to_start","to_end","edge_kind","from_cov","to_cov","bottleneck_cov"]


def make_chain_graph(copy: str, offset: int):
    """4-node linear-with-skip chain, all junctions, fixed shape."""
    nodes = []
    for i, (s, e, cov) in enumerate([(100,200,10),(400,500,11),(700,800,9),(1000,1100,12)]):
        nodes.append({
            "source":"rustle","chrom":copy,"bdstart":offset,"bdend":offset+1100,"strand":"+",
            "node_idx":i,"start":offset+s,"end":offset+e,"cov":cov,
            "hardstart":1 if i==0 else 0,"hardend":1 if i==3 else 0,
            "nchildren":2 if i==1 else (1 if i<3 else 0),
            "nparents":2 if i==3 else (1 if i>0 else 0),
        })
    edges = []
    for (a,b) in [(0,1),(1,2),(1,3),(2,3)]:
        edges.append({
            "source":"rustle","chrom":copy,"bdstart":offset,"bdend":offset+1100,"strand":"+",
            "from_idx":a,"to_idx":b,
            "from_start":nodes[a]["start"],"from_end":nodes[a]["end"],
            "to_start":nodes[b]["start"],"to_end":nodes[b]["end"],
            "edge_kind":"junction",
            "from_cov":nodes[a]["cov"],"to_cov":nodes[b]["cov"],
            "bottleneck_cov":min(nodes[a]["cov"], nodes[b]["cov"]),
        })
    return nodes, edges


def write(name, fields, rows):
    with (OUT/name).open("w") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)


# Iso triple
all_n, all_e = [], []
for copy in ("copy1","copy2","copy3"):
    n,e = make_chain_graph(copy, 0)
    all_n += n; all_e += e
write("iso_triple_nodes.tsv", NODE_FIELDS, all_n)
write("iso_triple_edges.tsv", EDGE_FIELDS, all_e)

# Non-iso pair: copyA = 4-node chain, copyB = 5-node linear (no skip)
n_a, e_a = make_chain_graph("copyA", 0)

# 5-node linear copyB
n_b = []
for i,(s,e,cov) in enumerate([(100,200,10),(300,400,11),(500,600,9),(700,800,10),(900,1000,12)]):
    n_b.append({
        "source":"rustle","chrom":"copyB","bdstart":0,"bdend":1000,"strand":"+",
        "node_idx":i,"start":s,"end":e,"cov":cov,
        "hardstart":1 if i==0 else 0,"hardend":1 if i==4 else 0,
        "nchildren":1 if i<4 else 0,
        "nparents":1 if i>0 else 0,
    })
e_b = []
for a,b in [(0,1),(1,2),(2,3),(3,4)]:
    e_b.append({
        "source":"rustle","chrom":"copyB","bdstart":0,"bdend":1000,"strand":"+",
        "from_idx":a,"to_idx":b,
        "from_start":n_b[a]["start"],"from_end":n_b[a]["end"],
        "to_start":n_b[b]["start"],"to_end":n_b[b]["end"],
        "edge_kind":"junction",
        "from_cov":n_b[a]["cov"],"to_cov":n_b[b]["cov"],
        "bottleneck_cov":min(n_b[a]["cov"], n_b[b]["cov"]),
    })
write("non_iso_nodes.tsv", NODE_FIELDS, n_a + n_b)
write("non_iso_edges.tsv", EDGE_FIELDS, e_a + e_b)
```

Run once. Commit the resulting TSVs.

- [ ] **Step 5.3 — Run test, verify fail**

`pytest tools/demo/tests/test_compare_graphs_isomorphism.py -v` — FAIL.

- [ ] **Step 5.4 — Implement `compare_graphs_isomorphism.py`**

```python
#!/usr/bin/env python3
"""Three-way splice-graph isomorphism comparator.

For each --locus given (≥2), build a NetworkX DiGraph from the matching nodes
and edges, then test pairwise isomorphism using VF2.

"Similar" definition: graph isomorphism preserving (a) node count, (b) edge
multiset by edge_kind, (c) source/sink topology (hardstart/hardend nodes).
Genomic coordinates normalized by offset within the locus before comparison.

Usage:
    compare_graphs_isomorphism.py \\
        --nodes-tsv ... --edges-tsv ... \\
        --locus chr19:104789647-104796276 \\
        --locus chr19:104830536-104837094 \\
        --locus chr19:104871356-104877901 \\
        --output-json out.json --output-png out.png
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from demo.common import (
    GraphEdge,
    GraphNode,
    group_by_bundle,
    read_edges_tsv,
    read_nodes_tsv,
)


def _parse_locus(s: str) -> tuple[str, int, int]:
    chrom, rng = s.split(":")
    start, end = rng.split("-")
    return (chrom, int(start), int(end))


def _build_nx(nodes: list[GraphNode], edges: list[GraphEdge]) -> nx.DiGraph:
    g = nx.DiGraph()
    if not nodes:
        return g
    bundles = list(group_by_bundle(nodes))
    (key, bnodes) = bundles[0]
    chrom, bs, be, _ = key
    bedges = [e for e in edges if e.chrom == chrom and e.bdstart == bs and e.bdend == be]
    # Normalize coordinates by subtracting bdstart so the same shape at different
    # genomic positions matches.
    for n in bnodes:
        g.add_node(
            n.node_idx,
            type="source" if n.hardstart else ("sink" if n.hardend else "internal"),
            length=n.end - n.start,
        )
    for e in bedges:
        g.add_edge(e.from_idx, e.to_idx, kind=e.edge_kind)
    return g


def _node_match(a: dict, b: dict) -> bool:
    return a.get("type") == b.get("type")


def _edge_match(a: dict, b: dict) -> bool:
    return a.get("kind") == b.get("kind")


def _compare_pair(g1: nx.DiGraph, g2: nx.DiGraph) -> dict:
    matcher = DiGraphMatcher(g1, g2, node_match=_node_match, edge_match=_edge_match)
    is_iso = matcher.is_isomorphic()
    out = {"isomorphic": is_iso}
    if is_iso:
        out["bijection"] = {str(k): v for k, v in matcher.mapping.items()}
    else:
        out["diff"] = {
            "node_counts": (g1.number_of_nodes(), g2.number_of_nodes()),
            "edge_counts": (g1.number_of_edges(), g2.number_of_edges()),
            "edge_kind_counts": (
                _kind_counts(g1), _kind_counts(g2),
            ),
        }
    return out


def _kind_counts(g: nx.DiGraph) -> dict:
    c: dict[str, int] = {}
    for _, _, d in g.edges(data=True):
        c[d.get("kind", "?")] = c.get(d.get("kind", "?"), 0) + 1
    return c


def _render_panel(loci, dots: list[str], output_png: Path, mappings: dict | None) -> None:
    """Compose N graphviz dot images into a horizontal panel via matplotlib."""
    import matplotlib.pyplot as plt
    from matplotlib.image import imread

    with TemporaryDirectory() as td:
        td_path = Path(td)
        pngs = []
        for i, dot in enumerate(dots):
            png = td_path / f"copy{i}.png"
            p = subprocess.run(
                ["dot", "-Tpng", "-o", str(png)],
                input=dot, text=True, capture_output=True,
            )
            if p.returncode != 0:
                raise SystemExit(f"graphviz failed for copy {i}: {p.stderr}")
            pngs.append(png)

        fig, axes = plt.subplots(1, len(pngs), figsize=(6 * len(pngs), 6))
        if len(pngs) == 1:
            axes = [axes]
        for ax, png, locus in zip(axes, pngs, loci):
            ax.imshow(imread(png))
            ax.set_title(f"{locus[0]}:{locus[1]}-{locus[2]}", fontsize=10)
            ax.axis("off")
        if mappings:
            note = "Bijection: " + "  /  ".join(
                f"{k}: " + ", ".join(f"{a}↔{b}" for a, b in sorted(v.items())[:5])
                + ("…" if len(v) > 5 else "")
                for k, v in mappings.items()
            )
            fig.suptitle(note, fontsize=10)
        fig.tight_layout()
        fig.savefig(output_png, dpi=150, bbox_inches="tight")
        plt.close(fig)


def _to_dot_for_graph(nodes: list[GraphNode], edges: list[GraphEdge], title: str) -> str:
    # Reuse render_splice_graph.to_dot for visual consistency.
    from demo.render_splice_graph import to_dot
    bundles = list(group_by_bundle(nodes))
    (key, bnodes) = bundles[0]
    chrom, bs, be, _ = key
    bedges = [e for e in edges if e.chrom == chrom and e.bdstart == bs and e.bdend == be]
    return to_dot(bnodes, bedges, title=title)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--nodes-tsv", required=True, type=Path)
    ap.add_argument("--edges-tsv", required=True, type=Path)
    ap.add_argument("--locus", action="append", required=True,
                    help="chrom:start-end (specify ≥2 times)")
    ap.add_argument("--output-json", required=True, type=Path)
    ap.add_argument("--output-png", required=True, type=Path)
    args = ap.parse_args(argv)

    if len(args.locus) < 2:
        raise SystemExit("need ≥2 --locus arguments")

    loci = [_parse_locus(s) for s in args.locus]
    graphs = []
    dots = []
    for locus in loci:
        nodes = read_nodes_tsv(args.nodes_tsv, locus)
        edges = read_edges_tsv(args.edges_tsv, locus)
        if not nodes:
            raise SystemExit(f"no nodes for {locus}")
        g = _build_nx(nodes, edges)
        graphs.append(g)
        dots.append(_to_dot_for_graph(nodes, edges, title=f"{locus[0]}:{locus[1]}-{locus[2]}"))

    # Compare all pairs (1↔2, 1↔3, etc.)
    pairwise = {}
    bijections = {}
    for i in range(len(graphs)):
        for j in range(i + 1, len(graphs)):
            label = f"{i+1}↔{j+1}"
            res = _compare_pair(graphs[i], graphs[j])
            pairwise[label] = res
            if res.get("isomorphic"):
                bijections[label] = res["bijection"]

    all_iso = all(p["isomorphic"] for p in pairwise.values())
    out = {
        "isomorphic": all_iso,
        "pairwise": pairwise,
        "bijection": bijections if all_iso else {},
        "graph_summary": [
            {"locus": f"{c}:{s}-{e}", "n_nodes": g.number_of_nodes(),
             "n_edges": g.number_of_edges(),
             "edge_kinds": _kind_counts(g)}
            for (c, s, e), g in zip(loci, graphs)
        ],
    }
    if not all_iso:
        out["diff"] = {k: v.get("diff") for k, v in pairwise.items() if not v["isomorphic"]}

    args.output_json.write_text(json.dumps(out, indent=2))
    _render_panel(loci, dots, args.output_png, bijections if all_iso else None)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
```

- [ ] **Step 5.5 — Run tests, verify pass**

Run: `python3 -m pytest tools/demo/tests/test_compare_graphs_isomorphism.py -v`

Expected: both tests PASS.

- [ ] **Step 5.6 — Commit**

```bash
git add tools/demo/compare_graphs_isomorphism.py tools/demo/tests/test_compare_graphs_isomorphism.py tools/demo/tests/fixtures/iso_triple_*.tsv tools/demo/tests/fixtures/non_iso_*.tsv
git commit -m "feat(demo): V2 three-way splice-graph isomorphism comparator (NetworkX VF2)"
```

---

### Task 6: V3 — Per-iteration flow decomposition renderer

**Files:**
- Create: `tools/demo/render_flow_iterations.py`
- Test: `tools/demo/tests/test_render_flow_iterations.py`

#### Step 6.1 — Write failing test

```python
# tools/demo/tests/test_render_flow_iterations.py
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TOOL = REPO / "tools" / "demo" / "render_flow_iterations.py"
FIXTURES = Path(__file__).parent / "fixtures"


def test_renders_iteration_grid(tmp_path):
    out_png = tmp_path / "flow.png"
    cp = subprocess.run([
        "python3", str(TOOL),
        "--flow-tsv", str(FIXTURES / "flow_3iter.tsv"),
        "--nodes-tsv", str(FIXTURES / "small_nodes.tsv"),
        "--edges-tsv", str(FIXTURES / "small_edges.tsv"),
        "--locus", "chr_test:0-20000",
        "--output-png", str(out_png),
    ], capture_output=True, text=True)
    assert cp.returncode == 0, cp.stderr
    assert out_png.exists() and out_png.stat().st_size > 1000
```

- [ ] **Step 6.2 — Create flow fixture**

```python
# in tools/demo/tests/fixtures/_create_flow_fixtures.py
import csv
from pathlib import Path

OUT = Path(__file__).parent

ROWS = [
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","bundle_run":1,"iter_idx":1,"path_len":4,"path_nodes":"S,0,1,3,T","bottleneck":12.5,"total_flow_after":12.5},
    {"source":"rustle","chrom":"chr_test","bdstart":1000,"bdend":4500,"strand":"-","bundle_run":1,"iter_idx":2,"path_len":5,"path_nodes":"S,0,1,2,3,T","bottleneck":2.0,"total_flow_after":14.5},
]

with (OUT/"flow_3iter.tsv").open("w") as f:
    w = csv.DictWriter(f, fieldnames=list(ROWS[0].keys()), delimiter="\t")
    w.writeheader(); w.writerows(ROWS)
```

Run once. Commit.

- [ ] **Step 6.3 — Run test, verify fail**

`pytest tools/demo/tests/test_render_flow_iterations.py -v` → FAIL.

- [ ] **Step 6.4 — Implement renderer**

```python
#!/usr/bin/env python3
"""Per-iteration flow decomposition renderer.

For each iteration row in the flow TSV, render the splice graph with the
augmenting path highlighted in orange, capacity-after-iteration shown on
edges. Compose into a horizontal grid.

Usage:
    render_flow_iterations.py --flow-tsv ... --nodes-tsv ... --edges-tsv ...
        --locus chrom:start-end --output-png out.png
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from demo.common import GraphEdge, GraphNode, group_by_bundle, read_edges_tsv, read_nodes_tsv, PALETTE
from demo.render_splice_graph import to_dot


@dataclass(frozen=True)
class FlowIter:
    chrom: str
    bdstart: int
    bdend: int
    bundle_run: int
    iter_idx: int
    path_nodes: list[str]
    bottleneck: float
    total_flow_after: float


def _parse_locus(s: str) -> tuple[str, int, int]:
    chrom, rng = s.split(":")
    start, end = rng.split("-")
    return (chrom, int(start), int(end))


def read_flow_tsv(path: Path, locus: tuple[str, int, int]) -> list[FlowIter]:
    out = []
    chrom_f, lstart, lend = locus
    with Path(path).open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["chrom"] != chrom_f:
                continue
            bs, be = int(row["bdstart"]), int(row["bdend"])
            if be < lstart or bs > lend:
                continue
            out.append(FlowIter(
                chrom=row["chrom"], bdstart=bs, bdend=be,
                bundle_run=int(row["bundle_run"]),
                iter_idx=int(row["iter_idx"]),
                path_nodes=row["path_nodes"].split(","),
                bottleneck=float(row["bottleneck"]),
                total_flow_after=float(row["total_flow_after"]),
            ))
    return sorted(out, key=lambda r: (r.bundle_run, r.iter_idx))


def _highlighted_dot(nodes: list[GraphNode], edges: list[GraphEdge], path_node_strs: list[str], iter_idx: int, bottleneck: float, title_extra: str) -> str:
    """Same as render_splice_graph.to_dot but highlights edges along the augmenting path in orange."""
    # Identify path edge pairs: consecutive node strs (skipping S/T sentinels).
    path_idx_pairs = set()
    cleaned = [p for p in path_node_strs if p not in ("S", "T")]
    for a, b in zip(cleaned, cleaned[1:]):
        try:
            path_idx_pairs.add((int(a), int(b)))
        except ValueError:
            continue

    # Reuse to_dot but tag highlighted edges by intercepting:
    base = to_dot(nodes, edges, title=f"iter {iter_idx} — bottleneck={bottleneck:.2f}{title_extra}")
    out_lines = []
    for line in base.splitlines():
        emitted = False
        for (a, b) in path_idx_pairs:
            tag = f"n{a} -> n{b}"
            if line.lstrip().startswith(tag):
                # rewrite color + penwidth
                line = line.replace(
                    'color="' + PALETTE["edge_collinear"] + '"',
                    f'color="{PALETTE["edge_augmenting"]}"'
                ).replace(
                    'color="' + PALETTE["edge_junction"] + '"',
                    f'color="{PALETTE["edge_augmenting"]}"'
                )
                # Bump pen width
                if "penwidth=" in line:
                    pre, _, post = line.rpartition("penwidth=")
                    val_str = post.split("]")[0]
                    line = pre + "penwidth=4.0" + post[len(val_str):]
                emitted = True
                break
        out_lines.append(line)
    return "\n".join(out_lines)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--flow-tsv", required=True, type=Path)
    ap.add_argument("--nodes-tsv", required=True, type=Path)
    ap.add_argument("--edges-tsv", required=True, type=Path)
    ap.add_argument("--locus", required=True)
    ap.add_argument("--output-png", required=True, type=Path)
    args = ap.parse_args(argv)

    locus = _parse_locus(args.locus)
    iters = read_flow_tsv(args.flow_tsv, locus)
    nodes = read_nodes_tsv(args.nodes_tsv, locus)
    edges = read_edges_tsv(args.edges_tsv, locus)
    if not iters:
        raise SystemExit(f"no flow iterations for locus {args.locus}")
    if not nodes:
        raise SystemExit(f"no nodes for locus {args.locus}")

    # Render one PNG per iteration, then stitch.
    import matplotlib.pyplot as plt
    from matplotlib.image import imread

    with TemporaryDirectory() as td:
        td_path = Path(td)
        pngs = []
        for it in iters:
            dot = _highlighted_dot(
                nodes, edges, it.path_nodes, it.iter_idx, it.bottleneck,
                title_extra=f"  total flow={it.total_flow_after:.2f}"
            )
            png = td_path / f"iter_{it.bundle_run:02d}_{it.iter_idx:02d}.png"
            p = subprocess.run(
                ["dot", "-Tpng", "-o", str(png)],
                input=dot, text=True, capture_output=True,
            )
            if p.returncode != 0:
                raise SystemExit(f"dot failed: {p.stderr}")
            pngs.append((it, png))

        cols = min(3, len(pngs))
        rows = (len(pngs) + cols - 1) // cols
        fig, axes = plt.subplots(rows, cols, figsize=(7 * cols, 5 * rows))
        if rows == 1 and cols == 1:
            axes = [[axes]]
        elif rows == 1:
            axes = [axes]
        elif cols == 1:
            axes = [[a] for a in axes]
        for k, (it, png) in enumerate(pngs):
            r, c = k // cols, k % cols
            ax = axes[r][c]
            ax.imshow(imread(png))
            ax.set_title(f"iter {it.iter_idx}: bn={it.bottleneck:.2f}, total={it.total_flow_after:.2f}", fontsize=10)
            ax.axis("off")
        # blank out unused subplots
        for k in range(len(pngs), rows * cols):
            r, c = k // cols, k % cols
            axes[r][c].axis("off")
        fig.suptitle(f"Edmonds-Karp iterations at {args.locus}", fontsize=12)
        fig.tight_layout()
        fig.savefig(args.output_png, dpi=150, bbox_inches="tight")
        plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
```

- [ ] **Step 6.5 — Run test, verify pass**

`python3 -m pytest tools/demo/tests/test_render_flow_iterations.py -v` → PASS.

- [ ] **Step 6.6 — Commit**

```bash
git add tools/demo/render_flow_iterations.py tools/demo/tests/test_render_flow_iterations.py tools/demo/tests/fixtures/flow_3iter.tsv
git commit -m "feat(demo): V3 per-iteration flow decomposition renderer"
```

---

### Task 7: V4 — Multi-mapper read pinning visualizer

**Files:**
- Create: `tools/demo/render_multimapper_pinning.py`
- Test: `tools/demo/tests/test_render_multimapper_pinning.py`

This visualizer parses Rustle's existing VG stderr trace (no new dump needed). It reads the trace, extracts per-multi-mapper EM weights before/after, picks one specific read crossing N copies, and produces a side-by-side bar chart: "StringTie 1/NH uniform" vs "Rustle EM after convergence".

#### Step 7.1 — Inspect VG stderr trace format

Run the assembler on synthetic_family with VG enabled and `RUSTLE_VG_TRACE=1`, capture stderr, identify the relevant lines.

```bash
RUSTLE_VG_TRACE=1 ./target/release/rustle -L --vg \
  -o /tmp/synth_vg.gtf test_data/synthetic_family/reads_sorted.bam \
  2> /tmp/vg_trace.log
grep -E "^\[VG\]" /tmp/vg_trace.log | head -40
```

Record the exact format. If the trace is sparse, fall back to `RUSTLE_TRACE_VGSNP=1` or another VG env var (per the audit table from earlier exploration).

#### Step 7.2 — Write failing test

```python
# tools/demo/tests/test_render_multimapper_pinning.py
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
TOOL = REPO / "tools" / "demo" / "render_multimapper_pinning.py"
FIXTURES = Path(__file__).parent / "fixtures"


def test_renders_pinning_for_one_read(tmp_path):
    out_png = tmp_path / "pinning.png"
    cp = subprocess.run([
        "python3", str(TOOL),
        "--vg-trace-log", str(FIXTURES / "vg_trace_sample.log"),
        "--read-id", "READ_42",
        "--output-png", str(out_png),
    ], capture_output=True, text=True)
    assert cp.returncode == 0, cp.stderr
    assert out_png.exists() and out_png.stat().st_size > 500
```

- [ ] **Step 7.3 — Create stub fixture**

```bash
cat > tools/demo/tests/fixtures/vg_trace_sample.log <<'EOF'
[VG] family_id=1 bundles=2 multimappers=4
[VG] em_iter=0 read_id=READ_42 copy=copy1 weight=0.333
[VG] em_iter=0 read_id=READ_42 copy=copy2 weight=0.333
[VG] em_iter=0 read_id=READ_42 copy=copy3 weight=0.334
[VG] em_iter=2 read_id=READ_42 copy=copy1 weight=0.95
[VG] em_iter=2 read_id=READ_42 copy=copy2 weight=0.03
[VG] em_iter=2 read_id=READ_42 copy=copy3 weight=0.02
EOF
```

> If the actual VG trace format observed in 7.1 differs, edit this fixture to match.

- [ ] **Step 7.4 — Run failing test**

- [ ] **Step 7.5 — Implement renderer**

```python
#!/usr/bin/env python3
"""Multi-mapper read pinning visualizer.

Parses [VG] stderr lines, picks one specified read, plots its weight at
each copy: BEFORE EM (uniform 1/NH) vs AFTER EM convergence.

Usage:
    render_multimapper_pinning.py --vg-trace-log path.log --read-id READ_X --output-png out.png
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt

LINE_RE = re.compile(
    r"^\[VG\]\s+em_iter=(\d+)\s+read_id=(\S+)\s+copy=(\S+)\s+weight=([0-9.eE+-]+)"
)


def parse_log(path: Path, read_id: str) -> dict[int, dict[str, float]]:
    by_iter: dict[int, dict[str, float]] = defaultdict(dict)
    with path.open() as f:
        for line in f:
            m = LINE_RE.match(line.strip())
            if not m:
                continue
            it, rid, copy, w = int(m[1]), m[2], m[3], float(m[4])
            if rid != read_id:
                continue
            by_iter[it][copy] = w
    return dict(by_iter)


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--vg-trace-log", required=True, type=Path)
    ap.add_argument("--read-id", required=True)
    ap.add_argument("--output-png", required=True, type=Path)
    args = ap.parse_args(argv)

    iters = parse_log(args.vg_trace_log, args.read_id)
    if 0 not in iters:
        raise SystemExit(f"no em_iter=0 row for read {args.read_id}")
    last_iter = max(iters.keys())

    pre = iters[0]
    post = iters[last_iter]
    copies = sorted(set(pre) | set(post))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, label, weights in [(axes[0], "Before EM\n(StringTie-equivalent: uniform 1/NH)", pre),
                                (axes[1], f"After EM (iter {last_iter})\nRustle reweighted", post)]:
        vals = [weights.get(c, 0.0) for c in copies]
        bars = ax.bar(copies, vals, color="#1F4E79")
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("read weight")
        ax.set_title(label, fontsize=10)
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.02, f"{v:.2f}",
                    ha="center", fontsize=9)
    fig.suptitle(f"Multi-mapper pinning: read = {args.read_id}", fontsize=12)
    fig.tight_layout()
    fig.savefig(args.output_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
```

- [ ] **Step 7.6 — Verify test passes**

- [ ] **Step 7.7 — Commit**

```bash
git add tools/demo/render_multimapper_pinning.py tools/demo/tests/test_render_multimapper_pinning.py tools/demo/tests/fixtures/vg_trace_sample.log
git commit -m "feat(demo): V4 multi-mapper read pinning visualizer (parses VG stderr trace)"
```

---

### Task 8: V5 — Phase 3 variation-graph sketch

**Files:**
- Create: `tools/demo/render_phase3_vg_sketch.py`

V5 is conceptual — no real data. Hand-coded 2-copy variation graph with bubble nodes at SNP positions, with one read aligning and settling on one path.

- [ ] **Step 8.1 — Implement (no test, manual smoke check is enough)**

```python
#!/usr/bin/env python3
"""Phase 3 conceptual sketch: variation graph for a 2-copy gene family.

Hand-coded illustration of where the project is going. Two copies share a
linear "spine" with bubbles at copy-discriminating SNP positions. A single
read carrying SNP variants is shown aligned to the graph, settling on
copy A's path.

Usage:
    render_phase3_vg_sketch.py --output-png out.png
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as patches


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-png", required=True, type=Path)
    args = ap.parse_args()

    fig, ax = plt.subplots(figsize=(11, 4.5))
    ax.set_xlim(0, 12)
    ax.set_ylim(-2.5, 3)
    ax.axis("off")

    # Spine nodes (shared)
    spine_y = 0
    spine_xs = [1, 3, 5.5, 8, 10.5]
    spine_labels = ["E1", "E2-pre", "E2-post", "E3", "E4"]
    for x, lbl in zip(spine_xs, spine_labels):
        c = patches.FancyBboxPatch((x - 0.45, spine_y - 0.3), 0.9, 0.6,
                                   boxstyle="round,pad=0.02", linewidth=1.2,
                                   edgecolor="#1F4E79", facecolor="#E8F1F8")
        ax.add_patch(c)
        ax.text(x, spine_y, lbl, ha="center", va="center", fontsize=9)

    # Spine edges
    for x1, x2 in zip(spine_xs, spine_xs[1:]):
        ax.annotate("", xy=(x2 - 0.5, 0), xytext=(x1 + 0.5, 0),
                    arrowprops=dict(arrowstyle="->", color="#1F4E79", lw=1.2))

    # Bubble at midpoint between E2-pre and E2-post: SNP G/T
    bx_mid = (spine_xs[1] + spine_xs[2]) / 2
    for dy, base, color, lbl in [(1.1, "G", "#2CA02C", "copy A"), (-1.1, "T", "#D62728", "copy B")]:
        c = patches.Circle((bx_mid, dy), 0.32, edgecolor=color, facecolor="white", linewidth=2)
        ax.add_patch(c)
        ax.text(bx_mid, dy, base, ha="center", va="center", fontsize=10, color=color, fontweight="bold")
        ax.text(bx_mid + 0.6, dy, lbl, ha="left", va="center", fontsize=8, color=color)
    # Bubble edges
    for dy in (1.1, -1.1):
        ax.annotate("", xy=(bx_mid - 0.32, dy * 0.7), xytext=(spine_xs[1] + 0.5, 0),
                    arrowprops=dict(arrowstyle="->", color="#888", lw=0.9))
        ax.annotate("", xy=(spine_xs[2] - 0.5, 0), xytext=(bx_mid + 0.32, dy * 0.7),
                    arrowprops=dict(arrowstyle="->", color="#888", lw=0.9))

    # Bubble at midpoint between E3 and E4: SNP C/A
    bx2 = (spine_xs[3] + spine_xs[4]) / 2
    for dy, base, color in [(1.1, "C", "#2CA02C"), (-1.1, "A", "#D62728")]:
        c = patches.Circle((bx2, dy), 0.32, edgecolor=color, facecolor="white", linewidth=2)
        ax.add_patch(c)
        ax.text(bx2, dy, base, ha="center", va="center", fontsize=10, color=color, fontweight="bold")
    for dy in (1.1, -1.1):
        ax.annotate("", xy=(bx2 - 0.32, dy * 0.7), xytext=(spine_xs[3] + 0.5, 0),
                    arrowprops=dict(arrowstyle="->", color="#888", lw=0.9))
        ax.annotate("", xy=(spine_xs[4] - 0.5, 0), xytext=(bx2 + 0.32, dy * 0.7),
                    arrowprops=dict(arrowstyle="->", color="#888", lw=0.9))

    # Read aligned on copy A path (top)
    read_y = 2.3
    ax.plot([spine_xs[0], spine_xs[-1]], [read_y, read_y], color="#E37222", lw=4, solid_capstyle="round")
    ax.text(spine_xs[0] - 0.6, read_y, "read", ha="right", va="center", fontsize=9, color="#E37222")
    # Mark read SNP calls
    ax.text(bx_mid, read_y - 0.3, "G", ha="center", color="#2CA02C", fontsize=9, fontweight="bold")
    ax.text(bx2, read_y - 0.3, "C", ha="center", color="#2CA02C", fontsize=9, fontweight="bold")
    # Annotate
    ax.text(6, read_y + 0.45, "SNPs G, C → copy A path", ha="center", fontsize=10, color="#2CA02C")

    ax.text(6, -2.2, "Future: align reads to family variation graph; resolve to specific copy by SNP support",
            ha="center", fontsize=10, style="italic", color="#555")

    fig.tight_layout()
    fig.savefig(args.output_png, dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
```

- [ ] **Step 8.2 — Smoke check**

Run: `python3 tools/demo/render_phase3_vg_sketch.py --output-png /tmp/v5_smoke.png && file /tmp/v5_smoke.png`

Expected: PNG image.

- [ ] **Step 8.3 — Commit**

```bash
git add tools/demo/render_phase3_vg_sketch.py
git commit -m "feat(demo): V5 Phase 3 variation-graph conceptual sketch"
```

---

## Phase 3 — Demo orchestration

### Task 9: Locus identification helper

**Files:**
- Create: `tools/demo/identify_locus_bundle.py`

Given a locus (chrom:start-end), this script reads `RUSTLE_PARITY_PARTITION_TSV` and reports the matching bundle ID(s). Used by the orchestrator to know which bundle the chr19 GOLGA6L triple lives in after a real Rustle run.

- [ ] **Step 9.1 — Implement**

```python
#!/usr/bin/env python3
"""Find the Rustle bundle containing a given genomic locus.

Reads RUSTLE_PARITY_PARTITION_TSV and reports bundles whose interval overlaps
the requested locus.

Usage:
    identify_locus_bundle.py --partition-tsv path.tsv --locus chrom:start-end
"""
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--partition-tsv", required=True, type=Path)
    ap.add_argument("--locus", required=True)
    args = ap.parse_args()

    chrom, rng = args.locus.split(":")
    lstart, lend = map(int, rng.split("-"))

    with args.partition_tsv.open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row.get("chrom") != chrom:
                continue
            bs, be = int(row["start"]), int(row["end"])
            if be >= lstart and bs <= lend:
                print(f"{row['chrom']}\t{bs}\t{be}\t{row.get('strand','.')}\t{row.get('signature','-')[:80]}")

if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 9.2 — Commit**

```bash
git add tools/demo/identify_locus_bundle.py
git commit -m "feat(demo): identify_locus_bundle helper"
```

---

### Task 10: Tier 1 orchestrator (synthetic family)

**Files:**
- Create: `tools/demo/run_demo_synthetic.sh`

- [ ] **Step 10.1 — Implement**

```bash
#!/usr/bin/env bash
# Tier 1 — synthetic_family walkthrough.
# Produces all dumps + V1 + V3 + V4 PNGs for the 95-read synthetic dataset.
set -euo pipefail

cd "$(dirname "$0")/../.."

OUT="tools/demo/out/synthetic"
mkdir -p "$OUT"

BAM="test_data/synthetic_family/reads_sorted.bam"
GTF_TRUTH="test_data/synthetic_family/truth.gtf"

GTF_BASELINE="$OUT/baseline.gtf"
GTF_VG="$OUT/vg.gtf"

NODES_TSV="$OUT/nodes.tsv"
EDGES_TSV="$OUT/edges.tsv"
FLOW_TSV="$OUT/flow.tsv"
READS_TSV="$OUT/reads.tsv"
PARTITION_TSV="$OUT/partition.tsv"
VG_LOG="$OUT/vg_trace.log"

# clean prior dumps so headers are fresh
rm -f "$NODES_TSV" "$EDGES_TSV" "$FLOW_TSV" "$READS_TSV" "$PARTITION_TSV"

echo "== Running Rustle baseline (no VG) =="
RUSTLE_PARITY_GRAPH_TSV="$NODES_TSV" \
RUSTLE_PARITY_GRAPH_EDGES_TSV="$EDGES_TSV" \
RUSTLE_FLOW_ITER_TSV="$FLOW_TSV" \
RUSTLE_PARITY_READ_DUMP_TSV="$READS_TSV" \
RUSTLE_PARITY_PARTITION_TSV="$PARTITION_TSV" \
./target/release/rustle -L -o "$GTF_BASELINE" "$BAM"

echo "== Running Rustle --vg (capturing trace) =="
RUSTLE_VG_TRACE=1 \
./target/release/rustle -L --vg -o "$GTF_VG" "$BAM" 2> "$VG_LOG"

echo "== Per-stage dump line counts =="
wc -l "$READS_TSV" "$PARTITION_TSV" "$NODES_TSV" "$EDGES_TSV" "$FLOW_TSV" "$GTF_BASELINE" "$GTF_VG"

echo "== Rendering V1 splice graph (first bundle) =="
FIRST_LOCUS=$(awk -F'\t' 'NR>1 {print $1":"$3"-"$4; exit}' "$NODES_TSV")
echo "  Locus = $FIRST_LOCUS"
python3 tools/demo/render_splice_graph.py \
  --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$FIRST_LOCUS" \
  --output-png "$OUT/v1_splice_graph.png" \
  --output-dot "$OUT/v1_splice_graph.dot" \
  --title "Tier 1 — synthetic_family bundle 1"

echo "== Rendering V3 flow iterations =="
python3 tools/demo/render_flow_iterations.py \
  --flow-tsv "$FLOW_TSV" --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$FIRST_LOCUS" \
  --output-png "$OUT/v3_flow_iterations.png"

echo "== Rendering V4 multi-mapper pinning (pick a multi-mapper) =="
# Read IDs in synthetic_family include "MULTI_*" pattern for cross-copy reads.
READ_ID=$(grep -m1 -oE 'read_id=MULTI_[A-Z0-9_]+' "$VG_LOG" | head -1 | cut -d= -f2 || echo "")
if [[ -n "$READ_ID" ]]; then
  python3 tools/demo/render_multimapper_pinning.py \
    --vg-trace-log "$VG_LOG" --read-id "$READ_ID" \
    --output-png "$OUT/v4_multimapper_${READ_ID}.png" || echo "WARN: V4 render failed"
else
  echo "WARN: no MULTI_* read found in VG trace; skipping V4"
fi

echo "== Tier 1 complete. Outputs in $OUT =="
ls -lh "$OUT"
```

Make executable.

```bash
chmod +x tools/demo/run_demo_synthetic.sh
```

- [ ] **Step 10.2 — Run end-to-end**

```bash
bash tools/demo/run_demo_synthetic.sh 2>&1 | tail -30
```

Expected: all sections complete; `tools/demo/out/synthetic/` populated with TSVs, GTFs, and 3 PNGs (V1, V3, V4).

- [ ] **Step 10.3 — Commit**

```bash
git add tools/demo/run_demo_synthetic.sh
git commit -m "feat(demo): Tier 1 synthetic-family orchestrator script"
```

---

### Task 11: Tier 2 orchestrator (chr19 GOLGA6L triple)

**Files:**
- Create: `tools/demo/run_demo_chr19_triple.sh`

- [ ] **Step 11.1 — Implement**

```bash
#!/usr/bin/env bash
# Tier 2 — chr19 GOLGA6L paralog triple from GGO_19.bam.
# Three loci on NC_073243.2: LOC115931294, LOC134757625, LOC101137218.
set -euo pipefail

cd "$(dirname "$0")/../.."

OUT="tools/demo/out/chr19_triple"
mkdir -p "$OUT"

BAM="/mnt/c/Users/jfris/Desktop/GGO_19.bam"
LOCUS_REGION="NC_073243.2:104780000-104900000"

LOCUS_A="NC_073243.2:104789647-104796276"   # LOC115931294
LOCUS_B="NC_073243.2:104830536-104837094"   # LOC134757625
LOCUS_C="NC_073243.2:104871356-104877901"   # LOC101137218

# slice BAM to the region for fast runs
REGION_BAM="$OUT/region.bam"
samtools view -bh "$BAM" "$LOCUS_REGION" > "$REGION_BAM"
samtools index "$REGION_BAM"

NODES_TSV="$OUT/nodes.tsv"
EDGES_TSV="$OUT/edges.tsv"
FLOW_TSV="$OUT/flow.tsv"
READS_TSV="$OUT/reads.tsv"
PARTITION_TSV="$OUT/partition.tsv"
VG_LOG="$OUT/vg_trace.log"

GTF_BASELINE="$OUT/baseline.gtf"
GTF_VG="$OUT/vg.gtf"

rm -f "$NODES_TSV" "$EDGES_TSV" "$FLOW_TSV" "$READS_TSV" "$PARTITION_TSV"

echo "== Rustle baseline on chr19 triple region =="
RUSTLE_PARITY_GRAPH_TSV="$NODES_TSV" \
RUSTLE_PARITY_GRAPH_EDGES_TSV="$EDGES_TSV" \
RUSTLE_FLOW_ITER_TSV="$FLOW_TSV" \
RUSTLE_PARITY_READ_DUMP_TSV="$READS_TSV" \
RUSTLE_PARITY_PARTITION_TSV="$PARTITION_TSV" \
./target/release/rustle -L -o "$GTF_BASELINE" "$REGION_BAM"

echo "== Rustle --vg on same region =="
RUSTLE_VG_TRACE=1 \
./target/release/rustle -L --vg -o "$GTF_VG" "$REGION_BAM" 2> "$VG_LOG"

echo "== Bundles overlapping each paralog locus =="
for L in "$LOCUS_A" "$LOCUS_B" "$LOCUS_C"; do
  echo "-- $L"
  python3 tools/demo/identify_locus_bundle.py --partition-tsv "$PARTITION_TSV" --locus "$L"
done

echo "== Rendering V1 individual splice graphs =="
python3 tools/demo/render_splice_graph.py --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$LOCUS_A" --output-png "$OUT/v1_locA.png" --output-dot "$OUT/v1_locA.dot" \
  --title "LOC115931294"
python3 tools/demo/render_splice_graph.py --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$LOCUS_B" --output-png "$OUT/v1_locB.png" --output-dot "$OUT/v1_locB.dot" \
  --title "LOC134757625"
python3 tools/demo/render_splice_graph.py --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$LOCUS_C" --output-png "$OUT/v1_locC.png" --output-dot "$OUT/v1_locC.dot" \
  --title "LOC101137218"

echo "== Rendering V2 three-way isomorphism =="
python3 tools/demo/compare_graphs_isomorphism.py \
  --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$LOCUS_A" --locus "$LOCUS_B" --locus "$LOCUS_C" \
  --output-json "$OUT/v2_isomorphism.json" \
  --output-png "$OUT/v2_isomorphism.png"

echo "  Iso result:"
python3 -c "import json,sys; d=json.load(open('$OUT/v2_isomorphism.json')); print('  isomorphic=%s' % d['isomorphic']); print('  pairwise=%s' % list(d['pairwise'].keys()))"

echo "== Rendering V3 flow iterations for locus A =="
python3 tools/demo/render_flow_iterations.py \
  --flow-tsv "$FLOW_TSV" --nodes-tsv "$NODES_TSV" --edges-tsv "$EDGES_TSV" \
  --locus "$LOCUS_A" --output-png "$OUT/v3_flow_locA.png"

echo "== Rendering V4 multi-mapper pinning (cross-copy read) =="
READ_ID=$(grep -m1 -oE 'read_id=\S+' "$VG_LOG" | head -1 | cut -d= -f2 || echo "")
if [[ -n "$READ_ID" ]]; then
  python3 tools/demo/render_multimapper_pinning.py \
    --vg-trace-log "$VG_LOG" --read-id "$READ_ID" \
    --output-png "$OUT/v4_multimapper_${READ_ID}.png" || echo "WARN: V4 render failed"
fi

echo "== Side-by-side baseline-vs-VG transcript counts =="
echo "Baseline transcripts:"
grep -cE '^[^#].*\ttranscript\t' "$GTF_BASELINE" || true
echo "VG transcripts:"
grep -cE '^[^#].*\ttranscript\t' "$GTF_VG" || true

echo "== Tier 2 complete. Outputs in $OUT =="
ls -lh "$OUT"
```

- [ ] **Step 11.2 — Run end-to-end**

```bash
chmod +x tools/demo/run_demo_chr19_triple.sh
bash tools/demo/run_demo_chr19_triple.sh 2>&1 | tail -40
```

Expected: TSVs populated; V1 (×3 individual paralog graphs), V2 (3-panel isomorphism), V3 (flow grid for locus A), V4 (multi-mapper) all produced.

> **If V2 reports `isomorphic: false`**, that's a real finding. Read `v2_isomorphism.json` to see which pair fails and why. We adjust the meeting narrative honestly: "graphs of LOC115931294 and LOC134757625 are isomorphic, but LOC101137218 has an extra junction at position X due to read evidence." This is still defensible — paralogs *should* have the same architecture but assemblies can drift if read coverage differs.

- [ ] **Step 11.3 — Commit**

```bash
git add tools/demo/run_demo_chr19_triple.sh
git commit -m "feat(demo): Tier 2 chr19 GOLGA6L triple orchestrator script"
```

---

## Phase 4 — Final polish for the meeting

### Task 12: Pre-rehearsal + fallback recording

- [ ] **Step 12.1 — Run both orchestrators end-to-end**

```bash
bash tools/demo/run_demo_synthetic.sh 2>&1 | tee /tmp/synth_run.log
bash tools/demo/run_demo_chr19_triple.sh 2>&1 | tee /tmp/chr19_run.log
```

- [ ] **Step 12.2 — Open each PNG and sanity-check**

For each of: `tools/demo/out/synthetic/{v1,v3,v4}_*.png` and `tools/demo/out/chr19_triple/{v1,v2,v3,v4}_*.png` — open in image viewer, confirm:
- Nodes are readable
- Edges have visible labels with capacities
- V2 panel shows 3 graphs side by side; if isomorphic, check that the title/caption mentions the bijection
- V3 grid shows multiple iterations with augmenting paths highlighted in orange
- V4 bar chart shows distinct before/after weight distributions

- [ ] **Step 12.3 — Write the runbook**

```markdown
<!-- tools/demo/runbook.md -->
# Meeting runbook

Total: ~25 minutes. Order matters — credibility first, results second, vision third.

## Setup (before he walks in)
- Terminal at `/mnt/c/Users/jfris/Desktop/Rustle/`
- Image viewer pointed at `tools/demo/out/`
- `git log --oneline -5` shows recent commits with timestamps (proves no last-minute hardcoding)

## Tier 1 — Synthetic walkthrough (5 min)
1. `cat test_data/synthetic_family/truth.gtf` — *"Ground truth, defined before reads existed."*
2. `head -30 test_data/synthetic_family/generate_isoseq.py` — *"Generator with `random.seed(42)`. The BAM is a deterministic function of these two files."*
3. `samtools view test_data/synthetic_family/reads_sorted.bam | head -5` — *"5 of 95 reads."*
4. `bash tools/demo/run_demo_synthetic.sh` — runs in seconds.
5. `cat tools/demo/out/synthetic/reads.tsv | column -t -s$'\t' | head -10` — *"Reads as parsed by Rustle."*
6. `cat tools/demo/out/synthetic/partition.tsv | head -5` — *"Bundles formed from the reads."*
7. Show V1 PNG — *"Splice graph from those bundles. Nodes = exons, edges = collinear or junctions."*
8. Show V3 PNG — *"Edmonds-Karp peels off augmenting paths one at a time. Each path becomes a transcript."*
9. `diff <(grep transcript tools/demo/out/synthetic/baseline.gtf | awk '{print $4,$5}') <(grep transcript test_data/synthetic_family/truth.gtf | awk '{print $4,$5}')` — *"Recovers all 3 truth transcripts."*

**Alternative splicing case:** Point out that MYG1_A's A2 transcript skips exon 3 — that creates the "skip" edge in the V1 graph (node 1 → node 3 directly). The flow decomposition (V3) extracts this path on iteration 2.

## Tier 2 — Real data: chr19 GOLGA6L triple (10 min)
1. `samtools view /mnt/c/Users/jfris/Desktop/GGO_19.bam NC_073243.2:104789647-104877901 | wc -l` — *"~N reads in this 88kb region."*
2. `bash tools/demo/run_demo_chr19_triple.sh`
3. Show V1 ×3 PNGs side by side — *"Three paralog splice graphs from real data. Same shape?"*
4. Show V2 PNG — *"VF2 isomorphism check confirms it. Bijection in v2_isomorphism.json."*
5. `cat tools/demo/out/chr19_triple/v2_isomorphism.json | head -30` — *"Pairwise isomorphism, with bijection mappings."*
6. Show V3 PNG (locus A) — *"Same Edmonds-Karp on real data."*
7. Show V4 PNG — *"This multi-mapper read crosses 2 of the 3 paralogs. StringTie distributes weight uniformly. Rustle's EM converges to a confident assignment."*
8. Compare baseline vs VG transcript counts at end of orchestrator log.

## Tier 3 — Closing slide (3 min)
Show the cross-chromosome GOLGA family table from `MULTI_COPY_FAMILY_PROOF.md`:
- 4 families, 52 loci, +3 exact / -4 wrong-strand / 1 fewer missing under VG
- Then 5 chromosomes spanned by the GOLGA family.

## Phase 3 — Vision (3 min)
Show V5 PNG. *"Where this goes: a single variation graph per family, with bubbles at SNP positions. Reads aligned to the graph resolve to a specific copy via SNP support. This handles unmapped reads, novel paralogs, and copy-level read assignment that's currently impossible because StringTie doesn't even use multi-mappers in long-read mode."*

## Q&A backup
- "Where's the assembly code?" — *"I'd rather show you what it does. Pick another locus and I'll re-run."*
- "Did you hardcode the result?" — *"The synthetic generator is deterministic from seed=42. Re-run it. Same BAM, same answer."*
- "Why don't the iso panels look identical?" — *"Different absolute positions, same topology. The bijection in the JSON shows the node correspondence."*
- "What's missing?" — *"Phase 3 isn't built. SNP-based copy assignment is scaffolded, not finished. Cross-chromosome novel-copy discovery is in scaffold."*

## Crash plan
If a live run fails:
- All TSVs and PNGs are pre-cached in `tools/demo/out/{synthetic,chr19_triple}/`. Show those instead.
- Pre-recorded terminal: `asciinema rec /tmp/demo.cast` while running each orchestrator the night before. Replay if needed.
```

- [ ] **Step 12.4 — Pre-record fallback terminal session**

```bash
asciinema rec tools/demo/out/synthetic_demo.cast --command "bash tools/demo/run_demo_synthetic.sh"
asciinema rec tools/demo/out/chr19_triple_demo.cast --command "bash tools/demo/run_demo_chr19_triple.sh"
```

> Skip if `asciinema` not installed; alternative: `script` command to capture terminal output.

- [ ] **Step 12.5 — Commit**

```bash
git add tools/demo/runbook.md tools/demo/out/
git commit -m "docs(demo): runbook + fallback recordings for advisor meeting"
```

---

## Self-Review

After writing this plan, I checked:

**Spec coverage:**
- ✅ Tier 1 synthetic walkthrough — Tasks 3-10
- ✅ Tier 2 chr19 GOLGA6L live demo — Tasks 9-11
- ✅ Tier 3 closing slide — Task 12.3 runbook references the existing `MULTI_COPY_FAMILY_PROOF.md` table
- ✅ V1 — Task 4
- ✅ V2 — Task 5
- ✅ V3 — Task 6
- ✅ V4 — Task 7
- ✅ V5 — Task 8
- ✅ Graph edges dump — Task 1
- ✅ Flow iter dump — Task 2
- ✅ Tools dir under `tools/demo/` — Task 3
- ✅ Alternative splicing handled — runbook Task 12.3 explicitly highlights MYG1_A's exon-3-skip
- ✅ "Similar" definition (iso preserving node count + edge multiset by kind + source/sink topology, normalized coords) — Task 5 implementation matches spec

**Placeholder scan:** none.

**Type consistency:** GraphNode/GraphEdge/FlowIter dataclasses are defined in common.py and imported by all renderers; field names match. Env-var names are consistent: `RUSTLE_PARITY_GRAPH_EDGES_TSV` and `RUSTLE_FLOW_ITER_TSV`. TSV column lists match between Rust emit code and Python read code.

**Time check:** Task 1 (~75 min) + Task 2 (~120 min) + Tasks 3-8 (~5-6 hours) + Tasks 9-11 (~2 hours) + Task 12 (~1.5 hours) = ~12 hours total. Tight but achievable. Drop priority if compressed: V5 → V4 → fallback recording → Tier 3 closing slide.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-04-26-rustle-multicopy-meeting-prep-plan.md`. Two execution options:

1. **Subagent-Driven (recommended)** — I dispatch a fresh subagent per task, review between tasks, fast iteration. Good for time-constrained work because parallel-eligible tasks (V4, V5, V8 fixtures) can run concurrently.

2. **Inline Execution** — Execute tasks in this session using executing-plans, batch execution with checkpoints. More linear; main thread does the work.

Which approach?

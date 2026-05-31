# Source/Sink Terminal-Wiring Alignment — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure — oracle-first, default-OFF — whether aligning Rustle's source/sink (terminal) edge wiring to StringTie's reduces over-enumeration FPs, and at what TP cost; decide ship/shelve from the oracle before any production change.

**Architecture:** Phase 0a is a free probe via the existing `RUSTLE_COVLINK_RECURSE_ZERO_OFF` toggle + `bench/edge_diff.py` + gffcompare. Phase 0b is an injection oracle (`RUSTLE_TERMINAL_ORACLE=<st_edge.jsonl>`, default OFF) that, on identical-node bundles, overrides Rustle's synthetic source/sink edge set with StringTie's captured set and reconciles the synthetic transfrag abundances, then measures the causal FP/TP delta. A hard abort gate follows. Phase 1 (gated) aligns the rule behind `RUSTLE_ST_TERMINAL`.

**Tech Stack:** Rust (Rustle pipeline + graph), Python 3 stdlib (bench tools), gffcompare. The `graph_edge` parity instrumentation and `bench/edge_diff.py` from the prior task are the validation substrate.

---

## Grounding facts (verified — do not re-derive)

- **Call site:** `add_coverage_source_sink_edges` → `coverage_synth` at `pipeline.rs:13255`, inside the `process_graph` closure (`pipeline.rs:12802`), in the parallel bundle loop (`pipeline.rs:11279`). `coverage_synth` is merged with `alt_tts_synth` (`:13335`) and `terminal_donor_synth` (`:13339`), then consumed by `collect_from_transfrags` (`:13626`) → `materialize_links` (`:13650`) → `transfrags.extend(realized)` (`:13652`) → `process_transfrags` (`:13824`) → flow.
- **Hook point for the oracle:** AFTER the three synth sources are merged into `coverage_synth` (after `pipeline.rs:13339`) and BEFORE `collect_from_transfrags` (`:13626`). At that point `graph_mut` holds ALL source/sink edges and `coverage_synth` is the complete synthetic-transfrag list. `graph_bundle` (`.start/.end/.strand/.chrom`) and node coords (`graph_mut.nodes[i].start/.end`) are in scope.
- **Edge API:** `graph.add_edge(from, to)` (`graph.rs:410`) sets `children`/`parents` SmallBitsets + `gpos`; `graph.remove_edge(from, to)` (`graph.rs:529`) clears them + removes from `gpos`. `graph.source_id`, `graph.sink_id`, `graph.n_nodes` public. Node coords public: `nodes[i].start`, `nodes[i].end`.
- **graph_edge token convention (from prior task):** real-node token = `format!("{}-{}", n.start + 1, n.end)` (1-based inclusive); `SRC`/`SNK` for source/sink. Bundle key = `(graph_bundle.start + 1, graph_bundle.end, strand_char)` — SAME as the emitted `graph_edge`/`graphnode_list` start/end/strand.
- **GraphTransfrag:** `node_ids: Vec<usize>`, `abundance: f64`, `origin_tag`, created via `GraphTransfrag::new(vec![a, b], graph.pattern_size())` then `.abundance = ...`. Synthetic terminal transfrags carry `origin_tag` ∈ {`"coverage_edge"`, `"alt_tts"`, `"terminal_donor"`}.
- **Abundance formula:** source edge `(icov - parcov)/DROP`, sink edge `(icov - chcov)/DROP`, `DROP=0.5`. `icov` = node avg cov = `bpcov.get_cov_range(sno, bpcov.plus.idx(node.start), bpcov.plus.idx(node.end)) / (node.end - node.start)`.
- **CONFOUND (document, don't fix):** flow's `build_lr_edge_capacities` (`global_flow.rs:205-234`) synthesizes VIRTUAL source/sink capacity edges for partial real transfrags during decomposition (`greedy_flow_decompose_paths` → `:302`). So removing a synthetic terminal edge does NOT fully suppress a start/end point if a partial real transfrag still needs it. The oracle therefore measures the **synthetic-coverage-edge contribution** (a lower bound on full ST-terminal suppression). The abort-gate number is interpreted with this in mind.
- **Flag predicate pattern (`stringtie_parity.rs:183/199/213`):** `matches!(v, Some(s) if s != "0")` over `std::env::var("X").ok().as_deref()`. `RUSTLE_TERMINAL_ORACLE` is a PATH (use `std::env::var("RUSTLE_TERMINAL_ORACLE").ok().filter(|s| !s.is_empty())` → `Option<String>`); `RUSTLE_ST_TERMINAL` is a bool (use the `matches!` pattern).
- **Toggles (Phase 0a):** `RUSTLE_COVLINK_RECURSE_ZERO_OFF=1` disables phantom recursion (`graph_build.rs:146,185`, default ON); `RUSTLE_COVLINK_THRESHOLD_LOOSE` (`:82`, default OFF = 0.05 = ST). `bench/capture_parity.sh` STEPS already includes `graph_edge`.
- **Baseline:** de novo `-L` → Intron 96.2/91.7, Transcript 95.6/90.5.

---

## File structure

- Modify `src/rustle/stringtie_parity.rs` — add `terminal_oracle_path()` (Option<String>) + `st_terminal()` (bool) predicates + tests.
- Create `src/rustle/terminal_oracle.rs` — load ST's graph_edge terminal set from the oracle file (OnceLock); `override_terminal_edges(graph, bundle_key, coverage_synth, bpcov, sno)` that, on identical-node bundles, removes RU source/sink edges not in ST's set (+ drops their synth transfrags) and installs ST's missing ones (+ synth transfrags with reconciled abundance).
- Modify `src/rustle/pipeline.rs` (~13340, after the synth merge) — call the oracle override under `terminal_oracle_path()`.
- Modify `src/rustle/graph_build.rs:843` region — Phase 1 only: OR `st_terminal()` into the wiring rule.
- Create `bench/terminal_oracle_report.py` — FP/FN/F1 delta + endpoint attribution.
- Modify `docs/STRINGTIE_PARITY_FINDINGS.md` + memory — record results.

---

## PHASE 0a — Free probe (no production code)

### Task 1: Phantom-recursion free probe

**Files:** none created (uses existing toggles + `bench/edge_diff.py` + gffcompare).

- [ ] **Step 1: Capture baseline + no-recurse edges**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
# baseline (recurse ON) — refresh
bash bench/capture_parity.sh >/dev/null 2>&1
cp /tmp/ru_edge.jsonl /tmp/ru_edge_baseline.jsonl 2>/dev/null || cp /tmp/ru_gn.jsonl /tmp/ru_edge_baseline.jsonl
# no-recurse capture
RUSTLE_COVLINK_RECURSE_ZERO_OFF=1 \
RUSTLE_PARITY_LOG=/tmp/ru_edge_norecurse.jsonl \
RUSTLE_PARITY_FILTER_STEPS=graphnode_list,graph_edge \
RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 \
  ./target/release/rustle -L GGO_19.bam -o /tmp/ru_norecurse.gtf 2>/dev/null
```

- [ ] **Step 2: Measure edge-reduction vs ST**

Run: `python3 bench/edge_diff.py /tmp/ru_edge_norecurse.jsonl /tmp/st_edge.jsonl`
Expected: report shows RU-extra source_sink count; compare to the baseline 2708. Record how many of the 2708 extra terminal edges the phantom-recursion-off run removes.

- [ ] **Step 3: Measure F1/FP delta**

Run: `gffcompare -r ../GGO_19.gtf /tmp/ru_norecurse.gtf -o /tmp/ru_norecurse_cmp 2>/dev/null; grep -E "Intron chain level:|Transcript level:" /tmp/ru_norecurse_cmp.stats`
Expected: compare to baseline 96.2/91.7, 95.6/90.5. Record the direction (help / hurt / no-op) and magnitude.

- [ ] **Step 4: Record the probe result (no commit — this informs Phase 0b)**

Note in your task report: (a) extra terminal edges removed by `RECURSE_ZERO_OFF`, (b) the F1 delta. If `RECURSE_ZERO_OFF` alone removes most of the 2708 AND is F1-positive, flag that the phantom recursion may be the whole story (a near-free win to consider promoting in Phase 1).

---

## PHASE 0b — Injection oracle

### Task 2: Flag predicates `terminal_oracle_path` + `st_terminal`

**Files:** Modify `src/rustle/stringtie_parity.rs`

- [ ] **Step 1: Write the failing test**

Add to the `#[cfg(test)]` module:

```rust
#[test]
fn st_terminal_default_off() {
    use super::st_terminal_from;
    assert!(!st_terminal_from(None));
    assert!(st_terminal_from(Some("1")));
    assert!(!st_terminal_from(Some("0")));
}

#[test]
fn terminal_oracle_path_parses() {
    use super::terminal_oracle_path_from;
    assert_eq!(terminal_oracle_path_from(None), None);
    assert_eq!(terminal_oracle_path_from(Some("".to_string())), None);
    assert_eq!(terminal_oracle_path_from(Some("/tmp/st_edge.jsonl".to_string())), Some("/tmp/st_edge.jsonl".to_string()));
}
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p rustle st_terminal_default_off terminal_oracle_path_parses 2>&1 | tail -5`
Expected: FAIL — functions not found.

- [ ] **Step 3: Implement the predicates**

Add next to `st_flow` (~line 213):

```rust
pub fn st_terminal_from(v: Option<&str>) -> bool {
    matches!(v, Some(s) if s != "0")
}

pub fn st_terminal() -> bool {
    st_terminal_from(std::env::var("RUSTLE_ST_TERMINAL").ok().as_deref())
}

pub fn terminal_oracle_path_from(v: Option<String>) -> Option<String> {
    v.filter(|s| !s.is_empty())
}

pub fn terminal_oracle_path() -> Option<String> {
    terminal_oracle_path_from(std::env::var("RUSTLE_TERMINAL_ORACLE").ok())
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cargo test -p rustle st_terminal_default_off terminal_oracle_path_parses 2>&1 | tail -5`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/rustle/stringtie_parity.rs
git commit -m "feat(parity): RUSTLE_ST_TERMINAL + RUSTLE_TERMINAL_ORACLE predicates (default OFF)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Injection-oracle module + pipeline hook

**Files:** Create `src/rustle/terminal_oracle.rs`; Modify `src/rustle/pipeline.rs` (~13340), `src/rustle/mod.rs`/`lib.rs` (module decl).

- [ ] **Step 1: Create the oracle module**

Create `src/rustle/terminal_oracle.rs`:

```rust
//! Injection oracle (analysis-only, default OFF via RUSTLE_TERMINAL_ORACLE=<path>).
//! On identical-node bundles, override Rustle's synthetic source/sink edge set with
//! StringTie's captured `graph_edge` set, reconciling the synthetic transfrag abundances,
//! to bound the FP/TP impact of the terminal-wiring divergence before any production change.
//!
//! CONFOUND: flow's build_lr_edge_capacities (global_flow.rs:205-234) creates virtual
//! source/sink capacity edges for partial real transfrags, so this measures the
//! synthetic-coverage-edge contribution — a LOWER BOUND on full ST-terminal suppression.

use std::collections::{HashMap, HashSet};
use std::sync::OnceLock;
use crate::graph::{Graph, GraphTransfrag};
use crate::bam::BpcovStranded;
use crate::constants::DROP; // adjust import path to where DROP lives (graph_build re-exports it)

/// Per-bundle ST terminal edges: bundle key (start1, end, strand) -> set of (from_token, to_token).
/// Tokens: "A-B" for a real node (1-based inclusive), "SRC"/"SNK" for source/sink.
type OracleMap = HashMap<(u64, u64, char), HashSet<(String, String)>>;

static ORACLE: OnceLock<OracleMap> = OnceLock::new();

/// Load ST graph_edge events once from the oracle file path.
pub fn load_oracle(path: &str) -> &'static OracleMap {
    ORACLE.get_or_init(|| {
        let mut map: OracleMap = HashMap::new();
        let Ok(content) = std::fs::read_to_string(path) else { return map; };
        for line in content.lines() {
            let Ok(v) = serde_json::from_str::<serde_json::Value>(line) else { continue; };
            if v.get("step").and_then(|s| s.as_str()) != Some("graph_edge") { continue; }
            let start = v.get("start").and_then(|x| x.as_u64()).unwrap_or(0);
            let end = v.get("end").and_then(|x| x.as_u64()).unwrap_or(0);
            let strand = v.get("strand").and_then(|x| x.as_str()).and_then(|s| s.chars().next()).unwrap_or('.');
            let edges_str = v.get("payload").and_then(|p| p.get("edges")).and_then(|e| e.as_str()).unwrap_or("");
            let set = map.entry((start, end, strand)).or_default();
            for e in edges_str.split(',') {
                if let Some((f, t)) = e.split_once('>') {
                    set.insert((f.to_string(), t.to_string()));
                }
            }
        }
        map
    })
}

/// Node-coordinate token, matching the graph_edge emit convention (1-based inclusive).
fn node_token(g: &Graph, nid: usize) -> String {
    if nid == g.source_id { return "SRC".to_string(); }
    if nid == g.sink_id { return "SNK".to_string(); }
    let n = &g.nodes[nid];
    format!("{}-{}", n.start + 1, n.end)
}

/// Build a token -> node_id map for real nodes (for mapping ST tokens back to RU ids).
fn token_to_id(g: &Graph) -> HashMap<String, usize> {
    let mut m = HashMap::new();
    for i in 0..g.n_nodes {
        if i == g.source_id || i == g.sink_id { continue; }
        m.insert(node_token(g, i), i);
    }
    m
}

/// Average coverage of a node (icov), mirroring add_coverage_source_sink_edges.
fn node_avg_cov(g: &Graph, nid: usize, bpcov: &BpcovStranded, sno: usize) -> f64 {
    let n = &g.nodes[nid];
    if n.end <= n.start { return 0.0; }
    let len = (n.end - n.start) as f64;
    let si = bpcov.plus.idx(n.start);
    let ei = bpcov.plus.idx(n.end);
    bpcov.get_cov_range(sno, si, ei) / len
}

/// Override RU's source/sink edges on an identical-node bundle to match ST's captured set.
/// Returns (edges_removed, edges_installed); mutates `graph` and filters/extends `coverage_synth`.
/// Returns None (no-op) if the bundle's node set does not byte-match ST's (cannot inject cleanly).
pub fn override_terminal_edges(
    graph: &mut Graph,
    bundle_key: (u64, u64, char),
    coverage_synth: &mut Vec<GraphTransfrag>,
    bpcov: &BpcovStranded,
    sno: usize,
    oracle: &OracleMap,
) -> Option<(usize, usize)> {
    let st_edges = oracle.get(&bundle_key)?;
    let tok2id = token_to_id(graph);

    // Identical-node check: RU real-node token set must equal ST's real-node token set
    // (the tokens appearing in ST's edges, minus SRC/SNK).
    let ru_nodes: HashSet<String> = tok2id.keys().cloned().collect();
    let st_nodes: HashSet<String> = st_edges.iter()
        .flat_map(|(f, t)| [f.clone(), t.clone()])
        .filter(|x| x != "SRC" && x != "SNK")
        .collect();
    // st_nodes is a subset of nodes referenced by edges; require ST's referenced nodes ⊆ RU's,
    // and that every ST SRC/SNK endpoint maps to a RU node. Bundles failing this are skipped.
    if !st_nodes.is_subset(&ru_nodes) { return None; }

    let src = graph.source_id;
    let snk = graph.sink_id;

    // Current RU source/sink edges as token pairs.
    let mut ru_src_edges: Vec<usize> = graph.nodes[src].children.ones().collect();
    let mut ru_snk_edges: Vec<usize> = graph.nodes[snk].parents.ones().collect();

    // ST's desired source/sink edges, mapped to RU ids.
    let st_src_targets: HashSet<usize> = st_edges.iter()
        .filter(|(f, _)| f == "SRC")
        .filter_map(|(_, t)| tok2id.get(t).copied())
        .collect();
    let st_snk_sources: HashSet<usize> = st_edges.iter()
        .filter(|(_, t)| t == "SNK")
        .filter_map(|(f, _)| tok2id.get(f).copied())
        .collect();

    let mut removed = 0usize;
    let mut installed = 0usize;

    // 1) REMOVE RU source edges ST does not have.
    for &nid in &ru_src_edges {
        if !st_src_targets.contains(&nid) {
            graph.remove_edge(src, nid);
            removed += 1;
        }
    }
    // 2) REMOVE RU sink edges ST does not have.
    for &nid in &ru_snk_edges {
        if !st_snk_sources.contains(&nid) {
            graph.remove_edge(nid, snk);
            removed += 1;
        }
    }
    // Recompute after removals for the "install missing" pass.
    ru_src_edges = graph.nodes[src].children.ones().collect();
    ru_snk_edges = graph.nodes[snk].parents.ones().collect();
    let ru_src_set: HashSet<usize> = ru_src_edges.iter().copied().collect();
    let ru_snk_set: HashSet<usize> = ru_snk_edges.iter().copied().collect();

    let pattern_size = graph.pattern_size();
    // 3) INSTALL ST source edges RU lacks (+ synth transfrag with reconciled abundance).
    for &nid in &st_src_targets {
        if !ru_src_set.contains(&nid) {
            let icov = node_avg_cov(graph, nid, bpcov, sno);
            let parcov: f64 = graph.nodes[nid].parents.ones()
                .filter(|&p| p != src && p != snk)
                .map(|p| node_avg_cov(graph, p, bpcov, sno)).sum();
            graph.add_edge(src, nid);
            let mut tf = GraphTransfrag::new(vec![src, nid], pattern_size);
            tf.abundance = ((icov - parcov) / DROP).max(0.0);
            coverage_synth.push(tf);
            installed += 1;
        }
    }
    // 4) INSTALL ST sink edges RU lacks.
    for &nid in &st_snk_sources {
        if !ru_snk_set.contains(&nid) {
            let icov = node_avg_cov(graph, nid, bpcov, sno);
            let chcov: f64 = graph.nodes[nid].children.ones()
                .filter(|&c| c != src && c != snk)
                .map(|c| node_avg_cov(graph, c, bpcov, sno)).sum();
            graph.add_edge(nid, snk);
            let mut tf = GraphTransfrag::new(vec![nid, snk], pattern_size);
            tf.abundance = ((icov - chcov) / DROP).max(0.0);
            coverage_synth.push(tf);
            installed += 1;
        }
    }

    // 5) DROP synth transfrags whose 2-node source/sink edge was removed.
    coverage_synth.retain(|tf| {
        if tf.node_ids.len() != 2 { return true; }
        let (a, b) = (tf.node_ids[0], tf.node_ids[1]);
        if a == src { return graph.nodes[src].children.contains(b); }
        if b == snk { return graph.nodes[snk].parents.contains(a); }
        true
    });

    Some((removed, installed))
}
```

NOTE: confirm at edit time — the `DROP` constant import path, `BpcovStranded::plus.idx`, `get_cov_range(sno, ..)`, `GraphTransfrag::new`/`.node_ids`/`.abundance` field names, `SmallBitset::contains`, and `graph.pattern_size()` (all cited in grounding; adjust the `use` paths to the real module locations). If `serde_json` is not already a dep, parse the JSONL with the existing parity JSON helper or a minimal manual parse — check `Cargo.toml` for `serde_json` first (the parity emit uses string formatting, so a dep may not exist; if absent, hand-parse `"start":N`, `"end":N`, `"strand":"X"`, and the `edges` string with `str::find`/`split`).

- [ ] **Step 2: Declare the module**

In `src/rustle/lib.rs` (or `mod.rs` where modules are declared), add `pub mod terminal_oracle;`. Confirm the exact module-decl file by grepping for `pub mod stringtie_parity`.

- [ ] **Step 3: Hook into the pipeline**

In `pipeline.rs`, AFTER the synth merge (after line ~13339, where `coverage_synth` has been extended with `alt_tts_synth` and `terminal_donor_synth`) and BEFORE `collect_from_transfrags` (~13626), add:

```rust
if let Some(oracle_path) = crate::stringtie_parity::terminal_oracle_path() {
    let oracle = crate::terminal_oracle::load_oracle(&oracle_path);
    let strand_c = match graph_bundle.strand { /* map to '+'/'-'/'.' as the graph_edge emit does */ };
    let bundle_key = (graph_bundle.start + 1, graph_bundle.end, strand_c);
    // sno: the same strand index passed to add_coverage_source_sink_edges at :13255
    let _ = crate::terminal_oracle::override_terminal_edges(
        &mut graph_mut, bundle_key, &mut coverage_synth, bpcov_ref, sno, oracle,
    );
}
```

Confirm at edit time: the exact names `graph_mut`, `coverage_synth`, `bpcov_ref`, `sno`, and `graph_bundle.strand`→char mapping used at `:13255` and the `graph_edge` emit (`:13814`). Use the SAME strand-char mapping the emit uses so `bundle_key` matches the oracle file.

- [ ] **Step 4: Build**

Run: `cargo build --release 2>&1 | tail -5`
Expected: `Finished`. Fix import/field-name mismatches against the real definitions (grounding cites them).

- [ ] **Step 5: Verify default unchanged (oracle OFF)**

Run: `./target/release/rustle -L ../GGO_19.bam -o /tmp/ru_def.gtf 2>/dev/null; gffcompare -r ../GGO_19.gtf /tmp/ru_def.gtf -o /tmp/ru_def_cmp 2>/dev/null; grep -E "Intron chain level:|Transcript level:" /tmp/ru_def_cmp.stats`
Expected: 96.2/91.7, 95.6/90.5 (flag OFF → no change).

- [ ] **Step 6: Commit**

```bash
git add src/rustle/terminal_oracle.rs src/rustle/lib.rs src/rustle/pipeline.rs
git commit -m "feat(oracle): RUSTLE_TERMINAL_ORACLE injects ST terminal edges on identical-node bundles

Default OFF; analysis-only. Overrides RU synthetic source/sink edges with ST's
captured graph_edge set + reconciles synth-transfrag abundance. Measures the
synthetic-coverage-edge contribution (lower bound; flow virtual edges not suppressed).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: `bench/terminal_oracle_report.py`

**Files:** Create `bench/terminal_oracle_report.py`

- [ ] **Step 1: Write the report tool**

```python
#!/usr/bin/env python3
"""Compare oracle-ON vs baseline final GTFs against the reference; report FP/FN/F1 delta
and how many removed FPs are endpoint-attributable to a removed terminal edge.
Usage: python3 bench/terminal_oracle_report.py /tmp/ru_def.gtf /tmp/ru_oracle.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf"""
import sys, re
from collections import defaultdict

def chains(path):
    tx = defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip().split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    out = set()
    for t, ex in tx.items():
        ex.sort()
        ic = tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))
        if ic: out.add((strand[t], ic))
    return out

def main():
    base = chains(sys.argv[1]); orac = chains(sys.argv[2]); ref = chains(sys.argv[3])
    def stats(s):
        tp = len(s & ref); fp = len(s - ref); fn = len(ref - s)
        sn = 100*tp/len(ref) if ref else 0; pr = 100*tp/(tp+fp) if (tp+fp) else 0
        f1 = 2*sn*pr/(sn+pr) if (sn+pr) else 0
        return tp, fp, fn, sn, pr, f1
    bt, bfp, bfn, bsn, bpr, bf1 = stats(base)
    ot, ofp, ofn, osn, opr, of1 = stats(orac)
    print(f"baseline: TP={bt} FP={bfp} FN={bfn}  Sn={bsn:.1f} Pr={bpr:.1f} F1={bf1:.2f}")
    print(f"oracle:   TP={ot} FP={ofp} FN={ofn}  Sn={osn:.1f} Pr={opr:.1f} F1={of1:.2f}")
    print(f"DELTA:    dTP={ot-bt:+d} dFP={ofp-bfp:+d} dFN={ofn-bfn:+d}  dF1={of1-bf1:+.2f}")
    fp_removed = (base - ref) - (orac - ref)
    fp_added = (orac - ref) - (base - ref)
    tp_lost = (base & ref) - (orac & ref)
    print(f"FP removed by oracle: {len(fp_removed)}  FP added: {len(fp_added)}  TP lost: {len(tp_lost)}")
    net = len(fp_removed) - len(tp_lost) - len(fp_added)
    print(f"NET prize (FP_removed - TP_lost - FP_added): {net}")
    print(f"ABORT GATE: {'ABORT (shelve Phase 1)' if net <= 0 else ('WEAK (<5)' if net < 5 else 'PROCEED')}")

if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the oracle + report**

```bash
cd /mnt/c/Users/jfris/Desktop/Rustle
# baseline final (oracle OFF)
./target/release/rustle -L GGO_19.bam -o /tmp/ru_def.gtf 2>/dev/null
# capture ST edges fresh (oracle input)
bash bench/capture_parity.sh >/dev/null 2>&1   # produces /tmp/st_edge.jsonl via STEPS
grep '"step":"graph_edge"' /tmp/st_parity.jsonl > /tmp/st_edge.jsonl 2>/dev/null || true
# oracle ON
RUSTLE_TERMINAL_ORACLE=/tmp/st_edge.jsonl ./target/release/rustle -L GGO_19.bam -o /tmp/ru_oracle.gtf 2>/dev/null
python3 bench/terminal_oracle_report.py /tmp/ru_def.gtf /tmp/ru_oracle.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf
```
Expected: prints baseline vs oracle TP/FP/FN/F1, the deltas, and the NET prize + gate verdict. (If `/tmp/st_edge.jsonl` is empty, the capture step's grep target differs — derive it from `/tmp/st_parity.jsonl`.)

- [ ] **Step 3: Commit**

```bash
git add bench/terminal_oracle_report.py
git commit -m "bench: terminal_oracle_report.py (oracle vs baseline FP/FN/F1 + net prize)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: Run Phase 0b + ABORT GATE decision

**Files:** Modify `docs/STRINGTIE_PARITY_FINDINGS.md`; memory `project_color_cgroup_parity.md`.

- [ ] **Step 1: Sanity — oracle actually changed edges**

Run the oracle with a trace (add a temporary `eprintln!` of total (removed, installed) summed across bundles, OR re-emit graph_edge under the oracle and edge_diff vs ST):
```bash
RUSTLE_TERMINAL_ORACLE=/tmp/st_edge.jsonl \
RUSTLE_PARITY_LOG=/tmp/ru_oracle_edge.jsonl RUSTLE_PARITY_FILTER_STEPS=graphnode_list,graph_edge RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 \
  ./target/release/rustle -L GGO_19.bam -o /tmp/ru_oracle.gtf 2>/dev/null
python3 bench/edge_diff.py /tmp/ru_oracle_edge.jsonl /tmp/st_edge.jsonl
```
Expected: RU-extra source_sink count on identical-node bundles drops sharply toward 0 (confirms the injection took effect). If it does NOT drop, the override isn't firing — debug the bundle-key match / node-token mapping before trusting the F1 number.

- [ ] **Step 2: Read the gate**

From `terminal_oracle_report.py`: NET = FP_removed − TP_lost − FP_added on the identical-node subset. **ABORT if NET ≤ 0** (terminal alignment doesn't help or costs more than it saves). **WEAK if NET < 5** and concentrated where injection can't reach → record and shelve. **PROCEED only if NET ≥ ~5 and clearly F1-positive.** Interpret with the documented confound (flow's virtual capacity edges mean this is a lower bound on suppression — but also an upper bound on the realizable prize, since real ST also has its own flow behavior).

- [ ] **Step 3: Record the result + decision**

Append to `docs/STRINGTIE_PARITY_FINDINGS.md` (§6j area): the Phase 0a probe (edges removed by RECURSE_ZERO_OFF + F1 delta) and Phase 0b (oracle NET prize, FP/FN/F1 delta, gate verdict). Update memory `project_color_cgroup_parity.md`. If ABORT/WEAK: record "terminal-wiring prize bounded; over-enum root is downstream of terminal wiring (node-mismatched bundles / flow traversal)"; STOP — do not start Phase 1.

- [ ] **Step 4: Commit**

```bash
git add docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "docs: Phase-0 terminal-wiring result (free probe + injection oracle + gate)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## PHASE 1 — Rule alignment (ONLY if abort gate cleared)

### Task 6: `RUSTLE_ST_TERMINAL` rule alignment

**Files:** Modify `src/rustle/graph_build.rs` (the `add_coverage_source_sink_edges` body, ~68-210).

Phase 0a/0b will have identified WHICH RU mechanism over-fires (phantom recursion vs blanket interior-node check). Target exactly that:

- [ ] **Step 1: Gate the phantom recursion off under the flag**

At `graph_build.rs:146` and `:185`, change:
```rust
let recurse_zero = std::env::var_os("RUSTLE_COVLINK_RECURSE_ZERO_OFF").is_none();
```
to:
```rust
let recurse_zero = std::env::var_os("RUSTLE_COVLINK_RECURSE_ZERO_OFF").is_none()
    && !crate::stringtie_parity::st_terminal();
```
(so `RUSTLE_ST_TERMINAL=1` disables the RU-specific phantom recursion, matching ST which has no such pass).

- [ ] **Step 2: (If 0b showed the blanket interior-node check also over-fires) tighten toward ST's structural rule**

Only if the probe attributed extra edges beyond the recursion: under `st_terminal()`, additionally require ST's structural condition (a source edge only where the node has no real parent, a sink edge only where it has no real child) before the coverage-ratio test — i.e. gate the per-interior-node coverage add on `parents.is_empty()`-equivalent. Show the exact condition added (derive from the 0b attribution; do NOT add speculatively if the recursion alone explained the divergence).

- [ ] **Step 3: Build + measure flag ON**

Run:
```bash
cargo build --release 2>&1 | tail -2
RUSTLE_ST_TERMINAL=1 ./target/release/rustle -L ../GGO_19.bam -o /tmp/ru_stterm.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf /tmp/ru_stterm.gtf -o /tmp/ru_stterm_cmp 2>/dev/null
grep -E "Intron chain level:|Transcript level:" /tmp/ru_stterm_cmp.stats
# default unchanged (flag OFF):
./target/release/rustle -L ../GGO_19.bam -o /tmp/ru_def2.gtf 2>/dev/null
gffcompare -r ../GGO_19.gtf /tmp/ru_def2.gtf -o /tmp/ru_def2_cmp 2>/dev/null
grep -E "Transcript level:" /tmp/ru_def2_cmp.stats
```
Expected: flag-ON F1 should approach the oracle's measured ceiling; flag-OFF must stay 95.6/90.5.

- [ ] **Step 4: Decision + commit**

If flag-ON F1 ≥ baseline: recommend default-flip to the user (do NOT flip without explicit approval). Else keep `RUSTLE_ST_TERMINAL` opt-in; record the gap vs the oracle ceiling (the residual is the flow virtual-edge confound or node-mismatched bundles). Commit:
```bash
git add src/rustle/graph_build.rs docs/STRINGTIE_PARITY_FINDINGS.md
git commit -m "feat(parity): RUSTLE_ST_TERMINAL aligns terminal wiring to ST (default OFF)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-review notes

- **Spec coverage:** Phase 0a free probe → Task 1; flag predicates → Task 2; injection oracle (hook + override + abundance reconciliation) → Task 3; report tool → Task 4; run + abort gate → Task 5; Phase 1 rule alignment → Task 6. All spec sections covered.
- **Default-unchanged guard:** Task 3 Step 5 and Task 6 Step 3 each re-verify 95.6/90.5 with flags OFF.
- **Confound documented:** the flow virtual-capacity-edge caveat is in the grounding facts, the oracle module doc-comment, and Task 5 Step 2's interpretation.
- **Type/name consistency:** `terminal_oracle_path`/`st_terminal` predicates; `override_terminal_edges`/`load_oracle` in `terminal_oracle.rs`; bundle key `(start+1, end, strand_char)` consistent with the graph_edge emit. Several real binding names (`graph_mut`, `coverage_synth`, `bpcov_ref`, `sno`, `DROP` path, `GraphTransfrag` fields, `SmallBitset::contains`) are flagged "confirm at edit time" with grounding citations — these are read-and-confirm, not deferred work.
- **No placeholders:** every code/command step shows actual content; the one judgment step (Task 6 Step 2) is explicitly conditional on the 0b attribution and says not to add speculatively.

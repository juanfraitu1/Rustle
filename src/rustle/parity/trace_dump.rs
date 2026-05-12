//! Optional **multi-kind** trace TSV for StringTie-style parity debugging.
//!
//! Enable with `RUSTLE_PARITY_TRACE_TSV=/path/trace.tsv`. Rows use a small fixed header and a
//! `payload` column (often `key=value;` pairs) so you can grep by `kind`. **`chrom` / `start` /
//! `end` match the parent pipeline bundle** (same key as `RUSTLE_PARITY_PARTITION_TSV`), even
//! when a row describes a per-subbundle graph (`run_tag` names the slice).
//!
//! | `kind` | When | Payload / columns |
//! |--------|------|---------------------|
//! | `cgroup_summary_v1` | After cgroup → bundlenodes for the pipeline bundle | `sub_bundles`, `bnodes`, `color_roots`, `assigned_reads` |
//! | `graph_topology_v1` | Per graph slice; `stage` column disambiguates | `post_graph_build`: after longtrim + prune + optional zero-collapse. `pre_read_map`: after coverage/terminal synth and pruning, **before** read→transfrag mapping. Payload: `run_tag`, `nodes`, `edges`, `edge_fp` |
//! | `read_map_summary_v1` | After `map_reads_to_graph*` | `run_tag`, `read_mapped_tfs`, `nodes`, `edges`, `edge_fp` (graph state used for mapping) |
//! | `partition_chain_detail_v1` | Optional (`RUSTLE_PARITY_CHAIN_DETAIL=1`) chain dump for partition analysis | `mode={subbundle|color_components}`, `chain_idx`, `segments` |
//! | `read_exons_v1` | Per read, only if `RUSTLE_PARITY_READ_EXONS=1` (or `true`) | `read_idx`, `exon_chain` (`start-end+…`, right end like partition emit) |
//!
//! StringTie does not yet emit matching `graph_topology_v1` / `cgroup_summary_v1` rows; add
//! parallel hooks in `rlink.cpp` when you need byte-for-byte graph parity. Until then, Rustle
//! rows anchor where divergence likely begins (cgroup counts vs bundlenode chain, then node/edge
//! counts vs topology fingerprint).

use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::{Mutex, OnceLock};

use anyhow::Result;

use crate::graph::Graph;

static WRITER: OnceLock<Mutex<Option<BufWriter<File>>>> = OnceLock::new();

fn writer_cell() -> &'static Mutex<Option<BufWriter<File>>> {
    WRITER.get_or_init(|| Mutex::new(None))
}

pub fn enabled() -> bool {
    std::env::var_os("RUSTLE_PARITY_TRACE_TSV").is_some()
}

fn trace_level() -> u32 {
    std::env::var("RUSTLE_PARITY_TRACE_LEVEL")
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(0)
}

pub fn event_enabled(level: u32) -> bool {
    enabled() && trace_level() >= level
}

pub fn read_exons_enabled() -> bool {
    match std::env::var("RUSTLE_PARITY_READ_EXONS") {
        Ok(s) => s == "1" || s.eq_ignore_ascii_case("true"),
        Err(_) => false,
    }
}

pub fn chain_detail_enabled() -> bool {
    match std::env::var("RUSTLE_PARITY_CHAIN_DETAIL") {
        Ok(s) => s == "1" || s.eq_ignore_ascii_case("true"),
        Err(_) => false,
    }
}

/// Exon intervals as `start-end` joined by `+`, with `end` emitted like partition signatures
/// (exclusive stored end → inclusive display via `end - 1`).
pub fn format_read_exon_chain(exons: &[(u64, u64)]) -> String {
    exons
        .iter()
        .map(|(s, e)| {
            let e_out = e.saturating_sub(1);
            format!("{}-{}", s, e_out)
        })
        .collect::<Vec<_>>()
        .join("+")
}

pub fn init() -> Result<()> {
    let path = match std::env::var_os("RUSTLE_PARITY_TRACE_TSV") {
        Some(p) => p,
        None => {
            *writer_cell().lock().expect("parity trace lock poisoned") = None;
            return Ok(());
        }
    };
    let path = Path::new(&path);
    let file = File::create(path)?;
    let mut w = BufWriter::new(file);
    writeln!(
        w,
        "kind\tsource\tstage\tchrom\tstart\tend\tstrand\tpayload"
    )?;
    *writer_cell().lock().expect("parity trace lock poisoned") = Some(w);
    Ok(())
}

fn graph_edge_fingerprint(g: &Graph) -> u64 {
    let mut pairs: Vec<(usize, usize)> = g
        .gpos
        .keys()
        .map(|&(a, b)| if a <= b { (a, b) } else { (b, a) })
        .collect();
    pairs.sort_unstable();
    let mut h = DefaultHasher::new();
    pairs.hash(&mut h);
    h.finish()
}

pub fn emit_cgroup_row(
    chrom: &str,
    start: u64,
    end: u64,
    strand: char,
    sub_bundles: usize,
    bnodes: usize,
    color_roots: usize,
    assigned_reads: usize,
) {
    if !enabled() {
        return;
    }
    let payload = format!(
        "sub_bundles={};bnodes={};color_roots={};assigned_reads={}",
        sub_bundles, bnodes, color_roots, assigned_reads
    );
    let mut guard = writer_cell().lock().expect("parity trace lock poisoned");
    let Some(w) = guard.as_mut() else {
        return;
    };
    let _ = writeln!(
        w,
        "cgroup_summary_v1\trustle\tpost_cgroup\t{}\t{}\t{}\t{}\t{}",
        chrom, start, end, strand, payload
    );
}

pub fn emit_graph_row(
    chrom: &str,
    start: u64,
    end: u64,
    strand: char,
    stage: &str,
    run_tag: &str,
    graph: &Graph,
) {
    if !enabled() {
        return;
    }
    let n_nodes = graph.nodes.len();
    let n_edges = graph.gpos.len();
    let fp = graph_edge_fingerprint(graph);
    let payload = format!(
        "run_tag={};nodes={};edges={};edge_fp={:016x}",
        run_tag, n_nodes, n_edges, fp
    );
    let mut guard = writer_cell().lock().expect("parity trace lock poisoned");
    let Some(w) = guard.as_mut() else {
        return;
    };
    let _ = writeln!(
        w,
        "graph_topology_v1\trustle\t{}\t{}\t{}\t{}\t{}\t{}",
        stage, chrom, start, end, strand, payload
    );
}

/// After read→transfrag mapping: how many read-backed transfrags and the graph fingerprint
/// at that point (per `run_tag` / graph slice).
pub fn emit_read_map_summary_row(
    chrom: &str,
    start: u64,
    end: u64,
    strand: char,
    run_tag: &str,
    read_mapped_transfrags: usize,
    graph: &Graph,
) {
    if !enabled() {
        return;
    }
    let n_nodes = graph.nodes.len();
    let n_edges = graph.gpos.len();
    let fp = graph_edge_fingerprint(graph);
    let payload = format!(
        "run_tag={};read_mapped_tfs={};nodes={};edges={};edge_fp={:016x}",
        run_tag, read_mapped_transfrags, n_nodes, n_edges, fp
    );
    let mut guard = writer_cell().lock().expect("parity trace lock poisoned");
    let Some(w) = guard.as_mut() else {
        return;
    };
    let _ = writeln!(
        w,
        "read_map_summary_v1\trustle\tpost_read_map\t{}\t{}\t{}\t{}\t{}",
        chrom, start, end, strand, payload
    );
}

pub fn emit_read_exon_row(chrom: &str, start: u64, end: u64, strand: char, read_idx: usize, exon_chain: &str) {
    if !enabled() || !read_exons_enabled() {
        return;
    }
    let payload = format!("read_idx={};exon_chain={}", read_idx, exon_chain);
    let mut guard = writer_cell().lock().expect("parity trace lock poisoned");
    let Some(w) = guard.as_mut() else {
        return;
    };
    let _ = writeln!(
        w,
        "read_exons_v1\trustle\tpost_cgroup\t{}\t{}\t{}\t{}\t{}",
        chrom, start, end, strand, payload
    );
}

pub fn emit_partition_chain_row(
    chrom: &str,
    start: u64,
    end: u64,
    strand: char,
    mode: &str,
    chain_idx: usize,
    segments: &str,
) {
    if !enabled() || !chain_detail_enabled() {
        return;
    }
    let payload = format!("mode={};chain_idx={};segments={}", mode, chain_idx, segments);
    let mut guard = writer_cell().lock().expect("parity trace lock poisoned");
    let Some(w) = guard.as_mut() else {
        return;
    };
    let _ = writeln!(
        w,
        "partition_chain_detail_v1\trustle\tpost_cgroup\t{}\t{}\t{}\t{}\t{}",
        chrom, start, end, strand, payload
    );
}

/// Generic high-granularity event row. Use level 1 for major function enter/exit,
/// level 2 for per-stage counters, level 3 for loop-level high-volume details.
pub fn emit_event_row(
    chrom: &str,
    start: u64,
    end: u64,
    strand: char,
    stage: &str,
    level: u32,
    payload: &str,
) {
    if !event_enabled(level) {
        return;
    }
    let mut guard = writer_cell().lock().expect("parity trace lock poisoned");
    let Some(w) = guard.as_mut() else {
        return;
    };
    let _ = writeln!(
        w,
        "trace_event_v1\trustle\t{}\t{}\t{}\t{}\t{}\tlevel={};{}",
        stage, chrom, start, end, strand, level, payload
    );
}

pub fn flush() {
    let mut guard = writer_cell().lock().expect("parity trace lock poisoned");
    if let Some(w) = guard.as_mut() {
        let _ = w.flush();
    }
}

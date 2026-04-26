//! Splice-graph **edges** dump for cross-tool parity / Python visualizers.
//!
//! Enable with `RUSTLE_PARITY_GRAPH_EDGES_TSV=/path/edges.tsv`. Complements
//! `RUSTLE_PARITY_GRAPH_TSV` (which dumps nodes only). Each row describes a
//! single directed edge between two non-source/non-sink nodes, classified as
//! `collinear` (adjacent base coordinates) or `junction` (anything else).
//!
//! Schema (header line):
//! `source	chrom	bdstart	bdend	strand	from_idx	to_idx	from_start	from_end	to_start	to_end	edge_kind	from_cov	to_cov	bottleneck_cov`
//!
//! Mirrors the inline emitter pattern used by `RUSTLE_PARITY_GRAPH_TSV` /
//! `RUSTLE_PARITY_TF_TSV` in `path_extract.rs` (verbose by design — do not DRY).

use std::io::Write;
use std::sync::{Mutex, OnceLock};

use crate::graph::Graph;

static WR: OnceLock<Mutex<Option<std::fs::File>>> = OnceLock::new();
static HDR: OnceLock<Mutex<bool>> = OnceLock::new();

/// Emit one row per directed edge between non-source/non-sink graph nodes.
///
/// Silently skips when the env var is unset or the file can't be opened.
pub fn emit(graph: &Graph, bundle_chrom: &str, bundle_strand: char, bundle_id: &str) {
    let gp = match std::env::var("RUSTLE_PARITY_GRAPH_EDGES_TSV") {
        Ok(p) if !p.is_empty() => p,
        _ => return,
    };
    let wr = WR.get_or_init(|| {
        Mutex::new(
            std::fs::OpenOptions::new()
                .create(true).append(true).open(&gp).ok(),
        )
    });
    let hdr = HDR.get_or_init(|| Mutex::new(false));
    let (bs, be) = bundle_id.split_once(':')
        .and_then(|(_, rng)| rng.split_once('-'))
        .map(|(s, e)| (
            s.trim().parse::<u64>().unwrap_or(0),
            e.trim().parse::<u64>().unwrap_or(0),
        ))
        .unwrap_or((0, 0));
    if let (Ok(mut f_opt), Ok(mut hdr_w)) = (wr.lock(), hdr.lock()) {
        if let Some(f) = f_opt.as_mut() {
            if !*hdr_w {
                let _ = writeln!(f, "source\tchrom\tbdstart\tbdend\tstrand\tfrom_idx\tto_idx\tfrom_start\tfrom_end\tto_start\tto_end\tedge_kind\tfrom_cov\tto_cov\tbottleneck_cov");
                *hdr_w = true;
            }
            for (idx, n) in graph.nodes.iter().enumerate() {
                if idx == graph.source_id || idx == graph.sink_id {
                    continue;
                }
                for to_idx in n.children.ones() {
                    if to_idx == graph.source_id || to_idx == graph.sink_id {
                        continue;
                    }
                    let to = match graph.nodes.get(to_idx) {
                        Some(t) => t,
                        None => continue,
                    };
                    let edge_kind = if to.start == n.end + 1 { "collinear" } else { "junction" };
                    let bottleneck = n.coverage.min(to.coverage);
                    let _ = writeln!(
                        f,
                        "rustle\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}",
                        bundle_chrom, bs, be, bundle_strand,
                        idx, to_idx,
                        n.start, n.end, to.start, to.end,
                        edge_kind,
                        n.coverage, to.coverage, bottleneck,
                    );
                }
            }
        }
    }
}

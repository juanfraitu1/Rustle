//! Per-bundle stage TSV (`--parity-stage-tsv` or env `RUSTLE_PARITY_STAGE_TSV`).
//!
//! In needy-locus runs (`--trace-reference`, `--debug-bundle` / `--only-debug-bundle`, or
//! `RUSTLE_TRACE_LOCUS`), rustle also writes `{output}.needy_bundle_summary.tsv` unless
//! `RUSTLE_DISABLE_AUTO_BUNDLE_SUMMARY` is set. Same schema as `--parity-stage-tsv`.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Mutex, OnceLock};

use anyhow::Result;

use crate::graph::Graph;
use crate::graph::GraphTransfrag;
use crate::path_extract::Transcript;
use crate::types::{Bundle, CBundlenode};

static STAGE_WRITER: OnceLock<Mutex<Option<BufWriter<File>>>> = OnceLock::new();
static STAGE_IO_ENABLED: AtomicBool = AtomicBool::new(false);

/// Fast path: `false` when no `--parity-stage-tsv` (avoids a mutex lock per `emit*` call).
#[inline]
pub fn is_enabled() -> bool {
    STAGE_IO_ENABLED.load(Ordering::Relaxed)
}

#[derive(Debug, Clone, Copy, Default)]
pub struct StageDetail {
    pub bundle_bnodes: usize,
    pub bundle_read_bnode_links: usize,
    pub bundle_color_count: usize,
    pub struct_hash: u64,
    pub graph_active_transfrags: usize,
    pub graph_long_transfrags: usize,
    pub graph_pattern_size: usize,
    pub graph_hard_boundary_nodes: usize,
    pub seed_unwitnessed: usize,
    pub seed_zero_flux: usize,
    pub seed_low_cov: usize,
    pub seed_back_fail: usize,
    pub seed_fwd_fail: usize,
    pub longrec_back_no_choice: usize,
    pub longrec_fwd_no_choice: usize,
    pub longrec_back_no_reach: usize,
    pub longrec_fwd_no_reach: usize,
    pub longrec_back_unreachable: usize,
    pub longrec_fwd_unreachable: usize,
    pub checktrf_total: usize,
    pub checktrf_pre_filter: usize,
    pub predcluster_before: usize,
    pub predcluster_removed: usize,
    pub pairwise_included_kill: usize,
    pub pairwise_secontained_kill: usize,
    pub pairwise_intronic_kill: usize,
    pub pairwise_bettercov_kill: usize,
    pub pairwise_other_kill: usize,
    pub isofrac_longunder_kill: usize,
    pub junction_support_removed: usize,
}

fn writer_cell() -> &'static Mutex<Option<BufWriter<File>>> {
    STAGE_WRITER.get_or_init(|| Mutex::new(None))
}

pub fn init<P: AsRef<Path>>(path: Option<P>) -> Result<()> {
    let mut guard = writer_cell().lock().expect("stage writer lock poisoned");
    if let Some(path) = path {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(
            writer,
            "stage_id\tchrom\tstart\tend\tstrand\tn_reads\tn_junctions\tn_nodes\tn_edges\tn_transfrags\tn_transcripts\tseed_total\tseed_stored\tseed_rescued\tbundle_bnodes\tbundle_read_bnode_links\tbundle_color_count\tstruct_hash\tgraph_active_transfrags\tgraph_long_transfrags\tgraph_pattern_size\tgraph_hard_boundary_nodes\tseed_unwitnessed\tseed_zero_flux\tseed_low_cov\tseed_back_fail\tseed_fwd_fail\tlongrec_back_no_choice\tlongrec_fwd_no_choice\tlongrec_back_no_reach\tlongrec_fwd_no_reach\tlongrec_back_unreachable\tlongrec_fwd_unreachable\tchecktrf_total\tchecktrf_pre_filter\tpredcluster_before\tpredcluster_removed\tpairwise_included_kill\tpairwise_secontained_kill\tpairwise_intronic_kill\tpairwise_bettercov_kill\tpairwise_other_kill\tisofrac_longunder_kill\tjunction_support_removed\tnote"
        )?;
        *guard = Some(writer);
        STAGE_IO_ENABLED.store(true, Ordering::Relaxed);
    } else {
        *guard = None;
        STAGE_IO_ENABLED.store(false, Ordering::Relaxed);
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn emit(
    stage_id: &str,
    chrom: &str,
    start: u64,
    end: u64,
    strand: char,
    n_reads: usize,
    n_junctions: usize,
    n_nodes: usize,
    n_edges: usize,
    n_transfrags: usize,
    n_transcripts: usize,
    seed_total: usize,
    seed_stored: usize,
    seed_rescued: usize,
    note: &str,
) {
    emit_with_detail(
        stage_id,
        chrom,
        start,
        end,
        strand,
        n_reads,
        n_junctions,
        n_nodes,
        n_edges,
        n_transfrags,
        n_transcripts,
        seed_total,
        seed_stored,
        seed_rescued,
        StageDetail::default(),
        note,
    );
}

#[allow(clippy::too_many_arguments)]
pub fn emit_with_detail(
    stage_id: &str,
    chrom: &str,
    start: u64,
    end: u64,
    strand: char,
    n_reads: usize,
    n_junctions: usize,
    n_nodes: usize,
    n_edges: usize,
    n_transfrags: usize,
    n_transcripts: usize,
    seed_total: usize,
    seed_stored: usize,
    seed_rescued: usize,
    detail: StageDetail,
    note: &str,
) {
    if !STAGE_IO_ENABLED.load(Ordering::Relaxed) {
        return;
    }
    let mut guard = writer_cell().lock().expect("stage writer lock poisoned");
    let Some(writer) = guard.as_mut() else {
        return;
    };
    let _ = writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        stage_id,
        chrom,
        start,
        end,
        strand,
        n_reads,
        n_junctions,
        n_nodes,
        n_edges,
        n_transfrags,
        n_transcripts,
        seed_total,
        seed_stored,
        seed_rescued,
        detail.bundle_bnodes,
        detail.bundle_read_bnode_links,
        detail.bundle_color_count,
        detail.struct_hash,
        detail.graph_active_transfrags,
        detail.graph_long_transfrags,
        detail.graph_pattern_size,
        detail.graph_hard_boundary_nodes,
        detail.seed_unwitnessed,
        detail.seed_zero_flux,
        detail.seed_low_cov,
        detail.seed_back_fail,
        detail.seed_fwd_fail,
        detail.longrec_back_no_choice,
        detail.longrec_fwd_no_choice,
        detail.longrec_back_no_reach,
        detail.longrec_fwd_no_reach,
        detail.longrec_back_unreachable,
        detail.longrec_fwd_unreachable,
        detail.checktrf_total,
        detail.checktrf_pre_filter,
        detail.predcluster_before,
        detail.predcluster_removed,
        detail.pairwise_included_kill,
        detail.pairwise_secontained_kill,
        detail.pairwise_intronic_kill,
        detail.pairwise_bettercov_kill,
        detail.pairwise_other_kill,
        detail.isofrac_longunder_kill,
        detail.junction_support_removed,
        note
    );
}

pub fn flush() {
    let mut guard = writer_cell().lock().expect("stage writer lock poisoned");
    if let Some(writer) = guard.as_mut() {
        let _ = writer.flush();
    }
}

pub fn graph_edge_count(graph: &Graph) -> usize {
    graph.edgeno
}

const HASH_OFFSET: u64 = 1469598103934665603;
const HASH_PRIME: u64 = 1099511628211;

fn hash_mix(hash: &mut u64, value: u64) {
    *hash ^= value;
    *hash = hash.wrapping_mul(HASH_PRIME);
}

fn count_bundlenodes(mut head: Option<&CBundlenode>) -> usize {
    let mut count = 0usize;
    while let Some(node) = head {
        count += 1;
        head = node.next.as_deref();
    }
    count
}

pub fn bundle_detail(bundle: &Bundle) -> StageDetail {
    let mut hash = HASH_OFFSET;
    hash_mix(&mut hash, bundle.start);
    hash_mix(&mut hash, bundle.end);
    hash_mix(&mut hash, bundle.strand as u64);
    hash_mix(&mut hash, bundle.reads.len() as u64);
    for read in &bundle.reads {
        hash_mix(&mut hash, read.ref_start);
        hash_mix(&mut hash, read.ref_end);
        hash_mix(&mut hash, read.exons.len() as u64);
    }
    let bundle_bnodes = count_bundlenodes(bundle.bundlenodes.as_ref());
    let bundle_read_bnode_links = bundle
        .read_bnodes
        .as_ref()
        .map(|rows| rows.iter().map(|row| row.len()).sum())
        .unwrap_or(0usize);
    let bundle_color_count = bundle
        .bnode_colors
        .as_ref()
        .map(|colors| colors.len())
        .unwrap_or(0usize);
    hash_mix(&mut hash, bundle_bnodes as u64);
    hash_mix(&mut hash, bundle_read_bnode_links as u64);
    hash_mix(&mut hash, bundle_color_count as u64);
    StageDetail {
        bundle_bnodes,
        bundle_read_bnode_links,
        bundle_color_count,
        struct_hash: hash,
        ..StageDetail::default()
    }
}

pub fn graph_detail(graph: &Graph, transfrags: &[GraphTransfrag]) -> StageDetail {
    let mut hash = HASH_OFFSET;
    hash_mix(&mut hash, graph.n_nodes as u64);
    hash_mix(&mut hash, graph.edgeno as u64);
    let mut hard_boundary_nodes = 0usize;
    for (idx, node) in graph.nodes.iter().enumerate() {
        if idx == graph.source_id || idx == graph.sink_id {
            continue;
        }
        hash_mix(&mut hash, node.start);
        hash_mix(&mut hash, node.end);
        hash_mix(&mut hash, node.children.count_ones() as u64);
        hash_mix(&mut hash, node.parents.count_ones() as u64);
        hash_mix(&mut hash, node.hardstart as u64);
        hash_mix(&mut hash, node.hardend as u64);
        if node.hardstart || node.hardend {
            hard_boundary_nodes += 1;
        }
    }
    let mut graph_active_transfrags = 0usize;
    let mut graph_long_transfrags = 0usize;
    for tf in transfrags {
        hash_mix(&mut hash, tf.node_ids.len() as u64);
        if let Some(&first) = tf.node_ids.first() {
            hash_mix(&mut hash, first as u64);
        }
        if let Some(&last) = tf.node_ids.last() {
            hash_mix(&mut hash, last as u64);
        }
        hash_mix(&mut hash, tf.longread as u64);
        hash_mix(&mut hash, tf.guide as u64);
        hash_mix(&mut hash, tf.trflong_seed as u64);
        if tf.abundance > 0.0 {
            graph_active_transfrags += 1;
        }
        if tf.longread {
            graph_long_transfrags += 1;
        }
    }
    StageDetail {
        struct_hash: hash,
        graph_active_transfrags,
        graph_long_transfrags,
        graph_pattern_size: graph.pattern_size(),
        graph_hard_boundary_nodes: hard_boundary_nodes,
        ..StageDetail::default()
    }
}

pub fn transcript_span(transcripts: &[Transcript]) -> Option<(u64, u64)> {
    let start = transcripts
        .iter()
        .flat_map(|tx| tx.exons.iter().map(|(s, _)| *s))
        .min()?;
    let end = transcripts
        .iter()
        .flat_map(|tx| tx.exons.iter().map(|(_, e)| *e))
        .max()?;
    Some((start, end))
}

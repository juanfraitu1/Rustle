//! Per-bundle JSONL snapshots for fast, 1:1 comparisons between:
//! - rustle pipeline internals (bundle -> graph -> transfrags -> transcripts)
//! - the reference assembler trace-derived “backward” packages (normalized with scripts)
//!
//! The intent is deterministic, compact-ish structure, not human readability.

use crate::graph::Graph;
use crate::path_extract::Transcript;
use crate::types::{Bundle, Junction};
use crate::GraphTransfrag;
use anyhow::{Context, Result};
use serde::Serialize;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SnapshotDetail {
    /// Only counts and spans (default unless forced).
    Summary,
    /// Include graph nodes/edges, transfrag paths, transcript exons.
    Full,
}

#[derive(Debug)]
pub struct SnapshotWriter {
    out: Option<BufWriter<File>>,
    detail: SnapshotDetail,
}

impl SnapshotWriter {
    pub fn new(path: Option<&str>) -> Result<Self> {
        let out = match path {
            None => None,
            Some(p) => {
                let f = File::create(Path::new(p))
                    .with_context(|| format!("create snapshot jsonl at {}", p))?;
                Some(BufWriter::new(f))
            }
        };
        // Default to summary. Full detail can be forced with env var to avoid
        // dumping huge JSON by accident on whole-genome runs.
        let detail = if std::env::var_os("RUSTLE_SNAPSHOT_FULL").is_some() {
            SnapshotDetail::Full
        } else {
            SnapshotDetail::Summary
        };
        Ok(Self { out, detail })
    }

    pub fn enabled(&self) -> bool {
        self.out.is_some()
    }

    pub fn detail(&self) -> SnapshotDetail {
        self.detail
    }

    pub fn set_detail(&mut self, detail: SnapshotDetail) {
        self.detail = detail;
    }

    pub fn emit<T: Serialize>(&mut self, item: &T) -> Result<()> {
        let Some(out) = self.out.as_mut() else {
            return Ok(());
        };
        serde_json::to_writer(&mut *out, item).context("serialize snapshot json")?;
        out.write_all(b"\n").context("write snapshot newline")?;
        Ok(())
    }
}

#[derive(Debug, Serialize)]
pub struct SnapshotEnvelope<'a, T: Serialize> {
    pub kind: &'static str,
    pub stage: &'a str,
    pub bundle: BundleKey<'a>,
    #[serde(flatten)]
    pub payload: T,
}

#[derive(Debug, Serialize)]
pub struct BundleKey<'a> {
    pub idx: usize,
    pub chrom: &'a str,
    pub start: u64,
    pub end: u64,
    pub strand: char,
}

#[derive(Debug, Serialize)]
pub struct BundleSummary {
    pub reads_total: usize,
    pub reads_multiexon: usize,
    pub junctions_total: usize,
}

#[derive(Debug, Serialize)]
pub struct BundleSummaryWithJunctions {
    pub summary: BundleSummary,
    pub junctions: Vec<JunctionSnap>,
}

#[derive(Debug, Serialize)]
pub struct JunctionSnap {
    pub donor: u64,
    pub acceptor: u64,
    pub strand: char,
}

#[derive(Debug, Serialize)]
pub struct GraphSummary {
    pub nodes: usize,
    pub edges: usize,
    pub source_id: usize,
    pub sink_id: usize,
}

#[derive(Debug, Serialize)]
pub struct GraphWithNodes {
    pub summary: GraphSummary,
    pub nodes: Vec<GraphNodeSnap>,
}

#[derive(Debug, Serialize)]
pub struct GraphNodeSnap {
    pub id: usize,
    pub start: u64,
    pub end: u64,
    pub cov: f64,
    pub hardstart: bool,
    pub hardend: bool,
}

#[derive(Debug, Serialize)]
pub struct TransfragSummary {
    pub transfrags: usize,
    pub longread: usize,
}

#[derive(Debug, Serialize)]
pub struct TransfragsFull {
    pub summary: TransfragSummary,
    pub transfrags: Vec<TransfragSnap>,
}

#[derive(Debug, Serialize)]
pub struct TransfragSnap {
    pub node_ids: Vec<usize>,
    pub abundance: f64,
    pub read_count: f64,
    pub srabund: f64,
    pub longread: bool,
    pub shortread: bool,
    pub longstart: u64,
    pub longend: u64,
    pub weak: u8,
    pub real: bool,
    pub guide: bool,
    pub killed_junction_orphan: bool,
    pub coverage_weak: bool,
    pub trflong_seed: bool,
    pub usepath: i32,
}

#[derive(Debug, Serialize)]
pub struct TranscriptSummary {
    pub transcripts: usize,
}

#[derive(Debug, Serialize)]
pub struct TranscriptsFull {
    pub summary: TranscriptSummary,
    pub transcripts: Vec<TranscriptSnap>,
}

#[derive(Debug, Serialize)]
pub struct TranscriptSnap {
    pub exons: Vec<(u64, u64)>,
    pub coverage: f64,
    pub longcov: f64,
    pub is_longread: bool,
    pub source: Option<String>,
    pub ref_transcript_id: Option<String>,
}

pub fn bundle_key<'a>(idx: usize, b: &'a Bundle) -> BundleKey<'a> {
    BundleKey {
        idx,
        chrom: &b.chrom,
        start: b.start,
        end: b.end,
        strand: b.strand,
    }
}

pub fn summarize_bundle(bundle: &Bundle) -> BundleSummary {
    let reads_total = bundle.reads.len();
    let reads_multiexon = bundle.reads.iter().filter(|r| r.exons.len() > 1).count();
    let junctions_total: usize = bundle.reads.iter().map(|r| r.junctions.len()).sum();
    BundleSummary {
        reads_total,
        reads_multiexon,
        junctions_total,
    }
}

pub fn snap_junctions(junctions: &[Junction], strand: char) -> Vec<JunctionSnap> {
    junctions
        .iter()
        .map(|j| JunctionSnap {
            donor: j.donor,
            acceptor: j.acceptor,
            strand,
        })
        .collect()
}

pub fn summarize_graph(graph: &Graph) -> GraphSummary {
    GraphSummary {
        nodes: graph.n_nodes,
        edges: graph.edgeno,
        source_id: graph.source_id,
        sink_id: graph.sink_id,
    }
}

pub fn snap_graph_nodes(graph: &Graph) -> Vec<GraphNodeSnap> {
    graph
        .nodes
        .iter()
        .enumerate()
        .map(|(id, n)| GraphNodeSnap {
            id,
            start: n.start,
            end: n.end,
            cov: n.coverage,
            hardstart: n.hardstart,
            hardend: n.hardend,
        })
        .collect()
}

pub fn summarize_transfrags(tfs: &[GraphTransfrag]) -> TransfragSummary {
    let transfrags = tfs.len();
    let longread = tfs.iter().filter(|t| t.longread).count();
    TransfragSummary { transfrags, longread }
}

pub fn snap_transfrags(tfs: &[GraphTransfrag]) -> Vec<TransfragSnap> {
    tfs.iter()
        .map(|t| TransfragSnap {
            node_ids: t.node_ids.clone(),
            abundance: t.abundance,
            read_count: t.read_count,
            srabund: t.srabund,
            longread: t.longread,
            shortread: t.shortread,
            longstart: t.longstart,
            longend: t.longend,
            weak: t.weak,
            real: t.real,
            guide: t.guide,
            killed_junction_orphan: t.killed_junction_orphan,
            coverage_weak: t.coverage_weak,
            trflong_seed: t.trflong_seed,
            usepath: t.usepath,
        })
        .collect()
}

pub fn summarize_transcripts(txs: &[Transcript]) -> TranscriptSummary {
    TranscriptSummary {
        transcripts: txs.len(),
    }
}

pub fn snap_transcripts(txs: &[Transcript]) -> Vec<TranscriptSnap> {
    txs.iter()
        .map(|t| TranscriptSnap {
            exons: t.exons.clone(),
            coverage: t.coverage,
            longcov: t.longcov,
            is_longread: t.is_longread,
            source: t.source.clone(),
            ref_transcript_id: t.ref_transcript_id.clone(),
        })
        .collect()
}

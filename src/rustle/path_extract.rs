//! Extract transcript paths from transfrags.
//! Path extension via extend_path_left/right (back_to_source_fast_long / fwd_to_sink_fast_long),
//! LR witness validation, poly-tail trimming, hardstart/hardend enforcement, and Edmonds-Karp flow.

use crate::bitset::NodeSet;
use crate::bitvec::GBitVec;
use crate::coord::{
    contiguous_half_open, len_half_open, overlap_len_half_open, overlaps_half_open,
};
use crate::graph::{Graph, GraphTransfrag};
use crate::max_flow::{
    guide_push_flow, long_max_flow_seeded_with_used_pathpat, push_guide_max_flow,
    push_max_flow_seeded_full, GuideFlowPriorRef,
};
use crate::types::{CPrediction, DetHashMap as HashMap, DetHashSet as HashSet, RunConfig};

#[inline]
fn node_start_or_zero(graph: &Graph, node_id: usize) -> u64 {
    graph.nodes.get(node_id).map(|n| n.start).unwrap_or(0)
}

#[inline]
fn node_end_or_zero(graph: &Graph, node_id: usize) -> u64 {
    graph.nodes.get(node_id).map(|n| n.end).unwrap_or(0)
}

/// `get_trf_long` / checktrf entry (~11133): guide OR abundance ≥ readthr.
/// After max-flow, `abundance` can drop below `readthr` while `read_count` still reflects evidence.
/// Use that only for **partially depleted** seeds (abundance above noise but below threshold) so we
/// do not open checktrf on zero-flux junk or destabilize pairwise competition (see STRG.251 panel).
#[inline]
fn checktrf_passes_abundance_gate(tf: &GraphTransfrag, readthr: f64) -> bool {
    // Guides always pass regardless of abundance (recovery of reference models).
    if tf.guide {
        return true;
    }
    // Standard abundance gate: transfrag must have enough remaining flow.
    // The partial-depletion relaxation (read_count >= readthr when abundance < readthr) was
    // removed because it allows reads already consumed by prior transcripts to generate
    // additional overlapping transcripts, inflating j/c class false positives.
    tf.abundance + EPS >= readthr
}

#[inline]
fn coord_min_node(graph: &Graph, a: usize, b: usize) -> usize {
    let (sa, sb) = (node_start_or_zero(graph, a), node_start_or_zero(graph, b));
    if sa < sb {
        a
    } else if sb < sa {
        b
    } else {
        a.min(b)
    }
}

#[inline]
fn coord_max_node(graph: &Graph, a: usize, b: usize) -> usize {
    let (sa, sb) = (node_start_or_zero(graph, a), node_start_or_zero(graph, b));
    if sa > sb {
        a
    } else if sb > sa {
        b
    } else {
        a.max(b)
    }
}

/// Parse RUSTLE_TRACE_LOCUS env var (format: "start-end") and check if node i's coordinates
/// overlap the trace range. Used for diagnostic tracing of back_to_source / fwd_to_sink.
fn trace_locus_active(graph: &Graph, i: usize) -> bool {
    let Ok(val) = std::env::var("RUSTLE_TRACE_LOCUS") else {
        return false;
    };
    let Some((s_str, e_str)) = val.split_once('-') else {
        return false;
    };
    let Ok(s) = s_str.trim().parse::<u64>() else {
        return false;
    };
    let Ok(e) = e_str.trim().parse::<u64>() else {
        return false;
    };
    graph
        .nodes
        .get(i)
        .map_or(false, |n| n.end >= s && n.start <= e)
}

/// Parse RUSTLE_TRACE_SEED_IDX and check whether the current seed index should emit
/// recursive path-extension debug even when the node itself is outside the traced locus
/// (e.g. source/sink or synthetic nodes).
fn trace_seed_active(seed_idx: usize) -> bool {
    std::env::var("RUSTLE_TRACE_SEED_IDX")
        .ok()
        .and_then(|v| v.trim().parse::<usize>().ok())
        == Some(seed_idx)
}

/// Global flag for detailed path extraction diagnostics.
/// Set RUSTLE_PATH_EXTRACT_DEBUG=1 for verbose per-seed diagnostics.
fn path_extract_debug() -> bool {
    std::env::var_os("RUSTLE_PATH_EXTRACT_DEBUG").is_some()
}

#[inline]
fn loop_trace_active() -> bool {
    std::env::var_os("RUSTLE_LOOP_TRACE").is_some()
        || std::env::var_os("RUSTLE_TRACE_LOG_STYLE").is_some()
}

#[inline]
fn same_node_coords(graph: &Graph, a: usize, b: usize) -> bool {
    match (graph.nodes.get(a), graph.nodes.get(b)) {
        (Some(na), Some(nb)) => na.start == nb.start && na.end == nb.end,
        _ => false,
    }
}

/// Check if detailed diagnostics should be emitted for this seed.
/// Either specific seed tracing is enabled, or global debug is on.
fn trace_seed_diagnostics(seed_idx: usize) -> bool {
    path_extract_debug() || trace_seed_active(seed_idx)
}

fn parse_trace_tf_ids() -> Vec<usize> {
    std::env::var("RUSTLE_TRACE_TF_IDS")
        .ok()
        .map(|v| {
            v.split(',')
                .filter_map(|s| s.trim().parse::<usize>().ok())
                .collect::<Vec<_>>()
        })
        .unwrap_or_default()
}

fn trace_tf_watch(
    seed_idx: usize,
    label: &str,
    watched: &[usize],
    transfrags: &[GraphTransfrag],
    graph: &Graph,
) {
    if watched.is_empty() || !trace_seed_diagnostics(seed_idx) {
        return;
    }
    eprintln!("[TRACE_TF_WATCH seed={}] {}", seed_idx, label);
    for &tid in watched {
        match transfrags.get(tid) {
            Some(tf) => {
                let first = tf.node_ids.first().copied().unwrap_or(usize::MAX);
                let last = tf.node_ids.last().copied().unwrap_or(usize::MAX);
                let first_coord = graph
                    .nodes
                    .get(first)
                    .map(|n| format!("{}-{}", n.start, n.end))
                    .unwrap_or_else(|| String::from("missing"));
                let last_coord = graph
                    .nodes
                    .get(last)
                    .map(|n| format!("{}-{}", n.start, n.end))
                    .unwrap_or_else(|| String::from("missing"));
                eprintln!(
                    "[TRACE_TF_WATCH seed={}]   t={} abund={:.4} weak={} longread={} nodes={} first={}({}) last={}({}) usepath={} longstart={} longend={}",
                    seed_idx,
                    tid,
                    tf.abundance,
                    tf.coverage_weak as u8,
                    tf.longread as u8,
                    tf.node_ids.len(),
                    first,
                    first_coord,
                    last,
                    last_coord,
                    tf.usepath,
                    tf.longstart,
                    tf.longend
                );
            }
            None => {
                eprintln!("[TRACE_TF_WATCH seed={}]   t={} missing", seed_idx, tid);
            }
        }
    }
}

fn trace_node_trf_watch(
    seed_idx: usize,
    label: &str,
    node_id: usize,
    graph: &Graph,
    watched: &[usize],
) {
    if watched.is_empty() || !trace_seed_diagnostics(seed_idx) {
        return;
    }
    let Some(node) = graph.nodes.get(node_id) else {
        return;
    };
    let mut hits = Vec::new();
    for &tid in watched {
        if node.trf_ids.contains(&tid) {
            hits.push(tid);
        }
    }
    eprintln!(
        "[TRACE_NODE_WATCH seed={}] {} node={}({}-{}) ntrf={} watched_hits={:?}",
        seed_idx,
        label,
        node_id,
        node.start,
        node.end,
        node.trf_ids.len(),
        hits
    );
}

/// Check whether a stored path contains the target intron (donor-acceptor pair).
/// `path_nodes` are graph node indices; an intron exists between consecutive nodes
/// where node_a.end < node_b.start.
fn path_has_intron(path_nodes: &[usize], graph: &Graph, target: (u64, u64)) -> bool {
    // Tracing helper: graph node boundaries may differ by a few bp due to junction snapping,
    // coordinate conventions, or local graph normalization. Use a small tolerance so the trace
    // remains informative.
    let tol: u64 = std::env::var("RUSTLE_TRACE_INTRON_TOL")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(5);
    for w in path_nodes.windows(2) {
        let (a, b) = (w[0], w[1]);
        if let (Some(na), Some(nb)) = (graph.nodes.get(a), graph.nodes.get(b)) {
            if na.end < nb.start
                && na.end.abs_diff(target.0) <= tol
                && nb.start.abs_diff(target.1) <= tol
            {
                return true;
            }
        }
    }
    false
}

/// Parse RUSTLE_TRACE_INTRON env var (format: "donor-acceptor").
fn parse_trace_intron() -> Option<(u64, u64)> {
    let val = std::env::var("RUSTLE_TRACE_INTRON").ok()?;
    let (a, b) = val.split_once('-')?;
    Some((a.trim().parse::<u64>().ok()?, b.trim().parse::<u64>().ok()?))
}

fn format_node_coord_list(graph: &Graph, nodes: &[usize]) -> String {
    let mut parts: Vec<String> = Vec::with_capacity(nodes.len());
    for &nid in nodes {
        match graph.nodes.get(nid) {
            Some(n) => parts.push(format!("{nid}({}-{})", n.start, n.end)),
            None => parts.push(format!("{nid}(missing)")),
        }
    }
    parts.join(" ")
}

fn format_exon_coord_list(exons: &[(u64, u64)]) -> String {
    exons
        .iter()
        .map(|(s, e)| format!("{s}-{e}"))
        .collect::<Vec<_>>()
        .join(",")
}

/// Optional first-exon left extension for source-linked starts.
/// Walk contiguous upstream parents (end == current.start) while preserving a simple
/// linear chain and return the earliest start coordinate.
fn source_contig_extension_start(graph: &Graph, start_nid: usize) -> Option<u64> {
    if start_nid >= graph.nodes.len() {
        return None;
    }
    let source_id = graph.source_id;
    let start_node = &graph.nodes[start_nid];
    if !start_node.parents.contains(source_id) {
        return None;
    }
    let mut cur = start_nid;
    let mut best = start_node.start;
    loop {
        let cnode = &graph.nodes[cur];
        let mut candidates: Vec<usize> = cnode
            .parents
            .ones()
            .filter(|&p| {
                p != source_id
                    && p < graph.nodes.len()
                    && graph.nodes[p].end == cnode.start
                    && graph.nodes[p].children.contains(cur)
            })
            .collect();
        if candidates.len() != 1 {
            break;
        }
        let p = candidates.pop().unwrap();
        if graph.nodes[p].children.count_ones() != 1 {
            break;
        }
        best = best.min(graph.nodes[p].start);
        cur = p;
    }
    Some(best)
}

/// Optional last-exon right extension for sink-linked ends.
/// Walk contiguous downstream children (start == current.end) while preserving a simple
/// linear chain and return the furthest end coordinate.
fn sink_contig_extension_end(graph: &Graph, end_nid: usize) -> Option<u64> {
    if end_nid >= graph.nodes.len() {
        return None;
    }
    let sink_id = graph.sink_id;
    let end_node = &graph.nodes[end_nid];
    if !end_node.children.contains(sink_id) {
        return None;
    }
    let mut cur = end_nid;
    let mut best = end_node.end;
    loop {
        let cnode = &graph.nodes[cur];
        let mut candidates: Vec<usize> = cnode
            .children
            .ones()
            .filter(|&c| {
                c != sink_id
                    && c < graph.nodes.len()
                    && graph.nodes[c].start == cnode.end
                    && graph.nodes[c].parents.contains(cur)
            })
            .collect();
        if candidates.len() != 1 {
            break;
        }
        let c = candidates.pop().unwrap();
        if graph.nodes[c].parents.count_ones() != 1 {
            break;
        }
        best = best.max(graph.nodes[c].end);
        cur = c;
    }
    Some(best)
}

fn node_can_reach(graph: &Graph, from: usize, to: usize) -> bool {
    if from == to {
        return true;
    }
    graph
        .nodes
        .get(from)
        .and_then(|n| n.childpat.as_ref())
        .map(|pat| pat.contains(to))
        .unwrap_or(false)
}

fn trace_checktrf_enqueue(
    t: usize,
    reason: &str,
    saved_abundance: f64,
    tf: &GraphTransfrag,
    tf_nodes: &[usize],
    graph: &Graph,
    path_nodes: Option<&[usize]>,
    back_ok: Option<bool>,
    fwd_ok: Option<bool>,
    flux: Option<f64>,
    coverage: Option<f64>,
    tocheck: Option<bool>,
) {
    let traced_tf = tf_nodes.iter().any(|&nid| trace_locus_active(graph, nid));
    let traced_path = path_nodes
        .map(|nodes| nodes.iter().any(|&nid| trace_locus_active(graph, nid)))
        .unwrap_or(false);
    if !traced_tf && !traced_path {
        return;
    }

    let path_nodes = path_nodes.unwrap_or(&[]);
    let back_str = back_ok.map(|v| if v { "1" } else { "0" }).unwrap_or("-");
    let fwd_str = fwd_ok.map(|v| if v { "1" } else { "0" }).unwrap_or("-");
    let flux_str = flux
        .map(|v| format!("{v:.4}"))
        .unwrap_or_else(|| String::from("-"));
    let cov_str = coverage
        .map(|v| format!("{v:.4}"))
        .unwrap_or_else(|| String::from("-"));
    let tocheck_str = tocheck.map(|v| if v { "1" } else { "0" }).unwrap_or("-");

    eprintln!(
        "[TRACE_CHECKTRF_ADD] t={} reason={} saved_abund={:.4} abund={:.4} guide={} longstart={} longend={} back={} fwd={} flux={} cov={} tocheck={} tf_nodes={} tf_coords={} path_nodes={} path_coords={}",
        t,
        reason,
        saved_abundance,
        tf.abundance,
        tf.guide,
        tf.longstart,
        tf.longend,
        back_str,
        fwd_str,
        flux_str,
        cov_str,
        tocheck_str,
        tf_nodes.len(),
        format_node_coord_list(graph, tf_nodes),
        path_nodes.len(),
        format_node_coord_list(graph, path_nodes)
    );
}

fn trace_checktrf_seed_snapshot(
    t: usize,
    tf: &GraphTransfrag,
    tf_nodes: &[usize],
    kept_paths: &[(Vec<usize>, f64, bool, usize)],
    graph: &Graph,
    out: &[Transcript],
) {
    if !tf_nodes.iter().any(|&nid| trace_locus_active(graph, nid)) {
        return;
    }
    eprintln!(
        "[TRACE_CHECKTRF_IN] t={} abund={:.4} guide={} longstart={} longend={} tf_nodes={} coords={}",
        t,
        tf.abundance,
        tf.guide,
        tf.longstart,
        tf.longend,
        tf_nodes.len(),
        format_node_coord_list(graph, tf_nodes)
    );
    eprintln!(
        "[TRACE_CHECKTRF_KEEPSET] t={} keeptrf={}",
        t,
        kept_paths.len()
    );
    for (kidx, (keep_nodes, keep_cov, keep_guide, out_idx)) in kept_paths.iter().enumerate() {
        let tx = out.get(*out_idx);
        let tx_cov = tx.map(|p| p.coverage).unwrap_or(0.0);
        let tx_exons = tx
            .map(|p| format_exon_coord_list(&p.exons))
            .unwrap_or_else(|| String::from("(missing_out)"));
        eprintln!(
            "[TRACE_CHECKTRF_KEEP] t={} k={} keep_cov={:.4} guide={} out_idx={} tx_cov={:.4} nodes={} exons={}",
            t,
            kidx,
            keep_cov,
            keep_guide,
            out_idx,
            tx_cov,
            format_node_coord_list(graph, keep_nodes),
            tx_exons
        );
    }
}

fn trace_checktrf_match_snapshot(
    t: usize,
    tmatch: &[usize],
    abundancesum: f64,
    kept_paths: &[(Vec<usize>, f64, bool, usize)],
    graph: &Graph,
    out: &[Transcript],
) {
    if tmatch.is_empty() {
        return;
    }
    let traced = tmatch.iter().copied().any(|k| {
        kept_paths
            .get(k)
            .map(|(nodes, _, _, _)| nodes.iter().any(|&nid| trace_locus_active(graph, nid)))
            .unwrap_or(false)
    });
    if !traced {
        return;
    }
    eprintln!(
        "[TRACE_CHECKTRF_OUT] t={} mode=matched kept_paths={:?} abundsum={:.4}",
        t, tmatch, abundancesum
    );
    for &k in tmatch {
        let Some((keep_nodes, keep_cov, keep_guide, out_idx)) = kept_paths.get(k) else {
            continue;
        };
        let tx = out.get(*out_idx);
        let tx_cov = tx.map(|p| p.coverage).unwrap_or(0.0);
        let tx_exons = tx
            .map(|p| format_exon_coord_list(&p.exons))
            .unwrap_or_else(|| String::from("(missing_out)"));
        eprintln!(
            "[TRACE_CHECKTRF_OUT_KEEP] t={} k={} keep_cov={:.4} guide={} out_idx={} tx_cov={:.4} nodes={} exons={}",
            t,
            k,
            keep_cov,
            keep_guide,
            out_idx,
            tx_cov,
            format_node_coord_list(graph, keep_nodes),
            tx_exons
        );
    }
}

/// Outcome of a single seed transfrag during extraction (for deep tracing).
#[derive(Debug, Clone)]
pub enum SeedOutcome {
    /// Skipped early: empty nodes, weak, not_seed, nascent_guide, short-read low support, empty real nodes.
    Skipped(&'static str),
    /// back_to_source_fast_long returned false → deferred to checktrf.
    BackToSourceFail,
    /// fwd_to_sink_fast_long returned false → deferred to checktrf.
    FwdToSinkFail,
    /// LR witness check rejected the path (unwitnessed consecutive splice pair).
    UnwitnessedSplice,
    /// Hard start/end boundaries don't match extended path + low abundance → checktrf.
    HardBoundaryMismatch,
    /// Flow returned 0 → deferred to checktrf.
    ZeroFlux,
    /// Coverage below EPSILON gate after flow.
    LowCoverage(f64),
    /// Eonly mode rejected non-guide prediction.
    EonlyNonGuide,
    /// Exons empty or transcript too short.
    TooShort,
    /// Checktrf: abundance below readthr gate.
    ChecktrfReadthr,
    /// Checktrf: skipped in eonly mode (non-guide transfrag).
    ChecktrfEonlySkip,
    /// Checktrf: redistributed to existing kept path (absorbed).
    ChecktrfRedistributed,
    /// Checktrf: rescued as independent transcript.
    ChecktrfRescued,
    /// Checktrf: incomplete path edges (rescue failed).
    ChecktrfIncomplete,
    /// Checktrf: rescue exons empty or too short or zero coverage.
    ChecktrfRescueFail,
    /// Successfully stored as transcript at output index.
    Stored(usize),
}

#[derive(Debug, Clone, Copy, Default)]
pub struct LongRecSummary {
    pub attempted: usize,
    pub succeeded: usize,
    pub fallback: usize,
    pub back_fail: usize,
    pub fwd_fail: usize,
    pub back_unreachable_minpath: usize,
    pub back_no_reach: usize,
    pub back_no_choice: usize,
    pub back_exclude_no_support: usize,
    pub fwd_unreachable_maxpath: usize,
    pub fwd_no_reach: usize,
    pub fwd_no_choice: usize,
    pub fwd_exclude_no_support: usize,
}

/// A single transcript: chrom, strand, exons (start,end), coverage, optional TPM/FPKM.
#[derive(Debug, Clone, Default)]
pub struct Transcript {
    pub chrom: String,
    pub strand: char,
    pub exons: Vec<(u64, u64)>,
    /// Final transcript coverage reported by assembly.
    ///
    /// For long-read extraction this is derived from depleted node usage
    /// (`nodeflux_abs * noderate`) and normalized by transcript length; for
    /// short-read extraction it is the analogous short-flow coverage.
    pub coverage: f64,
    pub exon_cov: Vec<f64>,
    /// TPM (transcripts per million); set by compute_tpm.
    pub tpm: f64,
    /// FPKM (fragments per kilobase million); set by compute_tpm.
    pub fpkm: f64,
/// Source tag for GTF and variant logic: `guide:*`, `flow` (long-read max-flow), `short_flow`,
    /// `checktrf_rescue`, `junction_path`, `ref_chain*`, rescues, or merge names from
    pub source: Option<String>,
    /// tlen<0 sentinel: true when prediction originates from long-read extraction.
    pub is_longread: bool,
    /// `p->longcov`: original transfrag abundance before flow depletion.
    pub longcov: f64,
    /// Strand-weighted bpcov sum over transcript exons (gene_abundance bpcov contribution).
    /// Set to 0.0 when bpcov is unavailable. Used in write_gene_abundance to approximate per-base coverage.
    pub bpcov_cov: f64,
    /// Forced transcript_id for GTF output. When set (eonly zero-cov guides), bypasses STRG.X.Y
    /// auto-numbering and uses the original guide transcript ID (guides[i]->getID()).
    pub transcript_id: Option<String>,
    /// Forced gene_id for GTF output. Used with transcript_id for eonly zero-cov guides.
    pub gene_id: Option<String>,
    /// reference_id GTF attribute: the guide transcript ID that this assembled transcript matches.
    pub ref_transcript_id: Option<String>,
    /// ref_gene_id GTF attribute: gene ID of the matched guide transcript.
    pub ref_gene_id: Option<String>,
    /// Parity: transcript terminates at a verified hard start node.
    pub hardstart: bool,
    /// Parity: transcript terminates at a verified hard end node.
    pub hardend: bool,
    // ── Variation graph fields (populated in --vg mode) ──────────────────────
    /// Family group ID (None = singleton, not part of a gene family).
    pub vg_family_id: Option<usize>,
    /// Copy index within the family group.
    pub vg_copy_id: Option<usize>,
    /// Number of copies in the family group.
    pub vg_family_size: Option<usize>,
}

#[inline]
fn gtf_source_long_flow(guide_tid: &Option<String>) -> Option<String> {
    guide_tid
        .as_ref()
        .map(|gid| format!("guide:{gid}"))
        .or_else(|| Some(String::from("flow")))
}

#[inline]
fn gtf_source_checktrf_rescue(guide_tid: &Option<String>) -> Option<String> {
    guide_tid
        .as_ref()
        .map(|gid| format!("guide:{gid}"))
        .or_else(|| Some(String::from("checktrf_rescue")))
}

#[inline]
fn gtf_source_short_flow(guide_tid: &Option<String>) -> Option<String> {
    guide_tid
        .as_ref()
        .map(|gid| format!("guide:{gid}"))
        .or_else(|| Some(String::from("short_flow")))
}

impl Transcript {
    /// Convert to `CPrediction` prediction object.
    pub fn to_cprediction(&self, geneno: i32) -> CPrediction {
        let start = self.exons.first().map(|e| e.0).unwrap_or(0);
        let end = self.exons.last().map(|e| e.1).unwrap_or(start);
        let exonic_len_u64 = self
            .exons
            .iter()
            .map(|(s, e)| e.saturating_sub(*s))
            .sum::<u64>();
        let exonic_len_i32 = exonic_len_u64.min(i32::MAX as u64) as i32;
        let tlen = if self.is_longread {
            -exonic_len_i32
        } else {
            exonic_len_i32
        };
        let exoncov = if self.exon_cov.len() == self.exons.len() {
            self.exon_cov.clone()
        } else {
            vec![self.coverage; self.exons.len()]
        };

        CPrediction {
            geneno,
            t_eq: self.ref_transcript_id.clone(),
            start,
            end,
            cov: self.coverage,
            longcov: self.longcov,
            strand: self.strand,
            tlen,
            flag: true,
            linkpred: None,
            exons: self.exons.clone(),
            exoncov,
            mergename: self.source.clone(),
            hardstart: self.hardstart,
            hardend: self.hardend,
        }
    }

    /// Build a transcript-shaped `CPrediction`.
    pub fn from_cprediction(chrom: String, pred: &CPrediction) -> Self {
        let exons = if pred.exons.is_empty() {
            vec![(pred.start, pred.end)]
        } else {
            pred.exons.clone()
        };
        let exon_cov = if pred.exoncov.len() == exons.len() {
            pred.exoncov.clone()
        } else {
            vec![pred.cov; exons.len()]
        };
        Self {
            chrom,
            strand: pred.strand,
            exons,
            coverage: pred.cov,
            exon_cov,
            tpm: 0.0,
            fpkm: 0.0,
            source: pred.mergename.clone(),
            is_longread: pred.tlen < 0 || pred.longcov > 0.0,
            longcov: pred.longcov,
            bpcov_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: pred.t_eq.clone(),
            ref_gene_id: None,
            hardstart: pred.hardstart,
            hardend: pred.hardend,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None,
        }
    }
}

/// Convert assembled transcripts to prediction objects.
pub fn transcripts_to_cpredictions(
    transcripts: &[Transcript],
    start_geneno: i32,
) -> Vec<CPrediction> {
    transcripts
        .iter()
        .enumerate()
        .map(|(i, tx)| tx.to_cprediction(start_geneno.saturating_add(i as i32)))
        .collect()
}

/// Convert prediction objects back to transcripts for Rust pipeline stages.
pub fn cpredictions_to_transcripts(chrom: &str, preds: &[CPrediction]) -> Vec<Transcript> {
    preds
        .iter()
        .map(|p| Transcript::from_cprediction(chrom.to_string(), p))
        .collect()
}

#[derive(Clone)]
struct GuideFlowState {
    pattern: GBitVec,
    abundance: f64,
    node_prior: HashMap<usize, f64>,
}

fn build_guide_node_prior(path: &[usize], nodeflux: &[f64], graph: &Graph) -> HashMap<usize, f64> {
    let mut out = Default::default();
    if path.len() < 3 {
        return out;
    }
    for (k, &nid) in path
        .iter()
        .enumerate()
        .skip(1)
        .take(path.len().saturating_sub(2))
    {
        let flux = nodeflux.get(k - 1).copied().unwrap_or(0.0);
        if flux <= EPS {
            continue;
        }
        let nlen = graph.nodes.get(nid).map(|n| n.length()).unwrap_or(0).max(1) as f64;
        out.insert(nid, flux * nlen);
    }
    out
}

fn intron_chain_from_exons(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
    if exons.len() < 2 {
        return Vec::new();
    }
    let mut out = Vec::with_capacity(exons.len().saturating_sub(1));
    for w in exons.windows(2) {
        let (_a0, a1) = w[0];
        let (b0, _b1) = w[1];
        if a1 < b0 {
            out.push((a1, b0));
        }
    }
    out
}

fn trace_intron_tol() -> u64 {
    std::env::var("RUSTLE_TRACE_INTRON_TOL")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(5)
}

fn intron_eq_tol(a: (u64, u64), b: (u64, u64), tol: u64) -> bool {
    a.0.abs_diff(b.0) <= tol && a.1.abs_diff(b.1) <= tol
}

fn intron_chains_match_tol(a: &[(u64, u64)], b: &[(u64, u64)], tol: u64) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter()
        .copied()
        .zip(b.iter().copied())
        .all(|(x, y)| intron_eq_tol(x, y, tol))
}

fn intron_chain_contains_all_tol(
    sup: &[(u64, u64)],
    sub: &[(u64, u64)],
    tol: u64,
) -> bool {
    sub.iter()
        .copied()
        .all(|needle| sup.iter().copied().any(|h| intron_eq_tol(h, needle, tol)))
}

/// Approximate strict assignment pre-EM step:
/// gather node-level prior mass from transfrags compatible with current guide pattern
/// and not claimed by previously processed guide patterns.
fn build_strict_guide_node_prior(
    path: &[usize],
    guide_pattern: &GBitVec,
    previous_guides: &[GuideFlowState],
    transfrags: &[GraphTransfrag],
    graph: &Graph,
) -> HashMap<usize, f64> {
    let mut out: HashMap<usize, f64> = Default::default();
    if path.len() < 3 {
        return out;
    }
    for &nid in path.iter().skip(1).take(path.len().saturating_sub(2)) {
        let Some(node) = graph.nodes.get(nid) else {
            continue;
        };
        let nlen = node.length().max(1) as f64;
        let mut mass = 0.0f64;
        for &tidx in &node.trf_ids {
            let Some(tf) = transfrags.get(tidx) else {
                continue;
            };
            if tf.abundance <= EPS {
                continue;
            }
            if !guide_pattern.contains_pattern(&tf.pattern) {
                continue;
            }
            let claimed_by_prev = previous_guides
                .iter()
                .any(|g| g.pattern.contains_pattern(&tf.pattern));
            if !claimed_by_prev {
                mass += tf.abundance;
            }
        }
        if mass > 0.0 {
            out.insert(nid, mass * nlen);
        }
    }
    out
}

/// guides_pushmaxflow pre-EM inspired node competition accounting for current guide.
/// Builds node priors from strict (unique-compatible) and loose (shared-compatible) transfrag counts.
fn build_cnodeguide_node_prior(
    path: &[usize],
    guide_pattern: &GBitVec,
    previous_guides: &[GuideFlowState],
    transfrags: &[GraphTransfrag],
    graph: &Graph,
) -> HashMap<usize, f64> {
    let mut out: HashMap<usize, f64> = Default::default();
    if path.len() < 3 {
        return out;
    }
    let first_real = path.get(1).copied().unwrap_or(0);
    let last_real = path
        .get(path.len().saturating_sub(2))
        .copied()
        .unwrap_or(first_real);

    for &nid in path.iter().skip(1).take(path.len().saturating_sub(2)) {
        let Some(node) = graph.nodes.get(nid) else {
            continue;
        };
        let nlen = node.length().max(1) as f64;

        let mut incount = 0.0f64;
        let mut outcount = 0.0f64;
        let mut guide_incount = 0.0f64;
        let mut guide_outcount = 0.0f64;
        let mut strict_in = 0.0f64;
        let mut strict_out = 0.0f64;
        let mut loose_in = 0.0f64;
        let mut loose_out = 0.0f64;
        let mut compat_seen = false;

        for tf in transfrags {
            if tf.abundance <= EPS || tf.node_ids.is_empty() {
                continue;
            }
            if !tf.node_id_set.contains(nid) {
                continue;
            }
            if !guide_pattern.contains_pattern(&tf.pattern) {
                continue;
            }
            compat_seen = true;
            let prev_compat = previous_guides
                .iter()
                .filter(|g| g.pattern.contains_pattern(&tf.pattern))
                .count();
            let comp_total = prev_compat + 1;

            let starts_here = tf.node_ids.first().copied() == Some(nid);
            let ends_here = tf.node_ids.last().copied() == Some(nid);
            let through_here = !starts_here && !ends_here && tf.pattern.get_bit(nid);

            if ends_here {
                incount += tf.abundance;
                if comp_total == 1 {
                    strict_in += tf.abundance;
                } else {
                    guide_incount += tf.abundance;
                    loose_in += tf.abundance;
                }
            } else if starts_here {
                outcount += tf.abundance;
                if comp_total == 1 {
                    strict_out += tf.abundance;
                } else {
                    guide_outcount += tf.abundance;
                    loose_out += tf.abundance;
                }
            } else if through_here {
                incount += tf.abundance;
                outcount += tf.abundance;
                if comp_total == 1 {
                    strict_in += tf.abundance;
                    strict_out += tf.abundance;
                } else {
                    guide_incount += tf.abundance;
                    guide_outcount += tf.abundance;
                    loose_in += tf.abundance;
                    loose_out += tf.abundance;
                }
            }
        }

        if !compat_seen {
            continue;
        }

        let rate = if incount > EPS && outcount > EPS {
            incount / outcount
        } else {
            1.0
        };

        let terminal_in = nid == first_real;
        let terminal_out = nid == last_real;
        let mut strict_mass = strict_in + rate * strict_out;
        if incount > EPS && outcount > EPS {
            if terminal_in {
                strict_mass += strict_in;
            }
            if terminal_out {
                strict_mass += rate * strict_out;
            }
        }

        let mut gcount = loose_in + rate * loose_out;
        if incount > EPS && outcount > EPS {
            if terminal_in {
                gcount += loose_in;
            }
            if terminal_out {
                gcount += rate * loose_out;
            }
        }

        let mut sumtrcount = guide_incount + rate * guide_outcount;
        let sumstrict = sumtrcount + strict_mass;
        if strict_mass > EPS {
            if sumstrict > EPS {
                sumtrcount *= sumtrcount / sumstrict;
            } else {
                sumtrcount = 0.0;
            }
        }

        let mut prior = strict_mass + gcount + sumtrcount;
        if prior <= EPS {
            prior = 1.0; // winner-takes-all fallback for compatible but uncovered competition.
        }
        out.insert(nid, prior * nlen);
    }
    out
}

#[derive(Clone)]
struct NodeTrPattern {
    guides: Vec<usize>,
    abund: f64,
}

#[derive(Clone)]
struct NodeGuideEmInfo {
    node_id: usize,
    node_len: f64,
    strict_current: f64,
    patterns: Vec<NodeTrPattern>,
}

/// CNodeGuide-style pattern-level EM for current guide vs known competitors (previous guides).
/// Returns (node_prior_for_current, refined_current_cov_mass).
fn cnodeguide_em_prior(
    path: &[usize],
    current_pattern: &GBitVec,
    current_cov_mass: f64,
    previous_guides: &[GuideFlowState],
    transfrags: &[GraphTransfrag],
    graph: &Graph,
) -> (HashMap<usize, f64>, f64) {
    if path.len() < 3 {
        return (Default::default(), current_cov_mass.max(0.0));
    }
    let mut competitors: Vec<&GBitVec> = Vec::with_capacity(previous_guides.len() + 1);
    competitors.push(current_pattern);
    for g in previous_guides {
        competitors.push(&g.pattern);
    }
    let m = competitors.len();
    let mut infos: Vec<NodeGuideEmInfo> = Vec::new();

    for &nid in path.iter().skip(1).take(path.len().saturating_sub(2)) {
        let Some(node) = graph.nodes.get(nid) else {
            continue;
        };
        let mut strict_current = 0.0f64;
        let mut patmap: HashMap<Vec<usize>, f64> = Default::default();
        for &tidx in &node.trf_ids {
            let Some(tf) = transfrags.get(tidx) else {
                continue;
            };
            if tf.abundance <= EPS || tf.node_ids.is_empty() {
                continue;
            }
            // Guide/transfrag compatibility used in CNodeGuide block.
            if !current_pattern.contains_pattern(&tf.pattern) {
                continue;
            }
            let mut comp: Vec<usize> = Vec::new();
            for (gi, gp) in competitors.iter().enumerate() {
                if gp.contains_pattern(&tf.pattern) {
                    comp.push(gi);
                }
            }
            if comp.is_empty() || !comp.contains(&0) {
                continue;
            }
            if comp.len() == 1 && comp[0] == 0 {
                strict_current += tf.abundance;
            } else {
                comp.sort_unstable();
                *patmap.entry(comp).or_insert(0.0) += tf.abundance;
            }
        }
        if strict_current <= EPS && patmap.is_empty() {
            continue;
        }
        let patterns = patmap
            .into_iter()
            .map(|(guides, abund)| NodeTrPattern { guides, abund })
            .collect();
        infos.push(NodeGuideEmInfo {
            node_id: nid,
            node_len: node.length().max(1) as f64,
            strict_current,
            patterns,
        });
    }
    if infos.is_empty() {
        return (Default::default(), current_cov_mass.max(0.0));
    }

    let mut cov = vec![0.0f64; m];
    cov[0] = current_cov_mass.max(1.0);
    for (i, g) in previous_guides.iter().enumerate() {
        cov[i + 1] = g.abundance.max(1.0);
    }

    let mut current_node_prior: HashMap<usize, f64> = Default::default();
    for _ in 0..10 {
        let mut next_cov = vec![0.0f64; m];
        current_node_prior.clear();
        for info in &infos {
            let mut alloc = vec![0.0f64; m];
            alloc[0] += info.strict_current;
            for p in &info.patterns {
                let mut denom = 0.0f64;
                for &gi in &p.guides {
                    denom += cov[gi].max(0.0);
                }
                if denom > EPS {
                    for &gi in &p.guides {
                        alloc[gi] += p.abund * cov[gi].max(0.0) / denom;
                    }
                } else {
                    // Winner-takes-all fallback (reference-style no prior mass).
                    alloc[p.guides[0]] += p.abund;
                }
            }
            for gi in 0..m {
                next_cov[gi] += alloc[gi] * info.node_len;
            }
            if alloc[0] > EPS {
                current_node_prior.insert(info.node_id, alloc[0] * info.node_len);
            }
        }
        let mut converged = true;
        for gi in 0..m {
            if (next_cov[gi] - cov[gi]).abs() > 1.0 {
                converged = false;
                break;
            }
        }
        cov = next_cov;
        if converged {
            break;
        }
    }
    (current_node_prior, cov[0].max(0.0))
}

/// Rawreads-mode extraction: emit transcript models directly from active transfrags
/// without max-flow subtraction (rawreads branch intent).
pub fn extract_rawreads_transcripts(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    bundle_chrom: &str,
    bundle_strand: char,
    _bundle_id: &str,
    config: &RunConfig,
) -> Vec<Transcript> {
    let mut out: Vec<Transcript> = Vec::new();
    let mut seen_paths: HashSet<Vec<usize>> = Default::default();
    for tf in transfrags {
        if tf.abundance <= EPS || tf.node_ids.len() < 2 || tf.weak != 0 {
            continue;
        }
        let key: Vec<usize> = tf
            .node_ids
            .iter()
            .copied()
            .filter(|&n| n != graph.source_id && n != graph.sink_id)
            .collect();
        if key.is_empty() {
            continue;
        }
        if !seen_paths.insert(key.clone()) {
            continue;
        }
        let exons = collect_path(
            &tf.node_ids,
            graph,
            (tf.longstart > 0).then_some(tf.longstart),
            (tf.longend > 0).then_some(tf.longend),
            true,
        );
        if exons.is_empty() {
            continue;
        }
        let tlen: u64 = exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
        if tlen < config.min_transcript_length {
            continue;
        }
        out.push(Transcript {
            chrom: bundle_chrom.to_string(),
            strand: bundle_strand,
            exons: exons.clone(),
            coverage: tf.abundance.max(0.0),
            exon_cov: vec![tf.abundance.max(0.0); exons.len()],
            tpm: 0.0,
            fpkm: 0.0,
            source: Some("rawreads".to_string()),
            is_longread: config.long_reads,
            longcov: 0.0,
            bpcov_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: None,
            ref_gene_id: None,
            hardstart: false,
            hardend: false,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None,
        });
    }
    out
}

/// Eonly (guide-only) transcript extraction (guides_pushmaxflow path).
/// Called when `config.eonly && !guides.is_empty()`. Delegates to `extract_transcripts`
/// with the graph/transfrags containing only guide transfrags.
pub fn extract_eonly_transcripts(
    graph: &mut Graph,
    transfrags: &mut Vec<GraphTransfrag>,
    bundle_chrom: &str,
    bundle_strand: char,
    bundle_id: &str,
    config: &RunConfig,
) -> Vec<Transcript> {
    // eonly: only guide-tagged transfrags participate.
    // For now, delegate to extract_transcripts (guides_pushmaxflow result is already in transfrags).
    extract_transcripts(
        graph,
        transfrags,
        bundle_chrom,
        bundle_strand,
        bundle_id,
        config,
        false,
        None,
        None,
    )
}

/// Short-read path extraction (parse_trf port): consume transfrags by short-read abundance order.
/// Called instead of extract_transcripts when running in short-read mode (!config.long_reads).
pub fn extract_shortread_transcripts(
    graph: &Graph,
    transfrags: &mut [GraphTransfrag],
    bundle_chrom: &str,
    bundle_strand: char,
    _bundle_id: &str,
    config: &RunConfig,
) -> Vec<Transcript> {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;

    // Build sorted index: short-read candidates sorted by (srabund + abundance) desc, node count desc.
    let mut indices: Vec<usize> = transfrags
        .iter()
        .enumerate()
        .filter(|(_, tf)| tf.srabund > EPS || !tf.longread)
        .map(|(i, _)| i)
        .collect();
    indices.sort_unstable_by(|&i, &j| {
        let a = &transfrags[i];
        let b = &transfrags[j];
        let eff_a = a.srabund + a.abundance;
        let eff_b = b.srabund + b.abundance;
        let ab = eff_b
            .partial_cmp(&eff_a)
            .unwrap_or(std::cmp::Ordering::Equal);
        if ab != std::cmp::Ordering::Equal {
            return ab;
        }
        b.node_ids.len().cmp(&a.node_ids.len())
    });

    let mut out: Vec<Transcript> = Vec::new();

    for idx in indices {
        if idx >= transfrags.len() {
            continue;
        }
        let eff = transfrags[idx].srabund + transfrags[idx].abundance;
        if eff <= EPS {
            continue;
        }

        // Build initial path from transfrag's real (non-source/sink) nodes.
        let mut path: Vec<usize> = transfrags[idx]
            .node_ids
            .iter()
            .copied()
            .filter(|&n| n != source_id && n != sink_id)
            .collect();
        if path.is_empty() {
            continue;
        }

        // Extend left (back_to_source_fast) and right (fwd_to_sink_fast).
        extend_path_left(&mut path, graph, transfrags);
        extend_path_right(&mut path, graph, transfrags);

        // Wrap path with source and sink.
        path.insert(0, source_id);
        path.push(sink_id);

        // Push flow and subtract from transfrags.
        let (flux, _nodeflux, _full) =
            push_max_flow_seeded_full(&path, transfrags, graph, false, Some(idx));

        if flux <= EPS {
            continue;
        }

        // Build exons from the path.
        let exons = collect_path(&path, graph, None, None, false);
        if exons.is_empty() {
            continue;
        }
        let tlen: u64 = exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
        if tlen < config.min_transcript_length {
            continue;
        }
        let n_exons = exons.len();
        let first_node = path[0];
        let last_node = *path.last().unwrap_or(&first_node);
        out.push(Transcript {
            chrom: bundle_chrom.to_string(),
            strand: bundle_strand,
            exons,
            coverage: flux.max(0.0),
            exon_cov: vec![flux.max(0.0); n_exons],
            tpm: 0.0,
            fpkm: 0.0,
            source: gtf_source_short_flow(&transfrags[idx].guide_tid),
            is_longread: false,
            longcov: 0.0,
            bpcov_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: None,
            ref_gene_id: None,
            hardstart: graph.nodes.get(first_node).map(|n| n.hardstart).unwrap_or(false),
            hardend: graph.nodes.get(last_node).map(|n| n.hardend).unwrap_or(false),
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None,
        });
    }

    out
}

// --- Constants header / ---
/// Min unaligned-tail reads to mark hardstart/hardend (POLY_TAIL_STOP_COUNT).
const POLY_TAIL_STOP_COUNT: u16 = 8;
/// CHI_THR: maximum tolerated unmatched internal distance for best_trf_match.
const CHI_THR_BP: i64 = 50;
/// Constants used by direct long-recursion port.
const DROP_FACTOR: f64 = 0.5;
const ERROR_PERC: f64 = 0.1;
const EPS: f64 = crate::constants::FLOW_EPSILON;
/// When checktrf rescues a failed transfrag, only redistribute it onto an existing kept path
/// if the intron chain agrees within this tolerance. This prevents collapsing distinct
/// alternative splice-site isoforms into a single dominant path.
#[allow(dead_code)]
#[allow(dead_code)]
const CHECKTRF_REDISTRIBUTE_INTRON_TOL: u64 = 5;

/// Build exons from a path of node IDs: merge overlapping/contiguous nodes (half-open).
/// optionally trim first/last exon by longstart/longend (collect_path; trim only, no extend).
pub fn collect_path(
    path: &[usize],
    graph: &Graph,
    longstart: Option<u64>,
    longend: Option<u64>,
    apply_longstart_longend: bool,
) -> Vec<(u64, u64)> {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let mut exons: Vec<(u64, u64)> = Vec::new();
    for (_i, &nid) in path.iter().enumerate() {
        if nid == source_id || nid == sink_id {
            continue;
        }
        let node = match graph.nodes.get(nid) {
            Some(n) => n,
            None => continue,
        };
        if node.end <= node.start {
            continue;
        }
        let (node_start, node_end) = (node.start, node.end);
        if let Some(last) = exons.last() {
            if node_start <= last.1 {
                let last = exons.last_mut().unwrap();
                last.1 = last.1.max(node_end);
                continue;
            }
        }
        exons.push((node_start, node_end));
    }
    if exons.is_empty() {
        return exons;
    }
    let orig_first_exon_start = exons[0].0;
    let orig_last_exon_end = exons.last().map(|e| e.1).unwrap_or(0);
    let mut source_extended = false;
    let mut sink_extended = false;
    if std::env::var_os("RUSTLE_SOURCE_CONTIG_EXTEND").is_some() {
        if let Some(first_real) = path
            .iter()
            .copied()
            .find(|&nid| nid != source_id && nid != sink_id)
        {
            if let Some(ext_start) = source_contig_extension_start(graph, first_real) {
                if ext_start < exons[0].0 {
                    exons[0].0 = ext_start;
                    source_extended = true;
                }
            }
        }
        if std::env::var_os("RUSTLE_SINK_CONTIG_EXTEND").is_some() {
            if let Some(last_real) = path
                .iter()
                .rev()
                .copied()
                .find(|&nid| nid != source_id && nid != sink_id)
            {
                if let Some(ext_end) = sink_contig_extension_end(graph, last_real) {
                    let li = exons.len() - 1;
                    if ext_end > exons[li].1 {
                        exons[li].1 = ext_end;
                        sink_extended = true;
                    }
                }
            }
        }
    }
    if apply_longstart_longend {
        if let Some(ls) = longstart {
            let (s, e) = exons[0];
            let left = if source_extended
                && std::env::var_os("RUSTLE_SOURCE_CONTIG_EXTEND").is_some()
                && ls >= orig_first_exon_start
            {
                // Keep source-contiguous extension even when longstart was trimmed to the
                // downstream node start by source-helper gating.
                s
            } else {
                s.max(ls).min(e)
            };
            if left < e {
                exons[0].0 = left;
            }
        }
        if let Some(le) = longend {
            let last_idx = exons.len() - 1;
            let (s, e) = exons[last_idx];
            let right = if sink_extended
                && std::env::var_os("RUSTLE_SINK_CONTIG_EXTEND").is_some()
                && le <= orig_last_exon_end
            {
                // Keep sink-contiguous extension even when longend was trimmed to a
                // pre-sink stub boundary by sink-helper gating.
                e
            } else {
                e.min(le).max(s)
            };
            if right > s {
                exons[last_idx].1 = right;
            }
        }
    }
    // Filter out any exons that became invalid after trimming
    exons.retain(|(s, e)| s < e);
    exons
}

/// Extend path toward source: pick parent with max supporting transfrag abundance.
/// path: first element is the leftmost real node (not source). Inserts nodes at front until source or no supported parent.
pub fn extend_path_left(path: &mut Vec<usize>, graph: &Graph, transfrags: &[GraphTransfrag]) {
    if path.is_empty() {
        return;
    }
    let unionpat = build_unionpat_active(transfrags, graph.pattern_size());
    let keept = init_keept_all_active(transfrags);
    if keept.is_empty() {
        return;
    }
    let nodecov: Vec<f64> = graph.nodes.iter().map(|n| n.nodecov.max(0.0)).collect();
    let mut leftpath = Vec::new();
    let left_seed = path[0];
    let mut visited: NodeSet = NodeSet::with_capacity(graph.n_nodes);
    extend_path_left_rec(
        left_seed,
        graph,
        transfrags,
        &unionpat,
        &keept,
        &mut leftpath,
        &nodecov,
        &mut visited,
    );
    if !leftpath.is_empty() {
        leftpath.reverse();
        let mut merged = leftpath;
        merged.extend(path.clone());
        merged.dedup();
        *path = merged;
    }
}

/// Extend path toward sink: pick child with max supporting transfrag abundance.
/// path: last element is the rightmost real node (not sink). Appends nodes until sink or no supported child.
pub fn extend_path_right(path: &mut Vec<usize>, graph: &Graph, transfrags: &[GraphTransfrag]) {
    if path.is_empty() {
        return;
    }
    let unionpat = build_unionpat_active(transfrags, graph.pattern_size());
    let keept = init_keept_all_active(transfrags);
    if keept.is_empty() {
        return;
    }
    let nodecov: Vec<f64> = graph.nodes.iter().map(|n| n.nodecov.max(0.0)).collect();
    let mut rightpath = Vec::new();
    let right_seed = *path.last().unwrap_or(&path[0]);
    let mut visited: NodeSet = NodeSet::with_capacity(graph.n_nodes);
    extend_path_right_rec(
        right_seed,
        graph,
        transfrags,
        &unionpat,
        &keept,
        &mut rightpath,
        &nodecov,
        &mut visited,
    );
    if !rightpath.is_empty() {
        path.extend(rightpath);
        path.dedup();
    }
}

fn first_last_real_nodes(
    tf: &GraphTransfrag,
    source_id: usize,
    sink_id: usize,
) -> Option<(usize, usize)> {
    let first = tf
        .node_ids
        .iter()
        .copied()
        .find(|&n| n != source_id && n != sink_id)?;
    let last = tf
        .node_ids
        .iter()
        .copied()
        .rfind(|&n| n != source_id && n != sink_id)?;
    Some((first, last))
}

fn edge_supported(tf: &GraphTransfrag, graph: &Graph, a: usize, b: usize) -> bool {
    graph
        .edge_bit_index(a, b)
        .map(|eid| tf.pattern.get_bit(eid))
        .unwrap_or(false)
}

fn build_unionpat_active(transfrags: &[GraphTransfrag], psize: usize) -> GBitVec {
    let mut unionpat = GBitVec::new(psize);
    for tf in transfrags {
        if tf.abundance > 0.0 && tf.weak == 0 {
            unionpat.or_assign(&tf.pattern);
        }
    }
    unionpat
}

fn init_keept_all_active(transfrags: &[GraphTransfrag]) -> Vec<usize> {
    transfrags
        .iter()
        .enumerate()
        .filter_map(|(i, tf)| {
            if tf.abundance > 0.0 && tf.weak == 0 {
                Some(i)
            } else {
                None
            }
        })
        .collect()
}

fn extend_path_left_rec(
    n: usize,
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    unionpat: &GBitVec,
    keept: &[usize],
    leftpath: &mut Vec<usize>,
    nodecov: &[f64],
    visited: &mut NodeSet,
) {
    if visited.contains(n) || leftpath.len() > graph.n_nodes {
        return;
    }
    visited.insert_grow(n);
    let source_id = graph.source_id;
    let Some(node) = graph.nodes.get(n) else {
        return;
    };

    let mut maxp: Option<usize> = None;
    let mut maxabund = 0.0f64;
    let mut tmax: Option<usize> = None; // best transfrag for current winner

    for p in node.parents.ones() {
        if p == source_id || !unionpat.get_bit(p) {
            continue;
        }
        let mut sumabund = 0.0f64;
        let mut tpar: Option<usize> = None; // best transfrag for this parent
        for &t in keept {
            let Some(tf) = transfrags.get(t) else {
                continue;
            };
            if tf.abundance <= 0.0 || !tf.pattern.get_bit(p) {
                continue;
            }
            if let Some((first_real, last_real)) =
                first_last_real_nodes(tf, source_id, graph.sink_id)
            {
                if last_real == p
                    || edge_supported(tf, graph, p, n)
                    || (first_real <= p && last_real > p)
                {
                    sumabund += tf.abundance;
                    replace_tmax(t, &mut tpar, transfrags, graph, nodecov);
                }
            }
        }
        // — weak-aware fork comparison
        if let Some(tp) = tpar {
            if sumabund > maxabund {
                let tp_weak = tf_weak(&transfrags[tp], graph, nodecov);
                if tmax.is_none()
                    || tmax
                        .map(|tm| tf_weak(&transfrags[tm], graph, nodecov))
                        .unwrap_or(true)
                    || !tp_weak
                {
                    maxabund = sumabund;
                    maxp = Some(p);
                    tmax = Some(tp);
                }
            } else if maxp.is_some() && (sumabund - maxabund).abs() < EPS {
                let tp_weak = tf_weak(&transfrags[tp], graph, nodecov);
                if tmax
                    .map(|tm| tf_weak(&transfrags[tm], graph, nodecov))
                    .unwrap_or(true)
                    || !tp_weak
                {
                    if nodecov.get(maxp.unwrap()).copied().unwrap_or(0.0)
                        < nodecov.get(p).copied().unwrap_or(0.0)
                    {
                        maxp = Some(p);
                        tmax = Some(tp);
                    }
                }
            }
        }
    }

    if let Some(p) = maxp {
        let mut ekeept: Vec<usize> = Vec::new();
        for &t in keept {
            let Some(tf) = transfrags.get(t) else {
                continue;
            };
            if tf.abundance <= 0.0 {
                continue;
            }
            let Some((first_real, last_real)) = first_last_real_nodes(tf, source_id, graph.sink_id)
            else {
                continue;
            };
            if tf.pattern.get_bit(p) {
                if last_real == p || edge_supported(tf, graph, p, n) {
                    ekeept.push(t);
                }
            } else if last_real < p || (first_real <= p && last_real > p) {
                ekeept.push(t);
            }
        }
        if !ekeept.is_empty() {
            leftpath.push(p);
            extend_path_left_rec(
                p, graph, transfrags, unionpat, &ekeept, leftpath, nodecov, visited,
            );
            return;
        }
    }

    let mut min_node = n;
    let mut fuzzy: Vec<usize> = Vec::new();
    for &t in keept {
        let Some(tf) = transfrags.get(t) else {
            continue;
        };
        let Some((first_real, last_real)) = first_last_real_nodes(tf, source_id, graph.sink_id)
        else {
            continue;
        };
        if last_real >= n && first_real < min_node {
            for &nid in &tf.node_ids {
                if nid == source_id || nid == graph.sink_id {
                    continue;
                }
                if nid == min_node {
                    break;
                }
                fuzzy.push(nid);
            }
            min_node = first_real;
        }
    }
    for &nid in fuzzy.iter().rev() {
        leftpath.push(nid);
    }
}

fn extend_path_right_rec(
    n: usize,
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    unionpat: &GBitVec,
    keept: &[usize],
    rightpath: &mut Vec<usize>,
    nodecov: &[f64],
    visited: &mut NodeSet,
) {
    if visited.contains(n) || rightpath.len() > graph.n_nodes {
        return;
    }
    visited.insert_grow(n);
    let sink_id = graph.sink_id;
    let Some(node) = graph.nodes.get(n) else {
        return;
    };

    let mut maxc: Option<usize> = None;
    let mut maxabund = 0.0f64;
    let mut tmax: Option<usize> = None; // best transfrag for current winner

    for c in node.children.ones() {
        if c == sink_id || !unionpat.get_bit(c) {
            continue;
        }
        let mut sumabund = 0.0f64;
        let mut tchild: Option<usize> = None; // best transfrag for this child
        for &t in keept {
            let Some(tf) = transfrags.get(t) else {
                continue;
            };
            if tf.abundance <= 0.0 || !tf.pattern.get_bit(c) {
                continue;
            }
            if let Some((first_real, last_real)) =
                first_last_real_nodes(tf, graph.source_id, sink_id)
            {
                if first_real == c
                    || edge_supported(tf, graph, n, c)
                    || (first_real < c && last_real >= c)
                {
                    sumabund += tf.abundance;
                    replace_tmax(t, &mut tchild, transfrags, graph, nodecov);
                }
            }
        }
        // — weak-aware fork comparison
        if let Some(tc) = tchild {
            if sumabund > maxabund {
                let tc_weak = tf_weak(&transfrags[tc], graph, nodecov);
                if tmax.is_none()
                    || tmax
                        .map(|tm| tf_weak(&transfrags[tm], graph, nodecov))
                        .unwrap_or(true)
                    || !tc_weak
                {
                    maxabund = sumabund;
                    maxc = Some(c);
                    tmax = Some(tc);
                }
            } else if maxc.is_some() && (sumabund - maxabund).abs() < EPS {
                let tc_weak = tf_weak(&transfrags[tc], graph, nodecov);
                if tmax
                    .map(|tm| tf_weak(&transfrags[tm], graph, nodecov))
                    .unwrap_or(true)
                    || !tc_weak
                {
                    if nodecov.get(maxc.unwrap()).copied().unwrap_or(0.0)
                        < nodecov.get(c).copied().unwrap_or(0.0)
                    {
                        maxc = Some(c);
                        tmax = Some(tc);
                    }
                }
            }
        }
    }

    if let Some(c) = maxc {
        let mut ekeept: Vec<usize> = Vec::new();
        for &t in keept {
            let Some(tf) = transfrags.get(t) else {
                continue;
            };
            if tf.abundance <= 0.0 {
                continue;
            }
            let Some((first_real, last_real)) = first_last_real_nodes(tf, graph.source_id, sink_id)
            else {
                continue;
            };
            if tf.pattern.get_bit(c) {
                if first_real == c || edge_supported(tf, graph, n, c) {
                    ekeept.push(t);
                }
            } else if first_real > c || (first_real < c && last_real >= c) {
                ekeept.push(t);
            }
        }
        if !ekeept.is_empty() {
            rightpath.push(c);
            extend_path_right_rec(
                c, graph, transfrags, unionpat, &ekeept, rightpath, nodecov, visited,
            );
            return;
        }
    }

    let mut max_node = n;
    let mut fuzzy: Vec<usize> = Vec::new();
    for &t in keept {
        let Some(tf) = transfrags.get(t) else {
            continue;
        };
        let Some((first_real, last_real)) = first_last_real_nodes(tf, graph.source_id, sink_id)
        else {
            continue;
        };
        if first_real <= n && last_real > max_node {
            let mut started = false;
            for &nid in &tf.node_ids {
                if nid == graph.source_id || nid == sink_id {
                    continue;
                }
                if nid == max_node {
                    started = true;
                    continue;
                }
                if started {
                    fuzzy.push(nid);
                }
            }
            max_node = last_real;
        }
    }
    rightpath.extend(fuzzy);
}

/// Add transfrag to path: OR pattern into pathpat, cov += abundance*len, update min/max node indices (add_transfrag_to_path).
/// Zeros tf.abundance. alltr receives the transfrag index for bookkeeping. Only adds if tf is on path (nodes ⊆ path).
pub fn add_transfrag_to_path(
    tf_idx: usize,
    tf: &mut GraphTransfrag,
    pathpat: &mut GBitVec,
    path: &[usize],
    graph: &Graph,
    min_node: &mut usize,
    max_node: &mut usize,
    cov: &mut f64,
    alltr: &mut Vec<usize>,
) {
    let path_set: HashSet<usize> = path.iter().copied().collect();
    let on_path = tf
        .node_ids
        .iter()
        .all(|&nid| nid == graph.source_id || nid == graph.sink_id || path_set.contains(&nid));
    if !on_path {
        return;
    }
    pathpat.or_assign(&tf.pattern);
    let len: u64 = tf
        .node_ids
        .iter()
        .filter(|&&nid| nid != graph.source_id && nid != graph.sink_id)
        .map(|&nid| graph.nodes.get(nid).map(|n| n.length()).unwrap_or(0))
        .sum();
    *cov += tf.abundance * (len as f64);
    alltr.push(tf_idx);
    if let (Some(&first), Some(&last)) = (tf.node_ids.first(), tf.node_ids.last()) {
        if first != graph.source_id && first < *min_node {
            *min_node = first;
        }
        if last != graph.sink_id && last > *max_node {
            *max_node = last;
        }
    }
    tf.abundance = 0.0;
}

/// Update path node coverage backward from transcript end (update_transcript_to_path_back).
/// path_incov[i] = inflow at node i, path_outcov[i] = outflow. Same length as path.
pub fn update_transcript_to_path_back(
    abundance: f64,
    trnode: &[usize],
    path: &[usize],
    path_incov: &mut [f64],
    path_outcov: &mut [f64],
) {
    if path.is_empty() || trnode.is_empty() {
        return;
    }
    let last_tr = *trnode.last().unwrap();
    let mut lastnode = last_tr;
    if path[0] < last_tr {
        lastnode = path[0];
    }
    if trnode.last() != path.last() {
        if let Some(last) = path_outcov.last_mut() {
            *last += abundance;
        }
    }
    for i in (0..path.len().saturating_sub(1)).rev() {
        if path[i] == lastnode {
            path_incov[i] += abundance;
            break;
        }
        path_incov[i] += abundance;
        path_outcov[i] += abundance;
    }
}

/// Update path node coverage forward from transcript start (update_transcript_to_path_fwd).
pub fn update_transcript_to_path_fwd(
    abundance: f64,
    trnode: &[usize],
    path: &[usize],
    path_incov: &mut [f64],
    path_outcov: &mut [f64],
) {
    if path.is_empty() || trnode.is_empty() {
        return;
    }
    if trnode[0] != path.last().copied().unwrap_or(0) {
        if let Some(last) = path_incov.last_mut() {
            *last += abundance;
        }
    }
    for i in (0..path.len().saturating_sub(1)).rev() {
        if path[i] == trnode[0] {
            path_outcov[i] += abundance;
            break;
        }
        path_incov[i] += abundance;
        path_outcov[i] += abundance;
    }
}

/// Update guide prediction exon coverage from path node coverage (update_guide_pred).
/// For each node in path, add nodeflux[i]*nodecov[node] to overlapping exons (by overlap length).
/// If nodeflux is None, uses 1.0 per node.
pub fn update_guide_pred(
    path: &[usize],
    graph: &Graph,
    nodecov: &[f64],
    exons: &[(u64, u64)],
    exoncov: &mut [f64],
    nodeflux: Option<&[f64]>,
) {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let gno = graph.n_nodes;
    if exoncov.len() != exons.len() {
        return;
    }
    let mut e = 0usize;
    let nex = exons.len();
    for (i, &nid) in path.iter().enumerate() {
        if nid == source_id || nid == sink_id || nid >= gno.saturating_sub(1) {
            continue;
        }
        let node = match graph.nodes.get(nid) {
            Some(n) => n,
            None => continue,
        };
        let flux = nodeflux.and_then(|f| f.get(i).copied()).unwrap_or(1.0);
        let cov = nodecov.get(nid).copied().unwrap_or(0.0);
        let addcov = flux * cov;
        while e < nex && exons[e].1 <= node.start {
            e += 1;
        }
        while e < nex && exons[e].0 < node.end {
            let (es, ee) = exons[e];
            let ovplen = overlap_len_half_open(es, ee, node.start, node.end) as f64;
            let exon_len = len_half_open(es, ee) as f64;
            if exon_len > 0.0 {
                exoncov[e] += addcov * ovplen / exon_len;
            }
            e += 1;
        }
    }
}

/// store_transcript-like exon coverage accumulation from path-node usage.
fn accumulate_exon_cov_from_path_usage(
    use_path: &[usize],
    use_start: usize,
    use_last: usize,
    graph: &Graph,
    exons: &[(u64, u64)],
    nodecov_before: &[f64],
    nodeflux: &[f64],
    nodeflux_is_proportion: bool,
    node_rates: Option<&[f64]>,
) -> (Vec<f64>, f64) {
    if exons.is_empty() || use_last < use_start {
        return (Vec::new(), 0.0);
    }
    let tlen: f64 = exons
        .iter()
        .map(|(s, e)| len_half_open(*s, *e) as f64)
        .sum::<f64>()
        .max(1.0);
    let mut exoncov = vec![0.0f64; exons.len()];
    let mut cov_total = 0.0f64;
    let mut k = 0usize;
    for p in use_start..=use_last {
        let nid = use_path[p];
        if nid == graph.source_id || nid == graph.sink_id || nid >= graph.nodes.len() {
            continue;
        }
        if k >= nodeflux.len() {
            break;
        }
        let Some(node) = graph.nodes.get(nid) else {
            k += 1;
            continue;
        };
        let addcov = if nodeflux_is_proportion {
            nodeflux[k] * nodecov_before.get(nid).copied().unwrap_or(0.0)
        } else {
            // ecov = nodeflux[j] * noderate[path[j]]
            let rate = node_rates
                .and_then(|r| r.get(nid).copied())
                .filter(|v| *v > 0.0)
                .unwrap_or(1.0);
            nodeflux[k] * rate
        };
        for (ei, &(es, ee)) in exons.iter().enumerate() {
            if !overlaps_half_open(es, ee, node.start, node.end) {
                continue;
            }
            let ovplen = overlap_len_half_open(es, ee, node.start, node.end) as f64;
            if ovplen <= 0.0 {
                continue;
            }
            let elen = len_half_open(es, ee) as f64;
            if elen > 0.0 {
                exoncov[ei] += addcov * ovplen / elen;
                cov_total += addcov * ovplen / tlen;
            }
        }
        k += 1;
    }
    (exoncov, cov_total)
}

// =============================================================================
// Path extraction helpers (parse_trflong)
// =============================================================================

/// Check if two consecutive graph nodes have a splice between them (gap > 0bp).
/// is_splice_between: a->end + 1 < b->start (1-based); 0-based: a->end < b->start.
#[inline]
fn is_splice_between(graph: &Graph, a: usize, b: usize) -> bool {
    if let (Some(na), Some(nb)) = (graph.nodes.get(a), graph.nodes.get(b)) {
        na.end < nb.start
    } else {
        false
    }
}

/// Check if two consecutive graph nodes are contiguous (no gap, no splice).
/// 0-based half-open: a.end == b.start.
#[inline]
fn nodes_are_contiguous(graph: &Graph, a: usize, b: usize) -> bool {
    if let (Some(na), Some(nb)) = (graph.nodes.get(a), graph.nodes.get(b)) {
        contiguous_half_open(na.end, nb.start)
    } else {
        false
    }
}

#[allow(dead_code)]
fn intron_chain_from_nodes(graph: &Graph, nodes: &[usize]) -> Vec<(u64, u64)> {
    let mut out: Vec<(u64, u64)> = Vec::new();
    for w in nodes.windows(2) {
        let (a, b) = (w[0], w[1]);
        let Some(na) = graph.nodes.get(a) else { continue };
        let Some(nb) = graph.nodes.get(b) else { continue };
        if is_splice_between(graph, a, b) {
            // Use the same coordinate convention everywhere in checktrf: donor=node.end, acceptor=node.start.
            out.push((na.end, nb.start));
        }
    }
    out
}

#[allow(dead_code)]
fn intron_chains_equal_tol(a: &[(u64, u64)], b: &[(u64, u64)], tol: u64) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter()
        .zip(b.iter())
        .all(|((d1, a1), (d2, a2))| d1.abs_diff(*d2) <= tol && a1.abs_diff(*a2) <= tol)
}

fn select_contiguous_parent_for_longstart(
    graph: &Graph,
    first: usize,
    target_start: u64,
) -> Option<usize> {
    let node = graph.nodes.get(first)?;
    let mut covering: Option<usize> = None;
    let mut covering_start = 0u64;
    let mut fallback: Option<usize> = None;
    let mut fallback_start = 0u64;
    for p in node.parents.ones() {
        if p == graph.source_id || p == graph.sink_id || !nodes_are_contiguous(graph, p, first) {
            continue;
        }
        let Some(pnode) = graph.nodes.get(p) else {
            continue;
        };
        if pnode.start <= target_start {
            if covering.is_none() || pnode.start > covering_start {
                covering = Some(p);
                covering_start = pnode.start;
            }
        } else if fallback.is_none() || pnode.start < fallback_start {
            fallback = Some(p);
            fallback_start = pnode.start;
        }
    }
    covering.or(fallback)
}

fn select_contiguous_child_for_longend(
    graph: &Graph,
    last: usize,
    target_end: u64,
) -> Option<usize> {
    let node = graph.nodes.get(last)?;
    let mut covering: Option<usize> = None;
    let mut covering_end = u64::MAX;
    let mut fallback: Option<usize> = None;
    let mut fallback_end = 0u64;
    for c in node.children.ones() {
        if c == graph.source_id || c == graph.sink_id || !nodes_are_contiguous(graph, last, c) {
            continue;
        }
        let Some(cnode) = graph.nodes.get(c) else {
            continue;
        };
        if cnode.end >= target_end {
            if covering.is_none() || cnode.end < covering_end {
                covering = Some(c);
                covering_end = cnode.end;
            }
        } else if fallback.is_none() || cnode.end > fallback_end {
            fallback = Some(c);
            fallback_end = cnode.end;
        }
    }
    covering.or(fallback)
}

fn materialize_longread_seed_nodes(
    graph: &Graph,
    base_nodes: &[usize],
    longstart: u64,
    longend: u64,
) -> Vec<usize> {
    if base_nodes.is_empty() {
        return Vec::new();
    }
    let mut out = base_nodes.to_vec();
    if longstart > 0 {
        loop {
            let first = out[0];
            let Some(first_node) = graph.nodes.get(first) else {
                break;
            };
            if first_node.start <= longstart {
                break;
            }
            let Some(parent) = select_contiguous_parent_for_longstart(graph, first, longstart)
            else {
                break;
            };
            if out[0] == parent {
                break;
            }
            out.insert(0, parent);
        }
    }
    if longend > 0 {
        loop {
            let last = *out.last().unwrap_or(&out[0]);
            let Some(last_node) = graph.nodes.get(last) else {
                break;
            };
            if last_node.end >= longend {
                break;
            }
            let Some(child) = select_contiguous_child_for_longend(graph, last, longend) else {
                break;
            };
            if out.last().copied() == Some(child) {
                break;
            }
            out.push(child);
        }
    }
    out
}

/// Check if any long-read transfrag contains the splice edge (la->ra) as a consecutive node pair.
///
/// check if any transfrag contains BOTH splice edges in order.
/// This matches the original algorithm's `has_lr_witness_two_splices(la, ra, lb, rb, transfrag)`.
fn has_lr_witness_two_splices(
    transfrags: &[GraphTransfrag],
    la: usize,
    ra: usize,
    lb: usize,
    rb: usize,
) -> bool {
    for tf in transfrags {
        // checks all transfrags (commented out longread filter)
        let nodes = &tf.node_ids;
        let mut found_first = false;
        for w in nodes.windows(2) {
            if !found_first {
                if w[0] == la && w[1] == ra {
                    found_first = true;
                }
            } else if w[0] == lb && w[1] == rb {
                return true;
            }
        }
    }
    false
}

/// Build set of individually-witnessed splice junctions for fallback checks.
/// This is a weaker witness condition than `has_lr_witness_two_splices` and allows stitching
/// long paths from partial reads as long as each individual splice junction is directly supported.
fn build_lr_witness_splice_junctions(
    graph: &Graph,
    transfrags: &[GraphTransfrag],
) -> HashSet<(u64, u64)> {
    // Keyed by (donor, acceptor) coordinates rather than node IDs. This makes the witness check
    // robust to path extension selecting different but coordinate-equivalent boundary nodes.
    let mut out: HashSet<(u64, u64)> = Default::default();
    for tf in transfrags {
        if !tf.longread {
            continue;
        }
        for w in tf.node_ids.windows(2) {
            let (a, b) = (w[0], w[1]);
            if !is_splice_between(graph, a, b) {
                continue;
            }
            let Some(na) = graph.nodes.get(a) else {
                continue;
            };
            let Some(nb) = graph.nodes.get(b) else {
                continue;
            };
            out.insert((na.end, nb.start));
        }
    }
    out
}

#[inline]
#[allow(dead_code)]
#[allow(dead_code)]
fn lr_witnesses_splice_edge(
    graph: &Graph,
    la: usize,
    ra: usize,
    lr_witness: &HashSet<(u64, u64)>,
) -> bool {
    let (Some(na), Some(nb)) = (graph.nodes.get(la), graph.nodes.get(ra)) else {
        return false;
    };
    lr_witness.contains(&(na.end, nb.start))
}

/// Extend path leftward through contiguous-only nodes (no new splices).
/// Walks graph parent edges where parent.end == current_first.start.
/// This preserves the intron chain while extending the first exon boundary.
fn extend_path_left_contiguous(path: &mut Vec<usize>, graph: &Graph) {
    if path.is_empty() {
        return;
    }
    loop {
        let left = path[0];
        let node = match graph.nodes.get(left) {
            Some(n) => n,
            None => break,
        };
        // Find best contiguous parent (highest coverage)
        let mut best: Option<usize> = None;
        let mut best_cov: f64 = -1.0;
        for p in node.parents.ones() {
            if p == graph.source_id || p == graph.sink_id {
                continue;
            }
            if nodes_are_contiguous(graph, p, left) {
                let pcov = graph.nodes.get(p).map(|n| n.coverage).unwrap_or(0.0);
                if pcov > best_cov {
                    best_cov = pcov;
                    best = Some(p);
                }
            }
        }
        match best {
            Some(p) => path.insert(0, p),
            None => break,
        }
    }
}

/// Extend path rightward through contiguous-only nodes (no new splices).
/// Walks graph child edges where current_last.end == child.start.
/// This preserves the intron chain while extending the last exon boundary.
fn extend_path_right_contiguous(path: &mut Vec<usize>, graph: &Graph) {
    if path.is_empty() {
        return;
    }
    loop {
        let right = *path.last().unwrap();
        let node = match graph.nodes.get(right) {
            Some(n) => n,
            None => break,
        };
        // Find best contiguous child (highest coverage)
        let mut best: Option<usize> = None;
        let mut best_cov: f64 = -1.0;
        for c in node.children.ones() {
            if c == graph.source_id || c == graph.sink_id {
                continue;
            }
            if nodes_are_contiguous(graph, right, c) {
                let ccov = graph.nodes.get(c).map(|n| n.coverage).unwrap_or(0.0);
                if ccov > best_cov {
                    best_cov = ccov;
                    best = Some(c);
                }
            }
        }
        match best {
            Some(c) => path.push(c),
            None => break,
        }
    }
}

/// Extend path rightward through unambiguous splice junctions.
/// Only extends when current node has exactly one splice child (non-contiguous, non-sink).
/// Then continues with contiguous extension on the new exon.
/// This is conservative: no chimera risk since there's only one possible continuation.
fn extend_path_right_unambiguous_splice(
    path: &mut Vec<usize>,
    graph: &Graph,
    transfrags: &[GraphTransfrag],
) {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    if path.is_empty() {
        return;
    }
    let mut visited_splices: HashSet<usize> = Default::default();
    let mut splice_hops = 0usize;
    const MAX_SPLICE_HOPS: usize = 1;
    const MIN_SPLICE_SUPPORT_COUNT: usize = 3;
    const MIN_SPLICE_SUPPORT_ABUND: f64 = 20.0;
    loop {
        if splice_hops >= MAX_SPLICE_HOPS {
            break;
        }
        // First extend contiguously as far as possible
        extend_path_right_contiguous(path, graph);
        let right = *path.last().unwrap();
        if visited_splices.contains(&right) {
            break;
        }
        visited_splices.insert(right);

        let node = match graph.nodes.get(right) {
            Some(n) => n,
            None => break,
        };

        // Collect splice children (non-contiguous, non-source/sink)
        let splice_children: Vec<usize> = node
            .children
            .ones()
            .filter(|&c| c != source_id && c != sink_id && !nodes_are_contiguous(graph, right, c))
            .collect();

        // Only extend if exactly one splice child
        if splice_children.len() != 1 {
            break;
        }
        let child = splice_children[0];

        // Require transfrag evidence: at least one transfrag spans this splice junction
        let mut support_count = 0usize;
        let mut support_abund = 0.0f64;
        let mut guide_support = 0usize;
        for tf in transfrags {
            if tf.abundance <= 0.0 {
                continue;
            }
            let ri = tf.node_ids.iter().position(|&x| x == right);
            let ci = tf.node_ids.iter().position(|&x| x == child);
            if matches!((ri, ci), (Some(r), Some(c)) if c == r + 1) {
                support_count += 1;
                support_abund += tf.abundance;
                if tf.guide {
                    guide_support += 1;
                }
            }
        }
        if guide_support == 0
            && (support_count < MIN_SPLICE_SUPPORT_COUNT
                || support_abund < MIN_SPLICE_SUPPORT_ABUND)
        {
            break;
        }

        // Skip if child is already on path
        if path.contains(&child) {
            break;
        }

        path.push(child);
        splice_hops += 1;
        // Continue loop to extend contiguously from the new node
    }
}

/// Extend path leftward through unambiguous splice junctions.
/// Only extends when current node has exactly one splice parent (non-contiguous, non-source).
fn extend_path_left_unambiguous_splice(
    path: &mut Vec<usize>,
    graph: &Graph,
    transfrags: &[GraphTransfrag],
) {
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    if path.is_empty() {
        return;
    }
    let mut visited_splices: HashSet<usize> = Default::default();
    let mut splice_hops = 0usize;
    const MAX_SPLICE_HOPS: usize = 1;
    const MIN_SPLICE_SUPPORT_COUNT: usize = 3;
    const MIN_SPLICE_SUPPORT_ABUND: f64 = 20.0;
    loop {
        if splice_hops >= MAX_SPLICE_HOPS {
            break;
        }
        // First extend contiguously as far as possible
        extend_path_left_contiguous(path, graph);
        let left = path[0];
        if visited_splices.contains(&left) {
            break;
        }
        visited_splices.insert(left);

        let node = match graph.nodes.get(left) {
            Some(n) => n,
            None => break,
        };

        // Collect splice parents (non-contiguous, non-source/sink)
        let splice_parents: Vec<usize> = node
            .parents
            .ones()
            .filter(|&p| p != source_id && p != sink_id && !nodes_are_contiguous(graph, p, left))
            .collect();

        // Only extend if exactly one splice parent
        if splice_parents.len() != 1 {
            break;
        }
        let parent = splice_parents[0];

        // Require transfrag evidence: at least one transfrag spans this splice junction
        let mut support_count = 0usize;
        let mut support_abund = 0.0f64;
        let mut guide_support = 0usize;
        for tf in transfrags {
            if tf.abundance <= 0.0 {
                continue;
            }
            let pi = tf.node_ids.iter().position(|&x| x == parent);
            let li = tf.node_ids.iter().position(|&x| x == left);
            if matches!((pi, li), (Some(p), Some(l)) if l == p + 1) {
                support_count += 1;
                support_abund += tf.abundance;
                if tf.guide {
                    guide_support += 1;
                }
            }
        }
        if guide_support == 0
            && (support_count < MIN_SPLICE_SUPPORT_COUNT
                || support_abund < MIN_SPLICE_SUPPORT_ABUND)
        {
            break;
        }

        // Skip if parent is already on path
        if path.contains(&parent) {
            break;
        }

        path.insert(0, parent);
        splice_hops += 1;
        // Continue loop to extend contiguously from the new node
    }
}

/// Port of best_trf_match for checktrf rescue.
/// Returns indices into keep_paths plus abundance sum of selected ties.
fn best_trf_match(
    tf_nodes: &[usize],
    keep_paths: &[(Vec<usize>, f64, bool, usize)],
    graph: &Graph,
    guide_only: bool,
    trace_tf: Option<usize>,
) -> (Vec<usize>, f64) {
    if tf_nodes.is_empty() || keep_paths.is_empty() {
        return (Vec::new(), 0.0);
    }
    let trace_btm =
        trace_tf.is_some() && tf_nodes.iter().any(|&nid| trace_locus_active(graph, nid));
    let tf_id_nonmono = tf_nodes.windows(2).any(|w| w[1] < w[0]);
    if trace_btm {
        eprintln!(
            "[TRACE_BTM] t={} guide_only={} tf_nodes={:?} id_nonmono={}",
            trace_tf.unwrap_or(usize::MAX),
            guide_only,
            tf_nodes,
            tf_id_nonmono
        );
    }
    let mut tmatch: Vec<usize> = Vec::new();
    let mut abundancesum = 0.0f64;
    let mut maxintersect: i64 = 0;
    let mut mininternaldist: i64 = i64::MAX / 4;
    let mut minflank: i64 = i64::MAX / 4;

    let nlen = |nid: usize| -> i64 { graph.nodes.get(nid).map(|n| n.length() as i64).unwrap_or(0) };
    let is_intron = |a: usize, b: usize| -> bool {
        match (graph.nodes.get(a), graph.nodes.get(b)) {
            (Some(na), Some(nb)) => nb.start > na.end,
            _ => false,
        }
    };

    for (kidx, (k_nodes, k_cov, k_guide, _out_idx)) in keep_paths.iter().enumerate() {
        let k_id_nonmono = k_nodes.windows(2).any(|w| w[1] < w[0]);
        if trace_btm {
            eprintln!(
                "[TRACE_BTM] t={} k={} keep_nodes={:?} guide={} keep_cov={:.4} id_nonmono={}",
                trace_tf.unwrap_or(usize::MAX),
                kidx,
                k_nodes,
                *k_guide,
                *k_cov,
                k_id_nonmono
            );
        }
        if *k_guide != guide_only {
            if trace_btm {
                eprintln!(
                    "[TRACE_BTM] t={} k={} skip=guide_mismatch",
                    trace_tf.unwrap_or(usize::MAX),
                    kidx
                );
            }
            continue;
        }
        if k_nodes.is_empty() {
            if trace_btm {
                eprintln!(
                    "[TRACE_BTM] t={} k={} skip=empty_keep",
                    trace_tf.unwrap_or(usize::MAX),
                    kidx
                );
            }
            continue;
        }
        if !(tf_nodes[0] <= *k_nodes.last().unwrap() && k_nodes[0] <= *tf_nodes.last().unwrap()) {
            if trace_btm {
                eprintln!(
                    "[TRACE_BTM] t={} k={} skip=coord_overlap raw_tf_first={} raw_tf_last={} raw_keep_first={} raw_keep_last={}",
                    trace_tf.unwrap_or(usize::MAX),
                    kidx,
                    tf_nodes[0],
                    *tf_nodes.last().unwrap(),
                    k_nodes[0],
                    *k_nodes.last().unwrap()
                );
            }
            continue;
        }

        let mut i = 0usize;
        let mut j = 0usize;
        while j < k_nodes.len() && tf_nodes[i] > k_nodes[j] {
            j += 1;
        }

        let mut leftdist: i64 = 0;
        let mut internaldist: i64 = 0;
        while i < tf_nodes.len() && j < k_nodes.len() && tf_nodes[i] < k_nodes[j] {
            leftdist += nlen(tf_nodes[i]);
            i += 1;
            if i < tf_nodes.len()
                && match (
                    graph.nodes.get(tf_nodes[i]),
                    graph.nodes.get(tf_nodes[i - 1]),
                ) {
                    (Some(curr), Some(prev)) => curr.end > prev.start.saturating_add(1),
                    _ => false,
                }
            {
                internaldist = leftdist;
                leftdist = 0;
                break;
            }
        }
        if internaldist > CHI_THR_BP {
            continue;
        }

        let mut intersect: i64 = 0;
        let mut intron = true;
        let mut rightdist: i64 = 0;

        while i < tf_nodes.len() {
            if j == k_nodes.len() {
                rightdist += nlen(tf_nodes[i]);
                if intron || (i > 0 && is_intron(tf_nodes[i - 1], tf_nodes[i])) {
                    internaldist += rightdist;
                    rightdist = 0;
                    intron = true;
                }
            } else if tf_nodes[i] == k_nodes[j] {
                intersect += nlen(tf_nodes[i]);
                internaldist += rightdist;
                rightdist = 0;
                intron = false;
                j += 1;
            } else {
                if !intron {
                    rightdist += nlen(tf_nodes[i]);
                    if i > 0 && is_intron(tf_nodes[i - 1], tf_nodes[i]) {
                        internaldist += rightdist;
                        rightdist = 0;
                        intron = true;
                    }
                } else {
                    internaldist += nlen(tf_nodes[i]);
                }
                if tf_nodes[i] > k_nodes[j] {
                    j += 1;
                }
            }
            i += 1;
        }

        if intersect <= 0 || internaldist > CHI_THR_BP {
            if trace_btm {
                eprintln!(
                    "[TRACE_BTM] t={} k={} reject intersect={} internaldist={} flank={}",
                    trace_tf.unwrap_or(usize::MAX),
                    kidx,
                    intersect,
                    internaldist,
                    leftdist + rightdist
                );
            }
            continue;
        }

        let flank = leftdist + rightdist;
        if trace_btm {
            eprintln!(
                "[TRACE_BTM] t={} k={} score intersect={} internaldist={} flank={} maxintersect={} mininternaldist={} minflank={}",
                trace_tf.unwrap_or(usize::MAX),
                kidx,
                intersect,
                internaldist,
                flank,
                maxintersect,
                mininternaldist,
                minflank
            );
        }
        if intersect > maxintersect {
            tmatch.clear();
            tmatch.push(kidx);
            abundancesum = *k_cov;
            mininternaldist = internaldist;
            minflank = flank;
            maxintersect = intersect;
        } else if intersect == maxintersect {
            if internaldist < mininternaldist {
                tmatch.clear();
                tmatch.push(kidx);
                mininternaldist = internaldist;
                abundancesum = *k_cov;
                minflank = flank;
            } else if flank == minflank {
                tmatch.push(kidx);
                abundancesum += *k_cov;
            } else if flank < minflank {
                tmatch.clear();
                tmatch.push(kidx);
                abundancesum = *k_cov;
                minflank = flank;
            }
        }
        if trace_btm {
            eprintln!(
                "[TRACE_BTM] t={} after k={} tmatch={:?} abundancesum={:.4}",
                trace_tf.unwrap_or(usize::MAX),
                kidx,
                tmatch,
                abundancesum
            );
        }
    }

    (tmatch, abundancesum)
}

/// Redistribute transfrag abundance into matched kept predictions, mirroring checktrf/second-pass update:
/// add per-node support to overlapping exons only when the node is present in kept path.
/// Rust stores per-base coverage, so bp-weighted additions are normalized by exon/tx lengths.
fn redistribute_transfrag_to_matches(
    tf: &GraphTransfrag,
    tmatch: &[usize],
    abundancesum: f64,
    kept_paths: &[(Vec<usize>, f64, bool, usize)],
    graph: &Graph,
    out: &mut [Transcript],
) -> bool {
    if tmatch.is_empty() || tf.node_ids.is_empty() || out.is_empty() {
        return false;
    }
    let mut any_updated = false;
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let first_real_idx = tf
        .node_ids
        .iter()
        .position(|&n| n != source_id && n != sink_id);
    let last_real_idx = tf
        .node_ids
        .iter()
        .rposition(|&n| n != source_id && n != sink_id);

    for (j, kmatch) in tmatch.iter().enumerate() {
        let Some((keep_nodes, keep_support, _keep_guide, out_idx)) = kept_paths.get(*kmatch) else {
            continue;
        };
        let oi = *out_idx;
        if oi >= out.len() {
            continue;
        }
        let abundprop = if abundancesum > 0.0 {
            tf.abundance * *keep_support / abundancesum
        } else if j == 0 {
            tf.abundance
        } else {
            0.0
        };
        if abundprop <= 0.0 {
            continue;
        }
        let tx_len: f64 = out[oi]
            .exons
            .iter()
            .map(|(s, e)| len_half_open(*s, *e).max(1) as f64)
            .sum::<f64>()
            .max(1.0);
        let mut p = 0usize;
        let mut i = 0usize;
        while i < tf.node_ids.len() && p < out[oi].exons.len() {
            let nid = tf.node_ids[i];
            if nid == source_id || nid == sink_id {
                i += 1;
                continue;
            }
            let Some(node) = graph.nodes.get(nid) else {
                i += 1;
                continue;
            };
            let (es, ee) = out[oi].exons[p];
            if node.end < es {
                i += 1;
                continue;
            }
            if ee < node.start {
                p += 1;
                continue;
            }
            if keep_nodes.contains(&nid) {
                let mut nodelen = node.length() as f64;
                if first_real_idx == Some(i) && tf.longstart > node.start {
                    nodelen -= tf.longstart.saturating_sub(node.start) as f64;
                }
                if last_real_idx == Some(i) && tf.longend > 0 && tf.longend < node.end {
                    nodelen -= node.end.saturating_sub(tf.longend) as f64;
                }
                if nodelen < 0.0 {
                    nodelen += node.length() as f64;
                }
                if nodelen > 0.0 {
                    let addcov_bp = nodelen * abundprop;
                    let exon_len = len_half_open(es, ee).max(1) as f64;
                    if p < out[oi].exon_cov.len() {
                        out[oi].exon_cov[p] += addcov_bp / exon_len;
                    }
                    out[oi].coverage += addcov_bp / tx_len;
                    any_updated = true;
                }
            }
            i += 1;
        }
    }
    any_updated
}

#[inline]
fn edge_set(pathpat: &mut GBitVec, graph: &Graph, a: usize, b: usize, val: bool) {
    if let Some(eid) = graph.edge_bit_index(a, b) {
        if val {
            pathpat.set_bit(eid);
        } else {
            pathpat.clear_bit(eid);
        }
    }
}

#[inline]
fn edge_bit(pathpat: &GBitVec, graph: &Graph, a: usize, b: usize) -> bool {
    graph
        .edge_bit_index(a, b)
        .map(|eid| pathpat.get_bit(eid))
        .unwrap_or(false)
}

fn trace_outgoing_edge_state(
    label: &str,
    seed_idx: usize,
    pivot: usize,
    pathpat: &GBitVec,
    graph: &Graph,
) {
    if !trace_locus_active(graph, pivot) {
        return;
    }
    let Some(node) = graph.nodes.get(pivot) else {
        return;
    };
    eprintln!(
        "[TRACE_PAT] stage={} idx={} pivot={} coord={}-{} node_bit={} pop={} nchildren={}",
        label,
        seed_idx,
        pivot,
        node.start,
        node.end,
        pathpat.get_bit(pivot),
        pathpat.len_bits(),
        node.children.count_ones()
    );
    for child in node.children.ones() {
        let ccoord = graph
            .nodes
            .get(child)
            .map(|n| (n.start, n.end))
            .unwrap_or((0, 0));
        eprintln!(
            "[TRACE_PAT]   child={} coord={}-{} edge_bit={} child_bit={} sink={}",
            child,
            ccoord.0,
            ccoord.1,
            edge_bit(pathpat, graph, pivot, child),
            pathpat.get_bit(child),
            child == graph.sink_id
        );
    }
}

fn pathpat_trace_active() -> bool {
    std::env::var_os("RUSTLE_TRACE_PATHPAT").is_some()
        || std::env::var_os("RUSTLE_TRACE_LOG_STYLE").is_some()
}

/// Summarize pathpat as (node_bits, edge_bits, set_bits_list) for trace output.
fn pathpat_summary(pat: &GBitVec, gno: usize) -> (usize, usize, String) {
    let mut node_bits = 0usize;
    let mut edge_bits = 0usize;
    let mut set_bits: Vec<usize> = Vec::new();
    // Pattern size = gno + gno*(gno-1)/2 edges typically, but just scan up to len_bits popcount
    // by iterating through potential positions.
    let max_pos = gno * gno; // upper bound on pattern size
    for i in 0..max_pos {
        if pat.get_bit(i) {
            if i < gno {
                node_bits += 1;
            } else {
                edge_bits += 1;
            }
            set_bits.push(i);
        }
    }
    let summary = if set_bits.len() <= 20 {
        format!("{:?}", set_bits)
    } else {
        format!("{:?}...+{}", &set_bits[..20], set_bits.len() - 20)
    };
    (node_bits, edge_bits, summary)
}

#[inline]
fn rebuild_flow_pathpat(path: &[usize], seed_pattern: &GBitVec, graph: &Graph) -> GBitVec {
    let mut pathpat = GBitVec::new(graph.pattern_size());
    for &nid in path {
        pathpat.set_bit(nid);
    }
    for w in path.windows(2) {
        edge_set(&mut pathpat, graph, w[0], w[1], true);
    }
    pathpat.or_assign(seed_pattern);
    pathpat
}

fn tf_weak(tf: &GraphTransfrag, graph: &Graph, nodecov: &[f64]) -> bool {
    // reuses `transfrag->weak` both for absorbed/grouped members set during
    // process_transfrags and for lazily computed coverage-weakness in
    // replace_transfrag/compute_weak. Rust stores the grouping marker eagerly in
    // `tf.weak`, so helper selection must honor it before the coverage-based test.
    if tf.weak != 0 {
        return true;
    }
    if tf.coverage_weak {
        return true;
    }
    let thr = DROP_FACTOR + ERROR_PERC;
    for w in tf.node_ids.windows(2) {
        let a = w[0];
        let b = w[1];
        if b != a + 1 {
            continue;
        }
        let (Some(na), Some(nb)) = (graph.nodes.get(a), graph.nodes.get(b)) else {
            continue;
        };
        if na.end != nb.start {
            continue;
        }
        let ca = nodecov.get(a).copied().unwrap_or(0.0);
        let cb = nodecov.get(b).copied().unwrap_or(0.0);
        if ca * thr > cb || cb * thr > ca {
            return true;
        }
    }
    false
}

/// Cached version of tf_weak: computes on first call, returns cached value thereafter.
/// Matches lazy-cache behavior where `transfrag[t]->weak` is computed once and reused.
fn tf_weak_cached(
    tidx: usize,
    transfrags: &[GraphTransfrag],
    graph: &Graph,
    nodecov: &[f64],
    cache: &mut Vec<Option<bool>>,
) -> bool {
    if let Some(cached) = cache.get(tidx).copied().flatten() {
        return cached;
    }
    let result = tf_weak(&transfrags[tidx], graph, nodecov);
    if tidx < cache.len() {
        cache[tidx] = Some(result);
    }
    result
}

fn replace_tmax(
    cand: usize,
    tmax: &mut Option<usize>,
    transfrags: &[GraphTransfrag],
    graph: &Graph,
    nodecov: &[f64],
) {
    match *tmax {
        None => *tmax = Some(cand),
        Some(cur) => {
            let cweak = tf_weak(&transfrags[cand], graph, nodecov);
            let curweak = tf_weak(&transfrags[cur], graph, nodecov);
            if (!curweak && !cweak && transfrags[cand].abundance > transfrags[cur].abundance)
                || (curweak && (!cweak || transfrags[cand].abundance > transfrags[cur].abundance))
            {
                *tmax = Some(cand);
            }
        }
    }
}

/// Find the next set bit in pathpattern that has a genomic coordinate greater than the current node.
/// This is needed because node IDs may not be in genomic order after graph remapping.
fn next_set_bit_by_coord(
    pathpat: &GBitVec,
    from_node: usize,
    max_node: usize,
    graph: &Graph,
) -> Option<usize> {
    let from_start = graph.nodes.get(from_node).map(|n| n.start).unwrap_or(0);
    let max_start = graph.nodes.get(max_node).map(|n| n.start).unwrap_or(u64::MAX);
    
    // Search all nodes in the graph for ones that:
    // 1. Have their bit set in pathpat
    // 2. Have start coordinate > from_start
    // 3. Have start coordinate <= max_start
    // Return the one with the smallest start coordinate
    graph.nodes
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != graph.source_id && *i != graph.sink_id)
        .filter(|(i, _)| pathpat.get_bit(*i))
        .filter(|(_, n)| n.start > from_start && n.start <= max_start)
        .min_by_key(|(_, n)| n.start)
        .map(|(i, _)| i)
}



#[inline]
fn first_last_raw(tf: &GraphTransfrag) -> Option<(usize, usize)> {
    Some((*tf.node_ids.first()?, *tf.node_ids.last()?))
}

fn format_path_nodes(path: &[usize], graph: &Graph, limit: usize) -> String {
    path.iter()
        .take(limit)
        .map(|&nid| {
            let coord = graph
                .nodes
                .get(nid)
                .map(|n| format!("{}-{}", n.start, n.end))
                .unwrap_or_else(|| "0-0".to_string());
            format!("{nid}({coord})")
        })
        .collect::<Vec<_>>()
        .join(" ")
}

#[derive(Default, Debug, Clone)]
struct LongRecDiag {
    back_unreachable_minpath: usize,
    back_no_reach: usize,
    back_no_choice: usize,
    back_exclude_no_support: usize,
    fwd_unreachable_maxpath: usize,
    fwd_no_reach: usize,
    fwd_no_choice: usize,
    fwd_exclude_no_support: usize,
    // Detailed diagnostics
    back_parents_checked: usize,
    back_parents_skipped_reachability: usize,
    _back_parents_skipped_edge: usize,
    _back_parents_skipped_coverage: usize,
    _back_parents_skipped_weak: usize,
    back_transfrags_checked: usize,
    back_transfrags_skipped_depleted: usize,
    back_transfrags_skipped_notlong: usize,
    back_transfrags_skipped_ends_at_sink: usize,
    _back_transfrags_skipped_pattern: usize,
    back_transfrags_skipped_onpath: usize,
}

fn onpath_long(
    trpattern: &GBitVec,
    trnode: &[usize],
    pathpattern: &GBitVec,
    minp: usize,
    maxp: usize,
    graph: &Graph,
) -> bool {
    if trnode.is_empty() {
        return false;
    }
    // compare based on genomic coordinates, not node IDs
    let minp_start = graph.nodes.get(minp).map(|n| n.start).unwrap_or(0);
    let maxp_start = graph.nodes.get(maxp).map(|n| n.start).unwrap_or(0);
    if minp_start > maxp_start {
        return false;
    }
    let source = graph.source_id;
    let sink = graph.sink_id;
    let mut j = 0usize;

    if edge_bit(pathpattern, graph, source, minp) {
        if trnode.len() >= 2 && trnode[0] == source {
            if trnode[1] != minp {
                return false;
            }
        } else {
            // compare based on genomic coordinates, not node IDs
            let trnode_0_start = graph.nodes.get(trnode[0]).map(|n| n.start).unwrap_or(0);
            let minp_start = graph.nodes.get(minp).map(|n| n.start).unwrap_or(0);
            if trnode_0_start < minp_start && !node_can_reach(graph, trnode[0], minp) {
                return false;
            }
        }
    }
    if edge_bit(pathpattern, graph, maxp, sink) {
        if trnode.len() >= 2 && *trnode.last().unwrap() == sink {
            if trnode[trnode.len() - 2] != maxp {
                return false;
            }
        } else {
            // compare based on genomic coordinates, not node IDs
            let trnode_last = *trnode.last().unwrap();
            let trnode_last_start = graph.nodes.get(trnode_last).map(|n| n.start).unwrap_or(0);
            let maxp_start = graph.nodes.get(maxp).map(|n| n.start).unwrap_or(0);
            if trnode_last_start > maxp_start && !node_can_reach(graph, maxp, trnode_last) {
                return false;
            }
        }
    }

    let mut prevp: Option<usize> = None;
    let mut p = minp;
    // Precompute node starts for coordinate-based comparisons
    let trnode_starts: Vec<u64> = trnode.iter().map(|n| graph.nodes.get(*n).map(|node| node.start).unwrap_or(0)).collect();
    let p_start = |node_id: usize| graph.nodes.get(node_id).map(|n| n.start).unwrap_or(0);
    loop {
        while j < trnode.len() && trnode_starts[j] < p_start(p) {
            if prevp == Some(source) || p == sink {
                return false;
            }
            if !graph
                .nodes
                .get(trnode[j])
                .and_then(|n| n.childpat.as_ref())
                .map(|pat| pat.contains(p))
                .unwrap_or(false)
            {
                return false;
            }
            j += 1;
        }
        if j == trnode.len() {
            return true;
        }
        if trnode_starts[j] > p_start(p) {
            if !graph
                .nodes
                .get(p)
                .and_then(|n| n.childpat.as_ref())
                .map(|pat| pat.contains(trnode[j]))
                .unwrap_or(false)
            {
                return false;
            }
            // Allow transfrag edges that aren't on path if:
            // 1. The path has a gap between these nodes (childpat reachability), OR
            // 2. The transfrag edge fills a gap in the path
            // This is needed for 5' terminal exon recovery where path has gaps
            if j > 0 && edge_bit(trpattern, graph, trnode[j - 1], trnode[j]) &&
                !edge_bit(pathpattern, graph, trnode[j - 1], trnode[j]) {
                // Check if path has reachability between these nodes via childpat
                let path_can_reach = graph
                    .nodes
                    .get(trnode[j - 1])
                    .and_then(|n| n.childpat.as_ref())
                    .map(|pat| pat.contains(trnode[j]))
                    .unwrap_or(false);
                if !path_can_reach {
                    return false;
                }
            }
        } else {
            j += 1;
        }

        if p == maxp {
            return true;
        }
        prevp = Some(p);
        let Some(np) = next_set_bit_by_coord(pathpattern, p, maxp, graph) else {
            return false;
        };
        p = np;
    }
}

fn onpath_long_reason(
    trpattern: &GBitVec,
    trnode: &[usize],
    pathpattern: &GBitVec,
    minp: usize,
    maxp: usize,
    graph: &Graph,
) -> &'static str {
    static ONPATH_DBG_RIGHT_GAP: std::sync::atomic::AtomicUsize =
        std::sync::atomic::AtomicUsize::new(0);
    static ONPATH_DBG_LEFT_GAP: std::sync::atomic::AtomicUsize =
        std::sync::atomic::AtomicUsize::new(0);
    static ONPATH_DBG_SOURCE_EDGE: std::sync::atomic::AtomicUsize =
        std::sync::atomic::AtomicUsize::new(0);
    static ONPATH_DBG_SINK_EDGE: std::sync::atomic::AtomicUsize =
        std::sync::atomic::AtomicUsize::new(0);
    const ONPATH_DBG_LIMIT: usize = 20;
    let trace_gate = std::env::var_os("RUSTLE_TRACE_LOCUS").is_some();

    if trnode.is_empty() {
        return "fail:empty_trnode";
    }
    // compare based on genomic coordinates, not node IDs
    let minp_start = node_start_or_zero(graph, minp);
    let maxp_start = node_start_or_zero(graph, maxp);
    if minp_start > maxp_start {
        return "fail:minp_gt_maxp";
    }

    let source = graph.source_id;
    let sink = graph.sink_id;
    let mut j = 0usize;

    if edge_bit(pathpattern, graph, source, minp) {
        if trnode.len() >= 2 && trnode[0] == source {
            if trnode[1] != minp {
                return "fail:source_edge_wrong_second";
            }
        } else {
            // compare based on genomic coordinates, not node IDs
            let trnode_0_start = node_start_or_zero(graph, trnode[0]);
            let minp_start = node_start_or_zero(graph, minp);
            if trnode_0_start < minp_start && !node_can_reach(graph, trnode[0], minp) {
                if std::env::var_os("RUSTLE_TRACE_ONPATH_DEBUG").is_some() {
                    let n = ONPATH_DBG_SOURCE_EDGE
                        .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if n < ONPATH_DBG_LIMIT {
                        if !trace_gate
                            || trace_locus_active(graph, trnode[0])
                            || trace_locus_active(graph, minp)
                        {
                        eprintln!(
                            "[ONPATH_DEBUG] reason=source_edge_starts_before_minp tr0={}({}-{}) minp={}({}-{}) can_reach={} path_has_source_edge={}",
                            trnode[0],
                            node_start_or_zero(graph, trnode[0]),
                            node_end_or_zero(graph, trnode[0]),
                            minp,
                            node_start_or_zero(graph, minp),
                            node_end_or_zero(graph, minp),
                            node_can_reach(graph, trnode[0], minp) as u8,
                            edge_bit(pathpattern, graph, source, minp) as u8
                        );
                        }
                    }
                }
                return "fail:source_edge_starts_before_minp";
            }
        }
    }
    if edge_bit(pathpattern, graph, maxp, sink) {
        if trnode.len() >= 2 && *trnode.last().unwrap() == sink {
            if trnode[trnode.len() - 2] != maxp {
                return "fail:sink_edge_wrong_penultimate";
            }
        } else {
            // compare based on genomic coordinates, not node IDs
            let trnode_last = *trnode.last().unwrap();
            let trnode_last_start = node_start_or_zero(graph, trnode_last);
            let maxp_start = node_start_or_zero(graph, maxp);
            if trnode_last_start > maxp_start && !node_can_reach(graph, maxp, trnode_last) {
                if std::env::var_os("RUSTLE_TRACE_ONPATH_DEBUG").is_some() {
                    let n = ONPATH_DBG_SINK_EDGE.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if n < ONPATH_DBG_LIMIT {
                        if !trace_gate
                            || trace_locus_active(graph, maxp)
                            || trace_locus_active(graph, trnode_last)
                        {
                        eprintln!(
                            "[ONPATH_DEBUG] reason=sink_edge_ends_after_maxp maxp={}({}-{}) tr_last={}({}-{}) can_reach={} path_has_sink_edge={}",
                            maxp,
                            node_start_or_zero(graph, maxp),
                            node_end_or_zero(graph, maxp),
                            trnode_last,
                            node_start_or_zero(graph, trnode_last),
                            node_end_or_zero(graph, trnode_last),
                            node_can_reach(graph, maxp, trnode_last) as u8,
                            edge_bit(pathpattern, graph, maxp, sink) as u8
                        );
                        }
                    }
                }
                return "fail:sink_edge_ends_after_maxp";
            }
        }
    }

    let mut prevp: Option<usize> = None;
    let mut p = minp;
    // Precompute node starts for coordinate-based comparisons
    let trnode_starts: Vec<u64> = trnode.iter().map(|&n| node_start_or_zero(graph, n)).collect();
    let p_start = |node_id: usize| node_start_or_zero(graph, node_id);
    loop {
        while j < trnode.len() && trnode_starts[j] < p_start(p) {
            if prevp == Some(source) {
                return "fail:left_gap_after_source";
            }
            if p == sink {
                return "fail:left_gap_hits_sink";
            }
            if !graph
                .nodes
                .get(trnode[j])
                .and_then(|n| n.childpat.as_ref())
                .map(|pat| pat.contains(p))
                .unwrap_or(false)
            {
                if std::env::var_os("RUSTLE_TRACE_ONPATH_DEBUG").is_some() {
                    let n = ONPATH_DBG_LEFT_GAP.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if n < ONPATH_DBG_LIMIT {
                        if !trace_gate
                            || trace_locus_active(graph, p)
                            || trace_locus_active(graph, trnode[j])
                        {
                        eprintln!(
                            "[ONPATH_DEBUG] reason=left_gap_unreachable p={}({}-{}) trj={}({}-{}) j={} tr_on_path={} tr_edge_prev={} prevp={:?}",
                            p,
                            node_start_or_zero(graph, p),
                            node_end_or_zero(graph, p),
                            trnode[j],
                            node_start_or_zero(graph, trnode[j]),
                            node_end_or_zero(graph, trnode[j]),
                            j,
                            pathpattern.get_bit(trnode[j]) as u8,
                            if j > 0 {
                                edge_bit(trpattern, graph, trnode[j - 1], trnode[j]) as u8
                            } else {
                                0
                            },
                            prevp
                        );
                        }
                    }
                }
                return "fail:left_gap_unreachable";
            }
            j += 1;
        }
        if j == trnode.len() {
            return "ok";
        }
        if trnode_starts[j] > p_start(p) {
            if !graph
                .nodes
                .get(p)
                .and_then(|n| n.childpat.as_ref())
                .map(|pat| pat.contains(trnode[j]))
                .unwrap_or(false)
            {
                if std::env::var_os("RUSTLE_TRACE_ONPATH_DEBUG").is_some() {
                    let n = ONPATH_DBG_RIGHT_GAP
                        .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    if n < ONPATH_DBG_LIMIT {
                        if !trace_gate
                            || trace_locus_active(graph, p)
                            || trace_locus_active(graph, trnode[j])
                        {
                        eprintln!(
                            "[ONPATH_DEBUG] reason=right_gap_unreachable p={}({}-{}) trj={}({}-{}) j={} tr_on_path={} p_childpat_present={} p_can_reach_trj={} tr_edge_prev={} prevp={:?}",
                            p,
                            node_start_or_zero(graph, p),
                            node_end_or_zero(graph, p),
                            trnode[j],
                            node_start_or_zero(graph, trnode[j]),
                            node_end_or_zero(graph, trnode[j]),
                            j,
                            pathpattern.get_bit(trnode[j]) as u8,
                            graph.nodes.get(p).and_then(|n| n.childpat.as_ref()).is_some() as u8,
                            node_can_reach(graph, p, trnode[j]) as u8,
                            if j > 0 {
                                edge_bit(trpattern, graph, trnode[j - 1], trnode[j]) as u8
                            } else {
                                0
                            },
                            prevp
                        );
                        }
                    }
                }
                return "fail:right_gap_unreachable";
            }
            // Allow transfrag edges that aren't on path if path has reachability
            // This is needed for 5' terminal exon recovery where path has gaps
            if j > 0 && edge_bit(trpattern, graph, trnode[j - 1], trnode[j]) &&
                !edge_bit(pathpattern, graph, trnode[j - 1], trnode[j]) {
                let path_can_reach = graph
                    .nodes
                    .get(trnode[j - 1])
                    .and_then(|n| n.childpat.as_ref())
                    .map(|pat| pat.contains(trnode[j]))
                    .unwrap_or(false);
                if !path_can_reach {
                    return "fail:tr_edge_crosses_gap";
                }
            }
        } else {
            j += 1;
        }

        if p == maxp {
            return "ok";
        }
        prevp = Some(p);
        let Some(np) = next_set_bit_by_coord(pathpattern, p, maxp, graph) else {
            return "fail:no_next_path_node";
        };
        p = np;
    }
}

fn onpath_short(
    trpattern: &GBitVec,
    trnode: &[usize],
    pathpattern: &GBitVec,
    minp: usize,
    maxp: usize,
    graph: &Graph,
) -> bool {
    if trnode.is_empty() {
        return false;
    }
    // compare based on genomic coordinates, not node IDs
    let minp_start = graph.nodes.get(minp).map(|n| n.start).unwrap_or(0);
    let maxp_start = graph.nodes.get(maxp).map(|n| n.start).unwrap_or(0);
    if minp_start > maxp_start {
        return false;
    }
    let mut j = 0usize;
    let mut prevp: Option<usize> = None;
    let mut p = minp;
    // Precompute node starts for coordinate-based comparisons
    let trnode_starts: Vec<u64> = trnode.iter().map(|n| graph.nodes.get(*n).map(|node| node.start).unwrap_or(0)).collect();
    let p_start = |node_id: usize| graph.nodes.get(node_id).map(|n| n.start).unwrap_or(0);
    loop {
        while j < trnode.len() && trnode_starts[j] < p_start(p) {
            let edge_on_path = prevp
                .and_then(|pp| graph.edge_bit_index(pp, p))
                .map(|eid| pathpattern.get_bit(eid))
                .unwrap_or(false);
            if edge_on_path {
                return false;
            }
            if !graph
                .nodes
                .get(trnode[j])
                .and_then(|n| n.childpat.as_ref())
                .map(|pat| pat.contains(p))
                .unwrap_or(false)
            {
                return false;
            }
            j += 1;
        }
        if j == trnode.len() {
            return true;
        }
        if trnode_starts[j] > p_start(p) {
            if !graph
                .nodes
                .get(p)
                .and_then(|n| n.childpat.as_ref())
                .map(|pat| pat.contains(trnode[j]))
                .unwrap_or(false)
            {
                return false;
            }
            // Allow transfrag edges that aren't on path if path has reachability
            if j > 0 && edge_bit(trpattern, graph, trnode[j - 1], trnode[j]) &&
                !edge_bit(pathpattern, graph, trnode[j - 1], trnode[j]) {
                let path_can_reach = graph
                    .nodes
                    .get(trnode[j - 1])
                    .and_then(|n| n.childpat.as_ref())
                    .map(|pat| pat.contains(trnode[j]))
                    .unwrap_or(false);
                if !path_can_reach {
                    return false;
                }
            }
        } else {
            j += 1;
        }
        if p == maxp {
            return true;
        }
        prevp = Some(p);
        let Some(np) = next_set_bit_by_coord(pathpattern, p, maxp, graph) else {
            return false;
        };
        p = np;
    }
}

fn fwd_to_sink_fast_long(
    seed_idx: usize,
    i: usize,
    path: &mut Vec<usize>,
    minpath: &mut usize,
    maxpath: &mut usize,
    pathpat: &mut GBitVec,
    transfrags: &[GraphTransfrag],
    graph: &mut Graph,
    nodecov: &[f64],
    require_longread: bool,
    diag: &mut LongRecDiag,
    visited: &mut HashSet<usize>,
    weak_cache: &mut Vec<Option<bool>>,
) -> bool {
    if !visited.insert(i) || path.len() > graph.n_nodes {
        diag.fwd_no_choice += 1;
        return false;
    }
    let gno = graph.n_nodes;
    let sink = graph.sink_id;
    let trace_fwd = trace_locus_active(graph, i) || trace_seed_active(seed_idx);
    let Some(inode) = graph.nodes.get(i) else {
        return false;
    };
    let has_sink_child = inode.children.contains(sink);
    if trace_fwd {
        eprintln!(
            "[TRACE_FWD] enter node i={} coord={}-{} minpath={} maxpath={} nchildren={} sink_child={}",
            i,
            inode.start,
            inode.end,
            *minpath,
            *maxpath,
            inode.children.count_ones(),
            has_sink_child
        );
        for c in inode.children.ones() {
            let ccoord = graph
                .nodes
                .get(c)
                .map(|n| (n.start, n.end))
                .unwrap_or((0, 0));
            eprintln!(
                "[TRACE_FWD]   path_edge child={} coord={}-{} edge_bit={} child_bit={} sink={}",
                c,
                ccoord.0,
                ccoord.1,
                edge_bit(pathpat, graph, i, c),
                pathpat.get_bit(c),
                c == sink
            );
        }
    }
    // Use genomic coordinates, not node IDs, for maxpath reachability check
    // Node IDs may not be ordered by coordinate (e.g., node 22 at 97417917 < node 9 at 97418573)
    let i_start = inode.start;
    let maxpath_start = node_start_or_zero(graph, *maxpath);
    // Sink node (no children) is the ultimate sink - any node can reach it
    let is_sink_node = inode.children.is_empty();
    if !is_sink_node
        && i_start < maxpath_start
        && !inode
            .childpat
            .as_ref()
            .map(|p| p.contains(*maxpath))
            .unwrap_or(false)
    {
        if trace_fwd {
            let max_coord = graph
                .nodes
                .get(*maxpath)
                .map(|n| (n.start, n.end))
                .unwrap_or((0, 0));
            eprintln!(
                "[TRACE_FWD]   FAIL: childpat does not contain maxpath={} coord={}-{} childpat_present={} i_before_maxpath={}",
                *maxpath,
                max_coord.0,
                max_coord.1,
                inode.childpat.is_some(),
                i_start < maxpath_start
            );
        }
        diag.fwd_unreachable_maxpath += 1;
        return false;
    }
    if inode.children.is_empty() {
        return true;
    }

    // Extract fields from inode before dropping the borrow so we can mutate graph.nodes
    // inside the child loop (for depleted-transfrag cleanup).
    let inode_end = inode.end;
    let children = inode.children.ones().collect::<Vec<_>>();

    let mut maxc: Option<usize> = None;
    let mut maxcov = 0.0;
    let mut tmax: Option<usize> = None;
    let mut exclude = false;
    let mut nextnode: Option<usize> = None;
    // Use genomic coordinates for reachability check, not node IDs
    let maxpath_coord = node_start_or_zero(graph, *maxpath);
    let i_coord = graph.nodes.get(i).map(|n| n.start).unwrap_or(0);
    let mut reach = maxpath_coord <= i_coord;

    // Fast-path, but only if `i+1` is an actual child edge. Node IDs are not guaranteed
    // coordinate-ordered, so we must not select `i+1` unless it is explicitly present in children.
    let next_id = i + 1;
    if next_id < gno && pathpat.get_bit(next_id) && children.iter().any(|&c| c == next_id) {
        maxc = Some(next_id);
        reach = true;
        if trace_fwd {
            eprintln!("[TRACE_FWD]   fast-path: i+1={} already on pathpat", next_id);
        }
    } else {
        for &c in &children {
            let childonpath = pathpat.get_bit(c);
            if edge_bit(pathpat, graph, i, c) {
                if trace_fwd {
                    eprintln!(
                        "[TRACE_FWD]   edge_bit fast-path: child={} coord={}-{}",
                        c,
                        graph.nodes.get(c).map(|n| n.start).unwrap_or(0),
                        graph.nodes.get(c).map(|n| n.end).unwrap_or(0)
                    );
                }
                maxc = Some(c);
                reach = true;
                tmax = None;
                break;
            }

            // Use genomic coordinates for reachability check
            if maxpath_coord > i_coord {
                if nextnode.is_none() {
                    nextnode = ((i + 2)..=*maxpath).find(|&j| pathpat.get_bit(j));
                }
                if let Some(nn) = nextnode {
                    if c != nn
                        && !graph
                            .nodes
                            .get(c)
                            .and_then(|n| n.childpat.as_ref())
                            .map(|p| p.contains(nn))
                            .unwrap_or(false)
                    {
                        continue;
                    }
                }
                reach = true;
            }

            // delete depleted transfrags from node's trf list
            if let Some(node) = graph.nodes.get_mut(c) {
                node.trf_ids.retain(|&t| {
                    transfrags
                        .get(t)
                        .map_or(false, |tf| tf.abundance >= EPS)
                });
            }
            let Some(cnode) = graph.nodes.get(c) else {
                continue;
            };
            let mut childcov = 0.0;
            let mut tchild: Option<usize> = None;
            let endpath = coord_max_node(graph, *maxpath, c);
            // : use minpath..endpath for onpath_long
            // compatibility, matching the original full-path-span check.
            // This prevents transfrags that are edge-compatible but path-incompatible
            // from inflating coverage, which causes over-extension.
            let compat_min = *minpath;
            let compat_max = endpath;
            // Coverage drop exclusion: if the child is coordinate-adjacent (touching)
            // and coverage drops significantly, skip it to avoid extending through
            // intergenic bridging nodes. StringTie uses `c == i+1` which works because
            // its node IDs are coordinate-ordered. After longtrim, Rustle's node IDs
            // are NOT coordinate-ordered, so we check coordinate adjacency directly.
            if inode_end == cnode.start
                && i < gno.saturating_sub(2)
                && c != sink
                && cnode.length() > 0
                && (nodecov.get(c).copied().unwrap_or(0.0) / cnode.length() as f64) < 1000.0
                && nodecov.get(i).copied().unwrap_or(0.0) * (DROP_FACTOR + ERROR_PERC)
                    > nodecov.get(c).copied().unwrap_or(0.0)
            {
                exclude = true;
            } else {
                pathpat.set_bit(c);
                edge_set(pathpat, graph, i, c, true);
                for &t in &cnode.trf_ids {
                    let Some(tf) = transfrags.get(t) else {
                        continue;
                    };
                    if tf.abundance < EPS {
                        if trace_fwd && c == sink {
                            eprintln!(
                                "[TRACE_FWD_SINK]   tf={} reject=depleted abundance={:.4}",
                                t, tf.abundance
                            );
                        }
                        continue;
                    }
                    if require_longread && !tf.longread {
                        if trace_fwd && c == sink {
                            eprintln!("[TRACE_FWD_SINK]   tf={} reject=not_longread", t);
                        }
                        continue;
                    }
                    if tf.node_ids.is_empty() {
                        if trace_fwd && c == sink {
                            eprintln!("[TRACE_FWD_SINK]   tf={} reject=empty_nodes", t);
                        }
                        continue;
                    }
                    // skip transfrags that start at source
                    let source = graph.source_id;
                    if require_longread && tf.node_ids.first() == Some(&source) {
                        if trace_fwd && c == sink {
                            eprintln!("[TRACE_FWD_SINK]   tf={} reject=starts_at_source", t);
                        }
                        continue;
                    }
                    let Some((first, last)) = first_last_raw(tf) else {
                        continue;
                    };
                    if c == sink {
                        // : For sink connections, accept only
                        // transfrags that start EXACTLY at current node and the path hasn't
                        // already extended forward. This prevents transfrags from earlier
                        // in the graph from inflating sink coverage, which causes
                        // over-extension (131 cases where Rustle adds extra exons).
                        let accepted = first == i && *maxpath <= i && last >= c;
                        if trace_fwd {
                            eprintln!(
                                "[TRACE_FWD_SINK]   tf={} first={} last={} abundance={:.4} maxpath={} i={} accepted={}",
                                t,
                                first,
                                last,
                                tf.abundance,
                                *maxpath,
                                i,
                                accepted
                            );
                        }
                        if accepted {
                            childcov += tf.abundance;
                            if tchild.is_none()
                                || tf.abundance > transfrags[tchild.unwrap()].abundance
                            {
                                tchild = Some(t);
                            }
                        }
                    } else if first <= i
                        && last >= c
                        && if require_longread {
                            onpath_long(
                                &tf.pattern,
                                &tf.node_ids,
                                pathpat,
                                compat_min,
                                compat_max,
                                graph,
                            )
                        } else {
                            onpath_short(
                                &tf.pattern,
                                &tf.node_ids,
                                pathpat,
                                compat_min,
                                compat_max,
                                graph,
                            )
                        }
                    {
                        if trace_fwd {
                            let weak = tf_weak_cached(t, transfrags, graph, nodecov, weak_cache);
                            let prev_tchild = tchild;
                            eprintln!(
                                "[FWD_TRACE seed={}]     trf t={} abund={:.1} first={} last={} weak={} ACCEPT pattern_pop={} nodes={:?}",
                                seed_idx,
                                t,
                                tf.abundance,
                                first,
                                last,
                                weak as u8,
                                tf.pattern.len_bits(),
                                &tf.node_ids[..tf.node_ids.len().min(8)]
                            );
                            replace_tmax(t, &mut tchild, transfrags, graph, nodecov);
                            eprintln!(
                                "[FWD_TRACE seed={}]       replace_tmax prev={:?} new={:?}",
                                seed_idx, prev_tchild, tchild
                            );
                        } else {
                            replace_tmax(t, &mut tchild, transfrags, graph, nodecov);
                        }
                        childcov += tf.abundance;
                    } else if trace_fwd && first <= i && last >= c {
                        let (onp, reason) = if require_longread {
                            let ok = onpath_long(
                                &tf.pattern,
                                &tf.node_ids,
                                pathpat,
                                compat_min,
                                compat_max,
                                graph,
                            );
                            let reason = onpath_long_reason(
                                &tf.pattern,
                                &tf.node_ids,
                                pathpat,
                                compat_min,
                                compat_max,
                                graph,
                            );
                            (ok, reason)
                        } else {
                            (false, "fail:short_mode_not_traced")
                        };
                        eprintln!(
                            "[FWD_TRACE seed={}]     trf t={} abund={:.1} first={} last={} FAIL onpath={} reason={} pattern_pop={} nodes={:?}",
                            seed_idx,
                            t,
                            tf.abundance,
                            first,
                            last,
                            onp,
                            reason,
                            tf.pattern.len_bits(),
                            &tf.node_ids[..tf.node_ids.len().min(8)]
                        );
                    }
                }
                if trace_fwd {
                    let ccoord = graph
                        .nodes
                        .get(c)
                        .map(|n| (n.start, n.end))
                        .unwrap_or((0, 0));
                    let ntrf = graph.nodes.get(c).map(|n| n.trf_ids.len()).unwrap_or(0);
                    eprintln!(
                        "[TRACE_FWD]   child c={} coord={}-{} ntrf={} childcov={:.1} tchild={:?}",
                        c, ccoord.0, ccoord.1, ntrf, childcov, tchild
                    );
                    if let Some(tid) = tchild {
                        if let Some(tf) = transfrags.get(tid) {
                            let first = tf.node_ids.first().copied().unwrap_or(usize::MAX);
                            let last = tf.node_ids.last().copied().unwrap_or(usize::MAX);
                            eprintln!(
                                "[TRACE_FWD]     support_tf={} first={} last={} nodes={} longstart={} longend={} tf_edge={} tf_sink_edge={}",
                                tid,
                                first,
                                last,
                                tf.node_ids.len(),
                                tf.longstart,
                                tf.longend,
                                edge_bit(&tf.pattern, graph, i, c),
                                edge_bit(&tf.pattern, graph, i, sink)
                            );
                        }
                    }
                }
                if childcov > maxcov {
                    let cur_tmax = tmax;
                    let cand_tchild = tchild;
                    let cur_weak = cur_tmax
                        .map(|tm| tf_weak_cached(tm, transfrags, graph, nodecov, weak_cache));
                    let cand_weak = cand_tchild
                        .map(|tc| tf_weak_cached(tc, transfrags, graph, nodecov, weak_cache));
                    let allow = tmax.is_none()
                        || tf_weak_cached(tmax.unwrap(), transfrags, graph, nodecov, weak_cache)
                        || tchild
                            .map(|tc| !tf_weak_cached(tc, transfrags, graph, nodecov, weak_cache))
                            .unwrap_or(false);
                    if trace_fwd {
                        eprintln!(
                            "[FWD_TRACE seed={}]   maxcov_update child={} childcov={:.4} prev_maxcov={:.4} prev_tmax={:?} prev_weak={:?} cand_tchild={:?} cand_weak={:?} allow={}",
                            seed_idx,
                            c,
                            childcov,
                            maxcov,
                            cur_tmax,
                            cur_weak,
                            cand_tchild,
                            cand_weak,
                            allow as u8
                        );
                    }
                    if allow {
                        maxcov = childcov;
                        maxc = Some(c);
                        tmax = tchild;
                    }
                } else if maxc.is_some() && childcov > maxcov - EPS {
                    let cur_tmax = tmax;
                    let cand_tchild = tchild;
                    let cur_weak = cur_tmax
                        .map(|tm| tf_weak_cached(tm, transfrags, graph, nodecov, weak_cache));
                    let cand_weak = cand_tchild
                        .map(|tc| tf_weak_cached(tc, transfrags, graph, nodecov, weak_cache));
                    let allow = tmax
                        .map(|tm| tf_weak_cached(tm, transfrags, graph, nodecov, weak_cache))
                        .unwrap_or(true)
                        || tchild
                            .map(|tc| !tf_weak_cached(tc, transfrags, graph, nodecov, weak_cache))
                            .unwrap_or(false);
                    let nodecov_old = nodecov.get(maxc.unwrap()).copied().unwrap_or(0.0);
                    let nodecov_new = nodecov.get(c).copied().unwrap_or(0.0);
                    if trace_fwd {
                        eprintln!(
                            "[FWD_TRACE seed={}]   tiebreak child={} childcov={:.4} maxcov={:.4} prev_maxc={:?} prev_tmax={:?} prev_weak={:?} cand_tchild={:?} cand_weak={:?} nodecov_old={:.4} nodecov_new={:.4} allow={}",
                            seed_idx,
                            c,
                            childcov,
                            maxcov,
                            maxc,
                            cur_tmax,
                            cur_weak,
                            cand_tchild,
                            cand_weak,
                            nodecov_old,
                            nodecov_new,
                            allow as u8
                        );
                    }
                    if allow && nodecov_old < nodecov_new {
                        maxc = Some(c);
                        tmax = tchild;
                    }
                }
                edge_set(pathpat, graph, i, c, false);
                if childonpath {
                    break;
                }
                pathpat.clear_bit(c);
            }
        }
    }
    if !reach {
        if trace_fwd {
            eprintln!(
                "[TRACE_FWD]   FAIL: no reach (maxpath={} > i={})",
                *maxpath, i
            );
        }
        diag.fwd_no_reach += 1;
        return false;
    }
    if maxc.is_none() {
        if exclude && i + 1 < gno && nodecov.get(i + 1).copied().unwrap_or(0.0) > 0.0 {
            let Some(cnode) = graph.nodes.get(i + 1) else {
                return false;
            };
            let mut childcov = 0.0;
            let mut tlocal: Option<usize> = None;
            let endpath = (*maxpath).max(i + 1);
            pathpat.set_bit(i + 1);
            edge_set(pathpat, graph, i, i + 1, true);
            for &t in &cnode.trf_ids {
                let Some(tf) = transfrags.get(t) else {
                    continue;
                };
                if tf.abundance < EPS
                    || (require_longread && !tf.longread)
                    || tf.node_ids.is_empty()
                {
                    continue;
                }
                let Some((first, last)) = first_last_raw(tf) else {
                    continue;
                };
                if first <= i
                    && last >= i + 1
                    && if require_longread {
                        onpath_long(
                            &tf.pattern,
                            &tf.node_ids,
                            pathpat,
                            *minpath,
                            endpath,
                            graph,
                        )
                    } else {
                        onpath_short(
                            &tf.pattern,
                            &tf.node_ids,
                            pathpat,
                            *minpath,
                            endpath,
                            graph,
                        )
                    }
                {
                    childcov += tf.abundance;
                    replace_tmax(t, &mut tlocal, transfrags, graph, nodecov);
                }
            }
            edge_set(pathpat, graph, i, i + 1, false);
            pathpat.clear_bit(i + 1);
            if childcov > 0.0 {
                maxc = Some(i + 1);
                tmax = tlocal;
            } else {
                diag.fwd_exclude_no_support += 1;
                return false;
            }
        } else {
            diag.fwd_no_choice += 1;
            return false;
        }
    }

    let c = maxc.unwrap();
    if trace_fwd {
        let ccoord = graph
            .nodes
            .get(c)
            .map(|n| (n.start, n.end))
            .unwrap_or((0, 0));
        eprintln!(
            "[TRACE_FWD]   CHOSE child c={} coord={}-{} maxcov={:.1} tmax={:?}",
            c, ccoord.0, ccoord.1, maxcov, tmax
        );
    }
    path.push(c);
    pathpat.set_bit(c);
    edge_set(pathpat, graph, i, c, true);
    if let Some(t) = tmax {
        pathpat.or_assign(&transfrags[t].pattern);
        // node IDs are coordinate-ordered; in Rust, explicitly
        // maintain min/max path endpoints by coordinate to match semantics.
        if let Some(&first) = transfrags[t].node_ids.first() {
            if node_start_or_zero(graph, first) < node_start_or_zero(graph, *minpath) {
                *minpath = first;
            }
        }
        if let Some(&last) = transfrags[t].node_ids.last() {
            if node_start_or_zero(graph, last) > node_start_or_zero(graph, *maxpath) {
                *maxpath = last;
            }
        }
    }
    fwd_to_sink_fast_long(
        seed_idx,
        c,
        path,
        minpath,
        maxpath,
        pathpat,
        transfrags,
        graph,
        nodecov,
        require_longread,
        diag,
        visited,
        weak_cache,
    )
}

fn back_to_source_fast_long(
    seed_idx: usize,
    i: usize,
    path: &mut Vec<usize>,
    minpath: &mut usize,
    maxpath: &mut usize,
    pathpat: &mut GBitVec,
    transfrags: &[GraphTransfrag],
    graph: &mut Graph,
    nodecov: &[f64],
    require_longread: bool,
    prefer_non_source: bool,
    diag: &mut LongRecDiag,
    visited: &mut HashSet<usize>,
    weak_cache: &mut Vec<Option<bool>>,
) -> bool {
    if !visited.insert(i) || path.len() > graph.n_nodes {
        diag.back_no_choice += 1;
        if trace_seed_diagnostics(seed_idx) {
            eprintln!(
                "[BACK_FAIL seed={}] node={} reason=cycle_or_overflow path_len={} n_nodes={}",
                seed_idx,
                i,
                path.len(),
                graph.n_nodes
            );
        }
        return false;
    }
    let gno = graph.n_nodes;
    let source = graph.source_id;
    let trace_back = trace_locus_active(graph, i) || trace_seed_diagnostics(seed_idx);
    let Some(inode) = graph.nodes.get(i) else {
        if trace_seed_diagnostics(seed_idx) {
            eprintln!(
                "[BACK_FAIL seed={}] node={} reason=node_not_found",
                seed_idx, i
            );
        }
        return false;
    };
    // Copy fields out of inode to avoid holding the immutable borrow across mutable graph access.
    let inode_start = inode.start;
    let inode_end = inode.end;
    let inode_parents = inode.parents.ones().collect::<Vec<_>>();
    // Use genomic coordinates, not node IDs, for minpath reachability check.
    let minpath_start = node_start_or_zero(graph, *minpath);
    // Source node (no parents) is the ultimate source - it can reach any node
    let is_source_node = inode_parents.is_empty();
    let inode_parentpat_contains_minpath = !is_source_node
        && minpath_start > inode_start
        && !inode
            .parentpat
            .as_ref()
            .map(|pp| pp.contains(*minpath))
            .unwrap_or(false);
    if trace_back {
        eprintln!(
            "[BACK_TRACE seed={}] enter node i={} coord={}-{} minpath={} maxpath={} nparents={} path_len={}",
            seed_idx, i, inode_start, inode_end, *minpath, *maxpath, inode_parents.len(), path.len()
        );
    }
    if inode_parentpat_contains_minpath {
        if trace_back {
            eprintln!(
                "[BACK_FAIL seed={}] node={} coord={}-{} reason=unreachable_minpath minpath={}",
                seed_idx, i, inode_start, inode_end, *minpath
            );
        }
        diag.back_unreachable_minpath += 1;
        diag.back_parents_skipped_reachability += 1;
        return false;
    }
    if inode_parents.is_empty() {
        if trace_back {
            eprintln!(
                "[BACK_TRACE seed={}] reached source-like terminal node={}",
                seed_idx, i
            );
        }
        return true;
    }

    let parents = inode_parents;
    let has_non_source_parent = parents.iter().any(|&p| p != source);

    let mut maxp: Option<usize> = None;
    let mut maxcov = 0.0;
    let mut tmax: Option<usize> = None;
    let mut exclude = false;
    let mut nextnode: Option<usize> = None;
    let mut reach = *minpath >= i;

    // Fast-path, but only if `i-1` is an actual parent edge. Node IDs are not guaranteed
    // coordinate-ordered, so we must not select `i-1` unless it is explicitly present in parents.
    let prev_id = i.saturating_sub(1);
    if i > 0 && pathpat.get_bit(prev_id) && parents.iter().any(|&p| p == prev_id) {
        maxp = Some(prev_id);
        reach = true;
        if trace_back {
            eprintln!("[TRACE_BACK]   fast-path: i-1={} already on pathpat", prev_id);
        }
    } else {
        for &p in &parents {
            // In max-sensitivity mode we often want to stitch partial long-read seeds upstream.
            // The code effectively gates source edges (keepsource/keepsink); our port can
            // otherwise over-prefer `source` due to large aggregated parentcov. Prefer real
            // (non-source) parents when they exist, and only fall back to source if needed.
            if prefer_non_source && has_non_source_parent && p == source {
                continue;
            }
            diag.back_parents_checked += 1;
            let parentonpath = pathpat.get_bit(p);
            if edge_bit(pathpat, graph, p, i) {
            // Accept unconditionally when the edge is already on the active path pattern.
                if trace_back {
                    eprintln!(
                        "[TRACE_BACK]   edge_bit fast-path: parent={} coord={}-{}",
                        p,
                        graph.nodes.get(p).map(|n| n.start).unwrap_or(0),
                        graph.nodes.get(p).map(|n| n.end).unwrap_or(0)
                    );
                }
                maxp = Some(p);
                reach = true;
                tmax = None;
                break;
            }

            if *minpath < i {
                if nextnode.is_none() && i >= 2 {
                    // the original algorithm's effective behavior here is to look for the closest
                    // already-on-path node *before* i, not to scan forward into the
                    // current node / downstream path. The forward scan blocks valid
                    // parent choices at loci like 23035526-23062635 by choosing the
                    // current sink-side node as the reachability target.
                    nextnode = ((*minpath)..=(i - 2)).rev().find(|&j| pathpat.get_bit(j));
                }
                if let Some(nn) = nextnode {
                    if p != nn
                        && !graph
                            .nodes
                            .get(p)
                            .and_then(|n| n.parentpat.as_ref())
                            .map(|pat| pat.contains(nn))
                            .unwrap_or(false)
                    {
                    continue;
                }
            }
                reach = true;
            }

            // Delete depleted transfrags from the node's support list before scoring parents.
            if let Some(node) = graph.nodes.get_mut(p) {
                node.trf_ids.retain(|&t| {
                    transfrags
                        .get(t)
                        .map_or(false, |tf| tf.abundance >= EPS)
                });
            }
            let Some(pnode) = graph.nodes.get(p) else {
                continue;
            };
            let mut parentcov = 0.0;
            let mut tpar: Option<usize> = None;
            // Match the original algorithm: compatibility is evaluated against the full current path span,
            // from the earliest node seen so far (minpath or parent) to the current global maxpath.
            let startpath = coord_min_node(graph, *minpath, p);
            let endpath = coord_max_node(graph, *maxpath, i);
            // Coverage drop exclusion for back_to_source: same coordinate-based fix
            // as fwd_to_sink. Check adjacency by coordinate, not by node ID.
            if inode_start == pnode.end
                && i > 1
                && p != source
                && pnode.length() > 0
                && (nodecov.get(p).copied().unwrap_or(0.0) / pnode.length() as f64) < 1000.0
                && nodecov.get(i).copied().unwrap_or(0.0) * (DROP_FACTOR + ERROR_PERC)
                    > nodecov.get(p).copied().unwrap_or(0.0)
            {
                exclude = true;
            } else {
                pathpat.set_bit(p);
                edge_set(pathpat, graph, p, i, true);
                for &t in &pnode.trf_ids {
                    diag.back_transfrags_checked += 1;
                    let Some(tf) = transfrags.get(t) else {
                        continue;
                    };
                    if tf.abundance < EPS {
                        diag.back_transfrags_skipped_depleted += 1;
                        if trace_seed_diagnostics(seed_idx) {
                            eprintln!("[BACK_TRACE seed={}]     skip trf={} reason=depleted abundance={:.4}",
                                seed_idx, t, tf.abundance);
                        }
                        continue;
                    }
                    if require_longread && !tf.longread {
                        diag.back_transfrags_skipped_notlong += 1;
                        if trace_seed_diagnostics(seed_idx) {
                            eprintln!(
                                "[BACK_TRACE seed={}]     skip trf={} reason=not_longread",
                                seed_idx, t
                            );
                        }
                        continue;
                    }
                    if tf.node_ids.is_empty() {
                        if trace_seed_diagnostics(seed_idx) {
                            eprintln!(
                                "[BACK_TRACE seed={}]     skip trf={} reason=empty_nodes",
                                seed_idx, t
                            );
                        }
                        continue;
                    }
                    // Skip transfrags that already terminate at sink.
                    let sink = graph.sink_id;
                    if require_longread && tf.node_ids.last() == Some(&sink) {
                        diag.back_transfrags_skipped_ends_at_sink += 1;
                        if trace_seed_diagnostics(seed_idx) {
                            eprintln!(
                                "[BACK_TRACE seed={}]     skip trf={} reason=ends_at_sink",
                                seed_idx, t
                            );
                        }
                        continue;
                    }
                    let Some((first, last)) = first_last_raw(tf) else {
                        continue;
                    };
                    if p == source {
                        // Only count source support when the transfrag ends exactly at `i`
                        // and `i` is at/left of the current minpath.
                        let min_start = node_start_or_zero(graph, *minpath);
                        let exact_source_match = last == i && min_start >= inode_start;
                        if trace_back {
                            eprintln!(
                                "[TRACE_BACK_SOURCE]   tf={} first={} last={} abundance={:.4} minpath={} exact={} accepted={}",
                                t, first, last, tf.abundance, *minpath, exact_source_match,
                                exact_source_match
                            );
                        }
                        if exact_source_match {
                            parentcov += tf.abundance;
                            if tpar.is_none() || tf.abundance > transfrags[tpar.unwrap()].abundance
                            {
                                tpar = Some(t);
                            }
                        }
                    } else if node_start_or_zero(graph, first) <= node_start_or_zero(graph, p)
                        && node_start_or_zero(graph, last) >= inode_start
                        && if require_longread {
                            onpath_long(
                                &tf.pattern,
                                &tf.node_ids,
                                pathpat,
                                startpath,
                                endpath,
                                graph,
                            )
                        } else {
                            onpath_short(
                                &tf.pattern,
                                &tf.node_ids,
                                pathpat,
                                startpath,
                                endpath,
                                graph,
                            )
                        }
                    {
                        if trace_back {
                            let weak = tf_weak_cached(t, transfrags, graph, nodecov, weak_cache);
                            let prev_tpar = tpar;
                            eprintln!(
                                "[BACK_TRACE seed={}]     trf t={} abund={:.1} first={} last={} weak={} ACCEPT pattern_pop={} nodes={:?}",
                                seed_idx,
                                t,
                                tf.abundance,
                                first,
                                last,
                                weak as u8,
                                tf.pattern.len_bits(),
                                &tf.node_ids[..tf.node_ids.len().min(8)]
                            );
                            replace_tmax(t, &mut tpar, transfrags, graph, nodecov);
                            eprintln!(
                                "[BACK_TRACE seed={}]       replace_tmax prev={:?} new={:?}",
                                seed_idx, prev_tpar, tpar
                            );
                        } else {
                            replace_tmax(t, &mut tpar, transfrags, graph, nodecov);
                        }
                        parentcov += tf.abundance;
                    } else if trace_back
                        && node_start_or_zero(graph, first) <= node_start_or_zero(graph, p)
                        && node_start_or_zero(graph, last) >= inode_start
                    {
                        diag.back_transfrags_skipped_onpath += 1;
                        let (onp, reason) = if require_longread {
                            let ok = onpath_long(
                                &tf.pattern,
                                &tf.node_ids,
                                pathpat,
                                startpath,
                                endpath,
                                graph,
                            );
                            let reason = onpath_long_reason(
                                &tf.pattern,
                                &tf.node_ids,
                                pathpat,
                                startpath,
                                endpath,
                                graph,
                            );
                            (ok, reason)
                        } else {
                            (false, "fail:short_mode_not_traced")
                        };
                        eprintln!("[BACK_TRACE seed={}]     trf t={} abund={:.1} first={} last={} FAIL onpath={} reason={} pattern_pop={} nodes={:?}",
                            seed_idx, t, tf.abundance, first, last, onp, reason, tf.pattern.len_bits(),
                            &tf.node_ids[..tf.node_ids.len().min(8)]);
                    }
                }
                if trace_back {
                    let pcoord = graph
                        .nodes
                        .get(p)
                        .map(|n| (n.start, n.end))
                        .unwrap_or((0, 0));
                    let ntrf = pnode.trf_ids.len();
                    let ndepleted = pnode
                        .trf_ids
                        .iter()
                        .filter(|&&t| {
                            transfrags
                                .get(t)
                                .map_or(true, |tf| tf.abundance < EPS)
                        })
                        .count();
                    let nskip_long = if require_longread {
                        pnode
                            .trf_ids
                            .iter()
                            .filter(|&&t| {
                                transfrags
                                    .get(t)
                                    .map_or(false, |tf| tf.abundance >= EPS && !tf.longread)
                            })
                            .count()
                    } else {
                        0
                    };
                    eprintln!("[TRACE_BACK]   parent p={} coord={}-{} ntrf={} depleted={} skip_notlong={} parentcov={:.1} tpar={:?}",
                        p, pcoord.0, pcoord.1, ntrf, ndepleted, nskip_long, parentcov, tpar);
                }
                if parentcov > maxcov {
                    let cur_tmax = tmax;
                    let cand_tpar = tpar;
                    let cur_weak = cur_tmax
                        .map(|tm| tf_weak_cached(tm, transfrags, graph, nodecov, weak_cache));
                    let cand_weak = cand_tpar
                        .map(|tp| tf_weak_cached(tp, transfrags, graph, nodecov, weak_cache));
                    let allow = tmax.is_none()
                        || tf_weak_cached(tmax.unwrap(), transfrags, graph, nodecov, weak_cache)
                        || tpar
                            .map(|tp| !tf_weak_cached(tp, transfrags, graph, nodecov, weak_cache))
                            .unwrap_or(false);
                    if trace_back {
                        eprintln!(
                            "[BACK_TRACE seed={}]   maxcov_update parent={} parentcov={:.4} prev_maxcov={:.4} prev_tmax={:?} prev_weak={:?} cand_tpar={:?} cand_weak={:?} allow={}",
                            seed_idx,
                            p,
                            parentcov,
                            maxcov,
                            cur_tmax,
                            cur_weak,
                            cand_tpar,
                            cand_weak,
                            allow as u8
                        );
                    }
                    if allow {
                        maxcov = parentcov;
                        maxp = Some(p);
                        tmax = tpar;
                    }
                } else if maxp.is_some() && parentcov > maxcov - EPS {
                    let cur_tmax = tmax;
                    let cand_tpar = tpar;
                    let cur_weak = cur_tmax
                        .map(|tm| tf_weak_cached(tm, transfrags, graph, nodecov, weak_cache));
                    let cand_weak = cand_tpar
                        .map(|tp| tf_weak_cached(tp, transfrags, graph, nodecov, weak_cache));
                    let allow = tmax
                        .map(|tm| tf_weak_cached(tm, transfrags, graph, nodecov, weak_cache))
                        .unwrap_or(true)
                        || tpar
                            .map(|tp| !tf_weak_cached(tp, transfrags, graph, nodecov, weak_cache))
                            .unwrap_or(false);
                    let nodecov_old = nodecov.get(maxp.unwrap()).copied().unwrap_or(0.0);
                    let nodecov_new = nodecov.get(p).copied().unwrap_or(0.0);
                    if trace_back {
                        eprintln!(
                            "[BACK_TRACE seed={}]   tiebreak parent={} parentcov={:.4} maxcov={:.4} prev_maxp={:?} prev_tmax={:?} prev_weak={:?} cand_tpar={:?} cand_weak={:?} nodecov_old={:.4} nodecov_new={:.4} allow={}",
                            seed_idx,
                            p,
                            parentcov,
                            maxcov,
                            maxp,
                            cur_tmax,
                            cur_weak,
                            cand_tpar,
                            cand_weak,
                            nodecov_old,
                            nodecov_new,
                            allow as u8
                        );
                    }
                    if allow && nodecov_old < nodecov_new {
                        maxp = Some(p);
                        tmax = tpar;
                    }
                }
                edge_set(pathpat, graph, p, i, false);
                if parentonpath {
                    break;
                }
                pathpat.clear_bit(p);
            }
        }
    }

    if !reach {
        if trace_back {
            eprintln!(
                "[BACK_FAIL seed={}] node={} coord={}-{} reason=no_reach minpath={} i={}",
                seed_idx, i, inode_start, inode_end, *minpath, i
            );
        }
        diag.back_no_reach += 1;
        if trace_seed_diagnostics(seed_idx) {
            eprintln!("[BACK_FAIL seed={}] SUMMARY: parents_checked={} transfrags_checked={} skipped_reachability={} skipped_onpath={}",
                seed_idx, diag.back_parents_checked, diag.back_transfrags_checked,
                diag.back_parents_skipped_reachability, diag.back_transfrags_skipped_onpath);
        }
        return false;
    }
    if maxp.is_none() {
        // A real source-terminal seed can fail here if explicit source support was
        // zeroed during keep-source gating even though the upstream path is coherent.
        // In that narrow case, accept the source transition directly.
        if parents.len() == 1
            && parents[0] == source
            && *minpath == i
        {
            if trace_back {
                eprintln!(
                    "[TRACE_BACK]   SOURCE_TERMINAL_FALLBACK node={} coord={}-{} minpath={} => allow",
                    i, inode_start, inode_end, *minpath
                );
            }
            return true;
        }
        if trace_back {
            eprintln!(
                "[BACK_FAIL seed={}] node={} coord={}-{} maxp=None exclude={} maxcov={:.1} parents_checked={}",
                seed_idx, i, inode_start, inode_end, exclude, maxcov, parents.len()
            );
        }
        if trace_seed_diagnostics(seed_idx) {
            eprintln!("[BACK_FAIL seed={}] SUMMARY: parents_checked={} transfrags_checked={} skipped_depleted={} skipped_notlong={} skipped_ends_at_sink={} skipped_onpath={}",
                seed_idx, diag.back_parents_checked, diag.back_transfrags_checked,
                diag.back_transfrags_skipped_depleted, diag.back_transfrags_skipped_notlong,
                diag.back_transfrags_skipped_ends_at_sink, diag.back_transfrags_skipped_onpath);
        }
        if exclude && i > 0 && nodecov.get(i - 1).copied().unwrap_or(0.0) > 0.0 {
            let Some(pnode) = graph.nodes.get(i - 1) else {
                return false;
            };
            let mut parentcov = 0.0;
            let mut tlocal: Option<usize> = None;
            let startpath = (*minpath).min(i - 1);
            pathpat.set_bit(i - 1);
            edge_set(pathpat, graph, i - 1, i, true);
            for &t in &pnode.trf_ids {
                let Some(tf) = transfrags.get(t) else {
                    continue;
                };
                if tf.abundance < EPS
                    || (require_longread && !tf.longread)
                    || tf.node_ids.is_empty()
                {
                    continue;
                }
                let Some((first, last)) = first_last_raw(tf) else {
                    continue;
                };
                if first <= i - 1
                    && last >= i
                    && if require_longread {
                        onpath_long(
                            &tf.pattern,
                            &tf.node_ids,
                            pathpat,
                            startpath,
                            *maxpath,
                            graph,
                        )
                    } else {
                        onpath_short(
                            &tf.pattern,
                            &tf.node_ids,
                            pathpat,
                            startpath,
                            *maxpath,
                            graph,
                        )
                    }
                {
                    parentcov += tf.abundance;
                    replace_tmax(t, &mut tlocal, transfrags, graph, nodecov);
                }
            }
            edge_set(pathpat, graph, i - 1, i, false);
            pathpat.clear_bit(i - 1);
            if parentcov > 0.0 {
                maxp = Some(i - 1);
                tmax = tlocal;
            } else {
                diag.back_exclude_no_support += 1;
                return false;
            }
        } else {
            diag.back_no_choice += 1;
            return false;
        }
    }

    let p = maxp.unwrap();
    if trace_back {
        let pcoord = graph
            .nodes
            .get(p)
            .map(|n| (n.start, n.end))
            .unwrap_or((0, 0));
        eprintln!(
            "[TRACE_BACK]   CHOSE parent p={} coord={}-{} maxcov={:.1} tmax={:?}",
            p, pcoord.0, pcoord.1, maxcov, tmax
        );
    }
    if p != source {
        path.push(p);
    }
    pathpat.set_bit(p);
    edge_set(pathpat, graph, p, i, true);
    if let Some(t) = tmax {
        if trace_back {
            let tf = &transfrags[t];
            let coords = tf
                .node_ids
                .iter()
                .take(10)
                .map(|&nid| {
                    graph
                        .nodes
                        .get(nid)
                        .map(|n| format!("{nid}({}-{})", n.start, n.end))
                        .unwrap_or_else(|| format!("{nid}(0-0)"))
                })
                .collect::<Vec<_>>()
                .join(" ");
            eprintln!(
                "[TRACE_BACK]   PATHPAT_OR via tmax={} abund={:.4} first={} last={} nodes={}",
                t,
                tf.abundance,
                tf.node_ids.first().copied().unwrap_or(usize::MAX),
                tf.node_ids.last().copied().unwrap_or(usize::MAX),
                coords
            );
        }
        pathpat.or_assign(&transfrags[t].pattern);
        // Keep min/max endpoints in coordinate space to match node-id ordering.
        if let Some(&first) = transfrags[t].node_ids.first() {
            if node_start_or_zero(graph, first) < node_start_or_zero(graph, *minpath) {
                *minpath = first;
            }
        }
        if let Some(&last) = transfrags[t].node_ids.last() {
            if node_start_or_zero(graph, last) > node_start_or_zero(graph, *maxpath) {
                *maxpath = last;
            }
        }
    }
    let _ = gno; // mirror; used in edge id operations.
    back_to_source_fast_long(
        seed_idx,
        p,
        path,
        minpath,
        maxpath,
        pathpat,
        transfrags,
        graph,
        nodecov,
        require_longread,
        prefer_non_source,
        diag,
        visited,
        weak_cache,
    )
}

/// Extract transcripts from transfrags with path extension (parse_trflong).
///
/// Process transfrags in abundance-descending order (longtrCmp):
/// 1. Build initial path from transfrag nodes
/// 2. Extend left/right through splice junctions using compatible transfrags
/// 3. Also extend contiguously at terminal exons
/// 4. Apply poly-tail trimming (strand-specific)
/// parse_trflong: seed ordering for long-read transcript extraction.
fn parse_trflong(transfrags: &[GraphTransfrag], _graph: &Graph) -> Vec<usize> {
    let mut seeded: Vec<usize> = transfrags
        .iter()
        .enumerate()
        .filter(|(_, tf)| {
            tf.trflong_seed && tf.weak == 0 && !tf.node_ids.is_empty() && tf.usepath >= 0
        })
        .map(|(i, _)| i)
        .collect();
    // ref: trflong order is determined by keeptrf insertion plus reverse iteration.
    // usepath encodes the vector order, so ascending usepath reproduces consumption.
    seeded.sort_unstable_by_key(|&i| transfrags[i].usepath);
    seeded
}

/// parse_trf: short-read ordering for transcript extraction (mixed mode).
fn parse_trf(_graph: &Graph, transfrags: &[GraphTransfrag]) -> Vec<usize> {
    let mut indices: Vec<usize> = transfrags
        .iter()
        .enumerate()
        .filter(|(_, tf)| tf.srabund > crate::constants::FLOW_EPSILON || !tf.longread)
        .map(|(i, _)| i)
        .collect();
    indices.sort_unstable_by(|&i, &j| {
        let a = &transfrags[i];
        let b = &transfrags[j];
        let eff_a = a.srabund + a.abundance;
        let eff_b = b.srabund + b.abundance;
        let ab = eff_b
            .partial_cmp(&eff_a)
            .unwrap_or(std::cmp::Ordering::Equal);
        if ab != std::cmp::Ordering::Equal {
            return ab;
        }
        b.node_ids.len().cmp(&a.node_ids.len())
    });
    indices
}

/// 5. Enforce hardstart/hardend constraints
/// 6. Build exons from extended path
/// 7. Zero main transfrag's abundance to prevent duplication
pub fn extract_transcripts(
    graph: &mut Graph,
    transfrags: &mut [GraphTransfrag],
    bundle_chrom: &str,
    bundle_strand: char,
    bundle_id: &str,
    config: &RunConfig,
    nascent: bool,
    mut seed_outcomes: Option<&mut Vec<(usize, SeedOutcome)>>,
    longrec_summary: Option<&mut LongRecSummary>,
) -> Vec<Transcript> {
    // Helper macro to record seed outcome when tracing is active.
    macro_rules! record_outcome {
        ($idx:expr, $outcome:expr) => {
            if let Some(ref mut outcomes) = seed_outcomes {
                outcomes.push(($idx, $outcome));
            }
        };
    }
    let audit_zero_flux = std::env::var_os("RUSTLE_AUDIT_ZERO_FLUX").is_some();
    let depletion_diag = std::env::var_os("RUSTLE_DEPLETION_DIAG").is_some();
    let plumb_debug = std::env::var_os("RUSTLE_PLUMB_DEBUG").is_some();
    if std::env::var("RUSTLE_DEBUG_EK").is_ok() {
        let first_real = graph
            .nodes
            .iter()
            .enumerate()
            .find(|(i, _)| *i != graph.source_id && *i != graph.sink_id)
            .map(|(_, n)| n.start)
            .unwrap_or(0);
        let last_real = graph
            .nodes
            .iter()
            .enumerate()
            .rev()
            .find(|(i, _)| *i != graph.source_id && *i != graph.sink_id)
            .map(|(_, n)| n.end)
            .unwrap_or(0);
        let n_seeds = transfrags.iter().filter(|tf| tf.trflong_seed).count();
        eprintln!(
            "[RUST_PARSE] extract_transcripts gno={} trflong={} locus={}-{}",
            graph.nodes.len(),
            n_seeds,
            first_real,
            last_real
        );
    }
    // Locus-level flow trace: set RUSTLE_TRACE_LOCUS=start-end to trace all seed processing.
    let trace_locus: Option<(u64, u64)> = std::env::var("RUSTLE_TRACE_LOCUS").ok().and_then(|v| {
        let parts: Vec<&str> = v.split('-').collect();
        if parts.len() == 2 {
            Some((parts[0].parse().ok()?, parts[1].parse().ok()?))
        } else {
            None
        }
    });
    let trace_chain_active = config.trace_intron_chain.is_some()
        && config
            .debug_bundle
            .as_ref()
            .map_or(true, |b| bundle_id.contains(b));
    let trace_chain = config.trace_intron_chain.as_ref();
    let mut _trace_chain_emitted = false;
    let trace_intron = parse_trace_intron();
    let mut out = Vec::new();
    let mut kept_paths: Vec<(Vec<usize>, f64, bool, usize)> = Vec::new(); // (inner nodes, support, guide, out_idx)
    let mut previous_guides: Vec<GuideFlowState> = Vec::new();
    let mut checktrf: Vec<usize> = Vec::new();
    // Default to strict behavior: failed direct long-rec seeds are deferred to checktrf.
    // Allow explicit opt-out only for diagnostics.
    let strict_longrec_checktrf_deferral =
        std::env::var_os("RUSTLE_DISABLE_STRICT_LONGREC_CHECKTRF_DEFERRAL").is_none();

    let mut longrec_attempted = 0usize;
    let mut longrec_succeeded = 0usize;
    let mut longrec_fallback = 0usize;
    let mut longrec_back_fail = 0usize;
    let mut longrec_fwd_fail = 0usize;
    let mut longrec_path_invalid = 0usize;
    let mut longrec_back_unreachable_minpath = 0usize;
    let mut longrec_back_no_reach = 0usize;
    let mut longrec_back_no_choice = 0usize;
    let mut longrec_back_exclude_no_support = 0usize;
    let mut longrec_fwd_unreachable_maxpath = 0usize;
    let mut longrec_fwd_no_reach = 0usize;
    let mut longrec_fwd_no_choice = 0usize;
    let mut longrec_fwd_exclude_no_support = 0usize;
    let mut zero_flux_set: HashSet<usize> = Default::default();
    let mut zero_flux_candidates = 0usize;
    let mut zero_flux_rescued = 0usize;
    let zero_flux_dropped_nomatch = 0usize;
    let mut zero_flux_dropped_lowabund = 0usize;
    let mut zero_flux_dropped_empty = 0usize;
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let mode = config.assembly_mode();
    let long_read_mode = mode.is_long_read();
    let mixed_mode = mode.is_mixed();
    let debug_ek = std::env::var("RUSTLE_DEBUG_EK").is_ok();
    let guided_mode = transfrags
        .iter()
        .any(|tf| tf.guide || tf.guide_tid.is_some());
    let _lr_witness_splice_junctions: Option<HashSet<(u64, u64)>> = if long_read_mode {
        Some(build_lr_witness_splice_junctions(graph, transfrags))
    } else {
        None
    };

    let mut order: Vec<usize> = Vec::new();
    if mixed_mode {
        // Mixed mode: long-pass by usepath/trflong, then short-pass by parse_trf ordering.
        let mut long_order: Vec<usize> = (0..transfrags.len()).collect();
        long_order.sort_unstable_by_key(|&a| transfrags[a].usepath);
        let mut seen: HashSet<usize> = Default::default();
        for i in long_order {
            if transfrags[i].usepath < -1 {
                seen.insert(i);
                order.push(i);
            }
        }
        let short_order = parse_trf(graph, transfrags);
        for i in short_order {
            if seen.insert(i) {
                order.push(i);
            }
        }
    } else if long_read_mode {
        order = parse_trflong(transfrags, graph);

        // Max-sensitivity diagnostic mode: also consider non-keeptrf long-read transfrags as seeds.
        // This helps recover transcripts that are stitchable from partials but never appear as a
        // keeptrf representative (common root cause of `not_extracted:junctions_present`).
        if config.max_sensitivity {
            let mut seen = vec![false; transfrags.len()];
            for &i in &order {
                if i < seen.len() {
                    seen[i] = true;
                }
            }
            let mut extra: Vec<usize> = (0..transfrags.len())
                .filter(|&i| {
                    let tf = &transfrags[i];
                    !seen[i]
                        && !tf.node_ids.is_empty()
                        && tf.longread
                        && tf.abundance > 0.0
                        && (tf.weak == 0
                            || (tf.real && !tf.killed_junction_orphan && !tf.coverage_weak))
                })
                .collect();
            extra.sort_unstable_by(|&a, &b| {
                let ta = &transfrags[a];
                let tb = &transfrags[b];
                let ab = tb
                    .abundance
                    .partial_cmp(&ta.abundance)
                    .unwrap_or(std::cmp::Ordering::Equal);
                if ab != std::cmp::Ordering::Equal {
                    return ab;
                }
                tb.node_ids.len().cmp(&ta.node_ids.len())
            });
            order.extend(extra);
        }
    } else {
        order = (0..transfrags.len()).collect();
        order.sort_unstable_by(|&a, &b| {
            transfrags[b]
                .abundance
                .partial_cmp(&transfrags[a].abundance)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
    }
    // Shared state across all paths in the bundle (nodecov, capacity from transfrag abundances).
    // Each path depletes local_nodecov and (via flow) transfrag abundances; no restore between paths.
    let mut local_nodecov: Vec<f64> = graph.nodes.iter().map(|n| n.nodecov.max(0.0)).collect();
    let local_noderate: Vec<f64> = graph
        .nodes
        .iter()
        .map(|n| if n.noderate > 0.0 { n.noderate } else { 1.0 })
        .collect();
    if loop_trace_active() {
        for _ in 0..graph.nodes.len() {
            eprintln!(
                "LOOP_find_transcripts: nodecov_compute gno={}",
                graph.nodes.len()
            );
        }
    }

    // seed ordering is fully determined by usepath (set in process_transfrags).
    // simply iterates trflong backwards — no second reordering.
    // usepath encodes: completes first (lowest usepath), incompletes second, within each group
    // highest-abundance first. parse_trflong's sort_by(usepath ascending) already matches this.

    let debug_detail = std::env::var_os("RUSTLE_DEBUG_DETAIL").is_some();
    let pathpat_trace = pathpat_trace_active();
    if debug_detail {
        eprintln!(
            "DEBUG_GRAPH gno={} source={} sink={} bundle={}",
            graph.nodes.len(),
            source_id,
            sink_id,
            bundle_id
        );
        for (t, tf) in transfrags.iter().enumerate() {
            if tf.longread && tf.abundance > 0.0 {
                eprint!("DEBUG_TF t={} abund={:.4} longstart={} longend={} seed={} usepath={} nodes={}:",
                    t, tf.abundance, tf.longstart, tf.longend, tf.trflong_seed, tf.usepath, tf.node_ids.len());
                for &nid in &tf.node_ids {
                    if let Some(n) = graph.nodes.get(nid) {
                        eprint!(" {}({}-{})", nid, n.start, n.end);
                    } else {
                        eprint!(" {}(?)", nid);
                    }
                }
                eprintln!();
            }
        }
        eprintln!(
            "DEBUG_ORDER len={} first_5={:?}",
            order.len(),
            &order[..order.len().min(5)]
        );
    }

    // cache tf_weak per-transfrag. computes transfrag[t]->weak lazily on
    // first access and caches the result. Subsequent checks use the cached value, reflecting
    // nodecov at the time of first evaluation (not current depleted nodecov).
    let mut weak_cache: Vec<Option<bool>> = vec![None; transfrags.len()];

    // Macro to emit DEBUG_SEED_PROC line at each seed exit point
    macro_rules! emit_debug_seed_proc {
        ($idx:expr, $graph:expr, $transfrags:expr, $real_nodes:expr,
         $back_ok:expr, $fwd_ok:expr, $flux:expr, $path_coords:expr, $outcome:expr) => {
            if debug_detail {
                let span_start = $real_nodes.first()
                    .and_then(|&n| $graph.nodes.get(n))
                    .map(|n| n.start).unwrap_or(0);
                let span_end = $real_nodes.last()
                    .and_then(|&n| $graph.nodes.get(n))
                    .map(|n| n.end).unwrap_or(0);
                let path_str = $path_coords.as_deref().unwrap_or("");
                eprintln!("DEBUG_SEED_PROC idx={} span={}-{} abund={:.4} nodes={} back={} fwd={} flux={:.4} path=[{}] outcome={}",
                    $idx, span_start, span_end, $transfrags[$idx].abundance,
                    $real_nodes.len(),
                    if $back_ok { "ok" } else { "fail" },
                    if $fwd_ok { "ok" } else { "fail" },
                    $flux, path_str, $outcome);
            }
        };
    }

    // Macro to emit DEBUG_SEED_DECISION trace (matches format)
    macro_rules! emit_debug_seed_decision {
        // Pattern with just idx, used_direct, branch (e.g., eonly_skip, HARD_BOUNDARY, UNWITNESSED)
        ($idx:expr, $used_direct:expr, $branch:expr) => {
            if debug_detail {
                eprintln!("DEBUG_SEED_DECISION t={} used_direct={} branch={}",
                    $idx, $used_direct, $branch);
            }
        };
        // Pattern with flux, cov, nexons (e.g., STORED, DIRECT_LOW_COV)
        ($idx:expr, $used_direct:expr, $branch:expr, $flux:expr, $cov:expr, $nexons:expr) => {
            if debug_detail {
                eprintln!("DEBUG_SEED_DECISION t={} used_direct={} branch={} flux={:.4} cov={:.4} nexons={}",
                    $idx, $used_direct, $branch, $flux, $cov, $nexons);
            }
        };
    }

    // Macro to emit DEBUG_CHK trace for checktrf outcomes (matches format)
    macro_rules! emit_debug_chk {
        ($t:expr, $abund:expr, $outcome:expr, $matched:expr, $abundsum:expr) => {
            if debug_detail {
                eprintln!("DEBUG_CHK t={} abund={:.4} outcome={} matched={} abundsum={:.4}",
                    $t, $abund, $outcome, $matched, $abundsum);
            }
        };
        ($t:expr, $abund:expr, $outcome:expr) => {
            if debug_detail {
                eprintln!("DEBUG_CHK t={} abund={:.4} outcome={}",
                    $t, $abund, $outcome);
            }
        };
    }

    let ntrflong = transfrags.iter().filter(|tf| tf.trflong_seed).count();
    for idx in order {
        if loop_trace_active() && long_read_mode && transfrags[idx].trflong_seed {
            eprintln!(
                "LOOP_parse_trflong: ntrflong={} nasc={}",
                ntrflong,
                if nascent { 1 } else { 0 }
            );
        }
        // Parity seed processing state — collected through the loop, emitted at exit points
        let mut debug_back_ok = false;
        let mut debug_fwd_ok = false;
        let mut debug_flux: f64 = 0.0;
        // Will be set once we have use_path
        let mut debug_path_coords: Option<String> = None;

        // Trace intron: report every seed that contains the target intron
        if let Some(ti) = trace_intron {
            if path_has_intron(&transfrags[idx].node_ids, graph, ti) {
                eprintln!(
                    "[TRACE_INTRON_TF] stage=seed_entry idx={} intron={}-{} abund={:.4} nodes={} seed={} weak={} longstart={} longend={} bundle={}",
                    idx, ti.0, ti.1, transfrags[idx].abundance,
                    transfrags[idx].node_ids.len(), transfrags[idx].trflong_seed,
                    transfrags[idx].weak, transfrags[idx].longstart, transfrags[idx].longend,
                    bundle_id
                );
            }
        }
        // Early trace: check if this transfrag overlaps the trace locus
        if let Some((lo, hi)) = trace_locus {
            let overlaps = transfrags[idx].node_ids.iter().any(|&nid| {
                graph
                    .nodes
                    .get(nid)
                    .map_or(false, |n| n.start <= hi && n.end >= lo)
            });
            if overlaps {
                let first_nid = transfrags[idx]
                    .node_ids
                    .first()
                    .and_then(|&nid| graph.nodes.get(nid));
                let last_nid = transfrags[idx]
                    .node_ids
                    .last()
                    .and_then(|&nid| graph.nodes.get(nid));
                let skip_reason = if transfrags[idx].node_ids.is_empty() {
                    "empty"
                } else if transfrags[idx].weak != 0 {
                    "weak"
                } else if !mixed_mode
                    && long_read_mode
                    && !config.max_sensitivity
                    && !transfrags[idx].trflong_seed
                {
                    "not_seed"
                } else {
                    "pass"
                };
                eprintln!(
                    "[TRACE_EARLY] idx={} abund={:.1} nodes={} seed={} weak={} skip={} span={}-{}",
                    idx,
                    transfrags[idx].abundance,
                    transfrags[idx].node_ids.len(),
                    transfrags[idx].trflong_seed,
                    transfrags[idx].weak,
                    skip_reason,
                    first_nid.map(|n| n.start).unwrap_or(0),
                    last_nid.map(|n| n.end).unwrap_or(0)
                );
            }
        }
        if transfrags[idx].node_ids.is_empty() {
            record_outcome!(idx, SeedOutcome::Skipped("empty"));
            continue;
        }
        if transfrags[idx].weak != 0
            && !(config.max_sensitivity
                && transfrags[idx].longread
                && transfrags[idx].real
                && !transfrags[idx].killed_junction_orphan
                && !transfrags[idx].coverage_weak)
        {
            record_outcome!(idx, SeedOutcome::Skipped("weak"));
            continue;
        }
        if !mixed_mode && long_read_mode && !config.max_sensitivity && !transfrags[idx].trflong_seed
        {
            record_outcome!(idx, SeedOutcome::Skipped("not_seed"));
            continue;
        }
        // nascent second pass processes only non-guide transfrags.
        if nascent && transfrags[idx].guide {
            record_outcome!(idx, SeedOutcome::Skipped("nascent_guide"));
            continue;
        }
        let is_single = transfrags[idx]
            .node_ids
            .iter()
            .filter(|&&n| n != source_id && n != sink_id)
            .count()
            <= 1;
        let thresh = if is_single {
            config.singlethr
        } else {
            config.readthr
        };
        let effective_support = if mixed_mode {
            transfrags[idx].abundance.max(transfrags[idx].srabund)
        } else {
            transfrags[idx].abundance
        };
        // parse_trflong does not apply readthr/singlethr as an early gate in long/mixed
        // long-read parsing; it defers low-support handling to flux/checktrf stages.
        if !long_read_mode && effective_support < thresh {
            record_outcome!(idx, SeedOutcome::Skipped("low_support"));
            continue;
        }

        let base_real_nodes: Vec<usize> = transfrags[idx]
            .node_ids
            .iter()
            .copied()
            .filter(|&n| n != source_id && n != sink_id)
            .collect();
        if base_real_nodes.is_empty() {
            record_outcome!(idx, SeedOutcome::Skipped("no_real_nodes"));
            continue;
        }
        let real_nodes: Vec<usize> = if long_read_mode && transfrags[idx].longread {
            materialize_longread_seed_nodes(
                graph,
                &base_real_nodes,
                transfrags[idx].longstart,
                transfrags[idx].longend,
            )
        } else {
            base_real_nodes.clone()
        };
        if real_nodes.is_empty() {
            record_outcome!(idx, SeedOutcome::Skipped("materialized_empty"));
            continue;
        }
        // single-node non-guide transfrags in long-read mode are never
        // stored as predictions. the original implementation zeroes their abundance and the cov computation
        // yields 0 for non-guides, so they always fail the store gate.
        // Skip them BEFORE path extension and flow extraction to prevent flow budget
        // depletion — this is the #1 precision/sensitivity fix, recovering ~500 TPs.
        if long_read_mode && real_nodes.len() == 1 && !transfrags[idx].guide {
            transfrags[idx].abundance = 0.0;
            record_outcome!(idx, SeedOutcome::Skipped("single_node_lr"));
            continue;
        }

        let debug_nodes: &[usize] = if transfrags[idx].longread {
            &base_real_nodes
        } else {
            &real_nodes
        };
        if trace_locus.map_or(false, |(lo, hi)| {
            real_nodes.iter().any(|&nid| {
                graph
                    .nodes
                    .get(nid)
                    .map_or(false, |n| n.start <= hi && n.end >= lo)
            })
        }) && real_nodes
            != transfrags[idx]
                .node_ids
                .iter()
                .copied()
                .filter(|&n| n != source_id && n != sink_id)
                .collect::<Vec<_>>()
        {
            let orig_last = transfrags[idx]
                .node_ids
                .iter()
                .rev()
                .copied()
                .find(|&n| n != source_id && n != sink_id)
                .unwrap_or(source_id);
            let new_last = *real_nodes.last().unwrap_or(&orig_last);
            eprintln!(
                "[TRACE_SEED_MAT] idx={} orig_last={} coord={}-{} new_last={} coord={}-{} longstart={} longend={}",
                idx,
                orig_last,
                graph.nodes.get(orig_last).map(|n| n.start).unwrap_or(0),
                graph.nodes.get(orig_last).map(|n| n.end).unwrap_or(0),
                new_last,
                graph.nodes.get(new_last).map(|n| n.start).unwrap_or(0),
                graph.nodes.get(new_last).map(|n| n.end).unwrap_or(0),
                transfrags[idx].longstart,
                transfrags[idx].longend
            );
        }

        if debug_ek {
            let tf = &transfrags[idx];
            eprintln!("[RUST_SEED] idx={} abund={:.6} guide={} nodes={} first={} last={} longstart={} longend={}",
                idx, tf.abundance, tf.guide, tf.node_ids.len(),
                tf.node_ids.first().unwrap_or(&0), tf.node_ids.last().unwrap_or(&0),
                tf.longstart, tf.longend);
        }

        // Snapshot values before mutable borrows
        let abundance = transfrags[idx].abundance;
        // read_count is the pre-depletion read mass (never modified after initialization).
        // Use this for longcov so it's stable regardless of how many competing transfrags
        // split the flow — this makes longcov comparable to the original algorithm's values.
        let read_count_snapshot = transfrags[idx].read_count;
        if depletion_diag && abundance <= 0.0 {
            eprintln!(
                "ZERO_AB_SEED idx={} nodes={}",
                idx,
                transfrags[idx].node_ids.len()
            );
        }
        let poly_start = transfrags[idx].poly_start_unaligned;
        let poly_end = transfrags[idx].poly_end_unaligned;
        let seed_pattern = transfrags[idx].pattern.clone();
        let mut tocheck = true;

        // Build path with optional strict direct long-recursion trial.
        let mut path: Vec<usize> = Vec::new();
        let mut flow_pathpat: Option<GBitVec> = None;
        let try_direct_longrec = long_read_mode && config.strict_longrec_port && !mixed_mode;
        let mut used_direct = false;
        let watched_tfs = parse_trace_tf_ids();
        if try_direct_longrec {
            longrec_attempted += 1;
            let mut diag = LongRecDiag::default();
            let mut pathpat = seed_pattern.clone();
            let seed_nodes = debug_nodes;
            // node IDs are coordinate-ordered; derive min/max endpoints by
            // coordinate explicitly since Rust node IDs are not guaranteed ordered.
            let mut minp = seed_nodes[0];
            let mut maxp = seed_nodes[0];
            for &nid in seed_nodes.iter().skip(1) {
                minp = coord_min_node(graph, minp, nid);
                maxp = coord_max_node(graph, maxp, nid);
            }
            let seed_minp = minp;
            let seed_maxp = maxp;
            trace_outgoing_edge_state("seed.raw", idx, maxp, &pathpat, graph);
            // Match helper traversal uses the stored transfrag nodes/pattern, not
            // the later materialized terminal-contiguous node expansion.
            for &nid in seed_nodes {
                pathpat.set_bit(nid);
            }
            trace_outgoing_edge_state("seed.after_nodes", idx, maxp, &pathpat, graph);
            // seed helper traversal sets source→minp and
            // maxp→sink edge bits so back_to_source / fwd_to_sink can take fast-path decisions.
            //
            // However, for partial long-read seeds this can prematurely "seal" the 5'/3' ends:
            // back_to_source immediately takes the source endcap and never stitches upstream
            // (observed at STRG.171 for junction-stitchable isoforms).
            //
            // In max-sensitivity mode, only set the endcaps when the seed endpoints are hard
            // boundaries (or guide-derived). This preserves standard behavior for trusted
            // endpoints while allowing stitching for partials.
            let minp_hardstart = graph
                .nodes
                .get(minp)
                .map(|n| n.hardstart)
                .unwrap_or(false);
            let maxp_hardend = graph
                .nodes
                .get(maxp)
                .map(|n| n.hardend)
                .unwrap_or(false);
            let minp_has_non_source_parent = graph
                .nodes
                .get(minp)
                .map(|n| n.parents.ones().any(|p| p != source_id))
                .unwrap_or(false);
            // Contiguous extension parent within the same exon (parent end == node start).
            // This is the specific case where the source-helper can wrongly dominate parentcov
            // and truncate the 5' end (e.g. STRG.171 seed 55: 33->105).
            // Only prefer non-source parent extension when the seed's first node
            // has a contiguous parent AND does NOT have source as parent.
            // If source is a parent, this node is a valid starting point (from a
            // junction or hardstart boundary) and back_to_source should stop here
            // rather than extending backward through contiguous exonic regions.
            let minp_has_source_parent = graph
                .nodes
                .get(minp)
                .map(|n| n.parents.contains(source_id))
                .unwrap_or(false);
            let _minp_has_contig_parent = !minp_has_source_parent
                && graph
                    .nodes
                    .get(minp)
                    .map(|n| {
                        n.parents.ones().any(|p| {
                            p != source_id
                                && graph
                                    .nodes
                                    .get(p)
                                    .map(|pn| pn.end == n.start)
                                    .unwrap_or(false)
                        })
                    })
                    .unwrap_or(false);
            let maxp_has_non_sink_child = graph
                .nodes
                .get(maxp)
                .map(|n| n.children.ones().any(|c| c != sink_id))
                .unwrap_or(false);
            if transfrags[idx].guide || (minp_hardstart && !minp_has_non_source_parent) {
                edge_set(&mut pathpat, graph, source_id, minp, true);
            }
            if transfrags[idx].guide || (maxp_hardend && !maxp_has_non_sink_child) {
                edge_set(&mut pathpat, graph, maxp, sink_id, true);
            }
            trace_outgoing_edge_state("seed.after_endcaps", idx, maxp, &pathpat, graph);
            if pathpat_trace {
                let (nb, eb, bits) = pathpat_summary(&pathpat, graph.n_nodes);
                eprintln!(
                    "parse_trflong: SEED idx={} abund={:.4} nnodes={} minp={} maxp={}",
                    idx,
                    abundance,
                    seed_nodes.len(),
                    minp,
                    maxp
                );
                eprintln!(
                    "parse_trflong: SEED_PATHPAT idx={} node_bits={} edge_bits={} bits={}",
                    idx, nb, eb, bits
                );
            }
            let maxi = minp;
            path.push(maxi);
            pathpat.set_bit(maxi);
            // Trace seed info for diagnostic locus analysis.
            if trace_locus_active(graph, minp) || trace_locus_active(graph, maxp) {
                let min_coord = graph
                    .nodes
                    .get(minp)
                    .map(|n| (n.start, n.end))
                    .unwrap_or((0, 0));
                let max_coord = graph
                    .nodes
                    .get(maxp)
                    .map(|n| (n.start, n.end))
                    .unwrap_or((0, 0));
                eprintln!("[TRACE_SEED] idx={} abund={:.1} nnodes={} minp={} ({}-{}) maxp={} ({}-{}) patpop={}",
                    idx, abundance, seed_nodes.len(), minp, min_coord.0, min_coord.1,
                    maxp, max_coord.0, max_coord.1, pathpat.len_bits());
            }
            trace_tf_watch(
                idx,
                "before_direct_longrec",
                &watched_tfs,
                transfrags,
                graph,
            );
            trace_node_trf_watch(idx, "source_before", source_id, graph, &watched_tfs);
            trace_node_trf_watch(idx, "minp_before", minp, graph, &watched_tfs);
            let mut visited_back: HashSet<usize> = Default::default();
            let mut visited_fwd: HashSet<usize> = Default::default();
            let back_ok = back_to_source_fast_long(
                idx,
                maxi,
                &mut path,
                &mut minp,
                &mut maxp,
                &mut pathpat,
                transfrags,
                graph,
                &local_nodecov,
                true,
                // Keep strict the original algorithm-style selection; we only apply a narrow internal override
                // when a contiguous exon-extension parent exists (see back_to_source_fast_long).
                // Never prefer non-source parents. Always let back_to_source stop
                // at source when available. This prevents over-extension that
                // consumes flow budget from other isoforms.
                std::env::var_os("RUSTLE_BACK_PREFER_NON_SOURCE").is_some(),
                &mut diag,
                &mut visited_back,
                &mut weak_cache,
            );
            let mut fwd_ok = false;
            if pathpat_trace || trace_seed_diagnostics(idx) || trace_locus_active(graph, seed_minp)
            {
                eprintln!(
                    "--- parse_trflong: back_to_source t={} result={} path_len={} path: {}",
                    idx,
                    if back_ok { "OK" } else { "FAIL" },
                    path.len(),
                    format_path_nodes(&path, graph, 30)
                );
            }
            if !back_ok {
                trace_tf_watch(idx, "after_back_fail", &watched_tfs, transfrags, graph);
                trace_node_trf_watch(
                    idx,
                    "source_after_back_fail",
                    source_id,
                    graph,
                    &watched_tfs,
                );
                trace_node_trf_watch(idx, "minp_after_back_fail", minp, graph, &watched_tfs);
            }
            if back_ok {
                path.push(source_id);
                path.reverse();
                fwd_ok = fwd_to_sink_fast_long(
                    idx,
                    maxi,
                    &mut path,
                    &mut minp,
                    &mut maxp,
                    &mut pathpat,
                    transfrags,
                    graph,
                    &local_nodecov,
                    true,
                    &mut diag,
                    &mut visited_fwd,
                    &mut weak_cache,
                );
                if pathpat_trace
                    || trace_seed_diagnostics(idx)
                    || trace_locus_active(graph, seed_maxp)
                {
                    eprintln!(
                        "--- parse_trflong: fwd_to_sink t={} result={} path_len={} path: {}",
                        idx,
                        if fwd_ok { "OK" } else { "FAIL" },
                        path.len(),
                        format_path_nodes(&path, graph, 30)
                    );
                }
                if !fwd_ok {
                    trace_tf_watch(idx, "after_fwd_fail", &watched_tfs, transfrags, graph);
                    trace_node_trf_watch(
                        idx,
                        "source_after_fwd_fail",
                        source_id,
                        graph,
                        &watched_tfs,
                    );
                    trace_node_trf_watch(idx, "maxp_after_fwd_fail", maxp, graph, &watched_tfs);
                }
                if fwd_ok {
                    // parse_trflong only checks back/fwd success.
                    // It does not enforce an extra source/sink endpoint validity gate here.
                    used_direct = true;
                }
                // NOTE: simple fwd fallback tested but causes regression (1596 vs 1606).
                // The fallback creates paths that produce wrong transcripts competing
                // with good ones. Checktrf rescue is a better outcome for most failed seeds.
            }
            debug_back_ok = back_ok;
            debug_fwd_ok = fwd_ok;
            if used_direct {
                longrec_succeeded += 1;
                flow_pathpat = Some(pathpat.clone());
                if pathpat_trace {
                    let (nb, eb, bits) = pathpat_summary(&pathpat, graph.n_nodes);
                    eprintln!(
                        "parse_trflong: FLOW_PATHPAT idx={} source=direct node_bits={} edge_bits={} bits={}",
                        idx, nb, eb, bits
                    );
                }
            } else {
                if let Some((lo, hi)) = trace_locus {
                    if real_nodes.iter().any(|&nid| {
                        graph
                            .nodes
                            .get(nid)
                            .map_or(false, |n| n.start <= hi && n.end >= lo)
                    }) {
                        let seed_min_coord = graph
                            .nodes
                            .get(seed_minp)
                            .map(|n| (n.start, n.end))
                            .unwrap_or((0, 0));
                        let seed_max_coord = graph
                            .nodes
                            .get(seed_maxp)
                            .map(|n| (n.start, n.end))
                            .unwrap_or((0, 0));
                        let cur_min_coord = graph
                            .nodes
                            .get(minp)
                            .map(|n| (n.start, n.end))
                            .unwrap_or((0, 0));
                        let cur_max_coord = graph
                            .nodes
                            .get(maxp)
                            .map(|n| (n.start, n.end))
                            .unwrap_or((0, 0));
                        let path_str = if path.is_empty() {
                            "-".to_string()
                        } else {
                            path.iter()
                                .map(|&nid| {
                                    let coord = graph
                                        .nodes
                                        .get(nid)
                                        .map(|n| format!("{}-{}", n.start, n.end))
                                        .unwrap_or_else(|| "0-0".to_string());
                                    format!("{nid}:{coord}")
                                })
                                .collect::<Vec<_>>()
                                .join(" ")
                        };
                        eprintln!(
                            "[TRACE_LONGREC_FAIL] idx={} back_ok={} fwd_ok={} seed_minp={} ({}-{}) seed_maxp={} ({}-{}) cur_minp={} ({}-{}) cur_maxp={} ({}-{}) path_len={} path={} diag={{back_unreachable_minpath:{} back_no_reach:{} back_no_choice:{} back_exclude_no_support:{} fwd_unreachable_maxpath:{} fwd_no_reach:{} fwd_no_choice:{} fwd_exclude_no_support:{}}}",
                            idx,
                            back_ok,
                            fwd_ok,
                            seed_minp,
                            seed_min_coord.0,
                            seed_min_coord.1,
                            seed_maxp,
                            seed_max_coord.0,
                            seed_max_coord.1,
                            minp,
                            cur_min_coord.0,
                            cur_min_coord.1,
                            maxp,
                            cur_max_coord.0,
                            cur_max_coord.1,
                            path.len(),
                            path_str,
                            diag.back_unreachable_minpath,
                            diag.back_no_reach,
                            diag.back_no_choice,
                            diag.back_exclude_no_support,
                            diag.fwd_unreachable_maxpath,
                            diag.fwd_no_reach,
                            diag.fwd_no_choice,
                            diag.fwd_exclude_no_support
                        );
                        trace_outgoing_edge_state(
                            "longrec.fail.seed_max",
                            idx,
                            seed_maxp,
                            &pathpat,
                            graph,
                        );
                        if maxp != seed_maxp {
                            trace_outgoing_edge_state(
                                "longrec.fail.cur_max",
                                idx,
                                maxp,
                                &pathpat,
                                graph,
                            );
                        }
                    }
                }
                if !back_ok {
                    longrec_back_fail += 1;
                } else if !fwd_ok {
                    longrec_fwd_fail += 1;
                } else {
                    longrec_path_invalid += 1;
                }
                path.clear();
                longrec_fallback += 1;
            }
            longrec_back_unreachable_minpath += diag.back_unreachable_minpath;
            longrec_back_no_reach += diag.back_no_reach;
            longrec_back_no_choice += diag.back_no_choice;
            longrec_back_exclude_no_support += diag.back_exclude_no_support;
            longrec_fwd_unreachable_maxpath += diag.fwd_unreachable_maxpath;
            longrec_fwd_no_reach += diag.fwd_no_reach;
            longrec_fwd_no_choice += diag.fwd_no_choice;
            longrec_fwd_exclude_no_support += diag.fwd_exclude_no_support;
            if config.verbose && !used_direct {
                eprintln!(
                    "      longrec_fail idx={} back_ok={} fwd_ok={} diag={{back_unreachable_minpath:{} back_no_reach:{} back_no_choice:{} back_exclude_no_support:{} fwd_unreachable_maxpath:{} fwd_no_reach:{} fwd_no_choice:{} fwd_exclude_no_support:{}}}",
                    idx,
                    back_ok,
                    fwd_ok,
                    diag.back_unreachable_minpath,
                    diag.back_no_reach,
                    diag.back_no_choice,
                    diag.back_exclude_no_support,
                    diag.fwd_unreachable_maxpath,
                    diag.fwd_no_reach,
                    diag.fwd_no_choice,
                    diag.fwd_exclude_no_support
                );
            }
        }

        if !used_direct {
            if !long_read_mode {
                let mut diag = LongRecDiag::default();
                let mut pathpat = seed_pattern.clone();
                // find minp/maxp by genomic coordinate
                let (minp_idx, maxp_idx) = real_nodes.iter().enumerate()
                    .map(|(i, &n)| (i, graph.nodes.get(n).map(|node| node.start).unwrap_or(0)))
                    .fold((0, 0), |(min_i, max_i), (i, start)| {
                        let min_start = graph.nodes.get(real_nodes[min_i]).map(|n| n.start).unwrap_or(0);
                        let max_start = graph.nodes.get(real_nodes[max_i]).map(|n| n.start).unwrap_or(0);
                        let new_min_i = if start < min_start { i } else { min_i };
                        let new_max_i = if start > max_start { i } else { max_i };
                        (new_min_i, new_max_i)
                    });
                let mut minp = real_nodes[minp_idx];
                let mut maxp = real_nodes[maxp_idx];
                // add edge bits for seed's consecutive node pairs.
                for &nid in &real_nodes {
                    pathpat.set_bit(nid);
                }
                for w in real_nodes.windows(2) {
                    edge_set(&mut pathpat, graph, w[0], w[1], true);
                }
                // source→minp and maxp→sink edge bits.
                edge_set(&mut pathpat, graph, source_id, minp, true);
                edge_set(&mut pathpat, graph, maxp, sink_id, true);
                let maxi = minp;
                path.push(maxi);
                pathpat.set_bit(maxi);
                let mut visited_back: HashSet<usize> = Default::default();
                let mut visited_fwd: HashSet<usize> = Default::default();
                if back_to_source_fast_long(
                    idx,
                    maxi,
                    &mut path,
                    &mut minp,
                    &mut maxp,
                    &mut pathpat,
                    transfrags,
                    graph,
                    &local_nodecov,
                    false,
                    false,
                    &mut diag,
                    &mut visited_back,
                    &mut weak_cache,
                ) {
                    debug_back_ok = true;
                    path.push(source_id);
                    path.reverse();
                    if fwd_to_sink_fast_long(
                        idx,
                        maxi,
                        &mut path,
                        &mut minp,
                        &mut maxp,
                        &mut pathpat,
                        transfrags,
                        graph,
                        &local_nodecov,
                        false,
                        &mut diag,
                        &mut visited_fwd,
                        &mut weak_cache,
                    ) {
                        debug_fwd_ok = true;
                        used_direct = true;
                    }
                }
                if !used_direct {
                    path.clear();
                }
            }
        }

        if !used_direct {
            if long_read_mode
                && !mixed_mode
                && try_direct_longrec
                && strict_longrec_checktrf_deferral
            {
                // when back/fwd fails in parse_trflong, defer
                // this seed to checktrf rescue and do not build a raw-node fallback path.
                if let Some((lo, hi)) = trace_locus {
                    if real_nodes.iter().any(|&nid| {
                        graph
                            .nodes
                            .get(nid)
                            .map_or(false, |n| n.start <= hi && n.end >= lo)
                    }) {
                        eprintln!(
                            "[TRACE_LOCUS] idx={} abund={:.1} LONGREC_FAIL → checktrf span={}-{}",
                            idx,
                            abundance,
                            graph
                                .nodes
                                .get(debug_nodes[0])
                                .map(|n| n.start)
                                .unwrap_or(0),
                            graph
                                .nodes
                                .get(*debug_nodes.last().unwrap())
                                .map(|n| n.end)
                                .unwrap_or(0)
                        );
                    }
                }
                trace_checktrf_enqueue(
                    idx,
                    "longrec_fail",
                    abundance,
                    &transfrags[idx],
                    debug_nodes,
                    graph,
                    None,
                    Some(debug_back_ok),
                    Some(debug_fwd_ok),
                    None,
                    None,
                    None,
                );
                checktrf.push(idx);
                emit_debug_seed_proc!(
                    idx,
                    graph,
                    transfrags,
                    debug_nodes,
                    debug_back_ok,
                    debug_fwd_ok,
                    debug_flux,
                    debug_path_coords,
                    "longrec_fail_checktrf"
                );
                if !debug_back_ok {
                    record_outcome!(idx, SeedOutcome::BackToSourceFail);
                } else {
                    record_outcome!(idx, SeedOutcome::FwdToSinkFail);
                }
                continue;
            }
            // Fallback: create path from raw nodes with extensions.
            let mut inner_path = real_nodes.clone();
            if config.enable_splice_path_extension
                && (transfrags[idx].guide || transfrags[idx].guide_tid.is_some())
            {
                extend_path_left_unambiguous_splice(&mut inner_path, graph, transfrags);
                extend_path_right_unambiguous_splice(&mut inner_path, graph, transfrags);
            }
            // In long-read mode, avoid unconditional contiguous graph extension:
            // it can create unsupported exon-chain combinations (gffcompare class 'j').
            if !long_read_mode {
                extend_path_left_contiguous(&mut inner_path, graph);
                extend_path_right_contiguous(&mut inner_path, graph);
            }

            path.push(source_id);
            path.extend_from_slice(&inner_path);
            path.push(sink_id);
        }

        let startnode = 1usize;
        let lastnode = path.len().saturating_sub(2);
        if lastnode < startnode || lastnode >= path.len() {
            record_outcome!(idx, SeedOutcome::Skipped("path_too_short"));
            continue;
        }

        // (parse_trflong): poly-tail trimming rewrites the concrete path,
        // not just [start,last] indices.
        let mut use_path = path;
        let mut use_start = 1usize;
        let mut use_last = use_path.len().saturating_sub(2);
        if use_last < use_start || use_last >= use_path.len() {
            record_outcome!(idx, SeedOutcome::Skipped("path_too_short"));
            continue;
        }

        // '+' strand: poly_end = 3' polyA → trim extended path back to tf's last real node
        if bundle_strand == '+' && poly_end >= POLY_TAIL_STOP_COUNT {
            let tf_last = *real_nodes.last().unwrap();
            if use_path[use_last] > tf_last {
                if let Some(cutpos) = (use_start..=use_last).find(|&pi| use_path[pi] == tf_last) {
                    if cutpos < use_last {
                        let mut newpath: Vec<usize> = Vec::with_capacity(cutpos + 2);
                        newpath.push(source_id);
                        newpath.extend_from_slice(&use_path[1..=cutpos]);
                        newpath.push(sink_id);
                        use_path = newpath;
                        if flow_pathpat.is_some() {
                            flow_pathpat =
                                Some(rebuild_flow_pathpat(&use_path, &seed_pattern, graph));
                        }
                        use_start = 1;
                        use_last = use_path.len().saturating_sub(2);
                    }
                }
            }
        }
        // '-' strand: poly_start = 3' polyT → trim path forward to tf's first real node
        if bundle_strand == '-' && poly_start >= POLY_TAIL_STOP_COUNT {
            let tf_first = real_nodes[0];
            if use_path[use_start] < tf_first {
                if let Some(keep_from) = (use_start..=use_last).find(|&pi| use_path[pi] == tf_first)
                {
                    if keep_from > use_start {
                        let mut newpath: Vec<usize> =
                            Vec::with_capacity(use_last.saturating_sub(keep_from) + 3);
                        newpath.push(source_id);
                        newpath.extend_from_slice(&use_path[keep_from..=use_last]);
                        newpath.push(sink_id);
                        use_path = newpath;
                        if flow_pathpat.is_some() {
                            flow_pathpat =
                                Some(rebuild_flow_pathpat(&use_path, &seed_pattern, graph));
                        }
                        use_start = 1;
                        use_last = use_path.len().saturating_sub(2);
                    }
                }
            }
        }

        // Hardstart/hardend enforcement (10277-10320):
        // If tf has both hard boundaries but extended path doesn't match, use tf's nodes directly
        let tf_first_node = real_nodes[0];
        let tf_last_node = *real_nodes.last().unwrap();
        let thardstart = graph
            .nodes
            .get(tf_first_node)
            .map(|n| n.hardstart)
            .unwrap_or(false);
        let thardend = graph
            .nodes
            .get(tf_last_node)
            .map(|n| n.hardend)
            .unwrap_or(false);

        let mut checkpath = true;
        if thardstart
            && thardend
            && (!same_node_coords(graph, use_path[use_start], tf_first_node)
                || !same_node_coords(graph, use_path[use_last], tf_last_node))
        {
            // if hard-boundary seed's path doesn't match
            // the extended path, only reconstruct from seed if abundance > CHI_WIN (100).
            // Low-abundance mismatches go to checktrf.
            if abundance > 100.0 {
                let mut newpath: Vec<usize> = Vec::with_capacity(real_nodes.len() + 2);
                newpath.push(source_id);
                newpath.extend_from_slice(&real_nodes);
                newpath.push(sink_id);
                use_path = newpath;
                if flow_pathpat.is_some() {
                    flow_pathpat =
                        Some(rebuild_flow_pathpat(&use_path, &seed_pattern, graph));
                }
                use_start = 1;
                use_last = use_path.len().saturating_sub(2);
                checkpath = false;
            } else {
                // Low-abundance hard-boundary mismatch: defer to checktrf
                trace_checktrf_enqueue(
                    idx,
                    "hard_boundary_low_abund",
                    abundance,
                    &transfrags[idx],
                    debug_nodes,
                    graph,
                    Some(&use_path[use_start..=use_last]),
                    Some(debug_back_ok),
                    Some(debug_fwd_ok),
                    None,
                    None,
                    None,
                );
                checktrf.push(idx);
                emit_debug_seed_proc!(
                    idx,
                    graph,
                    transfrags,
                    debug_nodes,
                    debug_back_ok,
                    debug_fwd_ok,
                    debug_flux,
                    debug_path_coords,
                    "hard_boundary"
                );
                record_outcome!(idx, SeedOutcome::Skipped("hard_boundary_low_abund"));
                continue;
            }
        }

        // Long-read: require LR witness across every consecutive PAIR of splice edges.
        // the original algorithm for each pair of consecutive splice edges in the
        // path, require at least one transfrag whose node list contains BOTH edges in order.
        // This prevents stitching novel junction combinations from separate reads.
        if checkpath && long_read_mode && use_last >= use_start && std::env::var_os("RUSTLE_ENABLE_WITNESS").is_some() {
            let mut splice_pos: Vec<usize> = Vec::new();
            for p in use_start..=use_last {
                if p == use_start {
                    continue;
                }
                let u = use_path[p - 1];
                let v = use_path[p];
                if is_splice_between(graph, u, v) {
                    splice_pos.push(p);
                }
            }
            let mut unwitnessed = false;
            // check consecutive PAIRS of splice edges, not individual edges.
            // For each pair (splice_k, splice_k+1), require one transfrag witnessing both.
            if splice_pos.len() >= 2 {
                for pair in splice_pos.windows(2) {
                    let (p_a, p_b) = (pair[0], pair[1]);
                    let la = use_path[p_a - 1];
                    let ra = use_path[p_a];
                    let lb = use_path[p_b - 1];
                    let rb = use_path[p_b];
                    if !has_lr_witness_two_splices(transfrags, la, ra, lb, rb) {
                        if trace_seed_active(idx) {
                            eprintln!(
                                "[TRACE_UNWIT_PAIR] idx={} edge1=({},{}) edge2=({},{}) coords=({}-{} -> {}-{}) ({}-{} -> {}-{})",
                                idx, la, ra, lb, rb,
                                graph.nodes.get(la).map(|n| n.start).unwrap_or(0),
                                graph.nodes.get(la).map(|n| n.end).unwrap_or(0),
                                graph.nodes.get(ra).map(|n| n.start).unwrap_or(0),
                                graph.nodes.get(ra).map(|n| n.end).unwrap_or(0),
                                graph.nodes.get(lb).map(|n| n.start).unwrap_or(0),
                                graph.nodes.get(lb).map(|n| n.end).unwrap_or(0),
                                graph.nodes.get(rb).map(|n| n.start).unwrap_or(0),
                                graph.nodes.get(rb).map(|n| n.end).unwrap_or(0),
                            );
                        }
                        unwitnessed = true;
                        break;
                    }
                }
            }
            if unwitnessed {
                // For hard-boundary seeds, keep compatible with defer to checktrf.
                // In max-sensitivity mode, also defer unwitnessed seeds so they can be
                // rediscovered by late passes without perturbing the greedy path extraction.
                if (thardstart && thardend) || config.max_sensitivity {
                    trace_checktrf_enqueue(
                        idx,
                        "unwitnessed",
                        abundance,
                        &transfrags[idx],
                        debug_nodes,
                        graph,
                        Some(&use_path[use_start..=use_last]),
                        Some(debug_back_ok),
                        Some(debug_fwd_ok),
                        None,
                        None,
                        None,
                    );
                    checktrf.push(idx);
                }
                // DEBUG: emit SEED_DECISION for unwitnessed branch
                emit_debug_seed_decision!(idx, 0, "UNWITNESSED");
                emit_debug_seed_proc!(
                    idx,
                    graph,
                    transfrags,
                    debug_nodes,
                    debug_back_ok,
                    debug_fwd_ok,
                    debug_flux,
                    debug_path_coords,
                    "unwitnessed"
                );
                record_outcome!(idx, SeedOutcome::UnwitnessedSplice);
                continue;
            }
        }

        // Capture path coordinates debug
        if debug_detail {
            let coords: Vec<String> = use_path[use_start..=use_last]
                .iter()
                .filter_map(|&n| {
                    graph
                        .nodes
                        .get(n)
                        .map(|nd| format!("{}-{}", nd.start, nd.end))
                })
                .collect();
            debug_path_coords = Some(coords.join(","));
        }
        if (long_read_mode || mixed_mode)
            && flow_pathpat.is_none()
            && std::env::var_os("RUSTLE_DISABLE_FLOW_PATHPAT_REBUILD").is_none()
        {
            flow_pathpat = Some(rebuild_flow_pathpat(&use_path, &seed_pattern, graph));
            if pathpat_trace {
                if let Some(ref fp) = flow_pathpat {
                    let (nb, eb, bits) = pathpat_summary(fp, graph.n_nodes);
                    eprintln!(
                        "parse_trflong: FLOW_PATHPAT idx={} source=rebuild node_bits={} edge_bits={} bits={}",
                        idx, nb, eb, bits
                    );
                }
            }
        }

        // Build exons from path[use_start..=use_last], merging contiguous nodes
        let mut exons: Vec<(u64, u64)> = Vec::new();
        let mut j = use_start;
        while j <= use_last {
            let nid = use_path[j];
            if nid >= graph.nodes.len() {
                j += 1;
                continue;
            }
            let node = &graph.nodes[nid];
            let nodestart = node.start;
            let mut nodeend = node.end;

            while j + 1 <= use_last {
                let next_nid = use_path[j + 1];
                if next_nid < graph.nodes.len()
                    && nodes_are_contiguous(graph, use_path[j], next_nid)
                {
                    j += 1;
                    nodeend = graph.nodes[next_nid].end;
                } else {
                    break;
                }
            }
            if nodeend > nodestart {
                exons.push((nodestart, nodeend));
            }
            j += 1;
        }

        if exons.is_empty() {
            emit_debug_seed_proc!(
                idx,
                graph,
                transfrags,
                debug_nodes,
                debug_back_ok,
                debug_fwd_ok,
                debug_flux,
                debug_path_coords,
                "too_short"
            );
            record_outcome!(idx, SeedOutcome::TooShort);
            continue;
        }

        // single-exon paths from multi-exon graphs are rejected
        // BEFORE flow extraction. the original implementation rejects all 2046 single-node paths
        // as "low_coverage" — they consume flow budget without producing valid
        // multi-exon predictions, starving downstream multi-exon path extraction.
        if long_read_mode && exons.len() == 1 && real_nodes.len() > 1 {
            // Single exon result from a multi-node seed: skip.
            record_outcome!(idx, SeedOutcome::Skipped("single_exon_from_multinode"));
            continue;
        }
        if long_read_mode && exons.len() == 1 && !transfrags[idx].guide {
            // Single-exon non-guide prediction: skip in LR mode (like 
            // These deplete flow budget without contributing to multi-exon assembly.
            transfrags[idx].abundance = 0.0;
            record_outcome!(idx, SeedOutcome::Skipped("single_exon_lr"));
            continue;
        }

        let mut source_ext_start: Option<u64> = None;
        let mut sink_ext_end: Option<u64> = None;
        if std::env::var_os("RUSTLE_SOURCE_CONTIG_EXTEND").is_some() {
            let start_nid = use_path[use_start];
            if let Some(ext_start) = source_contig_extension_start(graph, start_nid) {
                if !exons.is_empty() && ext_start < exons[0].0 {
                    exons[0].0 = ext_start;
                    source_ext_start = Some(ext_start);
                }
            }
            if std::env::var_os("RUSTLE_SINK_CONTIG_EXTEND").is_some() {
                let end_nid = use_path[use_last];
                if let Some(ext_end) = sink_contig_extension_end(graph, end_nid) {
                    if !exons.is_empty() {
                        let li = exons.len() - 1;
                        if ext_end > exons[li].1 {
                            exons[li].1 = ext_end;
                            sink_ext_end = Some(ext_end);
                        }
                    }
                }
            }
        }

        let length: u64 = exons.iter().map(|(s, e)| e - s).sum();
        if length < config.min_transcript_length {
            emit_debug_seed_proc!(
                idx,
                graph,
                transfrags,
                debug_nodes,
                debug_back_ok,
                debug_fwd_ok,
                debug_flux,
                debug_path_coords,
                "too_short"
            );
            record_outcome!(idx, SeedOutcome::TooShort);
            continue;
        }

        // Locus flow trace: check if any real node overlaps the trace range
        let trace_this = trace_locus.map_or(false, |(lo, hi)| {
            real_nodes.iter().any(|&nid| {
                graph
                    .nodes
                    .get(nid)
                    .map_or(false, |n| n.start <= hi && n.end >= lo)
            })
        });
        if trace_this {
            let first_rn = real_nodes.first().and_then(|&nid| graph.nodes.get(nid));
            let last_rn = real_nodes.last().and_then(|&nid| graph.nodes.get(nid));
            eprintln!("[TRACE_LOCUS] idx={} abund={:.1} nodes={} exons={} span={}-{} rnodes={}-{} complete={} seed={}",
                idx, abundance, real_nodes.len(), exons.len(),
                exons.first().map(|e| e.0).unwrap_or(0),
                exons.last().map(|e| e.1).unwrap_or(0),
                first_rn.map(|n| n.start).unwrap_or(0),
                last_rn.map(|n| n.end).unwrap_or(0),
                thardstart || thardend,
                transfrags[idx].trflong_seed);
        }

        // standard optional flow assignment: when enabled, push max-flow with subtraction
        // so overlapping paths compete for shared evidence.
        let mut out_exoncov: Vec<f64> = vec![0.0; exons.len()];
        let flow_flux;
        // True when (restore + nodecov-limited max_fl) is used for this path.
        // When true, nodecov depletion loop does NOT cap nflux to remaining nodecov.
        let long_read_nodecov_mode = false;

        // Debug: trace flow ON vs OFF differences
        let debug_bundle = config.debug_bundle.as_ref();
        let debug_flow = config.verbose && debug_bundle.map_or(false, |b| bundle_id.contains(b));
        // Per-transcript coverage debug for any transcript overlapping the trace locus
        // (env RUSTLE_COV_DEBUG=1 + optional RUSTLE_TRACE_LOCUS=start-end).
        let debug_cov =
            long_read_mode && std::env::var("RUSTLE_COV_DEBUG").is_ok() && trace_this;
        if debug_flow {
            eprintln!(
                "  [FLOW_DEBUG] idx={} guide={} longread={} abund={:.4} srabund={:.4} nodes={}",
                idx,
                transfrags[idx].guide,
                transfrags[idx].longread,
                transfrags[idx].abundance,
                transfrags[idx].srabund,
                transfrags[idx].node_ids.len()
            );
        }

        let coverage = {
            let mut flow_used: HashSet<usize> = Default::default();
            let mut nodeflux_is_proportion = false;
            let (flux, nodeflux) = if guided_mode
                && transfrags[idx].guide
                && long_read_mode
                && !mixed_mode
            {
                let guide_pattern = transfrags[idx].pattern.clone();
                let mut guide_abundance = push_guide_max_flow(&use_path, transfrags, graph, true);
                if guide_abundance <= 0.0 {
                    if debug_ek {
                        eprintln!("[RUST_GUIDE_ZERO] idx={} path_len={} abund=0", idx, use_path.len());
                    }
                    (0.0, Vec::new())
                } else {
                    // guides_pushmaxflow-style EM loop (nEM=10) with CNodeGuide-like trcount/gcount priors.
                    let mut cur_node_prior: HashMap<usize, f64> = Default::default();
                    let prev_refs_base: Vec<GuideFlowPriorRef<'_>> = previous_guides
                        .iter()
                        .map(|g| GuideFlowPriorRef {
                            pattern: &g.pattern,
                            abundance: g.abundance,
                            node_prior: &g.node_prior,
                        })
                        .collect();
                    for _ in 0..10 {
                        let (em_prior, em_cov) = cnodeguide_em_prior(
                            &use_path,
                            &guide_pattern,
                            guide_abundance,
                            &previous_guides,
                            transfrags,
                            graph,
                        );
                        if !em_prior.is_empty() {
                            cur_node_prior = em_prior;
                        } else if cur_node_prior.is_empty() {
                            cur_node_prior = build_cnodeguide_node_prior(
                                &use_path,
                                &guide_pattern,
                                &previous_guides,
                                transfrags,
                                graph,
                            );
                            if cur_node_prior.is_empty() {
                                cur_node_prior = build_strict_guide_node_prior(
                                    &use_path,
                                    &guide_pattern,
                                    &previous_guides,
                                    transfrags,
                                    graph,
                                );
                            }
                        }
                        if em_cov > 0.0 {
                            guide_abundance = em_cov;
                        }
                        let mut tf_clone = transfrags.to_vec();
                        let (next_flux, next_nodeflux) = guide_push_flow(
                            &use_path,
                            &guide_pattern,
                            guide_abundance,
                            &cur_node_prior,
                            &prev_refs_base,
                            &mut tf_clone,
                            graph,
                        );
                        if next_flux <= 0.0 {
                            break;
                        }
                        cur_node_prior = build_guide_node_prior(&use_path, &next_nodeflux, graph);
                        if (next_flux - guide_abundance).abs() < 1e-3 {
                            guide_abundance = next_flux;
                            break;
                        }
                        guide_abundance = next_flux;
                    }
                    let (gflux, gnodeflux) = guide_push_flow(
                        &use_path,
                        &guide_pattern,
                        guide_abundance,
                        &cur_node_prior,
                        &prev_refs_base,
                        transfrags,
                        graph,
                    );
                    if debug_ek {
                        eprintln!("[RUST_GUIDE_PUSH] idx={} gflux={:.4}", idx, gflux);
                    }
                    drop(prev_refs_base);
                    if gflux > 0.0 {
                        previous_guides.push(GuideFlowState {
                            pattern: guide_pattern,
                            abundance: gflux,
                            node_prior: build_guide_node_prior(&use_path, &gnodeflux, graph),
                        });
                    }
                    (gflux, gnodeflux)
                }
            } else if mixed_mode {
                if transfrags[idx].longread {
                    let (lf, lfpath, lused) = long_max_flow_seeded_with_used_pathpat(
                        &use_path,
                        transfrags,
                        graph,
                        false,
                        Some(idx),
                        flow_pathpat.as_ref(),
                    );
                    flow_used = lused.ones().collect();
                    let (sf, sfpath, sfull) =
                        push_max_flow_seeded_full(&use_path, transfrags, graph, true, Some(idx));
                    if sf > lf {
                        nodeflux_is_proportion = true;
                        let _ = sfull;
                        (sf, sfpath)
                    } else {
                        (lf, lfpath)
                    }
                } else {
                    let (sf, sfpath, sfull) =
                        push_max_flow_seeded_full(&use_path, transfrags, graph, true, Some(idx));
                    nodeflux_is_proportion = true;
                    let _ = sfull;
                    (sf, sfpath)
                }
            } else if long_read_mode {
                // Strict: one long_max_flow call per seed path.
                // No Rust-specific retry/capping branch.
                // long_max_flow: pass no_subtract=false so transfrag abundances
                // are depleted by used flow while nodecov is also depleted downstream.
                let (lf, lfpath, lused) = long_max_flow_seeded_with_used_pathpat(
                    &use_path,
                    transfrags,
                    graph,
                    false,
                    Some(idx),
                    flow_pathpat.as_ref(),
                );
                flow_used = lused.ones().collect();
                if debug_cov {
                    eprintln!("  [COV_DEBUG] long_flux={:.4} nodeflux={:?}", lf, &lfpath);
                }
                (lf, lfpath)
            } else {
                let (sf, sfpath, sfull) =
                    push_max_flow_seeded_full(&use_path, transfrags, graph, false, Some(idx));
                nodeflux_is_proportion = true;
                let _ = sfull;
                (sf, sfpath)
            };
            flow_flux = flux;
            debug_flux = flux;

            if depletion_diag {
                let ab_after = transfrags[idx].abundance;
                eprintln!(
                    "SEED_FLOW idx={} ab_before={:.6} flux={:.6} ab_after={:.6} path_len={} is_guide={}",
                    idx, abundance, flux, ab_after, use_path.len(), transfrags[idx].guide
                );
                if abundance > 0.0 && ab_after <= 0.0 && flux > 0.0 {
                    eprintln!("  SEED_FULLY_DEPLETED idx={} by_flux={:.6}", idx, flux);
                }
            }

            if debug_ek {
                eprintln!(
                    "[RUST_FLOW] idx={} flux={:.6} path_len={}",
                    idx,
                    flux,
                    use_path.len()
                );
                for nfi in use_start..=use_last {
                    let nid = use_path[nfi];
                    if nid == source_id || nid == sink_id {
                        continue;
                    }
                    let nc = local_nodecov.get(nid).copied().unwrap_or(0.0);
                    let nfx = if nfi >= use_start && (nfi - use_start) < nodeflux.len() {
                        nodeflux[nfi - use_start]
                    } else {
                        0.0
                    };
                    eprintln!("  node[{}]={} nflux={:.6} ncov={:.6}", nfi, nid, nfx, nc);
                }
            }

            if debug_flow && flux > 0.0 {
                eprintln!(
                    "  [FLOW_DEBUG] idx={} flux={:.4} nodes={} tocheck={}",
                    idx,
                    flux,
                    transfrags[idx].node_ids.len(),
                    tocheck
                );
            }

            // mixedMode parse_trflong start/end refinement using istranscript + longstart/longend.
            if mixed_mode && transfrags[idx].longread && use_start <= use_last && !exons.is_empty()
            {
                let start_nid = use_path[use_start];
                let end_nid = use_path[use_last];
                if let (Some(snode), Some(enode)) =
                    (graph.nodes.get(start_nid), graph.nodes.get(end_nid))
                {
                    let mut startpoint = snode.end;
                    let mut endpoint = enode.start;
                    if let Some(jnode) = graph.nodes.get(start_nid) {
                        for &u in &jnode.trf_ids {
                            if !flow_used.contains(&u) || u >= transfrags.len() {
                                continue;
                            }
                            let tfu = &transfrags[u];
                            if tfu.longread
                                && tfu.node_ids.first().copied().unwrap_or(source_id) != source_id
                                && tfu.node_ids.last().copied().unwrap_or(sink_id) != sink_id
                                && tfu.node_ids.first() == Some(&start_nid)
                                && tfu.longstart > 0
                                && tfu.longstart < startpoint
                            {
                                startpoint = tfu.longstart;
                            }
                        }
                    }
                    if let Some(jnode) = graph.nodes.get(end_nid) {
                        for &u in &jnode.trf_ids {
                            if !flow_used.contains(&u) || u >= transfrags.len() {
                                continue;
                            }
                            let tfu = &transfrags[u];
                            if tfu.longread
                                && tfu.node_ids.first().copied().unwrap_or(source_id) != source_id
                                && tfu.node_ids.last().copied().unwrap_or(sink_id) != sink_id
                                && tfu.node_ids.last() == Some(&end_nid)
                                && tfu.longend > 0
                                && tfu.longend > endpoint
                            {
                                endpoint = tfu.longend;
                            }
                        }
                    }
                    if startpoint == snode.end {
                        startpoint = snode.start;
                    }
                    if endpoint == enode.start {
                        endpoint = enode.end;
                    }
                    let (es, ee) = exons[0];
                    let mut left = es.max(startpoint).min(ee);
                    if source_ext_start.is_some()
                        && std::env::var_os("RUSTLE_SOURCE_CONTIG_EXTEND").is_some()
                        && startpoint == snode.start
                    {
                        // Keep source-contiguous extension when start refinement did not
                        // establish a stricter in-node start boundary.
                        left = es;
                    }
                    exons[0].0 = left;
                    let li = exons.len() - 1;
                    let (ls, le) = exons[li];
                    let mut right = le.min(endpoint).max(ls);
                    if sink_ext_end.is_some()
                        && std::env::var_os("RUSTLE_SINK_CONTIG_EXTEND").is_some()
                        && endpoint == enode.end
                    {
                        // Keep sink-contiguous extension when end refinement did not establish
                        // a stricter in-node end boundary.
                        right = le;
                    }
                    exons[li].1 = right;
                }
            }
            if !nodeflux.is_empty() {
                let (ecov, cov_from_nodes) = accumulate_exon_cov_from_path_usage(
                    &use_path,
                    use_start,
                    use_last,
                    graph,
                    &exons,
                    &local_nodecov,
                    &nodeflux,
                    nodeflux_is_proportion,
                    if nodeflux_is_proportion {
                        None
                    } else {
                        Some(&local_noderate)
                    },
                );
                if ecov.len() == exons.len() {
                    out_exoncov = ecov;
                }
                if nodeflux_is_proportion {
                    // store_transcript short-flow: nodecov *= (1-nodeflux).
                    let mut k = 0usize;
                    for p in use_start..=use_last {
                        let nid = use_path[p];
                        if nid == source_id || nid == sink_id || nid >= local_nodecov.len() {
                            continue;
                        }
                        if k >= nodeflux.len() {
                            break;
                        }
                        let nflux = nodeflux[k].clamp(0.0, 1.0);
                        local_nodecov[nid] = (local_nodecov[nid] * (1.0 - nflux)).max(0.0);
                        k += 1;
                    }
                    if cov_from_nodes > 0.0 {
                        cov_from_nodes
                    } else if flux > 0.0 {
                        flux
                    } else {
                        abundance.max(transfrags[idx].srabund)
                    }
                } else if long_read_mode {
                    // parse_trflong coverage formula
                    //   ecov = nodeflux_abs[j] * noderate[path[j]]
                    //   nodecov[path[j]] -= nodeflux_abs[j]
                    // long_max_flow_seeded returns ABSOLUTE nodecapacity values (same units as nodecov).
                    // So: ecov = abs_flow * (bpcov / nodecov) = per-base coverage contribution.
                    let mut cov_total = 0.0f64;
                    let mut lr_exoncov = vec![0.0f64; exons.len()];
                    let mut k = 0usize;
                    for p in use_start..=use_last {
                        let nid = use_path[p];
                        if nid == source_id || nid == sink_id || nid >= local_nodecov.len() {
                            continue;
                        }
                        if k >= nodeflux.len() {
                            break;
                        }
                        let nc_before = local_nodecov[nid].max(0.0);
                        // update_capacity skips the last node, so nodeflux[last]=0.
                        // Match exactly: use nodeflux[k] directly (0 for last node).
                        let raw_nflux = nodeflux[k];
                        // parse_trflong (10381, 10393): cap nodeflux by nodecov, then subtract.
                        let nflux = raw_nflux.min(local_nodecov[nid].max(0.0));
                        local_nodecov[nid] = (local_nodecov[nid] - nflux).max(0.0);
                        if local_nodecov[nid] < EPS {
                            local_nodecov[nid] = 0.0;
                        }
                        let ecov = nflux * local_noderate[nid];
                        if debug_cov {
                            let bpcov = graph.nodes.get(nid).map(|n| n.coverage).unwrap_or(0.0);
                            eprintln!("  [COV_DEBUG] node[{}] {}..{} nc_before={:.4} nodeflux_raw={:.4} nflux={:.4} rate={:.6} bpcov={:.4} contrib={:.6}",
                                nid, graph.nodes[nid].start, graph.nodes[nid].end,
                                nc_before, nodeflux[k], nflux, local_noderate[nid], bpcov, ecov);
                        }
                        cov_total += ecov;
                        // Distribute per-node ecov to overlapping exons
                        // accumulates raw ecov per exon, then divides by exon_len at line 10956.
                        if let Some(node) = graph.nodes.get(nid) {
                            for (ei, &(es, ee)) in exons.iter().enumerate() {
                                if overlaps_half_open(es, ee, node.start, node.end) {
                                    lr_exoncov[ei] += ecov;
                                }
                            }
                        }
                        k += 1;
                    }
                    // exoncov[i] /= exons[i].len()
                    for (ei, &(es, ee)) in exons.iter().enumerate() {
                        let elen = len_half_open(es, ee) as f64;
                        if elen > 0.0 {
                            lr_exoncov[ei] /= elen;
                        }
                    }
                    if lr_exoncov.iter().any(|v| *v > 0.0) {
                        out_exoncov = lr_exoncov;
                    }
                    if debug_cov {
                        let computed = if cov_total > 0.0 && length > 0 {
                            cov_total / length as f64
                        } else {
                            0.0
                        };
                        eprintln!("[COV_DEBUG] RESULT cov_total={:.6} length={} cov={:.6} long_read_nodecov_mode={}",
                            cov_total, length, computed, long_read_nodecov_mode);
                    }
                    let mut computed_cov = if cov_total > 0.0 && length > 0 {
                        cov_total / length as f64 // per-base: pred->cov /= abs(tlen) at 
                    } else if flux > 0.0 {
                        flux
                    } else {
                        abundance.max(transfrags[idx].srabund)
                    };
                    if transfrags[idx].trflong_seed
                        && !transfrags[idx].guide
                        && (thardstart || thardend)
                        && exons.len() >= 20
                    {
                        // Some complete long-read seed paths retain real long-read abundance but
                        // collapse to near-zero per-base coverage after node-rate normalization.
                        // Keep a modest seed-derived floor so they can still compete downstream.
                        let seed_support = abundance.max(flux).max(transfrags[idx].srabund);
                        let seed_floor = (4.0 * seed_support) / (exons.len().max(4) as f64);
                        if seed_floor > computed_cov {
                            computed_cov = seed_floor;
                        }
                    }
                    computed_cov
                } else if cov_from_nodes > 0.0 {
                    cov_from_nodes
                } else if flux > 0.0 {
                    flux
                } else {
                    abundance.max(transfrags[idx].srabund)
                }
            } else if flux > 0.0 {
                flux
            } else {
                abundance.max(transfrags[idx].srabund)
            }
        };
        
        if trace_this {
            eprintln!(
                "[TRACE_LOCUS]   → flux={:.4} cov={:.4} nodecov_mode={}",
                flow_flux, coverage, long_read_nodecov_mode
            );
        }
        if flow_flux > 0.0 || transfrags[idx].guide {
            tocheck = false;
        }

        // Debug: trace zero-flux decision
        if debug_flow && transfrags[idx].abundance > 0.0 {
            eprintln!(
                "  [FLOW_DEBUG] idx={} final_flux={:.4} tocheck={} guide={}",
                idx, flow_flux, tocheck, transfrags[idx].guide
            );
        }

        if long_read_mode && !mixed_mode && !transfrags[idx].guide && flow_flux <= 0.0 {
            // parse_trflong: zero-flux candidates are deferred to checktrf
            // without consuming their seed abundance.
            transfrags[idx].abundance = abundance;
            if tocheck {
                trace_checktrf_enqueue(
                    idx,
                    "zero_flux",
                    abundance,
                    &transfrags[idx],
                    debug_nodes,
                    graph,
                    Some(&use_path[use_start..=use_last]),
                    Some(debug_back_ok),
                    Some(debug_fwd_ok),
                    Some(flow_flux),
                    None,
                    Some(tocheck),
                );
                checktrf.push(idx);
                zero_flux_set.insert(idx);
                zero_flux_candidates += 1;
                if trace_this {
                    eprintln!("[TRACE_LOCUS]   → ZERO_FLUX checktrf idx={}", idx);
                }
                if debug_flow {
                    eprintln!("  [FLOW_DEBUG] idx={} ADDED TO CHECKTRF (zero-flux)", idx);
                }
                if audit_zero_flux {
                    let first = real_nodes.first().copied().unwrap_or(source_id);
                    let last = real_nodes.last().copied().unwrap_or(sink_id);
                    let c0 = graph.nodes.get(first).map(|n| n.start).unwrap_or(0);
                    let c1 = graph.nodes.get(last).map(|n| n.end).unwrap_or(0);
                    eprintln!(
                        "    [audit_zero_flux] gate idx={} nodes={} span={}-{} abund={:.4} polyS={} polyE={}",
                        idx,
                        real_nodes.len(),
                        c0,
                        c1,
                        abundance,
                        poly_start,
                        poly_end
                    );
                }
            }
            emit_debug_seed_proc!(
                idx,
                graph,
                transfrags,
                debug_nodes,
                debug_back_ok,
                debug_fwd_ok,
                debug_flux,
                debug_path_coords,
                "zero_flux"
            );
            record_outcome!(idx, SeedOutcome::ZeroFlux);
            continue;
        }
        // parse_trflong (line 10437):
        // Single-node transfrags are zeroed BEFORE the store check so that a single-node
        // transfrag that fails the coverage/length gate is still excluded from checktrf/post-pass.
        // Multi-node long-read transfrags are NOT zeroed here — their abundance is kept for
        // later checktrf and post-pass redistribution.
        if !long_read_mode || real_nodes.len() == 1 {
            transfrags[idx].abundance = 0.0;
            if mixed_mode {
                transfrags[idx].srabund = 0.0;
            }
        }

        // parse_trflong stores predictions with cov>epsilon (plus len checks).
        // The coverage gate uses epsilon for BOTH long-read and short-read modes.
        // guide predictions always stored: if(t || (cov && len>=mintranscriptlen))
        // store if: g || (!eonly && len>=mintranscriptlen && cov>epsilon)
        let is_guide_pred = transfrags[idx].guide || transfrags[idx].guide_tid.is_some();
        // uses epsilon (1e-6) for the coverage store gate. However, in bundle
        // graph mode, the unified graphs produce more thin flow paths than the original algorithm's
        // equivalent bundles (335 vs 64 seeds with flux < 1.0). Raise the gate to
        // filter ultra-low-coverage predictions that are mostly FPs.
        let min_cov_gate = if std::env::var_os("RUSTLE_BUNDLE_GRAPH").is_some() {
            0.5  // 0.5x coverage minimum in bundle mode (loses 3 TPs, removes ~22 FPs)
        } else {
            EPS
        };
        
        
        if coverage < min_cov_gate && !is_guide_pred {
            if trace_this {
                eprintln!(
                    "[TRACE_LOCUS]   → LOW_COV cov={:.4} < gate={:.4} → checktrf={}",
                    coverage, min_cov_gate, tocheck
                );
            }
            if debug_flow {
                eprintln!(
                    "  [FLOW_DEBUG] idx={} COVERAGE {} < MIN_GATE {} -> tocheck={}",
                    idx, coverage, min_cov_gate, tocheck
                );
            }
            if long_read_mode && tocheck {
                let source_sink_linked = if real_nodes.is_empty() {
                    false
                } else {
                    let firstn = real_nodes[0];
                    let lastn = *real_nodes.last().unwrap_or(&firstn);
                    let source_link = graph
                        .nodes
                        .get(firstn)
                        .map(|n| n.parents.contains(source_id))
                        .unwrap_or(false);
                    let sink_link = graph
                        .nodes
                        .get(lastn)
                        .map(|n| n.children.contains(sink_id))
                        .unwrap_or(false);
                    source_link && sink_link
                };
                if !mixed_mode || (!guided_mode || transfrags[idx].guide || source_sink_linked) {
                    trace_checktrf_enqueue(
                        idx,
                        "low_coverage",
                        abundance,
                        &transfrags[idx],
                        debug_nodes,
                        graph,
                        Some(&use_path[use_start..=use_last]),
                        Some(debug_back_ok),
                        Some(debug_fwd_ok),
                        Some(flow_flux),
                        Some(coverage),
                        Some(tocheck),
                    );
                    checktrf.push(idx);
                }
            }
            // DEBUG: emit SEED_DECISION for low_coverage branch
            emit_debug_seed_decision!(idx, 1, "DIRECT_LOW_COV", debug_flux, coverage, debug_nodes.len());
            // DEBUG: emit SEED_DECISION for low_coverage branch
            emit_debug_seed_decision!(idx, 1, "DIRECT_LOW_COV", debug_flux, coverage, 0);
            emit_debug_seed_proc!(
                idx,
                graph,
                transfrags,
                debug_nodes,
                debug_back_ok,
                debug_fwd_ok,
                debug_flux,
                debug_path_coords,
                "low_coverage"
            );
            // DEBUG: emit SEED_DECISION for low_coverage branch
            emit_debug_seed_decision!(idx, 1, "DIRECT_LOW_COV", debug_flux, coverage, 0);
            record_outcome!(idx, SeedOutcome::LowCoverage(coverage));
            continue;
        }

        // eonly stores only guide-matched predictions.
        // if(g || (!eonly && len>=mintranscriptlen && cov>epsilon)) — where g = guide flag.
        if config.eonly && !is_guide_pred {
            emit_debug_seed_decision!(idx, 1, "eonly_skip");
            emit_debug_seed_proc!(
                idx,
                graph,
                transfrags,
                debug_nodes,
                debug_back_ok,
                debug_fwd_ok,
                debug_flux,
                debug_path_coords,
                "eonly_skip"
            );
            record_outcome!(idx, SeedOutcome::EonlyNonGuide);
            continue;
        }

        // Emit DEBUG_SEED_DECISION for stored transcripts
        emit_debug_seed_decision!(idx, 1, "STORED", debug_flux, coverage, exons.len());

        if trace_this {
            let first_exon = exons.first().copied().unwrap_or((0, 0));
            let last_exon = exons.last().copied().unwrap_or((0, 0));
            eprintln!(
                "[TRACE_LOCUS]   → STORED idx={} cov={:.4} abund_after={:.4} exons={} first={}-{} last={}-{}",
                idx,
                coverage,
                transfrags[idx].abundance,
                exons.len(),
                first_exon.0,
                first_exon.1,
                last_exon.0,
                last_exon.1
            );
        }
        if trace_chain_active {
            if let Some(tc) = trace_chain {
                let introns = intron_chain_from_exons(&exons);
                let tol = trace_intron_tol();
                if intron_chains_match_tol(&introns, tc, tol) {
                    _trace_chain_emitted = true;
                    eprintln!(
                        "[TRACE_CHAIN] emit idx={} cov={:.4} exons={} bundle={}",
                        idx,
                        coverage,
                        exons.len(),
                        bundle_id
                    );
                } else if !tc.is_empty() && intron_chain_contains_all_tol(&introns, tc, tol) {
                    eprintln!(
                        "[TRACE_CHAIN] superset idx={} introns={} bundle={}",
                        idx,
                        introns.len(),
                        bundle_id
                    );
                }
            }
        }
        if pathpat_trace {
            let first_start = exons.first().map(|e| e.0).unwrap_or(0);
            let last_end = exons.last().map(|e| e.1).unwrap_or(0);
            eprintln!(
                "parse_trflong: PRED_STORED idx={} pred={} cov={:.4} exons={} span={}-{} longcov={:.4}",
                idx, out.len(), coverage, exons.len(), first_start, last_end, abundance
            );
        }
        // DEBUG: emit SEED_DECISION for stored transcripts
        emit_debug_seed_decision!(idx, 1, "STORED", flow_flux, coverage, exons.len());
        out.push(Transcript {
            chrom: bundle_chrom.to_string(),
            strand: bundle_strand,
            exons: exons.clone(),
            coverage,
            exon_cov: if out_exoncov.len() == exons.len() && out_exoncov.iter().any(|v| *v > 0.0) {
                out_exoncov.clone()
            } else {
                vec![coverage; exons.len()]
            },
            tpm: 0.0,
            fpkm: 0.0,
            source: gtf_source_long_flow(&transfrags[idx].guide_tid),
            is_longread: long_read_mode,
            longcov: read_count_snapshot, // pre-depletion read mass: stable, comparable to the original algorithm's longcov
            bpcov_cov: 0.0,
            transcript_id: None,
            gene_id: None,
            ref_transcript_id: None,
            ref_gene_id: None,
            hardstart: thardstart,
            hardend: thardend,
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None,
        });
        if debug_flow {
            let exons_str = exons
                .iter()
                .map(|(s, e)| format!("{}-{}", s, e))
                .collect::<Vec<_>>()
                .join(",");
            eprintln!(
                "  [FLOW_DEBUG] idx={} STORED exons={} cov={:.4}",
                idx, exons_str, coverage
            );
        }
        if plumb_debug {
            let nexons = exons.len();
            eprintln!(
                "[PLUMB] tx exons={} cov={:.6} longcov={:.6} flux={:.6}",
                nexons, coverage, abundance, flow_flux
            );
        }
        let out_idx = out.len() - 1;
        record_outcome!(idx, SeedOutcome::Stored(out_idx));
        if debug_ek && !exons.is_empty() {
            let elen: u64 = exons.iter().map(|(s, e)| e - s).sum();
            // Compute intron chain for tracing
            let introns: Vec<(u64, u64)> = exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
            let intron_str = introns.iter().map(|(d, a)| format!("{}-{}", d, a)).collect::<Vec<_>>().join(",");
            eprintln!(
                "[RUST_PRED] idx={} pred_idx={} cov_perbase={:.6} len={} start={} end={} nexons={} introns=[{}]",
                idx,
                out_idx,
                coverage,
                elen,
                exons[0].0,
                exons.last().unwrap().1,
                exons.len(),
                intron_str
            );
            eprint!("  exons:");
            for (es, ee) in &exons {
                eprint!(" {}-{}", es, ee);
            }
            eprintln!();
        }
        let kept_inner: Vec<usize> = use_path[use_start..=use_last].to_vec();
        if let Some(ti) = trace_intron {
            if path_has_intron(&use_path, graph, ti) {
                eprintln!(
                    "[TRACE_INTRON_TF] stage=extracted idx={} intron={}-{} cov={:.4} nodes={} exons={} bundle={}",
                    idx, ti.0, ti.1, coverage, use_path.len(), exons.len(), bundle_id
                );
            }
        }
        if debug_detail {
            emit_debug_seed_proc!(
                idx,
                graph,
                transfrags,
                debug_nodes,
                debug_back_ok,
                debug_fwd_ok,
                debug_flux,
                debug_path_coords,
                "stored"
            );
            eprint!(
                "DEBUG_PRED idx={} cov={:.4} nexons={} exons=",
                idx,
                coverage,
                exons.len()
            );
            for (i, (es, ee)) in exons.iter().enumerate() {
                if i > 0 {
                    eprint!(",");
                }
                eprint!("{}-{}", es, ee);
            }
            eprintln!();
        }
        kept_paths.push((kept_inner, coverage, transfrags[idx].guide, out_idx));
    }

    // checktrf rescue pass (parse_trflong follow-up):
    // 1) try to reassign skipped transfrag support to best kept paths
    // 2) if unmatched, store direct complete transfrag as an independent transcript
    let checktrf_enabled = std::env::var_os("RUSTLE_DISABLE_CHECKTRF").is_none();
    if long_read_mode && !checktrf.is_empty() && checktrf_enabled {
        // Dedup while preserving insertion order.
        let mut seen: HashSet<usize> = Default::default();
        let mut uniq_check = Vec::new();
        for t in checktrf {
            if seen.insert(t) {
                uniq_check.push(t);
            }
        }
        let has_guide_kept = kept_paths.iter().any(|(_, _, g, _)| *g);
        if debug_ek {
            eprintln!(
                "[RUST_CHECKTRF] count={} keeptrf={}",
                uniq_check.len(),
                kept_paths.len()
            );
        }
        for t in uniq_check {
            if t >= transfrags.len() {
                if zero_flux_set.contains(&t) {
                    zero_flux_dropped_empty += 1;
                }
                record_outcome!(t, SeedOutcome::ChecktrfRescueFail);
                continue;
            }
            if transfrags[t].node_ids.is_empty() {
                if zero_flux_set.contains(&t) {
                    zero_flux_dropped_empty += 1;
                }
                record_outcome!(t, SeedOutcome::ChecktrfRescueFail);
                continue;
            }
            // process checktrf only when guide OR abundance>=readthr.
            // The `(!eonly || guide)` logic is applied later for independent rescue.
            if !checktrf_passes_abundance_gate(&transfrags[t], config.readthr) {
                if zero_flux_set.contains(&t) {
                    zero_flux_dropped_lowabund += 1;
                }
                if let Some((lo, hi)) = trace_locus {
                    let in_range = transfrags[t].node_ids.iter().any(|&nid| {
                        graph
                            .nodes
                            .get(nid)
                            .map_or(false, |n| n.start <= hi && n.end >= lo)
                    });
                    if in_range {
                        eprintln!(
                            "[TRACE_CHECKTRF] t={} SKIPPED: abund={:.4} read_count={:.4} < readthr={:.4} guide={} nodes={:?}",
                            t,
                            transfrags[t].abundance,
                            transfrags[t].read_count,
                            config.readthr,
                            transfrags[t].guide,
                            &transfrags[t].node_ids[..transfrags[t].node_ids.len().min(12)]
                        );
                    }
                }
                record_outcome!(t, SeedOutcome::ChecktrfReadthr);
                continue;
            }
            if debug_ek {
                eprintln!(
                    "[RUST_CHECKTRF] t={} abund={:.6} guide={} nodes={}",
                    t,
                    transfrags[t].abundance,
                    transfrags[t].guide,
                    transfrags[t].node_ids.len()
                );
            }
            let tf_nodes: Vec<usize> = transfrags[t]
                .node_ids
                .iter()
                .copied()
                .filter(|&n| n != source_id && n != sink_id)
                .collect();
            if tf_nodes.is_empty() {
                continue;
            }
            trace_checktrf_seed_snapshot(t, &transfrags[t], &tf_nodes, &kept_paths, graph, &out);

            // checktrf (get_trf_long:11133-11297):
            // For multi-node longread: first try to redistribute to best-matched kept paths.
            // If no match found: fall through to independent rescue (get_trf_long:11180
            //   `else if(!eonly || guide)` — in non-eonly mode always rescues complete TFs).
            // Shortread / <=1-node: also attempt independent rescue (no prior redistribution).
            let is_shortread_tf = !transfrags[t].longread;
            if !is_shortread_tf && tf_nodes.len() > 1 {
                let mut tmatch: Vec<usize> = Vec::new();
                let mut abundancesum = 0.0;
                if has_guide_kept {
                    let (m, s) = best_trf_match(&tf_nodes, &kept_paths, graph, true, Some(t));
                    tmatch = m;
                    abundancesum = s;
                }
                if tmatch.is_empty() {
                    let (m, s) = best_trf_match(&tf_nodes, &kept_paths, graph, false, Some(t));
                    tmatch = m;
                    abundancesum = s;
                }
                if !tmatch.is_empty() {
                    // best_trf_match uses node-intersection to find overlapping kept paths.
                    // Two transcripts with distinct splice sites use different graph nodes, so
                    // best_trf_match naturally prevents cross-isoform redistribution.
                    // The secondary intron_chains_equal_tol filter was redundant and blocked
                    // valid redistributions (required same intron count, so different-length
                    // transfrag chains always cleared tmatch → independent rescue).
                    // the original algorithm redistributes based on best_trf_match result alone.
                    if tmatch.is_empty() {
                        // Treat as "no match" and fall through to independent rescue below.
                    } else {
                    {
                        if let Some((lo, hi)) = trace_locus {
                            let in_range = transfrags[t].node_ids.iter().any(|&nid| {
                                graph
                                    .nodes
                                    .get(nid)
                                    .map_or(false, |n| n.start <= hi && n.end >= lo)
                            });
                            if in_range {
                                eprintln!("[TRACE_CHECKTRF] t={} MATCHED to kept_paths={:?} abundsum={:.4} abund={:.4}",
                                t, &tmatch, abundancesum, transfrags[t].abundance);
                            }
                        }
                        if zero_flux_set.contains(&t) {
                            zero_flux_rescued += 1;
                            if audit_zero_flux {
                                eprintln!(
                                "    [audit_zero_flux] rescued idx={} matched_paths={} abundsum={:.4}",
                                t,
                                tmatch.len(),
                                abundancesum
                            );
                            }
                        }
                        if debug_detail {
                            eprintln!("DEBUG_CHK t={} abund={:.4} outcome=matched matched={:?} abundsum={:.4}",
                            t, transfrags[t].abundance, &tmatch, abundancesum);
                        }
                        trace_checktrf_match_snapshot(
                            t,
                            &tmatch,
                            abundancesum,
                            &kept_paths,
                            graph,
                            &out,
                        );
                        let updated = redistribute_transfrag_to_matches(
                            &transfrags[t],
                            &tmatch,
                            abundancesum,
                            &kept_paths,
                            graph,
                            &mut out,
                        );
                        if updated {
                            transfrags[t].abundance = 0.0;
                        }
                        // DEBUG: emit DEBUG_CHK for matched outcome
                        emit_debug_chk!(t, transfrags[t].abundance, "matched", tmatch.len(), abundancesum);
                        record_outcome!(t, SeedOutcome::ChecktrfRedistributed);
                        continue; // matched: done
                    }
                    }
                }
                // No match found: fall through to independent rescue below.
                if let Some((lo, hi)) = trace_locus {
                    let in_range = transfrags[t].node_ids.iter().any(|&nid| {
                        graph
                            .nodes
                            .get(nid)
                            .map_or(false, |n| n.start <= hi && n.end >= lo)
                    });
                    if in_range {
                        eprintln!("[TRACE_CHECKTRF] t={} NO MATCH → trying independent rescue abund={:.4}", t, transfrags[t].abundance);
                    }
                }
            }
            {
                // Independent rescue: store complete transfrag as its own prediction.
                // Applies to: shortread/<=1-node, AND multi-node longread with no kept-path match.
                // `else if(!eonly || guide)` — in non-eonly mode always rescues.
                // Non-contiguous jumps must have explicit transfrag edge pattern support.

                // in eonly mode, only guide transfrags are independently rescued.
                if config.eonly && !transfrags[t].guide {
                    record_outcome!(t, SeedOutcome::ChecktrfEonlySkip);
                    continue;
                }

                let rescue_nodes = &tf_nodes;
                let mut complete = true;
                for w in rescue_nodes.windows(2) {
                    let a = w[0];
                    let b = w[1];
                    let contiguous = nodes_are_contiguous(graph, a, b);
                    if contiguous {
                        continue;
                    }
                    let has_edge = graph
                        .nodes
                        .get(a)
                        .map(|n| n.children.contains(b))
                        .unwrap_or(false);
                    let has_pattern_edge = graph
                        .edge_bit_index(a, b)
                        .map(|eid| transfrags[t].pattern.get_bit(eid))
                        .unwrap_or(false);
                    if !has_edge || !has_pattern_edge {
                        complete = false;
                        break;
                    }
                }
                if !complete {
                    // else-branch candidates are consumed
                    // regardless of rescue success.
                    transfrags[t].abundance = 0.0;
                    record_outcome!(t, SeedOutcome::ChecktrfIncomplete);
                    continue;
                }
                let mut exons: Vec<(u64, u64)> = Vec::new();
                let mut exoncov: Vec<f64> = Vec::new();
                let mut cov_bp_total = 0.0f64;
                let source_ext_start = if std::env::var_os("RUSTLE_SOURCE_CONTIG_EXTEND").is_some()
                {
                    rescue_nodes
                        .first()
                        .and_then(|&nid| source_contig_extension_start(graph, nid))
                } else {
                    None
                };
                let mut j = 0usize;
                while j < rescue_nodes.len() {
                    let nid = rescue_nodes[j];
                    let Some(node) = graph.nodes.get(nid) else {
                        j += 1;
                        continue;
                    };
                    let start = node.start;
                    let mut end = node.end;
                    let mut exon_bp_cov = 0.0f64;
                    if nid < local_nodecov.len() {
                        let nlen = node.length() as f64;
                        if local_nodecov[nid] > EPS && nlen > 0.0 {
                            let mut addcov = transfrags[t].abundance * nlen;
                            let rate = local_noderate
                                .get(nid)
                                .copied()
                                .filter(|r| *r > 0.0)
                                .unwrap_or(1.0);
                            let mut newnodecov = local_nodecov[nid] - addcov / rate;
                            if newnodecov < 0.0 {
                                addcov = local_nodecov[nid] * rate;
                                newnodecov = 0.0;
                            }
                            local_nodecov[nid] = newnodecov;
                            exon_bp_cov += addcov;
                        }
                    }
                    while j + 1 < rescue_nodes.len()
                        && nodes_are_contiguous(graph, rescue_nodes[j], rescue_nodes[j + 1])
                    {
                        j += 1;
                        if let Some(nn) = graph.nodes.get(rescue_nodes[j]) {
                            end = nn.end;
                            let nid2 = rescue_nodes[j];
                            if nid2 < local_nodecov.len() {
                                let nlen = nn.length() as f64;
                                if local_nodecov[nid2] > EPS && nlen > 0.0 {
                                    let mut addcov = transfrags[t].abundance * nlen;
                                    let rate = local_noderate
                                        .get(nid2)
                                        .copied()
                                        .filter(|r| *r > 0.0)
                                        .unwrap_or(1.0);
                                    let mut newnodecov = local_nodecov[nid2] - addcov / rate;
                                    if newnodecov < 0.0 {
                                        addcov = local_nodecov[nid2] * rate;
                                        newnodecov = 0.0;
                                    }
                                    local_nodecov[nid2] = newnodecov;
                                    exon_bp_cov += addcov;
                                }
                            }
                        }
                    }
                    if end > start {
                        exons.push((start, end));
                        let elen = len_half_open(start, end).max(1) as f64;
                        exoncov.push(exon_bp_cov / elen);
                        cov_bp_total += exon_bp_cov;
                    }
                    j += 1;
                }
                if exons.is_empty() {
                    transfrags[t].abundance = 0.0;
                    record_outcome!(t, SeedOutcome::ChecktrfRescueFail);
                    continue;
                }
                if let Some(ext_start) = source_ext_start {
                    if ext_start < exons[0].0 {
                        exons[0].0 = ext_start;
                    }
                }
                if std::env::var_os("RUSTLE_SINK_CONTIG_EXTEND").is_some() {
                    if let Some(last_nid) = rescue_nodes.last().copied() {
                        if let Some(ext_end) = sink_contig_extension_end(graph, last_nid) {
                            let li = exons.len() - 1;
                            if ext_end > exons[li].1 {
                                exons[li].1 = ext_end;
                            }
                        }
                    }
                }
                // Compute coverage and apply gates on the full node span first.
                let full_length: u64 = exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();

                // Minimum transcript length gate (default 200bp).
                if !transfrags[t].guide && full_length < config.min_transcript_length {
                    transfrags[t].abundance = 0.0;
                    record_outcome!(t, SeedOutcome::ChecktrfRescueFail);
                    continue;
                }

                // Coverage gate on the full (unclipped) span: reject rescued
                // transcripts with negligible coverage vs already-kept paths.
                let full_cov = if full_length > 0 { cov_bp_total / (full_length as f64) } else { 0.0 };
                let kept_max_cov = kept_paths.iter()
                    .map(|(_, cov, _, _)| *cov)
                    .fold(0.0f64, |a, b| a.max(b));
                if kept_max_cov > 0.0 && full_cov < 0.03 * kept_max_cov && !transfrags[t].guide {
                    transfrags[t].abundance = 0.0;
                    record_outcome!(t, SeedOutcome::ChecktrfRescueFail);
                    continue;
                }

                // Clip terminal exons to longstart/longend (read boundaries)
                // AFTER gates have been applied on the full span.
                {
                    let longstart = transfrags[t].longstart;
                    let longend = transfrags[t].longend;
                    if longstart > 0 && !exons.is_empty() {
                        if longstart > exons[0].0 && longstart <= exons[0].1 {
                            exons[0].0 = longstart;
                        }
                    }
                    if longend > 0 && !exons.is_empty() {
                        let li = exons.len() - 1;
                        if longend >= exons[li].0 && longend < exons[li].1 {
                            exons[li].1 = longend;
                        }
                    }
                }

                // Final length after clipping for per-bp coverage.
                let length: u64 = exons.iter().map(|(s, e)| len_half_open(*s, *e)).sum();
                let coverage = if length > 0 {
                    cov_bp_total / (length as f64)
                } else {
                    0.0
                };
                if let Some((lo, hi)) = trace_locus {
                    let in_range = transfrags[t].node_ids.iter().any(|&nid| {
                        graph
                            .nodes
                            .get(nid)
                            .map_or(false, |n| n.start <= hi && n.end >= lo)
                    });
                    if in_range {
                        let first_exon = exons.first().copied().unwrap_or((0, 0));
                        let last_exon = exons.last().copied().unwrap_or((0, 0));
                        eprintln!(
                            "[TRACE_CHECKTRF] t={} RESCUED_INDEP cov={:.4} exons={} first={}-{} last={}-{} longstart={} longend={} rescue_nodes_first={} ({}-{}) rescue_nodes_last={} ({}-{})",
                            t,
                            coverage,
                            exons.len(),
                            first_exon.0,
                            first_exon.1,
                            last_exon.0,
                            last_exon.1,
                            transfrags[t].longstart,
                            transfrags[t].longend,
                            rescue_nodes.first().copied().unwrap_or(source_id),
                            rescue_nodes
                                .first()
                                .and_then(|&nid| graph.nodes.get(nid))
                                .map(|n| n.start)
                                .unwrap_or(0),
                            rescue_nodes
                                .first()
                                .and_then(|&nid| graph.nodes.get(nid))
                                .map(|n| n.end)
                                .unwrap_or(0),
                            rescue_nodes.last().copied().unwrap_or(source_id),
                            rescue_nodes
                                .last()
                                .and_then(|&nid| graph.nodes.get(nid))
                                .map(|n| n.start)
                                .unwrap_or(0),
                            rescue_nodes
                                .last()
                                .and_then(|&nid| graph.nodes.get(nid))
                                .map(|n| n.end)
                                .unwrap_or(0)
                        );
                        eprintln!(
                            "[TRACE_CHECKTRF_OUT] t={} mode=rescued cov={:.4} nodes={} coords={} exons={}",
                            t,
                            coverage,
                            rescue_nodes.len(),
                            format_node_coord_list(graph, rescue_nodes),
                            format_exon_coord_list(&exons)
                        );
                    }
                }
                let first_node = rescue_nodes[0];
                let last_node = *rescue_nodes.last().unwrap_or(&first_node);
                out.push(Transcript {
                    chrom: bundle_chrom.to_string(),
                    strand: bundle_strand,
                    exons: exons.clone(),
                    coverage,
                    exon_cov: if exoncov.len() == exons.len() && exoncov.iter().any(|v| *v > 0.0) {
                        exoncov
                    } else {
                        vec![coverage; exons.len()]
                    },
                    tpm: 0.0,
                    fpkm: 0.0,
                    source: gtf_source_checktrf_rescue(&transfrags[t].guide_tid),
                    is_longread: long_read_mode,
                    longcov: transfrags[t].abundance, // 
                    bpcov_cov: 0.0,
                    transcript_id: None,
                    gene_id: None,
                    ref_transcript_id: None,
                    ref_gene_id: None,
                    hardstart: graph.nodes.get(first_node).map(|n| n.hardstart).unwrap_or(false),
                    hardend: graph.nodes.get(last_node).map(|n| n.hardend).unwrap_or(false),
                    vg_family_id: None, vg_copy_id: None, vg_family_size: None,
                });
                let out_idx = out.len() - 1;
                record_outcome!(t, SeedOutcome::ChecktrfRescued);
                if debug_detail {
                    eprint!(
                        "DEBUG_CHK t={} abund={:.4} outcome=rescued cov={:.4} nexons={} exons=",
                        t,
                        transfrags[t].abundance,
                        coverage,
                        exons.len()
                    );
                    for (i, (es, ee)) in exons.iter().enumerate() {
                        if i > 0 {
                            eprint!(",");
                        }
                        eprint!("{}-{}", es, ee);
                    }
                    eprintln!();
                }
                // an independently rescued checktrf transcript
                // becomes a new keeptrf entry immediately, so later checktrf items and the post-pass
                // can match against it.
                if let Some(ti) = trace_intron {
                    if path_has_intron(&rescue_nodes, graph, ti) {
                        eprintln!(
                            "[TRACE_INTRON_TF] stage=checktrf_rescue t={} intron={}-{} cov={:.4} nodes={} exons={} bundle={}",
                            t, ti.0, ti.1, coverage, rescue_nodes.len(), exons.len(), bundle_id
                        );
                    }
                }
                kept_paths.push((rescue_nodes.clone(), coverage, transfrags[t].guide, out_idx));
                transfrags[t].abundance = 0.0;
            }
        }
    }
    // parse_trflong post-pass: for remaining long-read transfrags with internal paths,
    // redistribute support to best kept predictions (no independent transcript creation here).
    if long_read_mode && !out.is_empty() && !kept_paths.is_empty() {
        let has_guide_kept = kept_paths.iter().any(|(_, _, g, _)| *g);
        for t in 0..transfrags.len() {
            if !transfrags[t].longread || transfrags[t].abundance <= EPS {
                continue;
            }
            if transfrags[t].node_ids.len() <= 1 {
                continue;
            }
            let first_raw = transfrags[t].node_ids.first().copied().unwrap_or(source_id);
            let last_raw = transfrags[t].node_ids.last().copied().unwrap_or(sink_id);
            if first_raw == source_id || last_raw == sink_id {
                continue;
            }
            let tf_nodes: Vec<usize> = transfrags[t]
                .node_ids
                .iter()
                .copied()
                .filter(|&n| n != source_id && n != sink_id)
                .collect();
            if tf_nodes.len() <= 1 {
                continue;
            }
            let mut tmatch: Vec<usize> = Vec::new();
            let mut abundancesum = 0.0;
            if has_guide_kept {
                let (m, s) = best_trf_match(&tf_nodes, &kept_paths, graph, true, Some(t));
                tmatch = m;
                abundancesum = s;
            }
            if tmatch.is_empty() {
                let (m, s) = best_trf_match(&tf_nodes, &kept_paths, graph, false, Some(t));
                tmatch = m;
                abundancesum = s;
            }
            if tmatch.is_empty() {
                continue;
            }
            let updated = redistribute_transfrag_to_matches(
                &transfrags[t],
                &tmatch,
                abundancesum,
                &kept_paths,
                graph,
                &mut out,
            );
            if updated {
                if zero_flux_set.contains(&t) {
                    zero_flux_rescued += 1;
                    if audit_zero_flux {
                        eprintln!(
                            "    [audit_zero_flux] rescued_secondpass idx={} matched_paths={} abundsum={:.4}",
                            t,
                            tmatch.len(),
                            abundancesum
                        );
                    }
                }
                transfrags[t].abundance = 0.0;
            }
        }
    }
    // zero ALL longread transfrag abundances after attribution
    if long_read_mode {
        for tf in transfrags.iter_mut() {
            if tf.longread {
                tf.abundance = 0.0;
            }
        }
    }

    // keeptrf second pass re-estimate coverage using push_max_flow
    // on short-read transfrags. If new coverage > original, update the prediction.
    // DISABLED for long-read mode: causes coverage inflation due to incorrect formula
    if false && long_read_mode && !kept_paths.is_empty() {
        // Sort kept_paths indices by longtrCmp: guides first, then abundance desc, then node count desc
        let mut keeptrf_order: Vec<usize> = (0..kept_paths.len()).collect();
        keeptrf_order.sort_unstable_by(|&a, &b| {
            let (ref nodes_a, cov_a, guide_a, oidx_a) = kept_paths[a];
            let (ref nodes_b, cov_b, guide_b, oidx_b) = kept_paths[b];
            // guides first
            guide_b
                .cmp(&guide_a)
                .then_with(|| {
                    // then by abundance (coverage * length) descending
                    let len_a = out
                        .get(oidx_a)
                        .map(|t| t.exons.iter().map(|(s, e)| e - s).sum::<u64>())
                        .unwrap_or(1) as f64;
                    let len_b = out
                        .get(oidx_b)
                        .map(|t| t.exons.iter().map(|(s, e)| e - s).sum::<u64>())
                        .unwrap_or(1) as f64;
                    let abund_a = cov_a * len_a;
                    let abund_b = cov_b * len_b;
                    abund_b
                        .partial_cmp(&abund_a)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .then_with(|| nodes_b.len().cmp(&nodes_a.len()))
        });

        let source_id = graph.source_id;
        let sink_id = graph.sink_id;

        for &ki in &keeptrf_order {
            let (ref inner_nodes, _orig_cov, _is_guide, out_idx) = kept_paths[ki];
            if inner_nodes.is_empty() || out_idx >= out.len() {
                continue;
            }

            // Build full path: source + inner_nodes + sink
            let mut full_path = Vec::with_capacity(inner_nodes.len() + 2);
            if inner_nodes[0] != source_id {
                full_path.push(source_id);
            }
            full_path.extend_from_slice(inner_nodes);
            if *inner_nodes.last().unwrap() != sink_id {
                full_path.push(sink_id);
            }

            // Ensure source→first and last→sink edges exist
            graph.add_edge(full_path[0], full_path[1]);
            let last_idx = full_path.len() - 1;
            graph.add_edge(full_path[last_idx - 1], full_path[last_idx]);

            // Run push_max_flow on this path (uses remaining short-read transfrags)
            let (flux, nodeflux, _full) =
                push_max_flow_seeded_full(&full_path, transfrags, graph, false, None);

            if flux <= EPS {
                continue;
            }

            // Compute coverage from nodeflux and local_nodecov (store_transcript)
            let mut cov_total = 0.0f64;
            let mut len_total = 0u64;
            for (pi, &nid) in full_path.iter().enumerate() {
                if nid == source_id || nid == sink_id || nid >= local_nodecov.len() {
                    continue;
                }
                let Some(node) = graph.nodes.get(nid) else {
                    continue;
                };
                let node_len = node.end.saturating_sub(node.start) as f64; // uses end-start+1 (1-based)
                let nflux = if pi < nodeflux.len() {
                    nodeflux[pi]
                } else {
                    0.0
                };
                let usedcov = local_nodecov[nid] * nflux * node_len;
                cov_total += usedcov;
                len_total += node.end.saturating_sub(node.start);
                // Deplete local_nodecov
                local_nodecov[nid] = (local_nodecov[nid] * (1.0 - nflux)).max(0.0);
            }

            if len_total == 0 || cov_total <= 0.0 {
                continue;
            }
            let new_cov = cov_total / len_total as f64;

            // if new coverage > original, update prediction
            let orig_cov = out[out_idx].coverage;
            if new_cov > orig_cov {
                let orig_len = out[out_idx].exons.iter().map(|(s, e)| e - s).sum::<u64>();
                let new_len = len_total;
                if new_len < orig_len && orig_len > 0 {
                    // scale by length ratio
                    out[out_idx].coverage = new_cov * new_len as f64 / orig_len as f64;
                } else {
                    out[out_idx].coverage = new_cov;
                }
                // Update exon_cov proportionally
                if orig_cov > 0.0 {
                    let ratio = out[out_idx].coverage / orig_cov;
                    for ec in out[out_idx].exon_cov.iter_mut() {
                        *ec *= ratio;
                    }
                }
            }
        }
    }
    if config.verbose && longrec_attempted > 0 {
        eprintln!(
            "    strict_longrec_port {}:{}({}) attempted={} succeeded={} fallback={} back_fail={} fwd_fail={} path_invalid={} back_unreachable_minpath={} back_no_reach={} back_no_choice={} back_exclude_no_support={} fwd_unreachable_maxpath={} fwd_no_reach={} fwd_no_choice={} fwd_exclude_no_support={}",
            bundle_chrom,
            bundle_strand,
            graph.n_nodes,
            longrec_attempted,
            longrec_succeeded,
            longrec_fallback,
            longrec_back_fail,
            longrec_fwd_fail,
            longrec_path_invalid,
            longrec_back_unreachable_minpath,
            longrec_back_no_reach,
            longrec_back_no_choice,
            longrec_back_exclude_no_support,
            longrec_fwd_unreachable_maxpath,
            longrec_fwd_no_reach,
            longrec_fwd_no_choice,
            longrec_fwd_exclude_no_support
        );
    }
    if (config.verbose || audit_zero_flux) && zero_flux_candidates > 0 {
        eprintln!(
            "    [audit_zero_flux] summary candidates={} rescued={} drop_nomatch={} drop_lowabund={} drop_empty={}",
            zero_flux_candidates,
            zero_flux_rescued,
            zero_flux_dropped_nomatch,
            zero_flux_dropped_lowabund,
            zero_flux_dropped_empty
        );
    }
    if let Some(summary) = longrec_summary {
        *summary = LongRecSummary {
            attempted: longrec_attempted,
            succeeded: longrec_succeeded,
            fallback: longrec_fallback,
            back_fail: longrec_back_fail,
            fwd_fail: longrec_fwd_fail,
            back_unreachable_minpath: longrec_back_unreachable_minpath,
            back_no_reach: longrec_back_no_reach,
            back_no_choice: longrec_back_no_choice,
            back_exclude_no_support: longrec_back_exclude_no_support,
            fwd_unreachable_maxpath: longrec_fwd_unreachable_maxpath,
            fwd_no_reach: longrec_fwd_no_reach,
            fwd_no_choice: longrec_fwd_no_choice,
            fwd_exclude_no_support: longrec_fwd_exclude_no_support,
        };
    }

    // EM post-processing for guide abundance distribution.
    // guides_pushmaxflow runs when: eonly mode, OR (non-eonly AND non-mixedMode with guides).
    if config.eonly || (!config.eonly && !mixed_mode && guided_mode) {
        // Build guide_entries from guide transcripts in output
        let guide_entries: Vec<(Vec<usize>, usize)> = out
            .iter()
            .enumerate()
            .filter(|(_, tx)| {
                tx.source
                    .as_deref()
                    .map_or(false, |s| s.starts_with("guide:"))
            })
            .filter_map(|(idx, tx)| {
                // Map exon coordinates back to graph node IDs
                let mut node_ids: Vec<usize> = Vec::new();
                for (start, end) in &tx.exons {
                    // Find nodes that match these exon coordinates
                    for (nid, node) in graph.nodes.iter().enumerate() {
                        if node.start == *start && node.end == *end {
                            node_ids.push(nid);
                            break;
                        }
                    }
                }
                if !node_ids.is_empty() {
                    // Sort node_ids to ensure correct order
                    node_ids.sort();
                    Some((node_ids, idx))
                } else {
                    None
                }
            })
            .collect();

        if !guide_entries.is_empty() {
            // DISABLED: run_eonly_guide_em causes coverage inflation in long-read mode
            // The EM algorithm adds node_remaining * olen / nlen to gcov, but node_remaining
            // is in abundance units, causing massive inflation (e.g., 71 -> 2,664,062).
            // let local_nodecov: Vec<f64> = (0..graph.n_nodes)
            //     .map(|i| graph.nodes.get(i).map(|n| n.coverage).unwrap_or(0.0))
            //     .collect();
            // run_eonly_guide_em(&guide_entries, &mut out, &local_nodecov, graph, transfrags);
        }
    }

    out
}

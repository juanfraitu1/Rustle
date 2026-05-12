//! trace events for pairwise debugging with the original algorithm.
//!
//! Output format matches instrumentation so logs can be grepped,
//! counted, and diff'd.
//!
//! ## Environment variables
//!
//! | Variable | Events |
//! |----------|--------|
//! | `RUSTLE_TRACE_UPDATE_ABUND` | PATH_update_abund (LR_ENTRY, LR_FLAGS, SINK/SOURCE_*, FINAL_NODES, ADD_TO_TRF, NEW_TRF, SET_LONGSTART/END) |
//! | `RUSTLE_TRACE_READ_PATTERN` | PATH_read_pattern (INTERSECT, ADD_NODE, EXON_SKIP, UNITIG_TRIM_*) |
//! | `RUSTLE_TRACE_LONGTRIM` | LONGTRIM_BOUND, LONGTRIM_BNODE, LONGTRIM_SPLIT |
//! | `RUSTLE_TRACE_CREATE_GRAPH` | NEW_NODE, NODE_SHRINK_*, NODE_FINAL_END, SRCSINK_COVBUILD_SOURCE |
//! | `RUSTLE_TRACE_PREDCLUSTER` | pred[N] full dump |
//! | `RUSTLE_TRACE_PRED_FATE` | PRED_FATE, FINAL_FATE, FILTER, FINAL_KEEP, ENTER, BEFORE_PAIRWISE, etc. |
//! | `RUSTLE_TRACE_BUILD_GRAPHS` | JUNC_COLOR_BREAK, JUNC_DELETE, ROUTE_STRANDED/UNSTRANDED |
//! | `RUSTLE_TRACE_READ=N,M` | Only emit PATH traces for read indices N, M |
//! | `RUSTLE_TRACE_LOCUS=start-end` | Only emit traces for paths overlapping [start,end) |
//! | `RUSTLE_BUNDLE_COV_DUMP=1` | With `RUSTLE_TRACE_LOCUS`, emit `RUSTLE_BUNDLE_COV_TSV` rows at predcluster stages (`transcript_filter.rs`) |
//! | `RUSTLE_SOURCE_HISTOGRAM=1` | After assembly, print stderr table of transcript counts by GTF `source` (`flow`, `junction_path`, …) |
//! | `RUSTLE_DISABLE_TERMINAL_ALT_ACCEPTOR` | Skip terminal alt-acceptor rescue (fewer `terminal_alt_acceptor` transcripts) |
//! | `RUSTLE_DISABLE_MICRO_EXON_RESCUE` | Skip micro-exon insertion rescue |
//! | `RUSTLE_CHECKTRF_ABUNDANCE_FLOOR_FRAC` | Stricter partial-depletion floor for checktrf (default `0.05`) |
//! | `RUSTLE_READTHR_LONGCOV_FALLBACK` | Allow `longcov >= readthr` when `coverage` is below threshold (off by default; uses `pred->cov`) |

use crate::path_extract::Transcript;

// --- Gate helpers ---

#[inline]
fn trace_log_style_active() -> bool {
    std::env::var_os("RUSTLE_TRACE_LOG_STYLE").is_some()
}

pub fn path_update_abund_trace() -> bool {
    std::env::var_os("RUSTLE_TRACE_UPDATE_ABUND").is_some() || trace_log_style_active()
}

pub fn path_read_pattern_trace() -> bool {
    std::env::var_os("RUSTLE_TRACE_READ_PATTERN").is_some() || trace_log_style_active()
}

pub fn longtrim_trace() -> bool {
    std::env::var_os("RUSTLE_TRACE_LONGTRIM").is_some() || trace_log_style_active()
}

pub fn create_graph_trace() -> bool {
    std::env::var_os("RUSTLE_TRACE_CREATE_GRAPH").is_some() || trace_log_style_active()
}

fn span_overlaps_trace(start: u64, end: u64) -> bool {
    let Some((lo, hi)) = parse_trace_locus() else {
        return true;
    };
    // Rust spans are half-open [start,end); env locus is given as start-end (inclusive intent),
    // but we treat it as half-open overlap check for consistency with other tracing.
    start <= hi && end >= lo
}

#[inline]
fn to_1based_inclusive(start: u64, end: u64) -> (u64, u64) {
    // Rust graph nodes are half-open [start,end). traces use inclusive coordinates.
    // Convert: [s,e) -> [s+1, e] in 1-based inclusive space.
    (start.saturating_add(1), end)
}

fn any_node_overlaps_trace(nodes: &[(usize, u64, u64)]) -> bool {
    nodes.iter().any(|&(_, s, e)| span_overlaps_trace(s, e))
}

pub fn predcluster_trace() -> bool {
    std::env::var_os("RUSTLE_TRACE_PREDCLUSTER").is_some() || trace_log_style_active()
}

pub fn pred_fate_trace() -> bool {
    std::env::var_os("RUSTLE_TRACE_PRED_FATE").is_some() || trace_log_style_active()
}

pub fn build_graphs_trace() -> bool {
    std::env::var_os("RUSTLE_TRACE_BUILD_GRAPHS").is_some() || trace_log_style_active()
}

/// Parse RUSTLE_TRACE_READ=N,M to filter PATH traces to specific read indices.
pub fn trace_read_match(read_idx: usize) -> bool {
    let Ok(val) = std::env::var("RUSTLE_TRACE_READ") else {
        return true;
    };
    for part in val.split(',') {
        if let Ok(n) = part.trim().parse::<usize>() {
            if n == read_idx {
                return true;
            }
        }
    }
    false
}

/// Parse RUSTLE_TRACE_LOCUS=start-end for locus filtering.
pub fn parse_trace_locus() -> Option<(u64, u64)> {
    let val = std::env::var("RUSTLE_TRACE_LOCUS").ok()?;
    let (a, b) = val.split_once('-')?;
    let start = a.trim().parse::<u64>().ok()?;
    let end = b.trim().parse::<u64>().ok()?;
    Some((start.min(end), start.max(end)))
}

// --- PATH_update_abund events (568K events) ---

pub fn path_set_longstart(read_idx: usize, s: usize, g: usize, rstart: u64, longstart: u64) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    if let Some((lo, hi)) = parse_trace_locus() {
        if longstart < lo || longstart > hi {
            return;
        }
    }
    eprintln!(
        "PATH_update_abund: SET_LONGSTART read[{}] s={} g={} rstart={} longstart={}",
        read_idx, s, g, rstart, longstart
    );
}

pub fn path_set_longend(read_idx: usize, s: usize, g: usize, rend: u64, longend: u64) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    if let Some((lo, hi)) = parse_trace_locus() {
        if longend < lo || longend > hi {
            return;
        }
    }
    eprintln!(
        "PATH_update_abund: SET_LONGEND read[{}] s={} g={} rend={} longend={}",
        read_idx, s, g, rend, longend
    );
}

pub fn path_final_nodes(
    read_idx: usize,
    s: usize,
    g: usize,
    nnodes: usize,
    nodes: &[(usize, u64, u64)],
    abund: f64,
    is_lr: bool,
) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    if !any_node_overlaps_trace(nodes) {
        return;
    }
    let nodes_str: String = nodes
        .iter()
        .map(|(nid, st, en)| format!(" {}({}-{})", nid, st, en))
        .collect();
    eprintln!(
        "PATH_update_abund: FINAL_NODES read[{}] s={} g={} nnodes={} nodes={} abund={:.4} is_lr={}",
        read_idx,
        s,
        g,
        nnodes,
        nodes_str.trim(),
        abund,
        is_lr
    );
}

pub fn path_add_to_trf(
    read_idx: usize,
    s: usize,
    g: usize,
    nodes: &[usize],
    abund: f64,
    total: f64,
    is_lr: bool,
    is_sr: bool,
) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    // Locus gating best-effort: if caller also emits FINAL_NODES, this is redundant.
    let nodes_str: String = nodes.iter().map(|n| format!(" {}", n)).collect();
    eprintln!(
        "PATH_update_abund: ADD_TO_TRF read[{}] s={} g={} nodes={} abund={:.4} total={:.4} is_lr={} is_sr={}",
        read_idx, s, g, nodes_str.trim(), abund, total, is_lr, is_sr
    );
}

pub fn path_lr_entry(
    read_idx: usize,
    s: usize,
    g: usize,
    gno: usize,
    nnodes: usize,
    rstart: u64,
    rend: u64,
    nodes: &[(usize, u64, u64)],
) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    if !any_node_overlaps_trace(nodes) {
        return;
    }
    let nodes_str: String = nodes
        .iter()
        .map(|(nid, st, en)| format!(" {}({}-{})", nid, st, en))
        .collect();
    eprintln!(
        "PATH_update_abund: LR_ENTRY read[{}] s={} g={} gno={} nnodes={} rstart={} rend={} nodes={}",
        read_idx, s, g, gno, nnodes, rstart, rend, nodes_str.trim()
    );
}

pub fn path_lr_flags(
    read_idx: usize,
    s: usize,
    g: usize,
    is_source: bool,
    is_sink: bool,
    node0_start: u64,
    node0_hardstart: bool,
    node_last_end: u64,
    node_last_hardend: bool,
    rstart: u64,
    rend: u64,
) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    if !span_overlaps_trace(node0_start, node_last_end) {
        return;
    }
    eprintln!(
        "PATH_update_abund: LR_FLAGS read[{}] s={} g={} is_source={} is_sink={} node0_start={} node0_hardstart={} nodeLast_end={} nodeLast_hardend={} rstart={} rend={}",
        read_idx,
        s,
        g,
        is_source as i32,
        is_sink as i32,
        node0_start,
        node0_hardstart as i32,
        node_last_end,
        node_last_hardend as i32,
        rstart,
        rend
    );
}

pub fn path_sink_check_enter(read_idx: usize, s: usize, g: usize, last_i: usize) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    eprintln!(
        "PATH_update_abund: SINK_CHECK_ENTER read[{}] s={} g={} last_i={}",
        read_idx, s, g, last_i
    );
}

pub fn path_source_check_enter(read_idx: usize, s: usize, g: usize, nnodes: usize) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    eprintln!(
        "PATH_update_abund: SOURCE_CHECK_ENTER read[{}] s={} g={} nnodes={}",
        read_idx, s, g, nnodes
    );
}

pub fn path_sink_proximity(
    read_idx: usize,
    s: usize,
    g: usize,
    i: usize,
    dist: u64,
    longintronanchor: bool,
    prev_nid: usize,
    prev_start: u64,
    prev_end: u64,
    prev_cov: f64,
    prev_hardend: bool,
    nchild: usize,
    cur_nid: usize,
    cur_start: u64,
    cur_end: u64,
    cur_cov: f64,
    cur_len: u64,
    covdrop_test: bool,
) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    eprintln!(
        "PATH_update_abund: SINK_PROXIMITY read[{}] s={} g={} i={} dist={} longintronanchor={} node[i-1]={}({}-{},cov={:.1},hardend={},nchild={}) node[i]={}({}-{},cov={:.1},len={}) covdrop_test={}",
        read_idx,
        s,
        g,
        i,
        dist,
        longintronanchor as i32,
        prev_nid,
        prev_start,
        prev_end,
        prev_cov,
        prev_hardend as i32,
        nchild,
        cur_nid,
        cur_start,
        cur_end,
        cur_cov,
        cur_len,
        covdrop_test as i32
    );
}

pub fn path_source_proximity(
    read_idx: usize,
    s: usize,
    g: usize,
    i: usize,
    dist: u64,
    longintronanchor: bool,
    cur_nid: usize,
    cur_start: u64,
    cur_end: u64,
    cur_cov: f64,
    next_nid: usize,
    next_start: u64,
    next_end: u64,
    next_cov: f64,
    next_hardstart: bool,
    nparent: usize,
) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    eprintln!(
        "PATH_update_abund: SOURCE_PROXIMITY read[{}] s={} g={} i={} dist={} longintronanchor={} node[i]={}({}-{},cov={:.1}) node[i+1]={}({}-{},cov={:.1},hardstart={},nparent={})",
        read_idx,
        s,
        g,
        i,
        dist,
        longintronanchor as i32,
        cur_nid,
        cur_start,
        cur_end,
        cur_cov,
        next_nid,
        next_start,
        next_end,
        next_cov,
        next_hardstart as i32,
        nparent
    );
}

pub fn path_sink_trim(
    read_idx: usize,
    s: usize,
    g: usize,
    trim_type: u8,
    trimmed_to: usize,
    remaining: &[(usize, u64, u64)],
) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    if !any_node_overlaps_trace(remaining) {
        return;
    }
    let rem_str: String = remaining
        .iter()
        .map(|(nid, st, en)| format!(" {}({}-{})", nid, st, en))
        .collect();
    eprintln!(
        "PATH_update_abund: SINK_TRIM[{}] read[{}] s={} g={} trimmed_to={} remaining={}",
        trim_type,
        read_idx,
        s,
        g,
        trimmed_to,
        rem_str.trim()
    );
}

pub fn path_source_trim(
    read_idx: usize,
    s: usize,
    g: usize,
    trim_type: u8,
    shifted: usize,
    remaining: &[(usize, u64, u64)],
) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    if !any_node_overlaps_trace(remaining) {
        return;
    }
    let rem_str: String = remaining
        .iter()
        .map(|(nid, st, en)| format!(" {}({}-{})", nid, st, en))
        .collect();
    eprintln!(
        "PATH_update_abund: SOURCE_TRIM[{}] read[{}] s={} g={} shifted={} remaining={}",
        trim_type,
        read_idx,
        s,
        g,
        shifted,
        rem_str.trim()
    );
}

pub fn path_new_trf(
    read_idx: usize,
    s: usize,
    g: usize,
    trf_idx: usize,
    nodes: &[(usize, u64, u64)],
    abund: f64,
    is_lr: bool,
    is_sr: bool,
) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    if !any_node_overlaps_trace(nodes) {
        return;
    }
    let nodes_str: String = nodes
        .iter()
        .map(|(nid, st, en)| format!(" {}({}-{})", nid, st, en))
        .collect();
    eprintln!(
        "PATH_update_abund: NEW_TRF read[{}] s={} g={} trf_idx={} nodes={} abund={:.4} is_lr={} is_sr={}",
        read_idx, s, g, trf_idx, nodes_str.trim(), abund, is_lr, is_sr
    );
}

pub fn path_skip_single_node(
    read_idx: usize,
    s: usize,
    g: usize,
    node: usize,
    node_start: u64,
    node_end: u64,
    is_lr: bool,
) {
    if !path_update_abund_trace() || !trace_read_match(read_idx) {
        return;
    }
    if !span_overlaps_trace(node_start, node_end) {
        return;
    }
    eprintln!(
        "PATH_update_abund: SKIP_SINGLE_NODE read[{}] s={} g={} node={}({}-{}) is_lr={}",
        read_idx, s, g, node, node_start, node_end, is_lr as i32
    );
}

// --- PATH_read_pattern events (2.3M events) ---

pub fn path_read_pattern_intersect(
    read_idx: usize,
    s: usize,
    gnode: usize,
    node_start: u64,
    node_end: u64,
    ngraph: usize,
    bp: u64,
) {
    if !path_read_pattern_trace() || !trace_read_match(read_idx) {
        return;
    }
    if !span_overlaps_trace(node_start, node_end) {
        return;
    }
    eprintln!(
        "PATH_read_pattern: read[{}] s={} INTERSECT gnode={}({}-{}) ngraph={} bp={}",
        read_idx, s, gnode, node_start, node_end, ngraph, bp
    );
}

pub fn path_read_pattern_add_node(
    read_idx: usize,
    s: usize,
    g: usize,
    gnode: usize,
    node_start: u64,
    node_end: u64,
    total_nodes: usize,
) {
    if !path_read_pattern_trace() || !trace_read_match(read_idx) {
        return;
    }
    if !span_overlaps_trace(node_start, node_end) {
        return;
    }
    eprintln!(
        "PATH_read_pattern: read[{}] s={} ADD_NODE g={} gnode={}({}-{}) total_nodes={}",
        read_idx, s, g, gnode, node_start, node_end, total_nodes
    );
}

pub fn path_read_pattern_exon_skip(
    read_idx: usize,
    s: usize,
    bnode_gnode: usize,
    node_start: u64,
    node_end: u64,
    ngraph: usize,
    k: usize,
    ncoord: usize,
) {
    if !path_read_pattern_trace() || !trace_read_match(read_idx) {
        return;
    }
    if !span_overlaps_trace(node_start, node_end) {
        return;
    }
    eprintln!(
        "PATH_read_pattern: read[{}] s={} EXON_SKIP bnode_gnode={}({}-{}) ngraph={} k={} ncoord={}",
        read_idx, s, bnode_gnode, node_start, node_end, ngraph, k, ncoord
    );
}

pub fn path_read_pattern_unitig_trim_first(read_idx: usize, s: usize, g: usize) {
    if !path_read_pattern_trace() || !trace_read_match(read_idx) {
        return;
    }
    eprintln!(
        "PATH_read_pattern: read[{}] s={} UNITIG_TRIM_FIRST g={}",
        read_idx, s, g
    );
}

pub fn path_read_pattern_unitig_trim_last(read_idx: usize, s: usize, g: usize, gnode: usize) {
    if !path_read_pattern_trace() || !trace_read_match(read_idx) {
        return;
    }
    eprintln!(
        "PATH_read_pattern: read[{}] s={} UNITIG_TRIM_LAST g={} gnode={}",
        read_idx, s, g, gnode
    );
}

// --- LONGTRIM_* events (126K events) ---

pub fn longtrim_bound_start(s: usize, bnode_start: u64, bnode_end: u64, start_at: i32, cov: f64) {
    if !longtrim_trace() {
        return;
    }
    if !span_overlaps_trace(bnode_start, bnode_end) {
        return;
    }
    let (bs, be) = to_1based_inclusive(bnode_start, bnode_end);
    let start_1b = start_at.saturating_add(1);
    eprintln!(
        "LONGTRIM_BOUND s={} bnode={}-{} start_at={} cov={:.1}",
        s, bs, be, start_1b, cov
    );
}

pub fn longtrim_bound_end(s: usize, bnode_start: u64, bnode_end: u64, end_at: i32, cov: f64) {
    if !longtrim_trace() {
        return;
    }
    if !span_overlaps_trace(bnode_start, bnode_end) {
        return;
    }
    let (bs, be) = to_1based_inclusive(bnode_start, bnode_end);
    let end_1b = end_at.saturating_add(1);
    eprintln!(
        "LONGTRIM_BOUND s={} bnode={}-{} end_at={} cov={:.1}",
        s, bs, be, end_1b, cov
    );
}

pub fn longtrim_bnode(
    s: usize,
    bundle_start: u64,
    bundle_end: u64,
    bnode_start: u64,
    bnode_end: u64,
    refstart: u64,
) {
    if !longtrim_trace() {
        return;
    }
    if !span_overlaps_trace(bnode_start, bnode_end) {
        return;
    }
    let (bundle_cs, bundle_ce) = to_1based_inclusive(bundle_start, bundle_end);
    let (bs, be) = to_1based_inclusive(bnode_start, bnode_end);
    let refstart_1b = refstart.saturating_add(1);
    eprintln!(
        "LONGTRIM_BNODE s={} bundle={}-{} bnode={}-{} refstart={}",
        s, bundle_cs, bundle_ce, bs, be, refstart_1b
    );
}

pub fn longtrim_split(
    s: usize,
    split_type: &str,
    split_at: i32,
    tmpcov: f64,
    bnode_start: u64,
    bnode_end: u64,
) {
    if !longtrim_trace() {
        return;
    }
    if !span_overlaps_trace(bnode_start, bnode_end) {
        return;
    }
    let (bs, be) = to_1based_inclusive(bnode_start, bnode_end);
    let split_1b = split_at.saturating_add(1);
    let split_label = match split_type {
        "start" => "start_split_at",
        "end" => "end_split_at",
        other => other,
    };
    eprintln!(
        "LONGTRIM_SPLIT s={} {}={} tmpcov={:.2} bnode={}-{}",
        s, split_label, split_1b, tmpcov, bs, be
    );
}

// --- create_graph node events (38K events) ---

pub fn create_graph_new_node(s: usize, g: usize, nodeid: usize, start: u64, end: u64, note: &str) {
    if !create_graph_trace() {
        return;
    }
    if !span_overlaps_trace(start, end) {
        return;
    }
    let (cs, ce) = to_1based_inclusive(start, end);
    eprintln!(
        "--- create_graph: NEW_NODE s={} g={} nodeid={} {}-{} ({})",
        s, g, nodeid, cs, ce, note
    );
}

pub fn create_graph_node_shrink_jstart(
    s: usize,
    g: usize,
    nodeid: usize,
    start: u64,
    old_end: u64,
    new_end: u64,
    pos: u64,
) {
    if !create_graph_trace() {
        return;
    }
    if !span_overlaps_trace(start, old_end) && !span_overlaps_trace(start, new_end) {
        return;
    }
    let (ws, we) = to_1based_inclusive(start, old_end);
    let (ns, ne) = to_1based_inclusive(start, new_end);
    let pos_1b = pos.saturating_add(1);
    eprintln!(
        "--- create_graph: NODE_SHRINK_JSTART s={} g={} nodeid={} was={}-{} now={}-{} junc_start={}",
        s, g, nodeid, ws, we, ns, ne, pos_1b
    );
}

pub fn create_graph_node_shrink_jend(
    s: usize,
    g: usize,
    nodeid: usize,
    start: u64,
    old_end: u64,
    new_end: u64,
    pos: u64,
) {
    if !create_graph_trace() {
        return;
    }
    if !span_overlaps_trace(start, old_end) && !span_overlaps_trace(start, new_end) {
        return;
    }
    let (ws, we) = to_1based_inclusive(start, old_end);
    let (ns, ne) = to_1based_inclusive(start, new_end);
    let pos_1b = pos.saturating_add(1);
    eprintln!(
        "--- create_graph: NODE_SHRINK_JEND s={} g={} nodeid={} was={}-{} now={}-{} pos={}",
        s, g, nodeid, ws, we, ns, ne, pos_1b
    );
}

pub fn create_graph_node_final_end(s: usize, g: usize, nodeid: usize, start: u64, end: u64) {
    if !create_graph_trace() {
        return;
    }
    if !span_overlaps_trace(start, end) {
        return;
    }
    let (cs, ce) = to_1based_inclusive(start, end);
    eprintln!(
        "--- create_graph: NODE_FINAL_END s={} g={} nodeid={} {}-{}",
        s, g, nodeid, cs, ce
    );
}

pub fn create_graph_srcsink_covbuild_source(
    s: usize,
    created: bool,
    icov: f64,
    parcov: f64,
    threshold: f64,
    abundance: f64,
) {
    if !create_graph_trace() {
        return;
    }
    let status = if created { "CREATED" } else { "SUPPRESSED" };
    eprintln!(
        "--- create_graph: SRCSINK_COVBUILD_SOURCE {} s={} icov={:.4} parcov={:.4} threshold={:.4} abundance={:.4}",
        status, s, icov, parcov, threshold, abundance
    );
}

// --- print_predcluster fate events ---

pub fn predcluster_enter(npred: usize, geneno: usize, checkincomplete: bool) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- print_predcluster: ENTER npred={} geneno={} checkincomplete={}",
        npred, geneno, checkincomplete as i32
    );
}

pub fn predcluster_before_pairwise(npred: usize, alive: usize) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- print_predcluster: BEFORE_PAIRWISE npred={} alive={}",
        npred, alive
    );
}

pub fn predcluster_pred_fate_summary(npred: usize, kept: usize, killed: usize) {
    if !pred_fate_trace() {
        return;
    }
    eprintln!(
        "--- print_predcluster: PRED_FATE_SUMMARY npred={} kept={} killed={}",
        npred, kept, killed
    );
}

pub fn predcluster_pred_fate(
    n: usize,
    start: u64,
    end: u64,
    cov: f64,
    strand: char,
    exons: usize,
    guide: &str,
    fate: &str,
) {
    if !pred_fate_trace() {
        return;
    }
    eprintln!(
        "--- print_predcluster: PRED_FATE n={} {}-{} cov={:.4} strand={} exons={} guide={} fate={}",
        n, start, end, cov, strand, exons, guide, fate
    );
}

pub fn predcluster_final_survivors(npred: usize, alive: usize, killed: usize) {
    if !pred_fate_trace() {
        return;
    }
    eprintln!(
        "--- print_predcluster: FINAL_SURVIVORS npred={} alive={} killed={}",
        npred, alive, killed
    );
}

pub fn predcluster_final_fate(
    n: usize,
    start: u64,
    end: u64,
    cov: f64,
    strand: char,
    exons: usize,
    guide: &str,
    fate: &str,
) {
    if !pred_fate_trace() {
        return;
    }
    eprintln!(
        "--- print_predcluster: FINAL_FATE n={} {}-{} cov={:.4} strand={} exons={} guide={} fate={}",
        n, start, end, cov, strand, exons, guide, fate
    );
}

pub fn predcluster_filter(
    pred_idx: usize,
    start: u64,
    end: u64,
    cov: f64,
    strand: char,
    exons: usize,
    killed_by: Option<(usize, u64, u64, f64)>,
    reason: &str,
) {
    if !predcluster_trace() {
        return;
    }
    match killed_by {
        Some((killer_idx, k_start, k_end, k_cov)) => {
            eprintln!(
                "--- print_predcluster: FILTER pred[{}] {}-{} cov={:.4} strand={} exons={} KILLED_BY pred[{}] {}-{} cov={:.4} reason={}",
                pred_idx, start, end, cov, strand, exons, killer_idx, k_start, k_end, k_cov, reason
            );
        }
        None => {
            eprintln!(
                "--- print_predcluster: FILTER pred[{}] {}-{} cov={:.4} strand={} exons={} reason={}",
                pred_idx, start, end, cov, strand, exons, reason
            );
        }
    }
}

pub fn predcluster_final_keep(
    n: usize,
    start: u64,
    end: u64,
    cov: f64,
    strand: char,
    exons: usize,
    guide: &str,
) {
    if !pred_fate_trace() {
        return;
    }
    eprintln!(
        "--- print_predcluster: FINAL_KEEP pred[{}] {}-{} cov={:.4} strand={} exons={} guide={}",
        n, start, end, cov, strand, exons, guide
    );
}

pub fn predcluster_final_drop(
    n: usize,
    start: u64,
    end: u64,
    cov: f64,
    strand: char,
    exons: usize,
    guide: &str,
    reason: &str,
) {
    if !pred_fate_trace() {
        return;
    }
    eprintln!(
        "--- print_predcluster: FINAL_DROP pred[{}] {}-{} cov={:.4} strand={} exons={} guide={} reason={}",
        n, start, end, cov, strand, exons, guide, reason
    );
}

pub fn predcluster_final_ledger(
    npred: usize,
    flag_alive: usize,
    genes: usize,
    transcripts_printed: usize,
) {
    if !pred_fate_trace() {
        return;
    }
    eprintln!(
        "--- print_predcluster: FINAL_LEDGER npred={} flag_alive={} genes={} transcripts_printed={}",
        npred, flag_alive, genes, transcripts_printed
    );
}

/// Full pred[N] dump: coords, coverage, readcov, strand, exon structure (pred[N] 10K entries).
pub fn predcluster_pred_full(n: usize, tx: &Transcript) {
    if !predcluster_trace() {
        return;
    }
    let (start, end) = tx
        .exons
        .first()
        .zip(tx.exons.last())
        .map(|(f, l)| (f.0, l.1))
        .unwrap_or((0, 0));
    let exons_str: String = tx
        .exons
        .iter()
        .map(|(s, e)| format!("{}-{}", s, e))
        .collect::<Vec<_>>()
        .join(",");
    let guide = tx.source.as_deref().unwrap_or("(novel)");
    eprintln!(
        "--- print_predcluster: pred[{}] {}-{} cov={:.4} readcov={:.4} strand={} exons=[{}] guide={}",
        n,
        start,
        end,
        tx.coverage,
        tx.longcov,
        tx.strand,
        exons_str,
        guide
    );
}

// --- build_graphs junction events ---

pub fn junc_color_break(
    read_idx: usize,
    i: usize,
    junc_start: u64,
    junc_end: u64,
    strand_now: i8,
    mm: f64,
    changeleft: bool,
    changeright: bool,
    nreads: f64,
    nreads_good: f64,
    good_junc: i32,
) {
    if !build_graphs_trace() {
        return;
    }
    eprintln!(
        "--- build_graphs: JUNC_COLOR_BREAK read[{}] i={} junc={}-{} strand_now={} mm={:.1} changeleft={} changeright={} nreads={:.1} nreads_good={:.1} good_junc={}",
        read_idx,
        i,
        junc_start,
        junc_end,
        strand_now,
        mm,
        changeleft as i32,
        changeright as i32,
        nreads,
        nreads_good,
        good_junc
    );
}

pub fn junc_delete(junc_start: u64, junc_end: u64, reason: &str) {
    if !build_graphs_trace() {
        return;
    }
    eprintln!(
        "--- build_graphs: JUNC_DELETE junc={}-{} reason={}",
        junc_start, junc_end, reason
    );
}

pub fn junc_demote(junc_start: u64, junc_end: u64, reason: &str) {
    if !build_graphs_trace() {
        return;
    }
    eprintln!(
        "--- build_graphs: JUNC_DEMOTE junc={}-{} reason={}",
        junc_start, junc_end, reason
    );
}

pub fn route_stranded(read_idx: usize, group: usize, strand: i8) {
    if !build_graphs_trace() {
        return;
    }
    eprintln!(
        "--- build_graphs: ROUTE_STRANDED read[{}] group={} strand={}",
        read_idx, group, strand
    );
}

pub fn route_unstranded(read_idx: usize, group: usize) {
    if !build_graphs_trace() {
        return;
    }
    eprintln!(
        "--- build_graphs: ROUTE_UNSTRANDED read[{}] group={}",
        read_idx, group
    );
}

// --- parse_trflong seed details (11K events) ---

pub fn seed_pathpat(seed_idx: usize, pathpat_bits: usize, node_count: usize, edge_count: usize) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- parse_trflong: SEED_PATHPAT seed[{}] pathpat_bits={} node_count={} edge_count={}",
        seed_idx, pathpat_bits, node_count, edge_count
    );
}

pub fn seed_decision(seed_idx: usize, decision: &str, reason: Option<&str>) {
    if !predcluster_trace() {
        return;
    }
    match reason {
        Some(r) => eprintln!(
            "--- parse_trflong: SEED_DECISION seed[{}] {} reason={}",
            seed_idx, decision, r
        ),
        None => eprintln!(
            "--- parse_trflong: SEED_DECISION seed[{}] {}",
            seed_idx, decision
        ),
    }
}

pub fn flow_pathpat(seed_idx: usize, pathpat_bits: usize, flux: f64) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- parse_trflong: FLOW_PATHPAT seed[{}] pathpat_bits={} flux={:.4}",
        seed_idx, pathpat_bits, flux
    );
}

pub fn depletion_before(seed_idx: usize, nodecov_sum: f64) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- parse_trflong: DEPLETION_BEFORE seed[{}] nodecov_sum={:.4}",
        seed_idx, nodecov_sum
    );
}

pub fn depletion_after(seed_idx: usize, nodecov_sum: f64, depletion: f64) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- parse_trflong: DEPLETION_AFTER seed[{}] nodecov_sum={:.4} depletion={:.4}",
        seed_idx, nodecov_sum, depletion
    );
}

pub fn hard_boundary_check(seed_idx: usize, thardstart: bool, thardend: bool, decision: &str) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- parse_trflong: HARD_BOUNDARY_CHECK seed[{}] thardstart={} thardend={} decision={}",
        seed_idx, thardstart as i32, thardend as i32, decision
    );
}

// --- fwd_to_sink/back_to_source step-level events ---

pub fn flow_step(direction: &str, seed_idx: usize, step: usize, node: usize, action: &str) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- {}: step seed[{}] step={} node={} action={}",
        direction, seed_idx, step, node, action
    );
}

pub fn pathpat_or(seed_idx: usize, or_mask: usize, current_pathpat: usize) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- flow: PATHPAT_OR seed[{}] or_mask={} current_pathpat={}",
        seed_idx, or_mask, current_pathpat
    );
}

pub fn tiebreak(seed_idx: usize, node: usize, choice: usize, reason: &str) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- flow: TIEBREAK seed[{}] node={} choice={} reason={}",
        seed_idx, node, choice, reason
    );
}

// --- get_trf_long checktrf details ---

pub fn checktrf_gate(seed_idx: usize, passed: bool, reason: &str) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- checktrf: CHECKTRF_GATE seed[{}] passed={} reason={}",
        seed_idx, passed as i32, reason
    );
}

pub fn checktrf_gate_skip(seed_idx: usize, reason: &str) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- checktrf: CHECKTRF_GATE_SKIP seed[{}] reason={}",
        seed_idx, reason
    );
}

pub fn checktrf_rescue_start(seed_idx: usize, match_type: &str) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- checktrf: CHECKTRF_RESCUE_START seed[{}] match_type={}",
        seed_idx, match_type
    );
}

pub fn checktrf_complete(seed_idx: usize, pred_idx: usize, cov: f64, nexons: usize) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- checktrf: CHECKTRF_COMPLETE seed[{}] pred_idx={} cov={:.4} nexons={}",
        seed_idx, pred_idx, cov, nexons
    );
}

pub fn checktrf_reject(seed_idx: usize, reason: &str) {
    if !predcluster_trace() {
        return;
    }
    eprintln!(
        "--- checktrf: CHECKTRF_REJECT seed[{}] reason={}",
        seed_idx, reason
    );
}

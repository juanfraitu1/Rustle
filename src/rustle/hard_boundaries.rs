//! Annotate graph nodes with polyA/T evidence and hardstart/hardend.
//!
//! StringTie-parity port (rlink.cpp:5927-5953, 6162-6171):
//! - Aggregate per-transfrag poly_{start,end}_{aligned,unaligned} onto the
//!   first and last transfrag-endpoint nodes (sum of actual counts, not binary)
//! - POLY_TAIL_STOP_COUNT = 8 minimum reads to set hardstart/hardend
//! - For a `-` strand bundle, transcript 3' polyA lives at LOW ref coord, i.e.
//!   BAM 5' softclip polyT → `poly_start_unaligned` on the first node.
//! - For a `+` strand bundle, transcript 3' polyA lives at HIGH ref coord,
//!   i.e. BAM 3' softclip polyA → `poly_end_unaligned` on the last node.

use crate::graph::{Graph, GraphTransfrag};

/// Min polyA/T evidence (summed read count) to set hardstart/hardend on a node.
/// StringTie uses 8 (rlink.cpp:349 `POLY_TAIL_STOP_COUNT`).
pub const POLY_TAIL_STOP_COUNT: u32 = 8;

pub fn annotate_hard_boundaries(
    graph: &mut Graph,
    transfrags: &[GraphTransfrag],
    bundle_strand: char,
) {
    for tf in transfrags {
        if tf.node_ids.is_empty() {
            continue;
        }
        let first_nid = tf.node_ids[0];
        let last_nid = tf.node_ids[tf.node_ids.len() - 1];
        // Ignore source/sink-only boundary transfrags (StringTie rlink.cpp:5785).
        if first_nid == graph.source_id || last_nid == graph.sink_id {
            continue;
        }
        // First node: aggregate poly_start_{aligned,unaligned} counts.
        if let Some(n) = graph.node_mut(first_nid) {
            n.poly_start_aligned_total = n
                .poly_start_aligned_total
                .saturating_add(tf.poly_start_aligned as u32);
            n.poly_start_unaligned_total = n
                .poly_start_unaligned_total
                .saturating_add(tf.poly_start_unaligned as u32);
            // - strand: polyA at tx 3' = low ref coord = first node. Require
            // unaligned-tail evidence (softclip polyT) to meet the threshold.
            if bundle_strand == '-'
                && n.poly_start_unaligned_total >= POLY_TAIL_STOP_COUNT
            {
                n.hardstart = true;
            }
        }
        // Last node: aggregate poly_end_{aligned,unaligned} counts.
        if let Some(n) = graph.node_mut(last_nid) {
            n.poly_end_aligned_total = n
                .poly_end_aligned_total
                .saturating_add(tf.poly_end_aligned as u32);
            n.poly_end_unaligned_total = n
                .poly_end_unaligned_total
                .saturating_add(tf.poly_end_unaligned as u32);
            // + strand: polyA at tx 3' = high ref coord = last node.
            if bundle_strand == '+'
                && n.poly_end_unaligned_total >= POLY_TAIL_STOP_COUNT
            {
                n.hardend = true;
            }
        }
        // Per-transfrag strong-polyA boost (StringTie rlink.cpp:6162,6171):
        // if a single transfrag itself has >= threshold unaligned polyA count,
        // promote hardstart/hardend on its terminal node. This catches high-
        // support single-tx cases without waiting for aggregate accumulation.
        if bundle_strand == '-'
            && tf.poly_start_unaligned as u32 >= POLY_TAIL_STOP_COUNT
        {
            if let Some(n) = graph.node_mut(first_nid) {
                n.hardstart = true;
            }
        }
        if bundle_strand == '+'
            && tf.longread
            && tf.poly_end_unaligned as u32 >= POLY_TAIL_STOP_COUNT
        {
            if let Some(n) = graph.node_mut(last_nid) {
                n.hardend = true;
            }
        }
    }
}

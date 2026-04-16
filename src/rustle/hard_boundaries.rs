//! Annotate graph nodes with polyA/T evidence and hardstart/hardend (annotate_hard_boundaries).
//!  6046-6058; reference script:6745-6762.
//! POLY_TAIL_STOP_COUNT = 8: min unaligned-tail reads to mark hardstart/hardend.

use crate::graph::{Graph, GraphTransfrag};

/// Min polyA/T evidence (read count) to set hardstart/hardend on a node (header / 
pub const POLY_TAIL_STOP_COUNT: u32 = 3;

/// Aggregate polyA/T evidence per node from transfrags and set hard boundaries.
///
/// - plus strand (`+`): promote `hardend` from poly-end evidence only
/// - minus strand (`-`): promote `hardstart` from poly-start evidence only
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
        // 5785: ignore source/sink-only boundary transfrags for poly evidence promotion.
        if first_nid == graph.source_id || last_nid == graph.sink_id {
            continue;
        }
        if let Some(n) = graph.node_mut(first_nid) {
            // the original algorithm uses a read-count style threshold (POLY_TAIL_STOP_COUNT).
            // Treat any nonzero per-transfrag tail evidence as 1.
            //
            // Important: for '-' strand bundles, the transcript 3' polyA tail aligns to the
            // left genomic end, which corresponds to the read's 3' end (poly_end signal),
            // not the read's 5' end (poly_start / polyT signal).
            let inc = if bundle_strand == '-' {
                u32::from(tf.poly_end_unaligned > 0)
            } else {
                u32::from(tf.poly_start_unaligned > 0)
            };
            n.poly_start_unaligned_total = n.poly_start_unaligned_total.saturating_add(inc);
            if bundle_strand == '-'
                && n.poly_start_unaligned_total >= POLY_TAIL_STOP_COUNT
                && tf.longread
            {
                n.hardstart = true;
            }
        }
        if let Some(n) = graph.node_mut(last_nid) {
            let inc = u32::from(tf.poly_end_unaligned > 0);
            n.poly_end_unaligned_total = n.poly_end_unaligned_total.saturating_add(inc);
            if bundle_strand == '+'
                && n.poly_end_unaligned_total >= POLY_TAIL_STOP_COUNT
                && tf.longread
            {
                n.hardend = true;
            }
        }
    }
}

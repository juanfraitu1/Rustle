//! Nascent RNA transcript detection from intronic coverage (C++ reference).
//!
//! When `--nasc` is enabled, identifies partially-transcribed (nascent) RNA by looking
//! for coverage drops at intron boundaries. Creates nascent predictions linked to
//! their parent transcript.

use crate::bitvec::GBitVec;
use crate::bpcov::Bpcov;
use crate::graph::{Graph, GraphTransfrag};
use crate::path_extract::Transcript;

const DROP: f64 = 0.5;
const EPSILON: f64 = 1e-8;

/// A nascent transfrag candidate: partial path + coverage gap region.
pub struct NascentTf {
    pub node_ids: Vec<usize>,
    pub pattern: GBitVec,
    /// Start of the intronic gap (node index where coverage drops).
    pub longstart: usize,
    /// End of the intronic gap (node index where transcript resumes).
    pub longend: usize,
}

/// Identify nascent transcript candidates from coverage gaps in extracted path.
///
/// C++ create_nascent (C++ reference): for each intron gap in the path,
/// checks if coverage drops at the boundary. If so, creates a nascent transfrag
/// extending into the intron.
///
/// `sno`: 0 = reverse strand, >0 = forward strand.
pub fn create_nascent(
    path: &[usize],
    strand: char,
    graph: &Graph,
    bpcov: &Bpcov,
    bundle_start: u64,
) -> Vec<NascentTf> {
    let mut nascents = Vec::new();
    let sno = if strand == '-' { 0 } else { 1 };

    // C++ iterates from i=2 to path.len()-1 looking for gaps (path[i] > path[i-1]+1)
    for i in 2..path.len().saturating_sub(1) {
        if path[i] <= path[i - 1] + 1 {
            continue;
        }

        if sno > 0 {
            // Forward strand: check coverage drop at right boundary of path[i-1]
            let prev_node = match graph.nodes.get(path[i - 1]) {
                Some(n) => n,
                None => continue,
            };
            let next_node_id = path[i - 1] + 1;
            let next_node = match graph.nodes.get(next_node_id) {
                Some(n) => n,
                None => continue,
            };

            // Must be contiguous
            if prev_node.end + 1 != next_node.start {
                continue;
            }

            let left_idx = prev_node.end.saturating_sub(bundle_start) as usize;
            let right_idx = next_node.start.saturating_sub(bundle_start) as usize;
            let abund_left = bpcov.cov.get(left_idx).copied().unwrap_or(0.0);
            let abund_right = bpcov.cov.get(right_idx).copied().unwrap_or(0.0);
            if abund_right > abund_left * DROP {
                continue;
            }

            // Build nascent transfrag: path[0..i] + intronic nodes
            let psize = graph.pattern_size();
            let mut pat = GBitVec::new(psize);
            let mut nodes = Vec::new();
            for j in 0..i {
                nodes.push(path[j]);
                if path[j] < psize {
                    pat.set_bit(path[j]);
                }
                if j > 0 {
                    if let Some(eid) = graph.edge_bit_index(path[j - 1], path[j]) {
                        pat.set_bit(eid);
                    }
                }
            }
            // Extend into intron nodes
            let mut prev = path[i - 1];
            for j in (path[i - 1] + 1)..path[i] {
                let pn = match graph.nodes.get(prev) {
                    Some(n) => n,
                    None => break,
                };
                let cn = match graph.nodes.get(j) {
                    Some(n) => n,
                    None => break,
                };
                if pn.end + 1 != cn.start {
                    break;
                }
                nodes.push(j);
                if j < psize {
                    pat.set_bit(j);
                }
                if let Some(eid) = graph.edge_bit_index(prev, j) {
                    pat.set_bit(eid);
                }
                prev = j;
            }

            nascents.push(NascentTf {
                node_ids: nodes,
                pattern: pat,
                longstart: path[i - 1],
                longend: path[i],
            });
            if std::env::var_os("RUSTLE_PARITY_NASCENT").is_some() {
                let n = nascents.last().unwrap();
                let gap_start = graph.nodes.get(n.longstart).map(|n| n.start).unwrap_or(0);
                let gap_end = graph.nodes.get(n.longend).map(|n| n.end).unwrap_or(0);
                eprintln!("PARITY_NASCENT_CANDIDATE gap_start={} gap_end={} nodes={} longstart={} longend={}",
                    gap_start, gap_end, n.node_ids.len(), n.longstart, n.longend);
            }
        } else {
            // Reverse strand: check coverage drop at left boundary of path[i]
            let cur_node = match graph.nodes.get(path[i]) {
                Some(n) => n,
                None => continue,
            };
            let prev_node_id = path[i] - 1;
            let prev_node = match graph.nodes.get(prev_node_id) {
                Some(n) => n,
                None => continue,
            };

            if prev_node.end + 1 != cur_node.start {
                continue;
            }

            let left_idx = prev_node.end.saturating_sub(bundle_start) as usize;
            let right_idx = cur_node.start.saturating_sub(bundle_start) as usize;
            let abund_left = bpcov.cov.get(left_idx).copied().unwrap_or(0.0);
            let abund_right = bpcov.cov.get(right_idx).copied().unwrap_or(0.0);
            if abund_left > abund_right * DROP {
                continue;
            }

            // Build nascent: intronic nodes + path[i..]
            let psize = graph.pattern_size();
            let mut pat = GBitVec::new(psize);
            let mut nodes = Vec::new();

            // Find leftmost contiguous intronic node
            let mut j = path[i] as i64 - 1;
            while j > path[i - 1] as i64 {
                let pn = graph.nodes.get((j - 1) as usize);
                let cn = graph.nodes.get(j as usize);
                match (pn, cn) {
                    (Some(pn), Some(cn)) if pn.end + 1 == cn.start => {
                        j -= 1;
                    }
                    _ => break,
                }
            }

            // Add intronic nodes
            let start_j = j as usize;
            let mut prev_nid = start_j;
            for k in start_j..path[i] {
                nodes.push(k);
                if k < psize {
                    pat.set_bit(k);
                }
                if k > start_j {
                    if let Some(eid) = graph.edge_bit_index(k - 1, k) {
                        pat.set_bit(eid);
                    }
                }
                prev_nid = k;
            }

            // Add path[i..]
            for k in i..path.len() {
                nodes.push(path[k]);
                if path[k] < psize {
                    pat.set_bit(path[k]);
                }
                if k > i || !nodes.is_empty() {
                    let prev = if k == i { prev_nid } else { path[k - 1] };
                    if let Some(eid) = graph.edge_bit_index(prev, path[k]) {
                        pat.set_bit(eid);
                    }
                }
            }

            nascents.push(NascentTf {
                node_ids: nodes,
                pattern: pat,
                longstart: path[i - 1],
                longend: path[i],
            });
            if std::env::var_os("RUSTLE_PARITY_NASCENT").is_some() {
                let n = nascents.last().unwrap();
                let gap_start = graph.nodes.get(n.longstart).map(|n| n.start).unwrap_or(0);
                let gap_end = graph.nodes.get(n.longend).map(|n| n.end).unwrap_or(0);
                eprintln!("PARITY_NASCENT_CANDIDATE gap_start={} gap_end={} nodes={} longstart={} longend={}",
                    gap_start, gap_end, n.node_ids.len(), n.longstart, n.longend);
            }
        }
    }

    nascents
}

/// Compute max-flow for a nascent transfrag through the graph.
///
/// C++ nascent2max_flow (C++ reference). Simplified version:
/// computes forward/backward flow through nascent path using capacity ratios.
///
/// Returns (flow, nodeflux) where nodeflux[i] is the fraction of node coverage
/// to allocate to the nascent transcript.
pub fn nascent_max_flow(
    nascent: &NascentTf,
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    nodecov: &[f64],
) -> (f64, Vec<f64>) {
    let n = nascent.node_ids.len();
    let mut nodeflux = vec![0.0f64; n];

    if n == 0 {
        return (0.0, nodeflux);
    }

    // Compute capacity left/right and sum left/right for each node
    let mut cap_left = vec![0.0f64; n];
    let mut cap_right = vec![0.0f64; n];
    let mut sum_left = vec![0.0f64; n];
    let mut sum_right = vec![0.0f64; n];

    for (idx, &nid) in nascent.node_ids.iter().enumerate() {
        if nid == graph.source_id || nid == graph.sink_id {
            continue;
        }
        let node = match graph.nodes.get(nid) {
            Some(n) => n,
            None => continue,
        };

        for &tidx in &node.trf_ids {
            let tf = match transfrags.get(tidx) {
                Some(t) if t.abundance > EPSILON => t,
                _ => continue,
            };

            // Check if transfrag is compatible with nascent pattern
            let compatible = nascent.pattern.contains_pattern(&tf.pattern);

            if compatible {
                let tf_start = tf.node_ids.first().copied().unwrap_or(0);
                let tf_end = tf.node_ids.last().copied().unwrap_or(0);

                if tf_start < nid {
                    sum_left[idx] += tf.abundance;
                    cap_left[idx] += tf.abundance;
                }
                if tf_end > nid {
                    sum_right[idx] += tf.abundance;
                    cap_right[idx] += tf.abundance;
                }
            } else {
                let tf_start = tf.node_ids.first().copied().unwrap_or(0);
                let tf_end = tf.node_ids.last().copied().unwrap_or(0);
                if nid > tf_start {
                    sum_left[idx] += tf.abundance;
                }
                if nid < tf_end {
                    sum_right[idx] += tf.abundance;
                }
            }
        }

        if cap_left[idx] <= 0.0 || cap_right[idx] <= 0.0 {
            return (0.0, nodeflux);
        }
    }

    // Forward pass to compute flow
    let min_idx = 0;
    let max_idx = n;
    if min_idx >= max_idx {
        return (0.0, nodeflux);
    }

    let mut prev_flow = cap_left[min_idx];
    for idx in min_idx..max_idx {
        if sum_left[idx] <= 0.0 || sum_right[idx] <= 0.0 {
            return (0.0, nodeflux);
        }
        let perc_left = prev_flow / sum_left[idx];
        let mut perc_right = cap_right[idx] / sum_right[idx];
        if perc_right > perc_left {
            perc_right = perc_left;
        }
        prev_flow = perc_right * sum_right[idx];
    }

    if prev_flow <= 0.0 {
        return (0.0, nodeflux);
    }

    // Backward pass
    for idx in (min_idx..max_idx).rev() {
        if sum_right[idx] > 0.0 {
            nodeflux[idx] = prev_flow / sum_right[idx];
        }
        cap_right[idx] = prev_flow;
        if sum_left[idx] > 0.0 {
            prev_flow = nodeflux[idx] * sum_left[idx];
        }
    }

    if std::env::var_os("RUSTLE_PARITY_NASCENT").is_some() {
        eprintln!("PARITY_NASCENT_FLOW nodes={} flow={:.4}", n, nodeflux[min_idx]);
    }

    (nodeflux[min_idx], nodeflux)
}

/// Store a nascent transcript prediction from flow results.
///
/// C++ store_nascenttranscript (C++ reference).
/// Creates a Transcript from the nascent transfrag's nodes and flow.
pub fn store_nascent_transcript(
    nascent: &NascentTf,
    nodeflux: &[f64],
    nodecov: &mut [f64],
    graph: &Graph,
    strand: char,
    chrom: &str,
) -> Option<Transcript> {
    let mut cov = 0.0;
    let mut len: u64 = 0;
    let mut exons: Vec<(u64, u64)> = Vec::new();
    let mut exon_cov: Vec<f64> = Vec::new();
    let mut cur_excov = 0.0;

    for (idx, &nid) in nascent.node_ids.iter().enumerate() {
        let node = match graph.nodes.get(nid) {
            Some(n) => n,
            None => continue,
        };
        if nid == graph.source_id || nid == graph.sink_id {
            continue;
        }

        let flux = nodeflux.get(idx).copied().unwrap_or(0.0);
        let node_len = node.end.saturating_sub(node.start) as f64;
        let used_cov = nodecov.get(nid).copied().unwrap_or(0.0) * flux * node_len;

        // Extend or create exon
        if let Some(last) = exons.last_mut() {
            if node.start <= last.1 + 1 {
                // Contiguous: extend
                last.1 = node.end;
                cur_excov += used_cov;
            } else {
                // Gap: finish previous exon, start new one
                let prev_len = last.1.saturating_sub(last.0) + 1;
                exon_cov.push(if prev_len > 0 { cur_excov / prev_len as f64 } else { 0.0 });
                exons.push((node.start, node.end));
                cur_excov = used_cov;
            }
        } else {
            exons.push((node.start, node.end));
            cur_excov = used_cov;
        }

        let nlen = node.end.saturating_sub(node.start) + 1;
        len += nlen;
        cov += used_cov;

        // Deplete coverage for intronic nodes
        if nid > nascent.longstart && nid < nascent.longend {
            if let Some(nc) = nodecov.get_mut(nid) {
                *nc *= 1.0 - flux;
            }
        }
    }

    // Finalize last exon coverage
    if let Some(last) = exons.last() {
        let prev_len = last.1.saturating_sub(last.0) + 1;
        exon_cov.push(if prev_len > 0 { cur_excov / prev_len as f64 } else { 0.0 });
    }

    if len > 0 {
        cov /= len as f64;
    }

    if cov <= EPSILON {
        return None;
    }

    if std::env::var_os("RUSTLE_PARITY_NASCENT").is_some() {
        let introns: String = exons.windows(2)
            .map(|w| format!("{}-{}", w[0].1 + 1, w[1].0 - 1))
            .collect::<Vec<_>>()
            .join(",");
        eprintln!("PARITY_NASCENT_STORED nexons={} cov={:.4} introns={}", exons.len(), cov, introns);
    }
Some(Transcript {
    chrom: chrom.to_string(),
    strand,
    exons,
    coverage: cov,
    exon_cov,
    tpm: 0.0,
    fpkm: 0.0,
    source: Some("nascent".to_string()),
    is_longread: true,
    longcov: cov,
    bpcov_cov: 0.0,
    transcript_id: None,
    gene_id: None,
    ref_transcript_id: None,
    ref_gene_id: None,
    hardstart: true,
    hardend: true,
})
        nascent_chain: Vec::new(),
        mergename: Some('n'),
    })
}

/// Process nascent transcription for all extracted transcripts.
///
/// C++ remove_nascent_transcription (C++ reference).
/// For each nascent candidate, computes flow and stores as nascent prediction.
pub fn remove_nascent_transcription(
    nascents: &[NascentTf],
    graph: &Graph,
    transfrags: &[GraphTransfrag],
    nodecov: &mut [f64],
    strand: char,
    chrom: &str,
    parent: &mut Transcript,
) {
    for nascent in nascents {
        let (flux, nodeflux) = nascent_max_flow(nascent, graph, transfrags, nodecov);
        if flux <= EPSILON {
            continue;
        }
        if let Some(ntx) = store_nascent_transcript(nascent, &nodeflux, nodecov, graph, strand, chrom) {
            parent.nascent_chain.push(Box::new(ntx));
        }
    }
}

/// Check if prediction `p` is a nascent (partial) of prediction `n`.
///
/// C++ mark_if_nascent (C++ reference).
/// Returns true if p's intron chain is a subset of n's intron chain (strand-aware).
pub fn mark_if_nascent(parent: &Transcript, candidate: &mut Transcript) -> bool {
    if parent.strand != candidate.strand {
        return false;
    }

    let p_exons = &candidate.exons;
    let n_exons = &parent.exons;

    if p_exons.len() >= n_exons.len() {
        return false;
    }

    if parent.strand == '+' {
        // Forward: p must match first p_exons.len() introns of n
        if p_exons[0].0 > n_exons[0].1 {
            return false;
        }
        for i in 1..p_exons.len() {
            if p_exons[i - 1].1 != n_exons[i - 1].1 {
                return false; // Intron donor mismatch
            }
            if p_exons[i].0 != n_exons[i].0 {
                return false; // Intron acceptor mismatch
            }
        }
        let last_p = p_exons.len() - 1;
        if p_exons[last_p].1 < n_exons[last_p].1 {
            return false;
        }
        // Truncate candidate to end before next exon of parent
        let next_idx = p_exons.len();
        if next_idx < n_exons.len() {
            candidate.exons[last_p].1 = n_exons[next_idx].0 - 1;
        }
        // Extend start to match parent
        candidate.exons[0].0 = n_exons[0].0;
        true
    } else if parent.strand == '-' {
        // Reverse: p must match last p_exons.len() introns of n
        let p_last = p_exons.len() - 1;
        let n_last = n_exons.len() - 1;
        if p_exons[p_last].1 < n_exons[n_last].0 {
            return false;
        }
        let mut pi = p_last;
        let mut ni = n_last;
        while pi > 0 {
            if p_exons[pi - 1].1 != n_exons[ni - 1].1 {
                return false;
            }
            if p_exons[pi].0 != n_exons[ni].0 {
                return false;
            }
            pi -= 1;
            ni -= 1;
        }
        if p_exons[0].0 > n_exons[ni].0 {
            return false;
        }
        // Truncate candidate start to after previous intron of parent
        if ni > 0 {
            candidate.exons[0].0 = n_exons[ni - 1].1 + 1;
        }
        // Extend end to match parent
        candidate.exons[p_last].1 = n_exons[n_last].1;
        true
    } else {
        false
    }
}

/// Find and link nascent predictions to their parent transcripts.
///
/// C++ find_nascent_link (C++ reference).
/// For each prediction with fewer exons and lower coverage, checks if it's
/// a nascent of a higher-coverage prediction with more exons.
pub fn find_nascent_links(predictions: &mut [Transcript], start_idx: usize) {
    let n = predictions.len();
    for p in start_idx..n {
        if predictions[p].ref_transcript_id.is_some() {
            continue;
        }
        let p_nex = predictions[p].exons.len();
        for n_idx in start_idx..n {
            if n_idx == p {
                continue;
            }
            if predictions[n_idx].ref_transcript_id.is_some() {
                continue;
            }
            if predictions[n_idx].mergename == Some('n') {
                continue;
            }
            if predictions[n_idx].coverage <= predictions[p].coverage {
                continue;
            }
            if predictions[n_idx].exons.len() <= p_nex {
                continue;
            }
            if predictions[n_idx].strand != predictions[p].strand
                && predictions[p].strand != '.'
            {
                continue;
            }

            // Need to clone for borrow checker since mark_if_nascent mutates candidate
            let parent_clone = predictions[n_idx].clone();
            if mark_if_nascent(&parent_clone, &mut predictions[p]) {
                predictions[p].mergename = Some('n');
                if std::env::var_os("RUSTLE_PARITY_NASCENT").is_some() {
                    eprintln!("PARITY_NASCENT_LINK parent={} child={} parent_nexons={} child_nexons={}",
                        n_idx, p, predictions[n_idx].exons.len(), predictions[p].exons.len());
                }
                // Move p's nascent chain into n's chain
                let p_chain = std::mem::take(&mut predictions[p].nascent_chain);
                let p_clone = Box::new(predictions[p].clone());
                predictions[n_idx].nascent_chain.push(p_clone);
                predictions[n_idx].nascent_chain.extend(p_chain);
                break;
            }
        }
    }
}

/// Reconcile nascent chains when two predictions are merged.
///
/// C++ reconcile_nascents (C++ reference).
/// Adjusts endpoints of winner's nascents and merges loser's nascent chain.
pub fn reconcile_nascents(winner: &mut Transcript, loser: &Transcript) {
    // Adjust winner's existing nascents to match loser's endpoints if extended
    for nasc in &mut winner.nascent_chain {
        if nasc.exons.first().map(|e| e.0) == winner.exons.first().map(|e| e.0)
            && loser.exons.first().map(|e| e.0) != winner.exons.first().map(|e| e.0)
        {
            if let (Some(nasc_first), Some(loser_first)) =
                (nasc.exons.first_mut(), loser.exons.first())
            {
                nasc_first.0 = loser_first.0;
            }
        }
        if nasc.exons.last().map(|e| e.1) == winner.exons.last().map(|e| e.1)
            && loser.exons.last().map(|e| e.1) != winner.exons.last().map(|e| e.1)
        {
            if let (Some(nasc_last), Some(loser_last)) =
                (nasc.exons.last_mut(), loser.exons.last())
            {
                nasc_last.1 = loser_last.1;
            }
        }
    }

    // Merge loser's nascent chain into winner's, maintaining coordinate order
    for loser_nasc in &loser.nascent_chain {
        let loser_start = loser_nasc.exons.first().map(|e| e.0).unwrap_or(0);
        let loser_end = loser_nasc.exons.last().map(|e| e.1).unwrap_or(0);

        // Check if we can merge with an existing nascent
        let mut merged = false;
        for existing in &mut winner.nascent_chain {
            let ex_start = existing.exons.first().map(|e| e.0).unwrap_or(0);
            let ex_end = existing.exons.last().map(|e| e.1).unwrap_or(0);
            if ex_start == loser_start && ex_end == loser_end {
                // Same nascent: merge coverage
                existing.coverage += loser_nasc.coverage;
                for (i, ec) in existing.exon_cov.iter_mut().enumerate() {
                    if let Some(lc) = loser_nasc.exon_cov.get(i) {
                        *ec += lc;
                    }
                }
                merged = true;
                break;
            }
        }
        if !merged {
            winner.nascent_chain.push(loser_nasc.clone());
        }
    }
}

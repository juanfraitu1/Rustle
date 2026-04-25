//! Parse reference GTF for trace mode: load transcripts with exons (0-based start inclusive, end exclusive).
//! process_refguides: map guides to graph paths and set hardstart/hardend (process_refguides).

use crate::types::{DetHashMap as HashMap, DetHashSet as HashSet};
use anyhow::Result;
use std::path::Path;

use crate::graph::Graph;

/// Reference transcript from GTF: 0-based coordinates (exon end exclusive, matching our Transcript).
#[derive(Debug, Clone)]
pub struct RefTranscript {
    pub id: String,
    pub chrom: String,
    pub strand: char,
    /// Exons (start inclusive, end exclusive) 0-based
    pub exons: Vec<(u64, u64)>,
}

/// Parse GTF: collect exon lines by transcript_id, build RefTranscript per transcript.
/// GTF is 1-based inclusive; we convert to 0-based (start_incl, end_excl).
pub fn parse_reference_gtf<P: AsRef<Path>>(path: P) -> Result<Vec<RefTranscript>> {
    let s = std::fs::read_to_string(path)?;
    // transcript_id -> (chrom, strand, exons)
    let mut by_tid: HashMap<String, (String, char, Vec<(u64, u64)>)> = Default::default();
    for line in s.lines() {
        if line.starts_with('#') {
            continue;
        }
        let mut fields = line.splitn(9, '\t');
        let chrom = fields.next().unwrap_or("").to_string();
        let _source = fields.next().unwrap_or("");
        let feature = fields.next().unwrap_or("");
        let start_1b: u64 = fields.next().and_then(|x| x.parse().ok()).unwrap_or(0);
        let end_1b: u64 = fields.next().and_then(|x| x.parse().ok()).unwrap_or(0);
        let _score = fields.next();
        let strand_c = fields.next().unwrap_or(".").chars().next().unwrap_or('.');
        let _frame = fields.next();
        let attrs = fields.next().unwrap_or("");
        if chrom.is_empty() || feature.is_empty() {
            continue;
        }
        let strand = match strand_c {
            '-' => '-',
            '+' => '+',
            _ => '.',
        };
        let transcript_id = parse_gtf_attr(attrs, "transcript_id").unwrap_or_default();
        if transcript_id.is_empty() {
            continue;
        }
        if feature == "exon" {
            let start_0 = start_1b.saturating_sub(1);
            let end_0_excl = end_1b;
            by_tid
                .entry(transcript_id.clone())
                .or_insert_with(|| (chrom.clone(), strand, Vec::new()))
                .2
                .push((start_0, end_0_excl));
        }
    }
    let mut out = Vec::new();
    for (id, (chrom, strand, mut exons)) in by_tid {
        exons.sort_by_key(|e| e.0);
        if exons.is_empty() {
            continue;
        }
        out.push(RefTranscript {
            id,
            chrom,
            strand,
            exons,
        });
    }
    out.sort_by(|a, b| {
        (a.chrom.as_str(), a.exons.first().map(|e| e.0))
            .cmp(&(b.chrom.as_str(), b.exons.first().map(|e| e.0)))
    });
    Ok(out)
}

fn parse_gtf_attr(attrs: &str, key: &str) -> Option<String> {
    for part in attrs.split(';') {
        let part = part.trim();
        if part.starts_with(key) {
            let v = part
                .strip_prefix(key)?
                .trim()
                .trim_start_matches(' ')
                .strip_prefix('"')?;
            let end = v.find('"')?;
            return Some(v[..end].to_string());
        }
    }
    None
}

/// Result of mapping one guide to the graph (CGuide — guide mapping info).
#[derive(Debug, Clone)]
pub struct GuideInfo {
    pub node_ids: Vec<usize>,
    pub guide_index: usize,
    pub transcript_id: String,
    pub tx_start: u64,
    pub tx_end: u64,
}

/// Find a path of graph node IDs that matches the guide exons (find_guide_pat).
/// Direct port of.
/// Walks sequentially through graph nodes, following child edges for splice junctions
/// and contiguous sub-nodes within the same exon.
pub fn find_guide_pat(
    guide: &RefTranscript,
    graph: &Graph,
    bundle_start: u64,
    bundle_end: u64,
    _strand: char,
    ssdist: u64,
) -> Option<Vec<usize>> {
    let guide_start = guide.exons.first().map(|e| e.0).unwrap_or(0);
    let guide_end = guide.exons.last().map(|e| e.1).unwrap_or(0);
    if guide_end <= bundle_start || guide_start >= bundle_end {
        return None;
    }
    let source_id = graph.source_id;
    let sink_id = graph.sink_id;
    let gno = graph.nodes.len();
    let n_exons = guide.exons.len();
    if n_exons == 0 || gno < 3 {
        return None;
    }

    //: find start node — first interior node overlapping guide's first exon.
    // Walk forward past contiguous sub-nodes that also overlap first exon, choosing the last
    // contiguous chain node as start (the algorithm skips nodes linked by introns to later overlapping nodes).
    let first_exon_start = guide.exons[0].0;
    let first_exon_end = guide.exons[0].1;

    let mut start_i = None;
    let mut i = 1; // skip source (node 0)
    while i < gno.saturating_sub(1) {
        // skip source and sink
        if i == source_id || i == sink_id {
            i += 1;
            continue;
        }
        let node = &graph.nodes[i];
        if node.start < first_exon_end && node.end > first_exon_start {
            // This node overlaps the first guide exon.
            // Walk forward past contiguous sub-nodes that also overlap.
            let mut j = i + 1;
            while j < gno.saturating_sub(1) && j != source_id && j != sink_id {
                let nj = &graph.nodes[j];
                if !(nj.start < first_exon_end && nj.end > first_exon_start) {
                    break;
                }
                if nj.start != graph.nodes[j - 1].end
                    || !graph.nodes[j - 1].children.contains(j)
                    || graph.nodes[j - 1].children.is_empty()
                    || graph.nodes[j - 1]
                        .children
                        .ones()
                        .next()
                        .unwrap_or(usize::MAX)
                        != j
                {
                    // Not contiguous or not first child — this is an intron-linked node.
                    // i=j-1 then break (pick the last contiguous node before the jump).
                    start_i = Some(j.saturating_sub(1));
                    break;
                }
                j += 1;
            }
            if start_i.is_none() {
                start_i = Some(i);
            }
            break;
        }
        i += 1;
    }

    let mut i = start_i?;
    let mut nodes = vec![i];
    let mut j = 0usize; // guide exon index

    //: walk through guide exons, matching to graph nodes.
    while j < n_exons {
        let inode = &graph.nodes[i];
        let ex_end = guide.exons[j].1; // guide exon end (exclusive, 0-based)
        let end_diff = inode.end.abs_diff(ex_end);

        if ex_end + ssdist < inode.end {
            // Guide exon ends before node ends — guide exon is contained in this node.
            //: if last exon, create transfrag and return.
            if j == n_exons - 1 {
                return Some(nodes);
            }
            return Some(nodes); // guide fully consumed
        } else if end_diff <= ssdist {
            // Guide exon end matches node end exactly.
            //: move to next exon, find child matching next exon start.
            if j == n_exons - 1 {
                return Some(nodes);
            }
            j += 1;
            let next_ex_start = guide.exons[j].0;
            // Find child of current node whose start matches next exon start.
            let mut found_child = false;
            for child_id in inode.children.ones() {
                if child_id == source_id || child_id == sink_id {
                    continue;
                }
                if graph.nodes[child_id].start.abs_diff(next_ex_start) <= ssdist {
                    nodes.push(child_id);
                    i = child_id;
                    found_child = true;
                    break;
                }
            }
            if !found_child {
                break;
            }
        } else {
            // Guide exon end > node end — exon spans beyond this node.
            //: walk to contiguous next sub-node.
            // Find a contiguous child node (node.end == child.start, first child).
            let mut walked = false;
            for child_id in inode.children.ones() {
                if child_id == source_id || child_id == sink_id {
                    continue;
                }
                if graph.nodes[child_id].start == inode.end {
                    // Contiguous sub-node — walk into it.
                    nodes.push(child_id);
                    i = child_id;
                    walked = true;
                    break;
                }
            }
            if !walked {
                // No contiguous sub-node. If last exon, store; else fail.
                if j == n_exons - 1 {
                    return Some(nodes);
                }
                break;
            }
        }
    }

    if !nodes.is_empty() && j >= n_exons.saturating_sub(1) {
        return Some(nodes);
    }

    find_guide_pat_fallback(guide, graph, source_id, sink_id, ssdist)
}

fn find_guide_pat_fallback(
    guide: &RefTranscript,
    graph: &Graph,
    source_id: usize,
    sink_id: usize,
    ssdist: u64,
) -> Option<Vec<usize>> {
    fn overlaps_exon(node: &crate::graph::GraphNode, exon: (u64, u64), ssdist: u64) -> bool {
        node.start < exon.1.saturating_add(ssdist) && node.end.saturating_add(ssdist) > exon.0
    }

    fn search(
        graph: &Graph,
        guide: &RefTranscript,
        source_id: usize,
        sink_id: usize,
        ssdist: u64,
        node_id: usize,
        exon_idx: usize,
        path: &mut Vec<usize>,
        seen: &mut HashSet<(usize, usize)>,
    ) -> bool {
        if !seen.insert((node_id, exon_idx)) {
            return false;
        }
        let node = &graph.nodes[node_id];
        let exon = guide.exons[exon_idx];
        if !overlaps_exon(node, exon, ssdist) {
            return false;
        }
        path.push(node_id);
        let exon_end = exon.1;
        let end_matches = node.end.abs_diff(exon_end) <= ssdist || node.end >= exon_end;
        if end_matches {
            if exon_idx + 1 == guide.exons.len() {
                return true;
            }
            let next_exon = guide.exons[exon_idx + 1];
            for child_id in node.children.ones() {
                if child_id == source_id || child_id == sink_id {
                    continue;
                }
                let child = &graph.nodes[child_id];
                if overlaps_exon(child, next_exon, ssdist)
                    && search(
                        graph,
                        guide,
                        source_id,
                        sink_id,
                        ssdist,
                        child_id,
                        exon_idx + 1,
                        path,
                        seen,
                    )
                {
                    return true;
                }
            }
        }
        for child_id in node.children.ones() {
            if child_id == source_id || child_id == sink_id {
                continue;
            }
            let child = &graph.nodes[child_id];
            if child.start.abs_diff(node.end) <= ssdist
                && overlaps_exon(child, exon, ssdist)
                && search(
                    graph,
                    guide,
                    source_id,
                    sink_id,
                    ssdist,
                    child_id,
                    exon_idx,
                    path,
                    seen,
                )
            {
                return true;
            }
        }
        path.pop();
        false
    }

    let mut starts: Vec<usize> = graph
        .nodes
        .iter()
        .enumerate()
        .filter_map(|(nid, node)| {
            if nid == source_id || nid == sink_id {
                return None;
            }
            if overlaps_exon(node, guide.exons[0], ssdist) {
                Some(nid)
            } else {
                None
            }
        })
        .collect();
    starts.sort_by_key(|&nid| graph.nodes[nid].start);
    for start in starts {
        let mut path = Vec::new();
        let mut seen = HashSet::default();
        if search(
            graph,
            guide,
            source_id,
            sink_id,
            ssdist,
            start,
            0,
            &mut path,
            &mut seen,
        ) {
            return Some(path);
        }
    }
    None
}

/// Map reference guides to graph paths, set hardstart/hardend on guide terminal nodes (process_refguides).
pub fn process_refguides(
    graph: &mut Graph,
    guides: &[RefTranscript],
    bundle_start: u64,
    bundle_end: u64,
    strand: char,
    ssdist: u64,
    verbose: bool,
) -> Vec<GuideInfo> {
    let mut mapped = Vec::new();
    for (gi, guide) in guides.iter().enumerate() {
        if let Some(node_ids) =
            find_guide_pat(guide, graph, bundle_start, bundle_end, strand, ssdist)
        {
            let first_nid = *node_ids.first().unwrap_or(&0);
            let last_nid = *node_ids.last().unwrap_or(&0);
            if first_nid < graph.nodes.len() {
                crate::bump_hs!("reference_gtf.rs:397:hardstart");
                graph.nodes[first_nid].hardstart = true;
            }
            if last_nid < graph.nodes.len() {
                crate::bump_hs!("reference_gtf.rs:400:hardend");
                graph.nodes[last_nid].hardend = true;
            }
            mapped.push(GuideInfo {
                node_ids,
                guide_index: gi + 1,
                transcript_id: guide.id.clone(),
                tx_start: guide.exons.first().map(|e| e.0).unwrap_or(0),
                tx_end: guide.exons.last().map(|e| e.1).unwrap_or(0),
            });
        }
    }
    mapped.sort_by(|a, b| b.node_ids.len().cmp(&a.node_ids.len()));
    if verbose && !mapped.is_empty() {
        eprintln!(
            "      [Rustle] Guide mapping: {} guides mapped (bundle {}-{})",
            mapped.len(),
            bundle_start,
            bundle_end
        );
    }
    mapped
}

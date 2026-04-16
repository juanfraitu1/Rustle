//! Exon interval merging (merge_exons(CGene&, GVec<GSeg>&)).
//!
//! Merges a sorted list of (start, end) exons into an existing sorted list,
//! extending overlapping intervals and collapsing any that become redundant.

/// Merge new_exons into gene_exons in place (merge_exons(CGene&, GVec<GSeg>&)).
/// Both slices must be sorted by start. Overlapping or adjacent intervals in gene_exons
/// are extended to include new_exons and then merged so the result remains disjoint.
pub fn merge_exons_into(gene_exons: &mut Vec<(u64, u64)>, new_exons: &[(u64, u64)]) {
    let mut ig = 0usize;
    let mut ie = 0usize;
    while ie < new_exons.len() {
        let (ne_start, ne_end) = new_exons[ie];
        if ig >= gene_exons.len() || ne_end < gene_exons[ig].0 {
            gene_exons.insert(ig, (ne_start, ne_end));
            ie += 1;
            ig += 1;
            continue;
        }
        while ig < gene_exons.len() && ne_start > gene_exons[ig].1 {
            ig += 1;
        }
        if ig < gene_exons.len() {
            if ne_start <= gene_exons[ig].0 {
                gene_exons[ig].0 = ne_start;
            }
            if ne_end >= gene_exons[ig].1 {
                gene_exons[ig].1 = ne_end;
                ig += 1;
                while ig < gene_exons.len() && ne_end >= gene_exons[ig].0 {
                    if gene_exons[ig].1 > ne_end {
                        gene_exons[ig - 1].1 = gene_exons[ig].1;
                    }
                    gene_exons.remove(ig);
                }
            }
            ie += 1;
        }
    }
}

/// Merge overlapping intervals in a single sorted list; returns a new sorted disjoint list.
/// Useful when building gene exons from multiple transcripts (no prior gene list).
pub fn merge_exon_intervals(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
    if exons.is_empty() {
        return Vec::new();
    }
    let mut sorted: Vec<(u64, u64)> = exons.to_vec();
    sorted.sort_unstable_by_key(|e| e.0);
    let mut out = vec![sorted[0]];
    for &(s, e) in sorted.iter().skip(1) {
        let last = out.last_mut().unwrap();
        if s <= last.1 {
            if e > last.1 {
                last.1 = e;
            }
        } else {
            out.push((s, e));
        }
    }
    out
}

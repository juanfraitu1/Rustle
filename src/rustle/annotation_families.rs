//! -G2 annotation -> position-agnostic paralog families.
//! Clusters annotation transcripts into families by SHARED-EXON SEQUENCE similarity
//! (never coordinate proximity), gated by RunConfig.family_exon_similarity.

use anyhow::Result;
use std::path::Path;

/// One annotated transcript = one copy in a (potential) family.
#[derive(Clone, Debug)]
pub struct CopyStructure {
    pub copy_id: String,
    pub chrom: String,
    pub strand: char,
    /// 0-based half-open exons.
    pub exons: Vec<(u64, u64)>,
}

/// A paralog family: >=2 copies merged by shared-exon similarity, with the achieved similarity.
#[derive(Clone, Debug)]
pub struct AnnotationFamily {
    pub family_id: String,
    pub copies: Vec<CopyStructure>,
    pub achieved_min_sim: f64,
    pub achieved_mean_sim: f64,
}

/// Concatenated exon sequence of a copy, in genomic (exon-list) order.
///
/// Fetches each exon's bases from the genome using 0-based half-open coordinates
/// (matching `CopyStructure.exons`). Returns `None` if any exon is entirely out
/// of range for its chromosome. Strand is NOT considered here — all copies are
/// returned in the same forward-genomic convention so similarity is computed on
/// a consistent basis.
pub fn copy_sequence(
    copy: &CopyStructure,
    genome: &crate::genome::GenomeIndex,
) -> Option<Vec<u8>> {
    let mut seq = Vec::new();
    for &(s, e) in &copy.exons {
        let part = genome.fetch_sequence(&copy.chrom, s, e)?;
        seq.extend_from_slice(&part);
    }
    Some(seq)
}

// ── position-agnostic family clustering ──────────────────────────────────────

fn jaccard(a: &std::collections::HashSet<u64>, b: &std::collections::HashSet<u64>) -> f64 {
    if a.is_empty() && b.is_empty() {
        return 1.0;
    }
    let inter = a.intersection(b).count() as f64;
    let union = a.union(b).count() as f64;
    if union == 0.0 { 0.0 } else { inter / union }
}

/// Iterative path-halving union-find: find root of x.
fn uf_find(parent: &mut Vec<usize>, mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]]; // path halving
        x = parent[x];
    }
    x
}

/// Cluster copies into families purely by shared-exon minimizer-Jaccard >= threshold.
/// POSITION-AGNOSTIC: no chrom / coordinate / distance gate.
pub fn cluster_families(
    copies_with_seq: Vec<(CopyStructure, Vec<u8>)>,
    threshold: f64,
) -> Vec<AnnotationFamily> {
    let n = copies_with_seq.len();
    let mins: Vec<std::collections::HashSet<u64>> = copies_with_seq
        .iter()
        .map(|(_, s)| crate::vg_family::family_graph::minimizers(s, 15, 10))
        .collect();

    let mut parent: Vec<usize> = (0..n).collect();

    // Precompute pairwise Jaccard and union-find copies with sim >= threshold.
    let mut sim = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            let s = jaccard(&mins[i], &mins[j]);
            sim[i][j] = s;
            sim[j][i] = s;
            if s >= threshold {
                let ri = uf_find(&mut parent, i);
                let rj = uf_find(&mut parent, j);
                if ri != rj {
                    parent[ri] = rj;
                }
            }
        }
    }

    // Collect components.
    use std::collections::BTreeMap;
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = uf_find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }

    let mut out = Vec::new();
    for (gi, (_, members)) in groups.iter().enumerate() {
        // ≥2-copy invariant: singletons never form a paralog family.
        if members.len() < 2 { continue; }
        let mut pair_sims = Vec::new();
        for a in 0..members.len() {
            for b in (a + 1)..members.len() {
                pair_sims.push(sim[members[a]][members[b]]);
            }
        }
        // Singleton: achieved sims undefined; use 0.0 as sentinel.
        let (min_s, mean_s) = if pair_sims.is_empty() {
            (0.0, 0.0)
        } else {
            let mn = pair_sims.iter().cloned().fold(f64::INFINITY, f64::min);
            let me = pair_sims.iter().sum::<f64>() / pair_sims.len() as f64;
            (mn, me)
        };
        out.push(AnnotationFamily {
            family_id: format!("FAM{}", gi),
            copies: members.iter().map(|&m| copies_with_seq[m].0.clone()).collect(),
            achieved_min_sim: min_s,
            achieved_mean_sim: mean_s,
        });
    }
    out
}

/// Cluster within each strand group separately (same-strand is the only structural precondition,
/// since build_family_graph bails on mixed strands), then concat. Position-agnostic within strand.
pub fn families_from_grouped(
    copies_with_seq: Vec<(CopyStructure, Vec<u8>)>,
    threshold: f64,
) -> Vec<AnnotationFamily> {
    use std::collections::BTreeMap;
    let mut by_strand: BTreeMap<char, Vec<(CopyStructure, Vec<u8>)>> = BTreeMap::new();
    for cs in copies_with_seq {
        by_strand.entry(cs.0.strand).or_default().push(cs);
    }
    let mut out = Vec::new();
    for (_, group) in by_strand {
        out.extend(cluster_families(group, threshold));
    }
    // Reindex family_ids to be globally unique across all strand groups.
    // cluster_families restarts at "FAM0" per call, so + and - groups can collide.
    for (i, fam) in out.iter_mut().enumerate() {
        fam.family_id = format!("FAM{}", i);
    }
    out
}

/// Top-level: load -G2, fetch exon sequences from genome, cluster (position-agnostic, same-strand).
pub fn build_families<P: AsRef<std::path::Path>>(
    g2_path: P,
    genome: &crate::genome::GenomeIndex,
    threshold: f64,
) -> anyhow::Result<Vec<AnnotationFamily>> {
    let copies = load_copies(g2_path)?;
    let mut with_seq = Vec::new();
    for c in copies {
        if let Some(s) = copy_sequence(&c, genome) {
            with_seq.push((c, s));
        }
    }
    Ok(families_from_grouped(with_seq, threshold))
}

/// Load -G2 GTF into copy structures (one per transcript).
pub fn load_copies<P: AsRef<Path>>(path: P) -> Result<Vec<CopyStructure>> {
    let refs = crate::reference_gtf::parse_reference_gtf(path)?;
    Ok(refs
        .into_iter()
        .map(|r| CopyStructure {
            copy_id: r.id,
            chrom: r.chrom,
            strand: r.strand,
            exons: r.exons,
        })
        .collect())
}

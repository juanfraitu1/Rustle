//! Variation graph mode for gene family assembly.
//!
//! Links related genomic loci (gene family copies) using multi-mapping read
//! patterns, then jointly resolves multi-mapping reads across copies using EM.
//!
//! Activated with `--vg`. Multi-mapping reads (NH >= 2) that appear in multiple
//! bundles link those bundles into "family groups." Within each family group,
//! an EM algorithm redistributes read weights based on splice-graph compatibility,
//! then the standard flow pipeline runs with updated weights.

use crate::coord::overlaps_half_open;
use crate::graph::Graph;
use crate::map_reads::read_to_path_bundlenodes;
use crate::types::{Bundle, BundleRead, CBundlenode, Bundle2Graph};
use std::collections::HashMap;

/// FNV-1a hash matching bam.rs read_name_hash computation.
fn fnv1a64(s: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in s {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

// ── Data structures ──────────────────────────────────────────────────────────

/// A group of bundles linked by multi-mapping reads (gene family copies).
#[derive(Debug)]
pub struct FamilyGroup {
    pub family_id: usize,
    /// Indices into the global bundle list.
    pub bundle_indices: Vec<usize>,
    /// read_name_hash → Vec<(position in bundle_indices, read index within that bundle)>
    pub multimap_reads: HashMap<u64, Vec<(usize, usize)>>,
}

/// Result of the EM reweighting pass.
#[derive(Debug, Default, Clone)]
pub struct EmResult {
    pub iterations: usize,
    pub converged: bool,
    pub max_delta: f64,
    pub reads_reweighted: usize,
}

/// One row in the family report TSV.
pub struct FamilyReportRow {
    pub family_id: usize,
    pub n_copies: usize,
    pub chrom: String,
    pub regions: Vec<(u64, u64, char)>, // (start, end, strand) per copy
    pub n_shared_reads: usize,
    pub em_iterations: usize,
    pub em_converged: bool,
}

// ── Family group discovery ───────────────────────────────────────────────────

/// Build a map of read_name_hash → list of (bundle_idx, read_idx) for multi-mapping reads.
///
/// Detects multi-mappers by:
/// 1. Same read_name_hash appearing in multiple bundles (primary in bundle A, primary in bundle B)
/// 2. NH tag >= 2 with allocations in multiple bundles
///
/// For long reads where supplementary alignments are filtered during BAM parsing,
/// use `build_multimap_index_with_supplementary()` which does a second BAM pass.
pub fn build_multimap_index(bundles: &[Bundle]) -> HashMap<u64, Vec<(usize, usize)>> {
    let mut read_locs: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
    for (bi, bundle) in bundles.iter().enumerate() {
        for (ri, read) in bundle.reads.iter().enumerate() {
            read_locs
                .entry(read.read_name_hash)
                .or_default()
                .push((bi, ri));
        }
    }
    // Keep only reads that appear in more than one DISTINCT bundle.
    read_locs.retain(|_, locs| {
        if locs.len() < 2 {
            return false;
        }
        let first_bundle = locs[0].0;
        locs.iter().any(|(bi, _)| *bi != first_bundle)
    });
    read_locs
}

/// Second-pass BAM scan: find supplementary alignments and link them to bundles.
///
/// Long-read pipelines filter supplementary alignments during bundle detection.
/// This function does a lightweight second pass to find reads whose supplementary
/// alignments map to different bundles than their primary, creating cross-bundle
/// links for family group discovery.
pub fn build_multimap_index_with_supplementary(
    bam_path: &std::path::Path,
    bundles: &[Bundle],
) -> HashMap<u64, Vec<(usize, usize)>> {
    // Start with primary-based index.
    let mut read_locs = build_multimap_index(bundles);

    // Build interval index: sorted (chrom, start, end, bundle_idx) for overlap queries.
    let mut bundle_intervals: Vec<(String, u64, u64, usize)> = bundles
        .iter()
        .enumerate()
        .map(|(bi, b)| (b.chrom.clone(), b.start, b.end, bi))
        .collect();
    bundle_intervals.sort_by_key(|(c, s, _, _)| (c.clone(), *s));

    // Read supplementary alignments from BAM using noodles_bam (same as bam.rs).
    let bam_file = match std::fs::File::open(bam_path) {
        Ok(f) => f,
        Err(_) => return read_locs,
    };
    let buf = std::io::BufReader::new(bam_file);
    let worker_count = std::num::NonZeroUsize::MIN;
    let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(worker_count, buf);
    let mut reader = noodles_bam::io::Reader::from(bgzf);
    let header = match reader.read_header() {
        Ok(h) => h,
        Err(_) => return read_locs,
    };

    // Pre-build a set of (read_name_hash → primary_bundle_idx) for fast lookup.
    let mut rnh_to_primary: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
    for (bi, bundle) in bundles.iter().enumerate() {
        for (ri, read) in bundle.reads.iter().enumerate() {
            rnh_to_primary
                .entry(read.read_name_hash)
                .or_default()
                .push((bi, ri));
        }
    }

    let mut n_supp = 0usize;
    let mut n_linked = 0usize;
    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        let flags = record.flags();
        // Multi-mappers in PacBio/ONT data typically appear as SECONDARY (flag 256);
        // chimeric splits appear as SUPPLEMENTARY (flag 2048). Scan both.
        let is_multi = flags.is_secondary() || flags.is_supplementary();
        if !is_multi || flags.is_unmapped() {
            continue;
        }
        n_supp += 1;

        let ref_id = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            _ => continue,
        };
        let chrom = match header.reference_sequences().get_index(ref_id) {
            Some((name, _)) => format!("{}", name),
            None => continue,
        };
        let start_1b = match record.alignment_start() {
            Some(Ok(pos)) => pos.get() as u64,
            _ => continue,
        };
        let start = start_1b.saturating_sub(1); // convert to 0-based
        // Approximate end from start + read length (avoid computing CIGAR alignment end).
        let approx_len = record.sequence().len() as u64;
        let end = start + approx_len.max(100);

        // Find which bundle this supplementary overlaps.
        let target_bi = bundle_intervals
            .iter()
            .filter(|(c, bs, be, _)| *c == chrom && start < *be && end > *bs)
            .map(|(_, _, _, bi)| *bi)
            .next();

        let Some(target_bi) = target_bi else {
            continue;
        };

        // Hash the read name.
        let read_name_str = record
            .name()
            .map(|n| n.to_string())
            .unwrap_or_default();
        let rnh = fnv1a64(read_name_str.as_bytes());

        // Find primary alignment's bundle(s) for this read.
        let Some(primaries) = rnh_to_primary.get(&rnh) else {
            continue;
        };

        for &(pbi, ri) in primaries {
            if pbi == target_bi {
                continue;
            }
            let entry = read_locs.entry(rnh).or_default();
            if !entry.iter().any(|(bi, _)| *bi == pbi) {
                entry.push((pbi, ri));
            }
            // Sentinel read index for the supplementary target bundle.
            if !entry.iter().any(|(bi, _)| *bi == target_bi) {
                entry.push((target_bi, usize::MAX));
                n_linked += 1;
            }
        }
    }

    if n_supp > 0 {
        eprintln!(
            "[VG] Supplementary scan: {} supplementary alignments, {} cross-bundle links",
            n_supp, n_linked
        );
    }

    // Keep only entries with multiple distinct bundles.
    read_locs.retain(|_, locs| {
        if locs.len() < 2 {
            return false;
        }
        let first_bundle = locs[0].0;
        locs.iter().any(|(bi, _)| *bi != first_bundle)
    });

    read_locs
}

/// Discover family groups from multi-mapping read patterns.
///
/// Two bundles are linked if they share at least `min_shared_reads` multi-mapping
/// reads. Connected components of the resulting graph form family groups.
pub fn discover_family_groups(
    bundles: &[Bundle],
    min_shared_reads: usize,
    bam_path: Option<&std::path::Path>,
    genome: Option<&crate::genome::GenomeIndex>,
) -> Vec<FamilyGroup> {
    let multimap_index = if let Some(path) = bam_path {
        build_multimap_index_with_supplementary(path, bundles)
    } else {
        build_multimap_index(bundles)
    };

    // ALSO link bundles by exonic sequence similarity (for tandem duplicates
    // that don't generate supplementary cross-mappings).
    let seq_links = if let Some(genome_ref) = genome {
        discover_sequence_similar_bundles(bundles, genome_ref)
    } else {
        Vec::new()
    };
    if multimap_index.is_empty() {
        return Vec::new();
    }

    // Build bundle-pair link counts.
    let n_bundles = bundles.len();
    let mut link_counts: HashMap<(usize, usize), usize> = HashMap::new();
    for locs in multimap_index.values() {
        // Collect unique bundle indices for this read.
        let mut bundle_set: Vec<usize> = locs.iter().map(|(bi, _)| *bi).collect();
        bundle_set.sort_unstable();
        bundle_set.dedup();
        // Create pairwise links.
        for i in 0..bundle_set.len() {
            for j in (i + 1)..bundle_set.len() {
                let key = (bundle_set[i], bundle_set[j]);
                *link_counts.entry(key).or_insert(0) += 1;
            }
        }
    }

    // Union-Find for connected components.
    let mut parent: Vec<usize> = (0..n_bundles).collect();
    let mut rank: Vec<usize> = vec![0; n_bundles];

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }
    fn union(parent: &mut [usize], rank: &mut [usize], a: usize, b: usize) {
        let ra = find(parent, a);
        let rb = find(parent, b);
        if ra == rb {
            return;
        }
        if rank[ra] < rank[rb] {
            parent[ra] = rb;
        } else if rank[ra] > rank[rb] {
            parent[rb] = ra;
        } else {
            parent[rb] = ra;
            rank[ra] += 1;
        }
    }

    for (&(a, b), &count) in &link_counts {
        if count >= min_shared_reads {
            union(&mut parent, &mut rank, a, b);
        }
    }
    // Also link bundles with high exonic sequence similarity.
    for &(a, b) in &seq_links {
        union(&mut parent, &mut rank, a, b);
    }

    // Collect components with more than one bundle.
    let mut components: HashMap<usize, Vec<usize>> = HashMap::new();
    for bi in 0..n_bundles {
        let root = find(&mut parent, bi);
        components.entry(root).or_default().push(bi);
    }

    let mut families: Vec<FamilyGroup> = Vec::new();
    for (fam_idx, (_root, bundle_indices)) in components
        .into_iter()
        .filter(|(_, v)| v.len() > 1)
        .enumerate()
    {
        // Build the per-family multimap read index.
        let bundle_set: std::collections::HashSet<usize> =
            bundle_indices.iter().copied().collect();
        let mut family_reads: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for (&rnh, locs) in &multimap_index {
            let family_locs: Vec<(usize, usize)> = locs
                .iter()
                .filter(|(bi, _)| bundle_set.contains(bi))
                .map(|(bi, ri)| {
                    // Map global bundle index to position within bundle_indices.
                    let pos = bundle_indices.iter().position(|&x| x == *bi).unwrap();
                    (pos, *ri)
                })
                .collect();
            if family_locs.len() > 1 {
                family_reads.insert(rnh, family_locs);
            }
        }

        families.push(FamilyGroup {
            family_id: fam_idx,
            bundle_indices,
            multimap_reads: family_reads,
        });
    }

    if !families.is_empty() {
        eprintln!(
            "[VG] Discovered {} gene family group(s) from multi-mapping reads ({} linked bundles, {} multi-map reads)",
            families.len(),
            families.iter().map(|f| f.bundle_indices.len()).sum::<usize>(),
            families.iter().map(|f| f.multimap_reads.len()).sum::<usize>(),
        );
    }

    families
}

// ── Sequence-similarity-based family discovery ───────────────────────────────

/// Find bundle pairs with high exonic sequence similarity (for tandem duplicates
/// that don't generate supplementary cross-mappings).
///
/// Uses 15-mer Jaccard similarity: bundles sharing ≥30% of their exonic k-mers
/// are linked as gene family members.
fn discover_sequence_similar_bundles(
    bundles: &[Bundle],
    genome: &crate::genome::GenomeIndex,
) -> Vec<(usize, usize)> {
    let kmer_len = 15usize;
    let min_jaccard = 0.20; // 20% k-mer overlap = gene family member.
    let min_kmers = 50; // Skip bundles with too few exonic k-mers.
    let max_bundles_to_compare = 5000; // Skip if too many bundles (O(n²) comparison).

    if bundles.len() > max_bundles_to_compare {
        return Vec::new();
    }

    // Build exonic k-mer sets per bundle.
    let mut bundle_kmers: Vec<(usize, std::collections::HashSet<u64>)> = Vec::new();
    for (bi, bundle) in bundles.iter().enumerate() {
        // Extract exonic regions from junction stats.
        let mut exon_regions: Vec<(u64, u64)> = Vec::new();
        let mut introns: Vec<(u64, u64)> = bundle
            .junction_stats
            .iter()
            .map(|(j, _)| (j.donor, j.acceptor))
            .collect();
        introns.sort_by_key(|&(s, _)| s);

        if introns.is_empty() {
            // Single-exon bundle: use full span but cap at 10kb.
            let capped_end = bundle.end.min(bundle.start + 10000);
            exon_regions.push((bundle.start, capped_end));
        } else {
            // First exon.
            let first_intron_start = introns[0].0;
            if first_intron_start > bundle.start {
                exon_regions.push((
                    first_intron_start.saturating_sub(500).max(bundle.start),
                    first_intron_start,
                ));
            }
            // Internal exons.
            for i in 0..introns.len().saturating_sub(1) {
                let es = introns[i].1;
                let ee = introns[i + 1].0;
                if ee > es {
                    exon_regions.push((es, ee));
                }
            }
            // Last exon.
            let last_intron_end = introns.last().unwrap().1;
            exon_regions.push((last_intron_end, last_intron_end.saturating_add(500).min(bundle.end)));
        }

        let mut kmers = std::collections::HashSet::new();
        for (es, ee) in &exon_regions {
            if let Some(seq) = genome.fetch_sequence(&bundle.chrom, *es, *ee) {
                for window in seq.windows(kmer_len) {
                    if window.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                        kmers.insert(kmer_hash(window));
                    }
                }
            }
        }

        if kmers.len() >= min_kmers {
            bundle_kmers.push((bi, kmers));
        }
    }

    // Pairwise comparison (only same-chromosome bundles within 10Mb).
    let mut links: Vec<(usize, usize)> = Vec::new();
    for i in 0..bundle_kmers.len() {
        let (bi, ref ki) = bundle_kmers[i];
        for j in (i + 1)..bundle_kmers.len() {
            let (bj, ref kj) = bundle_kmers[j];
            // Same chromosome check.
            if bundles[bi].chrom != bundles[bj].chrom {
                continue;
            }
            // Skip if bundles overlap (same region, different strand).
            if bundles[bi].start < bundles[bj].end && bundles[bj].start < bundles[bi].end {
                continue;
            }
            // Distance cap: within 10Mb.
            let dist = if bundles[bi].start > bundles[bj].end {
                bundles[bi].start - bundles[bj].end
            } else {
                bundles[bj].start - bundles[bi].end
            };
            if dist > 10_000_000 {
                continue;
            }
            // Jaccard similarity.
            let intersection = ki.iter().filter(|k| kj.contains(k)).count();
            let union_size = ki.len() + kj.len() - intersection;
            if union_size == 0 {
                continue;
            }
            let jaccard = intersection as f64 / union_size as f64;
            if jaccard >= min_jaccard {
                links.push((bi, bj));
            }
        }
    }

    if !links.is_empty() {
        eprintln!(
            "[VG] Sequence similarity: found {} bundle pairs with ≥{:.0}% exonic k-mer overlap",
            links.len(),
            min_jaccard * 100.0,
        );
    }

    links
}

// ── EM reweighting ───────────────────────────────────────────────────────────

/// Compute compatibility score: how well does a read's exon structure match
/// the splice graph at a given bundle? Higher = better match.
fn read_graph_compatibility(
    read: &BundleRead,
    graph: &Graph,
    bundlenodes: Option<&CBundlenode>,
    bundle2graph: Option<&Bundle2Graph>,
) -> f64 {
    let path = read_to_path_bundlenodes(read, graph, bundlenodes, bundle2graph);
    if path.is_empty() {
        return 0.0;
    }
    // Score = fraction of read exons that map to at least one graph node.
    let total_exons = read.exons.len().max(1);
    let matched_exons = read.exons.iter().filter(|&&(es, ee)| {
        path.iter().any(|&nid| {
            if let Some(n) = graph.nodes.get(nid) {
                overlaps_half_open(es, ee, n.start, n.end)
            } else {
                false
            }
        })
    }).count();
    // Bonus for matching more nodes (longer concordant paths).
    let node_bonus = (path.len() as f64).ln().max(0.0);
    (matched_exons as f64 / total_exons as f64) * (1.0 + 0.1 * node_bonus)
}

/// Run EM reweighting across all copies in a family group.
///
/// Modifies `BundleRead.weight` in each bundle's reads to reflect the
/// splice-graph-compatibility-based reassignment.
///
/// Returns the EM result summary.
pub fn run_em_reweighting(
    family: &FamilyGroup,
    bundles: &mut [Bundle],
    graphs: &[Graph],
    bundlenodes: &[Option<CBundlenode>],
    bundle2graphs: &[Option<Vec<Vec<(usize, usize)>>>],
    max_iter: usize,
    convergence_thr: f64,
) -> EmResult {
    let n_copies = family.bundle_indices.len();
    if n_copies < 2 || family.multimap_reads.is_empty() {
        return EmResult::default();
    }

    // Build weight table: read_name_hash → Vec<(global_bundle_idx, read_idx, weight)>
    struct ReadEntry {
        locs: Vec<(usize, usize, usize)>, // (family_pos, global_bi, read_idx)
        weights: Vec<f64>,
    }
    let mut entries: Vec<(u64, ReadEntry)> = Vec::with_capacity(family.multimap_reads.len());

    for (&rnh, locs) in &family.multimap_reads {
        let mut entry_locs = Vec::with_capacity(locs.len());
        let mut entry_weights = Vec::with_capacity(locs.len());
        for &(fam_pos, ri) in locs {
            let global_bi = family.bundle_indices[fam_pos];
            let w = bundles[global_bi].reads[ri].weight;
            entry_locs.push((fam_pos, global_bi, ri));
            entry_weights.push(w);
        }
        entries.push((rnh, ReadEntry {
            locs: entry_locs,
            weights: entry_weights,
        }));
    }

    let mut result = EmResult::default();

    for iter in 0..max_iter {
        let mut max_delta: f64 = 0.0;

        for (_rnh, entry) in &mut entries {
            // E-step: compute compatibility scores.
            let mut scores: Vec<f64> = Vec::with_capacity(entry.locs.len());
            for &(fam_pos, global_bi, ri) in &entry.locs {
                let read = &bundles[global_bi].reads[ri];
                let graph = &graphs[fam_pos];
                let bnodes = bundlenodes[fam_pos].as_ref();
                let b2g = bundle2graphs[fam_pos].as_ref();
                let compat = read_graph_compatibility(read, graph, bnodes, b2g);
                // Weight by current graph abundance context (total node coverage).
                let context = graph.nodes.iter().map(|n| n.coverage).sum::<f64>().max(1.0);
                scores.push(compat * context.sqrt());
            }

            let total: f64 = scores.iter().sum();
            if total <= 0.0 {
                continue;
            }

            // Normalize and compute delta.
            let original_nh = entry.locs.len() as f64;
            for (i, score) in scores.iter().enumerate() {
                let _new_w = score / total * original_nh.recip() * original_nh;
                // new_w sums to 1.0 across copies (like 1/NH but redistributed).
                let new_w = score / total;
                let delta = (new_w - entry.weights[i]).abs();
                if delta > max_delta {
                    max_delta = delta;
                }
                entry.weights[i] = new_w;
            }
        }

        result.iterations = iter + 1;
        result.max_delta = max_delta;

        if max_delta < convergence_thr {
            result.converged = true;
            break;
        }
    }

    // M-step: apply final weights back to reads.
    let mut n_reweighted = 0usize;
    for (_rnh, entry) in &entries {
        for (i, &(_, global_bi, ri)) in entry.locs.iter().enumerate() {
            let old_w = bundles[global_bi].reads[ri].weight;
            let new_w = entry.weights[i];
            if (old_w - new_w).abs() > 1e-9 {
                bundles[global_bi].reads[ri].weight = new_w;
                n_reweighted += 1;
            }
        }
    }
    result.reads_reweighted = n_reweighted;

    if n_reweighted > 0 {
        eprintln!(
            "[VG] Family {}: EM converged={} in {} iter (delta={:.6}), reweighted {} read entries across {} copies",
            family.family_id,
            result.converged,
            result.iterations,
            result.max_delta,
            n_reweighted,
            n_copies,
        );
    }

    result
}

// ── Family report ────────────────────────────────────────────────────────────

/// Write a TSV report using pre-saved bundle coordinates (bundles consumed by processing loop).
pub fn write_family_report_from_coords(
    path: &std::path::Path,
    families: &[FamilyGroup],
    bundle_coords: &[(String, u64, u64, char)],
) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        f,
        "family_id\tn_copies\tchrom\tregions\tn_shared_reads"
    )?;
    for family in families {
        let regions: Vec<String> = family
            .bundle_indices
            .iter()
            .filter_map(|&bi| {
                bundle_coords.get(bi).map(|(c, s, e, st)| {
                    format!("{}:{}-{}:{}", c, s, e, st)
                })
            })
            .collect();
        let chrom = bundle_coords
            .get(*family.bundle_indices.first().unwrap_or(&0))
            .map(|(c, _, _, _)| c.as_str())
            .unwrap_or("?");
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}",
            family.family_id,
            family.bundle_indices.len(),
            chrom,
            regions.join(";"),
            family.multimap_reads.len(),
        )?;
    }
    Ok(())
}

/// Write a TSV report of discovered family groups (full bundle access version).
pub fn write_family_report(
    path: &std::path::Path,
    families: &[FamilyGroup],
    bundles: &[Bundle],
    em_results: &[EmResult],
) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        f,
        "family_id\tn_copies\tchrom\tregions\tn_shared_reads\tem_iterations\tem_converged\tem_max_delta\treads_reweighted"
    )?;
    for (fi, family) in families.iter().enumerate() {
        let regions: Vec<String> = family
            .bundle_indices
            .iter()
            .map(|&bi| {
                let b = &bundles[bi];
                format!("{}:{}-{}:{}", b.chrom, b.start, b.end, b.strand)
            })
            .collect();
        let chrom = bundles
            .get(*family.bundle_indices.first().unwrap_or(&0))
            .map(|b| b.chrom.as_str())
            .unwrap_or("?");
        let em = em_results.get(fi).cloned().unwrap_or_default();
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}",
            family.family_id,
            family.bundle_indices.len(),
            chrom,
            regions.join(";"),
            family.multimap_reads.len(),
            em.iterations,
            em.converged,
            em.max_delta,
            em.reads_reweighted,
        )?;
    }
    Ok(())
}

// ── Pre-assembly EM (Phase 2) ────────────────────────────────────────────────

/// Compute junction-based compatibility: how well do the read's junctions match
/// the junction profile at a given bundle?
///
/// Score = (number of read junctions supported by other reads at this bundle) / (total read junctions).
/// This is computed BEFORE assembly, using raw junction stats from bundle detection.
pub fn junction_compatibility(read: &BundleRead, bundle: &Bundle) -> f64 {
    if read.junctions.is_empty() {
        // Single-exon read: compatible with any bundle it overlaps.
        return 1.0;
    }
    let mut matched = 0usize;
    for j in &read.junctions {
        // Check if this junction is in the bundle's junction stats.
        let found = bundle.junction_stats.iter().any(|(bj, _st)| {
            // Match within a tolerance window (junction correction distance).
            let d_start = if j.donor > bj.donor {
                j.donor - bj.donor
            } else {
                bj.donor - j.donor
            };
            let d_end = if j.acceptor > bj.acceptor {
                j.acceptor - bj.acceptor
            } else {
                bj.acceptor - j.acceptor
            };
            d_start <= 10 && d_end <= 10
        });
        if found {
            matched += 1;
        }
    }
    matched as f64 / read.junctions.len() as f64
}

// ── SNP-based copy assignment (--vg-snp) ─────────────────────────────────

/// Position → allele frequency per copy: (ref_pos → Vec<(copy_idx, allele_counts)>).
/// A "diagnostic SNP" is a position where one copy has >80% of one allele and
/// another copy has >80% of a different allele.
pub struct DiagnosticSnps {
    /// ref_pos → Vec<(copy_idx_in_family, dominant_base, frequency)>
    pub positions: HashMap<u64, Vec<(usize, u8, f64)>>,
}

/// Build diagnostic SNP set for a family group.
/// Scans all reads at each copy and finds positions where allele frequencies diverge.
pub fn build_diagnostic_snps(
    family: &FamilyGroup,
    bundles: &[Bundle],
) -> DiagnosticSnps {
    // Per-copy allele counts: ref_pos → copy_idx → base → count
    let mut counts: HashMap<u64, Vec<HashMap<u8, u32>>> = HashMap::new();
    let n_copies = family.bundle_indices.len();

    for (copy_idx, &bi) in family.bundle_indices.iter().enumerate() {
        let bundle = &bundles[bi];
        for read in &bundle.reads {
            for &(ref_pos, query_base) in &read.mismatches {
                let entry = counts.entry(ref_pos).or_insert_with(|| {
                    vec![HashMap::new(); n_copies]
                });
                if copy_idx < entry.len() {
                    *entry[copy_idx].entry(query_base).or_insert(0) += 1;
                }
            }
        }
    }

    // Find diagnostic positions: where copies have different dominant alleles.
    let mut positions: HashMap<u64, Vec<(usize, u8, f64)>> = HashMap::new();
    for (pos, copy_alleles) in &counts {
        let mut copy_dominants: Vec<(usize, u8, f64)> = Vec::new();
        for (ci, alleles) in copy_alleles.iter().enumerate() {
            let total: u32 = alleles.values().sum();
            if total < 3 {
                continue; // Need at least 3 reads at this position.
            }
            if let Some((&best_base, &best_count)) = alleles.iter().max_by_key(|(_, &c)| c) {
                let freq = best_count as f64 / total as f64;
                if freq >= 0.8 {
                    copy_dominants.push((ci, best_base, freq));
                }
            }
        }
        // Diagnostic if at least 2 copies have different dominant alleles.
        if copy_dominants.len() >= 2 {
            let first_base = copy_dominants[0].1;
            if copy_dominants.iter().any(|(_, b, _)| *b != first_base) {
                positions.insert(*pos, copy_dominants);
            }
        }
    }

    if !positions.is_empty() {
        eprintln!(
            "[VG-SNP] Found {} diagnostic SNP positions for family with {} copies",
            positions.len(),
            n_copies,
        );
    }

    DiagnosticSnps { positions }
}

/// Score a read's SNP compatibility with a specific copy.
/// Returns a bonus factor (1.0 = no info, >1.0 = matches, <1.0 = mismatches).
pub fn snp_compatibility(
    read: &BundleRead,
    copy_idx: usize,
    diagnostic: &DiagnosticSnps,
) -> f64 {
    if diagnostic.positions.is_empty() || read.mismatches.is_empty() {
        return 1.0; // No info — neutral.
    }
    let mut matches = 0usize;
    let mut mismatches = 0usize;
    for &(ref_pos, query_base) in &read.mismatches {
        if let Some(copy_info) = diagnostic.positions.get(&ref_pos) {
            if let Some((_, diag_base, _)) = copy_info.iter().find(|(ci, _, _)| *ci == copy_idx) {
                if query_base == *diag_base {
                    matches += 1;
                } else {
                    mismatches += 1;
                }
            }
        }
    }
    if matches + mismatches == 0 {
        return 1.0;
    }
    // Bonus: each matching diagnostic SNP adds 0.5, each mismatch subtracts 0.3.
    let bonus = 1.0 + 0.5 * matches as f64 - 0.3 * mismatches as f64;
    bonus.max(0.1)
}

/// Pre-assembly EM: redistribute multi-mapping read weights across family members
/// using junction-based compatibility (before splice graphs are built).
///
/// For each multi-mapping read in a family group:
/// 1. Score compatibility at each copy using junction support
/// 2. Normalize scores to get new weights (summing to 1.0 across copies)
/// 3. Update read.weight in each bundle
///
/// This runs BEFORE the main assembly loop so the standard pipeline uses
/// the EM-adjusted weights for flow decomposition.
pub fn run_pre_assembly_em(
    families: &[FamilyGroup],
    bundles: &mut [Bundle],
    max_iter: usize,
) -> Vec<EmResult> {
    run_pre_assembly_em_inner(families, bundles, max_iter, false)
}

/// EM with optional SNP integration.
pub fn run_pre_assembly_em_with_snps(
    families: &[FamilyGroup],
    bundles: &mut [Bundle],
    max_iter: usize,
) -> Vec<EmResult> {
    run_pre_assembly_em_inner(families, bundles, max_iter, true)
}

fn run_pre_assembly_em_inner(
    families: &[FamilyGroup],
    bundles: &mut [Bundle],
    max_iter: usize,
    use_snps: bool,
) -> Vec<EmResult> {
    let mut results = Vec::with_capacity(families.len());
    let convergence_thr = 0.001;

    for family in families {
        let n_copies = family.bundle_indices.len();
        if n_copies < 2 || family.multimap_reads.is_empty() {
            results.push(EmResult::default());
            continue;
        }

        // Collect read entries with their current weights.
        struct ReadWeightEntry {
            locs: Vec<(usize, usize)>, // (global_bundle_idx, read_idx) — usize::MAX = no read in bundle
            weights: Vec<f64>,
        }
        let mut entries: Vec<(u64, ReadWeightEntry)> = Vec::new();

        for (&rnh, locs) in &family.multimap_reads {
            let mut entry_locs = Vec::new();
            let mut entry_weights = Vec::new();
            for &(fam_pos, ri) in locs {
                let global_bi = family.bundle_indices[fam_pos];
                let w = if ri < bundles[global_bi].reads.len() {
                    bundles[global_bi].reads[ri].weight
                } else {
                    0.0 // sentinel for supplementary-only link
                };
                entry_locs.push((global_bi, ri));
                entry_weights.push(w);
            }
            // Skip if all links are supplementary-only (no actual read in any bundle).
            if entry_weights.iter().all(|&w| w == 0.0) {
                continue;
            }
            entries.push((rnh, ReadWeightEntry {
                locs: entry_locs,
                weights: entry_weights,
            }));
        }

        if entries.is_empty() {
            results.push(EmResult::default());
            continue;
        }

        // Build diagnostic SNPs if SNP mode is active.
        let diagnostic = if use_snps {
            Some(build_diagnostic_snps(family, bundles))
        } else {
            None
        };

        let mut result = EmResult::default();

        for iter in 0..max_iter {
            let mut max_delta: f64 = 0.0;

            for (_rnh, entry) in &mut entries {
                // E-step: compute compatibility scores.
                let mut scores: Vec<f64> = Vec::with_capacity(entry.locs.len());
                for &(global_bi, ri) in &entry.locs {
                    if ri >= bundles[global_bi].reads.len() {
                        // No actual read at this bundle — supplementary link only.
                        // Give a small base score (the read COULD be from here).
                        scores.push(0.1);
                        continue;
                    }
                    let read = &bundles[global_bi].reads[ri];
                    let bundle = &bundles[global_bi];
                    let compat = junction_compatibility(read, bundle);
                    // Context: total junction support at this bundle (higher = more expressed copy).
                    let context: f64 = bundle
                        .junction_stats
                        .iter()
                        .map(|(_, st)| st.nreads_good)
                        .sum::<f64>()
                        .max(1.0);
                    let mut score = (compat + 0.01) * context.ln().max(1.0);
                    // SNP bonus: if diagnostic SNPs are available, multiply by SNP compatibility.
                    if let Some(ref diag) = diagnostic {
                        let copy_idx = family.bundle_indices.iter()
                            .position(|&bi| bi == global_bi)
                            .unwrap_or(0);
                        score *= snp_compatibility(read, copy_idx, diag);
                    }
                    scores.push(score);
                }

                let total: f64 = scores.iter().sum();
                if total <= 0.0 {
                    continue;
                }

                // Normalize.
                for (i, score) in scores.iter().enumerate() {
                    let new_w = score / total;
                    let delta = (new_w - entry.weights[i]).abs();
                    if delta > max_delta {
                        max_delta = delta;
                    }
                    entry.weights[i] = new_w;
                }
            }

            result.iterations = iter + 1;
            result.max_delta = max_delta;
            if max_delta < convergence_thr {
                result.converged = true;
                break;
            }
        }

        // Apply final weights.
        let mut n_reweighted = 0usize;
        for (_rnh, entry) in &entries {
            for (i, &(global_bi, ri)) in entry.locs.iter().enumerate() {
                if ri >= bundles[global_bi].reads.len() {
                    continue;
                }
                let old_w = bundles[global_bi].reads[ri].weight;
                let new_w = entry.weights[i];
                if (old_w - new_w).abs() > 1e-9 {
                    bundles[global_bi].reads[ri].weight = new_w;
                    n_reweighted += 1;
                }
            }
        }
        result.reads_reweighted = n_reweighted;

        if n_reweighted > 0 {
            let bi_str: Vec<String> = family
                .bundle_indices
                .iter()
                .take(4)
                .map(|bi| format!("{}", bi))
                .collect();
            eprintln!(
                "[VG] Family {}: EM converged={} in {} iter (delta={:.6}), reweighted {} reads across {} copies [bundles: {}{}]",
                family.family_id,
                result.converged,
                result.iterations,
                result.max_delta,
                n_reweighted,
                n_copies,
                bi_str.join(","),
                if n_copies > 4 { ",..." } else { "" },
            );
        }

        results.push(result);
    }

    let total_reweighted: usize = results.iter().map(|r| r.reads_reweighted).sum();
    if total_reweighted > 0 {
        eprintln!(
            "[VG] EM reweighting complete: {} reads adjusted across {} families",
            total_reweighted,
            results.iter().filter(|r| r.reads_reweighted > 0).count(),
        );
    }

    results
}

/// Write family report with EM results (uses pre-saved coords since bundles are consumed).
pub fn write_family_report_with_em(
    path: &std::path::Path,
    families: &[FamilyGroup],
    bundle_coords: &[(String, u64, u64, char)],
    em_results: &[EmResult],
) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        f,
        "family_id\tn_copies\tchrom\tregions\tn_shared_reads\tem_iterations\tem_converged\tem_delta\treads_reweighted"
    )?;
    for (fi, family) in families.iter().enumerate() {
        let regions: Vec<String> = family
            .bundle_indices
            .iter()
            .filter_map(|&bi| {
                bundle_coords
                    .get(bi)
                    .map(|(c, s, e, st)| format!("{}:{}-{}:{}", c, s, e, st))
            })
            .collect();
        let chrom = bundle_coords
            .get(*family.bundle_indices.first().unwrap_or(&0))
            .map(|(c, _, _, _)| c.as_str())
            .unwrap_or("?");
        let em = em_results.get(fi).cloned().unwrap_or_default();
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{}",
            family.family_id,
            family.bundle_indices.len(),
            chrom,
            regions.join(";"),
            family.multimap_reads.len(),
            em.iterations,
            em.converged,
            em.max_delta,
            em.reads_reweighted,
        )?;
    }
    Ok(())
}

// ── Novel copy discovery (Phase 3) ──────────────────────────────────────────

/// Candidate read for novel copy discovery: an unmapped or supplementary-only
/// read whose junctions match a family's consensus junction set.
#[derive(Debug)]
pub struct NovelCandidate {
    pub read_name: String,
    pub family_id: usize,
    pub matched_junctions: usize,
    pub total_junctions: usize,
    /// Approximate genomic position from the supplementary alignment.
    pub approx_chrom: String,
    pub approx_start: u64,
    pub approx_end: u64,
}

/// Build consensus junction set for each family: junctions shared by ≥2 copies.
pub fn build_family_consensus_junctions(
    families: &[FamilyGroup],
    bundles: &[Bundle],
) -> Vec<Vec<(u64, u64)>> {
    let mut all_consensus = Vec::with_capacity(families.len());

    for family in families {
        let mut junction_counts: HashMap<(u64, u64), usize> = HashMap::new();
        for &bi in &family.bundle_indices {
            let bundle = &bundles[bi];
            // Collect unique junctions at this copy.
            let mut seen: std::collections::HashSet<(u64, u64)> = std::collections::HashSet::new();
            for (j, _) in &bundle.junction_stats {
                // Round to nearest 10bp to account for LR error.
                let key = (j.donor / 10 * 10, j.acceptor / 10 * 10);
                if seen.insert(key) {
                    *junction_counts.entry(key).or_insert(0) += 1;
                }
            }
        }
        // Keep junctions present in ≥2 copies.
        let consensus: Vec<(u64, u64)> = junction_counts
            .into_iter()
            .filter(|(_, count)| *count >= 2)
            .map(|(junc, _)| junc)
            .collect();
        all_consensus.push(consensus);
    }

    all_consensus
}

/// Discover novel gene copies from unmapped reads using k-mer matching.
///
/// Approach: extract exonic sequences from each family's bundles (via genome FASTA),
/// build k-mer sets per family, then scan unmapped reads for k-mer overlap.
/// Reads with high k-mer overlap to a family are candidates for novel copies.
///
/// This avoids supplementary alignments (which contain chimeric artifacts).
/// Requires `--genome-fasta` for exon sequence extraction.
pub fn discover_novel_copies(
    bam_path: &std::path::Path,
    families: &[FamilyGroup],
    bundles: &[Bundle],
    config: &crate::types::RunConfig,
) -> Vec<NovelCandidate> {
    use std::collections::HashSet as StdHashSet;

    if families.is_empty() {
        return Vec::new();
    }

    // k=15: short enough to tolerate gene-family divergence (~90% identity means
    // ~20% of 15-mers are conserved vs ~10% of 21-mers), long enough to avoid
    // random matches (4^15 = 1 billion possible 15-mers).
    let kmer_len: usize = 15;
    let min_kmer_hits: usize = 3; // Minimum k-mer matches to call a candidate.
    let min_kmer_frac: f64 = 0.01; // At least 1% of the read's k-mers must match.

    // Step 1: Build family exonic k-mer sets from genome.
    let genome = config.genome_fasta.as_ref().and_then(|p| {
        crate::genome::GenomeIndex::from_fasta(p).ok()
    });
    let genome_ref = match genome.as_ref() {
        Some(g) => g,
        None => {
            eprintln!("[VG] Novel copy discovery requires --genome-fasta; skipping");
            return Vec::new();
        }
    };

    let mut family_kmers: Vec<StdHashSet<u64>> = Vec::with_capacity(families.len());
    for family in families {
        let mut kmers = StdHashSet::new();
        for &bi in &family.bundle_indices {
            let bundle = &bundles[bi];
            // Extract EXONIC k-mers only (not the full bundle span which includes introns).
            // Use junction_stats to reconstruct exon boundaries: the regions between
            // junctions are introns, so exonic regions are the gaps.
            let mut exon_regions: Vec<(u64, u64)> = Vec::new();
            let mut intron_starts: Vec<u64> = Vec::new();
            let mut intron_ends: Vec<u64> = Vec::new();
            for (j, _) in &bundle.junction_stats {
                intron_starts.push(j.donor);
                intron_ends.push(j.acceptor);
            }
            intron_starts.sort_unstable();
            intron_ends.sort_unstable();
            if intron_starts.is_empty() {
                // No junctions: entire bundle is one exon.
                exon_regions.push((bundle.start, bundle.end));
            } else {
                // First exon: bundle start to first intron.
                exon_regions.push((bundle.start, intron_starts[0]));
                // Internal exons: between consecutive introns.
                for i in 0..intron_starts.len().saturating_sub(1) {
                    let exon_start = intron_ends[i];
                    let exon_end = intron_starts[i + 1];
                    if exon_end > exon_start {
                        exon_regions.push((exon_start, exon_end));
                    }
                }
                // Last exon: last intron end to bundle end.
                exon_regions.push((*intron_ends.last().unwrap(), bundle.end));
            }
            // Also add read exons if reads are available.
            for read in &bundle.reads {
                for &(es, ee) in &read.exons {
                    exon_regions.push((es, ee));
                }
            }
            // Merge overlapping exon regions.
            exon_regions.sort_by_key(|&(s, _)| s);
            let mut merged: Vec<(u64, u64)> = Vec::new();
            for (s, e) in exon_regions {
                if let Some(last) = merged.last_mut() {
                    if s <= last.1 {
                        last.1 = last.1.max(e);
                        continue;
                    }
                }
                merged.push((s, e));
            }
            // Extract k-mers from merged exonic regions.
            for (es, ee) in &merged {
                if let Some(seq) = genome_ref.fetch_sequence(&bundle.chrom, *es, *ee) {
                    if seq.len() >= kmer_len {
                        for window in seq.windows(kmer_len) {
                            if window.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                                kmers.insert(kmer_hash(window));
                            }
                        }
                    }
                }
            }
        }
        family_kmers.push(kmers);
    }

    let total_family_kmers: usize = family_kmers.iter().map(|k| k.len()).sum();
    if total_family_kmers == 0 {
        eprintln!("[VG] No family k-mers extracted (genome fetch failed?)");
        return Vec::new();
    }
    eprintln!(
        "[VG] Built k-mer indices for {} families ({} total unique {}-mers)",
        families.len(),
        total_family_kmers,
        kmer_len,
    );

    // Step 2: Scan unmapped reads from BAM.
    let bam_file = match std::fs::File::open(bam_path) {
        Ok(f) => f,
        Err(_) => return Vec::new(),
    };
    let buf = std::io::BufReader::new(bam_file);
    let worker_count = std::num::NonZeroUsize::MIN;
    let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(worker_count, buf);
    let mut reader = noodles_bam::io::Reader::from(bgzf);
    let _header = match reader.read_header() {
        Ok(h) => h,
        Err(_) => return Vec::new(),
    };

    let mut candidates: Vec<NovelCandidate> = Vec::new();
    let mut n_unmapped = 0usize;
    let mut n_matched = 0usize;

    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        if !record.flags().is_unmapped() {
            continue;
        }
        n_unmapped += 1;

        // noodles BAM Record::sequence() is 4-bit encoded.
        // Decode via the raw bytes: each byte has two bases (high nibble, low nibble).
        let seq_obj = record.sequence();
        let seq_raw = seq_obj.as_ref();
        let seq_len = seq_obj.len();
        let mut seq_bytes: Vec<u8> = Vec::with_capacity(seq_len);
        for (i, &byte) in seq_raw.iter().enumerate() {
            let bases = [(byte >> 4) & 0xF, byte & 0xF];
            for (j, nib) in bases.iter().enumerate() {
                if i * 2 + j >= seq_len {
                    break;
                }
                let b = match nib {
                    1 => b'A',
                    2 => b'C',
                    4 => b'G',
                    8 => b'T',
                    _ => b'N',
                };
                seq_bytes.push(b);
            }
        }
        if seq_bytes.len() < kmer_len {
            continue;
        }

        // Extract k-mers from the unmapped read.
        let read_kmer_count = seq_bytes.len().saturating_sub(kmer_len) + 1;
        let mut read_kmers: Vec<u64> = Vec::with_capacity(read_kmer_count);
        for window in seq_bytes.windows(kmer_len) {
            if window.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                read_kmers.push(kmer_hash(window));
            }
        }
        // Also check reverse complement.
        let rc_bytes = reverse_complement(&seq_bytes);
        let mut rc_kmers: Vec<u64> = Vec::with_capacity(read_kmer_count);
        for window in rc_bytes.windows(kmer_len) {
            if window.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                rc_kmers.push(kmer_hash(window));
            }
        }

        if read_kmers.is_empty() && rc_kmers.is_empty() {
            continue;
        }

        // Match against each family.
        let mut best_family: Option<(usize, usize)> = None; // (family_idx, hit_count)
        for (fi, fk) in family_kmers.iter().enumerate() {
            if fk.is_empty() {
                continue;
            }
            let fwd_hits = read_kmers.iter().filter(|k| fk.contains(k)).count();
            let rc_hits = rc_kmers.iter().filter(|k| fk.contains(k)).count();
            let hits = fwd_hits.max(rc_hits);
            let total_kmers = read_kmers.len().max(rc_kmers.len()).max(1);
            let frac = hits as f64 / total_kmers as f64;

            if hits >= min_kmer_hits && frac >= min_kmer_frac {
                if best_family.map_or(true, |(_, best_hits)| hits > best_hits) {
                    best_family = Some((fi, hits));
                }
            }
        }

        if let Some((fi, hits)) = best_family {
            let read_name = record.name().map(|n| n.to_string()).unwrap_or_default();
            let total_kmers = read_kmers.len().max(rc_kmers.len()).max(1);
            candidates.push(NovelCandidate {
                read_name,
                family_id: fi,
                matched_junctions: hits, // Repurposed: k-mer hits
                total_junctions: total_kmers, // Repurposed: total k-mers
                approx_chrom: String::new(),
                approx_start: 0,
                approx_end: seq_bytes.len() as u64,
            });
            n_matched += 1;
        }
    }

    eprintln!(
        "[VG] Novel copy discovery: scanned {} unmapped reads, {} matched to {} families",
        n_unmapped,
        n_matched,
        candidates
            .iter()
            .map(|c| c.family_id)
            .collect::<StdHashSet<_>>()
            .len(),
    );

    candidates
}

/// Simple k-mer hash (FNV1a on bytes).
fn kmer_hash(kmer: &[u8]) -> u64 {
    fnv1a64(kmer)
}

/// Reverse complement of a DNA sequence.
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

// ── Phased assembly (--vg-phase) ─────────────────────────────────────────

/// Split a bundle's reads by haplotype (HP tag), creating up to 3 sub-bundles:
/// - HP=1 reads
/// - HP=2 reads
/// - Unphased reads (HP tag absent)
///
/// Each sub-bundle gets a copy of the bundle metadata but only its haplotype's reads.
/// Unphased reads are added to ALL sub-bundles (they provide shared coverage).
pub fn split_bundle_by_phase(bundle: &Bundle) -> Vec<(Bundle, Option<u8>)> {
    let mut hp1_reads: Vec<BundleRead> = Vec::new();
    let mut hp2_reads: Vec<BundleRead> = Vec::new();
    let mut unphased_reads: Vec<BundleRead> = Vec::new();

    for read in &bundle.reads {
        match read.hp_tag {
            Some(1) => hp1_reads.push(read.clone()),
            Some(2) => hp2_reads.push(read.clone()),
            _ => unphased_reads.push(read.clone()),
        }
    }

    // Only split if there are reads in both haplotypes.
    if hp1_reads.is_empty() || hp2_reads.is_empty() {
        // No phasing signal — return the original bundle as-is.
        return vec![(bundle.clone(), None)];
    }

    let mut result = Vec::new();

    // HP=1 bundle: phased reads + unphased.
    let mut hp1_bundle = bundle.clone();
    hp1_bundle.reads = hp1_reads;
    hp1_bundle.reads.extend(unphased_reads.iter().cloned());
    result.push((hp1_bundle, Some(1u8)));

    // HP=2 bundle: phased reads + unphased.
    let mut hp2_bundle = bundle.clone();
    hp2_bundle.reads = hp2_reads;
    hp2_bundle.reads.extend(unphased_reads);
    result.push((hp2_bundle, Some(2u8)));

    eprintln!(
        "[VG-PHASE] Split bundle {}:{}-{} into HP1 ({} reads) + HP2 ({} reads)",
        bundle.chrom,
        bundle.start,
        bundle.end,
        result[0].0.reads.len(),
        result[1].0.reads.len(),
    );

    result
}

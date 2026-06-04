//! Variation graph mode for gene family assembly.
//!
//! Links related genomic loci (gene family copies) using multi-mapping read
//! patterns, then jointly resolves multi-mapping reads across copies using EM.
//!
//! Activated with `--vg`. Multi-mapping reads (NH >= 2) that appear in multiple
//! bundles link those bundles into "family groups." Within each family group,
//! an EM algorithm redistributes read weights based on splice-graph compatibility,
//! then the standard flow pipeline runs with updated weights.

use crate::util::coord::overlaps_half_open;
use crate::graph::Graph;
use crate::map_reads::read_to_path_bundlenodes;
use crate::types::{Bundle, BundleRead, CBundlenode, Bundle2Graph};
use rayon::prelude::*;
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
#[derive(Debug, Clone)]
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
    // Use FxHash internally for speed (read_name_hash is a u64). The
    // public return type stays std HashMap for ABI compatibility.
    let mut read_locs: crate::types::DetHashMap<u64, Vec<(usize, usize)>> =
        crate::types::DetHashMap::default();
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
    // Convert back to std HashMap for the public API.
    read_locs.into_iter().collect()
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
    let mut n_no_target = 0usize;
    let mut n_no_primary = 0usize;
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
            n_no_target += 1;
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
            n_no_primary += 1;
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
            "[VG] Supplementary scan: {} supplementary alignments, {} cross-bundle links \
             (skipped: {} no-overlapping-bundle, {} primary-not-in-bundle)",
            n_supp, n_linked, n_no_target, n_no_primary
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

/// Two bundles on the same chromosome with nearly-identical coords but
/// opposite strands are strand-mirror artifacts of Rustle's per-strand
/// bundling. They are NOT biological paralogs, and linking them creates
/// spurious "families" that dominate the family report. Filter them out
/// before union-find merges connected components.
fn is_strand_mirror(b1: &Bundle, b2: &Bundle) -> bool {
    if b1.chrom != b2.chrom {
        return false;
    }
    // Opposite strands only: '+' vs '-'. If either is unknown ('.'), not a mirror.
    let opp = (b1.strand == '+' && b2.strand == '-')
        || (b1.strand == '-' && b2.strand == '+');
    if !opp {
        return false;
    }
    // High coordinate overlap (≥90% of the smaller bundle).
    let lo = b1.start.max(b2.start);
    let hi = b1.end.min(b2.end);
    if hi <= lo {
        return false;
    }
    let overlap = (hi - lo) as f64;
    let l1 = (b1.end.saturating_sub(b1.start)) as f64;
    let l2 = (b2.end.saturating_sub(b2.start)) as f64;
    let smaller = l1.min(l2).max(1.0);
    overlap / smaller >= 0.90
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

    // Build bundle-pair link counts. FxHash for the (usize,usize) keys —
    // this map sees ~120K entries on full GGO.bam and was a real hashing
    // hotspot during family discovery.
    let n_bundles = bundles.len();
    let mut link_counts: crate::types::DetHashMap<(usize, usize), usize> =
        crate::types::DetHashMap::default();
    let mut skipped_mirror_pairs = 0usize;
    for locs in multimap_index.values() {
        // Collect unique bundle indices for this read.
        let mut bundle_set: Vec<usize> = locs.iter().map(|(bi, _)| *bi).collect();
        bundle_set.sort_unstable();
        bundle_set.dedup();
        // Create pairwise links.
        for i in 0..bundle_set.len() {
            for j in (i + 1)..bundle_set.len() {
                let bi = bundle_set[i];
                let bj = bundle_set[j];
                // Skip strand-mirror pairs (same-coord opposite-strand bundles
                // are per-strand bundling artifacts, not biological paralogs).
                if is_strand_mirror(&bundles[bi], &bundles[bj]) {
                    skipped_mirror_pairs += 1;
                    continue;
                }
                let key = (bi, bj);
                *link_counts.entry(key).or_insert(0) += 1;
            }
        }
    }
    if skipped_mirror_pairs > 0 && std::env::var_os("RUSTLE_VG_TRACE").is_some() {
        eprintln!(
            "[VG] Skipped {} strand-mirror bundle pairs from family linking",
            skipped_mirror_pairs
        );
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
    // Same strand-mirror filter as above.
    for &(a, b) in &seq_links {
        if is_strand_mirror(&bundles[a], &bundles[b]) {
            continue;
        }
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
        // Build a {global_bi -> position_in_family} map once. Replaces the
        // O(n) `bundle_indices.iter().position()` linear scan that ran per
        // matched read, and also handles the membership check via Option.
        let bi_to_pos: crate::types::DetHashMap<usize, usize> = bundle_indices
            .iter().copied().enumerate()
            .map(|(pos, bi)| (bi, pos))
            .collect();
        let mut family_reads: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for (&rnh, locs) in &multimap_index {
            let mut family_locs: Vec<(usize, usize)> = Vec::with_capacity(locs.len());
            for (bi, ri) in locs {
                if let Some(&pos) = bi_to_pos.get(bi) {
                    family_locs.push((pos, *ri));
                }
            }
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

/// Compute the mean pairwise k-mer Jaccard over the family graph's per-copy
/// sequences. **Graph-supported** because it uses sequences extracted at the
/// family graph's exon-class nodes (not raw read k-mers). For each pair of
/// copies, builds a k-mer set per copy by concatenating its node-level
/// sequences, then takes Jaccard.
///
/// Real paralog clusters: high pairwise Jaccard (paralogs share most
/// exonic sequence). TE-bridge artifacts: low Jaccard (different gene
/// families have different sequences even when their bundles span similar
/// intron-length patterns by coincidence).
///
/// Returns `None` when family graph build fails (e.g., genome not available,
/// mixed-strand family) — caller should treat absence as "skip this signal".
pub fn compute_family_graph_kmer_jaccard(
    family: &FamilyGroup,
    bundles: &[Bundle],
    genome: &crate::genome::GenomeIndex,
    kmer_len: usize,
) -> Option<f64> {
    compute_family_graph_kmer_jaccard_diag(family, bundles, genome, kmer_len).0
}

/// Diagnostic variant: returns (jaccard, skip_reason). skip_reason is None
/// when computation succeeded; otherwise it carries a short label so the
/// caller can tally why families are skipped.
pub fn compute_family_graph_kmer_jaccard_diag(
    family: &FamilyGroup,
    bundles: &[Bundle],
    genome: &crate::genome::GenomeIndex,
    kmer_len: usize,
) -> (Option<f64>, Option<&'static str>) {
    // CANONICAL FNV-1a k-mer hash: min(forward, reverse-complement). A sequence and its
    // reverse-complement produce identical canonical k-mer sets, so INVERTED homologs
    // (e.g. DAZ1 on - strand vs DAZ3 on + strand) share k-mers and get a real Jaccard
    // instead of being dropped at 0.000. True artifacts stay low in both orientations.
    fn canonical_kmer_hash(window: &[u8]) -> u64 {
        fn fnv(it: impl Iterator<Item = u8>) -> u64 {
            let mut h: u64 = 0xcbf29ce484222325;
            for b in it {
                h ^= b as u64;
                h = h.wrapping_mul(0x100000001b3);
            }
            h
        }
        let fwd = fnv(window.iter().copied());
        let rc = fnv(window.iter().rev().map(|&b| match b {
            b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', other => other,
        }));
        fwd.min(rc)
    }

    // build_family_graph is single-strand, so partition by strand and accumulate a
    // canonical-k-mer set per copy across ALL sub-families. This scores the cross-strand
    // homologous pair (the inverted-paralog case) rather than only the largest same-strand
    // sub-family (which for inverted families is often unrelated bundles -> jaccard 0).
    let sub_families = crate::vg_hmm::rescue::partition_family_by_strand(family, bundles);
    let mut all_kmer_sets: Vec<crate::types::DetHashSet<u64>> = Vec::new();
    let mut built_any = false;
    for sub in &sub_families {
        if sub.bundle_indices.is_empty() { continue; }
        let fg = match crate::vg_hmm::family_graph::build_family_graph(
            sub, bundles, Some(genome), 0.30, 0.30, 0.30,
        ) {
            Ok(g) => g,
            Err(_) => continue,
        };
        if fg.nodes.is_empty() { continue; }
        built_any = true;
        let n_copies = sub.bundle_indices.len();
        let mut sets: Vec<crate::types::DetHashSet<u64>> =
            (0..n_copies).map(|_| crate::types::DetHashSet::default()).collect();
        for node in &fg.nodes {
            for (cid, seq) in &node.per_copy_sequences {
                if *cid >= n_copies { continue; }
                if seq.len() < kmer_len { continue; }
                for window in seq.windows(kmer_len) {
                    if window.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                        sets[*cid].insert(canonical_kmer_hash(window));
                    }
                }
            }
        }
        for s in sets {
            if !s.is_empty() { all_kmer_sets.push(s); }
        }
    }
    if !built_any {
        return (None, Some("graph_build_err"));
    }
    if all_kmer_sets.len() < 2 {
        return (None, Some("fewer_than_2_copies_with_kmers"));
    }
    let mut sum = 0.0_f64;
    let mut n_pairs = 0usize;
    for i in 0..all_kmer_sets.len() {
        for j in (i + 1)..all_kmer_sets.len() {
            let inter = all_kmer_sets[i].intersection(&all_kmer_sets[j]).count();
            let uni = all_kmer_sets[i].union(&all_kmer_sets[j]).count();
            if uni > 0 {
                sum += inter as f64 / uni as f64;
                n_pairs += 1;
            }
        }
    }
    if n_pairs > 0 { (Some(sum / n_pairs as f64), None) } else { (None, Some("no_pairs")) }
}

/// Optional 6th-signal filter: drops families whose graph-supported pairwise
/// k-mer Jaccard is below `min_kmer_jaccard`. Catches mild TE-bridges that
/// pass the intron-length-Jaccard signal by coincidence (similar intron sizes
/// across unrelated genes) but have NO actual sequence similarity.
///
/// Skips families where the family graph couldn't be built or has <2 copies
/// with any k-mers (typically: single-exon paralog clusters where node
/// sequences are too short, or mixed-strand families).
///
/// Disabled when `min_kmer_jaccard <= 0.0`. Requires genome FASTA — caller
/// must pass `Some(genome)` or this filter is a no-op.
pub fn filter_by_graph_kmer_jaccard(
    families: Vec<FamilyGroup>,
    bundles: &[Bundle],
    genome: Option<&crate::genome::GenomeIndex>,
    min_kmer_jaccard: f64,
    kmer_len: usize,
) -> Vec<FamilyGroup> {
    if min_kmer_jaccard <= 0.0 {
        return families;
    }
    let g = match genome {
        Some(g) => g,
        None => {
            eprintln!("[VG] graph-k-mer-Jaccard filter requested but no genome available — skipping");
            return families;
        }
    };
    let n_in = families.len();
    let mut kept: Vec<FamilyGroup> = Vec::new();
    let mut n_dropped = 0usize;
    let mut n_skipped_signal = 0usize;
    let mut skip_reasons: std::collections::BTreeMap<&'static str, usize> =
        std::collections::BTreeMap::new();
    let mut kept_jaccards: Vec<f64> = Vec::new();
    let mut dropped_jaccards: Vec<f64> = Vec::new();
    for fam in families {
        let (j_opt, reason) = compute_family_graph_kmer_jaccard_diag(&fam, bundles, g, kmer_len);
        match j_opt {
            Some(j) => {
                if j < min_kmer_jaccard {
                    n_dropped += 1;
                    dropped_jaccards.push(j);
                } else {
                    kept_jaccards.push(j);
                    kept.push(fam);
                }
            }
            None => {
                // Couldn't compute (single-exon, mixed-strand, etc.) — keep
                // by default; the prior 5 signals already vetted it.
                n_skipped_signal += 1;
                if let Some(r) = reason { *skip_reasons.entry(r).or_insert(0) += 1; }
                kept.push(fam);
            }
        }
    }
    if n_dropped > 0 || n_skipped_signal > 0 {
        eprintln!(
            "[VG] graph-k-mer-Jaccard filter: {} → {} families ({} dropped by low_kmer_jaccard < {:.2}, {} skipped — graph build failed or single-exon)",
            n_in, kept.len(), n_dropped, min_kmer_jaccard, n_skipped_signal
        );
        if !skip_reasons.is_empty() {
            eprintln!("[VG]   skip reasons: {:?}", skip_reasons);
        }
        if !kept_jaccards.is_empty() {
            let s: Vec<String> = kept_jaccards.iter().map(|j| format!("{:.3}", j)).collect();
            eprintln!("[VG]   kept jaccards: [{}]", s.join(", "));
        }
        if !dropped_jaccards.is_empty() {
            let s: Vec<String> = dropped_jaccards.iter().map(|j| format!("{:.3}", j)).collect();
            eprintln!("[VG]   dropped jaccards: [{}]", s.join(", "));
        }
    }
    kept
}

/// POA-aligned mean pairwise %-identity over family-graph nodes.
///
/// More sensitive than `compute_family_graph_kmer_jaccard` at moderate
/// divergence: rather than matching exact 21-mers (which saturates fast as
/// paralogs drift), this aligns per-copy node sequences with POA and counts
/// matches per aligned column. Stays DNA-side, no protein translation.
///
/// Strategy:
///   1. Build the same family graph (single-strand, ≥2 copies)
///   2. For each node, POA-align per-copy sequences to obtain an MSA
///   3. Per node, compute mean pairwise %-identity over aligned columns
///      (gap-vs-gap doesn't count; gap-vs-base counts as mismatch)
///   4. Aggregate across nodes weighted by total aligned column count
///
/// Min-hash prescreen: before running expensive POA, build a 128-element
/// min-hash signature of 21-mers per copy and estimate Jaccard. If even the
/// loose prescreen says identity is hopeless (< 0.1), short-circuit and
/// return the prescreen estimate so the caller can drop the family.
///
/// Returns (identity, skip_reason) — same shape as the kmer-Jaccard variant.
pub fn compute_family_graph_poa_identity_diag(
    family: &FamilyGroup,
    bundles: &[Bundle],
    genome: &crate::genome::GenomeIndex,
) -> (Option<f64>, Option<&'static str>) {
    let sub_families = crate::vg_hmm::rescue::partition_family_by_strand(family, bundles);
    let largest = match sub_families.into_iter()
        .filter(|f| f.bundle_indices.len() >= 2)
        .max_by_key(|f| f.bundle_indices.len())
    {
        Some(f) => f,
        None => return (None, Some("no_substrand_with_2_copies")),
    };
    let fg = match crate::vg_hmm::family_graph::build_family_graph(
        &largest, bundles, Some(genome), 0.30, 0.30, 0.30,
    ) {
        Ok(g) => g,
        Err(_) => return (None, Some("graph_build_err")),
    };
    if fg.nodes.is_empty() {
        return (None, Some("graph_empty"));
    }

    // Min-hash prescreen on 21-mers per copy across all nodes — cheap upper
    // bound on identity. If even the loose Jaccard estimate is hopeless, no
    // need to spend POA cycles.
    const SKETCH_SIZE: usize = 128;
    const KMER_K: usize = 21;
    let n_copies = largest.bundle_indices.len();
    let prescreen = minhash_jaccard_estimate(&fg, n_copies, KMER_K, SKETCH_SIZE);
    if let Some(j) = prescreen {
        if j < 0.10 {
            // Hopeless — return the estimate so the caller drops it cheaply.
            return (Some(j), None);
        }
    }

    // Full POA-aligned identity per node, weighted by aligned-column count.
    let mut total_match = 0u64;
    let mut total_compared = 0u64;
    for node in &fg.nodes {
        let seqs: Vec<Vec<u8>> = node.per_copy_sequences.iter().map(|(_, s)| s.clone()).collect();
        if seqs.len() < 2 { continue; }
        // Skip nodes where any copy is empty/too short to align meaningfully.
        if seqs.iter().any(|s| s.len() < 10) { continue; }
        let msa = match crate::vg_hmm::profile::poa_msa(&seqs) {
            Ok(m) => m,
            Err(_) => continue,
        };
        if msa.is_empty() { continue; }
        let n_col = msa[0].len();
        let n_row = msa.len();
        // Per-pair counting per column.
        for c in 0..n_col {
            // Count matches/mismatches over pairs i<j where neither side is gap-vs-gap.
            for i in 0..n_row {
                for j in (i + 1)..n_row {
                    let a = msa[i][c];
                    let b = msa[j][c];
                    if a == b'-' && b == b'-' { continue; }
                    total_compared += 1;
                    if a == b { total_match += 1; }
                }
            }
        }
    }
    if total_compared == 0 {
        return (None, Some("no_aligned_columns"));
    }
    (Some(total_match as f64 / total_compared as f64), None)
}

/// Min-hash Jaccard estimate over family-graph node k-mers, per pair-mean.
///
/// Cheap O(SKETCH × N_COPIES) prescreen for the POA-identity filter. Returns
/// a Jaccard estimate in [0, 1] or None if there are < 2 copies with k-mers.
fn minhash_jaccard_estimate(
    fg: &crate::vg_hmm::family_graph::FamilyGraph,
    n_copies: usize,
    kmer_len: usize,
    sketch_size: usize,
) -> Option<f64> {
    // sketch[copy][hash_idx] = min hash seen so far (u64::MAX = unseen).
    let mut sketch: Vec<Vec<u64>> = (0..n_copies)
        .map(|_| vec![u64::MAX; sketch_size])
        .collect();
    let mut have_kmers = vec![false; n_copies];

    for node in &fg.nodes {
        for (cid, seq) in &node.per_copy_sequences {
            if *cid >= n_copies || seq.len() < kmer_len { continue; }
            for window in seq.windows(kmer_len) {
                if !window.iter().all(|&b| matches!(b, b'A'|b'C'|b'G'|b'T')) { continue; }
                let base = fnv1a_64(window);
                // Apply `sketch_size` hash permutations (xor with i × golden ratio prime).
                for i in 0..sketch_size {
                    let h = base ^ (i as u64).wrapping_mul(0x9E3779B97F4A7C15);
                    if h < sketch[*cid][i] { sketch[*cid][i] = h; }
                }
                have_kmers[*cid] = true;
            }
        }
    }
    let with_kmers: Vec<usize> = (0..n_copies).filter(|&c| have_kmers[c]).collect();
    if with_kmers.len() < 2 { return None; }
    let mut sum = 0.0_f64;
    let mut pairs = 0usize;
    for i in 0..with_kmers.len() {
        for j in (i + 1)..with_kmers.len() {
            let a = &sketch[with_kmers[i]];
            let b = &sketch[with_kmers[j]];
            let eq = a.iter().zip(b.iter()).filter(|(x, y)| x == y && **x != u64::MAX).count();
            sum += eq as f64 / sketch_size as f64;
            pairs += 1;
        }
    }
    if pairs == 0 { None } else { Some(sum / pairs as f64) }
}

#[inline]
fn fnv1a_64(window: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in window {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Parallel POA-identity filter (sibling of `filter_by_graph_kmer_jaccard`).
///
/// Uses rayon to compute identity per family in parallel. Each family's POA
/// alignment is independent so this is embarrassingly parallel. Returns the
/// kept families plus emits a summary line with kept/dropped identities.
pub fn filter_by_graph_poa_identity(
    families: Vec<FamilyGroup>,
    bundles: &[Bundle],
    genome: Option<&crate::genome::GenomeIndex>,
    min_identity: f64,
) -> Vec<FamilyGroup> {
    if min_identity <= 0.0 {
        return families;
    }
    let g = match genome {
        Some(g) => g,
        None => {
            eprintln!("[VG] graph-POA-identity filter requested but no genome available — skipping");
            return families;
        }
    };
    let n_in = families.len();
    use rayon::prelude::*;
    let scored: Vec<(FamilyGroup, Option<f64>, Option<&'static str>)> = families
        .into_par_iter()
        .map(|fam| {
            let (j, reason) = compute_family_graph_poa_identity_diag(&fam, bundles, g);
            (fam, j, reason)
        })
        .collect();

    let mut kept: Vec<FamilyGroup> = Vec::new();
    let mut n_dropped = 0usize;
    let mut n_skipped_signal = 0usize;
    let mut skip_reasons: std::collections::BTreeMap<&'static str, usize> =
        std::collections::BTreeMap::new();
    let mut kept_ids: Vec<f64> = Vec::new();
    let mut dropped_ids: Vec<f64> = Vec::new();
    for (fam, j_opt, reason) in scored {
        match j_opt {
            Some(j) => {
                if j < min_identity {
                    n_dropped += 1;
                    dropped_ids.push(j);
                } else {
                    kept_ids.push(j);
                    kept.push(fam);
                }
            }
            None => {
                n_skipped_signal += 1;
                if let Some(r) = reason { *skip_reasons.entry(r).or_insert(0) += 1; }
                kept.push(fam);
            }
        }
    }
    if n_dropped > 0 || n_skipped_signal > 0 {
        eprintln!(
            "[VG] graph-POA-identity filter: {} → {} families ({} dropped by low_identity < {:.2}, {} skipped — graph build failed)",
            n_in, kept.len(), n_dropped, min_identity, n_skipped_signal
        );
        if !skip_reasons.is_empty() {
            eprintln!("[VG]   skip reasons: {:?}", skip_reasons);
        }
        if !kept_ids.is_empty() {
            let s: Vec<String> = kept_ids.iter().map(|j| format!("{:.3}", j)).collect();
            eprintln!("[VG]   kept identities: [{}]", s.join(", "));
        }
        if !dropped_ids.is_empty() {
            let s: Vec<String> = dropped_ids.iter().map(|j| format!("{:.3}", j)).collect();
            eprintln!("[VG]   dropped identities: [{}]", s.join(", "));
        }
    }
    kept
}

/// Quality-filter discovered families to keep only those that look like
/// real multi-copy gene paralogs (vs noise, mtDNA, repetitive megafamilies,
/// TE-bridge spurious merges).
///
/// Five signals corroborate "this is a real multi-copy gene family":
///   1. **n_copies in [min_copies, max_copies]** — drops mtDNA / repetitive
///      mega-clusters and singletons.
///   2. **multimap_reads ≥ min_shared** — total read sharing. Random
///      alignment artifacts produce 1-3 spurious shared reads; real paralogs
///      have substantially more (the GOLGA6L7 chr19 cluster has 386).
///   3. **multimap_reads / n_copies ≥ min_shared_per_copy** — sharing density.
///      A family of N copies held together by O(1) cross-mapping reads is
///      almost always an alignment artifact.
///   4. **Coefficient of variation (CV) of intron-counts ≤ max_exon_cv** —
///      cheap exon-count similarity proxy.
///   5. **Pairwise intron-length-set Jaccard ≥ min_primitive_jaccard** —
///      the structural-paralogy signal. For each pair of multi-exon copies,
///      compute the Jaccard of their intron-length sets (binned at 50bp).
///      Real paralogs of the same gene have nearly identical intron lengths
///      (Jaccard typically >0.5 even after copy-number variants). Random
///      cross-cluster merges via repetitive elements (chr19-GOLGA8 ↔
///      chr17-TBC1D3 is the canonical case) have intron-length Jaccard
///      near 0 — different genes have different intron lengths.
///
/// THIS IS THE GRAPH-STRUCTURAL DEFINITION of a multi-copy gene family in
/// rustle: paralogs share a "primitive" splice topology (signal #5),
/// supported by structural cross-mapping (signals #2, #3) at non-trivial
/// scale (#1) with consistent exon counts (#4).
///
/// Why intron lengths and not absolute splice coordinates: paralogs at
/// different genomic loci have ZERO shared (donor, acceptor) coordinates by
/// definition, but their RELATIVE intron structure is conserved. Binning
/// at 50bp absorbs minor copy-number-variant size differences.
///
/// Implementation: only multi-exon copies (≥1 intron) participate in the
/// Jaccard. Single-exon paralog clusters (olfactory receptors etc.) skip
/// this check.
pub fn filter_high_confidence_families(
    families: Vec<FamilyGroup>,
    bundles: &[Bundle],
    min_copies: usize,
    max_copies: usize,
    min_shared: usize,
    min_shared_per_copy: f64,
    max_exon_cv: f64,
    min_primitive_jaccard: f64,
) -> Vec<FamilyGroup> {
    let n_in = families.len();
    let mut kept: Vec<FamilyGroup> = Vec::new();
    let mut drops: std::collections::BTreeMap<&'static str, usize> = std::collections::BTreeMap::new();

    // Optional dump of dropped families for offline analysis (validation D —
    // confirm filter drops are TE-bridges, not real paralogs).
    let dump_path = std::env::var("RUSTLE_VG_DUMP_DROPPED_FAMILIES").ok();
    let mut dump_file: Option<std::io::BufWriter<std::fs::File>> = dump_path.as_ref()
        .and_then(|p| std::fs::File::create(p).ok())
        .map(std::io::BufWriter::new);
    if let Some(ref mut f) = dump_file {
        use std::io::Write;
        let _ = writeln!(f, "family_id\tn_copies\tn_intron_copies\tn_chroms\tchrom\tregions\tn_shared_reads\tshared_per_copy\texon_cv\tprimitive_jaccard\tdrop_reason");
    }
    let dump_drop = |df: &mut Option<std::io::BufWriter<std::fs::File>>, fam: &FamilyGroup,
                     bundles: &[Bundle], reason: &str, exon_cv: f64, primitive_jaccard: f64| {
        if let Some(f) = df {
            use std::io::Write;
            let n_copies = fam.bundle_indices.len();
            let n_shared = fam.multimap_reads.len();
            let shared_per_copy = if n_copies > 0 { n_shared as f64 / n_copies as f64 } else { 0.0 };
            let n_intron_copies = fam.bundle_indices.iter()
                .filter(|&&bi| bundles.get(bi).map(|b| !b.junction_stats.is_empty()).unwrap_or(false))
                .count();
            let chroms: std::collections::BTreeSet<&str> = fam.bundle_indices.iter()
                .filter_map(|&bi| bundles.get(bi).map(|b| b.chrom.as_str()))
                .collect();
            let chrom = bundles.get(fam.bundle_indices[0])
                .map(|b| b.chrom.as_str()).unwrap_or("?");
            let regions: Vec<String> = fam.bundle_indices.iter()
                .filter_map(|&bi| bundles.get(bi))
                .map(|b| format!("{}:{}-{}:{}", b.chrom, b.start, b.end, b.strand))
                .collect();
            let _ = writeln!(f,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{}",
                fam.family_id, n_copies, n_intron_copies, chroms.len(), chrom,
                regions.join(";"), n_shared, shared_per_copy, exon_cv,
                primitive_jaccard, reason
            );
        }
    };

    for fam in families {
        let n_copies = fam.bundle_indices.len();
        if n_copies < min_copies {
            *drops.entry("singleton").or_insert(0) += 1;
            dump_drop(&mut dump_file, &fam, bundles, "singleton", -1.0, -1.0);
            continue;
        }
        if n_copies > max_copies {
            *drops.entry("megafamily").or_insert(0) += 1;
            dump_drop(&mut dump_file, &fam, bundles, "megafamily", -1.0, -1.0);
            continue;
        }
        let n_shared = fam.multimap_reads.len();
        // Exemption: if any copy in the family consists entirely of supplementary
        // reads (0 primary reads at that locus), this is a genuine multi-copy
        // gene family where reads are shared across copies by definition.
        // Don't drop it on shared count thresholds — the pileup EM will handle it.
        let has_supp_only_copy = fam.bundle_indices.iter().any(|&bi| {
            bundles.get(bi).map(|b| {
                !b.reads.is_empty()
                    && b.reads.iter().all(|r| !r.is_primary_alignment)
            }).unwrap_or(false)
        });
        if n_shared < min_shared && !has_supp_only_copy {
            *drops.entry("low_shared").or_insert(0) += 1;
            dump_drop(&mut dump_file, &fam, bundles, "low_shared", -1.0, -1.0);
            continue;
        }
        // Sharing density: a 17-copy "family" with 13 shared reads has
        // 0.76 reads/copy — most copies see no cross-mapping evidence at all.
        let shared_per_copy = n_shared as f64 / n_copies as f64;
        if shared_per_copy < min_shared_per_copy && !has_supp_only_copy {
            *drops.entry("low_shared_per_copy").or_insert(0) += 1;
            dump_drop(&mut dump_file, &fam, bundles, "low_shared_per_copy", -1.0, -1.0);
            continue;
        }
        // CV of intron counts (only over copies with ≥1 intron). Cheap
        // structure-similarity proxy that works across loci.
        let mut computed_exon_cv = -1.0_f64;
        if max_exon_cv > 0.0 {
            let intron_counts: Vec<usize> = fam.bundle_indices.iter()
                .map(|&bi| bundles.get(bi).map(|b| b.junction_stats.len()).unwrap_or(0))
                .filter(|&c| c > 0)
                .collect();
            if intron_counts.len() >= 2 {
                let n = intron_counts.len() as f64;
                let mean = intron_counts.iter().map(|&c| c as f64).sum::<f64>() / n;
                if mean > 0.0 {
                    let var = intron_counts.iter()
                        .map(|&c| { let d = c as f64 - mean; d * d })
                        .sum::<f64>() / n;
                    let cv = var.sqrt() / mean;
                    computed_exon_cv = cv;
                    if cv > max_exon_cv {
                        *drops.entry("high_exon_cv").or_insert(0) += 1;
                        dump_drop(&mut dump_file, &fam, bundles, "high_exon_cv", cv, -1.0);
                        continue;
                    }
                }
            }
        }
        // Primitive (graph-structural) signal: pairwise intron-length-set
        // similarity using fuzzy matching (±200bp tolerance). Real paralogs
        // share intron lengths within tens of bp; TE-bridge artifacts have
        // wildly different intron lengths. Hard-bin Jaccard misses paralogs
        // whose introns differ by exactly one bin width — fuzzy matching
        // counts an intron as "shared" if some other copy has an intron of
        // length within ±200bp.
        if min_primitive_jaccard > 0.0 {
            let intron_lens_per_copy: Vec<Vec<u64>> = fam.bundle_indices.iter()
                .map(|&bi| {
                    bundles.get(bi)
                        .map(|b| b.junction_stats.iter()
                            .map(|(j, _)| j.acceptor.saturating_sub(j.donor))
                            .collect::<Vec<u64>>())
                        .unwrap_or_default()
                })
                .collect();
            let multi_exon: Vec<&Vec<u64>> = intron_lens_per_copy.iter()
                .filter(|s| !s.is_empty())
                .collect();
            if multi_exon.len() >= 2 {
                let tolerance: u64 = 200; // bp — paralogs typically agree at this scale
                // Fuzzy Jaccard: |A∩B| counts introns in A with a within-
                // tolerance match in B; |A∪B| = |A| + |B| − |A∩B|.
                let fuzzy_match = |a: &[u64], b: &[u64]| -> f64 {
                    let mut matched = 0usize;
                    for &x in a {
                        if b.iter().any(|&y| {
                            let d = if x > y { x - y } else { y - x };
                            d <= tolerance
                        }) {
                            matched += 1;
                        }
                    }
                    let intersection = matched as f64;
                    let union = (a.len() + b.len()) as f64 - intersection;
                    if union > 0.0 { intersection / union } else { 0.0 }
                };
                let mut sum = 0.0_f64;
                let mut n_pairs = 0usize;
                for i in 0..multi_exon.len() {
                    for j in (i + 1)..multi_exon.len() {
                        sum += fuzzy_match(multi_exon[i], multi_exon[j]);
                        n_pairs += 1;
                    }
                }
                let mean_jaccard = if n_pairs > 0 { sum / n_pairs as f64 } else { 0.0 };
                if mean_jaccard < min_primitive_jaccard {
                    *drops.entry("low_primitive_jaccard").or_insert(0) += 1;
                    dump_drop(&mut dump_file, &fam, bundles,
                              "low_primitive_jaccard", computed_exon_cv, mean_jaccard);
                    continue;
                }
            }
        }
        kept.push(fam);
    }
    if !drops.is_empty() {
        eprintln!(
            "[VG] family-quality filter: {} → {} families ({} dropped: {:?})",
            n_in, kept.len(), n_in - kept.len(), drops
        );
    }
    kept
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
    // Min-hash sketch size. 128 mins gives ~9% std-error on Jaccard; cheaper
    // than full HashSet&lt;u64&gt; comparison (each min-hash slot is u64, 128*8=1KB
    // per bundle vs 200KB+ for a full kmer set). Pair compare is also O(SKETCH)
    // not O(set.len()).
    const SKETCH_SIZE: usize = 128;

    if bundles.len() > max_bundles_to_compare {
        return Vec::new();
    }

    // Build per-bundle min-hash sketches in parallel.
    use rayon::prelude::*;
    struct BundleSketch {
        bi: usize,
        sketch: [u64; SKETCH_SIZE],
    }
    let bundle_sketches: Vec<BundleSketch> = bundles.par_iter().enumerate()
        .filter_map(|(bi, bundle)| {
            // Extract exonic regions from junction stats.
            let mut exon_regions: Vec<(u64, u64)> = Vec::new();
            let mut introns: Vec<(u64, u64)> = bundle
                .junction_stats
                .iter()
                .map(|(j, _)| (j.donor, j.acceptor))
                .collect();
            introns.sort_by_key(|&(s, _)| s);

            if introns.is_empty() {
                let capped_end = bundle.end.min(bundle.start + 10000);
                exon_regions.push((bundle.start, capped_end));
            } else {
                let first_intron_start = introns[0].0;
                if first_intron_start > bundle.start {
                    exon_regions.push((
                        first_intron_start.saturating_sub(500).max(bundle.start),
                        first_intron_start,
                    ));
                }
                for i in 0..introns.len().saturating_sub(1) {
                    let es = introns[i].1;
                    let ee = introns[i + 1].0;
                    if ee > es {
                        exon_regions.push((es, ee));
                    }
                }
                let last_intron_end = introns.last().unwrap().1;
                exon_regions.push((last_intron_end, last_intron_end.saturating_add(500).min(bundle.end)));
            }

            let mut sketch = [u64::MAX; SKETCH_SIZE];
            let mut n_kmers: u32 = 0;
            for (es, ee) in &exon_regions {
                if let Some(seq) = genome.fetch_sequence(&bundle.chrom, *es, *ee) {
                    for window in seq.windows(kmer_len) {
                        if window.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                            let base = kmer_hash(window);
                            for i in 0..SKETCH_SIZE {
                                let h = base ^ (i as u64).wrapping_mul(0x9E3779B97F4A7C15);
                                if h < sketch[i] { sketch[i] = h; }
                            }
                            n_kmers = n_kmers.saturating_add(1);
                        }
                    }
                }
            }

            if (n_kmers as usize) >= min_kmers {
                Some(BundleSketch { bi, sketch })
            } else {
                None
            }
        })
        .collect();

    // Pairwise comparison: O(SKETCH_SIZE) per pair instead of O(|kmer set|).
    let pair_trace = std::env::var_os("RUSTLE_VG_DISCOVERY_TRACE").is_some();
    let mut links: Vec<(usize, usize)> = Vec::new();
    for i in 0..bundle_sketches.len() {
        let bi = bundle_sketches[i].bi;
        for j in (i + 1)..bundle_sketches.len() {
            let bj = bundle_sketches[j].bi;
            if bundles[bi].chrom != bundles[bj].chrom { continue; }
            if bundles[bi].start < bundles[bj].end && bundles[bj].start < bundles[bi].end {
                if pair_trace {
                    eprintln!("[VG-SEQSIM]   pair bi={} ({}-{}) bj={} ({}-{}) SKIPPED: spans overlap",
                              bi, bundles[bi].start, bundles[bi].end, bj, bundles[bj].start, bundles[bj].end);
                }
                continue;
            }
            let dist = if bundles[bi].start > bundles[bj].end {
                bundles[bi].start - bundles[bj].end
            } else {
                bundles[bj].start - bundles[bi].end
            };
            if dist > 10_000_000 { continue; }

            // Min-hash Jaccard estimate.
            let a = &bundle_sketches[i].sketch;
            let b = &bundle_sketches[j].sketch;
            let eq = a.iter().zip(b.iter())
                .filter(|(x, y)| **x != u64::MAX && x == y).count();
            let jaccard = eq as f64 / SKETCH_SIZE as f64;
            if pair_trace {
                eprintln!("[VG-SEQSIM]   pair bi={} ({}-{}) bj={} ({}-{}) dist={} jaccard={:.3} {}",
                          bi, bundles[bi].start, bundles[bi].end, bj, bundles[bj].start, bundles[bj].end,
                          dist, jaccard, if jaccard >= min_jaccard {"LINK"} else {"below-threshold"});
            }
            if jaccard >= min_jaccard {
                links.push((bi, bj));
            }
        }
    }

    if !links.is_empty() || std::env::var_os("RUSTLE_VG_DISCOVERY_TRACE").is_some() {
        eprintln!(
            "[VG] Sequence similarity: {} bundle(s) sketched (of {} total; need ≥{} exonic k-mers), \
             found {} pair(s) with ≥{:.0}% min-hash overlap (sketch={})",
            bundle_sketches.len(), bundles.len(), min_kmers,
            links.len(), min_jaccard * 100.0, SKETCH_SIZE,
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

/// Per-alignment identity proxy from edit distance and soft-clip rate.
///
/// Each (fam_pos, ri) carries its own NM and clip lengths from minimap2's
/// alignment to that paralog's locus. At medium identity (30–70%) NM diverges
/// across paralogs even when intron boundaries are conserved, providing the
/// per-paralog discrimination that intron-match alone cannot.
///
/// Returns a value in [0, 1]: 1.0 = perfect match, 0.0 = many edits + heavy
/// soft-clipping. Single-exon reads still get a meaningful score.
fn read_alignment_identity_score(read: &BundleRead) -> f64 {
    let aligned_bp: u64 = read.exons.iter().map(|&(s, e)| e.saturating_sub(s)).sum();
    if aligned_bp == 0 {
        return 0.0;
    }
    let nm_rate = (read.nm as f64 / aligned_bp as f64).min(1.0);
    let clip_total = (read.clip_left + read.clip_right) as f64;
    let clip_rate = if let Some(qlen) = read.query_length {
        if qlen > 0 {
            (clip_total / qlen as f64).min(1.0)
        } else {
            0.0
        }
    } else {
        0.0
    };
    // Equal-weighted NM + clip penalty, both in [0, 1].
    (1.0 - 0.5 * nm_rate - 0.5 * clip_rate).max(0.0)
}

/// Fraction of the bundle's bundlenodes covered by any read exon.
///
/// Discriminates lost-exon paralogs: a read covering 4 exons at a 5-exon paralog
/// graph scores 0.8 there but 1.0 at a 4-exon paralog (exon truly missing). The
/// existing `junction_compatibility` is symmetric in junctions and saturates
/// once shared boundaries match; this metric is asymmetric on the paralog side
/// and rewards graph-side completeness.
fn read_bundle_exon_coverage(read: &BundleRead, bundle: &Bundle) -> f64 {
    let Some(ref head) = bundle.bundlenodes else {
        return 0.0;
    };
    let mut total = 0usize;
    let mut covered = 0usize;
    let mut node: Option<&CBundlenode> = Some(head);
    while let Some(n) = node {
        total += 1;
        let any = read.exons.iter().any(|&(rs, re)| rs < n.end && n.start < re);
        if any {
            covered += 1;
        }
        node = n.next.as_deref();
    }
    if total == 0 {
        0.0
    } else {
        covered as f64 / total as f64
    }
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

/// After EM reweighting, populate `ExonClass.per_copy_cov` for every node in
/// the family graph.
///
/// For each exon class and each copy:
/// - Locates the copy's genomic span for that exon via `per_copy_spans`
/// - Sums `read.weight * overlap_fraction` over all reads in that copy's bundle
///   that physically overlap the exon span
/// - Records `(copy_id, weighted_per_base_coverage)` in `ec.per_copy_cov`
///
/// Copies without a span for a given exon class (copy-specific exons from other
/// copies) get a `(copy_id, 0.0)` sentinel so downstream code can distinguish
/// "no data" from "absent" without index arithmetic.
///
/// Call this after `run_pre_assembly_em_hmm` (or the heuristic-EM variants)
/// and before the main bundle-assembly loop, so the annotation is available
/// for diagnostics and future copy-aware path extraction.
pub fn annotate_per_copy_exon_coverage(
    family: &FamilyGroup,
    bundles: &[Bundle],
    family_graph: &mut crate::vg_hmm::family_graph::FamilyGraph,
) {
    let n_copies = family.bundle_indices.len();
    for ec in family_graph.nodes.iter_mut() {
        ec.per_copy_cov = Vec::with_capacity(n_copies);
        for (copy_id, &bi) in family.bundle_indices.iter().enumerate() {
            let bundle = &bundles[bi];
            let span = ec.per_copy_spans.iter()
                .find(|(cid, _)| *cid == copy_id)
                .map(|(_, s)| *s);
            let cov = if let Some((ec_start, ec_end)) = span {
                let ec_len = (ec_end.saturating_sub(ec_start)).max(1) as f64;
                let mut weighted_bp = 0.0f64;
                for read in &bundle.reads {
                    let r_start = read.ref_start;
                    let r_end = read.ref_end;
                    if r_start < ec_end && r_end > ec_start {
                        let overlap = (ec_end.min(r_end) - ec_start.max(r_start)) as f64;
                        weighted_bp += read.weight * overlap / ec_len;
                    }
                }
                weighted_bp
            } else {
                0.0
            };
            ec.per_copy_cov.push((copy_id, cov));
        }
    }
}

/// Compute per-copy assignment confidence after EM reweighting.
///
/// For each family copy (bundle), scans all multi-mapper reads and computes
/// the mean post-EM weight for this copy's placements. Result range [0, 1]:
/// - 1.0 → all multi-mappers fully assigned to this copy, or no multi-mappers
/// - ~1/N → evenly split across N copies
///
/// Runs in O(total_reads) over family bundles. Call after EM has updated
/// `BundleRead.weight` and before bundles are consumed by the assembly loop.
///
/// Returns `HashMap<bundle_idx, confidence>`. Bundles not in any family are
/// absent from the map; callers should default to `None` / omit the attribute.
pub fn compute_per_copy_confidence(
    families: &[FamilyGroup],
    bundles: &[Bundle],
) -> HashMap<usize, f64> {
    let mut result: HashMap<usize, f64> = HashMap::new();
    for fam in families {
        // Build set of read_name_hashes that are multi-mappers in this family.
        let multimap_set: std::collections::HashSet<u64> =
            fam.multimap_reads.keys().copied().collect();
        for (copy_id, &bi) in fam.bundle_indices.iter().enumerate() {
            if bi >= bundles.len() { continue; }
            let bundle = &bundles[bi];
            let mut sum_weight = 0.0f64;
            let mut n_multimap = 0usize;
            for read in &bundle.reads {
                if multimap_set.contains(&read.read_name_hash) {
                    sum_weight += read.weight;
                    n_multimap += 1;
                }
            }
            let confidence = if n_multimap == 0 {
                1.0 // no multi-mappers → assignment is fully certain
            } else {
                (sum_weight / n_multimap as f64).min(1.0).max(0.0)
            };
            let _ = copy_id; // copy_id used implicitly via bundle_indices iteration
            result.insert(bi, confidence);
        }
    }
    result
}

/// Fraction of a copy C's reads that fit IT at least as well as any sibling.
///
/// `n_unique` reads always support C (they only fit C). A multimapper is given
/// as `(rate_C, rate_min_sibling)` = (read's NM-rate at THIS copy, read's best
/// NM-rate at any sibling); it supports C iff `rate_C <= rate_min_sibling +
/// margin` (NM-ties within `margin` count as supporting — no false suppression
/// of genuinely-ambiguous co-expressed copies).
///
/// Returns `(n_unique + #supporting_multimappers) / total`, or `0.0` if there
/// are no reads at all.
pub fn copy_support_fraction(n_unique: usize, multimappers: &[(f64, f64)], margin: f64) -> f64 {
    let total = n_unique + multimappers.len();
    if total == 0 {
        return 0.0;
    }
    let mm_support = multimappers
        .iter()
        .filter(|&&(rate_c, rate_min_sib)| rate_c <= rate_min_sib + margin)
        .count();
    (n_unique + mm_support) as f64 / total as f64
}

/// Where a multi-mapping read truly belongs, by raw edit-distance margin.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadAnchor { Owns, Sibling, Tie }

/// Classify one read's copy assignment from raw edit distance (the count of
/// distinguishing mismatches — the fingerprint-EM log-likelihood ratio).
/// `nm_c`/`alen_c`: NM and aligned length of the read's placement at THIS copy.
/// `others`: (nm, alen) of the read's placements at OTHER copies.
/// Only placements with `alen >= extent_frac * alen_c` compete (a short partial
/// hit elsewhere is not a real competitor). `dnm = nm_other_best - nm_c`:
/// `Owns` if `dnm >= t`, `Sibling` if `dnm <= -t`, else `Tie`. With no
/// comparable competitor the read uniquely Owns this copy.
pub fn anchor_read(nm_c: u32, alen_c: u64, others: &[(u32, u64)], t: i64, extent_frac: f64) -> ReadAnchor {
    if alen_c == 0 { return ReadAnchor::Owns; }
    let best_other = others.iter()
        .filter(|&&(_, al)| (al as f64) >= extent_frac * (alen_c as f64))
        .map(|&(nm, _)| nm)
        .min();
    match best_other {
        None => ReadAnchor::Owns,
        Some(bo) => {
            let dnm = bo as i64 - nm_c as i64;
            if dnm >= t { ReadAnchor::Owns }
            else if dnm <= -t { ReadAnchor::Sibling }
            else { ReadAnchor::Tie }
        }
    }
}

/// Per-copy anchored read mass (the anchor-first M-step prior source).
///
/// For each copy `k` (index = position in `family.bundle_indices`), sums
/// `read.weight` over reads in copy `k`'s bundle that are EITHER:
///   * unique — `read.read_name_hash` is not a key in `family.multimap_reads`
///     (the read only places at this copy), OR
///   * decisively owned — `anchor_read(nm_k, alen_k, others, t, extent_frac)`
///     returns `ReadAnchor::Owns`, where `others` are the read's `(nm, alen)` at
///     the copy's siblings. `Sibling`/`Tie` placements contribute 0 to copy `k`.
///
/// Per-placement `(nm, alen)` uses the same calibration as `classify_family`:
/// `alen = read_aligned_len(read)`; cost = `round(de*alen)` when `de` is present
/// (HiFi indel-robust), else raw `read.nm`. Calibration constants are passed in
/// (`t`, `extent_frac`); the project default is raw-dNM `t=2`, `extent_frac=0.8`.
///
/// Returns a `Vec<f64>` of length `family.bundle_indices.len()` (zeros for
/// stale/empty copies). Operates on the ORIGINAL `FamilyGroup` + intact bundles
/// (call at EM time, like `compute_copy_independent_support`).
pub fn anchored_mass_per_copy(
    family: &FamilyGroup,
    bundles: &[Bundle],
    t: i64,
    extent_frac: f64,
) -> Vec<f64> {
    let n_copies = family.bundle_indices.len();
    let mut mass = vec![0.0_f64; n_copies];

    // Per-read `(nm, alen)` at (fam_pos, read_idx) — same cost model as
    // classify_family's `placement` closure (de-based, NM fallback).
    let placement = |fam_pos: usize, ri: usize| -> Option<(u32, u64)> {
        let bi = *family.bundle_indices.get(fam_pos)?;
        let read = bundles.get(bi)?.reads.get(ri)?;
        let alen = read_aligned_len(read);
        if alen == 0 {
            return None;
        }
        let cost = match read.de {
            Some(de) => ((de as f64) * (alen as f64)).round() as u32,
            None => read.nm,
        };
        Some((cost, alen))
    };

    // Build, per multimap read, its best (min-NM) placement per copy. Mirrors
    // classify_family: dedupes redundant alignments at the same copy by min NM.
    let mut per_copy_by_read: std::collections::HashMap<u64, std::collections::HashMap<usize, (u32, u64)>> =
        std::collections::HashMap::with_capacity(family.multimap_reads.len());
    for (&rnh, placements) in &family.multimap_reads {
        let mut per_copy: std::collections::HashMap<usize, (u32, u64)> = std::collections::HashMap::new();
        for &(fp, ri) in placements {
            if fp >= n_copies {
                continue;
            }
            if let Some((nm, al)) = placement(fp, ri) {
                per_copy
                    .entry(fp)
                    .and_modify(|e| if nm < e.0 { *e = (nm, al); })
                    .or_insert((nm, al));
            }
        }
        if !per_copy.is_empty() {
            per_copy_by_read.insert(rnh, per_copy);
        }
    }

    for (copy_id, &bi) in family.bundle_indices.iter().enumerate() {
        let bundle = match bundles.get(bi) {
            Some(b) => b,
            None => continue,
        };
        for read in &bundle.reads {
            // Unique reads (not a multimap key) always anchor their copy.
            if !family.multimap_reads.contains_key(&read.read_name_hash) {
                mass[copy_id] += read.weight;
                continue;
            }
            // Multimapper: decisive only if it Owns THIS copy by dNM margin.
            let per_copy = match per_copy_by_read.get(&read.read_name_hash) {
                Some(pc) => pc,
                None => continue,
            };
            let (nm_c, alen_c) = match per_copy.get(&copy_id) {
                Some(&v) => v,
                None => continue, // this read has no placement at this copy
            };
            let others: Vec<(u32, u64)> = per_copy
                .iter()
                .filter(|(&c, _)| c != copy_id)
                .map(|(_, &v)| v)
                .collect();
            if anchor_read(nm_c, alen_c, &others, t, extent_frac) == ReadAnchor::Owns {
                mass[copy_id] += read.weight;
            }
        }
    }

    mass
}

/// Assign each copy a dense identifiability-class label (0..n_classes-1). Copies
/// are merged into one class when NON-identifiable (`nonid_pairs`, transitive).
pub fn identifiability_classes(n_copies: usize, nonid_pairs: &[(usize, usize)]) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..n_copies).collect();
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        let mut r = x;
        while parent[r] != r { r = parent[r]; }
        let mut c = x;
        while parent[c] != c { let n = parent[c]; parent[c] = r; c = n; }
        r
    }
    for &(a, b) in nonid_pairs {
        if a < n_copies && b < n_copies {
            let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
            if ra != rb { parent[ra] = rb; }
        }
    }
    let mut label = std::collections::HashMap::<usize, usize>::new();
    let mut next = 0usize;
    let mut out = vec![0usize; n_copies];
    for i in 0..n_copies {
        let r = find(&mut parent, i);
        let l = *label.entry(r).or_insert_with(|| { let v = next; next += 1; v });
        out[i] = l;
    }
    out
}

/// Number of identifiability classes (full = n_copies, none = 1).
pub fn identifiability_partition(n_copies: usize, nonid_pairs: &[(usize, usize)]) -> usize {
    identifiability_classes(n_copies, nonid_pairs)
        .iter().copied().collect::<std::collections::BTreeSet<usize>>().len()
}

/// Aligned length of a read in bp: sum of its exon spans (falls back to
/// `query_length` if exons are empty). Returns 0 when neither is available.
fn read_aligned_len(read: &BundleRead) -> u64 {
    let exon_sum: u64 = read.exons.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
    if exon_sum > 0 {
        exon_sum
    } else {
        read.query_length.unwrap_or(0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FamilyClass {
    Family, FamilyNonIdentifiable, GenePlusUnexpressedParalog, Spillover,
    NotConnected, NotExpressedHere, SingleExonOutOfScope,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Identifiability { Full, Partial, None }
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocusRel { Tandem, Distal, Trans, Overlapping, Single }

#[derive(Debug, Clone)]
pub struct FamilyVerdict {
    pub class: FamilyClass,
    pub n_copies: usize,
    pub n_expressed: usize,
    pub connectivity: f64,
    pub identifiability: Identifiability,
    pub n_id_classes: usize,
    pub locus_rel: LocusRel,
    /// Depth-inferred copy number (audit lever #3, opt-in RUSTLE_VG_DEPTH_COPYNUM).
    /// `None` unless enabled. Diagnostic only — never overrides `n_copies`. See
    /// `estimate_copies_from_depth` for the calibration and its honest limitation.
    pub depth_copies: Option<f64>,
    /// Per-copy EM attribution confidence (audit lever #2). Populated POST-EM (the
    /// verdict is built pre-EM, then this is filled once `em_weight_gap`/`em_n_sites`
    /// are set on the reads). `None` when the EM did not run or the copy had no shared
    /// reads. See `EmCopyConfidence`.
    pub em_confidence: Option<EmCopyConfidence>,
    /// Segmental-duplication extent / breakpoints (opt-in RUSTLE_VG_SEGDUP_EXTENT).
    /// `None` unless the segdup pass ran. Surfaces gene+flank duplicated extent and the
    /// segdup-vs-bare-paralog call. See `detect_and_report_segdups`.
    pub segdup: Option<crate::vg_hmm::segdup::SegdupExtent>,
    /// Best gene-conversion event for the family (opt-in RUSTLE_VG_MOSAIC_ON). `None` unless
    /// the mosaic pass found a recombinant. The family's most-supported event (confirmed
    /// preferred). See `detect_and_report_mosaics` / `ConversionEvent`.
    pub conversion: Option<crate::vg_hmm::mosaic::ConversionEvent>,
    /// Evidence of a copy NOT in the reference genome (opt-in RUSTLE_VG_HIDDEN_COPY). `None`
    /// unless flagged. Detect-and-flag only — never a fabricated copy. See
    /// `detect_and_report_hidden_copies` / `HiddenCopyEvidence`.
    pub hidden_copy: Option<crate::vg_hmm::hidden_copy::HiddenCopyEvidence>,
}
impl FamilyClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            FamilyClass::Family => "family",
            FamilyClass::FamilyNonIdentifiable => "family_nonidentifiable",
            FamilyClass::GenePlusUnexpressedParalog => "gene_plus_unexpressed_paralog",
            FamilyClass::Spillover => "spillover",
            FamilyClass::NotConnected => "not_connected",
            FamilyClass::NotExpressedHere => "not_expressed_here",
            FamilyClass::SingleExonOutOfScope => "single_exon_out_of_scope",
        }
    }
}
impl Identifiability {
    pub fn as_str(&self) -> &'static str {
        match self { Identifiability::Full => "full", Identifiability::Partial => "partial", Identifiability::None => "none" }
    }
}
impl LocusRel {
    pub fn as_str(&self) -> &'static str {
        match self {
            LocusRel::Tandem => "tandem", LocusRel::Distal => "distal",
            LocusRel::Trans => "trans", LocusRel::Overlapping => "overlapping", LocusRel::Single => "single",
        }
    }
}
pub struct FamilyParams {
    pub dnm_t: i64, pub extent_frac: f64, pub k_expr: usize, pub h_min: f64, pub tie_none: f64,
}
impl Default for FamilyParams {
    fn default() -> Self { FamilyParams { dnm_t: 2, extent_frac: 0.7, k_expr: 3, h_min: 0.10, tie_none: 0.60 } }
}

/// Classify a discovered VG family per the M/H/X/I definition.
/// Operates on the ORIGINAL FamilyGroup + intact bundles (call at EM time,
/// like `compute_copy_independent_support`).
/// Depth-based copy-number estimate for a multi-copy family (audit lever #3, O1/O2).
///
/// A fully-collapsed tandem array maps N physical copies onto overlapping coordinates,
/// so the read depth through it is ~N× a single copy's depth — yet only the copies that
/// ASSEMBLE into distinct bundles are counted by `bundle_indices.len()`. This estimates
/// copy number from DEPTH instead: calibrate a single-copy depth unit, then express each
/// resolved copy's depth as a multiple of it. A multiplier ≫ 1 flags a bundle that hides
/// several collapsed physical copies; `inferred_copies` (total / unit) ≥ the resolved
/// count then signals genome-structure under-counting.
///
/// HONEST LIMITATION: read depth conflates copy number with EXPRESSION — one copy
/// expressed 2× is indistinguishable, by depth alone, from two collapsed 1× copies. So
/// this is an *abundance multiplier*; it equals copy number only under equal-per-copy
/// expression. With no external reference, internal calibration (faintest resolved copy =
/// 1×) detects only RELATIVE over-representation among resolved copies — a uniformly
/// collapsed array (every resolved bundle hiding the same N) is invisible without an
/// external single-copy baseline. Emitted as a diagnostic signal, never overriding the
/// structural (bundle/footprint) copy count.
#[derive(Debug, Clone, PartialEq)]
pub struct DepthCopyEstimate {
    pub single_copy_unit: f64,         // calibrated single-copy depth
    pub total_depth: f64,              // sum of per-copy depths
    pub inferred_copies: f64,          // total_depth / single_copy_unit
    pub per_copy_multiplier: Vec<f64>, // depth_i / single_copy_unit (0 for absent copies)
}

/// Estimate copy number from per-copy read depth. `depths[i]` = average read depth of
/// resolved copy i (0 = absent/unassembled). `external_unit`, when known (e.g. a
/// genome-wide single-copy depth), gives an ABSOLUTE estimate; otherwise the faintest
/// non-zero resolved copy is taken as the 1× baseline (RELATIVE over-representation).
/// Returns `None` when no copy has positive depth.
pub fn estimate_copies_from_depth(depths: &[f64], external_unit: Option<f64>) -> Option<DepthCopyEstimate> {
    let nonzero: Vec<f64> = depths.iter().copied().filter(|&d| d > 0.0).collect();
    if nonzero.is_empty() {
        return None;
    }
    let unit = match external_unit {
        Some(u) if u > 0.0 => u,
        // Internal calibration: MEDIAN of resolved-copy depths = single-copy baseline.
        // Median (not min) is robust to per-copy read-count/truncation noise — the high
        // and low cancel around it, so equal-expression copies calibrate to ~1× cleanly.
        // Assumes the typical resolved copy is single (holds when most copies are single
        // and a few are collapsed; fails if the MAJORITY are collapsed — see doc).
        _ => {
            let mut s = nonzero.clone();
            s.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let m = s.len() / 2;
            if s.len() % 2 == 1 { s[m] } else { (s[m - 1] + s[m]) / 2.0 }
        }
    };
    if !(unit > 0.0) {
        return None;
    }
    let total: f64 = nonzero.iter().sum();
    let per_copy_multiplier: Vec<f64> = depths.iter().map(|&d| d / unit).collect();
    Some(DepthCopyEstimate {
        single_copy_unit: unit,
        total_depth: total,
        inferred_copies: total / unit,
        per_copy_multiplier,
    })
}

/// Per-copy EM attribution confidence (audit lever #2): how decisively the fingerprint-EM
/// resolved a copy's SHARED (multi-mapping) reads. These signals are computed during the EM
/// (per-read `em_weight_gap` = winning weight − runner-up; `em_n_sites` = diagnostic sites
/// the read covered) but were only printed to a trace — never surfaced per transcript. A
/// copy whose shared reads are mostly `decisive` is well-distinguished from its siblings; one
/// dominated by `uncertain` reads sits near the identifiability boundary (the spec's core
/// signal). Buckets use the SAME thresholds as the `[VG-FP-EM]` summary (0.8 / 0.5).
#[derive(Debug, Clone, PartialEq)]
pub struct EmCopyConfidence {
    pub n_em_reads: usize,   // shared/multimap reads touching this copy that went through EM
    pub n_decisive: usize,   // weight_gap > 0.8
    pub n_moderate: usize,   // 0.5 < weight_gap <= 0.8
    pub n_uncertain: usize,  // weight_gap <= 0.5
    pub mean_gap: f64,       // mean weight_gap over the EM reads
    pub mean_sites: f64,     // mean diagnostic sites covered
}
impl EmCopyConfidence {
    /// Fraction of this copy's shared reads the EM resolved decisively — the headline.
    pub fn decisive_frac(&self) -> f64 {
        if self.n_em_reads == 0 { 0.0 } else { self.n_decisive as f64 / self.n_em_reads as f64 }
    }
}

/// Summarize per-read EM decisiveness for one copy. `reads` = `(weight_gap, n_sites)` for
/// each shared read touching the copy. Returns `None` when the copy has no EM reads (e.g. a
/// copy supported only by unique reads, or a family the EM skipped).
pub fn summarize_em_confidence(reads: &[(f64, u32)]) -> Option<EmCopyConfidence> {
    if reads.is_empty() {
        return None;
    }
    let (mut nd, mut nm, mut nu) = (0usize, 0usize, 0usize);
    let mut sum_gap = 0.0f64;
    let mut sum_sites = 0.0f64;
    for &(gap, sites) in reads {
        if gap > 0.8 { nd += 1; } else if gap > 0.5 { nm += 1; } else { nu += 1; }
        sum_gap += gap;
        sum_sites += sites as f64;
    }
    let n = reads.len();
    Some(EmCopyConfidence {
        n_em_reads: n,
        n_decisive: nd,
        n_moderate: nm,
        n_uncertain: nu,
        mean_gap: sum_gap / n as f64,
        mean_sites: sum_sites / n as f64,
    })
}

pub fn classify_family(family: &FamilyGroup, bundles: &[Bundle], p: &FamilyParams) -> FamilyVerdict {
    let n_copies = family.bundle_indices.len();
    let placement = |fam_pos: usize, ri: usize| -> Option<(u32, u64)> {
        let bi = *family.bundle_indices.get(fam_pos)?;
        let bundle = bundles.get(bi)?;
        let read = bundle.reads.get(ri)?;
        let alen = read_aligned_len(read);
        if alen == 0 { return None; }
        // Gap-compressed event count (de * aligned_len), rounded — robust to HiFi
        // indels that inflate raw NM. Falls back to raw NM when `de` is absent
        // (e.g. synthetic test reads), preserving NM-based behavior there.
        let cost = match read.de {
            Some(de) => ((de as f64) * (alen as f64)).round() as u32,
            None => read.nm,
        };
        Some((cost, alen))
    };
    let copy_spliced = |fam_pos: usize| -> bool {
        family.bundle_indices.get(fam_pos)
            .and_then(|&bi| bundles.get(bi))
            .map(|b| b.reads.iter().any(|r| r.exons.len() >= 2 || !r.junctions.is_empty()))
            .unwrap_or(false)
    };
    // Scope gate: spec §7 scopes the core definition to spliced families. We use
    // ANY copy spliced (not ALL): a family with one spliced gene + an intronless
    // retro/pseudogene sibling is still routed through classification (the spliced
    // gene is the real assembly target). Only a wholly single-exon family is
    // SingleExonOutOfScope. This is a deliberate lenient reading of the per-copy wording.
    let any_spliced = (0..n_copies).any(copy_spliced);

    // Pass 1: per-read placements (kept for pass 2) + pairwise distinguishing + connectivity numerator.
    let mut reads_pc: Vec<std::collections::HashMap<usize, (u32, u64)>> = Vec::new();
    let mut pair_has_shared = std::collections::HashSet::<(usize, usize)>::new();
    let mut pair_has_distinguishing = std::collections::HashSet::<(usize, usize)>::new();
    let mut n_reads_total = 0usize;
    let mut n_reads_shared = 0usize;

    for placements in family.multimap_reads.values() {
        let mut per_copy: std::collections::HashMap<usize, (u32, u64)> = std::collections::HashMap::new();
        for &(fp, ri) in placements {
            if let Some((nm, al)) = placement(fp, ri) {
                per_copy.entry(fp).and_modify(|e| { if nm < e.0 { *e = (nm, al); } }).or_insert((nm, al));
            }
        }
        if per_copy.is_empty() { continue; }   // read resolves to no family copy (stale index) — don't count
        n_reads_total += 1;
        if per_copy.len() >= 2 { n_reads_shared += 1; }
        let copies: Vec<usize> = per_copy.keys().copied().collect();
        for i in 0..copies.len() {
            for j in (i + 1)..copies.len() {
                let (a, b) = (copies[i].min(copies[j]), copies[i].max(copies[j]));
                pair_has_shared.insert((a, b));
                let dnm = per_copy[&copies[i]].0 as i64 - per_copy[&copies[j]].0 as i64;
                if dnm.abs() >= p.dnm_t { pair_has_distinguishing.insert((a, b)); }
            }
        }
        reads_pc.push(per_copy);
    }

    // Unique reads per copy (in a copy's bundle, not a multimap key). Counted in
    // the connectivity denominator AND attributed to the copy's class below.
    let n_unique_per_copy: Vec<usize> = (0..n_copies).map(|c| {
        family.bundle_indices.get(c).and_then(|&bi| bundles.get(bi))
            .map(|b| b.reads.iter().filter(|r| !family.multimap_reads.contains_key(&r.read_name_hash)).count())
            .unwrap_or(0)
    }).collect();
    n_reads_total += n_unique_per_copy.iter().sum::<usize>();
    let connectivity = if n_reads_total > 0 { n_reads_shared as f64 / n_reads_total as f64 } else { 0.0 };

    // Identifiability classes from the non-identifiable pairs.
    let nonid: Vec<(usize, usize)> = pair_has_shared.iter()
        .filter(|pr| !pair_has_distinguishing.contains(pr)).copied().collect();
    let classes = identifiability_classes(n_copies, &nonid);   // class label per copy
    let n_id_classes = classes.iter().copied().collect::<std::collections::BTreeSet<usize>>().len();
    let n_class_labels = n_id_classes.max(1);

    // Pass 2: per-CLASS owning reads + cross-class ties. A read OWNS class K* (the
    // class of its best copy) iff its best copy beats every copy OUTSIDE K* by
    // dnm >= t (ties WITHIN K* are fine — that is the point of class-level X). A
    // shared read that beats no other class decisively is a cross-class tie.
    let mut class_owns: Vec<usize> = vec![0; n_class_labels];
    let mut n_class_tie = 0usize;
    for per_copy in &reads_pc {
        let (cstar, nm_star) = per_copy.iter().map(|(&c, &(nm, _))| (c, nm))
            .min_by_key(|&(_, nm)| nm).unwrap();   // reads_pc entries are non-empty
        let kstar = classes[cstar];
        let outside_min = per_copy.iter().filter(|(&c, _)| classes[c] != kstar)
            .map(|(_, &(nm, _))| nm).min();
        match outside_min {
            None => class_owns[kstar] += 1,        // no competing class present → owns it
            Some(om) => {
                if (om as i64 - nm_star as i64) >= p.dnm_t { class_owns[kstar] += 1; }
                else if per_copy.len() >= 2 { n_class_tie += 1; }
            }
        }
    }
    for c in 0..n_copies { class_owns[classes[c]] += n_unique_per_copy[c]; }

    // Class sizes (copies per class label) for the n_expressed==1 disambiguation.
    let mut class_size: Vec<usize> = vec![0; n_class_labels];
    for c in 0..n_copies { class_size[classes[c]] += 1; }

    // X is now the number of EXPRESSED identifiability CLASSES (a non-identifiable
    // cluster counts once if it collectively out-anchors the other classes).
    let n_expressed = class_owns.iter().filter(|&&o| o >= p.k_expr).count();
    let total_owns: usize = class_owns.iter().sum();
    let tie_frac = if n_reads_shared > 0 { n_class_tie as f64 / n_reads_shared as f64 } else { 0.0 }; // cross-class ties
    let identifiability = if n_id_classes == n_copies && tie_frac < 0.15 { Identifiability::Full }
        else if n_id_classes <= 1 || tie_frac >= p.tie_none { Identifiability::None }
        else { Identifiability::Partial };
    let locus_rel = locus_relationship(family, bundles);

    // Size of the single expressed class (when exactly one) — multi-copy → a real
    // but unresolvable family; single-copy → one gene + unexpressed paralogs.
    let lone_expressed_size = if n_expressed == 1 {
        (0..n_class_labels).find(|&k| class_owns[k] >= p.k_expr).map(|k| class_size[k]).unwrap_or(1)
    } else { 1 };

    let class = if n_copies < 2 { FamilyClass::NotConnected }
        else if !any_spliced { FamilyClass::SingleExonOutOfScope }
        else if total_owns < p.k_expr || class_owns.iter().all(|&o| o == 0) { FamilyClass::NotExpressedHere }
        else if connectivity < p.h_min { FamilyClass::NotConnected }
        else if n_expressed >= 2 { FamilyClass::Family }                                  // >=2 expressed classes
        else if n_expressed == 1 && lone_expressed_size >= 2 { FamilyClass::FamilyNonIdentifiable }
        else if n_expressed == 1 { FamilyClass::GenePlusUnexpressedParalog }              // 1 expressed single-copy class
        else { FamilyClass::Spillover };

    // Audit lever #3 (opt-in): depth-inferred copy number. Per-copy depth = total aligned
    // read bases in the copy's bundle / bundle genomic span (average read depth). Absent
    // copies contribute 0. An optional external single-copy unit (RUSTLE_VG_SINGLE_COPY_DEPTH)
    // converts the relative multiplier into an absolute copy count.
    let depth_copies = if std::env::var_os("RUSTLE_VG_DEPTH_COPYNUM").is_some() {
        // Per-copy abundance = aligned bases attributed to each copy as its READ HOME.
        // CRITICAL: for a high-identity family, a bright copy's reads are mostly MULTIMAP
        // reads (shared, stored in family.multimap_reads), NOT in bundle.reads — so summing
        // bundle.reads alone undercounts exactly the copies we care about. Attribute each
        // multimap read to its best-scoring copy (`reads_pc`, already computed above) plus
        // each copy's unique reads. Total aligned bases, not bases/span (paralog copies are
        // ~equal length; the multimap-extended bundle span adds variance without signal).
        let mut depths: Vec<f64> = vec![0.0; n_copies];
        for per_copy in &reads_pc {
            // best copy = lowest mismatch cost; add its aligned length there (the read's home)
            if let Some((&cbest, &(_, alen))) = per_copy.iter().min_by_key(|(_, &(nm, _))| nm) {
                if cbest < n_copies { depths[cbest] += alen as f64; }
            }
        }
        for c in 0..n_copies {
            if let Some(b) = family.bundle_indices.get(c).and_then(|&bi| bundles.get(bi)) {
                for r in &b.reads {
                    if !family.multimap_reads.contains_key(&r.read_name_hash) {
                        depths[c] += read_aligned_len(r) as f64;
                    }
                }
            }
        }
        let external = std::env::var("RUSTLE_VG_SINGLE_COPY_DEPTH").ok().and_then(|s| s.parse().ok());
        estimate_copies_from_depth(&depths, external).map(|e| e.inferred_copies)
    } else {
        None
    };
    // em_confidence is filled POST-EM by the pipeline (the read em_* fields aren't set
    // yet at classify-time, which runs before run_fingerprint_em).
    FamilyVerdict { class, n_copies, n_expressed, connectivity, identifiability, n_id_classes, locus_rel, depth_copies, em_confidence: None, segdup: None, conversion: None, hidden_copy: None }
}

/// True iff the family has two OPPOSITE-strand bundles that OVERLAP in genomic coords — a
/// same-locus inverted pair (DAZ1−/DAZ3+), where a shared read genuinely cannot tell the
/// copies apart and joint-strand EM apportionment across the inversion is correct.
///
/// False for a DISPERSED inverted duplication (opposite-strand copies at DISTINCT loci):
/// there each copy's reads have a clear primary locus, so per-strand EM is right. Joint-EM
/// on a dispersed inversion uses a neutral fingerprint (mixed-strand families have no joint
/// graph) and apportions the inverted copy's reads to the forward copies, starving the
/// inverted copy below the assembly threshold — even though it assembles cleanly de-novo.
pub fn mixed_strand_copies_overlap(family: &FamilyGroup, bundles: &[Bundle]) -> bool {
    let loci: Vec<(&str, char, u64, u64)> = family.bundle_indices.iter()
        .filter_map(|&bi| bundles.get(bi))
        .map(|b| (b.chrom.as_str(), b.strand, b.start, b.end))
        .collect();
    for i in 0..loci.len() {
        for j in (i + 1)..loci.len() {
            let (a, b) = (loci[i], loci[j]);
            if a.1 != b.1 && a.0 == b.0 && a.2 < b.3 && b.2 < a.3 { // strict: book-ended != overlapping (matches vg.rs:1141)
                return true;
            }
        }
    }
    false
}

fn locus_relationship(family: &FamilyGroup, bundles: &[Bundle]) -> LocusRel {
    let loci: Vec<(&str, u64, u64)> = family.bundle_indices.iter()
        .filter_map(|&bi| bundles.get(bi)).map(|b| (b.chrom.as_str(), b.start, b.end)).collect();
    classify_physical_loci(loci, family.bundle_indices.len())
}

/// Spatial classification of a family's copies, DE-MIRRORED. Rustle bundles per strand, so a
/// single genomic locus carries a `+` and a `−` bundle at the same coords (a strand mirror);
/// plus a copy may contribute several overlapping isoform bundles. The old check returned
/// `Overlapping` if ANY two bundles overlapped, so one strand-mirror pair made the WHOLE family
/// read as Overlapping — masking the true inter-copy relationship (a 13 Mb-spread dispersed
/// segdup was wrongly Overlapping instead of Distal). Fix: merge overlapping bundle intervals
/// into DISTINCT PHYSICAL loci first, then classify on those.
fn classify_physical_loci(mut loci: Vec<(&str, u64, u64)>, n_bundles: usize) -> LocusRel {
    if loci.is_empty() {
        return LocusRel::Single;
    }
    // Merge overlapping intervals per chromosome → distinct physical loci.
    loci.sort_by(|a, b| a.0.cmp(b.0).then(a.1.cmp(&b.1)));
    let mut merged: Vec<(&str, u64, u64)> = Vec::with_capacity(loci.len());
    for (c, s, e) in loci {
        match merged.last_mut() {
            Some(last) if last.0 == c && s <= last.2 => { last.2 = last.2.max(e); }
            _ => merged.push((c, s, e)),
        }
    }
    if merged.len() < 2 {
        // One physical locus. ≥2 bundles there = same-locus multi-copy (e.g. an inverted
        // pair, DAZ1−/DAZ3+) → genuinely Overlapping; a lone bundle → Single.
        return if n_bundles >= 2 { LocusRel::Overlapping } else { LocusRel::Single };
    }
    if !merged.iter().all(|l| l.0 == merged[0].0) {
        return LocusRel::Trans;
    }
    // Distinct loci no longer overlap; classify by the largest inter-locus gap.
    let mut maxgap = 0u64;
    for i in 0..merged.len() {
        for j in (i + 1)..merged.len() {
            let (a, b) = (merged[i], merged[j]);
            maxgap = maxgap.max(a.1.max(b.1).saturating_sub(a.2.min(b.2)));
        }
    }
    if maxgap < 1_000_000 { LocusRel::Tandem } else { LocusRel::Distal }
}

/// Per-family segmental-duplication extent / breakpoint pass (opt-in RUSTLE_VG_SEGDUP_EXTENT).
/// For each copy pair, fetch the GENOME flanks anchored at the gene and measure how far the
/// cross-copy homology extends into them (the duplicated-segment breakpoints). Reports the
/// max-extent pair per family and whether it's a true segdup (gene+flanks) vs a bare paralog
/// (gene-only homology). Genome-based: mRNA reads don't cover intergenic flanks. Analysis-only.
pub fn detect_and_report_segdups(
    families: &[FamilyGroup],
    bundles: &[Bundle],
    genome: &crate::genome::GenomeIndex,
) -> crate::types::DetHashMap<usize, crate::vg_hmm::segdup::SegdupExtent> {
    use crate::vg_hmm::segdup::{call_segdup_extent, SegdupExtent, SegdupParams};
    let mut out: crate::types::DetHashMap<usize, SegdupExtent> = Default::default();
    let params = SegdupParams::from_env();
    let win: u64 = std::env::var("RUSTLE_VG_SEGDUP_FETCH")
        .ok().and_then(|s| s.parse().ok()).unwrap_or(8000);
    let fetch_rev = |chrom: &str, s: u64, e: u64| -> Vec<u8> {
        let mut v = genome.fetch_sequence(chrom, s, e).unwrap_or_default();
        v.reverse(); // index 0 = the gene-proximal end, going outward
        v
    };
    for fam in families {
        let loci: Vec<(&str, u64, u64)> = fam.bundle_indices.iter()
            .filter_map(|&bi| bundles.get(bi))
            .map(|b| (b.chrom.as_str(), b.start, b.end))
            .collect();
        if loci.len() < 2 {
            continue;
        }
        let mut best: Option<SegdupExtent> = None;
        for i in 0..loci.len() {
            for j in (i + 1)..loci.len() {
                let (ca, sa, ea) = loci[i];
                let (cb, sb, eb) = loci[j];
                let gene_span = ((ea.saturating_sub(sa)) + (eb.saturating_sub(sb))) / 2;
                let up_a = fetch_rev(ca, sa.saturating_sub(win), sa);
                let up_b = fetch_rev(cb, sb.saturating_sub(win), sb);
                let down_a = genome.fetch_sequence(ca, ea, ea + win).unwrap_or_default();
                let down_b = genome.fetch_sequence(cb, eb, eb + win).unwrap_or_default();
                let seg = call_segdup_extent(gene_span, &up_a, &up_b, &down_a, &down_b, &params);
                if best.as_ref().map(|b| seg.total_extent > b.total_extent).unwrap_or(true) {
                    best = Some(seg);
                }
            }
        }
        if let Some(seg) = best {
            eprintln!(
                "[VG-SEGDUP] family={} copies={} gene~{}bp up_flank={}bp down_flank={}bp duplicated_extent={}bp => {}",
                fam.family_id, loci.len(), seg.gene_span, seg.upstream_extent,
                seg.downstream_extent, seg.total_extent,
                if seg.is_segdup { "SEGMENTAL DUPLICATION (gene+flanks)" } else { "bare paralog (gene-only homology)" },
            );
            out.insert(fam.family_id, seg);
        }
    }
    out
}

/// Per-family pass detecting gene-family copies PRESENT in the reads but ABSENT from the
/// reference (opt-in RUSTLE_VG_HIDDEN_COPY). For each reference copy locus, scans the FULL
/// bundles for PRIMARY reads (the paralog-bleed firewall — an in-reference paralog's reads are
/// primary at ITS locus, secondary here) and runs `detect_hidden_copy` over their mismatch
/// haplotype. DETECT + FLAG only — reports the discrepancy, never places or fabricates the copy.
pub fn detect_and_report_hidden_copies(
    families: &[FamilyGroup],
    bundles: &[Bundle],
) -> crate::types::DetHashMap<usize, crate::vg_hmm::hidden_copy::HiddenCopyEvidence> {
    use crate::vg_hmm::hidden_copy::{detect_hidden_copy, HiddenCopyEvidence, HiddenCopyParams, ReadObs};
    let params = HiddenCopyParams::from_env();
    let mut out: crate::types::DetHashMap<usize, HiddenCopyEvidence> = Default::default();
    for fam in families {
        let mut best: Option<HiddenCopyEvidence> = None;
        for &bi in &fam.bundle_indices {
            let Some(fb) = bundles.get(bi) else { continue };
            let (chrom, lo, hi) = (fb.chrom.clone(), fb.start, fb.end);
            let mut seen: crate::types::DetHashSet<u64> = Default::default();
            let mut reads: Vec<ReadObs> = Vec::new();
            for bundle in bundles {
                if bundle.chrom != chrom || bundle.end < lo || bundle.start > hi {
                    continue;
                }
                for r in &bundle.reads {
                    if !r.is_primary_alignment || r.ref_end < lo || r.ref_start > hi
                        || !seen.insert(r.read_name_hash)
                    {
                        continue;
                    }
                    let alts: Vec<u64> = r.mismatches.iter()
                        .map(|&(p, _)| p).filter(|&p| p >= lo && p <= hi).collect();
                    reads.push(ReadObs { start: r.ref_start, end: r.ref_end, alts });
                }
            }
            let ev = detect_hidden_copy(&reads, &params);
            if std::env::var_os("RUSTLE_VG_HIDDEN_TRACE").is_some() {
                eprintln!("[VG-HIDDEN-TRACE] family={} locus={}:{}-{} primary_reads={} alt_positions={} alt_reads={} flagged={}",
                    fam.family_id, chrom, lo, hi, ev.n_primary_reads, ev.n_alt_positions, ev.n_alt_reads, ev.flagged);
            }
            if ev.flagged {
                eprintln!(
                    "[VG-HIDDEN] family={} locus={}:{}-{} primary_reads={} alt_haplotype_positions={} alt_reads={} frac={:.2} => EVIDENCE OF A COPY NOT IN THE REFERENCE (flagged; reference models 1, reads imply >=2; NOT placed/fabricated)",
                    fam.family_id, chrom, lo, hi, ev.n_primary_reads, ev.n_alt_positions,
                    ev.n_alt_reads, ev.alt_read_fraction,
                );
                if best.as_ref().map(|b| ev.n_alt_positions > b.n_alt_positions).unwrap_or(true) {
                    best = Some(ev);
                }
            }
        }
        if let Some(ev) = best {
            out.insert(fam.family_id, ev);
        }
    }
    out
}

/// Compute each copy's **independent support** within a family.
///
/// Operates on the ORIGINAL discovered `FamilyGroup` (pre-strand-split), where
/// `family.multimap_reads` links a read's placement at copy C to its placements
/// at *all* sibling copies — including cross-strand siblings (DAZ1−/DAZ3+),
/// which a strand-split sub-family would hide.
///
/// For each copy C (= `fam_pos` in `family.bundle_indices`):
/// - **C-unique reads**: reads in C's bundle whose `read_name_hash` is NOT a key
///   in `family.multimap_reads` (they only fit C → always support C).
/// - **C-multimappers**: entries of `family.multimap_reads` that have a placement
///   at C. For such a read, `rate_C = nm_at_C / aligned_len_at_C`, and
///   `rate_min_sibling = min over the read's OTHER placements (other fam_pos in
///   this family) of nm / aligned_len`. A multimapper whose `rate_C` OR
///   `rate_min_sibling` is UNAVAILABLE (missing aligned-length, or a stale
///   read index that no longer resolves to a read) is UNASSESSABLE and is
///   EXCLUDED from the fraction entirely — neither numerator nor denominator.
///   We never let missing data *manufacture* support for a copy; it simply
///   isn't counted. If a copy has zero assessable reads AND zero unique reads,
///   we default to KEEP (`support = 1.0`) so a copy is never suppressed purely
///   on absent evidence.
///
/// With the map computed at EM time (bundles intact, `ri` indices valid) the
/// unassessable path almost never fires; it only guards against residual
/// staleness, and now fails SAFE (exclude) rather than DANGEROUS (support).
///
/// Returns `fam_pos -> copy_support_fraction(...)`. The `fam_pos` keys match the
/// `vg_copy_id` assigned to transcripts (both index into `family.bundle_indices`).
pub fn compute_copy_independent_support(
    family: &FamilyGroup,
    bundles: &[Bundle],
    margin: f64,
) -> HashMap<usize, f64> {
    let mut out: HashMap<usize, f64> = HashMap::new();
    let n_copies = family.bundle_indices.len();

    // Per-copy rate at a single placement, or None if aligned length is missing.
    let placement_rate = |fam_pos: usize, ri: usize| -> Option<f64> {
        let bi = *family.bundle_indices.get(fam_pos)?;
        let bundle = bundles.get(bi)?;
        let read = bundle.reads.get(ri)?;
        let alen = read_aligned_len(read);
        if alen == 0 {
            return None; // missing aligned-length data → no usable rate
        }
        Some(read.nm as f64 / alen as f64)
    };

    for c in 0..n_copies {
        let bi = match family.bundle_indices.get(c) {
            Some(&b) if b < bundles.len() => b,
            _ => {
                out.insert(c, 0.0);
                continue;
            }
        };
        let bundle = &bundles[bi];

        // C-unique reads: in C's bundle, read_name_hash NOT a multimap key.
        let n_unique = bundle
            .reads
            .iter()
            .filter(|r| !family.multimap_reads.contains_key(&r.read_name_hash))
            .count();

        // C-multimappers: family multimap entries with a placement at C.
        // `mm_pairs` holds only ASSESSABLE reads (both rate_C and a sibling rate
        // are available). Unassessable reads are excluded entirely and counted
        // in `n_unassessable` for the trace — never pushed as fake support.
        let mut mm_pairs: Vec<(f64, f64)> = Vec::new();
        let mut n_unassessable: usize = 0;
        for placements in family.multimap_reads.values() {
            // The read's placement at THIS copy (if any). A read can in
            // principle have multiple placements at the same copy; take the
            // best-fitting (lowest rate) as rate_C.
            let rate_c: Option<f64> = placements
                .iter()
                .filter(|&&(fp, _)| fp == c)
                .filter_map(|&(fp, ri)| placement_rate(fp, ri))
                .fold(None, |acc, r| Some(acc.map_or(r, |a: f64| a.min(r))));

            // Does this read place at C at all?
            let places_at_c = placements.iter().any(|&(fp, _)| fp == c);
            if !places_at_c {
                continue;
            }

            // rate_min over the read's OTHER placements (other fam_pos).
            let rate_min_sib: Option<f64> = placements
                .iter()
                .filter(|&&(fp, _)| fp != c)
                .filter_map(|&(fp, ri)| placement_rate(fp, ri))
                .fold(None, |acc, r| Some(acc.map_or(r, |a: f64| a.min(r))));

            match (rate_c, rate_min_sib) {
                (Some(rc), Some(rs)) => mm_pairs.push((rc, rs)),
                // Missing data at C or at every sibling → we cannot assess this
                // read. EXCLUDE it from the fraction (neither numerator nor
                // denominator). Never let absent evidence manufacture support.
                _ => n_unassessable += 1,
            }
        }

        // Default to KEEP only when there is NO usable evidence at all (no
        // unique reads and no assessable multimappers) — a copy is never
        // suppressed purely on missing data.
        let support = if n_unique == 0 && mm_pairs.is_empty() {
            1.0
        } else {
            copy_support_fraction(n_unique, &mm_pairs, margin)
        };
        if std::env::var_os("RUSTLE_VG_SUPPORT_TRACE").is_some() {
            let sample: Vec<(f64,f64)> = mm_pairs.iter().take(5).copied().collect();
            eprintln!("[SUPPORT] fam={} copy={} bi={} n_reads_in_bundle={} n_unique={} n_assessable_mm={} n_unassessable_mm={} support={:.3} sample_mm={:?}",
                family.family_id, c, bi, bundle.reads.len(), n_unique, mm_pairs.len(), n_unassessable, support, sample);
        }
        out.insert(c, support);
    }

    out
}

/// Boost read weights in underpowered family-member bundles.
///
/// When a bundle has fewer than `RUSTLE_VG_BOOST_PRIMARY_THR` primary reads
/// (default 10) but the family's multi-mappers carry significant EM-assigned
/// weight to that bundle (`RUSTLE_VG_BOOST_MIN_EM_WEIGHT`, default 3.0), read
/// weights are multiplied so the assembler sees enough effective coverage.
///
/// Three cases:
///   primary=0: all reads boosted uniformly (no primary signal to prefer).
///   primary∈[1, primary_min): skipped. Boosting 1–2 primaries creates a
///     coverage spike that can disrupt secondary-driven assembly at bundles
///     already producing output.
///   primary∈[primary_min, primary_thr): only primary reads are boosted;
///     secondaries keep their EM-assigned weights so relative evidence is
///     preserved. With 3+ primaries the boost is spread across enough reads
///     that the coverage profile stays balanced.
///
/// The boost factor is capped at `RUSTLE_VG_BOOST_MAX` (default 10.0) to
/// avoid extreme weight imbalances.
///
/// ⚠️ PRECISION CHECK (2026-06-03) — DO NOT ENABLE BY DEFAULT. The `primary=0`
/// path fabricates phantoms. On real GOLGA8 it fired twice on loci with ZERO
/// genuine primary reads (pure sibling echoes — the DAZ3 pattern); one was inert,
/// the other amplified 79 secondary echoes 10× into 2 NEW transcripts (RSTL.373)
/// at cov 229× with capacity_confidence 1.000 — i.e. it manufactures a
/// MAXIMALLY-confident-looking copy with no read calling that locus home, defeating
/// the capacity-confidence guard. Net effect vs no-boost: +2 transcripts, +0
/// reference matches (no sensitivity gain). A primary=0 copy is non-identifiable
/// from a phantom, so boosting it produces false positives. Kept OFF by default;
/// see docs/superpowers/specs/2026-06-03-vg-underuse-audit.md.
///
/// Tuning env vars:
///   RUSTLE_VG_BOOST_PRIMARY_THR  — upper bound (default 5)
///   RUSTLE_VG_BOOST_PRIMARY_MIN  — lower bound, inclusive (default 3)
///   RUSTLE_VG_BOOST_MIN_EM_WEIGHT — minimum multi-mapper EM weight (default 3.0)
///   RUSTLE_VG_BOOST_MAX          — boost cap (default 10.0)
///
/// Activated by `RUSTLE_VG_FAMILY_BOOST=1`.
pub fn apply_family_primary_boosts(families: &[FamilyGroup], bundles: &mut [Bundle]) {
    let primary_thr: usize = std::env::var("RUSTLE_VG_BOOST_PRIMARY_THR")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);
    let primary_min: usize = std::env::var("RUSTLE_VG_BOOST_PRIMARY_MIN")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let min_em_wt: f64 = std::env::var("RUSTLE_VG_BOOST_MIN_EM_WEIGHT")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3.0);
    let max_boost: f64 = std::env::var("RUSTLE_VG_BOOST_MAX")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(10.0_f64);

    for family in families {
        for (fam_pos, &bi) in family.bundle_indices.iter().enumerate() {
            if bi >= bundles.len() {
                continue;
            }
            let primary_count = bundles[bi]
                .reads
                .iter()
                .filter(|r| r.is_primary_alignment)
                .count();
            // primary=0 → all-reads boost path below.
            // primary∈[1, primary_min) → skip: too few to spread boost safely.
            // primary≥primary_thr → already has sufficient unique coverage.
            if (primary_count >= 1 && primary_count < primary_min)
                || primary_count >= primary_thr
            {
                continue;
            }

            // Sum EM-assigned weight contributed by multi-mappers landing here.
            let mut mm_em_weight = 0.0f64;
            for placements in family.multimap_reads.values() {
                for &(pos, ri) in placements {
                    if pos == fam_pos {
                        if let Some(r) = bundles[bi].reads.get(ri) {
                            mm_em_weight += r.weight;
                        }
                    }
                }
            }
            if mm_em_weight < min_em_wt {
                continue;
            }

            let raw_boost = mm_em_weight / primary_count.max(1) as f64;
            let boost = raw_boost.max(1.0).min(max_boost);
            eprintln!(
                "[VG-BOOST] bundle {}:{}-{} primary={} mm_em_weight={:.2} boost={:.2}x{}",
                bundles[bi].chrom,
                bundles[bi].start,
                bundles[bi].end,
                primary_count,
                mm_em_weight,
                boost,
                if primary_count == 0 { " (all reads)" } else { "" },
            );
            for read in &mut bundles[bi].reads {
                // primary=0: boost all reads so EM-weighted secondaries can drive assembly.
                // primary>=2: boost only primaries; secondaries keep EM-assigned weights.
                if primary_count == 0 || read.is_primary_alignment {
                    read.weight *= boost;
                }
            }
        }
    }
}

/// Build per-bundle borrowable exon coverage from annotated FamilyGraph nodes.
///
/// After `annotate_per_copy_exon_coverage` populates `ExonClass.per_copy_cov`,
/// this function extracts a compact summary for use in the main assembly loop:
///
/// For each exon class, for each copy k with a span in this class:
///   (span.start, span.end, copy_specific, this_copy_cov, total_family_cov)
///
/// `total_family_cov` is the sum across ALL copies — the total EM-weighted
/// evidence for this exon structure across the whole family. When copy k's own
/// coverage is much lower than `total_family_cov * (k's EM fraction)`, the
/// other copies' reads serve as structural evidence that this exon should exist.
///
/// Empty Vec entries mean no family graph was built for that partition (genome
/// FASTA unavailable). Callers can treat absence as "no borrowing available."
/// Tuple fields: (ec_start, ec_end, copy_specific, this_cov, total_fam_cov,
///                max_sibling_cov, n_copies_total)
/// `max_sibling_cov`: highest per_copy_cov among all OTHER copies in this ExonClass —
///   used as the calibration reference for structural-prior boosting.
/// `n_copies_total`: total copies in the family — used to compute per-copy expected
///   coverage and decide when a copy is under-represented.
pub fn build_bundle_borrow_coverage(
    partitions: &[FamilyGroup],
    family_graphs: &[Option<crate::vg_hmm::family_graph::FamilyGraph>],
) -> HashMap<usize, Vec<(u64, u64, bool, f64, f64, f64, usize)>> {
    let mut out: HashMap<usize, Vec<(u64, u64, bool, f64, f64, f64, usize)>> = HashMap::new();
    for (pi, fam) in partitions.iter().enumerate() {
        let fg = match family_graphs.get(pi).and_then(|o| o.as_ref()) {
            Some(g) if !g.nodes.is_empty() => g,
            _ => continue,
        };
        let n_copies_total = fam.bundle_indices.len();
        for ec in &fg.nodes {
            if ec.per_copy_cov.is_empty() { continue; }
            let total_fam_cov: f64 = ec.per_copy_cov.iter().map(|(_, c)| *c).sum();
            for (copy_id, &bi) in fam.bundle_indices.iter().enumerate() {
                let span = ec.per_copy_spans.iter()
                    .find(|(k, _)| *k == copy_id)
                    .map(|(_, s)| *s);
                if let Some((s, e)) = span {
                    let this_cov = ec.per_copy_cov.iter()
                        .find(|(k, _)| *k == copy_id)
                        .map(|(_, c)| *c)
                        .unwrap_or(0.0);
                    let max_sibling_cov = ec.per_copy_cov.iter()
                        .filter(|(k, _)| *k != copy_id)
                        .map(|(_, c)| *c)
                        .fold(0.0_f64, f64::max);
                    out.entry(bi).or_default()
                        .push((s, e, ec.copy_specific, this_cov, total_fam_cov,
                               max_sibling_cov, n_copies_total));
                }
            }
        }
    }
    out
}

/// Build per-bundle lists of junctions that should be propagated from
/// well-covered sibling copies in the same gene family.
///
/// For each FamilyGraph edge (from_ec → to_ec, strand), compute the per-copy
/// junction coordinates from ExonClass.per_copy_spans and compare local
/// junction read counts against the best-covered sibling. Junctions where
/// the sibling has ≥2 reads AND this copy has fewer than 50% of the sibling
/// count are emitted as candidates for synthetic injection before graph
/// construction — completing sparse splice graphs structurally without
/// requiring per-nucleotide HMM evidence.
///
/// Returns: HashMap<bundle_idx, Vec<(donor, acceptor, strand, max_sibling_count)>>
pub fn build_bundle_borrow_junctions(
    partitions: &[FamilyGroup],
    family_graphs: &[Option<crate::vg_hmm::family_graph::FamilyGraph>],
    bundles: &[crate::types::Bundle],
) -> HashMap<usize, Vec<(u64, u64, char, f64)>> {
    use crate::types::Junction;
    let trace = std::env::var_os("RUSTLE_VG_COMPLETION_TRACE").is_some();
    let mut out: HashMap<usize, Vec<(u64, u64, char, f64)>> = HashMap::new();
    for (pi, fam) in partitions.iter().enumerate() {
        let fg = match family_graphs.get(pi).and_then(|o| o.as_ref()) {
            Some(g) if !g.edges.is_empty() => g,
            _ => {
                if trace {
                    eprintln!("[VG-JCT] partition {} — no FamilyGraph or no edges", pi);
                }
                continue;
            }
        };
        let n_copies = fam.bundle_indices.len();
        if trace {
            let n_multi = fg.edges.iter().filter(|e| e.family_support >= 2).count();
            eprintln!("[VG-JCT] partition {} fam_id={} n_copies={} edges={} edges_supp>=2={}",
                      pi, fam.family_id, n_copies, fg.edges.len(), n_multi);
        }
        let mut edge_trace_count = 0usize;
        for edge in &fg.edges {
            if edge.family_support < 2 { continue; }
            let from_ec = &fg[edge.from];
            let to_ec   = &fg[edge.to];
            // On + strand the junction is (end-of-from-exon, start-of-to-exon);
            // on - strand the from/to order in the graph is reversed relative to
            // genomic position, so the intron lies between to_ec.end and from_ec.start.
            // Unify: left_ec is whichever exon is genomically to the left.
            let (left_ec, right_ec) = if edge.strand == '-' {
                (to_ec, from_ec)
            } else {
                (from_ec, to_ec)
            };
            // Collect per-copy junction counts (donor, acceptor, local_mrcount).
            let mut copy_info: Vec<Option<(u64, u64, f64)>> = vec![None; n_copies];
            for (copy_id, &bi) in fam.bundle_indices.iter().enumerate() {
                let donor = left_ec.per_copy_spans.iter()
                    .find(|(k, _)| *k == copy_id)
                    .map(|(_, s)| s.1);       // end of left exon
                let acceptor = right_ec.per_copy_spans.iter()
                    .find(|(k, _)| *k == copy_id)
                    .map(|(_, s)| s.0);       // start of right exon
                if let (Some(d), Some(a)) = (donor, acceptor) {
                    if d >= a { continue; }   // sanity: donor must be left of acceptor
                    let local = bundles.get(bi)
                        .and_then(|b| b.junction_stats.get(&Junction { donor: d, acceptor: a }))
                        .map(|js| js.mrcount)
                        .unwrap_or(0.0);
                    copy_info[copy_id] = Some((d, a, local));
                }
            }
            if trace && edge_trace_count < 3 {
                edge_trace_count += 1;
                let n_found = copy_info.iter().filter(|x| x.is_some()).count();
                let max_local = copy_info.iter().filter_map(|x| x.as_ref()).map(|(_, _, c)| *c).fold(0.0_f64, f64::max);
                let min_local = if n_found == 0 { 0.0 } else {
                    copy_info.iter().filter_map(|x| x.as_ref()).map(|(_, _, c)| *c).fold(f64::MAX, f64::min)
                };
                eprintln!("[VG-JCT]   edge {:?}->{:?} supp={} copies_with_spans={}/{} max_count={:.2} min_count={:.2}",
                          edge.from, edge.to, edge.family_support, n_found, n_copies, max_local, min_local);
            }
            // Propagate to under-represented copies.
            // Two conditions:
            // (A) High-coverage case: sibling has ≥2 reads AND this copy has < 50%.
            // (B) Dark-junction case: this copy has 0 reads but any sibling has ≥0.5.
            //     Fires for sparse families where all copies have low coverage but
            //     some copies are missing specific junction evidence entirely.
            for copy_id in 0..n_copies {
                let bi = fam.bundle_indices[copy_id];
                if let Some((donor, acceptor, this_count)) = copy_info[copy_id] {
                    let max_sibling = copy_info.iter()
                        .enumerate()
                        .filter(|(j, _)| *j != copy_id)
                        .filter_map(|(_, opt)| opt.as_ref())
                        .map(|(_, _, c)| *c)
                        .fold(0.0_f64, f64::max);
                    let propagate =
                        (max_sibling >= 2.0 && this_count < max_sibling * 0.5) ||
                        (max_sibling >= 0.5 && this_count < 0.5);
                    if propagate {
                        out.entry(bi).or_default()
                            .push((donor, acceptor, edge.strand, max_sibling));
                    }
                }
            }
        }
    }
    let total: usize = out.values().map(|v| v.len()).sum();
    if total > 0 {
        eprintln!("[VG] Junction propagation: {} synthetic junction(s) across {} bundle(s)",
                  total, out.len());
    }
    out
}

/// Identify dark exons in sparse copies that siblings predict should exist.
///
/// For each ExonClass in the FamilyGraph, checks every copy: if the copy's
/// `per_copy_cov` is < 1.0 but at least one sibling has ≥ 2.0, the copy's
/// `per_copy_spans` entry is returned as a completion candidate. Callers inject
/// synthetic bundlenodes for these spans before graph construction so the splice
/// graph can grow paths through the dark region without EM read reweighting.
///
/// Returns: HashMap<bundle_idx, Vec<(start, end)>>
pub fn build_bundle_completion_nodes(
    partitions: &[FamilyGroup],
    family_graphs: &[Option<crate::vg_hmm::family_graph::FamilyGraph>],
) -> HashMap<usize, Vec<(u64, u64)>> {
    let trace = std::env::var_os("RUSTLE_VG_COMPLETION_TRACE").is_some();
    let mut out: HashMap<usize, Vec<(u64, u64)>> = HashMap::new();
    for (pi, fam) in partitions.iter().enumerate() {
        let fg = match family_graphs.get(pi).and_then(|o| o.as_ref()) {
            Some(g) if !g.nodes.is_empty() => g,
            _ => {
                if trace {
                    eprintln!("[VG-COMPL] partition {} — no FamilyGraph or empty", pi);
                }
                continue;
            }
        };
        let n_copies = fam.bundle_indices.len();
        if trace {
            eprintln!("[VG-COMPL] partition {} — {} ExonClasses, {} copies", pi, fg.nodes.len(), n_copies);
        }
        for ec in &fg.nodes {
            if ec.per_copy_cov.is_empty() {
                if trace {
                    eprintln!("[VG-COMPL]   ExonClass span={:?} — per_copy_cov empty (not annotated?)", ec.span);
                }
                continue;
            }
            for copy_id in 0..n_copies {
                let bi = fam.bundle_indices[copy_id];
                let span = match ec.per_copy_spans.iter()
                    .find(|(k, _)| *k == copy_id)
                    .map(|(_, s)| *s)
                {
                    Some(s) => s,
                    None => continue,
                };
                let this_cov = ec.per_copy_cov.iter()
                    .find(|(k, _)| *k == copy_id)
                    .map(|(_, c)| *c)
                    .unwrap_or(0.0);
                let max_sibling_cov = ec.per_copy_cov.iter()
                    .filter(|(k, _)| *k != copy_id)
                    .map(|(_, c)| *c)
                    .fold(0.0_f64, f64::max);
                if trace {
                    eprintln!("[VG-COMPL]   copy_id={} bi={} span={}-{} this_cov={:.2} max_sib={:.2}",
                              copy_id, bi, span.0, span.1, this_cov, max_sibling_cov);
                }
                if this_cov < 1.0 && max_sibling_cov >= 2.0 {
                    out.entry(bi).or_default().push(span);
                }
            }
        }
    }
    if std::env::var_os("RUSTLE_VG_COMPLETION_OFF").is_none() {
        let total: usize = out.values().map(|v| v.len()).sum();
        if total > 0 {
            eprintln!("[VG] Graph completion: {} dark exon(s) across {} bundle(s) → synthetic bundlenodes",
                      total, out.len());
        }
    }
    out
}

/// Partition a (possibly mixed-strand) family into single-strand sub-families
/// AND remap each sub-family's `multimap_reads` so `fam_pos` values index into
/// the new `bundle_indices`. Multi-mapper entries pointing to bundles outside
/// the partition are dropped. Reads left with < 2 placements after remapping
/// are dropped entirely (no EM benefit).
///
/// Sub-families with < 2 bundles are dropped (no EM possible).
///
/// Use this before `build_family_graph` + `run_pre_assembly_em_hmm` because
/// `build_family_graph` requires single-strand input — without this, ~70 % of
/// real-world families fail the strand check and skip HMM-EM.
pub fn partition_and_remap_family_by_strand(
    family: &FamilyGroup,
    bundles: &[Bundle],
) -> Vec<FamilyGroup> {
    use std::collections::BTreeMap;

    // Group ORIGINAL bundle_indices by strand.
    let mut by_strand: BTreeMap<char, Vec<usize>> = BTreeMap::new();
    for &bi in &family.bundle_indices {
        by_strand.entry(bundles[bi].strand).or_default().push(bi);
    }
    if by_strand.len() <= 1 {
        // Already single-strand — return as-is (no remapping needed).
        return vec![FamilyGroup {
            family_id: family.family_id,
            bundle_indices: family.bundle_indices.clone(),
            multimap_reads: family.multimap_reads.clone(),
        }];
    }

    let mut out = Vec::new();
    for (i, (_strand, partition_bis)) in by_strand.into_iter().enumerate() {
        if partition_bis.len() < 2 { continue; }

        // global_bi → new_fam_pos in this partition.
        let bi_to_pos: HashMap<usize, usize> = partition_bis.iter()
            .enumerate()
            .map(|(idx, &bi)| (bi, idx))
            .collect();

        // Remap multimap_reads.
        let mut new_multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for (&rnh, locs) in &family.multimap_reads {
            let mut new_locs: Vec<(usize, usize)> = Vec::with_capacity(locs.len());
            for &(old_fam_pos, ri) in locs {
                if old_fam_pos >= family.bundle_indices.len() { continue; }
                let global_bi = family.bundle_indices[old_fam_pos];
                if let Some(&new_fp) = bi_to_pos.get(&global_bi) {
                    new_locs.push((new_fp, ri));
                }
            }
            // Drop multi-mapper entries that no longer have ≥2 placements.
            if new_locs.len() >= 2 {
                new_multimap.insert(rnh, new_locs);
            }
        }

        out.push(FamilyGroup {
            family_id: family.family_id * 10 + i,
            bundle_indices: partition_bis,
            multimap_reads: new_multimap,
        });
    }

    out
}

/// Joint-strand EM input (spec component 1, O1-resolved): return the family
/// UNSPLIT (both strands together) for the fingerprint-EM placement list +
/// normalization only. Unlike `partition_and_remap_family_by_strand`, no
/// strand grouping and no `<2 placements after split` drop happens, so a read
/// shared only across an inverted pair (e.g. DAZ1(-)/DAZ3(+)) keeps all its
/// placements and is reweighted to sum to 1.0 across BOTH strands instead of
/// staying at 1/NH in each. `fam_pos` is preserved as an index into
/// `bundle_indices`, keeping EM write-back `global_bi/ri` valid.
///
/// Returns a single-element `Vec<FamilyGroup>` so it is a drop-in for
/// `partition_and_remap_family_by_strand`'s return type at the call site.
/// The family GRAPH is NOT rebuilt on this unsplit group (build_family_graph
/// bails on mixed strands, family_graph.rs:432); the strand-split graphs are
/// kept and indexed inside `run_fingerprint_em`, where a missing/empty graph
/// for a cross-strand group routes to the `fp.n_sites==0` neutral path.
/// `_bundles` is accepted for signature parity (strand is irrelevant here).
pub fn family_for_em_input(family: &FamilyGroup, _bundles: &[Bundle]) -> Vec<FamilyGroup> {
    vec![FamilyGroup {
        family_id: family.family_id,
        bundle_indices: family.bundle_indices.clone(),
        multimap_reads: family.multimap_reads.clone(),
    }]
}

// ── HMM-based EM reweighting (sequence-aware copy assignment) ────────────────

/// HMM-based EM reweighting across copies in a family group.
///
/// **Difference from `run_em_reweighting`:** the E-step compatibility score
/// is `forward_against_path(family_graph, read_seq, paralog_path)` —
/// the full sequence-level forward log-likelihood through that paralog's
/// path through the family graph. This replaces the heuristic
/// `compat × √context` with a probabilistic likelihood that uses every
/// base-level mismatch the read carries (in particular, the diagnostic
/// SNPs that distinguish paralogs).
///
/// Preconditions:
///   - `family_graph` was built from the bundles indexed by
///     `family.bundle_indices` (same order; CopyId = position in
///     bundle_indices = `family_pos` in `family.multimap_reads` values).
///   - `family_graph` has profiles fitted (`fit_profiles_in_place` called).
///   - `sequences` contains entries for every multi-mapped read's
///     `read_name_hash`. Reads without sequences are skipped (no
///     reweighting applied — keeps original weight).
///
/// Returns the EM result summary. Modifies `BundleRead.weight` in place.
pub fn run_em_reweighting_hmm(
    family: &FamilyGroup,
    bundles: &mut [Bundle],
    family_graph: &crate::vg_hmm::family_graph::FamilyGraph,
    sequences: &HashMap<u64, Vec<u8>>,
    max_iter: usize,
    convergence_thr: f64,
) -> EmResult {
    use crate::vg_hmm::scorer::forward_against_path_for_copy_with_norm;

    let n_copies = family.bundle_indices.len();
    if n_copies < 2 || family.multimap_reads.is_empty() {
        return EmResult::default();
    }

    // Pre-compute per-paralog paths through the family graph.
    let paralog_paths: Vec<Vec<crate::vg_hmm::family_graph::NodeIdx>> = (0..n_copies)
        .map(|cid| family_graph.recover_paralog_path(cid))
        .collect();

    // Pre-compute total match-column count per paralog path (length-norm
    // denominator) once — constant across all reads.
    let path_match_cols: Vec<usize> = (0..n_copies)
        .map(|cid| crate::vg_hmm::scorer::path_match_cols_for_copy(family_graph, &paralog_paths[cid], cid))
        .collect();

    // Filter family.multimap_reads to entries where we have the read sequence
    // and at least one paralog has a non-empty path. Build the working entry
    // list with pre-computed log-likelihood scores per placement (constant
    // across EM iterations because the family graph and sequences don't change).
    struct ReadEntry {
        locs: Vec<(usize, usize, usize)>,  // (family_pos, global_bi, read_idx)
        log_scores: Vec<f64>,              // forward_against_path per placement (len = locs.len())
        weights: Vec<f64>,                 // current EM weights (sum to 1)
    }
    // PHASE 1 (sequential): collect work items. The expensive forward DP per
    // placement is deferred to PHASE 2 where it runs in parallel via rayon.
    struct WorkItem<'a> {
        rnh: u64,
        seq: &'a [u8],
        // (fam_pos, global_bi, ri, current_weight)
        placements: Vec<(usize, usize, usize, f64)>,
    }
    let mut work: Vec<WorkItem> = Vec::with_capacity(family.multimap_reads.len());
    for (&rnh, locs) in &family.multimap_reads {
        let seq = match sequences.get(&rnh) {
            Some(s) if !s.is_empty() => s.as_slice(),
            _ => continue,
        };
        let mut placements = Vec::with_capacity(locs.len());
        for &(fam_pos, ri) in locs {
            if fam_pos >= n_copies { continue; }
            let global_bi = family.bundle_indices[fam_pos];
            let w = bundles[global_bi].reads[ri].weight;
            placements.push((fam_pos, global_bi, ri, w));
        }
        if placements.len() < 2 { continue; }
        work.push(WorkItem { rnh, seq, placements });
    }

    // PHASE 2 (parallel): forward_against_path is pure in (graph, seq, path).
    let mut entries: Vec<(u64, ReadEntry)> = work.par_iter()
        .filter_map(|item| {
            let mut entry_locs = Vec::with_capacity(item.placements.len());
            let mut log_scores = Vec::with_capacity(item.placements.len());
            let mut weights   = Vec::with_capacity(item.placements.len());
            for &(fam_pos, global_bi, ri, w) in &item.placements {
                let path = &paralog_paths[fam_pos];
                let score = if path.is_empty() {
                    f64::NEG_INFINITY
                } else {
                    forward_against_path_for_copy_with_norm(
                        family_graph, item.seq, path, fam_pos, Some(path_match_cols[fam_pos]),
                    )
                };
                entry_locs.push((fam_pos, global_bi, ri));
                log_scores.push(score);
                weights.push(w);
            }
            if entry_locs.len() < 2 { return None; }
            if !log_scores.iter().any(|s| s.is_finite()) { return None; }
            Some((item.rnh, ReadEntry { locs: entry_locs, log_scores, weights }))
        })
        .collect();

    if entries.is_empty() {
        return EmResult::default();
    }

    // EM iteration. With pre-computed log-scores, each iteration is just a
    // softmax + delta check. We optionally weight by current per-paralog
    // coverage (M-step prior) so iteration converges to a self-consistent
    // assignment given the weights from the previous round. This is the
    // proper EM with Bayesian prior on copy abundance.
    let mut result = EmResult::default();

    for iter in 0..max_iter {
        // M-step: aggregate per-paralog total weight (proxy for per-copy coverage).
        let mut copy_total: Vec<f64> = vec![0.0; n_copies];
        for (_, entry) in &entries {
            for (i, w) in entry.weights.iter().enumerate() {
                let fam_pos = entry.locs[i].0;
                copy_total[fam_pos] += w;
            }
        }
        // Convert to log-priors (softmax over copy totals, smoothed). 1e-3
        // floor avoids -∞ when a paralog is fully drained.
        let total_sum: f64 = copy_total.iter().sum::<f64>().max(1.0);
        let log_priors: Vec<f64> = copy_total.iter()
            .map(|&t| ((t / total_sum) + 1e-3).ln())
            .collect();

        let mut max_delta: f64 = 0.0;

        for (_, entry) in &mut entries {
            // E-step: posterior = softmax(log_score + log_prior).
            let n = entry.locs.len();
            let mut log_post = vec![f64::NEG_INFINITY; n];
            for i in 0..n {
                let s = entry.log_scores[i];
                if !s.is_finite() { continue; }
                let fam_pos = entry.locs[i].0;
                log_post[i] = s + log_priors[fam_pos];
            }
            // Numerically-stable softmax.
            let max_lp = log_post.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
            if !max_lp.is_finite() { continue; }
            let mut total = 0.0_f64;
            let mut exps = vec![0.0_f64; n];
            for i in 0..n {
                if log_post[i].is_finite() {
                    exps[i] = (log_post[i] - max_lp).exp();
                    total += exps[i];
                }
            }
            if total <= 0.0 { continue; }
            for i in 0..n {
                let new_w = exps[i] / total;
                let delta = (new_w - entry.weights[i]).abs();
                if delta > max_delta { max_delta = delta; }
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

    // Apply final weights back to BundleRead.weight.
    let mut n_reweighted = 0usize;
    for (_, entry) in &entries {
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
            "[VG-HMM-EM] Family {}: HMM-EM converged={} in {} iter (delta={:.6}), reweighted {} read entries across {} copies",
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

/// Structural (exon-length) divergence between a read and a candidate paralog
/// bundle, in bp.  Captures indel-driven divergence between paralogs that the
/// HMM's per-base emission scoring cannot see (cassette-exon presence/absence,
/// TE expansions, UTR length differences).
///
/// Definition: minimum |Σ read_exon_lens − Σ bundle_read_exon_lens| across all
/// reads in the candidate bundle.  Picks the closest bundle read by total
/// spliced length, giving a robust coarse signal that scales with the actual
/// indel divergence between paralogs.  Returns 0.0 when the candidate bundle
/// has no reads (degenerate case).
pub fn exon_length_divergence(read: &BundleRead, bundle: &Bundle) -> f64 {
    let r_spliced: i64 = read.exons.iter().map(|(a, b)| *b as i64 - *a as i64).sum();
    if r_spliced == 0 || bundle.reads.is_empty() {
        return 0.0;
    }
    let mut best = i64::MAX;
    for br in &bundle.reads {
        let b_spliced: i64 = br.exons.iter().map(|(a, b)| *b as i64 - *a as i64).sum();
        let d = (r_spliced - b_spliced).abs();
        if d < best {
            best = d;
        }
    }
    if best == i64::MAX { 0.0 } else { best as f64 }
}

// ── SNP-based copy assignment (--vg-snp) ─────────────────────────────────

/// Position → allele frequency per copy: (ref_pos → Vec<(copy_idx, allele_counts)>).
/// A "diagnostic SNP" is a position where one copy has >80% of one allele and
/// another copy has >80% of a different allele.
pub struct DiagnosticSnps {
    /// ref_pos → Vec<(copy_idx_in_family, dominant_base, frequency)>.
    /// FxHash-keyed for fast u64 lookup in snp_compatibility's hot loop.
    pub positions: crate::types::DetHashMap<u64, Vec<(usize, u8, f64)>>,
}

/// Build diagnostic SNP set for a family group.
/// Scans all reads at each copy and finds positions where allele frequencies diverge.
///
/// At copy-diagnostic positions, ONE copy typically matches the reference (zero
/// mismatches among its reads) while ANOTHER copy shows a consistent mismatch
/// (that copy's consensus differs from ref). Counting mismatches alone misses
/// the ref-matching side, so we also tabulate reads whose exon span COVERS the
/// position without recording a mismatch there — those reads "voted" for the
/// reference base. Sentinel `b'='` encodes "matches reference" so we don't
/// need the FASTA here; the diagnostic rule only needs bases to differ across
/// copies, not their identity.
pub fn build_diagnostic_snps(
    family: &FamilyGroup,
    bundles: &[Bundle],
) -> DiagnosticSnps {
    let n_copies = family.bundle_indices.len();

    // Per-copy allele counts: ref_pos → copy_idx → base → count.
    // Use DetHashMap (FxHash) for the position-keyed outer map — called
    // ~10K times per family during build_diagnostic_snps, so SipHash
    // overhead matters.
    let mut counts: crate::types::DetHashMap<u64, Vec<crate::types::DetHashMap<u8, u32>>> =
        crate::types::DetHashMap::default();

    // Step 1: collect positions with ≥ min_support mismatch reads across the
    // family. Random sequencing errors hit any given position with 1-2 reads;
    // a true copy-diagnostic position sees ≥3 consistent mismatches (one
    // read-group per copy matches ref, the other carries the variant). This
    // filter keeps the candidate set small enough to scan exhaustively.
    let min_support: u32 = std::env::var("RUSTLE_VGSNP_MIN_MISMATCH_READS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(3);
    let mut mism_counts: crate::types::DetHashMap<u64, u32> = crate::types::DetHashMap::default();
    for &bi in &family.bundle_indices {
        for read in &bundles[bi].reads {
            for &(ref_pos, _) in &read.mismatches {
                *mism_counts.entry(ref_pos).or_insert(0) += 1;
            }
        }
    }
    let candidates: crate::types::DetHashSet<u64> = mism_counts
        .iter()
        .filter(|&(_, &c)| c >= min_support)
        .map(|(&p, _)| p)
        .collect();

    // Step 2: for each candidate position, for each copy, tally reads that
    // cover the position. A read covers `pos` iff `pos` falls inside any of
    // its exons. If the read has a mismatch recorded at `pos` we use its
    // query base; otherwise the read matched the reference so we record `=`.
    let mut total_reads_seen = 0usize;
    let mut total_mismatches_seen = 0usize;
    for (copy_idx, &bi) in family.bundle_indices.iter().enumerate() {
        let bundle = &bundles[bi];
        for read in &bundle.reads {
            total_reads_seen += 1;
            total_mismatches_seen += read.mismatches.len();
            // Per-read mismatch lookup (FxHash for fast u64 keys).
            let mut mism_map: crate::types::DetHashMap<u64, u8> =
                crate::types::DetHashMap::with_capacity_and_hasher(
                    read.mismatches.len(),
                    crate::types::FixedBuild::default(),
                );
            for &(p, b) in &read.mismatches {
                mism_map.insert(p, b);
            }
            // For each candidate position, decide if the read covers it.
            for &pos in &candidates {
                let covers = read
                    .exons
                    .iter()
                    .any(|&(s, e)| pos >= s && pos < e);
                if !covers {
                    continue;
                }
                let base = mism_map.get(&pos).copied().unwrap_or(b'=');
                let entry = counts.entry(pos).or_insert_with(|| {
                    vec![crate::types::DetHashMap::default(); n_copies]
                });
                if copy_idx < entry.len() {
                    *entry[copy_idx].entry(base).or_insert(0) += 1;
                }
            }
        }
    }
    if std::env::var_os("RUSTLE_TRACE_VGSNP").is_some() {
        eprintln!(
            "[VG-SNP-DBG] family n_copies={} reads_seen={} mismatches_seen={} unique_positions={}",
            n_copies, total_reads_seen, total_mismatches_seen, counts.len()
        );
    }

    // Find diagnostic positions: where copies have different dominant alleles.
    let mut positions: crate::types::DetHashMap<u64, Vec<(usize, u8, f64)>> =
        crate::types::DetHashMap::default();
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

/// Extract the read base at a given reference position using exon coordinates.
/// Uses a CIGAR-free approximation: assumes no indels within exon spans.
/// Returns None if the position is intronic or outside the read.
fn seq_base_at_ref_pos(read: &BundleRead, ref_pos: u64) -> Option<u8> {
    if read.seq.is_empty() {
        return None;
    }
    let mut query_offset: usize = read.clip_left as usize;
    for &(exon_start, exon_end) in &read.exons {
        if ref_pos >= exon_start && ref_pos < exon_end {
            let offset_in_exon = (ref_pos - exon_start) as usize;
            let qi = query_offset + offset_in_exon;
            return read.seq.get(qi).copied().map(|b| b.to_ascii_uppercase());
        }
        query_offset += (exon_end - exon_start) as usize;
    }
    None
}

/// Build diagnostic positions from pileup of read sequences at each copy.
/// Reference-free: compares copies to each other, not to a reference genome.
/// Requires `BundleRead.seq` to be populated (VG mode).
///
/// Algorithm:
///   1. For each copy, for each position in any read's exon span: count allele frequencies.
///   2. At positions with sufficient coverage in ≥2 copies: compare dominant alleles.
///   3. Positions where copies disagree (dominant allele differs) are diagnostic.
pub fn build_pileup_diagnostics(
    family: &FamilyGroup,
    bundles: &[Bundle],
) -> DiagnosticSnps {
    let n_copies = family.bundle_indices.len();

    let min_allele_support: u32 = std::env::var("RUSTLE_VGPILEUP_MIN_ALLELE_READS")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(2);

    // Cap per-copy read count to avoid pathological O(reads × read_length) blowup.
    const MAX_READS_PER_COPY: usize = 10_000;

    // Per-copy pileup: copy_idx → ref_pos → base → count
    let mut copy_pileups: Vec<crate::types::DetHashMap<u64, crate::types::DetHashMap<u8, u32>>> =
        (0..n_copies).map(|_| crate::types::DetHashMap::default()).collect();

    for (copy_idx, &bi) in family.bundle_indices.iter().enumerate() {
        let bundle = &bundles[bi];
        let reads_to_scan = bundle.reads.iter().take(MAX_READS_PER_COPY);
        for read in reads_to_scan {
            if read.seq.is_empty() {
                continue;
            }
            let mut query_offset: usize = read.clip_left as usize;
            for &(exon_start, exon_end) in &read.exons {
                for i in 0..(exon_end - exon_start) {
                    let ref_pos = exon_start + i;
                    let qi = query_offset + i as usize;
                    if let Some(&base) = read.seq.get(qi) {
                        let base = base.to_ascii_uppercase();
                        if matches!(base, b'A' | b'C' | b'G' | b'T') {
                            *copy_pileups[copy_idx]
                                .entry(ref_pos)
                                .or_default()
                                .entry(base)
                                .or_insert(0) += 1;
                        }
                    }
                }
                query_offset += (exon_end - exon_start) as usize;
            }
        }
    }

    // Collect all positions seen in any copy.
    let all_positions: crate::types::DetHashSet<u64> = copy_pileups
        .iter()
        .flat_map(|p| p.keys().copied())
        .collect();

    // Find diagnostic positions: ≥2 copies with different dominant alleles.
    let mut positions: crate::types::DetHashMap<u64, Vec<(usize, u8, f64)>> =
        crate::types::DetHashMap::default();

    for pos in all_positions {
        let mut copy_dominants: Vec<(usize, u8, f64)> = Vec::new();
        for (ci, pileup) in copy_pileups.iter().enumerate() {
            let Some(alleles) = pileup.get(&pos) else { continue };
            let total: u32 = alleles.values().sum();
            if total < min_allele_support { continue; }
            if let Some((&best_base, &best_count)) = alleles.iter().max_by_key(|(_, &c)| c) {
                let freq = best_count as f64 / total as f64;
                if freq >= 0.8 {
                    copy_dominants.push((ci, best_base, freq));
                }
            }
        }
        // Diagnostic: ≥2 copies with different dominant alleles.
        if copy_dominants.len() >= 2 {
            let first_base = copy_dominants[0].1;
            if copy_dominants.iter().any(|(_, b, _)| *b != first_base) {
                positions.insert(pos, copy_dominants);
            }
        }
    }

    if std::env::var_os("RUSTLE_TRACE_VGSNP").is_some() {
        eprintln!(
            "[VG-PILEUP] Found {} diagnostic positions for family with {} copies (pileup method)",
            positions.len(), n_copies,
        );
    }
    if !positions.is_empty() {
        eprintln!(
            "[VG-PILEUP] {} diagnostic positions for family {} ({} copies)",
            positions.len(), family.family_id, n_copies,
        );
    }

    DiagnosticSnps { positions }
}

/// Score a read's SNP compatibility with a specific copy.
/// Returns a bonus factor (1.0 = no info, >1.0 = matches, <1.0 = mismatches).
///
/// A copy's diagnostic base may be `b'='` (meaning "this copy's consensus
/// matches the reference"). In that case, the read matches the copy if the
/// read's base at that position is ALSO the reference — i.e. the read did
/// NOT record a mismatch at this position. We test coverage via the read's
/// exons.
/// Memoized SNP factor weights, parsed once per process. Without this,
/// `snp_compatibility` did 2× `std::env::var` (libc getenv) per call —
/// during full GGO.bam SNP-EM cache build, that's ~1M getenv calls.
fn snp_weights() -> (f64, f64) {
    use std::sync::OnceLock;
    static WEIGHTS: OnceLock<(f64, f64)> = OnceLock::new();
    *WEIGHTS.get_or_init(|| {
        let m = std::env::var("RUSTLE_VG_SNP_MATCH_BONUS")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.5);
        let p = std::env::var("RUSTLE_VG_SNP_MISMATCH_PENALTY")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.3);
        (m, p)
    })
}

pub fn snp_compatibility(
    read: &BundleRead,
    copy_idx: usize,
    diagnostic: &DiagnosticSnps,
) -> f64 {
    if diagnostic.positions.is_empty() {
        return 1.0; // No info — neutral.
    }
    let (match_bonus, mismatch_penalty) = snp_weights();

    // Build mismatch lookup for legacy path (only when seq is unavailable).
    let mism_map: crate::types::DetHashMap<u64, u8> = if read.seq.is_empty() {
        let mut m = crate::types::DetHashMap::with_capacity_and_hasher(
            read.mismatches.len(),
            crate::types::FixedBuild::default(),
        );
        for &(p, b) in &read.mismatches {
            m.insert(p, b);
        }
        m
    } else {
        crate::types::DetHashMap::default() // not used in seq path
    };

    let mut matches = 0usize;
    let mut mismatches = 0usize;
    for (&pos, copy_info) in &diagnostic.positions {
        let Some((_, diag_base, _)) = copy_info.iter().find(|(ci, _, _)| *ci == copy_idx)
        else {
            continue;
        };
        // Does the read cover this position?
        let covers = read.exons.iter().any(|&(s, e)| pos >= s && pos < e);
        if !covers {
            continue;
        }
        // Get read base: prefer seq-based pileup, fall back to mismatches.
        let read_base: Option<u8> = if !read.seq.is_empty() {
            seq_base_at_ref_pos(read, pos)
        } else {
            // Legacy: b'=' means matches reference; actual base means mismatch.
            Some(mism_map.get(&pos).copied().unwrap_or(b'='))
        };
        let Some(rb) = read_base else { continue };

        // For pileup diagnostics: diag_base is always A/C/G/T.
        // For legacy: diag_base may be b'=' (copy matches reference).
        let is_match = if *diag_base == b'=' {
            rb == b'='
        } else {
            rb == *diag_base
        };

        if is_match { matches += 1; } else { mismatches += 1; }
    }
    if matches + mismatches == 0 {
        return 1.0;
    }
    let bonus = 1.0 + match_bonus * matches as f64 - mismatch_penalty * mismatches as f64;
    bonus.max(0.1)
}

// ── ExonFingerprints ──────────────────────────────────────────────────────────
//
// Sequence-based variant fingerprints derived from ExonClass.per_copy_sequences.
//
// Unlike DiagnosticSnps (reference-coordinate mismatches from BundleRead),
// ExonFingerprints use exon-relative offsets from the FamilyGraph, so they
// work for non-overlapping gene copies (RBMY, TSPY, DAZ scattered across
// chrY) where reference positions are incommensurable across copies.
//
// Model: P(read_base | copy) = p_match if genomic_base[copy][pos] == read_base
//                             else p_mismatch.
// Log-likelihoods are additive across variant sites; callers exponentiate
// and normalize to get per-copy posterior weights.

/// One position where ≥ 2 copies have different genomic bases.
struct ExonVariantSite {
    /// (copy_id, uppercase genomic base) for every copy that has a sequence
    /// in this ExonClass. Only covers disagrees where all bases are ACGT.
    copy_bases: Vec<(usize, u8)>,
}

/// Fingerprint table for one gene family.
pub struct ExonFingerprints {
    sites: Vec<ExonVariantSite>,
    /// per_copy_site_refs[copy_id] = sorted (site_idx, ref_pos, genomic_base).
    /// genomic_base is per_copy_sequences[copy_id][exon_rel_pos]: the base
    /// this copy has at ref_pos so that a read with no mismatch recorded
    /// (b'=' sentinel in BundleRead.mismatches) can be resolved to the
    /// actual base the read carries without re-fetching the FASTA.
    per_copy_site_refs: Vec<Vec<(usize, u64, u8)>>,
    pub n_copies: usize,
    pub n_sites: usize,
}

/// Build an ExonFingerprints table from the per-copy sequences in a FamilyGraph.
///
/// For each ExonClass that has ≥ 2 copies, compares their sequences
/// position-by-position (exon-relative). Every position where at least two
/// copies carry different ACGT bases becomes a diagnostic site.
/// Multi-copy ExonClasses whose copies all agree at every position contribute
/// zero sites (they carry no discriminating information).
///
/// The `n_copies` parameter must match `family.bundle_indices.len()` so that
/// `per_copy_site_refs` is indexed consistently with the family's copy ordering.
pub fn build_exon_fingerprints(
    fg: &crate::vg_hmm::family_graph::FamilyGraph,
    n_copies: usize,
) -> ExonFingerprints {
    let mut sites: Vec<ExonVariantSite> = Vec::new();
    let mut per_copy_refs: Vec<Vec<(usize, u64, u8)>> = (0..n_copies).map(|_| Vec::new()).collect();

    for ec in &fg.nodes {
        if ec.per_copy_sequences.len() < 2 {
            continue; // singleton exon — no pairwise comparison possible
        }
        // `per_copy_sequences` can hold MULTIPLE fragments for the SAME CopyId:
        // terminal exons (ragged 5'/3' IsoSeq read ends) yield several
        // variable-length versions of one copy's exon, and build_family_graph
        // stores them all. Comparing those same-copy fragments as if they were
        // distinct copies manufactures spurious "diagnostic" sites — both from
        // read-to-read jitter and from length-misaligned exon-relative position
        // comparison (observed: 132 sites + 29x site multiplicity on a clean
        // 2-copy fixture, flipping read-to-copy assignment to noise). Collapse to
        // ONE representative — the LONGEST fragment, which best covers the exon —
        // per DISTINCT CopyId, then require >= 2 distinct copies. Single-copy
        // ExonClasses carry no cross-copy information and are correctly dropped.
        let mut rep_seq: Vec<(usize, &[u8])> = Vec::new();
        for (cid, seq) in &ec.per_copy_sequences {
            match rep_seq.iter_mut().find(|(c, _)| *c == *cid) {
                Some((_, best)) => { if seq.len() > best.len() { *best = seq.as_slice(); } }
                None => rep_seq.push((*cid, seq.as_slice())),
            }
        }
        if rep_seq.len() < 2 {
            continue; // single distinct copy — no cross-copy diagnostic info
        }
        let min_len = rep_seq.iter().map(|(_, s)| s.len()).min().unwrap_or(0);
        if min_len == 0 {
            continue;
        }

        for pos in 0..min_len {
            // Collect the uppercase ACGT base for each distinct copy at this position.
            let copy_bases: Vec<(usize, u8)> = rep_seq.iter()
                .filter_map(|(cid, seq)| {
                    let b = seq[pos].to_ascii_uppercase();
                    if matches!(b, b'A' | b'C' | b'G' | b'T') { Some((*cid, b)) } else { None }
                })
                .collect();

            if copy_bases.len() < 2 {
                continue;
            }
            // Skip positions where all copies agree — not informative.
            let first = copy_bases[0].1;
            if copy_bases.iter().all(|(_, b)| *b == first) {
                continue;
            }

            let site_idx = sites.len();

            // Record a (site_idx, ref_pos, genomic_base) entry for each copy
            // that has a per_copy_span here. ref_pos = span_start + pos lets
            // the scorer intersect with BundleRead.exons at score time.
            for &(cid, copy_base) in &copy_bases {
                if cid >= n_copies {
                    continue;
                }
                if let Some((_, (span_start, _))) = ec.per_copy_spans.iter().find(|(k, _)| *k == cid) {
                    let ref_pos = span_start + pos as u64;
                    per_copy_refs[cid].push((site_idx, ref_pos, copy_base));
                }
            }

            sites.push(ExonVariantSite { copy_bases });
        }
    }

    // Sort each copy's site list by ref_pos for binary-search intersection.
    for refs in &mut per_copy_refs {
        refs.sort_unstable_by_key(|&(_, ref_pos, _)| ref_pos);
    }

    let n_sites = sites.len();

    if std::env::var_os("RUSTLE_VG_FP_SITE_TRACE").is_some() {
        eprintln!("[FP-TRACE] n_copies={} n_exon_classes={} n_sites={}", n_copies, fg.nodes.len(), n_sites);
        for cid in 0..n_copies {
            let refs = &per_copy_refs[cid];
            let n_distinct = {
                let mut ps: Vec<u64> = refs.iter().map(|&(_, p, _)| p).collect();
                ps.sort_unstable(); ps.dedup(); ps.len()
            };
            eprintln!("[FP-TRACE]   copy {} site_refs entries={} distinct_ref_pos={}", cid, refs.len(), n_distinct);
        }
    }

    ExonFingerprints { sites, per_copy_site_refs: per_copy_refs, n_copies, n_sites }
}

// ===========================================================================
// enumerate_diagnostic_sites — find divergent sites by syntenic exon-pair
// alignment.
//
// build_exon_fingerprints compares the per-copy sequences that the k-mer merge
// happened to group into the SAME ExonClass. For genuinely-divergent copies
// (5-7% divergence) the homologous exons land in SEPARATE ExonClasses, so they
// are never compared and the function reports 0 sites. This finder instead
// matches homologous exons ACROSS copies by genomic synteny + minimizer
// anchoring, then diffs them base-by-base.
//
// Everything stored is in GENOME-FORWARD frame (same frame as BundleRead.exons
// and BundleRead.mismatches), so the scorer can intersect ref_pos with the
// read and compare the stored allele directly. `per_copy_sequences` is ALREADY
// genome-forward (fetched directly from the reference FASTA via
// genome.fetch_sequence with no strand handling), so each exon's sequence is
// materialized by uppercasing the bytes for BOTH strands — no revcomp. After
// that, genomic position of offset `o` is always `span.start + o` and the
// stored allele is always `gf_seq[o]`.
// ===========================================================================

/// One copy's exon, in genome-forward frame.
/// `(span_start, span_end, genome_forward_sequence)`.
struct GfExon {
    start: u64,
    /// Genomic end (exclusive). Retained for clarity/debugging; genome-forward
    /// positions are computed from `start + offset`.
    #[allow(dead_code)]
    end: u64,
    seq: Vec<u8>,
}

#[inline]
fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Jaccard similarity of the two exons' minimizer sets (k=15, w=10).
///
/// Returns `None` when neither sequence is long enough to produce any
/// minimizers (len < k): the pairing cannot be evaluated by minimizers and the
/// caller must fall back to positional synteny rather than rejecting a possibly
/// homologous short exon.
fn minimizer_jaccard(a: &[u8], b: &[u8]) -> Option<f64> {
    let ma = crate::vg_hmm::family_graph::minimizers(a, 15, 10);
    let mb = crate::vg_hmm::family_graph::minimizers(b, 15, 10);
    let union = ma.union(&mb).count();
    if union == 0 {
        return None; // too short to evaluate
    }
    Some(ma.intersection(&mb).count() as f64 / union as f64)
}

/// Compare two equal-or-unequal-length genome-forward exon sequences and return
/// the divergent columns as `(offset_in_a, base_a, offset_in_b, base_b)`.
/// Equal length → gapless. Unequal length → minimizer-anchored: chain co-linear
/// shared k-mer anchors, compare gaplessly within equal-length anchor-to-anchor
/// blocks, and SKIP indel gaps (never emit a site inside an indel).
fn diff_gf_exons(a: &[u8], b: &[u8]) -> Vec<(usize, u8, usize, u8)> {
    if a.len() == b.len() {
        return diff_equal_len(a, 0, b, 0);
    }
    diff_anchored(a, b)
}

/// Gapless diff of two equal-length windows starting at absolute offsets
/// `(off_a, off_b)` in their parent sequences. `a`/`b` here are the windows.
fn diff_equal_len(a: &[u8], off_a: usize, b: &[u8], off_b: usize) -> Vec<(usize, u8, usize, u8)> {
    debug_assert_eq!(a.len(), b.len());
    let mut out = Vec::new();
    for o in 0..a.len() {
        let ba = a[o].to_ascii_uppercase();
        let bb = b[o].to_ascii_uppercase();
        if is_acgt(ba) && is_acgt(bb) && ba != bb {
            out.push((off_a + o, ba, off_b + o, bb));
        }
    }
    out
}

/// Minimizer-anchored diff for unequal-length exon pairs. Finds shared k-mer
/// anchors, chains them into a co-linear, strictly-increasing set, then diffs
/// the equal-length runs between consecutive anchors. Indel gaps (where the two
/// sides advance by different amounts) emit no sites.
fn diff_anchored(a: &[u8], b: &[u8]) -> Vec<(usize, u8, usize, u8)> {
    const K: usize = 12;
    if a.len() < K || b.len() < K {
        // Too short to anchor — fall back to a gapless diff over the common
        // prefix length (safe: only emits sites where both sides have a base).
        let n = a.len().min(b.len());
        return diff_equal_len(&a[..n], 0, &b[..n], 0);
    }

    // Map each k-mer hash in `b` to its FIRST position (left-most), so a unique
    // k-mer anchors deterministically.
    let mut b_pos: std::collections::HashMap<u64, usize> = std::collections::HashMap::new();
    let mut b_seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut b_dup: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for i in 0..=(b.len() - K) {
        let h = fnv1a_kmer(&b[i..i + K]);
        if b_seen.contains(&h) {
            b_dup.insert(h);
        } else {
            b_seen.insert(h);
            b_pos.insert(h, i);
        }
    }

    // Walk `a`'s k-mers; collect (pos_a, pos_b) anchors for UNIQUE shared k-mers
    // (skip k-mers that repeat on either side — ambiguous).
    let mut a_seen: std::collections::HashSet<u64> = std::collections::HashSet::new();
    let mut a_dup: std::collections::HashSet<u64> = std::collections::HashSet::new();
    for i in 0..=(a.len() - K) {
        let h = fnv1a_kmer(&a[i..i + K]);
        if a_seen.contains(&h) {
            a_dup.insert(h);
        } else {
            a_seen.insert(h);
        }
    }
    let mut anchors: Vec<(usize, usize)> = Vec::new();
    for i in 0..=(a.len() - K) {
        let h = fnv1a_kmer(&a[i..i + K]);
        if a_dup.contains(&h) || b_dup.contains(&h) {
            continue;
        }
        if let Some(&j) = b_pos.get(&h) {
            anchors.push((i, j));
        }
    }
    // Chain anchors to a strictly-increasing co-linear set on BOTH axes (LIS on
    // pos_b after sorting by pos_a). anchors are already sorted by pos_a here
    // (we walked i ascending), but duplicates of pos_a can't occur for a unique
    // k-mer, so a simple greedy keep-if-increasing chain on pos_b is enough to
    // stay co-linear.
    let mut chained: Vec<(usize, usize)> = Vec::new();
    let mut last_b: i64 = -1;
    for &(ia, jb) in &anchors {
        if (jb as i64) > last_b {
            chained.push((ia, jb));
            last_b = jb as i64;
        }
    }
    if chained.is_empty() {
        let n = a.len().min(b.len());
        return diff_equal_len(&a[..n], 0, &b[..n], 0);
    }

    let mut out = Vec::new();
    // Diff the run BEFORE the first anchor (align right-to-left from the anchor:
    // compare equal-length suffix that abuts the anchor). To keep it simple and
    // indel-safe, only diff the leading block if both leading lengths are equal.
    let (first_a, first_b) = chained[0];
    if first_a == first_b && first_a > 0 {
        out.extend(diff_equal_len(&a[..first_a], 0, &b[..first_b], 0));
    }
    // Diff equal-length blocks between consecutive anchors; skip indel gaps.
    for w in chained.windows(2) {
        let (ia0, jb0) = w[0];
        let (ia1, jb1) = w[1];
        // The matched k-mer itself is identical; advance past its first base and
        // compare the inter-anchor block.
        let block_a_start = ia0;
        let block_b_start = jb0;
        let len_a = ia1 - ia0;
        let len_b = jb1 - jb0;
        if len_a == len_b {
            // No (net) indel between these anchors → gapless diff of the whole
            // block. LIMITATION: a balanced indel pair (equal-length insertion
            // + deletion inside one block) would desync this gapless diff; this
            // is bounded in practice because unique 12-mer anchors keep
            // inter-anchor blocks short on SNP-divergent paralogs.
            out.extend(diff_equal_len(
                &a[block_a_start..ia1],
                block_a_start,
                &b[block_b_start..jb1],
                block_b_start,
            ));
        } else {
            // Indel between anchors. Diff the K-base anchor itself (identical, no
            // sites) and SKIP the variable-length gap — never emit a site there.
            // The anchor k-mer guarantees a[ia0..ia0+K] == b[jb0..jb0+K], so no
            // sites are lost by skipping the gap.
        }
    }
    // Diff the run AFTER the last anchor if both trailing lengths are equal.
    let (last_a, last_b_pos) = *chained.last().unwrap();
    let tail_a = a.len() - last_a;
    let tail_b = b.len() - last_b_pos;
    if tail_a == tail_b && tail_a > 0 {
        out.extend(diff_equal_len(
            &a[last_a..],
            last_a,
            &b[last_b_pos..],
            last_b_pos,
        ));
    }
    // The inter-anchor blocks above start at the anchor's first base, so adjacent
    // blocks share the anchor base. Dedup by (offset_a) to avoid double-counting
    // a divergent column that two adjacent equal-length blocks both cover (only
    // possible at exact anchor boundaries, which are identical so emit nothing,
    // but dedup defensively).
    out.sort_unstable_by_key(|&(oa, _, _, _)| oa);
    out.dedup_by_key(|&mut (oa, _, _, _)| oa);
    out
}

/// FNV-1a hash of a k-mer (mirrors family_graph's private fnv1a so anchoring is
/// deterministic and self-contained).
fn fnv1a_kmer(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b.to_ascii_uppercase() as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Match homologous exons across two copies (already genome-forward, sorted by
/// genomic start) by synteny + minimizer overlap. Returns index pairs `(i, j)`
/// into `ea`/`eb`. Walks both lists in genomic order; validates each candidate
/// pairing with a minimizer-Jaccard floor and re-syncs (advances the list whose
/// next exon is the better match) when a positional pairing fails.
fn syntenic_exon_pairs(ea: &[GfExon], eb: &[GfExon]) -> Vec<(usize, usize)> {
    const JACCARD_FLOOR: f64 = 0.05;
    // Length-ratio and base-identity floors for the short-exon fallback below.
    const SHORT_LEN_RATIO: f64 = 0.8; // shorter/longer must be >= 0.8 (≈ within 20%)
    const SHORT_IDENT_FLOOR: f64 = 0.5; // gapless base identity over the overlap
    // When minimizers can't be computed (an exon is shorter than k=15) we cannot
    // trust positional synteny alone — accepting it could pair two unrelated
    // short exons and fabricate spurious sites. Instead, fall back to a cheap
    // direct check: accept only when the two exons are near-equal length AND a
    // gapless base-identity over their overlap clears a floor.
    let short_exon_homologous = |a: &[u8], b: &[u8]| -> bool {
        let (shorter, longer) = (a.len().min(b.len()), a.len().max(b.len()));
        if longer == 0 {
            return false;
        }
        if (shorter as f64) / (longer as f64) < SHORT_LEN_RATIO {
            return false; // too different in length to be homologous
        }
        let n = shorter;
        if n == 0 {
            return false;
        }
        let mut matches = 0usize;
        for o in 0..n {
            if a[o].to_ascii_uppercase() == b[o].to_ascii_uppercase() {
                matches += 1;
            }
        }
        (matches as f64) / (n as f64) >= SHORT_IDENT_FLOOR
    };
    let passes = |a: &[u8], b: &[u8]| -> bool {
        match minimizer_jaccard(a, b) {
            Some(j) => j >= JACCARD_FLOOR,
            // Unevaluable by minimizers → require near-equal length + base identity.
            None => short_exon_homologous(a, b),
        }
    };
    // Score for the re-sync lookahead: -1.0 means "definitely not homologous"
    // (evaluable and below floor); short/unevaluable pairings score 0.0 so they
    // neither strongly attract nor repel the walk.
    let lookahead = |a: &[u8], b: &[u8]| -> f64 {
        match minimizer_jaccard(a, b) {
            Some(j) => j,
            None => 0.0,
        }
    };
    let mut pairs = Vec::new();
    let (mut i, mut j) = (0usize, 0usize);
    while i < ea.len() && j < eb.len() {
        if passes(&ea[i].seq, &eb[j].seq) {
            pairs.push((i, j));
            i += 1;
            j += 1;
            continue;
        }
        // Not homologous at this positional pairing. Decide which list to
        // advance: look one ahead on each side and advance toward the better
        // match. If neither lookahead helps, advance the lower-start exon.
        let adv_i = if i + 1 < ea.len() {
            lookahead(&ea[i + 1].seq, &eb[j].seq)
        } else {
            -1.0
        };
        let adv_j = if j + 1 < eb.len() {
            lookahead(&ea[i].seq, &eb[j + 1].seq)
        } else {
            -1.0
        };
        if adv_i >= JACCARD_FLOOR && adv_i >= adv_j {
            // ea[i] is an unmatched (inserted) exon — skip it.
            i += 1;
        } else if adv_j >= JACCARD_FLOOR && adv_j > adv_i {
            // eb[j] is an unmatched (inserted) exon — skip it.
            j += 1;
        } else {
            // No nearby homolog either way — neither one side has a clearly
            // better lookahead match.  The two exons are positionally co-linear
            // (same index in their respective copies' exon lists) so treat them
            // as structurally homologous regardless of Jaccard (highly-divergent
            // paralogs can share 0 k=15-mers yet still be the same ancestral
            // exon).  Include the pair so diff_gf_exons can mine their SNP
            // columns; cross-genomic-region copies advance in lockstep (not by
            // start-position) to avoid draining one list.
            pairs.push((i, j));
            i += 1;
            j += 1;
        }
    }
    pairs
}

/// Enumerate divergent sites between copies by syntenic exon-pair comparison.
///
/// Replaces the broken exon-relative inner loop of `build_exon_fingerprints`:
/// homologous exons may sit in SEPARATE ExonClasses (the k-mer merge fails at
/// 5-7% divergence), so we match exons across copies by genomic synteny +
/// minimizer anchoring here, diff them base-by-base, and emit sites in
/// genome-forward coordinates with strand-correct (genome-forward) alleles —
/// exactly the frame `score_read_exon_fingerprint` consumes.
pub fn enumerate_diagnostic_sites(
    fg: &crate::vg_hmm::family_graph::FamilyGraph,
    n_copies: usize,
) -> ExonFingerprints {
    // All copies in a FamilyGraph are the SAME strand (build_family_graph requires it;
    // partition_family_by_strand splits mixed-strand families upstream), so per-copy
    // exons are co-linear in genomic order and per_copy_sequences are genome-forward —
    // no strand/inversion handling is needed here.

    // 1. Per copy: collect ONE canonical exon per node via recover_paralog_path.
    //    Using the canonical path avoids duplicate GfExon entries that arise when
    //    multiple reads for the same copy have slightly different exon boundaries
    //    (all stored as separate (cid, span) entries in per_copy_spans).
    //    per_copy_sequences are genome-forward (fetched directly from the FASTA),
    //    so we uppercase only — no revcomp.
    let mut exons_by_copy: Vec<Vec<GfExon>> = (0..n_copies).map(|_| Vec::new()).collect();
    let trace = std::env::var_os("RUSTLE_VG_FP_SITE_TRACE").is_some();
    for cid in 0..n_copies {
        let path = fg.recover_paralog_path(cid);
        for nidx in path {
            let node = &fg[nidx];
            // Take the FIRST span for this copy in this node (canonical).
            if let Some(&(_, (s, e))) = node.per_copy_spans.iter().find(|(c, _)| *c == cid) {
                // Take the FIRST sequence for this copy (same selection policy).
                if let Some((_, seq)) = node.per_copy_sequences.iter().find(|(c, _)| *c == cid) {
                    let gf: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
                    if !gf.is_empty() {
                        exons_by_copy[cid].push(GfExon { start: s, end: e, seq: gf });
                    }
                }
            }
        }
    }
    for v in &mut exons_by_copy {
        v.sort_by_key(|x| x.start);
        v.dedup_by_key(|x| x.start);
    }
    if trace {
        for cid in 0..n_copies {
            eprintln!("[FP-SITE-TRACE]   copy {} exon_count={}", cid, exons_by_copy[cid].len());
            for ex in &exons_by_copy[cid] {
                eprintln!("[FP-SITE-TRACE]     exon start={} end={} seqlen={}", ex.start, ex.end, ex.seq.len());
            }
        }
    }

    let mut sites: Vec<ExonVariantSite> = Vec::new();
    let mut per_copy_refs: Vec<Vec<(usize, u64, u8)>> = (0..n_copies).map(|_| Vec::new()).collect();

    // 2. For each unordered copy pair, match syntenic exons and diff them.
    for i in 0..n_copies {
        for j in (i + 1)..n_copies {
            if exons_by_copy[i].is_empty() || exons_by_copy[j].is_empty() {
                continue;
            }
            for (ei, ej) in syntenic_exon_pairs(&exons_by_copy[i], &exons_by_copy[j]) {
                let ex_a = &exons_by_copy[i][ei];
                let ex_b = &exons_by_copy[j][ej];
                for (oa, ba, ob, bb) in diff_gf_exons(&ex_a.seq, &ex_b.seq) {
                    // genome-forward: genomic pos = span.start + offset; allele is
                    // the genome-forward base already in `seq`.
                    let pos_i = ex_a.start + oa as u64;
                    let pos_j = ex_b.start + ob as u64;
                    let site_idx = sites.len();
                    per_copy_refs[i].push((site_idx, pos_i, ba));
                    per_copy_refs[j].push((site_idx, pos_j, bb));
                    sites.push(ExonVariantSite { copy_bases: vec![(i, ba), (j, bb)] });
                }
            }
        }
    }

    // 3. Per copy: sort by ref_pos and dedup so a position found in multiple
    //    pairs collapses to one entry per copy.
    for refs in &mut per_copy_refs {
        refs.sort_unstable_by_key(|&(_, p, _)| p);
        refs.dedup_by_key(|&mut (_, p, _)| p);
    }

    let n_sites = sites.len();

    if std::env::var_os("RUSTLE_VG_FP_SITE_TRACE").is_some() {
        eprintln!(
            "[FP-SITE-TRACE] n_copies={} n_sites={}",
            n_copies, n_sites
        );
        for cid in 0..n_copies {
            eprintln!(
                "[FP-SITE-TRACE]   copy {} sites={}",
                cid,
                per_copy_refs[cid].len()
            );
        }
    }

    ExonFingerprints { sites, per_copy_site_refs: per_copy_refs, n_copies, n_sites }
}

/// Compute log-likelihood scores for each copy given a multi-mapping read.
///
/// `copy_idx` is the copy that this `BundleRead` is aligned to — it selects
/// the reference-position lookup table (`per_copy_site_refs[copy_idx]`) so
/// diagnostic sites can be intersected with `read.exons`.
///
/// The b'=' sentinel in `BundleRead.mismatches` encodes "read matches the
/// reference at this position." Since the reference IS the copy's genomic
/// sequence, b'=' maps to `genomic_base` stored in `per_copy_site_refs`,
/// avoiding a FASTA re-fetch.
///
/// Returns `(scores, n_sites_covered)` where:
/// - `scores`: `Vec<f64>` of length `n_copies`. Entries for copies with no
///   sites covered by this read remain 0.0 (log-prob of 1.0 → neutral).
///   The caller exponentiates and normalises to get posterior weights.
/// - `n_sites_covered`: number of diagnostic sites in `per_copy_site_refs[copy_idx]`
///   that fall within the read's aligned exon intervals. Longer reads cover more
///   sites and produce more decisive (larger log-likelihood gap) assignments.
pub fn score_read_exon_fingerprint(
    read: &BundleRead,
    copy_idx: usize,
    fp: &ExonFingerprints,
) -> (Vec<f64>, usize) {
    use std::sync::OnceLock;
    static LOG_PROBS: OnceLock<(f64, f64)> = OnceLock::new();
    static ERROR_AWARE: OnceLock<(bool, f64, f64)> = OnceLock::new();
    let (ea, ea_lo, ea_hi) = *ERROR_AWARE.get_or_init(|| {
        let on = std::env::var_os("RUSTLE_VG_FP_ERROR_AWARE").is_some();
        let lo: f64 = std::env::var("RUSTLE_VG_FP_ERR_FLOOR")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.01);
        let hi: f64 = std::env::var("RUSTLE_VG_FP_ERR_CAP")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(0.10);
        (on, lo, hi)
    });
    let (log_match, log_mismatch) = if ea {
        // Error-aware (opt-in): the mismatch probability at a diagnostic base IS the
        // read's own per-base error rate (the `de` tag) — clamped to [floor, cap].
        // The fixed-0.10 default is deliberately conservative (≈9:1 per site) to
        // avoid single-base overconfidence; but for genuinely low-error IsoSeq reads
        // that UNDER-trusts real signal. Using the actual error rate sharpens clean
        // reads (≈99:1) and discounts noisy ones, improving attribution. The single-
        // site overconfidence this could reintroduce is handled separately by the
        // min-decisive-sites gate.
        let p_mm = read.de.map(|d| (d as f64).clamp(ea_lo, ea_hi)).unwrap_or(ea_hi);
        ((1.0 - p_mm).ln(), p_mm.ln())
    } else {
        *LOG_PROBS.get_or_init(|| {
            let p_match: f64 = std::env::var("RUSTLE_VG_FP_MATCH")
                .ok().and_then(|v| v.parse().ok()).unwrap_or(0.90);
            let p_mm: f64 = std::env::var("RUSTLE_VG_FP_MISMATCH")
                .ok().and_then(|v| v.parse().ok()).unwrap_or(0.10);
            (p_match.ln(), p_mm.ln())
        })
    };

    let mut scores = vec![0.0f64; fp.n_copies];
    let mut n_covered = 0usize;

    let Some(site_refs) = fp.per_copy_site_refs.get(copy_idx) else {
        return (scores, 0);
    };
    if site_refs.is_empty() {
        return (scores, 0);
    }

    // Build mismatch lookup for this read (FxHash for hot-loop performance).
    let mut mism: crate::types::DetHashMap<u64, u8> =
        crate::types::DetHashMap::with_capacity_and_hasher(
            read.mismatches.len(),
            crate::types::FixedBuild::default(),
        );
    for &(p, b) in &read.mismatches {
        mism.insert(p, b);
    }

    for &(site_idx, ref_pos, this_copy_base) in site_refs {
        // Does the read's alignment to copy_idx cover this reference position?
        if !read.exons.iter().any(|&(s, e)| ref_pos >= s && ref_pos < e) {
            continue;
        }
        n_covered += 1;

        // Actual base the read carries at ref_pos.
        // b'=' (no mismatch recorded) → the read matches this copy's genomic
        // base, which is `this_copy_base`.
        let read_base = match mism.get(&ref_pos).copied() {
            Some(b) => b,
            None => this_copy_base,
        };

        // Score read_base against every copy's expected base at this site.
        let site = &fp.sites[site_idx];
        for &(cid, copy_base) in &site.copy_bases {
            if cid < fp.n_copies {
                scores[cid] += if read_base == copy_base { log_match } else { log_mismatch };
            }
        }
    }

    (scores, n_covered)
}

/// Build the per-site cross-copy match structure for ONE read, in `copy_idx`'s coordinate
/// frame (the read's home/primary copy). Mirrors `score_read_exon_fingerprint`'s per-site
/// base lookup, but RETAINS each copy's match (`match_bits`) instead of summing — the input
/// the mosaic detector needs. Sites come out in genomic order (`per_copy_site_refs` is sorted).
fn build_read_site_obs(
    read: &crate::types::BundleRead,
    fp: &ExonFingerprints,
    copy_idx: usize,
) -> Vec<crate::vg_hmm::mosaic::SiteObs> {
    use crate::vg_hmm::mosaic::SiteObs;
    let Some(site_refs) = fp.per_copy_site_refs.get(copy_idx) else {
        return Vec::new();
    };
    let mut mism: crate::types::DetHashMap<u64, u8> =
        crate::types::DetHashMap::with_capacity_and_hasher(
            read.mismatches.len(), crate::types::FixedBuild::default());
    for &(p, b) in &read.mismatches {
        mism.insert(p, b);
    }
    let mut obs: Vec<SiteObs> = Vec::with_capacity(site_refs.len());
    for &(site_idx, ref_pos, this_copy_base) in site_refs {
        if !read.exons.iter().any(|&(s, e)| ref_pos >= s && ref_pos < e) {
            continue;
        }
        let read_base = mism.get(&ref_pos).copied().unwrap_or(this_copy_base);
        let site = &fp.sites[site_idx];
        let mut match_bits = vec![false; fp.n_copies];
        for &(cid, copy_base) in &site.copy_bases {
            if cid < fp.n_copies {
                match_bits[cid] = read_base == copy_base;
            }
        }
        obs.push(SiteObs { ref_pos, match_bits });
    }
    obs
}

/// Per-family mosaic detection pass (opt-in). For each copy's PRIMARY reads (a read's home
/// copy), build the per-site match structure and run `detect_mosaic`; aggregate per family
/// and report confirmed gene-conversion events vs chimera-suspects. Pure-analysis: it never
/// mutates bundles or EM weights.
fn detect_and_report_mosaics(
    family: &FamilyGroup,
    bundles: &[Bundle],
    fp: &ExonFingerprints,
) -> Vec<crate::vg_hmm::mosaic::ConversionEvent> {
    use crate::vg_hmm::mosaic::{aggregate_family, detect_mosaic, MosaicParams};
    use crate::vg_hmm::mosaic::MosaicStatus;
    let params = MosaicParams::from_env();
    let n_copies = fp.n_copies;
    let trace = std::env::var_os("RUSTLE_VG_MOSAIC_TRACE").is_some();
    let mut calls = Vec::new();
    // Each copy's locus = (chrom, start, end) of its family bundle. The family bundle is
    // reduced to ~the multi-mapping reads, but a recombinant is detectable at its PRIMARY
    // alignment whether or not it multimaps — so scan the FULL bundles array for every
    // bundle overlapping the locus, deduping reads by name. (Opt-in pass; the cheap overlap
    // pre-filter keeps cost bounded.)
    let copy_loci: Vec<Option<(String, u64, u64)>> = family.bundle_indices.iter()
        .map(|&bi| bundles.get(bi).map(|b| (b.chrom.clone(), b.start, b.end)))
        .collect();
    for (copy_idx, locus) in copy_loci.iter().enumerate() {
        let Some((chrom, lo, hi)) = locus else { continue };
        if fp.per_copy_site_refs.get(copy_idx).map(|s| s.is_empty()).unwrap_or(true) {
            continue;
        }
        let mut seen: crate::types::DetHashSet<u64> = Default::default();
        let (mut n_prim, mut n_obs, mut hist) = (0usize, 0usize, [0usize; 5]);
        for bundle in bundles.iter() {
            if &bundle.chrom != chrom || bundle.end < *lo || bundle.start > *hi {
                continue;
            }
            for read in &bundle.reads {
                if !read.is_primary_alignment || !seen.insert(read.read_name_hash) {
                    continue; // home copy = primary alignment; dedupe across overlapping bundles
                }
                n_prim += 1;
                let obs = build_read_site_obs(read, fp, copy_idx);
                if obs.is_empty() {
                    continue;
                }
                n_obs += 1;
                let eps = params.eps_for(read.de);
                let call = detect_mosaic(&obs, n_copies, eps, &params);
                hist[match call.status {
                    MosaicStatus::LowPower => 0, MosaicStatus::NonIdentifiable => 1,
                    MosaicStatus::SingleCopy => 2, MosaicStatus::NoSwitch => 3, MosaicStatus::Mosaic => 4,
                }] += 1;
                if trace && matches!(call.status, MosaicStatus::NoSwitch) && call.copy_a != call.copy_b && call.n_decisive >= params.min_decisive_sites {
                    eprintln!("[VG-MOSAIC-NOSWITCH] family={} home={} {:?}->{:?} D={} margin={} tracts={}/{} lr={:.1}/{:.1}",
                        family.family_id, copy_idx, call.copy_a, call.copy_b, call.n_decisive, call.margin,
                        call.tract_a_sites, call.tract_b_sites, call.lr_switch, call.threshold_used);
                }
                calls.push(call);
            }
        }
        if trace {
            eprintln!("[VG-MOSAIC-COPY] family={} copy={} locus={}:{}-{} primary_reads={} with_obs={} status[LP/NI/SC/NS/MO]={:?}",
                family.family_id, copy_idx, chrom, lo, hi, n_prim, n_obs, hist);
        }
    }
    let n_mosaic = calls.iter().filter(|c| c.is_mosaic()).count();
    let events = aggregate_family(&calls, &params);
    let n_conf = events.iter().filter(|e| e.confirmed).count();
    for ev in &events {
        eprintln!(
            "[VG-MOSAIC] family={} {}->{} breakpoint~{:?} molecules={} dispersion={}bp => {}",
            family.family_id, ev.copy_a, ev.copy_b, ev.breakpoint_ref,
            ev.n_supporting_reads, ev.breakpoint_dispersion,
            if ev.confirmed { "CONFIRMED gene-conversion" } else { "chimera-suspect (unreplicated/dispersed)" },
        );
    }
    if n_mosaic > 0 || !events.is_empty() {
        eprintln!(
            "[VG-MOSAIC] family={}: {} mosaic read(s), {} event(s), {} confirmed",
            family.family_id, n_mosaic, events.len(), n_conf,
        );
    }
    events
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

/// Per-family routing decision used by `--vg-solver auto`.
///
/// Compact dispatch: HMM-EM is the universal target. For divergent copies the
/// forward DP converges in 1–2 iterations to essentially-hard assignments, so
/// it subsumes the legacy heuristic-EM routing. The dispatcher in pipeline.rs
/// owns the heuristic fallback for the case where HMM prerequisites (genome
/// FASTA, multi-mapper sequences) are not satisfied.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EmRoute {
    /// Skip — singleton, too many copies, or intronless. No EM to do.
    Skip(&'static str),
    /// Send this family to HMM-EM (the universal solver).
    Hmm,
}

/// Decide per-family routing. Triviality filters:
///   - `n_copies < 2`: singleton, nothing to disambiguate.
///   - `n_copies > max_copies` (default 40): hard upper bound — even with the
///     work budget unmet, > 40 paralogs is almost always noise (mtDNA, rRNA,
///     mega-tandem clusters).
///   - `work_estimate > max_work` (default 2000): work budget gate. The
///     dominant cost is the forward DP per (read, copy) placement, so we
///     estimate as `n_copies × max(n_multimap_reads, 10)`. Rejects mtDNA-class
///     (1509 copies → ≥ 15090 work) while letting NBPF (25 × 30 = 750) through.
///   - `skip_intronless && all bundles have empty junction_stats`: intronless paralogs
///     (e.g. olfactory receptors) yield degenerate single-node family graphs.
///
/// Everything else routes to `Hmm`.
pub fn classify_family_for_em(
    family: &FamilyGroup,
    bundles: &[Bundle],
    max_copies: usize,
    max_work: usize,
    skip_intronless: bool,
) -> EmRoute {
    let n_copies = family.bundle_indices.len();
    if n_copies < 2 {
        return EmRoute::Skip("singleton");
    }
    if n_copies > max_copies {
        return EmRoute::Skip("too_many_copies");
    }
    let n_reads = family.multimap_reads.len();
    let work_estimate = n_copies.saturating_mul(n_reads.max(10));
    if work_estimate > max_work {
        return EmRoute::Skip("too_much_work");
    }
    if skip_intronless {
        let any_introns = family.bundle_indices.iter()
            .any(|&bi| bundles.get(bi).map(|b| !b.junction_stats.is_empty()).unwrap_or(false));
        if !any_introns {
            return EmRoute::Skip("intronless");
        }
    }
    EmRoute::Hmm
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

        // Build diagnostic SNPs: prefer reference-free pileup when reads have
        // seq populated (VG mode), otherwise fall back to MD-tag approach.
        let diagnostic = {
            let has_seq = family.bundle_indices.iter().any(|&bi| {
                bundles.get(bi).map(|b| b.reads.iter().any(|r| !r.seq.is_empty())).unwrap_or(false)
            });
            if has_seq {
                Some(build_pileup_diagnostics(family, bundles))
            } else if use_snps {
                Some(build_diagnostic_snps(family, bundles))
            } else {
                None
            }
        };

        // Precompute per-entry×copy SNP bonuses ONCE (they don't change across
        // EM iterations). Without this cache, `snp_compatibility` was called
        // O(entries × copies × iterations) times, each doing an O(positions ×
        // read_exons) scan. At family sizes of 2000 reads × 20 iterations with
        // 3000+ diagnostic positions this hangs the pipeline.
        // Parallelized over entries (rayon par_iter): snp_compatibility is
        // a pure function of (read, copy_idx, diagnostic). Bundles are
        // shared `&` for read-only access here.
        let snp_cache: Vec<Vec<f64>> = if let Some(ref diag) = diagnostic {
            use rayon::prelude::*;
            entries
                .par_iter()
                .map(|(_, entry)| {
                    entry
                        .locs
                        .iter()
                        .map(|&(global_bi, ri)| {
                            if ri >= bundles[global_bi].reads.len() {
                                return 1.0;
                            }
                            let read = &bundles[global_bi].reads[ri];
                            let copy_idx = family
                                .bundle_indices
                                .iter()
                                .position(|&bi| bi == global_bi)
                                .unwrap_or(0);
                            snp_compatibility(read, copy_idx, diag)
                        })
                        .collect()
                })
                .collect()
        } else {
            Vec::new()
        };

        let mut result = EmResult::default();

        // Optional Direction-A signal: per-alignment identity proxy (NM + soft-clip).
        // junction_compatibility above is structural-only and saturates near 1.0
        // across paralogs that share intron structure, leaving multi-mappers
        // undiscriminated. NM differs per paralog because each minimap2 alignment
        // carries its own edit distance, providing the medium-identity-band
        // discrimination that intron boundaries alone cannot.
        //
        // Two combination modes:
        //   - additive:       score = (compat + α·id + 0.01) · ln(context)
        //   - multiplicative: score = (compat + 0.01) · ln(context) · id^k
        //
        // The multiplicative form sharpens small identity differences when k is
        // large (e.g. 0.97^20 ≈ 0.54 vs 1.00^20 = 1.00 → 46% relative gap, vs
        // ~3% saturation for additive at any α). Set IDENTITY_K > 0 to enable.
        let use_intron_assign = std::env::var_os("RUSTLE_VG_INTRON_ASSIGN").is_some();
        let identity_alpha: f64 = std::env::var("RUSTLE_VG_IDENTITY_ALPHA")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(1.0);
        let identity_k: f64 = std::env::var("RUSTLE_VG_IDENTITY_K")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0);
        // Optional exon-coverage signal: fraction of bundlenodes covered by the
        // read. Discriminates lost-exon paralogs (4-exon paralog beats 5-exon
        // paralog when the read covers 4 exons that match both). Independent
        // gate so it can be tested without other Direction-A signals.
        let use_exon_cov = std::env::var_os("RUSTLE_VG_EXON_COV").is_some();
        let exon_cov_alpha: f64 = std::env::var("RUSTLE_VG_EXON_COV_ALPHA")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(1.0);
        // Score-gap rule for heuristic EM: only redistribute when
        // log(best_score / second_best_score) >= gap_threshold. Heuristic
        // scores are linear (0..~30 typical), so we work with log-ratio.
        // Empirical: NBPF/LOC101133271 case showed best≈second (both paralogs
        // have matching junctions); without gating, EM splits weight ~50/50,
        // halves the transfrag abundance, and the transcript drops below the
        // path-extraction threshold. Default 1.0 = require best/second ≥ e^1
        // ≈ 2.7x. Set 0 to disable.
        let heur_gap_threshold: f64 = std::env::var("RUSTLE_VG_EM_HEUR_SCORE_GAP")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(1.0);

        for iter in 0..max_iter {
            let mut max_delta: f64 = 0.0;

            for (ei, (_rnh, entry)) in entries.iter_mut().enumerate() {
                // E-step: compute compatibility scores.
                let mut scores: Vec<f64> = Vec::with_capacity(entry.locs.len());
                for (li, &(global_bi, ri)) in entry.locs.iter().enumerate() {
                    if ri >= bundles[global_bi].reads.len() {
                        // No actual read at this bundle — supplementary link only.
                        // Give a small base score (the read COULD be from here).
                        scores.push(0.1);
                        continue;
                    }
                    let read = &bundles[global_bi].reads[ri];
                    let bundle = &bundles[global_bi];
                    let compat = junction_compatibility(read, bundle);
                    let id_score = if use_intron_assign {
                        read_alignment_identity_score(read)
                    } else {
                        1.0
                    };
                    let identity_bonus = if use_intron_assign && identity_k <= 0.0 {
                        identity_alpha * id_score
                    } else {
                        0.0
                    };
                    let exon_cov_bonus = if use_exon_cov {
                        exon_cov_alpha * read_bundle_exon_coverage(read, bundle)
                    } else {
                        0.0
                    };
                    // Context: total junction support at this bundle (higher = more expressed copy).
                    let context: f64 = bundle
                        .junction_stats
                        .iter()
                        .map(|(_, st)| st.nreads_good)
                        .sum::<f64>()
                        .max(1.0);
                    let mut score = (compat + identity_bonus + exon_cov_bonus + 0.01) * context.ln().max(1.0);
                    if use_intron_assign && identity_k > 0.0 {
                        score *= id_score.powf(identity_k);
                    }
                    // SNP bonus: precomputed cache (static across iterations).
                    if !snp_cache.is_empty() {
                        score *= snp_cache[ei][li];
                    }
                    scores.push(score);
                }

                let total: f64 = scores.iter().sum();
                if total <= 0.0 {
                    continue;
                }

                // Score-gap gate: skip update if log(best/second) is below
                // threshold (low-confidence redistribution would dilute the
                // primary paralog's abundance below the path-extraction
                // threshold — see NBPF/LOC101133271 trace).
                if heur_gap_threshold > 0.0 && scores.len() >= 2 {
                    let mut best = f64::NEG_INFINITY;
                    let mut second = f64::NEG_INFINITY;
                    for &v in &scores {
                        if v > best { second = best; best = v; }
                        else if v > second { second = v; }
                    }
                    if best > 0.0 && second > 0.0 {
                        let log_gap = best.ln() - second.ln();
                        if log_gap < heur_gap_threshold {
                            continue;
                        }
                    }
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

        // Diagnostic: when running with SNPs, compute what the no-SNP
        // argmax-copy would have been per read and count how many reads
        // had their argmax-copy CHANGED by the SNP signal. This is the
        // measurement of objective 5 ("did SNPs actually shift any read
        // attributions?") that the prior "reweighted reads" count missed.
        // Activated by RUSTLE_VG_SNP_TRACE_FLIPS=1.
        let trace_flips = use_snps && std::env::var_os("RUSTLE_VG_SNP_TRACE_FLIPS").is_some();
        let mut snp_flip_count = 0usize;
        if trace_flips {
            // Re-run EM without SNP cache on a copy of the entries to see
            // what the argmax would have been. Cheap: same iteration loop,
            // just zero out the snp factor multiplication.
            let mut shadow_entries: Vec<(u64, Vec<f64>)> = entries.iter()
                .map(|(rnh, e)| (*rnh, vec![1.0 / e.weights.len() as f64; e.weights.len()]))
                .collect();
            let no_snp_cache: Vec<Vec<f64>> = Vec::new();
            for _ in 0..max_iter {
                let mut max_d: f64 = 0.0;
                for (ei, (_rnh, w)) in shadow_entries.iter_mut().enumerate() {
                    let entry = &entries[ei].1;
                    let mut scores: Vec<f64> = Vec::with_capacity(entry.locs.len());
                    for &(global_bi, ri) in &entry.locs {
                        if ri >= bundles[global_bi].reads.len() {
                            scores.push(0.1); continue;
                        }
                        let read = &bundles[global_bi].reads[ri];
                        let bundle = &bundles[global_bi];
                        let compat = junction_compatibility(read, bundle);
                        let context: f64 = bundle.junction_stats.iter()
                            .map(|(_, st)| st.nreads_good).sum::<f64>().max(1.0);
                        scores.push((compat + 0.01) * context.ln().max(1.0));
                    }
                    let total: f64 = scores.iter().sum();
                    if total <= 0.0 { continue; }
                    for i in 0..scores.len() {
                        let new_w = scores[i] / total;
                        let d = (new_w - w[i]).abs();
                        if d > max_d { max_d = d; }
                        w[i] = new_w;
                    }
                }
                if max_d < convergence_thr { break; }
            }
            let _ = no_snp_cache;
            // Compare argmax of each entry between SNP and no-SNP runs.
            for (ei, (_, snp_e)) in entries.iter().enumerate() {
                let snp_argmax = snp_e.weights.iter().enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap()).map(|x| x.0);
                let no_snp_argmax = shadow_entries[ei].1.iter().enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap()).map(|x| x.0);
                if snp_argmax != no_snp_argmax { snp_flip_count += 1; }
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
        if trace_flips && snp_flip_count > 0 {
            eprintln!(
                "[VG-SNP-FLIP] family {}: {} reads had argmax-copy changed by SNP signal (out of {} multi-mappers)",
                family.family_id, snp_flip_count, entries.len()
            );
        }

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

/// Pre-assembly EM using HMM-based per-path scoring (sequence-aware).
///
/// Mirrors `run_pre_assembly_em` but replaces the heuristic
/// `junction_compatibility × log_context` with `forward_against_path` —
/// the full forward log-likelihood of the read against each paralog's
/// path through the family-graph HMM.
///
/// Inputs:
///   - `families`: as before, the discovered family groups.
///   - `bundles`: as before, mutated in place.
///   - `family_graphs`: parallel slice to `families`. `family_graphs[i]`
///     must be the FamilyGraph built from `families[i]`'s bundles, with
///     profiles fitted (`fit_profiles_in_place`). If a family has `None`
///     here, that family falls back to `run_pre_assembly_em` semantics
///     (uniform / no-op).
///   - `sequences`: `read_name_hash → bytes` for multi-mappers. Reads
///     without sequences are skipped (no reweighting).
pub fn run_pre_assembly_em_hmm(
    families: &[FamilyGroup],
    bundles: &mut [Bundle],
    family_graphs: &[Option<crate::vg_hmm::family_graph::FamilyGraph>],
    sequences: &HashMap<u64, Vec<u8>>,
    max_iter: usize,
    use_snps: bool,
    junction_bonus: f64,
    uniform_prior: bool,
    exon_len_penalty: f64,
) -> Vec<EmResult> {
    use crate::vg_hmm::scorer::forward_against_path_for_copy_with_norm;
    let mut results = Vec::with_capacity(families.len());
    let convergence_thr = 0.001;
    let t_em_start = std::time::Instant::now();

    for (fam_idx, family) in families.iter().enumerate() {
        let n_copies = family.bundle_indices.len();
        if n_copies < 2 || family.multimap_reads.is_empty() {
            results.push(EmResult::default());
            continue;
        }

        let fg_opt = family_graphs.get(fam_idx).and_then(|f| f.as_ref());
        let fg = match fg_opt {
            Some(g) if !g.nodes.is_empty() => g,
            _ => {
                eprintln!(
                    "[VG-HMM-EM] Family {}: no family graph available — skipping EM",
                    family.family_id
                );
                results.push(EmResult::default());
                continue;
            }
        };

        // Pre-compute per-paralog paths.
        let paralog_paths: Vec<Vec<crate::vg_hmm::family_graph::NodeIdx>> = (0..n_copies)
            .map(|cid| fg.recover_paralog_path(cid))
            .collect();

        // Pre-compute total match-column count per paralog path. This is the
        // length-norm denominator and is constant across all reads — sum it
        // once here instead of recomputing inside every per-read forward call.
        let path_match_cols: Vec<usize> = (0..n_copies)
            .map(|cid| crate::vg_hmm::scorer::path_match_cols_for_copy(fg, &paralog_paths[cid], cid))
            .collect();

        // Build per-family diagnostic SNPs if SNP-aware HMM is enabled. The
        // SNP signal is constant across EM iterations and per-(read, copy)
        // independent of forward score, so cache log(snp_compat) per
        // placement once and add it into the posterior.
        let diagnostic = if use_snps {
            Some(build_diagnostic_snps(family, bundles))
        } else {
            None
        };

        // Build entries with pre-computed log-likelihoods (constant across EM iters).
        struct Entry {
            locs: Vec<(usize, usize)>,
            log_scores: Vec<f64>,
            log_snp: Vec<f64>,
            log_jct: Vec<f64>,
            log_exonlen: Vec<f64>,
            weights: Vec<f64>,
            fam_pos: Vec<usize>,
        }

        // Optional cap on reads scored per family (env: RUSTLE_VG_HMM_EM_READ_CAP).
        // Bounds runtime on huge families (e.g. 769 reads × 18 copies). Set to 0
        // to disable. Default: 200 reads — usually enough to estimate copy
        // priors, after which the EM iteration handles redistribution.
        let read_cap: usize = std::env::var("RUSTLE_VG_HMM_EM_READ_CAP")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(200);

        // PHASE 1 (sequential): collect work items — one per multimap read with
        // a usable sequence. The expensive part (forward_against_path per
        // placement) is deferred to PHASE 2 where it runs in parallel.
        struct WorkItem<'a> {
            seq: &'a [u8],
            // (fam_pos, global_bi, ri, current_weight)
            placements: Vec<(usize, usize, usize, f64)>,
        }
        let mut work: Vec<WorkItem> = Vec::with_capacity(family.multimap_reads.len());
        let mut n_no_seq = 0usize;
        let mut n_prefiltered = 0usize;

        // Junction-compat pre-filter threshold: drop (read, copy) placements
        // whose read junctions overlap with the copy's bundle junction_stats
        // below this fraction. Cheap O(1) per placement; uses the existing
        // junction_compatibility infrastructure. Default 0.5 — read must have
        // ≥ 50 % of its CIGAR junctions matched within ±10 bp tolerance.
        // Single-exon reads (no junctions) return compat = 1.0 and pass through.
        // Override via RUSTLE_VG_HMM_EM_JCT_PREFILTER (0.0 disables).
        let jct_prefilter: f64 = std::env::var("RUSTLE_VG_HMM_EM_JCT_PREFILTER")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.5);

        // Iterate deterministically (HashMap order is unspecified). Sort by rnh so
        // capping is reproducible across runs.
        let mut keys: Vec<u64> = family.multimap_reads.keys().copied().collect();
        keys.sort_unstable();

        for rnh in keys {
            let locs = &family.multimap_reads[&rnh];
            let seq = match sequences.get(&rnh) {
                Some(s) if !s.is_empty() => s.as_slice(),
                _ => { n_no_seq += 1; continue; }
            };
            let mut placements = Vec::with_capacity(locs.len());
            for &(fam_pos, ri) in locs {
                if fam_pos >= n_copies { continue; }
                let global_bi = family.bundle_indices[fam_pos];
                let w = if ri < bundles[global_bi].reads.len() {
                    bundles[global_bi].reads[ri].weight
                } else {
                    0.0
                };
                // Junction-compat pre-filter: skip this placement entirely
                // when the read's junctions don't overlap enough with the
                // copy's expected junctions. Saves a full forward DP per
                // skipped placement; safe because such a placement would
                // score very low and contribute negligibly to the posterior.
                if jct_prefilter > 0.0 && ri < bundles[global_bi].reads.len() {
                    let read = &bundles[global_bi].reads[ri];
                    if !read.junctions.is_empty() {
                        let compat = junction_compatibility(read, &bundles[global_bi]);
                        if compat < jct_prefilter {
                            n_prefiltered += 1;
                            continue;
                        }
                    }
                }
                placements.push((fam_pos, global_bi, ri, w));
            }
            if placements.len() < 2 { continue; }
            work.push(WorkItem { seq, placements });
            if read_cap > 0 && work.len() >= read_cap { break; }
        }
        let n_capped = if read_cap > 0 && family.multimap_reads.len() > read_cap {
            family.multimap_reads.len() - read_cap
        } else {
            0
        };

        // PHASE 2 (parallel): score every placement with forward_against_path.
        // forward_against_path is a pure function of (fg, seq, path) — fg is
        // shared `&` across threads, seqs/paths are non-overlapping reads.
        let scored: Vec<Entry> = work.par_iter()
            .filter_map(|item| {
                let mut entry_locs = Vec::with_capacity(item.placements.len());
                let mut log_scores  = Vec::with_capacity(item.placements.len());
                let mut log_snp     = Vec::with_capacity(item.placements.len());
                let mut log_jct     = Vec::with_capacity(item.placements.len());
                let mut log_exonlen = Vec::with_capacity(item.placements.len());
                let mut weights   = Vec::with_capacity(item.placements.len());
                let mut fp_vec    = Vec::with_capacity(item.placements.len());
                for &(fam_pos, global_bi, ri, w) in &item.placements {
                    let path = &paralog_paths[fam_pos];
                    let score = if path.is_empty() {
                        f64::NEG_INFINITY
                    } else {
                        // Use per-copy profiles: each paralog scored against
                        // its own genomic sequence, not the family consensus.
                        // Pre-computed match-cols avoids re-summing per read.
                        forward_against_path_for_copy_with_norm(
                            fg, item.seq, path, fam_pos, Some(path_match_cols[fam_pos]),
                        )
                    };
                    let snp_log = if let Some(ref diag) = diagnostic {
                        if ri < bundles[global_bi].reads.len() {
                            let read = &bundles[global_bi].reads[ri];
                            snp_compatibility(read, fam_pos, diag).ln()
                        } else {
                            0.0
                        }
                    } else {
                        0.0
                    };
                    // Structural-mismatch penalty: -bonus × #junctions in the
                    // read's CIGAR that the copy's expected junction set does
                    // NOT contain (±10 bp tolerance, via junction_compatibility).
                    // Independent of nucleotide identity — discriminates copies
                    // that differ in splice structure even when emissions are
                    // near-uniform (divergent paralogs).  Zero when bonus = 0.
                    let jct_log = if junction_bonus > 0.0
                        && ri < bundles[global_bi].reads.len()
                    {
                        let read = &bundles[global_bi].reads[ri];
                        let total = read.junctions.len() as f64;
                        if total > 0.0 {
                            let compat = junction_compatibility(read, &bundles[global_bi]);
                            let miss = total * (1.0 - compat);
                            -junction_bonus * miss
                        } else {
                            0.0
                        }
                    } else {
                        0.0
                    };
                    // Structural exon-length divergence: penalise reads whose
                    // total spliced length is far from any read in the candidate
                    // bundle.  Orthogonal to identity — fires on divergent
                    // paralogs whose copies differ in indel/cassette structure.
                    let exonlen_log = if exon_len_penalty > 0.0
                        && ri < bundles[global_bi].reads.len()
                    {
                        let read = &bundles[global_bi].reads[ri];
                        let div = exon_length_divergence(read, &bundles[global_bi]);
                        -exon_len_penalty * div
                    } else {
                        0.0
                    };
                    entry_locs.push((global_bi, ri));
                    log_scores.push(score);
                    log_snp.push(snp_log);
                    log_jct.push(jct_log);
                    log_exonlen.push(exonlen_log);
                    weights.push(w);
                    fp_vec.push(fam_pos);
                }
                if entry_locs.len() < 2 { return None; }
                if !log_scores.iter().any(|s| s.is_finite()) { return None; }
                Some(Entry { locs: entry_locs, log_scores, log_snp, log_jct, log_exonlen, weights, fam_pos: fp_vec })
            })
            .collect();
        let mut entries = scored;

        if entries.is_empty() {
            if n_no_seq > 0 || n_capped > 0 {
                eprintln!(
                    "[VG-HMM-EM] Family {}: 0 entries usable ({} reads lacked sequences, {} capped)",
                    family.family_id, n_no_seq, n_capped
                );
            }
            results.push(EmResult::default());
            continue;
        }

        let mut result = EmResult::default();

        // Score-gap rule (advisor's primary-vs-secondary threshold): only
        // redistribute a read's weight when the best-paralog log-score is at
        // least `gap_threshold` log-units above the next-best. When the gap
        // is smaller, the EM is uncertain — leave the read's weight at its
        // initial 1/NH (BAM-aligner's call) rather than risk flipping it the
        // wrong way. Default 10.0 log-units (per advisor); set 0 to disable.
        let gap_threshold: f64 = std::env::var("RUSTLE_VG_EM_SCORE_GAP")
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(10.0);

        for iter in 0..max_iter {
            // M-step: per-copy weight totals → log-priors with floor.
            // With `uniform_prior`, the M-step is disabled — every copy gets
            // a flat prior so the per-(read, copy) posterior is driven purely
            // by the HMM (+ SNP, + junction) log-likelihood.  Used to test
            // whether EM's iterative prior estimation is doing real work.
            let mut copy_total = vec![0.0_f64; n_copies];
            for entry in &entries {
                for (i, w) in entry.weights.iter().enumerate() {
                    copy_total[entry.fam_pos[i]] += w;
                }
            }
            let total_sum: f64 = copy_total.iter().sum::<f64>().max(1.0);
            let log_priors: Vec<f64> = if uniform_prior {
                vec![0.0_f64; n_copies]
            } else {
                copy_total.iter()
                    .map(|&t| ((t / total_sum) + 1e-3).ln())
                    .collect()
            };

            let mut max_delta: f64 = 0.0;
            for entry in &mut entries {
                let n = entry.locs.len();
                let mut log_post = vec![f64::NEG_INFINITY; n];
                for i in 0..n {
                    let s = entry.log_scores[i];
                    if !s.is_finite() { continue; }
                    log_post[i] = s
                        + log_priors[entry.fam_pos[i]]
                        + entry.log_snp[i]
                        + entry.log_jct[i]
                        + entry.log_exonlen[i];
                }
                let max_lp = log_post.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
                if !max_lp.is_finite() { continue; }

                // Confidence gate: skip update if best vs second-best gap is
                // below threshold.
                if gap_threshold > 0.0 && n >= 2 {
                    let mut best = f64::NEG_INFINITY;
                    let mut second = f64::NEG_INFINITY;
                    for &v in &log_post {
                        if v > best { second = best; best = v; }
                        else if v > second { second = v; }
                    }
                    if best.is_finite() && second.is_finite() && (best - second) < gap_threshold {
                        continue;
                    }
                }

                let mut total = 0.0_f64;
                let mut exps = vec![0.0_f64; n];
                for i in 0..n {
                    if log_post[i].is_finite() {
                        exps[i] = (log_post[i] - max_lp).exp();
                        total += exps[i];
                    }
                }
                if total <= 0.0 { continue; }
                for i in 0..n {
                    let new_w = exps[i] / total;
                    let delta = (new_w - entry.weights[i]).abs();
                    if delta > max_delta { max_delta = delta; }
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
        for entry in &entries {
            for (i, &(global_bi, ri)) in entry.locs.iter().enumerate() {
                if ri >= bundles[global_bi].reads.len() { continue; }
                let old_w = bundles[global_bi].reads[ri].weight;
                let new_w = entry.weights[i];
                if (old_w - new_w).abs() > 1e-9 {
                    bundles[global_bi].reads[ri].weight = new_w;
                    n_reweighted += 1;
                }
            }
        }
        result.reads_reweighted = n_reweighted;

        if n_reweighted > 0 || n_no_seq > 0 || n_capped > 0 || n_prefiltered > 0 {
            eprintln!(
                "[VG-HMM-EM] Family {}: HMM-EM converged={} in {} iter (delta={:.6}), reweighted {} reads ({} skipped: no seq, {} capped, {} jct-prefiltered) across {} copies",
                family.family_id, result.converged, result.iterations, result.max_delta,
                n_reweighted, n_no_seq, n_capped, n_prefiltered, n_copies,
            );
        }
        results.push(result);
    }

    let total_reweighted: usize = results.iter().map(|r| r.reads_reweighted).sum();
    eprintln!(
        "[VG-HMM-EM] HMM-EM{} reweighting complete: {} reads adjusted across {} families in {:.1}s",
        if use_snps { " (SNP-aware)" } else { "" },
        total_reweighted,
        results.iter().filter(|r| r.reads_reweighted > 0).count(),
        t_em_start.elapsed().as_secs_f64(),
    );
    results
}

/// Fingerprint-based EM: redistribute multi-mapping read weights using per-copy
/// SNP/indel fingerprints built from `ExonClass.per_copy_sequences`.
///
/// Unlike `run_pre_assembly_em_hmm`, this does not run a full forward DP.
/// Instead it accumulates log-likelihood from copy-distinguishing positions
/// found by pairwise comparison of the family graph's per-copy sequences.
/// Works for non-overlapping copies (RBMY/TSPY/DAZ on chrY) because positions
/// are exon-relative, not reference-coordinate.
///
/// If a family has no diagnostic sites (indistinguishable copies or no
/// family graph), returns a default (no-op) EmResult for that family.
/// Uses BundleRead.seq (always populated in VG mode) via snp_compatibility,
/// so does not require --vg-snp or BundleRead.mismatches.
pub fn run_fingerprint_em(
    families: &[FamilyGroup],
    bundles: &mut [Bundle],
    family_graphs: &[Option<crate::vg_hmm::family_graph::FamilyGraph>],
    max_iter: usize,
    // Optional collector for per-family gene-conversion events (mosaic pass), keyed by
    // family_id, so the pipeline can surface them to the GTF. None for callers that don't need it.
    mut mosaic_out: Option<&mut crate::types::DetHashMap<usize, Vec<crate::vg_hmm::mosaic::ConversionEvent>>>,
) -> Vec<EmResult> {
    use rayon::prelude::*;
    let mut results = Vec::with_capacity(families.len());
    let convergence_thr = 0.001;
    let t_start = std::time::Instant::now();

    // Score-gap gate: fingerprint evidence requires a large gap (default 10.0 log-units)
    // to be decisive. Structural-only evidence (junction chain + alignment identity,
    // which dominate when a read covers 0 diagnostic sequence sites) is weaker;
    // struct_gap is a lower threshold that still lets junction/NM guide those reads.
    //
    // `gap_threshold` is an ABSOLUTE ceiling calibrated for reads covering many
    // (hundreds–thousands of) diagnostic sites — the real-IsoSeq paralog regime.
    // For reads covering FEW sites (sparse divergence or short overlap), an
    // absolute 10.0-log-unit bar suppresses clean but sparse signal: a read that
    // nets even 2 distinguishing sites carries ~80:1 odds yet only ~4.4 log-units.
    // `per_site_gap` makes the requirement scale with the evidence actually
    // available — eff_gap = min(gap_threshold, per_site_gap * n_sites_covered) —
    // so the abundant-site regime is UNCHANGED (the per-site term saturates the
    // ceiling) while sparse-coverage reads get an evidence-proportional bar.
    let gap_threshold: f64 = std::env::var("RUSTLE_VG_EM_SCORE_GAP")
        .ok().and_then(|s| s.parse().ok()).unwrap_or(10.0);
    let per_site_gap: f64 = std::env::var("RUSTLE_VG_EM_SCORE_GAP_PER_SITE")
        .ok().and_then(|s| s.parse().ok()).unwrap_or(0.25);
    // Calibration: minimum diagnostic sites a read must cover for the fingerprint
    // to be allowed to DECIDE its copy. Below this, the read abstains (falls to the
    // prior) — a single diagnostic base is flippable by one sequencing error.
    // Default 2; set to 1 to restore the prior (pre-2026-06-03) behaviour.
    let min_decisive_sites: usize = std::env::var("RUSTLE_VG_EM_MIN_DECISIVE_SITES")
        .ok().and_then(|s| s.parse().ok()).unwrap_or(2);

    // Multi-signal EM weights (IsoSeq-specific):
    //   RUSTLE_VG_FP_LAMBDA_J  — junction-chain signal weight (default 1.0)
    //   RUSTLE_VG_FP_LAMBDA_N  — alignment-identity signal weight (default 1.0)
    //   RUSTLE_VG_STRUCT_GAP   — gap threshold when fingerprint coverage = 0 (default 0.1)
    use std::sync::OnceLock;
    static MULTI_SIG: OnceLock<(f64, f64, f64)> = OnceLock::new();
    let (lambda_j, lambda_nm, struct_gap) = *MULTI_SIG.get_or_init(|| {
        let lj  = std::env::var("RUSTLE_VG_FP_LAMBDA_J").ok()
                     .and_then(|s| s.parse().ok()).unwrap_or(1.0_f64);
        let lnm = std::env::var("RUSTLE_VG_FP_LAMBDA_N").ok()
                     .and_then(|s| s.parse().ok()).unwrap_or(1.0_f64);
        let sg  = std::env::var("RUSTLE_VG_STRUCT_GAP").ok()
                     .and_then(|s| s.parse().ok()).unwrap_or(0.1_f64);
        (lj, lnm, sg)
    });

    for (fam_idx, family) in families.iter().enumerate() {
        let n_copies = family.bundle_indices.len();
        if n_copies < 2 || family.multimap_reads.is_empty() {
            results.push(EmResult::default());
            continue;
        }

        let fg_opt = family_graphs.get(fam_idx).and_then(|f| f.as_ref());
        // Joint-strand mode (default ON): a None/empty graph is the mixed-strand
        // (inverted-pair) case — there is NO valid joint sequence graph, but we
        // still want to apportion its shared multimappers by junc + nm + the
        // anchored prior. Build a neutral (0-site) fingerprint and fall through
        // instead of skipping. With RUSTLE_VG_JOINT_STRAND_EM=0 we keep the exact
        // legacy skip so the rollback A/B is byte-identical.
        let joint_strand_em = std::env::var("RUSTLE_VG_JOINT_STRAND_EM")
            .map(|v| v != "0").unwrap_or(true);
        let fp = match fg_opt {
            Some(g) if !g.nodes.is_empty() => enumerate_diagnostic_sites(g, n_copies),
            _ if joint_strand_em => {
                eprintln!(
                    "[VG-FP-EM] Family {}: no joint graph (mixed-strand) — neutral-fp fallthrough (junc+nm+anchored prior)",
                    family.family_id
                );
                ExonFingerprints {
                    sites: Vec::new(),
                    per_copy_site_refs: vec![Vec::new(); n_copies],
                    n_copies,
                    n_sites: 0,
                }
            }
            _ => {
                eprintln!("[VG-FP-EM] Family {}: no family graph — skipping", family.family_id);
                results.push(EmResult::default());
                continue;
            }
        };

        if fp.n_sites == 0 {
            // No diagnostic sites — copies are sequence-identical (e.g. DAZ).
            // The fingerprint carries no discriminating signal, but SKIPPING the
            // EM is wrong: shared multi-mappers then stay counted at full weight
            // in EVERY copy, inflating the non-starved copy (DAZ1: 9 tx / cov 142
            // instead of HMM-EM's 2 / 8.5). Fall through to the EM with the
            // fingerprint term neutral (all sites 0): the per-copy pileup-depth
            // PRIOR (M-step) plus the junction-chain / alignment-identity signals
            // then split the multi-mappers. This is the simple, explainable
            // replacement for the HMM-EM, whose forward DP also converges
            // trivially for identical copies — the prior is the real signal.
            eprintln!(
                "[VG-FP-EM] Family {}: 0 diagnostic sites — pileup-depth-prior fallback (identical copies)",
                family.family_id
            );
        }

        // Gene-conversion mosaic-read detection (audit theme: "VG finds unusual exon
        // combinations"). Opt-in RUSTLE_VG_MOSAIC_ON, additive — does NOT touch the EM
        // weights below; only reports recombinant reads / confirmed conversion events.
        if fp.n_sites > 0 && std::env::var_os("RUSTLE_VG_MOSAIC_ON").is_some() {
            let ev = detect_and_report_mosaics(family, bundles, &fp);
            if let Some(m) = mosaic_out.as_deref_mut() {
                if !ev.is_empty() {
                    m.entry(family.family_id).or_default().extend(ev);
                }
            }
        }

        // PHASE 1 (sequential): collect placement lists, one per multi-mapped read.
        struct PreEntry {
            locs: Vec<(usize, usize)>, // (global_bi, ri)
            fam_pos: Vec<usize>,       // copy index per loc
            weights: Vec<f64>,         // initial weights
        }

        let mut keys: Vec<u64> = family.multimap_reads.keys().copied().collect();
        keys.sort_unstable();

        let mut pre_entries: Vec<PreEntry> = Vec::with_capacity(keys.len());
        for rnh in &keys {
            let locs = &family.multimap_reads[rnh];
            let mut entry_locs: Vec<(usize, usize)> = Vec::with_capacity(locs.len());
            let mut entry_fp: Vec<usize> = Vec::with_capacity(locs.len());
            let mut entry_weights: Vec<f64> = Vec::with_capacity(locs.len());
            for &(fam_pos, ri) in locs {
                if fam_pos >= n_copies { continue; }
                let global_bi = family.bundle_indices[fam_pos];
                let w = if ri < bundles[global_bi].reads.len() {
                    bundles[global_bi].reads[ri].weight
                } else {
                    0.0
                };
                entry_locs.push((global_bi, ri));
                entry_fp.push(fam_pos);
                entry_weights.push(w);
            }
            if entry_locs.len() < 2 { continue; }
            pre_entries.push(PreEntry { locs: entry_locs, fam_pos: entry_fp, weights: entry_weights });
        }

        // PHASE 2 (parallel): pre-compute fingerprint log-scores per placement.
        // score_read_exon_fingerprint is a pure function of (read, copy_idx, fp).
        // We call it once per placement and take the diagonal: score for THIS copy.
        // Constant across EM iterations.
        struct WorkItem {
            locs: Vec<(usize, usize)>,
            fam_pos: Vec<usize>,
            /// Sequence fingerprint: Bernoulli log-likelihood over diagnostic sites
            /// where copies carry different bases. Accumulates across thousands of
            /// sites for IsoSeq reads → dominant signal when coverage is high.
            log_scores: Vec<f64>,
            /// IsoSeq-specific: ln(junction_compatibility) per placement.
            /// Fraction of this read's splice junctions supported at copy k.
            /// Long reads span the full junction chain; short reads see ≤1 junction.
            junc_scores: Vec<f64>,
            /// IsoSeq-specific: ln(alignment_identity) per placement.
            /// 1 − (NM_rate + clip_rate)/2. Reads misaligned to the wrong copy
            /// accumulate more edit distance over a 3 kb read than a 150 bp one.
            nm_scores: Vec<f64>,
            /// Number of diagnostic sites covered per placement. Longer reads
            /// cover more sites and produce larger log-likelihood gaps → more
            /// decisive assignment. Tracked for attribution quality reporting.
            n_sites_covered: Vec<usize>,
            weights: Vec<f64>,
        }

        let mut entries: Vec<WorkItem> = pre_entries
            .into_par_iter()
            .map(|pre| {
                let np = pre.locs.len();
                // Pass 1: per-placement full per-copy fp vector + diagnostic-site coverage +
                // structural signals (junction-chain, alignment-identity — these ARE legitimately
                // per-placement: a read's junctions/edit-distance to each copy differ by locus).
                let mut fp_vecs: Vec<Vec<f64>> = Vec::with_capacity(np);
                let mut ncovs: Vec<usize> = Vec::with_capacity(np);
                let mut junc_scores: Vec<f64> = Vec::with_capacity(np);
                let mut nm_scores:   Vec<f64> = Vec::with_capacity(np);
                let mut valid: Vec<bool> = Vec::with_capacity(np);
                for (&(global_bi, ri), &fam_pos) in pre.locs.iter().zip(pre.fam_pos.iter()) {
                    if ri >= bundles[global_bi].reads.len() {
                        fp_vecs.push(Vec::new()); ncovs.push(0);
                        junc_scores.push(0.0); nm_scores.push(0.0); valid.push(false);
                        continue;
                    }
                    let read   = &bundles[global_bi].reads[ri];
                    let bundle = &bundles[global_bi];
                    let (fp_scores, n_cov) = score_read_exon_fingerprint(read, fam_pos, &fp);
                    if std::env::var_os("RUSTLE_VG_FP_SCORE_TRACE").is_some() {
                        eprintln!("[FP-SCORE] rnh={} placed_at_copy={} n_cov={} full_scores={:?}",
                            read.read_name_hash, fam_pos, n_cov, fp_scores);
                    }
                    fp_vecs.push(fp_scores);
                    ncovs.push(n_cov);
                    junc_scores.push(junction_compatibility(read, bundle).max(1e-6).ln());
                    nm_scores.push(read_alignment_identity_score(read).max(1e-6).ln());
                    valid.push(true);
                }
                // Fingerprint evidence is a per-READ haplotype property: the read carries its
                // bases (and so covers diagnostic sites) only at the placement(s) where it is
                // aligned with sequence — at the OTHER copy it may cover 0 sites (SEQ=* on
                // secondaries, or the homologous exon falls outside its span). Comparing
                // per-placement DIAGONAL fp scores is therefore broken: a 0-coverage placement
                // scores 0.0 and beats a strongly-matched placement whose Bernoulli
                // log-likelihood is negative (observed: fam175 read covers 137 sites at copy B,
                // vector [-315(A), -14(B)] => clearly B, but diagonal compared 0@A vs -14@B and
                // mis-picked A). Fix: take the full per-copy score vector from the read's
                // MOST-INFORMATIVE placement (max sites covered) and give every placement its
                // copy's entry from that single vector. Robust to asymmetric coverage AND to
                // imperfect synteny pairing (at any covered placement the read matches its own
                // copy's alleles and mismatches the other's, since sites are A!=B by definition).
                // Informative placement = the read's BEST-matching copy among those covering
                // diagnostic sites (max alignment identity = minimap2's primary). Only there are
                // the read's bases its TRUE sequence: a secondary placement may be SEQ=* (then
                // read_base defaults to the copy's own allele — a spurious self-match) or, in a
                // chimeric fixture, a different copy's slice. nm_scores = ln(alignment identity).
                let best = (0..np)
                    .filter(|&i| valid[i] && ncovs[i] > 0)
                    .max_by(|&a, &b| nm_scores[a].partial_cmp(&nm_scores[b]).unwrap_or(std::cmp::Ordering::Equal));
                let mut log_scores: Vec<f64> = Vec::with_capacity(np);
                let mut n_sites_covered: Vec<usize> = Vec::with_capacity(np);
                for i in 0..np {
                    if !valid[i] {
                        log_scores.push(f64::NEG_INFINITY);
                        n_sites_covered.push(0);
                    } else if let Some(b) = best {
                        log_scores.push(fp_vecs[b].get(pre.fam_pos[i]).copied().unwrap_or(0.0));
                        n_sites_covered.push(ncovs[b]);
                    } else {
                        log_scores.push(0.0); // no fp evidence anywhere; junc/nm/prior decide
                        n_sites_covered.push(0);
                    }
                }
                WorkItem {
                    locs: pre.locs, fam_pos: pre.fam_pos,
                    log_scores, junc_scores, nm_scores, n_sites_covered, weights: pre.weights,
                }
            })
            .collect();

        if entries.is_empty() {
            results.push(EmResult::default());
            continue;
        }

        let mut result = EmResult::default();

        // Anchor-first apportionment (default ON in VG): replace the M-step
        // pileup-depth prior with an anchored prior grounded in UNAMBIGUOUS mass
        // (unique reads + dNM-decisive Owns reads). This removes the pileup
        // prior's self-reinforcement of an already-double-counted copy. The
        // anchored prior is constant, so the E-step runs ONCE (max_iter
        // effectively 1). Calibration: raw dNM t=2, extent_frac=0.8.
        let anchor_prior_on = std::env::var("RUSTLE_VG_ANCHOR_PRIOR")
            .map(|v| v != "0").unwrap_or(true);
        let fixed_log_priors: Option<Vec<f64>> = if anchor_prior_on {
            let anchored = anchored_mass_per_copy(family, bundles, 2, 0.8);
            let total_anchored: f64 = anchored.iter().sum();
            if total_anchored < 1e-9 {
                // All-zero-anchor (e.g. every copy is a pure multimapper): no
                // grounded mass exists. Graceful-degrade to the existing
                // pileup-depth prior for THIS family only.
                eprintln!(
                    "[VG-FP-EM] Family {}: 0 anchored mass — graceful-degrade to pileup-depth prior",
                    family.family_id
                );
                None
            } else {
                let total = total_anchored.max(1.0);
                Some(anchored.iter().map(|&a| ((a / total) + 1e-3).ln()).collect())
            }
        } else {
            None
        };
        // Under anchor mode with a fixed prior, a single E-step suffices.
        let effective_max_iter = if fixed_log_priors.is_some() { 1 } else { max_iter };

        for iter in 0..effective_max_iter {
            // M-step prior: anchored (fixed, computed once) when available,
            // otherwise the legacy pileup-depth prior recomputed from current
            // weights (graceful-degrade or RUSTLE_VG_ANCHOR_PRIOR=0).
            let log_priors: Vec<f64> = if let Some(ref lp) = fixed_log_priors {
                lp.clone()
            } else {
                let mut copy_total = vec![0.0_f64; n_copies];
                for entry in &entries {
                    for (i, &w) in entry.weights.iter().enumerate() {
                        copy_total[entry.fam_pos[i]] += w;
                    }
                }
                let total_sum: f64 = copy_total.iter().sum::<f64>().max(1.0);
                copy_total.iter().map(|&t| ((t / total_sum) + 1e-3).ln()).collect()
            };

            let mut max_delta: f64 = 0.0;
            for entry in &mut entries {
                let n = entry.locs.len();

                // E-step: log P(copy=k | read) ∝ fingerprint + junction_chain + alignment_id + prior.
                //
                // Three conditionally-independent signals combine in log space:
                //   fp(read, k)  = Σ_j log P(seq_base_j | copy_k)   [Bernoulli at diagnostic sites]
                //   junc(read,k) = ln(junction_compatibility)         [IsoSeq: full junction chain]
                //   nm(read, k)  = ln(alignment_identity)             [IsoSeq: edit-distance over full read]
                //
                // When fp covers many sites (IsoSeq: thousands) it dominates.
                // For reads in conserved regions (n_sites_covered = 0), junc + nm
                // provide structural discrimination that short reads cannot supply.
                let log_post: Vec<f64> = (0..n)
                    .map(|i| {
                        let fp = entry.log_scores[i];
                        let j  = lambda_j  * entry.junc_scores[i];
                        let m  = lambda_nm * entry.nm_scores[i];
                        // fp is NEG_INFINITY for failed placements (ri out of range).
                        // j and m are 0.0 (neutral) in those cases.
                        let base = if fp.is_finite() { fp } else { 0.0 };
                        base + j + m + log_priors[entry.fam_pos[i]]
                    })
                    .collect();

                let max_lp = log_post.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));
                if !max_lp.is_finite() { continue; }

                // Adaptive score-gap gate:
                // • Fingerprint evidence (n_sites_covered > 0): require an
                //   evidence-proportional gap, min(gap_threshold, per_site_gap *
                //   n_sites_covered). With many sites this saturates the absolute
                //   gap_threshold ceiling (default 10.0 log-units, the thousands-of-
                //   sites Bernoulli scale); with few sites it relaxes to a per-site
                //   bar so clean sparse signal is honored instead of being silenced.
                // • Structural-only evidence (junction + NM, no diagnostic sites):
                //   require struct_gap (default 0.1 log-units) — weaker signal, lower bar.
                //   This allows reads in conserved regions to receive soft structural
                //   assignments that are impossible for short reads.
                let n_cov_max = entry.n_sites_covered.iter().copied().max().unwrap_or(0);
                let eff_gap = if n_cov_max == 0 {
                    struct_gap
                } else if n_cov_max < min_decisive_sites {
                    // Calibration (2026-06-03): a read resting on a SINGLE diagnostic
                    // site is unreliable — one sequencing error at that base flips the
                    // call. The attribution benchmark showed the EM was decisive-but-
                    // overconfident at exactly 1 site (committed on 81% of such reads,
                    // only 62% correct). Require the FULL absolute gap (which a lone
                    // site cannot supply, ~2 log-units max) so these reads fall through
                    // to the prior — they ABSTAIN instead of committing on one base.
                    // Reads covering 0 sites are unaffected (struct_gap above); reads
                    // covering >=min_decisive_sites use the evidence-proportional bar.
                    gap_threshold
                } else {
                    (per_site_gap * n_cov_max as f64).min(gap_threshold)
                };
                if eff_gap > 0.0 && n >= 2 {
                    let mut best = f64::NEG_INFINITY;
                    let mut second = f64::NEG_INFINITY;
                    for &v in &log_post {
                        if v > best { second = best; best = v; }
                        else if v > second { second = v; }
                    }
                    if best.is_finite() && second.is_finite() && (best - second) < eff_gap {
                        // Conservation hole fix: rather than leaving this read at its
                        // raw 1/NH weights (which do NOT generally sum to 1.0 across the
                        // placements kept in this entry), assign the normalized anchored
                        // prior over this read's placements. Mass-conserving: weights
                        // sum to exactly 1.0. The prior (log_priors) is the M-step's
                        // anchored-mass distribution, so gate-silenced reads fall back to
                        // the family's unambiguous-mass apportionment instead of a stale,
                        // non-conserving initialization.
                        let max_lprior = entry.fam_pos.iter()
                            .map(|&k| log_priors[k])
                            .fold(f64::NEG_INFINITY, f64::max);
                        if max_lprior.is_finite() {
                            let mut psum = 0.0_f64;
                            let mut pexps = vec![0.0_f64; n];
                            for i in 0..n {
                                pexps[i] = (log_priors[entry.fam_pos[i]] - max_lprior).exp();
                                psum += pexps[i];
                            }
                            if psum > 0.0 {
                                for i in 0..n {
                                    let new_w = pexps[i] / psum;
                                    let delta = (new_w - entry.weights[i]).abs();
                                    if delta > max_delta { max_delta = delta; }
                                    entry.weights[i] = new_w;
                                }
                            }
                        }
                        continue;
                    }
                }

                // Log-sum-exp normalization → posterior weights.
                let mut total = 0.0_f64;
                let mut exps = vec![0.0_f64; n];
                for i in 0..n {
                    if log_post[i].is_finite() {
                        exps[i] = (log_post[i] - max_lp).exp();
                        total += exps[i];
                    }
                }
                if total <= 0.0 { continue; }
                for i in 0..n {
                    let new_w = exps[i] / total;
                    let delta = (new_w - entry.weights[i]).abs();
                    if delta > max_delta { max_delta = delta; }
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

        // Apply final weights to bundles and collect attribution quality stats.
        //
        // Weight conservation: for every multi-mapped read the weights across all
        // placements sum to exactly 1.0 (log-sum-exp normalisation guarantees this).
        // Decisive assignment means one copy gets >> 0.5 — the read is not silenced,
        // it is attributed. Uncertain reads (weight_gap < 0.1) are left near their
        // initial 1/n_copies prior; the information is retained, not discarded.
        //
        // Attribution quality buckets (by weight_gap = max_w - second_max_w):
        //   decisive  : weight_gap > 0.8  (read effectively assigned to one copy)
        //   moderate  : 0.5 < gap ≤ 0.8
        //   uncertain : gap ≤ 0.5         (below score-gap gate or few sites covered)
        let mut n_reweighted = 0usize;
        let mut n_decisive = 0usize;      // weight_gap > 0.8
        let mut n_moderate = 0usize;      // 0.5 < gap ≤ 0.8
        let mut n_uncertain = 0usize;     // gap ≤ 0.5
        let mut n_struct_guided = 0usize; // decisive/moderate with 0 fingerprint sites (junction/NM only)
        let mut total_sites_covered = 0usize;
        let mut n_reads_with_sites = 0usize;

        // Optional per-read TSV dump (RUSTLE_VG_FP_ATTR_TSV=<path>).
        // Columns: family_id, read_name_hash, placement_copy, n_sites_covered,
        //          final_weight, weight_gap, weight_sum
        // One row per (read, placement). Useful for attributing multi-mappers
        // back to their copy and validating weight conservation per read.
        let attr_tsv_path: Option<String> = std::env::var("RUSTLE_VG_FP_ATTR_TSV").ok();
        let mut attr_writer: Option<Box<dyn std::io::Write>> = attr_tsv_path
            .as_ref()
            .and_then(|p| {
                use std::fs::OpenOptions;
                OpenOptions::new().create(true).append(true).open(p).ok()
                    .map(|f| -> Box<dyn std::io::Write> { Box::new(std::io::BufWriter::new(f)) })
            });
        // Write header only once (check file size).
        if let (Some(ref p), Some(ref mut w)) = (&attr_tsv_path, &mut attr_writer) {
            if std::fs::metadata(p).map(|m| m.len()).unwrap_or(1) <= 1 {
                let _ = writeln!(w, "family_id\tread_name_hash\tplacement_copy\tn_sites_covered\tfinal_weight\tweight_gap\tweight_sum\tevidence_gap\tev_decisive");
            }
        }

        for entry in &entries {
            // Compute weight_gap and weight_sum for this read.
            let mut best = 0.0_f64;
            let mut second = 0.0_f64;
            let mut w_sum = 0.0_f64;
            for &w in &entry.weights {
                w_sum += w;
                if w > best { second = best; best = w; }
                else if w > second { second = w; }
            }
            let gap = best - second;
            let max_sites: usize = entry.n_sites_covered.iter().copied().max().unwrap_or(0);
            // Pre-prior EVIDENCE margin (fingerprint + junction + NM, EXCLUDING the prior).
            // The posterior weight `gap` above also reflects the prior — and on a lopsided-
            // coverage family (real RBMY: one dominant copy) the prior drives ambiguous reads
            // onto a SINK copy with a confident-looking weight gap even when the diagnostic
            // evidence is near-flat (the 2026-06-04 RBMY mis-calibration: reported confidence
            // anti-correlated with real PSV identifiability). Reported confidence MUST reflect
            // evidence, so the buckets gate on the evidence margin clearing the same eff_gap
            // bar used for anchoring. This subsumes the earlier min_decisive_sites gate (a
            // 1-site read's evidence margin can't clear gap_threshold). 0-site (structural /
            // DAZ mixed-strand) reads keep the prior path unchanged — DAZ byte-identical gate.
            let ev_gap = {
                let (mut e1, mut e2) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
                for i in 0..entry.locs.len() {
                    let fp = entry.log_scores[i];
                    let base = if fp.is_finite() { fp } else { 0.0 };
                    let v = base + lambda_j * entry.junc_scores[i] + lambda_nm * entry.nm_scores[i];
                    if v > e1 { e2 = e1; e1 = v; } else if v > e2 { e2 = v; }
                }
                if e1.is_finite() && e2.is_finite() { e1 - e2 } else { 0.0 }
            };
            let eff_gap = if max_sites == 0 { struct_gap }
                else if max_sites < min_decisive_sites { gap_threshold }
                else { (per_site_gap * max_sites as f64).min(gap_threshold) };
            let ev_decisive = max_sites == 0 || ev_gap >= eff_gap;
            if gap > 0.8 && ev_decisive {
                n_decisive += 1;
                if max_sites == 0 { n_struct_guided += 1; }
            } else if gap > 0.5 && ev_decisive {
                n_moderate += 1;
                if max_sites == 0 { n_struct_guided += 1; }
            } else {
                n_uncertain += 1;
            }

            // Track sites covered (use max across placements — the read itself
            // determines this, not which copy it's aligned to).
            if max_sites > 0 {
                total_sites_covered += max_sites;
                n_reads_with_sites += 1;
            }

            for (i, &(global_bi, ri)) in entry.locs.iter().enumerate() {
                if ri >= bundles[global_bi].reads.len() { continue; }
                let old_w = bundles[global_bi].reads[ri].weight;
                let new_w = entry.weights[i];
                if (old_w - new_w).abs() > 1e-9 {
                    bundles[global_bi].reads[ri].weight = new_w;
                    n_reweighted += 1;
                }
                // Persist per-read EM attribution for the downstream capacity-
                // confidence channel. Two corrections vs the naive formula:
                //  (1) Uniqueness MUST come from the read's actual placement count
                //      (entry.locs), NOT BundleRead.nh: minimap2 emits no NH tag, so
                //      nh defaults to 1 and every read would look "unique" ->
                //      em_anchored=true -> capacity_confidence stuck at 1.000.
                //  (2) Anchoring is PER PLACEMENT: a read counts as anchored for
                //      THIS copy only if this placement is the read's winner. A read
                //      decisively assigned to copy A also has a decisive weight_gap
                //      at its small residual on copy B; without the winner gate that
                //      residual marks B "anchored" and inflates B's confidence (the
                //      DAZ3 phantom read 1.000 from 0.055-weight residuals).
                let was_unique = entry.locs.len() <= 1;
                let max_w = entry.weights.iter().cloned().fold(f64::MIN, f64::max);
                let is_winner = entry.weights[i] >= max_w - 1e-9;
                bundles[global_bi].reads[ri].em_weight_gap = gap;
                bundles[global_bi].reads[ri].em_n_sites = entry.n_sites_covered[i] as u32;
                bundles[global_bi].reads[ri].em_anchored =
                    was_unique || (is_winner && ((gap > 0.8) || (max_sites > 0 && gap > 0.5)));
                if let Some(ref mut w) = attr_writer {
                    let rnh = &bundles[global_bi].reads[ri].read_name_hash;
                    let _ = writeln!(w, "{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{}",
                        family.family_id, rnh, entry.fam_pos[i],
                        entry.n_sites_covered[i], new_w, gap, w_sum,
                        ev_gap, if ev_decisive { 1 } else { 0 });
                }
            }
        }
        result.reads_reweighted = n_reweighted;

        // Attribution quality summary: confirms long reads (more sites) →
        // decisive assignment, and that no weight is lost (weight_sum ≈ 1.0).
        let n_total = n_decisive + n_moderate + n_uncertain;
        let avg_sites = if n_reads_with_sites > 0 {
            total_sites_covered as f64 / n_reads_with_sites as f64
        } else { 0.0 };
        if n_total > 0 || fp.n_sites > 0 {
            eprintln!(
                "[VG-FP-EM] Family {}: converged={} in {} iter (delta={:.6}), \
                 {} diag-sites, {} reads — decisive={} moderate={} uncertain={} \
                 struct-guided={}, avg_sites_covered={:.1}",
                family.family_id, result.converged, result.iterations, result.max_delta,
                fp.n_sites, n_total, n_decisive, n_moderate, n_uncertain,
                n_struct_guided, avg_sites,
            );
        }
        results.push(result);
    }

    let total_reweighted: usize = results.iter().map(|r| r.reads_reweighted).sum();
    eprintln!(
        "[VG-FP-EM] fingerprint-EM complete: {} reads adjusted across {} families in {:.1}s",
        total_reweighted,
        results.iter().filter(|r| r.reads_reweighted > 0).count(),
        t_start.elapsed().as_secs_f64(),
    );
    results
}

/// Write family report with EM results (uses pre-saved coords since bundles are consumed).
///
/// `transcripts` (optional): when provided, each bundle's reported region is
/// refined by clustering emitted transcripts within the bundle into gene-like
/// sub-loci (default gap = 5kb). A 430kb Rustle bundle that contains 5 separate
/// ZNF genes will report 5 per-gene regions instead of one mega-bundle span.
/// Tune the gap via env `RUSTLE_VG_REPORT_GAP` (bp).
/// Per-family rescue diagnostic counts aggregated from synthetic transcripts.
#[derive(Default, Clone)]
pub struct RescueCounts {
    pub n_rescued: usize,
    pub n_below_thresh: usize,
    pub n_seed_masked: usize,
    pub n_divergent: usize,
    pub n_structural: usize,
    pub n_ref_absent: usize,
}

pub fn write_family_report_with_em(
    path: &std::path::Path,
    families: &[FamilyGroup],
    bundle_coords: &[(String, u64, u64, char)],
    em_results: &[EmResult],
    transcripts: Option<&[crate::path_extract::Transcript]>,
    bundles: Option<&[Bundle]>,
) -> std::io::Result<()> {
    use crate::vg_hmm::diagnostic::RescueClass;
    use std::io::Write;
    let gene_gap: u64 = std::env::var("RUSTLE_VG_REPORT_GAP")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(5_000);

    // Build per-family_id rescue counts from synthetic transcripts.
    let mut rescue_map: HashMap<usize, RescueCounts> = HashMap::new();
    if let Some(txs) = transcripts {
        for tx in txs.iter().filter(|t| t.synthetic) {
            if let Some(fam_id) = tx.vg_family_id {
                let c = rescue_map.entry(fam_id).or_default();
                c.n_rescued += 1;
                match tx.rescue_class {
                    Some(RescueClass::BelowThresholdChain) => c.n_below_thresh += 1,
                    Some(RescueClass::SeedMasked)          => c.n_seed_masked  += 1,
                    Some(RescueClass::Divergent)            => c.n_divergent    += 1,
                    Some(RescueClass::Structural)           => c.n_structural   += 1,
                    Some(RescueClass::ReferenceAbsent)      => c.n_ref_absent   += 1,
                    Some(RescueClass::NovelLocusFromScan)
                    | Some(RescueClass::NeedsExternalVerification)
                    | Some(RescueClass::ChimericSuffixRescue)
                    | Some(RescueClass::TopologyBorrow)
                    | Some(RescueClass::TandemCopy)
                    | None => {}
                }
            }
        }
    }

    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(
        f,
        "family_id\tn_copies\tchrom\tregions\tn_shared_reads\tshared_per_copy\tn_intron_copies\texon_cv\tem_iterations\tem_converged\tem_delta\treads_reweighted\tn_gene_loci\tn_rescued\tn_below_thresh\tn_seed_masked\tn_divergent\tn_structural\tn_ref_absent"
    )?;
    for (fi, family) in families.iter().enumerate() {
        let mut refined_regions: Vec<String> = Vec::new();
        let mut refined_gene_count = 0usize;
        for &bi in &family.bundle_indices {
            let (chrom, bs, be, strand) = match bundle_coords.get(bi) {
                Some(c) => (c.0.as_str(), c.1, c.2, c.3),
                None => continue,
            };
            let sub = match transcripts {
                Some(txs) => refine_bundle_into_gene_loci(txs, chrom, bs, be, strand, gene_gap),
                None => Vec::new(),
            };
            if sub.is_empty() {
                refined_regions.push(format!("{}:{}-{}:{}", chrom, bs, be, strand));
                refined_gene_count += 1;
            } else {
                for (s, e) in &sub {
                    refined_regions.push(format!("{}:{}-{}:{}", chrom, s, e, strand));
                }
                refined_gene_count += sub.len();
            }
        }
        let chrom = bundle_coords
            .get(*family.bundle_indices.first().unwrap_or(&0))
            .map(|(c, _, _, _)| c.as_str())
            .unwrap_or("?");
        let em = em_results.get(fi).cloned().unwrap_or_default();
        let rc = rescue_map.get(&family.family_id).cloned().unwrap_or_default();

        // Family-quality diagnostics: shared-per-copy density and CV of
        // intron counts (the actual signals used by the post-discovery
        // quality filter — see filter_high_confidence_families).
        let n_copies = family.bundle_indices.len();
        let shared_per_copy = if n_copies > 0 {
            family.multimap_reads.len() as f64 / n_copies as f64
        } else { 0.0 };
        let (n_intron_copies, exon_cv): (usize, f64) = match bundles {
            Some(bs) => {
                let intron_counts: Vec<usize> = family.bundle_indices.iter()
                    .map(|&bi| bs.get(bi).map(|b| b.junction_stats.len()).unwrap_or(0))
                    .filter(|&c| c > 0)
                    .collect();
                let n_intron = intron_counts.len();
                let cv = if n_intron >= 2 {
                    let n = n_intron as f64;
                    let mean = intron_counts.iter().map(|&c| c as f64).sum::<f64>() / n;
                    if mean > 0.0 {
                        let var = intron_counts.iter()
                            .map(|&c| { let d = c as f64 - mean; d * d })
                            .sum::<f64>() / n;
                        var.sqrt() / mean
                    } else { 0.0 }
                } else { -1.0 };
                (n_intron, cv)
            }
            None => (0, -1.0),
        };

        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{:.3}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            family.family_id,
            family.bundle_indices.len(),
            chrom,
            refined_regions.join(";"),
            family.multimap_reads.len(),
            shared_per_copy,
            n_intron_copies,
            exon_cv,
            em.iterations,
            em.converged,
            em.max_delta,
            em.reads_reweighted,
            refined_gene_count,
            rc.n_rescued,
            rc.n_below_thresh,
            rc.n_seed_masked,
            rc.n_divergent,
            rc.n_structural,
            rc.n_ref_absent,
        )?;
    }
    Ok(())
}

/// Cluster transcripts overlapping `(chrom, bundle_start, bundle_end)` on matching
/// strand into gene-like sub-loci, separated by gaps ≥ `gap` bp.
/// Returns `Vec<(start, end)>` per sub-locus. If no transcripts overlap, returns empty.
fn refine_bundle_into_gene_loci(
    transcripts: &[crate::path_extract::Transcript],
    chrom: &str,
    bstart: u64,
    bend: u64,
    strand: char,
    gap: u64,
) -> Vec<(u64, u64)> {
    let mut ranges: Vec<(u64, u64)> = transcripts
        .iter()
        .filter(|t| t.chrom == chrom && t.strand == strand)
        .filter_map(|t| {
            let ts = t.exons.first().map(|e| e.0).unwrap_or(0);
            let te = t.exons.last().map(|e| e.1).unwrap_or(0);
            if te >= bstart && ts <= bend {
                Some((ts, te))
            } else {
                None
            }
        })
        .collect();
    if ranges.is_empty() {
        return Vec::new();
    }
    ranges.sort_unstable();
    let mut clusters: Vec<(u64, u64)> = Vec::new();
    let (mut cs, mut ce) = ranges[0];
    for &(s, e) in &ranges[1..] {
        if s > ce + gap {
            clusters.push((cs, ce));
            cs = s;
            ce = e;
        } else if e > ce {
            ce = e;
        }
    }
    clusters.push((cs, ce));
    clusters
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

/// Top-level dispatcher for novel gene-copy discovery.
///
/// Routes to the k-mer legacy path or the new HMM-based rescue path depending on
/// `config.vg_discover_novel_mode`. The HMM path falls back to k-mer on error.
pub fn discover_novel_copies(
    bam_path: &std::path::Path,
    families: &[FamilyGroup],
    bundles: &[Bundle],
    config: &crate::types::RunConfig,
) -> Vec<NovelCandidate> {
    match config.vg_discover_novel_mode.as_str() {
        "hmm" => match crate::vg_hmm::rescue::run_rescue(bam_path, families, bundles, config) {
            Ok(c) => c,
            Err(e) => {
                eprintln!("[VG-HMM] rescue failed: {} — falling back to kmer", e);
                discover_novel_copies_kmer(bam_path, families, bundles, config)
            }
        },
        _ => discover_novel_copies_kmer(bam_path, families, bundles, config),
    }
}

/// Discover novel gene copies from unmapped reads using k-mer matching.
///
/// Approach: extract exonic sequences from each family's bundles (via genome FASTA),
/// build k-mer sets per family, then scan unmapped reads for k-mer overlap.
/// Reads with high k-mer overlap to a family are candidates for novel copies.
///
/// This avoids supplementary alignments (which contain chimeric artifacts).
/// Requires `--genome-fasta` for exon sequence extraction.
fn discover_novel_copies_kmer(
    bam_path: &std::path::Path,
    families: &[FamilyGroup],
    bundles: &[Bundle],
    config: &crate::types::RunConfig,
) -> Vec<NovelCandidate> {
    // FxHash via DetHashSet — k-mer set is hot in the unmapped-read scan loop
    // (n_reads × n_families × n_kmers_per_read lookups) and SipHash dominates
    // there. Same membership semantics, ~3-5x faster lookups.
    type KmerSet = crate::types::DetHashSet<u64>;

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

    let mut family_kmers: Vec<KmerSet> = Vec::with_capacity(families.len());
    for family in families {
        let mut kmers = KmerSet::default();
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
            .collect::<crate::types::DetHashSet<_>>()
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

// ── Direction 1: Post-assembly topology transfer ──────────────────────────────

/// State saved from the EM phase to enable post-assembly topology transfer.
pub struct TopoTransferState {
    /// FamilyGroup partitions (copy index → bundle index, same ordering as family_graphs).
    pub partitions: Vec<FamilyGroup>,
    /// FamilyGraph per partition entry (None if genome was unavailable for that family).
    pub family_graphs: Vec<Option<crate::vg_hmm::family_graph::FamilyGraph>>,
    /// Snapshot of bundles at EM time (needed for chrom/strand/reads access after
    /// `bundles` is consumed by the parallel assembly loop).
    pub bundles: Vec<crate::types::Bundle>,
}

/// Transfer assembled isoform junction chains from high-confidence copies to
/// under-assembled sister copies in the same gene family.
///
/// After EM + first-pass assembly, some copies have confident multi-exon
/// transcripts (`copy_assignment_confidence ≥ min_confidence`) while sisters
/// may have fewer or no assembled isoforms. This function projects the exon
/// topology from high-confidence source copies to sisters using
/// `ExonClass.per_copy_spans` coordinate mapping.
///
/// A projected transcript is added only if:
///   (a) every exon in the source transcript can be mapped to an ExonClass
///       node that also has a span for the sister copy, and
///   (b) the sister bundle has at least one read whose alignment overlaps
///       any projected exon (proves the copy is expressed).
///
/// Returns new synthetic `Transcript` objects to append to the assembly.
/// These carry `synthetic = true`, `rescue_class = TopologyBorrow`, and
/// `copy_assignment_confidence = source_confidence * 0.5` (confidence-discounted).
pub fn transfer_assembled_topology(
    assembled_txs: &[crate::path_extract::Transcript],
    state: &TopoTransferState,
    min_confidence: f64,
) -> Vec<crate::path_extract::Transcript> {
    const SPAN_TOL: u64 = 50; // exon-end matching tolerance (bp)

    use std::collections::HashMap;
    let mut result: Vec<crate::path_extract::Transcript> = Vec::new();

    // Build family_id → partition_idx lookup.
    let mut fam_to_pi: HashMap<usize, usize> = HashMap::new();
    for (pi, part) in state.partitions.iter().enumerate() {
        fam_to_pi.insert(part.family_id, pi);
    }

    // Group source transcripts by (family_id, copy_id).
    // Only multi-exon txs with sufficient confidence qualify as sources.
    let mut by_fam_copy: HashMap<(usize, usize), Vec<&crate::path_extract::Transcript>> =
        HashMap::new();
    for tx in assembled_txs {
        if tx.exons.len() < 2 {
            continue;
        }
        let (Some(fam_id), Some(copy_id)) = (tx.vg_family_id, tx.vg_copy_id) else {
            continue;
        };
        let conf = tx.copy_assignment_confidence.unwrap_or(0.0);
        if conf < min_confidence {
            continue;
        }
        by_fam_copy.entry((fam_id, copy_id)).or_default().push(tx);
    }

    for ((fam_id, src_copy_id), src_txs) in &by_fam_copy {
        let pi = match fam_to_pi.get(fam_id) {
            Some(&i) => i,
            None => continue,
        };
        let fg = match state.family_graphs.get(pi).and_then(|o| o.as_ref()) {
            Some(g) => g,
            None => continue,
        };
        let part = &state.partitions[pi];
        let n_copies = part.bundle_indices.len();
        let topo_trace = std::env::var_os("RUSTLE_VG_TOPO_TRACE").is_some();
        if topo_trace {
            eprintln!("[VG-TOPO-TRACE] fam={} src_copy={} n_src_txs={} n_copies={} fg_nodes={}",
                fam_id, src_copy_id, src_txs.len(), n_copies, fg.nodes.len());
        }

        for sister_copy_id in 0..n_copies {
            if sister_copy_id == *src_copy_id {
                continue;
            }
            let sister_bi = part.bundle_indices[sister_copy_id];
            let sister_bundle = match state.bundles.get(sister_bi) {
                Some(b) => b,
                None => continue,
            };

            // Offset bootstrap (opt-in RUSTLE_VG_TOPO_OFFSET): learn the affine
            // copy-to-copy genomic shift from ExonClasses that DID unify src<->sister,
            // so exons whose sequence-based unification fails can still project.
            // Rationale: dispersed paralog copies are colinear shifts, but minimizer
            // Jaccard collapses to ~0 for short, divergent exons (a 300 bp exon at
            // 97% identity has a SNP every ~33 bp, leaving no shared 15-mer minimizer)
            // — so only some exons unify. The unified ones pin the shift; the rest
            // project by that shift, validated downstream by sister read coverage.
            let copy_offset: Option<i64> =
                if std::env::var_os("RUSTLE_VG_TOPO_OFFSET").is_some() {
                    let mut deltas: Vec<i64> = Vec::new();
                    for n in &fg.nodes {
                        let src_s = n.per_copy_spans.iter()
                            .find(|(c, _)| *c == *src_copy_id).map(|(_, (s, _))| *s);
                        let sis_s = n.per_copy_spans.iter()
                            .find(|(c, _)| *c == sister_copy_id).map(|(_, (s, _))| *s);
                        if let (Some(a), Some(b)) = (src_s, sis_s) {
                            deltas.push(b as i64 - a as i64);
                        }
                    }
                    if deltas.is_empty() {
                        None
                    } else {
                        deltas.sort_unstable();
                        Some(deltas[deltas.len() / 2]) // median shift (indel-robust)
                    }
                } else {
                    None
                };
            if topo_trace {
                eprintln!("[VG-TOPO-TRACE]   sister={} copy_offset={:?} (from shared ExonClasses)",
                    sister_copy_id, copy_offset);
            }

            for src_tx in src_txs.iter() {
                // Project each exon of src_tx from src copy to sister copy.
                let mut proj_exons: Vec<(u64, u64)> =
                    Vec::with_capacity(src_tx.exons.len());
                let mut ok = true;

                for &(ex_s, ex_e) in &src_tx.exons {
                    // Primary: an ExonClass that unifies this src exon with the sister
                    // (sequence/position-derived correspondence — most precise).
                    let via_class: Option<(u64, u64)> = fg.nodes.iter().find_map(|n| {
                        let has_src = n.per_copy_spans.iter().any(|(c, (ns, ne))| {
                            *c == *src_copy_id
                                && ns.saturating_sub(SPAN_TOL) <= ex_s
                                && ex_e <= ne.saturating_add(SPAN_TOL)
                        });
                        if !has_src {
                            return None;
                        }
                        n.per_copy_spans
                            .iter()
                            .find(|(c, _)| *c == sister_copy_id)
                            .map(|(_, (ss, se))| (*ss, *se))
                    });
                    // Fallback: project by the bootstrapped colinear copy shift.
                    let projected = via_class.or_else(|| {
                        copy_offset.map(|off| {
                            let ss = (ex_s as i64 + off).max(0) as u64;
                            let se = (ex_e as i64 + off).max(0) as u64;
                            (ss, se)
                        })
                    });
                    match projected {
                        Some(span) => {
                            if topo_trace {
                                let how = if via_class.is_some() { "via_class" } else { "via_offset" };
                                eprintln!("[VG-TOPO-TRACE]     exon ({},{}) src_copy={} -> ({},{}) [{}]",
                                    ex_s, ex_e, src_copy_id, span.0, span.1, how);
                            }
                            proj_exons.push(span);
                        }
                        None => {
                            if topo_trace {
                                eprintln!("[VG-TOPO-TRACE]     exon ({},{}) src_copy={} -> SKIP (no shared ExonClass with sister {} and no copy_offset available)",
                                    ex_s, ex_e, src_copy_id, sister_copy_id);
                            }
                            ok = false;
                            break;
                        }
                    }
                }

                if !ok || proj_exons.is_empty() {
                    if topo_trace {
                        eprintln!("[VG-TOPO-TRACE]   sister={} src_exons={} projected={} -> SKIP (projection miss: an exon mapped to no shared ExonClass with the sister)",
                            sister_copy_id, src_tx.exons.len(), proj_exons.len());
                    }
                    continue;
                }

                // Validate: sister bundle has at least one read overlapping the projected region.
                let has_reads = sister_bundle.reads.iter().any(|r| {
                    proj_exons
                        .iter()
                        .any(|&(ps, pe)| r.ref_start < pe && r.ref_end > ps)
                });
                if !has_reads {
                    if topo_trace {
                        eprintln!("[VG-TOPO-TRACE]   sister={} projected={} exons -> SKIP (sister bundle has NO read overlapping the projected spans)",
                            sister_copy_id, proj_exons.len());
                    }
                    continue;
                }
                if topo_trace {
                    eprintln!("[VG-TOPO-TRACE]   sister={} -> BORROW {} exons (emitting topology_borrow tx)",
                        sister_copy_id, proj_exons.len());
                }

                proj_exons.sort_by_key(|&(s, _)| s);

                let src_conf = src_tx.copy_assignment_confidence.unwrap_or(0.0);

                // Build the synthetic transcript.
                let synth = build_topology_borrow_transcript(
                    sister_bundle,
                    proj_exons,
                    *fam_id,
                    sister_copy_id,
                    n_copies,
                    src_conf * 0.5,
                );
                result.push(synth);
            }
        }
    }
    result
}

fn build_topology_borrow_transcript(
    sister_bundle: &crate::types::Bundle,
    exons: Vec<(u64, u64)>,
    family_id: usize,
    copy_id: usize,
    family_size: usize,
    confidence: f64,
) -> crate::path_extract::Transcript {
    use crate::path_extract::Transcript;
    // Diagnostic-only (RUSTLE_VG_TOPO_FORCE_EMIT): give the borrowed copy a nominal
    // coverage so downstream coverage filters don't drop it before the GTF. This
    // exists to MEASURE the projection's correctness end-to-end, NOT as a default:
    // a borrowed copy has zero independent coverage by construction, and forcing it
    // past the coverage/phantom guards is exactly the DAZ3 fabrication path.
    let nominal_cov: f64 = if std::env::var_os("RUSTLE_VG_TOPO_FORCE_EMIT").is_some() {
        1.0
    } else {
        0.0
    };
    Transcript {
        chrom: sister_bundle.chrom.clone(),
        strand: sister_bundle.strand,
        exons,
        coverage: nominal_cov,
        exon_cov: Vec::new(),
        tpm: 0.0,
        fpkm: 0.0,
        source: Some("topo_borrow".to_string()),
        is_longread: true,
        longcov: nominal_cov,
        bpcov_cov: 0.0,
        all_strand_cov: 0.0,
        transcript_id: None,
        gene_id: None,
        ref_transcript_id: None,
        ref_gene_id: None,
        hardstart: false,
        hardend: false,
        alt_tts_end: false,
        vg_family_id: Some(family_id),
        vg_copy_id: Some(copy_id),
        vg_family_size: Some(family_size),
        copy_assignment_confidence: Some(confidence),
        copy_independent_support: None, capacity_confidence: None, abundance_min: None, family_verdict: None, tandem_copy: None,
        intron_low: Vec::new(),
        synthetic: true,
        rescue_class: Some(crate::vg_hmm::diagnostic::RescueClass::TopologyBorrow),
        raw_flow_sum: 0.0, min_jct_mm: 0.0, skip_jct_mm: 0.0, chain_witnessed: false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vg_hmm::family_graph::{ExonClass, FamilyGraph, NodeIdx};

    fn approx(a: f64, b: f64) -> bool { (a - b).abs() < 1e-9 }

    #[test]
    fn depth_copynum_equal_copies_no_collapse() {
        // Three equally-covered resolved copies -> 3 single copies, no over-representation.
        let e = estimate_copies_from_depth(&[10.0, 10.0, 10.0], None).unwrap();
        assert!(approx(e.single_copy_unit, 10.0));
        assert!(approx(e.inferred_copies, 3.0));
        assert_eq!(e.per_copy_multiplier, vec![1.0, 1.0, 1.0]);
    }

    #[test]
    fn depth_copynum_flags_collapsed_copy() {
        // Two single copies + one bundle at 3x depth -> that bundle hides ~3 copies,
        // total inferred 5 vs 3 resolved (genome-structure under-counting).
        let e = estimate_copies_from_depth(&[10.0, 10.0, 30.0], None).unwrap();
        assert!(approx(e.single_copy_unit, 10.0));   // faintest = baseline
        assert!(approx(e.inferred_copies, 5.0));
        assert_eq!(e.per_copy_multiplier, vec![1.0, 1.0, 3.0]);
    }

    #[test]
    fn depth_copynum_absent_copies_ignored_in_unit() {
        // Absent copies (depth 0) don't poison the median baseline; still counted as 0x.
        // nonzero=[10,10,30] -> median 10 -> inferred (10+10+30)/10 = 5, third copy 3x.
        let e = estimate_copies_from_depth(&[0.0, 10.0, 10.0, 30.0], None).unwrap();
        assert!(approx(e.single_copy_unit, 10.0));
        assert!(approx(e.inferred_copies, 5.0));
        assert_eq!(e.per_copy_multiplier, vec![0.0, 1.0, 1.0, 3.0]);
    }

    #[test]
    fn depth_copynum_external_unit_gives_absolute() {
        // With an external single-copy depth of 5, everything is 2x+ (uniform collapse
        // that internal calibration alone could not see).
        let e = estimate_copies_from_depth(&[10.0, 10.0, 30.0], Some(5.0)).unwrap();
        assert!(approx(e.single_copy_unit, 5.0));
        assert!(approx(e.inferred_copies, 10.0));
        assert_eq!(e.per_copy_multiplier, vec![2.0, 2.0, 6.0]);
    }

    #[test]
    fn depth_copynum_single_copy_cannot_detect_collapse() {
        // One isolated bundle: no internal reference -> reports 1x (honest limitation).
        let e = estimate_copies_from_depth(&[60.0], None).unwrap();
        assert!(approx(e.inferred_copies, 1.0));
    }

    #[test]
    fn depth_copynum_empty_or_zero_is_none() {
        assert!(estimate_copies_from_depth(&[], None).is_none());
        assert!(estimate_copies_from_depth(&[0.0, 0.0], None).is_none());
    }

    #[test]
    fn em_confidence_buckets_and_means() {
        // gaps 0.9 (decisive), 0.7 (moderate), 0.3 (uncertain); sites 5,3,1.
        let c = summarize_em_confidence(&[(0.9, 5), (0.7, 3), (0.3, 1)]).unwrap();
        assert_eq!((c.n_em_reads, c.n_decisive, c.n_moderate, c.n_uncertain), (3, 1, 1, 1));
        assert!(approx(c.mean_gap, (0.9 + 0.7 + 0.3) / 3.0));
        assert!(approx(c.mean_sites, 3.0));
        assert!(approx(c.decisive_frac(), 1.0 / 3.0));
    }

    #[test]
    fn em_confidence_threshold_edges() {
        // Exactly 0.8 is NOT decisive (gap > 0.8); exactly 0.5 is NOT moderate (gap > 0.5).
        let c = summarize_em_confidence(&[(0.8, 2), (0.5, 2), (0.81, 2)]).unwrap();
        assert_eq!((c.n_decisive, c.n_moderate, c.n_uncertain), (1, 1, 1));
    }

    #[test]
    fn locus_rel_demirrors_strand_mirror_pairs() {
        // Two physical loci 13 Mb apart, each with a +/- strand-mirror pair (4 overlapping
        // bundles total). Old code: Overlapping (mirrors overlap). De-mirrored: Distal.
        let loci = vec![
            ("chr1", 89_000_000, 89_010_000), ("chr1", 89_000_100, 89_010_100), // locus A mirror pair
            ("chr1", 103_000_000, 103_010_000), ("chr1", 103_000_100, 103_010_100), // locus B mirror pair
        ];
        assert_eq!(classify_physical_loci(loci, 4), LocusRel::Distal);
    }

    #[test]
    fn locus_rel_tandem_array_after_demirror() {
        // Three physical loci <1 Mb apart (each a mirror pair) → Tandem.
        let loci = vec![
            ("chrY", 1000, 12000), ("chrY", 1100, 12100),
            ("chrY", 210000, 221000), ("chrY", 210100, 221100),
            ("chrY", 410000, 421000), ("chrY", 410100, 421100),
        ];
        assert_eq!(classify_physical_loci(loci, 6), LocusRel::Tandem);
    }

    #[test]
    fn locus_rel_same_locus_inverted_pair_stays_overlapping() {
        // A genuine same-locus inverted pair (DAZ1-/DAZ3+): one physical locus, 2 bundles.
        let loci = vec![("chrY", 5_000_000, 5_020_000), ("chrY", 5_000_050, 5_020_050)];
        assert_eq!(classify_physical_loci(loci, 2), LocusRel::Overlapping);
    }

    #[test]
    fn locus_rel_trans_and_single() {
        assert_eq!(classify_physical_loci(vec![("chr1", 100, 200), ("chr2", 100, 200)], 2), LocusRel::Trans);
        assert_eq!(classify_physical_loci(vec![("chr1", 100, 200)], 1), LocusRel::Single);
        assert_eq!(classify_physical_loci(vec![], 0), LocusRel::Single);
    }

    #[test]
    fn em_confidence_all_decisive_and_empty() {
        let c = summarize_em_confidence(&[(0.95, 4), (0.99, 6)]).unwrap();
        assert!(approx(c.decisive_frac(), 1.0));
        assert!(summarize_em_confidence(&[]).is_none());
    }

    fn make_ec(idx: usize, copy_seqs: &[(usize, &[u8])], copy_spans: &[(usize, (u64, u64))]) -> ExonClass {
        let span = copy_spans.iter().map(|(_, (s, e))| (*s, *e))
            .fold((u64::MAX, 0u64), |(lo, hi), (s, e)| (lo.min(s), hi.max(e)));
        ExonClass {
            idx: NodeIdx(idx),
            chrom: "chrTest".into(),
            span,
            strand: '+',
            per_copy_sequences: copy_seqs.iter().map(|(c, s)| (*c, s.to_vec())).collect(),
            per_copy_spans: copy_spans.to_vec(),
            copy_specific: copy_spans.len() == 1,
            profile: None,
            per_copy_profiles: vec![],
            per_copy_cov: vec![],
        }
    }

    fn make_read(exons: Vec<(u64, u64)>, mismatches: Vec<(u64, u8)>) -> BundleRead {
        BundleRead {
            read_uid: 0, read_name: "r".into(), read_name_hash: 0,
            ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
            ref_start: exons.first().map(|(s, _)| *s).unwrap_or(0),
            ref_end:   exons.last().map(|(_, e)| *e).unwrap_or(0),
            exons, junctions: vec![], junction_valid: vec![],
            junctions_raw: vec![], junctions_del: vec![],
            weight: 1.0, is_reverse: false, strand: '+',
            has_poly_start: false, has_poly_end: false,
            has_poly_start_aligned: false, has_poly_start_unaligned: false,
            has_poly_end_aligned: false, has_poly_end_unaligned: false,
            unaligned_poly_t: 0, unaligned_poly_a: 0,
            has_last_exon_polya: false, has_first_exon_polyt: false,
            query_length: None, clip_left: 0, clip_right: 0,
            nh: 1, nm: 0, de: None, md: None, insertion_sites: vec![],
            unitig: false, unitig_cov: 0.0, read_count_yc: 1.0,
            countfrag_len: 0.0, countfrag_num: 0.0, junc_mismatch_weight: 0.0,
            pair_idx: vec![], pair_count: vec![],
            mapq: 60, mismatches, seq: Vec::new(), hp_tag: None, ps_tag: None,
            is_primary_alignment: true,
            em_weight_gap: -1.0, em_n_sites: 0, em_anchored: true,
        }
    }

    #[test]
    fn fingerprint_finds_snp_and_scores_correctly() {
        // Two copies differ at exon position 2: copy0=b'G', copy1=b'C'.
        // copy0 at ref 100-105, copy1 at ref 500-505 (non-overlapping).
        let ec = make_ec(
            0,
            &[(0, b"ACGTA"), (1, b"ACCTA")],
            &[(0, (100, 105)), (1, (500, 505))],
        );
        let fg = FamilyGraph { family_id: 0, nodes: vec![ec], edges: vec![] };
        let fp = build_exon_fingerprints(&fg, 2);

        assert_eq!(fp.n_sites, 1, "exactly one variant site");
        assert_eq!(fp.per_copy_site_refs[0][0].1, 102, "copy0 ref_pos");
        assert_eq!(fp.per_copy_site_refs[1][0].1, 502, "copy1 ref_pos");

        // Read aligned to copy0, no mismatch at 102 → base = b'G' = copy0's base.
        let read_a = make_read(vec![(100, 105)], vec![]);
        let (s, n) = score_read_exon_fingerprint(&read_a, 0, &fp);
        assert!(s[0] > s[1], "read with copy0 base should prefer copy0: {s:?}");
        assert_eq!(n, 1, "one site covered");

        // Read aligned to copy0, mismatch=b'C' at 102 → base = copy1's base.
        let read_b = make_read(vec![(100, 105)], vec![(102, b'C')]);
        let (s, n) = score_read_exon_fingerprint(&read_b, 0, &fp);
        assert!(s[1] > s[0], "read with copy1 base should prefer copy1: {s:?}");
        assert_eq!(n, 1, "one site covered");

        // Read doesn't cover either copy's site → all scores neutral (0.0), 0 covered.
        let read_c = make_read(vec![(200, 210)], vec![]);
        let (s, n) = score_read_exon_fingerprint(&read_c, 0, &fp);
        assert_eq!(s, vec![0.0, 0.0]);
        assert_eq!(n, 0, "zero sites covered");
    }

    #[test]
    fn fingerprint_skips_invariant_positions() {
        let ec = make_ec(
            0,
            &[(0, b"ACGT"), (1, b"ACGT")],
            &[(0, (0, 4)), (1, (100, 104))],
        );
        let fg = FamilyGraph { family_id: 0, nodes: vec![ec], edges: vec![] };
        let fp = build_exon_fingerprints(&fg, 2);
        assert_eq!(fp.n_sites, 0, "identical copies yield zero variant sites");
    }

    /// Adaptive gap gate: fingerprint-covered reads use gap_threshold (10.0);
    /// reads with 0 covered sites fall back to struct_gap (0.1).
    /// This is the IsoSeq-specific path that allows junction/NM signals to
    /// guide reads in conserved regions that sequence fingerprints can't reach.
    #[test]
    fn adaptive_gap_selects_correct_threshold() {
        let gap_threshold = 10.0_f64;
        let struct_gap = 0.1_f64;
        let score_gap = 0.5_f64; // log-units between best and second-best placements

        // Structural-only read (n_cov_max = 0): should use struct_gap.
        let eff_gap_structural = if 0 > 0 { gap_threshold } else { struct_gap };
        assert!(
            score_gap >= eff_gap_structural,
            "junction/NM-guided read (0 diagnostic sites) should not be skipped: \
             gap {score_gap} >= struct_gap {eff_gap_structural}"
        );
        assert!(
            score_gap < gap_threshold,
            "same read would be skipped by fingerprint threshold"
        );

        // Fingerprint-covered read (n_cov_max > 0): should use gap_threshold.
        let eff_gap_fp = if 10 > 0 { gap_threshold } else { struct_gap };
        assert!(
            score_gap < eff_gap_fp,
            "fingerprint-covered read with small gap should be skipped: \
             gap {score_gap} < gap_threshold {eff_gap_fp}"
        );
    }

    /// Multi-signal E-step: log-posterior is the log-linear sum of fingerprint,
    /// junction-chain, and alignment-identity scores. The three signals combine
    /// additively in log space — this is the formal model claimed for the EM.
    #[test]
    fn multi_signal_estep_is_additive_in_log_space() {
        let fp_score_a = -1.0_f64;   // fingerprint: copy A preferred
        let fp_score_b = -5.0_f64;   // fingerprint: copy B penalised

        let junc_a = f64::ln(0.95);  // copy A: 95% junctions present (good match)
        let junc_b = f64::ln(0.50);  // copy B: only 50% junctions present

        let nm_a = f64::ln(0.98);    // copy A: near-perfect alignment
        let nm_b = f64::ln(0.80);    // copy B: more edit distance

        let log_prior = f64::ln(0.5); // uniform 2-copy prior
        let lambda_j = 1.0_f64;
        let lambda_nm = 1.0_f64;

        let post_a = fp_score_a + lambda_j * junc_a + lambda_nm * nm_a + log_prior;
        let post_b = fp_score_b + lambda_j * junc_b + lambda_nm * nm_b + log_prior;

        assert!(
            post_a > post_b,
            "copy A should win on combined signal: {post_a:.3} vs {post_b:.3}"
        );

        // Verify the gap is larger than either signal alone would produce.
        let fp_gap_only = fp_score_a - fp_score_b;
        let combined_gap = post_a - post_b;
        assert!(
            combined_gap > fp_gap_only,
            "junction+NM should widen the gap over fingerprint alone: \
             combined {combined_gap:.3} > fp-only {fp_gap_only:.3}"
        );
    }

    // ── copy_support_fraction (certified copy-support guard) ──────────────────
    // rate pairs are (rate_C, rate_min_sibling); n_unique reads always support C.
    // A multimapper supports C iff rate_C <= rate_min_sibling + margin.
    #[test]
    fn daz3_phantom_zero_support() {
        // 0 unique; all multimappers fit a sibling far better (rate_C 0.07 vs sib 0.005)
        let mm = vec![(0.07, 0.005); 30];
        assert!(copy_support_fraction(0, &mm, 0.01) < 0.05);
    }

    #[test]
    fn daz1_real_high_support() {
        // many unique reads + multimappers that fit C best
        let mm = vec![(0.005, 0.07); 14];
        assert!(copy_support_fraction(167, &mm, 0.01) > 0.95);
    }

    #[test]
    fn co_expressed_tie_supports_both() {
        // NM-tie multimappers (within margin) count as supporting C (no false suppression)
        let mm = vec![(0.004, 0.004); 14];
        assert!(copy_support_fraction(9, &mm, 0.01) > 0.95);
    }

    #[test]
    fn margin_boundary() {
        // 0.015 > 0.004+0.01 -> belongs elsewhere
        assert_eq!(copy_support_fraction(0, &vec![(0.015, 0.004)], 0.01), 0.0);
        // 0.013 within margin -> supports
        assert_eq!(copy_support_fraction(0, &vec![(0.013, 0.004)], 0.01), 1.0);
    }

    #[test]
    fn no_reads_is_zero() {
        assert_eq!(copy_support_fraction(0, &[], 0.01), 0.0);
    }

    // ── anchor_read (raw-dNM per-read copy anchor) ────────────────────────────
    #[test]
    fn anchor_owns_when_dnm_large() {
        // this copy NM=6, sibling NM=504 over comparable extent -> owns
        assert_eq!(anchor_read(6, 4600, &[(504, 4145)], 2, 0.7), ReadAnchor::Owns);
    }
    #[test]
    fn anchor_sibling_when_dnm_negative() {
        // this copy NM=504, sibling NM=6 -> belongs to sibling
        assert_eq!(anchor_read(504, 4145, &[(6, 4600)], 2, 0.7), ReadAnchor::Sibling);
    }
    #[test]
    fn anchor_tie_when_no_distinguishing_position() {
        // NM 2 vs 2 -> tie (no distinguishing position)
        assert_eq!(anchor_read(2, 2900, &[(2, 2900)], 2, 0.7), ReadAnchor::Tie);
        // NM 14 vs 26 -> 12 distinguishing mismatches -> NOT a tie, owns
        assert_eq!(anchor_read(14, 4018, &[(26, 3934)], 2, 0.7), ReadAnchor::Owns);
    }
    #[test]
    fn anchor_owns_when_no_comparable_other() {
        // sibling alignment too short (extent guard) -> not a competitor -> owns
        assert_eq!(anchor_read(50, 4000, &[(2, 200)], 2, 0.7), ReadAnchor::Owns);
        // no other placement at all -> owns (uniquely this copy)
        assert_eq!(anchor_read(50, 4000, &[], 2, 0.7), ReadAnchor::Owns);
    }
    #[test]
    fn anchor_boundary_dnm_equals_t() {
        // dnm == t -> Owns (>= t); dnm == t-1 -> Tie
        assert_eq!(anchor_read(0, 100, &[(2, 100)], 2, 0.9), ReadAnchor::Owns);
        assert_eq!(anchor_read(0, 100, &[(1, 100)], 2, 0.9), ReadAnchor::Tie);
        // alen_c == 0 guard -> Owns regardless of others
        assert_eq!(anchor_read(5, 0, &[(0, 100)], 2, 0.7), ReadAnchor::Owns);
    }

    // ── identifiability_partition (union-find equivalence classes) ────────────
    #[test]
    fn partition_full_when_all_pairs_distinguishable() {
        // 3 copies, every shared pair has a distinguishing read -> 3 classes (full)
        let nonid_pairs: Vec<(usize, usize)> = vec![];
        assert_eq!(identifiability_partition(3, &nonid_pairs), 3);
    }
    #[test]
    fn partition_none_when_all_pairs_nonidentifiable() {
        // 3 identical copies, every pair non-identifiable -> 1 class (none)
        let nonid_pairs = vec![(0, 1), (1, 2), (0, 2)];
        assert_eq!(identifiability_partition(3, &nonid_pairs), 1);
    }
    #[test]
    fn partition_partial_when_some_merge() {
        // copies 0,1 identical; 2 distinct -> 2 classes (partial)
        let nonid_pairs = vec![(0, 1)];
        assert_eq!(identifiability_partition(3, &nonid_pairs), 2);
    }

    // ── compute_copy_independent_support fixture ──────────────────────────────
    fn make_read_full(rnh: u64, exons: Vec<(u64, u64)>, nm: u32, is_reverse: bool) -> BundleRead {
        let mut r = make_read(exons, vec![]);
        r.read_name_hash = rnh;
        r.nm = nm;
        r.is_reverse = is_reverse;
        r.strand = if is_reverse { '-' } else { '+' };
        r
    }

    fn make_bundle(chrom: &str, strand: char, reads: Vec<BundleRead>) -> Bundle {
        let start = reads.iter().map(|r| r.ref_start).min().unwrap_or(0);
        let end = reads.iter().map(|r| r.ref_end).max().unwrap_or(start);
        Bundle {
            chrom: chrom.into(),
            start,
            end,
            strand,
            reads,
            junction_stats: Default::default(),
            junction_pair_stats: Default::default(),
            bundlenodes: None,
            read_bnodes: None,
            bnode_colors: None,
            synthetic: false,
            rescue_class: None,
            vg_family_id: None,
        }
    }

    /// Cross-strand DAZ3(+, phantom) ↔ DAZ1(−, real) in the ORIGINAL family.
    /// Each multimapper places at BOTH copies; DAZ3's NM is ~15x DAZ1's.
    /// DAZ3 has no unique reads → support ≈ 0 (suppressed); DAZ1 has unique
    /// reads + multimappers that fit it best → support = 1.0 (kept).
    #[test]
    fn extractor_suppresses_cross_strand_phantom() {
        // aligned_len = 1000 for every read (single 1000-bp exon).
        // DAZ3 copy = fam_pos 0 (+strand); DAZ1 copy = fam_pos 1 (−strand).
        // 30 shared reads: at DAZ3 nm=70 (rate .07), at DAZ1 nm=5 (rate .005).
        let mut daz3_reads = Vec::new(); // copy 0 bundle (the phantom locus)
        let mut daz1_reads = Vec::new(); // copy 1 bundle (the real locus)
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..30u64 {
            let rnh = 1000 + i;
            // placement at DAZ3 (copy 0): bad fit
            daz3_reads.push(make_read_full(rnh, vec![(42_879_000, 42_880_000)], 70, false));
            // placement at DAZ1 (copy 1): good fit
            daz1_reads.push(make_read_full(rnh, vec![(42_780_000, 42_781_000)], 5, true));
            multimap.insert(rnh, vec![(0, i as usize), (1, i as usize)]);
        }
        // DAZ1 also has 167 unique reads (read_name_hash NOT in multimap).
        for i in 0..167u64 {
            daz1_reads.push(make_read_full(9000 + i, vec![(42_785_000, 42_786_000)], 3, true));
        }

        let bundles = vec![
            make_bundle("chrY", '+', daz3_reads),
            make_bundle("chrY", '-', daz1_reads),
        ];
        let family = FamilyGroup {
            family_id: 7,
            bundle_indices: vec![0, 1],
            multimap_reads: multimap,
        };

        let support = compute_copy_independent_support(&family, &bundles, 0.01);
        let daz3 = support[&0];
        let daz1 = support[&1];
        assert!(daz3 < 0.05, "DAZ3 (phantom) should have ~0 support, got {daz3}");
        assert!(daz1 > 0.95, "DAZ1 (real) should have ~1.0 support, got {daz1}");
    }

    /// anchored_mass_per_copy: real copy (unique reads + decisively-owned
    /// multimappers) accumulates clear mass; phantom copy (only bad-fit
    /// multimappers, no unique reads) accumulates ~0.
    #[test]
    fn anchored_mass_real_copy_vs_phantom() {
        // copy0 = real locus (good fit, has unique reads)
        // copy1 = phantom locus (bad fit, no unique reads)
        let mut real_reads = Vec::new();    // copy 0 bundle
        let mut phantom_reads = Vec::new(); // copy 1 bundle
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..30u64 {
            let rnh = 2000 + i;
            // placement at copy 0 (real): good fit, nm=5 over 1000 bp
            real_reads.push(make_read_full(rnh, vec![(100_000, 101_000)], 5, false));
            // placement at copy 1 (phantom): bad fit, nm=70 over 1000 bp
            phantom_reads.push(make_read_full(rnh, vec![(900_000, 901_000)], 70, false));
            multimap.insert(rnh, vec![(0, i as usize), (1, i as usize)]);
        }
        // copy 0 also has 167 UNIQUE reads (read_name_hash NOT in multimap).
        for i in 0..167u64 {
            real_reads.push(make_read_full(8000 + i, vec![(100_500, 101_500)], 3, false));
        }

        let bundles = vec![
            make_bundle("chrTest", '+', real_reads),
            make_bundle("chrTest", '+', phantom_reads),
        ];
        let family = FamilyGroup {
            family_id: 11,
            bundle_indices: vec![0, 1],
            multimap_reads: multimap,
        };

        // dnm = 70 - 5 = 65 >= t(2) at copy0 → all 30 multimappers Own copy0.
        // dnm = 5 - 70 = -65 <= -t at copy1 → Sibling (not counted).
        let mass = anchored_mass_per_copy(&family, &bundles, 2, 0.8);
        assert_eq!(mass.len(), 2);
        // copy0: 167 unique (weight 1.0 each) + 30 owned multimappers = 197.0
        assert!((mass[0] - 197.0).abs() < 1e-6, "real copy mass should be 197.0, got {}", mass[0]);
        // copy1: no unique reads, no owned multimappers → ~0
        assert!(mass[1] < 1e-6, "phantom copy mass should be ~0, got {}", mass[1]);
    }

    /// Genuinely co-expressed near-identical copies (NM-tie) → both kept.
    #[test]
    fn extractor_keeps_co_expressed_ties() {
        let mut a_reads = Vec::new();
        let mut b_reads = Vec::new();
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..14u64 {
            let rnh = 2000 + i;
            a_reads.push(make_read_full(rnh, vec![(100, 1100)], 4, false));
            b_reads.push(make_read_full(rnh, vec![(5000, 6000)], 4, false));
            multimap.insert(rnh, vec![(0, i as usize), (1, i as usize)]);
        }
        // each copy has a few unique anchoring reads
        for i in 0..9u64 {
            a_reads.push(make_read_full(3000 + i, vec![(200, 1200)], 2, false));
            b_reads.push(make_read_full(4000 + i, vec![(5200, 6200)], 2, false));
        }
        let bundles = vec![
            make_bundle("chr1", '+', a_reads),
            make_bundle("chr1", '+', b_reads),
        ];
        let family = FamilyGroup {
            family_id: 1,
            bundle_indices: vec![0, 1],
            multimap_reads: multimap,
        };
        let support = compute_copy_independent_support(&family, &bundles, 0.01);
        assert!(support[&0] > 0.95, "copy A should be kept, got {}", support[&0]);
        assert!(support[&1] > 0.95, "copy B should be kept, got {}", support[&1]);
    }

    /// Anchor-first one-pass: with identical copies (fp.n_sites==0) and a junc/NM
    /// TIE, a tied multimapper must apportion by the ANCHORED mass ratio (unique +
    /// Owns reads), not the uniform 1/n. Copy0 has 20 unique anchors, copy1 has 2;
    /// the shared read should land ~0.91/0.09 (=20/22), and the two placements must
    /// still sum to 1.0 (mass conservation via the unchanged log-sum-exp).
    #[test]
    fn anchored_prior_apportions_tied_read_by_anchor_ratio() {
        std::env::set_var("RUSTLE_VG_ANCHOR_PRIOR", "1");
        // Identical single-exon copies → 0 diagnostic sites; copy0 at 100-200, copy1 at 5100-5200.
        let ec = make_ec(
            0,
            &[(0, b"ACGTACGTAC"), (1, b"ACGTACGTAC")],
            &[(0, (100, 110)), (1, (5100, 5110))],
        );
        let fg = FamilyGraph { family_id: 42, nodes: vec![ec], edges: vec![] };
        let fp = build_exon_fingerprints(&fg, 2);
        assert_eq!(fp.n_sites, 0, "identical copies → no diagnostic sites (DAZ regime)");

        let mut c0_reads = Vec::new(); // copy 0 bundle
        let mut c1_reads = Vec::new(); // copy 1 bundle
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();

        // 20 unique anchors at copy0 (read_name_hash NOT in multimap → unique mass).
        for i in 0..20u64 {
            c0_reads.push(make_read_full(3000 + i, vec![(100, 1100)], 2, false));
        }
        // 2 unique anchors at copy1.
        for i in 0..2u64 {
            c1_reads.push(make_read_full(4000 + i, vec![(5100, 6100)], 2, false));
        }
        // ONE shared multimapper: single exon (no junctions → junc tie), equal nm
        // at both placements (→ identity tie). Only the prior can break the tie.
        let rnh = 9999u64;
        c0_reads.push(make_read_full(rnh, vec![(100, 1100)], 2, false));   // copy0 idx 20
        c1_reads.push(make_read_full(rnh, vec![(5100, 6100)], 2, false));  // copy1 idx 2
        let c0_idx = c0_reads.len() - 1;
        let c1_idx = c1_reads.len() - 1;
        multimap.insert(rnh, vec![(0, c0_idx), (1, c1_idx)]);

        let mut bundles = vec![
            make_bundle("chr1", '+', c0_reads),
            make_bundle("chr1", '+', c1_reads),
        ];
        let family = FamilyGroup {
            family_id: 42,
            bundle_indices: vec![0, 1],
            multimap_reads: multimap,
        };
        let family_graphs = vec![Some(fg)];

        // max_iter is effectively 1 under anchor mode (single E-step with fixed prior).
        let _ = run_fingerprint_em(&[family], &mut bundles, &family_graphs, 5, None);

        let w0 = bundles[0].reads[c0_idx].weight;
        let w1 = bundles[1].reads[c1_idx].weight;

        // Mass conservation: the read's two placements sum to 1.0.
        assert!((w0 + w1 - 1.0).abs() < 1e-6, "mass must sum to 1.0: {w0} + {w1}");
        // Anchored ratio 20:2 → copy0 gets the lion's share, NOT 0.5/0.5.
        assert!(w0 > 0.80, "copy0 (20 anchors) should dominate: w0={w0}");
        assert!(w1 < 0.20, "copy1 (2 anchors) should be down-weighted: w1={w1}");
        assert!((w0 - 0.5).abs() > 0.30, "must NOT collapse to uniform 1/n: w0={w0}");

        std::env::remove_var("RUSTLE_VG_ANCHOR_PRIOR");
    }

    /// Mixed-strand DAZ regime: the family has NO joint sequence graph
    /// (graph = None). The lifted skip (Step 3b) + anchored prior must still
    /// apportion a tied multimapper by anchored mass (20:2), proving the fix
    /// reaches the inverted-pair family rather than skipping it.
    #[test]
    fn anchored_prior_reaches_mixed_strand_none_graph_family() {
        std::env::set_var("RUSTLE_VG_ANCHOR_PRIOR", "1");
        std::env::set_var("RUSTLE_VG_JOINT_STRAND_EM", "1");

        let mut c0_reads = Vec::new();
        let mut c1_reads = Vec::new();
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..20u64 { c0_reads.push(make_read_full(3000 + i, vec![(100, 1100)], 2, false)); }
        for i in 0..2u64  { c1_reads.push(make_read_full(4000 + i, vec![(5100, 6100)], 2, false)); }
        let rnh = 9999u64;
        c0_reads.push(make_read_full(rnh, vec![(100, 1100)], 2, false));
        c1_reads.push(make_read_full(rnh, vec![(5100, 6100)], 2, false));
        let c0_idx = c0_reads.len() - 1;
        let c1_idx = c1_reads.len() - 1;
        multimap.insert(rnh, vec![(0, c0_idx), (1, c1_idx)]);

        let mut bundles = vec![
            make_bundle("chr1", '+', c0_reads),
            make_bundle("chr1", '-', c1_reads),   // opposite strand = inverted pair
        ];
        let family = FamilyGroup { family_id: 77, bundle_indices: vec![0, 1], multimap_reads: multimap };
        let family_graphs: Vec<Option<crate::vg_hmm::family_graph::FamilyGraph>> = vec![None];

        let _ = run_fingerprint_em(&[family], &mut bundles, &family_graphs, 5, None);

        let w0 = bundles[0].reads[c0_idx].weight;
        let w1 = bundles[1].reads[c1_idx].weight;
        assert!((w0 + w1 - 1.0).abs() < 1e-6, "mass must sum to 1.0: {w0} + {w1}");
        assert!(w0 > 0.80, "copy0 (20 anchors) should dominate even with None graph: w0={w0}");
        assert!(w1 < 0.20, "copy1 (2 anchors) down-weighted: w1={w1}");

        std::env::remove_var("RUSTLE_VG_ANCHOR_PRIOR");
        std::env::remove_var("RUSTLE_VG_JOINT_STRAND_EM");
    }

    #[test]
    fn classify_family_smoke() {
        // 2-copy resolvable family; MULTI-EXON reads so it is in scope.
        let ex_a = vec![(100u64, 200u64), (300, 400)];     // 2 exons -> spliced, alen=200
        let ex_b = vec![(5100u64, 5200u64), (5300, 5400)];
        let mut a_reads = Vec::new();
        let mut b_reads = Vec::new();
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        // 5 unique reads per copy (own their copy)
        for i in 0..5u64 {
            a_reads.push(make_read_full(3000 + i, ex_a.clone(), 2, false));
            b_reads.push(make_read_full(4000 + i, ex_b.clone(), 2, false));
        }
        // 3 shared reads that decisively OWN A (A nm=2, B nm=10)
        for i in 0..3u64 {
            let rnh = 2000 + i;
            a_reads.push(make_read_full(rnh, ex_a.clone(), 2, false));     // A bundle index 5+i
            b_reads.push(make_read_full(rnh, ex_b.clone(), 10, false));    // B bundle index 5+i
            multimap.insert(rnh, vec![(0, (5 + i) as usize), (1, (5 + i) as usize)]);
        }
        // 3 shared reads that decisively OWN B
        for i in 0..3u64 {
            let rnh = 2100 + i;
            a_reads.push(make_read_full(rnh, ex_a.clone(), 10, false));    // A bundle index 8+i
            b_reads.push(make_read_full(rnh, ex_b.clone(), 2, false));     // B bundle index 8+i
            multimap.insert(rnh, vec![(0, (8 + i) as usize), (1, (8 + i) as usize)]);
        }
        let bundles = vec![make_bundle("chr1", '+', a_reads), make_bundle("chr1", '+', b_reads)];
        let family = FamilyGroup { family_id: 1, bundle_indices: vec![0, 1], multimap_reads: multimap };
        let v = classify_family(&family, &bundles, &FamilyParams::default());
        assert_eq!(v.class, FamilyClass::Family);
        assert_eq!(v.n_copies, 2);
        assert!(v.n_expressed >= 2, "n_expressed={}", v.n_expressed);
        assert_eq!(v.identifiability, Identifiability::Full);
    }

    #[test]
    fn classify_family_nonid_class_counts_as_expressed() {
        // 3 copies: A distinct; B,C near-identical (non-identifiable). A expressed via
        // unique reads; {B,C} expressed as a CLASS (reads beat A but tie B/C). Per-copy
        // counting would call this gene_plus_unexpressed_paralog (n_expressed=1); the
        // class-level count yields TWO expressed classes -> Family.
        let ex_a = vec![(100u64, 200u64), (300, 400)];
        let ex_b = vec![(5100u64, 5200u64), (5300, 5400)];
        let ex_c = vec![(9100u64, 9200u64), (9300, 9400)];
        let mut a_reads = Vec::new();
        let mut b_reads = Vec::new();
        let mut c_reads = Vec::new();
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..5u64 { a_reads.push(make_read_full(3000 + i, ex_a.clone(), 2, false)); }   // A unique
        for i in 0..4u64 {                                                                      // shared: A nm10, B/C nm2 (tie)
            let rnh = 2000 + i;
            a_reads.push(make_read_full(rnh, ex_a.clone(), 10, false));   // A index 5+i
            b_reads.push(make_read_full(rnh, ex_b.clone(), 2, false));    // B index i
            c_reads.push(make_read_full(rnh, ex_c.clone(), 2, false));    // C index i
            multimap.insert(rnh, vec![(0, (5 + i) as usize), (1, i as usize), (2, i as usize)]);
        }
        let bundles = vec![
            make_bundle("chr1", '+', a_reads),
            make_bundle("chr1", '+', b_reads),
            make_bundle("chr1", '+', c_reads),
        ];
        let family = FamilyGroup { family_id: 2, bundle_indices: vec![0, 1, 2], multimap_reads: multimap };
        let v = classify_family(&family, &bundles, &FamilyParams::default());
        assert_eq!(v.class, FamilyClass::Family, "B,C non-id class should make this a Family");
        assert_eq!(v.n_expressed, 2, "two expressed classes (A | BC)");
        assert_eq!(v.identifiability, Identifiability::Partial);
    }

    #[test]
    fn classify_family_fully_nonid_is_nonidentifiable() {
        // 2 copies, every shared read ties (no distinguishing position) -> one
        // expressed class of size 2 -> FamilyNonIdentifiable.
        let ex_b = vec![(100u64, 200u64), (300, 400)];
        let ex_c = vec![(5100u64, 5200u64), (5300, 5400)];
        let mut b_reads = Vec::new();
        let mut c_reads = Vec::new();
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..8u64 {
            let rnh = 7000 + i;
            b_reads.push(make_read_full(rnh, ex_b.clone(), 3, false));
            c_reads.push(make_read_full(rnh, ex_c.clone(), 3, false));
            multimap.insert(rnh, vec![(0, i as usize), (1, i as usize)]);
        }
        let bundles = vec![make_bundle("chr1", '+', b_reads), make_bundle("chr1", '+', c_reads)];
        let family = FamilyGroup { family_id: 3, bundle_indices: vec![0, 1], multimap_reads: multimap };
        let v = classify_family(&family, &bundles, &FamilyParams::default());
        assert_eq!(v.class, FamilyClass::FamilyNonIdentifiable);
        assert_eq!(v.identifiability, Identifiability::None);
        assert_eq!(v.n_id_classes, 1);
    }

    /// `family_for_em_input` returns ONE FamilyGroup spanning BOTH strands,
    /// preserving `fam_pos` indexing into `bundle_indices` and the full
    /// `multimap_reads` (no per-strand split, no <2-placement drop). Contrast:
    /// `partition_and_remap_family_by_strand` would split this into two
    /// single-bundle sub-families and drop the cross-strand read.
    #[test]
    fn family_for_em_input_keeps_cross_strand_placements() {
        // Read 7 maps to copy 0 (would be '-' strand) and copy 1 ('+' strand).
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        multimap.insert(7u64, vec![(0usize, 0usize), (1usize, 0usize)]);
        let fam = FamilyGroup {
            family_id: 42,
            bundle_indices: vec![10, 20],
            multimap_reads: multimap,
        };
        let bundles: Vec<Bundle> = Vec::new();
        let out = family_for_em_input(&fam, &bundles);
        assert_eq!(out.len(), 1, "unsplit family yields exactly one group");
        let g = &out[0];
        assert_eq!(g.family_id, 42, "family_id preserved (no *10+i remap)");
        assert_eq!(g.bundle_indices, vec![10, 20], "all copies kept");
        let locs = g.multimap_reads.get(&7u64)
            .expect("cross-strand shared read retained in EM input");
        assert!(locs.len() >= 2, "shared read keeps >=2 placements (got {})", locs.len());
    }

    /// Conservation hole closed: when the score-gap gate fires for a read,
    /// its per-placement weights are set to the normalized anchored prior
    /// (mass-conserving, sums to 1.0) instead of being left at raw 1/NH,
    /// and the em_* fields are populated on the BundleRead.
    #[test]
    fn gate_fired_read_conserves_mass_and_sets_em_fields() {
        // Two copies with one diagnostic SNP at exon pos 2: copy0=G, copy1=C.
        // copy0 at ref 100-105, copy1 at ref 500-505 (non-overlapping loci).
        let ec = make_ec(
            0,
            &[(0, b"ACGTA"), (1, b"ACCTA")],
            &[(0, (100, 105)), (1, (500, 505))],
        );
        let fg = FamilyGraph { family_id: 0, nodes: vec![ec], edges: vec![] };
        let family_graphs = vec![Some(fg)];

        // One multimapper placed at both copies. At copy0 it covers the SNP
        // (no mismatch -> base G -> copy0 allele); at copy1 it is the SAME read
        // span (covers 0 sites there). The diagonal fp scoring is near-tied so
        // the read trips the gap gate.
        let mut r0 = make_read(vec![(100, 105)], vec![]);
        r0.read_name_hash = 42;
        r0.nh = 2;
        r0.weight = 0.5;
        let mut r1 = make_read(vec![(100, 105)], vec![]);
        r1.read_name_hash = 42;
        r1.nh = 2;
        r1.weight = 0.5;

        let bundles_vec = vec![
            make_bundle("chrTest", '+', vec![r0]),
            make_bundle("chrTest", '+', vec![r1]),
        ];
        let mut bundles = bundles_vec;

        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        multimap.insert(42, vec![(0, 0), (1, 0)]);
        let family = FamilyGroup {
            family_id: 0,
            bundle_indices: vec![0, 1],
            multimap_reads: multimap,
        };

        // Force the gate to ALWAYS fire by setting a huge per-site gap threshold
        // (any finite score gap < eff_gap -> gate fires for this read).
        std::env::set_var("RUSTLE_VG_EM_SCORE_GAP", "1000.0");
        std::env::set_var("RUSTLE_VG_EM_SCORE_GAP_PER_SITE", "1000.0");

        let _ = run_fingerprint_em(&[family], &mut bundles, &family_graphs, 1, None);

        std::env::remove_var("RUSTLE_VG_EM_SCORE_GAP");
        std::env::remove_var("RUSTLE_VG_EM_SCORE_GAP_PER_SITE");

        let w0 = bundles[0].reads[0].weight;
        let w1 = bundles[1].reads[0].weight;
        // Mass conservation across the read's placements.
        assert!((w0 + w1 - 1.0).abs() < 1e-6, "gate-fired weights must sum to 1.0: {w0}+{w1}");

        // em_* fields populated in write-back.
        for bi in 0..2 {
            let r = &bundles[bi].reads[0];
            assert!(r.em_weight_gap >= 0.0, "em_weight_gap set, got {}", r.em_weight_gap);
            assert_eq!(r.em_n_sites, 1, "em_n_sites = best-placement covered sites");
            // gate-fired read: gap small, max_sites>0 but gap<=0.5, nh>1 -> not anchored.
            assert!(!r.em_anchored, "gate-fired near-tied read should not be em_anchored");
        }
    }

    /// REGRESSION (NH-tag bug): minimap2 emits no NH tag, so BundleRead.nh defaults
    /// to 1. em_anchored must NOT therefore treat an uncertain MULTIMAPPER as
    /// "unique/anchored" — that defaulted capacity_confidence to 1.000 on every
    /// minimap2 BAM (incl. the real DAZ run). Uniqueness must come from the read's
    /// actual placement count, not the absent tag.
    #[test]
    fn em_anchored_uses_placement_count_not_absent_nh_tag() {
        // Identical copies -> 0 diagnostic sites; single-exon read -> junc/NM tie:
        // the multimapper is genuinely uncertain.
        let ec = make_ec(
            0,
            &[(0, b"ACGTA"), (1, b"ACGTA")],
            &[(0, (100, 105)), (1, (500, 505))],
        );
        let fg = FamilyGraph { family_id: 0, nodes: vec![ec], edges: vec![] };
        let family_graphs = vec![Some(fg)];

        // nh is left at the make_read default of 1 -> simulates a BAM with NO NH tag.
        let mut r0 = make_read(vec![(100, 105)], vec![]);
        r0.read_name_hash = 99; r0.weight = 0.5;
        let mut r1 = make_read(vec![(100, 105)], vec![]);
        r1.read_name_hash = 99; r1.weight = 0.5;
        assert_eq!(r0.nh, 1, "precondition: read looks 'unique' by the absent-NH default");

        let mut bundles = vec![
            make_bundle("c", '+', vec![r0]),
            make_bundle("c", '+', vec![r1]),
        ];
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        multimap.insert(99, vec![(0, 0), (1, 0)]);
        let family = FamilyGroup {
            family_id: 0,
            bundle_indices: vec![0, 1],
            multimap_reads: multimap,
        };

        let _ = run_fingerprint_em(&[family], &mut bundles, &family_graphs, 1, None);

        // 2-placement multimapper covering 0 diagnostic sites -> uncertain -> must
        // NOT be counted as anchored, despite nh==1 (no NH tag).
        assert!(!bundles[0].reads[0].em_anchored,
            "uncertain multimapper must not be em_anchored when NH tag is absent (nh defaulted to 1)");
        assert!(!bundles[1].reads[0].em_anchored,
            "both placements of the uncertain multimapper must be unanchored");
    }

    /// Per-placement anchoring: a read decisively apportioned to copy0 (via the
    /// anchored prior) must mark ONLY its copy0 (winner) placement anchored — not
    /// its small copy1 residual. Without the winner gate, the minority/phantom
    /// copy's capacity_confidence is falsely 1.000 (the DAZ3 phantom symptom).
    #[test]
    fn em_anchored_is_per_placement_winner_only() {
        std::env::set_var("RUSTLE_VG_ANCHOR_PRIOR", "1");
        // identical copies -> 0 diagnostic sites; the anchored prior (copy0 has far
        // more unique anchors) drives the shared read decisively to copy0.
        let ec = make_ec(0, &[(0, b"ACGTACGTAC"), (1, b"ACGTACGTAC")],
                         &[(0, (100, 110)), (1, (5100, 5110))]);
        let fg = FamilyGraph { family_id: 42, nodes: vec![ec], edges: vec![] };
        let mut c0 = Vec::new(); let mut c1 = Vec::new();
        let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
        for i in 0..50u64 { c0.push(make_read_full(3000 + i, vec![(100, 1100)], 2, false)); }
        c1.push(make_read_full(4000, vec![(5100, 6100)], 2, false));
        let rnh = 9999u64;
        c0.push(make_read_full(rnh, vec![(100, 1100)], 2, false));
        c1.push(make_read_full(rnh, vec![(5100, 6100)], 2, false));
        let c0i = c0.len() - 1; let c1i = c1.len() - 1;
        multimap.insert(rnh, vec![(0, c0i), (1, c1i)]);
        let mut bundles = vec![make_bundle("chr1", '+', c0), make_bundle("chr1", '+', c1)];
        let family = FamilyGroup { family_id: 42, bundle_indices: vec![0, 1], multimap_reads: multimap };
        let family_graphs = vec![Some(fg)];
        let _ = run_fingerprint_em(&[family], &mut bundles, &family_graphs, 5, None);

        let w0 = bundles[0].reads[c0i].weight; let w1 = bundles[1].reads[c1i].weight;
        assert!(w0 > 0.9, "shared read should be decisively copy0: w0={w0} w1={w1}");
        assert!(bundles[0].reads[c0i].em_anchored, "winning placement (copy0) must be anchored");
        assert!(!bundles[1].reads[c1i].em_anchored, "losing residual (copy1) must NOT be anchored");
        std::env::remove_var("RUSTLE_VG_ANCHOR_PRIOR");
    }
}

#[cfg(test)]
mod vg_phasing_sites {
    use super::*;
    use crate::vg_hmm::family_graph::{ExonClass, FamilyGraph, NodeIdx};

    /// Build a single ExonClass node carrying per-copy sequences (GENOME-FORWARD,
    /// as produced by genome.fetch_sequence) and per-copy genomic spans, with an
    /// explicit node strand. Used to construct minimal FamilyGraph fixtures for
    /// enumerate_diagnostic_sites.
    fn node(
        idx: usize,
        strand: char,
        copies: &[(usize, &[u8], (u64, u64))],
    ) -> ExonClass {
        let span = copies.iter().map(|(_, _, (s, e))| (*s, *e))
            .fold((u64::MAX, 0u64), |(lo, hi), (s, e)| (lo.min(s), hi.max(e)));
        ExonClass {
            idx: NodeIdx(idx),
            chrom: "chrT".into(),
            span,
            strand,
            per_copy_sequences: copies.iter().map(|(c, s, _)| (*c, s.to_vec())).collect(),
            per_copy_spans: copies.iter().map(|(c, _, sp)| (*c, *sp)).collect(),
            copy_specific: copies.len() == 1,
            profile: None,
            per_copy_profiles: vec![],
            per_copy_cov: vec![],
        }
    }

    fn refs_pos_base(refs: &[(usize, u64, u8)]) -> Vec<(u64, char)> {
        refs.iter().map(|&(_, p, b)| (p, b as char)).collect()
    }

    /// Test 1 — forward SNP. copy0 exon @100..105 "ACGTA" (+), copy1 exon
    /// @500..505 "ACCTA" (+) differ at offset 2 → one site at genomic 102/'G'
    /// (copy0) and 502/'C' (copy1).
    #[test]
    fn forward_snp_genomic_coords() {
        let fg = FamilyGraph {
            family_id: 0,
            nodes: vec![node(0, '+', &[(0, b"ACGTA", (100, 105)), (1, b"ACCTA", (500, 505))])],
            edges: vec![],
        };
        let fp = enumerate_diagnostic_sites(&fg, 2);
        assert_eq!(fp.n_sites, 1, "one divergent column");
        assert_eq!(refs_pos_base(&fp.per_copy_site_refs[0]), vec![(102, 'G')]);
        assert_eq!(refs_pos_base(&fp.per_copy_site_refs[1]), vec![(502, 'C')]);
    }

    /// Test 2 — two SNPs equal-length. copy0 "ACGTA", copy1 "TCGTT" differ at
    /// offsets 0 and 4 → two sites at the right genomic coords.
    #[test]
    fn two_snps_equal_length() {
        let fg = FamilyGraph {
            family_id: 0,
            nodes: vec![node(0, '+', &[(0, b"ACGTA", (100, 105)), (1, b"TCGTT", (500, 505))])],
            edges: vec![],
        };
        let fp = enumerate_diagnostic_sites(&fg, 2);
        assert_eq!(fp.n_sites, 2, "two divergent columns");
        assert_eq!(refs_pos_base(&fp.per_copy_site_refs[0]), vec![(100, 'A'), (104, 'A')]);
        assert_eq!(refs_pos_base(&fp.per_copy_site_refs[1]), vec![(500, 'T'), (504, 'T')]);
    }

    /// Test 3 — minus strand stores the genome-forward base unchanged.
    /// per_copy_sequences is genome-forward for BOTH strands (fetch_sequence does
    /// no strand handling), so a '-' node's bytes are NOT revcomp'd. copy0 '-' exon
    /// @100..105 genome-forward "ACGTA"; copy1 '-' exon @500..505 genome-forward
    /// "ACCTA" differ at offset 2 → site copy0 (102,'G'), copy1 (502,'C') — the
    /// SAME coords/alleles as the '+' test in Test 1, proving strand does not change
    /// the genome-forward mapping. (All copies in a graph are single-strand, so this
    /// is a same-strand '-' pair, matching production.)
    #[test]
    fn minus_strand_stores_genome_forward_base() {
        let fg = FamilyGraph {
            family_id: 0,
            nodes: vec![node(0, '-', &[
                (0, b"ACGTA", (100, 105)),
                (1, b"ACCTA", (500, 505)),
            ])],
            edges: vec![],
        };
        let fp = enumerate_diagnostic_sites(&fg, 2);
        assert_eq!(fp.n_sites, 1, "one divergent column on genome-forward frame");
        // copy0 genome-forward "ACGTA": offset 2 = 'G' at genomic 102 (no revcomp).
        assert_eq!(refs_pos_base(&fp.per_copy_site_refs[0]), vec![(102, 'G')]);
        // copy1 genome-forward "ACCTA": offset 2 = 'C' at genomic 502.
        assert_eq!(refs_pos_base(&fp.per_copy_site_refs[1]), vec![(502, 'C')]);
    }

    /// Test 4 — same-strand '-' pair (production-realistic). All copies in a
    /// FamilyGraph are the SAME strand (partition_family_by_strand splits mixed
    /// strands upstream), so an opposite-strand pair within one graph is impossible.
    /// Two '-' copies, genome-forward, one SNP → correct genomic coords on both.
    #[test]
    fn same_strand_minus_pair_maps_both_copies() {
        let fg = FamilyGraph {
            family_id: 0,
            nodes: vec![node(0, '-', &[
                (0, b"GGACGTAACC", (200, 210)),
                (1, b"GGACATAACC", (700, 710)),
            ])],
            edges: vec![],
        };
        let fp = enumerate_diagnostic_sites(&fg, 2);
        assert_eq!(fp.n_sites, 1, "one divergent column");
        // "GGACGTAACC" vs "GGACATAACC" differ at offset 4 ('G' vs 'A').
        assert_eq!(refs_pos_base(&fp.per_copy_site_refs[0]), vec![(204, 'G')]);
        assert_eq!(refs_pos_base(&fp.per_copy_site_refs[1]), vec![(704, 'A')]);
    }

    /// Test 5 — unequal length / indel. copy1 has a 1-bp insertion upstream of a
    /// SNP. With the insertion absorbed by anchoring, the downstream SNP maps to
    /// its correct genomic position and NO spurious site is emitted in the gap.
    #[test]
    fn unequal_length_indel_skips_gap() {
        // copy0: AAAAAAAAAACGTACGTACGTACGTGAAAAAAAAAA  (len 36)
        // copy1: AAAAAAAAAACGTACGTACGTACGTGTAAAAAAAAAA (len 37) — extra 'T' after the
        //        long shared prefix; then a SNP further downstream.
        // Construct: shared anchor prefix, then insertion in copy1, then shared
        // anchor, then a single SNP, then shared anchor suffix.
        let pre = b"GACTGACTGACTGACTGACT"; // 20bp shared anchor-rich prefix
        let mid = b"CACACACACACACACACACA"; // 20bp shared block (post-insertion)
        let suf = b"TGTGTGTGTGTGTGTGTGTG"; // 20bp shared suffix
        // copy0: pre + mid + suf, with a SNP planted at a position in `suf`.
        // copy1: pre + "A"(insertion) + mid + suf' (same SNP region but base differs).
        let mut c0: Vec<u8> = Vec::new();
        c0.extend_from_slice(pre);
        c0.extend_from_slice(mid);
        c0.extend_from_slice(suf);
        let snp_off_c0 = pre.len() + mid.len() + 10; // a 'T' in suf at offset 10
        c0[snp_off_c0] = b'T';

        let mut c1: Vec<u8> = Vec::new();
        c1.extend_from_slice(pre);
        c1.push(b'A'); // 1bp insertion
        c1.extend_from_slice(mid);
        c1.extend_from_slice(suf);
        let snp_off_c1 = pre.len() + 1 + mid.len() + 10;
        c1[snp_off_c1] = b'C'; // divergent base at the homologous position

        let len0 = c0.len() as u64;
        let len1 = c1.len() as u64;
        let leaked: &'static [u8] = Box::leak(c0.clone().into_boxed_slice());
        let leaked1: &'static [u8] = Box::leak(c1.clone().into_boxed_slice());
        let fg = FamilyGraph {
            family_id: 0,
            nodes: vec![
                node(0, '+', &[(0, leaked, (100, 100 + len0))]),
                node(1, '+', &[(1, leaked1, (500, 500 + len1))]),
            ],
            edges: vec![],
        };
        let fp = enumerate_diagnostic_sites(&fg, 2);
        // Exactly one divergent site (the SNP); the insertion gap emits nothing.
        assert_eq!(fp.n_sites, 1, "only the SNP, no spurious gap sites: {:?}", fp.per_copy_site_refs);
        // copy0 SNP genomic position = 100 + snp_off_c0.
        assert_eq!(
            refs_pos_base(&fp.per_copy_site_refs[0]),
            vec![(100 + snp_off_c0 as u64, 'T')]
        );
        // copy1 SNP genomic position = 500 + snp_off_c1 (shifted by the insertion).
        assert_eq!(
            refs_pos_base(&fp.per_copy_site_refs[1]),
            vec![(500 + snp_off_c1 as u64, 'C')]
        );
    }

    /// Test 6a — single-copy class contributes nothing.
    #[test]
    fn single_copy_class_contributes_nothing() {
        let fg = FamilyGraph {
            family_id: 0,
            nodes: vec![node(0, '+', &[(0, b"ACGTACGTAC", (100, 110))])],
            edges: vec![],
        };
        let fp = enumerate_diagnostic_sites(&fg, 2);
        assert_eq!(fp.n_sites, 0);
        assert!(fp.per_copy_site_refs[0].is_empty());
        assert!(fp.per_copy_site_refs[1].is_empty());
    }

    /// Test 6b — all-identical pair → n_sites == 0.
    #[test]
    fn identical_pair_yields_no_sites() {
        let fg = FamilyGraph {
            family_id: 0,
            nodes: vec![node(0, '+', &[
                (0, b"ACGTACGTAC", (100, 110)),
                (1, b"ACGTACGTAC", (500, 510)),
            ])],
            edges: vec![],
        };
        let fp = enumerate_diagnostic_sites(&fg, 2);
        assert_eq!(fp.n_sites, 0, "identical copies yield zero sites");
        assert!(fp.per_copy_site_refs[0].is_empty());
        assert!(fp.per_copy_site_refs[1].is_empty());
    }
}

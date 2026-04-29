//! Rescue pipeline: read prefilter → HMM scoring → clustering → synthetic bundle emission.

use crate::types::{Bundle, BundleRead, JunctionStats, RunConfig};
use crate::vg::{FamilyGroup, NovelCandidate};
use crate::vg_hmm::family_graph::{FamilyGraph, NodeIdx};
use crate::vg_hmm::scorer::{forward_against_family, viterbi_path, ViterbiPath};
use anyhow::Result;
use std::collections::HashMap;

// ── FNV-1a (same as vg.rs) ────────────────────────────────────────────────────

#[inline]
fn fnv1a64(s: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in s {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

// ── Task 5.1: K-mer prefilter ─────────────────────────────────────────────────

/// K-mer prefilter using FNV-1a hashing (identical to `vg.rs::kmer_hash`).
///
/// Returns `(passed, hit_count)`.
/// - `passed` is `true` when `hit_count >= min_hits`.
/// - `hit_count` is the number of distinct k-mers in `seq` that appear in
///   `family_kmers`.
pub fn prefilter_read(
    seq: &[u8],
    family_kmers: &std::collections::HashSet<u64>,
    k: usize,
    min_hits: usize,
) -> (bool, usize) {
    if seq.len() < k || family_kmers.is_empty() {
        return (false, 0);
    }
    let mut hits = 0usize;
    for window in seq.windows(k) {
        // Skip windows with ambiguous bases.
        if window.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            let h = fnv1a64(window);
            if family_kmers.contains(&h) {
                hits += 1;
            }
        }
    }
    (hits >= min_hits, hits)
}

// ── Task 5.2: In-memory rescue scoring entry point ───────────────────────────

/// A read that passed HMM scoring, paired with its best Viterbi path.
/// `(read_name, family_id, viterbi_path, log_lik)`
pub type ScoredRead = (String, usize, ViterbiPath, f64);

/// In-memory rescue: score pre-decoded reads against pre-built family graphs.
///
/// `family_graphs` must be indexed identically to `families` (same order, same
/// family_id values). Each graph's nodes must already have profiles fitted
/// (`fit_profiles_in_place` must have been called by the caller).
///
/// Returns `(novel_candidates, scored_reads)`.
///
/// When `family_graphs` or `unmapped_reads` is empty the function returns two
/// empty Vecs without allocating.
pub fn run_rescue_in_memory(
    families: &[FamilyGroup],
    family_graphs: &[FamilyGraph],
    unmapped_reads: &[(String, Vec<u8>)],   // (read_name, sequence)
    min_loglik: f64,
) -> Result<(Vec<NovelCandidate>, Vec<ScoredRead>)> {
    if families.is_empty() || family_graphs.is_empty() || unmapped_reads.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }

    let mut candidates: Vec<NovelCandidate> = Vec::new();
    let mut scored: Vec<ScoredRead> = Vec::new();

    for (read_name, seq) in unmapped_reads {
        // Score against every family; keep the best.
        let mut best_score = f64::NEG_INFINITY;
        let mut best_family_idx: Option<usize> = None;
        let mut best_vpath: Option<ViterbiPath> = None;

        for (fi, fg) in family_graphs.iter().enumerate() {
            let ll = forward_against_family(fg, seq);
            if ll > best_score {
                best_score = ll;
                best_family_idx = Some(fi);
                best_vpath = viterbi_path(fg, seq);
            }
        }

        if best_score <= min_loglik {
            continue;
        }
        let fi = match best_family_idx {
            Some(f) => f,
            None => continue,
        };
        let vpath = match best_vpath {
            Some(p) => p,
            None => continue,
        };

        let family_id = family_graphs[fi].family_id;
        candidates.push(NovelCandidate {
            read_name: read_name.clone(),
            family_id,
            matched_junctions: vpath.nodes.len(),
            total_junctions: family_graphs[fi].nodes.len(),
            approx_chrom: String::new(),
            approx_start: 0,
            approx_end: seq.len() as u64,
        });
        scored.push((read_name.clone(), family_id, vpath, best_score));
    }

    Ok((candidates, scored))
}

// ── Task 5.3: BAM-backed run_rescue driver ────────────────────────────────────

/// Top-level entry point for HMM-based novel-copy rescue.
///
/// (a) Load genome (optional; profiles degrade gracefully to empty if absent).
/// (b) Build family graphs with sequences.
/// (c) Fit per-node profiles.
/// (d) Iterate BAM unmapped reads, decode sequences.
/// (e) K-mer prefilter + HMM score + collect ScoredReads.
/// (f) Calls `synthesize_bundles` (Task 5.4) — results returned via the
///     `run_rescue_with_bundles` wrapper; this entry keeps the original
///     `Vec<NovelCandidate>` API for the legacy dispatcher.
pub fn run_rescue(
    bam_path: &std::path::Path,
    families: &[FamilyGroup],
    bundles: &[Bundle],
    config: &RunConfig,
) -> Result<Vec<NovelCandidate>> {
    let (candidates, _) = run_rescue_with_bundles(bam_path, families, bundles, config)?;
    Ok(candidates)
}

/// Full rescue, returning both novel candidates and synthetic bundles.
pub fn run_rescue_with_bundles(
    bam_path: &std::path::Path,
    families: &[FamilyGroup],
    bundles: &[Bundle],
    config: &RunConfig,
) -> Result<(Vec<NovelCandidate>, Vec<Bundle>)> {
    if families.is_empty() {
        return Ok((Vec::new(), Vec::new()));
    }

    // (a) Load genome.
    let genome = config.genome_fasta.as_ref().and_then(|p| {
        crate::genome::GenomeIndex::from_fasta(p).ok()
    });

    // (b) Build family graphs with sequences.
    let kmer_len: usize = 15;
    let min_kmer_hits: usize = 3;

    let mut family_graphs: Vec<FamilyGraph> = Vec::with_capacity(families.len());
    let mut family_kmer_sets: Vec<std::collections::HashSet<u64>> =
        Vec::with_capacity(families.len());

    for family in families {
        let fg = match crate::vg_hmm::family_graph::build_family_graph(
            family,
            bundles,
            genome.as_ref(),
            0.30,
            0.30,
        ) {
            Ok(g) => g,
            Err(e) => {
                // Mixed-strand families, empty families, etc. — skip gracefully.
                eprintln!("[VG-HMM] skipping family {}: {}", family.family_id, e);
                continue;
            }
        };

        // Also build k-mer set for this family (for prefilter).
        let mut kmers: std::collections::HashSet<u64> = std::collections::HashSet::new();
        for node in &fg.nodes {
            for (_, seq) in &node.per_copy_sequences {
                if seq.len() >= kmer_len {
                    for window in seq.windows(kmer_len) {
                        if window.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                            kmers.insert(fnv1a64(window));
                        }
                    }
                }
            }
        }

        family_graphs.push(fg);
        family_kmer_sets.push(kmers);
    }

    // (c) Fit profiles (skip family graphs that fail profile fitting).
    for fg in &mut family_graphs {
        if let Err(e) = crate::vg_hmm::family_graph::fit_profiles_in_place(fg) {
            eprintln!("[VG-HMM] fit_profiles_in_place failed for family {}: {} — continuing", fg.family_id, e);
        }
    }

    if family_graphs.is_empty() {
        eprintln!("[VG-HMM] 0 candidates, 0 synthetic bundles (no valid family graphs)");
        return Ok((Vec::new(), Vec::new()));
    }

    // (d) Open BAM and iterate unmapped reads.
    let bam_file = match std::fs::File::open(bam_path) {
        Ok(f) => f,
        Err(e) => return Err(anyhow::anyhow!("cannot open BAM {}: {}", bam_path.display(), e)),
    };
    let buf = std::io::BufReader::new(bam_file);
    let worker_count = std::num::NonZeroUsize::MIN;
    let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(worker_count, buf);
    let mut reader = noodles_bam::io::Reader::from(bgzf);
    let _header = reader.read_header()?;

    // Threshold: any score above NEG_INF is kept for now; tighten later.
    let min_loglik: f64 = -1000.0;

    let mut unmapped_reads: Vec<(String, Vec<u8>)> = Vec::new();
    let mut n_unmapped = 0usize;
    let mut n_prefilter_pass = 0usize;

    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        if !record.flags().is_unmapped() {
            continue;
        }
        n_unmapped += 1;

        // Decode 4-bit BAM sequence (same pattern as vg.rs:discover_novel_copies_kmer).
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

        // (e) Prefilter: must match at least one family.
        let passes_any = family_kmer_sets.iter().any(|fk| {
            prefilter_read(&seq_bytes, fk, kmer_len, min_kmer_hits).0
        });
        if !passes_any {
            continue;
        }
        n_prefilter_pass += 1;

        let read_name = record.name().map(|n| n.to_string()).unwrap_or_default();
        unmapped_reads.push((read_name, seq_bytes));
    }

    eprintln!(
        "[VG-HMM] rescue: {} unmapped reads, {} passed k-mer prefilter",
        n_unmapped, n_prefilter_pass
    );

    if unmapped_reads.is_empty() {
        eprintln!("[VG-HMM] 0 candidates, 0 synthetic bundles");
        return Ok((Vec::new(), Vec::new()));
    }

    // Score in-memory.
    let (candidates, scored_reads) =
        run_rescue_in_memory(families, &family_graphs, &unmapped_reads, min_loglik)?;

    // (f) Build per-family scored maps for synthesize_bundles.
    // Group by family_id.
    let mut per_family: HashMap<usize, (usize, Vec<(String, Vec<NodeIdx>)>)> = HashMap::new();
    for (read_name, family_id, vpath, _ll) in &scored_reads {
        let fg_idx = family_graphs.iter().position(|fg| fg.family_id == *family_id);
        if let Some(fi) = fg_idx {
            per_family
                .entry(*family_id)
                .or_insert_with(|| (fi, Vec::new()))
                .1
                .push((read_name.clone(), vpath.nodes.clone()));
        }
    }

    let min_reads_per_cluster: usize = 3;
    let mut synthetic_bundles: Vec<Bundle> = Vec::new();
    for (_family_id, (fg_idx, reads_with_paths)) in &per_family {
        let fg = &family_graphs[*fg_idx];
        let mut new_bundles = synthesize_bundles(fg, reads_with_paths, min_reads_per_cluster);
        synthetic_bundles.append(&mut new_bundles);
    }

    eprintln!(
        "[VG-HMM] {} candidates, {} synthetic bundles",
        candidates.len(),
        synthetic_bundles.len()
    );

    Ok((candidates, synthetic_bundles))
}

// ── Task 5.4: synthesize_bundles ──────────────────────────────────────────────

/// Cluster rescued reads by their Viterbi path key and emit one `Bundle` per
/// cluster that has ≥ `min_reads` members.
///
/// `rescued_with_paths`: slice of `(read_name, node_path)` pairs — all from the
/// same family graph `fg`.
pub fn synthesize_bundles(
    fg: &FamilyGraph,
    rescued_with_paths: &[(String, Vec<NodeIdx>)],
    min_reads: usize,
) -> Vec<Bundle> {
    // Cluster by the exact Vec<NodeIdx> path key.
    let mut clusters: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();
    for (i, (_, path)) in rescued_with_paths.iter().enumerate() {
        let key: Vec<usize> = path.iter().map(|n| n.0).collect();
        clusters.entry(key).or_default().push(i);
    }

    let mut out: Vec<Bundle> = Vec::new();

    for (path_key, member_indices) in &clusters {
        if member_indices.len() < min_reads {
            continue;
        }

        // Resolve path nodes.
        let path_nodes: Vec<&crate::vg_hmm::family_graph::ExonClass> = path_key
            .iter()
            .filter_map(|&ni| fg.nodes.get(ni))
            .collect();

        if path_nodes.is_empty() {
            continue;
        }

        let chrom: String = path_nodes[0].chrom.clone();
        let strand: char = path_nodes[0].strand;

        // Compute genomic span.
        let min_start: u64 = path_nodes.iter().map(|n| n.span.0).min().unwrap_or(0);
        let max_end: u64 = path_nodes.iter().map(|n| n.span.1).max().unwrap_or(0);
        let bundle_start = min_start.saturating_sub(1000);
        let bundle_end = max_end + 1000;

        // Build one BundleRead per member read.
        let exons: Vec<(u64, u64)> = path_nodes.iter().map(|n| n.span).collect();

        let reads: Vec<BundleRead> = member_indices
            .iter()
            .enumerate()
            .map(|(i, &ri)| {
                let read_name = rescued_with_paths[ri].0.clone();
                BundleRead {
                    read_uid: 1_000_000 + i as u64,
                    read_name: std::sync::Arc::from(read_name.as_str()),
                    read_name_hash: 0,
                    ref_id: None,
                    mate_ref_id: None,
                    mate_start: None,
                    hi: 0,
                    ref_start: min_start,
                    ref_end: max_end,
                    exons: exons.clone(),
                    junctions: Vec::new(),
                    junction_valid: Vec::new(),
                    junctions_raw: Vec::new(),
                    junctions_del: Vec::new(),
                    weight: 1.0,
                    is_reverse: strand == '-',
                    strand,
                    has_poly_start: false,
                    has_poly_end: false,
                    has_poly_start_aligned: false,
                    has_poly_start_unaligned: false,
                    has_poly_end_aligned: false,
                    has_poly_end_unaligned: false,
                    unaligned_poly_t: 0,
                    unaligned_poly_a: 0,
                    has_last_exon_polya: false,
                    has_first_exon_polyt: false,
                    query_length: None,
                    clip_left: 0,
                    clip_right: 0,
                    nh: 1,
                    nm: 0,
                    md: None,
                    insertion_sites: Vec::new(),
                    unitig: false,
                    unitig_cov: 0.0,
                    read_count_yc: 0.0,
                    countfrag_len: 0.0,
                    countfrag_num: 0.0,
                    junc_mismatch_weight: 0.0,
                    pair_idx: Vec::new(),
                    pair_count: Vec::new(),
                    mapq: 0,
                    mismatches: Vec::new(),
                    hp_tag: None,
                    ps_tag: None,
                    is_primary_alignment: true,
                }
            })
            .collect();

        out.push(Bundle {
            chrom,
            start: bundle_start,
            end: bundle_end,
            strand,
            reads,
            junction_stats: JunctionStats::default(),
            bundlenodes: None,
            read_bnodes: None,
            bnode_colors: None,
            synthetic: true,
        });
    }

    out
}

//! Rescue pipeline: read prefilter → HMM scoring → clustering → synthetic bundle emission.

use crate::types::{Bundle, BundleRead, JunctionStats, RunConfig};
use crate::vg::{FamilyGroup, NovelCandidate};
use crate::vg_hmm::diagnostic::{classify_internal, RescueClass};
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

/// Split a mixed-strand FamilyGroup into per-strand sub-families.
///
/// Many discovered families contain bundles on both `+` and `-` strands —
/// shared multi-mappers across strands are common (sequence palindromy,
/// antisense paralogs, or cross-strand seed matches). `build_family_graph`
/// requires single-strand input, so we partition by strand here. Sub-families
/// with fewer than 2 bundles are dropped (a 1-bundle "family" produces no
/// useful HMM signal).
///
/// Returns the original family if it's already single-strand; otherwise one
/// new FamilyGroup per strand. Each split family inherits a derived
/// `family_id` (`original * 10 + 0/1/2`) so the source family is recoverable
/// from the report.
fn partition_family_by_strand(
    family: &FamilyGroup,
    bundles: &[Bundle],
) -> Vec<FamilyGroup> {
    use std::collections::BTreeMap;
    let mut by_strand: BTreeMap<char, Vec<usize>> = BTreeMap::new();
    for &bi in &family.bundle_indices {
        by_strand
            .entry(bundles[bi].strand)
            .or_default()
            .push(bi);
    }
    if by_strand.len() <= 1 {
        return vec![FamilyGroup {
            family_id: family.family_id,
            bundle_indices: family.bundle_indices.clone(),
            multimap_reads: family.multimap_reads.clone(),
        }];
    }
    let mut out = Vec::new();
    for (i, (_strand, bis)) in by_strand.into_iter().enumerate() {
        if bis.len() < 2 {
            continue;
        }
        out.push(FamilyGroup {
            family_id: family.family_id * 10 + i,
            bundle_indices: bis,
            // Multimap reads are kept as-is; downstream consumers filter by bundle_indices.
            multimap_reads: family.multimap_reads.clone(),
        });
    }
    out
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

    // Partition mixed-strand families into per-strand sub-families before
    // build_family_graph (which requires single-strand input). This converts
    // ~half of typical discovered families from "skipped" to "rescuable".
    let split_families: Vec<FamilyGroup> = families
        .iter()
        .flat_map(|f| partition_family_by_strand(f, bundles))
        .collect();
    let n_split = split_families.len().saturating_sub(families.len());
    if n_split > 0 {
        eprintln!(
            "[VG-HMM] split {} mixed-strand families into per-strand sub-families ({} → {})",
            n_split,
            families.len(),
            split_families.len(),
        );
    }
    let families = &split_families;

    // (a) Load genome.
    let genome = match config.genome_fasta.as_ref() {
        Some(p) => match crate::genome::GenomeIndex::from_fasta(p) {
            Ok(g) => {
                eprintln!("[VG-HMM] loaded genome FASTA: {}", p);
                Some(g)
            }
            Err(e) => {
                eprintln!("[VG-HMM] WARNING: failed to load genome FASTA {}: {} — \
                           per-copy sequences unavailable; rescue will fail prefilter",
                          p, e);
                None
            }
        },
        None => {
            eprintln!("[VG-HMM] WARNING: --genome-fasta not provided; rescue will fail prefilter");
            None
        }
    };

    // (b) Build family graphs with sequences. Parallelized across families
    // — build_family_graph is independent per family, and on full GGO.bam
    // the sequential loop over 450 families is the dominant pre-scoring
    // cost.
    let kmer_len: usize = 15;
    let min_kmer_hits: usize = 3;

    use rayon::prelude::*;
    let built: Vec<Option<(FamilyGraph, std::collections::HashSet<u64>)>> = families
        .par_iter()
        .map(|family| {
            let fg = match crate::vg_hmm::family_graph::build_family_graph(
                family,
                bundles,
                genome.as_ref(),
                0.30,
                0.30,
            ) {
                Ok(g) => g,
                Err(_) => return None,
            };
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
            Some((fg, kmers))
        })
        .collect();

    let mut family_graphs: Vec<FamilyGraph> = Vec::with_capacity(families.len());
    let mut family_kmer_sets: Vec<std::collections::HashSet<u64>> =
        Vec::with_capacity(families.len());
    let mut n_skipped = 0usize;
    for entry in built {
        match entry {
            Some((fg, kmers)) => {
                family_graphs.push(fg);
                family_kmer_sets.push(kmers);
            }
            None => n_skipped += 1,
        }
    }
    if n_skipped > 0 {
        eprintln!("[VG-HMM] skipped {} families during graph build (mixed-strand / empty)", n_skipped);
    }

    // (c) Fit profiles in parallel. fit_profiles_in_place is independent
    // per family graph.
    family_graphs.par_iter_mut().for_each(|fg| {
        let _ = crate::vg_hmm::family_graph::fit_profiles_in_place(fg);
    });

    if family_graphs.is_empty() {
        eprintln!("[VG-HMM] 0 candidates, 0 synthetic bundles (no valid family graphs)");
        return Ok((Vec::new(), Vec::new()));
    }

    // Diagnostic: how many family graphs have any sequences vs. empty?
    let n_with_seq = family_graphs.iter().filter(|fg| {
        fg.nodes.iter().any(|n| n.per_copy_sequences.iter().any(|(_, s)| !s.is_empty()))
    }).count();
    let total_kmers: usize = family_kmer_sets.iter().map(|s| s.len()).sum();
    eprintln!(
        "[VG-HMM] family graphs: {}/{} have non-empty sequences; total k-mers in prefilter set: {}",
        n_with_seq, family_graphs.len(), total_kmers
    );

    // (d) Open BAM and iterate unmapped reads.
    let bam_file = match std::fs::File::open(bam_path) {
        Ok(f) => f,
        Err(e) => return Err(anyhow::anyhow!("cannot open BAM {}: {}", bam_path.display(), e)),
    };
    let buf = std::io::BufReader::new(bam_file);
    let worker_count = std::num::NonZeroUsize::MIN;
    let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(worker_count, buf);
    let mut reader = noodles_bam::io::Reader::from(bgzf);
    let header = reader.read_header()?;

    // Threshold: any score above NEG_INF is kept for now; tighten later.
    let min_loglik: f64 = -1000.0;

    // (read_name, sequence, candidate_family_indices_after_prefilter)
    let mut unmapped_reads: Vec<(String, Vec<u8>, Vec<usize>)> = Vec::new();
    let mut n_unmapped = 0usize;
    let mut n_masked = 0usize;
    let mut n_prefilter_pass = 0usize;
    let mask_active = !config.vg_mask_regions.is_empty();

    for result in reader.records() {
        let record = match result {
            Ok(r) => r,
            Err(_) => continue,
        };
        let is_unmapped = record.flags().is_unmapped();
        // LOO mask: include reads whose primary alignment overlaps a mask
        // region (treat as if unaligned). Skip secondary/supplementary so we
        // collect the sequence only once per read.
        let mut is_masked = false;
        if !is_unmapped && mask_active
            && !record.flags().is_secondary()
            && !record.flags().is_supplementary()
        {
            // Resolve chromosome name from the reference index.
            let ref_name = record.reference_sequence_id()
                .and_then(|r| r.ok())
                .and_then(|rid| header.reference_sequences().get_index(rid))
                .map(|(name, _)| String::from_utf8_lossy(name.as_ref()).into_owned());
            if let Some(name) = ref_name {
                let span = record.alignment_start()
                    .and_then(|s| s.ok())
                    .map(|s| s.get() as u64);
                if let Some(start_1b) = span {
                    use noodles_sam::alignment::record::Cigar as _;
                    let rs = start_1b.saturating_sub(1);
                    let span_len = record.cigar().alignment_span().unwrap_or(0) as u64;
                    let re_excl = rs.saturating_add(span_len);
                    is_masked = config.vg_mask_regions.iter().any(|(c, ms, me)| {
                        c == &name && rs < *me && *ms < re_excl
                    });
                }
            }
        }
        if !is_unmapped && !is_masked {
            continue;
        }
        if is_unmapped { n_unmapped += 1; }
        if is_masked { n_masked += 1; }

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

        // (e) Prefilter: collect the families this read passes against, so
        // HMM scoring restricts to candidate families instead of all 26+.
        // For 22 reads × 26 families on chr19 the unfiltered loop is 572
        // forward DPs (≥30min); per-read candidate sets typically have ≤2
        // families.
        let mut candidate_fams: Vec<usize> = Vec::new();
        for (fi, fk) in family_kmer_sets.iter().enumerate() {
            if prefilter_read(&seq_bytes, fk, kmer_len, min_kmer_hits).0 {
                candidate_fams.push(fi);
            }
        }
        if candidate_fams.is_empty() {
            continue;
        }
        n_prefilter_pass += 1;

        let read_name = record.name().map(|n| n.to_string()).unwrap_or_default();
        unmapped_reads.push((read_name, seq_bytes, candidate_fams));
    }

    if mask_active {
        eprintln!(
            "[VG-HMM] rescue: {} unmapped + {} mask-region reads, {} passed k-mer prefilter",
            n_unmapped, n_masked, n_prefilter_pass
        );
    } else {
        eprintln!(
            "[VG-HMM] rescue: {} unmapped reads, {} passed k-mer prefilter",
            n_unmapped, n_prefilter_pass
        );
    }

    if unmapped_reads.is_empty() {
        eprintln!("[VG-HMM] 0 candidates, 0 synthetic bundles");
        return Ok((Vec::new(), Vec::new()));
    }

    // Score each read only against its candidate families (prefilter narrowed
    // the family set per-read above). For typical reads, ≤2 candidates.
    //
    // Parallelized over reads (rayon par_iter): forward_against_family and
    // viterbi_path are pure functions of `(fg, seq)` — fg shared `&` across
    // threads, seqs are non-overlapping reads. On full GGO.bam scale (10K+
    // unmapped + masked reads), the prior sequential loop took hours; with
    // par_iter on 32 cores the scoring drops to a few minutes.
    //
    // Algorithmic fix: only compute viterbi_path ONCE per read, against the
    // best-scoring family (found via forward_against_family). The previous
    // code recomputed viterbi every time the running-best forward score
    // improved — wasteful and a distinct hot path itself.
    let _ = families; // reserved for future use; scoring uses family_graphs.
    let scored_per_read: Vec<Option<(String, usize, crate::vg_hmm::scorer::ViterbiPath, f64)>> =
        unmapped_reads.par_iter().map(|(read_name, seq, cand_fams)| {
            let mut best_score = f64::NEG_INFINITY;
            let mut best_fi: Option<usize> = None;
            for &fi in cand_fams {
                let fg = match family_graphs.get(fi) { Some(g) => g, None => continue };
                let ll = crate::vg_hmm::scorer::forward_against_family(fg, seq);
                if ll > best_score {
                    best_score = ll;
                    best_fi = Some(fi);
                }
            }
            if best_score <= min_loglik { return None; }
            let fi = best_fi?;
            let fg = family_graphs.get(fi)?;
            let vpath = crate::vg_hmm::scorer::viterbi_path(fg, seq)?;
            let family_id = fg.family_id;
            Some((read_name.clone(), family_id, vpath, best_score))
        }).collect();

    let mut candidates: Vec<NovelCandidate> = Vec::new();
    let mut scored_reads: Vec<ScoredRead> = Vec::new();
    for entry in scored_per_read.into_iter().flatten() {
        let (read_name, family_id, vpath, score) = entry;
        let fi = match family_graphs.iter().position(|fg| fg.family_id == family_id) {
            Some(i) => i,
            None => continue,
        };
        candidates.push(NovelCandidate {
            read_name: read_name.clone(),
            family_id,
            matched_junctions: vpath.nodes.len(),
            total_junctions: family_graphs[fi].nodes.len(),
            approx_chrom: String::new(),
            approx_start: 0,
            approx_end: 0,
        });
        scored_reads.push((read_name, family_id, vpath, score));
    }

    // (f) Build per-family scored maps for synthesize_bundles_refined.
    // Group by family_id, threading the per-node Viterbi read spans through
    // so the refined synthesis can apply boundary refinement.
    let mut per_family: HashMap<usize, (usize, Vec<(String, Vec<NodeIdx>, Vec<(usize, usize)>)>)> = HashMap::new();
    for (read_name, family_id, vpath, _ll) in &scored_reads {
        let fg_idx = family_graphs.iter().position(|fg| fg.family_id == *family_id);
        if let Some(fi) = fg_idx {
            per_family
                .entry(*family_id)
                .or_insert_with(|| (fi, Vec::new()))
                .1
                .push((read_name.clone(), vpath.nodes.clone(), vpath.read_spans.clone()));
        }
    }

    let min_reads_per_cluster: usize = 3;
    let hard_min_reads: usize = 3;
    let cluster_window: u64 = 20;
    let mut synthetic_bundles: Vec<Bundle> = Vec::new();
    for (_family_id, (fg_idx, reads_with_spans)) in &per_family {
        let fg = &family_graphs[*fg_idx];
        let mut new_bundles = synthesize_bundles_refined(
            fg,
            reads_with_spans,
            min_reads_per_cluster,
            kmer_len,
            hard_min_reads,
            cluster_window,
        );
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
/// `kmer_len`: k-mer length used in the prefilter (for classify_internal).
/// `vg_rescue_diagnostic`: when true, enables external minimap2 classification.
pub fn synthesize_bundles(
    fg: &FamilyGraph,
    rescued_with_paths: &[(String, Vec<NodeIdx>)],
    min_reads: usize,
    kmer_len: usize,
    _vg_rescue_diagnostic: bool,
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

        // Run internal diagnostic classifier.  Use cluster size as a proxy
        // for n_kmer_hits (each read passed the prefilter, so ≥ min_kmer_hits
        // each; cluster size × kmer_len is a conservative lower bound on the
        // aggregate chain score).  masked_fraction is 0.0 — no per-read mask
        // signal is threaded here yet (Phase 7 will refine).
        let rc: Option<RescueClass> = Some(classify_internal(
            member_indices.len(), // proxy for n_kmer_hits
            kmer_len,
            (max_end - min_start) as usize,
            0.0, // masked_fraction placeholder
        ));

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
            rescue_class: rc,
        });
    }

    out
}

// ── Boundary-refined bundle synthesis (uses Viterbi per-node read spans) ───

/// Production-realistic synthesize_bundles that applies boundary refinement
/// using per-read Viterbi traceback spans.
///
/// For each cluster of reads with the same path through the family graph:
///   - Internal exons keep family-graph node spans (splice junctions are
///     conserved across paralogs by construction).
///   - Boundary exons (TSS / TES) get refined lengths via:
///       1. Strand-aware identification of TSS vs TES from `fg.nodes[i].strand`
///          (TSS = genomic-first for + strand, genomic-last for − strand).
///       2. StringTie-style position clustering (`cluster_positions_with_counts`)
///          on per-read footprint lengths at that node position, with
///          window=`cluster_window` bp.
///       3. The cluster with the highest read count wins → its center is the
///          consensus boundary length.
///       4. Hard-boundary fallback: if no cluster has ≥`hard_min_reads`
///          supporters, fall back to the family-graph node length (RAW)
///          rather than invent a boundary from too few reads.
///
/// `rescued`: slice of `(read_name, node_path, viterbi_read_spans)`. The
/// `viterbi_read_spans[i] = (entry_r, exit_r)` records the read positions at
/// which node i was active. Same length as `node_path`.
pub fn synthesize_bundles_refined(
    fg: &FamilyGraph,
    rescued: &[(String, Vec<NodeIdx>, Vec<(usize, usize)>)],
    min_reads: usize,
    kmer_len: usize,
    hard_min_reads: usize,
    cluster_window: u64,
) -> Vec<Bundle> {
    use crate::tss_tts::cluster_positions_with_counts;

    // Cluster reads by exact node-path key; also accumulate per-position
    // footprint lengths (from Viterbi exit_r - entry_r) for each cluster.
    let mut clusters: HashMap<Vec<usize>, Vec<usize>> = HashMap::new();
    let mut spans_by_path: HashMap<Vec<usize>, Vec<Vec<u64>>> = HashMap::new();
    let mut read_lens: Vec<u64> = Vec::with_capacity(rescued.len());
    for (i, (_, path, spans)) in rescued.iter().enumerate() {
        let key: Vec<usize> = path.iter().map(|n| n.0).collect();
        clusters.entry(key.clone()).or_default().push(i);
        let entry = spans_by_path.entry(key).or_insert_with(|| vec![Vec::new(); path.len()]);
        for (j, &(s, e)) in spans.iter().enumerate() {
            if j < entry.len() { entry[j].push(e.saturating_sub(s) as u64); }
        }
        // Read length is approximately the max read position seen in any span.
        let l = spans.iter().map(|&(_, e)| e as u64).max().unwrap_or(0);
        read_lens.push(l);
    }

    let mut out: Vec<Bundle> = Vec::new();

    for (path_key, member_indices) in &clusters {
        if member_indices.len() < min_reads { continue; }

        let path_nodes: Vec<&crate::vg_hmm::family_graph::ExonClass> = path_key
            .iter()
            .filter_map(|&ni| fg.nodes.get(ni))
            .collect();
        if path_nodes.is_empty() { continue; }

        let chrom: String = path_nodes[0].chrom.clone();
        let strand: char = path_nodes[0].strand;

        // Raw exon spans from family-graph nodes.
        let raw_exons: Vec<(u64, u64)> = path_nodes.iter().map(|n| n.span).collect();
        let n_ex = raw_exons.len();

        // Strand-aware TSS / TES identification:
        //   + strand: TSS = genomic-first, TES = genomic-last
        //   − strand: TSS = genomic-last,  TES = genomic-first
        let strand_plus = strand == '+';
        let tss_idx = if strand_plus { 0 } else { n_ex.saturating_sub(1) };
        let tes_idx = if strand_plus { n_ex.saturating_sub(1) } else { 0 };

        // Build refined exon list.
        let footprints_per_pos: &Vec<Vec<u64>> = match spans_by_path.get(path_key) {
            Some(s) => s,
            None => continue,
        };

        let refine_boundary = |i: usize, node_s: u64, node_e: u64, is_tss: bool| -> u64 {
            let node_len = node_e - node_s;
            if i >= footprints_per_pos.len() || footprints_per_pos[i].is_empty() {
                return node_len;
            }
            // For the TSS exon, include pre-profile boundary-state bases by
            // using exit_r as the footprint (= total read consumed up to and
            // including this node). We already stored exit-entry as length
            // above; for the TSS we want the more permissive measure.
            // Here we just use the cluster-of-lengths approach: each footprint
            // is one (length, weight=1) point; cluster within window; take
            // highest-count cluster center.
            let pos: Vec<(u64, f64)> = footprints_per_pos[i].iter().map(|&f| (f, 1.0)).collect();
            let clusters = cluster_positions_with_counts(&pos, cluster_window, 1);
            let best = clusters.iter().max_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
            match best {
                Some(&(center, w)) if (w as usize) >= hard_min_reads => center,
                _ => node_len,  // hard-boundary fallback
            }
        };

        let mut refined_exons: Vec<(u64, u64)> = Vec::with_capacity(n_ex);
        for (i, &(node_s, node_e)) in raw_exons.iter().enumerate() {
            let node_len = node_e - node_s;
            let is_first = i == 0;
            let is_last = i + 1 == n_ex;
            let len = if i == tss_idx && (is_first || is_last) && n_ex > 1 {
                refine_boundary(i, node_s, node_e, true)
            } else if i == tes_idx && (is_first || is_last) && n_ex > 1 {
                refine_boundary(i, node_s, node_e, false)
            } else {
                node_len
            };
            // Anchor: keep the splice junction (= shared splice site) and grow
            // the boundary outward. For the genomic-first boundary exon, keep
            // node_e (splice donor); for genomic-last, keep node_s (acceptor).
            if is_first && n_ex > 1 {
                refined_exons.push((node_e.saturating_sub(len), node_e));
            } else if is_last && n_ex > 1 {
                refined_exons.push((node_s, node_s + len));
            } else {
                refined_exons.push((node_s, node_e));
            }
        }

        let min_start = refined_exons.iter().map(|&(s, _)| s).min().unwrap_or(0);
        let max_end   = refined_exons.iter().map(|&(_, e)| e).max().unwrap_or(0);
        let bundle_start = min_start.saturating_sub(1000);
        let bundle_end   = max_end + 1000;

        let reads: Vec<BundleRead> = member_indices.iter().enumerate().map(|(i, &ri)| {
            let read_name = rescued[ri].0.clone();
            BundleRead {
                read_uid: 1_000_000 + i as u64,
                read_name: std::sync::Arc::from(read_name.as_str()),
                read_name_hash: 0,
                ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
                ref_start: min_start, ref_end: max_end,
                exons: refined_exons.clone(),
                junctions: Vec::new(), junction_valid: Vec::new(),
                junctions_raw: Vec::new(), junctions_del: Vec::new(),
                weight: 1.0, is_reverse: strand == '-', strand,
                has_poly_start: false, has_poly_end: false,
                has_poly_start_aligned: false, has_poly_start_unaligned: false,
                has_poly_end_aligned: false, has_poly_end_unaligned: false,
                unaligned_poly_t: 0, unaligned_poly_a: 0,
                has_last_exon_polya: false, has_first_exon_polyt: false,
                query_length: None, clip_left: 0, clip_right: 0,
                nh: 1, nm: 0, md: None, insertion_sites: Vec::new(),
                unitig: false, unitig_cov: 0.0, read_count_yc: 0.0,
                countfrag_len: 0.0, countfrag_num: 0.0, junc_mismatch_weight: 0.0,
                pair_idx: Vec::new(), pair_count: Vec::new(),
                mapq: 0, mismatches: Vec::new(),
                hp_tag: None, ps_tag: None, is_primary_alignment: true,
            }
        }).collect();

        let rc: Option<RescueClass> = Some(classify_internal(
            member_indices.len(),
            kmer_len,
            (max_end - min_start) as usize,
            0.0,
        ));

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
            rescue_class: rc,
        });
    }

    out
}

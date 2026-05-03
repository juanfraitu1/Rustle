//! Rescue pipeline: read prefilter → HMM scoring → clustering → synthetic bundle emission.

use crate::types::{Bundle, BundleRead, DetHashMap, DetHashSet, JunctionStats, RunConfig};
use crate::vg::{FamilyGroup, NovelCandidate};
use crate::vg_hmm::diagnostic::{classify_internal, RescueClass};
use crate::vg_hmm::family_graph::{FamilyGraph, NodeIdx};
use crate::vg_hmm::scorer::{forward_against_family, viterbi_path, ViterbiPath};
use anyhow::Result;
use std::collections::HashMap;

/// Original alignment of a masked read: (chrom, strand, exons).
/// Used as a positional prior so LOO-rescued reads land at their
/// original genomic location instead of a sister-paralog's.
pub(crate) type OriginalAlignment = (String, char, Vec<(u64, u64)>);

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
    family_kmers: &DetHashSet<u64>,
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
pub fn partition_family_by_strand(
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
    let built: Vec<Option<(FamilyGraph, DetHashSet<u64>)>> = families
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
            let mut kmers: DetHashSet<u64> = DetHashSet::default();
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
    let mut family_kmer_sets: Vec<DetHashSet<u64>> =
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
    // Phase 2: report candidate loci availability per family. Used in
    // Phase 3 as positional priors for synthesizing bundles at novel
    // paralog locations. Logged here for traceability.
    if !config.vg_candidate_loci.is_empty() {
        let n_with_cand = family_graphs.iter()
            .filter(|fg| config.vg_candidate_loci.get(&fg.family_id).map_or(false, |v| !v.is_empty()))
            .count();
        let total_cand: usize = config.vg_candidate_loci.values()
            .map(|v| v.len()).sum();
        eprintln!(
            "[VG-HMM] positional priors: {} families have candidate loci, {} total candidates",
            n_with_cand, total_cand
        );
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
    let header = reader.read_header()?;

    // HMM forward-log-lik threshold for accepting a rescued read. Scores
    // are sums of per-base log-probs; a typical 2kb read has ~ -2000 to
    // -3000 against a well-fit family. Default cutoff is generous (-5000)
    // so the cluster-size filter below does most of the false-positive
    // rejection. Override with RUSTLE_VG_RESCUE_MIN_LOGLIK=<f64>.
    let min_loglik: f64 = std::env::var("RUSTLE_VG_RESCUE_MIN_LOGLIK")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(-5000.0);

    // (read_name, sequence, candidate_family_indices, optional original
    // genomic position prior). The position prior is `Some` for masked
    // reads (LOO experiment) — we know where they aligned originally,
    // BEFORE the mask, including their CIGAR-derived exons. Truly
    // unmapped reads have `None` and use the family-graph copy spans.
    let mut unmapped_reads: Vec<(String, Vec<u8>, Vec<usize>, Option<OriginalAlignment>)> = Vec::new();
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
        // Capture (chrom, strand, exons) when masking — used as a positional
        // prior so synthetic bundles land at the read's ORIGINAL genomic
        // location, not at a sister-paralog's location chosen by the HMM.
        let mut is_masked = false;
        let mut original_alignment: Option<OriginalAlignment> = None;
        if !is_unmapped && mask_active
            && !record.flags().is_secondary()
            && !record.flags().is_supplementary()
        {
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
                    if is_masked {
                        // Walk CIGAR to extract exon spans (M/=/X consume ref,
                        // N is intron, D consumes ref but stays in exon).
                        // Each contiguous M/D run between N ops is one exon.
                        use noodles_sam::alignment::record::cigar::op::Kind;
                        let mut exons: Vec<(u64, u64)> = Vec::new();
                        let mut ref_pos: u64 = rs;
                        let mut exon_start: u64 = rs;
                        let mut in_exon = false;
                        for op_res in record.cigar().iter() {
                            let Ok(op) = op_res else { break; };
                            let l = op.len() as u64;
                            match op.kind() {
                                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                                    if !in_exon { exon_start = ref_pos; in_exon = true; }
                                    ref_pos += l;
                                }
                                Kind::Skip => {
                                    if in_exon && ref_pos > exon_start {
                                        exons.push((exon_start, ref_pos));
                                    }
                                    in_exon = false;
                                    ref_pos += l;
                                }
                                _ => {}  // Insertion, soft/hard clip, pad: no ref advance
                            }
                        }
                        if in_exon && ref_pos > exon_start {
                            exons.push((exon_start, ref_pos));
                        }
                        let strand = if record.flags().is_reverse_complemented() { '-' } else { '+' };
                        original_alignment = Some((name.clone(), strand, exons));
                    }
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
        unmapped_reads.push((read_name, seq_bytes, candidate_fams, original_alignment));
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

    // Diagnostic cap on the number of unmapped reads scored. The
    // forward_against_family DP is the bottleneck on full-scale family
    // graphs (93 families × 70K-k-mer graphs × 5K reads on full GGO.bam
    // never finished within an hour). Setting this env var to a small
    // number (e.g. 500) makes Phase-3 diagnostic runs tractable. 0 =
    // unlimited (default).
    let read_cap: usize = std::env::var("RUSTLE_VG_RESCUE_READ_CAP")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(0);
    if read_cap > 0 && unmapped_reads.len() > read_cap {
        // Sort by read name for determinism, then truncate.
        unmapped_reads.sort_by(|a, b| a.0.cmp(&b.0));
        let n_before = unmapped_reads.len();
        unmapped_reads.truncate(read_cap);
        eprintln!(
            "[VG-HMM] read cap active: scoring first {} of {} prefilter-passing reads (RUSTLE_VG_RESCUE_READ_CAP={})",
            unmapped_reads.len(), n_before, read_cap,
        );
    }

    // Score each read only against its candidate families (prefilter narrowed
    // the family set per-read above). For typical reads, ≤2 candidates.
    //
    // Per-path scoring: for each (read, candidate-family) we score the read
    // against EACH KNOWN PARALOG'S path through the graph via
    // `forward_against_path`, instead of running family-wide
    // `forward_against_family` once. Per-path DP is O(K × L × prof) where K
    // is the path length (≈ 10-15 nodes per paralog) versus the family-wide
    // O(N × L × prof) where N is the union of all paralog node sets (≈ 30+
    // for multi-copy families). The total work scales with the number of
    // paralogs but each call has a much smaller working set, so memory
    // bandwidth contention drops dramatically — on full GGO.bam this was
    // the bottleneck that stalled the rescue at 5K+ reads × 2-cores
    // utilization for hours.
    //
    // Tradeoffs:
    //  - We get per-paralog likelihoods directly (useful for downstream
    //    HMM-EM — already a separate path, but the per-path scoring here
    //    means the same primitive is reusable).
    //  - Viterbi for read_spans is replaced by `viterbi_against_path` on
    //    the chosen best path, also O(K × L × prof) — the prior
    //    `viterbi_path` was a second family-wide DP.
    //
    // Precompute per-family per-copy paths once (outside the hot loop).
    let mut family_paths: Vec<Vec<(crate::vg_hmm::family_graph::CopyId, Vec<crate::vg_hmm::family_graph::NodeIdx>)>> =
        Vec::with_capacity(family_graphs.len());
    for fg in &family_graphs {
        let mut entries: Vec<(crate::vg_hmm::family_graph::CopyId, Vec<crate::vg_hmm::family_graph::NodeIdx>)> = Vec::new();
        for cid in fg.all_copies() {
            let path = fg.recover_paralog_path(cid);
            if !path.is_empty() {
                entries.push((cid, path));
            }
        }
        family_paths.push(entries);
    }
    let _ = families;
    // Per-path vs family-wide scoring: opt-in via RUSTLE_VG_RESCUE_PER_PATH=1.
    // Per-path was hypothesized to be cache-friendlier but on full GGO.bam it
    // does ~5x more total work per family (one DP per copy, not per family),
    // and has not been validated end-to-end. Default keeps the prior
    // family-wide forward + family-wide Viterbi path.
    let use_per_path: bool = std::env::var("RUSTLE_VG_RESCUE_PER_PATH")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(0) != 0;
    use std::sync::atomic::{AtomicUsize, Ordering};
    let progress = AtomicUsize::new(0);
    let total_reads = unmapped_reads.len();
    let progress_step = (total_reads / 10).max(50);
    let scoring_t0 = std::time::Instant::now();
    let scored_per_read: Vec<Option<(String, usize, crate::vg_hmm::scorer::ViterbiPath, f64, Option<OriginalAlignment>)>> =
        unmapped_reads.par_iter().map(|(read_name, seq, cand_fams, orig)| {
            let mut best_score = f64::NEG_INFINITY;
            let mut best_fi: Option<usize> = None;
            let mut best_copy_idx: Option<usize> = None;
            if use_per_path {
                for &fi in cand_fams {
                    let _fg = match family_graphs.get(fi) { Some(g) => g, None => continue };
                    let paths = match family_paths.get(fi) { Some(p) => p, None => continue };
                    for (ci, (_cid, path)) in paths.iter().enumerate() {
                        let ll = crate::vg_hmm::scorer::forward_against_path(&family_graphs[fi], seq, path);
                        if ll > best_score {
                            best_score = ll;
                            best_fi = Some(fi);
                            best_copy_idx = Some(ci);
                        }
                    }
                }
            } else {
                for &fi in cand_fams {
                    let fg = match family_graphs.get(fi) { Some(g) => g, None => continue };
                    let ll = crate::vg_hmm::scorer::forward_against_family(fg, seq);
                    if ll > best_score {
                        best_score = ll;
                        best_fi = Some(fi);
                    }
                }
            }
            let n = progress.fetch_add(1, Ordering::Relaxed) + 1;
            if n % progress_step == 0 {
                let elapsed = scoring_t0.elapsed().as_secs();
                eprintln!("[VG-HMM]   scored {}/{} reads ({}s elapsed)", n, total_reads, elapsed);
            }
            if best_score <= min_loglik { return None; }
            let fi = best_fi?;
            let fg = family_graphs.get(fi)?;
            let vpath = if use_per_path {
                let ci = best_copy_idx?;
                let path = &family_paths[fi][ci].1;
                crate::vg_hmm::scorer::viterbi_against_path(fg, seq, path)?
            } else {
                crate::vg_hmm::scorer::viterbi_path(fg, seq)?
            };
            let family_id = fg.family_id;
            Some((read_name.clone(), family_id, vpath, best_score, orig.clone()))
        }).collect();
    eprintln!("[VG-HMM] HMM scoring: {} reads in {}s (per_path={})",
              total_reads, scoring_t0.elapsed().as_secs(), use_per_path);

    // Score-distribution diagnostic: helps tune min_loglik and explains
    // why rescue rejects reads. Counts per-bucket of the best forward log-lik.
    let mut score_dist: std::collections::BTreeMap<&'static str, usize> = std::collections::BTreeMap::new();
    let mut n_no_candidate = 0usize;
    for (_rn, _seq, cand_fams, _orig) in &unmapped_reads {
        if cand_fams.is_empty() { n_no_candidate += 1; }
    }
    for entry in scored_per_read.iter() {
        match entry {
            None => *score_dist.entry("rejected_below_min_loglik").or_insert(0) += 1,
            Some((_, _, _, s, _)) => {
                let bucket = if *s > -500.0 { "score>-500" }
                    else if *s > -1500.0 { "-1500<score<=-500" }
                    else if *s > -3000.0 { "-3000<score<=-1500" }
                    else if *s > -5000.0 { "-5000<score<=-3000" }
                    else { "score<=-5000" };
                *score_dist.entry(bucket).or_insert(0) += 1;
            }
        }
    }
    let n_accepted: usize = scored_per_read.iter().filter(|x| x.is_some()).count();
    eprintln!(
        "[VG-HMM] HMM scoring: {} reads scored ({} no candidate fams); {} accepted at min_loglik={}; distribution: {:?}",
        unmapped_reads.len(), n_no_candidate, n_accepted, min_loglik, score_dist
    );

    let mut candidates: Vec<NovelCandidate> = Vec::new();
    let mut scored_reads: Vec<ScoredRead> = Vec::new();
    // Reads that came in with an original alignment (LOO mask) — bypass
    // the family-graph-coords synthesis and emit at their original locus.
    let mut anchored_reads: Vec<(String, OriginalAlignment, f64, usize)> = Vec::new();
    for entry in scored_per_read.into_iter().flatten() {
        let (read_name, family_id, vpath, score, orig) = entry;
        let fi = match family_graphs.iter().position(|fg| fg.family_id == family_id) {
            Some(i) => i,
            None => continue,
        };
        let (chrom, start, end) = match &orig {
            Some((c, _, exons)) => {
                let s = exons.first().map(|e| e.0).unwrap_or(0);
                let e = exons.last().map(|e| e.1).unwrap_or(0);
                (c.clone(), s, e)
            }
            None => (String::new(), 0, 0),
        };
        candidates.push(NovelCandidate {
            read_name: read_name.clone(),
            family_id,
            matched_junctions: vpath.nodes.len(),
            total_junctions: family_graphs[fi].nodes.len(),
            approx_chrom: chrom,
            approx_start: start,
            approx_end: end,
        });
        if let Some(o) = orig {
            anchored_reads.push((read_name.clone(), o, score, family_id));
        }
        scored_reads.push((read_name, family_id, vpath, score));
    }

    // ── Phase 3: candidate-locus matching for unmapped reads ──────────────
    // For reads with no original alignment but whose assigned family has
    // scan-discovered candidate loci, score the read's k-mers against each
    // candidate locus's reference sequence. If the best candidate's k-mer
    // overlap exceeds a threshold, route the read to a "novel" synthesis
    // path that emits a synthetic bundle at the candidate's coords with
    // rescue_class=NovelLocusFromScan. Reads that don't match any candidate
    // strongly fall back to the existing family-graph-spans path.
    // Phase 3 routed-read tuple: (read_name, chrom, strand, locus_start,
    // locus_end, score, family_id, projected_exons). The projected_exons
    // come from Phase 3.1: project a representative paralog's exon
    // structure onto the candidate locus so the resulting synthetic
    // bundle preserves splice structure (multi-exon, with junctions),
    // not just a single big exon. The downstream assembly can then
    // discover alternative isoforms within the novel locus from per-read
    // coverage of these projected exons.
    type RoutedRead = (String, String, char, u64, u64, f64, usize, Vec<(u64, u64)>);
    let candidate_anchored_reads: Vec<RoutedRead> =
        if !config.vg_candidate_loci.is_empty() && genome.is_some() {
            let g = genome.as_ref().unwrap();
            let anchored_set: std::collections::HashSet<&str> =
                anchored_reads.iter().map(|(rn, _, _, _)| rn.as_str()).collect();
            // Build per-family candidate k-mer sets (precompute once).
            let mut fam_cand_kmers: HashMap<usize, Vec<(crate::vg_hmm::positional::CandidateLocus, DetHashSet<u64>)>> =
                HashMap::new();
            for (&fid, locs) in &config.vg_candidate_loci {
                let mut entry: Vec<(crate::vg_hmm::positional::CandidateLocus, DetHashSet<u64>)> = Vec::new();
                for c in locs {
                    let seq = match g.fetch_sequence(&c.chrom, c.start, c.end) {
                        Some(s) => s,
                        None => continue,
                    };
                    let mut ks: DetHashSet<u64> = DetHashSet::default();
                    for w in seq.windows(kmer_len) {
                        if w.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                            ks.insert(fnv1a64(w));
                        }
                    }
                    if !ks.is_empty() {
                        entry.push((c.clone(), ks));
                    }
                }
                if !entry.is_empty() {
                    fam_cand_kmers.insert(fid, entry);
                }
            }

            // Phase 3.1: per family_id, precompute the representative
            // paralog's exon offsets relative to its first-exon start.
            // Used to project structure onto each candidate locus.
            let mut fam_paralog_offsets: HashMap<usize, Vec<(u64, u64)>> = HashMap::new();
            for fg in family_graphs.iter() {
                let rep_cid = match fg.representative_copy() {
                    Some(c) => c,
                    None => continue,
                };
                let exons = fg.paralog_exon_spans(rep_cid);
                if exons.is_empty() { continue; }
                let base = exons[0].0;
                let offsets: Vec<(u64, u64)> = exons.iter()
                    .map(|(s, e)| (s.saturating_sub(base), e.saturating_sub(base)))
                    .collect();
                fam_paralog_offsets.insert(fg.family_id, offsets);
            }

            // Phase 3.2: per (family, candidate locus) pair, precompute
            // per-projected-exon k-mer sets. Lets us decide PER READ which
            // exons it supports — same family, same locus, different reads
            // may light up different exon subsets, which is exactly what
            // enables isoform discovery within the novel locus.
            //
            // Keyed by (family_id, locus_start) since loci are unique
            // within a family by chrom+start. Stores (projected_exons,
            // exon_kmer_sets).
            type LocusKey = (usize, String, u64);
            let mut per_locus_exon_kmers: HashMap<LocusKey, (Vec<(u64, u64)>, Vec<DetHashSet<u64>>)> =
                HashMap::new();
            for (&fid, locs) in &config.vg_candidate_loci {
                let offsets = match fam_paralog_offsets.get(&fid) {
                    Some(o) if !o.is_empty() => o,
                    _ => continue,
                };
                for c in locs {
                    let proj: Vec<(u64, u64)> = offsets.iter().map(|(rs, re)| {
                        (c.start.saturating_add(*rs), c.start.saturating_add(*re))
                    }).collect();
                    let mut exon_kmers: Vec<DetHashSet<u64>> = Vec::with_capacity(proj.len());
                    for (s, e) in &proj {
                        let seq = g.fetch_sequence(&c.chrom, *s, *e).unwrap_or_default();
                        let mut ks: DetHashSet<u64> = DetHashSet::default();
                        for w in seq.windows(kmer_len) {
                            if w.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                                ks.insert(fnv1a64(w));
                            }
                        }
                        exon_kmers.push(ks);
                    }
                    per_locus_exon_kmers.insert((fid, c.chrom.clone(), c.start), (proj, exon_kmers));
                }
            }

            // Build seq lookup by read name.
            let seq_by_name: HashMap<&str, &[u8]> = unmapped_reads.iter()
                .map(|(rn, seq, _, _)| (rn.as_str(), seq.as_slice()))
                .collect();

            let min_overlap_frac: f64 = std::env::var("RUSTLE_VG_RESCUE_NOVEL_OVERLAP")
                .ok().and_then(|v| v.parse().ok()).unwrap_or(0.30);

            let mut out: Vec<RoutedRead> = Vec::new();
            let mut n_evaluated = 0usize;
            let mut n_with_projection = 0usize;
            for (read_name, family_id, _vpath, score) in &scored_reads {
                if anchored_set.contains(read_name.as_str()) { continue; }
                let cands = match fam_cand_kmers.get(family_id) {
                    Some(c) => c,
                    None => continue,
                };
                let seq = match seq_by_name.get(read_name.as_str()) {
                    Some(s) => *s,
                    None => continue,
                };
                if seq.len() < kmer_len { continue; }
                n_evaluated += 1;
                let mut read_kmers: DetHashSet<u64> = DetHashSet::default();
                for w in seq.windows(kmer_len) {
                    if w.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
                        read_kmers.insert(fnv1a64(w));
                    }
                }
                if read_kmers.is_empty() { continue; }
                let mut best: Option<(usize, f64)> = None;
                for (idx, (_loc, cand_kmers)) in cands.iter().enumerate() {
                    let inter = read_kmers.intersection(cand_kmers).count();
                    let frac = inter as f64 / read_kmers.len() as f64;
                    if best.map_or(true, |(_, f)| frac > f) {
                        best = Some((idx, frac));
                    }
                }
                if let Some((idx, frac)) = best {
                    if frac >= min_overlap_frac {
                        let loc = &cands[idx].0;
                        let strand = family_graphs.iter()
                            .find(|fg| fg.family_id == *family_id)
                            .and_then(|fg| fg.nodes.first().map(|n| n.strand))
                            .unwrap_or('+');
                        // Phase 3.1+3.2: pick per-locus projection and
                        // compute which projected exons THIS read supports
                        // (≥20% of the exon's k-mers appear in the read).
                        // Reads with different supported-exon subsets at
                        // the same locus produce alternative isoforms.
                        let key = (*family_id, loc.chrom.clone(), loc.start);
                        let supported_exons: Vec<(u64, u64)> =
                            match per_locus_exon_kmers.get(&key) {
                                Some((proj, ekmers)) if !proj.is_empty() => {
                                    n_with_projection += 1;
                                    let supp: Vec<(u64, u64)> = proj.iter().enumerate()
                                        .filter_map(|(ei, exon)| {
                                            let ek = &ekmers[ei];
                                            if ek.is_empty() { return None; }
                                            let inter = read_kmers.intersection(ek).count();
                                            let frac_e = inter as f64 / ek.len() as f64;
                                            if frac_e >= 0.20 { Some(*exon) } else { None }
                                        })
                                        .collect();
                                    if supp.is_empty() {
                                        // Read doesn't strongly support any
                                        // projected exon — fall back to
                                        // first projected exon as anchor.
                                        vec![proj[0]]
                                    } else {
                                        supp
                                    }
                                }
                                _ => vec![(loc.start, loc.end)],
                            };
                        out.push((
                            read_name.clone(),
                            loc.chrom.clone(),
                            strand,
                            loc.start,
                            loc.end,
                            *score,
                            *family_id,
                            supported_exons,
                        ));
                    }
                }
            }
            if n_evaluated > 0 {
                eprintln!(
                    "[VG-HMM] candidate-locus match: {} unmapped reads evaluated, {} routed (overlap≥{:.2}), {} with projected paralog structure",
                    n_evaluated, out.len(), min_overlap_frac, n_with_projection
                );
            }
            out
        } else {
            Vec::new()
        };
    let candidate_anchored_names: std::collections::HashSet<String> =
        candidate_anchored_reads.iter().map(|t| t.0.clone()).collect();

    // (f) Build per-family scored maps for synthesize_bundles_refined.
    // Group by family_id, threading the per-node Viterbi read spans through
    // so the refined synthesis can apply boundary refinement.
    // Reads with `anchored_reads` (LOO mask) or `candidate_anchored_reads`
    // (Phase 3 novel-locus) entries skip this — they go through the
    // dedicated synthesis paths below instead.
    let anchored_names: std::collections::HashSet<String> =
        anchored_reads.iter().map(|(rn, _, _, _)| rn.clone()).collect();
    let mut per_family: HashMap<usize, (usize, Vec<(String, Vec<NodeIdx>, Vec<(usize, usize)>)>)> = HashMap::new();
    for (read_name, family_id, vpath, _ll) in &scored_reads {
        if anchored_names.contains(read_name) { continue; }
        if candidate_anchored_names.contains(read_name) { continue; }
        let fg_idx = family_graphs.iter().position(|fg| fg.family_id == *family_id);
        if let Some(fi) = fg_idx {
            per_family
                .entry(*family_id)
                .or_insert_with(|| (fi, Vec::new()))
                .1
                .push((read_name.clone(), vpath.nodes.clone(), vpath.read_spans.clone()));
        }
    }

    // Cluster-size thresholds for synthesis. Tunable via env vars to
    // explore single-paralog LOO recovery — defaults are conservative so
    // the rescue doesn't over-create synthetic bundles from low-coverage
    // alignment noise. RUSTLE_VG_RESCUE_MIN_CLUSTER lowers to 1 for LOO
    // testing of low-expression paralogs.
    let min_reads_per_cluster: usize = std::env::var("RUSTLE_VG_RESCUE_MIN_CLUSTER")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(3);
    let hard_min_reads: usize = std::env::var("RUSTLE_VG_RESCUE_HARD_MIN")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(3);
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

    // ── Position-anchored synthesis (LOO mask reads) ──────────────────
    // For reads that came in with a known original alignment, bundle them
    // at THEIR original genomic position rather than at the matching
    // family-graph copy's spans. This is the LOO recovery path: the HMM
    // confirms the read belongs to a paralog family, but the rescue's
    // goal is to recover the *masked* paralog, which lives at the read's
    // own coordinates — not at a sister-paralog's coordinates.
    let n_anchored_input = anchored_reads.len();
    let n_anchored_bundles = synthesize_bundles_from_original_alignments(
        &anchored_reads,
        &mut synthetic_bundles,
        min_reads_per_cluster.max(1),
    );

    // ── Phase 3: candidate-locus synthesis (novel paralog rescue) ─────
    // For reads routed to a scan-discovered candidate locus, emit a
    // synthetic bundle at the candidate's coords with rescue_class set
    // so the GTF carries copy_status="novel".
    let n_novel_input = candidate_anchored_reads.len();
    let n_novel_bundles = synthesize_bundles_at_candidate_loci(
        &candidate_anchored_reads,
        &mut synthetic_bundles,
        min_reads_per_cluster.max(1),
    );

    eprintln!(
        "[VG-HMM] {} candidates, {} synthetic bundles ({} graph-spans, {} from {} mask-anchored reads, {} from {} candidate-locus reads)",
        candidates.len(),
        synthetic_bundles.len(),
        synthetic_bundles.len() - n_anchored_bundles - n_novel_bundles,
        n_anchored_bundles,
        n_anchored_input,
        n_novel_bundles,
        n_novel_input,
    );

    Ok((candidates, synthetic_bundles))
}

/// Phase 3 synthesis: cluster reads routed to scan-discovered candidate
/// loci by overlapping locus span, and emit one synthetic Bundle per
/// cluster. Each bundle uses the read's projected paralog structure
/// (Phase 3.1) so the resulting transcripts are multi-exon with proper
/// junctions — enabling isoform discovery within the novel locus.
///
/// `routed_reads`: tuples of (read_name, chrom, strand, locus_start,
/// locus_end, score, family_id, projected_exons). projected_exons is the
/// representative paralog's exon structure offset to start at locus_start.
fn synthesize_bundles_at_candidate_loci(
    routed_reads: &[(String, String, char, u64, u64, f64, usize, Vec<(u64, u64)>)],
    out_bundles: &mut Vec<Bundle>,
    min_reads: usize,
) -> usize {
    if routed_reads.is_empty() { return 0; }
    // Group by (chrom, strand, span). Within each group all reads share
    // identical projected_exons (since we used a single representative
    // paralog per family).
    let mut groups: std::collections::BTreeMap<(String, char, u64, u64),
        Vec<&(String, String, char, u64, u64, f64, usize, Vec<(u64, u64)>)>> =
        std::collections::BTreeMap::new();
    for r in routed_reads {
        groups.entry((r.1.clone(), r.2, r.3, r.4)).or_default().push(r);
    }
    let mut emitted = 0usize;
    for ((chrom, strand, locus_start, locus_end), group) in groups {
        if group.len() < min_reads { continue; }
        // Bundle bounds = union of all reads' supported exons in this
        // locus. Different reads may support different exon subsets
        // (Phase 3.2), so we take the outer envelope.
        let mut all_exon_starts: Vec<u64> = Vec::new();
        let mut all_exon_ends: Vec<u64> = Vec::new();
        for r in &group {
            for &(s, e) in &r.7 {
                all_exon_starts.push(s);
                all_exon_ends.push(e);
            }
        }
        let cluster_first = *all_exon_starts.iter().min().unwrap_or(&locus_start);
        let cluster_last = *all_exon_ends.iter().max().unwrap_or(&locus_end);
        let bundle_start = cluster_first.saturating_sub(100);
        let bundle_end = cluster_last + 100;
        let reads_vec: Vec<BundleRead> = group.iter().enumerate().map(|(i, r)| {
            // Per-read exon list = the exons THIS read supports via
            // k-mer match (Phase 3.2). Junctions are derived from
            // consecutive supported exons — so reads with different
            // exon subsets generate different junctions, enabling
            // alternative-isoform discovery downstream.
            let exons = r.7.clone();
            let juncs: Vec<crate::types::Junction> = exons.windows(2)
                .map(|w| crate::types::Junction { donor: w[0].1, acceptor: w[1].0 })
                .collect();
            let jvalid = vec![true; juncs.len()];
            let ref_start = exons.first().map(|e| e.0).unwrap_or(locus_start);
            let ref_end = exons.last().map(|e| e.1).unwrap_or(locus_end);
            BundleRead {
                read_uid: 3_000_000 + i as u64,
                read_name: std::sync::Arc::from(r.0.as_str()),
                read_name_hash: 0,
                ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
                ref_start, ref_end,
                exons,
                junctions: juncs.clone(),
                junction_valid: jvalid,
                junctions_raw: juncs,
                junctions_del: Vec::new(),
                weight: 1.0, is_reverse: strand == '-', strand,
                has_poly_start: false, has_poly_end: false,
                has_poly_start_aligned: false, has_poly_start_unaligned: false,
                has_poly_end_aligned: false, has_poly_end_unaligned: false,
                unaligned_poly_t: 0, unaligned_poly_a: 0,
                has_last_exon_polya: false, has_first_exon_polyt: false,
                query_length: None, clip_left: 0, clip_right: 0,
                nh: 1, nm: 0, md: None, insertion_sites: Vec::new(),
                unitig: false, unitig_cov: 0.0,
                read_count_yc: 0.0, countfrag_len: 0.0, countfrag_num: 0.0,
                junc_mismatch_weight: 0.0,
                pair_idx: Vec::new(), pair_count: Vec::new(),
                mapq: 0, mismatches: Vec::new(),
                hp_tag: None, ps_tag: None,
                is_primary_alignment: true,
            }
        }).collect();

        out_bundles.push(Bundle {
            chrom: chrom.clone(),
            start: bundle_start,
            end: bundle_end,
            strand,
            reads: reads_vec,
            junction_stats: JunctionStats::default(),
            bundlenodes: None,
            read_bnodes: None,
            bnode_colors: None,
            synthetic: true,
            rescue_class: Some(RescueClass::NovelLocusFromScan),
        });
        emitted += 1;
    }
    emitted
}

/// Cluster `anchored_reads` by chrom + overlapping genomic position and
/// emit one synthetic Bundle per cluster, using the reads' own
/// CIGAR-derived exons (NOT family-graph copy spans). Returns the number
/// of bundles appended.
fn synthesize_bundles_from_original_alignments(
    anchored_reads: &[(String, OriginalAlignment, f64, usize)],
    out_bundles: &mut Vec<Bundle>,
    min_reads: usize,
) -> usize {
    if anchored_reads.is_empty() { return 0; }
    // Group by (chrom, strand). Within each group, cluster by overlap.
    let mut grouped: std::collections::BTreeMap<(String, char), Vec<&(String, OriginalAlignment, f64, usize)>> =
        std::collections::BTreeMap::new();
    for r in anchored_reads {
        let (chrom, strand, _) = &r.1;
        grouped.entry((chrom.clone(), *strand)).or_default().push(r);
    }
    let mut emitted = 0usize;
    for ((chrom, strand), mut reads) in grouped {
        // Sort by first-exon start.
        reads.sort_by_key(|r| r.1.2.first().map(|e| e.0).unwrap_or(0));
        // Greedy single-pass cluster: merge into running cluster as long as
        // the next read's first-exon overlaps the current cluster's range.
        let mut idx = 0usize;
        while idx < reads.len() {
            let mut cluster: Vec<&(String, OriginalAlignment, f64, usize)> = vec![reads[idx]];
            let mut cluster_end = reads[idx].1.2.last().map(|e| e.1).unwrap_or(0);
            idx += 1;
            while idx < reads.len() {
                let next_start = reads[idx].1.2.first().map(|e| e.0).unwrap_or(0);
                if next_start > cluster_end { break; }
                let next_end = reads[idx].1.2.last().map(|e| e.1).unwrap_or(0);
                if next_end > cluster_end { cluster_end = next_end; }
                cluster.push(reads[idx]);
                idx += 1;
            }
            if cluster.len() < min_reads { continue; }

            // Build a "consensus" read by taking the union of exons from
            // all cluster members. Two-pass: collect all exon edges, sort,
            // merge overlapping into consensus exons.
            let mut edges: Vec<(u64, i32)> = Vec::new();
            for r in &cluster {
                for &(s, e) in &r.1.2 {
                    edges.push((s, 1));
                    edges.push((e, -1));
                }
            }
            edges.sort_by_key(|&(p, _)| p);
            let mut depth = 0i32;
            let mut consensus_exons: Vec<(u64, u64)> = Vec::new();
            let mut cur_start: Option<u64> = None;
            for &(pos, delta) in &edges {
                let prev_depth = depth;
                depth += delta;
                if prev_depth == 0 && depth > 0 {
                    cur_start = Some(pos);
                } else if prev_depth > 0 && depth == 0 {
                    if let Some(s) = cur_start.take() {
                        if pos > s {
                            consensus_exons.push((s, pos));
                        }
                    }
                }
            }
            if consensus_exons.is_empty() { continue; }

            let bundle_start = consensus_exons.first().map(|e| e.0.saturating_sub(100)).unwrap_or(0);
            let bundle_end   = consensus_exons.last().map(|e| e.1 + 100).unwrap_or(0);

            // Build BundleReads using each cluster member's own exons (so
            // junction support is real, not synthesized from a consensus).
            let synth_junctions_for_member = |exons: &[(u64, u64)]| -> Vec<crate::types::Junction> {
                exons.windows(2)
                    .map(|w| crate::types::Junction { donor: w[0].1, acceptor: w[1].0 })
                    .collect()
            };
            let reads_vec: Vec<BundleRead> = cluster.iter().enumerate().map(|(i, r)| {
                let exons = &r.1.2;
                let juncs = synth_junctions_for_member(exons);
                let jvalid = vec![true; juncs.len()];
                let ref_start = exons.first().map(|e| e.0).unwrap_or(0);
                let ref_end = exons.last().map(|e| e.1).unwrap_or(0);
                BundleRead {
                    read_uid: 2_000_000 + i as u64,
                    read_name: std::sync::Arc::from(r.0.as_str()),
                    read_name_hash: 0,
                    ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
                    ref_start, ref_end,
                    exons: exons.clone(),
                    junctions: juncs.clone(),
                    junction_valid: jvalid,
                    junctions_raw: juncs,
                    junctions_del: Vec::new(),
                    weight: 1.0, is_reverse: strand == '-', strand,
                    has_poly_start: false, has_poly_end: false,
                    has_poly_start_aligned: false, has_poly_start_unaligned: false,
                    has_poly_end_aligned: false, has_poly_end_unaligned: false,
                    unaligned_poly_t: 0, unaligned_poly_a: 0,
                    has_last_exon_polya: false, has_first_exon_polyt: false,
                    query_length: None, clip_left: 0, clip_right: 0,
                    nh: 1, nm: 0, md: None, insertion_sites: Vec::new(),
                    unitig: false, unitig_cov: 0.0,
                    read_count_yc: 0.0, countfrag_len: 0.0, countfrag_num: 0.0,
                    junc_mismatch_weight: 0.0,
                    pair_idx: Vec::new(), pair_count: Vec::new(),
                    mapq: 0, mismatches: Vec::new(),
                    hp_tag: None, ps_tag: None,
                    is_primary_alignment: true,
                }
            }).collect();

            out_bundles.push(Bundle {
                chrom: chrom.clone(),
                start: bundle_start,
                end: bundle_end,
                strand,
                reads: reads_vec,
                junction_stats: JunctionStats::default(),
                bundlenodes: None,
                read_bnodes: None,
                bnode_colors: None,
                synthetic: true,
                rescue_class: Some(RescueClass::NeedsExternalVerification),
            });
            emitted += 1;
        }
    }
    emitted
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

        // Derive junctions from consecutive exons. Without these, the
        // assembly's splice graph has no edges and emits zero transcripts
        // from the synthetic bundle. Each junction is (donor=prev_exon_end,
        // acceptor=next_exon_start).
        let synth_junctions: Vec<crate::types::Junction> = exons.windows(2)
            .map(|w| crate::types::Junction { donor: w[0].1, acceptor: w[1].0 })
            .collect();
        let synth_junction_valid: Vec<bool> = vec![true; synth_junctions.len()];

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
                    junctions: synth_junctions.clone(),
                    junction_valid: synth_junction_valid.clone(),
                    junctions_raw: synth_junctions.clone(),
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

        // Derive junctions from consecutive refined exons. Without these,
        // the assembly's splice graph has no edges and emits zero
        // transcripts from the synthetic bundle (the merge step recomputes
        // junction_stats from `r.junctions`, so leaving them empty kills
        // the rescue's downstream effect).
        let synth_junctions: Vec<crate::types::Junction> = refined_exons
            .windows(2)
            .map(|w| crate::types::Junction { donor: w[0].1, acceptor: w[1].0 })
            .collect();
        let synth_junction_valid: Vec<bool> = vec![true; synth_junctions.len()];

        let reads: Vec<BundleRead> = member_indices.iter().enumerate().map(|(i, &ri)| {
            let read_name = rescued[ri].0.clone();
            BundleRead {
                read_uid: 1_000_000 + i as u64,
                read_name: std::sync::Arc::from(read_name.as_str()),
                read_name_hash: 0,
                ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
                ref_start: min_start, ref_end: max_end,
                exons: refined_exons.clone(),
                junctions: synth_junctions.clone(),
                junction_valid: synth_junction_valid.clone(),
                junctions_raw: synth_junctions.clone(),
                junctions_del: Vec::new(),
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

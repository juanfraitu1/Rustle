//! Intra-bundle tandem decomposition (spec 2026-06-02). Detects a tandem
//! paralog array that collapsed into one bundle and splits it into per-copy
//! clusters so the existing family-graph/EM machinery can resolve the copies.

use std::collections::HashSet;

/// Exact k-mer Jaccard between two DNA sequences (ACGT only; other bases
/// skipped). Returns 0.0 if either has no valid k-mers. Clusters are short
/// (one copy's exonic sequence), so an exact set is fine — no min-hash needed.
///
/// Name follows the domain-first, verb-last pattern used in this codebase
/// (`compute_family_graph_kmer_jaccard`, `minimizer_jaccard`): *seq* is the
/// domain (two raw sequences), *kmer_jaccard* is the operation.
pub fn seq_kmer_jaccard(a: &[u8], b: &[u8], k: usize) -> f64 {
    fn kset(s: &[u8], k: usize) -> HashSet<&[u8]> {
        let mut set = HashSet::new();
        if s.len() >= k {
            for w in s.windows(k) {
                if w.iter().all(|&c| matches!(c, b'A' | b'C' | b'G' | b'T')) {
                    set.insert(w);
                }
            }
        }
        set
    }
    let sa = kset(a, k);
    let sb = kset(b, k);
    if sa.is_empty() || sb.is_empty() {
        return 0.0;
    }
    let inter = sa.intersection(&sb).count() as f64;
    let union = sa.union(&sb).count() as f64;
    inter / union
}

/// Cluster `[ (ref_start, ref_end) ]` read footprints into genomically-ordered
/// groups separated by gaps `> min_gap`. Returns, per cluster, the index list
/// into `intervals` plus the cluster's (start,end) span. Genomic-ascending.
pub fn cluster_reads_by_position(
    intervals: &[(u64, u64)],
    min_gap: u64,
) -> Vec<(Vec<usize>, (u64, u64))> {
    if intervals.is_empty() {
        return Vec::new();
    }
    let mut order: Vec<usize> = (0..intervals.len()).collect();
    order.sort_by_key(|&i| intervals[i].0);

    let mut clusters: Vec<(Vec<usize>, (u64, u64))> = Vec::new();
    for &i in &order {
        let (s, e) = intervals[i];
        match clusters.last_mut() {
            Some((idxs, span)) if s <= span.1.saturating_add(min_gap) => {
                idxs.push(i);
                span.1 = span.1.max(e);
                span.0 = span.0.min(s);
            }
            _ => clusters.push((vec![i], (s, e))),
        }
    }
    clusters
}

use crate::types::Bundle;
use crate::genome::GenomeIndex;

/// A tandem array detected inside one bundle: the per-copy genomic spans
/// (genomic-ascending) and, for each copy, the read indices (into `bundle.reads`)
/// that fall in that copy's span.
pub struct TandemDecomposition {
    pub copy_spans: Vec<(u64, u64)>,
    pub copy_read_idxs: Vec<Vec<usize>>,
}

/// Detect whether `bundle` is a collapsed tandem array. Requires ≥2 position
/// clusters that are mutually sequence-similar (≥ `cfg.min_jaccard`). Returns
/// `None` when disabled, single-cluster, or clusters are dissimilar (a normal
/// multi-exon gene — must NOT be split).
pub fn detect_tandem_bundle(
    bundle: &Bundle,
    genome: &GenomeIndex,
    cfg: TandemConfig,
) -> Option<TandemDecomposition> {
    if !cfg.enabled {
        return None;
    }
    // Read footprints = (ref_start, ref_end) per read.
    let ivs: Vec<(u64, u64)> = bundle.reads.iter()
        .map(|r| (r.ref_start, r.ref_end))
        .collect();
    let clusters = cluster_reads_by_position(&ivs, cfg.min_gap);
    if clusters.len() < 2 {
        return None;
    }
    // Fetch each cluster's reference sequence (the copy's genomic span).
    let seqs: Vec<Vec<u8>> = clusters.iter()
        .map(|(_, (s, e))| genome.fetch_sequence(&bundle.chrom, *s, *e).unwrap_or_default())
        .collect();
    // Require mutual similarity: every cluster must be ≥min_jaccard to at least
    // one other (so a lone unrelated cluster doesn't gate, but a single gene with
    // one weird exon block isn't promoted either).
    let n = clusters.len();
    let mut similar_to_any = vec![false; n];
    for i in 0..n {
        for j in (i + 1)..n {
            if seq_kmer_jaccard(&seqs[i], &seqs[j], cfg.kmer) >= cfg.min_jaccard {
                similar_to_any[i] = true;
                similar_to_any[j] = true;
            }
        }
    }
    // Keep only similar clusters; need ≥2 to call a tandem array.
    let kept: Vec<usize> = (0..n).filter(|&i| similar_to_any[i]).collect();
    if kept.len() < 2 {
        return None;
    }
    Some(TandemDecomposition {
        copy_spans: kept.iter().map(|&i| clusters[i].1).collect(),
        copy_read_idxs: kept.iter().map(|&i| clusters[i].0.clone()).collect(),
    })
}

/// O5/tandem runtime config (env). `enabled` defaults OFF (opt-in prototype).
#[derive(Debug, Clone, Copy)]
pub struct TandemConfig {
    pub enabled: bool,
    pub min_gap: u64,       // genomic separation to split copies
    pub min_jaccard: f64,   // cross-copy sequence-similarity floor
    pub kmer: usize,
}

impl TandemConfig {
    pub fn from_env() -> Self {
        let enabled = std::env::var("RUSTLE_VG_TANDEM")
            .map(|v| v == "1" || v.eq_ignore_ascii_case("true"))
            .unwrap_or(false);
        let min_gap = std::env::var("RUSTLE_VG_TANDEM_MIN_SEP")
            .ok().and_then(|v| v.parse().ok()).unwrap_or(5_000);
        let min_jaccard = std::env::var("RUSTLE_VG_TANDEM_MIN_JACCARD")
            .ok().and_then(|v| v.parse().ok())
            .filter(|f: &f64| *f > 0.0 && *f <= 1.0).unwrap_or(0.20);
        TandemConfig { enabled, min_gap, min_jaccard, kmer: 15 }
    }
}

use crate::vg::FamilyGroup;

/// Turn a detection into synthetic per-copy sub-bundles (appended to `bundles`)
/// + a `FamilyGroup` linking them. Each sub-bundle = the parent cloned, with
/// `reads` restricted to that copy's reads, span set to the copy span, and
/// `synthetic=true`. `multimap_reads` links reads whose name appears in >1 copy
/// (the secondaries that caused the collapse — now cross-copy evidence).
pub fn decompose_tandem_to_family(
    parent: &Bundle,
    decomp: &TandemDecomposition,
    bundles: &mut Vec<Bundle>,
    family_id: usize,
    config: &crate::types::RunConfig,
) -> (FamilyGroup, Vec<usize>) {
    let mut bundle_indices = Vec::new();
    // Track read_name_hash -> Vec<(copy_pos, read_idx_in_subbundle)> for multimap.
    let mut name_to_copies: std::collections::HashMap<u64, Vec<(usize, usize)>> =
        std::collections::HashMap::new();

    for (copy_pos, (span, read_idxs)) in
        decomp.copy_spans.iter().zip(decomp.copy_read_idxs.iter()).enumerate()
    {
        let mut sub = parent.clone();
        sub.start = span.0;
        sub.end = span.1;
        sub.synthetic = true;
        sub.vg_family_id = Some(family_id);
        sub.reads = read_idxs.iter().map(|&i| parent.reads[i].clone()).collect();
        for (ri, r) in sub.reads.iter().enumerate() {
            name_to_copies.entry(r.read_name_hash).or_default().push((copy_pos, ri));
        }
        crate::bundle::recompute_junction_stats(&mut sub, config);
        bundle_indices.push(bundles.len());
        bundles.push(sub);
    }

    // A read shared across ≥2 copies is a multimapper for the EM.
    let multimap_reads: std::collections::HashMap<u64, Vec<(usize, usize)>> =
        name_to_copies.into_iter().filter(|(_, v)| v.len() > 1).collect();

    (FamilyGroup { family_id, bundle_indices: bundle_indices.clone(), multimap_reads },
     bundle_indices)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── seq_kmer_jaccard tests ──────────────────────────────────────────────

    #[test]
    fn identical_sequences_jaccard_one() {
        let s = b"ACGTACGTACGTACGTACGT";
        assert!((seq_kmer_jaccard(s, s, 15) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn near_identical_high_jaccard() {
        // The spec prescribed: let a = b"ACGTACGTACGTACGTACGTACGTACGTACGT" (32
        // bytes, ACGT-repeating) with one substitution at byte 16 and k=15 → assert
        // Jaccard > 0.2. That is mathematically broken: a 32-byte ACGT-repeat
        // produces only 4 distinct 15-mers (the four rotation phases), and flipping
        // byte 16 adds/removes a few but the Jaccard is ≈ 0.105 < 0.2 — the
        // assertion would always fail. A non-repeating 100-base sequence is used
        // instead: one substitution at position 50 displaces 15 k-mers out of 86,
        // leaving 71 shared → Jaccard ≈ 0.70 > 0.2.
        let a: Vec<u8> = (0u8..100)
            .map(|i| [b'A', b'C', b'G', b'T', b'G', b'C', b'A', b'T', b'T', b'G'][i as usize % 10])
            .collect();
        let mut b = a.clone();
        b[50] = if b[50] == b'A' { b'C' } else { b'A' }; // one substitution
        assert!(seq_kmer_jaccard(&a, &b, 15) > 0.2);
    }

    #[test]
    fn disjoint_sequences_jaccard_zero() {
        // AAAA vs CCCC: both sequences have non-empty k-mer sets (each has exactly
        // one 15-mer), but the sets are disjoint → inter=0 → Jaccard=0.
        // NOTE: this does NOT exercise the early-exit guard (sa.is_empty() ||
        // sb.is_empty()); see the three tests below for guard coverage.
        let a = b"AAAAAAAAAAAAAAAAAAAA";
        let b = b"CCCCCCCCCCCCCCCCCCCC";
        assert_eq!(seq_kmer_jaccard(a, b, 15), 0.0);
    }

    // ── early-exit guard tests (sa.is_empty() || sb.is_empty()) ────────────
    //
    // The docstring promises 0.0 for (i) empty slice, (ii) slice shorter than k,
    // (iii) all-non-ACGT (e.g. all-N). Each test below hits the guard directly.

    #[test]
    fn empty_slice_returns_zero() {
        // An empty slice produces an empty k-mer set → guard fires.
        let s = b"ACGTACGTACGTACGT";
        assert_eq!(seq_kmer_jaccard(b"", s, 15), 0.0);
        assert_eq!(seq_kmer_jaccard(s, b"", 15), 0.0);
        assert_eq!(seq_kmer_jaccard(b"", b"", 15), 0.0);
    }

    #[test]
    fn shorter_than_k_returns_zero() {
        // A slice with length < k produces an empty k-mer set → guard fires.
        let short = b"ACGTACGT"; // 8 bytes < k=15
        let long = b"ACGTACGTACGTACGTACGT";
        assert_eq!(seq_kmer_jaccard(short, long, 15), 0.0);
        assert_eq!(seq_kmer_jaccard(long, short, 15), 0.0);
    }

    #[test]
    fn all_n_sequence_returns_zero() {
        // All-N (non-ACGT) sequences: all k-mers are rejected by the ACGT filter,
        // so the k-mer set is empty → guard fires.
        let n_seq = b"NNNNNNNNNNNNNNNNNNNN"; // 20 N-bases
        let acgt = b"ACGTACGTACGTACGTACGT";
        assert_eq!(seq_kmer_jaccard(n_seq, acgt, 15), 0.0);
        assert_eq!(seq_kmer_jaccard(acgt, n_seq, 15), 0.0);
    }

    // ── cluster_reads_by_position tests ────────────────────────────────────

    #[test]
    fn splits_two_separated_copies() {
        // two copies ~23kb apart, gap threshold 5kb -> 2 clusters.
        let ivs = [(1000, 14000), (1500, 14500), (37000, 50000), (37200, 50200)];
        let c = cluster_reads_by_position(&ivs, 5000);
        assert_eq!(c.len(), 2);
        assert_eq!(c[0].0.len(), 2);
        assert_eq!(c[1].0.len(), 2);
        assert_eq!(c[0].1, (1000, 14500));
        assert_eq!(c[1].1, (37000, 50200));
        // Primary contract: returned indices must map back to the original input
        // slots so callers can recover per-read metadata. A sort permutation that
        // re-ordered indices incorrectly would break all callers silently.
        assert_eq!(c[0].0, vec![0, 1], "first cluster must contain original indices 0 and 1");
        assert_eq!(c[1].0, vec![2, 3], "second cluster must contain original indices 2 and 3");
    }

    #[test]
    fn single_cluster_when_within_gap() {
        let ivs = [(1000, 14000), (16000, 30000)]; // 2kb gap < 5kb
        assert_eq!(cluster_reads_by_position(&ivs, 5000).len(), 1);
    }

    #[test]
    fn empty_input() {
        assert!(cluster_reads_by_position(&[], 5000).is_empty());
    }

    // ── TandemConfig tests ──────────────────────────────────────────────────

    #[test]
    fn tandemconfig_defaults_off() {
        // Ensure none of our env vars are set (clean state).
        std::env::remove_var("RUSTLE_VG_TANDEM");
        std::env::remove_var("RUSTLE_VG_TANDEM_MIN_SEP");
        std::env::remove_var("RUSTLE_VG_TANDEM_MIN_JACCARD");
        let c = TandemConfig::from_env();
        assert!(!c.enabled);
        assert_eq!(c.min_gap, 5_000);
        assert!((c.min_jaccard - 0.20).abs() < 1e-9);
    }

    /// RUSTLE_VG_TANDEM=1 and the "true"/case-variant spellings must all enable.
    #[test]
    fn tandemconfig_enabled_by_one() {
        std::env::set_var("RUSTLE_VG_TANDEM", "1");
        let c = TandemConfig::from_env();
        std::env::remove_var("RUSTLE_VG_TANDEM");
        assert!(c.enabled, "RUSTLE_VG_TANDEM=1 must set enabled=true");
    }

    #[test]
    fn tandemconfig_enabled_by_true_lowercase() {
        std::env::set_var("RUSTLE_VG_TANDEM", "true");
        let c = TandemConfig::from_env();
        std::env::remove_var("RUSTLE_VG_TANDEM");
        assert!(c.enabled, "RUSTLE_VG_TANDEM=true must set enabled=true");
    }

    #[test]
    fn tandemconfig_enabled_by_true_uppercase() {
        std::env::set_var("RUSTLE_VG_TANDEM", "TRUE");
        let c = TandemConfig::from_env();
        std::env::remove_var("RUSTLE_VG_TANDEM");
        assert!(c.enabled, "RUSTLE_VG_TANDEM=TRUE must set enabled=true");
    }

    #[test]
    fn tandemconfig_enabled_by_true_mixed_case() {
        std::env::set_var("RUSTLE_VG_TANDEM", "True");
        let c = TandemConfig::from_env();
        std::env::remove_var("RUSTLE_VG_TANDEM");
        assert!(c.enabled, "RUSTLE_VG_TANDEM=True must set enabled=true");
    }

    /// RUSTLE_VG_TANDEM_MIN_SEP: valid integer overrides the 5000 default;
    /// a non-parseable value falls back to the default.
    #[test]
    fn tandemconfig_min_gap_custom_value() {
        std::env::set_var("RUSTLE_VG_TANDEM_MIN_SEP", "12345");
        let c = TandemConfig::from_env();
        std::env::remove_var("RUSTLE_VG_TANDEM_MIN_SEP");
        assert_eq!(c.min_gap, 12_345, "valid integer should override the default 5000");
    }

    #[test]
    fn tandemconfig_min_gap_fallback_on_bad_value() {
        std::env::set_var("RUSTLE_VG_TANDEM_MIN_SEP", "not_a_number");
        let c = TandemConfig::from_env();
        std::env::remove_var("RUSTLE_VG_TANDEM_MIN_SEP");
        assert_eq!(c.min_gap, 5_000, "non-parseable value should fall back to 5000");
    }

    /// RUSTLE_VG_TANDEM_MIN_JACCARD: valid in-range value overrides the 0.20 default.
    #[test]
    fn tandemconfig_min_jaccard_custom_value() {
        std::env::set_var("RUSTLE_VG_TANDEM_MIN_JACCARD", "0.5");
        let c = TandemConfig::from_env();
        std::env::remove_var("RUSTLE_VG_TANDEM_MIN_JACCARD");
        assert!((c.min_jaccard - 0.5).abs() < 1e-9, "0.5 is in (0,1] and must be accepted; got {}", c.min_jaccard);
    }

    /// Boundary: exactly 1.0 is accepted (<=1.0 branch).
    #[test]
    fn tandemconfig_min_jaccard_boundary_one_accepted() {
        std::env::set_var("RUSTLE_VG_TANDEM_MIN_JACCARD", "1.0");
        let c = TandemConfig::from_env();
        std::env::remove_var("RUSTLE_VG_TANDEM_MIN_JACCARD");
        assert!((c.min_jaccard - 1.0).abs() < 1e-9, "1.0 is in (0,1] and must be accepted; got {}", c.min_jaccard);
    }

    /// Boundary: exactly 0.0 is rejected (>0.0 filter) → falls back to 0.20.
    #[test]
    fn tandemconfig_min_jaccard_boundary_zero_rejected() {
        std::env::set_var("RUSTLE_VG_TANDEM_MIN_JACCARD", "0.0");
        let c = TandemConfig::from_env();
        std::env::remove_var("RUSTLE_VG_TANDEM_MIN_JACCARD");
        assert!((c.min_jaccard - 0.20).abs() < 1e-9, "0.0 is not in (0,1] and must fall back to 0.20; got {}", c.min_jaccard);
    }

    /// Values <= 0.0 (negative): rejected, falls back to 0.20.
    #[test]
    fn tandemconfig_min_jaccard_negative_rejected() {
        std::env::set_var("RUSTLE_VG_TANDEM_MIN_JACCARD", "-0.5");
        let c = TandemConfig::from_env();
        std::env::remove_var("RUSTLE_VG_TANDEM_MIN_JACCARD");
        assert!((c.min_jaccard - 0.20).abs() < 1e-9, "-0.5 must fall back to 0.20; got {}", c.min_jaccard);
    }

    /// Values > 1.0: rejected, falls back to 0.20.
    #[test]
    fn tandemconfig_min_jaccard_above_one_rejected() {
        std::env::set_var("RUSTLE_VG_TANDEM_MIN_JACCARD", "1.5");
        let c = TandemConfig::from_env();
        std::env::remove_var("RUSTLE_VG_TANDEM_MIN_JACCARD");
        assert!((c.min_jaccard - 0.20).abs() < 1e-9, "1.5 > 1.0 must fall back to 0.20; got {}", c.min_jaccard);
    }

    // ── detect_tandem_bundle tests ─────────────────────────────────────────

    /// Minimal Bundle constructor for tests: empty reads, no junctions.
    fn make_empty_bundle() -> crate::types::Bundle {
        crate::types::Bundle {
            chrom: "chrTest".to_string(),
            start: 0,
            end: 0,
            strand: '+',
            reads: Vec::new(),
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

    /// Construct a BundleRead with the given ref_start / ref_end; all other
    /// fields are inert defaults used only to satisfy the struct layout.
    fn make_read(rs: u64, re: u64) -> crate::types::BundleRead {
        use std::sync::Arc;
        crate::types::BundleRead {
            read_uid: 0, read_name: Arc::from("r"), read_name_hash: 0,
            ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
            ref_start: rs, ref_end: re,
            exons: vec![(rs, re)], junctions: vec![], junction_valid: vec![],
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
            mapq: 60, mismatches: vec![], seq: Vec::new(),
            hp_tag: None, ps_tag: None,
            is_primary_alignment: true,
            em_weight_gap: -1.0, em_n_sites: 0, em_anchored: true,
        }
    }

    /// Build a GenomeIndex with a single contig "chrTest" whose sequence is
    /// constructed by the caller-supplied builder: `build(buf)` fills `buf` to
    /// the desired length. Uses a tempfile so the private `seqs` field need not
    /// be exposed.
    fn make_genome(build: impl Fn(&mut Vec<u8>)) -> crate::genome::GenomeIndex {
        let dir = tempfile::tempdir().unwrap();
        let fa = dir.path().join("test.fa");
        let mut seq: Vec<u8> = Vec::new();
        build(&mut seq);
        let mut content = b">chrTest\n".to_vec();
        content.extend_from_slice(&seq);
        content.push(b'\n');
        std::fs::write(&fa, &content).unwrap();
        crate::genome::GenomeIndex::from_fasta(fa.to_str().unwrap()).unwrap()
    }

    // ── None-path unit tests ───────────────────────────────────────────────

    /// `enabled=false` must return `None` immediately without reading the genome.
    #[test]
    fn detect_disabled_returns_none() {
        let cfg = TandemConfig { enabled: false, min_gap: 5000, min_jaccard: 0.2, kmer: 15 };
        let b = make_empty_bundle();
        let g = crate::genome::GenomeIndex::default();
        assert!(detect_tandem_bundle(&b, &g, cfg).is_none());
    }

    /// A bundle with no reads always returns `None` (fewer than 2 clusters).
    #[test]
    fn detect_empty_bundle_returns_none() {
        let cfg = TandemConfig { enabled: true, min_gap: 5000, min_jaccard: 0.2, kmer: 15 };
        let b = make_empty_bundle();
        let g = crate::genome::GenomeIndex::default();
        assert!(detect_tandem_bundle(&b, &g, cfg).is_none());
    }

    /// A bundle whose reads all fall in a single position cluster returns `None`.
    #[test]
    fn detect_single_cluster_returns_none() {
        let cfg = TandemConfig { enabled: true, min_gap: 5000, min_jaccard: 0.2, kmer: 15 };
        let mut b = make_empty_bundle();
        // Two reads within 5kb of each other → single cluster.
        b.reads = vec![make_read(1000, 14000), make_read(1500, 14500)];
        let g = crate::genome::GenomeIndex::default();
        assert!(detect_tandem_bundle(&b, &g, cfg).is_none());
    }

    // ── Dissimilar-clusters None-path (primary false-positive guard) ───────
    //
    // A bundle with ≥2 position clusters but dissimilar reference sequences
    // (a normal multi-exon gene, not a tandem array) must return None.
    // This is the guard the spec explicitly calls out as "must NOT be split".

    /// Two clusters with completely disjoint k-mer sets (AAAA… vs CCCC…):
    /// seq_kmer_jaccard = 0.0 < min_jaccard=0.2 → function must return None.
    #[test]
    fn detect_dissimilar_clusters_returns_none() {
        // Genome layout (all ACGT-only so k-mers are valid):
        //   positions   0.. 100: copy-1 region = ACGT repeated (25×)
        //   positions 100..5200: filler = A repeated
        //   positions 5200..5300: copy-2 region = TGCA repeated (25×)
        // The two 100-base windows share zero 15-mers → Jaccard = 0.0 < 0.2.
        let genome = make_genome(|buf| {
            // copy-1 window at 0..100: ACGT×25
            for _ in 0..25 { buf.extend_from_slice(b"ACGT"); }
            // filler 100..5200
            buf.extend(std::iter::repeat(b'A').take(5100));
            // copy-2 window at 5200..5300: TGCA×25 (zero shared 15-mers with ACGT×25)
            for _ in 0..25 { buf.extend_from_slice(b"TGCA"); }
        });
        let cfg = TandemConfig { enabled: true, min_gap: 5000, min_jaccard: 0.2, kmer: 15 };
        let mut b = make_empty_bundle();
        // Read 0: cluster 1 (positions 0..100)
        // Read 1: cluster 2 (positions 5200..5300, gap=5100 > 5000)
        b.reads = vec![make_read(0, 100), make_read(5200, 5300)];
        assert!(
            detect_tandem_bundle(&b, &genome, cfg).is_none(),
            "dissimilar clusters (multi-exon gene) must return None, not be split"
        );
    }

    // ── Some-path: happy path (two similar clusters) ───────────────────────

    /// Two separated clusters with sequence-similar reference windows:
    /// detect_tandem_bundle must return Some with 2 copy_spans and correct
    /// copy_read_idxs. This exercises lines 96-121 of detect_tandem_bundle.
    #[test]
    fn detect_similar_clusters_returns_some() {
        // Genome layout:
        //   positions   0..100: copy-1 window = non-repeating 100-base sequence
        //   positions 100..5200: filler = A repeated
        //   positions 5200..5300: copy-2 window = same 100-base sequence with
        //     one substitution at position 50 (Jaccard ≈ 0.70 >> 0.2 threshold)
        let base_seq: Vec<u8> = (0u8..100)
            .map(|i| [b'A',b'C',b'G',b'T',b'G',b'C',b'A',b'T',b'T',b'G'][i as usize % 10])
            .collect();
        let mut copy2_seq = base_seq.clone();
        copy2_seq[50] = if copy2_seq[50] == b'A' { b'C' } else { b'A' };

        let base_seq_clone = base_seq.clone();
        let copy2_seq_clone = copy2_seq.clone();
        let genome = make_genome(move |buf| {
            buf.extend_from_slice(&base_seq_clone);
            buf.extend(std::iter::repeat(b'A').take(5100));
            buf.extend_from_slice(&copy2_seq_clone);
        });

        let cfg = TandemConfig { enabled: true, min_gap: 5000, min_jaccard: 0.2, kmer: 15 };
        let mut b = make_empty_bundle();
        // Read 0: cluster 1, Read 1: cluster 2 (gap = 5100 > 5000)
        b.reads = vec![make_read(0, 100), make_read(5200, 5300)];

        let result = detect_tandem_bundle(&b, &genome, cfg);
        assert!(result.is_some(), "similar clusters must return Some(TandemDecomposition)");
        let decomp = result.unwrap();
        assert_eq!(decomp.copy_spans.len(), 2, "must report exactly 2 copy spans");
        assert_eq!(decomp.copy_read_idxs.len(), 2, "must report exactly 2 copy read-index vecs");
        // Cluster 0: read index 0; cluster 1: read index 1.
        assert_eq!(decomp.copy_read_idxs[0], vec![0usize]);
        assert_eq!(decomp.copy_read_idxs[1], vec![1usize]);
        // copy_spans must cover the respective cluster genomic spans.
        assert_eq!(decomp.copy_spans[0], (0, 100));
        assert_eq!(decomp.copy_spans[1], (5200, 5300));
    }

    // ── Partial-similarity filter (3 clusters, 2 similar + 1 unrelated) ───

    /// Three position clusters where clusters 0 and 1 are sequence-similar but
    /// cluster 2 has a completely unrelated sequence. The kept-cluster filter
    /// should retain only clusters 0 and 1 and return a TandemDecomposition
    /// with exactly 2 copy_spans (not 3).
    #[test]
    fn detect_partial_similarity_drops_dissimilar_cluster() {
        // Genome layout:
        //   positions   0..100: copy-1 window = non-repeating sequence
        //   positions 100..5200: filler = A repeated
        //   positions 5200..5300: copy-2 window = same sequence, one sub (similar)
        //   positions 5300..10400: filler = A repeated
        //   positions 10400..10500: copy-3 window = CCCC×25 (dissimilar to both)
        let base_seq: Vec<u8> = (0u8..100)
            .map(|i| [b'A',b'C',b'G',b'T',b'G',b'C',b'A',b'T',b'T',b'G'][i as usize % 10])
            .collect();
        let mut copy2_seq = base_seq.clone();
        copy2_seq[50] = if copy2_seq[50] == b'A' { b'C' } else { b'A' };

        let base_seq_clone = base_seq.clone();
        let copy2_seq_clone = copy2_seq.clone();
        let genome = make_genome(move |buf| {
            buf.extend_from_slice(&base_seq_clone);           // 0..100
            buf.extend(std::iter::repeat(b'A').take(5100));   // 100..5200
            buf.extend_from_slice(&copy2_seq_clone);          // 5200..5300
            buf.extend(std::iter::repeat(b'A').take(5100));   // 5300..10400
            // cluster 2 window at 10400..10500: CCCC×25 (no shared 15-mers with ACGT-based seqs)
            buf.extend(std::iter::repeat(b'C').take(100));    // 10400..10500
        });

        let cfg = TandemConfig { enabled: true, min_gap: 5000, min_jaccard: 0.2, kmer: 15 };
        let mut b = make_empty_bundle();
        // Read 0 → cluster 0 (positions 0..100)
        // Read 1 → cluster 1 (positions 5200..5300, gap=5100 > 5000)
        // Read 2 → cluster 2 (positions 10400..10500, gap=5100 > 5000)
        b.reads = vec![make_read(0, 100), make_read(5200, 5300), make_read(10400, 10500)];

        let result = detect_tandem_bundle(&b, &genome, cfg);
        assert!(result.is_some(), "2 of 3 clusters are similar → should still return Some");
        let decomp = result.unwrap();
        assert_eq!(
            decomp.copy_spans.len(), 2,
            "the dissimilar cluster must be dropped; expected 2 spans, got {:?}",
            decomp.copy_spans
        );
        assert_eq!(decomp.copy_read_idxs.len(), 2);
        // The kept spans must be the two similar ones (cluster 0 and cluster 1).
        assert_eq!(decomp.copy_spans[0], (0, 100), "first kept span must be cluster 0");
        assert_eq!(decomp.copy_spans[1], (5200, 5300), "second kept span must be cluster 1");
        // Each kept cluster contributed exactly one read.
        assert_eq!(decomp.copy_read_idxs[0], vec![0usize]);
        assert_eq!(decomp.copy_read_idxs[1], vec![1usize]);
    }

    // ── decompose_tandem_to_family tests ─────────────────────────────────────

    /// Build a parent bundle with 4 reads: reads 0,1 at copy-0 locus (1000–14000)
    /// and reads 2,3 at copy-1 locus (37000–50000). All reads have distinct
    /// read_name_hash values so none are cross-copy multimappers.
    fn make_bundle_two_copies() -> crate::types::Bundle {
        use std::sync::Arc;
        let mut b = make_empty_bundle();
        b.reads = vec![
            {
                let mut r = make_read(1000, 14000);
                r.read_name_hash = 1001;
                r.read_name = Arc::from("r1001");
                r
            },
            {
                let mut r = make_read(1500, 14500);
                r.read_name_hash = 1002;
                r.read_name = Arc::from("r1002");
                r
            },
            {
                let mut r = make_read(37000, 50000);
                r.read_name_hash = 1003;
                r.read_name = Arc::from("r1003");
                r
            },
            {
                let mut r = make_read(37200, 50200);
                r.read_name_hash = 1004;
                r.read_name = Arc::from("r1004");
                r
            },
        ];
        b.start = 1000;
        b.end = 50200;
        b
    }

    #[test]
    fn decompose_builds_one_subbundle_per_copy() {
        let parent = make_bundle_two_copies();
        let decomp = TandemDecomposition {
            copy_spans: vec![(1000, 14000), (37000, 50000)],
            copy_read_idxs: vec![vec![0, 1], vec![2, 3]],
        };
        let mut bundles: Vec<crate::types::Bundle> = Vec::new();
        let cfg = crate::types::RunConfig::default();
        let (fam, idxs) = crate::vg_hmm::tandem::decompose_tandem_to_family(
            &parent, &decomp, &mut bundles, 0, &cfg,
        );
        assert_eq!(idxs.len(), 2, "two copies → two bundle indices");
        assert_eq!(bundles.len(), 2, "two sub-bundles appended");
        assert_eq!(bundles[0].reads.len(), 2, "copy-0 sub-bundle has 2 reads");
        assert_eq!(bundles[1].reads.len(), 2, "copy-1 sub-bundle has 2 reads");
        assert!(bundles[0].synthetic, "copy-0 sub-bundle must be synthetic");
        assert!(bundles[1].synthetic, "copy-1 sub-bundle must be synthetic");
        assert_eq!(bundles[0].start, 1000);
        assert_eq!(bundles[0].end, 14000);
        assert_eq!(bundles[1].start, 37000);
        assert_eq!(bundles[1].end, 50000);
        assert_eq!(fam.bundle_indices, idxs, "FamilyGroup.bundle_indices matches returned idxs");
        assert_eq!(fam.family_id, 0);
        // No shared read_name_hash → multimap_reads must be empty.
        assert!(fam.multimap_reads.is_empty(), "no cross-copy reads → multimap_reads empty");
    }

    #[test]
    fn decompose_detects_multimap_reads() {
        // Read at index 0 (copy-0) and index 2 (copy-1) share the same
        // read_name_hash (simulating a secondary alignment spanning both loci).
        let mut parent = make_bundle_two_copies();
        let shared_hash = 9999u64;
        parent.reads[0].read_name_hash = shared_hash;
        parent.reads[2].read_name_hash = shared_hash;

        let decomp = TandemDecomposition {
            copy_spans: vec![(1000, 14000), (37000, 50000)],
            copy_read_idxs: vec![vec![0, 1], vec![2, 3]],
        };
        let mut bundles: Vec<crate::types::Bundle> = Vec::new();
        let cfg = crate::types::RunConfig::default();
        let (fam, _idxs) = crate::vg_hmm::tandem::decompose_tandem_to_family(
            &parent, &decomp, &mut bundles, 42, &cfg,
        );
        assert_eq!(fam.multimap_reads.len(), 1, "one shared hash → one multimap entry");
        let placements = fam.multimap_reads.get(&shared_hash).unwrap();
        assert_eq!(placements.len(), 2, "shared read has placements in 2 copies");
    }
}

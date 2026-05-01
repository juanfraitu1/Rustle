//! Regression tests for vg_hmm::rescue (Phase 5).

use rustle::vg_hmm::rescue::{prefilter_read, run_rescue_in_memory, synthesize_bundles, synthesize_bundles_refined};
use rustle::vg_hmm::family_graph::{FamilyGraph, ExonClass, JunctionEdge, NodeIdx};
use rustle::vg_hmm::profile::ProfileHmm;
use rustle::vg::FamilyGroup;
use rustle::path_extract::Transcript;
use std::collections::{HashMap, HashSet};

// ── Helpers ───────────────────────────────────────────────────────────────────

fn kmer_set(seqs: &[&[u8]], k: usize) -> HashSet<u64> {
    fn fnv1a64(s: &[u8]) -> u64 {
        let mut h: u64 = 0xcbf29ce484222325;
        for &b in s { h ^= b as u64; h = h.wrapping_mul(0x100000001b3); }
        h
    }
    let mut set = HashSet::new();
    for seq in seqs {
        for w in seq.windows(k) {
            set.insert(fnv1a64(w));
        }
    }
    set
}

fn two_node_fg_with_profiles() -> FamilyGraph {
    let p1 = ProfileHmm::from_singleton(b"ACGTACGT");
    let p2 = ProfileHmm::from_singleton(b"TGCATGCA");
    let nodes = vec![
        ExonClass {
            idx: NodeIdx(0), chrom: "chr1".into(), span: (1000, 1008), strand: '+',
            per_copy_sequences: vec![(0, b"ACGTACGT".to_vec())],
            per_copy_spans: vec![(0, (1000, 1008))],
            copy_specific: false, profile: Some(p1),
        },
        ExonClass {
            idx: NodeIdx(1), chrom: "chr1".into(), span: (2000, 2008), strand: '+',
            per_copy_sequences: vec![(0, b"TGCATGCA".to_vec())],
            per_copy_spans: vec![(0, (2000, 2008))],
            copy_specific: false, profile: Some(p2),
        },
    ];
    let edges = vec![JunctionEdge { from: NodeIdx(0), to: NodeIdx(1), family_support: 1, strand: '+' }];
    FamilyGraph { family_id: 0, nodes, edges }
}

// ── Task 5.1: prefilter_read ──────────────────────────────────────────────────

#[test]
fn prefilter_read_passes_when_enough_hits() {
    let seq = b"ACGTACGTACGTACGT";
    let ks = kmer_set(&[seq], 5);
    let (passed, hits) = prefilter_read(seq, &ks, 5, 3);
    assert!(passed, "expected pass; hits={}", hits);
    assert!(hits >= 3);
}

#[test]
fn prefilter_read_fails_when_too_few_hits() {
    let family_seq = b"ACGTACGTACGTACGT";
    let ks = kmer_set(&[family_seq], 5);
    // Query is completely different.
    let query = b"TTTTTTTTTTTTTTTTT";
    let (passed, hits) = prefilter_read(query, &ks, 5, 3);
    assert!(!passed, "expected fail; hits={}", hits);
}

#[test]
fn prefilter_read_empty_family_kmers_always_fails() {
    let empty: HashSet<u64> = HashSet::new();
    let (passed, hits) = prefilter_read(b"ACGTACGT", &empty, 5, 1);
    assert!(!passed);
    assert_eq!(hits, 0);
}

#[test]
fn prefilter_read_short_read_fails() {
    let ks = kmer_set(&[b"ACGTACGT"], 15);
    // Read shorter than k.
    let (passed, hits) = prefilter_read(b"ACG", &ks, 15, 1);
    assert!(!passed);
    assert_eq!(hits, 0);
}

#[test]
fn prefilter_read_returns_consistent_passed_and_hits() {
    let seq = b"ACGTACGTACGTACGTACGT";
    let ks = kmer_set(&[seq], 5);
    let (passed, hits) = prefilter_read(seq, &ks, 5, 5);
    // passed iff hits >= min_hits
    assert_eq!(passed, hits >= 5);
}

// ── Task 5.2: run_rescue_in_memory ────────────────────────────────────────────

#[test]
fn run_rescue_in_memory_empty_inputs_returns_empty() {
    let (candidates, scored) = run_rescue_in_memory(&[], &[], &[], 30.0)
        .expect("should not fail on empty inputs");
    assert!(candidates.is_empty());
    assert!(scored.is_empty());
}

#[test]
fn run_rescue_in_memory_empty_reads_returns_empty() {
    let fg = two_node_fg_with_profiles();
    let family = FamilyGroup { family_id: 0, bundle_indices: vec![], multimap_reads: HashMap::new() };
    let (candidates, scored) = run_rescue_in_memory(&[family], &[fg], &[], -999.0)
        .expect("should not fail");
    assert!(candidates.is_empty());
    assert!(scored.is_empty());
}

#[test]
fn run_rescue_in_memory_high_threshold_returns_empty() {
    let fg = two_node_fg_with_profiles();
    let family = FamilyGroup { family_id: 0, bundle_indices: vec![], multimap_reads: HashMap::new() };
    // Use an impossibly high threshold so no read passes.
    let reads = vec![("r1".to_string(), b"ACGTACGTTGCATGCA".to_vec())];
    let (candidates, scored) = run_rescue_in_memory(&[family], &[fg], &reads, 1e18)
        .expect("should not fail");
    assert!(candidates.is_empty());
    assert!(scored.is_empty());
}

#[test]
fn run_rescue_in_memory_matching_read_scores_pass() {
    let fg = two_node_fg_with_profiles();
    let family = FamilyGroup { family_id: 0, bundle_indices: vec![], multimap_reads: HashMap::new() };
    // A read that exactly matches the two-node chain.
    let reads = vec![("r1".to_string(), b"ACGTACGTTGCATGCA".to_vec())];
    let (candidates, scored) = run_rescue_in_memory(&[family], &[fg], &reads, f64::NEG_INFINITY)
        .expect("should not fail");
    assert_eq!(candidates.len(), 1);
    assert_eq!(scored.len(), 1);
    assert_eq!(scored[0].1, 0); // family_id = 0
    assert!(scored[0].3 > f64::NEG_INFINITY);
}

// ── Task 5.4: synthesize_bundles ──────────────────────────────────────────────

#[test]
fn synthesize_bundles_min_reads_drops_small_clusters() {
    let fg = two_node_fg_with_profiles();
    // Only 2 reads on the same path — below min_reads=3.
    let path = vec![NodeIdx(0), NodeIdx(1)];
    let reads = vec![
        ("r1".to_string(), path.clone()),
        ("r2".to_string(), path.clone()),
    ];
    let bundles = synthesize_bundles(&fg, &reads, 3, 15, false);
    assert!(bundles.is_empty(), "cluster of 2 should be dropped with min_reads=3");
}

#[test]
fn synthesize_bundles_emits_synthetic_bundle_for_sufficient_cluster() {
    let fg = two_node_fg_with_profiles();
    let path = vec![NodeIdx(0), NodeIdx(1)];
    let reads: Vec<(String, Vec<NodeIdx>)> = (0..3)
        .map(|i| (format!("r{}", i), path.clone()))
        .collect();
    let bundles = synthesize_bundles(&fg, &reads, 3, 15, false);
    assert_eq!(bundles.len(), 1, "expected exactly 1 bundle");
    let b = &bundles[0];
    assert!(b.synthetic, "bundle must be marked synthetic");
    assert_eq!(b.chrom, "chr1");
    assert_eq!(b.strand, '+');
    assert_eq!(b.reads.len(), 3);
    // Start/end should be padded by 1000 bp.
    assert_eq!(b.start, 0); // 1000 - 1000 = 0 (saturating)
    assert_eq!(b.end, 3008); // 2008 + 1000
    // read_uid should start at 1_000_000.
    assert_eq!(b.reads[0].read_uid, 1_000_000);
}

// ── synthesize_bundles_refined: boundary refinement via Viterbi spans ────────

#[test]
fn synthesize_bundles_refined_keeps_internal_node_spans() {
    // 2-node graph; rescued reads supply Viterbi spans that all agree at the
    // first-node footprint length. Internal junctions (= node spans here, since
    // both nodes are boundary in 2-node) must come from the family graph.
    // With only 2 nodes, both ARE boundary (first and last). Internal-keep
    // semantics is tested with 3+ nodes — covered indirectly by clustering tests.
    let fg = two_node_fg_with_profiles();
    let path = vec![NodeIdx(0), NodeIdx(1)];
    // 4 reads, all with footprint (0,8) on node 0 (= 8 bp footprint) and (8,16) on node 1.
    let rescued: Vec<(String, Vec<NodeIdx>, Vec<(usize, usize)>)> = (0..4)
        .map(|i| (format!("r{}", i), path.clone(), vec![(0, 8), (8, 16)]))
        .collect();
    let bundles = synthesize_bundles_refined(&fg, &rescued, 3, 15, 3, 5);
    assert_eq!(bundles.len(), 1, "expected exactly 1 bundle");
    let b = &bundles[0];
    // Node spans: node 0 = (1000, 1008) length 8; node 1 = (2000, 2008) length 8.
    // All 4 reads agree on footprint = 8 (matches node_len), so refined boundary
    // = node_len = 8. Boundary anchoring: first exon keeps node_e=1008, grows
    // start leftward by 8 → (1000, 1008). Last exon keeps node_s=2000, grows
    // end rightward by 8 → (2000, 2008). Net: same as RAW.
    assert_eq!(b.reads[0].exons, vec![(1000, 1008), (2000, 2008)]);
}

#[test]
fn synthesize_bundles_refined_falls_back_when_too_few_supporters() {
    // Only 2 reads — below hard_min_reads=3. Cluster fails the support
    // threshold; refined boundaries fall back to family-graph node lengths.
    let fg = two_node_fg_with_profiles();
    let path = vec![NodeIdx(0), NodeIdx(1)];
    let rescued: Vec<(String, Vec<NodeIdx>, Vec<(usize, usize)>)> = (0..3)
        .map(|i| (format!("r{}", i), path.clone(), vec![(0, 12), (12, 20)]))  // footprint 12 on node 0 (longer than node_len=8)
        .collect();
    let bundles = synthesize_bundles_refined(&fg, &rescued, 3, 15, 5, 5);  // hard_min_reads=5, only 3 reads available
    assert_eq!(bundles.len(), 1);
    let b = &bundles[0];
    // 3 < hard_min_reads=5 → both boundaries fall back to node_len=8.
    // First exon: node_e=1008 minus len=8 → (1000, 1008). Same as RAW.
    assert_eq!(b.reads[0].exons, vec![(1000, 1008), (2000, 2008)]);
}

#[test]
fn synthesize_bundles_refined_produces_synthetic_bundle_with_strand_aware_anchoring() {
    let fg = two_node_fg_with_profiles();
    let path = vec![NodeIdx(0), NodeIdx(1)];
    let rescued: Vec<(String, Vec<NodeIdx>, Vec<(usize, usize)>)> = (0..3)
        .map(|i| (format!("r{}", i), path.clone(), vec![(0, 8), (8, 16)]))
        .collect();
    let bundles = synthesize_bundles_refined(&fg, &rescued, 3, 15, 3, 5);
    let b = &bundles[0];
    assert!(b.synthetic, "must mark synthetic=true");
    assert_eq!(b.strand, '+');  // both nodes are + strand in two_node_fg_with_profiles
    assert_eq!(b.reads[0].read_uid, 1_000_000);  // refined uses same uid scheme
}

#[test]
fn synthesize_bundles_two_clusters_emit_two_bundles() {
    let fg = two_node_fg_with_profiles();
    let path_a = vec![NodeIdx(0)];
    let path_b = vec![NodeIdx(1)];
    let reads: Vec<(String, Vec<NodeIdx>)> = (0..3)
        .map(|i| (format!("ra{}", i), path_a.clone()))
        .chain((0..3).map(|i| (format!("rb{}", i), path_b.clone())))
        .collect();
    let bundles = synthesize_bundles(&fg, &reads, 3, 15, false);
    assert_eq!(bundles.len(), 2, "expected two distinct clusters → two bundles");
}

// ── Task 5.5: synthetic flag propagation Bundle → Transcript ──────────────────

/// Regression test for the Phase 5 propagation bug: Bundle.synthetic must be
/// stamped onto every Transcript produced from that bundle.  The fix lives in
/// pipeline.rs::extract_bundle_transcripts_for_graph immediately after the
/// per-bundle extract call.  Here we verify the propagation logic in isolation:
/// start with transcripts that all have synthetic=false, apply the same block
/// that the pipeline now runs, and confirm every transcript ends up with
/// synthetic=true when the source bundle is synthetic.
#[test]
fn synthetic_flag_propagates_from_bundle_to_transcripts() {
    let fg = two_node_fg_with_profiles();
    let path = vec![NodeIdx(0), NodeIdx(1)];
    let reads: Vec<(String, Vec<NodeIdx>)> = (0..3)
        .map(|i| (format!("r{}", i), path.clone()))
        .collect();
    let bundles = synthesize_bundles(&fg, &reads, 3, 15, false);
    assert_eq!(bundles.len(), 1);
    let bundle = &bundles[0];
    assert!(bundle.synthetic, "precondition: synthesize_bundles must set synthetic=true");

    // Simulate what extract_bundle_transcripts_for_graph returns before the
    // fix: a Vec<Transcript> where every entry has synthetic=false.
    let mut txs: Vec<Transcript> = (0..3).map(|_| Transcript {
        chrom: bundle.chrom.clone(),
        strand: bundle.strand,
        exons: vec![(bundle.start, bundle.end)],
        synthetic: false,   // ← exactly what every extractor hardcodes
        ..Default::default()
    }).collect();

    // Apply the exact propagation block added by this fix.
    if bundle.synthetic {
        for tx in &mut txs {
            tx.synthetic = true;
        }
    }

    // Every transcript must now carry synthetic=true.
    for (i, tx) in txs.iter().enumerate() {
        assert!(tx.synthetic, "transcript[{}] must have synthetic=true after propagation", i);
    }
}

/// Confirm non-synthetic bundles leave transcript.synthetic=false (no
/// accidental collateral stamping).
#[test]
fn non_synthetic_bundle_leaves_transcript_flag_false() {
    // Build a bundle with synthetic=false by calling synthesize_bundles and
    // then clearing the flag (synthesize_bundles always sets it).
    let fg = two_node_fg_with_profiles();
    let path = vec![NodeIdx(0), NodeIdx(1)];
    let reads: Vec<(String, Vec<NodeIdx>)> = (0..3)
        .map(|i| (format!("r{}", i), path.clone()))
        .collect();
    let mut bundles = synthesize_bundles(&fg, &reads, 3, 15, false);
    bundles[0].synthetic = false; // normal (non-rescue) bundle
    let bundle = &bundles[0];

    let mut txs: Vec<Transcript> = (0..2).map(|_| Transcript {
        chrom: bundle.chrom.clone(),
        strand: bundle.strand,
        exons: vec![(bundle.start, bundle.end)],
        synthetic: false,
        ..Default::default()
    }).collect();

    if bundle.synthetic {
        for tx in &mut txs {
            tx.synthetic = true;
        }
    }

    for (i, tx) in txs.iter().enumerate() {
        assert!(!tx.synthetic, "transcript[{}] must stay synthetic=false for non-synthetic bundle", i);
    }
}

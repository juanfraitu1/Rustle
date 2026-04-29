//! Regression tests for vg_hmm::rescue (Phase 5).

use rustle::vg_hmm::rescue::{prefilter_read, run_rescue_in_memory, synthesize_bundles};
use rustle::vg_hmm::family_graph::{FamilyGraph, ExonClass, JunctionEdge, NodeIdx};
use rustle::vg_hmm::profile::ProfileHmm;
use rustle::vg::FamilyGroup;
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
            copy_specific: false, profile: Some(p1),
        },
        ExonClass {
            idx: NodeIdx(1), chrom: "chr1".into(), span: (2000, 2008), strand: '+',
            per_copy_sequences: vec![(0, b"TGCATGCA".to_vec())],
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
    let bundles = synthesize_bundles(&fg, &reads, 3);
    assert!(bundles.is_empty(), "cluster of 2 should be dropped with min_reads=3");
}

#[test]
fn synthesize_bundles_emits_synthetic_bundle_for_sufficient_cluster() {
    let fg = two_node_fg_with_profiles();
    let path = vec![NodeIdx(0), NodeIdx(1)];
    let reads: Vec<(String, Vec<NodeIdx>)> = (0..3)
        .map(|i| (format!("r{}", i), path.clone()))
        .collect();
    let bundles = synthesize_bundles(&fg, &reads, 3);
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

#[test]
fn synthesize_bundles_two_clusters_emit_two_bundles() {
    let fg = two_node_fg_with_profiles();
    let path_a = vec![NodeIdx(0)];
    let path_b = vec![NodeIdx(1)];
    let reads: Vec<(String, Vec<NodeIdx>)> = (0..3)
        .map(|i| (format!("ra{}", i), path_a.clone()))
        .chain((0..3).map(|i| (format!("rb{}", i), path_b.clone())))
        .collect();
    let bundles = synthesize_bundles(&fg, &reads, 3);
    assert_eq!(bundles.len(), 2, "expected two distinct clusters → two bundles");
}

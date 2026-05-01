//! Regression tests for HMM-based EM reweighting (`run_em_reweighting_hmm`).

use std::collections::HashMap;
use std::sync::Arc;

use rustle::types::{Bundle, BundleRead, JunctionStats};
use rustle::vg::{run_em_reweighting_hmm, FamilyGroup};
use rustle::vg_hmm::family_graph::{ExonClass, FamilyGraph, JunctionEdge, NodeIdx};
use rustle::vg_hmm::profile::ProfileHmm;

// ── Helpers ────────────────────────────────────────────────────────────────

fn make_bundle_read(name: &str, name_hash: u64, weight: f64) -> BundleRead {
    BundleRead {
        read_uid: name_hash,
        read_name: Arc::from(name),
        read_name_hash: name_hash,
        ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
        ref_start: 0, ref_end: 100,
        exons: Vec::new(),
        junctions: Vec::new(), junction_valid: Vec::new(),
        junctions_raw: Vec::new(), junctions_del: Vec::new(),
        weight,
        is_reverse: false, strand: '+',
        has_poly_start: false, has_poly_end: false,
        has_poly_start_aligned: false, has_poly_start_unaligned: false,
        has_poly_end_aligned: false, has_poly_end_unaligned: false,
        unaligned_poly_t: 0, unaligned_poly_a: 0,
        has_last_exon_polya: false, has_first_exon_polyt: false,
        query_length: None, clip_left: 0, clip_right: 0,
        nh: 2, nm: 0, md: None, insertion_sites: Vec::new(),
        unitig: false, unitig_cov: 0.0, read_count_yc: 0.0,
        countfrag_len: 0.0, countfrag_num: 0.0, junc_mismatch_weight: 0.0,
        pair_idx: Vec::new(), pair_count: Vec::new(),
        mapq: 0, mismatches: Vec::new(),
        hp_tag: None, ps_tag: None, is_primary_alignment: true,
    }
}

fn empty_bundle(chrom: &str) -> Bundle {
    Bundle {
        chrom: chrom.to_string(),
        start: 0, end: 100, strand: '+',
        reads: Vec::new(),
        junction_stats: JunctionStats::default(),
        bundlenodes: None, read_bnodes: None, bnode_colors: None,
        synthetic: false, rescue_class: None,
    }
}

/// Build a 2-paralog FamilyGraph where each paralog has a single
/// **copy-specific** node — paralog 0's node has profile fitted from
/// `seq0`, paralog 1's node from `seq1`. A read that matches `seq0`
/// scores higher against paralog 0's path than paralog 1's, and
/// vice versa.
fn two_paralog_disjoint_fg(seq0: &[u8], seq1: &[u8]) -> FamilyGraph {
    let p0 = ProfileHmm::from_singleton(seq0);
    let p1 = ProfileHmm::from_singleton(seq1);
    let nodes = vec![
        ExonClass {
            idx: NodeIdx(0), chrom: "x".into(), span: (0, seq0.len() as u64), strand: '+',
            per_copy_sequences: vec![(0, seq0.to_vec())],
            per_copy_spans: vec![(0, (0, seq0.len() as u64))],
            copy_specific: true, profile: Some(p0),
        },
        ExonClass {
            idx: NodeIdx(1), chrom: "x".into(), span: (1000, 1000 + seq1.len() as u64), strand: '+',
            per_copy_sequences: vec![(1, seq1.to_vec())],
            per_copy_spans: vec![(1, (1000, 1000 + seq1.len() as u64))],
            copy_specific: true, profile: Some(p1),
        },
    ];
    FamilyGraph { family_id: 0, nodes, edges: Vec::new() }
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[test]
fn em_no_reweighting_when_family_too_small() {
    let fg = two_paralog_disjoint_fg(b"ACGTACGT", b"TGCATGCA");
    let mut bundles = vec![empty_bundle("x")];
    let family = FamilyGroup {
        family_id: 0,
        bundle_indices: vec![0],
        multimap_reads: HashMap::new(),
    };
    let result = run_em_reweighting_hmm(&family, &mut bundles, &fg, &HashMap::new(), 20, 1e-6);
    assert_eq!(result.iterations, 0);
    assert_eq!(result.reads_reweighted, 0);
}

#[test]
fn em_no_reweighting_when_no_multimappers() {
    let fg = two_paralog_disjoint_fg(b"ACGTACGT", b"TGCATGCA");
    let mut bundles = vec![empty_bundle("x"), empty_bundle("x")];
    let family = FamilyGroup {
        family_id: 0,
        bundle_indices: vec![0, 1],
        multimap_reads: HashMap::new(),  // empty — no multi-mappers
    };
    let result = run_em_reweighting_hmm(&family, &mut bundles, &fg, &HashMap::new(), 20, 1e-6);
    assert_eq!(result.iterations, 0);
    assert_eq!(result.reads_reweighted, 0);
}

#[test]
fn em_assigns_paralog_0_read_to_paralog_0() {
    // Paralog 0 sequence = ACGTACGT; paralog 1 sequence = TGCATGCA.
    // Multi-mapped read with sequence "ACGTACGT" should converge to weight ~1
    // for paralog 0 and ~0 for paralog 1.
    let fg = two_paralog_disjoint_fg(b"ACGTACGT", b"TGCATGCA");

    let mut bundles = vec![empty_bundle("x"), empty_bundle("x")];
    let rnh = 0xDEAD_BEEFu64;
    bundles[0].reads.push(make_bundle_read("r1", rnh, 0.5));
    bundles[1].reads.push(make_bundle_read("r1", rnh, 0.5));

    let mut multimap = HashMap::new();
    multimap.insert(rnh, vec![(0, 0), (1, 0)]);  // (fam_pos, read_idx)
    let family = FamilyGroup {
        family_id: 0,
        bundle_indices: vec![0, 1],
        multimap_reads: multimap,
    };

    let mut sequences = HashMap::new();
    sequences.insert(rnh, b"ACGTACGT".to_vec());

    let result = run_em_reweighting_hmm(&family, &mut bundles, &fg, &sequences, 20, 1e-6);
    assert!(result.iterations >= 1, "expected ≥1 EM iteration");
    assert!(result.reads_reweighted >= 1, "expected ≥1 reweight (read mass moved)");

    // After EM, paralog 0's bundle should hold the mass, paralog 1's near 0.
    let w0 = bundles[0].reads[0].weight;
    let w1 = bundles[1].reads[0].weight;
    assert!(w0 > w1, "paralog 0 (matching seq) should outweigh paralog 1: w0={} w1={}", w0, w1);
    assert!(w0 > 0.7, "paralog 0 weight should be dominant (>0.7), got {}", w0);
    assert!(w1 < 0.3, "paralog 1 weight should be small (<0.3), got {}", w1);
    // Weights sum to approximately 1 (softmax normalization).
    assert!((w0 + w1 - 1.0).abs() < 1e-6, "weights should sum to ~1: w0+w1={}", w0+w1);
}

#[test]
fn em_assigns_paralog_1_read_to_paralog_1() {
    // Mirror test: read matches paralog 1.
    let fg = two_paralog_disjoint_fg(b"ACGTACGT", b"TGCATGCA");

    let mut bundles = vec![empty_bundle("x"), empty_bundle("x")];
    let rnh = 0xCAFE_BABEu64;
    bundles[0].reads.push(make_bundle_read("r1", rnh, 0.5));
    bundles[1].reads.push(make_bundle_read("r1", rnh, 0.5));

    let mut multimap = HashMap::new();
    multimap.insert(rnh, vec![(0, 0), (1, 0)]);
    let family = FamilyGroup {
        family_id: 0,
        bundle_indices: vec![0, 1],
        multimap_reads: multimap,
    };

    let mut sequences = HashMap::new();
    sequences.insert(rnh, b"TGCATGCA".to_vec());

    run_em_reweighting_hmm(&family, &mut bundles, &fg, &sequences, 20, 1e-6);

    let w0 = bundles[0].reads[0].weight;
    let w1 = bundles[1].reads[0].weight;
    assert!(w1 > w0, "paralog 1 should win: w0={} w1={}", w0, w1);
    assert!(w1 > 0.7, "paralog 1 weight should be dominant, got {}", w1);
}

#[test]
fn em_skips_read_with_no_sequence_provided() {
    // Without a sequence in the store, EM cannot score the read; weight stays at 0.5.
    let fg = two_paralog_disjoint_fg(b"ACGTACGT", b"TGCATGCA");
    let mut bundles = vec![empty_bundle("x"), empty_bundle("x")];
    let rnh = 0xFEEDu64;
    bundles[0].reads.push(make_bundle_read("r1", rnh, 0.5));
    bundles[1].reads.push(make_bundle_read("r1", rnh, 0.5));

    let mut multimap = HashMap::new();
    multimap.insert(rnh, vec![(0, 0), (1, 0)]);
    let family = FamilyGroup {
        family_id: 0,
        bundle_indices: vec![0, 1],
        multimap_reads: multimap,
    };

    // Empty sequence store.
    let result = run_em_reweighting_hmm(&family, &mut bundles, &fg, &HashMap::new(), 20, 1e-6);
    assert_eq!(result.reads_reweighted, 0);
    // Weights unchanged.
    assert_eq!(bundles[0].reads[0].weight, 0.5);
    assert_eq!(bundles[1].reads[0].weight, 0.5);
}

#[test]
fn em_converges_with_many_reads_split_across_paralogs() {
    // 10 reads, 5 from paralog 0 (ACGTACGT) and 5 from paralog 1 (TGCATGCA).
    // Each is multi-mapped between the two bundles. Final per-bundle coverage
    // should approximately reflect the true split (5:5).
    let fg = two_paralog_disjoint_fg(b"ACGTACGT", b"TGCATGCA");

    let mut bundles = vec![empty_bundle("x"), empty_bundle("x")];
    let mut multimap: HashMap<u64, Vec<(usize, usize)>> = HashMap::new();
    let mut sequences: HashMap<u64, Vec<u8>> = HashMap::new();

    for i in 0..10u64 {
        let rnh = 0xA000 + i;
        let seq: &[u8] = if i < 5 { b"ACGTACGT" } else { b"TGCATGCA" };
        bundles[0].reads.push(make_bundle_read(&format!("r{}", i), rnh, 0.5));
        bundles[1].reads.push(make_bundle_read(&format!("r{}", i), rnh, 0.5));
        let read_idx = bundles[0].reads.len() - 1;
        multimap.insert(rnh, vec![(0, read_idx), (1, read_idx)]);
        sequences.insert(rnh, seq.to_vec());
    }

    let family = FamilyGroup {
        family_id: 0,
        bundle_indices: vec![0, 1],
        multimap_reads: multimap,
    };
    let result = run_em_reweighting_hmm(&family, &mut bundles, &fg, &sequences, 20, 1e-6);
    assert!(result.converged, "EM should converge in <20 iter");

    // Total weight in bundle 0 should ≈ 5 (the 5 ACGT reads), bundle 1 ≈ 5.
    let total0: f64 = bundles[0].reads.iter().map(|r| r.weight).sum();
    let total1: f64 = bundles[1].reads.iter().map(|r| r.weight).sum();
    assert!((total0 - 5.0).abs() < 0.5, "expected total0 ≈ 5, got {}", total0);
    assert!((total1 - 5.0).abs() < 0.5, "expected total1 ≈ 5, got {}", total1);
}

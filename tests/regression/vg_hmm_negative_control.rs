//! Negative-control test: low-complexity reads must not pass the HMM rescue threshold.
//!
//! Adapts the Phase 7.3 plan to the actual Phase 5 `run_rescue_in_memory` signature:
//!
//!   run_rescue_in_memory(
//!       families:      &[FamilyGroup],
//!       family_graphs: &[FamilyGraph],
//!       unmapped_reads: &[(String, Vec<u8>)],
//!       min_loglik:     f64,
//!   ) -> Result<(Vec<NovelCandidate>, Vec<ScoredRead>)>
//!
//! The function already has a unit test for the "empty inputs → empty output" case.
//! Here we cover two additional negative-control scenarios:
//!
//!   1. 50 poly-A reads scored against a two-node family graph with a very low
//!      threshold (-1000.0).  Expected: zero NovelCandidates, because the forward
//!      algorithm assigns essentially 0 probability to AAAA… against an ACGT model.
//!
//!   2. 50 random-sequence reads scored against the same family graph but with
//!      an impossibly high threshold (1e18). Expected: zero NovelCandidates regardless.

use rustle::vg_hmm::rescue::run_rescue_in_memory;
use rustle::vg_hmm::family_graph::{FamilyGraph, ExonClass, JunctionEdge, NodeIdx};
use rustle::vg_hmm::profile::ProfileHmm;
use rustle::vg::FamilyGroup;
use std::collections::HashMap;

fn make_family_graph_for_negative_control() -> FamilyGraph {
    // Two-node graph: exon 0 = "ACGTACGT", exon 1 = "TGCATGCA".
    // The profile HMM is trained on these sequences; poly-A reads should score
    // poorly (low log-likelihood) against this model.
    let p1 = ProfileHmm::from_singleton(b"ACGTACGT");
    let p2 = ProfileHmm::from_singleton(b"TGCATGCA");
    FamilyGraph {
        family_id: 0,
        nodes: vec![
            ExonClass {
                idx: NodeIdx(0),
                chrom: "chrX".into(),
                span: (1000, 1008),
                strand: '+',
                per_copy_sequences: vec![(0, b"ACGTACGT".to_vec())],
                copy_specific: false,
                profile: Some(p1),
            },
            ExonClass {
                idx: NodeIdx(1),
                chrom: "chrX".into(),
                span: (2000, 2008),
                strand: '+',
                per_copy_sequences: vec![(0, b"TGCATGCA".to_vec())],
                copy_specific: false,
                profile: Some(p2),
            },
        ],
        edges: vec![JunctionEdge {
            from: NodeIdx(0),
            to: NodeIdx(1),
            family_support: 1,
            strand: '+',
        }],
    }
}

/// Poly-A reads (length 1000) with a very permissive threshold still yield zero
/// rescues, because the forward algorithm assigns essentially 0 probability to
/// 1000 consecutive A's against an ACGT-model family.
#[test]
fn poly_a_reads_yield_no_rescue() {
    let fg = make_family_graph_for_negative_control();
    let family = FamilyGroup {
        family_id: 0,
        bundle_indices: vec![],
        multimap_reads: HashMap::new(),
    };

    // 50 poly-A reads of length 1000 — maximally low-complexity.
    let reads: Vec<(String, Vec<u8>)> = (0..50)
        .map(|i| (format!("poly_a_{}", i), vec![b'A'; 1000]))
        .collect();

    let (candidates, scored) = run_rescue_in_memory(
        &[family],
        &[fg],
        &reads,
        // Threshold: any score above -1000 is kept.
        // A 1000-base poly-A read should score far below -1000 against an ACGT
        // profile HMM — each position has ~0 emission probability.
        -1000.0,
    )
    .expect("run_rescue_in_memory must not fail");

    assert!(
        candidates.is_empty(),
        "expected 0 NovelCandidates for poly-A reads; got {}",
        candidates.len()
    );
    assert!(
        scored.is_empty(),
        "expected 0 ScoredReads for poly-A reads; got {}",
        scored.len()
    );
}

/// Any reads (even realistic sequences) score 0 rescues when the threshold is
/// set impossibly high (1e18).  This verifies the threshold gate works correctly.
#[test]
fn any_reads_yield_no_rescue_at_high_threshold() {
    let fg = make_family_graph_for_negative_control();
    let family = FamilyGroup {
        family_id: 0,
        bundle_indices: vec![],
        multimap_reads: HashMap::new(),
    };

    // Mix of sequences including one that matches the family exactly.
    let reads: Vec<(String, Vec<u8>)> = vec![
        ("r0".into(), b"ACGTACGTTGCATGCA".to_vec()), // exact family match
        ("r1".into(), b"AAAAAAAAAAAAAAAA".to_vec()), // poly-A
        ("r2".into(), b"GCGCGCGCGCGCGCGC".to_vec()), // low-complexity
    ];

    let (candidates, scored) = run_rescue_in_memory(
        &[family],
        &[fg],
        &reads,
        1e18, // impossibly high threshold
    )
    .expect("run_rescue_in_memory must not fail");

    assert!(
        candidates.is_empty(),
        "expected 0 NovelCandidates with threshold=1e18; got {}",
        candidates.len()
    );
    assert!(
        scored.is_empty(),
        "expected 0 ScoredReads with threshold=1e18; got {}",
        scored.len()
    );
}

/// Empty families list → zero rescues regardless of reads.
#[test]
fn empty_families_yield_no_rescue() {
    let reads: Vec<(String, Vec<u8>)> = (0..50)
        .map(|i| (format!("r{}", i), vec![b'A'; 1000]))
        .collect();

    let (candidates, scored) = run_rescue_in_memory(&[], &[], &reads, -999.0)
        .expect("run_rescue_in_memory must not fail on empty families");

    assert!(candidates.is_empty(), "expected 0 candidates; got {}", candidates.len());
    assert!(scored.is_empty(), "expected 0 scored reads; got {}", scored.len());
}

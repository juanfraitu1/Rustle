use rustle::vg_hmm::profile::ProfileHmm;
use rustle::vg_hmm::scorer::forward_against_profile;

#[test]
fn perfectly_matching_read_scores_higher_than_random() {
    let p = ProfileHmm::from_singleton(b"ACGTACGT");
    let good = forward_against_profile(&p, b"ACGTACGT");
    let bad  = forward_against_profile(&p, b"TTTTTTTT");
    assert!(good > bad, "good={} bad={}", good, bad);
}

use rustle::vg_hmm::family_graph::{ExonClass, FamilyGraph, JunctionEdge, NodeIdx};
use rustle::vg_hmm::scorer::forward_against_family;

fn two_node_chain() -> FamilyGraph {
    let p1 = ProfileHmm::from_singleton(b"ACGT");
    let p2 = ProfileHmm::from_singleton(b"TGCA");
    let nodes = vec![
        ExonClass { idx: NodeIdx(0), chrom: "x".into(), span: (0, 4), strand: '+',
            per_copy_sequences: vec![(0, b"ACGT".to_vec())], copy_specific: true, profile: Some(p1) },
        ExonClass { idx: NodeIdx(1), chrom: "x".into(), span: (10, 14), strand: '+',
            per_copy_sequences: vec![(0, b"TGCA".to_vec())], copy_specific: true, profile: Some(p2) },
    ];
    let edges = vec![JunctionEdge { from: NodeIdx(0), to: NodeIdx(1), family_support: 1, strand: '+' }];
    FamilyGraph { family_id: 0, nodes, edges }
}

#[test]
fn family_forward_path_through_two_nodes() {
    let fg = two_node_chain();
    let on_path = forward_against_family(&fg, b"ACGTTGCA");
    let off_path = forward_against_family(&fg, b"AAAAAAAA");
    assert!(on_path > off_path, "on_path={} off_path={}", on_path, off_path);
}

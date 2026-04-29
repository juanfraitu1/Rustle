use rustle::vg_hmm::family_graph::{ExonClass, FamilyGraph, JunctionEdge, NodeIdx};
use rustle::types::{Bundle, BundleRead, JunctionStats};
use std::sync::Arc;

fn mk_bundle_with_reads(start: u64, end: u64, exons: Vec<(u64, u64)>) -> Bundle {
    let read = BundleRead {
        read_uid: 1, read_name: Arc::from("r1"), read_name_hash: 0,
        ref_id: None, mate_ref_id: None, mate_start: None, hi: 0,
        ref_start: start, ref_end: end, exons: exons.clone(),
        junctions: Vec::new(), junction_valid: Vec::new(),
        junctions_raw: Vec::new(), junctions_del: Vec::new(),
        weight: 1.0, is_reverse: false, strand: '+',
        has_poly_start: false, has_poly_end: false,
        has_poly_start_aligned: false, has_poly_start_unaligned: false,
        has_poly_end_aligned: false, has_poly_end_unaligned: false,
        unaligned_poly_t: 0, unaligned_poly_a: 0,
        has_last_exon_polya: false, has_first_exon_polyt: false,
        query_length: None, clip_left: 0, clip_right: 0, nh: 1, nm: 0, md: None,
        insertion_sites: Vec::new(),
        unitig: false, unitig_cov: 0.0, read_count_yc: 0.0,
        countfrag_len: 0.0, countfrag_num: 0.0, junc_mismatch_weight: 0.0,
        pair_idx: Vec::new(), pair_count: Vec::new(),
        mapq: 0, mismatches: Vec::new(), hp_tag: None, ps_tag: None,
        is_primary_alignment: true,
    };
    Bundle {
        chrom: "chr1".into(), start, end, strand: '+',
        reads: vec![read],
        junction_stats: JunctionStats::default(),
        bundlenodes: None, read_bnodes: None, bnode_colors: None,
    }
}

#[test]
fn empty_family_graph_constructs() {
    let fg = FamilyGraph {
        family_id: 0,
        nodes: Vec::<ExonClass>::new(),
        edges: Vec::<JunctionEdge>::new(),
    };
    assert_eq!(fg.family_id, 0);
    assert!(fg.nodes.is_empty());
    assert!(fg.edges.is_empty());
}

#[test]
fn exon_class_carries_per_copy_sequences() {
    let cls = ExonClass {
        idx: NodeIdx(0),
        chrom: "chrX".into(),
        span: (1000, 1500),
        strand: '+',
        per_copy_sequences: vec![(0, b"ACGT".to_vec()), (1, b"ACAT".to_vec())],
        copy_specific: false,
    };
    assert_eq!(cls.per_copy_sequences.len(), 2);
    assert!(!cls.copy_specific);
}

#[test]
fn exon_extraction_returns_unique_sorted_intervals() {
    use rustle::vg_hmm::family_graph::extract_copy_exons;
    let b = mk_bundle_with_reads(100, 500, vec![(100, 200), (300, 400), (450, 500)]);
    let exons = extract_copy_exons(&b);
    assert_eq!(exons, vec![(100, 200), (300, 400), (450, 500)]);
}

#[test]
fn build_family_graph_two_copy_smoke() {
    use rustle::vg_hmm::family_graph::build_family_graph;
    use rustle::vg::FamilyGroup;
    use std::collections::HashMap;
    // Two bundles, identical exon coordinates.
    let b0 = mk_bundle_with_reads(100, 500, vec![(100, 200), (300, 400)]);
    let b1 = mk_bundle_with_reads(100, 500, vec![(100, 200), (300, 400)]);
    let bundles = vec![b0, b1];
    let family = FamilyGroup { family_id: 7, bundle_indices: vec![0, 1], multimap_reads: HashMap::new() };

    // No genome FASTA available in this test — pass `None`; build should still
    // produce nodes (with empty per-copy sequences) and edges.
    let fg = build_family_graph(&family, &bundles, None, 0.30, 0.30).unwrap();
    assert_eq!(fg.family_id, 7);
    assert_eq!(fg.nodes.len(), 2, "two shared exon classes expected");
    assert!(fg.nodes.iter().all(|n| !n.copy_specific));
}

#[test]
fn junction_edges_count_unique_copies_supporting_each_edge() {
    use rustle::vg_hmm::family_graph::collect_family_junctions;
    // Three copies, two of which share a junction at (1000,1100); one has a private (2000,2100).
    let per_copy = vec![
        ('+', vec![(1000u64, 1100u64), (3000, 3100)]),
        ('+', vec![(1000, 1100)]),
        ('+', vec![(2000, 2100)]),
    ];
    let edges = collect_family_junctions(&per_copy);
    // Expect 3 distinct junctions with supports 2, 1, 1.
    let supp = |donor, accept| edges.iter()
        .find(|e| e.donor == donor && e.acceptor == accept)
        .map(|e| e.family_support).unwrap_or(0);
    assert_eq!(supp(1000, 1100), 2);
    assert_eq!(supp(2000, 2100), 1);
    assert_eq!(supp(3000, 3100), 1);
}

#[test]
fn minimizer_jaccard_splits_position_cluster_when_sequences_diverge() {
    use rustle::vg_hmm::family_graph::refine_by_minimizer_jaccard;
    // Two "exons" at the same position cluster but with no shared k-mers.
    let cluster = vec![(0usize, b"AAAAAAAAAAAAAAAA".to_vec()),
                       (1usize, b"GGGGGGGGGGGGGGGG".to_vec())];
    let split = refine_by_minimizer_jaccard(&cluster, 0.30, 15, 10);
    assert_eq!(split.len(), 2, "fully divergent should split");
}

#[test]
fn minimizer_jaccard_keeps_similar_sequences_together() {
    use rustle::vg_hmm::family_graph::refine_by_minimizer_jaccard;
    let s = b"ACGTACGTACGTACGTACGTACGTACGT".to_vec();
    let cluster = vec![(0usize, s.clone()), (1usize, s.clone())];
    let split = refine_by_minimizer_jaccard(&cluster, 0.30, 15, 10);
    assert_eq!(split.len(), 1, "identical should stay together");
}

#[test]
fn position_overlap_clusters_partition_exons() {
    use rustle::vg_hmm::family_graph::cluster_by_position;
    // copy 0: exons at 100-200 and 300-400
    // copy 1: exons at 110-210 (overlaps copy0[0]) and 500-600 (no overlap)
    let copy0 = vec![(100u64, 200u64), (300, 400)];
    let copy1 = vec![(110u64, 210u64), (500, 600)];
    let clusters = cluster_by_position(&[("chrA", '+', copy0), ("chrA", '+', copy1)], 0.30);
    // Expect 3 clusters: {(c0,0),(c1,0)}, {(c0,1)}, {(c1,1)}
    assert_eq!(clusters.len(), 3);
    assert!(clusters.iter().any(|c| c.len() == 2)); // the overlapping pair
}

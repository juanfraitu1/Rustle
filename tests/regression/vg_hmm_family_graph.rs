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

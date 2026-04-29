use rustle::vg_hmm::family_graph::{ExonClass, FamilyGraph, JunctionEdge, NodeIdx};

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

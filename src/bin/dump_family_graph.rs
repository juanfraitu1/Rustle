/// Smoke checkpoint: Task 1.6.5 of the vg-hmm novel-copy plan.
///
/// Constructs synthetic Bundles from hard-coded GOLGA6L7 paralog exon
/// coordinates (extracted from GGO_genomic.gff), calls build_family_graph,
/// and prints node/edge counts + structure to stdout.
///
/// Expected: 27 nodes (9 per copy × 3, all copy_specific=true) because the
/// three GOLGA6L7 copies are tandem paralogs at non-overlapping genomic
/// positions on NC_073243.2.
///
/// Delete after Phase 7.

use std::collections::HashMap;
use std::sync::Arc;

use rustle::types::{Bundle, BundleRead, Junction, JunctionStat, JunctionStats};
use rustle::vg::FamilyGroup;
use rustle::vg_hmm::family_graph::build_family_graph;

// ── Hard-coded GOLGA6L7 exon coordinates (0-based, half-open) ──────────────
// Extracted from GGO_genomic.gff: $4-1 gives 0-based start, $5 is end.
//
// LOC115931294 / XM_055367547.2 (copy 0)
const COPY0_EXONS: [(u64, u64); 9] = [
    (104789646, 104789918),
    (104791223, 104791352),
    (104792166, 104792196),
    (104792465, 104792516),
    (104792604, 104792698),
    (104792780, 104792887),
    (104794129, 104794188),
    (104794517, 104794659),
    (104794953, 104796276),
];

// LOC134757625 / XM_063700849.1 (copy 1)
const COPY1_EXONS: [(u64, u64); 9] = [
    (104830535, 104830737),
    (104832042, 104832171),
    (104832985, 104833015),
    (104833284, 104833335),
    (104833423, 104833517),
    (104833599, 104833706),
    (104834952, 104835011),
    (104835340, 104835482),
    (104835776, 104837094),
];

// LOC101137218 / XM_055367548.2 (copy 2)
const COPY2_EXONS: [(u64, u64); 9] = [
    (104871355, 104871545),
    (104872850, 104872979),
    (104873793, 104873823),
    (104874092, 104874143),
    (104874231, 104874325),
    (104874407, 104874514),
    (104875758, 104875817),
    (104876146, 104876288),
    (104876582, 104877901),
];

fn build_junction_stats(exons: &[(u64, u64)]) -> JunctionStats {
    let mut stats: JunctionStats = HashMap::with_hasher(fxhash::FxBuildHasher::default());
    // 8 junctions for 9 exons: donor = exon[i].end, acceptor = exon[i+1].start
    for i in 0..exons.len() - 1 {
        let junc = Junction { donor: exons[i].1, acceptor: exons[i + 1].0 };
        stats.insert(junc, JunctionStat::default());
    }
    stats
}

fn make_bundle(exons: &[(u64, u64)]) -> Bundle {
    let start = exons.iter().map(|e| e.0).min().unwrap();
    let end = exons.iter().map(|e| e.1).max().unwrap();

    let read = BundleRead {
        read_uid: 1,
        read_name: Arc::from("synthetic_r1"),
        read_name_hash: 0,
        ref_id: None,
        mate_ref_id: None,
        mate_start: None,
        hi: 0,
        ref_start: start,
        ref_end: end,
        exons: exons.to_vec(),
        junctions: Vec::new(),
        junction_valid: Vec::new(),
        junctions_raw: Vec::new(),
        junctions_del: Vec::new(),
        weight: 1.0,
        is_reverse: false,
        strand: '+',
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
    };

    Bundle {
        chrom: "NC_073243.2".to_string(),
        start,
        end,
        strand: '+',
        reads: vec![read],
        junction_stats: build_junction_stats(exons),
        bundlenodes: None,
        read_bnodes: None,
        bnode_colors: None,
    }
}

fn main() {
    let bundles = vec![
        make_bundle(&COPY0_EXONS),
        make_bundle(&COPY1_EXONS),
        make_bundle(&COPY2_EXONS),
    ];

    let family = FamilyGroup {
        family_id: 0,
        bundle_indices: vec![0, 1, 2],
        multimap_reads: HashMap::new(),
    };

    let fg = build_family_graph(&family, &bundles, None, 0.30, 0.30)
        .expect("build_family_graph failed");

    println!(
        "[smoke] family graph: {} nodes, {} edges",
        fg.n_nodes(),
        fg.n_edges()
    );
    for n in &fg.nodes {
        println!(
            "  node {:>2}: {} {}-{} ({} copies, copy_specific={})",
            n.idx.0,
            n.chrom,
            n.span.0,
            n.span.1,
            n.per_copy_sequences.len(),
            n.copy_specific
        );
    }
    for e in &fg.edges {
        println!(
            "  edge: {:>2} -> {:>2} (support={}, strand={})",
            e.from.0, e.to.0, e.family_support, e.strand
        );
    }
}

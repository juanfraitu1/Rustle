//! Layer-2 orchestrator. Implements C2–C3 (family discovery + variation-graph
//! merge); C4–C6 (secondary amendment, constrained flow-decomposition, emission)
//! are deferred to later milestones. Reads Layer-1 splice graphs and the C1
//! side-index; never re-bundles, never re-splits, never shifts a coordinate.
//! Default-off behind `config.vg_layer2`.

use crate::genome::GenomeIndex;
use crate::graph::Graph;
use crate::path_extract::Transcript;
use crate::types::{DetHashMap, DetHashSet};
use crate::vg::FamilyGroup;
use crate::vg_family::family_graph::FamilyGraph;
use crate::vg_family::secondary_index::SecondaryIndex;
use anyhow::Result;

/// A candidate junction-edge folded into a family graph from a secondary.
#[derive(Debug, Clone, PartialEq)]
pub struct AmendCandidate {
    pub intron: (u64, u64),
    pub support: f64,
}

/// Result of amending a family graph with secondaries: candidate edges PLUS the
/// (node_idx, copy_id) memberships the secondaries imply.
#[derive(Debug, Default)]
pub struct GraphAmendment {
    pub edges: Vec<AmendCandidate>,
    pub copy_membership: Vec<(usize, usize)>,
}

/// (C4 amend) Fold side-index secondaries into candidate edges + copy memberships.
/// A junction is accepted only if its donor/acceptor already exist as Layer-1 node
/// boundaries in the merged graph (no invented coordinates — additivity). Each
/// accepted junction registers the secondary's copy on the donor+acceptor nodes.
/// `copy_of_locus[side_index_locus] = family copy_id`.
pub fn amend_family_graph(
    fg: &FamilyGraph,
    secondaries: &[crate::vg_family::secondary_index::SecondaryAlignment],
    copy_of_locus: &DetHashMap<usize, usize>,
) -> GraphAmendment {
    let mut end_to_node: DetHashMap<u64, usize> = DetHashMap::default();
    let mut start_to_node: DetHashMap<u64, usize> = DetHashMap::default();
    for (i, n) in fg.nodes.iter().enumerate() {
        end_to_node.insert(n.span.1, i);
        start_to_node.insert(n.span.0, i);
    }
    let mut support: DetHashMap<(u64, u64), f64> = DetHashMap::default();
    let mut membership: DetHashSet<(usize, usize)> = DetHashSet::default();
    for s in secondaries {
        let copy = s.locus.and_then(|l| copy_of_locus.get(&l).copied());
        for &(d, a) in &s.introns {
            if let (Some(&u), Some(&v)) = (end_to_node.get(&d), start_to_node.get(&a)) {
                *support.entry((d, a)).or_default() += 1.0;
                if let Some(c) = copy {
                    membership.insert((u, c));
                    membership.insert((v, c));
                }
            }
        }
    }
    let mut edges: Vec<AmendCandidate> = support
        .into_iter()
        .map(|(intron, support)| AmendCandidate { intron, support })
        .collect();
    edges.sort_by(|x, y| x.intron.cmp(&y.intron));
    let mut copy_membership: Vec<(usize, usize)> = membership.into_iter().collect();
    copy_membership.sort_unstable();
    GraphAmendment { edges, copy_membership }
}

/// A recovered copy/isoform path (exon chain in genomic coordinates). A path that
/// violates allele-linkage (a copy claiming a diagnostic node it has no allele for)
/// is DROPPED during decomposition — never emitted.
#[derive(Debug, Clone, PartialEq)]
pub struct FamilyPath {
    pub exons: Vec<(u64, u64)>,
    pub flow: f64,
    pub copy_id: usize,
}

/// (C4 decompose) Recover copies (paths) and isoforms (sub-paths) by CONSTRAINED
/// flow-decomposition of the amended family graph. `enumerate_diagnostic_sites`
/// yields PSV columns; any candidate path requiring alleles from two different
/// copies at a diagnostic node is rejected as chimeric (thesis "constrained
/// flow-decomposition under allele-linkage").
pub fn decompose_family_paths(
    fg: &FamilyGraph,
    amendment: &GraphAmendment,
    min_flow: f64,
) -> Vec<FamilyPath> {
    let n_copies = fg
        .nodes
        .iter()
        .flat_map(|n| n.per_copy_sequences.iter().map(|(c, _)| *c))
        .max()
        .map(|m| m + 1)
        .unwrap_or(0);
    if n_copies == 0 {
        return Vec::new();
    }
    // PSV-driven allele-linkage constraint. `enumerate_diagnostic_sites` tells us
    // whether the family is identifiable at all (n_sites > 0); we enforce linkage
    // only then. A node is DIAGNOSTIC when its present copies carry >= 2 distinct
    // alleles (a PSV column lives there). The constraint: a copy may borrow
    // strength on shared / non-diagnostic backbone nodes, but may NOT claim a
    // diagnostic node where it has no native allele — that would assign it another
    // copy's distinguishing sequence (a chimeric, allele-mixing path).
    let fingerprints = crate::vg::enumerate_diagnostic_sites(fg, n_copies);
    let enforce_linkage = fingerprints.n_sites > 0;
    let mut diagnostic_node: DetHashSet<usize> = DetHashSet::default();
    for (ni, n) in fg.nodes.iter().enumerate() {
        let mut distinct: DetHashSet<&[u8]> = DetHashSet::default();
        for (_, s) in &n.per_copy_sequences {
            distinct.insert(s.as_slice());
        }
        if distinct.len() >= 2 {
            diagnostic_node.insert(ni);
        }
    }
    let mut adj: DetHashMap<usize, Vec<usize>> = DetHashMap::default();
    let mut edge_flow: DetHashMap<(usize, usize), f64> = DetHashMap::default();
    let mut end_to_node: DetHashMap<u64, usize> = DetHashMap::default();
    let mut start_to_node: DetHashMap<u64, usize> = DetHashMap::default();
    for (i, n) in fg.nodes.iter().enumerate() {
        end_to_node.insert(n.span.1, i);
        start_to_node.insert(n.span.0, i);
    }
    for e in &fg.edges {
        let (u, v) = (e.from.0, e.to.0);
        adj.entry(u).or_default().push(v);
        *edge_flow.entry((u, v)).or_default() += e.family_support as f64;
    }
    for c in &amendment.edges {
        if let (Some(&u), Some(&v)) =
            (end_to_node.get(&c.intron.0), start_to_node.get(&c.intron.1))
        {
            adj.entry(u).or_default().push(v);
            *edge_flow.entry((u, v)).or_default() += c.support;
        }
    }
    for v in adj.values_mut() {
        v.sort_unstable();
        v.dedup();
    }
    let mut node_copies: DetHashMap<usize, DetHashSet<usize>> = DetHashMap::default();
    for (ni, n) in fg.nodes.iter().enumerate() {
        let set = node_copies.entry(ni).or_default();
        for (c, _) in &n.per_copy_sequences {
            set.insert(*c);
        }
    }
    for (ni, c) in &amendment.copy_membership {
        node_copies.entry(*ni).or_default().insert(*c);
    }
    let mut paths: Vec<FamilyPath> = Vec::new();
    for copy in 0..n_copies {
        let mut nodes_of_copy: Vec<usize> = (0..fg.nodes.len())
            .filter(|ni| node_copies.get(ni).map(|s| s.contains(&copy)).unwrap_or(false))
            .collect();
        nodes_of_copy.sort_by_key(|&i| fg.nodes[i].span.0);
        if nodes_of_copy.len() < 2 {
            continue;
        }
        let connected = nodes_of_copy.windows(2).all(|w| {
            edge_flow.get(&(w[0], w[1])).map(|f| *f >= min_flow).unwrap_or(false)
        });
        if !connected {
            continue;
        }
        // Allele-linkage: reject this copy's path if it traverses a diagnostic
        // node where the copy has no native allele (it would borrow another copy's
        // distinguishing sequence = chimeric). Shared / non-diagnostic nodes are
        // freely borrowable (the starved-copy strength-borrowing case).
        let linkage_ok = !enforce_linkage
            || nodes_of_copy.iter().all(|&ni| {
                !(diagnostic_node.contains(&ni)
                    && !fg.nodes[ni]
                        .per_copy_sequences
                        .iter()
                        .any(|(cc, _)| *cc == copy))
            });
        if !linkage_ok {
            continue;
        }
        let exons: Vec<(u64, u64)> =
            nodes_of_copy.iter().map(|&i| fg.nodes[i].span).collect();
        let flow = nodes_of_copy
            .windows(2)
            .map(|w| edge_flow[&(w[0], w[1])])
            .fold(f64::INFINITY, f64::min);
        paths.push(FamilyPath { exons, flow, copy_id: copy });
    }
    paths.sort_by(|a, b| a.exons.cmp(&b.exons).then(a.copy_id.cmp(&b.copy_id)));
    paths
}

/// Clustering thresholds passed through to the family-graph merge. Layer-2
/// families are already screened by `min_similarity` at discovery, so positional
/// overlap is left permissive (0.0 — accept any) rather than the bundle path's
/// stricter `min_pos_recip`; minimizer-refinement is likewise left to the
/// `merge_min_jaccard` (= `min_similarity`) stage. (NOT the bundle-path values.)
const FG_MIN_POS_RECIP: f64 = 0.0;
const FG_REFINE_MIN_JACCARD: f64 = 0.0;

/// Result of a Layer-2 pass for one chromosome.
#[derive(Debug, Default)]
pub struct Layer2Output {
    pub families: Vec<FamilyGroup>,
    pub family_graphs: Vec<Option<FamilyGraph>>,
    pub novel_transcripts: Vec<Transcript>,
}

/// One per-locus input to Layer 2: the Layer-1 splice graph + its chrom/strand.
pub struct Layer1Locus<'a> {
    pub chrom: String,
    pub strand: char,
    pub graph: &'a Graph,
}

/// Run Layer 2 over one chromosome's Layer-1 loci.
///
/// `primary_locus[read_name_hash] = bundle index of the read's primary`. The
/// side-index is consumed (pruned to family-candidate loci, capped per locus).
///
/// CALLER CONTRACT (M4.3 wiring): every value in `primary_locus` and every
/// `SecondaryAlignment::locus` in `side_index` MUST be a valid index into `loci`
/// (i.e. < `loci.len()`) — they are used to index `loci` directly. The pipeline
/// builds both from the same bundle list `loci` is derived from, so this holds by
/// construction; a `debug_assert!` below catches violations in test/debug builds.
///
/// With `genome = None`, similarity is never computed (`similarity` stays empty);
/// families then form on link count alone iff `min_similarity <= 0.0`.
#[allow(clippy::too_many_arguments)]
pub fn run_layer2(
    loci: &[Layer1Locus<'_>],
    mut side_index: SecondaryIndex,
    primary_locus: &DetHashMap<u64, usize>,
    genome: Option<&GenomeIndex>,
    min_link: u32,
    min_similarity: f64,
    per_locus_cap: usize,
    kmer_len: usize,
    new_copies: bool,
) -> Result<Layer2Output> {
    // (C2) candidate links → compute similarity only for linked, same-chrom pairs.
    let links = crate::vg::build_multimap_index_from_secondary_index(&side_index, primary_locus);
    debug_assert!(
        links.iter().all(|((a, b), _)| *a < loci.len() && *b < loci.len()),
        "run_layer2: link locus index out of range of loci slice (caller contract)"
    );
    let mut similarity: DetHashMap<(usize, usize), f64> = DetHashMap::default();
    if let Some(g) = genome {
        for ((a, b), count) in &links {
            if *count < min_link {
                continue;
            }
            let (la, lb) = (&loci[*a], &loci[*b]);
            if la.chrom == lb.chrom {
                let sim = crate::vg::exon_kmer_similarity_between_graphs(
                    la.graph, lb.graph, &la.chrom, g, kmer_len,
                );
                similarity.insert((*a, *b), sim);
            }
        }
    }

    let families = crate::vg::discover_family_groups_layer2(
        &side_index, primary_locus, &similarity, min_link, min_similarity,
    );

    // (C1 prune) keep only family-candidate loci, cap per locus (logged drops).
    let keep: DetHashSet<usize> = families
        .iter()
        .flat_map(|f| f.bundle_indices.iter().copied())
        .collect();
    let pruned = side_index.prune_to_loci(&keep);
    let capped = side_index.cap_per_locus(per_locus_cap);
    if pruned + capped > 0 {
        eprintln!(
            "[layer2] side-index pruned {pruned} (non-candidate) + capped {capped} (per-locus>{per_locus_cap})"
        );
    }

    // (C3) merge each family's Layer-1 graphs into a variation graph.
    let mut family_graphs: Vec<Option<FamilyGraph>> = Vec::with_capacity(families.len());
    for fam in &families {
        let copies: Vec<(String, char, &Graph)> = fam
            .bundle_indices
            .iter()
            .map(|&i| (loci[i].chrom.clone(), loci[i].strand, loci[i].graph))
            .collect();
        let fg = crate::vg_family::family_graph::build_family_graph_from_layer1_graphs(
            fam.family_id,
            &copies,
            genome,
            FG_MIN_POS_RECIP,
            min_similarity,
            FG_REFINE_MIN_JACCARD,
        )
        .ok();
        family_graphs.push(fg);
    }

    // (M4.3c) Materialize recovered copy/isoform paths into emittable Transcripts.
    // For each family: amend its merged graph with the family's pruned secondaries,
    // then constrained-flow-decompose into per-copy paths. Each surviving path becomes
    // a UnionBaseline-tagged Transcript so the pipeline's union-by-chain stage emits it
    // only if its intron chain is absent from the VG output (strictly additive).
    let mut novel_transcripts: Vec<Transcript> = Vec::new();
    for (fi, fam) in families.iter().enumerate() {
        let Some(fg) = family_graphs[fi].as_ref() else { continue };
        let mut copy_of_locus: DetHashMap<usize, usize> = DetHashMap::default();
        for (copy, &loc) in fam.bundle_indices.iter().enumerate() {
            copy_of_locus.insert(loc, copy);
        }
        let member_set: DetHashSet<usize> = fam.bundle_indices.iter().copied().collect();
        let secs: Vec<crate::vg_family::secondary_index::SecondaryAlignment> = side_index
            .alignments().iter()
            .filter(|a| a.locus.map(|l| member_set.contains(&l)).unwrap_or(false))
            .cloned().collect();
        let am = amend_family_graph(fg, &secs, &copy_of_locus);
        let paths = decompose_family_paths(fg, &am, 1.0);
        for p in paths {
            let chrom = loci[fam.bundle_indices[0]].chrom.clone();
            let strand = loci[fam.bundle_indices[0]].strand;
            let mut t = Transcript::default();
            t.chrom = chrom;
            t.strand = strand;
            t.exons = p.exons;
            t.coverage = p.flow;
            t.longcov = p.flow;
            t.is_longread = true;
            t.synthetic = true;
            t.vg_family_id = Some(fam.family_id);
            t.vg_copy_id = Some(p.copy_id);
            t.rescue_class = Some(crate::vg_family::diagnostic::RescueClass::UnionBaseline);
            novel_transcripts.push(t);
        }
    }
    let _ = new_copies; // C5 (M5)
    Ok(Layer2Output {
        families,
        family_graphs,
        novel_transcripts,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sec_aln(
        h: u64, start: u64, end: u64, introns: &[(u64, u64)], locus: usize,
    ) -> crate::vg_family::secondary_index::SecondaryAlignment {
        crate::vg_family::secondary_index::SecondaryAlignment {
            read_name_hash: h, chrom: "chrT".into(), ref_start: start, ref_end: end,
            introns: introns.to_vec(), nm: 0, is_supplementary: false, locus: Some(locus),
        }
    }

    #[test]
    fn amend_adds_candidate_edges_from_secondaries() {
        let g0 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160), (300, 360), (500, 560)],
        );
        let g1 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160)],
        );
        let genome = crate::vg_family::family_graph::tests_support::make_two_copy_genome_3exon();
        let fg = crate::vg_family::family_graph::build_family_graph_from_layer1_graphs(
            0, &[("chrT".to_string(), '+', &g0), ("chrT".to_string(), '+', &g1)],
            Some(&genome), 0.0, 0.5, 0.0,
        ).unwrap();
        let secs = vec![
            sec_aln(7, 100, 360, &[(160, 300)], 1),
            sec_aln(8, 300, 560, &[(360, 500)], 1),
        ];
        let mut copy_of_locus: DetHashMap<usize, usize> = DetHashMap::default();
        copy_of_locus.insert(0, 0);
        copy_of_locus.insert(1, 1);
        let am = amend_family_graph(&fg, &secs, &copy_of_locus);
        assert_eq!(am.edges.len(), 2, "two candidate junction-edges added");
        assert!(am.edges.iter().any(|c| c.intron == (160, 300)));
        assert!(am.copy_membership.iter().any(|(_, c)| *c == 1),
            "starved copy registered as contributor to traversed nodes");
    }

    #[test]
    fn decompose_recovers_starved_copy_path() {
        let g0 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160), (300, 360), (500, 560)],
        );
        let g1 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160)],
        );
        let genome = crate::vg_family::family_graph::tests_support::make_two_copy_genome_3exon();
        let fg = crate::vg_family::family_graph::build_family_graph_from_layer1_graphs(
            0, &[("chrT".to_string(), '+', &g0), ("chrT".to_string(), '+', &g1)],
            Some(&genome), 0.0, 0.5, 0.0,
        ).unwrap();
        let secs = vec![
            sec_aln(7, 100, 360, &[(160, 300)], 1),
            sec_aln(8, 300, 560, &[(360, 500)], 1),
        ];
        let mut copy_of_locus: DetHashMap<usize, usize> = DetHashMap::default();
        copy_of_locus.insert(0, 0);
        copy_of_locus.insert(1, 1);
        let am = amend_family_graph(&fg, &secs, &copy_of_locus);
        let paths = decompose_family_paths(&fg, &am, 1.0);
        assert!(
            paths.iter().any(|p| p.exons == vec![(100, 160), (300, 360), (500, 560)]),
            "starved copy's full chain recovered: {paths:?}"
        );
    }

    #[test]
    fn decompose_forbids_allele_mixing_path() {
        // 3-copy family. Middle node M (idx 1) is DIAGNOSTIC: copy0 and copy1 carry
        // distinct alleles (AAAAAA vs CCCCCC); copy2 has NO native allele at M.
        // copy2 is registered on M only via an AMENDMENT (a secondary). copy2's
        // walk A→M→C is edge-connected, but claiming the diagnostic node M with no
        // allele = chimeric → it MUST be rejected, while copy0/copy1 (native M
        // alleles) survive. This genuinely exercises the linkage guard: deleting
        // the guard would let copy2's chimeric path leak (test would fail).
        let fg = build_psv_three_copy_fg();
        let am = GraphAmendment {
            edges: Vec::new(),
            copy_membership: vec![(1, 2)], // amendment puts copy 2 onto diagnostic node M
        };
        let paths = decompose_family_paths(&fg, &am, 0.0);
        assert!(
            paths.iter().any(|p| p.copy_id == 0 && p.exons.len() == 3),
            "copy 0 (native M allele) recovered: {paths:?}"
        );
        assert!(
            paths.iter().any(|p| p.copy_id == 1 && p.exons.len() == 3),
            "copy 1 (native M allele) recovered: {paths:?}"
        );
        assert!(
            !paths.iter().any(|p| p.copy_id == 2),
            "copy 2 has no allele at diagnostic M → its chimeric path is rejected: {paths:?}"
        );
    }

    /// 3 copies; node M (idx 1) shared-coordinate but PSV-distinguishable between
    /// copy0/copy1 (distinct alleles), with copy2 ABSENT at M. Nodes A and C are
    /// shared identical across all three copies.
    fn build_psv_three_copy_fg() -> FamilyGraph {
        use crate::vg_family::family_graph::{ExonClass, JunctionEdge, NodeIdx};
        let mk_node = |idx: usize, span: (u64, u64), seqs: &[(usize, &[u8])]| ExonClass {
            idx: NodeIdx(idx),
            chrom: "chrT".into(),
            span,
            strand: '+',
            per_copy_sequences: seqs.iter().map(|(c, s)| (*c, s.to_vec())).collect(),
            per_copy_spans: seqs.iter().map(|(c, _)| (*c, span)).collect(),
            copy_specific: seqs.len() == 1,
            per_copy_cov: Vec::new(),
        };
        let nodes = vec![
            mk_node(0, (100, 160), &[(0, b"ACGTAC"), (1, b"ACGTAC"), (2, b"ACGTAC")]),
            // diagnostic: copy0 vs copy1 differ; copy2 absent (recovered only via amendment).
            mk_node(1, (300, 360), &[(0, b"AAAAAA"), (1, b"CCCCCC")]),
            mk_node(2, (500, 560), &[(0, b"TTGGCC"), (1, b"TTGGCC"), (2, b"TTGGCC")]),
        ];
        let edges = vec![
            JunctionEdge { from: NodeIdx(0), to: NodeIdx(1), family_support: 5, strand: '+' },
            JunctionEdge { from: NodeIdx(1), to: NodeIdx(2), family_support: 5, strand: '+' },
        ];
        FamilyGraph { family_id: 0, nodes, edges }
    }

    #[test]
    fn run_layer2_forms_family_and_merges_graph() {
        let g0 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160), (300, 360)],
        );
        let g1 = crate::vg_family::family_graph::tests_support::make_layer1_graph(
            "chrT", '+', &[(100, 160), (5300, 5360)],
        );
        let genome = crate::vg_family::family_graph::tests_support::make_two_copy_genome();

        let loci = vec![
            Layer1Locus { chrom: "chrT".into(), strand: '+', graph: &g0 },
            Layer1Locus { chrom: "chrT".into(), strand: '+', graph: &g1 },
        ];

        let mut si = SecondaryIndex::new();
        si.push(crate::vg_family::secondary_index::SecondaryAlignment {
            read_name_hash: 7,
            chrom: "chrT".into(),
            ref_start: 100,
            ref_end: 160,
            introns: vec![],
            nm: 0,
            is_supplementary: false,
            locus: Some(1),
        });
        let mut primary_locus: DetHashMap<u64, usize> = DetHashMap::default();
        primary_locus.insert(7, 0);

        let out = run_layer2(
            &loci, si, &primary_locus, Some(&genome),
            1, 0.3, 1000, 15, false,
        )
        .unwrap();

        assert_eq!(out.families.len(), 1, "one family from C2");
        assert!(out.family_graphs[0].is_some(), "family graph merged (C3)");
        // M4.3c: emission is now wired. This fixture carries no amending secondaries
        // (the lone secondary has no introns), so any emitted paths come purely from the
        // merged family graph's own edges; assert the emitted transcripts (if any) are
        // well-formed and carry the UnionBaseline tag that drives union-by-chain.
        for t in &out.novel_transcripts {
            assert_eq!(
                t.rescue_class,
                Some(crate::vg_family::diagnostic::RescueClass::UnionBaseline),
                "emitted Layer-2 transcripts are tagged UnionBaseline"
            );
            assert_eq!(t.vg_family_id, Some(0), "tagged with the family id");
            assert!(t.exons.len() >= 2, "emitted paths are multi-exon");
        }
    }
}

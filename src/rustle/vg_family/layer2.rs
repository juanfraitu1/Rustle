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
        // Emit this copy's OWN coordinates per node — NEVER the node's `span` (the
        // UNION of all contributing copies' spans, which would fabricate coordinates).
        //   - native node (this copy has a `per_copy_spans` entry) → that coordinate.
        //   - amendment-BORROWED node (no native entry): placeable ONLY when the node
        //     carries >=2 per_copy_spans entries that AGREE — i.e. multiple copies
        //     independently confirm the SAME coordinate (genuine co-location). A
        //     single foreign entry is a DIFFERENT-LOCUS paralog exon (the common
        //     cross-mapped-secondary case): borrowing it would assign this copy the
        //     OTHER copy's genomic coordinates — a wrong-coordinate false positive.
        //     Placing a cross-locus borrowed copy correctly requires homology
        //     coordinate-transfer (the thesis "isoform transfer via homology" — not
        //     yet implemented) → DROP the path rather than emit fabricated coords.
        let exons: Option<Vec<(u64, u64)>> = nodes_of_copy
            .iter()
            .map(|&i| {
                let n = &fg.nodes[i];
                if let Some((_, sp)) = n.per_copy_spans.iter().find(|(c, _)| *c == copy) {
                    return Some(*sp);
                }
                if n.per_copy_spans.len() >= 2 {
                    let f = n.per_copy_spans[0].1;
                    if n.per_copy_spans.iter().all(|(_, s)| *s == f) {
                        return Some(f);
                    }
                }
                None
            })
            .collect();
        let Some(exons) = exons else {
            continue; // cross-locus borrow → coords unknowable here → don't emit
        };
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

    if std::env::var_os("RUSTLE_LAYER2_DEBUG").is_some() {
        eprintln!(
            "[layer2-dbg] {} cross-map link(s) over {} loci; min_link={min_link} min_sim={min_similarity}",
            links.len(),
            loci.len()
        );
        for ((a, b), count) in &links {
            let sim = similarity.get(&(*a, *b)).copied();
            let na = loci.get(*a).map(|l| l.graph.nodes.len()).unwrap_or(0);
            let nb = loci.get(*b).map(|l| l.graph.nodes.len()).unwrap_or(0);
            eprintln!(
                "[layer2-dbg]   link ({a},{b}) count={count} sim={sim:?} graph_nodes a={na} b={nb}"
            );
        }
        eprintln!("[layer2-dbg] {} families formed", families.len());
    }

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
            t.vg_family_size = Some(fam.bundle_indices.len());
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

// =============================================================================
// Homology coordinate-transfer (C5 building blocks, PURE — no wiring yet).
//
// Recover a STARVED paralog copy's exon coordinates from its OWN secondary
// alignments (those that mapped to its native region but whose PRIMARY went to
// the well-expressed paralog). NO coordinate is ever invented: every emitted
// coordinate is a real secondary's exon, anchored by the starved copy's
// surviving native exon(s).
// =============================================================================

/// Reconstruct each secondary's exon chain from `(ref_start, introns, ref_end)`
/// and return the CONSENSUS exon chain for this group of secondaries (all from
/// ONE copy's locus). Deterministic. Returns empty if there is no supported
/// consensus (never fabricate from a single noisy read).
///
/// Junctions are crisp splice sites (group key); exon ends jitter (median).
/// The winning chain needs >=2 supporters; interior boundaries are the chosen
/// junction coords; the 5' start and 3' end are the median over supporters
/// (even count → lower-middle element).
fn consensus_exon_chain(
    secondaries: &[&crate::vg_family::secondary_index::SecondaryAlignment],
) -> Vec<(u64, u64)> {
    if secondaries.is_empty() {
        return Vec::new();
    }
    // Group by full junction-chain identity (the introns tuple).
    let mut groups: DetHashMap<Vec<(u64, u64)>, Vec<&crate::vg_family::secondary_index::SecondaryAlignment>> =
        DetHashMap::default();
    for s in secondaries {
        groups.entry(s.introns.clone()).or_default().push(s);
    }
    // Pick the MODE junction-chain. Deterministic tie-break:
    //   (a) most support, (b) longest chain (more junctions),
    //   (c) lexicographically smallest intron tuple.
    let mut best: Option<&Vec<(u64, u64)>> = None;
    let mut best_support = 0usize;
    for (chain, members) in &groups {
        let support = members.len();
        let take = match best {
            None => true,
            Some(b) => {
                let b_support = groups[b].len();
                (support, chain.len(), std::cmp::Reverse(chain))
                    > (b_support, b.len(), std::cmp::Reverse(b))
            }
        };
        if take {
            best = Some(chain);
            best_support = support;
        }
    }
    let chain = match best {
        Some(c) if best_support >= 2 => c.clone(),
        _ => return Vec::new(), // no consensus (single read / unsupported)
    };
    let supporters = &groups[&chain];

    // Terminal exon coords vary across supporters → MEDIAN (lower-middle).
    let median = |mut v: Vec<u64>| -> u64 {
        v.sort_unstable();
        v[(v.len() - 1) / 2]
    };
    let start = median(supporters.iter().map(|s| s.ref_start).collect());
    let end = median(supporters.iter().map(|s| s.ref_end).collect());

    // Build exons: start .. (donor) , (acceptor) .. (next donor) , ... , (last acceptor) .. end.
    // Interior boundaries are exactly the junction coords (crisp).
    let mut exons = Vec::with_capacity(chain.len() + 1);
    let mut cur_start = start;
    for &(donor, acceptor) in &chain {
        exons.push((cur_start, donor));
        cur_start = acceptor;
    }
    exons.push((cur_start, end));
    exons
}

/// Map a starved copy's consensus exon chain onto the family backbone nodes by
/// CO-LINEAR position, pinned by the copy's surviving NATIVE exon(s) (its
/// `per_copy_spans` entries). Returns `(node_idx, transferred_exon_span)` for
/// each backbone node the copy is placed on. Requires >=1 native anchor (else
/// the frame is unpinnable → returns empty; never guesses the frame). Returns
/// empty on any frame inconsistency rather than emitting a wrong-frame mapping.
fn positional_mapping(
    fg: &crate::vg_family::family_graph::FamilyGraph,
    starved_copy_id: usize,
    consensus: &[(u64, u64)],
) -> Vec<(usize, (u64, u64))> {
    if consensus.is_empty() {
        return Vec::new();
    }
    // 1. Backbone order: node indices sorted by genomic start. Stable tie-break
    //    on (span, idx) keeps it total-order deterministic.
    let mut order: Vec<usize> = (0..fg.nodes.len()).collect();
    order.sort_by(|&a, &b| {
        fg.nodes[a]
            .span
            .cmp(&fg.nodes[b].span)
            .then(a.cmp(&b))
    });
    if order.is_empty() {
        return Vec::new();
    }

    // 2. Starved copy's NATIVE coords per backbone rank.
    let mut native_at_rank: DetHashMap<usize, (u64, u64)> = DetHashMap::default();
    for (rank, &node_idx) in order.iter().enumerate() {
        for &(c, sp) in &fg.nodes[node_idx].per_copy_spans {
            if c == starved_copy_id {
                native_at_rank.insert(rank, sp);
            }
        }
    }
    if native_at_rank.is_empty() {
        return Vec::new(); // C5: refuse to guess the frame.
    }

    // 3. Anchor the consensus to the backbone. For every (rank, native span)
    //    find a consensus exon that overlaps it; the implied offset is
    //    `rank - ci`. Use the FIRST/lowest anchor deterministically; any
    //    second anchor that implies a DIFFERENT offset = frame inconsistency.
    let overlaps =
        |a: (u64, u64), b: (u64, u64)| a.0 < b.1 && b.0 < a.1;
    let mut ranks: Vec<usize> = native_at_rank.keys().copied().collect();
    ranks.sort_unstable();
    let mut offset: Option<i64> = None;
    for &rank in &ranks {
        let native = native_at_rank[&rank];
        for (ci, &cons) in consensus.iter().enumerate() {
            if overlaps(cons, native) {
                let off = rank as i64 - ci as i64;
                match offset {
                    None => offset = Some(off),
                    Some(o) if o != off => return Vec::new(), // inconsistent frame
                    Some(_) => {}
                }
            }
        }
    }
    let offset = match offset {
        Some(o) => o,
        None => return Vec::new(), // native exists but no consensus exon overlaps it
    };

    // 4. Map each consensus exon ci → backbone rank r = ci + offset → order[r].
    //    Validate every exon lands in range (contiguous co-linear by construction).
    let mut out: Vec<(usize, (u64, u64))> = Vec::with_capacity(consensus.len());
    for (ci, &cons) in consensus.iter().enumerate() {
        let r = ci as i64 + offset;
        if r < 0 || r as usize >= order.len() {
            return Vec::new(); // maps out of range → frame inconsistent.
        }
        out.push((order[r as usize], cons));
    }

    // 6. Determinism: sort by node_idx.
    out.sort_by_key(|&(node_idx, _)| node_idx);
    out
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
    fn decompose_recovers_native_copy_defers_crosslocus_borrow() {
        // copy 0 is fully native (all 3 exons survived Layer 1); copy 1 is starved
        // (only the first exon survived) and would borrow exons 2+3 via secondaries.
        // The well copy MUST recover its full chain at ITS OWN coordinates. The
        // starved copy's borrowed nodes are single-entry (copy-0-only), so its true
        // coordinates are unknowable without homology coordinate-transfer (future
        // work) — it is conservatively DROPPED rather than emitted at copy 0's coords.
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
        // The native copy recovers its full chain at its own coordinates.
        assert!(
            paths.iter().any(|p| p.copy_id == 0
                && p.exons == vec![(100, 160), (300, 360), (500, 560)]),
            "native copy 0 full chain recovered: {paths:?}"
        );
        // No path is emitted with FABRICATED coordinates for the starved copy: every
        // emitted exon is a real per-copy span, never a node union or a borrowed
        // foreign coordinate. (Cross-locus starved recovery is deferred.)
        assert!(
            paths.iter().all(|p| p.exons.iter().all(|&(s, e)| e - s == 60)),
            "no fabricated/union-span exon (all exons are real 60bp spans): {paths:?}"
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

    // -------------------------------------------------------------------------
    // T1 — consensus_exon_chain
    // -------------------------------------------------------------------------

    fn sec_owned(
        start: u64, end: u64, introns: &[(u64, u64)],
    ) -> crate::vg_family::secondary_index::SecondaryAlignment {
        sec_aln(0, start, end, introns, 0)
    }

    #[test]
    fn consensus_full_agreement_returns_exact_chain() {
        let s = [
            sec_owned(30000, 30660, &[(30060, 30300), (30360, 30600)]),
            sec_owned(30000, 30660, &[(30060, 30300), (30360, 30600)]),
            sec_owned(30000, 30660, &[(30060, 30300), (30360, 30600)]),
        ];
        let refs: Vec<&_> = s.iter().collect();
        let got = consensus_exon_chain(&refs);
        assert_eq!(
            got,
            vec![(30000, 30060), (30300, 30360), (30600, 30660)],
            "3 identical secondaries reconstruct that exact exon chain"
        );
    }

    #[test]
    fn consensus_jitter_uses_median_terminals_exact_junctions() {
        // Same introns; ref_start/ref_end jitter a few bp. Interior boundaries are
        // junction-anchored (crisp); terminals are the MEDIAN of the supporters.
        let s = [
            sec_owned(29998, 30658, &[(30060, 30300), (30360, 30600)]),
            sec_owned(30000, 30660, &[(30060, 30300), (30360, 30600)]),
            sec_owned(30003, 30662, &[(30060, 30300), (30360, 30600)]),
        ];
        let refs: Vec<&_> = s.iter().collect();
        let got = consensus_exon_chain(&refs);
        assert_eq!(
            got,
            vec![(30000, 30060), (30300, 30360), (30600, 30660)],
            "median start=30000, median end=30660; interior junctions exact"
        );
    }

    #[test]
    fn consensus_disagreement_returns_mode() {
        // 3 of chain X, 2 of chain Y (different introns) → X wins on support.
        let s = [
            sec_owned(100, 560, &[(160, 300), (360, 500)]),
            sec_owned(100, 560, &[(160, 300), (360, 500)]),
            sec_owned(100, 560, &[(160, 300), (360, 500)]),
            sec_owned(100, 560, &[(160, 320), (380, 500)]),
            sec_owned(100, 560, &[(160, 320), (380, 500)]),
        ];
        let refs: Vec<&_> = s.iter().collect();
        let got = consensus_exon_chain(&refs);
        assert_eq!(
            got,
            vec![(100, 160), (300, 360), (500, 560)],
            "mode chain X recovered, not Y"
        );
    }

    #[test]
    fn consensus_single_read_returns_empty() {
        let s = [sec_owned(100, 560, &[(160, 300), (360, 500)])];
        let refs: Vec<&_> = s.iter().collect();
        assert!(
            consensus_exon_chain(&refs).is_empty(),
            "support<2 → no consensus (never fabricate from one read)"
        );
    }

    #[test]
    fn consensus_empty_input_returns_empty() {
        assert!(consensus_exon_chain(&[]).is_empty());
    }

    // -------------------------------------------------------------------------
    // T2 — positional_mapping
    // -------------------------------------------------------------------------

    /// Backbone of `spans.len()` nodes; for each (rank, span) the listed copies
    /// are recorded in per_copy_spans at the GIVEN per-copy coordinate. `well`
    /// supplies the node's union `span` (the well copy's coordinate).
    fn mk_backbone(
        well: &[(u64, u64)],
        per_copy: &[Vec<(usize, (u64, u64))>],
    ) -> FamilyGraph {
        use crate::vg_family::family_graph::{ExonClass, NodeIdx};
        let nodes = well
            .iter()
            .enumerate()
            .map(|(i, &span)| ExonClass {
                idx: NodeIdx(i),
                chrom: "chrT".into(),
                span,
                strand: '+',
                per_copy_sequences: per_copy[i].iter().map(|(c, _)| (*c, Vec::new())).collect(),
                per_copy_spans: per_copy[i].clone(),
                copy_specific: per_copy[i].len() == 1,
                per_copy_cov: Vec::new(),
            })
            .collect();
        FamilyGraph { family_id: 0, nodes, edges: Vec::new() }
    }

    #[test]
    fn positional_cross_locus_returns_starved_copy_coords() {
        // Well copy (0) at backbone coords; starved copy (1) native ONLY at node 0
        // (its real coord 30000-30060). Consensus = the starved copy's own 3-exon
        // chain. Mapping must place exons 2+3 onto node1/node2 at the STARVED
        // copy's coords, NOT the well copy's backbone coords.
        let fg = mk_backbone(
            &[(100, 160), (300, 360), (500, 560)],
            &[
                vec![(0, (100, 160)), (1, (30000, 30060))],
                vec![(0, (300, 360))],
                vec![(0, (500, 560))],
            ],
        );
        let consensus = [(30000, 30060), (30300, 30360), (30600, 30660)];
        let got = positional_mapping(&fg, 1, &consensus);
        assert_eq!(
            got,
            vec![
                (0, (30000, 30060)),
                (1, (30300, 30360)),
                (2, (30600, 30660)),
            ],
            "transferred spans are the STARVED copy's coords, not the well copy's"
        );
    }

    #[test]
    fn positional_no_anchor_returns_empty() {
        // Starved copy (1) has zero per_copy_spans entries anywhere.
        let fg = mk_backbone(
            &[(100, 160), (300, 360), (500, 560)],
            &[vec![(0, (100, 160))], vec![(0, (300, 360))], vec![(0, (500, 560))]],
        );
        let consensus = [(30000, 30060), (30300, 30360), (30600, 30660)];
        assert!(
            positional_mapping(&fg, 1, &consensus).is_empty(),
            "no native anchor → frame unpinnable → empty (never guess)"
        );
    }

    #[test]
    fn positional_overflow_returns_empty() {
        // Anchor at node 0 fixes offset 0, but the consensus has more exons than
        // the backbone can hold → frame inconsistent → empty.
        let fg = mk_backbone(
            &[(100, 160), (300, 360)],
            &[vec![(0, (100, 160)), (1, (30000, 30060))], vec![(0, (300, 360))]],
        );
        let consensus = [(30000, 30060), (30300, 30360), (30600, 30660)];
        assert!(
            positional_mapping(&fg, 1, &consensus).is_empty(),
            "consensus longer than backbone from the anchor → empty"
        );
    }

    #[test]
    fn positional_inconsistent_offsets_return_empty() {
        // Starved copy native at node 0 AND node 2, but the consensus overlaps
        // those native spans at offsets that disagree → frame inconsistent → empty.
        // node0 native (30000,30060) overlaps consensus[0] → offset 0.
        // node2 native (30600,30660) overlaps consensus[1] (30550,30660) → offset 1.
        let fg = mk_backbone(
            &[(100, 160), (300, 360), (500, 560)],
            &[
                vec![(0, (100, 160)), (1, (30000, 30060))],
                vec![(0, (300, 360))],
                vec![(0, (500, 560)), (1, (30600, 30660))],
            ],
        );
        let consensus = [(30000, 30060), (30550, 30660), (30900, 30960)];
        assert!(
            positional_mapping(&fg, 1, &consensus).is_empty(),
            "anchors imply inconsistent offsets → empty (no wrong-frame map)"
        );
    }
}

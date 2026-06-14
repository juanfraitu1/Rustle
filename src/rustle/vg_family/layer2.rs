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
    /// Homology-transferred per-copy coordinates: `(node_idx, copy_id, exon_span)`.
    /// Populated when a starved copy's exon coordinates are recovered from its OWN
    /// secondaries' consensus and positionally mapped onto the backbone. Each span
    /// is a REAL secondary's exon at the STARVED copy's coordinates — never the well
    /// copy's backbone span. Consumed by `decompose_family_paths` to emit the copy.
    pub transferred_coords: Vec<(usize, usize, (u64, u64))>,
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
    // Homology coordinate-transfer (C5): for each family copy, reconstruct its OWN
    // exon coordinates from the consensus of its secondaries and positionally map
    // them onto the backbone. This recovers a STARVED copy whose Layer-1 graph was
    // empty — its real coordinates live only in its secondaries. Every transferred
    // span is a real secondary's exon; positional_mapping refuses ambiguous frames.
    let mut by_copy: DetHashMap<usize, Vec<&crate::vg_family::secondary_index::SecondaryAlignment>> =
        DetHashMap::default();
    for s in secondaries {
        if let Some(copy) = s.locus.and_then(|l| copy_of_locus.get(&l).copied()) {
            by_copy.entry(copy).or_default().push(s);
        }
    }
    let mut copies: Vec<usize> = by_copy.keys().copied().collect();
    copies.sort_unstable();
    let mut transferred_coords: Vec<(usize, usize, (u64, u64))> = Vec::new();
    for copy in copies {
        let secs = &by_copy[&copy];
        let consensus = consensus_exon_chain(secs);
        if consensus.is_empty() {
            continue;
        }
        for (node_idx, span) in positional_mapping(fg, copy, &consensus) {
            transferred_coords.push((node_idx, copy, span));
            membership.insert((node_idx, copy));
        }
    }
    transferred_coords.sort_unstable();

    let mut edges: Vec<AmendCandidate> = support
        .into_iter()
        .map(|(intron, support)| AmendCandidate { intron, support })
        .collect();
    edges.sort_by(|x, y| x.intron.cmp(&y.intron));
    let mut copy_membership: Vec<(usize, usize)> = membership.into_iter().collect();
    copy_membership.sort_unstable();
    GraphAmendment { edges, copy_membership, transferred_coords }
}

/// Where a recovered isoform came from. `Native` = the copy's own reads (base
/// decompose path or its own enumerated alt-splice chains). `Transferred` = a
/// donor copy's isoform mapped onto this copy and confirmed by the recipient's
/// own per-junction evidence ("transfer proposes, recipient disposes").
#[derive(Debug, Clone, PartialEq)]
pub enum IsoformSource {
    Native,
    Transferred { donor_copy: usize },
}

/// A recovered copy/isoform path (exon chain in genomic coordinates). A path that
/// violates allele-linkage (a copy claiming a diagnostic node it has no allele for)
/// is DROPPED during decomposition — never emitted.
#[derive(Debug, Clone, PartialEq)]
pub struct FamilyPath {
    pub exons: Vec<(u64, u64)>,
    pub flow: f64,
    pub copy_id: usize,
    /// (prev_exon.end, next_exon.start) per junction — always derived from
    /// `exons` via `junctions_of`. Carried so dedup / union-by-chain can key on
    /// the intron chain without recomputing it.
    pub junction_chain: Vec<(u64, u64)>,
    pub source: IsoformSource,
}

/// The intron chain (junctions) implied by an exon chain: one `(donor, acceptor)`
/// per adjacent exon pair. Single-exon chains have no junctions (empty). This is
/// the ONE canonical derivation — every FamilyPath's `junction_chain` comes from
/// here so the field can never drift out of sync with `exons`.
fn junctions_of(exons: &[(u64, u64)]) -> Vec<(u64, u64)> {
    exons.windows(2).map(|w| (w[0].1, w[1].0)).collect()
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
    // Copy count must include copies introduced ONLY via the amendment — a STARVED
    // copy whose Layer-1 graph was empty contributes no `per_copy_sequences` to any
    // node, so it would be invisible if we counted graph nodes alone. Its homology-
    // transferred coordinates / memberships are what re-introduce it.
    let n_copies = fg
        .nodes
        .iter()
        .flat_map(|n| n.per_copy_sequences.iter().map(|(c, _)| *c))
        .chain(amendment.copy_membership.iter().map(|(_, c)| *c))
        .chain(amendment.transferred_coords.iter().map(|(_, c, _)| *c))
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
    // Homology-transferred per-copy coordinates: (node_idx, copy_id) → recovered
    // exon span (a real secondary's exon at the STARVED copy's own coordinates).
    let mut transferred: DetHashMap<(usize, usize), (u64, u64)> = DetHashMap::default();
    for &(ni, c, sp) in &amendment.transferred_coords {
        transferred.insert((ni, c), sp);
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
                // Homology-transferred coordinate (NEW): the starved copy's OWN exon
                // span, recovered from its secondaries' consensus and positionally
                // mapped onto this backbone node. Takes precedence over the >=2-agree
                // co-located borrow because it carries the copy's TRUE coordinates.
                if let Some(sp) = transferred.get(&(i, copy)) {
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
        let junction_chain = junctions_of(&exons);
        paths.push(FamilyPath {
            exons,
            flow,
            copy_id: copy,
            junction_chain,
            source: IsoformSource::Native,
        });
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
    // Secondary-enriched isoform discovery (T4 wiring). Part A (per-copy alt-splice
    // from own secondaries) is ON by default under --vg-layer2 (additive; emit only
    // unions chains absent from VG output). Opt out for bisection. Part B (cross-copy
    // transfer) stays OFF until real-data gffcompare proves it net-positive.
    let isoform_min_support: usize = std::env::var("RUSTLE_VG_LAYER2_MIN_ISOFORM_SUPPORT")
        .ok().and_then(|v| v.parse::<usize>().ok()).filter(|&k| k >= 1).unwrap_or(2);
    let enable_multi_isoform = std::env::var_os("RUSTLE_VG_LAYER2_NO_MULTI_ISOFORM").is_none();
    let enable_transfer = std::env::var_os("RUSTLE_VG_LAYER2_ISOFORM_TRANSFER").is_some();

    // (C2) candidate links → compute similarity only for linked, same-chrom pairs.
    let links = crate::vg::build_multimap_index_from_secondary_index(&side_index, primary_locus);
    debug_assert!(
        links.iter().all(|((a, b), _)| *a < loci.len() && *b < loci.len()),
        "run_layer2: link locus index out of range of loci slice (caller contract)"
    );
    // Derive a locus's exon SEQUENCES with a starved-copy fallback: prefer the
    // Layer-1 graph nodes; if the graph is empty (the starved copy whose primaries
    // cross-mapped away → 0 exon nodes), reconstruct the copy's exon seqs from the
    // CONSENSUS of its OWN secondaries' coordinates. This keeps the SAME k-mer-
    // Jaccard threshold — it only sources the starved copy's real exon sequence from
    // its reads instead of its (empty) graph, so families can still form.
    let locus_exon_seqs = |locus_idx: usize, g: &GenomeIndex| -> Vec<Vec<u8>> {
        let l = &loci[locus_idx];
        let mut seqs = crate::vg::graph_exon_seqs(l.graph, &l.chrom, g);
        if seqs.is_empty() {
            let secs = side_index.secondaries_for_locus(locus_idx);
            let consensus = consensus_exon_chain(&secs);
            seqs = consensus
                .iter()
                .filter_map(|&(s, e)| g.fetch_sequence(&l.chrom, s, e))
                .collect();
        }
        seqs
    };

    let mut similarity: DetHashMap<(usize, usize), f64> = DetHashMap::default();
    if let Some(g) = genome {
        for ((a, b), count) in &links {
            if *count < min_link {
                continue;
            }
            // Defensive: a cross-map link index must be in range of `loci`
            // (the caller re-syncs side-index loci to the loci snapshot, but never
            // index out of range even if that contract is violated).
            let (Some(la), Some(lb)) = (loci.get(*a), loci.get(*b)) else {
                continue;
            };
            if la.chrom == lb.chrom {
                let sim = crate::vg::exon_kmer_similarity(
                    &locus_exon_seqs(*a, g),
                    &locus_exon_seqs(*b, g),
                    kmer_len,
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
            .filter_map(|&i| loci.get(i).map(|l| (l.chrom.clone(), l.strand, l.graph)))
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
        // Group secondaries by the copy of their ALIGNMENT locus — the SAME grouping
        // amend_family_graph uses internally (by_copy). Satisfies emit's "own secondaries
        // only" precondition by construction: each bucket holds exactly the secondaries
        // that align at that copy's locus.
        let mut secondaries_by_copy: DetHashMap<usize, Vec<&crate::vg_family::secondary_index::SecondaryAlignment>> = DetHashMap::default();
        for s in &secs {
            if let Some(copy) = s.locus.and_then(|l| copy_of_locus.get(&l).copied()) {
                secondaries_by_copy.entry(copy).or_default().push(s);
            }
        }
        // min_flow stays 1.0 so the base decompose behavior is unchanged (preserves
        // starved-copy + baseline-floor recovery — legs 6/7). Part A is independently
        // gated at >=2 support by K + the G7 isofrac floor (max(2.0, 0.01*copy_max)),
        // which already subsumes a 2.0 min_flow, so 1.0 here loses no Part-A precision.
        let paths = emit_family_isoforms(
            fg, &am, &secondaries_by_copy, 1.0,
            isoform_min_support, enable_multi_isoform, enable_transfer,
        );
        if std::env::var_os("RUSTLE_LAYER2_DEBUG").is_some() {
            eprintln!(
                "[layer2-dbg] fam {} fg_nodes={} secs={} transferred={} paths={}",
                fam.family_id, fg.nodes.len(), secs.len(),
                am.transferred_coords.len(), paths.len()
            );
            for (ni, c, sp) in &am.transferred_coords {
                eprintln!("[layer2-dbg]   transferred node={ni} copy={c} span={sp:?}");
            }
            for p in &paths {
                eprintln!(
                    "[layer2-dbg]   path copy={} flow={} source={:?} exons={:?}",
                    p.copy_id, p.flow, p.source, p.exons
                );
            }
        }
        for p in paths {
            // Use the recovered copy's OWN bundle for chrom/strand (not the family's
            // first), falling back to the first member; guarded against any stale idx.
            let anchor = fam
                .bundle_indices
                .get(p.copy_id)
                .copied()
                .or_else(|| fam.bundle_indices.first().copied())
                .and_then(|i| loci.get(i));
            let Some(anchor) = anchor else { continue };
            let chrom = anchor.chrom.clone();
            let strand = anchor.strand;
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
    build_exons_from_chain(supporters, &chain)
}

/// Build exon spans for a known intron chain from its supporting secondaries:
/// interior boundaries are the exact junctions; terminal exon start/end are the
/// MEDIAN of supporters' ref_start/ref_end (lower-middle, deterministic).
///
/// Factored out of `consensus_exon_chain` (behavior-preserving): given the chain
/// and the secondaries that traversed it, produce the consensus exon coordinates.
fn build_exons_from_chain(
    supporters: &[&crate::vg_family::secondary_index::SecondaryAlignment],
    chain: &[(u64, u64)],
) -> Vec<(u64, u64)> {
    if supporters.is_empty() {
        return Vec::new();
    }
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
    for &(donor, acceptor) in chain {
        exons.push((cur_start, donor));
        cur_start = acceptor;
    }
    exons.push((cur_start, end));
    exons
}

/// Return EVERY distinct intron-chain carried by >= K of these secondaries (one
/// copy's), each with its consensus exons + support count. Deterministic:
/// sorted by (support DESC, chain.len() DESC, chain ASC). The single most
/// important false-positive defense is that we only ever return chains a real
/// molecule fully traversed >= K times — never a synthesized/graph-walk chain.
fn enumerate_secondary_chains(
    secondaries: &[&crate::vg_family::secondary_index::SecondaryAlignment],
    k: usize,
) -> Vec<(Vec<(u64, u64)> /*exons*/, Vec<(u64, u64)> /*intron chain*/, usize /*support*/)> {
    debug_assert!(k >= 1, "enumerate_secondary_chains: K must be >= 1");
    if secondaries.is_empty() {
        return Vec::new();
    }
    // Group by full junction-chain identity (the introns tuple) — the SAME
    // grouping consensus_exon_chain uses.
    let mut groups: DetHashMap<
        Vec<(u64, u64)>,
        Vec<&crate::vg_family::secondary_index::SecondaryAlignment>,
    > = DetHashMap::default();
    for s in secondaries {
        groups.entry(s.introns.clone()).or_default().push(s);
    }
    let mut out: Vec<(Vec<(u64, u64)>, Vec<(u64, u64)>, usize)> = Vec::new();
    for (chain, members) in &groups {
        let support = members.len();
        if support < k {
            continue; // never return a chain a real molecule didn't traverse >= K times
        }
        let exons = build_exons_from_chain(members, chain);
        out.push((exons, chain.clone(), support));
    }
    // Deterministic order: support DESC, chain.len() DESC, chain ASC.
    out.sort_by(|a, b| {
        b.2.cmp(&a.2)
            .then(b.1.len().cmp(&a.1.len()))
            .then(a.1.cmp(&b.1))
    });
    out
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
        // No native anchor. Only one case is unambiguous: the consensus spans the
        // WHOLE backbone (exact length match). Then co-linear order pins the frame
        // with no anchor needed — consensus exon `i` → `order[i]` at the STARVED
        // copy's OWN consensus coordinate (never the well copy's backbone span).
        // Any other length mismatch leaves the frame ambiguous → refuse.
        if consensus.len() == fg.nodes.len() {
            let mut out: Vec<(usize, (u64, u64))> = consensus
                .iter()
                .enumerate()
                .map(|(ci, &cons)| (order[ci], cons))
                .collect();
            out.sort_by_key(|&(node_idx, _)| node_idx);
            return out;
        }
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

// =============================================================================
// Isoform transfer across copies (C5 building blocks, PURE — no wiring yet).
//
// "Transfer proposes, recipient disposes." A donor copy's isoform (a real exon
// chain) is mapped onto a recipient copy's OWN coordinates via the co-linear
// backbone, then admitted ONLY if the recipient has its own per-junction
// evidence. We never emit the donor's coordinates for the recipient, and never
// claim a node the recipient cannot legitimately place.
// =============================================================================

/// Backbone node indices in canonical co-linear order (by span, then idx) — the
/// SAME order positional_mapping uses. Rank = position in this Vec.
fn canonical_node_order(fg: &FamilyGraph) -> Vec<usize> {
    let mut order: Vec<usize> = (0..fg.nodes.len()).collect();
    order.sort_by(|&a, &b| fg.nodes[a].span.cmp(&fg.nodes[b].span).then(a.cmp(&b)));
    order
}

/// Map donor copy A's isoform (given as A's exon spans) onto recipient copy B's
/// COORDINATES via the co-linear backbone. For each of A's exons: find its
/// backbone node (A's per_copy_spans match), take that node's rank, and emit B's
/// coordinate at that node — B's own per_copy_spans entry, ELSE a >=2-agree
/// co-located shared span (the borrow rule from decompose), ELSE abort (return
/// empty: B can't legitimately claim a copy-A-unique node). Result exons must be
/// strictly genomically ordered (co-linear) else abort. NEVER emits A's coords.
fn map_isoform_across_copies(
    fg: &FamilyGraph,
    donor_copy: usize,
    donor_exons: &[(u64, u64)],
    recipient_copy: usize,
) -> Vec<(u64, u64)> {
    debug_assert!(
        donor_copy != recipient_copy,
        "map_isoform_across_copies: donor and recipient must differ"
    );
    if donor_exons.is_empty() {
        return Vec::new();
    }
    let order = canonical_node_order(fg);
    // Rank of each backbone node (position in canonical order).
    let mut rank_of_node: DetHashMap<usize, usize> = DetHashMap::default();
    for (rank, &node_idx) in order.iter().enumerate() {
        rank_of_node.insert(node_idx, rank);
    }
    let mut placed: Vec<(usize /*rank*/, (u64, u64) /*recipient coord*/)> =
        Vec::with_capacity(donor_exons.len());
    for &dex in donor_exons {
        // 1. Locate A's backbone node: the node whose per_copy_spans entry for donor A
        //    OVERLAPS this donor exon. Exact equality is too strict — the donor isoform's
        //    exons are secondary-consensus coords (median terminal ends + crisp junctions),
        //    which essentially never equal the Layer-1-assembled node span to the base. The
        //    node a donor exon OVERLAPS is the node it belongs to; the recipient still
        //    supplies its OWN coordinate below (we never emit A's coords). Pick the largest
        //    overlap so a donor exon grazing a neighbor maps to its true node.
        let overlaps = |a: (u64, u64), b: (u64, u64)| a.0 < b.1 && b.0 < a.1;
        let node_idx = fg
            .nodes
            .iter()
            .enumerate()
            .filter_map(|(i, n)| {
                n.per_copy_spans
                    .iter()
                    .find(|&&(c, _)| c == donor_copy)
                    .map(|&(_, sp)| (i, sp))
            })
            .filter(|&(_, sp)| overlaps(sp, dex))
            .max_by_key(|&(_, sp)| dex.1.min(sp.1).saturating_sub(dex.0.max(sp.0)))
            .map(|(i, _)| i);
        let Some(node_idx) = node_idx else {
            return Vec::new(); // donor exon overlaps no backbone node for A → abort
        };
        let Some(&rank) = rank_of_node.get(&node_idx) else {
            return Vec::new();
        };
        let n = &fg.nodes[node_idx];
        // 2. Recipient B's coordinate at this node:
        //    (a) B's own per_copy_spans entry → that coordinate.
        //    (b) ELSE a >=2-agree co-located shared span (genuine co-location).
        //    (c) ELSE abort: copy-A-unique node, B cannot legitimately claim it.
        let recip_coord = if let Some((_, sp)) =
            n.per_copy_spans.iter().find(|(c, _)| *c == recipient_copy)
        {
            *sp
        } else if n.per_copy_spans.len() >= 2 {
            let f = n.per_copy_spans[0].1;
            if n.per_copy_spans.iter().all(|(_, s)| *s == f) {
                f
            } else {
                return Vec::new();
            }
        } else {
            return Vec::new();
        };
        placed.push((rank, recip_coord));
    }
    // 3. Co-linearity: emit in backbone-rank order; recipient coords must be
    //    strictly genomically ordered (non-overlapping, ascending) else abort.
    placed.sort_by_key(|&(rank, _)| rank);
    let exons: Vec<(u64, u64)> = placed.iter().map(|&(_, sp)| sp).collect();
    for w in exons.windows(2) {
        // each exon's end must not exceed the next exon's start (strict order)
        if w[0].1 > w[1].0 {
            return Vec::new();
        }
    }
    // Each exon must itself be a valid (start < end) span.
    if exons.iter().any(|&(s, e)| s >= e) {
        return Vec::new();
    }
    exons
}

/// "Transfer proposes, recipient disposes": a transferred isoform is admissible
/// ONLY if the recipient B has its OWN evidence for EVERY junction of the
/// recipient chain — i.e. each (donor,acceptor) junction of `recipient_exons`
/// appears (within +-5bp jitter, matching collect_family_junctions) among B's own
/// secondaries' introns OR among fg.edges between B-member nodes. Returns true iff
/// all junctions are recipient-supported with >= K support.
fn verify_recipient_support(
    fg: &FamilyGraph,
    recipient_copy: usize,
    recipient_exons: &[(u64, u64)],
    recipient_secondaries: &[&crate::vg_family::secondary_index::SecondaryAlignment],
    k: usize,
) -> bool {
    // The recipient chain's junctions: (prev_exon.end, next_exon.start).
    let junctions: Vec<(u64, u64)> = recipient_exons
        .windows(2)
        .map(|w| (w[0].1, w[1].0))
        .collect();
    if junctions.is_empty() {
        return false; // a single-exon "chain" has no junction to verify
    }
    // Splice sites are crisp; allow the small (<=5 bp) jitter typical of noisy
    // long-read soft-clip/alignment tails when matching a recipient junction.
    let close = |a: u64, b: u64| -> bool { a.max(b) - a.min(b) <= 5 };
    // Nodes that recipient B is a member of (for the fg.edges path).
    let mut b_member_node: DetHashSet<usize> = DetHashSet::default();
    for (ni, n) in fg.nodes.iter().enumerate() {
        if n.per_copy_spans.iter().any(|(c, _)| *c == recipient_copy)
            || n.per_copy_sequences.iter().any(|(c, _)| *c == recipient_copy)
        {
            b_member_node.insert(ni);
        }
    }
    for &(jd, ja) in &junctions {
        // (A) Count B's own secondaries whose intron chain carries this junction.
        let mut support = 0usize;
        for s in recipient_secondaries {
            if s.introns.iter().any(|&(d, a)| close(d, jd) && close(a, ja)) {
                support += 1;
            }
        }
        // (B) ELSE fg.edges between B-native nodes whose donor/acceptor match — at
        // the RECIPIENT's OWN coordinate (per_copy_spans), NOT the node's union
        // `span` (which, for dispersed paralogs at different loci, is genomically
        // nonsensical and would never match the recipient's junction). An edge
        // contributes only when B is natively placed at BOTH endpoints.
        if support < k {
            for e in &fg.edges {
                if !b_member_node.contains(&e.from.0) || !b_member_node.contains(&e.to.0) {
                    continue;
                }
                let de = fg.nodes[e.from.0]
                    .per_copy_spans
                    .iter()
                    .find(|(c, _)| *c == recipient_copy)
                    .map(|(_, sp)| sp.1);
                let ae = fg.nodes[e.to.0]
                    .per_copy_spans
                    .iter()
                    .find(|(c, _)| *c == recipient_copy)
                    .map(|(_, sp)| sp.0);
                if let (Some(de), Some(ae)) = (de, ae) {
                    if close(de, jd) && close(ae, ja) {
                        support += e.family_support as usize;
                    }
                }
            }
        }
        if support < k {
            return false; // recipient lacks its own evidence for this junction
        }
    }
    true
}

/// Compose the validated per-copy base (`decompose_family_paths`) with two additive
/// discovery channels, under strict false-positive gates:
///   Part A (enable_multi_isoform): per-copy ALT-SPLICE isoforms from the copy's OWN
///     secondaries (enumerate_secondary_chains). More isoforms, never a synthesized
///     chain — only chains a real molecule traversed >= k times.
///   Part B (enable_transfer): cross-copy TRANSFER — a donor copy's isoform mapped
///     onto a recipient's own coordinates (map_isoform_across_copies) and admitted
///     ONLY if the recipient has its own per-junction evidence (verify_recipient_support).
///
/// REGRESSION ANCHOR: with both flags false, returns exactly decompose_family_paths(..)
/// (now carrying source=Native + junction_chain). All additions are gated off.
///
/// # Preconditions
/// - `secondaries_by_copy` keys are family `copy_id`s — the SAME copy indices used in
///   `fg.nodes[].per_copy_spans` / `per_copy_sequences`. Part A/B index the graph by
///   these ids, so a mismatched keyspace would map a copy's reads onto another copy's
///   coordinates.
/// - Each value MUST contain ONLY that copy's OWN secondaries. Part A labels every
///   enumerated chain `Native` with NO allele-linkage guard (a chain from a copy's own
///   reads is, by assumption, a real molecule at that copy's locus). So if a foreign /
///   cross-mapped secondary leaks into a copy's bucket, its chain is emitted as that
///   copy's native isoform — a SILENT false positive with no downstream gate to catch
///   it. This is a known historical hazard in this codebase: cross-mapped secondaries
///   contaminating per-copy graphs. The caller owns this partition; we trust it here.
#[allow(clippy::too_many_arguments)]
pub fn emit_family_isoforms(
    fg: &FamilyGraph,
    amendment: &GraphAmendment,
    secondaries_by_copy: &DetHashMap<
        usize,
        Vec<&crate::vg_family::secondary_index::SecondaryAlignment>,
    >,
    min_flow: f64,
    k: usize,
    enable_multi_isoform: bool,
    enable_transfer: bool,
) -> Vec<FamilyPath> {
    debug_assert!(k >= 1, "emit_family_isoforms: K must be >= 1");

    // 1. The VALIDATED base. These are Native and pushed FIRST so they win every
    //    dedup tie (a recipient's native isoform always beats an identical transfer).
    let mut out = decompose_family_paths(fg, amendment, min_flow);

    // Canonical copy iteration order — sort the keys once so Part A / Part B are
    // deterministic regardless of map insertion order.
    let mut copies: Vec<usize> = secondaries_by_copy.keys().copied().collect();
    copies.sort_unstable();

    // --- Part A: per-copy alt-splice isoforms from the copy's OWN secondaries ----
    // No allele-linkage guard here: a chain enumerated from the copy's OWN reads is a
    // real molecule at the copy's own locus, so decompose's diagnostic-node guard
    // (which defends against BORROWING another copy's node) does not apply. The
    // chains-only >= k rule is the false-positive defense.
    if enable_multi_isoform {
        for &copy in &copies {
            let secs = &secondaries_by_copy[&copy];
            let chains = enumerate_secondary_chains(secs, k);
            // Isoform-fraction floor is relative to THIS copy's strongest chain.
            let copy_max = chains.iter().map(|(_, _, sup)| *sup).max().unwrap_or(0) as f64;
            for (exons, chain, support) in &chains {
                let flow = *support as f64;
                // G6 (min_flow): respect the same flow floor the base decompose uses.
                if flow < min_flow {
                    continue;
                }
                // G7 (isofrac floor): an alt isoform must clear BOTH an absolute floor
                // (>= 2 supporters — never a near-singleton) and 1% of the copy's
                // dominant chain (suppresses dribble alongside a high-coverage copy).
                if flow < (2.0_f64).max(0.01 * copy_max) {
                    continue;
                }
                // Single-exon chains carry no junction → no informative alt-splice
                // signal here; they are noise in this channel.
                if exons.len() < 2 {
                    continue;
                }
                out.push(FamilyPath {
                    exons: exons.clone(),
                    flow,
                    copy_id: copy,
                    junction_chain: chain.clone(),
                    source: IsoformSource::Native,
                });
            }
        }
    }

    // --- Part B: cross-copy isoform transfer ("transfer proposes, recipient disposes")
    if enable_transfer {
        // Fallback for a recipient with no secondaries — borrows nothing donor- or
        // recipient-specific, so hoist it once instead of rebuilding it per recipient.
        let empty: Vec<&crate::vg_family::secondary_index::SecondaryAlignment> = Vec::new();
        // RUSTLE_LAYER2_DEBUG funnel: pinpoints WHERE transfers die (donor chains →
        // exact backbone-match of donor exons → map success → recipient verify). On
        // dispersed real paralogs the exact `sp == dex` match in map_isoform_across_copies
        // is the suspected bottleneck; this counts it directly.
        let dbg = std::env::var_os("RUSTLE_LAYER2_DEBUG").is_some();
        let (mut n_chains, mut n_chain_unmatched, mut n_map_attempt, mut n_map_ok, mut n_verify_ok) =
            (0usize, 0usize, 0usize, 0usize, 0usize);
        for &donor in &copies {
            let donor_secs = &secondaries_by_copy[&donor];
            // The donor's REAL isoforms (only chains a molecule traversed >= k times).
            let donor_chains = enumerate_secondary_chains(donor_secs, k);
            if dbg {
                n_chains += donor_chains.len();
                // Recipient-independent probe: does EVERY donor exon of each chain have
                // an EXACT backbone node (donor per_copy_spans sp == dex)? If not, map
                // returns empty for ALL recipients (the exact-match hypothesis).
                for (donor_exons, _c, _s) in &donor_chains {
                    let all_matched = donor_exons.iter().all(|&dex| {
                        fg.nodes.iter().any(|n| {
                            n.per_copy_spans.iter().any(|&(c, sp)| c == donor && sp == dex)
                        })
                    });
                    if !all_matched {
                        n_chain_unmatched += 1;
                    }
                }
            }
            for &recipient in &copies {
                if recipient == donor {
                    continue;
                }
                // Recipient's own secondaries (empty slice if it has none → the gate
                // can still pass via fg.edges, but never via the donor's reads).
                let recip_secs = secondaries_by_copy.get(&recipient).unwrap_or(&empty);
                for (donor_exons, _donor_chain, _support) in &donor_chains {
                    if dbg {
                        n_map_attempt += 1;
                    }
                    // Map donor's isoform onto the recipient's OWN coordinates.
                    let recip_exons = map_isoform_across_copies(fg, donor, donor_exons, recipient);
                    if recip_exons.is_empty() {
                        continue; // recipient cannot legitimately place this chain
                    }
                    if recip_exons.len() < 2 {
                        continue; // single-exon → nothing to verify / transfer
                    }
                    if dbg {
                        n_map_ok += 1;
                    }
                    // GATE: the recipient must have its OWN per-junction evidence.
                    if !verify_recipient_support(fg, recipient, &recip_exons, recip_secs, k) {
                        continue;
                    }
                    if dbg {
                        n_verify_ok += 1;
                    }
                    // Conservative simplification: the isoform is RECOVERED (verified
                    // to >= k support), not flow-quantified — record the verified
                    // support floor as its flow rather than a fabricated abundance.
                    // NOTE: this is a RECOVERY FLOOR, not an abundance estimate. It is
                    // NOT comparable to native `flow` values (real per-chain support);
                    // downstream code must not rank/filter across native vs transferred
                    // paths by `flow` as if they were on the same scale.
                    let flow = k as f64;
                    // G6 (min_flow): respect the base decompose flow floor.
                    if flow < min_flow {
                        continue;
                    }
                    let jc = junctions_of(&recip_exons);
                    out.push(FamilyPath {
                        exons: recip_exons,
                        flow,
                        copy_id: recipient,
                        junction_chain: jc,
                        source: IsoformSource::Transferred { donor_copy: donor },
                    });
                }
            }
        }
        if dbg {
            eprintln!(
                "[layer2-dbg] Part B funnel fam={}: donor_chains={} chains_with_unmatched_donor_exon={} map_attempts={} map_ok={} verify_ok(emitted)={}",
                fg.family_id, n_chains, n_chain_unmatched, n_map_attempt, n_map_ok, n_verify_ok
            );
        }
    }

    // --- G8 dedup by (copy_id, junction_chain), keeping the FIRST occurrence -------
    // First-wins preserves base > Part A > Part B priority (a transferred isoform
    // identical to a recipient's native isoform is dropped in favor of the native).
    let mut seen: DetHashSet<(usize, Vec<(u64, u64)>)> = DetHashSet::default();
    out.retain(|p| seen.insert((p.copy_id, p.junction_chain.clone())));

    // --- Final deterministic total order: (exons, copy_id, source_rank) ------------
    // After the G8 dedup, no two surviving paths share (copy_id, junction_chain), and
    // since junction_chain is derived from exons no two share (exons, copy_id) either —
    // so (exons, copy_id) already totally orders the output and the source_rank tiebreak
    // never actually decides ordering. It is a total-order safety net only. (With every
    // path Native this also reduces to decompose's own (exons, copy_id) ordering, keeping
    // the both-flags-off regression anchor byte-identical.)
    let source_rank = |s: &IsoformSource| -> (usize, usize) {
        match s {
            IsoformSource::Native => (0, 0),
            IsoformSource::Transferred { donor_copy } => (1, *donor_copy),
        }
    };
    out.sort_by(|a, b| {
        a.exons
            .cmp(&b.exons)
            .then(a.copy_id.cmp(&b.copy_id))
            .then(source_rank(&a.source).cmp(&source_rank(&b.source)))
    });

    out
}

// =============================================================================
// All-secondary new-copy emission (C5 / M5.2, PURE — not yet wired).
//
// Layer 1 drops a cluster of secondaries when NO primary bundle exists at its
// region (an unexpressed paralog copy whose reads all placed primarily on the
// well-expressed sibling). Task A's `all_secondary_regions` clusters them; this
// emitter turns each cluster into candidate NEW-COPY transcripts. As everywhere
// in Layer 2, no coordinate is invented — exons come from `build_exons_from_chain`
// (crisp junctions + median terminal ends).
// =============================================================================

/// Emit candidate NEW-COPY transcripts from all-secondary regions (clusters of
/// secondaries Layer 1 dropped because no primary bundle exists there). Each
/// cluster may yield several isoforms: we group its alignments by intron chain,
/// and emit one transcript per distinct chain supported by >= `min_reads`
/// DISTINCT reads (read_name_hash — never raw alignment count, since one molecule
/// can have multiple secondary/supplementary placements in the cluster) AND with
/// >= 1 intron (>= 2 exons). Single-exon clusters are skipped (unvalidatable by
/// chain, FP-prone). Coordinates come from `build_exons_from_chain` (median
/// terminal ends + crisp junctions) — never invented. Strand = majority strand of
/// the chain's supporters (tie -> '+'). Gated: returns empty unless `new_copies`.
/// These are tagged `NovelLocusFromScan` and flow through the additive
/// union-by-chain stage downstream (they only ADD chains absent from VG output).
fn candidate_new_copy_transcripts(
    regions: &[Vec<&crate::vg_family::secondary_index::SecondaryAlignment>],
    new_copies: bool,
    min_reads: usize,
) -> Vec<Transcript> {
    // 1. The gate. Default-off path returns nothing.
    if !new_copies {
        return Vec::new();
    }

    let mut out: Vec<Transcript> = Vec::new();
    for cluster in regions {
        // 2. Group this cluster's alignments by their full intron-chain identity.
        let mut groups: DetHashMap<
            Vec<(u64, u64)>,
            Vec<&crate::vg_family::secondary_index::SecondaryAlignment>,
        > = DetHashMap::default();
        for s in cluster {
            groups.entry(s.introns.clone()).or_default().push(s);
        }

        // Collect surviving (distinct_reads, chain, supporters) for deterministic
        // ordering before we touch the output (never iterate a map into output).
        let mut survivors: Vec<(usize, Vec<(u64, u64)>, Vec<&crate::vg_family::secondary_index::SecondaryAlignment>)> =
            Vec::new();
        for (chain, supporters) in &groups {
            // Single-exon (no junction) → unvalidatable by chain; skip.
            if chain.is_empty() {
                continue;
            }
            // DISTINCT-read support — count read_name_hash, NOT alignment count: one
            // molecule can place secondary+supplementary in the same cluster.
            let mut reads: DetHashSet<u64> = DetHashSet::default();
            for s in supporters {
                reads.insert(s.read_name_hash);
            }
            let distinct = reads.len();
            if distinct < min_reads {
                continue;
            }
            survivors.push((distinct, chain.clone(), supporters.clone()));
        }
        // Deterministic per-cluster order: distinct DESC, chain.len() DESC, chain ASC
        // (mirrors enumerate_secondary_chains).
        survivors.sort_by(|a, b| {
            b.0.cmp(&a.0)
                .then(b.1.len().cmp(&a.1.len()))
                .then(a.1.cmp(&b.1))
        });

        for (distinct, chain, supporters) in survivors {
            let exons = build_exons_from_chain(&supporters, &chain);
            if exons.len() < 2 {
                continue; // defensive: a junction chain must yield >= 2 exons.
            }
            // Majority strand of the chain's supporters; tie / all-equal -> '+'.
            let plus = supporters.iter().filter(|s| s.strand == '+').count();
            let minus = supporters.iter().filter(|s| s.strand == '-').count();
            let strand = if minus > plus { '-' } else { '+' };

            let mut t = Transcript::default();
            // All alignments in a cluster share a chromosome (clusters are per-chrom).
            t.chrom = supporters[0].chrom.clone();
            t.strand = strand;
            t.exons = exons;
            t.coverage = distinct as f64;
            t.longcov = distinct as f64;
            t.is_longread = true;
            t.synthetic = true;
            t.rescue_class = Some(crate::vg_family::diagnostic::RescueClass::NovelLocusFromScan);
            // Standalone scan candidates — not family members; leave vg_* default None.
            out.push(t);
        }
    }

    // 3. Final deterministic output order: (chrom, exons) so emission never
    //    depends on map iteration order across clusters.
    out.sort_by(|a, b| a.chrom.cmp(&b.chrom).then(a.exons.cmp(&b.exons)));
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
            introns: introns.to_vec(), nm: 0, strand: '+', is_supplementary: false, locus: Some(locus),
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
            ..Default::default()
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
            strand: '+',
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
    // T1 — enumerate_secondary_chains
    // -------------------------------------------------------------------------

    #[test]
    fn enumerate_returns_both_chains_above_k_sorted() {
        // Chain X: 3 supporters (2 introns). Chain Y: 2 supporters (1 intron).
        // Both >= K=2 → both returned. Sorted by support DESC → X first.
        let s = [
            sec_owned(100, 560, &[(160, 300), (360, 500)]),
            sec_owned(100, 560, &[(160, 300), (360, 500)]),
            sec_owned(100, 560, &[(160, 300), (360, 500)]),
            sec_owned(100, 560, &[(160, 500)]),
            sec_owned(100, 560, &[(160, 500)]),
        ];
        let refs: Vec<&_> = s.iter().collect();
        let got = enumerate_secondary_chains(&refs, 2);
        assert_eq!(got.len(), 2, "both >=K chains returned: {got:?}");
        // X first (3 support), then Y (2 support).
        assert_eq!(got[0].1, vec![(160, 300), (360, 500)], "X (mode) first");
        assert_eq!(got[0].2, 3);
        assert_eq!(got[0].0, vec![(100, 160), (300, 360), (500, 560)], "X exons");
        assert_eq!(got[1].1, vec![(160, 500)], "Y second");
        assert_eq!(got[1].2, 2);
        assert_eq!(got[1].0, vec![(100, 160), (500, 560)], "Y exons");
    }

    #[test]
    fn enumerate_filters_chain_below_k() {
        // Chain X: 2 supporters; chain Y: 1 supporter (< K=2) → only X returned.
        let s = [
            sec_owned(100, 560, &[(160, 300), (360, 500)]),
            sec_owned(100, 560, &[(160, 300), (360, 500)]),
            sec_owned(100, 560, &[(160, 320), (380, 500)]),
        ];
        let refs: Vec<&_> = s.iter().collect();
        let got = enumerate_secondary_chains(&refs, 2);
        assert_eq!(got.len(), 1, "the <K chain is filtered: {got:?}");
        assert_eq!(got[0].1, vec![(160, 300), (360, 500)]);
        assert_eq!(got[0].2, 2);
    }

    #[test]
    fn enumerate_single_chain_matches_consensus_exons() {
        // One distinct chain at >=K support → exactly one result, and its exons
        // match consensus_exon_chain's single mode chain (factor consistency).
        let s = [
            sec_owned(29998, 30658, &[(30060, 30300), (30360, 30600)]),
            sec_owned(30000, 30660, &[(30060, 30300), (30360, 30600)]),
            sec_owned(30003, 30662, &[(30060, 30300), (30360, 30600)]),
        ];
        let refs: Vec<&_> = s.iter().collect();
        let got = enumerate_secondary_chains(&refs, 2);
        assert_eq!(got.len(), 1, "one chain returned");
        assert_eq!(got[0].0, consensus_exon_chain(&refs), "exons match consensus");
        assert_eq!(got[0].0, vec![(30000, 30060), (30300, 30360), (30600, 30660)]);
    }

    #[test]
    fn enumerate_empty_input_returns_empty() {
        assert!(enumerate_secondary_chains(&[], 2).is_empty());
    }

    #[test]
    fn enumerate_k_default_two_filters_singletons() {
        // K=2 default behavior: a lone read forms a singleton group → excluded.
        let s = [sec_owned(100, 560, &[(160, 300), (360, 500)])];
        let refs: Vec<&_> = s.iter().collect();
        assert!(
            enumerate_secondary_chains(&refs, 2).is_empty(),
            "single read (<K=2) yields no chain"
        );
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
    fn positional_no_anchor_partial_consensus_returns_empty() {
        // Starved copy (1) has zero per_copy_spans entries anywhere AND the consensus
        // does NOT span the whole backbone (2 exons vs 4 nodes) → frame unpinnable.
        // (The exact-length, whole-backbone case is covered by the anchor-free test;
        // a length mismatch with no anchor stays ambiguous → empty, never guess.)
        let fg = mk_backbone(
            &[(100, 160), (300, 360), (500, 560), (700, 760)],
            &[
                vec![(0, (100, 160))],
                vec![(0, (300, 360))],
                vec![(0, (500, 560))],
                vec![(0, (700, 760))],
            ],
        );
        let consensus = [(30000, 30060), (30300, 30360)];
        assert!(
            positional_mapping(&fg, 1, &consensus).is_empty(),
            "no native anchor + partial-span consensus → frame unpinnable → empty (never guess)"
        );
    }

    #[test]
    fn positional_anchor_free_full_span_maps_starved_coords() {
        // Starved copy (1) has ZERO native per_copy_spans entries anywhere, but the
        // consensus spans the WHOLE backbone (len 3 == 3 nodes). Co-linear order pins
        // the frame with no anchor needed: consensus exon i → order[i], at the STARVED
        // copy's OWN 30k coords (NOT the well copy's backbone coords).
        let fg = mk_backbone(
            &[(100, 160), (300, 360), (500, 560)],
            &[vec![(0, (100, 160))], vec![(0, (300, 360))], vec![(0, (500, 560))]],
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
            "anchor-free full-span map yields the STARVED copy's own coords, not the well copy's"
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

    // -------------------------------------------------------------------------
    // T2 — map_isoform_across_copies / verify_recipient_support
    // -------------------------------------------------------------------------

    #[test]
    fn map_isoform_returns_recipient_coords_not_donor() {
        // 3-node backbone. Donor A (copy 0) has its own coords at every node;
        // recipient B (copy 1) has its OWN distinct coords (the 30k locus) at
        // every node. A's SKIP isoform uses nodes 0 and 2 only. The mapped
        // recipient chain must be B's coords for nodes 0 and 2 — never A's.
        let fg = mk_backbone(
            &[(100, 160), (300, 360), (500, 560)],
            &[
                vec![(0, (100, 160)), (1, (30000, 30060))],
                vec![(0, (300, 360)), (1, (30300, 30360))],
                vec![(0, (500, 560)), (1, (30600, 30660))],
            ],
        );
        // Donor A's skip isoform: exon 0 (node 0) then exon 2 (node 2).
        let donor_exons = [(100, 160), (500, 560)];
        let got = map_isoform_across_copies(&fg, 0, &donor_exons, 1);
        assert_eq!(
            got,
            vec![(30000, 30060), (30600, 30660)],
            "recipient B's own coords for the shared nodes, NOT donor A's"
        );
        // Explicitly assert none of A's coordinates leaked.
        assert!(
            !got.iter().any(|&sp| donor_exons.contains(&sp)),
            "no donor coordinate appears in the transferred chain"
        );
    }

    #[test]
    fn map_isoform_aborts_on_copy_a_unique_node() {
        // Node 1 is copy-A-UNIQUE (only copy 0 present, single entry → not a
        // >=2-agree co-located borrow). A's full chain touches node 1, which B
        // cannot legitimately claim → abort (empty).
        let fg = mk_backbone(
            &[(100, 160), (300, 360), (500, 560)],
            &[
                vec![(0, (100, 160)), (1, (30000, 30060))],
                vec![(0, (300, 360))], // copy-0-unique node → B cannot claim it
                vec![(0, (500, 560)), (1, (30600, 30660))],
            ],
        );
        let donor_exons = [(100, 160), (300, 360), (500, 560)];
        assert!(
            map_isoform_across_copies(&fg, 0, &donor_exons, 1).is_empty(),
            "recipient cannot claim a copy-A-unique node → abort"
        );
    }

    #[test]
    fn map_isoform_overlap_not_exact() {
        // Co-located paralog: donor (copy 0) and recipient (copy 1) share the SAME
        // backbone nodes at the SAME coordinates. The donor isoform's exons are
        // secondary-consensus coords (median terminals + crisp junctions) that
        // OVERLAP — but never exactly equal — the Layer-1-assembled node spans.
        // Exact-equality lookup finds no node and returns empty (the real-chrY bug);
        // overlap matching maps each donor exon to the node it overlaps and emits the
        // RECIPIENT's own coordinates.
        let fg = mk_backbone(
            &[(100, 160), (300, 360), (500, 560)],
            &[
                vec![(0, (100, 160)), (1, (100, 160))],
                vec![(0, (300, 360)), (1, (300, 360))],
                vec![(0, (500, 560)), (1, (500, 560))],
            ],
        );
        // Donor exons overlap each node but equal none (off-by-1 terminals/junctions).
        let donor_exons = [(101, 161), (301, 361), (501, 559)];
        let got = map_isoform_across_copies(&fg, 0, &donor_exons, 1);
        assert_eq!(
            got,
            vec![(100, 160), (300, 360), (500, 560)],
            "overlap matching maps donor exons to their nodes and emits recipient coords"
        );
        // Negative: a donor exon overlapping NO backbone node → abort (empty).
        assert!(
            map_isoform_across_copies(&fg, 0, &[(5000, 5060)], 1).is_empty(),
            "donor exon overlapping no backbone node → abort"
        );
    }

    /// Build a 3-node backbone graph that ALSO carries junction edges, so
    /// verify_recipient_support can be exercised on its fg.edges path.
    fn mk_backbone_with_edges(
        well: &[(u64, u64)],
        per_copy: &[Vec<(usize, (u64, u64))>],
        edges: &[(usize, usize)],
    ) -> FamilyGraph {
        use crate::vg_family::family_graph::{JunctionEdge, NodeIdx};
        let mut fg = mk_backbone(well, per_copy);
        fg.edges = edges
            .iter()
            .map(|&(f, t)| JunctionEdge {
                from: NodeIdx(f),
                to: NodeIdx(t),
                family_support: 5,
                strand: '+',
            })
            .collect();
        fg
    }

    #[test]
    fn verify_true_when_recipient_secondaries_cover_junction() {
        // Recipient B (copy 1) at the 30k locus. Transferred chain has one
        // junction (30060 -> 30300). B has 2 of its OWN secondaries carrying it
        // → >= K=2 → admissible.
        let fg = mk_backbone_with_edges(
            &[(100, 160), (300, 360)],
            &[
                vec![(0, (100, 160)), (1, (30000, 30060))],
                vec![(0, (300, 360)), (1, (30300, 30360))],
            ],
            &[],
        );
        let recipient_exons = [(30000, 30060), (30300, 30360)];
        let s = [
            sec_owned(30000, 30360, &[(30060, 30300)]),
            sec_owned(30000, 30360, &[(30060, 30300)]),
        ];
        let refs: Vec<&_> = s.iter().collect();
        assert!(
            verify_recipient_support(&fg, 1, &recipient_exons, &refs, 2),
            "B's own secondaries cover the transferred junction → admissible"
        );
    }

    #[test]
    fn verify_false_when_recipient_has_no_evidence() {
        // Same transferred chain, but B has ZERO supporting secondaries and no
        // fg.edge between its member nodes → not admissible.
        let fg = mk_backbone_with_edges(
            &[(100, 160), (300, 360)],
            &[
                vec![(0, (100, 160)), (1, (30000, 30060))],
                vec![(0, (300, 360)), (1, (30300, 30360))],
            ],
            &[],
        );
        let recipient_exons = [(30000, 30060), (30300, 30360)];
        assert!(
            !verify_recipient_support(&fg, 1, &recipient_exons, &[], 2),
            "recipient has 0 evidence → not admissible (recipient disposes)"
        );
    }

    #[test]
    fn verify_false_when_recipient_support_below_k() {
        // B has exactly ONE supporting secondary (< K=2) → not admissible.
        let fg = mk_backbone_with_edges(
            &[(100, 160), (300, 360)],
            &[
                vec![(0, (100, 160)), (1, (30000, 30060))],
                vec![(0, (300, 360)), (1, (30300, 30360))],
            ],
            &[],
        );
        let recipient_exons = [(30000, 30060), (30300, 30360)];
        let s = [sec_owned(30000, 30360, &[(30060, 30300)])];
        let refs: Vec<&_> = s.iter().collect();
        assert!(
            !verify_recipient_support(&fg, 1, &recipient_exons, &refs, 2),
            "single supporter (<K=2) → not admissible"
        );
    }

    #[test]
    fn verify_true_via_fg_edge_between_recipient_native_nodes() {
        // Recipient B (copy 1) at the 30k locus, NO secondaries — but an fg.edge
        // (node 0 -> node 1) between B-native nodes carries the junction. The edge
        // coordinate MUST be derived from B's OWN per_copy_spans (30060 -> 30300),
        // NOT the node union span (which, with copy-0 at 100-360, would be a
        // nonsensical 30060 -> 300). This exercises the fg.edges fallback path and
        // is a regression guard for the union-span bug.
        let fg = mk_backbone_with_edges(
            &[(100, 160), (300, 360)],
            &[
                vec![(0, (100, 160)), (1, (30000, 30060))],
                vec![(0, (300, 360)), (1, (30300, 30360))],
            ],
            &[(0, 1)],
        );
        let recipient_exons = [(30000, 30060), (30300, 30360)];
        assert!(
            verify_recipient_support(&fg, 1, &recipient_exons, &[], 2),
            "fg.edge between B-native nodes supports the junction at B's own coords"
        );
    }

    // -------------------------------------------------------------------------
    // T3 — emit_family_isoforms (orchestrator)
    // -------------------------------------------------------------------------

    /// Empty secondaries-by-copy map (the both-flags-off regression-anchor input).
    fn empty_secs_map<'a>(
    ) -> DetHashMap<usize, Vec<&'a crate::vg_family::secondary_index::SecondaryAlignment>> {
        DetHashMap::default()
    }

    #[test]
    fn emit_both_off_equals_decompose() {
        // REGRESSION ANCHOR: both channels gated off → emit == decompose byte-for-byte
        // (including the new junction_chain + source=Native fields). Built from the
        // same fixture as decompose_recovers_native_copy_defers_crosslocus_borrow.
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
        let empty = empty_secs_map();
        let emitted = emit_family_isoforms(&fg, &am, &empty, 1.0, 2, false, false);
        let base = decompose_family_paths(&fg, &am, 1.0);
        assert_eq!(emitted, base, "both flags off → exactly decompose's output");
    }

    #[test]
    fn emit_part_a_discovers_alt_splice_isoform() {
        // ONE copy whose secondaries carry TWO distinct intron chains, each with 3
        // supporters: a full chain [(160,300),(360,500)] and an exon-skipping chain
        // [(160,500)]. With enable_multi_isoform=true,k=2 emit returns BOTH chains
        // for that copy (both Native). With the flag off, the alt (skip) is absent.
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
        let am = GraphAmendment::default();
        // copy 0's own secondaries: 3 full-chain + 3 skip-chain.
        let secs: Vec<_> = vec![
            sec_aln(1, 100, 560, &[(160, 300), (360, 500)], 0),
            sec_aln(2, 100, 560, &[(160, 300), (360, 500)], 0),
            sec_aln(3, 100, 560, &[(160, 300), (360, 500)], 0),
            sec_aln(4, 100, 560, &[(160, 500)], 0),
            sec_aln(5, 100, 560, &[(160, 500)], 0),
            sec_aln(6, 100, 560, &[(160, 500)], 0),
        ];
        let refs: Vec<&_> = secs.iter().collect();
        let mut by_copy: DetHashMap<usize, Vec<&crate::vg_family::secondary_index::SecondaryAlignment>> =
            DetHashMap::default();
        by_copy.insert(0, refs);

        let on = emit_family_isoforms(&fg, &am, &by_copy, 1.0, 2, true, false);
        let full = vec![(160, 300), (360, 500)];
        let skip = vec![(160, 500)];
        assert!(
            on.iter().any(|p| p.copy_id == 0 && p.junction_chain == full && p.source == IsoformSource::Native),
            "full alt-splice chain emitted as Native: {on:?}"
        );
        assert!(
            on.iter().any(|p| p.copy_id == 0 && p.junction_chain == skip && p.source == IsoformSource::Native),
            "exon-skipping chain emitted as Native: {on:?}"
        );

        let off = emit_family_isoforms(&fg, &am, &by_copy, 1.0, 2, false, false);
        assert!(
            !off.iter().any(|p| p.junction_chain == skip),
            "with multi-isoform OFF the alt (skip) chain is absent: {off:?}"
        );
    }

    #[test]
    fn emit_part_a_respects_support_floor() {
        // This case covers the K-FILTER inside enumerate_secondary_chains (NOT G7): a
        // lone read (support 1 < k=2) never even survives enumeration, so it is gone
        // before Part A's G7 isofrac floor runs. G7 is exercised independently by
        // emit_part_a_g7_floor_drops_minor_isoform below (a chain that CLEARS k yet
        // falls under the isofrac floor). Here: a high-support dominant + a sub-k chain.
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
        let am = GraphAmendment::default();
        // 4 supporters of the dominant full chain; a 1-supporter (sub-k) skip chain.
        // The skip chain is dropped by enumerate's support<k filter — it never reaches
        // the G7 floor. (G7's distinct behavior is proven in the sibling test.)
        let mut secs: Vec<crate::vg_family::secondary_index::SecondaryAlignment> = Vec::new();
        for h in 0..4u64 {
            secs.push(sec_aln(h, 100, 560, &[(160, 300), (360, 500)], 0));
        }
        // a lone read on a different chain (skip) — support 1, below k=2.
        secs.push(sec_aln(100, 100, 560, &[(160, 500)], 0));
        let refs: Vec<&_> = secs.iter().collect();
        let mut by_copy: DetHashMap<usize, Vec<&crate::vg_family::secondary_index::SecondaryAlignment>> =
            DetHashMap::default();
        by_copy.insert(0, refs);

        let on = emit_family_isoforms(&fg, &am, &by_copy, 1.0, 2, true, false);
        let dominant = vec![(160, 300), (360, 500)];
        let skip = vec![(160, 500)];
        assert!(
            on.iter().any(|p| p.copy_id == 0 && p.junction_chain == dominant),
            "dominant chain (support 4 >= floor) emitted: {on:?}"
        );
        assert!(
            !on.iter().any(|p| p.junction_chain == skip),
            "sub-k minor chain NOT emitted (dropped by enumerate's support<k filter): {on:?}"
        );
    }

    #[test]
    fn emit_part_a_g7_floor_drops_minor_isoform() {
        // G7 isofrac-floor coverage that is INDEPENDENT of enumerate's support<k filter.
        // One copy carries a DOMINANT chain with ~300 supporters and a MINOR alt chain
        // with support == 2 (== k, so it CLEARS the k-filter and survives enumeration).
        // With copy_max=300 the G7 floor is max(2.0, 0.01*300)=3.0, so the support-2 alt
        // is dropped by G7 — NOT by the k-filter. This is the case the existing
        // emit_part_a_respects_support_floor test does NOT reach: deleting the G7 skip
        // line in Part A makes ONLY this test fail (the alt would then be emitted).
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
        let am = GraphAmendment::default();
        let dominant_chain = [(160, 300), (360, 500)];
        let minor_chain = [(160, 500)];
        // ~300 supporters of the dominant chain (copy_max=300 → G7 floor = 3.0) ...
        let mut secs: Vec<crate::vg_family::secondary_index::SecondaryAlignment> = Vec::new();
        for i in 0..300u64 {
            secs.push(sec_aln(i, 100, 560, &dominant_chain, 0));
        }
        // ... plus EXACTLY 2 supporters of the minor alt chain: support==k=2 clears the
        // k-filter, but 2 < 3.0 floor → G7 must drop it.
        secs.push(sec_aln(1000, 100, 560, &minor_chain, 0));
        secs.push(sec_aln(1001, 100, 560, &minor_chain, 0));
        let refs: Vec<&_> = secs.iter().collect();
        let mut by_copy: DetHashMap<usize, Vec<&crate::vg_family::secondary_index::SecondaryAlignment>> =
            DetHashMap::default();
        by_copy.insert(0, refs);

        let on = emit_family_isoforms(&fg, &am, &by_copy, 1.0, 2, true, false);
        let dominant = dominant_chain.to_vec();
        let minor = minor_chain.to_vec();
        assert!(
            on.iter().any(|p| p.copy_id == 0 && p.junction_chain == dominant),
            "dominant chain (support 300 >= floor) emitted: {on:?}"
        );
        assert!(
            !on.iter().any(|p| p.junction_chain == minor),
            "minor chain (support 2 >= k but < 3.0 isofrac floor) dropped by G7: {on:?}"
        );
    }

    #[test]
    fn emit_part_b_transfers_with_recipient_evidence() {
        // Two co-located 3-exon copies (per_copy_spans line up on the backbone).
        // Donor copy 0 has a full isoform (its enumerated chain). Recipient copy 1
        // has per_copy_spans on the same nodes AND its own secondaries carrying the
        // SAME junctions → transfer admitted at copy 1's OWN coords.
        let fg = mk_backbone_with_edges(
            &[(100, 160), (300, 360), (500, 560)],
            &[
                vec![(0, (100, 160)), (1, (30000, 30060))],
                vec![(0, (300, 360)), (1, (30300, 30360))],
                vec![(0, (500, 560)), (1, (30600, 30660))],
            ],
            &[],
        );
        let am = GraphAmendment::default();
        // Donor 0's own secondaries → full 3-exon chain at copy 0's coords.
        let donor_secs: Vec<_> = vec![
            sec_aln(1, 100, 560, &[(160, 300), (360, 500)], 0),
            sec_aln(2, 100, 560, &[(160, 300), (360, 500)], 0),
        ];
        // Recipient 1's own secondaries → the SAME junctions at copy 1's coords.
        let recip_secs: Vec<_> = vec![
            sec_aln(3, 30000, 30660, &[(30060, 30300), (30360, 30600)], 1),
            sec_aln(4, 30000, 30660, &[(30060, 30300), (30360, 30600)], 1),
        ];
        let dref: Vec<&_> = donor_secs.iter().collect();
        let rref: Vec<&_> = recip_secs.iter().collect();
        let mut by_copy: DetHashMap<usize, Vec<&crate::vg_family::secondary_index::SecondaryAlignment>> =
            DetHashMap::default();
        by_copy.insert(0, dref);
        by_copy.insert(1, rref);

        let on = emit_family_isoforms(&fg, &am, &by_copy, 1.0, 2, false, true);
        let transferred: Vec<&FamilyPath> = on
            .iter()
            .filter(|p| p.copy_id == 1 && p.source == IsoformSource::Transferred { donor_copy: 0 })
            .collect();
        assert!(!transferred.is_empty(), "a transferred isoform for copy 1 is emitted: {on:?}");
        // It must be at copy 1's OWN coordinates (the 30k locus), never copy 0's.
        assert!(
            transferred.iter().any(|p| p.exons == vec![(30000, 30060), (30300, 30360), (30600, 30660)]),
            "transferred chain is at recipient copy 1's own coords: {transferred:?}"
        );
        assert!(
            transferred.iter().all(|p| p.exons.iter().all(|&(s, _)| s >= 30000)),
            "no donor (copy 0, <1000) coordinate leaked into the transferred chain: {transferred:?}"
        );
    }

    #[test]
    fn emit_part_b_rejects_without_recipient_evidence() {
        // Same donor isoform, but recipient copy 1 has NO secondaries carrying the
        // junctions and no fg.edge support → verify_recipient_support fails → no
        // Transferred path. Also: with enable_transfer=false, none appears.
        let fg = mk_backbone_with_edges(
            &[(100, 160), (300, 360), (500, 560)],
            &[
                vec![(0, (100, 160)), (1, (30000, 30060))],
                vec![(0, (300, 360)), (1, (30300, 30360))],
                vec![(0, (500, 560)), (1, (30600, 30660))],
            ],
            &[],
        );
        let am = GraphAmendment::default();
        let donor_secs: Vec<_> = vec![
            sec_aln(1, 100, 560, &[(160, 300), (360, 500)], 0),
            sec_aln(2, 100, 560, &[(160, 300), (360, 500)], 0),
        ];
        let dref: Vec<&_> = donor_secs.iter().collect();
        let mut by_copy: DetHashMap<usize, Vec<&crate::vg_family::secondary_index::SecondaryAlignment>> =
            DetHashMap::default();
        by_copy.insert(0, dref);
        // copy 1 deliberately has NO secondaries entry.

        let on = emit_family_isoforms(&fg, &am, &by_copy, 1.0, 2, false, true);
        assert!(
            !on.iter().any(|p| p.copy_id == 1 && matches!(p.source, IsoformSource::Transferred { .. })),
            "no recipient evidence → no Transferred path for copy 1: {on:?}"
        );
        let off = emit_family_isoforms(&fg, &am, &by_copy, 1.0, 2, false, false);
        assert!(
            !off.iter().any(|p| matches!(p.source, IsoformSource::Transferred { .. })),
            "with transfer OFF no Transferred path appears: {off:?}"
        );
    }

    #[test]
    fn emit_dedup_prefers_native_over_transfer() {
        // copy 1 has its OWN full isoform (Part A Native) AND copy 0's identical
        // isoform would transfer onto copy 1 with the SAME (copy_id, junction_chain).
        // Dedup (base/Part-A first) keeps the Native one; the Transferred dup drops.
        let fg = mk_backbone_with_edges(
            &[(100, 160), (300, 360), (500, 560)],
            &[
                vec![(0, (100, 160)), (1, (30000, 30060))],
                vec![(0, (300, 360)), (1, (30300, 30360))],
                vec![(0, (500, 560)), (1, (30600, 30660))],
            ],
            &[],
        );
        let am = GraphAmendment::default();
        let donor_secs: Vec<_> = vec![
            sec_aln(1, 100, 560, &[(160, 300), (360, 500)], 0),
            sec_aln(2, 100, 560, &[(160, 300), (360, 500)], 0),
        ];
        // Recipient 1 has its OWN secondaries for the same chain → Part A emits it
        // as Native; Part B would also propose it (identical key) → dedup to one.
        let recip_secs: Vec<_> = vec![
            sec_aln(3, 30000, 30660, &[(30060, 30300), (30360, 30600)], 1),
            sec_aln(4, 30000, 30660, &[(30060, 30300), (30360, 30600)], 1),
        ];
        let dref: Vec<&_> = donor_secs.iter().collect();
        let rref: Vec<&_> = recip_secs.iter().collect();
        let mut by_copy: DetHashMap<usize, Vec<&crate::vg_family::secondary_index::SecondaryAlignment>> =
            DetHashMap::default();
        by_copy.insert(0, dref);
        by_copy.insert(1, rref);

        let on = emit_family_isoforms(&fg, &am, &by_copy, 1.0, 2, true, true);
        let recip_chain = vec![(30060, 30300), (30360, 30600)];
        let matching: Vec<&FamilyPath> = on
            .iter()
            .filter(|p| p.copy_id == 1 && p.junction_chain == recip_chain)
            .collect();
        assert_eq!(matching.len(), 1, "exactly one path survives dedup for that key: {matching:?}");
        assert_eq!(
            matching[0].source,
            IsoformSource::Native,
            "the survivor is the Native isoform (base/Part-A wins over transfer)"
        );
    }

    // -------------------------------------------------------------------------
    // M5.2 candidate_new_copy_transcripts — spliced consensus emission for
    // all-secondary regions (PURE; not yet wired into run_layer2).
    // -------------------------------------------------------------------------

    // `sec_aln` defaults strand to '+'; this variant lets a test pick the strand
    // (needed for the strand-consensus test).
    fn sec_aln_s(
        h: u64, start: u64, end: u64, introns: &[(u64, u64)], locus: usize, strand: char,
    ) -> crate::vg_family::secondary_index::SecondaryAlignment {
        let mut s = sec_aln(h, start, end, introns, locus);
        s.strand = strand;
        s
    }

    #[test]
    fn new_copy_gated_off_returns_empty() {
        // Non-empty regions but the gate is off → nothing emitted (proves the
        // default-off path).
        let regions_owned = vec![vec![
            sec_aln(1, 1000, 1600, &[(1200, 1400)], 0),
            sec_aln(2, 1000, 1600, &[(1200, 1400)], 0),
            sec_aln(3, 1000, 1600, &[(1200, 1400)], 0),
        ]];
        let regions: Vec<Vec<&_>> = regions_owned.iter().map(|c| c.iter().collect()).collect();
        let out = candidate_new_copy_transcripts(&regions, false, 3);
        assert!(out.is_empty(), "gate off → empty: {out:?}");
    }

    #[test]
    fn new_copy_emits_multiexon_consensus() {
        // One cluster, 3 alignments with DISTINCT reads sharing intron chain
        // [(1200,1400)] → exons (1000,1200),(1400,1600). Expect one transcript.
        let regions_owned = vec![vec![
            sec_aln(1, 1000, 1600, &[(1200, 1400)], 0),
            sec_aln(2, 1000, 1600, &[(1200, 1400)], 0),
            sec_aln(3, 1000, 1600, &[(1200, 1400)], 0),
        ]];
        let regions: Vec<Vec<&_>> = regions_owned.iter().map(|c| c.iter().collect()).collect();
        let out = candidate_new_copy_transcripts(&regions, true, 3);
        assert_eq!(out.len(), 1, "exactly one transcript: {out:?}");
        let t = &out[0];
        assert_eq!(t.exons, vec![(1000, 1200), (1400, 1600)], "consensus exons");
        assert_eq!(t.coverage, 3.0, "coverage == distinct read count");
        assert_eq!(t.longcov, 3.0, "longcov == distinct read count");
        assert!(t.synthetic, "synthetic candidate");
        assert!(t.is_longread, "is_longread");
        assert_eq!(
            t.rescue_class,
            Some(crate::vg_family::diagnostic::RescueClass::NovelLocusFromScan),
            "tagged NovelLocusFromScan"
        );
        assert_eq!(t.vg_family_id, None, "standalone scan candidate, no family");
        assert_eq!(t.vg_copy_id, None, "standalone scan candidate, no copy id");
    }

    #[test]
    fn new_copy_distinct_read_gate() {
        // Cluster A: 3 ALIGNMENTS but only 2 DISTINCT reads (read 1 appears twice)
        // → must NOT be emitted with min_reads=3 (proves distinct-read counting,
        // NOT alignment counting). Cluster B: 3 DISTINCT reads → emitted. This
        // brackets the gate: same min_reads, only the distinct count differs.
        let region_a = vec![
            sec_aln(1, 1000, 1600, &[(1200, 1400)], 0),
            sec_aln(1, 1000, 1600, &[(1200, 1400)], 0), // same read again
            sec_aln(2, 1000, 1600, &[(1200, 1400)], 0),
        ];
        let region_b = vec![
            sec_aln(10, 5000, 5600, &[(5200, 5400)], 1),
            sec_aln(11, 5000, 5600, &[(5200, 5400)], 1),
            sec_aln(12, 5000, 5600, &[(5200, 5400)], 1),
        ];
        let regions: Vec<Vec<&_>> = vec![region_a.iter().collect(), region_b.iter().collect()];
        let out = candidate_new_copy_transcripts(&regions, true, 3);
        assert_eq!(out.len(), 1, "only the 3-distinct-read cluster emits: {out:?}");
        assert_eq!(
            out[0].exons,
            vec![(5000, 5200), (5400, 5600)],
            "the emitted transcript is cluster B (3 distinct reads)"
        );
    }

    #[test]
    fn new_copy_skips_single_exon() {
        // Single-exon alignments (empty intron chain) are unvalidatable by chain
        // and FP-prone → skipped entirely, even with plenty of distinct reads.
        let region = vec![
            sec_aln(1, 1000, 1600, &[], 0),
            sec_aln(2, 1000, 1600, &[], 0),
            sec_aln(3, 1000, 1600, &[], 0),
        ];
        let regions: Vec<Vec<&_>> = vec![region.iter().collect()];
        let out = candidate_new_copy_transcripts(&regions, true, 3);
        assert!(out.is_empty(), "single-exon cluster skipped: {out:?}");
    }

    #[test]
    fn new_copy_strand_consensus() {
        // Majority of supporters are '-' → emitted strand is '-'.
        let region = vec![
            sec_aln_s(1, 1000, 1600, &[(1200, 1400)], 0, '-'),
            sec_aln_s(2, 1000, 1600, &[(1200, 1400)], 0, '-'),
            sec_aln_s(3, 1000, 1600, &[(1200, 1400)], 0, '+'),
        ];
        let regions: Vec<Vec<&_>> = vec![region.iter().collect()];
        let out = candidate_new_copy_transcripts(&regions, true, 3);
        assert_eq!(out.len(), 1, "one transcript: {out:?}");
        assert_eq!(out[0].strand, '-', "majority strand '-'");
    }

    #[test]
    fn new_copy_multiple_isoforms_one_cluster() {
        // One cluster containing two distinct intron chains, each with >= min_reads
        // distinct reads → two transcripts (distinct exon chains).
        let region = vec![
            // chain X: [(1200,1400)] → exons (1000,1200),(1400,1600)
            sec_aln(1, 1000, 1600, &[(1200, 1400)], 0),
            sec_aln(2, 1000, 1600, &[(1200, 1400)], 0),
            sec_aln(3, 1000, 1600, &[(1200, 1400)], 0),
            // chain Y: [(1200,1500)] → exons (1000,1200),(1500,1700)
            sec_aln(4, 1000, 1700, &[(1200, 1500)], 0),
            sec_aln(5, 1000, 1700, &[(1200, 1500)], 0),
            sec_aln(6, 1000, 1700, &[(1200, 1500)], 0),
        ];
        let regions: Vec<Vec<&_>> = vec![region.iter().collect()];
        let out = candidate_new_copy_transcripts(&regions, true, 3);
        assert_eq!(out.len(), 2, "two distinct isoforms: {out:?}");
        let mut chains: Vec<Vec<(u64, u64)>> = out.iter().map(|t| t.exons.clone()).collect();
        chains.sort();
        assert_eq!(
            chains,
            vec![
                vec![(1000, 1200), (1400, 1600)],
                vec![(1000, 1200), (1500, 1700)],
            ],
            "both chains emitted"
        );
    }

    #[test]
    fn new_copy_deterministic() {
        // The SAME regions built with alignments in two different orders must give
        // identical output (count + each transcript's (chrom, strand, exons)).
        let mk = |order: &[usize]| -> Vec<crate::vg_family::secondary_index::SecondaryAlignment> {
            let pool = vec![
                sec_aln(1, 1000, 1600, &[(1200, 1400)], 0),
                sec_aln(2, 1000, 1600, &[(1200, 1400)], 0),
                sec_aln(3, 1000, 1600, &[(1200, 1400)], 0),
                sec_aln(4, 1000, 1700, &[(1200, 1500)], 0),
                sec_aln(5, 1000, 1700, &[(1200, 1500)], 0),
                sec_aln(6, 1000, 1700, &[(1200, 1500)], 0),
            ];
            order.iter().map(|&i| pool[i].clone()).collect()
        };
        let a = mk(&[0, 1, 2, 3, 4, 5]);
        let b = mk(&[5, 0, 4, 1, 3, 2]);
        let ra: Vec<Vec<&_>> = vec![a.iter().collect()];
        let rb: Vec<Vec<&_>> = vec![b.iter().collect()];
        let oa = candidate_new_copy_transcripts(&ra, true, 3);
        let ob = candidate_new_copy_transcripts(&rb, true, 3);
        assert_eq!(oa.len(), ob.len(), "same transcript count regardless of input order");
        let key = |ts: &[Transcript]| -> Vec<(String, char, Vec<(u64, u64)>)> {
            ts.iter().map(|t| (t.chrom.clone(), t.strand, t.exons.clone())).collect()
        };
        assert_eq!(key(&oa), key(&ob), "byte-identical (chrom, strand, exons) ordering");
    }
}

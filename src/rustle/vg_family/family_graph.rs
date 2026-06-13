//! Family graph: union of per-copy splice graphs across a gene family.
//!
//! Nodes are exon-equivalence classes (groups of exons across copies that
//! we treat as the same exon). Edges are junctions observed in any copy.
//! Copy-specific exons become singleton nodes — that is what makes a
//! family graph different from the per-copy splice graphs it is built
//! from.

use std::ops::Index;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct NodeIdx(pub usize);

pub type CopyId = usize;

#[derive(Debug, Clone)]
pub struct ExonClass {
    pub idx: NodeIdx,
    pub chrom: String,
    pub span: (u64, u64),
    pub strand: char,
    /// (CopyId, raw genome-forward (forward-strand) exon sequence fetched from
    /// the reference FASTA).
    pub per_copy_sequences: Vec<(CopyId, Vec<u8>)>,
    /// (CopyId, original genomic span) for each contributing copy. The node's
    /// `span` field is the union; `per_copy_spans` preserves the exact
    /// per-copy coordinates so a transcript can be reconstructed losslessly.
    pub per_copy_spans: Vec<(CopyId, (u64, u64))>,
    /// True iff exactly one copy contributes — i.e. a bubble branch.
    pub copy_specific: bool,
    /// Per-copy weighted coverage for this exon class, populated after EM by
    /// `annotate_per_copy_exon_coverage`. Entry `(k, cov)` gives the
    /// EM-weighted per-base coverage from copy k's reads. Copies that don't
    /// contribute to this exon class have cov=0.0 as a sentinel entry.
    /// Empty until annotation is called.
    pub per_copy_cov: Vec<(CopyId, f64)>,
}

#[derive(Debug, Clone, Copy)]
pub struct JunctionEdge {
    pub from: NodeIdx,
    pub to: NodeIdx,
    pub family_support: u32,
    pub strand: char,
}

#[derive(Debug, Clone)]
pub struct FamilyGraph {
    pub family_id: usize,
    pub nodes: Vec<ExonClass>,
    pub edges: Vec<JunctionEdge>,
}

impl FamilyGraph {
    pub fn n_nodes(&self) -> usize { self.nodes.len() }
    pub fn n_edges(&self) -> usize { self.edges.len() }

    /// Recover the path of `NodeIdx` that paralog `cid` traverses through this
    /// graph, in genomic-ascending order. Walks `nodes`, keeps only those
    /// where `per_copy_sequences` contains an entry for `cid`. Used as the
    /// E-step of HMM-based EM: the per-paralog path is what each in-family
    /// copy "looks like" through the family graph.
    ///
    /// Returns an empty Vec if `cid` does not contribute to any node.
    pub fn recover_paralog_path(&self, cid: CopyId) -> Vec<NodeIdx> {
        let mut nodes_with_pos: Vec<(u64, NodeIdx)> = self.nodes.iter()
            .filter(|n| n.per_copy_sequences.iter().any(|(c, _)| *c == cid))
            .map(|n| (n.span.0, n.idx))
            .collect();
        nodes_with_pos.sort_by_key(|(s, _)| *s);
        nodes_with_pos.into_iter().map(|(_, idx)| idx).collect()
    }

    /// Set of all CopyIds that contribute to at least one exon class.
    /// Sorted ascending (deterministic for downstream iteration).
    pub fn all_copies(&self) -> Vec<CopyId> {
        let mut s: std::collections::BTreeSet<CopyId> = std::collections::BTreeSet::new();
        for n in &self.nodes {
            for (cid, _) in &n.per_copy_sequences {
                s.insert(*cid);
            }
        }
        s.into_iter().collect()
    }

    /// Pick the "representative" copy whose path traverses the most exon
    /// classes — i.e., the most complete known paralog. Used for Phase 3.1
    /// projection of exon structure onto a candidate locus.
    pub fn representative_copy(&self) -> Option<CopyId> {
        let mut counts: std::collections::BTreeMap<CopyId, usize> = std::collections::BTreeMap::new();
        for n in &self.nodes {
            for (cid, _) in &n.per_copy_sequences {
                *counts.entry(*cid).or_insert(0) += 1;
            }
        }
        counts.into_iter().max_by_key(|(_, c)| *c).map(|(cid, _)| cid)
    }

    /// Genomic-ordered (start, end) spans for the given copy's exons. Returns
    /// empty Vec if the copy doesn't contribute. Used by Phase 3.1 to
    /// project paralog structure onto a novel candidate locus.
    pub fn paralog_exon_spans(&self, cid: CopyId) -> Vec<(u64, u64)> {
        let path = self.recover_paralog_path(cid);
        path.into_iter().filter_map(|nidx| {
            self.nodes[nidx.0].per_copy_spans.iter()
                .find(|(c, _)| *c == cid)
                .map(|(_, sp)| *sp)
        }).collect()
    }
}

impl Index<NodeIdx> for FamilyGraph {
    type Output = ExonClass;
    fn index(&self, idx: NodeIdx) -> &ExonClass { &self.nodes[idx.0] }
}

use crate::types::Bundle;

/// Collect, dedup, and sort exonic intervals from a bundle's reads.
/// Returns half-open (start, end) pairs on the bundle's chromosome.
pub fn extract_copy_exons(bundle: &Bundle) -> Vec<(u64, u64)> {
    // Primary path: reads carry exact aligned-exon spans.
    let mut all: Vec<(u64, u64)> = bundle.reads.iter()
        .flat_map(|r| r.exons.iter().copied())
        .collect();
    if all.is_empty() {
        // Fallback: lightweight bundle clones (e.g., the slim copy passed to
        // discover_novel_copies — pipeline.rs:8137) strip reads to save
        // memory. Reconstruct exons from junction_stats: the regions between
        // consecutive (donor, acceptor) pairs are introns, so exonic regions
        // are the gaps. First/last exon use bundle.{start,end}.
        let mut donors: Vec<u64> = Vec::new();
        let mut acceptors: Vec<u64> = Vec::new();
        for (j, _) in &bundle.junction_stats {
            donors.push(j.donor);
            acceptors.push(j.acceptor);
        }
        donors.sort_unstable();
        acceptors.sort_unstable();
        if donors.is_empty() {
            return vec![(bundle.start, bundle.end)];
        }
        let n = donors.len();
        let mut out = Vec::with_capacity(n + 1);
        out.push((bundle.start, donors[0]));
        for i in 0..n.saturating_sub(1) {
            let s = acceptors[i];
            let e = donors[i + 1];
            if e > s {
                out.push((s, e));
            }
        }
        out.push((acceptors[n - 1], bundle.end));
        out.sort_unstable();
        out.dedup();
        return out;
    }
    all.sort_unstable();
    all.dedup();
    all
}

/// One member of an exon cluster: (copy id, exon index within copy).
pub type ExonRef = (CopyId, usize);

/// Cluster exons across copies by reciprocal-overlap fraction on the same
/// (chrom, strand). `min_recip` is the minimum reciprocal-overlap (e.g. 0.30
/// = 30%) for two exons to join the same cluster. Single-copy exons emerge
/// as singleton clusters.
pub fn cluster_by_position(
    copies: &[(&str, char, Vec<(u64, u64)>)],
    min_recip: f64,
) -> Vec<Vec<ExonRef>> {
    // Flat list of (copy_id, exon_idx, chrom, strand, start, end).
    let mut all: Vec<(CopyId, usize, &str, char, u64, u64)> = Vec::new();
    for (cid, (chrom, strand, exons)) in copies.iter().enumerate() {
        for (ei, &(s, e)) in exons.iter().enumerate() {
            all.push((cid, ei, chrom, *strand, s, e));
        }
    }

    // Union-Find over the flat list.
    let n = all.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        if p[x] == x { x } else { let r = find(p, p[x]); p[x] = r; r }
    }
    fn union(p: &mut [usize], a: usize, b: usize) {
        let ra = find(p, a); let rb = find(p, b);
        if ra != rb { p[ra] = rb; }
    }

    for i in 0..n {
        for j in (i + 1)..n {
            let (_, _, ci, si, ai, bi) = all[i];
            let (_, _, cj, sj, aj, bj) = all[j];
            if ci != cj || si != sj { continue; }
            // Reciprocal overlap fraction.
            let inter_s = ai.max(aj); let inter_e = bi.min(bj);
            if inter_e <= inter_s { continue; }
            let inter = (inter_e - inter_s) as f64;
            let li = (bi - ai) as f64; let lj = (bj - aj) as f64;
            let r = (inter / li).min(inter / lj);
            if r >= min_recip { union(&mut parent, i, j); }
        }
    }

    // Group flat indices by root.
    use std::collections::BTreeMap;
    let mut groups: BTreeMap<usize, Vec<ExonRef>> = BTreeMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push((all[i].0, all[i].1));
    }
    groups.into_values().collect()
}

use std::collections::HashSet;

fn fnv1a(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes { h ^= b as u64; h = h.wrapping_mul(0x100000001b3); }
    h
}

pub(crate) fn minimizers(seq: &[u8], k: usize, w: usize) -> HashSet<u64> {
    let mut out = HashSet::new();
    if seq.len() < k { return out; }
    let n = seq.len() - k + 1;
    for win_start in 0..n.saturating_sub(w).max(1) {
        let win_end = (win_start + w).min(n);
        let mut best: Option<u64> = None;
        for i in win_start..win_end {
            let h = fnv1a(&seq[i..i + k]);
            best = Some(best.map_or(h, |b| b.min(h)));
        }
        if let Some(h) = best { out.insert(h); }
    }
    out
}

/// Default similarity bar for unifying homologous exons across paralog copies
/// (both the Stage-1b cross-copy merge and the Stage-2 within-cluster refine).
/// Overridable per-run via `RUSTLE_VG_FAMILY_MERGE_JACCARD` (applied to both
/// gates by `family_merge_jaccard()`).
///
/// 0.30 from an indel-inclusive TDD study + assembly benchmark (2026-06-06; see
/// the `merge_*` tests below + bench/multi_copy_eval/MERGE_JACCARD_THRESHOLD.md):
///
/// The main VG assembly path (pipeline.rs em_partitions) previously HARDCODED the
/// cross-copy merge bar at **0.05** — far too low. On the flagship dispersed
/// paralog (RBMY / LOC129530243) recovered RefSeq transcripts are MONOTONE in the
/// bar: **0.05→5, 0.30→7, 0.40→8**, because a low bar OVER-merges near-identical-
/// but-DISTINCT copies (union-find folds a copy on a single above-bar edge and
/// chains transitively; the completion gate at vg.rs:2803 then fabricates phantom
/// exons on coverage asymmetry with no similarity-magnitude guard), collapsing
/// copy-distinct isoforms. Nothing on the 15-locus win-loci panel regresses when
/// raising 0.05→0.30. Adversarial constructions confirm false-merge modes that
/// 0.30 refuses but a low bar accepts: embedded shared TE/domain (J≈0.22), poly-A/
/// low-complexity cores (J≈0.17). A pure sequence-SEPARATION view (unrelated exons
/// share zero k=15 minimizers, ceiling 0.000) would even tolerate a much lower
/// bar — but that view ignores the assembly over-merge cost, which is the binding
/// constraint. 0.30 is the benchmark-coherent operating point (0.40 is a candidate
/// pending genome-wide validation). The remaining lever for starved-copy recovery
/// is coverage/flow + a similarity-magnitude guard on the completion gate.
pub const FAMILY_MERGE_JACCARD_DEFAULT: f64 = 0.30;

/// Resolves the merge bar from `RUSTLE_VG_FAMILY_MERGE_JACCARD` or the default,
/// applied CONSISTENTLY to both the Stage-1b cross-copy merge and the Stage-2
/// within-cluster refine (so a single env knob controls the whole merge decision
/// — required for a clean threshold sweep).
pub fn family_merge_jaccard() -> f64 {
    std::env::var("RUSTLE_VG_FAMILY_MERGE_JACCARD")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(FAMILY_MERGE_JACCARD_DEFAULT)
}

/// Within a position-overlap cluster, split by minimizer-Jaccard similarity.
pub fn refine_by_minimizer_jaccard(
    cluster: &[(CopyId, Vec<u8>)],
    min_jaccard: f64,
    k: usize,
    w: usize,
) -> Vec<Vec<CopyId>> {
    let n = cluster.len();
    let mins: Vec<HashSet<u64>> = cluster.iter().map(|(_, s)| minimizers(s, k, w)).collect();

    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], x: usize) -> usize {
        if p[x] == x { x } else { let r = find(p, p[x]); p[x] = r; r }
    }
    fn union(p: &mut [usize], a: usize, b: usize) {
        let ra = find(p, a); let rb = find(p, b);
        if ra != rb { p[ra] = rb; }
    }

    for i in 0..n {
        for j in (i + 1)..n {
            let inter = mins[i].intersection(&mins[j]).count() as f64;
            let union_sz = mins[i].union(&mins[j]).count() as f64;
            if union_sz == 0.0 { continue; }
            if inter / union_sz >= min_jaccard { union(&mut parent, i, j); }
        }
    }
    use std::collections::BTreeMap;
    let mut g: BTreeMap<usize, Vec<CopyId>> = BTreeMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        g.entry(r).or_default().push(cluster[i].0);
    }
    g.into_values().collect()
}

/// Augment position clusters by merging singleton exons across copies via
/// sequence similarity.
///
/// Position clusters correctly group overlapping copies (tandem arrays) but
/// leave non-overlapping copies (e.g. RBMY copies scattered across the Y
/// chromosome) as singletons. This function:
///   1. Keeps multi-copy position clusters unchanged.
///   2. Collects singleton exons (one per cluster).
///   3. Runs cross-copy minimizer-Jaccard clustering on those singletons.
///   4. Merges singleton clusters that are sequence-similar.
///
/// Same-copy exons are never merged (a copy's exon 3 cannot be in the same
/// ExonClass as its exon 7, even if they look similar).
///
/// Returns the same `Vec<Vec<ExonRef>>` format as `cluster_by_position`.
pub fn merge_singletons_by_sequence(
    pos_clusters: Vec<Vec<ExonRef>>,
    copies: &[(&str, char, Vec<(u64, u64)>)],
    genome: &crate::genome::GenomeIndex,
    min_jaccard: f64,
) -> Vec<Vec<ExonRef>> {
    use std::collections::BTreeMap;
    // Env overrides for the cross-copy exon-unification (the lever that gates ALL
    // structure-sharing: TOPO_BORROW / completion / borrow project via the SHARED
    // ExonClasses this merge produces). The default 0.30 minimizer-Jaccard bar is
    // too strict for 95-99%-identical dispersed paralog pairs (k=15 minimizers are
    // SNP-sensitive; ~9 SNPs in a 300 bp exon push Jaccard to ~0.3), so homologous
    // exons stay single-copy and the borrow machinery is inert (validated
    // 2026-06-03). Unrelated exons share ~0 minimizers, so a much lower bar still
    // cleanly separates homologous from unrelated. `_K`/`_W` retune the minimizer
    // scheme for robustness; `_JACCARD` overrides the merge bar.
    let min_jaccard = std::env::var("RUSTLE_VG_FAMILY_MERGE_JACCARD")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(min_jaccard);
    let mk: usize = std::env::var("RUSTLE_VG_FAMILY_MERGE_K")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(15);
    let mw: usize = std::env::var("RUSTLE_VG_FAMILY_MERGE_W")
        .ok().and_then(|v| v.parse().ok()).unwrap_or(10);

    // Split into genuinely-multi-COPY clusters (keep as-is) and single-copy
    // clusters (fed to the cross-copy sequence merge below).
    //
    // BUG FIX (RUSTLE_VG_FAMILY_MERGE_SAME_COPY, opt-in): the original code split on
    // cluster LENGTH (`len > 1` = "multi, keep"). But for DISPERSED paralog copies,
    // a length>1 cluster is usually SAME-COPY isoform variants grouped by position —
    // a SINGLE-copy cluster that was wrongly treated as "already unified" and so
    // bypassed the cross-copy merge. Its homologous core exon then never unified
    // with the other copies → single-copy ExonClasses → the entire structure-sharing
    // machinery (TOPO_BORROW / completion / borrow) is inert at moderate identity
    // (validated 2026-06-03). The fix: route a cluster by the number of distinct
    // COPIES it spans, not its exon count, so a same-copy variant cluster enters the
    // merge (represented by one exon) and can unify with homologous clusters in
    // sibling copies. Default OFF preserves byte-identical behaviour.
    let route_by_copy = std::env::var_os("RUSTLE_VG_FAMILY_MERGE_SAME_COPY").is_some();
    let mut multi: Vec<Vec<ExonRef>> = Vec::new();
    let mut singletons: Vec<ExonRef> = Vec::new(); // one representative ExonRef per single-copy cluster
    for cluster in pos_clusters {
        let n_copies_in: usize = {
            let mut cs: Vec<usize> = cluster.iter().map(|&(c, _)| c).collect();
            cs.sort_unstable();
            cs.dedup();
            cs.len()
        };
        let is_single_copy = if route_by_copy { n_copies_in <= 1 } else { cluster.len() <= 1 };
        if !is_single_copy {
            multi.push(cluster);
        } else if let Some(&er) = cluster.first() {
            // Representative exon for the merge; the rest of the cluster's
            // (same-copy variant) exons travel with it.
            singletons.push(er);
            for &extra in cluster.iter().skip(1) {
                multi.push(vec![extra]); // keep variants as their own (copy-private) classes
            }
        }
    }

    if singletons.is_empty() {
        // No singletons to merge — early exit (all clusters already multi-copy).
        multi.extend(singletons.into_iter().map(|er| vec![er]));
        return multi;
    }

    // Compute minimizer sets for all singleton exons.
    let seqs_and_refs: Vec<(ExonRef, Vec<u8>)> = singletons.iter()
        .map(|&(cid, ei)| {
            let (chrom, _, exons) = &copies[cid];
            let (s, e) = exons[ei];
            let seq = genome.fetch_sequence(chrom, s, e).unwrap_or_default();
            ((cid, ei), seq)
        })
        .collect();

    let n = seqs_and_refs.len();
    let mins: Vec<HashSet<u64>> = seqs_and_refs.iter()
        .map(|(_, s)| minimizers(s, mk, mw))
        .collect();

    let mut parent: Vec<usize> = (0..n).collect();
    fn find_s(p: &mut [usize], x: usize) -> usize {
        if p[x] == x { x } else { let r = find_s(p, p[x]); p[x] = r; r }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            let (cid_i, _) = seqs_and_refs[i].0;
            let (cid_j, _) = seqs_and_refs[j].0;
            if cid_i == cid_j { continue; } // never merge same-copy exons
            let inter = mins[i].intersection(&mins[j]).count() as f64;
            let union_sz = mins[i].union(&mins[j]).count() as f64;
            if union_sz == 0.0 { continue; }
            if inter / union_sz >= min_jaccard {
                let ri = find_s(&mut parent, i);
                let rj = find_s(&mut parent, j);
                if ri != rj { parent[ri] = rj; }
            }
        }
    }

    let mut merged: BTreeMap<usize, Vec<ExonRef>> = BTreeMap::new();
    for i in 0..n {
        let r = find_s(&mut parent, i);
        merged.entry(r).or_default().push(seqs_and_refs[i].0);
    }

    // Combine multi-copy position clusters with the merged singleton groups.
    let mut result = multi;
    result.extend(merged.into_values());
    result
}

/// A junction observation before binding to graph nodes.
#[derive(Debug, Clone, Copy)]
pub struct RawJunction {
    pub donor: u64,
    pub acceptor: u64,
    pub strand: char,
    pub family_support: u32,
}

/// Aggregate per-copy junctions into family junctions with copy-support counts.
pub fn collect_family_junctions(per_copy: &[(char, Vec<(u64, u64)>)]) -> Vec<RawJunction> {
    use std::collections::BTreeMap;
    // Round to nearest 10 bp to absorb long-read jitter (matches existing
    // build_family_consensus_junctions in vg.rs).
    let mut counts: BTreeMap<(char, u64, u64), u32> = BTreeMap::new();
    for (strand, juncs) in per_copy {
        let mut seen: HashSet<(u64, u64)> = HashSet::new();
        for &(d, a) in juncs {
            let key = (d / 10 * 10, a / 10 * 10);
            if seen.insert(key) {
                *counts.entry((*strand, key.0, key.1)).or_insert(0) += 1;
            }
        }
    }
    counts.into_iter()
        .map(|((strand, d, a), n)| RawJunction { donor: d, acceptor: a, strand, family_support: n })
        .collect()
}

use crate::genome::GenomeIndex;
use crate::vg::FamilyGroup;
use anyhow::Result;

/// Shared cluster→FamilyGraph assembly used by BOTH `build_family_graph` (bundle
/// path) and `build_family_graph_from_layer1_graphs` (Layer-1-graph path).
/// per_copy_exons[c] = (chrom, strand, sorted exon spans); per_copy_juncs[c] =
/// (strand, junctions). The three thresholds pass through to the existing stages.
fn assemble_family_graph_from_copy_exons(
    family_id: usize,
    per_copy_exons: &[(String, char, Vec<(u64, u64)>)],
    per_copy_juncs: &[(char, Vec<(u64, u64)>)],
    genome: Option<&GenomeIndex>,
    min_pos_recip: f64,
    merge_min_jaccard: f64,
    refine_min_jaccard: f64,
) -> Result<FamilyGraph> {
    // 1. Collect (chrom, strand, exons) per copy. Borrow the owned chrom String
    //    as &str so the rest of the body is verbatim from build_family_graph.
    let copies: Vec<(&str, char, Vec<(u64, u64)>)> = per_copy_exons.iter()
        .map(|(c, s, e)| (c.as_str(), *s, e.clone()))
        .collect();
    if copies.is_empty() {
        return Ok(FamilyGraph { family_id, nodes: Vec::new(), edges: Vec::new() });
    }

    // Anyhow-error if strands are mixed — caller should have filtered this.
    let strand0 = copies[0].1;
    if copies.iter().any(|(_, s, _)| *s != strand0) {
        anyhow::bail!("family {} has mixed strands; build_family_graph requires single-strand families", family_id);
    }

    // 2. Stage 1: position-overlap clusters.
    let pos_clusters = cluster_by_position(&copies, min_pos_recip);

    // Stage 1b: if position clustering produced all singletons (copies don't
    // overlap genomically, e.g. RBMY / TSPY scattered across a chromosome),
    // fall back to cross-copy sequence similarity clustering. This groups
    // structurally equivalent exons from different loci so ExonClasses span
    // multiple copies — enabling coverage borrowing and graph completion.
    // Guard: skip if total exons is very large (O(n²) cost).
    let trace_compl = std::env::var_os("RUSTLE_VG_COMPLETION_TRACE").is_some();
    // Stage 1b: merge singleton position clusters by sequence similarity.
    // Extends ExonClass grouping to non-overlapping copies (e.g. RBMY / TSPY
    // scattered across a chromosome). Multi-copy position clusters are kept
    // intact; only un-grouped singleton exons are re-clustered by sequence.
    // Guard: skip if too many exons (O(n²) cost on singletons).
    let effective_clusters: Vec<Vec<ExonRef>> = {
        let total_exons: usize = copies.iter().map(|(_, _, e)| e.len()).sum();
        let n_singletons = pos_clusters.iter().filter(|c| c.len() == 1).count();
        if trace_compl {
            eprintln!("[FG] family={} total_exons={} pos_clusters={} singletons={} has_genome={}",
                family_id, total_exons, pos_clusters.len(), n_singletons, genome.is_some());
        }
        if n_singletons >= 2 && total_exons <= 2000 {
            if let Some(g) = genome {
                let merged = merge_singletons_by_sequence(pos_clusters, &copies, g, merge_min_jaccard);
                if trace_compl {
                    let n_multi = merged.iter().filter(|c| c.len() > 1).count();
                    eprintln!("[FG]   after merge: clusters={} multi-copy={}", merged.len(), n_multi);
                }
                merged
            } else {
                pos_clusters
            }
        } else {
            pos_clusters
        }
    };

    // 3. Stage 2: refine by minimizer Jaccard (requires sequences). If no genome,
    //    skip refinement — every position cluster becomes one ExonClass.
    let mut nodes: Vec<ExonClass> = Vec::new();
    for cluster in &effective_clusters {
        // Recover sequences (or empty stubs if no genome) and original spans.
        let with_seq: Vec<(CopyId, Vec<u8>)> = cluster.iter().map(|&(cid, ei)| {
            let (chrom, _, exons) = &copies[cid];
            let (s, e) = exons[ei];
            let seq = match genome {
                Some(g) => g.fetch_sequence(chrom, s, e).unwrap_or_default(),
                None => Vec::new(),
            };
            (cid, seq)
        }).collect();
        // If we have sequences, refine; otherwise skip.
        let groups: Vec<Vec<CopyId>> = if with_seq.iter().all(|(_, s)| !s.is_empty()) {
            refine_by_minimizer_jaccard(&with_seq, refine_min_jaccard, 15, 10)
        } else {
            vec![with_seq.iter().map(|(c, _)| *c).collect()]
        };
        for cids in groups {
            let mut per_copy_sequences: Vec<(CopyId, Vec<u8>)> = with_seq.iter()
                .filter(|(c, _)| cids.contains(c))
                .cloned().collect();
            per_copy_sequences.sort_by_key(|(c, _)| *c);
            // Original per-copy spans (lossless preservation of input coords).
            let mut per_copy_spans: Vec<(CopyId, (u64, u64))> = Vec::with_capacity(cids.len());
            // Representative span = union of contributing exon spans.
            let mut s_min = u64::MAX; let mut e_max = 0u64;
            for &cid in &cids {
                let pos_in_cluster = cluster.iter().find(|&&(c, _)| c == cid)
                    .ok_or_else(|| anyhow::anyhow!(
                        "copy {} produced by minimizer-Jaccard refinement not found in originating position cluster (invariant violated)", cid))?;
                let (_, _, exs) = &copies[pos_in_cluster.0];
                let (s, e) = exs[pos_in_cluster.1];
                per_copy_spans.push((cid, (s, e)));
                s_min = s_min.min(s); e_max = e_max.max(e);
            }
            per_copy_spans.sort_by_key(|(c, _)| *c);
            nodes.push(ExonClass {
                idx: NodeIdx(nodes.len()),
                chrom: copies[0].0.to_string(),
                span: (s_min, e_max),
                strand: strand0,
                per_copy_sequences,
                per_copy_spans,
                copy_specific: cids.len() == 1,
                per_copy_cov: Vec::new(),
            });
        }
    }

    // 4. Junction edges: collect junctions from each copy (passed in by the
    //    caller), bind to nodes by mapping (donor, acceptor) → nodes whose span
    //    ends at donor / starts at acceptor.
    let raw = collect_family_junctions(per_copy_juncs);

    // r.donor / r.acceptor are rounded to nearest 10 bp by collect_family_junctions
    // (jitter absorption). Search per_copy_spans (not span) so that multi-copy
    // ExonClasses built by sequence clustering — where span is the union and may
    // cover the whole chromosome — still bind junctions correctly.
    // Group by ExonClass pair (f, t) and sum family_support so that
    // non-overlapping copies at different absolute positions (e.g. RBMY)
    // contribute to the same structural edge rather than creating N separate
    // edges each with support=1.
    use std::collections::HashMap as HMap;
    let mut edge_map: HMap<(NodeIdx, NodeIdx), (u32, char)> = HMap::new();
    for r in &raw {
        let from = nodes.iter().find(|n| {
            n.per_copy_spans.iter().any(|(_, (_, e))| {
                e / 10 * 10 == r.donor || *e == r.donor
            })
        }).map(|n| n.idx);
        let to = nodes.iter().find(|n| {
            n.per_copy_spans.iter().any(|(_, (s, _))| {
                s / 10 * 10 == r.acceptor || *s == r.acceptor
            })
        }).map(|n| n.idx);
        if let (Some(f), Some(t)) = (from, to) {
            let e = edge_map.entry((f, t)).or_insert((0, r.strand));
            e.0 += r.family_support;
        }
    }
    let edges: Vec<JunctionEdge> = edge_map.into_iter()
        .map(|((f, t), (support, strand))| JunctionEdge {
            from: f, to: t, family_support: support, strand,
        })
        .collect();

    if trace_compl {
        let multi_edges = edges.iter().filter(|e| e.family_support >= 2).count();
        eprintln!("[FG] family={} edges={} edges_support>=2={}",
                  family_id, edges.len(), multi_edges);
        for e in edges.iter().filter(|e| e.family_support >= 2).take(5) {
            eprintln!("[FG]   edge {:?}->{:?} support={}", e.from, e.to, e.family_support);
        }
    }

    Ok(FamilyGraph { family_id, nodes, edges })
}

pub fn build_family_graph(
    family: &FamilyGroup,
    bundles: &[Bundle],
    genome: Option<&GenomeIndex>,
    min_pos_recip: f64,
    merge_min_jaccard: f64,
    refine_min_jaccard: f64,
) -> Result<FamilyGraph> {
    // 1. Collect (chrom, strand, exons) per copy from the family's bundles.
    let per_copy_exons: Vec<(String, char, Vec<(u64, u64)>)> = family.bundle_indices.iter()
        .map(|&bi| {
            let b = &bundles[bi];
            (b.chrom.clone(), b.strand, extract_copy_exons(b))
        })
        .collect();
    // 2. Per-copy junctions from each bundle's junction_stats.
    let per_copy_juncs: Vec<(char, Vec<(u64, u64)>)> = family.bundle_indices.iter()
        .map(|&bi| {
            let b = &bundles[bi];
            let mut js = Vec::new();
            for (j, _) in &b.junction_stats {
                js.push((j.donor, j.acceptor));
            }
            (b.strand, js)
        })
        .collect();
    assemble_family_graph_from_copy_exons(
        family.family_id,
        &per_copy_exons,
        &per_copy_juncs,
        genome,
        min_pos_recip,
        merge_min_jaccard,
        refine_min_jaccard,
    )
}

/// Layer-2 (C3) family-graph merge sourced from Layer-1 splice graphs.
/// Each graph contributes its exon-node spans as one copy and its node adjacency
/// as that copy's junctions. Delegates to the shared routine (fetches sequences,
/// binds edges). No coordinate is invented.
pub fn build_family_graph_from_layer1_graphs(
    family_id: usize,
    copies: &[(String, char, &crate::graph::Graph)],
    genome: Option<&GenomeIndex>,
    min_pos_recip: f64,
    merge_min_jaccard: f64,
    refine_min_jaccard: f64,
) -> Result<FamilyGraph> {
    let per_copy_exons: Vec<(String, char, Vec<(u64, u64)>)> = copies
        .iter()
        .map(|(chrom, strand, g)| {
            let mut exons: Vec<(u64, u64)> = g
                .nodes
                .iter()
                .filter(|n| n.node_id != g.source_id && n.node_id != g.sink_id && n.end > n.start)
                .map(|n| (n.start, n.end))
                .collect();
            exons.sort_unstable();
            exons.dedup();
            (chrom.clone(), *strand, exons)
        })
        .collect();
    let per_copy_juncs: Vec<(char, Vec<(u64, u64)>)> = copies
        .iter()
        .map(|(_, strand, g)| {
            // node_id == vec index (graph.rs add_node sets node_id = nodes.len()),
            // so g.nodes[child] indexes the child node directly.
            let mut js: Vec<(u64, u64)> = Vec::new();
            for n in &g.nodes {
                if n.node_id == g.source_id || n.node_id == g.sink_id || !(n.end > n.start) {
                    continue;
                }
                for child in n.children.ones() {
                    if child == g.sink_id {
                        continue;
                    }
                    let c = &g.nodes[child];
                    if c.start > n.end {
                        js.push((n.end, c.start));
                    }
                }
            }
            js.sort_unstable();
            js.dedup();
            (*strand, js)
        })
        .collect();
    assemble_family_graph_from_copy_exons(
        family_id,
        &per_copy_exons,
        &per_copy_juncs,
        genome,
        min_pos_recip,
        merge_min_jaccard,
        refine_min_jaccard,
    )
}

/// Multiple sequence alignment via poasta POA graph traversal.
/// Returns rows of equal length, padded with b'-' for gaps.
/// One row per input sequence, in the same order.
/// Errors if fewer than 2 sequences provided.
///
/// Pure POA-MSA utility used by the graph-POA-identity family filter
/// (`compute_family_graph_poa_identity_diag`). No HMM/profile dependency.
pub fn poa_msa(seqs: &[Vec<u8>]) -> Result<Vec<Vec<u8>>> {
    use anyhow::anyhow;
    use poasta::graphs::poa::{POAGraph, POANodeIndex};
    use poasta::aligner::PoastaAligner;
    use poasta::aligner::config::AffineDijkstra;
    use poasta::aligner::scoring::{GapAffine, AlignmentType};
    use poasta::io::fasta::poa_graph_to_fasta;

    if seqs.len() < 2 {
        return Err(anyhow!("poa_msa requires at least 2 sequences"));
    }

    // Build a POA graph by aligning each sequence progressively.
    // Uses affine-gap Dijkstra (no A* heuristic) with mismatch=1, gap_extend=1, gap_open=2.
    let gap_costs = GapAffine::new(1, 1, 2);
    let aligner: PoastaAligner<'_, AffineDijkstra> =
        PoastaAligner::new(AffineDijkstra(gap_costs), AlignmentType::Global);

    let mut graph: POAGraph<u32> = POAGraph::new();
    let unit_weights: Vec<usize> = vec![1; seqs.iter().map(|s| s.len()).max().unwrap_or(0)];

    for (i, seq) in seqs.iter().enumerate() {
        let w = &unit_weights[..seq.len()];
        if graph.is_empty() {
            // First sequence: add unaligned (no existing graph to align to).
            graph
                .add_alignment_with_weights(&format!("seq{i}"), seq, None, w)
                .map_err(|e| anyhow!("poasta add_alignment_with_weights failed: {e}"))?;
        } else {
            // Subsequent sequences: align against the current graph.
            let aln_result = aligner.align::<u32, _>(&graph, seq);
            let alignment = if aln_result.alignment.is_empty() {
                None
            } else {
                Some(&aln_result.alignment as &poasta::aligner::Alignment<POANodeIndex<u32>>)
            };
            graph
                .add_alignment_with_weights(&format!("seq{i}"), seq, alignment, w)
                .map_err(|e| anyhow!("poasta add_alignment_with_weights failed: {e}"))?;
        }
    }

    // Extract the MSA rows by using poasta's public poa_graph_to_fasta function,
    // which walks the graph in topological order and writes aligned FASTA records.
    // We write into a Vec<u8> buffer and parse the resulting FASTA lines.
    let mut fasta_buf: Vec<u8> = Vec::new();
    poa_graph_to_fasta(&graph, &mut fasta_buf)
        .map_err(|e| anyhow!("poa_graph_to_fasta failed: {e}"))?;

    // Parse the FASTA output: each sequence record is a set of lines starting with '>'.
    // We collect the sequences in order and convert them to Vec<u8> rows.
    let mut msa: Vec<Vec<u8>> = Vec::with_capacity(seqs.len());
    let mut current_seq: Option<Vec<u8>> = None;
    for line in fasta_buf.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if let Some(seq) = current_seq.take() {
                msa.push(seq);
            }
            current_seq = Some(Vec::new());
        } else if let Some(ref mut seq) = current_seq {
            seq.extend_from_slice(line);
        }
    }
    if let Some(seq) = current_seq {
        msa.push(seq);
    }

    if msa.len() != seqs.len() {
        return Err(anyhow!(
            "poa_msa: expected {} rows from graph but got {}",
            seqs.len(),
            msa.len()
        ));
    }

    // poasta's fasta_aln_for_seq has an off-by-one in the trailing-gap fill for
    // sequences whose path ends before max_col, so rows can occasionally come
    // back one column short. Pad with trailing gaps to the max row length —
    // this is content-neutral for an MSA (gaps are not part of the ungapped
    // sequence) and the round-trip ungap(row) ≡ original input is preserved.
    let n_cols = msa.iter().map(|r| r.len()).max().unwrap_or(0);
    for row in msa.iter_mut() {
        if row.len() < n_cols {
            row.extend(std::iter::repeat(b'-').take(n_cols - row.len()));
        }
    }
    Ok(msa)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refine_zero_threshold_keeps_all_merged() {
        let seq_a: Vec<u8> = b"ATCGATCGATCGATCGATCGATCGATCGATCG"
            .iter().cycle().take(200).copied().collect();
        let seq_b: Vec<u8> = b"GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
            .iter().cycle().take(200).copied().collect();

        let cluster = vec![(0usize, seq_a.clone()), (1usize, seq_b.clone())];

        let groups_zero = refine_by_minimizer_jaccard(&cluster, 0.0, 15, 10);
        assert_eq!(groups_zero.len(), 1,
            "threshold=0.0 should keep all sequences in one ExonClass");
        assert_eq!(groups_zero[0].len(), 2);

        let groups_high = refine_by_minimizer_jaccard(&cluster, 0.30, 15, 10);
        assert_eq!(groups_high.len(), 2,
            "sequences share no k-mers (Jaccard=0.0); any threshold > 0.0 splits them");
    }

    #[test]
    fn refine_identical_sequences_always_merged() {
        let seq: Vec<u8> = b"ATCGATCGATCGATCGATCG".iter().cycle().take(200).copied().collect();
        let cluster = vec![(0usize, seq.clone()), (1usize, seq.clone())];

        for threshold in [0.0_f64, 0.05, 0.30, 0.99] {
            let groups = refine_by_minimizer_jaccard(&cluster, threshold, 15, 10);
            assert_eq!(groups.len(), 1,
                "identical sequences should always stay merged (threshold={threshold})");
        }
    }

    // --- Merge-Jaccard sensitivity sweep (analysis harness, not a pass/fail test) ---
    //
    // Characterises the EXACT metric the family-merge step uses (`minimizers`
    // k=15/w=10 + Jaccard, threshold env RUSTLE_VG_FAMILY_MERGE_JACCARD, default
    // 0.30). Produces (1) a controlled SNP-divergence curve on real exon
    // sequences and (2) an unrelated-exon baseline, so the separation margin and
    // a recommended threshold can be read off directly.
    //
    // Guarded by env MERGE_SWEEP_DATA (path to a `id<TAB>seq` TSV) and #[ignore],
    // so it never runs in the normal suite. Run with:
    //   MERGE_SWEEP_DATA=/tmp/merge_sweep_exons.tsv \
    //     cargo test -p rustle merge_jaccard_sensitivity_sweep -- --ignored --nocapture
    struct SplitMix64(u64);
    impl SplitMix64 {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        }
        fn unit(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn jaccard(a: &HashSet<u64>, b: &HashSet<u64>) -> f64 {
        let inter = a.intersection(b).count() as f64;
        let uni = a.union(b).count() as f64;
        if uni == 0.0 { 0.0 } else { inter / uni }
    }

    /// Substitution-only mutation at per-base rate `r`; returns (mutated, #SNPs).
    fn mutate_snps(seq: &[u8], r: f64, rng: &mut SplitMix64) -> (Vec<u8>, usize) {
        const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let mut out = seq.to_vec();
        let mut n = 0;
        for x in out.iter_mut() {
            if rng.unit() < r {
                let cur = *x;
                // pick a different base
                let mut nb = B[(rng.next_u64() % 4) as usize];
                while nb == cur { nb = B[(rng.next_u64() % 4) as usize]; }
                *x = nb;
                n += 1;
            }
        }
        (out, n)
    }

    /// COHERENCE GUARD (permanent, runs in the normal suite).
    ///
    /// The family-merge default bar must BRACKET the real homology signal, once
    /// indels are modelled, on BOTH sides:
    ///   (lower bound) genuinely near-identical paralog copies — the real regime,
    ///     ≤~1% substitution + small indels, e.g. RBMY at 98-99.8% identity —
    ///     must unify into ONE `ExonClass` (else the borrow/completion machinery
    ///     is inert); guards against the bar being too HIGH.
    ///   (upper bound — separate adversarial guards below) unrelated / TE-sharing
    ///     / low-complexity exons must NEVER merge; guards against too LOW.
    ///
    /// NOTE the divergence here is ≤1% (NOT 2-3%): the indel-inclusive sweep +
    /// the RBMY assembly benchmark showed that forcing the 2-3% band to merge
    /// (which a ~0.10 bar would do) OVER-merges distinct copies and REGRESSES
    /// recovered transcripts. So the coherent property is "merge the near-
    /// identical, leave the moderately-diverged distinct" — exactly what 0.30
    /// does. Uses the production metric (minimizer-Jaccard, k=15, w=10) at the
    /// production default. Hermetic panel: bench/multi_copy_eval/merge_sweep_exons.tsv.
    #[test]
    fn family_merge_default_merges_near_identical_paralogs() {
        let panel = load_panel();
        let fam: Vec<&(String, Vec<u8>)> =
            panel.iter().filter(|(id, _)| id.starts_with("LRPAP1")).collect();
        let unrel: Vec<&(String, Vec<u8>)> =
            panel.iter().filter(|(id, _)| id.starts_with("UNREL")).collect();
        let thr = FAMILY_MERGE_JACCARD_DEFAULT;

        // Near-identical (≤1% div + small indels) copy-vs-copy recall at default bar.
        let (mut merged, mut total) = (0usize, 0usize);
        for (ei, exon) in fam.iter().enumerate() {
            let mut rng = SplitMix64(0x00C0_FFEE_u64 ^ ei as u64);
            for _ in 0..40 {
                let (a, _, _) = mutate_indels(&exon.1, 0.008, 0.002, 3, &mut rng);
                let (b, _, _) = mutate_indels(&exon.1, 0.008, 0.002, 3, &mut rng);
                let j = jaccard(&minimizers(&a, 15, 10), &minimizers(&b, 15, 10));
                if j >= thr { merged += 1; }
                total += 1;
            }
        }
        let recall = merged as f64 / total as f64;
        assert!(recall >= 0.85,
            "default merge bar {thr} fails to unify NEAR-IDENTICAL paralog copies \
             (only {:.0}% merge; need >=85%) — bar is too HIGH", recall * 100.0);

        // Unrelated exons must never merge at the same bar (specificity).
        let mins: Vec<HashSet<u64>> = unrel.iter().map(|(_, s)| minimizers(s, 15, 10)).collect();
        for i in 0..mins.len() {
            for j in (i + 1)..mins.len() {
                let s = jaccard(&mins[i], &mins[j]);
                assert!(s < thr,
                    "unrelated exons would falsely merge at bar {thr} (sim={s:.3})");
            }
        }
    }

    fn rand_seq(n: usize, rng: &mut SplitMix64) -> Vec<u8> {
        const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
        (0..n).map(|_| B[(rng.next_u64() % 4) as usize]).collect()
    }

    /// ADVERSARIAL GUARD: a short exon-fragment fully embedded inside a larger,
    /// otherwise-unrelated exon must NOT cause a merge at the default bar. This is
    /// the "small-set-in-big-set" failure mode of a containment metric: the small
    /// fragment's minimizers ⊆ the big exon's, so containment ≈ 1.0 (would merge),
    /// but these are NOT homologous whole exons. The production Jaccard metric
    /// (union denominator) correctly refuses. This test is WHY the default stays
    /// Jaccard rather than switching to containment.
    #[test]
    fn fragment_in_larger_exon_does_not_merge_under_jaccard_default() {
        let mut rng = SplitMix64(0x0BAD_F00D);
        let small = rand_seq(45, &mut rng);            // a short exon fragment
        let mut big = rand_seq(300, &mut rng);
        big.extend_from_slice(&small);                 // fragment literally inside big
        big.extend(rand_seq(300, &mut rng));           // unrelated flanks
        let (ms, mb) = (minimizers(&small, 15, 10), minimizers(&big, 15, 10));

        let cont = containment(&ms, &mb);
        let jacc = jaccard(&ms, &mb);
        // containment WOULD merge (this documents the risk we are avoiding)…
        assert!(cont >= 0.90,
            "sanity: fragment is contained, containment should be ~1 (got {cont:.3})");
        // …but the production Jaccard default must NOT merge them.
        assert!(jacc < FAMILY_MERGE_JACCARD_DEFAULT,
            "Jaccard default {} must not merge a fragment into an unrelated larger exon \
             (jacc={jacc:.3})", FAMILY_MERGE_JACCARD_DEFAULT);
    }

    /// ADVERSARIAL GUARD: two genuinely unrelated exons that happen to share a
    /// low-complexity tract (e.g. a (CAG)_n microsatellite or a common short
    /// motif) must NOT merge at the default bar — the shared repeat contributes
    /// few distinct minimizers, so Jaccard stays well below the bar.
    #[test]
    fn shared_low_complexity_tract_does_not_merge_at_default() {
        let mut rng = SplitMix64(0xFEED_BEEF);
        let repeat: Vec<u8> = b"CAG".iter().cycle().take(60).copied().collect();
        let mut a = rand_seq(280, &mut rng); a.extend_from_slice(&repeat); a.extend(rand_seq(20, &mut rng));
        let mut b = rand_seq(20, &mut rng);  b.extend_from_slice(&repeat); b.extend(rand_seq(280, &mut rng));
        let j = jaccard(&minimizers(&a, 15, 10), &minimizers(&b, 15, 10));
        assert!(j < FAMILY_MERGE_JACCARD_DEFAULT,
            "unrelated exons sharing only a low-complexity tract must not merge at bar {} (jacc={j:.3})",
            FAMILY_MERGE_JACCARD_DEFAULT);
    }

    /// ADVERSARIAL GUARD (from the merge-metric adversarial workflow): two
    /// otherwise-unrelated exons that share ONE embedded high-complexity block
    /// (stand-in for a common Alu/L1 fragment or protein domain) must NOT merge
    /// at the default bar. Measured Jaccard ≈ 0.2 — this is EXACTLY a case a 0.10
    /// bar would wrongly merge (fabricating a chimeric ExonClass across genes),
    /// and a key reason the default stays 0.30.
    #[test]
    fn embedded_shared_element_does_not_merge_at_default() {
        let element = { let mut r = SplitMix64(0xABCD_0001); rand_seq(110, &mut r) };
        let mut rng = SplitMix64(0x1111_2222);
        let mut a = rand_seq(95, &mut rng); a.extend_from_slice(&element); a.extend(rand_seq(95, &mut rng));
        let mut b = rand_seq(95, &mut rng); b.extend_from_slice(&element); b.extend(rand_seq(95, &mut rng));
        let j = jaccard(&minimizers(&a, 15, 10), &minimizers(&b, 15, 10));
        assert!(j > 0.05, "sanity: shared element should give non-trivial overlap (j={j:.3})");
        assert!(j < FAMILY_MERGE_JACCARD_DEFAULT,
            "unrelated exons sharing an embedded element must not merge at default bar {} (j={j:.3})",
            FAMILY_MERGE_JACCARD_DEFAULT);
    }

    /// ADVERSARIAL GUARD: two unrelated exons that share a long poly-A / AT-rich
    /// core but have distinct ends must NOT merge at the default. Measured Jaccard
    /// ≈ 0.17 — again inside the danger zone of a ~0.10 bar, refused by 0.30.
    #[test]
    fn poly_a_core_does_not_merge_at_default() {
        let core: Vec<u8> = std::iter::repeat(b'A').take(600).collect();
        let mut a = core.clone(); a.extend_from_slice(b"CGTCGTCGTC");
        let mut b = core.clone(); b.extend_from_slice(b"GCTGCTGCTG");
        let j = jaccard(&minimizers(&a, 15, 10), &minimizers(&b, 15, 10));
        assert!(j < FAMILY_MERGE_JACCARD_DEFAULT,
            "unrelated exons sharing only a poly-A core must not merge at default bar {} (j={j:.3})",
            FAMILY_MERGE_JACCARD_DEFAULT);
    }

    /// CHARACTERIZATION (documents a KNOWN metric limitation, not threshold-
    /// fixable): two unrelated exons each DOMINATED by the same simple repeat
    /// collapse to a tiny minimizer set and false-merge at Jaccard ≈ 1.0 — at ANY
    /// threshold including 0.30. Real exons are rarely pure-satellite, and this is
    /// mitigated upstream by repeat-masking (diagnostic.rs `SeedMasked`,
    /// masked_fraction>0.30). Encoded so the limitation is tracked; a true fix
    /// would low-complexity-mask before minimizing.
    #[test]
    fn pure_repeat_exons_collapse_known_metric_limitation() {
        let a: Vec<u8> = b"AC".iter().cycle().take(146).copied().collect();
        let mut b: Vec<u8> = vec![b'G'];
        b.extend(b"AC".iter().cycle().take(146).copied());
        b.push(b'T');
        let j = jaccard(&minimizers(&a, 15, 10), &minimizers(&b, 15, 10));
        assert!(j >= 0.90,
            "documents the known low-complexity collapse (pure-repeat exons false-merge); got j={j:.3}. \
             If this dropped, a low-complexity mask was added — update the metric docs.");
    }

    fn containment(a: &HashSet<u64>, b: &HashSet<u64>) -> f64 {
        let inter = a.intersection(b).count() as f64;
        let denom = a.len().min(b.len()) as f64;
        if denom == 0.0 { 0.0 } else { inter / denom }
    }

    /// SNP + indel mutation. `snp` = per-base substitution rate, `indel_rate` =
    /// per-position indel-event rate, `mean_indel_len` = geometric mean length.
    /// Insertions add random bases; deletions remove a run. Returns
    /// (mutated, n_snp, n_indel_events). Length-changing — this is what makes the
    /// metric choice (Jaccard vs containment) matter.
    fn mutate_indels(seq: &[u8], snp: f64, indel_rate: f64, mean_indel_len: usize,
                     rng: &mut SplitMix64) -> (Vec<u8>, usize, usize) {
        const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
        let mut out = Vec::with_capacity(seq.len() + 8);
        let (mut nsnp, mut nindel) = (0usize, 0usize);
        let mut i = 0usize;
        while i < seq.len() {
            if indel_rate > 0.0 && rng.unit() < indel_rate {
                nindel += 1;
                let mut len = 1usize;
                let cont = 1.0 - 1.0 / (mean_indel_len.max(1) as f64);
                while rng.unit() < cont { len += 1; if len >= 20 { break; } }
                if rng.next_u64() & 1 == 0 {
                    for _ in 0..len { out.push(B[(rng.next_u64() % 4) as usize]); } // insertion
                    continue; // original base processed next iteration
                } else {
                    i += len; // deletion of a run
                    continue;
                }
            }
            let base = seq[i];
            if rng.unit() < snp {
                let mut nb = B[(rng.next_u64() % 4) as usize];
                while nb == base { nb = B[(rng.next_u64() % 4) as usize]; }
                out.push(nb); nsnp += 1;
            } else {
                out.push(base);
            }
            i += 1;
        }
        (out, nsnp, nindel)
    }

    fn load_panel() -> Vec<(String, Vec<u8>)> {
        let path = concat!(env!("CARGO_MANIFEST_DIR"),
                           "/bench/multi_copy_eval/merge_sweep_exons.tsv");
        let raw = std::fs::read_to_string(path).expect("read committed exon panel");
        raw.lines().filter_map(|l| {
            let mut it = l.splitn(2, '\t');
            let id = it.next()?.to_string();
            let seq = it.next()?.trim().as_bytes().to_vec();
            if seq.len() < 15 { return None; }
            Some((id, seq))
        }).collect()
    }

    /// `n` paralog copies, each an INDEPENDENT mutation of `ancestor` (so pairwise
    /// divergence ≈ 2× per-copy — models real paralog pairs descending from a
    /// common ancestor).
    fn homolog_cluster(ancestor: &[u8], n: usize, snp: f64, indel: f64, ilen: usize,
                       seed: u64) -> Vec<(CopyId, Vec<u8>)> {
        let mut rng = SplitMix64(seed ^ 0x5DEE_CE66_u64);
        (0..n).map(|c| {
            let (m, _, _) = mutate_indels(ancestor, snp, indel, ilen, &mut rng);
            (c as CopyId, m)
        }).collect()
    }

    // Indel-inclusive separation sweep: how SNPs+indels move minimizer similarity,
    // and which (metric, k, w, threshold) operating point is most COHERENT —
    // homologs unify, unrelated stay separate, max margin. Deterministic, hermetic
    // (reads the committed exon panel). Run:
    //   cargo test -p rustle merge_metric_indel_sweep -- --ignored --nocapture
    #[test]
    #[ignore]
    fn merge_metric_indel_sweep() {
        let panel = load_panel();
        let fam: Vec<&(String, Vec<u8>)> =
            panel.iter().filter(|(id, _)| id.starts_with("LRPAP1")).collect();
        let unrel: Vec<&(String, Vec<u8>)> =
            panel.iter().filter(|(id, _)| id.starts_with("UNREL")).collect();
        const TRIALS: usize = 40;

        // Pairwise homolog similarity at a given (metric,k,w) and divergence: for
        // each real exon, mutate twice (two independent copies), measure copy-vs-copy
        // similarity (the quantity union-find thresholds). Returns the mean and the
        // fraction that would MERGE at each candidate threshold.
        let eval = |use_cont: bool, k: usize, w: usize, snp: f64, indel: f64|
            -> (f64, Vec<(f64, f64)>) {
            let thrs = [0.30_f64, 0.25, 0.20, 0.15, 0.10, 0.07, 0.05];
            let mut sims: Vec<f64> = Vec::new();
            for (ei, exon) in fam.iter().enumerate() {
                let seed = 0xD1B5_4A32_u64.wrapping_mul(ei as u64 + 1)
                    .wrapping_add((snp * 1e6) as u64 + (indel * 1e5) as u64);
                let mut rng = SplitMix64(seed);
                for _ in 0..TRIALS {
                    let (a, _, _) = mutate_indels(&exon.1, snp, indel, 3, &mut rng);
                    let (b, _, _) = mutate_indels(&exon.1, snp, indel, 3, &mut rng);
                    let (ma, mb) = (minimizers(&a, k, w), minimizers(&b, k, w));
                    sims.push(if use_cont { containment(&ma, &mb) } else { jaccard(&ma, &mb) });
                }
            }
            let mean = sims.iter().sum::<f64>() / sims.len() as f64;
            let merge_frac: Vec<(f64, f64)> = thrs.iter().map(|&t| {
                (t, sims.iter().filter(|&&s| s >= t).count() as f64 / sims.len() as f64)
            }).collect();
            (mean, merge_frac)
        };

        // unrelated false-merge ceiling for a (metric,k,w)
        let neg_ceiling = |use_cont: bool, k: usize, w: usize| -> f64 {
            let mins_un: Vec<HashSet<u64>> = unrel.iter().map(|(_, s)| minimizers(s, k, w)).collect();
            let mins_fa: Vec<HashSet<u64>> = fam.iter().map(|(_, s)| minimizers(s, k, w)).collect();
            let mut mx = 0.0_f64;
            for i in 0..mins_un.len() {
                for j in (i + 1)..mins_un.len() {
                    let s = if use_cont { containment(&mins_un[i], &mins_un[j]) }
                            else { jaccard(&mins_un[i], &mins_un[j]) };
                    mx = mx.max(s);
                }
            }
            for a in &mins_fa { for b in &mins_un {
                let s = if use_cont { containment(a, b) } else { jaccard(a, b) };
                mx = mx.max(s);
            }}
            mx
        };

        eprintln!("\n===== A. indel-inclusive divergence curve (Jaccard vs Containment, k=15,w=10) =====");
        eprintln!("model: snp=div, indel=div/3 (mean len 3); copy-vs-copy similarity, {} trials/exon", TRIALS);
        eprintln!("{:>6}  {:>18}  {:>18}", "div%", "Jaccard mean", "Containment mean");
        for div in [0.0_f64, 0.01, 0.02, 0.03, 0.05, 0.08] {
            let (jm, _) = eval(false, 15, 10, div, div / 3.0);
            let (cm, _) = eval(true, 15, 10, div, div / 3.0);
            eprintln!("{:>5.0}%  {:>18.3}  {:>18.3}", div * 100.0, jm, cm);
        }

        eprintln!("\n===== B. operating-point grid @ realistic paralog divergence (2.5% snp + 0.8% indel) =====");
        eprintln!("recall = % copy-pairs merging; ceil = max unrelated sim (must be < threshold)");
        eprintln!("{:>5} {:>2} {:>2}  {:>8} {:>7} {:>7} {:>7} {:>7} {:>7}",
            "metric", "k", "w", "ceil", "R@.30", "R@.20", "R@.15", "R@.10", "R@.05");
        let (snp_r, ind_r) = (0.025_f64, 0.008_f64);
        for &use_cont in &[false, true] {
            for &k in &[11usize, 13, 15] {
                for &w in &[10usize] {
                    let (_mean, mf) = eval(use_cont, k, w, snp_r, ind_r);
                    let ceil = neg_ceiling(use_cont, k, w);
                    let g = |t: f64| mf.iter().find(|(tt, _)| (*tt - t).abs() < 1e-9).map(|(_, f)| *f).unwrap_or(0.0);
                    eprintln!("{:>5} {:>2} {:>2}  {:>8.4} {:>6.0}% {:>6.0}% {:>6.0}% {:>6.0}% {:>6.0}%",
                        if use_cont { "cont" } else { "jacc" }, k, w, ceil,
                        100.0 * g(0.30), 100.0 * g(0.20), 100.0 * g(0.15), 100.0 * g(0.10), 100.0 * g(0.05));
                }
            }
        }

        eprintln!("\n===== C. recommended operating point: jaccard,k=15,w=10 — recall vs divergence =====");
        let ceil = neg_ceiling(false, 15, 10);
        eprintln!("unrelated ceiling = {:.4}", ceil);
        eprintln!("{:>6} {:>7} {:>7} {:>7} {:>7}", "thr", "1%", "2%", "3%", "5%");
        for &t in &[0.30_f64, 0.20, 0.15, 0.10, 0.07] {
            let r: Vec<f64> = [0.01_f64, 0.02, 0.03, 0.05].iter().map(|&d| {
                let (_m, mf) = eval(false, 15, 10, d, d / 3.0);
                mf.iter().find(|(tt, _)| (*tt - t).abs() < 1e-9).map(|(_, f)| *f).unwrap_or(0.0)
            }).collect();
            let fm = if ceil >= t { "FALSE-MERGE" } else { "clean" };
            eprintln!("{:>6.2} {:>6.0}% {:>6.0}% {:>6.0}% {:>6.0}%   {}", t,
                100.0 * r[0], 100.0 * r[1], 100.0 * r[2], 100.0 * r[3], fm);
        }
        eprintln!();
    }

    #[test]
    #[ignore]
    fn merge_jaccard_sensitivity_sweep() {
        let path = match std::env::var("MERGE_SWEEP_DATA") {
            Ok(p) => p,
            Err(_) => { eprintln!("MERGE_SWEEP_DATA unset; skipping sweep"); return; }
        };
        let raw = std::fs::read_to_string(&path).expect("read MERGE_SWEEP_DATA");
        let exons: Vec<(String, Vec<u8>)> = raw.lines().filter_map(|l| {
            let mut it = l.splitn(2, '\t');
            let id = it.next()?.to_string();
            let seq = it.next()?.trim().as_bytes().to_vec();
            if seq.len() < 15 { return None; }
            Some((id, seq))
        }).collect();

        let fam: Vec<&(String, Vec<u8>)> =
            exons.iter().filter(|(id, _)| id.starts_with("LRPAP1")).collect();
        let unrel: Vec<&(String, Vec<u8>)> =
            exons.iter().filter(|(id, _)| id.starts_with("UNREL")).collect();
        assert!(!fam.is_empty() && !unrel.is_empty(), "need both family and unrel exons");

        const K: usize = 15;
        const W: usize = 10;
        const TRIALS: usize = 40;
        let rates = [0.0_f64, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.10];

        eprintln!("\n===== PART 1: controlled SNP-divergence curve (k={K}, w={W}) =====");
        eprintln!("real exons: {} (LRPAP1), {} trials/rate", fam.len(), TRIALS);
        eprintln!("{:>7} {:>9} {:>8} {:>8} {:>8} {:>10} {:>12}",
            "div%", "meanJacc", "minJacc", "maxJacc", "meanSNP", "merge@.30", "merge@.10");
        // per-rate aggregate of mean Jaccard for the crossing computation
        let mut rate_meanj: Vec<(f64, f64)> = Vec::new();
        for &r in &rates {
            let (mut sj, mut mn, mut mx, mut ssnp, mut merge30, mut merge10, mut cnt) =
                (0.0_f64, 1.0_f64, 0.0_f64, 0.0_f64, 0usize, 0usize, 0usize);
            for (ei, exon) in fam.iter().enumerate() {
                let seq = &exon.1;
                let orig = minimizers(seq, K, W);
                // deterministic seed per (exon, rate)
                let seed = 0xD1B5_4A32_u64
                    .wrapping_mul(ei as u64 + 1)
                    .wrapping_add((r * 1e6) as u64);
                let mut rng = SplitMix64(seed ^ 0x1234_5678_9ABC_DEF0);
                for _ in 0..TRIALS {
                    let (mt, nsnp) = mutate_snps(seq, r, &mut rng);
                    let mm = minimizers(&mt, K, W);
                    let j = jaccard(&orig, &mm);
                    sj += j; mn = mn.min(j); mx = mx.max(j); ssnp += nsnp as f64;
                    if j >= 0.30 { merge30 += 1; }
                    if j >= 0.10 { merge10 += 1; }
                    cnt += 1;
                }
            }
            let meanj = sj / cnt as f64;
            rate_meanj.push((r, meanj));
            eprintln!("{:>7.1} {:>9.3} {:>8.3} {:>8.3} {:>8.1} {:>9.0}% {:>11.0}%",
                r * 100.0, meanj, mn, mx, ssnp / cnt as f64,
                100.0 * merge30 as f64 / cnt as f64,
                100.0 * merge10 as f64 / cnt as f64);
        }

        eprintln!("\n===== PART 2: unrelated-exon baseline (false-merge ceiling) =====");
        // all distinct pairs among unrel, plus cross fam×unrel (truly non-homologous)
        let mut neg: Vec<(String, f64)> = Vec::new();
        let mins_un: Vec<(String, HashSet<u64>)> =
            unrel.iter().map(|(id, s)| (id.clone(), minimizers(s, K, W))).collect();
        let mins_fa: Vec<(String, HashSet<u64>)> =
            fam.iter().map(|(id, s)| (id.clone(), minimizers(s, K, W))).collect();
        for i in 0..mins_un.len() {
            for j in (i + 1)..mins_un.len() {
                neg.push((format!("{}|{}", mins_un[i].0, mins_un[j].0),
                          jaccard(&mins_un[i].1, &mins_un[j].1)));
            }
        }
        for a in &mins_fa {
            for b in &mins_un {
                neg.push((format!("{}|{}", a.0, b.0), jaccard(&a.1, &b.1)));
            }
        }
        let max_neg = neg.iter().map(|x| x.1).fold(0.0_f64, f64::max);
        let mean_neg = neg.iter().map(|x| x.1).sum::<f64>() / neg.len() as f64;
        eprintln!("non-homologous pairs: {}", neg.len());
        eprintln!("  mean Jaccard = {:.4}, MAX Jaccard = {:.4}", mean_neg, max_neg);
        let mut sorted = neg.clone();
        sorted.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        eprintln!("  top-3 non-homologous: {:?}",
            sorted.iter().take(3).map(|(k,v)| format!("{k}={v:.4}")).collect::<Vec<_>>());

        eprintln!("\n===== PART 3: threshold sweep (recall vs false-merge) =====");
        eprintln!("recall = % of mutated-copy trials that still MERGE at each divergence");
        eprintln!("{:>7} {:>9} {:>9} {:>9} {:>9}   {:>10}",
            "thresh", "R@0.5%", "R@1%", "R@2%", "R@3%", "false-merge");
        let report_rates = [0.005_f64, 0.01, 0.02, 0.03];
        // recompute recall per (threshold, rate) using same deterministic mutation
        let recall = |thr: f64, rr: f64| -> f64 {
            let (mut hit, mut cnt) = (0usize, 0usize);
            for (ei, exon) in fam.iter().enumerate() {
                let seq = &exon.1;
                let orig = minimizers(seq, K, W);
                let seed = 0xD1B5_4A32_u64
                    .wrapping_mul(ei as u64 + 1)
                    .wrapping_add((rr * 1e6) as u64);
                let mut rng = SplitMix64(seed ^ 0x1234_5678_9ABC_DEF0);
                for _ in 0..TRIALS {
                    let (mt, _n) = mutate_snps(seq, rr, &mut rng);
                    let j = jaccard(&orig, &minimizers(&mt, K, W));
                    if j >= thr { hit += 1; }
                    cnt += 1;
                }
            }
            100.0 * hit as f64 / cnt as f64
        };
        for &thr in &[0.30_f64, 0.25, 0.20, 0.15, 0.10, 0.07, 0.05, 0.03, 0.02] {
            let r: Vec<f64> = report_rates.iter().map(|&rr| recall(thr, rr)).collect();
            let fm = if max_neg >= thr { "FALSE-MERGE" } else { "clean" };
            eprintln!("{:>7.2} {:>8.0}% {:>8.0}% {:>8.0}% {:>8.0}%   {:>10} (max_neg={:.3})",
                thr, r[0], r[1], r[2], r[3], fm, max_neg);
        }
        eprintln!("\nrecommended bar = just above max non-homologous Jaccard ({:.3}); ", max_neg);
        eprintln!("any threshold in ({:.3}, R_min] keeps unrelated exons separate while", max_neg);
        eprintln!("recovering near-identical paralog exons the 0.30 default rejects.\n");
    }

    #[test]
    fn merge_from_layer1_graphs_shares_homologous_exon_and_has_edges() {
        let g0 = tests_support::make_layer1_graph("chrT", '+', &[(100, 160), (300, 360)]);
        let g1 = tests_support::make_layer1_graph("chrT", '+', &[(100, 160), (5300, 5360)]);
        let genome = tests_support::make_two_copy_genome();
        let fg = build_family_graph_from_layer1_graphs(
            0,
            &[("chrT".to_string(), '+', &g0), ("chrT".to_string(), '+', &g1)],
            Some(&genome),
            0.0, 0.5, 0.0,
        ).expect("merge two layer-1 graphs");
        let shared = fg.nodes.iter().find(|n| n.span == (100, 160)).expect("homologous exon present");
        assert_eq!(shared.per_copy_sequences.len(), 2, "both copies contribute → shared node");
        assert!(!shared.per_copy_sequences[0].1.is_empty(), "sequence channel populated");
        assert!(!shared.copy_specific, "shared exon is not copy-specific");
        let n_private = fg.nodes.iter().filter(|n| n.copy_specific).count();
        assert_eq!(n_private, 2, "one private exon per copy");
        assert!(!fg.edges.is_empty(), "junction edges derived from Layer-1 adjacency");
    }
}

#[cfg(test)]
pub(crate) mod tests_support {
    use super::*;

    /// Minimal real Layer-1 splice graph: source(0), one node per exon (chained), sink.
    pub fn make_layer1_graph(_chrom: &str, _strand: char, exons: &[(u64, u64)]) -> crate::graph::Graph {
        let mut g = crate::graph::Graph::new();
        let source = g.add_node(0, 0).node_id;
        let mut exon_node_ids = Vec::new();
        for &(s, e) in exons {
            exon_node_ids.push(g.add_node(s, e).node_id);
        }
        let sink = g.add_node(u64::MAX, u64::MAX).node_id;
        g.source_id = source;
        g.sink_id = sink;
        g.n_nodes = g.nodes.len();
        let mut chain = vec![source];
        chain.extend(exon_node_ids.iter().copied());
        chain.push(sink);
        for w in chain.windows(2) {
            g.nodes[w[0]].children.insert(w[1]);
            g.nodes[w[1]].parents.insert(w[0]);
        }
        g
    }

    /// Same deterministic genome as make_two_copy_genome (covers 3-exon spans up to 560).
    pub fn make_two_copy_genome_3exon() -> GenomeIndex {
        make_two_copy_genome()
    }

    /// GenomeIndex over chrT (>=60kb) with deterministic pseudo-sequence so shared
    /// coords read identical bytes across copies and distinct spans differ.
    ///
    /// The base at each position is a non-periodic hash of its ABSOLUTE index
    /// (SplitMix64 finalizer), NOT a short cyclic pattern: a low-period pattern
    /// makes every window share the same minimizer set (Jaccard 1.0 everywhere),
    /// which would collapse distinct exons. With per-index hashing, identical
    /// coordinates always read identical bytes (so homologous exons still merge)
    /// while distinct spans get distinct minimizers (so private exons stay split).
    pub fn make_two_copy_genome() -> GenomeIndex {
        use std::io::Write;
        let dir = Box::leak(Box::new(tempfile::tempdir().unwrap()));
        let fa = dir.path().join("chrT.fa");
        let len = 60_000usize;
        let mut seq = vec![b'A'; len];
        for (i, b) in seq.iter_mut().enumerate() {
            let mut z = (i as u64).wrapping_add(0x9E3779B97F4A7C15);
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^= z >> 31;
            *b = match z % 4 { 0 => b'A', 1 => b'C', 2 => b'G', _ => b'T' };
        }
        let mut f = std::fs::File::create(&fa).unwrap();
        writeln!(f, ">chrT").unwrap();
        for chunk in seq.chunks(70) {
            f.write_all(chunk).unwrap();
            f.write_all(b"\n").unwrap();
        }
        drop(f);
        GenomeIndex::from_fasta(fa.to_str().unwrap()).expect("load test genome")
    }
}


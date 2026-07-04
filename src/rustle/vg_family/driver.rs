//! O1 family-definition DRIVER -- native Rust port of `bench/family_rna_refine.py`.
//!
//! CAPSTONE of the O1 migration (`docs/RETIREMENT_AND_MIGRATION.md` Tier-1 driver):
//! the ORCHESTRATION that turns the raw de-novo family graph into the shipped RNA
//! multi-copy family catalog. Every building block it needs is ALREADY ported (and
//! byte-parity tested) in this crate; this module only wires them together in the
//! EXACT order of the Python `build_catalog` / `write_outputs`, reproducing the
//! shipped catalog `bench/family_rna_refine.tsv` byte-for-byte (md5 dca64cbd).
//!
//! REUSED (nothing re-derived here -- the Python source is the spec, line refs are
//! `bench/family_rna_refine.py`):
//!   * `super::family_loaders` : `load_meta` / `load_annot` / `build_genes_dict`
//!     (gene projection) / `load_raw_families` / `load_edges` / `components_from_edges`;
//!   * `super::edge_oracle`    : `load_pair_core` / `load_universal_aln` / `load_allele`
//!     + `demote_gene` (the copy-vs-allele decision) + `gene_pair_key`;
//!   * `super::family_definition` : `refine_families` (gamma-quasi-clique refinement,
//!     GAMMA=0.20, deterministic splitter) + `distinct_loci` (>=2-loci predicate);
//!   * `super::recombinant_split` : `split_block` (recombinant/mosaic split gate);
//!   * `super::multi_repeat_bridge` : `split_families_repeat_bridge` (multi-repeat-bridge
//!     gate, T=8/C=2) + `load_node_mult`;
//!   * `super::bridge_detector` : `load_skeletons` (exon-list r_skel) + `load_strand`;
//!   * `super::repeat_catalog` : `load_skeletons` (introns-string vskel) + `IndexedFasta`.
//!
//! ONE remaining not-yet-standalone input: `load_repeat_mult` LOADS the precomputed
//! per-edge `min_shared_mult` column from the "# SECTION edges" block of the shipped
//! `bench/vg_repeat_catalog.tsv` (exactly as the Python driver LOADS it -- the
//! multiplicity COMPUTATION of that section was deferred in the `repeat_catalog` port;
//! the driver consumes it as a precomputed input, like it consumes precomputed
//! core/aln). See the module report.
//!
//! NOT ported (deliberately out of driver scope; does NOT affect the catalog TSV):
//!   * `validate()` (the DNA diploid-oracle SCORING that only feeds the JSON summary),
//!   * the `family_rna_refine.json` summary + the `_report` stdout block.
//! The acceptance is the CATALOG TSV md5; the JSON is a validation-only side output.

use std::collections::{BTreeSet, HashMap};
use std::io::Write;

use super::bridge_detector::ExonFetcher;
use super::edge_oracle::{self, AlleleSignal};
use super::family_definition::{self, Gene};
use super::family_loaders::{self, GeneMeta};
use super::multi_repeat_bridge::{self, StageLoaders};
use super::recombinant_split;
use super::repeat_catalog::{self, IndexedFasta};

// --------------------------------------------------------------------------- constants
// RECALL-PRESERVING deploy point (family_rna_refine.py:153-164; RNA_ONLY_EDGE_ORACLE.md
// sec 2; do NOT re-derive):
/// `core_recip` KEEP threshold (`family_rna_refine.py:153  CORE_MIN`).
pub const CORE_MIN: f64 = 0.19;
/// `aln_frac` KEEP threshold (`family_rna_refine.py:154  ALN_MIN`).
pub const ALN_MIN: f64 = 0.24;
/// Recall-preserving gamma-quasi-clique cohesion (`family_rna_refine.py:155  GAMMA`
/// == `genome_family_def.GAMMA` == `family_definition::GAMMA`).
pub const GAMMA: f64 = family_definition::GAMMA;
/// HIGH-PRECISION cohesion (`family_rna_refine.py:202  HIGH_PRECISION_GAMMA`;
/// PRECISION_RECALL_FRONTIER.md). Swaps ONLY gamma; everything else unchanged.
pub const HIGH_PRECISION_GAMMA: f64 = 0.40;
/// REPEAT-HUB gate cut (`family_rna_refine.py:164  REPEAT_MULT_MIN`; VG_REPEAT_CATALOG.md
/// tail): a cross-gene edge whose shared sequence is ONLY an extreme repeat
/// (`min_shared_mult >= 20`) is CUT even when core/aln pass.
pub const REPEAT_MULT_MIN: i64 = 20;
/// MULTI-REPEAT-BRIDGE gate T (`family_rna_refine.py:178` == `multi_repeat_bridge::WIRED_MULT_MIN`).
pub const REPEAT_BRIDGE_MULT_MIN: usize = multi_repeat_bridge::WIRED_MULT_MIN;
/// MULTI-REPEAT-BRIDGE gate C (`family_rna_refine.py:179` == `multi_repeat_bridge::WIRED_COUNT_MIN`).
pub const REPEAT_BRIDGE_COUNT_MIN: usize = multi_repeat_bridge::WIRED_COUNT_MIN;
/// allele DEMOTE `balanced_frac >=` conjunct (`family_rna_refine.py:182` ==
/// `edge_oracle::DEMOTE_BAL_MIN`).
pub const DEMOTE_BAL_MIN: f64 = edge_oracle::DEMOTE_BAL_MIN;
/// allele DEMOTE `copy_like <=` conjunct (`family_rna_refine.py:183` ==
/// `edge_oracle::DEMOTE_COPY_MAX`).
pub const DEMOTE_COPY_MAX: f64 = edge_oracle::DEMOTE_COPY_MAX;

// --------------------------------------------------------------------------- inputs
/// The input file paths + gate/gamma switches for one catalog build. Defaults mirror
/// the Python module constants (`family_er_pr` / `rna_only_edge_oracle` /
/// `recombination_bridge_detector` / `vg_repeat_catalog`).
#[derive(Clone, Debug)]
pub struct Inputs {
    /// `family_er_pr.META` = denovo transcript meta (`dn -> chrom/start/end`).
    pub meta: String,
    /// `family_er_pr.ANNOT` = annotation intervals.
    pub annot: String,
    /// `family_er_pr.DENOVO_FAMS` = raw de-novo families.
    pub raw_fams: String,
    /// `family_er_pr.DENOVO_EDGES` == `rna_only_edge_oracle.EDGES_TSV` (edges + core_recip).
    pub edges: String,
    /// `rna_only_edge_oracle.SHAREDLEN_UNIV_TSV` (universal aln_frac).
    pub univ_aln: String,
    /// `rna_only_edge_oracle.A1_TSV` (read-consensus allele signal).
    pub allele: String,
    /// `family_rna_refine.VG_REPEAT_TSV` (per-edge min_shared_mult + node multiplicity).
    pub vg_repeat: String,
    /// `recombination_bridge_detector.SKEL_TSV` (de-novo skeletons; both exon-list + introns-string).
    pub skeletons: String,
    /// `recombination_bridge_detector.GENOME` == `vg_repeat_catalog.GENOME` (soft-masked GGO.fasta).
    pub genome: String,
}

impl Default for Inputs {
    fn default() -> Self {
        const SCRATCH: &str = "/home/juanfra/winloci_scratch";
        const BENCH: &str = "bench";
        Inputs {
            meta: format!("{SCRATCH}/denovo_transcripts.meta.tsv"),
            annot: format!("{SCRATCH}/annot_intervals.tsv"),
            raw_fams: format!("{BENCH}/denovo_families.tsv"),
            edges: format!("{BENCH}/denovo_family_edges.tsv"),
            univ_aln: format!("{BENCH}/ri_sharedlen_universal.tsv"),
            allele: format!("{BENCH}/a1_read_consensus_o1.tsv"),
            vg_repeat: format!("{BENCH}/vg_repeat_catalog.tsv"),
            skeletons: format!("{SCRATCH}/denovo_skeletons.tsv"),
            genome: format!("{SCRATCH}/GGO.fasta"),
        }
    }
}

/// The gate/gamma switches (all gates DEFAULT-ON, matching the Python default).
#[derive(Clone, Copy, Debug)]
pub struct Options {
    pub repeat_gate: bool,
    pub split_recombinants: bool,
    pub repeat_bridge_gate: bool,
    pub gamma: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            repeat_gate: true,
            split_recombinants: true,
            repeat_bridge_gate: true,
            gamma: GAMMA,
        }
    }
}

// --------------------------------------------------------------------------- load_repeat_mult
/// Load the per-edge VG `min_shared_mult`, keyed by gene-pair, from the "# SECTION edges"
/// block of `vg_repeat_catalog.tsv`. Faithful port of `load_repeat_mult`
/// (`family_rna_refine.py:273`): only the `min_shared_mult` column of the per-edge section
/// is read (the soft-mask column is NEVER read); a blank cell is omitted (no repeat cut).
///
/// This is the ONE precomputed input the driver LOADS rather than recomputes -- the
/// per-edge multiplicity COMPUTATION was deferred in the `repeat_catalog` port; the driver
/// consumes it exactly as the Python driver does (like it consumes precomputed core/aln).
pub fn load_repeat_mult(text: &str) -> HashMap<(String, String), i64> {
    let mut out: HashMap<(String, String), i64> = HashMap::new();
    let mut in_edges = false;
    let mut ix: Option<HashMap<String, usize>> = None;
    for ln in text.lines() {
        // Python iterates the whole file; enters the section at "# SECTION edges".
        if ln.starts_with("# SECTION edges") {
            in_edges = true;
            ix = None;
            continue;
        }
        if !in_edges {
            continue;
        }
        if ln.starts_with("gene_a\t") {
            let hdr: Vec<&str> = ln.split('\t').collect();
            let m: HashMap<String, usize> =
                hdr.iter().enumerate().map(|(i, h)| (h.to_string(), i)).collect();
            assert!(
                m.contains_key("min_shared_mult"),
                "vg_repeat_catalog.tsv missing min_shared_mult column"
            );
            ix = Some(m);
            continue;
        }
        let idx = match &ix {
            None => continue,
            Some(m) => m,
        };
        let f: Vec<&str> = ln.split('\t').collect();
        let msm = f[idx["min_shared_mult"]];
        if msm.is_empty() {
            continue;
        }
        let key = edge_oracle::gene_pair_key(f[idx["gene_a"]], f[idx["gene_b"]]);
        out.insert(key, msm.parse().expect("min_shared_mult int"));
    }
    out
}

// --------------------------------------------------------------------------- build result
/// Everything the writer + the report need from one `build_catalog` run.
pub struct BuildResult {
    /// The multi-copy catalog: one block per family, each block a DN-sorted `Vec<String>`.
    pub catalog: Vec<Vec<String>>,
    /// `dn -> gene metadata` (the Python `genes` dict).
    pub genes: HashMap<String, GeneMeta>,
    /// `dn -> gene` (floored projection; the Python `gene_of_dn2`).
    pub gene_of_dn: HashMap<String, String>,
    pub gamma: f64,
    pub repeat_gate: bool,
    pub split_recombinants: bool,
    pub repeat_bridge_gate: bool,
    pub n_families_split: usize,
    pub n_families_repeat_bridge_split: usize,
    pub n_demoted: usize,
    // edge bookkeeping (report only; matches the Python counters).
    pub n_dn_within: usize,
    pub n_dn_cross_kept: usize,
    pub n_dn_cross_cut: usize,
    pub n_dn_cross_cut_repeat: usize,
}

/// Number of distinct physical loci among a DN block (>=2-loci multi-copy predicate).
/// Port of `genome_family_def.distinct_loci(block, genes)` over the block's DN coords.
fn distinct_loci_dn(block: &[String], genes: &HashMap<String, GeneMeta>) -> usize {
    let coords: Vec<(String, i64, i64)> = block
        .iter()
        .map(|dn| {
            let g = &genes[dn];
            (g.chrom.clone(), g.start, g.end)
        })
        .collect();
    family_definition::distinct_loci(&coords)
}

/// `collections.Counter(seq).most_common(1)[0]` -- highest count, ties broken by
/// FIRST-SEEN in iteration order. `None` iff `seq` is empty.
fn most_common<'a>(seq: &[&'a String]) -> Option<(&'a String, usize)> {
    let mut order: Vec<&String> = Vec::new();
    let mut counts: HashMap<&String, usize> = HashMap::new();
    for &g in seq {
        if !counts.contains_key(g) {
            order.push(g);
        }
        *counts.entry(g).or_insert(0) += 1;
    }
    let mut best: Option<(&String, usize)> = None;
    for g in order {
        let c = counts[g];
        match best {
            None => best = Some((g, c)),
            Some((_, bc)) => {
                if c > bc {
                    best = Some((g, c));
                }
            }
        }
    }
    best
}

/// Recombinant-split wrapper (`recombinant_split.split_families`, :210): apply `split_block`
/// to every block, extend the catalog by its parts, count the families that split. The
/// Python `split_info` detail is JSON-only; the driver only needs the split COUNT.
fn split_families_recombinant<F: ExonFetcher>(
    refined: &[Vec<String>],
    genes: &HashMap<String, GeneMeta>,
    gene_of_dn: &HashMap<String, String>,
    skel: &HashMap<(String, i64, i64), Vec<(i64, i64)>>,
    strand: &HashMap<String, String>,
    fa: &F,
) -> (Vec<Vec<String>>, usize) {
    let mut new_refined: Vec<Vec<String>> = Vec::new();
    let mut n_split = 0usize;
    for b in refined {
        let parts = recombinant_split::split_block(b, genes, gene_of_dn, skel, strand, fa);
        if parts.len() > 1 {
            n_split += 1;
        }
        new_refined.extend(parts);
    }
    (new_refined, n_split)
}

// --------------------------------------------------------------------------- build_catalog
/// Apply the recall-preserving RNA-only gate + repeat-hub gate + shipped gamma refinement +
/// recombinant-split gate + multi-repeat-bridge gate + allele demote. Faithful port of
/// `build_catalog` (`family_rna_refine.py:309`); reproduces every step, guard, ordering,
/// and gate on/off wiring. No DNA read.
pub fn build_catalog(inp: &Inputs, opt: Options) -> std::io::Result<BuildResult> {
    // ---- shipped graph context + RNA features (exact ported loaders) ----
    let meta = family_loaders::load_meta(&inp.meta)?;
    let annot = family_loaders::load_annot(&inp.annot)?;
    let raw_fams = family_loaders::load_raw_families(&inp.raw_fams)?;
    let edges_text = std::fs::read_to_string(&inp.edges)?;
    let edge_pairs = family_loaders::load_edges_str(&edges_text);

    let bg = family_loaders::build_genes_dict(&raw_fams, &meta, &annot);
    let genes = bg.genes; // dn -> GeneMeta (the Python `genes`)
    let gene_of_dn = bg.gene_of_dn; // dn -> gene (floored; == RO.load_gene_of_dn())

    // pair_core reads the SAME edges file; gene_of_dn == RO.load_gene_of_dn() (both are the
    // floored build_genes_dict projection), so we reuse it here.
    let pair_core = edge_oracle::load_pair_core_str(&gene_of_dn, &edges_text);
    let univ_aln = edge_oracle::load_universal_aln(&inp.univ_aln)?;
    let allele: HashMap<String, AlleleSignal> = edge_oracle::load_allele(&inp.allele)?;
    let pair_repeat_mult: HashMap<(String, String), i64> = if opt.repeat_gate {
        load_repeat_mult(&std::fs::read_to_string(&inp.vg_repeat)?)
    } else {
        HashMap::new()
    };

    // ---- node universe = union of raw-family members (Python `all_nodes`) ----
    let mut all_nodes: BTreeSet<String> = BTreeSet::new();
    for f in &raw_fams {
        for dn in f {
            all_nodes.insert(dn.clone());
        }
    }

    // ---- alignment KEEP/CUT + repeat-hub gate on unique undirected cross-gene DN edges ----
    let core_aln_keep = |k: &(String, String)| -> bool {
        let c = pair_core.get(k).copied().unwrap_or(0.0);
        let a = univ_aln.get(k).copied().flatten().unwrap_or(0.0); // absent OR empty-cell -> 0.0
        c >= CORE_MIN && a >= ALN_MIN
    };
    let repeat_hub = |k: &(String, String)| -> bool {
        match pair_repeat_mult.get(k) {
            Some(&m) => m >= REPEAT_MULT_MIN,
            None => false, // absent -> no repeat cut (fall through to core+aln)
        }
    };

    // Unique undirected DN edges among nodes in the universe (== networkx Graph.edges()).
    let mut dn_edges: BTreeSet<(String, String)> = BTreeSet::new();
    for (a, b) in &edge_pairs {
        if a != b && all_nodes.contains(a) && all_nodes.contains(b) {
            dn_edges.insert(edge_oracle::gene_pair_key(a, b)); // sorted unordered pair
        }
    }

    let mut kept: Vec<(String, String)> = Vec::new();
    let (mut n_dn_within, mut n_dn_cross_kept, mut n_dn_cross_cut, mut n_dn_cross_cut_repeat) =
        (0usize, 0usize, 0usize, 0usize);
    for (u, v) in &dn_edges {
        let ga = gene_of_dn.get(u);
        let gb = gene_of_dn.get(v);
        match (ga, gb) {
            (Some(ga), Some(gb)) if ga != gb => {
                let k = edge_oracle::gene_pair_key(ga, gb);
                let keep_ca = core_aln_keep(&k);
                let repeat_only = opt.repeat_gate && repeat_hub(&k);
                if keep_ca && !repeat_only {
                    kept.push((u.clone(), v.clone()));
                    n_dn_cross_kept += 1;
                } else {
                    n_dn_cross_cut += 1;
                    if keep_ca && repeat_only {
                        n_dn_cross_cut_repeat += 1;
                    }
                }
            }
            // within-gene / unannotated: never an over-merge, never repeat-gated -> always kept.
            _ => {
                kept.push((u.clone(), v.clone()));
                n_dn_within += 1;
            }
        }
    }

    // ---- components over the kept edges + shipped gamma-quasi-clique refinement ----
    let all_nodes_vec: Vec<String> = all_nodes.iter().cloned().collect();
    let comps = family_loaders::components_from_edges(&all_nodes_vec, &kept);

    // Assign node indices in SORTED-DN order (the universe = all_nodes); this makes each
    // refined block come out ascending-index == DN-sorted (matching the Python `sorted(block)`).
    let dn_to_idx: HashMap<&str, usize> =
        all_nodes_vec.iter().enumerate().map(|(i, s)| (s.as_str(), i)).collect();
    let genes_vec: Vec<Gene> = all_nodes_vec.iter().map(|dn| genes[dn].to_gene()).collect();
    let fams_idx: Vec<Vec<usize>> = comps
        .iter()
        .map(|c| c.iter().map(|dn| dn_to_idx[dn.as_str()]).collect())
        .collect();
    let edges_idx: Vec<(usize, usize)> = kept
        .iter()
        .map(|(u, v)| (dn_to_idx[u.as_str()], dn_to_idx[v.as_str()]))
        .collect();
    let refined_idx = family_definition::refine_families(&fams_idx, &edges_idx, &genes_vec, opt.gamma);
    // Map back to DN. Ascending index == DN-sorted (index order == sorted-DN order).
    let mut refined: Vec<Vec<String>> = refined_idx
        .iter()
        .map(|blk| blk.iter().map(|&i| all_nodes_vec[i].clone()).collect())
        .collect();
    // refine_families already applies the >=2-distinct-loci filter internally (matches Python).

    // Loaders shared by the two split gates (loaded once, lazily, only if needed).
    let need_recs = opt.split_recombinants || opt.repeat_bridge_gate;
    let r_skel = if need_recs {
        Some(super::bridge_detector::load_skeletons(&inp.skeletons)?)
    } else {
        None
    };
    let strand = if need_recs {
        Some(super::bridge_detector::load_strand(&inp.meta)?)
    } else {
        None
    };
    let fa = if need_recs {
        Some(IndexedFasta::open(&inp.genome)?)
    } else {
        None
    };

    // ---- recombinant-split gate (VG path-colinearity; DEFAULT-ON) ----
    let mut n_families_split = 0usize;
    if opt.split_recombinants {
        let (r, n) = split_families_recombinant(
            &refined,
            &genes,
            &gene_of_dn,
            r_skel.as_ref().unwrap(),
            strand.as_ref().unwrap(),
            fa.as_ref().unwrap(),
        );
        n_families_split = n;
        refined = r
            .into_iter()
            .filter(|b| distinct_loci_dn(b, &genes) >= 2)
            .collect();
    }

    // ---- multi-repeat-bridge gate (3rd VG-native gate; DEFAULT-ON) ----
    let mut n_families_repeat_bridge_split = 0usize;
    if opt.repeat_bridge_gate {
        let vskel = repeat_catalog::load_skeletons(&inp.skeletons)?;
        let (mult, cyclic) = multi_repeat_bridge::load_node_mult(&inp.vg_repeat)?;
        let stage = StageLoaders {
            r_skel: r_skel.as_ref().unwrap(),
            strand: strand.as_ref().unwrap(),
            vskel: &vskel,
            meta: &meta,
            fa: fa.as_ref().unwrap(),
            mult: &mult,
            cyclic: &cyclic,
        };
        let (r, info) = multi_repeat_bridge::split_families_repeat_bridge(
            &refined,
            &genes,
            &gene_of_dn,
            REPEAT_BRIDGE_MULT_MIN,
            REPEAT_BRIDGE_COUNT_MIN,
            &stage,
        );
        n_families_repeat_bridge_split = info.len();
        refined = r
            .into_iter()
            .filter(|b| distinct_loci_dn(b, &genes) >= 2)
            .collect();
    }

    // ---- allele DEMOTE (RNA read signal only; exact oracle decision) ----
    let mut catalog: Vec<Vec<String>> = Vec::new();
    let mut n_demoted = 0usize;
    for b in &refined {
        // gs = [gene_of_dn[dn] for dn in b if dn in gene_of_dn]  (iterate b in its order)
        let gs: Vec<&String> = b.iter().filter_map(|dn| gene_of_dn.get(dn)).collect();
        if !gs.is_empty() {
            let (dom, cnt) = most_common(&gs).expect("gs is non-empty");
            let homog = cnt as f64 / gs.len() as f64;
            if edge_oracle::demote_gene(&allele, dom)
                && homog >= 0.5
                && distinct_loci_dn(b, &genes) >= 2
            {
                n_demoted += 1;
                continue; // alleles -> singletons, dropped from the catalog
            }
        }
        let mut sb = b.clone();
        sb.sort(); // catalog.append(sorted(b))
        catalog.push(sb);
    }

    Ok(BuildResult {
        catalog,
        genes,
        gene_of_dn,
        gamma: opt.gamma,
        repeat_gate: opt.repeat_gate,
        split_recombinants: opt.split_recombinants,
        repeat_bridge_gate: opt.repeat_bridge_gate,
        n_families_split,
        n_families_repeat_bridge_split,
        n_demoted,
        n_dn_within,
        n_dn_cross_kept,
        n_dn_cross_cut,
        n_dn_cross_cut_repeat,
    })
}

// --------------------------------------------------------------------------- write_outputs
/// Write the `family_rna_refine.tsv` catalog byte-for-byte. Faithful port of the TSV half
/// of `write_outputs` (`family_rna_refine.py:506`): sort families by their sorted member
/// tuple, then long-format rows `family_id  n_loci  dominant_gene  member_dn  member_gene
/// chrom  start  end`. (`dominant_gene` = the most-common non-NA gene, first-seen tie-break.)
pub fn write_outputs<W: Write>(built: &BuildResult, out: &mut W) -> std::io::Result<()> {
    // fams_sorted = sorted(catalog, key=lambda b: tuple(sorted(b))). Each block is already
    // DN-sorted, so sort the LIST of blocks lexicographically by their member sequence.
    let mut fams_sorted: Vec<&Vec<String>> = built.catalog.iter().collect();
    fams_sorted.sort();

    writeln!(out, "family_id\tn_loci\tdominant_gene\tmember_dn\tmember_gene\tchrom\tstart\tend")?;
    for (fid, b) in fams_sorted.iter().enumerate() {
        // dom = Counter(g for g in [gene_of_dn.get(dn,"NA") for dn in b] if g != "NA").most_common(1)
        let non_na: Vec<&String> = b.iter().filter_map(|dn| built.gene_of_dn.get(dn)).collect();
        let dom: String = most_common(&non_na)
            .map(|(g, _)| g.clone())
            .unwrap_or_else(|| "NA".to_string());
        let nl = distinct_loci_dn(b, &built.genes);
        for dn in b.iter() {
            let g = &built.genes[dn];
            let gene = built.gene_of_dn.get(dn).map(|s| s.as_str()).unwrap_or("NA");
            writeln!(
                out,
                "{fid}\t{nl}\t{dom}\t{dn}\t{gene}\t{}\t{}\t{}",
                g.chrom, g.start, g.end
            )?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `load_repeat_mult` reads ONLY the `min_shared_mult` column of the "# SECTION edges"
    /// block, keyed by the sorted gene-pair, skipping blank cells and all node-section rows.
    #[test]
    fn load_repeat_mult_edges_section_only() {
        let text = "\
# SECTION nodes (mult>=2)
node_id\tmult\tn_tx\tcyclic
12345\t7\t3\t0
# SECTION edges (eval population)
gene_a\tgene_b\tcls\tn_shared_nodes\tmin_shared_mult\tmax_shared_mult
GB\tGA\tgenuine\t3\t38\t99
GC\tGD\ttruthbar\t1\t\t4
GE\tGF\tgenuine\t2\t9\t9";
        let m = load_repeat_mult(text);
        // node-section rows never leak in; blank min_shared_mult omitted.
        assert_eq!(m.len(), 2);
        // key is the SORTED gene-pair (frozenset), value = min_shared_mult int.
        assert_eq!(m.get(&("GA".to_string(), "GB".to_string())).copied(), Some(38));
        assert_eq!(m.get(&("GE".to_string(), "GF".to_string())).copied(), Some(9));
        assert!(!m.contains_key(&("GC".to_string(), "GD".to_string()))); // blank -> omitted
        assert!(!m.contains_key(&("12345".to_string(), "7".to_string()))); // node row ignored
    }

    /// `most_common` = `Counter(seq).most_common(1)[0]`: highest count, ties broken by
    /// FIRST-SEEN insertion order (matches CPython `Counter.most_common(1)`).
    #[test]
    fn most_common_first_seen_tie_break() {
        let s = |x: &str| x.to_string();
        let (a, b, c) = (s("A"), s("B"), s("C"));
        // B appears twice -> winner
        let seq = vec![&a, &b, &c, &b];
        assert_eq!(most_common(&seq), Some((&b, 2)));
        // all singletons -> first-seen (A) wins the tie
        let seq2 = vec![&c, &a, &b];
        assert_eq!(most_common(&seq2), Some((&c, 1)));
        // empty -> None
        let empty: Vec<&String> = Vec::new();
        assert_eq!(most_common(&empty), None);
    }
}

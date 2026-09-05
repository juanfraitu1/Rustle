//! Multi-copy gene families defined by MCL clustering of the ANNOTATION, **by sequence**, then
//! corroborated by RNA.
//!
//! ⭐ **THE CONTRIBUTION IS NOT THE CLUSTERING — IT IS THE CORROBORATION.** Clustering annotated gene
//! sequences is what OrthoFinder does, and `docs/NEGATIVE_RESULTS_REGISTER.md:474` already diagnosed the
//! gap it fills: transitive closure over the same homology graph gives superfamilies of 145 and 114 genes
//! by chaining subfamilies through domain hubs — *"It is OrthoFinder without MCL, and skipping MCL is
//! why."* What this module adds is the second stage: **DNA proposes, RNA disposes.** Measured genome-wide
//! on gorilla (ledger §6de), a repeat clique is a PERFECT clique and a real family is not —
//! median DNA density **1.000** for the 391/1,101 clusters with zero RNA support versus **0.700** for the
//! corroborated ones, **at identical median size 4.0**, so the separation is not the size confound.
//!
//! ⭐ **NO GENE SYMBOL EVER ENTERS THE DEFINITION.** `RFPL` is 9/9 LOC-named and `GOLGA6` 14/14, and
//! genome-wide **95.2%** of members across 1,021 product-defined families are `LOC*`-named with **67.9%**
//! of families entirely so — a `gene=NPIP*` grep recovers **1 of 44** members of the NPIP cluster. Symbols
//! are attached only when REPORTING (see [`Cluster::members`], which carries coordinates, not names).
//!
//! ⚠⚠ **THE COVERAGE DENOMINATOR IS EXONIC, NOT THE GENOMIC SPAN** (ledger §6dc). The median gorilla gene
//! is **23.15% exon**, so a `cov_longer >= 0.30` floor measured against the SPAN demands more coverage
//! than the median gene *has exonic sequence at all*: it removed **1,525/2,428 = 62.8%** of genes that
//! already passed identity and length. Switching the denominator to exonic bases took the pilot graph from
//! **903 to 1,920 nodes with 0 lost**.
//!
//! ⚠ **`cov_longer`, NOT `cov_shorter`.** A ~300 bp Alu covers most of a short fragment and almost none of
//! a real gene, so a shorter-side weight lets shared repeat drive the clustering — every one of the 22
//! adjudicated NPIP false merges was Alu-mediated (§6cr). Weighting by the LONGER side asks how much of
//! BOTH objects is shared, which is what paralogy means and what a shared repeat cannot satisfy.
//!
//! **STATUS:** OTHER-BINARY  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use std::collections::{BTreeMap, BTreeSet};

/// Admission thresholds for an edge of the homology graph. Defaults are the values every measurement in
/// §6da–§6de used; changing one changes the graph, so they are recorded in the params certificate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GraphParams {
    /// Minimum alignment identity (`nmatch / blocklen`).
    pub min_identity: f64,
    /// Minimum fraction of the LONGER sequence covered by the alignment (see the module note).
    pub min_cov_longer: f64,
    /// Minimum alignment block length in bp.
    pub min_bp: u64,
    /// ⭐ Charge `cov_longer`'s NUMERATOR in EXONIC bases too, instead of aligned genomic span.
    ///
    /// The span numerator and the exonic denominator are in DIFFERENT UNITS: the node's sequence is the
    /// whole gene span (introns included), so `max(qe-qs, te-ts)` counts intronic bases and divides them
    /// by an exon-union length. Measured on the pilot: **6,388/24,286 = 26.3%** of admitted edges have
    /// `aln >= denominator`, so the floor is vacuous for them, and **702/4,276 = 16.4%** of within-cluster
    /// edges are driven by an alignment that is <5% exonic. Default OFF ⟹ byte-identical.
    pub exonic_overlap: bool,
    /// ⭐ Reject a pair whose ANNOTATION INTERVALS overlap on the same contig.
    ///
    /// The only self-comparison guard is `q == t` — string equality of the FASTA header — so two distinct
    /// gene records whose intervals overlap align to each other as if they were paralogs. Measured:
    /// **108/4,276** intra-cluster edges join fully NESTED intervals at identity exactly 1.000, and
    /// 54/139 clusters contain at least one. Default OFF ⟹ byte-identical.
    pub reject_overlapping: bool,
    /// ⭐ Minimum ABSOLUTE exonic bases the pair's alignments must jointly cover on the longer gene.
    ///
    /// An ADDITIVE guard, not a replacement for `cov_longer`. ⚠ Replacing the coverage measure with an
    /// exonic one (`exonic_overlap`) over-restricts: a recent segmental duplication copies INTRONS too,
    /// so a real paralog's alignment is legitimately mostly intronic — measured, it shattered the NPIP
    /// cluster 43 -> 19/14/4/4. This clause instead demands that an edge rest on SOME exonic evidence,
    /// which is what the repeat-driven merges lack (MCL2's 33 unrelated genes: <5% exonic over 1.5-2.5 kb
    /// alignments ⟹ <125 exonic bp). 0 = off ⟹ byte-identical.
    ///
    /// ⭐⭐ **SET IT TO 1, NOT 300** (§6dt). The distribution is a WALL AT ZERO, not a gradient:
    /// **53,305/63,361 = 84.1%** of candidate pairs share **exactly zero** exonic bases, only 97 pairs
    /// fall below 32 bp, and moving the floor 1 -> 300 changes the surviving set by **1.6 points**. So
    /// the threshold is not doing the work — the zero/non-zero boundary is — and the clause is properly
    /// a STRUCTURAL requirement with no free number. At 1 bp it also beats 300 on the adjudicated set
    /// (MCL3 27/27 vs 22; the repeat clique 0/33 vs 1/33). ⚠ 300 was anchored to `min_bp` and that
    /// anchor is REFUTED: 300 sits at the EDGE of the 1-200 plateau and costs MCL3 five members.
    pub min_exonic_bp: u64,
}

impl Default for GraphParams {
    fn default() -> Self {
        Self {
            min_identity: 0.70,
            min_cov_longer: 0.30,
            min_bp: 300,
            exonic_overlap: false,
            reject_overlapping: false,
            min_exonic_bp: 0,
        }
    }
}

/// A gene, addressed the way the FASTA headers of an all-vs-all address it.
///
/// ⚠ Coordinates are **GFF 1-based, verbatim** — `NC_073241.2:31346-41669` is exactly GFF `31346 41669`.
/// A `start - 1` key joins **0/4,477** against the annotation and silently falls back to the span
/// denominator, producing a byte-identical graph that reads as "the fix is inert" (§6dd).
/// ⭐ **Always report the join rate.**
pub type GeneKey = (String, u64, u64);

/// One weighted, undirected homology graph over annotated genes.
#[derive(Debug, Default)]
pub struct HomologyGraph {
    /// Node index -> gene, in the order the nodes were first seen (stable for a given PAF).
    pub genes: Vec<GeneKey>,
    /// `(i, j)` with `i < j` -> weight `identity * cov_longer`, capped at 1.0.
    pub edges: BTreeMap<(usize, usize), f64>,
    /// Genes whose exonic length was unknown, so the span was used. Reported, never silent.
    pub missing_exonic: usize,
    /// Edges whose numerator was computed from EXON BLOCKS (`exonic_overlap`). ⚠ **Always report this**:
    /// a zero join is the signature of the §6dd coordinate bug, which produced a byte-identical graph
    /// that read as "the fix is inert".
    pub exonic_overlap_joined: usize,
    /// Edges that reached the numerator but had NO exon blocks, so the span numerator was used.
    pub exonic_overlap_missing: usize,
    /// Pairs dropped by `reject_overlapping`. Reported, never silent.
    pub rejected_overlapping: usize,
    /// Pairs dropped by `min_exonic_bp` — the edge rested on no exonic evidence. Reported, never silent.
    pub rejected_no_exonic: usize,
    /// PAF records between two annotations of ONE locus (`LocusMap`) — a locus aligned to itself —
    /// skipped. Reported, never silent.
    pub same_locus_records: usize,
}

/// ⭐ A LOCUS, not an annotation record, is the node (§6ee).
///
/// Two annotation records whose EXON-UNIONS share ≥1 bp on one contig are the SAME transcribed DNA — a
/// lncRNA model drawn over a gene's exons, two Gnomon models of one transcription unit, an antisense
/// model over the same exons. Left as separate nodes they align to each other at identity 1.000 and each
/// collects the same paralogue edges, so a family counts one locus twice and O2 ties by construction
/// (NPIP: 11 such pairs, 607/1,221 tied reads). Components of that relation are one locus, represented by
/// the record with the greatest exon-union length (ties: lowest start, then lowest end).
/// ⚠ GENOMIC overlap is NOT the criterion: a gene inside another gene's intron is a distinct locus, and
/// merging by genomic overlap folded 7,575 intronic genes genome-wide and dissolved MCL6 (22/23 members
/// sit in one host's introns). Strand is ignored on purpose: for DNA homology an antisense model over the
/// same exons is the same locus.
#[derive(Debug, Default)]
pub struct LocusMap {
    /// Every annotation in a multi-record locus -> its representative (representatives map to themselves).
    pub rep_of: BTreeMap<GeneKey, GeneKey>,
    /// Number of loci that hold more than one annotation record.
    pub n_multi: usize,
}

impl LocusMap {
    pub fn representative<'a>(&'a self, k: &'a GeneKey) -> &'a GeneKey {
        self.rep_of.get(k).unwrap_or(k)
    }
    pub fn is_representative(&self, k: &GeneKey) -> bool {
        self.rep_of.get(k).map_or(true, |r| r == k)
    }
    /// Annotation records folded into another record's locus.
    pub fn n_merged(&self) -> usize {
        self.rep_of.iter().filter(|(k, r)| k != r).count()
    }
}

/// Build the locus map from every gene's merged exon blocks (absolute coordinates, half-open, as
/// `mcl_families::exonic_blocks` emits them). A gene with no blocks is its own locus.
pub fn loci_from_exon_blocks(exon_blocks: &BTreeMap<GeneKey, Vec<(u64, u64)>>) -> LocusMap {
    let genes: Vec<&GeneKey> = exon_blocks.keys().collect();
    let idx: BTreeMap<&GeneKey, usize> = genes.iter().enumerate().map(|(i, g)| (*g, i)).collect();
    let mut parent: Vec<usize> = (0..genes.len()).collect();
    fn find(p: &mut Vec<usize>, mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    // Sweep every exon block per contig; a block starting before the running max end overlaps the
    // block that set it, so the two genes are joined.
    let mut by_contig: BTreeMap<&str, Vec<(u64, u64, usize)>> = BTreeMap::new();
    for (g, blocks) in exon_blocks {
        for &(s, e) in blocks {
            by_contig.entry(g.0.as_str()).or_default().push((s, e, idx[g]));
        }
    }
    for (_, mut v) in by_contig {
        v.sort_unstable();
        let mut max_end = 0u64;
        let mut owner = usize::MAX;
        for (s, e, gi) in v {
            if owner != usize::MAX && s < max_end {
                let (a, b) = (find(&mut parent, gi), find(&mut parent, owner));
                if a != b {
                    parent[a.max(b)] = a.min(b);
                }
            }
            if e > max_end || owner == usize::MAX {
                max_end = max_end.max(e);
                owner = gi;
            }
        }
    }
    let mut comps: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..genes.len() {
        let r = find(&mut parent, i);
        comps.entry(r).or_default().push(i);
    }
    let exlen = |i: usize| exon_blocks[genes[i]].iter().map(|(s, e)| e - s).sum::<u64>();
    let mut m = LocusMap::default();
    for (_, members) in comps {
        if members.len() < 2 {
            continue;
        }
        let rep = *members
            .iter()
            .max_by(|&&a, &&b| exlen(a).cmp(&exlen(b)).then_with(|| (genes[b].1, genes[b].2).cmp(&(genes[a].1, genes[a].2))))
            .unwrap();
        for &i in &members {
            m.rep_of.insert(genes[i].clone(), genes[rep].clone());
        }
        m.n_multi += 1;
    }
    m
}

impl HomologyGraph {
    pub fn n_nodes(&self) -> usize {
        self.genes.len()
    }
    pub fn n_edges(&self) -> usize {
        self.edges.len()
    }
}

/// Parse one gene key out of a FASTA/PAF name of the form `CONTIG:START-END`.
pub fn parse_gene_key(name: &str) -> Option<GeneKey> {
    let (c, r) = name.rsplit_once(':')?;
    let (a, b) = r.split_once('-')?;
    Some((c.to_string(), a.parse().ok()?, b.parse().ok()?))
}

/// Build the weighted graph from an all-vs-all PAF.
///
/// `exonic_len` maps a gene to its EXON-UNION length; a gene absent from the map falls back to the PAF's
/// own sequence length (the genomic span) and is counted in [`HomologyGraph::missing_exonic`].
pub fn graph_from_paf(
    paf: &str,
    exonic_len: &BTreeMap<GeneKey, u64>,
    exon_blocks: &BTreeMap<GeneKey, Vec<(u64, u64)>>,
    p: &GraphParams,
) -> HomologyGraph {
    graph_from_paf_loci(paf, exonic_len, exon_blocks, p, None)
}

/// [`graph_from_paf`] with an optional [`LocusMap`]: when given, the NODE of every record is its locus
/// representative, while the record's own lengths and exon blocks still decide admission — so an edge
/// admitted for ANY annotation of a locus is an edge of the locus (a folded lncRNA model may carry an
/// alignment its host record lacks: NPIP's 4,743-read copy had only such edges). Records between two
/// annotations of one locus are skipped (`same_locus_records`). Coordinates are never remapped.
pub fn graph_from_paf_loci(
    paf: &str,
    exonic_len: &BTreeMap<GeneKey, u64>,
    exon_blocks: &BTreeMap<GeneKey, Vec<(u64, u64)>>,
    p: &GraphParams,
    loci: Option<&LocusMap>,
) -> HomologyGraph {
    let mut g = HomologyGraph::default();
    let mut idx: BTreeMap<GeneKey, usize> = BTreeMap::new();
    let mut missing: BTreeSet<GeneKey> = BTreeSet::new();

    let mut node_of = |k: GeneKey, genes: &mut Vec<GeneKey>| -> usize {
        if let Some(&i) = idx.get(&k) {
            return i;
        }
        let i = genes.len();
        genes.push(k.clone());
        idx.insert(k, i);
        i
    };

    // ⭐ Under `exonic_overlap`, coverage is a property of the PAIR, not of one PAF record: minimap2
    // splits one paralogous alignment across many records, so a per-record floor rejects a real family
    // whose exons are jointly covered. Accumulate the pair's intervals, threshold ONCE at the end.
    // (Measured: per-record thresholding shattered the NPIP cluster 43 -> 17/15/3/3/1.)
    type PairAcc = (Vec<(u64, u64)>, Vec<(u64, u64)>, u64, u64);
    let mut acc: BTreeMap<(GeneKey, GeneKey), PairAcc> = BTreeMap::new();

    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 11 {
            continue;
        }
        let (q, t) = (f[0], f[5]);
        if q == t {
            continue;
        }
        let (Some(qk), Some(tk)) = (parse_gene_key(q), parse_gene_key(t)) else { continue };
        if let Some(l) = loci {
            if l.representative(&qk) == l.representative(&tk) {
                g.same_locus_records += 1;
                continue;
            }
        }
        // Two DISTINCT gene records whose 1-based inclusive intervals intersect are not paralogs; the
        // `q == t` guard above only catches byte-identical headers.
        if p.reject_overlapping && qk.0 == tk.0 && qk.1 <= tk.2 && tk.1 <= qk.2 {
            g.rejected_overlapping += 1;
            continue;
        }
        let (Ok(ql), Ok(qs), Ok(qe)) = (f[1].parse::<u64>(), f[2].parse::<u64>(), f[3].parse::<u64>())
        else {
            continue;
        };
        let (Ok(tl), Ok(ts), Ok(te)) = (f[6].parse::<u64>(), f[7].parse::<u64>(), f[8].parse::<u64>())
        else {
            continue;
        };
        let (Ok(nmatch), Ok(blocklen)) = (f[9].parse::<u64>(), f[10].parse::<u64>()) else { continue };
        if blocklen < p.min_bp {
            continue;
        }
        let identity = nmatch as f64 / blocklen.max(1) as f64;
        if identity < p.min_identity {
            continue;
        }
        let aln = (qe.saturating_sub(qs)).max(te.saturating_sub(ts));
        let dq = match exonic_len.get(&qk) {
            Some(&v) => v,
            None => {
                missing.insert(qk.clone());
                ql
            }
        };
        let dt = match exonic_len.get(&tk) {
            Some(&v) => v,
            None => {
                missing.insert(tk.clone());
                tl
            }
        };
        if p.exonic_overlap || p.min_exonic_bp > 0 {
            // Defer: bank this record's intervals on each side and decide once the pair is complete.
            let (k, qf) = if qk <= tk {
                ((qk.clone(), tk.clone()), true)
            } else {
                ((tk.clone(), qk.clone()), false)
            };
            let e = acc.entry(k).or_insert_with(|| (Vec::new(), Vec::new(), 0, 0));
            if qf {
                e.0.push((qs, qe));
                e.1.push((ts, te));
            } else {
                e.0.push((ts, te));
                e.1.push((qs, qe));
            }
            e.2 += nmatch;
            e.3 += blocklen;
            continue;
        }
        let cov_longer = (aln as f64 / dq.max(dt).max(1) as f64).min(1.0);
        if cov_longer < p.min_cov_longer {
            continue;
        }
        let w = identity * cov_longer;
        let (qn, tn) = match loci {
            Some(l) => (l.representative(&qk).clone(), l.representative(&tk).clone()),
            None => (qk, tk),
        };
        let (a, b) = (node_of(qn, &mut g.genes), node_of(tn, &mut g.genes));
        let key = if a < b { (a, b) } else { (b, a) };
        let e = g.edges.entry(key).or_insert(0.0);
        if w > *e {
            *e = w;
        }
    }
    // Resolve the deferred pairs: merge each side's intervals, charge the LONGER gene's own exonic bases
    // against its own exonic denominator, and threshold once.
    for ((ak, bk), (aiv, biv, nmatch, blocklen)) in acc {
        let da = exonic_len.get(&ak).copied().unwrap_or(0);
        let db = exonic_len.get(&bk).copied().unwrap_or(0);
        let (gk, iv, den) = if da >= db { (&ak, aiv, da) } else { (&bk, biv, db) };
        let Some(blocks) = exon_blocks.get(gk) else {
            g.exonic_overlap_missing += 1;
            continue;
        };
        g.exonic_overlap_joined += 1;
        let merged = merge_intervals(iv);
        let covered =
            merged.iter().map(|&(s0, e0)| exonic_bases_in(blocks, gk.1, s0, e0)).sum::<u64>();
        if covered < p.min_exonic_bp {
            g.rejected_no_exonic += 1;
            continue;
        }
        // `exonic_overlap` REPLACES the numerator; otherwise keep the span numerator (unioned over the
        // pair's records, so a split alignment is not penalised for being split).
        let numer = if p.exonic_overlap {
            covered
        } else {
            merged.iter().map(|&(s0, e0)| e0 - s0).sum::<u64>()
        };
        let cov_longer = (numer as f64 / den.max(1) as f64).min(1.0);
        if cov_longer < p.min_cov_longer {
            continue;
        }
        let identity = nmatch as f64 / blocklen.max(1) as f64;
        let w = identity * cov_longer;
        let (an, bn) = match loci {
            Some(l) => (l.representative(&ak).clone(), l.representative(&bk).clone()),
            None => (ak, bk),
        };
        let (a, b) = (node_of(an, &mut g.genes), node_of(bn, &mut g.genes));
        let key = if a < b { (a, b) } else { (b, a) };
        let e = g.edges.entry(key).or_insert(0.0);
        if w > *e {
            *e = w;
        }
    }
    g.missing_exonic = missing.len();
    g
}

// ---------------------------------------------------------------------------------------------------
// ⭐ DUPLICON-FIRST CORE REFINEMENT (§6eh, pre-registered adj/core/PREREG.md)
//
// A gene family is the set of loci sharing a CORE — the segment linked by segmental-duplication pairs
// to most of the family — not the set of annotation models that happen to align. On NPIP the core is the
// ~23 kb LCR16a duplicon: it sits whole inside every 125–308 kb chimeric model (ABCC1+NPIP, SORL1+NPIP),
// EIF3C carries 808 bp of it, the ABCC1-region records none, and no SMG1P (LCR16u) record is linked to
// even half of NPIP. One constant, "half", used three times; no other number.
// ---------------------------------------------------------------------------------------------------

/// SEDEF pairs (BED, 0-based half-open, cols 1–3 and 4–6), indexed by the contig of EITHER side.
#[derive(Debug, Default)]
pub struct SdPairs {
    /// contig -> (start, end, other_contig, other_start, other_end), sorted by start.
    by_contig: BTreeMap<String, Vec<(u64, u64, String, u64, u64)>>,
    max_len: BTreeMap<String, u64>,
}

impl SdPairs {
    pub fn from_bed_str(text: &str) -> SdPairs {
        let mut sd = SdPairs::default();
        for line in text.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 6 {
                continue;
            }
            let (Ok(s1), Ok(e1), Ok(s2), Ok(e2)) =
                (f[1].parse::<u64>(), f[2].parse::<u64>(), f[4].parse::<u64>(), f[5].parse::<u64>())
            else {
                continue;
            };
            sd.push(f[0], s1, e1, f[3], s2, e2);
            sd.push(f[3], s2, e2, f[0], s1, e1);
        }
        for v in sd.by_contig.values_mut() {
            v.sort_unstable();
        }
        sd
    }
    fn push(&mut self, c: &str, s: u64, e: u64, oc: &str, os: u64, oe: u64) {
        self.by_contig.entry(c.to_string()).or_default().push((s, e, oc.to_string(), os, oe));
        let m = self.max_len.entry(c.to_string()).or_insert(0);
        *m = (*m).max(e.saturating_sub(s));
    }
    pub fn n_pairs(&self) -> usize {
        self.by_contig.values().map(|v| v.len()).sum::<usize>() / 2
    }
    /// Pairs whose THIS side overlaps `[s, e)` on `contig`.
    fn overlapping(&self, contig: &str, s: u64, e: u64) -> impl Iterator<Item = &(u64, u64, String, u64, u64)> {
        let v = self.by_contig.get(contig).map(|v| v.as_slice()).unwrap_or(&[]);
        let lo = s.saturating_sub(self.max_len.get(contig).copied().unwrap_or(0));
        let i = v.partition_point(|x| x.0 < lo);
        v[i..].iter().take_while(move |x| x.0 < e).filter(move |x| x.1 > s)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoreStatus {
    /// The cluster failed the depth gate (no SD evidence): every member is left exactly as it was.
    Untouched,
    /// `core_bp >= span/2`: a full or partial copy of the duplicon.
    KeptFull,
    /// `core_bp < span/2` but `>= median_core/2`: a chimeric model carrying a full core; node := core hull.
    KeptTrimmed,
    /// Neither: a record carrying at most a sliver of the family's core.
    Dropped,
}

#[derive(Debug, Clone)]
pub struct CoreRecord {
    pub member: GeneKey,
    pub gate_passed: bool,
    /// Maximum over the member's bases of the number of distinct other members linked there.
    pub max_depth: usize,
    pub core_bp: u64,
    pub span: u64,
    pub median_core: u64,
    pub status: CoreStatus,
    /// 1-based inclusive hull of the core segments, when a core exists.
    pub hull: Option<(u64, u64)>,
}

/// Apply the pre-registered core rule to one cluster. Members are GFF 1-based inclusive.
pub fn refine_cluster_cores(members: &[GeneKey], sd: &SdPairs) -> Vec<CoreRecord> {
    let n = members.len();
    if n < 2 {
        return members
            .iter()
            .map(|m| CoreRecord {
                member: m.clone(),
                gate_passed: false,
                max_depth: 0,
                core_bp: 0,
                span: m.2.saturating_sub(m.1) + 1,
                median_core: 0,
                status: CoreStatus::Untouched,
                hull: None,
            })
            .collect();
    }
    let half_others = (n - 1) as f64 / 2.0;
    // Per member: depth profile from every SD pair whose one side overlaps the member and whose other
    // side overlaps a DIFFERENT member; depth counts distinct partner members.
    let mut recs: Vec<(usize, u64, Vec<(u64, u64)>)> = Vec::with_capacity(n); // (max_depth, core_bp, core segments)
    for (ri, r) in members.iter().enumerate() {
        let (rs, re) = (r.1.saturating_sub(1), r.2);
        let mut events: Vec<(u64, i32, usize)> = Vec::new();
        for (s, e, oc, os, oe) in sd.overlapping(&r.0, rs, re) {
            for (oi, o) in members.iter().enumerate() {
                if oi == ri || o.0 != *oc {
                    continue;
                }
                if *oe > o.1.saturating_sub(1) && *os < o.2 {
                    events.push((rs.max(*s), 1, oi));
                    events.push((re.min(*e), -1, oi));
                }
            }
        }
        events.sort_unstable();
        let mut count = vec![0i32; n];
        let mut depth = 0usize;
        let mut max_depth = 0usize;
        let mut last: Option<u64> = None;
        let mut segs: Vec<(u64, u64)> = Vec::new();
        for (pos, d, oi) in events {
            if let Some(l) = last {
                if pos > l && depth as f64 >= half_others {
                    if let Some(t) = segs.last_mut().filter(|t| t.1 == l) {
                        t.1 = pos;
                    } else {
                        segs.push((l, pos));
                    }
                }
            }
            if d > 0 {
                if count[oi] == 0 {
                    depth += 1;
                }
                count[oi] += 1;
            } else {
                count[oi] -= 1;
                if count[oi] == 0 {
                    depth -= 1;
                }
            }
            max_depth = max_depth.max(depth);
            last = Some(pos);
        }
        let core_bp = segs.iter().map(|(a, b)| b - a).sum::<u64>();
        recs.push((max_depth, core_bp, segs));
    }
    let mut ratios: Vec<f64> = recs.iter().map(|(md, _, _)| *md as f64 / (n - 1) as f64).collect();
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let gate = median_f(&ratios) >= 0.5;
    let mut cores: Vec<u64> = recs.iter().map(|(_, c, _)| *c).collect();
    cores.sort_unstable();
    let median_core = cores[cores.len() / 2];
    members
        .iter()
        .zip(recs)
        .map(|(m, (max_depth, core_bp, segs))| {
            let span = m.2.saturating_sub(m.1) + 1;
            let hull = if segs.is_empty() { None } else { Some((segs[0].0 + 1, segs.last().unwrap().1)) };
            let status = if !gate {
                CoreStatus::Untouched
            } else if core_bp as f64 >= span as f64 / 2.0 {
                CoreStatus::KeptFull
            } else if core_bp as f64 >= median_core as f64 / 2.0 && core_bp > 0 {
                CoreStatus::KeptTrimmed
            } else {
                CoreStatus::Dropped
            };
            CoreRecord { member: m.clone(), gate_passed: gate, max_depth, core_bp, span, median_core, status, hull }
        })
        .collect()
}

fn median_f(sorted: &[f64]) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    sorted[sorted.len() / 2]
}

/// Sort and merge half-open intervals so overlapping alignment records are not double-counted.
fn merge_intervals(mut v: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    if v.is_empty() {
        return v;
    }
    v.sort_unstable();
    let mut out = vec![v[0]];
    for &(s, e) in &v[1..] {
        let last = out.last_mut().expect("non-empty");
        if s <= last.1 {
            last.1 = last.1.max(e);
        } else {
            out.push((s, e));
        }
    }
    out
}

/// Bases of the local, 0-based half-open interval `[s, e)` that fall inside an exon.
///
/// `blocks` are ABSOLUTE 0-based half-open merged exon intervals; `gene_start` is the GFF **1-based**
/// gene start, so an absolute base `b` sits at local `b - (gene_start - 1)`. ⚠ Getting that conversion
/// wrong is the §6dd bug: a `start - 1` key joined 0/4,477 and silently fell back to the span.
fn exonic_bases_in(blocks: &[(u64, u64)], gene_start: u64, s: u64, e: u64) -> u64 {
    let off = gene_start.saturating_sub(1);
    let mut acc = 0u64;
    for &(bs, be) in blocks {
        let (lo, hi) = (bs.saturating_sub(off).max(s), be.saturating_sub(off).min(e));
        if hi > lo {
            acc += hi - lo;
        }
    }
    acc
}

/// Markov clustering over a sparse symmetric graph.
///
/// ⭐ **INFLATION SITS ON A PLATEAU FOR NPIP — AND ONLY FOR NPIP** (§6dd, caveated §6dx). `I = 2.8..3.2`
/// gives identical NPIP recovery (31/31) — but the real 48-member tandem clique has best-Jaccard **0.00 at
/// I=3.2**, and low-cohesion SD-embedded subfamilies are stable (1.00) at every inflation. Do not generalise
/// the plateau beyond NPIP. Measured on the pilot, `I = 2.8..3.2` gives
/// identical NPIP recovery (31/31), identical principal-cluster size (43) and identical truth
/// concentration (26), with **zero clusters over 100 members**; the cliff is at 3.6 (NPIP 31 -> 18). A
/// constant on a plateau is a choice; a constant on a cliff is indefensible. The default is the middle.
/// ⚠ §6ec (2026-09-04): that cliff is a PRUNE artefact. With the absolute post-inflation `prune` (1e-5) a
/// near-uniform clique of n nodes empties once n+1 > 10^(5/I) (61 at I=2.8, 37 at 3.2); at prune 1e-9
/// NPIP is 44/44 from I=2.0 to 4.0 and the 84+22-copy tandem array, dissolved genome-wide at 1e-5,
/// clusters. The default is kept for byte-identity; `mcl_families --prune` exposes it.
///
/// Determinism: the graph is a `BTreeMap`, columns are normalised in index order, and ties keep the
/// lowest index — the same PAF and parameters always yield the same clustering.
pub fn mcl(g: &HomologyGraph, inflation: f64, max_iter: usize, prune: f64) -> Vec<Vec<usize>> {
    let n = g.n_nodes();
    if n == 0 {
        return Vec::new();
    }
    // Column-major sparse matrix with self-loops, as MCL requires.
    let mut col: Vec<BTreeMap<usize, f64>> = vec![BTreeMap::new(); n];
    for (&(a, b), &w) in &g.edges {
        col[a].insert(b, w);
        col[b].insert(a, w);
    }
    for (i, c) in col.iter_mut().enumerate() {
        c.insert(i, 1.0);
    }
    normalize(&mut col);

    for _ in 0..max_iter {
        let next = expand(&col);
        let mut next = next;
        for c in next.iter_mut() {
            for v in c.values_mut() {
                *v = v.powf(inflation);
            }
            c.retain(|_, v| *v >= prune);
        }
        normalize(&mut next);
        let done = converged(&col, &next);
        col = next;
        if done {
            break;
        }
    }
    attractors(&col, n)
}

fn normalize(col: &mut [BTreeMap<usize, f64>]) {
    for c in col.iter_mut() {
        let s: f64 = c.values().sum();
        if s > 0.0 {
            for v in c.values_mut() {
                *v /= s;
            }
        }
    }
}

/// One MCL expansion step: `M <- M * M`, column-major.
fn expand(col: &[BTreeMap<usize, f64>]) -> Vec<BTreeMap<usize, f64>> {
    col.iter()
        .map(|c| {
            let mut out: BTreeMap<usize, f64> = BTreeMap::new();
            for (&k, &wk) in c {
                for (&i, &wi) in &col[k] {
                    *out.entry(i).or_insert(0.0) += wi * wk;
                }
            }
            out
        })
        .collect()
}

fn converged(a: &[BTreeMap<usize, f64>], b: &[BTreeMap<usize, f64>]) -> bool {
    a.iter().zip(b.iter()).all(|(x, y)| {
        x.len() == y.len()
            && x.iter().zip(y.iter()).all(|((i, u), (j, v))| i == j && (u - v).abs() < 1e-7)
    })
}

/// Read clusters off the converged matrix: every column's heaviest row is its attractor, and columns
/// sharing an attractor (transitively) are one cluster. Ties keep the lowest index.
fn attractors(col: &[BTreeMap<usize, f64>], n: usize) -> Vec<Vec<usize>> {
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut Vec<usize>, mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    for (j, c) in col.iter().enumerate() {
        let Some((&r, _)) = c.iter().max_by(|a, b| {
            a.1.partial_cmp(b.1).unwrap_or(std::cmp::Ordering::Equal).then(b.0.cmp(a.0))
        }) else {
            continue;
        };
        let (ra, rb) = (find(&mut parent, j), find(&mut parent, r));
        if ra != rb {
            parent[ra.max(rb)] = ra.min(rb);
        }
    }
    let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    out.sort_by(|a, b| b.len().cmp(&a.len()).then(a[0].cmp(&b[0])));
    out
}

/// A cluster as reported: member genes plus the two statistics that decide whether it is a family.
#[derive(Debug, Clone)]
pub struct Cluster {
    pub members: Vec<GeneKey>,
    /// Fraction of the possible member pairs that carry a homology edge. ⭐ A repeat clique is a PERFECT
    /// clique (median **1.000**); a real family is not (**0.700**) — §6de, size-controlled.
    pub density: f64,
    /// ⭐ **COHESION CERTIFICATE**: `e_in / (e_in + e_out)` — the share of this cluster's incident edges
    /// that stay inside it. ⚠ Distinct from `density`, which sees only internal pairs and is therefore
    /// blind to a slice cut out of a hairball.
    ///
    /// Calibrated (§6du, admitted edges at `min_exonic_bp = 1`): adjudicated REAL families measure
    /// **0.970** (NPIP, n=43) and **0.923** (n=27); the adjudicated repeat clique measures **0.283** and
    /// the unadjudicated 101-member cluster **0.557**. Over the pilot, **73/142** clusters sit in a clean
    /// mode at >= 0.9, while **20/46** clusters with n >= 5 fall below 0.75, covering 253/620 members.
    /// ⚠ Calibrated on 3 real + 2 artefact clusters — REPORTED, never used as a filter.
    pub frac_in: f64,
    /// Fraction of members with RNA read support, when a BAM was supplied. `None` means NOT MEASURED —
    /// never conflate it with 0.0, which is the repeat-clique signature.
    pub corroborated: Option<f64>,
}

/// Assemble reportable clusters. `min_size` drops singletons (and, at 3, the pairs that carry no
/// density signal). `corroborated_members` is the caller's RNA verdict per gene, or `None` if no BAM.
pub fn build_clusters(
    g: &HomologyGraph,
    parts: &[Vec<usize>],
    min_size: usize,
    corroborated: Option<&dyn Fn(&GeneKey) -> bool>,
) -> Vec<Cluster> {
    parts
        .iter()
        .filter(|p| p.len() >= min_size)
        .map(|p| {
            let n = p.len();
            let possible = n * (n - 1) / 2;
            let mut have = 0usize;
            for (a, x) in p.iter().enumerate() {
                for y in p.iter().skip(a + 1) {
                    let k = if x < y { (*x, *y) } else { (*y, *x) };
                    if g.edges.contains_key(&k) {
                        have += 1;
                    }
                }
            }
            // Cohesion: count every edge incident on a member, then split inside/outside. `density`
            // cannot see e_out, which is exactly how a hairball slice passes as a family.
            let inside: BTreeSet<usize> = p.iter().copied().collect();
            let (mut e_in, mut e_out) = (0usize, 0usize);
            for &(x, y) in g.edges.keys() {
                match (inside.contains(&x), inside.contains(&y)) {
                    (true, true) => e_in += 1,
                    (true, false) | (false, true) => e_out += 1,
                    _ => {}
                }
            }
            let members: Vec<GeneKey> = p.iter().map(|&i| g.genes[i].clone()).collect();
            let corr = corroborated.map(|f| {
                members.iter().filter(|m| f(m)).count() as f64 / n.max(1) as f64
            });
            Cluster {
                members,
                density: if possible == 0 { 0.0 } else { have as f64 / possible as f64 },
                frac_in: if e_in + e_out == 0 { 0.0 } else { e_in as f64 / (e_in + e_out) as f64 },
                corroborated: corr,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn paf_line(q: &str, ql: u64, qs: u64, qe: u64, t: &str, tl: u64, ts: u64, te: u64, nm: u64, bl: u64) -> String {
        format!("{q}\t{ql}\t{qs}\t{qe}\t+\t{t}\t{tl}\t{ts}\t{te}\t{nm}\t{bl}\t60")
    }

    #[test]
    fn gene_key_parses_gff_one_based_coordinates_verbatim() {
        // ⚠ The headers carry GFF 1-based coordinates AS WRITTEN. A `start - 1` key joins nothing and
        // silently falls back to the span denominator (§6dd) — pin the convention.
        assert_eq!(
            parse_gene_key("NC_073241.2:31346-41669"),
            Some(("NC_073241.2".to_string(), 31346, 41669))
        );
        // A contig name containing ':' must still split on the LAST one.
        assert_eq!(parse_gene_key("a:b:10-20"), Some(("a:b".to_string(), 10, 20)));
        assert_eq!(parse_gene_key("no-colon"), None);
    }

    #[test]
    fn coverage_denominator_is_exonic_not_span() {
        // Two genes with 10 kb spans but only 1 kb of exon each, sharing a 600 bp alignment.
        // Against the SPAN: 600/10000 = 0.06, rejected. Against EXON: 600/1000 = 0.60, admitted.
        // This is the §6dc defect that excluded 62.8% of genes.
        let paf = paf_line("c:1-10000", 10_000, 0, 600, "c:20001-30000", 10_000, 0, 600, 570, 600);
        let p = GraphParams::default();

        let empty = BTreeMap::new();
        let span_graph = graph_from_paf(&paf, &empty, &BTreeMap::new(), &p);
        assert_eq!(span_graph.n_edges(), 0, "span denominator must reject this pair");
        assert_eq!(span_graph.missing_exonic, 2, "both genes fell back, and that must be counted");

        let mut exonic = BTreeMap::new();
        exonic.insert(("c".to_string(), 1, 10_000), 1_000);
        exonic.insert(("c".to_string(), 20_001, 30_000), 1_000);
        let exon_graph = graph_from_paf(&paf, &exonic, &BTreeMap::new(), &p);
        assert_eq!(exon_graph.n_edges(), 1, "exonic denominator must admit it");
        assert_eq!(exon_graph.missing_exonic, 0);
    }

    #[test]
    fn weight_uses_the_longer_side_so_a_short_fragment_cannot_drive_clustering() {
        // A 300 bp fragment fully covered by a 300 bp alignment to a 10 kb gene. On the SHORTER side that
        // is coverage 1.0; on the LONGER side 0.03, below the floor. Every one of the 22 adjudicated NPIP
        // false merges was exactly this shape — a shared Alu (§6cr).
        let paf = paf_line("c:1-300", 300, 0, 300, "c:1001-11000", 10_000, 0, 300, 290, 300);
        let mut exonic = BTreeMap::new();
        exonic.insert(("c".to_string(), 1, 300), 300);
        exonic.insert(("c".to_string(), 1001, 11_000), 10_000);
        let g = graph_from_paf(&paf, &exonic, &BTreeMap::new(), &GraphParams::default());
        assert_eq!(g.n_edges(), 0, "a fragment-sized alignment must not become an edge");
    }

    #[test]
    fn mcl_splits_two_cliques_joined_by_a_single_bridge() {
        // The superfamily failure this module exists to prevent: transitive closure chains subfamilies
        // through a hub (register:474 — 145- and 114-gene superfamilies). MCL must cut the bridge.
        let mut g = HomologyGraph::default();
        for i in 0..6 {
            g.genes.push(("c".to_string(), i as u64 * 100 + 1, i as u64 * 100 + 50));
        }
        for (a, b) in [(0, 1), (0, 2), (1, 2), (3, 4), (3, 5), (4, 5)] {
            g.edges.insert((a, b), 0.95);
        }
        g.edges.insert((2, 3), 0.72); // the single weak bridge

        let parts = mcl(&g, 2.8, 100, 1e-5);
        let big: Vec<&Vec<usize>> = parts.iter().filter(|p| p.len() >= 2).collect();
        assert_eq!(big.len(), 2, "expected the two cliques to separate, got {parts:?}");
        assert!(big.iter().all(|p| p.len() == 3), "each clique keeps its three members: {parts:?}");
    }

    #[test]
    fn absolute_prune_empties_large_uniform_cliques_and_a_size_safe_prune_does_not() {
        // §6ec: after inflation every entry of a near-uniform n-clique is ~(1/(n+1))^I, so an absolute
        // prune p empties the whole column once n+1 > p^(-1/I) — 61 nodes at I=2.8, p=1e-5. That is how
        // the anchored 84+22-copy tandem array dissolved genome-wide. A 100-clique must survive at 1e-9.
        let mut g = HomologyGraph::default();
        for i in 0..100 {
            g.genes.push(("c".to_string(), i as u64 * 100 + 1, i as u64 * 100 + 50));
        }
        // Weights vary deterministically in [0.90, 0.99]: a perfectly uniform clique is a numerical
        // knife-edge (every column's self entry ties its neighbours) and is not what a family looks like.
        for a in 0..100 {
            for b in (a + 1)..100 {
                g.edges.insert((a, b), 0.90 + 0.09 * (((a * 7 + b * 13) % 17) as f64 / 16.0));
            }
        }
        let old = mcl(&g, 2.8, 100, 1e-5);
        let largest_old = old.iter().map(|p| p.len()).max().unwrap_or(0);
        assert!(largest_old < 50, "1e-5 must shatter a 100-clique at I=2.8 (largest {largest_old}); if it no longer does, the prune semantics changed and §6ec must be re-measured");
        let new = mcl(&g, 2.8, 100, 1e-9);
        let largest_new = new.iter().map(|p| p.len()).max().unwrap_or(0);
        assert!(largest_new >= 95, "1e-9 must keep the 100-clique (largest {largest_new}): {new:?}");
    }

    #[test]
    fn loci_are_exon_overlap_components_not_genomic_overlap() {
        let k = |s: u64, e: u64| ("c".to_string(), s, e);
        let mut bl: BTreeMap<GeneKey, Vec<(u64, u64)>> = BTreeMap::new();
        bl.insert(k(100, 1000), vec![(100, 200), (900, 1000)]); // host: two exons, big intron
        bl.insert(k(300, 400), vec![(300, 400)]); // INTRONIC gene inside the host: a distinct locus
        bl.insert(k(150, 250), vec![(150, 250)]); // overlaps the host's first exon: same locus
        bl.insert(k(950, 1500), vec![(950, 1050), (1400, 1500)]); // overlaps the host's last exon: same locus, longest exon-union? 200 vs host 200 -> tie -> lowest start = host
        bl.insert(k(1450, 1600), vec![(1450, 1600)]); // overlaps the previous only: chained in
        bl.insert(k(2000, 2500), vec![(2000, 2500)]); // separate
        bl.insert(("d".to_string(), 100, 1000), vec![(100, 200)]); // other contig
        let m = loci_from_exon_blocks(&bl);
        assert_eq!(m.n_multi, 1, "{:?}", m.rep_of);
        assert_eq!(m.n_merged(), 3);
        assert!(m.is_representative(&k(300, 400)), "an intronic gene is its own locus");
        assert!(m.is_representative(&k(2000, 2500)));
        assert!(m.is_representative(&("d".to_string(), 100, 1000)));
        // exon-union lengths: host 200, (150,250) 100, (950,1500) 200, (1450,1600) 150 -> tie host vs (950,1500) -> lowest start wins
        for g in [k(150, 250), k(950, 1500), k(1450, 1600)] {
            assert_eq!(m.representative(&g), &k(100, 1000), "{g:?}");
        }
    }

    #[test]
    fn a_nested_annotation_is_not_a_second_node() {
        // Host H (c:1001-3000) contains lncRNA N (c:1201-2200). H aligns to paralogue P, N to paralogue Q
        // (each at full coverage, so both edges pass the floors). Without the locus map the family counts
        // H's locus twice: 4 nodes, 2 edges. With it N is folded into H and N's edge becomes H's:
        // 3 nodes {H, P, Q}, 2 edges, and N is not a node.
        let mut paf = paf_line("c:1001-3000", 2000, 0, 2000, "c:9001-11000", 2000, 0, 2000, 1900, 2000);
        paf.push('\n');
        paf.push_str(&paf_line("c:1201-2200", 1000, 0, 1000, "c:20001-21000", 1000, 0, 1000, 950, 1000));
        let mut ex = BTreeMap::new();
        ex.insert(("c".to_string(), 1001, 3000), 2000);
        ex.insert(("c".to_string(), 1201, 2200), 1000);
        ex.insert(("c".to_string(), 9001, 11000), 2000);
        ex.insert(("c".to_string(), 20001, 21000), 1000);
        let p = GraphParams { min_bp: 100, ..GraphParams::default() };
        let plain = graph_from_paf(&paf, &ex, &BTreeMap::new(), &p);
        assert_eq!((plain.n_nodes(), plain.n_edges()), (4, 2), "{:?}", plain.genes);
        let mut bl: BTreeMap<GeneKey, Vec<(u64, u64)>> = BTreeMap::new();
        bl.insert(("c".to_string(), 1001, 3000), vec![(1000, 3000)]);
        bl.insert(("c".to_string(), 1201, 2200), vec![(1200, 2200)]); // over the host's exon
        bl.insert(("c".to_string(), 9001, 11000), vec![(9000, 11000)]);
        bl.insert(("c".to_string(), 20001, 21000), vec![(20000, 21000)]);
        let m = loci_from_exon_blocks(&bl);
        assert_eq!(m.n_merged(), 1);
        let g = graph_from_paf_loci(&paf, &ex, &BTreeMap::new(), &p, Some(&m));
        assert_eq!((g.n_nodes(), g.n_edges()), (3, 2), "{:?}", g.genes);
        assert_eq!(g.same_locus_records, 0);
        assert!(g.genes.contains(&("c".to_string(), 1001, 3000)));
        assert!(!g.genes.contains(&("c".to_string(), 1201, 2200)));
        // and a record between H and N (the locus aligned to itself) is skipped, not an edge
        paf.push('\n');
        paf.push_str(&paf_line("c:1201-2200", 1000, 0, 1000, "c:1001-3000", 2000, 200, 1200, 1000, 1000));
        let g2 = graph_from_paf_loci(&paf, &ex, &BTreeMap::new(), &p, Some(&m));
        assert_eq!((g2.n_nodes(), g2.n_edges(), g2.same_locus_records), (3, 2, 1));
    }

    #[test]
    fn core_refinement_keeps_copies_trims_chimeras_drops_slivers_and_respects_the_gate() {
        // Family of 5 on contig c: four 20 kb copies at 100k, 200k, 300k, 400k, all pairwise SD-linked over
        // their whole span; one chimeric 150 kb model at 600k-750k whose first 20 kb is linked to all four;
        // one 40 kb record at 900k with a 500 bp sliver linked to all four (EIF3C-shaped).
        let k = |s: u64, e: u64| ("c".to_string(), s, e);
        let copies = [k(100_001, 120_000), k(200_001, 220_000), k(300_001, 320_000), k(400_001, 420_000)];
        let chimera = k(600_001, 750_000);
        let sliver = k(900_001, 940_000);
        let mut bed = String::new();
        let mut pair = |a: (u64, u64), b: (u64, u64)| bed.push_str(&format!("c\t{}\t{}\tc\t{}\t{}\n", a.0, a.1, b.0, b.1));
        for i in 0..4 {
            for j in (i + 1)..4 {
                pair((copies[i].1 - 1, copies[i].2), (copies[j].1 - 1, copies[j].2));
            }
            pair((copies[i].1 - 1, copies[i].2), (600_000, 620_000)); // the chimera's core
            pair((copies[i].1 - 1, copies[i].2), (900_000, 900_500)); // the sliver
        }
        let sd = SdPairs::from_bed_str(&bed);
        assert_eq!(sd.n_pairs(), 6 + 4 + 4);
        let mut members: Vec<GeneKey> = copies.to_vec();
        members.push(chimera.clone());
        members.push(sliver.clone());
        let recs = refine_cluster_cores(&members, &sd);
        assert!(recs.iter().all(|r| r.gate_passed), "{recs:?}");
        for r in &recs[..4] {
            assert_eq!(r.status, CoreStatus::KeptFull, "{r:?}");
            assert_eq!(r.core_bp, 20_000);
        }
        let ch = &recs[4];
        assert_eq!(ch.status, CoreStatus::KeptTrimmed, "{ch:?}");
        assert_eq!(ch.hull, Some((600_001, 620_000)));
        let sl = &recs[5];
        assert_eq!(sl.status, CoreStatus::Dropped, "{sl:?}");
        assert_eq!(sl.core_bp, 500);
        // An old family with NO SD pairs fails the gate and is untouched.
        let none = refine_cluster_cores(&members, &SdPairs::from_bed_str(""));
        assert!(none.iter().all(|r| !r.gate_passed && r.status == CoreStatus::Untouched));
    }

    #[test]
    fn mcl_is_deterministic() {
        let mut g = HomologyGraph::default();
        for i in 0..8 {
            g.genes.push(("c".to_string(), i as u64 + 1, i as u64 + 2));
        }
        for (a, b) in [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (5, 6), (6, 7)] {
            g.edges.insert((a, b), 0.9);
        }
        let a = mcl(&g, 2.8, 100, 1e-5);
        let b = mcl(&g, 2.8, 100, 1e-5);
        assert_eq!(a, b, "same graph and parameters must give the same clustering");
    }

    #[test]
    fn density_separates_a_perfect_clique_from_a_real_family() {
        // §6de, genome-wide and size-controlled: zero-corroboration clusters have median density 1.000,
        // corroborated ones 0.700, at identical median size 4.0.
        let mut g = HomologyGraph::default();
        for i in 0..8 {
            g.genes.push(("c".to_string(), i as u64 + 1, i as u64 + 2));
        }
        for (a, b) in [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
            g.edges.insert((a, b), 0.9); // a perfect 4-clique
        }
        for (a, b) in [(4, 5), (5, 6), (6, 7), (4, 7)] {
            g.edges.insert((a, b), 0.9); // a 4-cycle: 4 of 6 possible pairs
        }
        let parts = vec![vec![0, 1, 2, 3], vec![4, 5, 6, 7]];
        let cs = build_clusters(&g, &parts, 3, None);
        assert_eq!(cs[0].density, 1.0, "the clique is perfect");
        assert!((cs[1].density - 4.0 / 6.0).abs() < 1e-9, "the family is not");
        assert!(cs[0].corroborated.is_none(), "no BAM => NOT MEASURED, never 0.0");
    }

    #[test]
    fn corroboration_is_none_not_zero_when_unmeasured() {
        // ⚠ `None` (no BAM) and `Some(0.0)` (measured, no support — the repeat-clique signature) are
        // different claims and must never be conflated.
        let mut g = HomologyGraph::default();
        for i in 0..3 {
            g.genes.push(("c".to_string(), i as u64 + 1, i as u64 + 2));
        }
        g.edges.insert((0, 1), 0.9);
        g.edges.insert((1, 2), 0.9);
        g.edges.insert((0, 2), 0.9);
        let parts = vec![vec![0, 1, 2]];
        assert!(build_clusters(&g, &parts, 3, None)[0].corroborated.is_none());
        let never = |_: &GeneKey| false;
        assert_eq!(build_clusters(&g, &parts, 3, Some(&never))[0].corroborated, Some(0.0));
    }

    /// ⭐ RED-BEFORE: the MCL0 mechanism. A 794 bp bridge that is 1.1% exonic saturates `cov_longer`
    /// under the span numerator (794 non-exonic bases / a 679 bp exonic denominator) and is REJECTED
    /// once the numerator is charged in exonic bases too.
    #[test]
    fn exonic_overlap_rejects_a_non_exonic_bridge_the_span_numerator_admits() {
        let q: GeneKey = ("NC_1".to_string(), 1, 5000);
        let t: GeneKey = ("NC_2".to_string(), 1, 5000);
        // The alignment sits at local [2400, 3194) — 794 bp, entirely INTRONIC on both genes.
        let paf = "NC_1:1-5000\t5000\t2400\t3194\t+\tNC_2:1-5000\t5000\t2400\t3194\t790\t794\n";
        let exonic: BTreeMap<GeneKey, u64> =
            [(q.clone(), 679u64), (t.clone(), 679u64)].into_iter().collect();
        // Exons live at absolute 0-based [0, 679) — disjoint from the alignment.
        let blocks: BTreeMap<GeneKey, Vec<(u64, u64)>> =
            [(q.clone(), vec![(0u64, 679u64)]), (t.clone(), vec![(0u64, 679u64)])].into_iter().collect();

        let span = graph_from_paf(paf, &exonic, &blocks, &GraphParams::default());
        assert_eq!(span.n_edges(), 1, "span numerator admits the intronic bridge (794/679 -> capped 1.0)");

        let p = GraphParams { exonic_overlap: true, ..GraphParams::default() };
        let exon = graph_from_paf(paf, &exonic, &blocks, &p);
        assert_eq!(exon.n_edges(), 0, "exonic numerator sees 0 shared exonic bases and rejects it");
        assert_eq!(exon.exonic_overlap_joined, 1, "the join MUST fire - a 0 join is the §6dd bug");
        assert_eq!(exon.exonic_overlap_missing, 0);
    }

    /// A genuinely exonic alignment must still be admitted — the fix is a restriction, not a wall.
    #[test]
    fn exonic_overlap_keeps_a_real_exonic_alignment() {
        let q: GeneKey = ("NC_1".to_string(), 1, 5000);
        let t: GeneKey = ("NC_2".to_string(), 1, 5000);
        let paf = "NC_1:1-5000\t5000\t0\t679\t+\tNC_2:1-5000\t5000\t0\t679\t670\t679\n";
        let exonic: BTreeMap<GeneKey, u64> =
            [(q.clone(), 679u64), (t.clone(), 679u64)].into_iter().collect();
        let blocks: BTreeMap<GeneKey, Vec<(u64, u64)>> =
            [(q.clone(), vec![(0u64, 679u64)]), (t.clone(), vec![(0u64, 679u64)])].into_iter().collect();
        let p = GraphParams { exonic_overlap: true, ..GraphParams::default() };
        let g = graph_from_paf(paf, &exonic, &blocks, &p);
        assert_eq!(g.n_edges(), 1, "fully exonic alignment survives");
        assert_eq!(g.exonic_overlap_joined, 1);
    }

    /// ⚠ The §6dd coordinate trap, pinned: blocks are ABSOLUTE, the gene start is GFF 1-based.
    #[test]
    fn exonic_bases_in_converts_absolute_blocks_to_the_local_frame() {
        // Gene at GFF 1-based 1001..2000; exon at absolute 0-based [1000, 1100) == local [0, 100).
        assert_eq!(exonic_bases_in(&[(1000, 1100)], 1001, 0, 100), 100);
        assert_eq!(exonic_bases_in(&[(1000, 1100)], 1001, 50, 150), 50);
        assert_eq!(exonic_bases_in(&[(1000, 1100)], 1001, 200, 300), 0);
    }

    /// The nesting artefact: 108 intra-cluster edges joined fully nested intervals at identity 1.000.
    #[test]
    fn reject_overlapping_drops_a_nested_annotation_pair() {
        let outer: GeneKey = ("NC_1".to_string(), 1000, 9000);
        let inner: GeneKey = ("NC_1".to_string(), 2000, 3000);
        let paf = "NC_1:1000-9000\t8001\t1000\t2001\t+\tNC_1:2000-3000\t1001\t0\t1001\t1001\t1001\n";
        let exonic: BTreeMap<GeneKey, u64> =
            [(outer.clone(), 1001u64), (inner.clone(), 1001u64)].into_iter().collect();
        let b = BTreeMap::new();
        assert_eq!(graph_from_paf(paf, &exonic, &b, &GraphParams::default()).n_edges(), 1);
        let p = GraphParams { reject_overlapping: true, ..GraphParams::default() };
        let g = graph_from_paf(paf, &exonic, &b, &p);
        assert_eq!(g.n_edges(), 0, "a gene nested inside another is not its paralog");
        assert_eq!(g.rejected_overlapping, 1);
    }

    /// Both new clauses are OFF by default, so every prior measurement stays byte-identical.
    #[test]
    fn the_new_clauses_are_off_by_default() {
        let d = GraphParams::default();
        assert!(!d.exonic_overlap, "flipping this is a THESIS EDIT, not a code edit");
        assert!(!d.reject_overlapping);
    }


    /// ⭐ The shipped remedy: an ADDITIVE exonic floor keeps a real, intron-spanning paralogy edge and
    /// drops one that rests on no exonic sequence. Measured on the pilot: NPIP 43/43 intact, the
    /// 33-gene repeat clique reduced to 1, MCL0 split at its adjudicated cut vertex.
    #[test]
    fn min_exonic_bp_drops_a_non_exonic_edge_and_keeps_an_intron_spanning_one() {
        let q: GeneKey = ("NC_1".to_string(), 1, 20000);
        let t: GeneKey = ("NC_2".to_string(), 1, 20000);
        let exonic: BTreeMap<GeneKey, u64> =
            [(q.clone(), 4000u64), (t.clone(), 4000u64)].into_iter().collect();
        // Exons at absolute [0,2000) and [15000,17000); introns everywhere else.
        let blocks: BTreeMap<GeneKey, Vec<(u64, u64)>> = [
            (q.clone(), vec![(0u64, 2000u64), (15000u64, 17000u64)]),
            (t.clone(), vec![(0u64, 2000u64), (15000u64, 17000u64)]),
        ]
        .into_iter()
        .collect();
        let p = GraphParams { min_exonic_bp: 300, ..GraphParams::default() };

        // A whole-locus duplication: 16 kb spanning introns AND both exons. Must SURVIVE.
        let real = "NC_1:1-20000\t20000\t500\t16500\t+\tNC_2:1-20000\t20000\t500\t16500\t15400\t16000\n";
        let g = graph_from_paf(real, &exonic, &blocks, &p);
        assert_eq!(g.n_edges(), 1, "an intron-spanning real paralogy edge must survive");
        assert_eq!(g.rejected_no_exonic, 0);

        // A purely intronic repeat bridge of the same length. Must DIE.
        let rep = "NC_1:1-20000\t20000\t3000\t13000\t+\tNC_2:1-20000\t20000\t3000\t13000\t8500\t10000\n";
        let g2 = graph_from_paf(rep, &exonic, &blocks, &p);
        assert_eq!(g2.n_edges(), 0, "an edge resting on zero exonic bases must be rejected");
        assert_eq!(g2.rejected_no_exonic, 1, "and the rejection must be COUNTED, never silent");

        // Off by default: the same repeat bridge is admitted on the default path.
        assert_eq!(graph_from_paf(rep, &exonic, &blocks, &GraphParams::default()).n_edges(), 1);
        assert_eq!(GraphParams::default().min_exonic_bp, 0);
    }


    /// ⭐ `frac_in` sees a hairball slice that `density` cannot: a 4-clique cut out of a larger blob has
    /// density 1.000 and cohesion 0.5. Calibrated on the pilot at real 0.923-0.970 vs artefact 0.283.
    #[test]
    fn frac_in_catches_a_slice_that_density_calls_perfect() {
        let mut g = HomologyGraph::default();
        for i in 0..8u64 {
            g.genes.push(("NC_1".to_string(), i * 1000 + 1, i * 1000 + 500));
        }
        // 0..4 is a perfect clique; each of its members also leaks one edge to the outside block.
        for a in 0..4 {
            for b in (a + 1)..4 {
                g.edges.insert((a, b), 1.0);
            }
            g.edges.insert((a, 4 + a), 1.0);
        }
        let c = build_clusters(&g, &[vec![0, 1, 2, 3]], 3, None);
        assert_eq!(c.len(), 1);
        assert!((c[0].density - 1.0).abs() < 1e-9, "density is blind: a perfect internal clique");
        assert!(
            (c[0].frac_in - 0.6).abs() < 1e-9,
            "cohesion sees the leak: 6 internal / (6 internal + 4 leaving) = 0.6, got {}",
            c[0].frac_in
        );
    }

}

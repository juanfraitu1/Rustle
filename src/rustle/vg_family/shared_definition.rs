//! The shared multi-copy family definition (`docs/seeded_family_definition.md` §0★★), instantiated on the
//! de novo homology catalog. OPT-IN: `RUSTLE_SHARED_DEFINITION=1`; unset leaves the catalog byte-identical.
//!
//! Port of the pre-registered Python prototype `bench/denovo_shared_def.py` (ledger §6jz-§6kd, prereg Addenda
//! AB-AF). Construction, in order:
//! 1. NODES — the catalog's read-supported reps, consolidated to gene level: exon chains cut at introns longer
//!    than [`MAX_INTRON`] (P99.9 of annotated GGO introns), pieces under [`MIN_PIECE`] bp dropped, same-strand
//!    pieces with overlapping exons merged (Addendum AB).
//! 2. READ-LOCUS NODES — loci the pipeline dropped: same-strand primary MAPQ >= 1 reads grouped by exon
//!    overlap, exons = bases at read depth >= 2, added when they overlap no existing node (Addendum AC). With
//!    `RUSTLE_SD_READ_LOCUS_SPLIT=1`, a chained read group is first split into sub-loci linked by >= 2 reads
//!    (Addendum AF-3).
//! 3. EDGES — the guided finders on genomic DNA: a node's spliced transcript (`minimap2 -x splice -uf`) hitting
//!    another node's exons at identity >= 0.80 over >= 0.50 of the transcript, or its gene body (`-x asm20`)
//!    chaining onto another node's exons at identity >= 0.80 over >= 0.50 of the shorter body.
//! 4. FAMILIES — triangle-supported leader neighbourhoods (confirmed on a fresh gorilla substrate, §6kd).
//!
//! Every tie-break mirrors the prototype so the two produce the same families on the same inputs (AF-1 gate).
//!
//! **STATUS:** OPT-IN  (docs/MODULE_STATUS.md; reached only when `RUSTLE_SHARED_DEFINITION` is set)

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::io::Write;

use anyhow::{Context, Result};

use crate::genome::GenomeIndex;
use crate::vg_family::denovo_assemble::BamRead;
use crate::vg_family::family_detect::DenovoTranscript;

/// Introns longer than this cut a node's exon chain (P99.9 of 1,092,233 annotated GGO_genomic.gff introns).
pub const MAX_INTRON: u64 = 271_359;
/// Pieces and read loci with fewer exonic bases than this are not nodes.
pub const MIN_PIECE: u64 = 100;
const MIN_ID: f64 = 0.80;
const MIN_COV: f64 = 0.50;
const MIN_LOCUS_READS: usize = 3;
const FLAGS: [&str; 5] = ["-c", "-N", "50", "-p", "0.1"];

fn env_on(name: &str) -> bool {
    matches!(std::env::var(name), Ok(v) if !v.is_empty() && v != "0")
}

/// `RUSTLE_SHARED_DEFINITION=1`: build the homology catalog with the shared definition.
pub fn enabled() -> bool {
    env_on("RUSTLE_SHARED_DEFINITION")
}

/// `RUSTLE_SD_READ_LOCUS_SPLIT=1`: split chained read groups before adding read-locus nodes (Addendum AF-3).
pub fn split_enabled() -> bool {
    env_on("RUSTLE_SD_READ_LOCUS_SPLIT")
}

/// Default junction-support floor for read-isoform widening (§6m0: k = 5 keeps FAMILY R at the annotated
/// ceiling 0.963 and lifts FAMILY F strict 0.450 -> 0.515 on the development substrate; lower k buys more
/// isoforms at the cost of precision, higher k the reverse).
pub const ISOFORM_MIN_READS: u64 = 5;

/// Read-isoform widening is ON by default; `RUSTLE_SD_READ_ISOFORM=0` restores the single-representative node.
pub fn isoform_enabled() -> bool {
    !matches!(std::env::var("RUSTLE_SD_READ_ISOFORM"), Ok(v) if v == "0")
}

/// `RUSTLE_SD_ISOFORM_K`: junction-support floor k, default [`ISOFORM_MIN_READS`].
pub fn isoform_k() -> u64 {
    std::env::var("RUSTLE_SD_ISOFORM_K").ok().and_then(|v| v.parse().ok()).unwrap_or(ISOFORM_MIN_READS)
}

/// One node of the copy graph: a gene-level locus with its exons and the exons of its representative transcript.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SdNode {
    pub chrom: String,
    pub strand: char,
    pub n_reads: u64,
    pub exons: Vec<(u64, u64)>,
    pub rep_exons: Vec<(u64, u64)>,
    /// Spliced query chains for this node: the representative chain, plus every read-supported isoform
    /// admitted by [`widen_with_read_isoforms`] (§6m0). Always non-empty; `[rep_exons]` when widening is off.
    pub tx_chains: Vec<Vec<(u64, u64)>>,
}

impl SdNode {
    pub fn start(&self) -> u64 {
        self.exons.first().map(|e| e.0).unwrap_or(0)
    }
    pub fn end(&self) -> u64 {
        self.exons.last().map(|e| e.1).unwrap_or(0)
    }
}

/// Sort and merge intervals; touching intervals (`s <= last.end`) merge, as the prototype's `merge`.
pub fn merge(mut iv: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    iv.sort();
    let mut out: Vec<(u64, u64)> = Vec::with_capacity(iv.len());
    for (s, e) in iv {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

fn exon_sum(ex: &[(u64, u64)]) -> u64 {
    ex.iter().map(|(s, e)| e - s).sum()
}

struct UnionFind(Vec<usize>);

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind((0..n).collect())
    }
    fn find(&mut self, mut x: usize) -> usize {
        while self.0[x] != x {
            self.0[x] = self.0[self.0[x]];
            x = self.0[x];
        }
        x
    }
    /// `parent[find(a)] = find(b)`, the prototype's orientation.
    fn attach(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        self.0[ra] = rb;
    }
}

/// Connected components of intervals that overlap by >= 1 bp, per key, by the prototype's sweep; groups are
/// returned in order of their first member index.
fn overlap_groups<K: Ord + Clone>(items: &[(K, Vec<(u64, u64)>)]) -> Vec<Vec<usize>> {
    let mut uf = UnionFind::new(items.len());
    let mut by: BTreeMap<K, Vec<(u64, u64, usize)>> = BTreeMap::new();
    for (i, (k, blocks)) in items.iter().enumerate() {
        for &(s, e) in blocks {
            by.entry(k.clone()).or_default().push((s, e, i));
        }
    }
    for iv in by.values_mut() {
        iv.sort();
        let (mut end, mut owner): (i64, Option<usize>) = (-1, None);
        for &(s, e, i) in iv.iter() {
            if let Some(o) = owner {
                if (s as i64) < end {
                    uf.attach(i, o);
                }
            }
            if (e as i64) > end {
                end = e as i64;
                owner = Some(i);
            }
        }
    }
    let mut order: Vec<usize> = Vec::new();
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..items.len() {
        let r = uf.find(i);
        groups.entry(r).or_insert_with(|| {
            order.push(r);
            Vec::new()
        }).push(i);
    }
    order.into_iter().map(|r| groups.remove(&r).unwrap_or_default()).collect()
}

/// Nodes as the catalog emits them: one per rep, exons from its intron chain (the node dump's `exons`).
pub fn nodes_from_reps(reps: &[DenovoTranscript]) -> Vec<SdNode> {
    reps.iter()
        .map(|r| {
            let mut exons = crate::vg_family::catalog_input::exon_blocks(r.start, r.end, &r.introns);
            exons.sort();
            SdNode { chrom: r.chrom.clone(), strand: r.strand, n_reads: r.n_reads as u64, tx_chains: vec![exons.clone()], rep_exons: exons.clone(), exons }
        })
        .collect()
}

/// Addendum AB consolidation to gene-level loci.
pub fn consolidate(ab1: &[SdNode]) -> Vec<SdNode> {
    let mut pieces: Vec<SdNode> = Vec::new();
    for n in ab1 {
        if n.exons.is_empty() {
            continue;
        }
        let mut cur = vec![n.exons[0]];
        for &b in &n.exons[1..] {
            if b.0 > cur.last().unwrap().1 && b.0 - cur.last().unwrap().1 > MAX_INTRON {
                pieces.push(SdNode { exons: cur.clone(), rep_exons: cur.clone(), tx_chains: vec![cur.clone()], ..n.clone() });
                cur = vec![b];
            } else {
                cur.push(b);
            }
        }
        pieces.push(SdNode { rep_exons: cur.clone(), tx_chains: vec![cur.clone()], exons: cur, ..n.clone() });
    }
    pieces.retain(|p| exon_sum(&p.exons) >= MIN_PIECE);
    let items: Vec<((String, char), Vec<(u64, u64)>)> =
        pieces.iter().map(|p| ((p.chrom.clone(), p.strand), p.exons.clone())).collect();
    let mut out: Vec<SdNode> = overlap_groups(&items)
        .into_iter()
        .map(|g| {
            let mut rep = &pieces[g[0]];
            for &i in &g[1..] {
                let p = &pieces[i];
                if (p.n_reads, exon_sum(&p.exons)) > (rep.n_reads, exon_sum(&rep.exons)) {
                    rep = p;
                }
            }
            SdNode {
                chrom: rep.chrom.clone(),
                strand: rep.strand,
                n_reads: g.iter().map(|&i| pieces[i].n_reads).sum(),
                exons: merge(g.iter().flat_map(|&i| pieces[i].exons.iter().copied()).collect()),
                tx_chains: vec![rep.exons.clone()],
                rep_exons: rep.exons.clone(),
            }
        })
        .collect();
    out.sort_by(|a, b| (a.chrom.as_str(), a.start(), a.strand).cmp(&(b.chrom.as_str(), b.start(), b.strand)));
    out
}

/// One read's exon blocks (CIGAR segments between `N`), chromosome and transcript strand.
#[derive(Clone, Debug)]
pub struct ReadBlocks {
    pub chrom: String,
    pub strand: char,
    pub blocks: Vec<(u64, u64)>,
}

/// Primary, MAPQ >= 1 reads (the expression gate's reads); strand = read orientation, flipped by `ts:A:-`.
pub fn read_blocks(reads: &[BamRead]) -> Vec<ReadBlocks> {
    reads
        .iter()
        .filter(|r| !r.is_secondary && !r.is_supplementary && r.mapq >= 1)
        .map(|r| {
            let mut strand = if r.reverse { '-' } else { '+' };
            if r.ts == Some('-') {
                strand = if strand == '-' { '+' } else { '-' };
            }
            let (mut blocks, mut pos, mut cur) = (Vec::new(), r.read.ref_start, r.read.ref_start);
            for &(op, n) in &r.read.cigar {
                match op {
                    'M' | 'D' | '=' | 'X' => pos += n,
                    'N' => {
                        if pos > cur {
                            blocks.push((cur, pos));
                        }
                        pos += n;
                        cur = pos;
                    }
                    _ => {}
                }
            }
            if pos > cur {
                blocks.push((cur, pos));
            }
            ReadBlocks { chrom: r.chrom.clone(), strand, blocks }
        })
        .collect()
}

/// Bases covered by >= 2 reads' blocks, merged (at equal coordinates block ends sort before starts).
pub fn depth2_exons(block_lists: &[&[(u64, u64)]]) -> Vec<(u64, u64)> {
    let mut ev: Vec<(u64, i8)> = Vec::new();
    for blocks in block_lists {
        for &(s, e) in blocks.iter() {
            ev.push((s, 1));
            ev.push((e, -1));
        }
    }
    ev.sort();
    let (mut depth, mut start, mut ex) = (0i64, None::<u64>, Vec::new());
    for (x, d) in ev {
        let prev = depth;
        depth += d as i64;
        if prev < 2 && depth >= 2 {
            start = Some(x);
        } else if prev >= 2 && depth < 2 {
            if let Some(s) = start {
                if x > s {
                    ex.push((s, x));
                }
            }
        }
    }
    merge(ex)
}

/// Addendum AF-3: components of segments linked by >= 2 reads overlapping both, with their supporting read counts,
/// in order of each component's first segment.
pub fn split_linked(segments: &[(u64, u64)], block_lists: &[&[(u64, u64)]]) -> Vec<(Vec<usize>, usize)> {
    let touch: Vec<Vec<usize>> = block_lists
        .iter()
        .map(|blocks| {
            let set: BTreeSet<usize> = segments
                .iter()
                .enumerate()
                .filter(|(_, &(s, e))| blocks.iter().any(|&(b0, b1)| b0 < e && s < b1))
                .map(|(k, _)| k)
                .collect();
            set.into_iter().collect()
        })
        .collect();
    let mut link: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    for t in &touch {
        for i in 0..t.len() {
            for j in (i + 1)..t.len() {
                *link.entry((t[i], t[j])).or_insert(0) += 1;
            }
        }
    }
    let mut uf = UnionFind::new(segments.len());
    for (&(i, j), &n) in &link {
        if n >= 2 {
            uf.attach(j, i);
        }
    }
    let mut comps: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for k in 0..segments.len() {
        let r = uf.find(k);
        comps.entry(r).or_default().push(k);
    }
    let mut out: Vec<Vec<usize>> = comps.into_values().collect();
    out.sort_by_key(|v| v[0]);
    out.into_iter()
        .map(|ks| {
            let sup = touch.iter().filter(|t| t.iter().any(|k| ks.contains(k))).count();
            (ks, sup)
        })
        .collect()
}

/// Exon intervals of a node set, for overlap queries.
pub struct ExonIndex {
    by: HashMap<String, Vec<(u64, u64, usize)>>,
}

impl ExonIndex {
    pub fn new(nodes: &[SdNode]) -> Self {
        let mut by: HashMap<String, Vec<(u64, u64, usize)>> = HashMap::new();
        for (i, n) in nodes.iter().enumerate() {
            for &(s, e) in &n.exons {
                by.entry(n.chrom.clone()).or_default().push((s, e, i));
            }
        }
        for v in by.values_mut() {
            v.sort();
        }
        ExonIndex { by }
    }
    /// Nodes with an exon overlapping `[s, e)` by >= 1 bp.
    pub fn hits(&self, chrom: &str, s: u64, e: u64) -> BTreeSet<usize> {
        let mut out = BTreeSet::new();
        if let Some(v) = self.by.get(chrom) {
            let hi = v.partition_point(|x| x.0 < e);
            for &(x0, x1, i) in &v[..hi] {
                if x1 > s && x0 < e {
                    out.insert(i);
                }
            }
        }
        out
    }
}

/// Addendum AC read-locus nodes (AF-3 split when `split`), merged into `base` and sorted as the prototype.
pub fn with_read_locus_nodes(base: &[SdNode], reads: &[ReadBlocks], split: bool) -> (Vec<SdNode>, usize) {
    let items: Vec<((String, char), Vec<(u64, u64)>)> =
        reads.iter().map(|r| ((r.chrom.clone(), r.strand), r.blocks.clone())).collect();
    let bidx = ExonIndex::new(base);
    let mut added: Vec<SdNode> = Vec::new();
    for members in overlap_groups(&items) {
        if members.len() < MIN_LOCUS_READS {
            continue;
        }
        let (chrom, strand) = (reads[members[0]].chrom.clone(), reads[members[0]].strand);
        let lists: Vec<&[(u64, u64)]> = members.iter().map(|&i| reads[i].blocks.as_slice()).collect();
        let ex = depth2_exons(&lists);
        let consider = |sub: Vec<(u64, u64)>, n: usize, added: &mut Vec<SdNode>| {
            if n < MIN_LOCUS_READS || exon_sum(&sub) < MIN_PIECE {
                return;
            }
            if sub.iter().any(|&(s, e)| !bidx.hits(&chrom, s, e).is_empty()) {
                return;
            }
            added.push(SdNode { chrom: chrom.clone(), strand, n_reads: n as u64, tx_chains: vec![sub.clone()], rep_exons: sub.clone(), exons: sub });
        };
        if split {
            for (ks, sup) in split_linked(&ex, &lists) {
                consider(ks.iter().map(|&k| ex[k]).collect(), sup, &mut added);
            }
        } else {
            consider(ex, members.len(), &mut added);
        }
    }
    let n_added = added.len();
    let mut nodes: Vec<SdNode> = base.to_vec();
    nodes.extend(added);
    nodes.sort_by(|a, b| (a.chrom.as_str(), a.start(), a.strand).cmp(&(b.chrom.as_str(), b.start(), b.strand)));
    (nodes, n_added)
}

/// Read-isoform widening (§6m0, ledger 2026-09-17). For each node, every intron chain observed in its own
/// same-strand primary reads whose junctions are EACH carried by >= `k` reads becomes a spliced query, and its
/// blocks are merged into the node's exon union. A node can only widen: its previous exons and representative
/// chain are always kept, and no node is created, removed or merged here.
///
/// Why: the catalog emits one representative per locus before [`consolidate`] runs, so a node's exon set is one
/// isoform (705/710 nodes on the human ideal substrate). Widening lifts full-length NPIP copies from 5/27 to
/// 17/27 at k = 5 with FAMILY R held at the annotated arm's own 0.963.
pub fn widen_with_read_isoforms(nodes: &[SdNode], reads: &[ReadBlocks], k: u64) -> (Vec<SdNode>, usize) {
    let idx = ExonIndex::new(nodes);
    // Per node: chain (intron vector) -> (support, min start, max end).
    let mut per_node: Vec<BTreeMap<Vec<(u64, u64)>, (u64, u64, u64)>> =
        vec![BTreeMap::new(); nodes.len()];
    for r in reads {
        if r.blocks.len() < 2 {
            continue;
        }
        let introns: Vec<(u64, u64)> = r.blocks.windows(2).map(|w| (w[0].1, w[1].0)).collect();
        let (first, last) = (r.blocks[0].0, r.blocks[r.blocks.len() - 1].1);
        let mut seen: BTreeSet<usize> = BTreeSet::new();
        for &(s, e) in &r.blocks {
            for v in idx.hits(&r.chrom, s, e) {
                if nodes[v].strand == r.strand {
                    seen.insert(v);
                }
            }
        }
        for v in seen {
            let ent = per_node[v].entry(introns.clone()).or_insert((0, u64::MAX, 0));
            ent.0 += 1;
            ent.1 = ent.1.min(first);
            ent.2 = ent.2.max(last);
        }
    }
    let mut widened = 0usize;
    let out: Vec<SdNode> = nodes
        .iter()
        .enumerate()
        .map(|(v, n)| {
            // A junction is supported when >= k reads of this node carry it, counted over all of the node's chains.
            let mut junc: BTreeMap<(u64, u64), u64> = BTreeMap::new();
            for (chain, &(sup, _, _)) in &per_node[v] {
                for &j in chain {
                    *junc.entry(j).or_insert(0) += sup;
                }
            }
            let mut chains: Vec<Vec<(u64, u64)>> = n.tx_chains.clone();
            let mut exons = n.exons.clone();
            let mut added = false;
            for (chain, &(_, first, last)) in &per_node[v] {
                if chain.iter().any(|j| junc.get(j).copied().unwrap_or(0) < k) {
                    continue;
                }
                let blocks = crate::vg_family::catalog_input::exon_blocks(first, last, chain);
                if blocks.is_empty() || chains.iter().any(|c| *c == blocks) {
                    continue;
                }
                exons.extend(blocks.iter().copied());
                chains.push(blocks);
                added = true;
            }
            if !added {
                return n.clone();
            }
            widened += 1;
            chains.sort();
            chains.dedup();
            SdNode { exons: merge(exons), tx_chains: chains, ..n.clone() }
        })
        .collect();
    (out, widened)
}

/// One PAF record (the fields the finders use).
#[derive(Clone, Debug)]
pub struct PafRec {
    pub q: String,
    pub qlen: u64,
    pub qs: u64,
    pub qe: u64,
    pub strand: char,
    pub chrom: String,
    pub clen: u64,
    pub ts: u64,
    pub te: u64,
    pub nm: u64,
    pub bl: u64,
    pub cg: String,
}

pub fn parse_paf(text: &str) -> Vec<PafRec> {
    text.lines()
        .filter_map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 12 {
                return None;
            }
            Some(PafRec {
                q: f[0].to_string(),
                qlen: f[1].parse().ok()?,
                qs: f[2].parse().ok()?,
                qe: f[3].parse().ok()?,
                strand: f[4].chars().next()?,
                chrom: f[5].to_string(),
                clen: f[6].parse().ok()?,
                ts: f[7].parse().ok()?,
                te: f[8].parse().ok()?,
                nm: f[9].parse().ok()?,
                bl: f[10].parse().ok()?,
                cg: f[12..].iter().find_map(|x| x.strip_prefix("cg:Z:")).unwrap_or("").to_string(),
            })
        })
        .collect()
}

fn cigar_ops(cg: &str) -> Vec<(u64, char)> {
    let mut out = Vec::new();
    let mut n: u64 = 0;
    for c in cg.chars() {
        if let Some(d) = c.to_digit(10) {
            n = n * 10 + d as u64;
        } else {
            out.push((n, c));
            n = 0;
        }
    }
    out
}

/// Spliced hits passing identity >= 0.80 (matches / block) and query coverage >= 0.50.
pub fn transcript_hits(recs: &[PafRec]) -> Vec<PafRec> {
    recs.iter()
        .filter(|h| h.bl > 0 && h.qlen > 0 && h.nm as f64 / h.bl as f64 >= MIN_ID && (h.qe - h.qs) as f64 / h.qlen as f64 >= MIN_COV)
        .cloned()
        .collect()
}

/// Genomic exon blocks of a spliced alignment (split at `N`).
pub fn tx_exon_blocks(h: &PafRec) -> Vec<(u64, u64)> {
    let (mut blocks, mut pos, mut cur) = (Vec::new(), h.ts, h.ts);
    for (n, op) in cigar_ops(&h.cg) {
        match op {
            'M' | '=' | 'X' | 'D' => pos += n,
            'N' => {
                if pos > cur {
                    blocks.push((cur, pos));
                }
                pos += n;
                cur = pos;
            }
            _ => {}
        }
    }
    if pos > cur {
        blocks.push((cur, pos));
    }
    blocks
}

/// A gene-body chain: query, chromosome, strand and its target interval.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Chain {
    pub q: String,
    pub chrom: String,
    pub strand: char,
    pub s: u64,
    pub e: u64,
}

/// The prototype's `gene_body_chains`: records of one (query, chromosome, strand) chained in target order while
/// the gap stays <= the query length, the span <= twice it and the query advances; a chain is kept at identity
/// >= 0.80 with aligned query >= 0.50 of min(query length, extrapolated target span).
pub fn gene_body_chains(recs: &[PafRec]) -> Vec<Chain> {
    let mut order: Vec<(String, String, char)> = Vec::new();
    let mut by: HashMap<(String, String, char), Vec<&PafRec>> = HashMap::new();
    for r in recs {
        let k = (r.q.clone(), r.chrom.clone(), r.strand);
        by.entry(k.clone()).or_insert_with(|| {
            order.push(k);
            Vec::new()
        }).push(r);
    }
    let mut chains = Vec::new();
    for k in order {
        let mut rs = by.remove(&k).unwrap_or_default();
        rs.sort_by_key(|r| (r.ts, r.te));
        let lq = rs[0].qlen as i64;
        let strand = k.2;
        let mut cur: Vec<&PafRec> = Vec::new();
        let emit = |cur: &Vec<&PafRec>, chains: &mut Vec<Chain>| {
            let nm: u64 = cur.iter().map(|x| x.nm).sum();
            let bl: u64 = cur.iter().map(|x| x.bl).sum();
            let qiv = merge(cur.iter().map(|x| (x.qs, x.qe)).collect());
            let aligned: i64 = qiv.iter().map(|(s, e)| (e - s) as i64).sum();
            let ts = cur.iter().map(|x| x.ts).min().unwrap() as i64;
            let te = cur.iter().map(|x| x.te).max().unwrap() as i64;
            let (q0, q1) = (qiv[0].0 as i64, qiv.last().unwrap().1 as i64);
            let (xs, xe) = if strand == '+' { (ts - q0, te + (lq - q1)) } else { (ts - (lq - q1), te + q0) };
            let (xs, xe) = (xs.max(0), xe.min(cur[0].clen as i64));
            if bl > 0 && nm as f64 / bl as f64 >= MIN_ID && aligned as f64 >= MIN_COV * (lq.min(xe - xs)) as f64 {
                chains.push(Chain { q: k.0.clone(), chrom: k.1.clone(), strand, s: ts as u64, e: te as u64 });
            }
        };
        for r in rs {
            if !cur.is_empty() {
                let gap = r.ts as i64 - cur.iter().map(|x| x.te).max().unwrap() as i64;
                let span = r.te as i64 - cur.iter().map(|x| x.ts).min().unwrap() as i64;
                let last_qs = cur.last().unwrap().qs;
                let ordered = if strand == '+' { r.qs >= last_qs } else { r.qs <= last_qs };
                if !(gap <= lq && span <= 2 * lq && ordered) {
                    emit(&cur, &mut chains);
                    cur.clear();
                }
            }
            cur.push(r);
        }
        if !cur.is_empty() {
            emit(&cur, &mut chains);
        }
    }
    chains
}

/// Query keys shared by nodes with the same representative transcript / gene body (the prototype's dedupe keys).
pub fn tx_key(n: &SdNode) -> String {
    tx_key_for(n, &n.rep_exons)
}

/// The query key of one spliced chain of a node (the same form as [`tx_key`], which is this for `rep_exons`).
pub fn tx_key_for(n: &SdNode, chain: &[(u64, u64)]) -> String {
    let ex: Vec<String> = chain.iter().map(|(a, b)| format!("{a}-{b}")).collect();
    format!("{}|{}|{}", n.chrom, n.strand, ex.join(","))
}

pub fn body_key(n: &SdNode) -> String {
    format!("{}|{}|{}", n.chrom, n.start(), n.end())
}

/// Symmetric node pairs joined by an exon edge or a gene-body edge, with counts per kind.
pub fn edges(
    nodes: &[SdNode],
    tx_by_key: &HashMap<String, Vec<PafRec>>,
    chains_by_key: &HashMap<String, Vec<Chain>>,
) -> (BTreeSet<(usize, usize)>, usize, usize) {
    let idx = ExonIndex::new(nodes);
    let spliced: Vec<bool> = nodes.iter().map(|n| n.exons.len() >= 2).collect();
    let mut kinds: BTreeSet<((usize, usize), u8)> = BTreeSet::new();
    for (u, nu) in nodes.iter().enumerate() {
        let (us, ue) = (nu.start(), nu.end());
        // (kind, chrom, s, e, hit_s, hit_e, orientation on the target)
        let mut cand: Vec<(u8, &str, u64, u64, u64, u64, char)> = Vec::new();
        for chain in &nu.tx_chains {
            for h in tx_by_key.get(&tx_key_for(nu, chain)).map(|v| v.as_slice()).unwrap_or(&[]) {
                for (s, e) in tx_exon_blocks(h) {
                    cand.push((0, &h.chrom, s, e, h.ts, h.te, h.strand));
                }
            }
        }
        for c in chains_by_key.get(&body_key(nu)).map(|v| v.as_slice()).unwrap_or(&[]) {
            let orient = if c.strand == '+' { nu.strand } else if nu.strand == '+' { '-' } else if nu.strand == '-' { '+' } else { nu.strand };
            cand.push((1, &c.chrom, c.s, c.e, c.s, c.e, orient));
        }
        for (kind, chrom, s, e, hs, he, orient) in cand {
            if chrom == nu.chrom && hs < ue && us < he {
                continue;
            }
            for v in idx.hits(chrom, s, e) {
                if v == u || (spliced[u] && spliced[v] && nodes[v].strand != orient) {
                    continue;
                }
                kinds.insert(((u.min(v), u.max(v)), kind));
            }
        }
    }
    let n_exon = kinds.iter().filter(|(_, k)| *k == 0).count();
    let n_body = kinds.len() - n_exon;
    (kinds.into_iter().map(|(p, _)| p).collect(), n_exon, n_body)
}

/// Triangle-supported leader neighbourhoods: nodes in order of reads (desc), degree (desc), chromosome, start and
/// index; an unassigned node with unassigned neighbours leads a family of itself, those neighbours, and unassigned
/// nodes adjacent to >= 2 members of that star.
pub fn triangle_leaders(nodes: &[SdNode], pairs: &BTreeSet<(usize, usize)>) -> Vec<Vec<usize>> {
    let mut adj: Vec<BTreeSet<usize>> = vec![BTreeSet::new(); nodes.len()];
    for &(a, b) in pairs {
        adj[a].insert(b);
        adj[b].insert(a);
    }
    let mut order: Vec<usize> = (0..nodes.len()).collect();
    order.sort_by(|&a, &b| {
        (std::cmp::Reverse(nodes[a].n_reads), std::cmp::Reverse(adj[a].len()), nodes[a].chrom.as_str(), nodes[a].start(), a)
            .cmp(&(std::cmp::Reverse(nodes[b].n_reads), std::cmp::Reverse(adj[b].len()), nodes[b].chrom.as_str(), nodes[b].start(), b))
    });
    let mut seen = vec![false; nodes.len()];
    let mut fams = Vec::new();
    for i in order {
        if seen[i] {
            continue;
        }
        let free: Vec<usize> = adj[i].iter().copied().filter(|&y| !seen[y]).collect();
        if free.is_empty() {
            continue;
        }
        let mut star: BTreeSet<usize> = free.iter().copied().collect();
        star.insert(i);
        let second: BTreeSet<usize> = star
            .iter()
            .flat_map(|&y| adj[y].iter().copied())
            .filter(|&z| !seen[z] && !star.contains(&z) && adj[z].intersection(&star).count() >= 2)
            .collect();
        let fam: BTreeSet<usize> = star.union(&second).copied().collect();
        for &x in &fam {
            seen[x] = true;
        }
        fams.push(fam.into_iter().collect::<Vec<usize>>());
    }
    fams.retain(|f| f.len() >= 2);
    fams.sort_by(|a, b| (std::cmp::Reverse(a.len()), a[0]).cmp(&(std::cmp::Reverse(b.len()), b[0])));
    fams
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

fn spliced_seq(genome: &GenomeIndex, chrom: &str, exons: &[(u64, u64)], strand: char) -> Vec<u8> {
    let mut s: Vec<u8> = Vec::new();
    for &(a, b) in exons {
        if let Some(x) = genome.fetch_sequence(chrom, a, b) {
            s.extend(x.iter().map(|c| c.to_ascii_uppercase()));
        }
    }
    if strand == '-' { revcomp(&s) } else { s }
}

fn run_minimap2(minimap2: &str, preset: &[&str], threads: usize, target: &std::path::Path, query: &std::path::Path) -> Result<String> {
    let out = std::process::Command::new(minimap2)
        .args(FLAGS)
        .args(preset)
        .arg("-t")
        .arg(threads.max(1).to_string())
        .arg(target)
        .arg(query)
        .output()
        .with_context(|| format!("running {minimap2}"))?;
    anyhow::ensure!(out.status.success(), "minimap2 {:?} failed: {}", preset, String::from_utf8_lossy(&out.stderr));
    Ok(String::from_utf8_lossy(&out.stdout).into_owned())
}

/// The full construction on a catalog's reps and reads. Returns the nodes, the families (node indices) and the
/// edge pairs.
pub fn build(
    reps: &[DenovoTranscript],
    reads: &[ReadBlocks],
    genome: &GenomeIndex,
    minimap2: &str,
    threads: usize,
) -> Result<(Vec<SdNode>, Vec<Vec<usize>>, BTreeSet<(usize, usize)>)> {
    let ab2 = consolidate(&nodes_from_reps(reps));
    let (nodes, n_added) = with_read_locus_nodes(&ab2, reads, split_enabled());
    eprintln!(
        "[shared-definition] {} reps -> {} gene-level loci + {} read-locus nodes{} = {} nodes",
        reps.len(), ab2.len(), n_added, if split_enabled() { " (split)" } else { "" }, nodes.len()
    );
    let nodes = if isoform_enabled() {
        let k = isoform_k();
        let (w, n_widened) = widen_with_read_isoforms(&nodes, reads, k);
        let queries: usize = w.iter().map(|n| n.tx_chains.len()).sum();
        eprintln!(
            "[shared-definition] read-isoform widening k={k}: {n_widened} of {} nodes widened, {queries} spliced queries",
            w.len()
        );
        w
    } else {
        eprintln!("[shared-definition] read-isoform widening OFF (RUSTLE_SD_READ_ISOFORM=0)");
        nodes
    };
    let dir = std::env::temp_dir().join(format!("rustle_sd_{}_{}", std::process::id(), reps.len()));
    std::fs::create_dir_all(&dir)?;
    let (target, txfa, bodyfa) = (dir.join("target.fa"), dir.join("tx.fa"), dir.join("body.fa"));
    {
        // Target contig order: `RUSTLE_SD_TARGET_ORDER` (comma list) when set, otherwise sorted names.
        let mut names: Vec<String> = genome.chroms().map(|(c, _)| c.to_string()).collect();
        names.sort();
        if let Ok(v) = std::env::var("RUSTLE_SD_TARGET_ORDER") {
            let want: Vec<String> = v.split(',').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect();
            let mut rest: Vec<String> = names.iter().filter(|n| !want.contains(n)).cloned().collect();
            names = want.into_iter().filter(|w| genome.chrom_len(w) > 0).collect();
            names.append(&mut rest);
        }
        let mut fh = std::io::BufWriter::new(std::fs::File::create(&target)?);
        for c in &names {
            if let Some(seq) = genome.fetch_sequence(c, 0, genome.chrom_len(c)) {
                writeln!(fh, ">{c}")?;
                fh.write_all(&seq)?;
                writeln!(fh)?;
            }
        }
    }
    let (mut tx_seen, mut body_seen) = (BTreeSet::new(), BTreeSet::new());
    {
        let mut ft = std::io::BufWriter::new(std::fs::File::create(&txfa)?);
        let mut fb = std::io::BufWriter::new(std::fs::File::create(&bodyfa)?);
        for n in &nodes {
            let bk = body_key(n);
            for chain in &n.tx_chains {
                let tk = tx_key_for(n, chain);
                if tx_seen.insert(tk.clone()) {
                    writeln!(ft, ">{tk}")?;
                    ft.write_all(&spliced_seq(genome, &n.chrom, chain, n.strand))?;
                    writeln!(ft)?;
                }
            }
            if body_seen.insert(bk.clone()) {
                writeln!(fb, ">{bk}")?;
                let body = genome.fetch_sequence(&n.chrom, n.start(), n.end()).unwrap_or_default();
                fb.write_all(&body.iter().map(|c| c.to_ascii_uppercase()).collect::<Vec<u8>>())?;
                writeln!(fb)?;
            }
        }
    }
    let tx_paf = run_minimap2(minimap2, &["-x", "splice", "-uf"], threads, &target, &txfa)?;
    let body_paf = run_minimap2(minimap2, &["-x", "asm20"], threads, &target, &bodyfa)?;
    let _ = std::fs::remove_dir_all(&dir);
    let mut tx_by_key: HashMap<String, Vec<PafRec>> = HashMap::new();
    for h in transcript_hits(&parse_paf(&tx_paf)) {
        tx_by_key.entry(h.q.clone()).or_default().push(h);
    }
    let mut chains_by_key: HashMap<String, Vec<Chain>> = HashMap::new();
    for c in gene_body_chains(&parse_paf(&body_paf)) {
        chains_by_key.entry(c.q.clone()).or_default().push(c);
    }
    let (pairs, n_exon, n_body) = edges(&nodes, &tx_by_key, &chains_by_key);
    let fams = triangle_leaders(&nodes, &pairs);
    eprintln!(
        "[shared-definition] edges: {n_exon} exon + {n_body} gene-body ({} pairs); {} triangle-supported families holding {} loci",
        pairs.len(), fams.len(), fams.iter().map(|f| f.len()).sum::<usize>()
    );
    Ok((nodes, fams, pairs))
}

/// A node as a catalog copy (sequence = its exons in transcript orientation).
pub fn node_transcript(n: &SdNode, genome: &GenomeIndex) -> DenovoTranscript {
    DenovoTranscript {
        tid: format!("SD~{}_{}_{}", n.chrom, n.start(), n.end()),
        chrom: n.chrom.clone(),
        start: n.start(),
        end: n.end(),
        n_reads: n.n_reads.min(u32::MAX as u64) as u32,
        strand: n.strand,
        introns: n.exons.windows(2).map(|w| (w[0].1, w[1].0)).collect(),
        seq: spliced_seq(genome, &n.chrom, &n.exons, n.strand),
        distinguishing_uniq: 0,
        core_bp: 0,
        stub: false,
        tes: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn node(chrom: &str, strand: char, reads: u64, exons: &[(u64, u64)]) -> SdNode {
        SdNode { chrom: chrom.into(), strand, n_reads: reads, exons: exons.to_vec(), rep_exons: exons.to_vec(), tx_chains: vec![exons.to_vec()] }
    }

    #[test]
    fn read_isoform_widening_admits_supported_chains_and_only_widens() {
        // One node with a 2-exon representative; reads carry a second isoform with an extra exon.
        let n = node("chr1", '+', 10, &[(100, 200), (500, 600)]);
        let rb = |blocks: &[(u64, u64)]| ReadBlocks { chrom: "chr1".into(), strand: '+', blocks: blocks.to_vec() };
        let mut reads: Vec<ReadBlocks> = Vec::new();
        for _ in 0..5 {
            reads.push(rb(&[(100, 200), (300, 350), (500, 600)]));
        }
        // A third chain seen twice only: below k, so it must not be admitted.
        for _ in 0..2 {
            reads.push(rb(&[(100, 200), (700, 800)]));
        }
        let (out, widened) = widen_with_read_isoforms(&[n.clone()], &reads, 5);
        assert_eq!(widened, 1);
        let w = &out[0];
        // Only widens: every original exon block survives, and the rep chain stays a query.
        for e in &n.exons {
            assert!(w.exons.iter().any(|x| x.0 <= e.0 && e.1 <= x.1), "{e:?} lost");
        }
        assert!(w.tx_chains.contains(&n.rep_exons));
        assert_eq!(w.rep_exons, n.rep_exons);
        // The supported isoform is now a query and its exon is in the union.
        assert_eq!(w.tx_chains.len(), 2);
        assert!(w.exons.iter().any(|&(a, b)| a <= 300 && 350 <= b));
        // The 2-read chain contributed nothing.
        assert!(!w.exons.iter().any(|&(a, b)| a <= 700 && 800 <= b));
    }

    #[test]
    fn read_isoform_widening_is_inert_without_support_and_on_unspliced_reads() {
        let n = node("chr1", '-', 4, &[(100, 200), (500, 600)]);
        let unspliced = ReadBlocks { chrom: "chr1".into(), strand: '-', blocks: vec![(100, 600)] };
        let wrong_strand = ReadBlocks { chrom: "chr1".into(), strand: '+', blocks: vec![(100, 200), (300, 350), (500, 600)] };
        let reads: Vec<ReadBlocks> = std::iter::repeat(unspliced).take(9)
            .chain(std::iter::repeat(wrong_strand).take(9)).collect();
        let (out, widened) = widen_with_read_isoforms(&[n.clone()], &reads, 5);
        assert_eq!(widened, 0);
        assert_eq!(out[0], n);
    }

    #[test]
    fn isoform_knobs_default_to_on_at_k_five() {
        if std::env::var("RUSTLE_SD_READ_ISOFORM").is_err() {
            assert!(isoform_enabled());
        }
        if std::env::var("RUSTLE_SD_ISOFORM_K").is_err() {
            assert_eq!(isoform_k(), 5);
            assert_eq!(ISOFORM_MIN_READS, 5);
        }
    }

    #[test]
    fn shared_definition_is_off_by_default() {
        if std::env::var("RUSTLE_SHARED_DEFINITION").is_err() {
            assert!(!enabled());
        }
    }

    #[test]
    fn merge_joins_touching_intervals() {
        assert_eq!(merge(vec![(5, 10), (0, 5), (20, 30), (25, 26)]), vec![(0, 10), (20, 30)]);
    }

    #[test]
    fn consolidate_cuts_long_introns_drops_small_pieces_and_merges_overlaps() {
        let a = node("c", '+', 5, &[(0, 200), (200 + MAX_INTRON + 1, 200 + MAX_INTRON + 50)]);
        let b = node("c", '+', 9, &[(150, 400)]);
        let c = node("c", '-', 1, &[(150, 400)]);
        let out = consolidate(&[a, b, c]);
        // the 49 bp tail piece is dropped; a's first piece and b merge on '+'; c stays on '-'
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].exons, vec![(0, 400)]);
        assert_eq!(out[0].n_reads, 14);
        assert_eq!(out[0].rep_exons, vec![(150, 400)]);
        assert_eq!(out[1].strand, '-');
    }

    #[test]
    fn depth_and_split_match_the_prototype_toy() {
        let a: Vec<(u64, u64)> = vec![(0, 100), (200, 300)];
        let b: Vec<(u64, u64)> = vec![(1000, 1100)];
        let rt: Vec<(u64, u64)> = vec![(250, 300), (1000, 1050)];
        let lists: Vec<&[(u64, u64)]> = vec![&a, &a, &a, &b, &b, &b, &rt];
        let ex = depth2_exons(&lists);
        assert_eq!(ex, vec![(0, 100), (200, 300), (1000, 1100)]);
        assert_eq!(split_linked(&ex, &lists), vec![(vec![0, 1], 4), (vec![2], 4)]);
    }

    #[test]
    fn read_locus_nodes_skip_groups_touching_existing_nodes_unless_split() {
        let base = vec![node("c", '+', 10, &[(1000, 1100)])];
        let mk = |b: Vec<(u64, u64)>| ReadBlocks { chrom: "c".into(), strand: '+', blocks: b };
        let mut reads = vec![mk(vec![(0, 100), (200, 300)]); 3];
        reads.extend(vec![mk(vec![(1000, 1100)]); 3]);
        reads.push(mk(vec![(250, 300), (1000, 1050)]));
        let (whole, n_whole) = with_read_locus_nodes(&base, &reads, false);
        assert_eq!(n_whole, 0);
        assert_eq!(whole.len(), 1);
        let (split, n_split) = with_read_locus_nodes(&base, &reads, true);
        assert_eq!(n_split, 1);
        assert_eq!(split[0].exons, vec![(0, 100), (200, 300)]);
        assert_eq!(split[0].n_reads, 4);
    }

    #[test]
    fn transcript_blocks_and_hits() {
        let h = PafRec {
            q: "q".into(), qlen: 830, qs: 0, qe: 830, strand: '+', chrom: "c".into(), clen: 5000, ts: 2000, te: 4230,
            nm: 830, bl: 830, cg: "300M600N250M800N280M".into(),
        };
        assert_eq!(tx_exon_blocks(&h), vec![(2000, 2300), (2900, 3150), (3950, 4230)]);
        assert_eq!(transcript_hits(&[h.clone()]).len(), 1);
        let low = PafRec { nm: 600, ..h };
        assert!(transcript_hits(&[low]).is_empty());
    }

    #[test]
    fn gene_body_chain_joins_ordered_records_and_applies_floors() {
        let r = |qs, qe, ts, te| PafRec {
            q: "b".into(), qlen: 1000, qs, qe, strand: '+', chrom: "c".into(), clen: 100_000, ts, te,
            nm: 95, bl: 100, cg: String::new(),
        };
        let chains = gene_body_chains(&[r(500, 1000, 10_600, 11_100), r(0, 400, 10_000, 10_400)]);
        assert_eq!(chains, vec![Chain { q: "b".into(), chrom: "c".into(), strand: '+', s: 10_000, e: 11_100 }]);
        // a lone 200 bp record covers < 0.50 of the 1,000 bp body
        assert!(gene_body_chains(&[r(0, 200, 50_000, 50_200)]).is_empty());
    }

    #[test]
    fn triangle_leaders_stop_at_one_hop_unless_two_star_members_support() {
        // leader 0 (most reads) - {1,2}; 3 touches 1 and 2 (joins); 4 touches only 2 (does not); 4-5 pair
        let nodes: Vec<SdNode> = (0..6).map(|i| node("c", '+', if i == 0 { 50 } else { 5 }, &[(i * 1000, i * 1000 + 500)])).collect();
        let pairs: BTreeSet<(usize, usize)> = [(0, 1), (0, 2), (1, 3), (2, 3), (2, 4), (4, 5)].into_iter().collect();
        assert_eq!(triangle_leaders(&nodes, &pairs), vec![vec![0, 1, 2, 3], vec![4, 5]]);
    }

    #[test]
    fn edges_ignore_own_locus_and_wrong_orientation() {
        let nodes = vec![
            node("c", '+', 5, &[(0, 100), (200, 300)]),
            node("c", '+', 5, &[(10_000, 10_100), (10_200, 10_300)]),
            node("c", '-', 5, &[(20_000, 20_100), (20_200, 20_300)]),
        ];
        let hit = |ts: u64, te: u64, strand: char| PafRec {
            q: tx_key(&nodes[0]), qlen: 200, qs: 0, qe: 200, strand, chrom: "c".into(), clen: 100_000, ts, te,
            nm: 200, bl: 200, cg: "100M100N100M".into(),
        };
        let mut tx: HashMap<String, Vec<PafRec>> = HashMap::new();
        tx.insert(tx_key(&nodes[0]), vec![hit(0, 300, '+'), hit(10_000, 10_300, '+'), hit(20_000, 20_300, '+')]);
        let (pairs, n_exon, n_body) = edges(&nodes, &tx, &HashMap::new());
        assert_eq!(pairs.into_iter().collect::<Vec<_>>(), vec![(0, 1)]);
        assert_eq!((n_exon, n_body), (1, 0));
    }
}

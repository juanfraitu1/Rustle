//! O3 candidate copies (spec docs/superpowers/specs/2026-10-02-o3-candidates-design.md): the reference-absent-copy chain of
//! PREREG_rna_allele_haplotype_count Amendments 7-11 without IsoCon. Pure functions here, plus the cached minimap2 runner and the output
//! writers; the BAM passes and the flow of the stage are in the binary.
//!
//! **STATUS:** OTHER-BINARY — reached only from the `o3_candidates` binary (`src/bin/o3_candidates.rs`; the driver's `candidates` stage runs that binary, an OPT-IN stage since ruling R14: Amendment 12 failed, its re-run A13 passed (docs/O3_CANDIDATES_ACCEPTANCE_A13_2026-10-03.md); and Amendment 14 failed, docs/O3_CANDIDATES_CONTROL_A14_2026-10-03.md)  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use crate::vg_family::run_cache as rc;
use anyhow::Context;
use std::collections::{BTreeMap, HashMap};
use std::io::Write;
use std::path::Path;

pub const KMER_K: usize = 31;
pub const SKETCH_W: usize = 5;
/// Reads of the attribution set shorter than this are never aligned nor attributed (prereg Amendment 13: unmapped reads >= 300 bp only;
/// Amendment 13c, ruling R19: the floor holds for the poorly placed reads as well).
pub const MIN_UNMAPPED_LEN: usize = 300;
/// `attribute_by_hits` (prereg Amendment 13b): the best hit must cover at least this fraction of the READ (`(qe - qs) / qlen`, whatever the
/// target's length) ...
pub const ATTRIB_MIN_READ_COV: f64 = 0.5;
/// ... and its gap-compressed divergence `de` must be at most this: the family definition's identity floor (0.80) in read space.
pub const ATTRIB_MAX_DE: f64 = 0.20;

/// `is_poorly_placed` (prereg Amendment 13b): a primary record whose `de` exceeds this places its read poorly.
pub const POORLY_PLACED_DE: f32 = 0.02;

/// `structural_scores` (prereg Amendment 13e): a member is eligible to be its cluster's template when it is aligned to at least
/// min(0.5 x (n - 1), `ELIGIBLE_PARTNER_CAP`) of the n - 1 other members. The cap is half of the all-vs-all's 100 hits per query (`MM2_AVA`'s
/// `-N 100`, with `--dual=no`), which keeps most members of a cluster of several hundred reads under half of the others.
pub const ELIGIBLE_PARTNER_CAP: usize = 50;

/// Prereg Amendment 13b: a mapped read joins the attribution set ("poorly placed") when no record of it lies on a copy of the run's families
/// (`netted` false: it is in no net of pass A; ruling R18) and its PRIMARY record has `de > POORLY_PLACED_DE` or MAPQ 0. The caller passes the
/// primary's `de` (absent -> 0.0: MAPQ alone decides) and its MAPQ (unavailable, 255, is not 0).
pub fn is_poorly_placed(de: f32, mapq: u8, netted: bool) -> bool { !netted && (de > POORLY_PLACED_DE || mapq == 0) }

fn code(b: u8) -> Option<u64> { match b { b'A' | b'a' => Some(0), b'C' | b'c' => Some(1), b'G' | b'g' => Some(2), b'T' | b't' => Some(3), _ => None } }

/// Canonical (min of forward / reverse-complement) 2-bit k-mers, 1 <= k <= 32 (asserted); windows with N are skipped.
pub fn canonical_kmers(seq: &[u8], k: usize) -> Vec<u64> {
    assert!((1..=32).contains(&k), "canonical_kmers: k = {k} must lie in 1..=32 (2-bit codes in a u64)");
    let mask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    let mut out = Vec::with_capacity(seq.len().saturating_sub(k) + 1);
    for &b in seq {
        match code(b) {
            Some(c) => { fw = ((fw << 2) | c) & mask; rv = (rv >> 2) | ((3 - c) << (2 * (k - 1))); valid += 1; }
            None => { valid = 0; fw = 0; rv = 0; }
        }
        if valid >= k { out.push(fw.min(rv)); }
    }
    out
}

fn mix(x: u64) -> u64 { let mut z = x.wrapping_add(0x9E3779B97F4A7C15); z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }

/// (k, w) minimizers over the canonical k-mers: the smallest hashed k-mer of every window of w consecutive k-mers, deduplicated
/// consecutively; sorted and deduplicated so sketches compare by merge. w >= 1 (asserted; k as in `canonical_kmers`).
pub fn minimizer_sketch(seq: &[u8], k: usize, w: usize) -> Vec<u64> {
    assert!(w >= 1, "minimizer_sketch: the window w = {w} must be at least 1");
    let km: Vec<u64> = canonical_kmers(seq, k).into_iter().map(mix).collect();
    let mut out = Vec::new();
    if km.len() < w { out.extend(km.iter().copied().min()); }
    for win in km.windows(w) { let m = *win.iter().min().unwrap(); if out.last() != Some(&m) { out.push(m); } }
    out.sort_unstable(); out.dedup(); out
}

/// Shared minimizers / size of the smaller sketch (both sorted, deduplicated).
pub fn sketch_share(a: &[u64], b: &[u64]) -> f64 {
    let (mut i, mut j, mut shared) = (0, 0, 0usize);
    while i < a.len() && j < b.len() { if a[i] == b[j] { shared += 1; i += 1; j += 1; } else if a[i] < b[j] { i += 1; } else { j += 1; } }
    shared as f64 / a.len().min(b.len()).max(1) as f64
}

/// One minimap2 PAF alignment: columns 1-11 (mapq is not kept) plus the `de:f:` and `cs:Z:` tags. Coordinates are 0-based half-open as in
/// PAF; `strand` is `b'+'` or `b'-'`; `de` is minimap2's gap-compressed divergence (1.0 when the tag is absent); `cs` is the raw string.
#[derive(Clone, Debug, PartialEq)]
pub struct PafHit {
    pub q: String, pub qlen: usize, pub qs: usize, pub qe: usize, pub strand: u8,
    pub t: String, pub tlen: usize, pub ts: usize, pub te: usize,
    pub matches: usize, pub block: usize, pub de: f64, pub cs: Option<String>,
}

/// Parses minimap2 PAF text into one `PafHit` per alignment line. Lenient, like the repo's other PAF readers
/// (`shared_definition::parse_paf`): a line with fewer than 12 tab-separated columns, a non-numeric fixed column or a strand other than
/// `+`/`-` (blank lines, a truncated last line) is skipped. Tags are looked up after column 12: `de:f:` (absent or unreadable -> 1.0, the
/// worst divergence: such a hit never clusters and never merges; linking is decided by d, not `de`) and `cs:Z:` (absent -> `None`).
pub fn parse_paf(text: &str) -> Vec<PafHit> { text.lines().filter_map(parse_paf_line).collect() }

fn parse_paf_line(line: &str) -> Option<PafHit> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 12 { return None; }
    let tag = |prefix: &str| f[12..].iter().find_map(|x| x.strip_prefix(prefix));
    Some(PafHit {
        q: f[0].to_string(), qlen: f[1].parse().ok()?, qs: f[2].parse().ok()?, qe: f[3].parse().ok()?,
        strand: match f[4] { "+" => b'+', "-" => b'-', _ => return None },
        t: f[5].to_string(), tlen: f[6].parse().ok()?, ts: f[7].parse().ok()?, te: f[8].parse().ok()?,
        matches: f[9].parse().ok()?, block: f[10].parse().ok()?,
        de: tag("de:f:").and_then(|v| v.parse::<f64>().ok()).unwrap_or(1.0),
        cs: tag("cs:Z:").map(str::to_string),
    })
}

/// Per query name, the hit with the most matching bases; a tie keeps the first encountered (stable).
pub fn best_by_matches(hits: &[PafHit]) -> HashMap<String, PafHit> {
    let mut best: HashMap<String, PafHit> = HashMap::new();
    for h in hits { if best.get(&h.q).map_or(true, |b| h.matches > b.matches) { best.insert(h.q.clone(), h.clone()); } }
    best
}

/// Per query name, the hit with the highest identity x coverage (`id_cov`); a tie keeps the first encountered (stable). This chooses the GENOME
/// hit that `classify` judges (ruling R4), as `best_hits` of `bench/rna_allele/link_test.py` does; pairwise (all-vs-all) hits keep
/// `best_by_matches`, as `best_pairs` of `bench/rna_allele/merge_test.py` does.
pub fn best_by_id_cov(hits: &[PafHit]) -> HashMap<String, PafHit> {
    let mut best: HashMap<String, PafHit> = HashMap::new();
    for h in hits { if best.get(&h.q).map_or(true, |b| id_cov(h) > id_cov(b)) { best.insert(h.q.clone(), h.clone()); } }
    best
}

/// The family each read of the attribution set joins (prereg Amendments 13 / 13b, which retire spec §5.2's k-mer index), from its
/// `MM2_ATTRIB` hits on the targets (this run's net reads, `{family}|{read}`, and every family's copy sequences): per read, its best hit by
/// matches (the first encountered on a tie, as `best_by_matches`) among the hits whose target `family_of_target` names, and the read joins
/// that target's family iff the hit covers >= `ATTRIB_MIN_READ_COV` of the READ (the query's aligned span over its length,
/// `(qe - qs) / qlen`; NOT `shorter_cov`, which would judge the target's covered fraction when the read is the longer sequence) and its `de`
/// is <= `ATTRIB_MAX_DE`. A best hit that fails either joins the read to no family: a lesser hit never stands in (a read whose best hit is a
/// 40% shared exon of another family's copy stays out). A target absent from the map (a partner row, a record the copies table does not
/// hold) is ignored. Returns read -> family.
pub fn attribute_by_hits(hits: &[PafHit], family_of_target: &HashMap<String, String>) -> HashMap<String, String> {
    let mut best: HashMap<&str, (&PafHit, &String)> = HashMap::new();
    for h in hits {
        let Some(family) = family_of_target.get(&h.t) else { continue };
        if best.get(h.q.as_str()).map_or(true, |(b, _)| h.matches > b.matches) {
            best.insert(h.q.as_str(), (h, family));
        }
    }
    // the read's covered fraction, guarded as the quantities below are (a zero length counts as 1, an end before its start as no span)
    let read_cov = |h: &PafHit| h.qe.saturating_sub(h.qs) as f64 / h.qlen.max(1) as f64;
    best.into_iter()
        .filter(|(_, (h, _))| read_cov(h) >= ATTRIB_MIN_READ_COV && h.de <= ATTRIB_MAX_DE)
        .map(|(read, (_, family))| (read.to_string(), family.clone()))
        .collect()
}

// The three quantities below divide by lengths from the PAF; a zero length counts as 1 (no NaN / inf) and an end before its start as an
// empty span (no usize wrap-around).

/// Identity x coverage = `matches / block x (qe - qs) / qlen`, in the chain's order of evaluation (`bench/rna_allele/link_test.py:141`).
/// `>= 0.999` means the sequence is already in the reference; below it the sequence is flagged.
pub fn id_cov(h: &PafHit) -> f64 { h.matches as f64 / h.block.max(1) as f64 * h.qe.saturating_sub(h.qs) as f64 / h.qlen.max(1) as f64 }

/// Whole-length divergence `1 - matches / qlen` (`bench/rna_allele/yag_test.py:245`): the quantity the link rule compares with delta.
pub fn whole_length_d(h: &PafHit) -> f64 { 1.0 - h.matches as f64 / h.qlen.max(1) as f64 }

/// The aligned span on the SHORTER of query / target over that sequence's length; the query on a tie (`qlen <= tlen`, as `best_pairs` in
/// `bench/rna_allele/merge_test.py`). The merge and clustering rules join a pair when this is >= 0.5.
pub fn shorter_cov(h: &PafHit) -> f64 {
    if h.qlen <= h.tlen { h.qe.saturating_sub(h.qs) as f64 / h.qlen.max(1) as f64 } else { h.te.saturating_sub(h.ts) as f64 / h.tlen.max(1) as f64 }
}

/// One operation of a minimap2 `cs:Z:` string (the short form written by `--cs`; the long form of `--cs=long` is not supported).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CsOp {
    /// `:N` — N identical bases.
    Eq(usize),
    /// `*xy` — a substitution, reference base x -> query base y, bytes as written (lowercase from minimap2).
    Sub(u8, u8),
    /// `+seq` — bases in the query only (an insertion in the query).
    Ins(Vec<u8>),
    /// `-seq` — bases in the reference only (a deletion from the reference).
    Del(Vec<u8>),
    /// `~xxNNyy` — an intron of NN reference bases (the donor `xx` and acceptor `yy` dinucleotides are dropped).
    Intron(usize),
}

fn is_cs_base(b: u8) -> bool { matches!(b.to_ascii_lowercase(), b'a' | b'c' | b'g' | b't' | b'n') }

/// The decimal number at `b[i..]` and the index after it; `None` when there is no digit or the value overflows `usize`.
fn cs_number(b: &[u8], mut i: usize) -> Option<(usize, usize)> {
    let start = i;
    let mut n = 0usize;
    while i < b.len() && b[i].is_ascii_digit() { n = n.checked_mul(10)?.checked_add((b[i] - b'0') as usize)?; i += 1; }
    (i > start).then_some((n, i))
}

fn cs_err(b: &[u8], at: usize, what: &str) -> anyhow::Error {
    anyhow::anyhow!("malformed cs string: {what} at byte {at}, near `{}`", String::from_utf8_lossy(&b[at.saturating_sub(12)..(at + 12).min(b.len())]))
}

/// Parses a minimap2-produced `cs:Z:` short string into its operations (`:N` `*xy` `+seq` `-seq` `~xxNNyy`; bases are `acgtn`, either case).
/// Anything else is a hard error naming the byte offset and its neighbourhood, never a partial list and never a panic: an unknown operation
/// character, `:` or the intron length without digits, a number that overflows `usize`, `*` or `~` without its two bases, `+`/`-` without bases.
pub fn parse_cs(cs: &str) -> anyhow::Result<Vec<CsOp>> {
    let b = cs.as_bytes();
    let (mut ops, mut i) = (Vec::new(), 0usize);
    let two_bases = |i: usize| i + 1 < b.len() && is_cs_base(b[i]) && is_cs_base(b[i + 1]);
    while i < b.len() {
        let at = i;
        i += 1;
        match b[at] {
            b':' => {
                let (n, j) = cs_number(b, i).ok_or_else(|| cs_err(b, at, "`:` without a length"))?;
                ops.push(CsOp::Eq(n)); i = j;
            }
            b'*' => {
                if !two_bases(i) { return Err(cs_err(b, at, "`*` without two bases")); }
                ops.push(CsOp::Sub(b[i], b[i + 1])); i += 2;
            }
            op @ (b'+' | b'-') => {
                let j = i + b[i..].iter().take_while(|&&x| is_cs_base(x)).count();
                if j == i { return Err(cs_err(b, at, "`+`/`-` without bases")); }
                let seq = b[i..j].to_vec();
                ops.push(if op == b'+' { CsOp::Ins(seq) } else { CsOp::Del(seq) }); i = j;
            }
            b'~' => {
                if !two_bases(i) { return Err(cs_err(b, at, "`~` without a two-base donor")); }
                let (n, j) = cs_number(b, i + 2).ok_or_else(|| cs_err(b, at, "`~` without an intron length"))?;
                if !two_bases(j) { return Err(cs_err(b, at, "`~` without a two-base acceptor")); }
                ops.push(CsOp::Intron(n)); i = j + 2;
            }
            _ => return Err(cs_err(b, at, "unknown operation")),
        }
    }
    Ok(ops)
}

/// Indels of at least this many bases are STRUCTURE (isoforms, ruling R2), not sequencing errors; shorter ones follow the 50% majority.
const STRUCT_MIN_INDEL: usize = 20;
/// Members that must cover a column before the vote (or a < 20 bp indel majority) may change the template there.
const VOTE_MIN_COVER: usize = 3;
/// Members that must carry an insertion of >= 20 bp before it may enter the consensus (prereg Amendment 15: the floor of 3,
/// AND a majority of the members covering that column — `2 x count >= covering`; the pre-15 "whatever its share" rule duplicated
/// one biological insertion at every column a minority carrier happened to land it, see docs/O3_CANDIDATES_CONSENSUS_DEFECT_2026-10-03.md).
const STRUCT_MIN_SUPPORT: usize = 3;
/// Members that must cover a template column at either end of the consensus for it to survive the end trim (spec §5.4, ruling R5).
const END_MIN_COVER: usize = 2;

/// Clusters one family's reads (spec §5.3, with minimap2's all-vs-all hits in place of the greedy pass; ruling §9b): union-find over the reads,
/// joining a pair when its best hit covers >= 50% of the shorter read (`shorter_cov`) with `de <= delta`. The pair's best hit is the one with
/// the most matching bases, the first on a tie, whichever read is `q` (`best_pairs` of `bench/rna_allele/merge_test.py`, which this
/// reproduces). Self hits and hits naming a read that is not in `names` are ignored; the strand is not consulted (`de` and the spans do not
/// depend on orientation). `names` must be distinct. Returns the clusters as ascending index lists ordered by (size descending, first
/// index). The clusters are a function of the set of joined pairs, so neither the hash order of the pair table nor the order of `ava`
/// (beyond the first-wins tie rule) reaches the output.
pub fn cluster_reads(names: &[String], ava: &[PafHit], delta: f64) -> Vec<Vec<usize>> {
    let idx: HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut best: HashMap<(usize, usize), &PafHit> = HashMap::new();
    for h in ava {
        let (Some(&a), Some(&b)) = (idx.get(h.q.as_str()), idx.get(h.t.as_str())) else { continue };
        if a == b { continue; }
        let key = (a.min(b), a.max(b));
        if best.get(&key).map_or(true, |o| h.matches > o.matches) { best.insert(key, h); }
    }
    fn find(p: &mut [usize], mut x: usize) -> usize { while p[x] != x { p[x] = p[p[x]]; x = p[x]; } x }
    let mut par: Vec<usize> = (0..names.len()).collect();
    for (&(a, b), h) in &best {
        if shorter_cov(h) >= 0.5 && h.de <= delta {
            let (ra, rb) = (find(&mut par, a), find(&mut par, b));
            if ra != rb { par[ra.max(rb)] = ra.min(rb); }
        }
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..names.len() { let r = find(&mut par, i); groups.entry(r).or_default().push(i); }   // pushed in index order: each group is ascending
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    out.sort_by(|a, b| b.len().cmp(&a.len()).then(a[0].cmp(&b[0])));                                // first indices are distinct: a total order
    out
}

/// The bases of the indels of at least `STRUCT_MIN_INDEL` (20) bp in one pairwise `cs`: insertions (`+`), deletions (`-`) and introns (`~`,
/// a stretch of the target the query skips) alike, so the count is the same whichever of the two sequences is the query.
fn big_indel_bases(ops: &[CsOp]) -> u64 {
    ops.iter()
        .map(|op| match op {
            CsOp::Ins(s) | CsOp::Del(s) if s.len() >= STRUCT_MIN_INDEL => s.len() as u64,
            CsOp::Intron(n) if *n >= STRUCT_MIN_INDEL => *n as u64,
            _ => 0,
        })
        .sum()
}

/// The terminal bases of one sequence that a pairwise alignment leaves uncovered (prereg Amendment 13d): the `start` bases before its aligned
/// span and the `len - end` after it, each end counted only from `STRUCT_MIN_INDEL` (20) bp. Both ends count alike, so the strand of the hit
/// (which end is the read's 5' end) does not matter.
fn uncovered_ends(start: usize, end: usize, len: usize) -> u64 {
    [start, len.saturating_sub(end)].into_iter().filter(|&bases| bases >= STRUCT_MIN_INDEL).map(|bases| bases as u64).sum()
}

/// One member's place in its cluster under the structural distance of prereg Amendment 13d (`structural_scores`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct StructuralScore {
    /// The member: an index into the net's names.
    pub member: usize,
    /// Its aligned partners: the other members it has a hit with in the all-vs-all, in either direction.
    pub aligned: usize,
    /// The sum of d(member, p) over its aligned partners.
    pub sum_d: u64,
    /// Aligned to >= min(0.5 x (n - 1), `ELIGIBLE_PARTNER_CAP`) of the cluster's n - 1 other members (prereg Amendment 13e).
    pub eligible: bool,
}
impl StructuralScore {
    /// The mean d over its aligned partners; `None` without one.
    pub fn mean_d(&self) -> Option<f64> { (self.aligned > 0).then(|| self.sum_d as f64 / self.aligned as f64) }
}

/// How a cluster's template was chosen (prereg Amendments 13 / 13d / 13e).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TemplateChoice {
    /// The medoid: the eligible member with the lowest mean d over its aligned partners (`structural_template`); also the only member of a
    /// cluster of one.
    Medoid(usize),
    /// No eligible member (none is aligned to min(half of the others, `ELIGIBLE_PARTNER_CAP`) of them), and some member has an aligned
    /// partner: the longest member WITH an aligned partner (a mean), the smaller name on equal lengths. A partner-less member is passed over
    /// however long it is (Amendment 13e: it is never chosen while another member has a mean). The binary logs it.
    LongestAligned(usize),
    /// No member has an aligned partner (no mean anywhere): the longest member, the smaller name on equal lengths (Amendment 13e). The binary
    /// logs it apart from `LongestAligned`.
    LongestUnaligned(usize),
    /// `refined_template` only: the refinement kept the cluster's template, which stays.
    Kept(usize),
}
impl TemplateChoice {
    /// The chosen member.
    pub fn member(self) -> usize {
        match self {
            TemplateChoice::Medoid(m) | TemplateChoice::LongestAligned(m) | TemplateChoice::LongestUnaligned(m) | TemplateChoice::Kept(m) => m,
        }
    }
}

/// Per member of one cluster, its aligned partners, its summed structural distance to them and its eligibility (prereg Amendment 13d).
/// A pair's alignment is its best hit in the net's all-vs-all `ava` (`MM2_AVA`): the most matches, the first on a tie, whichever of the two
/// reads is the query (the pair rule of `cluster_reads`); a pair with no hit is not aligned and adds nothing. For a member m and an aligned
/// partner p, d(m, p) = the bases of the indels of >= 20 bp in that alignment (`big_indel_bases`: insertions, deletions and `~` introns
/// alike, the same from either side) PLUS p's terminal bases the alignment leaves uncovered, each end counted from 20 bp (`uncovered_ends`:
/// `ts` and `tlen - te` when p is the target, `qs` and `qlen - qe` when it is the query). d is not symmetric: a fragment inside its partner
/// pays the partner's uncovered ends, the partner pays nothing for the fragment, and a member's own unaligned end costs its partners, not
/// it. A best hit without a `cs` aligns its pair with the uncovered ends as d. A member is eligible when it is aligned to >= min(0.5 x
/// (n - 1), `ELIGIBLE_PARTNER_CAP`) of the cluster's n - 1 other members (Amendment 13e: one partner in a cluster of 2 or 3, half of the
/// others up to 101 members, 50 beyond, where the all-vs-all's 100 hits per query keep most members under half). Hits naming a read outside
/// `members` and self hits are ignored. `members` index `names` (the net's reads, distinct) and are distinct. Returned in `members` order; a malformed `cs` of a pair's
/// best hit is an error naming the pair (never a silent 0). No hash order reaches the result (pairs are visited in index order).
pub fn structural_scores(members: &[usize], names: &[String], ava: &[PafHit]) -> anyhow::Result<Vec<StructuralScore>> {
    let idx: HashMap<&str, usize> = members.iter().enumerate().map(|(i, &m)| (names[m].as_str(), i)).collect();
    let mut best: BTreeMap<(usize, usize), &PafHit> = BTreeMap::new();
    for h in ava {
        let (Some(&a), Some(&b)) = (idx.get(h.q.as_str()), idx.get(h.t.as_str())) else { continue };
        if a == b { continue; }
        let key = (a.min(b), a.max(b));
        if best.get(&key).map_or(true, |o| h.matches > o.matches) { best.insert(key, h); }
    }
    let n = members.len();
    let (mut aligned, mut sum_d) = (vec![0usize; n], vec![0u64; n]);
    for h in best.into_values() {
        let indel = match h.cs.as_deref() {
            Some(cs) => big_indel_bases(&parse_cs(cs).with_context(|| format!("the structural template: {} against {}", h.q, h.t))?),
            None => 0,
        };
        let (q, t) = (idx[h.q.as_str()], idx[h.t.as_str()]);
        // the query pays the target's uncovered ends, the target the query's
        sum_d[q] += indel + uncovered_ends(h.ts, h.te, h.tlen);
        sum_d[t] += indel + uncovered_ends(h.qs, h.qe, h.qlen);
        aligned[q] += 1;
        aligned[t] += 1;
    }
    // aligned >= min(0.5 x (n - 1), cap), in integers: 2 x aligned >= min(n - 1, 2 x cap)
    let need_twice = n.saturating_sub(1).min(2 * ELIGIBLE_PARTNER_CAP);
    Ok((0..n).map(|i| StructuralScore { member: members[i], aligned: aligned[i], sum_d: sum_d[i], eligible: 2 * aligned[i] >= need_twice }).collect())
}

/// The template of a cluster (prereg Amendments 13d / 13e, the medoid under the structural distance, in place of Amendment 13's lowest total of
/// indel bases and of the longest read before it): among the eligible members with an aligned partner (`structural_scores`), the one with
/// the LOWEST mean d over its aligned partners (compared exactly, as cross-multiplied sums), ties to the longest (`lens`, indexed like
/// `names`), then to the smallest name: `Medoid`. A read that retains an intron or skips an exon pays the indel in every pair, a fragment pays
/// its partners' uncovered ends, and a member aligned to fewer than min(half of the others, 50) is not eligible (Amendment 13e). A member
/// with no aligned partner has no mean and is never chosen while another member has one, the fallbacks included: without an eligible member
/// the template is the longest member that has an aligned partner (`LongestAligned`), and only when no member has one the longest member
/// (`LongestUnaligned`); the binary logs both. A cluster of one member is its own template; an empty one is an error, and so is a malformed
/// `cs`. A function of the cluster as a set, the hits and the names: neither the order of `members` nor (beyond
/// the first-wins tie between two hits of one pair) the order of `ava` reaches it.
pub fn structural_template(members: &[usize], names: &[String], ava: &[PafHit], lens: &[usize]) -> anyhow::Result<TemplateChoice> {
    anyhow::ensure!(!members.is_empty(), "the structural template of a cluster without members");
    if let [only] = members { return Ok(TemplateChoice::Medoid(*only)); }
    let scores = structural_scores(members, names, ava)?;
    let medoid = scores.iter().filter(|s| s.eligible && s.aligned > 0).min_by(|a, b| {
        (u128::from(a.sum_d) * b.aligned as u128)
            .cmp(&(u128::from(b.sum_d) * a.aligned as u128))
            .then(lens[b.member].cmp(&lens[a.member]))
            .then_with(|| names[a.member].cmp(&names[b.member]))
    });
    Ok(match medoid {
        Some(s) => TemplateChoice::Medoid(s.member),
        // Amendment 13e: a member with no aligned partner has no mean and is never chosen while another member has one
        None => match longest_of(scores.iter().filter(|s| s.aligned > 0).map(|s| s.member), names, lens) {
            Some(m) => TemplateChoice::LongestAligned(m),
            None => TemplateChoice::LongestUnaligned(longest_of(members.iter().copied(), names, lens).expect("a non-empty cluster has a longest member")),
        },
    })
}

/// The longest of `pool` (`lens`, indexed like `names`), the smaller name on equal lengths; `None` for an empty pool.
fn longest_of(pool: impl Iterator<Item = usize>, names: &[String], lens: &[usize]) -> Option<usize> {
    pool.max_by(|&a, &b| lens[a].cmp(&lens[b]).then_with(|| names[b].cmp(&names[a])))
}

/// The template of a refined cluster's kept set (prereg Amendment 13): `Kept(template)` while the refinement keeps the template, else (the
/// refinement split it off) the `structural_template` of `kept`, over the pairs of the kept set only.
pub fn refined_template(template: usize, kept: &[usize], names: &[String], ava: &[PafHit], lens: &[usize]) -> anyhow::Result<TemplateChoice> {
    if kept.contains(&template) { Ok(TemplateChoice::Kept(template)) } else { structural_template(kept, names, ava, lens) }
}

/// Template-and-vote consensus of one cluster (spec §5.4, rulings R2 and R5). `template` is the sequence the members were aligned to; each
/// member's hit carries its `cs` against it (`ts..te` on the template). The member sequences are not read (the `cs` holds every base the vote
/// needs) and the strand is not consulted (a `cs` always lies along the target's forward strand). A member without a `cs` abstains; a
/// malformed `cs` is an error naming the member. The template votes (and covers) only if it is itself among the members, as a self hit.
/// Per template column, over the members that COVER it: `Eq`, `Sub` and `< 20 bp` `Del` columns cover; an `Intron` that a member splices out
/// and a `>= 20 bp` `Del` (an isoform's absence) do not: such a member carries no base there, so it is neither a vote for the template base
/// nor a vote to delete (R5).
/// * a column covered by >= 3 members takes the substituted base (the query base, upper case) with the most votes when that count exceeds the
///   number of covering members that carry the template base (no substitution, no deletion there); a tie keeps the template; two variants
///   with equal support take the smaller base; fewer than 3 covering members keep the template base;
/// * the insertions before a column (after the last column when `t == n`) are voted by size class, the long class first (prereg Amendment
///   13), and at most ONE of them is inserted. First the insertions >= 20 bp, which are STRUCTURE, not errors (R2): the most frequent one
///   (ties to the smaller sequence) is inserted when >= 3 members carry it AND its carriers are a majority of the members covering that
///   column (`2 x count >= covering`, prereg Amendment 15 — the pre-15 "whatever its share" rule is what duplicated one biological
///   insertion at six columns in GWFAM37:c1, reproduced byte for byte in docs/O3_CANDIDATES_CONSENSUS_DEFECT_2026-10-03.md). Only when no
///   long insertion passes, the insertions < 20 bp: the most frequent one (ties to the smaller) is inserted when >= 50% of >= 3 covering
///   members carry it. One winner, because two insertions before one column would make a sequence no member carries (a member's `cs` holds
///   at most one insertion there); the long class first, because 3 reads carrying an exon are evidence of an isoform while a short
///   insertion, even a majority one, is error-sized: 3 carriers of a 24-bp insertion win over 4 of 8 carrying a 1-bp one at that column
///   (when the 3 are a majority of the covering members);
/// * before the vote, each member's `cs` is normalised (prereg Amendment 15b): an insertion that begins with the template bases its
///   adjacent skip removes (`+X·E ~|X|`, minimap2 splice:hq's placement of an exon the template lacks, 107 of 134 insertion-plus-skip pairs
///   in the defective clusters) is shortened by them and the skip becomes matches, so the one event lands at one column; and each distinct
///   long insertion enters the consensus at most ONCE, at the column with the most carriers (Amendment 15b's one-copy rule: carriers are
///   not re-counted across columns)
/// * a deletion < 20 bp removes the column when >= 50% of >= 3 covering members delete it (a plain majority of the covering members, also
///   inside an exon that other members skip); a deletion >= 20 bp is STRUCTURE (R2) and never applied (the template's exon stays), so the
///   consensus is the exon union of the cluster's isoforms with SNV-level majority voting.
///
/// End trim (spec §5.4, R5): after the vote, the leading and trailing template columns covered by fewer than 2 members are dropped. An
/// insertion goes with the column it precedes (the one after the last column, with the last column). When no column is covered by 2 members
/// nothing survives: the consensus is empty (as it is for an empty member list). The output is upper case. A `cs` that runs past the template
/// end is clamped (nothing beyond the end is voted or appended); the frame of the `cs` (`ts` and the template) is trusted, not checked.
pub fn consensus_from_template(template: &[u8], member_hits: &[(&[u8], &PafHit)]) -> anyhow::Result<Vec<u8>> {
    let n = template.len();
    let mut cover = vec![0usize; n];
    let mut subs: Vec<HashMap<u8, usize>> = vec![HashMap::new(); n];
    let mut dels = vec![0usize; n];
    let mut ins: Vec<HashMap<Vec<u8>, usize>> = vec![HashMap::new(); n + 1];
    for (_, h) in member_hits {
        let Some(cs) = h.cs.as_deref() else { continue };
        let ops = parse_cs(cs).map_err(|e| anyhow::anyhow!("member {} against {}: {e}", h.q, h.t))?;
        // Amendment 15b (a): `+X·E ~|X|` -> `:|X| +E` before the vote (the insertion begins with the skipped template bases)
        let ops = normalize_cs_ops(ops, template, h.ts);
        let mut t = h.ts;
        for op in ops {
            match op {
                CsOp::Eq(len) => { for p in t..t.saturating_add(len).min(n) { cover[p] += 1; } t = t.saturating_add(len); }
                CsOp::Sub(_, qb) => { if t < n { cover[t] += 1; *subs[t].entry(qb.to_ascii_uppercase()).or_insert(0) += 1; } t = t.saturating_add(1); }
                // a < 20 bp deletion is an error-sized gap: the member covers the columns and votes to delete them; a >= 20 bp one is an isoform's
                // absence (R5, as for an `Intron`): the member carries no base there, so it neither covers nor votes
                CsOp::Del(seq) => {
                    if seq.len() < STRUCT_MIN_INDEL { for p in t..t.saturating_add(seq.len()).min(n) { cover[p] += 1; dels[p] += 1; } }
                    t = t.saturating_add(seq.len());
                }
                CsOp::Ins(mut seq) => { if t <= n { seq.make_ascii_uppercase(); *ins[t].entry(seq).or_insert(0) += 1; } }
                CsOp::Intron(len) => { t = t.saturating_add(len); }
            }
        }
    }
    // Amendment 15b (b): one copy per identical long insertion per consensus. The same biological piece lands at different columns on
    // different carriers (minimap2 places it relative to where the read starts), so without this rule the vote inserts it once per
    // column. Aggregated across columns, it is placed at the column with the most carriers (ties to the earliest) and removed elsewhere;
    // the count at the surviving column stays column-local (Amendment 15's majority test reads carriers against the members covering
    // THAT column).
    for p in 0..=n {
        let seqs: Vec<Vec<u8>> = ins[p].keys().filter(|s| s.len() >= STRUCT_MIN_INDEL).cloned().collect();
        for seq in seqs {
            let at = |q: usize| ins[q].get(&seq).copied().unwrap_or(0);
            let columns: Vec<usize> = (0..=n).filter(|&q| at(q) > 0).collect();
            if columns.len() > 1 {
                let keep = columns.iter().copied().max_by_key(|&q| (at(q), std::cmp::Reverse(q))).unwrap();
                for q in columns { if q != keep { ins[q].remove(&seq); } }
            }
        }
    }
    // spec §5.4 / R5 end trim: columns lo..=hi survive, the leading and trailing ones covered by fewer than END_MIN_COVER members do not; an
    // insertion goes with the column it precedes, so the one after the last column survives only with the last column. No covered column
    // (an empty template, no member) returns here, so n >= 1 below.
    let (Some(lo), Some(hi)) = (cover.iter().position(|&c| c >= END_MIN_COVER), cover.iter().rposition(|&c| c >= END_MIN_COVER)) else { return Ok(Vec::new()); };
    let last = if hi + 1 == n { n } else { hi };
    let mut out = Vec::with_capacity(hi - lo + 64);
    for p in lo..=last {
        // A13: the most frequent insertion of one size class before column p (the larger count, then the smaller sequence: a total order, no
        // hash order reaches the output); the long class first, the short class only when no long insertion passes Amendment 15's test
        let most = |long: bool| ins[p].iter().filter(|(s, _)| (s.len() >= STRUCT_MIN_INDEL) == long).max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0)));
        let covering = cover[p.min(n - 1)];
        // Amendment 15 (c): a long insertion needs the floor of 3 carriers AND a majority of the members covering the column
        if let Some((seq, _)) = most(true).filter(|&(_, &cnt)| cnt >= STRUCT_MIN_SUPPORT && 2 * cnt >= covering) {
            out.extend_from_slice(seq);
        } else if let Some((seq, _)) = most(false).filter(|&(_, &cnt)| covering >= VOTE_MIN_COVER && 2 * cnt >= covering) {
            out.extend_from_slice(seq);
        }
        if p == n { break; }
        if cover[p] >= VOTE_MIN_COVER && 2 * dels[p] >= cover[p] { continue; }
        let mut base = template[p].to_ascii_uppercase();
        if cover[p] >= VOTE_MIN_COVER {
            let same = cover[p].saturating_sub(subs[p].values().sum::<usize>() + dels[p]);
            if let Some((&b, &c)) = subs[p].iter().max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0))) { if c > same { base = b; } }
        }
        out.push(base);
    }
    Ok(out)
}

/// Prereg Amendment 15b (a): normalise one member's `cs` before the vote. minimap2 `splice:hq` writes an exon the template lacks as
/// `+X·E ~|X|` — the insertion begins with the template bases the adjacent skip removes (107 of 134 insertion-plus-skip pairs in the
/// defective clusters of docs/O3_CANDIDATES_CONSENSUS_DEFECT_2026-10-03.md), and the insertion's column then depends on where the read
/// starts. When an insertion is immediately followed by a skip over template bases `X` and the insertion begins with exactly those bases,
/// the pair becomes matches over `X` and the insertion shortened by them (`+X·E ~|X|` -> `:|X| +E`), so one biological event lands at one
/// column. Chained skips repeat; anything else passes through untouched. `ts` is the hit's start on the template.
fn normalize_cs_ops(ops: Vec<CsOp>, template: &[u8], ts: usize) -> Vec<CsOp> {
    let n = template.len();
    let (mut out, mut t) = (Vec::with_capacity(ops.len()), ts);
    let mut iter = ops.into_iter().peekable();
    while let Some(op) = iter.next() {
        match op {
            CsOp::Ins(seq) => {
                let mut seq = seq;
                while let Some(&CsOp::Intron(skip)) = iter.peek() {
                    let begins_with_the_skipped = seq.len() >= skip
                        && t + skip <= n
                        && seq[..skip].eq_ignore_ascii_case(&template[t..t + skip]);
                    if !begins_with_the_skipped { break; }
                    out.push(CsOp::Eq(skip));
                    t += skip;
                    seq = seq[skip..].to_vec();
                    iter.next();
                    if seq.is_empty() { break; }
                }
                if !seq.is_empty() { out.push(CsOp::Ins(seq)); }
            }
            other => {
                match &other {
                    CsOp::Eq(len) => t = t.saturating_add(*len),
                    CsOp::Sub(..) => t = t.saturating_add(1),
                    CsOp::Del(s) => t = t.saturating_add(s.len()),
                    CsOp::Intron(len) => t = t.saturating_add(*len),
                    CsOp::Ins(_) => {}
                }
                out.push(other);
            }
        }
    }
    out
}

/// The number of substitutions in a pairwise `cs`: spec §5.5's k, the columns that distinguish two consensus sequences. Indels are not counted.
pub fn distinguishing_columns(ops: &[CsOp]) -> usize { ops.iter().filter(|o| matches!(o, CsOp::Sub(..))).count() }

/// P(X >= n_small) for X ~ Binomial(n_small + n_large, eps^k), in log space with a log-factorial table; the variant is real iff it is < alpha
/// (spec §5.5: n_small = the reads of the smaller cluster, n_large = the larger's, k = `distinguishing_columns`). No distinguishing column or
/// no reads in the small cluster is never real; an `eps^k` that underflows to 0 gives P = 0: real.
pub fn variant_is_real(n_small: usize, n_large: usize, k: usize, eps: f64, alpha: f64) -> bool {
    if k == 0 || n_small == 0 { return false; }
    let n = n_small + n_large;
    let p = eps.powi(i32::try_from(k).unwrap_or(i32::MAX));
    let (lp, lq) = (p.ln(), (1.0 - p).ln());
    let mut lf = vec![0f64; n + 1];
    for i in 1..=n { lf[i] = lf[i - 1] + (i as f64).ln(); }
    let tail: f64 = (n_small..=n).map(|x| (lf[n] - lf[x] - lf[n - x] + x as f64 * lp + (n - x) as f64 * lq).exp()).sum();
    tail < alpha
}

/// One refinement pass over a cluster (spec §5.4's one polishing pass, then the split). `members` pairs each member read with its hit against
/// the consensus; `template` is that consensus. Neither it nor the read sequences are consulted (the rule reads the hits only; the parameters
/// keep the signature parallel to `consensus_from_template`). Returns the ascending indices of the members that fit the consensus
/// (`de <= delta` and `shorter_cov >= 0.5`) and of the rest; the caller makes the rest a new cluster when it holds >= `--min-cluster` reads
/// and places a member that has no hit at all.
pub fn refine_cluster(_template: &[u8], members: &[(&[u8], &PafHit)], delta: f64) -> (Vec<usize>, Vec<usize>) {
    (0..members.len()).partition(|&i| { let h = members[i].1; h.de <= delta && shorter_cov(h) >= 0.5 })
}

// ---- flag / link / merge fates, the flag floor and the exon-union representative (spec §5.6-§5.7) ------------------------------------

/// Identity x coverage (`id_cov`) from which a consensus counts as already in the reference (`bench/rna_allele/link_test.py:141`).
const IN_REFERENCE_ID_COV: f64 = 0.999;

/// One cluster consensus of a family, the unit the flag / link / merge chain judges: `id` is its sequence name in the consensus FASTA (the
/// `q` of its PAF hits), `n_reads` the reads behind it after refinement and merging, `seq` the consensus.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ClusterSeq { pub family: String, pub id: String, pub n_reads: usize, pub seq: Vec<u8> }

/// The chain's verdict on one consensus (spec §5.6; Amendments 7-9). Loci are `chrom:start-end` of the best genome hit in PAF coordinates
/// (0-based, half-open).
#[derive(Clone, Debug, PartialEq)]
pub enum Fate {
    /// Identity x coverage >= 0.999: the sequence is already in the reference, so it is no copy of anything and is dropped.
    InReference,
    /// Within delta of a reference locus (`d` = `whole_length_d` of the hit): an allele of that locus, never a candidate.
    Linked { locus: String, d: f64 },
    /// Beyond delta of every reference locus: a new copy, merged by `components` and kept by `is_flagged`. `nearest` is the best hit's
    /// locus, `"none"` (with d = 1.0) when the genome gave no hit.
    NewCopy { nearest: String, d: f64 },
}

/// The verdict on one consensus from its best genome hit (`best_by_id_cov`; `None` = no hit at all). Identity x coverage (`id_cov`) >= 0.999
/// is `InReference`, tested first as `contigs` of `bench/rna_allele/link_test.py` does; otherwise the whole-length divergence d
/// (`whole_length_d`) <= `delta` is `Linked` to the hit's locus, and anything else a `NewCopy` with that locus as `nearest`. No hit is a
/// `NewCopy` with `nearest = "none"` and d = 1.0 whatever `delta`. The cluster is not consulted (its hit carries every quantity used; the
/// parameter is the plan's interface).
pub fn classify(_c: &ClusterSeq, best_genome_hit: Option<&PafHit>, delta: f64) -> Fate {
    let Some(h) = best_genome_hit else { return Fate::NewCopy { nearest: "none".into(), d: 1.0 } };
    if id_cov(h) >= IN_REFERENCE_ID_COV { return Fate::InReference; }
    let (locus, d) = (format!("{}:{}-{}", h.t, h.ts, h.te), whole_length_d(h));
    if d <= delta { Fate::Linked { locus, d } } else { Fate::NewCopy { nearest: locus, d } }
}

/// The merge step (spec §5.6.4, Amendment 8): the components of a family's NEW-COPY consensus sequences `ids` under their all-vs-all hits
/// `ava`. The rule is the read clustering's, so this is `cluster_reads` (a pair joins on its best hit when that covers >= 50% of the shorter
/// sequence with `de <= delta`), with its ordering: ascending index lists by (size descending, first index).
pub fn components(ids: &[String], ava: &[PafHit], delta: f64) -> Vec<Vec<usize>> { cluster_reads(ids, ava, delta) }

/// The flag floor (ruling R1): a component of new-copy clusters is a flagged candidate iff the reads behind its clusters sum to
/// >= `min_support` (`--min-support`, 6 = twice IsoCon's 3-read transcript minimum). It counts reads, not clusters: at delta a copy's
/// isoforms merge into one cluster, so the cluster count is not IsoCon's transcript count. An empty component sums to 0.
pub fn is_flagged(component_clusters: &[&ClusterSeq], min_support: usize) -> bool {
    component_clusters.iter().map(|c| c.n_reads).sum::<usize>() >= min_support
}

/// What `union_sequence_with_note` set aside, counted per member after the backbone.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct UnionNote {
    /// Members for which the closure found no alignment to the current union: nothing is taken from them.
    pub no_hit: usize,
    /// Members whose best hit is not on the `+` strand: skipped whole (see `union_sequence`).
    pub skipped_minus: usize,
}

/// The exon-union representative of a component (spec §5.7). `members_longest_first` are the component's consensus sequences, the longest
/// first (the caller sorts; lengths are not looked at). The backbone is the first member (upper-cased); each later member is passed with the
/// CURRENT union to `hits_vs_current`, which returns its best hit with the member as query and the union as target, so an exon that two
/// members share is taken once: the second one finds it already there. That hit must skip a long stretch of the union in ONE alignment (a
/// member that skips an exon): minimap2 `-x splice:hq -uf -c --cs` does (the skip is a `~` intron of the `cs`); `-x asm20 -c --cs` does not,
/// a skip of some hundred bp costs more than the shorter side earns, the alignment is cut there, and what lies beyond the cut comes back as an
/// 'unaligned' prefix or suffix and is inserted again (ruling R6, spec §9b, records the measurement behind the splice preset). From the
/// hit, in member order:
/// * the unaligned prefix `member[..qs]`, when `qs >= 20`, goes in front of target position `ts`, and the unaligned suffix `member[qe..]`,
///   when `qlen - qe >= 20`, in front of position `te` (so with `ts == 0` / `te == tlen` they are prepended / appended);
/// * every insertion of the `cs` of >= 20 bp goes in front of the target base it precedes. An insertion < 20 bp is an error or a microindel,
///   not an exon (the 20 bp of the consensus' R2); deletions (the union has what the member lacks) and substitutions (the backbone's allele
///   stays) change nothing.
///
/// The pieces of one member are all placed against the union as it was when the member was aligned: they are applied from the highest target
/// position down, so the offsets of those still to be applied stay valid, and pieces at one position keep their member order. Everything
/// inserted is upper case (`cs` writes lower case).
///
/// A member with no hit adds nothing. Nor does one whose hit is not on the `+` strand (ruled): inserted `cs` bases already lie along the
/// target, but the unaligned ends of a `-` member would have to be reverse-complemented, which is not done. `union_sequence_with_note` counts
/// both kinds. A hit without a `cs` contributes its unaligned ends only. A hit that does not fit its sequences (query or target length,
/// spans, a `cs` that does not end at `te`) and a malformed `cs` are errors naming the member.
///
/// What one linear sequence cannot express is exon ORDER. Members that order exons inconsistently (an exon the union holds before another,
/// a later member holds after it; this also arises when an insertion and a deletion sit side by side in a `cs` and the aligner chose their
/// order) cannot all be read from the union: the out-of-order exon is inserted again, a duplicate. Alternative first (last) exons end up side
/// by side at `ts` (`te`). A piece that the aligner left unaligned through divergence rather than structure is inserted although the union
/// has it, and a long unaligned end of a member (junk, a primer remnant) is taken like an exon: the rule is length only.
pub fn union_sequence(members_longest_first: &[Vec<u8>], hits_vs_current: impl FnMut(&[u8], &[u8]) -> Option<PafHit>) -> anyhow::Result<Vec<u8>> {
    union_sequence_with_note(members_longest_first, hits_vs_current).map(|(union, _)| union)
}

/// `union_sequence` with the count of the members it set aside (no hit; hit not on the `+` strand).
pub fn union_sequence_with_note(
    members_longest_first: &[Vec<u8>],
    mut hits_vs_current: impl FnMut(&[u8], &[u8]) -> Option<PafHit>,
) -> anyhow::Result<(Vec<u8>, UnionNote)> {
    let mut note = UnionNote::default();
    let Some((backbone, rest)) = members_longest_first.split_first() else { return Ok((Vec::new(), note)) };
    let mut union = backbone.to_ascii_uppercase();
    for (k, member) in rest.iter().enumerate() {
        let Some(hit) = hits_vs_current(member.as_slice(), union.as_slice()) else { note.no_hit += 1; continue };
        if hit.strand != b'+' { note.skipped_minus += 1; continue; }
        let pieces = union_pieces(member, &union, &hit).map_err(|e| anyhow::anyhow!("member {}: {e}", k + 1))?;
        // the pieces are in member order, so their positions never decrease: taken from the last one back, each goes in at the highest position
        // still to do (the offsets of the pieces before it stay valid), and pieces at one position are applied last-first, so that they end up in member order
        for (pos, seq) in pieces.into_iter().rev() { union.splice(pos..pos, seq); }
    }
    Ok((union, note))
}

/// The pieces of `member` that the union lacks (see `union_sequence`), in member order, each with the union position it goes in front of;
/// the positions never decrease along the list. An error when the hit does not fit the two sequences or its `cs` is malformed.
fn union_pieces(member: &[u8], union: &[u8], h: &PafHit) -> anyhow::Result<Vec<(usize, Vec<u8>)>> {
    anyhow::ensure!(h.qlen == member.len(), "the hit's query length {} is not the member's {} bp", h.qlen, member.len());
    anyhow::ensure!(h.tlen == union.len(), "the hit's target length {} is not the current union's {} bp", h.tlen, union.len());
    anyhow::ensure!(h.qs <= h.qe && h.qe <= h.qlen, "query span {}..{} lies outside the {} bp member", h.qs, h.qe, h.qlen);
    anyhow::ensure!(h.ts <= h.te && h.te <= h.tlen, "target span {}..{} lies outside the {} bp union", h.ts, h.te, h.tlen);
    let mut pieces = Vec::new();
    if h.qs >= STRUCT_MIN_INDEL { pieces.push((h.ts, member[..h.qs].to_ascii_uppercase())); }
    if let Some(cs) = h.cs.as_deref() {
        let mut t = h.ts;                                                   // the target position the cs has reached
        for op in parse_cs(cs)? {
            match op {
                CsOp::Eq(n) | CsOp::Intron(n) => t = t.saturating_add(n),
                CsOp::Sub(..) => t = t.saturating_add(1),
                CsOp::Del(s) => t = t.saturating_add(s.len()),
                CsOp::Ins(mut s) => if s.len() >= STRUCT_MIN_INDEL { s.make_ascii_uppercase(); pieces.push((t, s)); },
            }
        }
        anyhow::ensure!(t == h.te, "the cs ends at target position {t}, not at te = {}", h.te);
    }
    if h.qlen - h.qe >= STRUCT_MIN_INDEL { pieces.push((h.te, member[h.qe..].to_ascii_uppercase())); }
    // ts <= every position of the cs walk <= te (the walk ends at te): the order that `union_sequence_with_note` relies on
    debug_assert!(pieces.windows(2).all(|w| w[0].0 <= w[1].0), "pieces out of order: {:?}", pieces.iter().map(|p| p.0).collect::<Vec<_>>());
    Ok(pieces)
}

// ---- the minimap2 runner (cached through run_cache) and the output writers (spec §5.8; plan task 7) ---------------------------------------

/// All-vs-all, one direction per pair (`--dual=no`; the self hits it lets through are ignored by `cluster_reads`) with up to 100 secondary
/// hits per sequence (`-N 100 -p 0.1`): the reads of a net (read clustering) and the consensus sequences of a family (the significance merge,
/// the components of new-copy consensus sequences). `--cs` carries the substitution columns.
/// Ruling R11: `--dual=no`, not `-X`: with `-X` minimap2 keeps every chain and ignores `-N`/`-p`.
pub const MM2_AVA: &[&str] = &["-x", "asm20", "-c", "--cs", "--dual=no", "-N", "100", "-p", "0.1", "--secondary=yes"];
/// Members against their template (the consensus vote) and against their cluster's consensus (`refine_cluster`). A splice preset (prereg
/// Amendment 13), the same arguments as `MM2_UNION`: a member that skips an exon of its template is aligned across the skip in ONE alignment
/// (the skip is a `~` of its `cs`, which neither covers nor votes), where `asm20` cut the alignment at the skip (ruling R6's measurement) and
/// the member then voted on one side of it only.
pub const MM2_MEMBERS: &[&str] = &["-x", "splice:hq", "-uf", "-c", "--cs", "-N", "5", "-p", "0.5"];
/// A member against the CURRENT union (`union_sequence`; ruling R6, spec §9b, records the measurement that chose it). A splice preset: it
/// skips an exon the member lacks in ONE alignment, where `asm20` cuts the alignment at the skip and the far side comes back as an unaligned end
/// and is inserted again.
pub const MM2_UNION: &[&str] = &["-x", "splice:hq", "-uf", "-c", "--cs", "-N", "5", "-p", "0.5"];
/// Consensus sequences against the primary genome's splice index (`--index`): the genome hits that `classify` judges (spec §5.6).
pub const MM2_GENOME: &[&str] = &["-x", "splice:hq", "-uf", "-c", "--eqx", "-N", "20"];
/// The attribution set (the unmapped reads >= `MIN_UNMAPPED_LEN` and the poorly placed reads: the query) against this run's net reads and
/// every family's copy sequences (`--copies-fa`: the target), in one call per run (prereg Amendment 13b): the hits `attribute_by_hits`
/// judges. `map-ont`: read against read, the family definition's edge alignment in read space; up to 5 secondary hits within half the best
/// score; no `--cs` (the rule reads matches, spans and `de` only).
pub const MM2_ATTRIB: &[&str] = &["-x", "map-ont", "-c", "-N", "5", "-p", "0.5"];

/// The name of the k-th candidate of `family` (k = the component's order within the family): `cand_<family>_<k>`.
pub fn candidate_id(family: &str, k: usize) -> String { format!("cand_{family}_{k}") }

/// Runs `minimap2 <args> -t <threads> <target> <query>` and leaves its PAF at `out_paf`; `target` may be a `.mmi` index. The binary is
/// `RUSTLE_MINIMAP2` (default `minimap2`) and its stderr is discarded. A run that fails is an error naming the whole command (an unstartable
/// binary names `RUSTLE_MINIMAP2` too) and leaves no stale or partial `out_paf` behind.
///
/// With `cache = Some(root)` (the run_cache root, `run_cache::cache_root()`) the PAF is a pinned `paf` entry keyed on the command line, the
/// minimap2 build and the content hash of every byte of `target` and `query` (`paf_key`; paths and mtimes do not count): a hit replays the
/// cached PAF as `out_paf` (a hard link, else a copy), a miss runs minimap2 and links the product into a new entry; only a successful run is
/// stored. `-t <threads>` is appended after the key is made (ruling R10): the thread count never enters it, so a run with other threads
/// replays the same entry (minimap2's output does not depend on it). Every call reads both inputs in full to hash them, which for a multi-GB
/// `.mmi` takes about as long as loading it: key such a target with `minimap2_keyed`. An input that cannot be read is an error naming it.
/// `out_paf` is unlinked before it is written, never truncated: it may be a hard link to a cache payload.
pub fn minimap2(args: &[&str], target: &Path, query: &Path, out_paf: &Path, cache: Option<&Path>, threads: usize) -> anyhow::Result<()> {
    run_minimap2(&minimap2_binary(), args, target, None, query, out_paf, cache, threads)
}

/// `minimap2` with the TARGET named in the cache key by `target_key` instead of the content hash of its bytes (ruling R7): the genome's
/// splice index (`--index`, 12-14 GB) is keyed by `run_cache::file_fingerprint` (canonical path, size, mtime), as the binaries key their
/// BAM and FASTA, and is never read for the key. The query is still hashed in full.
pub fn minimap2_keyed(args: &[&str], target: &Path, target_key: &str, query: &Path, out_paf: &Path, cache: Option<&Path>, threads: usize) -> anyhow::Result<()> {
    run_minimap2(&minimap2_binary(), args, target, Some(target_key), query, out_paf, cache, threads)
}

/// `RUSTLE_MINIMAP2`, else `minimap2`.
pub fn minimap2_binary() -> String { std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string()) }

/// `minimap2` / `minimap2_keyed` with the binary given (the tests' seam: no process-wide variable to set).
#[allow(clippy::too_many_arguments)]
fn run_minimap2(mm2: &str, args: &[&str], target: &Path, target_key: Option<&str>, query: &Path, out_paf: &Path, cache: Option<&Path>, threads: usize) -> anyhow::Result<()> {
    let entry = cache.map(|root| paf_entry(root, mm2, args, target, target_key, query)).transpose()?;
    if let Some(e) = entry.as_ref().filter(|e| e.is_hit()) {
        // a replay that fails (another run replacing the entry) falls through to running minimap2
        if let Ok(linked) = e.replay("out.paf", out_paf) {
            eprintln!("[cache] o3_candidates: minimap2 PAF replayed from {} ({}; minimap2 skipped)", e.dir.display(), if linked { "hard link" } else { "copy" });
            return Ok(());
        }
    }
    // R10: the thread count joins the command only now, after the key
    let threads = threads.max(1).to_string();
    let shown = format!("{mm2} {} -t {threads} {} {}", args.join(" "), target.display(), query.display());
    unlink_if_present(out_paf)?;
    let ran = std::fs::File::create(out_paf).with_context(|| format!("creating {}", out_paf.display())).and_then(|paf| {
        let status = std::process::Command::new(mm2).args(args).args(["-t", threads.as_str()]).arg(target).arg(query)
            .stdout(paf).stderr(std::process::Stdio::null()).status()
            .with_context(|| format!("running `{shown}` (RUSTLE_MINIMAP2 names the minimap2 binary)"))?;
        anyhow::ensure!(status.success(), "`{shown}` failed ({status})");
        Ok(())
    });
    if let Err(e) = ran {
        let _ = std::fs::remove_file(out_paf);
        return Err(e);
    }
    if let Some(e) = entry.as_ref() {
        let stored = e.staging().and_then(|st| { e.stage_link(&st, "out.paf", out_paf)?; e.commit(&st) });
        if let Err(err) = stored { eprintln!("[cache] could not store the PAF ({err:#}); continuing"); }
    }
    Ok(())
}

/// How the cache key of one minimap2 call names its target.
enum TargetId {
    /// The content hash of every byte of the target file (the default).
    Content(rc::ContentHash),
    /// A caller-given key, e.g. `run_cache::file_fingerprint` of a multi-GB `.mmi` (`minimap2_keyed`; ruling R7).
    Key(String),
}

/// The pinned `paf` cache entry of one call: its key (`paf_key`) hashes the query in full, and the target too unless `target_key` names it.
fn paf_entry(root: &Path, mm2: &str, args: &[&str], target: &Path, target_key: Option<&str>, query: &Path) -> anyhow::Result<rc::Entry> {
    let hash = |p: &Path| rc::ContentHash::of_file(p).with_context(|| format!("hashing {} for the minimap2 cache key", p.display()));
    let target = match target_key { Some(k) => TargetId::Key(k.to_string()), None => TargetId::Content(hash(target)?) };
    let key = paf_key(&format!("{mm2} {}", args.join(" ")), &rc::minimap2_version(mm2), &target, &hash(query)?);
    Ok(rc::Entry::new(root, "paf", key).pinned())
}

/// The key text of one minimap2 call (`paf/<fnv(key)>/key.tsv`, which a hit must equal byte for byte): the command line (binary and
/// arguments without `-t`; not the paths, the content stands for them), the minimap2 build, the target (content hash and byte length, or the
/// caller's key under its own label) and the query (content hash and byte length), in that order, so the two files in each other's roles are
/// another key.
fn paf_key(cmd: &str, minimap2_version: &str, target: &TargetId, query: &rc::ContentHash) -> String {
    let target = match target {
        TargetId::Content(h) => format!("target_hash\tcontent128:{}\ntarget_bytes\t{}\n", h.hex(), h.len()),
        TargetId::Key(k) => format!("target_key\t{k}\n"),
    };
    format!("rustle o3 minimap2 v1\ncmd\t{cmd}\nminimap2\t{minimap2_version}\n{target}query_hash\tcontent128:{}\nquery_bytes\t{}\n", query.hex(), query.len())
}

/// One candidate copy of a family: a component of its new-copy cluster consensus sequences (`components`), judged by `is_flagged` and
/// represented by the exon union (`union_sequence`). The caller builds it; `write_outputs` writes it.
#[derive(Clone, Debug, PartialEq)]
pub struct Candidate {
    pub family: String,
    /// The candidate's name `cand_<family>_<k>` (`candidate_id`): its contig's FASTA name, the `candidate` column of both tables.
    pub id: String,
    /// The component's clusters, all of this family (`write_outputs` refuses another family's).
    pub clusters: Vec<ClusterSeq>,
    /// The exon-union representative; written as the contig when `flagged`.
    pub union: Vec<u8>,
    /// `is_flagged` of the component's clusters (reads >= `--min-support`), computed by the caller: the `flagged` column, and only a flagged
    /// candidate gets a contig.
    pub flagged: bool,
    /// The nearest reference locus (`chrom:start-end` of the best genome hit, `"none"` without one) and the whole-length divergence d to it
    /// (`Fate::NewCopy`'s fields). `clusters.tsv` shows this d on each of the candidate's member rows too: a member's own d is not carried.
    pub nearest: String,
    pub d: f64,
    /// The family's net: its reads before the `--max-reads` cap and the reads used after it (a capped net is sampled, never truncated silently:
    /// spec §8). PER FAMILY, so the caller sets the same two numbers on every candidate of the family and each row repeats them.
    pub n_net: usize,
    pub n_used: usize,
}

const CANDIDATES_HEADER: &str = "family\tcandidate\tn_clusters\tn_reads\tflagged\tunion_len\tnearest_locus\td\tn_net\tn_used";
const CLUSTERS_HEADER: &str = "family\tcluster\tcandidate\tn_reads\tconsensus_len\tfate\tlinked_to\td";

/// Writes the stage's four products (spec §4, §5.8), every float with 5 decimals, deterministic: rows are sorted by family, then candidate
/// (names ending in a number sort by it: `cand_F_2` before `cand_F_10`, `MCL2` before `MCL10`), whatever the input order.
/// * `<prefix>.candidates.tsv` (`family candidate n_clusters n_reads flagged union_len nearest_locus d n_net n_used`): one row per candidate,
///   flagged or not; `n_reads` sums its clusters' reads, `flagged` is `1` / `0` (the caller's `is_flagged`: reads >= `--min-support`).
/// * `<prefix>.clusters.tsv` (`family cluster candidate n_reads consensus_len fate linked_to d`): one row per cluster. A candidate's clusters
///   have its id, fate `new_copy`, `linked_to` `-` and the candidate's d; the `linked` clusters (cluster, locus, d) have candidate `-`, fate
///   `linked`, their locus and their own d and come last in their family. `InReference` clusters are not written (the caller drops them).
/// * `<prefix>.contigs.fa`: `>cand_<family>_<k>` and the union on one line, for the flagged candidates only.
/// * `<prefix>.nets.fa`: `>name` and the read on one line for every read of `nets_for_patch` (family, reads), in the given order: the input
///   of the patch realignment, whose caller picks the families (those with a flagged candidate).
///
/// A flagged candidate with an empty union and a cluster of another family than its candidate's are refused before anything is written. Each
/// file is unlinked and created anew, never truncated in place (run_cache replays products by hard link).
pub fn write_outputs(prefix: &str, cands: &[Candidate], linked: &[(ClusterSeq, String, f64)], nets_for_patch: &[(String, Vec<(String, Vec<u8>)>)]) -> anyhow::Result<()> {
    for c in cands {
        anyhow::ensure!(!c.flagged || !c.union.is_empty(), "candidate {} is flagged but its union sequence is empty: there is no contig to write", c.id);
        if let Some(k) = c.clusters.iter().find(|k| k.family != c.family) {
            anyhow::bail!("candidate {} (family {}) holds cluster {} of family {}", c.id, c.family, k.id, k.family);
        }
    }
    let mut sorted: Vec<&Candidate> = cands.iter().collect();
    sorted.sort_by(|a, b| name_order(&a.family, &b.family).then_with(|| name_order(&a.id, &b.id)));
    write_file(&format!("{prefix}.candidates.tsv"), |w| {
        writeln!(w, "{CANDIDATES_HEADER}")?;
        for c in &sorted {
            let n_reads: usize = c.clusters.iter().map(|k| k.n_reads).sum();
            writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.5}\t{}\t{}", c.family, c.id, c.clusters.len(), n_reads, c.flagged as u8, c.union.len(), c.nearest, c.d, c.n_net, c.n_used)?;
        }
        Ok(())
    })?;
    // `locus` is Some for a linked cluster
    struct Row<'a> { family: &'a str, cluster: &'a str, candidate: &'a str, n_reads: usize, len: usize, locus: Option<&'a str>, d: f64 }
    let mut rows: Vec<Row> = Vec::new();
    for c in &sorted {
        for k in &c.clusters { rows.push(Row { family: &c.family, cluster: &k.id, candidate: &c.id, n_reads: k.n_reads, len: k.seq.len(), locus: None, d: c.d }); }
    }
    for (k, locus, d) in linked { rows.push(Row { family: &k.family, cluster: &k.id, candidate: "-", n_reads: k.n_reads, len: k.seq.len(), locus: Some(locus), d: *d }); }
    rows.sort_by(|a, b| name_order(a.family, b.family)
        .then(a.locus.is_some().cmp(&b.locus.is_some()))
        .then_with(|| name_order(a.candidate, b.candidate))
        .then_with(|| name_order(a.cluster, b.cluster)));
    write_file(&format!("{prefix}.clusters.tsv"), |w| {
        writeln!(w, "{CLUSTERS_HEADER}")?;
        for r in &rows {
            writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.5}", r.family, r.cluster, r.candidate, r.n_reads, r.len, if r.locus.is_some() { "linked" } else { "new_copy" }, r.locus.unwrap_or("-"), r.d)?;
        }
        Ok(())
    })?;
    write_file(&format!("{prefix}.contigs.fa"), |w| {
        for c in sorted.iter().filter(|c| c.flagged) {
            writeln!(w, ">{}", c.id)?;
            w.write_all(&c.union)?;
            writeln!(w)?;
        }
        Ok(())
    })?;
    write_file(&format!("{prefix}.nets.fa"), |w| {
        for (name, seq) in nets_for_patch.iter().flat_map(|(_family, reads)| reads) {
            writeln!(w, ">{name}")?;
            w.write_all(seq)?;
            writeln!(w)?;
        }
        Ok(())
    })
}

/// The `--max-reads` cap of one family's net (spec §5.1.4): the names (distinct) sorted, and when there are more than `cap` of them a sample
/// of `cap`: the sorted names shuffled by Fisher-Yates driven by splitmix64 seeded with 1, the first `cap` kept. Returned sorted. A function of
/// the SET of names and `cap` only (no hash order, no clock, no thread), so the same on every run and machine.
pub fn sample_net(names: &[String], cap: usize) -> Vec<String> {
    let mut v = names.to_vec();
    v.sort_unstable();
    if v.len() > cap {
        let mut state = 1u64;                                               // splitmix64 seed 1: outputs mix(1), mix(1 + g), ...
        let mut next = || { let r = mix(state); state = state.wrapping_add(0x9E37_79B9_7F4A_7C15); r };
        for i in (1..v.len()).rev() { let j = (next() % (i as u64 + 1)) as usize; v.swap(i, j); }
        v.truncate(cap);
        v.sort_unstable();
    }
    v
}

const FAMILIES_HEADER: &str = "family\tn_net\tn_used\tn_clusters\tn_in_reference\tn_linked\tn_new\tn_candidates\tn_flagged";

/// One family's row of `<prefix>.families.tsv` (ruling R8): every family the stage was given gets one, zeros where nothing happened, so a
/// family without a candidate is visible with where its reads went (no net, no cluster, all clusters already in the reference or linked).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FamilyCounts {
    pub family: String,
    /// The net before the `--max-reads` cap and the reads used after it (the numbers `Candidate` carries).
    pub n_net: usize,
    pub n_used: usize,
    /// The clusters the chain judged (after the `--min-cluster` floor, the refinement and the significance merge) and their fates.
    pub n_clusters: usize,
    pub n_in_reference: usize,
    pub n_linked: usize,
    pub n_new: usize,
    /// The components of the new-copy clusters (the candidates) and the flagged ones among them.
    pub n_candidates: usize,
    pub n_flagged: usize,
}

/// Writes `<prefix>.families.tsv` (`family n_net n_used n_clusters n_in_reference n_linked n_new n_candidates n_flagged`), one row per
/// given family, sorted by family as `write_outputs` sorts (`MCL2` before `MCL10`); the file is unlinked and created anew.
pub fn write_family_table(prefix: &str, rows: &[FamilyCounts]) -> anyhow::Result<()> {
    let mut sorted: Vec<&FamilyCounts> = rows.iter().collect();
    sorted.sort_by(|a, b| name_order(&a.family, &b.family));
    write_file(&format!("{prefix}.families.tsv"), |w| {
        writeln!(w, "{FAMILIES_HEADER}")?;
        for r in sorted {
            writeln!(w, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", r.family, r.n_net, r.n_used, r.n_clusters, r.n_in_reference, r.n_linked, r.n_new, r.n_candidates, r.n_flagged)?;
        }
        Ok(())
    })
}

const READS_HEADER: &str = "read\tfamily\tcluster";

/// The members and consensus sequences of the clusters written to `clusters.tsv` (the caller passes those: the new-copy and the linked ones):
/// `<prefix>.reads.tsv` (`read family cluster`, spec §4: read -> cluster, diagnostics that O2 does not read) and `<prefix>.clusters.fa`
/// (`>cluster` and its consensus on one line), the inputs of the representative measure (Amendment 12, A12-2: each read against its union
/// and its cluster consensus). Clusters in `clusters.tsv` order (family, then cluster id, `name_order`), reads sorted within a cluster; each
/// file unlinked and created anew.
pub fn write_cluster_members(prefix: &str, clusters: &[(&ClusterSeq, &[String])]) -> anyhow::Result<()> {
    let mut sorted: Vec<&(&ClusterSeq, &[String])> = clusters.iter().collect();
    sorted.sort_by(|a, b| name_order(&a.0.family, &b.0.family).then_with(|| name_order(&a.0.id, &b.0.id)));
    write_file(&format!("{prefix}.reads.tsv"), |w| {
        writeln!(w, "{READS_HEADER}")?;
        for (c, reads) in &sorted {
            let mut names: Vec<&String> = reads.iter().collect();
            names.sort();
            for r in names { writeln!(w, "{r}\t{}\t{}", c.family, c.id)?; }
        }
        Ok(())
    })?;
    write_file(&format!("{prefix}.clusters.fa"), |w| {
        for (c, _) in &sorted {
            writeln!(w, ">{}", c.id)?;
            w.write_all(&c.seq)?;
            writeln!(w)?;
        }
        Ok(())
    })
}

/// The order of names that end in a number: the text before the trailing digits, then the number (by digit count, then digits: it cannot
/// overflow), then the whole name. `MCL2` < `MCL10`; `cand_F_2` < `cand_F_10`.
fn name_order(a: &str, b: &str) -> std::cmp::Ordering {
    fn key(s: &str) -> (&str, usize, &str) {
        let stem = s.trim_end_matches(|c: char| c.is_ascii_digit());
        let digits = s[stem.len()..].trim_start_matches('0');
        (stem, digits.len(), digits)
    }
    key(a).cmp(&key(b)).then_with(|| a.cmp(b))
}

/// Unlinks `path` when it exists. A product is unlinked before it is rewritten, never truncated in place: run_cache replays products by
/// hard link, so an old product may share its inode with a cache payload.
fn unlink_if_present(path: &Path) -> anyhow::Result<()> {
    match std::fs::remove_file(path) {
        Ok(()) => Ok(()),
        Err(e) if e.kind() == std::io::ErrorKind::NotFound => Ok(()),
        Err(e) => Err(e).with_context(|| format!("removing {}", path.display())),
    }
}

/// A new file at `path` (any old one unlinked first) holding what `fill` writes, flushed; an error names the path.
fn write_file(path: &str, fill: impl FnOnce(&mut dyn Write) -> std::io::Result<()>) -> anyhow::Result<()> {
    unlink_if_present(Path::new(path))?;
    let mut w = std::io::BufWriter::new(std::fs::File::create(path).with_context(|| format!("creating {path}"))?);
    fill(&mut w).and_then(|()| w.flush()).with_context(|| format!("writing {path}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rand_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed; (0..n).map(|_| { x ^= x << 13; x ^= x >> 7; x ^= x << 17; b"ACGT"[(x % 4) as usize] }).collect()
    }
    fn mutate(s: &[u8], rate_per_kb: usize, seed: u64) -> Vec<u8> {
        // substitutes every (1000 / rate_per_kb)-th base; the new base is always different from the original
        let mut v = s.to_vec(); let mut x = seed.max(1);
        for i in (0..v.len()).step_by(1000 / rate_per_kb.max(1)) {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            let orig = b"ACGT".iter().position(|&c| c == v[i]).expect("ACGT input");
            v[i] = b"ACGT"[(orig + 1 + (x % 3) as usize) % 4];   // index + 1..=3 (mod 4) never returns the original index
        }
        v
    }
    #[test]
    fn canonical_kmers_are_strand_free() {
        let s = rand_seq(500, 7);
        let rc = crate::vg_family::seq_utils::reverse_complement(&s);
        let mut a = canonical_kmers(&s, 31); let mut b = canonical_kmers(&rc, 31);
        a.sort(); b.sort();
        assert_eq!(a, b);
    }
    #[test]
    fn sketch_share_separates_copies_from_errors() {
        let copy_a = rand_seq(3000, 1);
        let same_copy_read = mutate(&copy_a, 2, 3);        // 0.2%: HiFi-like (every 500th base: 6 substitutions in 3 kb)
        let copy_b = mutate(&copy_a, 25, 5);               // 2.5% diverged paralog (every 40th base: 75 substitutions in 3 kb)
        let real_subs = |x: &[u8]| x.iter().zip(&copy_a).filter(|(p, q)| p != q).count();
        assert_eq!((real_subs(&same_copy_read), real_subs(&copy_b)), (6, 75), "every substitution must be real: exactly 0.2% and 2.5% of 3000 bp");
        let sa = minimizer_sketch(&copy_a, 31, 5);
        assert!(sketch_share(&sa, &minimizer_sketch(&same_copy_read, 31, 5)) > 0.8);
        assert!(sketch_share(&sa, &minimizer_sketch(&copy_b, 31, 5)) < 0.6);
    }

    fn paf_hit(q: &str, t: &str, matches: usize) -> PafHit {
        PafHit { q: q.into(), qlen: 1000, qs: 0, qe: 1000, strand: b'+', t: t.into(), tlen: 1000, ts: 0, te: 1000, matches, block: 1000, de: 0.01, cs: None }
    }
    #[test]
    fn paf_quantities_match_the_chain() {
        let line = "r1\t1000\t10\t990\t+\tc1\t5000\t100\t1100\t960\t1000\t60\tNM:i:40\tde:f:0.0150\tcs:Z::500*ac:300+tt:179-g:0";
        let h = &parse_paf(line)[0];
        assert!((whole_length_d(h) - 0.04).abs() < 1e-9);          // 1 - 960/1000
        assert!((id_cov(h) - 0.96 * 0.98).abs() < 1e-9);
        assert!((h.de - 0.015).abs() < 1e-9);
        let ops = parse_cs(h.cs.as_deref().unwrap()).unwrap();
        assert!(matches!(ops[1], CsOp::Sub(b'a', b'c')));
        assert!(matches!(&ops[3], CsOp::Ins(v) if v == b"tt"));
        assert!(matches!(&ops[5], CsOp::Del(v) if v == b"g"));
    }
    #[test]
    fn parse_paf_reads_many_lines_and_defaults_absent_tags() {
        // a CRLF-terminated '-' hit with both tags, a blank line, a hit with another tag only, three unusable lines (truncated to 4
        // columns, strand '.', non-numeric qlen: all skipped), and a hit with exactly the 12 fixed columns and no tag at all
        let text = [
            "a\t120\t3\t117\t-\tb\t260\t50\t170\t105\t118\t60\tde:f:0.01\tcs:Z::100\r",
            "",
            "c\t300\t10\t290\t+\td\t350\t5\t285\t270\t281\t0\ttp:A:P",
            "e\t100\t0\t5",
            "x\t100\t0\t100\t.\ty\t100\t0\t100\t100\t100\t60\tde:f:0.0",
            "z\tNA\t0\t100\t+\ty\t100\t0\t100\t100\t100\t60\tde:f:0.0",
            "g\t50\t1\t49\t+\th\t60\t2\t52\t40\t48\t30",
        ].join("\n") + "\n";
        let hits = parse_paf(&text);
        assert_eq!(hits.len(), 3);
        let (a, c, g) = (&hits[0], &hits[1], &hits[2]);
        // every fixed column lands in its own field (all values distinct so a swapped column cannot pass)
        assert_eq!((a.q.as_str(), a.qlen, a.qs, a.qe, a.strand, a.t.as_str(), a.tlen, a.ts, a.te, a.matches, a.block), ("a", 120, 3, 117, b'-', "b", 260, 50, 170, 105, 118));
        assert!((a.de - 0.01).abs() < 1e-12);
        assert_eq!(a.cs.as_deref(), Some(":100"));              // the CR of CRLF is not part of the cs string
        assert_eq!((c.q.as_str(), c.qlen, c.qs, c.qe, c.strand, c.t.as_str(), c.tlen, c.ts, c.te, c.matches, c.block), ("c", 300, 10, 290, b'+', "d", 350, 5, 285, 270, 281));
        assert_eq!((c.de, c.cs.as_deref()), (1.0, None));       // absent de -> 1.0 (never clusters, never merges), absent cs -> None
        assert_eq!((g.q.as_str(), g.t.as_str(), g.matches, g.de, g.cs.as_deref()), ("g", "h", 40, 1.0, None));
        assert!(parse_paf("").is_empty() && parse_paf("\n\n").is_empty());
    }
    #[test]
    fn parse_cs_walks_every_operation_and_rejects_malformed_strings() {
        let ops = parse_cs(":10*ag~gt1500ag+acg-t:5*nc").unwrap();
        assert_eq!(ops, vec![CsOp::Eq(10), CsOp::Sub(b'a', b'g'), CsOp::Intron(1500), CsOp::Ins(b"acg".to_vec()), CsOp::Del(b"t".to_vec()), CsOp::Eq(5), CsOp::Sub(b'n', b'c')]);
        assert_eq!(parse_cs("").unwrap(), vec![]);
        assert_eq!(parse_cs(":0").unwrap(), vec![CsOp::Eq(0)]);                       // a zero-length match is kept, not dropped
        assert_eq!(parse_cs("*AG+ACG").unwrap(), vec![CsOp::Sub(b'A', b'G'), CsOp::Ins(b"ACG".to_vec())]);   // either case, bytes kept as written
        // every malformed string is an error, including ones whose prefix is a valid operation list (never a partial result)
        for bad in [":", ":x", "*a", "*ac*", "*a-", "+", "-:3", "~", "~gt", "~gt15", "~gt15a", "~gt15ag7", "~g1ag", "~1215ag", "~gtag", ":10;", "10", " :3", "=ACGT", ":5\n", "*aé", ":99999999999999999999999"] {
            assert!(parse_cs(bad).is_err(), "{bad:?} must be rejected");
        }
        let msg = parse_cs(":10;").unwrap_err().to_string();
        assert!(msg.contains("byte 3") && msg.contains(":10;"), "the error names the offset and the neighbourhood: {msg}");
    }
    #[test]
    fn best_by_matches_keeps_the_hit_with_most_matches_per_query() {
        let hits = vec![
            paf_hit("r1", "t_low", 900), paf_hit("r2", "t_only", 500), paf_hit("r1", "t_high", 950), paf_hit("r1", "t_later_low", 940),
            paf_hit("r3", "t_first", 700), paf_hit("r3", "t_tie", 700),
        ];
        let best = best_by_matches(&hits);
        assert_eq!(best.len(), 3);
        assert_eq!(best["r1"].t, "t_high");                     // the higher-matches hit wins whatever its position
        assert_eq!(best["r2"].t, "t_only");
        assert_eq!(best["r3"].t, "t_first");                    // a tie keeps the first encountered
        assert!(best_by_matches(&[]).is_empty());
    }
    #[test]
    fn shorter_cov_uses_the_span_on_the_shorter_sequence() {
        let mut h = paf_hit("q", "t", 700);
        (h.qlen, h.qs, h.qe) = (3000, 100, 1000);               // query: 900 aligned bases of 3000
        (h.tlen, h.ts, h.te) = (1000, 150, 950);                // target: 800 aligned bases of 1000, the shorter sequence
        assert!((shorter_cov(&h) - 0.8).abs() < 1e-12);         // 800 / 1000, not 900 / 3000
        (h.qlen, h.qs, h.qe, h.tlen, h.ts, h.te) = (1000, 150, 950, 3000, 100, 1000);   // roles swapped: the query is the shorter
        assert!((shorter_cov(&h) - 0.8).abs() < 1e-12);
        (h.qlen, h.qs, h.qe, h.tlen, h.ts, h.te) = (1000, 0, 900, 1000, 0, 800);        // equal lengths: the query's span (the chain's qlen <= tlen)
        assert!((shorter_cov(&h) - 0.9).abs() < 1e-12);
    }
    #[test]
    fn quantities_divide_by_the_right_lengths() {
        // distinct qlen / tlen / block, so that a swapped denominator changes the value (the brief's hit has qlen == block)
        let mut h = paf_hit("q", "t", 950);
        (h.qlen, h.qs, h.qe, h.tlen, h.ts, h.te, h.block) = (1000, 20, 980, 4000, 100, 1100, 1100);
        assert!((id_cov(&h) - 950.0 / 1100.0 * 960.0 / 1000.0).abs() < 1e-12);
        assert!((whole_length_d(&h) - 0.05).abs() < 1e-12);               // 1 - 950 / qlen: neither / block nor / tlen
    }
    #[test]
    fn quantities_stay_finite_on_empty_or_inverted_hits() {
        let mut h = paf_hit("q", "t", 0);
        (h.qlen, h.qs, h.qe, h.tlen, h.ts, h.te, h.block) = (0, 0, 0, 0, 0, 0, 0);
        assert_eq!((id_cov(&h), whole_length_d(&h), shorter_cov(&h)), (0.0, 1.0, 0.0));
        h.qlen = 10;                                                                  // an empty target: the target-side division is guarded too
        assert_eq!(shorter_cov(&h), 0.0);
        (h.qlen, h.qs, h.qe, h.tlen, h.ts, h.te, h.block, h.matches) = (100, 60, 40, 100, 60, 40, 10, 5);   // end before start: zero span, not a wrapped one
        assert_eq!((id_cov(&h), shorter_cov(&h)), (0.0, 0.0));
    }

    // ---- attribution of the unmapped and poorly placed reads to the nets and copies (prereg Amendments 13 / 13b) --------------------------

    /// A hit of a read of the attribution set (`qlen` bases) on a target of 5,000 bp: `span` read bases aligned (read coverage `span / qlen`,
    /// the fraction the rule reads), `matches` of them identical, divergence `de`. A target shorter than the read is built field by field.
    fn read_hit(read: &str, copy: &str, qlen: usize, span: usize, matches: usize, de: f64) -> PafHit {
        PafHit { q: read.into(), qlen, qs: 0, qe: span, strand: b'+', t: copy.into(), tlen: 5000, ts: 100, te: 100 + span, matches, block: span, de, cs: None }
    }
    /// The `--copies-fa` record names of two families' copies -> their family (F1 has two copies).
    fn copy_families() -> HashMap<String, String> {
        [("F1|0|chr1:100-5100|+|nexon=1", "F1"), ("F1|1|chr1:9000-14000|+|nexon=1", "F1"), ("F2|0|chr2:100-5100|+|nexon=1", "F2")]
            .into_iter()
            .map(|(t, f)| (t.to_string(), f.to_string()))
            .collect()
    }
    const F1_A: &str = "F1|0|chr1:100-5100|+|nexon=1";
    const F1_B: &str = "F1|1|chr1:9000-14000|+|nexon=1";
    const F2_A: &str = "F2|0|chr2:100-5100|+|nexon=1";
    fn attributed(hits: &[PafHit]) -> Vec<(String, String)> {
        let mut v: Vec<(String, String)> = attribute_by_hits(hits, &copy_families()).into_iter().collect();
        v.sort();
        v
    }
    fn pairs(v: &[(&str, &str)]) -> Vec<(String, String)> { v.iter().map(|(r, f)| (r.to_string(), f.to_string())).collect() }

    #[test]
    fn attribution_gives_a_read_to_the_family_of_a_best_hit_covering_most_of_it() {
        // (a) the best hit covers 90% of the read at de 0.05: the read joins that copy's family; each read is judged on its own hits
        let hits = vec![read_hit("r1", F1_B, 1000, 900, 880, 0.05), read_hit("r2", F2_A, 2000, 1900, 1890, 0.01)];
        assert_eq!(attributed(&hits), pairs(&[("r1", "F1"), ("r2", "F2")]));
        assert!(attributed(&[]).is_empty());
    }
    #[test]
    fn attribution_refuses_a_best_hit_covering_under_half_of_the_read() {
        // (b) a best hit covering 40% of the read attributes nothing
        assert!(attributed(&[read_hit("r", F1_A, 1000, 400, 398, 0.01)]).is_empty());
        // Review Focus 1: a read on copies of two families (a shared exon). Its best hit (most matches) covers 40%: no family, and no fall-back
        // to the other family's lesser hit although that one covers 60% (a 200-bp insertion inside it: fewer matches, one gap at de 0.01)
        let shared = vec![read_hit("r", F1_A, 1000, 400, 398, 0.01), read_hit("r", F2_A, 1000, 600, 390, 0.01)];
        assert!(attributed(&shared).is_empty());
        // the bound is inclusive: a best hit covering exactly half of the read attributes, one covering 0.4999 of it does not
        assert_eq!(attributed(&[read_hit("r", F2_A, 1000, 500, 495, 0.01)]), pairs(&[("r", "F2")]));
        assert!(attributed(&[read_hit("r", F2_A, 10_000, 4999, 4990, 0.01)]).is_empty());
    }
    #[test]
    fn attribution_measures_coverage_on_the_read_even_when_the_copy_is_shorter() {
        // Amendment 13: the best hit must cover >= 50% of the READ. A 3,000-bp read whose best hit spans 899 bp of a 900-bp copy covers 30% of
        // the read and is refused, although it covers the shorter sequence (the copy) almost whole
        let long = PafHit { qlen: 3000, qs: 1000, qe: 1899, tlen: 900, ts: 0, te: 899, ..read_hit("r", F1_A, 3000, 899, 895, 0.01) };
        assert!((shorter_cov(&long) - 899.0 / 900.0).abs() < 1e-12, "the shorter-sequence fraction, which the rule must not read");
        assert!(attributed(&[long]).is_empty());
        // a 1,200-bp read with a 900-bp hit on a 900-bp copy covers 75% of the read: attributed
        let fits = PafHit { qlen: 1200, qs: 150, qe: 1050, tlen: 900, ts: 0, te: 900, ..read_hit("r", F2_A, 1200, 900, 897, 0.01) };
        assert_eq!(attributed(&[fits]), pairs(&[("r", "F2")]));
    }
    #[test]
    fn attribution_refuses_a_divergent_best_hit() {
        // (c) Amendment 13b: de 0.2001 > 0.20 attributes nothing, however well the hit covers the read; the bound is inclusive (de 0.20 attributes)
        assert!(attributed(&[read_hit("r", F1_A, 1000, 1000, 799, 0.2001)]).is_empty());
        assert_eq!(attributed(&[read_hit("r", F1_A, 1000, 1000, 800, 0.20)]), pairs(&[("r", "F1")]));
        // a divergent best hit is not replaced by a closer lesser one either
        assert!(attributed(&[read_hit("r", F1_A, 1000, 1000, 780, 0.25), read_hit("r", F2_A, 1000, 800, 775, 0.005)]).is_empty());
    }
    #[test]
    fn attribution_prefers_a_net_read_of_one_family_over_a_lesser_copy_of_another() {
        // Amendment 13b: the targets are the run's net reads (`{family}|{read}`) and the copies; a read whose best hit is a net read of F1 joins
        // F1 although a lesser hit lies on a copy of F2, in either order of the hits
        let mut map = copy_families();
        map.insert("F1|m17".to_string(), "F1".to_string());
        let (net, copy) = (read_hit("r", "F1|m17", 1000, 950, 940, 0.12), read_hit("r", F2_A, 1000, 900, 880, 0.01));
        let got = |hits: &[PafHit]| attribute_by_hits(hits, &map).into_iter().collect::<Vec<_>>();
        assert_eq!(got(&[net.clone(), copy.clone()]), pairs(&[("r", "F1")]));
        assert_eq!(got(&[copy, net]), pairs(&[("r", "F1")]));
    }
    #[test]
    fn poorly_placed_reads_are_unnetted_with_a_divergent_or_mapq0_primary() {
        // Amendment 13b: a read in a net of the run (ruling R18) is never poorly placed, whatever its primary record
        assert!(!is_poorly_placed(0.5, 0, true));
        // de above 0.02 or MAPQ 0 makes an un-netted read poorly placed; de exactly 0.02 at a good MAPQ does not
        assert!(is_poorly_placed(0.021, 60, false));
        assert!(!is_poorly_placed(0.02, 60, false));
        assert!(is_poorly_placed(0.0, 0, false));
        assert!(!is_poorly_placed(0.002, 1, false));
        // MAPQ 255 (unavailable) is not MAPQ 0
        assert!(!is_poorly_placed(0.0, 255, false));
    }
    #[test]
    fn attribution_follows_the_hit_with_most_matches_not_the_lowest_de() {
        // (d) two hits that both pass: the one with more matches decides although the other has the lower de, in either order
        let a = read_hit("r", F1_A, 1000, 1000, 920, 0.08);
        let b = read_hit("r", F2_A, 1000, 900, 899, 0.001);
        assert_eq!(attributed(&[a.clone(), b.clone()]), pairs(&[("r", "F1")]));
        assert_eq!(attributed(&[b, a]), pairs(&[("r", "F1")]));
    }
    #[test]
    fn attribution_ignores_targets_outside_the_map_and_keeps_the_first_of_tied_hits() {
        // a target the map does not name (a partner row, a record the copies table lacks) is ignored: the best mapped hit decides
        let partner = read_hit("r", "F9|3|chr9:1-6000|+|nexon=1", 1000, 1000, 990, 0.01);
        assert_eq!(attributed(&[partner.clone(), read_hit("r", F2_A, 1000, 900, 880, 0.02)]), pairs(&[("r", "F2")]));
        assert!(attributed(&[partner]).is_empty());
        // equal matches: the first encountered (the `best_by_matches` rule)
        let (f1, f2) = (read_hit("r", F1_A, 1000, 900, 870, 0.03), read_hit("r", F2_A, 1000, 900, 870, 0.03));
        assert_eq!(attributed(&[f1.clone(), f2.clone()]), pairs(&[("r", "F1")]));
        assert_eq!(attributed(&[f2, f1]), pairs(&[("r", "F2")]));
    }

    // ---- read clustering, template-and-vote consensus, the real-vs-error test, one refinement pass --------------------------------------

    /// A hit between two 3000 bp reads covering `cov` of the shorter one (both are 3000 bp), for the clustering tests.
    fn ava_hit(a: &str, b: &str, matches: usize, cov: f64, de: f64) -> PafHit {
        let span = (3000.0 * cov).round() as usize;
        PafHit { q: a.into(), qlen: 3000, qs: 0, qe: span, strand: b'+', t: b.into(), tlen: 3000, ts: 0, te: span, matches, block: 3000, de, cs: None }
    }
    fn names_of(s: &str) -> Vec<String> { s.split_whitespace().map(str::to_string).collect() }

    /// An edit of a synthetic member relative to the template, in template coordinates (ascending, non-overlapping).
    enum Edit<'a> { Sub(usize, u8), Ins(usize, &'a [u8]), Del(usize, usize) }
    impl Edit<'_> { fn at(&self) -> usize { match self { Edit::Sub(p, _) | Edit::Ins(p, _) | Edit::Del(p, _) => *p } } }
    /// A synthetic member: its sequence (the template with `edits` applied) and its full-length `+` hit against the template, with the short
    /// `cs` minimap2 would write (bases lower case). `matches`, `block` and `de` are fillers: the vote reads `ts` and `cs` only.
    fn member_of(template: &[u8], edits: &[Edit]) -> (Vec<u8>, PafHit) {
        let lower = |s: &[u8]| String::from_utf8(s.to_ascii_lowercase()).unwrap();
        let (mut seq, mut cs, mut t) = (Vec::new(), String::new(), 0usize);
        for e in edits {
            assert!(e.at() >= t, "edits must ascend and not overlap");
            if e.at() > t { cs += &format!(":{}", e.at() - t); seq.extend_from_slice(&template[t..e.at()]); }
            t = e.at();
            match e {
                Edit::Sub(p, b) => { cs += &format!("*{}{}", lower(&template[*p..*p + 1]), lower(&[*b])); seq.push(*b); t = p + 1; }
                Edit::Ins(_, s) => { cs += &format!("+{}", lower(s)); seq.extend_from_slice(s); }
                Edit::Del(p, len) => { cs += &format!("-{}", lower(&template[*p..p + len])); t = p + len; }
            }
        }
        if t < template.len() { cs += &format!(":{}", template.len() - t); seq.extend_from_slice(&template[t..]); }
        let hit = PafHit { q: "m".into(), qlen: seq.len(), qs: 0, qe: seq.len(), strand: b'+', t: "t".into(), tlen: template.len(), ts: 0, te: template.len(), matches: template.len(), block: template.len(), de: 0.0, cs: Some(cs) };
        (seq, hit)
    }
    /// A synthetic member that aligns to the window `template[ts..te]` only: `edits` are in the window's coordinates, the hit starts at `ts`.
    fn member_span(template: &[u8], ts: usize, te: usize, edits: &[Edit]) -> (Vec<u8>, PafHit) {
        let (seq, mut h) = member_of(&template[ts..te], edits);
        (h.ts, h.te, h.tlen) = (ts, te, template.len());
        (seq, h)
    }
    /// `total` members, the first `carriers` of them with `edits` and the rest equal to the template.
    fn cluster_of(template: &[u8], carriers: usize, total: usize, edits: &[Edit]) -> Vec<(Vec<u8>, PafHit)> {
        let none: &[Edit] = &[];
        (0..total).map(|i| member_of(template, if i < carriers { edits } else { none })).collect()
    }
    fn consensus_of(template: &[u8], members: &[(Vec<u8>, PafHit)]) -> Vec<u8> {
        let refs: Vec<(&[u8], &PafHit)> = members.iter().map(|(s, h)| (s.as_slice(), h)).collect();
        consensus_from_template(template, &refs).unwrap()
    }
    /// The k-th base different from `b`.
    fn other_base(b: u8, k: usize) -> u8 { *b"ACGT".iter().filter(|&&c| c != b.to_ascii_uppercase()).nth(k).unwrap() }

    #[test]
    fn clustering_splits_two_percent_and_joins_point_two_percent() {
        // 6 reads: r0..r2 from copy A (0.2% errors), r3..r5 from copy B (2% from A); hits written as the minimap2 ava would give them
        let names: Vec<String> = (0..6).map(|i| format!("r{i}")).collect();
        let hit = |a: usize, b: usize, de: f64| PafHit { q: format!("r{a}"), qlen: 3000, qs: 0, qe: 3000, strand: b'+', t: format!("r{b}"), tlen: 3000, ts: 0, te: 3000, matches: 2900, block: 3000, de, cs: None };
        let ava = vec![hit(0,1,0.002), hit(1,2,0.003), hit(0,2,0.002), hit(3,4,0.002), hit(4,5,0.002), hit(0,3,0.021), hit(2,5,0.019)];
        let cl = cluster_reads(&names, &ava, 0.00958);
        assert_eq!(cl, vec![vec![0,1,2], vec![3,4,5]]);
    }
    #[test]
    fn clustering_decides_on_the_pairs_best_hit_and_needs_half_of_the_shorter_read() {
        let names = names_of("a b c d e f g h i j k l m n o p q r");           // indices 0..=17
        let d = 0.00958;
        let ava = vec![
            ava_hit("a", "b", 2900, 1.0, 0.002), ava_hit("b", "a", 2950, 1.0, 0.05),    // the hit with the MOST MATCHES decides, and it is too divergent: not joined
            ava_hit("c", "d", 2900, 1.0, 0.002), ava_hit("d", "c", 2900, 1.0, 0.05),    // equal matches: the first hit is kept, here a good one: joined
            ava_hit("e", "f", 2900, 1.0, 0.05), ava_hit("f", "e", 2900, 1.0, 0.002),    // ... and here the first is the divergent one: not joined
            ava_hit("g", "h", 2900, 1.0, d),                                            // de == delta joins
            ava_hit("i", "j", 2900, 0.5, 0.001),                                        // exactly half of the shorter read joins
            ava_hit("k", "l", 2900, 0.49, 0.001),                                       // just under half does not
            PafHit { strand: b'-', ..ava_hit("n", "m", 2900, 1.0, 0.002) },             // an opposite-strand pair, named the other way round: joined
            ava_hit("o", "o", 3000, 1.0, 0.0), ava_hit("o", "zz", 3000, 1.0, 0.0),      // a self hit and a hit to an unknown read are ignored
            ava_hit("p", "q", 2900, 1.0, 0.002), ava_hit("q", "r", 2900, 1.0, 0.002),    // p-q-r is one cluster with no p-r hit at all
        ];
        let want: Vec<Vec<usize>> = vec![vec![15, 16, 17], vec![2, 3], vec![6, 7], vec![8, 9], vec![12, 13], vec![0], vec![1], vec![4], vec![5], vec![10], vec![11], vec![14]];
        assert_eq!(cluster_reads(&names, &ava, d), want);                               // by size descending, then first index
    }
    #[test]
    fn clustering_of_nothing_and_of_unlinked_reads() {
        assert!(cluster_reads(&[], &[], 0.01).is_empty());
        let names = names_of("x y z");
        assert_eq!(cluster_reads(&names, &[], 0.01), vec![vec![0], vec![1], vec![2]]);   // singletons, in index order
        assert_eq!(cluster_reads(&names, &[ava_hit("z", "x", 2900, 1.0, 0.9)], 0.01), vec![vec![0], vec![1], vec![2]]);
    }
    #[test]
    fn clustering_matches_a_brute_force_reference_on_random_hit_lists() {
        // reproducible xorshift; 40 reads, 50 random hits (either direction, repeated pairs, self hits) so that components stay small and varied
        let mut x = 0x2545F4914F6CDD1Du64;
        let mut next = move || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let (n, d) = (40usize, 0.00958);
        let names: Vec<String> = (0..n).map(|i| format!("r{i}")).collect();
        for round in 0..40 {
            let ava: Vec<PafHit> = (0..50).map(|_| {
                let (a, b) = ((next() % n as u64) as usize, (next() % n as u64) as usize);
                let (m, c, e) = (2800 + (next() % 200) as usize, [0.3, 0.45, 0.5, 0.75, 1.0][(next() % 5) as usize], [0.001, 0.005, d, 0.012, 0.03][(next() % 5) as usize]);
                ava_hit(&names[a], &names[b], m, c, e)
            }).collect();
            // reference: per unordered pair the first hit with the most matches; joined iff cov >= 0.5 && de <= delta; components by label propagation
            let mut joined = vec![];
            for a in 0..n { for b in a + 1..n {
                let mut best: Option<&PafHit> = None;
                for h in &ava { if (h.q == names[a] && h.t == names[b]) || (h.q == names[b] && h.t == names[a]) { if best.map_or(true, |o| h.matches > o.matches) { best = Some(h); } } }
                if let Some(h) = best { if shorter_cov(h) >= 0.5 && h.de <= d { joined.push((a, b)); } }
            } }
            let mut label: Vec<usize> = (0..n).collect();
            loop {
                let mut changed = false;
                for &(a, b) in &joined { let m = label[a].min(label[b]); if label[a] != m || label[b] != m { label[a] = m; label[b] = m; changed = true; } }
                if !changed { break; }
            }
            let mut groups: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
            for i in 0..n { groups.entry(label[i]).or_default().push(i); }
            let mut want: Vec<Vec<usize>> = groups.into_values().collect();
            want.sort_by(|p, q| q.len().cmp(&p.len()).then(p[0].cmp(&q[0])));
            assert_eq!(cluster_reads(&names, &ava, d), want, "round {round}");
        }
    }

    // ---- the template: the medoid under a structural distance (prereg Amendments 13 / 13d) ------------------------------------------------

    /// A hand-written all-vs-all hit: read `q` (`qlen` bases, `qs..qe` aligned) on read `t` (`tlen` bases, `ts..te` aligned) with the
    /// `cs` given; `matches` decides between two hits of one pair, `block` and `de` are fillers.
    #[allow(clippy::too_many_arguments)]
    fn pair_hit(q: &str, qlen: usize, (qs, qe): (usize, usize), t: &str, tlen: usize, (ts, te): (usize, usize), matches: usize, cs: &str) -> PafHit {
        PafHit { q: q.into(), qlen, qs, qe, strand: b'+', t: t.into(), tlen, ts, te, matches, block: (qe - qs).max(te - ts), de: 0.002, cs: Some(cs.into()) }
    }
    /// `pair_hit` of two reads aligned end to end.
    fn whole_hit(q: &str, qlen: usize, t: &str, tlen: usize, matches: usize, cs: &str) -> PafHit { pair_hit(q, qlen, (0, qlen), t, tlen, (0, tlen), matches, cs) }
    /// `structural_scores` as `(member, aligned partners, sum of d, eligible)`, sorted by member.
    fn sorted_scores(members: &[usize], names: &[String], ava: &[PafHit]) -> Vec<(usize, usize, u64, bool)> {
        let mut v: Vec<(usize, usize, u64, bool)> = structural_scores(members, names, ava).unwrap().iter().map(|s| (s.member, s.aligned, s.sum_d, s.eligible)).collect();
        v.sort();
        v
    }
    fn template(members: &[usize], names: &[String], ava: &[PafHit], lens: &[usize]) -> TemplateChoice { structural_template(members, names, ava, lens).unwrap() }

    #[test]
    fn retained_intron_read_is_not_the_template() {
        // Review Focus 2. A = a clean 900-bp read; B = A with a 300-bp intron retained after base 400 (1,200 bp: the longest member, the old
        // template); C = A with 0.2% substitutions (2 in 900). Every hit covers both reads end to end, so d is the indel length: B's pairs carry
        // the 300 bases as a `-` (B the target) or a `+` (B the query), A-C none. Mean d: A 150, B 300, C 150; A and C tie on the mean and the
        // length, the name gives A
        let intron = lc(&rand_seq(300, 211));
        let names = names_of("A B C");
        let lens = [900, 1200, 900];
        let ava = vec![
            whole_hit("A", 900, "B", 1200, 900, &format!(":400-{intron}:500")),
            whole_hit("B", 1200, "C", 900, 898, &format!(":400+{intron}:50*ag:349*tc:99")),
            whole_hit("A", 900, "C", 900, 898, ":450*ag:349*tc:99"),
        ];
        assert_eq!(sorted_scores(&[0, 1, 2], &names, &ava), vec![(0, 2, 300, true), (1, 2, 600, true), (2, 2, 300, true)]);
        assert_eq!(template(&[0, 1, 2], &names, &ava, &lens), TemplateChoice::Medoid(0));
        // the cluster's order and the PAF's order do not matter
        let rev: Vec<PafHit> = ava.iter().rev().cloned().collect();
        assert_eq!(template(&[2, 1, 0], &names, &rev, &lens), TemplateChoice::Medoid(0));
        // the retained-intron read loses whatever its name: renamed to sort first, B is still not the template
        let renamed = names_of("x 0 y");
        let ava0: Vec<PafHit> = ava.iter().map(|h| PafHit { q: rename(&h.q), t: rename(&h.t), ..h.clone() }).collect();
        fn rename(s: &str) -> String { match s { "A" => "x", "B" => "0", _ => "y" }.to_string() }
        assert_eq!(template(&[0, 1, 2], &renamed, &ava0, &lens), TemplateChoice::Medoid(0));
    }
    #[test]
    fn skipping_read_is_not_the_template() {
        // r1 lacks a 200-bp exon of the others (after base 300): 700 bp. Its hits cover both reads end to end, so d is the indel length, and
        // carry the exon as a 200-bp `-` (r1 the query), a 200-bp `+` (r1 the target) and, as a splice preset writes it, a 200-bp `~` intron;
        // every other pair is clean. Mean d: r1 200, the others 200 / 3 each; the name gives r0
        let exon = lc(&rand_seq(200, 223));
        let names = names_of("r0 r1 r2 r3");
        let lens = [900, 700, 900, 900];
        let ava = vec![
            whole_hit("r1", 700, "r0", 900, 700, &format!(":300-{exon}:400")),
            whole_hit("r2", 900, "r1", 700, 700, &format!(":300+{exon}:400")),
            whole_hit("r1", 700, "r3", 900, 700, ":300~gt200ag:400"),
            whole_hit("r0", 900, "r2", 900, 900, ":900"),
            whole_hit("r0", 900, "r3", 900, 900, ":900"),
            whole_hit("r2", 900, "r3", 900, 900, ":900"),
        ];
        assert_eq!(sorted_scores(&[0, 1, 2, 3], &names, &ava), vec![(0, 3, 200, true), (1, 3, 600, true), (2, 3, 200, true), (3, 3, 200, true)]);
        assert_eq!(template(&[0, 1, 2, 3], &names, &ava, &lens), TemplateChoice::Medoid(0));
        // without r0 the skipping read is still not chosen: r2 and r3 tie at a mean of 100 (r1: 200), the name gives r2
        assert_eq!(template(&[1, 2, 3], &names, &ava, &lens), TemplateChoice::Medoid(2));
    }
    #[test]
    fn fragment_is_not_the_template() {
        // Amendment 13d. F, a 300-bp read, lies wholly inside the full-length reads A (925 bp: a 25-bp exon of its own), B and C (900 bp). F's
        // alignments carry no indel but leave >= 600 bases of each partner uncovered (both ends >= 20 bp), so d(F, p) >= 600 while d(p, F) = 0;
        // A-B and A-C carry A's 25-bp exon. Sums of d over 3 partners: A 50, B 25, C 25, F 625 + 600 + 600; B and C tie on the mean (8.3) and the
        // length, the name gives B. (Amendment 13's total of indel bases chose F: 0 against A 50, B 25, C 25.)
        let x25 = lc(&rand_seq(25, 257));
        let names = names_of("A B C F");
        let lens = [925, 900, 900, 300];
        let ava = vec![
            whole_hit("A", 925, "B", 900, 900, &format!(":400+{x25}:500")),
            whole_hit("C", 900, "A", 925, 900, &format!(":400-{x25}:500")),
            whole_hit("B", 900, "C", 900, 900, ":900"),
            pair_hit("F", 300, (0, 300), "A", 925, (500, 800), 300, ":300"),        // A uncovered: 500 before, 125 after
            pair_hit("B", 900, (300, 600), "F", 300, (0, 300), 300, ":300"),        // F the target: B (the query) uncovered 300 + 300
            pair_hit("F", 300, (0, 300), "C", 900, (300, 600), 300, ":300"),
        ];
        assert_eq!(sorted_scores(&[0, 1, 2, 3], &names, &ava), vec![(0, 3, 50, true), (1, 3, 25, true), (2, 3, 25, true), (3, 3, 1825, true)]);
        assert_eq!(template(&[0, 1, 2, 3], &names, &ava, &lens), TemplateChoice::Medoid(1));
        // with the fragment and one full-length read only, the full-length read is the template (d 0 against F's 625)
        assert_eq!(template(&[0, 3], &names, &ava, &lens), TemplateChoice::Medoid(0));
    }
    #[test]
    fn the_distance_adds_the_partners_uncovered_ends_of_20_bp_and_more() {
        // Amendment 13d: d(m, p) = the >= 20 bp indels of the pair's best hit + p's terminal bases that the hit leaves uncovered, each end only
        // from 20 bp, from whichever side m sits on. a (the query, 144 bp, 100 aligned) on b (161 bp, 122 aligned): b has 19 bases uncovered
        // before ts (not counted) and 20 after te (counted); a has 25 before qs (counted) and 19 after qe (not); the cs carries a 22-bp deletion
        let x22 = lc(&rand_seq(22, 271));
        let names = names_of("a b");
        let h = pair_hit("a", 144, (25, 125), "b", 161, (19, 141), 100, &format!(":50-{x22}:50"));
        assert_eq!(sorted_scores(&[0, 1], &names, &[h]), vec![(0, 1, 22 + 20, true), (1, 1, 22 + 25, true)]);
        // the same alignment written the other way round (b the query): the same distances
        let flipped = pair_hit("b", 161, (19, 141), "a", 144, (25, 125), 100, &format!(":50+{x22}:50"));
        assert_eq!(sorted_scores(&[0, 1], &names, &[flipped]), vec![(0, 1, 42, true), (1, 1, 47, true)]);
        // the uncovered end counts against the member that leaves it, not against its owner: r1 skips a 200-bp exon but carries a 250-bp 5'
        // extension of its own that no partner covers, so d(r1, p) = 200 and d(p, r1) = 450; in a cluster of three it is the medoid (mean 200
        // against 225), in a larger one the partners' means fall (450 / (n - 1)) and it is not
        let exon = lc(&rand_seq(200, 223));
        let three = names_of("r1 r2 r3");
        let ava = vec![
            pair_hit("r2", 900, (0, 900), "r1", 950, (250, 950), 700, &format!(":300+{exon}:400")),
            pair_hit("r1", 950, (250, 950), "r3", 900, (0, 900), 700, ":300~gt200ag:400"),
            whole_hit("r2", 900, "r3", 900, 900, ":900"),
        ];
        assert_eq!(sorted_scores(&[0, 1, 2], &three, &ava), vec![(0, 2, 400, true), (1, 2, 450, true), (2, 2, 450, true)]);
        assert_eq!(template(&[0, 1, 2], &three, &ava, &[950, 900, 900]), TemplateChoice::Medoid(0));
    }
    #[test]
    fn sparse_member_is_ineligible() {
        // Amendment 13d: a member is eligible when it is aligned to >= 50% of the cluster's other members. m6 is aligned to m0 only (1 of 6),
        // d 0 both ways; m0..m5 are aligned pairwise, each pair carrying a 30-bp exon swap (d 60). Mean d: m6 0 over 1 partner, m0 50, m1..m5
        // 60: m6 is skipped and m0 is the template
        let (xa, xb) = (lc(&rand_seq(30, 263)), lc(&rand_seq(30, 269)));
        let names = names_of("m0 m1 m2 m3 m4 m5 m6");
        let lens = [900; 7];
        let all: Vec<usize> = (0..7).collect();
        let swap = format!(":200+{xa}:200-{xb}:470");
        let mut ava: Vec<PafHit> = Vec::new();
        for i in 0..6 {
            for j in i + 1..6 { ava.push(whole_hit(&names[i], 900, &names[j], 900, 840, &swap)); }
        }
        ava.push(whole_hit("m6", 900, "m0", 900, 900, ":900"));
        let s = sorted_scores(&all, &names, &ava);
        assert_eq!((s[0], s[6]), ((0, 6, 300, true), (6, 1, 0, false)));
        assert!(s[1..6].iter().all(|&(_, aligned, sum_d, eligible)| (aligned, sum_d, eligible) == (5, 300, true)), "{s:?}");
        assert_eq!(template(&all, &names, &ava, &lens), TemplateChoice::Medoid(0));
        // aligned to 2 of 6 it is still under half (m1 now ties m0 at a mean of 50: the name keeps m0) ...
        ava.push(whole_hit("m1", 900, "m6", 900, 900, ":900"));
        assert_eq!(sorted_scores(&all, &names, &ava)[6], (6, 2, 0, false));
        assert_eq!(template(&all, &names, &ava, &lens), TemplateChoice::Medoid(0));
        // ... aligned to 3 of 6 (exactly half) it is eligible and, at a mean of 0, the template
        ava.push(whole_hit("m6", 900, "m2", 900, 900, ":900"));
        assert_eq!(sorted_scores(&all, &names, &ava)[6], (6, 3, 0, true));
        assert_eq!(template(&all, &names, &ava, &lens), TemplateChoice::Medoid(6));
    }
    #[test]
    fn small_cluster_eligibility_needs_one_aligned_partner() {
        // Amendment 13e (in place of 13d's "every member of a cluster of < 4 is eligible"): eligible = aligned to >= min(0.5 x (n - 1), 50)
        // others, which in clusters of 2 or 3 is one aligned partner. Of p, q, r only p-q is aligned (d 0): r, with no partner, is not eligible
        // (it would have no mean either) and is not the template although it is the longest; p and q tie, the name gives p
        let names = names_of("p q r s");
        let ava = vec![whole_hit("p", 900, "q", 900, 900, ":900")];
        let lens = [900, 900, 1000, 900];
        assert_eq!(sorted_scores(&[0, 1, 2], &names, &ava), vec![(0, 1, 0, true), (1, 1, 0, true), (2, 0, 0, false)]);
        assert_eq!(template(&[0, 1, 2], &names, &ava, &lens), TemplateChoice::Medoid(0));
        // p, q, s with p-q carrying a 25-bp exon swap (d 50) and s aligned to p only (d 0): s is eligible and the template
        let (xa, xb) = (lc(&rand_seq(25, 277)), lc(&rand_seq(25, 281)));
        let swap = whole_hit("p", 900, "q", 900, 850, &format!(":200+{xa}:200-{xb}:475"));
        let ava = vec![swap.clone(), whole_hit("s", 900, "p", 900, 900, ":900")];
        assert_eq!(sorted_scores(&[0, 1, 3], &names, &ava), vec![(0, 2, 50, true), (1, 1, 50, true), (3, 1, 0, true)]);
        assert_eq!(template(&[0, 1, 3], &names, &ava, &lens), TemplateChoice::Medoid(3));
        // with 4 members 2 partners are needed (min(1.5, 50) rounds up): with r (here 900 bp) aligned to p and q (d 50), s aligned to 1 of 3 is
        // not eligible and p (mean 100 / 3) is the template
        let swap_r = |q: &str| whole_hit(q, 900, "r", 900, 850, &format!(":200+{xa}:200-{xb}:475"));
        let ava = vec![swap, whole_hit("s", 900, "p", 900, 900, ":900"), swap_r("p"), swap_r("q")];
        assert_eq!(sorted_scores(&[0, 1, 2, 3], &names, &ava), vec![(0, 3, 100, true), (1, 2, 100, true), (2, 2, 100, true), (3, 1, 0, false)]);
        assert_eq!(template(&[0, 1, 2, 3], &names, &ava, &[900; 4]), TemplateChoice::Medoid(0));
    }
    #[test]
    fn large_cluster_eligibility_is_cap_aware() {
        // Amendment 13e: eligible = aligned to >= min(0.5 x (n - 1), 50) other members (the all-vs-all keeps at most 100 hits per query, so in a
        // cluster of hundreds of reads half of the others is out of reach). n members on a ring, each aligned to its k nearest on either side
        // (2k partners; clean end-to-end pairs, d 0 everywhere); m017 is the longest (the lengths decide the ties only)
        let ring = |n: usize, k: usize| -> (Vec<String>, Vec<PafHit>, Vec<usize>) {
            let names: Vec<String> = (0..n).map(|i| format!("m{i:03}")).collect();
            let ava = (0..n).flat_map(|i| (1..=k).map(move |j| (i, (i + j) % n))).map(|(a, b)| whole_hit(&names[a], 900, &names[b], 900, 900, ":900")).collect();
            let mut lens = vec![900; n];
            lens[17] = 950;
            (names, ava, lens)
        };
        let check = |n: usize, k: usize, eligible: bool| {
            let (names, ava, lens) = ring(n, k);
            let all: Vec<usize> = (0..n).collect();
            let s = structural_scores(&all, &names, &ava).unwrap();
            assert!(s.iter().all(|x| (x.aligned, x.eligible) == (2 * k, eligible)), "n = {n}, k = {k}");
            let want = if eligible { TemplateChoice::Medoid(17) } else { TemplateChoice::LongestAligned(17) };
            assert_eq!(template(&all, &names, &ava, &lens), want, "n = {n}, k = {k}");
        };
        // 300 members: 60 partners each -> all eligible (60 >= min(149.5, 50)), the medoid a tie broken by length; exactly 50 -> all; 40 -> none,
        // and the template is the longest member
        check(300, 30, true);
        check(300, 25, true);
        check(300, 20, false);
        // below the cap the half rules: 61 members need 30 partners (min(30, 50)); 100 members need 50 (49.5, rounded up by the integer count)
        check(61, 15, true);
        check(61, 14, false);
        check(100, 25, true);
        check(100, 24, false);
    }
    #[test]
    fn no_eligible_member_takes_the_longest() {
        // Amendments 13d / 13e: 5 members each aligned to at most 1 of the 4 others (pairs m0-m1 and m2-m3, m4 none): no member is eligible, so
        // the template is the longest member that has an aligned partner (`LongestAligned`: the binary logs it), equal lengths to the smaller
        // name. m4, the longest of all, has no partner and no mean: it is never chosen while another member has one
        let names = names_of("m0 m1 m2 m3 m4");
        let ava = vec![whole_hit("m0", 900, "m1", 900, 900, ":900"), whole_hit("m2", 900, "m3", 900, 900, ":900")];
        let all = [0, 1, 2, 3, 4];
        assert!(sorted_scores(&all, &names, &ava).iter().all(|&(_, _, _, eligible)| !eligible));
        assert_eq!(template(&all, &names, &ava, &[900, 900, 950, 900, 1000]), TemplateChoice::LongestAligned(2));
        assert_eq!(template(&all, &names, &ava, &[900, 980, 980, 900, 700]), TemplateChoice::LongestAligned(1));
        // a 4-read set as a refinement split-off or a merge leaves it: m0-m1 is its only pair, so nobody has the 2 partners needed; the two
        // partner-less members are the longest (2,000 and 1,500 bp) and are passed over for m1 (950 bp)
        assert_eq!(template(&[0, 1, 3, 4], &names, &ava[..1], &[900, 950, 970, 2000, 1500]), TemplateChoice::LongestAligned(1));
        // no member with an aligned partner at all (no mean anywhere): the longest member (`LongestUnaligned`, logged apart)
        assert_eq!(template(&[0, 1, 4], &names, &[], &[900, 980, 970, 900, 700]), TemplateChoice::LongestUnaligned(1));
        assert_eq!(template(&[2, 3, 4], &names, &ava[..1], &[900, 980, 970, 900, 700]), TemplateChoice::LongestUnaligned(2));
        // one member is its own template, not a fallback; no member is an error
        assert_eq!(template(&[3], &names, &[], &[900; 5]), TemplateChoice::Medoid(3));
        assert!(structural_template(&[], &names, &[], &[900; 5]).is_err());
    }
    #[test]
    fn tie_breaks_by_length_then_name() {
        // Review Focus 3: two members with identical structure (a clean pair: d 0 both ways) tie on the mean: the longer one, else the smaller
        // name, whatever the order of the cluster, of the PAF and of the net
        let names = names_of("b a");
        let ava = vec![whole_hit("b", 900, "a", 900, 900, ":900")];
        let pick = |members: &[usize], lens: &[usize]| template(members, &names, &ava, lens);
        assert_eq!(pick(&[0, 1], &[900, 900]), TemplateChoice::Medoid(1));                        // equal lengths: "a"
        assert_eq!(pick(&[1, 0], &[900, 900]), TemplateChoice::Medoid(1));
        assert_eq!(pick(&[0, 1], &[900, 901]), TemplateChoice::Medoid(1));                        // "a" is longer
        assert_eq!(pick(&[0, 1], &[901, 900]), TemplateChoice::Medoid(0));                        // "b" is longer: length before name
        // the mean comes before the length: q (920 bp) carries a 20-bp insertion against p and r, so p (900 bp) is the template ...
        let three = names_of("p q r");
        let ins20 = lc(&rand_seq(20, 229));
        let ins = |n: usize| vec![
            whole_hit("q", 900 + n, "p", 900, 900, &format!(":450+{}:450", &ins20[..n])),
            whole_hit("r", 900, "q", 900 + n, 900, &format!(":450-{}:450", &ins20[..n])),
            whole_hit("p", 900, "r", 900, 900, ":900"),
        ];
        assert_eq!(sorted_scores(&[0, 1, 2], &three, &ins(20)), vec![(0, 2, 20, true), (1, 2, 40, true), (2, 2, 20, true)]);
        assert_eq!(template(&[0, 1, 2], &three, &ins(20), &[900, 920, 900]), TemplateChoice::Medoid(0));
        // ... while a 19-bp insertion is no structure: every d is 0 and the longest member (q) is the template
        assert_eq!(sorted_scores(&[0, 1, 2], &three, &ins(19)), vec![(0, 2, 0, true), (1, 2, 0, true), (2, 2, 0, true)]);
        assert_eq!(template(&[0, 1, 2], &three, &ins(19), &[900, 919, 900]), TemplateChoice::Medoid(1));
    }
    #[test]
    fn the_structural_score_reads_each_member_pairs_best_hit_only() {
        let big = lc(&rand_seq(300, 233));
        let names = names_of("m0 m1 m2 out");
        // the pair's best hit by matches decides, in either order of the PAF: a secondary hit with fewer matches (a 300-bp `+`, 300 bases of m0
        // uncovered) does not count
        let clean = whole_hit("m0", 900, "m1", 900, 899, ":900");
        let secondary = pair_hit("m1", 900, (0, 900), "m0", 900, (300, 900), 600, &format!(":400+{big}:200"));
        for ava in [vec![clean.clone(), secondary.clone()], vec![secondary.clone(), clean.clone()]] {
            assert_eq!(sorted_scores(&[0, 1, 2], &names, &ava), vec![(0, 1, 0, true), (1, 1, 0, true), (2, 0, 0, false)]);
        }
        // equal matches: the first hit encountered (as `cluster_reads` picks a pair's hit); d(m0, m1) = 300, d(m1, m0) = 300 + m0's 300 uncovered
        let tied = PafHit { matches: 899, ..secondary.clone() };
        assert_eq!(sorted_scores(&[0, 1, 2], &names, &[clean.clone(), tied.clone()]), vec![(0, 1, 0, true), (1, 1, 0, true), (2, 0, 0, false)]);
        assert_eq!(sorted_scores(&[0, 1, 2], &names, &[tied, clean.clone()]), vec![(0, 1, 300, true), (1, 1, 600, true), (2, 0, 0, false)]);
        // hits to a read outside the cluster and self hits are ignored; a best hit without a cs aligns its pair with the uncovered ends as d
        let ignored = vec![
            whole_hit("m2", 900, "out", 1200, 890, &format!(":400+{big}:500")),
            whole_hit("m2", 900, "m2", 900, 900, &format!(":400+{big}:500")),
            PafHit { cs: None, ..pair_hit("m0", 900, (0, 900), "m2", 950, (0, 900), 900, "") },
        ];
        assert_eq!(sorted_scores(&[0, 1, 2], &names, &ignored), vec![(0, 1, 50, true), (1, 0, 0, false), (2, 1, 0, true)]);
        // the scores come back in the cluster's order, with the mean over the aligned partners
        let pair = vec![whole_hit("m2", 600, "m0", 900, 590, &format!(":400-{big}:200"))];
        let got = structural_scores(&[2, 1, 0], &names, &pair).unwrap();
        assert_eq!(got.iter().map(|s| (s.member, s.sum_d, s.mean_d())).collect::<Vec<_>>(), vec![(2, 300, Some(300.0)), (1, 0, None), (0, 300, Some(300.0))]);
        // a malformed cs of a best hit is an error naming the pair, never a silent 0
        let bad = vec![whole_hit("m0", 900, "m1", 900, 900, ":400+")];
        let msg = format!("{:#}", structural_scores(&[0, 1, 2], &names, &bad).unwrap_err());
        assert!(msg.contains("m0") && msg.contains("m1") && msg.contains("malformed cs"), "{msg}");
        assert!(structural_template(&[0, 1, 2], &names, &bad, &[900, 900, 900]).is_err());
    }
    #[test]
    fn refine_re_templates_the_kept_set_when_its_template_was_split_off() {
        // prereg Amendment 13: X (the cluster's template) is split off by the refinement; the kept set {A, B, C} is re-templated by the same
        // rule over ITS pairs only. Over the whole cluster X's pairs count and C is the medoid (mean d: A 800/3, B 300, C 700/3, X 300: A-X
        // carries a 500-bp indel and leaves 100 bases of A uncovered, C-X leaves 400 of X); over the kept pairs A and C tie at 150 (B, a
        // retained intron, 300) and the name gives A
        let (intron, other) = (lc(&rand_seq(300, 241)), lc(&rand_seq(500, 251)));
        let names = names_of("A B C X");
        let lens = [900, 1200, 900, 1300];
        let ava = vec![
            whole_hit("A", 900, "B", 1200, 900, &format!(":400-{intron}:500")),
            whole_hit("B", 1200, "C", 900, 900, &format!(":400+{intron}:500")),
            whole_hit("A", 900, "C", 900, 900, ":900"),
            pair_hit("A", 900, (50, 850), "X", 1300, (0, 1300), 800, &format!(":400-{other}:400")),
            pair_hit("C", 900, (0, 900), "X", 1300, (0, 900), 900, ":900"),
        ];
        assert_eq!(template(&[0, 1, 2, 3], &names, &ava, &lens), TemplateChoice::Medoid(2), "over the whole cluster, X's pairs count");
        assert_eq!(refined_template(3, &[0, 1, 2], &names, &ava, &lens).unwrap(), TemplateChoice::Medoid(0));
        // a template that the refinement keeps stays the template, although the kept set's medoid would be another member
        assert_eq!(refined_template(3, &[0, 1, 2, 3], &names, &ava, &lens).unwrap(), TemplateChoice::Kept(3));
        assert_eq!(refined_template(1, &[0, 1, 2], &names, &ava, &lens).unwrap(), TemplateChoice::Kept(1));
    }

    #[test]
    fn consensus_vote_fixes_errors_and_keeps_majority_indels() {
        let template = b"ACGTACGTACGTTTTTACGTACGT".to_vec();          // template carries a 1-base error at index 4 (A instead of G)
        let truth    = b"ACGTGCGTACGTTTTTACGTACGT".to_vec();
        // three members equal to truth: cs vs template = ":4*ag:19"  (4 eq, sub a->g, 19 eq)
        let h = PafHit { q: "m".into(), qlen: 24, qs: 0, qe: 24, strand: b'+', t: "t".into(), tlen: 24, ts: 0, te: 24, matches: 23, block: 24, de: 0.04, cs: Some(":4*ag:19".into()) };
        let members: Vec<(&[u8], &PafHit)> = vec![(&truth[..], &h), (&truth[..], &h), (&truth[..], &h)];
        assert_eq!(consensus_from_template(&template, &members).unwrap(), truth);
    }
    #[test]
    fn consensus_keeps_an_exon_that_two_of_five_members_skip_and_inserts_a_24bp_insertion_carried_by_three() {
        // ruling R2: the template carries an exon (40..64, 24 bp) that 2 of 5 members skip (a 24 bp `-` in their cs); 3 of 5 carry a 24 bp insertion at 100
        let template = rand_seq(120, 41);
        let ins: &[u8] = b"GATTACAGATTACACCGGTTAACC";
        assert_eq!(ins.len(), 24);
        let carrier = || member_of(&template, &[Edit::Ins(100, ins)]);
        let skipper = || member_of(&template, &[Edit::Del(40, 24)]);
        let members = vec![carrier(), carrier(), carrier(), skipper(), skipper()];
        let mut want = template.clone(); want.splice(100..100, ins.iter().copied());
        let got = consensus_of(&template, &members);
        assert_eq!(got.len(), 120 + 24);                                  // the exon stayed and the insertion went in
        assert_eq!(got, want);
    }
    #[test]
    fn consensus_treats_indels_of_20_bp_and_more_as_structure_under_the_amendment_15_majority() {
        let template = rand_seq(120, 43);
        let ins24: &[u8] = b"GATTACAGATTACACCGGTTAACC";
        // deletions: 19 bp follows the 50% majority (4 of 5 -> applied) but 20 bp never does, even when every member carries it
        let mut minus19 = template.clone(); minus19.drain(40..59);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 4, 5, &[Edit::Del(40, 19)])), minus19);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 4, 5, &[Edit::Del(40, 20)])), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 5, 5, &[Edit::Del(40, 24)])), template);
        // insertions (prereg Amendment 15, the consensus-defect fix): a >= 20 bp insertion needs the floor of 3 carriers AND a
        // majority of the members covering the column. 3 of 8 (37.5%) — inserted under the pre-15 rule, the defect's core case — no longer
        // enters; 2 carriers never did (the floor), whatever their share ...
        let with = |ins: &[u8]| { let mut v = template.clone(); v.splice(100..100, ins.iter().copied()); v };
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 8, &[Edit::Ins(100, &ins24[..20])])), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 8, &[Edit::Ins(100, &ins24[..20])])), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 2, &[Edit::Ins(100, &ins24[..20])])), template);   // 100%, but only 2 members
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 3, &[Edit::Ins(100, &ins24[..20])])), template);   // 67% of 3: 20 bp needs 3 carriers
        // ... a majority does: 4 of 8 (exactly 50%), 5 of 8, and all of 3 ...
        assert_eq!(consensus_of(&template, &cluster_of(&template, 4, 8, &[Edit::Ins(100, &ins24[..20])])), with(&ins24[..20]));
        assert_eq!(consensus_of(&template, &cluster_of(&template, 5, 8, &[Edit::Ins(100, &ins24[..20])])), with(&ins24[..20]));
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 3, &[Edit::Ins(100, &ins24[..20])])), with(&ins24[..20]));
        // ... while 19 bp needs 50% of the covering members: 3 of 8 no, 4 of 8 yes
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 8, &[Edit::Ins(100, &ins24[..19])])), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 4, 8, &[Edit::Ins(100, &ins24[..19])])), with(&ins24[..19]));
    }
    #[test]
    fn consensus_votes_the_insertions_of_20_bp_and_more_before_the_short_majority_at_one_position() {
        let template = rand_seq(120, 47);
        let ins24: &[u8] = b"GATTACAGATTACACCGGTTAACC";
        let with = |ins: &[u8]| { let mut v = template.clone(); v.splice(100..100, ins.iter().copied()); v };
        // prereg Amendments 13 + 15: at one position the insertions >= 20 bp are considered first — but only when their carriers are a
        // majority of the covering members. 8 members cover column 100: 3 carry a 24-bp insertion there, 4 carry a 1-bp insertion there
        // (50% of the covering members) and 1 carries none. The 24-bp one is NOT a majority (2x3 = 6 < 8), so the 1-bp one is inserted: one
        // winner per position (a member's cs holds one insertion before a column, so no member supports both)
        let mut members = cluster_of(&template, 3, 3, &[Edit::Ins(100, ins24)]);
        members.extend(cluster_of(&template, 4, 5, &[Edit::Ins(100, b"T")]));
        assert_eq!(members.len(), 8);
        assert_eq!(consensus_of(&template, &members), with(b"T"));
        // when the long insertion IS a majority it still wins over the short one: 5 of 9 carry the 24-bp one, 4 of 9 the 1-bp one
        let mut members = cluster_of(&template, 5, 5, &[Edit::Ins(100, ins24)]);
        members.extend(cluster_of(&template, 4, 4, &[Edit::Ins(100, b"T")]));
        assert_eq!(members.len(), 9);
        assert_eq!(consensus_of(&template, &members), with(ins24));
        // a long insertion under 3 carriers leaves the position to the short rule: 2 carriers of the 24-bp one, 4 of 8 with the 1-bp one
        let mut members = cluster_of(&template, 2, 2, &[Edit::Ins(100, ins24)]);
        members.extend(cluster_of(&template, 4, 6, &[Edit::Ins(100, b"T")]));
        assert_eq!(consensus_of(&template, &members), with(b"T"));
        // ... and the short rule still needs its 50%: 2 carriers of the 24-bp one, 3 of 8 with the 1-bp one: nothing is inserted
        let mut members = cluster_of(&template, 2, 2, &[Edit::Ins(100, ins24)]);
        members.extend(cluster_of(&template, 3, 6, &[Edit::Ins(100, b"T")]));
        assert_eq!(consensus_of(&template, &members), template);
        // two long insertions: the more frequent one wins (3 against 4 carriers), whichever sorts first
        let other: &[u8] = b"CCCCGGGGAAAATTTTCCCCGGGGAA";
        let mut members = cluster_of(&template, 3, 3, &[Edit::Ins(100, ins24)]);
        members.extend(cluster_of(&template, 4, 4, &[Edit::Ins(100, other)]));
        assert_eq!(consensus_of(&template, &members), with(other));
    }
    #[test]
    fn consensus_small_indels_follow_the_fifty_percent_majority_of_three_or_more_covering_members() {
        let template = rand_seq(80, 81);
        let (mut minus2, mut plus3) = (template.clone(), template.clone());
        minus2.drain(30..32); plus3.splice(50..50, b"CCA".iter().copied());
        let (del2, ins3) = ([Edit::Del(30, 2)], [Edit::Ins(50, b"CCA")]);
        // 3 of 5 and exactly 2 of 4 (>= 50%) apply; 2 of 5 does not; two covering members are below the floor of 3
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 5, &del2)), minus2);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 4, &del2)), minus2);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 5, &del2)), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 2, &del2)), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 5, &ins3)), plus3);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 4, &ins3)), plus3);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 5, &ins3)), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 2, &ins3)), template);
        // an insertion before the first column and one after the last are placed there
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 4, &[Edit::Ins(0, b"GG")])), [&b"GG"[..], &template].concat());
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 4, &[Edit::Ins(80, b"GG")])), [&template[..], &b"GG"[..]].concat());
        // the most frequent insertion sequence wins over a rarer one, whichever sorts first (3 TT + 1 AC) ...
        let mut skew = cluster_of(&template, 3, 3, &[Edit::Ins(50, b"TT")]);
        skew.extend(cluster_of(&template, 1, 1, &[Edit::Ins(50, b"AC")]));
        let mut plus_tt = template.clone(); plus_tt.splice(50..50, b"TT".iter().copied());
        assert_eq!(consensus_of(&template, &skew), plus_tt);
        // ... and equally frequent ones tie to the smaller sequence (3 + 3 of 6, each at 50%), whatever the hash order: every call builds its own
        // HashMaps with their own random seed, so repeating the call exposes a missing tie-break
        let mut mixed = cluster_of(&template, 3, 3, &[Edit::Ins(50, b"TT")]);
        mixed.extend(cluster_of(&template, 3, 3, &[Edit::Ins(50, b"AC")]));
        let mut plus_ac = template.clone(); plus_ac.splice(50..50, b"AC".iter().copied());
        for _ in 0..40 { assert_eq!(consensus_of(&template, &mixed), plus_ac); }
    }
    #[test]
    fn consensus_substitution_needs_three_covering_members_and_more_votes_than_the_template_base() {
        let template = rand_seq(50, 71);
        let alt = |k: usize| other_base(template[20], k);
        let variant = |k: usize| member_of(&template, &[Edit::Sub(20, alt(k))]);
        let plain = || member_of(&template, &[]);
        let mut with_alt = template.clone(); with_alt[20] = alt(0);
        assert_eq!(consensus_of(&template, &[variant(0), variant(0)]), template);                              // 2 covering members: below the floor
        assert_eq!(consensus_of(&template, &[variant(0), variant(0), variant(0)]), with_alt);                  // 3 of 3
        assert_eq!(consensus_of(&template, &[variant(0), variant(0), plain()]), with_alt);                     // 2 votes > 1 template base
        assert_eq!(consensus_of(&template, &[variant(0), variant(0), plain(), plain()]), template);            // a tie keeps the template
        assert_eq!(consensus_of(&template, &[variant(0), variant(0), variant(0), plain(), plain()]), with_alt);
        // two variants with the same support: the smaller base, provided it beats the template base (whatever the hash order: repeated calls)
        let mut with_min = template.clone(); with_min[20] = alt(0).min(alt(1));
        for _ in 0..40 {
            assert_eq!(consensus_of(&template, &[variant(0), variant(0), variant(1), variant(1), plain()]), with_min);
            assert_eq!(consensus_of(&template, &[variant(1), variant(1), variant(0), variant(0), plain()]), with_min);
        }
        // the better supported variant wins over a rarer one, also when it is the larger base (3 x alt(1) against 1 x alt(0))
        let mut with_big = template.clone(); with_big[20] = alt(1);
        assert!(alt(1) > alt(0));
        assert_eq!(consensus_of(&template, &[variant(1), variant(1), variant(1), variant(0)]), with_big);
        assert_eq!(consensus_of(&template, &[variant(1), variant(1), variant(1), variant(0), plain()]), with_big);
        // members that delete the column carry no base: they are neither template carriers nor votes (2 variants, 2 deleters, 1 template base)
        let del1 = || member_of(&template, &[Edit::Del(20, 1)]);
        assert_eq!(consensus_of(&template, &[variant(0), variant(0), del1(), del1(), plain()]), with_alt);
    }
    #[test]
    fn consensus_recovers_the_truth_from_ten_reads_with_private_errors_and_a_homopolymer_indel() {
        let mut truth = rand_seq(200, 17);
        truth[99] = b'A'; truth[100..106].fill(b'T'); truth[106] = b'C';                    // a 6-T homopolymer at 100..106
        let mut template = truth.clone(); template.insert(100, b'T');                       // the template has 7 T's (an insertion error) ...
        template[151] = if truth[150] == b'G' { b'A' } else { b'G' };                       // ... and a substitution error (truth 150 = template 151)
        let fix = truth[150];
        let mut members = vec![];
        for p in [10usize, 30, 60, 90, 120, 140, 170] {                                      // 7 members: their own error, the homopolymer fixed, the template error fixed
            let mut e = vec![Edit::Sub(p, other_base(template[p], 0)), Edit::Del(100, 1), Edit::Sub(151, fix)];
            e.sort_by_key(|x| x.at()); members.push(member_of(&template, &e));
        }
        members.push(member_of(&template, &[Edit::Del(100, 2), Edit::Sub(151, fix)]));       // one that also lost a second T (5 T's)
        members.push(member_of(&template, &[Edit::Sub(151, fix)]));                          // one that shares the template's 7 T's
        members.push(member_of(&template, &[Edit::Del(100, 1), Edit::Sub(151, fix), Edit::Sub(190, other_base(template[190], 1))]));
        assert_eq!(members.len(), 10);
        assert_eq!(consensus_of(&template, &members), truth);
    }
    #[test]
    fn consensus_intron_and_long_deletion_members_do_not_cover_but_a_short_deletion_member_does() {
        // ruling R5. 8 members; 4 lack columns 40..64; of the 4 that have them, 2 carry a 5 bp insertion inside (at 50)
        let template = rand_seq(100, 61);
        let mut plus5 = template.clone(); plus5.splice(50..50, b"ACGTA".iter().copied());
        let (carrier, plain) = (|| member_of(&template, &[Edit::Ins(50, b"ACGTA")]), || member_of(&template, &[]));
        let eight = |absent: &dyn Fn() -> (Vec<u8>, PafHit)| vec![carrier(), carrier(), plain(), plain(), absent(), absent(), absent(), absent()];
        // an intron, a 24 bp `-` and a 20 bp `-` all leave the 4 members out of the coverage: 2 of the 4 covering members carry the insertion (50%)
        let spliced = || { let (s, mut h) = member_of(&template, &[Edit::Del(40, 24)]); h.cs = Some(":40~gt24ag:36".into()); (s, h) };
        assert_eq!(consensus_of(&template, &eight(&spliced)), plus5);
        assert_eq!(consensus_of(&template, &eight(&|| member_of(&template, &[Edit::Del(40, 24)]))), plus5);
        assert_eq!(consensus_of(&template, &eight(&|| member_of(&template, &[Edit::Del(40, 20)]))), plus5);
        // a 19 bp deletion is error-sized: those 4 members cover the columns and delete them (4 of 8 = 50% removes 40..59), and the insertion
        // has only 2 of 8 (25%)
        let mut minus19 = template.clone(); minus19.drain(40..59);
        assert_eq!(consensus_of(&template, &eight(&|| member_of(&template, &[Edit::Del(40, 19)]))), minus19);
    }
    #[test]
    fn consensus_members_that_skip_an_exon_do_not_vote_inside_it_and_the_carriers_decide() {
        // ruling R5: a member with a >= 20 bp deletion has no base there (an isoform's absence): it neither covers the columns nor counts as
        // agreeing with the template, so the vote inside the exon is among the members that carry it and needs 3 of them
        let template = rand_seq(120, 53);                                  // exon 40..64, variant column 50
        let alt = other_base(template[50], 0);
        let mut with_alt = template.clone(); with_alt[50] = alt;
        let (skip, variant, plain) = (|| member_of(&template, &[Edit::Del(40, 24)]), || member_of(&template, &[Edit::Sub(50, alt)]), || member_of(&template, &[]));
        assert_eq!(consensus_of(&template, &[variant(), variant(), skip(), skip(), skip()]), template);                        // 3 of 5 skip: 2 carriers are below the floor of 3
        assert_eq!(consensus_of(&template, &[variant(), variant(), plain(), skip(), skip(), skip()]), with_alt);               // 3 of 6 skip: 2 against 1 among the carriers flips it
        assert_eq!(consensus_of(&template, &[variant(), variant(), variant(), skip(), skip(), skip()]), with_alt);             // 3 against 0
        assert_eq!(consensus_of(&template, &[variant(), variant(), plain(), plain(), skip(), skip(), skip()]), template);      // 2 against 2: the template stays
    }
    #[test]
    fn consensus_small_deletions_are_a_plain_majority_of_the_covering_members_even_inside_a_skipped_exon() {
        // ruling R5: no guard for columns that a member skips with a >= 20 bp deletion, and the skippers do not dilute the share
        let template = rand_seq(120, 53);
        let mut minus1 = template.clone(); minus1.remove(50);
        let (skip, del1, plain) = (|| member_of(&template, &[Edit::Del(40, 24)]), || member_of(&template, &[Edit::Del(50, 1)]), || member_of(&template, &[]));
        assert_eq!(consensus_of(&template, &[del1(), del1(), del1(), skip(), skip()]), minus1);          // 3 of the 3 carriers delete the column
        assert_eq!(consensus_of(&template, &[del1(), del1(), plain(), skip(), skip()]), minus1);         // 2 of 3
        assert_eq!(consensus_of(&template, &[del1(), plain(), plain(), skip(), skip()]), template);      // 1 of 3
        assert_eq!(consensus_of(&template, &[del1(), del1(), skip(), skip(), skip()]), template);        // 2 carriers: below the floor of 3
    }
    #[test]
    fn consensus_trims_template_ends_covered_by_fewer_than_two_members() {
        // spec §5.4 as ruled in R5: leading and trailing template columns with coverage < 2 are dropped (the template counts only as a member: here a self hit)
        let template = rand_seq(150, 33);
        let selfhit = || member_of(&template, &[]);
        let window = |ts: usize, te: usize| member_span(&template, ts, te, &[]);
        // the template is 30 bp longer at the 5' end than every other member (they start at 30 and at 45) and 20 bp longer at the 3' end:
        // the consensus starts where the second member starts and stops where the others stop
        assert_eq!(consensus_of(&template, &[selfhit(), window(30, 130), window(45, 130)]), template[30..130].to_vec());
        // coverage of exactly 2 survives, also at the ends ...
        assert_eq!(consensus_of(&template, &[selfhit(), selfhit()]), template);
        assert_eq!(consensus_of(&template, &[selfhit(), window(10, 140)]), template[10..140].to_vec());
        // ... one member covers nothing twice and no member covers nothing: an empty consensus
        assert!(consensus_of(&template, &[selfhit()]).is_empty());
        assert!(consensus_of(&template, &[]).is_empty());
        // only the ENDS are trimmed: a stretch in the middle that a single member covers stays
        assert_eq!(consensus_of(&template, &[selfhit(), window(0, 60), window(90, 150)]), template);
        // an insertion goes with the column it precedes: one before the first surviving column stays, one after the last surviving column (it
        // precedes a trimmed one) goes. (Alignments that begin or end with an insertion: minimap2 does not write them, a vote over `cs` must handle them.)
        let gg = |ts: usize, te: usize, at: usize| member_span(&template, ts, te, &[Edit::Ins(at, b"GG")]);
        let mut want = b"GG".to_vec(); want.extend_from_slice(&template[30..130]);
        assert_eq!(consensus_of(&template, &[selfhit(), gg(30, 130, 0), gg(30, 130, 0), gg(30, 130, 0)]), want);
        assert_eq!(consensus_of(&template, &[gg(20, 120, 100), gg(20, 120, 100), gg(20, 120, 100)]), template[20..120].to_vec());
    }
    #[test]
    fn consensus_places_a_member_that_starts_inside_the_template_at_its_ts() {
        let template = rand_seq(60, 5);
        let alt = other_base(template[30], 0);
        let mut with_alt = template.clone(); with_alt[30] = alt;
        let partial = || { let (s, mut h) = member_of(&template[25..], &[Edit::Sub(5, alt)]); h.ts = 25; h.te = 60; h.tlen = 60; (s, h) };   // cs `:5*xy:29` from ts = 25
        // two full-length members keep every column covered twice; the three carriers of the variant at 30 are one of them and the two partial members
        assert_eq!(consensus_of(&template, &[member_of(&template, &[Edit::Sub(30, alt)]), member_of(&template, &[]), partial(), partial()]), with_alt);
        // the two partial members alone cover 25..60 twice (below the vote's 3): the consensus is that stretch, unvoted, starting at their ts
        assert_eq!(consensus_of(&template, &[partial(), partial()]), template[25..].to_vec());
    }
    #[test]
    fn consensus_ignores_members_without_a_cs_and_rejects_a_malformed_one() {
        let template = rand_seq(40, 91);
        let alt = other_base(template[10], 0);
        let mut with_alt = template.clone(); with_alt[10] = alt;
        let (s, h) = member_of(&template, &[Edit::Sub(10, alt)]);
        let none = PafHit { cs: None, ..h.clone() };
        let read = &s[..];
        let mut three: Vec<(&[u8], &PafHit)> = vec![(read, &h), (read, &h), (read, &h)];
        three.extend([(read, &none), (read, &none), (read, &none), (read, &none), (read, &none)]);                  // five members with no alignment string abstain
        assert_eq!(consensus_from_template(&template, &three).unwrap(), with_alt);
        let two: Vec<(&[u8], &PafHit)> = vec![(read, &h), (read, &h), (read, &none), (read, &none), (read, &none)];
        assert_eq!(consensus_from_template(&template, &two).unwrap(), template);                                    // only 2 members count
        let bad = PafHit { q: "read7".into(), cs: Some(":10*a".into()), ..h.clone() };
        let err = consensus_from_template(&template, &[(read, &h), (read, &bad)]).unwrap_err().to_string();
        assert!(err.contains("read7") && err.contains("malformed cs"), "the error names the member and the cs problem: {err}");
    }
    #[test]
    fn consensus_survives_an_empty_template_and_a_cs_longer_than_the_template_and_returns_upper_case() {
        let hit = |cs: &str, ts: usize| PafHit { q: "m".into(), qlen: 0, qs: 0, qe: 0, strand: b'+', t: "t".into(), tlen: 0, ts, te: 0, matches: 0, block: 0, de: 0.0, cs: Some(cs.into()) };
        let read: &[u8] = b"";
        // an empty template: an insertion-only cs from three members must not panic (a bare cover[n - 1] would underflow)
        let ins = hit("+acgtacgtacgtacgtacgtacgt", 0);
        assert_eq!(consensus_from_template(b"", &[(read, &ins), (read, &ins), (read, &ins)]).unwrap(), b"".to_vec());
        assert_eq!(consensus_from_template(b"", &[]).unwrap(), b"".to_vec());
        // a cs that runs past the end of the template (garbage in) is clamped: nothing beyond the end is voted or appended
        let template = rand_seq(24, 5);
        let over = hit(":50*ac:10-acg+tt", 0);
        assert_eq!(consensus_from_template(&template, &[(read, &over), (read, &over), (read, &over)]).unwrap(), template);
        let after = hit(":5+aaaa", 30);                                                                      // members that start past the end cover nothing ...
        assert!(consensus_from_template(&template, &[(read, &after), (read, &after), (read, &after)]).unwrap().is_empty());   // ... and nothing survives the end trim
        let whole = hit(":9", 0);
        assert_eq!(consensus_from_template(b"acgtnacgt", &[(read, &whole), (read, &whole)]).unwrap(), b"ACGTNACGT".to_vec());   // lower-case template: upper-case consensus
    }

    #[test]
    fn consensus_a15b_shortens_an_insertion_that_begins_with_the_skipped_template_bases() {
        // prereg Amendment 15b (a), regression for the `+X·E ~|X|` cs pattern: minimap2 splice:hq writes an exon E the template lacks
        // with the skipped template bases X prepended, at a column that depends on where the read starts; normalised (`:|X| +E`), the
        // event lands at ONE column as E alone. 3 carriers of 5 covering: a majority, so E enters — but as 25 bp at column 120, not the
        // raw 45 bp piece at column 100 the pre-fix vote would have inserted
        let template = rand_seq(200, 97);
        let e: &[u8] = b"ACGTTGCAACGTTGCAACGTTGCAA";
        assert_eq!(e.len(), 25);
        let x = &template[100..120];
        let mut piece = x.to_vec();
        piece.extend_from_slice(e);
        let carrier = |q: &str| {
            (
                Vec::new(),
                PafHit { q: q.into(), qlen: 160 + piece.len(), qs: 0, qe: 160 + piece.len(), strand: b'+', t: "t".into(), tlen: 200,
                    ts: 0, te: 200, matches: 160 + e.len(), block: 160 + piece.len(), de: 0.0,
                    cs: Some(format!(":100+{}~gt20ag:60", lc(&piece))) },
            )
        };
        let plain = || member_of(&template, &[]);
        let members = vec![carrier("c1"), carrier("c2"), carrier("c3"), plain(), plain()];
        let mut want = template.clone();
        want.splice(120..120, e.iter().copied());
        assert_eq!(consensus_of(&template, &members), want);
        // chained skips: `+X'·X·E ~10 ~20` — the insertion keeps feeding the next skip and still lands as E at column 120
        let x2 = &template[90..100];
        let mut piece2 = x2.to_vec();
        piece2.extend_from_slice(&piece);
        let carrier2 = |q: &str| {
            (
                Vec::new(),
                PafHit { q: q.into(), qlen: 150 + piece2.len(), qs: 0, qe: 150 + piece2.len(), strand: b'+', t: "t".into(), tlen: 200,
                    ts: 0, te: 200, matches: 150 + e.len(), block: 150 + piece2.len(), de: 0.0,
                    cs: Some(format!(":90+{}~gt10ag~ct20ac:60", lc(&piece2))) },
            )
        };
        let members2 = vec![carrier2("d1"), carrier2("d2"), carrier2("d3"), plain(), plain()];
        assert_eq!(consensus_of(&template, &members2), want);
    }
    #[test]
    fn consensus_a15b_leaves_an_insertion_that_does_not_begin_with_the_skipped_bases() {
        // the other cs pattern: `+E ~|X|` where E does NOT begin with the skipped template bases — no normalisation, the insertion
        // stays where the aligner put it and answers only to the Amendment 15 majority test
        let template = rand_seq(200, 98);
        let e: &[u8] = b"TTGCAAGCTTGCAAGCTTGCAAGCT";
        assert_eq!(e.len(), 25);
        let carrier = |q: &str| {
            (
                Vec::new(),
                PafHit { q: q.into(), qlen: 160 + e.len(), qs: 0, qe: 160 + e.len(), strand: b'+', t: "t".into(), tlen: 200,
                    ts: 0, te: 200, matches: 160 + e.len(), block: 160 + e.len(), de: 0.0,
                    cs: Some(format!(":100+{}~gt20ag:60", lc(e))) },
            )
        };
        let plain = || member_of(&template, &[]);
        let members = vec![carrier("c1"), carrier("c2"), carrier("c3"), carrier("c4"), carrier("c5"), plain(), plain(), plain()];
        let mut want = template.clone();
        want.splice(100..100, e.iter().copied());
        assert_eq!(consensus_of(&template, &members), want); // 5 of 8 covering: a majority, inserted at column 100 as-is
    }
    #[test]
    fn consensus_a15b_inserts_one_identical_long_insertion_once_at_its_best_supported_column() {
        // prereg Amendment 15b (b): one copy per identical long insertion per consensus — the same 24-bp piece carried at two columns
        // (as when carriers' alignments place it differently) enters ONCE, at the column with the most carriers (ties to the earliest).
        // 5 carriers at each of columns 100 and 150: both pass the majority test on their own (2x5 = 10 >= 10 covering), so the
        // pre-(b) vote inserted it twice; the one-copy rule keeps column 100 alone
        let template = rand_seq(200, 99);
        let ins24: &[u8] = b"GATTACAGATTACACCGGTTAACC";
        let carrier = |q: &str, at: usize| {
            (
                Vec::new(),
                PafHit { q: q.into(), qlen: 224, qs: 0, qe: 224, strand: b'+', t: "t".into(), tlen: 200, ts: 0, te: 200, matches: 200,
                    block: 224, de: 0.0, cs: Some(format!(":{at}+{}:{}", lc(ins24), 200 - at)) },
            )
        };
        let mut members = Vec::new();
        for i in 0..5 { members.push(carrier(&format!("a{i}"), 100)); }
        for i in 0..5 { members.push(carrier(&format!("b{i}"), 150)); }
        let mut want = template.clone();
        want.splice(100..100, ins24.iter().copied());
        let got = consensus_of(&template, &members);
        assert_eq!(got.len(), 224, "one 24-bp copy, not two");
        assert_eq!(got, want);
    }

    #[test]
    fn significance_keeps_supported_variants_and_merges_singletons() {
        assert!(variant_is_real(5, 200, 3, 0.001, 0.01));     // 5 reads carrying 3 distinguishing bases: not error
        assert!(!variant_is_real(2, 200, 1, 0.001, 0.01));    // 2 reads, 1 column: P(X>=2) with p=0.001, n=202 ~ 0.018 >= 0.01 -> merge
        assert!(!variant_is_real(3, 10, 0, 0.001, 0.01));
    }
    #[test]
    fn significance_is_the_binomial_tail_of_the_small_cluster() {
        // exact P(X >= n_small), X ~ Binomial(n_small + n_large, eps^k), eps = 0.001, from rational arithmetic outside this file: the alpha pairs
        // bracket each value, so the p-value is pinned to about three digits (including n = 1000, where the log-factorial table matters)
        let real = |ns: usize, nl: usize, k: usize, alpha: f64| variant_is_real(ns, nl, k, 0.001, alpha);
        assert!(!real(2, 200, 1, 0.0177) && real(2, 200, 1, 0.0178));          // 0.0177860
        assert!(!real(1, 200, 1, 0.1821) && real(1, 200, 1, 0.1822));          // 0.1821698
        assert!(!real(3, 200, 1, 0.00118) && real(3, 200, 1, 0.00119));        // 0.00118318
        assert!(!real(2, 200, 2, 2.02e-8) && real(2, 200, 2, 2.04e-8));        // eps^2 = 1e-6: 2.02983e-8
        assert!(!real(30, 970, 1, 9.4e-34) && real(30, 970, 1, 9.6e-34));      // n = 1000: 9.5031e-34
        // no distinguishing column or no reads is never a real variant; a vanishing p (underflow of eps^k) is real; alpha = 0 admits nothing
        assert!(!real(5, 200, 0, 0.5) && !real(0, 200, 3, 0.5));
        assert!(real(2, 200, 400, 0.01) && !real(2, 200, 400, 0.0));
        // more supporting reads or more columns only make a variant more real
        assert!(!real(2, 200, 1, 0.01) && real(3, 200, 1, 0.01) && real(2, 200, 2, 0.01));
    }
    #[test]
    fn distinguishing_columns_count_substitutions_only() {
        assert_eq!(distinguishing_columns(&parse_cs(":10*ag:5+acgt:3-tt~gt100ag*ct:2").unwrap()), 2);
        assert_eq!(distinguishing_columns(&parse_cs(":50+aaaa-cc:20").unwrap()), 0);
        assert_eq!(distinguishing_columns(&[]), 0);
    }

    #[test]
    fn refine_splits_the_members_that_do_not_fit_the_consensus() {
        let consensus = rand_seq(3000, 3);
        let hit = |de: f64, span: usize| PafHit { q: "m".into(), qlen: 3000, qs: 0, qe: span, strand: b'+', t: "cons".into(), tlen: 3000, ts: 0, te: span, matches: span, block: span, de, cs: None };
        let d = 0.00958;
        // fits; de exactly delta fits; de above delta splits (also between delta and 2 delta); half of the shorter read fits; under half splits
        let hits = [hit(0.002, 3000), hit(d, 3000), hit(0.021, 3000), hit(0.001, 1500), hit(0.001, 1499), hit(0.015, 3000)];
        let read: &[u8] = b"";
        let members: Vec<(&[u8], &PafHit)> = hits.iter().map(|h| (read, h)).collect();
        assert_eq!(refine_cluster(&consensus, &members, d), (vec![0, 1, 3], vec![2, 4, 5]));
        assert_eq!(refine_cluster(&consensus, &members[..2], d), (vec![0, 1], vec![]));              // nothing to split
        assert_eq!(refine_cluster(&consensus, &[], d), (vec![], vec![]));
    }

    // ---- flag / link / merge fates, the flag floor and the exon-union representative (spec §5.6-§5.7) ---------------------------------------

    fn cluster(n_reads: usize) -> ClusterSeq { ClusterSeq { family: "F".into(), id: "F:c0".into(), n_reads, seq: vec![b'A'; 1000] } }
    /// A genome hit of the 1000 bp consensus `F:c0` on chr7:500-1500 with the given matches, block and query span.
    fn genome_hit(matches: usize, block: usize, qs: usize, qe: usize) -> PafHit {
        PafHit { q: "F:c0".into(), qlen: 1000, qs, qe, strand: b'+', t: "chr7".into(), tlen: 2_000_000, ts: 500, te: 1500, matches, block, de: 0.01, cs: None }
    }

    #[test]
    fn linked_cluster_never_enters_a_component() {
        let c = ClusterSeq { family: "F".into(), id: "F:c0".into(), n_reads: 5, seq: vec![b'A'; 1000] };
        let h = PafHit { q: "F:c0".into(), qlen: 1000, qs: 0, qe: 1000, strand: b'+', t: "chr1".into(), tlen: 1_000_000, ts: 10, te: 1010, matches: 995, block: 1000, de: 0.005, cs: None };
        assert!(matches!(classify(&c, Some(&h), 0.00958), Fate::Linked { .. }));
        let h2 = PafHit { matches: 950, de: 0.05, ..h.clone() };
        assert!(matches!(classify(&c, Some(&h2), 0.00958), Fate::NewCopy { .. }));
        assert!(matches!(classify(&c, None, 0.00958), Fate::NewCopy { .. }));
    }
    #[test]
    fn classify_tests_the_reference_first_then_the_link_distance_and_names_the_nearest_locus() {
        let (c, delta) = (cluster(5), 0.00958);
        let linked = |f: Fate| match f { Fate::Linked { locus, d } => Some((locus, d)), _ => None };
        // identity x coverage = matches / block x span / qlen; exactly 0.999 is already in the reference (>=), although its d = 0.001 would link it
        assert_eq!(classify(&c, Some(&genome_hit(999, 1000, 0, 1000)), delta), Fate::InReference);
        let (locus, d) = linked(classify(&c, Some(&genome_hit(998, 1000, 0, 1000)), delta)).expect("0.998 is not in the reference: linked");
        assert_eq!(locus, "chr7:500-1500");                                  // `chrom:ts-te` of the hit
        assert!((d - 0.002).abs() < 1e-12);                                   // 1 - matches / qlen
        // identity and coverage both count for the reference test: 999 matches in a 1100 bp block, or a clean hit over 995 of 1000 bases, are not in it
        assert!(linked(classify(&c, Some(&genome_hit(999, 1100, 0, 1000)), delta)).is_some());
        assert!(linked(classify(&c, Some(&genome_hit(995, 995, 5, 1000)), delta)).is_some());
        // d is 1 - matches / qlen whatever the block: 600 clean matches over 600 bases is 40% divergent, a new copy beside that locus
        match classify(&c, Some(&genome_hit(600, 600, 0, 600)), delta) {
            Fate::NewCopy { nearest, d } => { assert_eq!(nearest, "chr7:500-1500"); assert!((d - 0.4).abs() < 1e-12); }
            f => panic!("{f:?}"),
        }
        // d == delta links, a hair below it does not
        let h = genome_hit(990, 1000, 0, 1000);
        let d = whole_length_d(&h);
        assert!(matches!(classify(&c, Some(&h), d), Fate::Linked { .. }));
        assert!(matches!(classify(&c, Some(&h), d - 1e-9), Fate::NewCopy { .. }));
        // no genome hit at all: a new copy with no nearest locus and d = 1, even when delta would admit 1
        assert_eq!(classify(&c, None, delta), Fate::NewCopy { nearest: "none".into(), d: 1.0 });
        assert_eq!(classify(&c, None, 5.0), Fate::NewCopy { nearest: "none".into(), d: 1.0 });
    }
    #[test]
    fn components_join_new_copy_consensus_sequences_by_the_clustering_rule() {
        let (ids, d) = (names_of("n0 n1 n2 n3 n4"), 0.00958);
        let ava = vec![
            ava_hit("n0", "n1", 2900, 1.0, 0.002),                 // joined
            ava_hit("n1", "n2", 2900, 1.0, 0.015),                 // too divergent: between delta and 2 delta
            ava_hit("n2", "n3", 2900, 0.5, d),                     // exactly half of the shorter sequence and de == delta: joined
            ava_hit("n3", "n4", 2900, 0.49, 0.001),                // just under half
        ];
        assert_eq!(components(&ids, &ava, d), vec![vec![0, 1], vec![2, 3], vec![4]]);
        assert_eq!(components(&ids, &[], d), vec![vec![0], vec![1], vec![2], vec![3], vec![4]]);
        assert!(components(&[], &[], d).is_empty());
    }
    #[test]
    fn best_by_id_cov_prefers_identity_times_coverage_over_matches() {
        let hit = |q: &str, t: &str, matches: usize, block: usize, qs: usize, qe: usize| PafHit { q: q.into(), t: t.into(), matches, block, qs, qe, ..paf_hit(q, t, matches) };
        let hits = vec![
            hit("q", "wide", 900, 1100, 0, 1000),                  // the longer hit (1000 bases, 900 matches) but gappy: 900 / 1100 x 1.00 = 0.818
            hit("q", "clean", 850, 860, 40, 900),                  // the shorter, cleaner one: 850 / 860 x 0.86 = 0.850
            hit("q", "tiny", 300, 300, 0, 300),                    // perfect identity over 30% of the query: 0.300, so identity alone is not the rule either
            hit("r", "first", 800, 800, 0, 800), hit("r", "second", 800, 800, 0, 800),   // equal identity x coverage
            hit("s", "only", 10, 100, 0, 100),
        ];
        let best = best_by_id_cov(&hits);
        assert_eq!(best.len(), 3);
        assert_eq!(best["q"].t, "clean");                          // identity x coverage, not the number of matches ...
        assert_eq!(best_by_matches(&hits)["q"].t, "wide");         // ... which is what the pairwise rule keeps
        assert_eq!(best["r"].t, "first");                          // a tie keeps the first encountered
        assert_eq!(best["s"].t, "only");
        assert!(best_by_id_cov(&[]).is_empty());
    }
    #[test]
    fn is_flagged_sums_the_reads_of_the_component_against_the_floor() {
        let (a, b, c) = (cluster(3), cluster(3), cluster(2));
        assert!(is_flagged(&[&a, &b], 6));                         // 3 + 3 = 6: exactly the floor
        assert!(!is_flagged(&[&a, &c], 6));                        // 3 + 2 = 5
        assert!(is_flagged(&[&cluster(6)], 6) && !is_flagged(&[&cluster(5)], 6));   // the floor counts reads, not clusters: one cluster of 6 is enough
        assert!(!is_flagged(&[&a, &b, &c], 9) && is_flagged(&[&a, &b, &c], 8));
        assert!(!is_flagged(&[], 6));
        assert!(is_flagged(&[&c], 0));
    }

    /// Concatenation of byte slices.
    fn cat(parts: &[&[u8]]) -> Vec<u8> { parts.concat() }
    /// A hit of a member (query) on a union (target) as minimap2 `-x splice:hq -uf -c --cs` (`MM2_UNION`) writes it: `+` strand; matches, block
    /// and de are fillers.
    fn uhit(qlen: usize, qs: usize, qe: usize, tlen: usize, ts: usize, te: usize, cs: &str) -> PafHit {
        PafHit { q: "m".into(), qlen, qs, qe, strand: b'+', t: "u".into(), tlen, ts, te, matches: te.saturating_sub(ts), block: qe.saturating_sub(qs), de: 0.0, cs: Some(cs.into()) }
    }
    /// The lower-case string minimap2 writes inside a `cs`.
    fn lc(s: &[u8]) -> String { String::from_utf8(s.to_ascii_lowercase()).unwrap() }

    #[test]
    fn union_contains_each_exon_once() {
        let e1 = b"ACGTACGTACGTACGTACGTACGT".to_vec(); let e2 = b"TTGACCATGACCATGACCATGACC".to_vec(); let e3 = b"GGCATTGGCATTGGCATTGGCATT".to_vec();
        let iso_a: Vec<u8> = [e1.clone(), e3.clone()].concat();                 // skips e2
        let iso_b: Vec<u8> = [e1.clone(), e2.clone(), e3.clone()].concat();
        // the closure plays minimap2: iso_b vs union(iso_a) has a 24-bp insertion after e1
        let u = union_sequence(&[iso_b.clone(), iso_a.clone()], |_m, _u| None).unwrap();  // longest first: nothing to add from iso_a
        assert_eq!(u, iso_b);
        let u2 = union_sequence(&[iso_a.clone(), iso_b.clone()], |_m, _u| Some(PafHit { q: "b".into(), qlen: 72, qs: 0, qe: 72, strand: b'+', t: "u".into(), tlen: 48, ts: 0, te: 48, matches: 48, block: 72, de: 0.0, cs: Some(format!(":24+{}:24", String::from_utf8_lossy(&e2).to_lowercase())) })).unwrap();
        assert_eq!(u2, iso_b);
    }
    #[test]
    fn union_reads_a_hit_parsed_from_paf_text() {
        let (e1, e2, e3) = (rand_seq(24, 7), rand_seq(24, 8), rand_seq(24, 9));
        let (a, b) = (cat(&[&e1, &e3]), cat(&[&e1, &e2, &e3]));
        let line = format!("b\t72\t0\t72\t+\tu\t48\t0\t48\t48\t72\t60\tNM:i:24\tde:f:0.0\tcs:Z::24+{}:24\n", lc(&e2));
        let hit = parse_paf(&line).remove(0);
        assert_eq!(union_sequence(&[a, b.clone()], |_, _| Some(hit.clone())).unwrap(), b);
    }
    #[test]
    fn union_takes_unaligned_ends_in_front_of_ts_and_te() {
        let (a, core, c) = (rand_seq(30, 1), rand_seq(100, 2), rand_seq(30, 3));
        let (pre, suf) = (rand_seq(25, 4), rand_seq(22, 5));
        let one = |hit: PafHit, backbone: &[u8], member: &[u8]| union_sequence(&[backbone.to_vec(), member.to_vec()], |_, _| Some(hit.clone())).unwrap();
        // the member aligns inside the union (ts = 30, te = 130 of 160): its prefix goes in front of ts and its suffix in front of te, not to the ends
        let (u, m) = (cat(&[&a, &core, &c]), cat(&[&pre, &core, &suf]));
        assert_eq!(one(uhit(147, 25, 125, 160, 30, 130, ":100"), &u, &m), cat(&[&a, &pre, &core, &suf, &c]));
        // from the union's own ends (ts = 0, te = tlen) they are prepended and appended
        let (u, m) = (cat(&[&core, &c]), cat(&[&pre, &core]));
        assert_eq!(one(uhit(125, 25, 125, 130, 0, 100, ":100"), &u, &m), cat(&[&pre, &core, &c]));
        let (u, m) = (cat(&[&a, &core]), cat(&[&core, &suf]));
        assert_eq!(one(uhit(122, 0, 100, 130, 30, 130, ":100"), &u, &m), cat(&[&a, &core, &suf]));
        // both at once
        let (u, m) = (core.clone(), cat(&[&pre, &core, &suf]));
        assert_eq!(one(uhit(147, 25, 125, 100, 0, 100, ":100"), &u, &m), m);
        // the output is upper case whatever the case of the backbone or of the member the ends come from
        let lower = core.to_ascii_lowercase();
        assert_eq!(one(uhit(125, 25, 125, 100, 0, 100, ":100"), &lower, &cat(&[&pre, &core])), cat(&[&pre, &core]));
        assert_eq!(one(uhit(147, 25, 125, 100, 0, 100, ":100"), &u, &m.to_ascii_lowercase()), m);
    }
    #[test]
    fn union_takes_a_piece_only_from_20_bp() {
        let (core, e1) = (rand_seq(100, 2), 40usize);
        for len in [19usize, 20, 21] {
            let (x, take) = (rand_seq(len, 60 + len as u64), len >= 20);
            let run = |hit: PafHit, member: &[u8]| union_sequence(&[core.clone(), member.to_vec()], |_, _| Some(hit.clone())).unwrap();
            // an insertion of the cs
            let m = cat(&[&core[..e1], &x, &core[e1..]]);
            assert_eq!(run(uhit(100 + len, 0, 100 + len, 100, 0, 100, &format!(":{e1}+{}:{}", lc(&x), 100 - e1)), &m), if take { m.clone() } else { core.clone() }, "insertion of {len} bp");
            // an unaligned prefix, an unaligned suffix
            let m = cat(&[&x, &core]);
            assert_eq!(run(uhit(len + 100, len, len + 100, 100, 0, 100, ":100"), &m), if take { m.clone() } else { core.clone() }, "prefix of {len} bp");
            let m = cat(&[&core, &x]);
            assert_eq!(run(uhit(100 + len, 0, 100, 100, 0, 100, ":100"), &m), if take { m.clone() } else { core.clone() }, "suffix of {len} bp");
        }
    }
    #[test]
    fn union_applies_the_pieces_of_one_member_from_the_highest_position_and_keeps_their_order_at_one_position() {
        let ex: Vec<Vec<u8>> = (0..5).map(|i| rand_seq(24 + 3 * i, 20 + i as u64)).collect();       // 24, 27, 30, 33, 36 bp
        let (fl, fr) = (rand_seq(30, 11), rand_seq(20, 12));
        // union = flank e0 e2 e4 flank, member e0 e1 e2 e3 e4: two insertions, at 30 + 24 and 30 + 24 + 30 of the union; the alignment starts at ts = 30
        let u = cat(&[&fl, &ex[0], &ex[2], &ex[4], &fr]);
        let m = cat(&[&ex[0], &ex[1], &ex[2], &ex[3], &ex[4]]);
        let cs = format!(":{}+{}:{}+{}:{}", ex[0].len(), lc(&ex[1]), ex[2].len(), lc(&ex[3]), ex[4].len());
        let hit = uhit(m.len(), 0, m.len(), u.len(), fl.len(), u.len() - fr.len(), &cs);
        assert_eq!(union_sequence(&[u, m.clone()], |_, _| Some(hit.clone())).unwrap(), cat(&[&fl, &m, &fr]));
        // two pieces at ONE position stay in member order: an unaligned prefix then a cs that starts with an insertion; an insertion that ends
        // the cs then an unaligned suffix (minimap2 does not write such alignments; the rule must still be deterministic)
        let (p, x, y, s, core) = (rand_seq(21, 31), rand_seq(22, 32), rand_seq(23, 33), rand_seq(24, 34), rand_seq(60, 35));
        let m = cat(&[&p, &x, &core, &y, &s]);
        let hit = uhit(m.len(), p.len(), p.len() + x.len() + 60 + y.len(), 60, 0, 60, &format!("+{}:60+{}", lc(&x), lc(&y)));
        assert_eq!(union_sequence(&[core, m.clone()], |_, _| Some(hit.clone())).unwrap(), m);
    }
    #[test]
    fn union_tracks_the_target_position_through_every_cs_operation() {
        // a 158 bp union; the member has, in order: a substitution at 10, a 2 bp deletion at 16..18, an intron of 100 bases at 26..126 and a 25 bp
        // insertion in front of 138. The walk lands on 138 only if every operation advances the target by its own length (and starts at ts).
        let (u, x) = (rand_seq(158, 77), rand_seq(25, 78));
        let alt = other_base(u[10], 0);
        let m = cat(&[&u[..10], &[alt], &u[11..16], &u[18..26], &u[126..138], &x, &u[138..]]);
        assert_eq!(m.len(), 10 + 1 + 5 + 8 + 12 + 25 + 20);
        let cs = format!(":10*{}{}:5-{}:8~gt100ag:12+{}:20", lc(&u[10..11]), lc(&[alt]), lc(&u[16..18]), lc(&x));
        let want = cat(&[&u[..138], &x, &u[138..]]);
        assert_eq!(union_sequence(&[u.clone(), m.clone()], |_, _| Some(uhit(m.len(), 0, m.len(), 158, 0, 158, &cs))).unwrap(), want);
        // the same member against a union that has a 12 bp flank on the left (ts = 12) and a 9 bp one on the right (te = tlen - 9)
        let (fl, fr) = (rand_seq(12, 79), rand_seq(9, 80));
        let uf = cat(&[&fl, &u, &fr]);
        assert_eq!(union_sequence(&[uf.clone(), m.clone()], |_, _| Some(uhit(m.len(), 0, m.len(), uf.len(), 12, 12 + 158, &cs))).unwrap(), cat(&[&fl, &want, &fr]));
    }
    #[test]
    fn union_aligns_each_member_to_the_union_as_it_has_grown() {
        let ex: Vec<Vec<u8>> = (0..4).map(|i| rand_seq(25 + i, 40 + i as u64)).collect();           // e0..e3, 25..28 bp
        let (b, m1, m2) = (cat(&[&ex[0], &ex[3]]), cat(&[&ex[0], &ex[1], &ex[3]]), cat(&[&ex[0], &ex[2], &ex[3]]));
        let mut seen = vec![];
        let got = union_sequence(&[b.clone(), m1.clone(), m2.clone()], |m, u| {
            seen.push(u.len());
            // m1 vs e0 e3: an insertion of e1; m2 vs e0 e1 e3 (the grown union): the cs deletes e1 and inserts e2
            let cs = if m == &m1[..] { format!(":{}+{}:{}", ex[0].len(), lc(&ex[1]), ex[3].len()) } else { format!(":{}-{}+{}:{}", ex[0].len(), lc(&ex[1]), lc(&ex[2]), ex[3].len()) };
            Some(uhit(m.len(), 0, m.len(), u.len(), 0, u.len(), &cs))
        }).unwrap();
        assert_eq!(seen, vec![b.len(), b.len() + ex[1].len()]);                                      // the second member met the union with e1 in it
        assert_eq!(got, cat(&[&ex[0], &ex[1], &ex[2], &ex[3]]));
    }
    #[test]
    fn union_skips_members_without_a_plus_hit_and_counts_them() {
        let (core, x) = (rand_seq(60, 3), rand_seq(25, 4));
        let m = cat(&[&core[..30], &x, &core[30..]]);                                                // the union with a 25 bp exon inside
        let hit = |strand: u8| PafHit { strand, ..uhit(85, 0, 85, 60, 0, 60, &format!(":30+{}:30", lc(&x))) };
        // a '-' hit, no hit, then a '+' hit: only the last one adds the exon, and the first two are counted
        let mut results = vec![Some(hit(b'-')), None, Some(hit(b'+'))].into_iter();
        let members = vec![core.clone(), m.clone(), m.clone(), m.clone()];
        let (u, note) = union_sequence_with_note(&members, |_, _| results.next().unwrap()).unwrap();
        assert_eq!((u, note), (m.clone(), UnionNote { no_hit: 1, skipped_minus: 1 }));
        // the '-' hit alone adds nothing, though the same hit on the '+' strand would add the exon
        let (u, note) = union_sequence_with_note(&members[..2], |_, _| Some(hit(b'-'))).unwrap();
        assert_eq!((u, note), (core.clone(), UnionNote { no_hit: 0, skipped_minus: 1 }));
        assert_eq!(union_sequence_with_note(&members[..2], |_, _| Some(hit(b'.'))).unwrap().1, UnionNote { no_hit: 0, skipped_minus: 1 });   // anything but '+'
        assert_eq!(union_sequence(&members[..2], |_, _| Some(hit(b'+'))).unwrap(), m);
        assert_eq!(union_sequence_with_note(&members[..1], |_, _| unreachable!()).unwrap(), (core, UnionNote::default()));
    }
    #[test]
    fn union_takes_the_ends_of_a_hit_without_a_cs_and_nothing_inside_it() {
        let (core, x, suf) = (rand_seq(100, 2), rand_seq(30, 6), rand_seq(22, 5));
        let m = cat(&[&core[..40], &x, &core[40..], &suf]);                                          // an exon inside AND an unaligned suffix
        let hit = PafHit { cs: None, ..uhit(m.len(), 0, 130, 100, 0, 100, "") };
        assert_eq!(union_sequence(&[core.clone(), m], |_, _| Some(hit.clone())).unwrap(), cat(&[&core, &suf]));
    }
    #[test]
    fn union_of_nothing_and_of_one_member_never_calls_the_closure() {
        let never = |_: &[u8], _: &[u8]| -> Option<PafHit> { panic!("no member after the backbone: nothing to align") };
        assert_eq!(union_sequence(&[], never).unwrap(), Vec::<u8>::new());
        assert_eq!(union_sequence(&[b"acgtNacgt".to_vec()], never).unwrap(), b"ACGTNACGT".to_vec());   // upper case, like the consensus
    }
    #[test]
    fn union_rejects_a_hit_that_does_not_fit_its_sequences_and_a_malformed_cs() {
        let (core, x) = (rand_seq(60, 3), rand_seq(20, 4));
        let m = cat(&[&core, &x]);                                                                   // 80 bp: the union plus a 20 bp exon at its end
        let ok = uhit(80, 0, 80, 60, 0, 60, &format!(":60+{}", lc(&x)));
        assert_eq!(union_sequence(&[core.clone(), m.clone()], |_, _| Some(ok.clone())).unwrap(), m);   // the control: this hit is fine
        let bad = [
            ("query length", PafHit { qlen: 81, ..ok.clone() }, "query length"),
            ("target length", PafHit { tlen: 61, ..ok.clone() }, "target length"),
            ("qe past qlen", PafHit { qe: 90, ..ok.clone() }, "query span"),
            ("qs after qe", PafHit { qs: 50, qe: 40, ..ok.clone() }, "query span"),
            ("te past tlen", PafHit { te: 70, ..ok.clone() }, "target span"),
            ("ts after te", PafHit { ts: 50, te: 40, ..ok.clone() }, "target span"),
            ("cs short of te", PafHit { cs: Some(format!(":59+{}", lc(&x))), ..ok.clone() }, "cs ends"),
            ("cs past te", PafHit { cs: Some(format!(":61+{}", lc(&x))), ..ok.clone() }, "cs ends"),
            ("malformed cs", PafHit { cs: Some(":60+".into()), ..ok.clone() }, "malformed cs"),
        ];
        for (what, hit, needle) in bad {
            let err = union_sequence(&[core.clone(), m.clone()], |_, _| Some(hit.clone())).expect_err(what).to_string();
            assert!(err.contains("member 1") && err.contains(needle), "{what}: {err}");
        }
    }

    /// A stand-in for minimap2 on exon-structured sequences (the property test): both sequences are concatenations of distinct exons of `pool`.
    /// A member exon found in the union after the previous aligned one is a match (`:len`; the union stretch between two aligned exons is a
    /// `-`), one the union lacks between two aligned exons is a `+`, and member exons before the first / after the last aligned one are the
    /// unaligned query ends (`qs`, `qlen - qe`). No shared exon: no hit.
    fn exon_hit(pool: &[Vec<u8>], member: &[u8], union: &[u8]) -> Option<PafHit> {
        let find = |hay: &[u8], e: &[u8]| hay.windows(e.len()).position(|w| w == e);
        let mut in_member: Vec<(usize, &Vec<u8>)> = pool.iter().filter_map(|e| find(member, e).map(|p| (p, e))).collect();
        in_member.sort_by_key(|x| x.0);
        let (mut cs, mut pending, mut first, mut t, mut qe) = (String::new(), String::new(), None::<(usize, usize)>, 0usize, 0usize);   // first = (qs, ts)
        for (q, e) in in_member {
            match find(union, e) {
                Some(p) if first.is_none() || p >= t => {
                    if first.is_none() { first = Some((q, p)); } else { cs += &pending; if p > t { cs += &format!("-{}", lc(&union[t..p])); } }
                    pending.clear();
                    cs += &format!(":{}", e.len()); t = p + e.len(); qe = q + e.len();
                }
                _ if first.is_some() => pending += &format!("+{}", lc(e)),
                _ => {}
            }
        }
        let (qs, ts) = first?;
        Some(uhit(member.len(), qs, qe, union.len(), ts, t, &cs))
    }
    #[test]
    fn union_of_random_collinear_isoforms_is_the_exon_pool_in_slot_order() {
        // slots in genomic order: a 5' exon, c0, o0, c1, o1, c2, o2, c3, a 3' exon. The c's are in every isoform, every other exon in a random half
        let sizes = [31usize, 40, 24, 52, 27, 45, 33, 60, 22];
        let pool: Vec<Vec<u8>> = sizes.iter().enumerate().map(|(i, &n)| rand_seq(n, 100 + i as u64)).collect();
        let mut x = 0x9E3779B97F4A7C15u64;
        let mut coin = move || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x & 8 != 0 };
        for trial in 0..200usize {
            let isoforms: Vec<Vec<usize>> = (0..2 + trial % 5).map(|_| (0..pool.len()).filter(|i| i % 2 == 1 || coin()).collect()).collect();
            let mut seqs: Vec<Vec<u8>> = isoforms.iter().map(|iso| iso.iter().flat_map(|&i| pool[i].iter().copied()).collect()).collect();
            seqs.sort_by(|a, b| b.len().cmp(&a.len()));                                              // longest first (a stable sort)
            let used: Vec<usize> = (0..pool.len()).filter(|i| isoforms.iter().any(|iso| iso.contains(i))).collect();
            let want: Vec<u8> = used.iter().flat_map(|&i| pool[i].iter().copied()).collect();
            let got = union_sequence(&seqs, |m, u| exon_hit(&pool, m, u)).unwrap();
            assert_eq!(got, want, "trial {trial}: isoforms {isoforms:?}");
        }
    }

    // ---- the minimap2 runner and the output writers (plan task 7) -----------------------------------------------------------------------

    fn cl(family: &str, id: &str, n_reads: usize, len: usize) -> ClusterSeq { ClusterSeq { family: family.into(), id: id.into(), n_reads, seq: vec![b'A'; len] } }
    /// A candidate of a net of 120 reads, 100 of them used.
    fn cand(family: &str, id: &str, clusters: Vec<ClusterSeq>, union: &[u8], flagged: bool, nearest: &str, d: f64) -> Candidate {
        Candidate { family: family.into(), id: id.into(), clusters, union: union.to_vec(), flagged, nearest: nearest.into(), d, n_net: 120, n_used: 100 }
    }
    fn prefix_in(dir: &tempfile::TempDir) -> String { dir.path().join("t.cand").to_str().unwrap().to_string() }
    fn text(prefix: &str, ext: &str) -> String { std::fs::read_to_string(format!("{prefix}.{ext}")).unwrap() }
    const OUTPUT_EXTS: [&str; 4] = ["candidates.tsv", "clusters.tsv", "contigs.fa", "nets.fa"];
    const CANDIDATES_HEADER_TEXT: &str = "family\tcandidate\tn_clusters\tn_reads\tflagged\tunion_len\tnearest_locus\td\tn_net\tn_used";
    const CLUSTERS_HEADER_TEXT: &str = "family\tcluster\tcandidate\tn_reads\tconsensus_len\tfate\tlinked_to\td";

    #[test]
    fn write_outputs_writes_the_four_files_with_the_registered_headers_and_a_contig_for_the_flagged_candidate_only() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = prefix_in(&dir);
        let cands = vec![
            cand("MCL0", "cand_MCL0_0", vec![cl("MCL0", "MCL0:c0", 4, 100), cl("MCL0", "MCL0:c1", 3, 101)], b"ACGTACGTAC", true, "chr1:100-1100", 0.0123456),
            cand("MCL0", "cand_MCL0_1", vec![cl("MCL0", "MCL0:c3", 2, 60)], b"GGGGCC", false, "none", 1.0),
        ];
        let linked = vec![(cl("MCL0", "MCL0:c2", 5, 90), "chr1:2000-3000".to_string(), 0.00321)];
        let nets = vec![
            ("MCL0".to_string(), vec![("r1".to_string(), b"ACGT".to_vec()), ("r2".to_string(), b"GGCC".to_vec())]),
            ("MCL3".to_string(), vec![("r9".to_string(), b"TTTT".to_vec())]),
        ];
        write_outputs(&prefix, &cands, &linked, &nets).unwrap();
        // one row per candidate, flagged or not; n_reads sums its clusters' reads; the flag is 1 / 0; floats with 5 decimals; n_net / n_used as given
        assert_eq!(text(&prefix, "candidates.tsv"), format!("{CANDIDATES_HEADER_TEXT}\n\
            MCL0\tcand_MCL0_0\t2\t7\t1\t10\tchr1:100-1100\t0.01235\t120\t100\n\
            MCL0\tcand_MCL0_1\t1\t2\t0\t6\tnone\t1.00000\t120\t100\n"));
        // one row per cluster: the candidates' members (fate new_copy, the candidate's d), then the linked clusters (candidate -, their own locus and d)
        assert_eq!(text(&prefix, "clusters.tsv"), format!("{CLUSTERS_HEADER_TEXT}\n\
            MCL0\tMCL0:c0\tcand_MCL0_0\t4\t100\tnew_copy\t-\t0.01235\n\
            MCL0\tMCL0:c1\tcand_MCL0_0\t3\t101\tnew_copy\t-\t0.01235\n\
            MCL0\tMCL0:c3\tcand_MCL0_1\t2\t60\tnew_copy\t-\t1.00000\n\
            MCL0\tMCL0:c2\t-\t5\t90\tlinked\tchr1:2000-3000\t0.00321\n"));
        // the contigs: the flagged candidate only, its union on one line; the nets: every read of every family, in the given order
        assert_eq!(text(&prefix, "contigs.fa"), ">cand_MCL0_0\nACGTACGTAC\n");
        assert_eq!(text(&prefix, "nets.fa"), ">r1\nACGT\n>r2\nGGCC\n>r9\nTTTT\n");
    }

    #[test]
    fn write_outputs_of_nothing_is_the_two_headers_and_two_empty_fastas() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = prefix_in(&dir);
        write_outputs(&prefix, &[], &[], &[]).unwrap();
        assert_eq!(text(&prefix, "candidates.tsv"), format!("{CANDIDATES_HEADER_TEXT}\n"));
        assert_eq!(text(&prefix, "clusters.tsv"), format!("{CLUSTERS_HEADER_TEXT}\n"));
        assert_eq!((text(&prefix, "contigs.fa"), text(&prefix, "nets.fa")), (String::new(), String::new()));
    }

    #[test]
    fn write_outputs_sorts_by_family_then_candidate_with_trailing_numbers_in_numeric_order_whatever_the_input_order() {
        assert_eq!(candidate_id("MCL0", 3), "cand_MCL0_3");
        let mk = |f: &str, k: usize| cand(f, &candidate_id(f, k), vec![cl(f, &format!("{f}:c{k}"), 6, 50), cl(f, &format!("{f}:b{k}"), 1, 40)], b"ACGTAC", true, "none", 1.0);
        let sorted = vec![mk("MCL1", 0), mk("MCL1", 2), mk("MCL1", 10), mk("MCL2", 0), mk("MCL10", 0)];
        let shuffled: Vec<Candidate> = [4usize, 2, 0, 3, 1].iter().map(|&i| sorted[i].clone()).collect();
        let linked = vec![(cl("MCL1", "MCL1:y", 2, 10), "chr1:3-4".to_string(), 0.002), (cl("MCL2", "MCL2:z", 2, 10), "chr1:1-2".to_string(), 0.001)];
        let reversed: Vec<(ClusterSeq, String, f64)> = linked.iter().rev().cloned().collect();
        let dir = tempfile::tempdir().unwrap();
        let (pa, pb) = (dir.path().join("a.cand").to_str().unwrap().to_string(), dir.path().join("b.cand").to_str().unwrap().to_string());
        write_outputs(&pa, &sorted, &linked, &[]).unwrap();
        write_outputs(&pb, &shuffled, &reversed, &[]).unwrap();
        for ext in OUTPUT_EXTS { assert_eq!(text(&pa, ext), text(&pb, ext), "{ext}: the input order must not show"); }
        let col = |s: &str, c: usize| -> Vec<String> { s.lines().skip(1).map(|l| l.split('\t').nth(c).unwrap().to_string()).collect() };
        let (cands_tsv, clusters_tsv, contigs) = (text(&pa, "candidates.tsv"), text(&pa, "clusters.tsv"), text(&pa, "contigs.fa"));
        assert_eq!(col(&cands_tsv, 0), ["MCL1", "MCL1", "MCL1", "MCL2", "MCL10"]);
        assert_eq!(col(&cands_tsv, 1), ["cand_MCL1_0", "cand_MCL1_2", "cand_MCL1_10", "cand_MCL2_0", "cand_MCL10_0"]);
        // clusters: by family, then the candidates in order with their clusters by id, a family's linked clusters last
        assert_eq!(col(&clusters_tsv, 1), ["MCL1:b0", "MCL1:c0", "MCL1:b2", "MCL1:c2", "MCL1:b10", "MCL1:c10", "MCL1:y", "MCL2:b0", "MCL2:c0", "MCL2:z", "MCL10:b0", "MCL10:c0"]);
        assert_eq!(contigs.lines().filter(|l| l.starts_with('>')).collect::<Vec<_>>(), [">cand_MCL1_0", ">cand_MCL1_2", ">cand_MCL1_10", ">cand_MCL2_0", ">cand_MCL10_0"]);
    }

    #[cfg(unix)]
    #[test]
    fn write_outputs_replaces_its_files_instead_of_writing_through_a_hard_link() {
        // run_cache replays products by hard link: a writer that truncated its old product in place would rewrite the cache entry's payload
        let dir = tempfile::tempdir().unwrap();
        let prefix = prefix_in(&dir);
        for ext in OUTPUT_EXTS {
            let payload = dir.path().join(format!("payload.{ext}"));
            std::fs::write(&payload, b"cached payload\n").unwrap();
            std::fs::hard_link(&payload, format!("{prefix}.{ext}")).unwrap();
        }
        write_outputs(&prefix, &[], &[], &[]).unwrap();
        for ext in OUTPUT_EXTS {
            assert_eq!(std::fs::read(dir.path().join(format!("payload.{ext}"))).unwrap(), b"cached payload\n", "{ext}: the linked payload must be untouched");
        }
        assert_eq!(text(&prefix, "candidates.tsv"), format!("{CANDIDATES_HEADER_TEXT}\n"));
    }

    #[test]
    fn write_outputs_refuses_inconsistent_input_before_writing_anything() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = prefix_in(&dir);
        let nothing_written = |p: &str| OUTPUT_EXTS.iter().all(|e| !Path::new(&format!("{p}.{e}")).exists());
        // a flagged candidate needs a union: it is its contig
        let no_union = cand("F", "cand_F_0", vec![cl("F", "F:c0", 6, 50)], b"", true, "none", 1.0);
        let msg = format!("{:#}", write_outputs(&prefix, &[no_union.clone()], &[], &[]).unwrap_err());
        assert!(msg.contains("cand_F_0") && msg.contains("union"), "{msg}");
        // a candidate holds the clusters of its own family only
        let mixed = cand("F", "cand_F_0", vec![cl("F", "F:c0", 6, 50), cl("G", "G:c0", 6, 50)], b"ACGT", true, "none", 1.0);
        let msg = format!("{:#}", write_outputs(&prefix, &[mixed], &[], &[]).unwrap_err());
        assert!(msg.contains("cand_F_0") && msg.contains("G:c0"), "{msg}");
        assert!(nothing_written(&prefix), "a refused input leaves no partial output");
        // an unflagged candidate has no contig, so an empty union is fine
        write_outputs(&prefix, &[Candidate { flagged: false, ..no_union }], &[], &[]).unwrap();
        assert_eq!(text(&prefix, "contigs.fa"), "");
        // an output that cannot be created names its path
        let bad = dir.path().join("no_such_dir").join("t.cand");
        let msg = format!("{:#}", write_outputs(bad.to_str().unwrap(), &[], &[], &[]).unwrap_err());
        assert!(msg.contains("no_such_dir"), "{msg}");
    }

    #[test]
    fn the_minimap2_argument_sets_are_the_registered_ones() {
        assert_eq!(MM2_AVA, ["-x", "asm20", "-c", "--cs", "--dual=no", "-N", "100", "-p", "0.1", "--secondary=yes"]);   // R11: not -X
        // prereg Amendment 13: the vote and refinement alignments use the splice preset (asm20 cut alignments at exon skips, as R6 found)
        assert_eq!(MM2_MEMBERS, ["-x", "splice:hq", "-uf", "-c", "--cs", "-N", "5", "-p", "0.5"]);
        assert_eq!(MM2_UNION, ["-x", "splice:hq", "-uf", "-c", "--cs", "-N", "5", "-p", "0.5"]);
        assert_eq!(MM2_GENOME, ["-x", "splice:hq", "-uf", "-c", "--eqx", "-N", "20"]);
        // prereg Amendment 13b: the attribution preset and the rule's values
        assert_eq!(MM2_ATTRIB, ["-x", "map-ont", "-c", "-N", "5", "-p", "0.5"]);
        assert_eq!((MIN_UNMAPPED_LEN, ATTRIB_MIN_READ_COV, ATTRIB_MAX_DE, POORLY_PLACED_DE), (300, 0.5, 0.20, 0.02));
        // prereg Amendment 13e: the template's eligibility caps the aligned partners needed at 50
        assert_eq!(ELIGIBLE_PARTNER_CAP, 50);
    }

    #[test]
    fn the_paf_key_covers_every_byte_of_the_target_and_the_query_the_command_and_the_build_and_nothing_else() {
        let dir = tempfile::tempdir().unwrap();
        let file = |name: &str, bytes: &[u8]| { let p = dir.path().join(name); std::fs::write(&p, bytes).unwrap(); p };
        let root = dir.path().join("cache");
        let entry = |t: &Path, q: &Path, args: &[&str]| paf_entry(&root, "/bin/false", args, t, None, q).unwrap();
        let (target, query) = (b">t\nACGTACGTACGTACGTACGT\n".to_vec(), b">q\nGGGGCCCCAAAATTTTGGGG\n".to_vec());
        let base = entry(&file("t.fa", &target), &file("q.fa", &query), MM2_AVA);
        assert!(base.key.starts_with("rustle o3 minimap2 v1\ncmd\t/bin/false -x asm20 -c --cs --dual=no -N 100 -p 0.1 --secondary=yes\nminimap2\t"), "{}", base.key);
        assert!(base.pin && base.dir.starts_with(root.join("paf")), "a pinned entry of kind paf: {}", base.dir.display());
        // the same bytes written again later (a new mtime) and under other names: the same key, so the same entry
        std::thread::sleep(std::time::Duration::from_millis(20));
        let again = entry(&file("t.fa", &target), &file("q.fa", &query), MM2_AVA);
        let renamed = entry(&file("t2.fa", &target), &file("q2.fa", &query), MM2_AVA);
        assert_eq!((&again.key, &again.dir), (&base.key, &base.dir));
        assert_eq!((&renamed.key, &renamed.dir), (&base.key, &base.dir));
        // any one byte of the query or of the target, the size unchanged: another key and another entry
        let flip = |b: &[u8], i: usize| { let mut v = b.to_vec(); v[i] = if v[i] == b'A' { b'C' } else { b'A' }; v };
        for i in 0..query.len() {
            let e = entry(&file("t.fa", &target), &file("qx.fa", &flip(&query, i)), MM2_AVA);
            assert!(e.key != base.key && e.dir != base.dir, "query byte {i}");
        }
        for i in 0..target.len() {
            let e = entry(&file("tx.fa", &flip(&target, i)), &file("q.fa", &query), MM2_AVA);
            assert!(e.key != base.key && e.dir != base.dir, "target byte {i}");
        }
        // one byte more or fewer; the two files in each other's roles; other arguments; another binary
        let longer: Vec<u8> = query.iter().copied().chain(*b"A").collect();
        assert_ne!(entry(&file("t.fa", &target), &file("ql.fa", &longer), MM2_AVA).key, base.key);
        assert_ne!(entry(&file("t.fa", &target), &file("qs.fa", &query[..query.len() - 1]), MM2_AVA).key, base.key);
        assert_ne!(entry(&file("q.fa", &query), &file("t.fa", &target), MM2_AVA).key, base.key);
        assert_ne!(entry(&file("t.fa", &target), &file("q.fa", &query), MM2_MEMBERS).key, base.key);
        assert_ne!(paf_entry(&root, "/bin/true", MM2_AVA, &file("t.fa", &target), None, &file("q.fa", &query)).unwrap().key, base.key);
        // the minimap2 build is part of the key
        let hash = |b: &[u8]| { let mut h = rc::ContentHash::default(); h.update(b); h };
        let (ht, hq) = (TargetId::Content(hash(&target)), hash(&query));
        assert_eq!(paf_key("minimap2 -x asm20", "2.30-r1287", &ht, &hq), paf_key("minimap2 -x asm20", "2.30-r1287", &TargetId::Content(hash(&target)), &hash(&query)));
        assert_ne!(paf_key("minimap2 -x asm20", "2.30-r1287", &ht, &hq), paf_key("minimap2 -x asm20", "2.28-r1209", &ht, &hq));
        // R7: a target named by the caller's key (the genome index's file fingerprint) is never read: an absent file is no error; the key text
        // carries the given key under its own label instead of a content hash; another key is another entry; the query still counts
        let absent = dir.path().join("absent.mmi");
        let keyed = |k: &str, q: &Path| paf_entry(&root, "/bin/false", MM2_GENOME, &absent, Some(k), q).unwrap();
        let k1 = keyed("/x/G.mmi\t13600000000\t1", &file("q.fa", &query));
        assert!(k1.key.contains("\ntarget_key\t/x/G.mmi\t13600000000\t1\nquery_hash\tcontent128:") && !k1.key.contains("target_hash"), "{}", k1.key);
        assert_eq!(keyed("/x/G.mmi\t13600000000\t1", &file("q2.fa", &query)).key, k1.key);
        assert_ne!(keyed("/x/G.mmi\t13600000000\t2", &file("q.fa", &query)).key, k1.key);
        assert_ne!(keyed("/x/G.mmi\t13600000000\t1", &file("ql.fa", &longer)).key, k1.key);
    }

    #[cfg(unix)]
    #[test]
    fn a_failing_minimap2_names_its_command_and_leaves_neither_a_paf_nor_a_cache_entry() {
        let dir = tempfile::tempdir().unwrap();
        let (t, q, out, root) = (dir.path().join("t.fa"), dir.path().join("q.fa"), dir.path().join("out.paf"), dir.path().join("cache"));
        std::fs::write(&t, b">t\nACGT\n").unwrap();
        std::fs::write(&q, b">q\nACGT\n").unwrap();
        for cache in [None, Some(root.as_path())] {
            std::fs::write(&out, b"an earlier product\n").unwrap();
            let msg = format!("{:#}", run_minimap2("/bin/false", MM2_AVA, &t, None, &q, &out, cache, 3).unwrap_err());
            for want in ["/bin/false -x asm20 -c --cs --dual=no -N 100 -p 0.1 --secondary=yes -t 3", t.to_str().unwrap(), q.to_str().unwrap()] { assert!(msg.contains(want), "{msg}"); }
            assert!(!out.exists(), "a failed run leaves no stale or partial PAF");
            assert!(!root.join("paf").exists(), "a failed run commits nothing to the cache");
        }
        // a binary that cannot be started: the error names it and the variable that overrides it
        let msg = format!("{:#}", run_minimap2("/nonexistent/minimap2", MM2_GENOME, &t, None, &q, &out, None, 1).unwrap_err());
        assert!(msg.contains("/nonexistent/minimap2") && msg.contains("RUSTLE_MINIMAP2"), "{msg}");
        assert!(!out.exists());
        // an input the cache key cannot hash is named too
        let msg = format!("{:#}", run_minimap2("/bin/false", MM2_GENOME, &dir.path().join("absent.mmi"), None, &q, &out, Some(&root), 1).unwrap_err());
        assert!(msg.contains("absent.mmi"), "{msg}");
    }

    #[cfg(unix)]
    #[test]
    fn minimap2_replays_a_hit_by_content_and_reruns_on_any_change_without_touching_the_cached_paf() {
        // `/bin/echo` stands for minimap2: it prints its command line, so a PAF says which files the run that wrote it was given
        let dir = tempfile::tempdir().unwrap();
        let file = |name: &str, bytes: &[u8]| { let p = dir.path().join(name); std::fs::write(&p, bytes).unwrap(); p };
        let (target, query) = (b">t\nACGTACGTACGT\n".to_vec(), b">q\nGGGGCCCCAAAA\n".to_vec());
        let mut query3 = query.clone();
        query3[5] = b'T';
        let (t1, q1, t2, q2, q3) = (file("t1.fa", &target), file("q1.fa", &query), file("t2.fa", &target), file("q2.fa", &query), file("q3.fa", &query3));
        let (out, out2, root) = (dir.path().join("out.paf"), dir.path().join("out2.paf"), dir.path().join("cache"));
        let run = |args: &[&str], t: &Path, q: &Path, o: &Path, cached: bool| run_minimap2("/bin/echo", args, t, None, q, o, cached.then_some(root.as_path()), 1).unwrap();
        let said = |o: &Path| std::fs::read_to_string(o).unwrap();
        let line = |args: &[&str], t: &Path, q: &Path| format!("{} -t 1 {} {}\n", args.join(" "), t.display(), q.display());
        // a miss runs the program, and its output is the PAF
        run(MM2_AVA, &t1, &q1, &out, true);
        assert_eq!(said(&out), line(MM2_AVA, &t1, &q1));
        // the same bytes under other names: a hit, replayed from the first run (its output names t1 and q1): the key is the content, not the paths
        run(MM2_AVA, &t2, &q2, &out2, true);
        assert_eq!(said(&out2), line(MM2_AVA, &t1, &q1));
        // one query byte changed: a miss, run again into the very path the first run's PAF was linked to ...
        run(MM2_AVA, &t2, &q3, &out, true);
        assert_eq!(said(&out), line(MM2_AVA, &t2, &q3));
        // ... which left the first run's cached PAF as it was
        run(MM2_AVA, &t2, &q2, &out, true);
        assert_eq!(said(&out), line(MM2_AVA, &t1, &q1));
        // other arguments, or the two files in each other's roles: misses; without a cache root nothing is replayed
        run(MM2_MEMBERS, &t2, &q2, &out, true);
        assert_eq!(said(&out), line(MM2_MEMBERS, &t2, &q2));
        run(MM2_AVA, &q2, &t2, &out, true);
        assert_eq!(said(&out), line(MM2_AVA, &q2, &t2));
        run(MM2_AVA, &t2, &q2, &out, false);
        assert_eq!(said(&out), line(MM2_AVA, &t2, &q2));
        // R10: the thread count is on the command line but not in the key: 4 threads replay the 1-thread run's entry; uncached, `-t 4` is passed
        run_minimap2("/bin/echo", MM2_AVA, &t2, None, &q2, &out, Some(&root), 4).unwrap();
        assert_eq!(said(&out), line(MM2_AVA, &t1, &q1));
        run_minimap2("/bin/echo", MM2_AVA, &t2, None, &q2, &out, None, 4).unwrap();
        assert_eq!(said(&out), format!("{} -t 4 {} {}\n", MM2_AVA.join(" "), t2.display(), q2.display()));
    }

    #[cfg(unix)]
    #[test]
    fn a_keyed_target_is_replayed_by_its_key_and_never_read() {
        // R7: the target (here a path that does not exist: `/bin/echo` never opens it) is named in the key by the caller's string
        let dir = tempfile::tempdir().unwrap();
        let (q, out, root) = (dir.path().join("q.fa"), dir.path().join("out.paf"), dir.path().join("cache"));
        std::fs::write(&q, b">q\nACGTACGT\n").unwrap();
        let (g1, g2) = (dir.path().join("G1.mmi"), dir.path().join("G2.mmi"));
        let run = |t: &Path, key: &str| { run_minimap2("/bin/echo", MM2_GENOME, t, Some(key), &q, &out, Some(&root), 2).unwrap(); std::fs::read_to_string(&out).unwrap() };
        let line = |t: &Path| format!("{} -t 2 {} {}\n", MM2_GENOME.join(" "), t.display(), q.display());
        assert_eq!(run(&g1, "fp-1"), line(&g1));                                                      // a miss runs
        assert_eq!(run(&g2, "fp-1"), line(&g1));                                                      // the same key: replayed, whatever the path
        assert_eq!(run(&g1, "fp-2"), line(&g1));                                                      // another key (a rebuilt index): a miss ...
        assert_eq!(run(&g2, "fp-2"), line(&g1));                                                      // ... stored under its own key
        assert_eq!(run(&g2, "fp-3"), line(&g2));
    }

    // ---- the binary's helpers: the net sample, the family and cluster-member tables (plan task 8) -------------------------------------

    #[test]
    fn sampling_is_deterministic() {
        // Review Focus 2: the sample of an over-cap net is the same across runs and machines. Pinned from an independent Python implementation of
        // splitmix64 (seed 1) + Fisher-Yates over the sorted names (j = next() % (i + 1) for i = n - 1 down to 1), not from this code's output.
        let names: Vec<String> = (0..40).map(|i| format!("read{i:02}")).collect();
        let pinned = sample_net(&names, 6);
        assert_eq!(pinned, ["read02", "read11", "read16", "read18", "read26", "read28"].map(String::from).to_vec());
        // the input order does not matter (the names are sorted first), nor does a repeat call
        let mut shuffled = names.clone();
        shuffled.reverse();
        shuffled.swap(3, 17);
        assert_eq!(sample_net(&shuffled, 6), pinned);
        for _ in 0..5 { assert_eq!(sample_net(&names, 6), pinned); }
        // at or under the cap: every name, sorted; a cap of 0: nothing
        assert_eq!(sample_net(&shuffled, 40), names);
        assert_eq!(sample_net(&shuffled[..5], 6), { let mut v = shuffled[..5].to_vec(); v.sort(); v });
        assert!(sample_net(&names, 0).is_empty());
        // a larger cap keeps another fixed set (same reference)
        let big: Vec<String> = [2, 4, 5, 6, 9, 10, 11, 12, 13, 15, 16, 17, 18, 20, 22, 24, 26, 28, 29, 30, 32, 34, 35, 37, 38].iter().map(|i| format!("read{i:02}")).collect();
        assert_eq!(sample_net(&names, 25), big);
    }

    #[test]
    fn write_family_table_writes_one_row_per_family_sorted_with_the_ruled_header() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = prefix_in(&dir);
        let row = |f: &str, n: usize| FamilyCounts { family: f.into(), n_net: n, n_used: n.min(1000), n_clusters: 3, n_in_reference: 1, n_linked: 1, n_new: 1, n_candidates: 1, n_flagged: 1 };
        write_family_table(&prefix, &[row("MCL10", 1500), FamilyCounts { family: "MCL2".into(), ..Default::default() }, row("MCL1", 40)]).unwrap();
        assert_eq!(text(&prefix, "families.tsv"), "family\tn_net\tn_used\tn_clusters\tn_in_reference\tn_linked\tn_new\tn_candidates\tn_flagged\n\
            MCL1\t40\t40\t3\t1\t1\t1\t1\t1\n\
            MCL2\t0\t0\t0\t0\t0\t0\t0\t0\n\
            MCL10\t1500\t1000\t3\t1\t1\t1\t1\t1\n");
        write_family_table(&prefix, &[]).unwrap();
        assert_eq!(text(&prefix, "families.tsv"), format!("{FAMILIES_HEADER}\n"));
    }

    #[test]
    fn write_cluster_members_lists_each_clusters_reads_and_its_consensus() {
        let dir = tempfile::tempdir().unwrap();
        let prefix = prefix_in(&dir);
        let (c10, c2, g0) = (ClusterSeq { seq: b"ACGT".to_vec(), ..cl("MCL0", "MCL0:c10", 2, 4) }, ClusterSeq { seq: b"GGA".to_vec(), ..cl("MCL0", "MCL0:c2", 2, 3) }, ClusterSeq { seq: b"T".to_vec(), ..cl("MCL1", "MCL1:c0", 1, 1) });
        let (r10, r2, r0) = (vec!["r9".to_string(), "r1".to_string()], vec!["r5".to_string(), "r3".to_string()], vec!["r1".to_string()]);
        write_cluster_members(&prefix, &[(&g0, &r0), (&c10, &r10), (&c2, &r2)]).unwrap();
        // clusters by family then id (c2 before c10), reads sorted inside a cluster; a read may be in two families' clusters
        assert_eq!(text(&prefix, "reads.tsv"), "read\tfamily\tcluster\nr3\tMCL0\tMCL0:c2\nr5\tMCL0\tMCL0:c2\nr1\tMCL0\tMCL0:c10\nr9\tMCL0\tMCL0:c10\nr1\tMCL1\tMCL1:c0\n");
        assert_eq!(text(&prefix, "clusters.fa"), ">MCL0:c2\nGGA\n>MCL0:c10\nACGT\n>MCL1:c0\nT\n");
    }

    // ---- canonical k-mers and sketches against independent references (Task 2's deferred items) ------------------------------------------

    /// A reproducible generator for the reference tests below (not the module's code paths).
    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 { let mut x = self.0; x ^= x << 13; x ^= x >> 7; x ^= x << 17; self.0 = x; x }
        fn below(&mut self, n: u64) -> u64 { (self.next() >> 11) % n }
        fn acgt(&mut self, n: usize) -> Vec<u8> { (0..n).map(|_| b"ACGT"[self.below(4) as usize]).collect() }
    }
    /// Canonical k-mers by explicit window slicing and an explicit reverse complement, with u128 accumulators: no rolling state.
    fn brute_canonical(seq: &[u8], k: usize) -> Vec<u64> {
        let code = |b: u8| match b.to_ascii_uppercase() { b'A' => Some(0u128), b'C' => Some(1), b'G' => Some(2), b'T' => Some(3), _ => None };
        (0..(seq.len() + 1).saturating_sub(k)).filter_map(|i| {
            let c: Vec<u128> = seq[i..i + k].iter().map(|&b| code(b)).collect::<Option<_>>()?;
            let fw = c.iter().fold(0u128, |a, &x| (a << 2) | x);
            let rv = c.iter().rev().fold(0u128, |a, &x| (a << 2) | (3 - x));
            Some(fw.min(rv) as u64)
        }).collect()
    }

    #[test]
    fn canonical_kmers_equal_a_brute_force_reference_with_n_lower_case_and_iupac_for_every_k() {
        let mut rng = Rng(0x1234_5678_9abc_def1);
        let (mut nonempty, mut with_other_bytes) = (0usize, 0usize);
        for k in [1usize, 2, 3, 7, 16, 21, 30, 31, 32] {
            for _ in 0..400 {
                // N runs, IUPAC codes and lower case mixed into random sequence of 0..160 bases
                let len = rng.below(160) as usize;
                let mut s = Vec::with_capacity(len);
                while s.len() < len {
                    match rng.below(100) {
                        0..=1 => { for _ in 0..1 + rng.below(40) { if s.len() < len { s.push(b'N'); } } }
                        2..=3 => s.push(b"RYKMSWnBDHV"[rng.below(11) as usize]),
                        _ => { let b = b"ACGT"[rng.below(4) as usize]; s.push(if rng.below(2) == 0 { b.to_ascii_lowercase() } else { b }); }
                    }
                }
                let got = canonical_kmers(&s, k);
                assert_eq!(got, brute_canonical(&s, k), "k = {k}, seq {}", String::from_utf8_lossy(&s));
                if !got.is_empty() { nonempty += 1; if s.iter().any(|b| !b"ACGTacgt".contains(b)) { with_other_bytes += 1; } }
            }
        }
        assert!(nonempty > 1000 && with_other_bytes > 200, "the cases must exercise the reference: {nonempty} / {with_other_bytes}");
    }

    #[test]
    fn canonical_kmers_edge_cases_and_the_parameter_asserts() {
        let mut rng = Rng(5);
        assert!(canonical_kmers(b"", 31).is_empty() && canonical_kmers(b"ACGT", 31).is_empty());
        assert!(canonical_kmers(&rng.acgt(30), 31).is_empty());
        assert_eq!(canonical_kmers(&rng.acgt(31), 31).len(), 1);
        let s32 = rng.acgt(32);
        assert_eq!(canonical_kmers(&s32, 32), brute_canonical(&s32, 32));                            // k = 32: the full-width mask
        assert_eq!(canonical_kmers(&s32, 32).len(), 1);
        // an N in the middle of a 71-mer leaves 5 windows on each side, none across it; an N run before a clean tail resets the window
        let mut v = rng.acgt(71); v[35] = b'N';
        assert_eq!(canonical_kmers(&v, 31).len(), 5 + 5);
        let mut t = b"NNN".to_vec(); t.extend(rng.acgt(100));
        assert_eq!(canonical_kmers(&t, 31).len(), 100 - 31 + 1);
        assert!(canonical_kmers(&[b'N'; 500], 31).is_empty());
        let mut sparse = rng.acgt(500); for i in (10..500).step_by(20) { sparse[i] = b'N'; }        // no 31-window free of N
        assert!(canonical_kmers(&sparse, 31).is_empty());
        // k outside 1..=32 and w = 0 are refused, not computed into garbage
        for bad in [0usize, 33, 64] { assert!(std::panic::catch_unwind(|| canonical_kmers(b"ACGTACGTAC", bad)).is_err(), "k = {bad}"); }
        assert!(std::panic::catch_unwind(|| minimizer_sketch(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGT", 31, 0)).is_err(), "w = 0");
        assert_eq!(minimizer_sketch(&rng.acgt(34), 31, 5).len(), 1);                                 // fewer k-mers than w: one minimizer
    }
}

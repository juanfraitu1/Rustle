//! Missing copies from RNA alone (thesis objective O3; §6ze, `docs/PREREG_o3_rna_only_2026-09-23.md`): everything that can be said about a possible
//! reference-absent copy from one BAM, stopping where only DNA can go (copy number).
//!
//! **STATUS:** OTHER-BINARY  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)
//!
//! Per locus: the per-read `de` divergence mixture (the held-out-validated S2 statistic, 08-14) splits the
//! pile into a HOST cluster and a divergent SUB-PILE; the sub-pile's mismatches are tested for CONSISTENCY
//! (a hidden copy's reads share the same sites — its PSVs — while error, editing and somatic drift do not);
//! a spliced consensus of the hidden copy is built by patching the reference at those sites; the consensus
//! is aligned to the whole primary genome (a better home elsewhere = an unannotated paralogue, not a missing
//! copy) and, optionally, to confirmation genomes (here the parental haplotypes); the screens that killed
//! every earlier candidate (immunoglobulin hypermutation, run-exclusive contamination, RNA editing) are
//! applied; the verdict is the first screen that fires, else `reference_absent_candidate`.
//!
//! Nothing here uses read depth to decide — expression is not dosage (RESULTS_DETECTOR.md); the expected
//! DNA depth ratio is reported as what DNA would have to show.

use crate::types::{DetHashMap, DetHashSet};
use std::collections::{BTreeMap};

/// One primary read of a locus pile: `de:f`, 0-based reference start, `--eqx` CIGAR ops, sequence
/// (soft-clips kept, as `AlignedRead::seq`), name.
#[derive(Clone, Debug)]
pub struct PileRead {
    pub name: String,
    pub de: f64,
    pub ref_start: u64,
    pub ops: Vec<(char, u64)>,
    pub seq: Vec<u8>,
}

/// Deterministic 1-D 2-means on a divergence vector — a line-for-line port of `detector.py::two_means`.
/// Returns `(m_high, median_high, delta, mid)`: the high cluster's fraction, its median, the median gap to
/// the low cluster, and the split threshold (`x <= mid` is low). `None` below 10 values.
pub fn two_means(v: &[f64]) -> Option<(f64, f64, f64, f64)> {
    if v.len() < 10 {
        return None;
    }
    let mut x: Vec<f64> = v.to_vec();
    x.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let (lo, hi) = (x[0], x[x.len() - 1]);
    if hi - lo < 1e-12 {
        return Some((0.0, x[x.len() - 1], 0.0, hi));
    }
    let (mut a, mut b) = (lo, hi);
    for _ in 0..50 {
        let mid = (a + b) / 2.0;
        let l: Vec<f64> = x.iter().copied().filter(|&t| t <= mid).collect();
        let h: Vec<f64> = x.iter().copied().filter(|&t| t > mid).collect();
        if l.is_empty() || h.is_empty() {
            break;
        }
        let (na, nb) = (mean(&l), mean(&h));
        if (na - a).abs() < 1e-12 && (nb - b).abs() < 1e-12 {
            break;
        }
        a = na;
        b = nb;
    }
    let mid = (a + b) / 2.0;
    let l: Vec<f64> = x.iter().copied().filter(|&t| t <= mid).collect();
    let h: Vec<f64> = x.iter().copied().filter(|&t| t > mid).collect();
    if h.is_empty() || l.is_empty() {
        return Some((0.0, x[x.len() - 1], 0.0, mid));
    }
    let (mh, ml) = (median_sorted(&h), median_sorted(&l));
    Some((h.len() as f64 / x.len() as f64, mh, mh - ml, mid))
}

fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}

/// Median of an already-sorted slice (Python `statistics.median` semantics).
fn median_sorted(v: &[f64]) -> f64 {
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

/// The mixture split of a pile. `sub`/`host` are indices into the pile.
#[derive(Clone, Debug)]
pub struct Split {
    pub m: f64,
    pub d_high: f64,
    pub delta: f64,
    pub sub: Vec<usize>,
    pub host: Vec<usize>,
}

/// Fire rule of the S2 detector: `m_min <= m <= 0.5`, `delta >= delta_min`, sub-pile >= `min_sub` reads.
pub fn split_pile(reads: &[PileRead], m_min: f64, delta_min: f64, min_sub: usize) -> Option<Split> {
    let de: Vec<f64> = reads.iter().map(|r| r.de).collect();
    let (m, d_high, delta, mid) = two_means(&de)?;
    if !(m >= m_min && m <= 0.5 && delta >= delta_min) {
        return None;
    }
    let sub: Vec<usize> = (0..reads.len()).filter(|&i| reads[i].de > mid).collect();
    let host: Vec<usize> = (0..reads.len()).filter(|&i| reads[i].de <= mid).collect();
    if sub.len() < min_sub {
        return None;
    }
    Some(Split { m, d_high, delta, sub, host })
}

/// Per-position tallies from `--eqx` CIGARs: coverage (aligned `=`/`X` bases), mismatches, and the
/// mismatching read base per position (majority taken later).
#[derive(Default)]
struct Tally {
    cov: DetHashMap<u64, u32>,
    mism: DetHashMap<u64, u32>,
    bases: DetHashMap<u64, DetHashMap<u8, u32>>,
    total_mism: u64,
}

fn tally(reads: &[&PileRead], keep_bases: bool) -> Tally {
    let mut t = Tally::default();
    for r in reads {
        let mut rp = r.ref_start;
        let mut qp = 0usize;
        for &(op, n) in &r.ops {
            match op {
                '=' => {
                    for k in 0..n {
                        *t.cov.entry(rp + k).or_insert(0) += 1;
                    }
                    rp += n;
                    qp += n as usize;
                }
                'X' => {
                    for k in 0..n {
                        *t.cov.entry(rp + k).or_insert(0) += 1;
                        *t.mism.entry(rp + k).or_insert(0) += 1;
                        if keep_bases {
                            if let Some(&b) = r.seq.get(qp + k as usize) {
                                *t.bases.entry(rp + k).or_default().entry(b.to_ascii_uppercase()).or_insert(0) += 1;
                            }
                        }
                    }
                    t.total_mism += n;
                    rp += n;
                    qp += n as usize;
                }
                'M' => {
                    // non-eqx CIGAR: counts as coverage only; such a BAM cannot yield PSV sites
                    for k in 0..n {
                        *t.cov.entry(rp + k).or_insert(0) += 1;
                    }
                    rp += n;
                    qp += n as usize;
                }
                'D' | 'N' => rp += n,
                'I' | 'S' => qp += n as usize,
                _ => {}
            }
        }
    }
    t
}

/// A PSV site: reference position, sub-pile majority base, and the tallies behind it.
#[derive(Clone, Debug, PartialEq)]
pub struct PsvSite {
    pub pos: u64,
    pub base: u8,
    pub n_sub_cov: u32,
    pub n_sub_mism: u32,
    pub n_host_cov: u32,
    pub n_host_mism: u32,
}

#[derive(Clone, Debug, Default)]
pub struct Consistency {
    pub sites: Vec<PsvSite>,
    /// sub-pile mismatches on PSV sites / all sub-pile mismatches
    pub shared_frac: f64,
    /// fraction of PSV substitutions that are A>G or T>C on the reference strand (editing-type)
    pub editing_frac: f64,
}

/// PSV sites as pre-registered: sub-pile coverage >= 3 with >= 80% mismatched, host coverage >= 3 with
/// <= 20% mismatched. `ref_seq` is the locus reference sequence starting at `ref_offset` (0-based).
pub fn consistency(sub: &[&PileRead], host: &[&PileRead], ref_seq: &[u8], ref_offset: u64) -> Consistency {
    let ts = tally(sub, true);
    let th = tally(host, false);
    let mut sites = Vec::new();
    let mut on_sites = 0u64;
    let mut editing = 0usize;
    for (&pos, &nm) in &ts.mism {
        let cs = *ts.cov.get(&pos).unwrap_or(&0);
        let ch = *th.cov.get(&pos).unwrap_or(&0);
        let mh = *th.mism.get(&pos).unwrap_or(&0);
        if cs >= 3 && (nm as f64) >= 0.8 * cs as f64 && ch >= 3 && (mh as f64) <= 0.2 * ch as f64 {
            let base = ts
                .bases
                .get(&pos)
                .and_then(|m| m.iter().max_by_key(|(b, c)| (**c, std::cmp::Reverse(**b))).map(|(b, _)| *b))
                .unwrap_or(b'N');
            on_sites += nm as u64;
            let rb = pos
                .checked_sub(ref_offset)
                .and_then(|i| ref_seq.get(i as usize))
                .map(|b| b.to_ascii_uppercase())
                .unwrap_or(b'N');
            if (rb == b'A' && base == b'G') || (rb == b'T' && base == b'C') {
                editing += 1;
            }
            sites.push(PsvSite { pos, base, n_sub_cov: cs, n_sub_mism: nm, n_host_cov: ch, n_host_mism: mh });
        }
    }
    sites.sort_by_key(|s| s.pos);
    let shared_frac = if ts.total_mism > 0 { on_sites as f64 / ts.total_mism as f64 } else { 0.0 };
    let editing_frac = if sites.is_empty() { 0.0 } else { editing as f64 / sites.len() as f64 };
    Consistency { sites, shared_frac, editing_frac }
}

/// The spliced patched consensus: the template read's `=`/`X` blocks taken from the reference and patched at
/// PSV sites. Returns `(sequence, blocks)`; blocks are 0-based half-open reference intervals.
pub fn patched_consensus(template: &PileRead, sites: &[PsvSite], ref_seq: &[u8], ref_offset: u64) -> (Vec<u8>, Vec<(u64, u64)>) {
    let patch: DetHashMap<u64, u8> = sites.iter().map(|s| (s.pos, s.base)).collect();
    let mut seq = Vec::new();
    let mut blocks: Vec<(u64, u64)> = Vec::new();
    let mut rp = template.ref_start;
    for &(op, n) in &template.ops {
        match op {
            '=' | 'X' | 'M' => {
                for k in 0..n {
                    let pos = rp + k;
                    let b = match patch.get(&pos) {
                        Some(&b) => b,
                        None => pos
                            .checked_sub(ref_offset)
                            .and_then(|i| ref_seq.get(i as usize))
                            .map(|b| b.to_ascii_uppercase())
                            .unwrap_or(b'N'),
                    };
                    seq.push(b);
                }
                match blocks.last_mut() {
                    Some(last) if last.1 == rp => last.1 = rp + n,
                    _ => blocks.push((rp, rp + n)),
                }
                rp += n;
            }
            'D' => {
                // a deletion in the template read: keep the reference bases (the consensus follows the reference
                // where the template is not informative) — extend the block, do not add sequence gaps
                for k in 0..n {
                    let pos = rp + k;
                    let b = pos
                        .checked_sub(ref_offset)
                        .and_then(|i| ref_seq.get(i as usize))
                        .map(|b| b.to_ascii_uppercase())
                        .unwrap_or(b'N');
                    seq.push(b);
                }
                match blocks.last_mut() {
                    Some(last) if last.1 == rp => last.1 = rp + n,
                    _ => blocks.push((rp, rp + n)),
                }
                rp += n;
            }
            'N' => rp += n,
            _ => {}
        }
    }
    (seq, blocks)
}

/// One read carrying a reference exon it also skips: minimap2 keeps the longest colinear exon run and encodes a
/// displaced exon as an INSERTION (addendum 2). `exon` is the matched stretch inside one of the read's intron gaps.
#[derive(Clone, Debug, PartialEq)]
pub struct Rearrangement {
    pub name: String,
    /// reference position of the insertion (0-based, the base after which the read inserts)
    pub ins_ref_pos: u64,
    pub ins_len: usize,
    pub exon: (u64, u64),
    pub identity: f64,
    /// the matched exon is ALSO covered by the read's aligned blocks: the exon appears twice in the read
    /// (tandem / rolling-circle duplication, the circRNA class), not a rearrangement
    pub duplicated: bool,
}

/// Best ungapped placement of `ins` inside `target`: 12-mer seed votes on the offset (refined ± 16), then the
/// longest stretch of the overlap with <= 10% mismatches. Returns `(target_start, matched_len, identity)` of that
/// stretch, or `None` when no stretch reaches `min_len`.
fn place_insertion(ins: &[u8], target: &[u8], min_len: usize) -> Option<(usize, usize, f64)> {
    const K: usize = 12;
    if ins.len() < K || target.len() < K {
        return None;
    }
    let mut index: DetHashMap<&[u8], Vec<usize>> = DetHashMap::default();
    for i in 0..=target.len() - K {
        index.entry(&target[i..i + K]).or_default().push(i);
    }
    let mut votes: DetHashMap<i64, u32> = DetHashMap::default();
    let mut i = 0;
    while i + K <= ins.len() {
        if let Some(ps) = index.get(&ins[i..i + K]) {
            for &p in ps.iter().take(8) {
                *votes.entry(p as i64 - i as i64).or_insert(0) += 1;
            }
        }
        i += 4;
    }
    let (&off0, _) = votes.iter().max_by_key(|(o, c)| (**c, std::cmp::Reverse(**o)))?;
    let mut best: Option<(usize, usize, f64)> = None;
    for off in off0 - 16..=off0 + 16 {
        let start = off.max(0) as usize;
        let ins_start = (start as i64 - off) as usize;
        if ins_start >= ins.len() || start >= target.len() {
            continue;
        }
        let n = (ins.len() - ins_start).min(target.len() - start);
        // mismatch prefix sums, then the longest window with mismatches <= 10% of its length
        let mut pre = vec![0u32; n + 1];
        for k in 0..n {
            pre[k + 1] = pre[k] + (!ins[ins_start + k].eq_ignore_ascii_case(&target[start + k])) as u32;
        }
        let (mut lo, mut run): (usize, Option<(usize, usize)>) = (0, None);
        for hi in 1..=n {
            while lo < hi && (pre[hi] - pre[lo]) as f64 > 0.10 * (hi - lo) as f64 {
                lo += 1;
            }
            if run.map_or(true, |(a, b)| hi - lo > b - a) {
                run = Some((lo, hi));
            }
        }
        if let Some((a, b)) = run {
            let len = b - a;
            let idn = 1.0 - (pre[b] - pre[a]) as f64 / len as f64;
            if len >= min_len && best.map_or(true, |x| len > x.1) {
                best = Some((start + a, len, idn));
            }
        }
    }
    best
}

/// Exon-order rearrangements among `reads` (addendum 2): every insertion >= `min_ins` is searched in the locus
/// reference (`ref_seq`, starting at `ref_offset`); a stretch of >= 50 bp at >= 0.90 identity that is not at the
/// insertion point itself is reference sequence the read carries OUT OF ORDER. If that stretch is also covered by
/// the read's aligned blocks the exon appears twice (`duplicated`: tandem / rolling-circle); otherwise the read
/// skips it (`rearranged`).
pub fn find_rearrangements(reads: &[&PileRead], ref_seq: &[u8], ref_offset: u64, min_ins: usize) -> Vec<Rearrangement> {
    let mut out = Vec::new();
    for r in reads {
        let (mut rp, mut qp) = (r.ref_start, 0usize);
        let mut blocks: Vec<(u64, u64)> = Vec::new();
        let mut inserts: Vec<(u64, usize, usize)> = Vec::new(); // (ref pos, query pos, len)
        for &(op, n) in &r.ops {
            match op {
                '=' | 'X' | 'M' | 'D' => {
                    match blocks.last_mut() {
                        Some(b) if b.1 == rp => b.1 = rp + n,
                        _ => blocks.push((rp, rp + n)),
                    }
                    rp += n;
                    if op != 'D' {
                        qp += n as usize;
                    }
                }
                'N' => rp += n,
                'I' => {
                    if n as usize >= min_ins {
                        inserts.push((rp, qp, n as usize));
                    }
                    qp += n as usize;
                }
                'S' => qp += n as usize,
                _ => {}
            }
        }
        for (ins_ref_pos, q, len) in inserts {
            let Some(ins) = r.seq.get(q..q + len) else { continue };
            let Some((start, mlen, idn)) = place_insertion(ins, ref_seq, 50) else { continue };
            let exon = (ref_offset + start as u64, ref_offset + (start + mlen) as u64);
            let covered: i64 = blocks.iter().map(|&(bs, be)| (be.min(exon.1) as i64 - bs.max(exon.0) as i64).max(0)).sum();
            let duplicated = covered >= (mlen as i64) / 2;
            out.push(Rearrangement { name: r.name.clone(), ins_ref_pos, ins_len: len, exon, identity: idn, duplicated });
        }
    }
    out
}

/// Clusters of `rearranged` (not duplicated) reads sharing the same matched exon (± `tol` bp) and insertion
/// site (± `tol` bp); returns clusters with >= `min_reads` reads, largest first.
pub fn rearrangement_clusters(rs: &[Rearrangement], tol: u64, min_reads: usize) -> Vec<Vec<Rearrangement>> {
    let mut clusters: Vec<Vec<Rearrangement>> = Vec::new();
    for r in rs.iter().filter(|r| !r.duplicated) {
        let near = |a: u64, b: u64| a.abs_diff(b) <= tol;
        match clusters.iter_mut().find(|c| near(c[0].exon.0, r.exon.0) && near(c[0].exon.1, r.exon.1) && near(c[0].ins_ref_pos, r.ins_ref_pos)) {
            Some(c) => c.push(r.clone()),
            None => clusters.push(vec![r.clone()]),
        }
    }
    let mut out: Vec<Vec<Rearrangement>> = clusters.into_iter().filter(|c| c.len() >= min_reads).collect();
    out.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a[0].exon.0.cmp(&b[0].exon.0)));
    out
}

/// The sub-pile read with the longest aligned reference span (deterministic tie-break on name).
pub fn template_read<'a>(sub: &[&'a PileRead]) -> Option<&'a PileRead> {
    sub.iter()
        .copied()
        .max_by(|a, b| aligned_span(a).cmp(&aligned_span(b)).then_with(|| b.name.cmp(&a.name)))
}

fn aligned_span(r: &PileRead) -> u64 {
    r.ops.iter().filter(|(op, _)| matches!(op, '=' | 'X' | 'M' | 'D')).map(|(_, n)| *n).sum()
}

/// Sequencing-run label of a read name: the prefix before the first `.` or `/` when it looks like a run id
/// (SRA `SRR27438212.n`, PacBio movie `m64076_221110_210557/zmw/ccs`: one to four letters then >= 3 digits);
/// otherwise `None` — a simulated or renamed BAM has no run structure and the screen must stay inert.
pub fn run_of(name: &str) -> Option<&str> {
    let cut = name.find(|c| c == '.' || c == '/').unwrap_or(name.len());
    let p = &name[..cut];
    let letters = p.bytes().take_while(|b| b.is_ascii_alphabetic()).count();
    let digits = p.bytes().skip(letters).take_while(|b| b.is_ascii_digit()).count();
    if (1..=4).contains(&letters) && digits >= 3 {
        Some(p)
    } else {
        None
    }
}

/// Two-sided Fisher exact test on [[a, b], [c, d]] via the hypergeometric distribution (log-factorials).
pub fn fisher_exact(a: u64, b: u64, c: u64, d: u64) -> f64 {
    let n = a + b + c + d;
    if n == 0 {
        return 1.0;
    }
    let (r1, c1) = (a + b, a + c);
    let lf = |k: u64| -> f64 { (1..=k).map(|i| (i as f64).ln()).sum() };
    let lp = |x: u64| -> f64 {
        // P(X = x) for X ~ Hypergeom(n, r1, c1)
        let (bb, cc) = (r1 - x, c1 - x);
        let dd = n - r1 - c1 + x;
        lf(r1) + lf(n - r1) + lf(c1) + lf(n - c1) - lf(n) - lf(x) - lf(bb) - lf(cc) - lf(dd)
    };
    let lo = r1.saturating_sub(n - c1);
    let hi = r1.min(c1);
    let p_obs = lp(a);
    let mut p = 0.0;
    for x in lo..=hi {
        let px = lp(x);
        if px <= p_obs + 1e-9 {
            p += px.exp();
        }
    }
    p.min(1.0)
}

/// Run-exclusivity screen: (p, top run, fraction of the sub-pile on that run, n distinct runs in the pile).
pub fn run_screen(sub: &[&PileRead], host: &[&PileRead]) -> (f64, String, f64, usize) {
    let mut runs: BTreeMap<&str, (u64, u64)> = BTreeMap::new();
    let mut unlabelled = 0usize;
    for r in sub {
        match run_of(&r.name) {
            Some(k) => runs.entry(k).or_default().0 += 1,
            None => unlabelled += 1,
        }
    }
    for r in host {
        match run_of(&r.name) {
            Some(k) => runs.entry(k).or_default().1 += 1,
            None => unlabelled += 1,
        }
    }
    let n_runs = runs.len();
    if n_runs < 2 || sub.is_empty() || unlabelled > 0 {
        return (1.0, String::new(), 0.0, n_runs);
    }
    let (top, (a, c)) = runs.iter().max_by_key(|(k, (s, _))| (*s, std::cmp::Reverse(**k))).map(|(k, v)| (*k, *v)).unwrap();
    let b = sub.len() as u64 - a;
    let d = host.len() as u64 - c;
    (fisher_exact(a, b, c, d), top.to_string(), a as f64 / sub.len() as f64, n_runs)
}

/// One PAF hit of a consensus against a genome.
#[derive(Clone, Debug)]
pub struct Hit {
    pub query: String,
    pub qcov: f64,
    pub identity: f64,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

/// Parse `minimap2 -c` PAF: identity = matches / alignment block length (cols 10/11), coverage = aligned
/// query span / query length.
pub fn parse_paf_hits(paf: &str) -> Vec<Hit> {
    let mut out = Vec::new();
    for line in paf.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let qlen: f64 = f[1].parse().unwrap_or(0.0);
        let (qs, qe): (f64, f64) = (f[2].parse().unwrap_or(0.0), f[3].parse().unwrap_or(0.0));
        let matches: f64 = f[9].parse().unwrap_or(0.0);
        let blen: f64 = f[10].parse().unwrap_or(0.0);
        if qlen <= 0.0 || blen <= 0.0 {
            continue;
        }
        out.push(Hit {
            query: f[0].to_string(),
            qcov: (qe - qs) / qlen,
            identity: matches / blen,
            chrom: f[5].to_string(),
            start: f[7].parse().unwrap_or(0),
            end: f[8].parse().unwrap_or(0),
        });
    }
    out
}

/// Best hit at the locus (overlapping `chrom:[start,end)`) and best hit elsewhere, among hits with
/// coverage >= `min_cov`.
pub fn home(hits: &[Hit], chrom: &str, start: u64, end: u64, min_cov: f64) -> (Option<Hit>, Option<Hit>) {
    let mut at: Option<Hit> = None;
    let mut other: Option<Hit> = None;
    for h in hits.iter().filter(|h| h.qcov >= min_cov) {
        let overlaps = h.chrom == chrom && h.start < end && h.end > start;
        let slot = if overlaps { &mut at } else { &mut other };
        if slot.as_ref().map_or(true, |b| h.identity > b.identity) {
            *slot = Some(h.clone());
        }
    }
    (at, other)
}

/// The pre-registered verdict order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Verdict {
    Contamination,
    ForeignSpecies,
    Hypermutation,
    RnaEditing,
    Scattered,
    UnannotatedParalogue,
    ReferenceAbsentCandidate,
}

impl Verdict {
    pub fn as_str(&self) -> &'static str {
        match self {
            Verdict::Contamination => "contamination",
            Verdict::ForeignSpecies => "foreign_species",
            Verdict::Hypermutation => "hypermutation",
            Verdict::RnaEditing => "rna_editing",
            Verdict::Scattered => "scattered",
            Verdict::UnannotatedParalogue => "unannotated_paralogue",
            Verdict::ReferenceAbsentCandidate => "reference_absent_candidate",
        }
    }
}

pub struct VerdictInput {
    pub run_p: f64,
    pub run_top_frac: f64,
    pub is_ig_tr: bool,
    pub n_psv: usize,
    pub shared_frac: f64,
    pub editing_frac: f64,
    pub host_identity: Option<f64>,
    pub other_identity: Option<f64>,
    /// best identity of the consensus in a FOREIGN genome (`--foreign`, e.g. human for a gorilla library)
    pub foreign_identity: Option<f64>,
    pub delta: f64,
    /// the locus fired on the exon-order rearrangement detector only (no divergence mixture): the PSV
    /// consistency and editing rules do not apply (addendum 2)
    pub structural_only: bool,
}

pub fn verdict(v: &VerdictInput) -> Verdict {
    if v.run_p < 1e-3 && v.run_top_frac >= 0.95 {
        return Verdict::Contamination;
    }
    // addendum 1: the consensus is a near-perfect match in another species' genome while `delta`-far from
    // its own host — cross-species contamination (the same shape as the confirmation rule, other genome)
    if confirmed(v.foreign_identity, v.host_identity, v.delta) && v.foreign_identity.map_or(false, |f| f >= 0.995) {
        return Verdict::ForeignSpecies;
    }
    if v.is_ig_tr {
        return Verdict::Hypermutation;
    }
    if !v.structural_only {
        if v.n_psv >= 5 && v.editing_frac >= 0.8 {
            return Verdict::RnaEditing;
        }
        if v.n_psv < 3 || v.shared_frac < 0.5 {
            return Verdict::Scattered;
        }
    }
    if let (Some(h), Some(o)) = (v.host_identity, v.other_identity) {
        if o > h {
            return Verdict::UnannotatedParalogue;
        }
    }
    Verdict::ReferenceAbsentCandidate
}

/// Confirmation rule: a near-perfect home in the confirm genome, `delta/2` closer than the primary host.
pub fn confirmed(conf_identity: Option<f64>, host_identity: Option<f64>, delta: f64) -> bool {
    match (conf_identity, host_identity) {
        (Some(c), Some(h)) => c >= 0.99 && c - h >= delta / 2.0,
        _ => false,
    }
}

/// Immunoglobulin / T-cell-receptor test on a gene name + description.
pub fn is_ig_tr(name: &str, description: &str) -> bool {
    let n = name.to_ascii_uppercase();
    let d = description.to_ascii_lowercase();
    let prefix = ["IGH", "IGK", "IGL", "TRA", "TRB", "TRD", "TRG"];
    let by_name = prefix.iter().any(|p| n.starts_with(p) && n.len() > 3 && matches!(n.as_bytes()[3], b'V' | b'D' | b'J' | b'C'));
    by_name || d.contains("immunoglobulin") || d.contains("t cell receptor") || d.contains("t-cell receptor")
}

/// Loci from a GTF (`gene_id` spans), a GFF (`gene`/`pseudogene` features; name from `Name=`/`gene=`),
/// or a BED. Returns `(id, chrom, start0, end, name)` in file order.
pub fn load_loci(path: &str) -> anyhow::Result<Vec<(String, String, u64, u64, String)>> {
    use std::io::BufRead;
    let f = std::fs::File::open(path).map_err(|e| anyhow::anyhow!("opening {path}: {e}"))?;
    let is_bed = path.ends_with(".bed");
    let is_gff = path.ends_with(".gff") || path.ends_with(".gff3") || path.ends_with(".gff.gz");
    let mut order: Vec<String> = Vec::new();
    let mut spans: DetHashMap<String, (String, u64, u64, String)> = DetHashMap::default();
    for line in std::io::BufReader::new(f).lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let fs: Vec<&str> = line.split('\t').collect();
        if is_bed {
            if fs.len() < 3 {
                continue;
            }
            let id = fs.get(3).map(|s| s.to_string()).unwrap_or_else(|| format!("{}:{}-{}", fs[0], fs[1], fs[2]));
            let (s, e) = (fs[1].parse::<u64>()?, fs[2].parse::<u64>()?);
            if spans.insert(id.clone(), (fs[0].to_string(), s, e, id.clone())).is_none() {
                order.push(id);
            }
            continue;
        }
        if fs.len() < 9 {
            continue;
        }
        let (s, e) = match (fs[3].parse::<u64>(), fs[4].parse::<u64>()) {
            (Ok(s), Ok(e)) => (s.saturating_sub(1), e),
            _ => continue,
        };
        if is_gff {
            if fs[2] != "gene" && fs[2] != "pseudogene" {
                continue;
            }
            let attrs: DetHashMap<&str, &str> = fs[8].split(';').filter_map(|kv| kv.split_once('=')).collect();
            let id = attrs.get("ID").map(|s| s.to_string()).unwrap_or_else(|| format!("{}:{}-{}", fs[0], s, e));
            let name = attrs.get("Name").or(attrs.get("gene")).map(|s| s.to_string()).unwrap_or_else(|| id.clone());
            if spans.insert(id.clone(), (fs[0].to_string(), s, e, name)).is_none() {
                order.push(id);
            }
        } else {
            if fs[2] != "exon" && fs[2] != "transcript" {
                continue;
            }
            let Some(g) = gtf_attr(fs[8], "gene_id") else { continue };
            if g.is_empty() {
                continue;
            }
            let key = format!("{}\t{}", fs[0], g);
            match spans.get_mut(&key) {
                Some(v) => {
                    v.1 = v.1.min(s);
                    v.2 = v.2.max(e);
                }
                None => {
                    spans.insert(key.clone(), (fs[0].to_string(), s, e, g.to_string()));
                    order.push(key);
                }
            }
        }
    }
    Ok(order
        .into_iter()
        .map(|k| {
            let (c, s, e, n) = spans.remove(&k).unwrap();
            let id = if k.contains('\t') { n.clone() } else { k };
            (id, c, s, e, n)
        })
        .collect())
}

fn gtf_attr<'a>(s: &'a str, key: &str) -> Option<&'a str> {
    let pat = format!("{key} \"");
    let i = s.find(&pat)? + pat.len();
    let j = s[i..].find('"')? + i;
    Some(&s[i..j])
}

/// Annotation genes for the hypermutation screen: per chrom, `(start0, end, is_ig_tr)`.
pub fn load_ig_tr(gff: &str) -> anyhow::Result<DetHashMap<String, Vec<(u64, u64)>>> {
    use std::io::BufRead;
    let f = std::fs::File::open(gff).map_err(|e| anyhow::anyhow!("opening {gff}: {e}"))?;
    let mut out: DetHashMap<String, Vec<(u64, u64)>> = DetHashMap::default();
    for line in std::io::BufReader::new(f).lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fs: Vec<&str> = line.split('\t').collect();
        if fs.len() < 9 || (fs[2] != "gene" && fs[2] != "pseudogene") {
            continue;
        }
        let attrs: DetHashMap<&str, &str> = fs[8].split(';').filter_map(|kv| kv.split_once('=')).collect();
        let name = attrs.get("Name").or(attrs.get("gene")).copied().unwrap_or("");
        let desc = attrs.get("description").copied().unwrap_or("");
        if is_ig_tr(name, desc) {
            if let (Ok(s), Ok(e)) = (fs[3].parse::<u64>(), fs[4].parse::<u64>()) {
                out.entry(fs[0].to_string()).or_default().push((s.saturating_sub(1), e));
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use crate::types::{DetHashMap, DetHashSet};
    use super::*;

    fn read(name: &str, de: f64, start: u64, cigar: &str, seq: &str) -> PileRead {
        let mut ops = Vec::new();
        let mut n = 0u64;
        for ch in cigar.chars() {
            if ch.is_ascii_digit() {
                n = n * 10 + (ch as u64 - '0' as u64);
            } else {
                ops.push((ch, n));
                n = 0;
            }
        }
        PileRead { name: name.into(), de, ref_start: start, ops, seq: seq.as_bytes().to_vec() }
    }

    #[test]
    fn two_means_matches_the_python_reference() {
        // detector.py::two_means on a 60 host + 20 sub-pile vector (seeded in Python): (0.25, 0.0285, 0.02645)
        let v = [
            0.0018, 0.0024, 0.0018, 0.0017, 0.0013, 0.0018, 0.0029, 0.0023, 0.0028, 0.0022, 0.0023, 0.0021, 0.0007, 0.0027, 0.0024,
            0.0024, 0.0006, 0.0006, 0.0013, 0.0016, 0.0022, 0.002, 0.0024, 0.0015, 0.0022, 0.0023, 0.0015, 0.0034, 0.0024, 0.003,
            0.0015, 0.0014, 0.0017, 0.0019, 0.0025, 0.0022, 0.0016, 0.0012, 0.0016, 0.003, 0.0014, 0.0022, 0.0023, 0.0008, 0.002,
            0.003, 0.0004, 0.0017, 0.0019, 0.0013, 0.0024, 0.002, 0.0008, 0.0027, 0.0025, 0.0028, 0.0032, 0.0023, 0.0021, 0.001,
            0.0318, 0.0282, 0.0286, 0.0262, 0.0271, 0.0284, 0.0339, 0.0239, 0.0256, 0.0307, 0.0343, 0.0317, 0.0243, 0.0224, 0.0311,
            0.0278, 0.0266, 0.0329, 0.0333, 0.0305,
        ];
        let (m, dh, delta, _) = two_means(&v).unwrap();
        assert!((m - 0.25).abs() < 1e-12);
        assert!((dh - 0.0285).abs() < 1e-9);
        assert!((delta - 0.02645).abs() < 1e-9);
        // host only: (0.55, 0.0024, 0.0009) — fires nothing because m > 0.5
        let (m2, dh2, d2, _) = two_means(&v[..60]).unwrap();
        assert!((m2 - 0.55).abs() < 1e-12 && (dh2 - 0.0024).abs() < 1e-9 && (d2 - 0.0009).abs() < 1e-9);
        // constant vector
        assert_eq!(two_means(&[0.001; 12]).unwrap().0, 0.0);
        assert!(two_means(&[0.001; 9]).is_none());
    }

    #[test]
    fn split_pile_fires_on_the_mixture_and_not_on_the_host_alone() {
        let mut reads: Vec<PileRead> = (0..30).map(|i| read(&format!("h{i}"), 0.002 + 0.0001 * (i % 5) as f64, 0, "10=", "ACGTACGTAC")).collect();
        assert!(split_pile(&reads, 0.10, 0.01, 3).is_none());
        reads.extend((0..8).map(|i| read(&format!("s{i}"), 0.03 + 0.001 * i as f64, 0, "10=", "ACGTACGTAC")));
        let s = split_pile(&reads, 0.10, 0.01, 3).unwrap();
        assert_eq!(s.sub.len(), 8);
        assert_eq!(s.host.len(), 30);
        assert!(s.delta > 0.02);
        // a sub-pile of 2 does not fire (min_sub)
        let small: Vec<PileRead> = reads[..32].to_vec();
        assert!(split_pile(&small, 0.10, 0.01, 3).is_none());
    }

    #[test]
    fn psv_sites_require_sharing_across_the_subpile_and_absence_in_the_host() {
        // reference 20 bp; sub-pile reads all mismatch at pos 5 (ref A -> G) and pos 12 (ref C -> T);
        // one sub read also has a private error at pos 17; host reads are clean except one at pos 12.
        let r = b"ACGTAACGTACGCTAGCTAG";
        let sub: Vec<PileRead> = (0..4)
            .map(|i| {
                let mut s = r.to_vec();
                s[5] = b'G';
                s[12] = b'T';
                let cigar = if i == 0 {
                    s[17] = b'A';
                    "5=1X6=1X4=1X2="
                } else {
                    "5=1X6=1X7="
                };
                read(&format!("s{i}"), 0.03, 0, cigar, std::str::from_utf8(&s).unwrap())
            })
            .collect();
        let mut host: Vec<PileRead> = (0..6).map(|i| read(&format!("h{i}"), 0.002, 0, "20=", std::str::from_utf8(r).unwrap())).collect();
        let mut h6 = r.to_vec();
        h6[12] = b'T';
        host.push(read("h6", 0.003, 0, "12=1X7=", std::str::from_utf8(&h6).unwrap()));
        let subr: Vec<&PileRead> = sub.iter().collect();
        let hostr: Vec<&PileRead> = host.iter().collect();
        let c = consistency(&subr, &hostr, r, 0);
        assert_eq!(c.sites.len(), 2);
        assert_eq!((c.sites[0].pos, c.sites[0].base), (5, b'G'));
        assert_eq!((c.sites[1].pos, c.sites[1].base), (12, b'T'));
        // 8 of 9 sub-pile mismatches fall on the two sites
        assert!((c.shared_frac - 8.0 / 9.0).abs() < 1e-12);
        // one of the two PSVs is A>G (editing-type)
        assert!((c.editing_frac - 0.5).abs() < 1e-12);
        // the consensus is the reference patched at both sites, over the template's blocks
        let t = template_read(&subr).unwrap();
        let (cons, blocks) = patched_consensus(t, &c.sites, r, 0);
        let mut want = r.to_vec();
        want[5] = b'G';
        want[12] = b'T';
        assert_eq!(cons, want);
        assert_eq!(blocks, vec![(0, 20)]);
    }

    #[test]
    fn spliced_template_yields_concatenated_blocks() {
        let r = b"ACGTAACGTACGCTAGCTAGGGGGGTTTTT";
        let t = read("t", 0.02, 2, "5=10N8=", "GTAAC" /* unused past block */);
        let (cons, blocks) = patched_consensus(&t, &[], r, 0);
        assert_eq!(blocks, vec![(2, 7), (17, 25)]);
        assert_eq!(cons, b"GTAACTAGGGGGG".to_vec());
    }

    #[test]
    fn fisher_and_run_screen() {
        // 10/10 sub-pile on run A, host 5 A / 95 B -> tiny p
        let p = fisher_exact(10, 0, 5, 95);
        assert!(p < 1e-6, "{p}");
        assert!((fisher_exact(5, 5, 5, 5) - 1.0).abs() < 1e-9);
        let sub: Vec<PileRead> = (0..10).map(|i| read(&format!("SRR100.{i}"), 0.03, 0, "5=", "ACGTA")).collect();
        let host: Vec<PileRead> = (0..100).map(|i| read(&format!("SRR{}.{i}", if i < 5 { 100 } else { 200 }), 0.002, 0, "5=", "ACGTA")).collect();
        let (p, top, frac, n) = run_screen(&sub.iter().collect::<Vec<_>>(), &host.iter().collect::<Vec<_>>());
        assert_eq!((top.as_str(), n), ("SRR100", 2));
        assert!(frac == 1.0 && p < 1e-3);
        // single-run BAM: screen inert
        let one: Vec<PileRead> = (0..10).map(|i| read(&format!("m64076/{i}/ccs"), 0.03, 0, "5=", "ACGTA")).collect();
        let (p1, _, _, n1) = run_screen(&one.iter().collect::<Vec<_>>(), &one.iter().collect::<Vec<_>>());
        assert_eq!((p1, n1), (1.0, 1));
        // names without run structure (a simulation) never make a "run": screen inert
        assert_eq!(run_of("extra|gene-NTSR1|rna-NM_002531.3|0"), None);
        assert_eq!(run_of("SRR27438212.581327"), Some("SRR27438212"));
        assert_eq!(run_of("m64076_221110_210557/78250013/ccs"), Some("m64076_221110_210557"));
        let sim: Vec<PileRead> = (0..10).map(|i| read(&format!("extra|g|t|{i}"), 0.03, 0, "5=", "ACGTA")).collect();
        let simh: Vec<PileRead> = (0..30).map(|i| read(&format!("template|g|t|{i}"), 0.002, 0, "5=", "ACGTA")).collect();
        let (ps, _, _, ns) = run_screen(&sim.iter().collect::<Vec<_>>(), &simh.iter().collect::<Vec<_>>());
        assert_eq!((ps, ns), (1.0, 0));
    }

    #[test]
    fn verdict_order_and_confirmation() {
        let base = VerdictInput { run_p: 1.0, run_top_frac: 0.0, is_ig_tr: false, n_psv: 6, shared_frac: 0.9, editing_frac: 0.0, host_identity: Some(0.98), other_identity: Some(0.95), foreign_identity: None, delta: 0.02, structural_only: false };
        assert_eq!(verdict(&base), Verdict::ReferenceAbsentCandidate);
        assert_eq!(verdict(&VerdictInput { n_psv: 0, shared_frac: 0.0, structural_only: true, ..base }), Verdict::ReferenceAbsentCandidate);
        assert_eq!(verdict(&VerdictInput { n_psv: 0, structural_only: true, other_identity: Some(0.999), ..base }), Verdict::UnannotatedParalogue);
        assert_eq!(verdict(&VerdictInput { foreign_identity: Some(0.999), ..base }), Verdict::ForeignSpecies);
        assert_eq!(verdict(&VerdictInput { foreign_identity: Some(0.985), ..base }), Verdict::ReferenceAbsentCandidate); // not near-perfect
        assert_eq!(verdict(&VerdictInput { foreign_identity: Some(0.999), is_ig_tr: true, ..base }), Verdict::ForeignSpecies); // before the IG screen
        assert_eq!(verdict(&VerdictInput { other_identity: Some(0.995), ..base }), Verdict::UnannotatedParalogue);
        assert_eq!(verdict(&VerdictInput { n_psv: 2, ..base }), Verdict::Scattered);
        assert_eq!(verdict(&VerdictInput { shared_frac: 0.3, ..base }), Verdict::Scattered);
        assert_eq!(verdict(&VerdictInput { editing_frac: 0.9, ..base }), Verdict::RnaEditing);
        assert_eq!(verdict(&VerdictInput { is_ig_tr: true, editing_frac: 0.9, ..base }), Verdict::Hypermutation);
        assert_eq!(verdict(&VerdictInput { run_p: 1e-5, run_top_frac: 1.0, is_ig_tr: true, ..base }), Verdict::Contamination);
        assert!(confirmed(Some(0.998), Some(0.98), 0.02));
        assert!(!confirmed(Some(0.998), Some(0.995), 0.02)); // no margin: the consensus was barely patched
        assert!(!confirmed(Some(0.985), Some(0.96), 0.02)); // not a near-perfect home
        assert!(!confirmed(None, Some(0.98), 0.02));
    }

    #[test]
    fn rearranged_exon_is_found_as_an_insertion_matching_a_skipped_exon() {
        // reference: E1 (0..60) intron (60..160) E2 (160..220) intron (220..320) E3 (320..380) intron E4 (480..540)
        let mut rf = vec![b'A'; 540];
        let gen = |seed: u64| -> Vec<u8> {
            let mut x = seed;
            (0..60).map(|_| { x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); b"ACGT"[((x >> 33) % 4) as usize] }).collect()
        };
        let (e1, e2, e3, e4) = (gen(1), gen(2), gen(3), gen(4));
        rf[0..60].copy_from_slice(&e1); rf[160..220].copy_from_slice(&e2); rf[320..380].copy_from_slice(&e3); rf[480..540].copy_from_slice(&e4);
        // read with exons 2 and 3 swapped: E1 E3 E2 E4. minimap2-style: E1 aligned, N to E3, E3 aligned, E2 inserted, N to E4.
        let mut seq = Vec::new(); seq.extend(&e1); seq.extend(&e3); seq.extend(&e2); seq.extend(&e4);
        let read = PileRead { name: "r1".into(), de: 0.002, ref_start: 0, ops: vec![('=', 60), ('N', 260), ('=', 60), ('I', 60), ('N', 100), ('=', 60)], seq };
        let rs = find_rearrangements(&[&read], &rf, 0, 50);
        assert_eq!(rs.len(), 1);
        assert_eq!(rs[0].exon, (160, 220));
        assert!(rs[0].identity > 0.99 && !rs[0].duplicated);
        assert_eq!(rs[0].ins_ref_pos, 380);
        // a read that carries E2 twice (in place AND inserted) is a duplication, not a rearrangement
        let mut seq2 = Vec::new(); seq2.extend(&e1); seq2.extend(&e2); seq2.extend(&e2); seq2.extend(&e3);
        let read2 = PileRead { name: "r2".into(), de: 0.002, ref_start: 0, ops: vec![('=', 60), ('N', 100), ('=', 60), ('I', 60), ('N', 100), ('=', 60)], seq: seq2 };
        let rs2 = find_rearrangements(&[&read2], &rf, 0, 50);
        assert_eq!(rs2.len(), 1);
        assert!(rs2[0].duplicated);
        // an insertion of random sequence matches no gap
        let read3 = PileRead { name: "r3".into(), de: 0.002, ref_start: 0, ops: vec![('=', 60), ('I', 60), ('N', 100), ('=', 60)], seq: [e1.clone(), vec![b'C'; 60], e2.clone()].concat() };
        assert!(find_rearrangements(&[&read3], &rf, 0, 50).is_empty());
        // clustering: three reads with the same rearrangement fire, a lone one does not
        let three: Vec<Rearrangement> = (0..3).map(|i| Rearrangement { name: format!("x{i}"), ins_ref_pos: 380 + i, ins_len: 60, exon: (160 + i, 220 + i), identity: 0.98, duplicated: false }).collect();
        let mut all = three.clone();
        all.push(Rearrangement { name: "lone".into(), ins_ref_pos: 100, ins_len: 70, exon: (900, 970), identity: 0.95, duplicated: false });
        let cl = rearrangement_clusters(&all, 20, 3);
        assert_eq!(cl.len(), 1);
        assert_eq!(cl[0].len(), 3);
    }

    #[test]
    fn ig_tr_names_and_descriptions() {
        assert!(is_ig_tr("IGLV3-1", ""));
        assert!(is_ig_tr("TRBC2", ""));
        assert!(is_ig_tr("LOC1", "immunoglobulin lambda variable 3-1-like"));
        assert!(!is_ig_tr("IGF1", "insulin like growth factor 1"));
        assert!(!is_ig_tr("TRAF3", "TNF receptor associated factor 3"));
    }

    #[test]
    fn paf_hits_and_home() {
        let paf = "c1\t1000\t0\t1000\t+\tchr1\t5000000\t100000\t101000\t980\t1000\t60\n\
                   c1\t1000\t0\t950\t+\tchr2\t5000000\t7000\t7950\t940\t950\t0\n\
                   c1\t1000\t0\t400\t+\tchr3\t5000000\t1\t401\t400\t400\t0\n";
        let hits = parse_paf_hits(paf);
        assert_eq!(hits.len(), 3);
        let (at, other) = home(&hits, "chr1", 100500, 100600, 0.8);
        assert!((at.unwrap().identity - 0.98).abs() < 1e-12);
        let o = other.unwrap();
        assert_eq!(o.chrom, "chr2"); // the chr3 hit fails coverage
        assert!((o.identity - 940.0 / 950.0).abs() < 1e-12);
    }

    #[test]
    fn loci_from_gtf_gff_bed() {
        let dir = std::env::temp_dir().join(format!("o3rna_loci_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let gtf = dir.join("a.gtf");
        std::fs::write(&gtf, "c1\tx\ttranscript\t101\t200\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\nc1\tx\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\nc1\tx\texon\t301\t400\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";\nc2\tx\texon\t1\t10\t.\t-\t.\tgene_id \"G2\"; transcript_id \"T3\";\n").unwrap();
        let l = load_loci(gtf.to_str().unwrap()).unwrap();
        assert_eq!(l, vec![("G1".into(), "c1".into(), 100, 400, "G1".into()), ("G2".into(), "c2".into(), 0, 10, "G2".into())]);
        let gff = dir.join("b.gff");
        std::fs::write(&gff, "c1\tR\tgene\t11\t20\t.\t+\t.\tID=gene-A;Name=A;description=x\nc1\tR\tmRNA\t11\t20\t.\t+\t.\tID=rna-A;Parent=gene-A\n").unwrap();
        assert_eq!(load_loci(gff.to_str().unwrap()).unwrap(), vec![("gene-A".into(), "c1".into(), 10, 20, "A".into())]);
        let bed = dir.join("c.bed");
        std::fs::write(&bed, "c1\t5\t9\tL1\n").unwrap();
        assert_eq!(load_loci(bed.to_str().unwrap()).unwrap(), vec![("L1".into(), "c1".into(), 5, 9, "L1".into())]);
        let _ = std::fs::remove_dir_all(&dir);
    }
}

// ---- merged 2026-10-05: was `vg_family/missing_copy_flag_pass.rs`, now the inline module below (one component) ----
#[allow(clippy::all)]
pub mod missing_copy_flag_pass {
//! O3 flag-pass detector (`docs/superpowers/specs/2026-09-10-o3-flag-pass-integration-design.md`):
//! ports `bench/missing_copy_flag_pass.py`'s missing-copy detector natively into `copy_assign`.
//!
//! **STATUS:** OPT-IN  (reachable via `copy_assign --flag-missing-copies`, src/bin/copy_assign.rs; default off)
    use crate::types::{DetHashMap, DetHashSet};


// `lgamma` lived in the retired ASJ module (dropped objective, tag notebook-2026-09-23c); kept here verbatim.
/// Lanczos approximation of `ln Γ(x)`.
pub(crate) fn lgamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_1,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        std::f64::consts::PI.ln() - (std::f64::consts::PI * x).sin().ln() - lgamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let t = x + G + 0.5;
        let mut a = C[0];
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    Divergent,
    Structural,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LocusClass {
    OtherFamily,
    AnnotatedNoUnit,
    Unannotated,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Flag {
    MissingCopy,
    Untestable,
    NoFlag,
}

#[derive(Clone, Debug)]
pub struct RawPair {
    pub family_id: String,
    pub copy_idx: String,
    pub is_partner: bool,
    pub n_rejected: usize,
    pub n_aligned: usize,
    pub covered_kb: f64,
    pub n_sites: usize,
    pub ctl_n: usize,
    pub ctl_covered_kb: f64,
    pub ctl_n_sites: usize,
    pub p_uncorrected: Option<f64>,
    pub class: Class,
    pub med_mismatch: f64,
    pub med_unaligned: f64,
}

#[derive(Clone, Debug)]
pub struct OrphanLocus {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_reads: usize,
    pub n_orphans: usize,
    pub class: LocusClass,
    pub n_genes_overlapping: usize,
    pub other_family_units: Vec<String>,
}

#[derive(Clone, Debug)]
pub struct FlaggedPair {
    pub pair: RawPair,
    pub flag: Flag,
    pub p_corrected_threshold: f64,
}

#[derive(Clone, Debug)]
pub struct O3Params {
    pub alpha: f64,
    pub max_reads: usize,
    pub min_reads: usize,
}

impl Default for O3Params {
    fn default() -> Self {
        O3Params { alpha: 0.001, max_reads: 500, min_reads: 3 }
    }
}

/// P(X >= k) for X ~ Poisson(lam). Mirrors `bench/missing_copy_flag_pass.py`'s `poisson_tail` exactly (same
/// closed-form sum, same edge cases): k==0 -> 1.0 always true; lam<=0 with k>=1 -> 0.0 (a Poisson(0) is
/// a point mass at 0).
pub fn poisson_tail(k: usize, lam: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if lam <= 0.0 {
        return 0.0;
    }
    let mut s = 0.0_f64;
    for i in 0..k {
        s += (-lam + (i as f64) * lam.ln() - lgamma((i + 1) as f64)).exp();
    }
    (1.0 - s).clamp(0.0, 1.0)
}

/// Genome-wide Bonferroni aggregation (Phase 2): n_pairs = every RawPair with a computed p-value, from
/// the WHOLE run, not just one family. threshold = alpha / max(n_pairs, 1). Pure, no I/O.
pub fn finalize_flags(all_pairs: &[RawPair], alpha: f64) -> Vec<FlaggedPair> {
    let n_pairs = all_pairs.iter().filter(|p| p.p_uncorrected.is_some()).count();
    let threshold = alpha / (n_pairs.max(1) as f64);
    all_pairs
        .iter()
        .map(|p| {
            let flag = match p.p_uncorrected {
                Some(pv) if pv < threshold => Flag::MissingCopy,
                Some(_) => Flag::NoFlag,
                None => Flag::Untestable,
            };
            FlaggedPair { pair: p.clone(), flag, p_corrected_threshold: threshold }
        })
        .collect()
}

pub(crate) struct AlignmentSummary {
    pub covered_kb: f64,
    pub n_sites: usize,
    pub per_read: DetHashMap<String, (usize, i64, usize)>,
}

/// Parses `minimap2 -x splice -c --eqx -N 1` PAF output. Target-position coverage/mismatch tallies (no
/// query-sequence lookup needed — see the module doc comment on why allele identity is dropped).
/// PAF columns used: [0]=query name [1]=query len [2]=query start [3]=query end [7]=target start,
/// [12..]=tags (the `cg:Z:` CIGAR tag).
pub(crate) fn parse_paf_consistency(paf_text: &str) -> AlignmentSummary {
    let mut cov: DetHashMap<u64, u32> = DetHashMap::default();
    let mut mism: DetHashMap<u64, u32> = DetHashMap::default();
    let mut per_read: DetHashMap<String, (usize, i64, usize)> = DetHashMap::default();
    for line in paf_text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 13 {
            continue;
        }
        let cg = match f[12..].iter().find(|t| t.starts_with("cg:Z:")) {
            Some(t) => &t[5..],
            None => continue,
        };
        let qlen: i64 = f[1].parse().unwrap_or(0);
        let qstart: i64 = f[2].parse().unwrap_or(0);
        let qend: i64 = f[3].parse().unwrap_or(0);
        let mut t: u64 = f[7].parse().unwrap_or(0);
        let mut nx: usize = 0;
        for (num, op) in crate::family::shared_definition::cigar_ops(cg) {
            match op {
                '=' => {
                    for k in 0..num {
                        *cov.entry(t + k).or_insert(0) += 1;
                    }
                    t += num;
                }
                'X' => {
                    nx += num as usize;
                    for k in 0..num {
                        *cov.entry(t + k).or_insert(0) += 1;
                        *mism.entry(t + k).or_insert(0) += 1;
                    }
                    t += num;
                }
                'D' | 'N' => t += num,
                _ => {} // 'I'/'S' advance only the query position, which this function never tracks
            }
        }
        let unal = qlen - (qend - qstart);
        let aligned_len = (qend - qstart) as usize;
        let name = f[0].to_string();
        let better = per_read.get(&name).map_or(true, |&(_, _, prev)| aligned_len > prev);
        if better {
            per_read.insert(name, (nx, unal, aligned_len));
        }
    }
    let n_sites = mism
        .iter()
        .filter(|(p, &c)| c >= 3 && (c as f64) >= 0.5 * (*cov.get(p).unwrap_or(&0) as f64))
        .count();
    let covered_kb = cov.values().filter(|&&c| c >= 3).count() as f64 / 1000.0;
    AlignmentSummary { covered_kb, n_sites, per_read }
}

#[cfg(test)]
mod paf_tests {
    use super::*;

    #[test]
    fn single_clean_read_no_mismatches() {
        // read "r1", qlen 100, aligned 0..100 on target starting at 1000, 100 matched bases
        let paf = "r1\t100\t0\t100\t+\tY\t2000\t1000\t1100\t100\t100\t60\tcg:Z:100=";
        let s = parse_paf_consistency(paf);
        assert_eq!(s.n_sites, 0);
        assert_eq!(s.per_read.get("r1"), Some(&(0, 0, 100)));
        assert!((s.covered_kb - 0.0).abs() < 1e-9); // 100 positions >= coverage 3? NO - coverage is 1 here
    }

    #[test]
    fn covered_kb_requires_coverage_at_least_three() {
        // three reads all covering the same 10bp target window -> those 10 positions reach coverage 3
        let paf = concat!(
            "r1\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "r2\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "r3\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=",
        );
        let s = parse_paf_consistency(paf);
        assert!((s.covered_kb - 0.010).abs() < 1e-9, "10 positions at cov>=3 => 0.010 kb, got {}", s.covered_kb);
    }

    #[test]
    fn consistent_mismatch_site_needs_majority_and_floor_of_three() {
        // 4 reads: 3 mismatch at target pos 5, 1 matches -> count 3 >= 3 AND 3 >= 0.5*4 (cov=4) -> consistent
        let paf = concat!(
            "a\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "b\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "c\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "d\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=",
        );
        let s = parse_paf_consistency(paf);
        assert_eq!(s.n_sites, 1);
    }

    #[test]
    fn mismatch_below_majority_is_not_consistent() {
        // 4 reads, only 1 mismatches at pos 5 (cov 4, mismatch count 1 < 0.5*4=2) -> not consistent
        let paf = concat!(
            "a\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:5=1X4=\n",
            "b\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "c\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n",
            "d\t10\t0\t10\t+\tY\t100\t0\t10\t10\t10\t60\tcg:Z:10=",
        );
        let s = parse_paf_consistency(paf);
        assert_eq!(s.n_sites, 0);
    }

    #[test]
    fn keeps_the_longer_alignment_when_a_read_has_two_paf_lines() {
        let paf = concat!(
            "r1\t100\t0\t40\t+\tY\t1000\t0\t40\t40\t40\t60\tcg:Z:40=\n",
            "r1\t100\t0\t90\t+\tY\t1000\t100\t190\t90\t90\t60\tcg:Z:90=",
        );
        let s = parse_paf_consistency(paf);
        assert_eq!(s.per_read.get("r1").map(|&(_, _, len)| len), Some(90));
    }

    #[test]
    fn unaligned_bases_is_query_len_minus_aligned_span() {
        // qlen 100, aligned query 10..90 -> unal = 100 - (90-10) = 20
        let paf = "r1\t100\t10\t90\t+\tY\t1000\t0\t80\t80\t80\t60\tcg:Z:80=";
        let s = parse_paf_consistency(paf);
        assert_eq!(s.per_read.get("r1"), Some(&(0, 20, 80)));
    }
}

#[cfg(test)]
mod tests {
    use crate::types::{DetHashMap, DetHashSet};
    use super::*;

    fn pair(p: Option<f64>) -> RawPair {
        RawPair {
            family_id: "F".into(), copy_idx: "0".into(), is_partner: false,
            n_rejected: 5, n_aligned: 5, covered_kb: 1.0, n_sites: 1,
            ctl_n: 5, ctl_covered_kb: 1.0, ctl_n_sites: 1,
            p_uncorrected: p, class: Class::Divergent, med_mismatch: 1.0, med_unaligned: 0.0,
        }
    }

    #[test]
    fn poisson_tail_k_zero_is_always_one() {
        assert_eq!(poisson_tail(0, 5.0), 1.0);
        assert_eq!(poisson_tail(0, 0.0), 1.0);
    }

    #[test]
    fn poisson_tail_zero_lambda_with_positive_k_is_zero() {
        assert_eq!(poisson_tail(3, 0.0), 0.0);
    }

    #[test]
    fn poisson_tail_matches_hand_computed_value() {
        // P(X >= 1) for Poisson(1) = 1 - P(X=0) = 1 - e^-1 = 0.6321205588...
        let p = poisson_tail(1, 1.0);
        assert!((p - 0.6321205588).abs() < 1e-6, "got {p}");
        // P(X >= 5) for Poisson(1) is small; hand value from scipy.stats.poisson.sf(4, 1) = 0.003659846827...
        let p2 = poisson_tail(5, 1.0);
        assert!((p2 - 0.003659846827).abs() < 1e-6, "got {p2}");
    }

    #[test]
    fn finalize_flags_labels_below_threshold_as_missing_copy() {
        // 2 testable pairs -> threshold = 0.001 / 2 = 0.0005
        let pairs = vec![pair(Some(0.0001)), pair(Some(0.9))];
        let flagged = finalize_flags(&pairs, 0.001);
        assert_eq!(flagged[0].flag, Flag::MissingCopy);
        assert_eq!(flagged[1].flag, Flag::NoFlag);
        assert!((flagged[0].p_corrected_threshold - 0.0005).abs() < 1e-12);
    }

    #[test]
    fn finalize_flags_untestable_when_p_is_none() {
        let pairs = vec![pair(None), pair(Some(0.0001))];
        let flagged = finalize_flags(&pairs, 0.001);
        assert_eq!(flagged[0].flag, Flag::Untestable);
        // n_pairs counts only the testable one, so threshold = 0.001 / 1 = 0.001
        assert!((flagged[1].p_corrected_threshold - 0.001).abs() < 1e-12);
    }

    #[test]
    fn finalize_flags_on_empty_input_does_not_divide_by_zero() {
        let flagged = finalize_flags(&[], 0.001);
        assert!(flagged.is_empty());
    }
}

/// Orphan-locus classification (§ "Scope note" in the design doc): OtherFamily > AnnotatedNoUnit >
/// Unannotated precedence. `all_units_by_chrom`/`genes_by_chrom` built ONCE per run, shared read-only.
/// `pub` (not `pub(crate)`): Task 5 (`src/bin/copy_assign.rs`) is a separate binary crate and calls this
/// directly per the integration plan's "Consumes" list -- `pub(crate)` here was invisible to it (E0603).
/// Pure visibility widening, no signature/behavior change; every other item Task 5 needs was already `pub`.
///
/// `own_family_ids` (final whole-branch-review fix round, follow-up): a SET, not a single id --
/// `bench/missing_copy_flag_pass.py`'s own exclusion (`u[2] != cp[next(iter(cp))]['family_id'].split('_')[0]`) is
/// trivially always correct because each Python sweep is ISOLATED to exactly one catalog family. The
/// caller here is a locally co-located physical sweep group that can legitimately bundle MORE THAN ONE
/// true catalog family (the same fact Fix 1 addressed for the pair detector) -- passing a single id (the
/// group's own, possibly arbitrary, local label) would incorrectly classify an orphan cluster's OWN
/// bundled sibling family's real unit as "OtherFamily" whenever that sibling isn't the one label chosen to
/// represent the group. `own_family_ids` is every TRUE catalog family id actually present among the
/// calling group's own copies; a unit belongs to "this same group" (excluded from `OtherFamily`) iff its
/// `fid` is a member of this set.
pub fn classify_orphan_locus(
    chrom: &str,
    start: u64,
    end: u64,
    own_family_ids: &DetHashSet<String>,
    all_units_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64, String, String)>>,
    genes_by_chrom: &std::collections::BTreeMap<String, Vec<(u64, u64)>>,
) -> (LocusClass, usize, Vec<String>) {
    let overlaps = |a0: u64, a1: u64, b0: u64, b1: u64| a0 < b1 && b0 < a1;
    let other_family_units: Vec<String> = all_units_by_chrom
        .get(chrom)
        .map(|v| {
            v.iter()
                .filter(|(s, e, fid, _)| overlaps(start, end, *s, *e) && !own_family_ids.contains(fid))
                .take(3)
                .map(|(_, _, fid, cidx)| format!("{fid}:{cidx}"))
                .collect()
        })
        .unwrap_or_default();
    let n_genes_overlapping = genes_by_chrom
        .get(chrom)
        .map(|v| v.iter().filter(|(s, e)| overlaps(start, end, *s, *e)).count())
        .unwrap_or(0);
    let class = if !other_family_units.is_empty() {
        LocusClass::OtherFamily
    } else if n_genes_overlapping > 0 {
        LocusClass::AnnotatedNoUnit
    } else {
        LocusClass::Unannotated
    };
    (class, n_genes_overlapping, other_family_units)
}

#[cfg(test)]
mod locus_tests {
    use super::*;
    use std::collections::BTreeMap;

    fn units() -> BTreeMap<String, Vec<(u64, u64, String, String)>> {
        let mut m = BTreeMap::new();
        m.insert("chr1".to_string(), vec![(100, 200, "OTHERFAM".to_string(), "3".to_string())]);
        m
    }

    fn genes() -> BTreeMap<String, Vec<(u64, u64)>> {
        let mut m = BTreeMap::new();
        m.insert("chr1".to_string(), vec![(500, 600)]);
        m
    }

    fn myfam() -> DetHashSet<String> {
        ["MYFAM".to_string()].into_iter().collect()
    }

    #[test]
    fn overlapping_another_familys_unit_wins_other_family() {
        let (class, _, units_hit) = classify_orphan_locus("chr1", 150, 250, &myfam(), &units(), &genes());
        assert_eq!(class, LocusClass::OtherFamily);
        assert_eq!(units_hit, vec!["OTHERFAM:3".to_string()]);
    }

    #[test]
    fn own_family_unit_does_not_count_as_other_family() {
        let mut u = units();
        u.get_mut("chr1").unwrap().push((150, 250, "MYFAM".to_string(), "0".to_string()));
        let (class, _, units_hit) = classify_orphan_locus("chr1", 150, 250, &myfam(), &u, &genes());
        // still hits OTHERFAM's unit at 100-200 too, so still OtherFamily -- verifies own-family rows
        // are excluded, not that the whole overlap set is
        assert_eq!(class, LocusClass::OtherFamily);
        assert!(!units_hit.iter().any(|s| s.starts_with("MYFAM:")));
    }

    #[test]
    fn no_other_family_unit_but_gene_overlap_is_annotated_no_unit() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chr1", 550, 650, &myfam(), &units(), &genes());
        assert_eq!(class, LocusClass::AnnotatedNoUnit);
        assert_eq!(n_genes, 1);
        assert!(units_hit.is_empty());
    }

    #[test]
    fn nothing_overlapping_is_unannotated() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chr1", 9000, 9100, &myfam(), &units(), &genes());
        assert_eq!(class, LocusClass::Unannotated);
        assert_eq!(n_genes, 0);
        assert!(units_hit.is_empty());
    }

    #[test]
    fn other_family_units_capped_at_three() {
        let mut u = BTreeMap::new();
        u.insert(
            "chr1".to_string(),
            (0..5).map(|i| (100u64, 200u64, format!("FAM{i}"), "0".to_string())).collect(),
        );
        let (_, _, units_hit) = classify_orphan_locus("chr1", 100, 200, &myfam(), &u, &BTreeMap::new());
        assert_eq!(units_hit.len(), 3);
    }

    #[test]
    fn unknown_chrom_is_unannotated() {
        let (class, n_genes, units_hit) = classify_orphan_locus("chrZZZ", 0, 10, &myfam(), &units(), &genes());
        assert_eq!(class, LocusClass::Unannotated);
        assert_eq!(n_genes, 0);
        assert!(units_hit.is_empty());
    }

    #[test]
    fn own_family_ids_set_excludes_every_member_not_just_the_first() {
        // Final whole-branch-review fix round, follow-up: a caller that bundles TWO true catalog
        // families (e.g. a locally co-located `fa` spanning more than one catalog family, the same fact
        // Fix 1 addressed for the pair detector) must have BOTH of its own family ids excluded from
        // `OtherFamily`, not just whichever single label happened to be chosen to represent the group.
        let mut u = BTreeMap::new();
        u.insert(
            "chr1".to_string(),
            vec![
                (100, 200, "SIBLING_B".to_string(), "0".to_string()), // bundled INTO the same caller group
                (100, 200, "TRULY_OTHER".to_string(), "0".to_string()), // NOT part of the caller's group
            ],
        );
        let own = ["SIBLING_A".to_string(), "SIBLING_B".to_string()].into_iter().collect();
        let (class, _, units_hit) = classify_orphan_locus("chr1", 100, 200, &own, &u, &BTreeMap::new());
        assert_eq!(class, LocusClass::OtherFamily, "TRULY_OTHER's unit is real external evidence");
        assert_eq!(
            units_hit,
            vec!["TRULY_OTHER:0".to_string()],
            "SIBLING_B must be excluded as self (it is in own_family_ids) even though it is not the \
             single id a bare-string comparison would have used"
        );
    }

    #[test]
    fn own_family_ids_set_with_no_external_unit_is_not_other_family() {
        // The bundle's OWN two sibling families both overlap the locus; with no truly external unit
        // present, the class must fall through to AnnotatedNoUnit/Unannotated, not OtherFamily -- the
        // single-id bug (comparing against only one arbitrary label) would have wrongly reported
        // OtherFamily here for whichever sibling wasn't the chosen label.
        let mut u = BTreeMap::new();
        u.insert(
            "chr1".to_string(),
            vec![
                (100, 200, "SIBLING_A".to_string(), "0".to_string()),
                (100, 200, "SIBLING_B".to_string(), "0".to_string()),
            ],
        );
        let own = ["SIBLING_A".to_string(), "SIBLING_B".to_string()].into_iter().collect();
        let (class, _, units_hit) = classify_orphan_locus("chr1", 100, 200, &own, &u, &genes());
        assert_eq!(class, LocusClass::Unannotated);
        assert!(units_hit.is_empty());
    }
}

/// One candidate copy Y's test batch (its own origin-rejected reads) and control batch (its own
/// certificate-accepted reads), both as `(name, sequence)` pairs ready to realign.
pub struct PairInput {
    pub copy_idx: String,
    pub is_partner: bool,
    pub rejected: Vec<(String, Vec<u8>)>,
    pub accepted: Vec<(String, Vec<u8>)>,
}

/// Realigns one batch of reads to a target window via `minimap2 -x splice -c --eqx -N 1`, mirroring
/// `copy_assign_pipeline.rs`'s `minimap2_msa_pair` temp-file convention (pid+atomic-nonce names,
/// `RUSTLE_MINIMAP2` env override, `Drop`-based cleanup) so concurrent region-parallel workers never
/// collide on the same path.
fn realign_batch(target_seq: &[u8], reads: &[(String, Vec<u8>)]) -> anyhow::Result<AlignmentSummary> {
    use std::io::Write;
    if reads.is_empty() {
        return Ok(AlignmentSummary { covered_kb: 0.0, n_sites: 0, per_read: DetHashMap::default() });
    }
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    use std::sync::atomic::{AtomicUsize, Ordering};
    static O3_NONCE: AtomicUsize = AtomicUsize::new(0);
    let nonce = O3_NONCE.fetch_add(1, Ordering::Relaxed);
    let ypath = dir.join(format!("rustle_o3_y_{pid}_{nonce}.fa"));
    let rpath = dir.join(format!("rustle_o3_reads_{pid}_{nonce}.fa"));
    struct Cleanup(std::path::PathBuf, std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
            let _ = std::fs::remove_file(&self.1);
        }
    }
    let _cl = Cleanup(ypath.clone(), rpath.clone());
    {
        let mut y = std::fs::File::create(&ypath)?;
        y.write_all(b">Y\n")?;
        y.write_all(target_seq)?;
        y.write_all(b"\n")?;
        let mut r = std::fs::File::create(&rpath)?;
        for (name, seq) in reads {
            writeln!(r, ">{name}")?;
            r.write_all(seq)?;
            r.write_all(b"\n")?;
        }
    }
    // Minor (final whole-branch review): `-t 4` matches the Python reference's own hardcoded thread count,
    // but this codebase's existing `minimap2_msa_pair` pattern this function is modeled on uses `-t 1`.
    // Under `copy_assign --region-parallel N` this multiplies to up to 4xN minimap2 threads concurrently
    // -- on this project's 5-core WSL2 box (see the crash-rule notes), combining `--region-parallel` with
    // `--flag-missing-copies` could overload the machine. Left at `-t 4` (not changed here -- this is a
    // documentation-only note), but `--flag-missing-copies` runs on this machine should stay serial
    // (no `--region-parallel`) until/unless this is revisited.
    let out = std::process::Command::new(&mm2)
        .args(["-x", "splice", "-c", "--eqx", "-N", "1", "-t", "4"])
        .arg(&ypath)
        .arg(&rpath)
        .output()?;
    if !out.status.success() {
        anyhow::bail!("minimap2 failed: {}", String::from_utf8_lossy(&out.stderr));
    }
    Ok(parse_paf_consistency(&String::from_utf8_lossy(&out.stdout)))
}

fn median(mut v: Vec<f64>) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(v[v.len() / 2])
}

/// Realignment target window for one Y-copy, extracted as a pure function so Task 7's window-selection
/// fix (using the catalog's L2 locus extent when present, matching `bench/missing_copy_flag_pass.py`'s `detector()`
/// lines 45-48) is directly unit-testable without a real `minimap2`/`GenomeIndex`. `locus` is the
/// catalog's `(locus_start, locus_end)` for this copy when one is recorded; `longest_rejected` is the
/// length of the longest rejected read, used only in the `None` fallback to pad the bare copy span (so a
/// read realigned against the window cannot hang off its edge).
fn locus_or_padded_window(s: u64, e: u64, locus: Option<(u64, u64)>, longest_rejected: u64) -> (u64, u64) {
    match locus {
        Some((locus_s, locus_e)) => (locus_s.min(s), locus_e.max(e)),
        None => (s.saturating_sub(longest_rejected), e + longest_rejected),
    }
}

/// Phase 1: one family's raw (uncorrected) missing-copy pair statistics. Skips any Y with fewer than
/// `params.min_reads` rejected reads (matches the Python's `len(names) < 3: continue` -- not emitted as
/// Untestable, simply absent from the output). A `minimap2` failure for one Y is logged to stderr and
/// that pair is skipped -- see the design doc's Error Handling section for why this must not abort.
///
/// `copy_span_by_catalog_idx` maps to `(chrom, start, end, locus_extent)`: the realignment TARGET window
/// mirrors `bench/missing_copy_flag_pass.py`'s `detector()` exactly (lines 45-48) -- `(min(locus_start, start),
/// max(locus_end, end))` when the catalog carries an L2 locus extent for this copy, else the copy's own
/// span padded by the longest rejected read on each side (see [`locus_or_padded_window`], which does this
/// computation and is unit-tested directly). Using the bare copy span unconditionally (the pre-fix
/// behaviour) mis-sizes the window whenever a copy's locus extent differs from its own span -- found
/// during Task 7's reproduction gate: `covered_kb`/`n_sites`/`rate_per_kb`/`p` all move even when
/// `n_rejected` matches exactly (e.g. MCL117_073244 copy 1: same 4 rejected reads, Python rate 4.27/kb
/// vs the unfixed Rust 59.81/kb).
pub fn detect_missing_copy_pairs(
    family_id: &str,
    copy_span_by_catalog_idx: &DetHashMap<String, (String, u64, u64, Option<(u64, u64)>)>,
    genome: &crate::genome::GenomeIndex,
    inputs: &[PairInput],
    params: &O3Params,
) -> Vec<RawPair> {
    let mut out = Vec::new();
    for input in inputs {
        if input.rejected.len() < params.min_reads {
            continue;
        }
        let Some((chrom, s, e, locus)) = copy_span_by_catalog_idx.get(&input.copy_idx) else {
            continue;
        };
        // Task 7 fix (reproduction-gate finding): the Python caps each side at `max_reads` by NAME order
        // (`sorted(names)[:max_reads]` / `sorted(n for n in truth if ...)[:max_reads]`), not by collection
        // order -- irrelevant when a family stays under the cap (most pairs), but for a copy with more
        // than `max_reads` eligible reads (e.g. MCL121_073244 copy 0's 500-capped control pool) the two
        // orders sample a DIFFERENT subset, changing `covered_kb`/`n_sites`/`p` even when the uncapped
        // rejected side matches exactly. Sorting first makes the cap deterministic and Python-identical.
        let mut rejected_sorted = input.rejected.clone();
        rejected_sorted.sort_by(|a, b| a.0.cmp(&b.0));
        let rejected: Vec<_> = rejected_sorted.into_iter().take(params.max_reads).collect();
        let pad = rejected.iter().map(|(_, seq)| seq.len() as u64).max().unwrap_or(0);
        let (ls, le) = locus_or_padded_window(*s, *e, *locus, pad);
        let target = match genome.fetch_sequence(chrom, ls, le) {
            Some(t) => t,
            None => continue,
        };
        let test = match realign_batch(&target, &rejected) {
            Ok(s) => s,
            Err(err) => {
                eprintln!("[o3-flag-pass] realignment failed for {family_id}:{}: {err}", input.copy_idx);
                continue;
            }
        };
        let mut accepted_sorted = input.accepted.clone(); // same name-order cap as `rejected` above
        accepted_sorted.sort_by(|a, b| a.0.cmp(&b.0));
        let accepted: Vec<_> = accepted_sorted.into_iter().take(params.max_reads).collect();
        let ctl = if accepted.len() >= params.min_reads {
            match realign_batch(&target, &accepted) {
                Ok(s) => s,
                Err(err) => {
                    eprintln!("[o3-flag-pass] control realignment failed for {family_id}:{}: {err}", input.copy_idx);
                    AlignmentSummary { covered_kb: 0.0, n_sites: 0, per_read: DetHashMap::default() }
                }
            }
        } else {
            AlignmentSummary { covered_kb: 0.0, n_sites: 0, per_read: DetHashMap::default() }
        };
        // rate floored at 1 site over the control's covered kb, so a control with zero observed sites
        // never claims a zero rate (which would trivially pass every test) -- design doc, RawPair section.
        let p_uncorrected = if test.covered_kb > 0.0 && ctl.covered_kb > 0.0 {
            let rate = (ctl.n_sites.max(1) as f64) / ctl.covered_kb;
            Some(poisson_tail(test.n_sites, rate * test.covered_kb))
        } else {
            None
        };
        let mismatches: Vec<f64> = test.per_read.values().map(|&(nx, _, _)| nx as f64).collect();
        let unaligned: Vec<f64> = test.per_read.values().map(|&(_, unal, _)| unal as f64).collect();
        let med_mismatch = median(mismatches).unwrap_or(f64::NAN);
        let med_unaligned = median(unaligned).unwrap_or(f64::NAN);
        let class = if !med_mismatch.is_nan() && med_mismatch > med_unaligned {
            Class::Divergent
        } else {
            Class::Structural
        };
        out.push(RawPair {
            family_id: family_id.to_string(),
            copy_idx: input.copy_idx.clone(),
            is_partner: input.is_partner,
            n_rejected: input.rejected.len(),
            n_aligned: rejected.len(),
            covered_kb: test.covered_kb,
            n_sites: test.n_sites,
            ctl_n: accepted.len(),
            ctl_covered_kb: ctl.covered_kb,
            ctl_n_sites: ctl.n_sites,
            p_uncorrected,
            class,
            med_mismatch,
            med_unaligned,
        });
    }
    out
}

#[cfg(test)]
mod pair_detector_tests {
    use super::*;
    use crate::genome::GenomeIndex;

    #[test]
    fn fewer_than_min_reads_is_skipped_entirely() {
        let genome = GenomeIndex::from_seqs(&[("chrT", &[b'A'; 100])]);
        let mut spans = DetHashMap::default();
        spans.insert("0".to_string(), ("chrT".to_string(), 0u64, 100u64, None));
        let inputs = vec![PairInput {
            copy_idx: "0".to_string(),
            is_partner: false,
            rejected: vec![("r1".to_string(), vec![b'A'; 50])], // only 1, min_reads default is 3
            accepted: vec![],
        }];
        let pairs = detect_missing_copy_pairs("F", &spans, &genome, &inputs, &O3Params::default());
        assert!(pairs.is_empty());
    }

    #[test]
    fn missing_copy_span_is_skipped_not_errored() {
        let genome = GenomeIndex::from_seqs(&[("chrT", &[b'A'; 100])]);
        let spans = DetHashMap::default(); // no span for "0"
        let inputs = vec![PairInput {
            copy_idx: "0".to_string(),
            is_partner: false,
            rejected: vec![
                ("r1".to_string(), vec![b'A'; 50]),
                ("r2".to_string(), vec![b'A'; 50]),
                ("r3".to_string(), vec![b'A'; 50]),
            ],
            accepted: vec![],
        }];
        let pairs = detect_missing_copy_pairs("F", &spans, &genome, &inputs, &O3Params::default());
        assert!(pairs.is_empty());
    }

    #[test]
    fn locus_extent_wider_than_the_copy_span_wins() {
        // Locus extent (40,160) strictly contains the copy span (50,120): the window must widen to the
        // locus, matching `bench/missing_copy_flag_pass.py`'s `(min(locus_start, start), max(locus_end, end))`.
        assert_eq!(locus_or_padded_window(50, 120, Some((40, 160)), 999), (40, 160));
    }

    #[test]
    fn locus_extent_narrower_than_the_copy_span_does_not_shrink_it() {
        // Locus extent (60,110) sits INSIDE the copy span (50,120): min/max never shrinks below the
        // copy's own span, so the window collapses to the bare span, not the narrower locus.
        assert_eq!(locus_or_padded_window(50, 120, Some((60, 110)), 999), (50, 120));
    }

    #[test]
    fn no_locus_falls_back_to_padding_by_the_longest_rejected_read() {
        // No catalog locus extent recorded (`None`): pad the bare copy span by the longest rejected
        // read's length on each side.
        assert_eq!(locus_or_padded_window(1000, 2000, None, 250), (750, 2250));
        // saturating_sub must not underflow when the pad exceeds the start coordinate.
        assert_eq!(locus_or_padded_window(100, 200, None, 500), (0, 700));
    }

    #[test]
    fn min_reads_guard_below_threshold_returns_empty() {
        // Smoke test for detect_missing_copy_pairs when input is below the min_reads threshold.
        // Verifies the function returns early without crashing when given insufficient rejected reads.
        // The Some(locus) destructure and locus_or_padded_window call logic are covered separately
        // by the three locus_or_padded_window unit tests directly above.
        let genome = GenomeIndex::from_seqs(&[("chrT", &[b'A'; 300])]);
        let mut spans = DetHashMap::default();
        // locus extent (0,300) is wider than the bare copy span (100,200).
        spans.insert("0".to_string(), ("chrT".to_string(), 100u64, 200u64, Some((0u64, 300u64))));
        let inputs = vec![PairInput {
            copy_idx: "0".to_string(),
            is_partner: false,
            rejected: vec![("r1".to_string(), vec![b'A'; 50])], // 1 < default min_reads (3)
            accepted: vec![],
        }];
        let pairs = detect_missing_copy_pairs("F", &spans, &genome, &inputs, &O3Params::default());
        assert!(pairs.is_empty(), "below min_reads, skipped before realign_batch as usual");
    }
}
}

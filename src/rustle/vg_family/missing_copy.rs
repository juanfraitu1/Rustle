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

use std::collections::{BTreeMap, HashMap};

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
    cov: HashMap<u64, u32>,
    mism: HashMap<u64, u32>,
    bases: HashMap<u64, HashMap<u8, u32>>,
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
    let patch: HashMap<u64, u8> = sites.iter().map(|s| (s.pos, s.base)).collect();
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
    let mut index: HashMap<&[u8], Vec<usize>> = HashMap::new();
    for i in 0..=target.len() - K {
        index.entry(&target[i..i + K]).or_default().push(i);
    }
    let mut votes: HashMap<i64, u32> = HashMap::new();
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
    let mut spans: HashMap<String, (String, u64, u64, String)> = HashMap::new();
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
            let attrs: HashMap<&str, &str> = fs[8].split(';').filter_map(|kv| kv.split_once('=')).collect();
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
pub fn load_ig_tr(gff: &str) -> anyhow::Result<HashMap<String, Vec<(u64, u64)>>> {
    use std::io::BufRead;
    let f = std::fs::File::open(gff).map_err(|e| anyhow::anyhow!("opening {gff}: {e}"))?;
    let mut out: HashMap<String, Vec<(u64, u64)>> = HashMap::new();
    for line in std::io::BufReader::new(f).lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let fs: Vec<&str> = line.split('\t').collect();
        if fs.len() < 9 || (fs[2] != "gene" && fs[2] != "pseudogene") {
            continue;
        }
        let attrs: HashMap<&str, &str> = fs[8].split(';').filter_map(|kv| kv.split_once('=')).collect();
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

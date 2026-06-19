//! Allele-specific junctions (ASJ) — advisor interest #3. An ASJ is a splice junction whose USAGE depends
//! on the allele a molecule carries. Long HiFi reads observe BOTH the heterozygous-SNP allele AND the read's
//! introns on the SAME molecule, so allele→junction linkage is PER-MOLECULE — no statistical phasing (the
//! short-read sQTL problem is sidestepped). The het loci that confound copy-detection are exactly the
//! substrate (heterozygosity = the phasing signal: confound → feature).
//!
//! Port of `bench/allele_specific_junctions.py`, reusing `copy_split`'s per-read allele + intron extraction.
//! Per gene: call a balanced het ANCHOR SNP, partition reads by the allele they carry, and for each junction
//! compute PSI = used/spanning within each allele; a 2×2 Fisher exact + |ΔPSI| effect size gives the call.
//! Genome-wide BH-FDR + a |ΔPSI| floor (applied by the caller) select the ASJs.

use std::collections::BTreeMap;

use super::copy_split::{allele_at, intron_chain_of, AlignedRead};

/// ASJ scan thresholds (defaults mirror the python).
#[derive(Clone, Copy, Debug)]
pub struct AsjParams {
    pub d_min: u32,    // anchor SNP min depth
    pub c_min: u32,    // anchor SNP min reads per allele
    pub bal_lo: f64,   // anchor balanced window: minor >= bal_lo * depth
    pub bal_hi: f64,   // anchor balanced window: major <= bal_hi * depth
    pub min_j: u32,    // junction observed >= this many times to be tested
    pub min_span: u32, // >= this many spanning reads PER ALLELE to test a junction
    pub dpsi: f64,     // report |ΔPSI| >= this (effect size floor; applied by the caller)
}
impl Default for AsjParams {
    fn default() -> Self {
        AsjParams { d_min: 12, c_min: 5, bal_lo: 0.35, bal_hi: 0.65, min_j: 3, min_span: 5, dpsi: 0.30 }
    }
}

fn base_idx(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn ref_end(r: &AlignedRead) -> u64 {
    let mut rpos = r.ref_start;
    for &(op, len) in &r.cigar {
        if matches!(op, 'M' | '=' | 'X' | 'D' | 'N') {
            rpos += len;
        }
    }
    rpos
}

/// Per-ref-position ACGT counts over `[start, end)`, one CIGAR walk per read.
fn pileup(reads: &[AlignedRead], start: u64, end: u64) -> BTreeMap<u64, [u32; 4]> {
    let mut pile: BTreeMap<u64, [u32; 4]> = BTreeMap::new();
    for r in reads {
        let mut qpos = 0usize;
        let mut rpos = r.ref_start;
        for &(op, len) in &r.cigar {
            match op {
                'M' | '=' | 'X' => {
                    for k in 0..len {
                        let rp = rpos + k;
                        if rp >= start && rp < end {
                            if let Some(bi) = r.seq.get(qpos + k as usize).copied().and_then(base_idx) {
                                pile.entry(rp).or_insert([0; 4])[bi] += 1;
                            }
                        }
                    }
                    qpos += len as usize;
                    rpos += len;
                }
                'D' | 'N' => rpos += len,
                'I' | 'S' => qpos += len as usize,
                _ => {}
            }
        }
    }
    pile
}

/// het ANCHOR: the most-balanced, best-supported biallelic het SNP in `[start, end)`, as `(pos, major, minor)`
/// bases. `None` if no balanced het column qualifies. Deterministic (BTreeMap ascending → first max wins,
/// matching the python's first-max-over-columns).
pub fn call_anchor(reads: &[AlignedRead], start: u64, end: u64, p: &AsjParams) -> Option<(u64, u8, u8)> {
    const AC: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut best: Option<(u32, u64, u8, u8)> = None; // (minor_count, pos, major, minor)
    for (&pos, counts) in pileup(reads, start, end).iter() {
        let depth: u32 = counts.iter().sum();
        let mut idx = [0usize, 1, 2, 3];
        idx.sort_by(|&a, &b| counts[b].cmp(&counts[a])); // stable desc
        let (major, minor) = (counts[idx[0]], counts[idx[1]]);
        if depth >= p.d_min
            && minor >= p.c_min
            && (minor as f64) >= p.bal_lo * depth as f64
            && (major as f64) <= p.bal_hi * depth as f64
            && best.map_or(true, |b| minor > b.0)
        {
            best = Some((minor, pos, AC[idx[0]], AC[idx[1]]));
        }
    }
    best.map(|(_, pos, x, y)| (pos, x, y))
}

/// One tested junction at an anchor: per-allele used/spanning, PSI, |ΔPSI|, Fisher p.
#[derive(Clone, Debug)]
pub struct AsjCall {
    pub anchor_pos: u64,
    pub ax: u8,
    pub ay: u8,
    pub donor: u64,
    pub acceptor: u64,
    pub used_x: u32,
    pub span_x: u32,
    pub used_y: u32,
    pub span_y: u32,
    pub psi_x: f64,
    pub psi_y: f64,
    pub dpsi: f64,
    pub p: f64,
}

/// Scan a gene `[start, end)`: call the anchor, partition reads by allele, and test every junction seen
/// `>= min_j` times. Returns ALL tested calls with their Fisher p + |ΔPSI| (FDR + the ΔPSI floor are applied
/// genome-wide by the caller). Constitutive-in-both-alleles junctions are skipped (focuses the test burden).
pub fn scan_gene_asj(reads: &[AlignedRead], start: u64, end: u64, p: &AsjParams) -> Vec<AsjCall> {
    let (pos, ax, ay) = match call_anchor(reads, start, end, p) {
        Some(a) => a,
        None => return Vec::new(),
    };
    // per read covering the anchor: (allele, intron set, ref_start, ref_end)
    let mut parts: Vec<(u8, Vec<(u64, u64)>, u64, u64)> = Vec::new();
    for r in reads {
        let re = ref_end(r);
        if r.ref_start > pos || re <= pos {
            continue;
        }
        match allele_at(r, pos) {
            Some(b) if b == ax || b == ay => parts.push((b, intron_chain_of(r), r.ref_start, re)),
            _ => {}
        }
    }
    let nx = parts.iter().filter(|(a, ..)| *a == ax).count() as u32;
    let ny = parts.iter().filter(|(a, ..)| *a == ay).count() as u32;
    if nx < p.c_min || ny < p.c_min {
        return Vec::new();
    }
    // candidate junctions: introns seen >= min_j across the partitioned reads.
    let mut jc: BTreeMap<(u64, u64), u32> = BTreeMap::new();
    for (_, js, _, _) in &parts {
        for &j in js {
            *jc.entry(j).or_insert(0) += 1;
        }
    }
    let mut out = Vec::new();
    for (&(d, a), &n) in jc.iter() {
        if n < p.min_j {
            continue;
        }
        let (mut ux, mut sx, mut uy, mut sy) = (0u32, 0u32, 0u32, 0u32);
        for (al, js, rs, re) in &parts {
            if *rs <= d && *re >= a {
                let uses = js.contains(&(d, a)) as u32;
                if *al == ax {
                    sx += 1;
                    ux += uses;
                } else {
                    sy += 1;
                    uy += uses;
                }
            }
        }
        if sx < p.min_span || sy < p.min_span {
            continue;
        }
        let (psi_x, psi_y) = (ux as f64 / sx as f64, uy as f64 / sy as f64);
        // skip junctions constitutive in BOTH alleles (can't be allele-specific).
        if (psi_x >= 0.98 && psi_y >= 0.98) || (psi_x <= 0.02 && psi_y <= 0.02) {
            continue;
        }
        out.push(AsjCall {
            anchor_pos: pos,
            ax,
            ay,
            donor: d,
            acceptor: a,
            used_x: ux,
            span_x: sx,
            used_y: uy,
            span_y: sy,
            psi_x,
            psi_y,
            dpsi: (psi_x - psi_y).abs(),
            p: fisher_exact_2x2(ux, sx - ux, uy, sy - uy),
        });
    }
    out
}

/// A TRANSVERSION (purine↔pyrimidine: A/C, A/T, G/C, G/T) is unambiguously GENETIC; a transition (A/G or C/T)
/// could be an RNA-editing site. The anchor/PSV base-change type is the genetic-vs-edit confound control.
pub fn is_transversion(a: u8, b: u8) -> bool {
    let purine = |x: u8| matches!(x.to_ascii_uppercase(), b'A' | b'G');
    purine(a) != purine(b)
}

/// Anchor→junction distance = `min(|anchor − donor|, |anchor − acceptor|)`. `<= ~100bp` (splice-proximal) makes
/// a splice-disrupting-variant mechanism plausible; a distal anchor fully flipping a junction is more likely
/// two distinct transcript populations than allele-modulated splicing of one gene.
pub fn anchor_junction_dist(anchor: u64, donor: u64, acceptor: u64) -> u64 {
    let d = anchor.abs_diff(donor);
    let a = anchor.abs_diff(acceptor);
    d.min(a)
}

/// Benjamini-Hochberg q-values for `pvals` (returned in input order). Genome-wide multiple-testing control.
pub fn bh_fdr(pvals: &[f64]) -> Vec<f64> {
    let m = pvals.len();
    if m == 0 {
        return Vec::new();
    }
    let mut idx: Vec<usize> = (0..m).collect();
    idx.sort_by(|&i, &j| pvals[i].partial_cmp(&pvals[j]).unwrap_or(std::cmp::Ordering::Equal));
    let mut q = vec![0.0; m];
    let mut prev = 1.0_f64;
    for rank in (0..m).rev() {
        let i = idx[rank];
        let v = (pvals[i] * m as f64 / (rank + 1) as f64).min(prev);
        q[i] = v;
        prev = v;
    }
    q
}

/// Lanczos approximation of `ln Γ(x)`.
fn lgamma(x: f64) -> f64 {
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

fn lchoose(n: u32, k: u32) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    lgamma(n as f64 + 1.0) - lgamma(k as f64 + 1.0) - lgamma((n - k) as f64 + 1.0)
}

/// Two-sided Fisher exact test on the 2×2 table `[[a, b], [c, d]]` (matches `scipy.stats.fisher_exact`):
/// sum the hypergeometric probabilities of all tables (same margins) no more likely than the observed.
pub fn fisher_exact_2x2(a: u32, b: u32, c: u32, d: u32) -> f64 {
    let (r1, r2, col1, n) = (a + b, c + d, a + c, a + b + c + d);
    if n == 0 {
        return 1.0;
    }
    let lcn = lchoose(n, col1);
    let pmf = |k: u32| (lchoose(r1, k) + lchoose(r2, col1 - k) - lcn).exp();
    let p_obs = pmf(a);
    let (lo, hi) = (col1.saturating_sub(r2), col1.min(r1));
    let mut total = 0.0;
    for k in lo..=hi {
        let pk = pmf(k);
        if pk <= p_obs * (1.0 + 1e-7) {
            total += pk;
        }
    }
    total.min(1.0)
}

// ============================ MULTI-SNP HAPLOTYPE PHASING ============================
// Single-anchor phasing needs a read to cover ONE specific het SNP. Multi-SNP phasing assigns each read to
// one of the two DIPLOID haplotypes using ALL het SNPs (read-backed 2-means), so reads covering any subset
// of SNPs get phased → more reach, and the phasing QUALITY (how cleanly reads bipartition) is itself a
// confound check (clean 2-way = diploid het; messy/>2 = paralog/artifact). Port of the python multi-SNP scan.

const MAX_SNP: usize = 40; // cap het SNPs per gene (highest-coverage)
const PHASE_ITERS: usize = 6;
const AGREE: f64 = 0.80; // a read is confidently phased if it agrees with one hap at >= this frac

/// All balanced het SNP columns in `[start, end)` as `(pos, major, minor)`, highest-coverage first, capped at
/// `MAX_SNP`. (The multi-SNP scan uses a slightly looser balance window than the single-anchor `call_anchor`.)
pub fn het_snps(reads: &[AlignedRead], start: u64, end: u64, p: &AsjParams) -> Vec<(u64, u8, u8)> {
    const AC: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut snps: Vec<(u32, u64, u8, u8)> = Vec::new(); // (minor, pos, major, minor)
    for (&pos, counts) in pileup(reads, start, end).iter() {
        let depth: u32 = counts.iter().sum();
        let mut idx = [0usize, 1, 2, 3];
        idx.sort_by(|&a, &b| counts[b].cmp(&counts[a]));
        let (major, minor) = (counts[idx[0]], counts[idx[1]]);
        if depth >= p.d_min
            && minor >= p.c_min
            && (minor as f64) >= p.bal_lo * depth as f64
            && (major as f64) <= p.bal_hi * depth as f64
        {
            snps.push((minor, pos, AC[idx[0]], AC[idx[1]]));
        }
    }
    snps.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1))); // minor desc, then pos asc (deterministic)
    snps.truncate(MAX_SNP);
    snps.into_iter().map(|(_, pos, a, b)| (pos, a, b)).collect()
}

fn read_alleles(read: &AlignedRead, snps: &[(u64, u8, u8)]) -> Vec<Option<u8>> {
    snps.iter().map(|&(pos, _, _)| allele_at(read, pos)).collect()
}

/// Deterministic 2-means haplotype phasing. `reads[r][i]` = read r's base at SNP i (or None). Seed = SNP 0's
/// two alleles; iterate (assign reads to the closer hap, re-estimate each hap's consensus per SNP). Returns
/// `hap[r]` in {Some(0), Some(1), None}` (None = not confidently phased) and `n_diff` = SNPs where the two
/// haplotype consensuses actually differ (the phasing's information content).
pub fn phase_multisnp(reads: &[Vec<Option<u8>>], snps: &[(u64, u8, u8)]) -> (Vec<Option<usize>>, usize) {
    let k = snps.len();
    let mut h0: Vec<Option<u8>> = vec![None; k];
    let mut h1: Vec<Option<u8>> = vec![None; k];
    h0[0] = Some(snps[0].1);
    h1[0] = Some(snps[0].2);
    let mut assign: Vec<Option<usize>> = vec![None; reads.len()];
    for _ in 0..PHASE_ITERS {
        for (ri, al) in reads.iter().enumerate() {
            let (mut a0, mut a1, mut cov) = (0u32, 0u32, 0u32);
            for i in 0..k {
                if let (Some(b), Some(hb0)) = (al[i], h0[i]) {
                    cov += 1;
                    if b == hb0 {
                        a0 += 1;
                    }
                    if h1[i] == Some(b) {
                        a1 += 1;
                    }
                }
            }
            assign[ri] = if cov == 0 { None } else if a0 >= a1 { Some(0) } else { Some(1) };
        }
        for i in 0..k {
            for grp in 0..2 {
                let (mut t_a, mut t_b) = (0u32, 0u32);
                for (ri, al) in reads.iter().enumerate() {
                    if assign[ri] == Some(grp) {
                        match al[i] {
                            Some(b) if b == snps[i].1 => t_a += 1,
                            Some(b) if b == snps[i].2 => t_b += 1,
                            _ => {}
                        }
                    }
                }
                let cons = if t_a == 0 && t_b == 0 {
                    None
                } else if t_a >= t_b {
                    Some(snps[i].1)
                } else {
                    Some(snps[i].2)
                };
                if grp == 0 {
                    h0[i] = cons;
                } else {
                    h1[i] = cons;
                }
            }
        }
    }
    let n_diff = (0..k).filter(|&i| h0[i].is_some() && h1[i].is_some() && h0[i] != h1[i]).count();
    let mut hap = vec![None; reads.len()];
    for (ri, al) in reads.iter().enumerate() {
        let (mut a0, mut a1, mut cov) = (0u32, 0u32, 0u32);
        for i in 0..k {
            if h0[i].is_none() || h1[i].is_none() || h0[i] == h1[i] {
                continue;
            }
            if let Some(b) = al[i] {
                cov += 1;
                if Some(b) == h0[i] {
                    a0 += 1;
                } else if Some(b) == h1[i] {
                    a1 += 1;
                }
            }
        }
        if cov == 0 {
            continue;
        }
        let (best, frac) = if a0 >= a1 { (0, a0 as f64 / cov as f64) } else { (1, a1 as f64 / cov as f64) };
        if frac >= AGREE {
            hap[ri] = Some(best);
        }
    }
    (hap, n_diff)
}

/// A multi-SNP ASJ call: the gene's het-SNP count + phase quality (confound metric) + the per-haplotype test.
#[derive(Clone, Debug)]
pub struct AsjMultiCall {
    pub n_snps: usize,
    pub phase_qual: f64,
    pub donor: u64,
    pub acceptor: u64,
    pub used0: u32,
    pub span0: u32,
    pub used1: u32,
    pub span1: u32,
    pub psi0: f64,
    pub psi1: f64,
    pub dpsi: f64,
    pub p: f64,
}

/// Multi-SNP ASJ scan: het SNPs → 2-means phasing → per-haplotype junction test. A superset of the
/// single-anchor scan (with one het SNP it degenerates to it) — more reach, plus `phase_qual` as a confound.
pub fn scan_gene_asj_multisnp(reads: &[AlignedRead], start: u64, end: u64, p: &AsjParams) -> Vec<AsjMultiCall> {
    phase_and_test(reads, &het_snps(reads, start, end, p), p)
}

/// Copy-distinguishing PSV columns (like [`het_snps`] but WITHOUT the diploid-balance window — paralog copies
/// have unequal expression, so the two copy-alleles need not be 50/50; only `>= c_min` reads each). Returned
/// `(pos, allele_a, allele_b)`, highest-minor first, capped at `MAX_SNP`.
pub fn discover_psv_alleles(reads: &[AlignedRead], start: u64, end: u64, p: &AsjParams) -> Vec<(u64, u8, u8)> {
    const AC: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut psvs: Vec<(u32, u64, u8, u8)> = Vec::new();
    for (&pos, counts) in pileup(reads, start, end).iter() {
        let mut idx = [0usize, 1, 2, 3];
        idx.sort_by(|&a, &b| counts[b].cmp(&counts[a]));
        let minor = counts[idx[1]];
        if minor >= p.c_min {
            psvs.push((minor, pos, AC[idx[0]], AC[idx[1]]));
        }
    }
    psvs.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    psvs.truncate(MAX_SNP);
    psvs.into_iter().map(|(_, pos, a, b)| (pos, a, b)).collect()
}

/// CONNECT ASJ TO COPIES: a junction whose usage depends on which COPY a read belongs to — the "allele" is
/// the copy identity. Same phasing + test as the het-SNP scan, but the variants are copy-distinguishing PSVs,
/// so the two phased groups are the two paralog COPIES and a call is a COPY-SPECIFIC junction. (This is the
/// thesis join: copy DETECTION/ASSIGNMENT (PSVs) meets junction MODELLING (ASJ).) `phase_qual` again flags a
/// clean 2-copy split vs a messy/>2 one. Pairwise on the two dominant copies.
pub fn scan_gene_copy_specific_junctions(reads: &[AlignedRead], start: u64, end: u64, p: &AsjParams) -> Vec<AsjMultiCall> {
    phase_and_test(reads, &discover_psv_alleles(reads, start, end, p), p)
}

/// Shared core: phase reads into 2 groups by their multi-VARIANT haplotype (het SNPs → diploid haplotypes;
/// PSVs → paralog copies) and test each junction's PSI per group. Used by both the multi-SNP ASJ scan and the
/// copy-specific-junction scan.
fn phase_and_test(reads: &[AlignedRead], variants: &[(u64, u8, u8)], p: &AsjParams) -> Vec<AsjMultiCall> {
    let snps = variants;
    if snps.is_empty() {
        return Vec::new();
    }
    let mut ral: Vec<Vec<Option<u8>>> = Vec::new();
    let mut rmeta: Vec<(Vec<(u64, u64)>, u64, u64)> = Vec::new();
    for r in reads {
        let al = read_alleles(r, snps);
        if al.iter().all(Option::is_none) {
            continue;
        }
        ral.push(al);
        rmeta.push((intron_chain_of(r), r.ref_start, ref_end(r)));
    }
    if ral.len() < 2 * p.c_min as usize {
        return Vec::new();
    }
    let (hap, n_diff) = phase_multisnp(&ral, snps);
    if n_diff < 1 {
        return Vec::new();
    }
    let phased: Vec<usize> = (0..hap.len()).filter(|&i| hap[i].is_some()).collect();
    let snp_covering = ral.iter().filter(|al| al.iter().any(Option::is_some)).count();
    if snp_covering == 0 {
        return Vec::new();
    }
    let phase_qual = phased.len() as f64 / snp_covering as f64;
    let n0 = phased.iter().filter(|&&i| hap[i] == Some(0)).count() as u32;
    let n1 = phased.iter().filter(|&&i| hap[i] == Some(1)).count() as u32;
    if n0 < p.c_min || n1 < p.c_min {
        return Vec::new();
    }
    let mut jc: BTreeMap<(u64, u64), u32> = BTreeMap::new();
    for &i in &phased {
        for &j in &rmeta[i].0 {
            *jc.entry(j).or_insert(0) += 1;
        }
    }
    let mut out = Vec::new();
    for (&(d, a), &nj) in jc.iter() {
        if nj < p.min_j {
            continue;
        }
        let mut cnt = [[0u32; 2]; 2]; // [hap][used, spanning]
        for &i in &phased {
            let (js, rs, re) = (&rmeta[i].0, rmeta[i].1, rmeta[i].2);
            if rs <= d && re >= a {
                let h = hap[i].unwrap();
                cnt[h][1] += 1;
                if js.contains(&(d, a)) {
                    cnt[h][0] += 1;
                }
            }
        }
        let (u0, s0, u1, s1) = (cnt[0][0], cnt[0][1], cnt[1][0], cnt[1][1]);
        if s0 < p.min_span || s1 < p.min_span {
            continue;
        }
        let (psi0, psi1) = (u0 as f64 / s0 as f64, u1 as f64 / s1 as f64);
        if (psi0 >= 0.98 && psi1 >= 0.98) || (psi0 <= 0.02 && psi1 <= 0.02) {
            continue;
        }
        out.push(AsjMultiCall {
            n_snps: snps.len(),
            phase_qual,
            donor: d,
            acceptor: a,
            used0: u0,
            span0: s0,
            used1: u1,
            span1: s1,
            psi0,
            psi1,
            dpsi: (psi0 - psi1).abs(),
            p: fisher_exact_2x2(u0, s0 - u0, u1, s1 - u1),
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fisher_matches_scipy() {
        let approx = |a: f64, b: f64| (a - b).abs() < 1e-5;
        assert!(approx(fisher_exact_2x2(8, 2, 1, 9), 0.005_477_49), "{}", fisher_exact_2x2(8, 2, 1, 9));
        assert!(approx(fisher_exact_2x2(5, 5, 5, 5), 1.0));
        assert!(approx(fisher_exact_2x2(10, 5, 3, 12), 0.025_327_69));
        assert!(approx(fisher_exact_2x2(3, 7, 8, 2), 0.069_778_52));
        assert!(fisher_exact_2x2(14, 0, 0, 18) < 1e-6);
    }

    #[test]
    fn transversion_and_distance() {
        assert!(is_transversion(b'A', b'C')); // purine vs pyrimidine
        assert!(is_transversion(b'G', b'T'));
        assert!(!is_transversion(b'A', b'G')); // both purine -> transition
        assert!(!is_transversion(b'C', b'T')); // both pyrimidine -> transition
        assert_eq!(anchor_junction_dist(100, 95, 300), 5); // nearest boundary
        assert_eq!(anchor_junction_dist(290, 95, 300), 10);
    }

    #[test]
    fn bh_fdr_monotone() {
        let q = bh_fdr(&[0.001, 0.5, 0.01, 0.9]);
        assert!(q[0] <= q[2] && q[2] <= q[1]); // smaller p -> smaller q
        assert!(q.iter().all(|&v| (0.0..=1.0).contains(&v)));
    }

    /// Build a read: exon up to `intron.0`, intron, exon to `end`; or no intron (retained). `base` sits at
    /// ref position `anchor` (within the first exon, which starts at ref 0).
    fn read(anchor: u64, base: u8, end: u64, intron: Option<(u64, u64)>) -> AlignedRead {
        let (cigar, seqlen) = match intron {
            Some((d, a)) => (
                vec![('M', d), ('N', a - d), ('M', end - a)],
                (d + (end - a)) as usize,
            ),
            None => (vec![('M', end)], end as usize),
        };
        let mut seq = vec![b'A'; seqlen];
        seq[anchor as usize] = base; // anchor is in the first exon (ref==query offset there)
        AlignedRead { ref_start: 0, cigar, seq }
    }

    #[test]
    fn scan_finds_allele_specific_junction() {
        // anchor at ref 50: 6 reads carry G + USE junction (60,160); 6 carry T + RETAIN it. Both alleles span
        // the junction -> a clean ASJ: PSI(G)=1, PSI(T)=0, |dPSI|=1.
        let mut reads = Vec::new();
        for _ in 0..6 {
            reads.push(read(50, b'G', 260, Some((60, 160))));
        }
        for _ in 0..6 {
            reads.push(read(50, b'T', 260, None));
        }
        let calls = scan_gene_asj(&reads, 0, 260, &AsjParams::default());
        assert_eq!(calls.len(), 1, "one tested ASJ junction");
        let c = &calls[0];
        assert_eq!((c.donor, c.acceptor), (60, 160));
        assert!((c.dpsi - 1.0).abs() < 1e-9, "full allele switch");
        assert!(c.p < 0.01, "significant: {}", c.p);
        assert_eq!(c.anchor_pos, 50);
    }

    #[test]
    fn no_anchor_no_calls() {
        // homozygous locus (all G at 50) -> no balanced het anchor -> no calls.
        let reads: Vec<AlignedRead> = (0..12).map(|_| read(50, b'G', 260, Some((60, 160)))).collect();
        assert!(scan_gene_asj(&reads, 0, 260, &AsjParams::default()).is_empty());
    }

    #[test]
    fn multisnp_phases_two_snps_and_finds_asj() {
        // two het SNPs (30, 50). haplotype A: G at both + USES junction (60,160); haplotype B: T at both +
        // RETAINS it. The 2-means phasing must bipartition cleanly and recover the full-switch ASJ.
        let mk = |b: u8, intron: bool| {
            let (cigar, slen) = if intron {
                (vec![('M', 60u64), ('N', 100), ('M', 100)], 160usize)
            } else {
                (vec![('M', 260u64)], 260usize)
            };
            let mut seq = vec![b'A'; slen];
            seq[30] = b;
            seq[50] = b;
            AlignedRead { ref_start: 0, cigar, seq }
        };
        let mut reads = Vec::new();
        for _ in 0..7 {
            reads.push(mk(b'G', true));
        }
        for _ in 0..7 {
            reads.push(mk(b'T', false));
        }
        let calls = scan_gene_asj_multisnp(&reads, 0, 260, &AsjParams::default());
        assert_eq!(calls.len(), 1, "one ASJ via multi-SNP phasing");
        assert_eq!((calls[0].donor, calls[0].acceptor), (60, 160));
        assert!((calls[0].dpsi - 1.0).abs() < 1e-9, "full allele switch");
        assert_eq!(calls[0].n_snps, 2, "both het SNPs used");
        assert!(calls[0].phase_qual > 0.9, "clean diploid bipartition");
    }

    #[test]
    fn copy_specific_junction_from_unbalanced_psvs() {
        // two COPIES at unequal expression (10 vs 5 reads) distinguished by 2 PSVs; copy0 USES junction
        // (60,160), copy1 RETAINS. The het scan rejects this (10:5 is not a balanced diploid het), but the
        // copy-specific scan accepts it (PSVs need no balance) -> a copy-specific junction.
        let mk = |b: u8, intron: bool| {
            let (cigar, slen) = if intron {
                (vec![('M', 60u64), ('N', 100), ('M', 100)], 160usize)
            } else {
                (vec![('M', 260u64)], 260usize)
            };
            let mut seq = vec![b'A'; slen];
            seq[30] = b;
            seq[50] = b;
            AlignedRead { ref_start: 0, cigar, seq }
        };
        let mut reads = Vec::new();
        for _ in 0..10 {
            reads.push(mk(b'G', true));
        }
        for _ in 0..5 {
            reads.push(mk(b'T', false));
        }
        assert!(
            scan_gene_asj_multisnp(&reads, 0, 260, &AsjParams::default()).is_empty(),
            "10:5 is unbalanced -> not a diploid het anchor"
        );
        let calls = scan_gene_copy_specific_junctions(&reads, 0, 260, &AsjParams::default());
        assert_eq!(calls.len(), 1, "copy-specific junction recovered from unbalanced PSVs");
        assert!((calls[0].dpsi - 1.0).abs() < 1e-9, "full copy switch");
    }
}

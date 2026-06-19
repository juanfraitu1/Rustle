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
}

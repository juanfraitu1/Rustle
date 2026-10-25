//! Statistics and columns (o3port.py `poisson_tail`, `_small_tail`, `err`, `informative`, `patterns`, `compatible`), with the two
//! CPython numerics they depend on: `math.lgamma` (CPython's own Lanczos, not libm) and `sum()` of floats (Neumaier since 3.12).

use super::{Cols, Obs, ObsX, Pattern, ReadId, ALPHA, MIN_Q, MIN_READS};
use std::collections::{BTreeMap, BTreeSet};

/// CPython `math.lgamma` (Modules/mathmodule.c `m_lgamma`, finite positive arguments: the only ones used here).
pub fn lgamma_cpython(x: f64) -> f64 {
    if x == x.floor() && x <= 2.0 {
        return if x <= 0.0 { f64::INFINITY } else { 0.0 }; // lgamma(1) = lgamma(2) = 0.0
    }
    let absx = x.abs();
    if absx < 1e-20 {
        return -absx.ln();
    }
    let mut r = lanczos_sum(absx).ln() - LANCZOS_G;
    r += (absx - 0.5) * ((absx + LANCZOS_G - 0.5).ln() - 1.0);
    assert!(x > 0.0, "lgamma_cpython: negative arguments are not ported");
    r
}

const LANCZOS_G: f64 = 6.024680040776729583740234375;
const LANCZOS_NUM: [f64; 13] = [
    23531376880.410759688572007674451636754734846804940,
    42919803642.649098768957899047001988850926355848959,
    35711959237.355668049440185451547166705960488635843,
    17921034426.037209699919755754458931112671403265390,
    6039542586.3520280050642916443072979210699388420708,
    1439720407.3117216736632230727949123939715485786772,
    248874557.86205415651146038641322942321632125127801,
    31426415.585400194380614231628318205362874684987640,
    2876370.6289353724412254090516208496135991145378768,
    186056.26539522349504029498971604569928220784236328,
    8071.6720023658162106380029022722506138218516325024,
    210.82427775157934587250973392071336271166969580291,
    2.5066282746310002701649081771338373386264310793408,
];
const LANCZOS_DEN: [f64; 13] = [
    0.0, 39916800.0, 120543840.0, 150917976.0, 105258076.0, 45995730.0, 13339535.0, 2637558.0, 357423.0, 32670.0, 1925.0, 66.0,
    1.0,
];

fn lanczos_sum(x: f64) -> f64 {
    let (mut num, mut den) = (0.0f64, 0.0f64);
    if x < 5.0 {
        for i in (0..13).rev() {
            num = num * x + LANCZOS_NUM[i];
            den = den * x + LANCZOS_DEN[i];
        }
    } else {
        for i in 0..13 {
            num = num / x + LANCZOS_NUM[i];
            den = den / x + LANCZOS_DEN[i];
        }
    }
    num / den
}

/// CPython `sum()` of floats starting from the int 0 (Neumaier compensation, bltinmodule.c).
pub fn py_sum(xs: impl IntoIterator<Item = f64>) -> f64 {
    let (mut f, mut c) = (0.0f64, 0.0f64);
    for x in xs {
        let t = f + x;
        if f.abs() >= x.abs() {
            c += (f - t) + x;
        } else {
            c += (x - t) + f;
        }
        f = t;
    }
    if c != 0.0 && c.is_finite() {
        f += c;
    }
    f
}

/// P(X >= k) for X ~ Poisson(lam).
pub fn poisson_tail(k: usize, lam: f64) -> f64 {
    if k == 0 {
        return 1.0;
    }
    if lam <= 0.0 {
        return 0.0;
    }
    let mut term = (-lam).exp();
    let mut cdf = term;
    for j in 1..k {
        term *= lam / j as f64;
        cdf += term;
    }
    if cdf < 1.0 - 1e-12 {
        (1.0 - cdf).max(0.0)
    } else {
        small_tail(k, lam)
    }
}

/// The tail by its leading terms when 1 - cdf underflows (lam small); log space, so a large k cannot overflow.
fn small_tail(k: usize, lam: f64) -> f64 {
    let mut lt = -lam + k as f64 * lam.ln() - lgamma_cpython((k + 1) as f64);
    let (mut s, mut j) = (0.0f64, k);
    while lt > -700.0 && j < k + 50 {
        s += lt.exp();
        j += 1;
        lt += lam.ln() - (j as f64).ln();
    }
    s
}

/// A read's per-base error probability at a column: its empirical rate when known, else its base quality.
pub fn err(x: &ObsX) -> f64 {
    x.eps.unwrap_or_else(|| 10f64.powf(-(x.qual as f64) / 10.0))
}

/// `informative({c: [x for x in xs if keep(x.read)] for c, xs in obs.items()}, n_tests)`.
pub fn informative(obs: &Obs, keep: impl Fn(ReadId) -> bool, n_tests: Option<usize>) -> Cols {
    let mut cand: Vec<(i64, Vec<&ObsX>, Vec<(u8, usize)>, Vec<u8>)> = Vec::new();
    for (&col, xs) in obs {
        let good: Vec<&ObsX> = xs.iter().filter(|x| keep(x.read) && x.qual >= MIN_Q && b"ACGT".contains(&x.base)).collect();
        let cnt = counter(good.iter().map(|x| x.base));
        let alleles: Vec<u8> = cnt.iter().filter(|(_, n)| *n >= MIN_READS).map(|(b, _)| *b).collect();
        if alleles.len() >= 2 {
            cand.push((col, good, cnt, alleles));
        }
    }
    let m = n_tests.unwrap_or(cand.len().max(1));
    let mut out = Cols::new();
    for (col, good, cnt, alleles) in cand {
        let lam = py_sum(good.iter().map(|x| err(x) / 3.0));
        let n_of = |b: u8| cnt.iter().find(|(c, _)| *c == b).map_or(0, |(_, n)| *n);
        let major = *alleles.iter().max_by_key(|&&b| (n_of(b), b)).unwrap();
        let mut real = BTreeSet::from([major]);
        for &b in &alleles {
            if b != major && poisson_tail(n_of(b), lam) < ALPHA / m as f64 {
                real.insert(b);
            }
        }
        if real.len() >= 2 {
            out.insert(col, real);
        }
    }
    out
}

/// `collections.Counter` of bytes: (item, count) in first-appearance order.
pub fn counter(xs: impl IntoIterator<Item = u8>) -> Vec<(u8, usize)> {
    let mut out: Vec<(u8, usize)> = Vec::new();
    for b in xs {
        match out.iter_mut().find(|(c, _)| *c == b) {
            Some(e) => e.1 += 1,
            None => out.push((b, 1)),
        }
    }
    out
}

/// `Counter.most_common(1)[0]`: the highest count, the first inserted on a tie.
pub fn most_common<T: Copy>(cnt: &[(T, usize)]) -> Option<(T, usize)> {
    let mut best: Option<(T, usize)> = None;
    for &(b, n) in cnt {
        if best.map_or(true, |(_, m)| n > m) {
            best = Some((b, n));
        }
    }
    best
}

/// `{read: {col: allele}}` restricted to informative columns, a base that is no supported allele counted as missing.
pub fn patterns(reads: &BTreeSet<ReadId>, obs: &Obs, cols: &Cols) -> BTreeMap<ReadId, Pattern> {
    let mut pat: BTreeMap<ReadId, Pattern> = reads.iter().map(|&r| (r, Pattern::new())).collect();
    for (col, real) in cols {
        for x in obs.get(col).map_or(&[][..], |v| &v[..]) {
            if x.qual >= MIN_Q && real.contains(&x.base) {
                if let Some(p) = pat.get_mut(&x.read) {
                    p.insert(*col, x.base);
                }
            }
        }
    }
    pat
}

pub fn compatible(a: &Pattern, b: &Pattern) -> bool {
    let mut shared = false;
    for (c, x) in a {
        if let Some(y) = b.get(c) {
            if x != y {
                return false;
            }
            shared = true;
        }
    }
    shared
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hex(x: f64) -> String {
        x.to_le_bytes().iter().map(|b| format!("{b:02x}")).collect()
    }

    #[test]
    fn lgamma_matches_cpython() {
        // struct.pack('<d', math.lgamma(x)).hex() under Python 3.14.4
        for (x, h) in [
            (0.5, "c0a1e748d050e23f"),
            (1.0, "0000000000000000"),
            (2.0, "0000000000000000"),
            (3.7, "889b355e62d9f63f"),
            (10.0, "cef7a721899a2940"),
            (25.5, "2238cf3ed0314c40"),
            (171.2, "4130bcb2ce1c8640"),
            (3.0, "ec39fafe422ee63f"),
            (7.0, "cb01cf3a27511a40"),
            (51.0, "321febdd498f6240"),
        ] {
            assert_eq!(hex(lgamma_cpython(x)), h, "lgamma({x})");
        }
    }

    #[test]
    fn py_sum_is_compensated() {
        // sum([0.1] * 10) == 1.0 in Python >= 3.12 (plain addition gives 0.9999999999999999)
        assert_eq!(py_sum(std::iter::repeat(0.1).take(10)), 1.0);
        assert_eq!(py_sum([1e100, 1.0, -1e100]), 1.0);
        assert_eq!(py_sum(std::iter::empty()), 0.0);
    }

    #[test]
    fn poisson_tail_small_and_large() {
        assert!((poisson_tail(1, 0.5) - (1.0 - 2.718281828f64.powf(-0.5))).abs() < 1e-6);
        let t = poisson_tail(5, 1e-6);
        assert!(t > 0.0 && t < 1e-28);
        let t = poisson_tail(400, 1e-3);
        assert!(t == 0.0 || t < 1e-300);
        assert_eq!(poisson_tail(0, 1.0), 1.0);
        assert_eq!(poisson_tail(3, 0.0), 0.0);
    }

    fn x(read: ReadId, base: u8, qual: u8) -> ObsX {
        ObsX { read, base, qual, eps: None }
    }

    #[test]
    fn informative_needs_two_carriers_and_significance() {
        // 10 A / 10 G at Q30 -> informative; a single G -> not; 2 G among 1,000 reads at Q10 -> explained by error
        let mut obs = Obs::new();
        obs.insert(1, (0..10).map(|i| x(i, b'A', 30)).chain((100..110).map(|i| x(i, b'G', 30))).collect());
        obs.insert(2, (0..19).map(|i| x(i, b'A', 30)).chain([x(100, b'G', 30)]).collect());
        obs.insert(3, (0..998).map(|i| x(i, b'A', 10)).chain([x(2000, b'G', 10), x(2001, b'G', 10)]).collect());
        let cols = informative(&obs, |_| true, None);
        assert_eq!(cols.len(), 1);
        assert_eq!(cols[&1], BTreeSet::from([b'A', b'G']));
        // the read filter: without the G carriers column 1 is not informative
        assert!(informative(&obs, |r| r < 100, None).is_empty());
    }

    #[test]
    fn patterns_drop_unsupported_and_low_quality_bases() {
        let mut obs = Obs::new();
        obs.insert(5, vec![x(0, b'A', 30), x(1, b'G', 30), x(2, b'T', 30), x(3, b'A', 5)]);
        let mut cols = Cols::new();
        cols.insert(5, BTreeSet::from([b'A', b'G']));
        let p = patterns(&BTreeSet::from([0, 1, 2, 3]), &obs, &cols);
        assert_eq!(p[&0], Pattern::from([(5, b'A')]));
        assert_eq!(p[&1], Pattern::from([(5, b'G')]));
        assert!(p[&2].is_empty() && p[&3].is_empty());
        assert!(compatible(&p[&0], &p[&0]));
        assert!(!compatible(&p[&0], &p[&1]));
        assert!(!compatible(&p[&0], &p[&2])); // nothing shared
    }
}

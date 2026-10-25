//! Putative copies: slot holders whose reads behave like another copy's (o3port.py `junction_usage`, `foreign_splicing`,
//! `foreign_clipping`, `foreign_coverage`, `putative_copies`).

use super::pile::{Pile, OP_D, OP_EQ, OP_M, OP_N, OP_S, OP_X};
use super::slots::Matching;
use super::{ReadId, CANONICAL, CLIP_BIN, CLIP_BP, MIN_READS, MIN_SUPPORT};
use std::collections::{BTreeMap, BTreeSet, HashMap};

/// A junction (intron) `(start, end)` on the reference.
pub type Junction = (i64, i64);
/// {read: (alignment start, alignment end, its introns)}.
pub type Usage = BTreeMap<ReadId, (i64, i64, BTreeSet<Junction>)>;

pub fn junction_usage(pile: &Pile) -> Usage {
    let mut out = Usage::new();
    for (&n, rec) in pile {
        let (mut pos, mut js) = (rec.start, BTreeSet::new());
        for &(op, k) in &rec.cigar {
            if op == OP_N {
                js.insert((pos, pos + k as i64));
            }
            if [OP_M, OP_D, OP_N, OP_EQ, OP_X].contains(&op) {
                pos += k as i64;
            }
        }
        out.insert(n, (rec.start, rec.end, js));
    }
    out
}

/// >= MIN_READS of haplotype i's own junctions on the slot copy are non-canonical on the copy's DNA; `site(j)` = (donor, acceptor).
pub fn foreign_splicing(reads_i: &BTreeSet<ReadId>, usage: &Usage, site: impl Fn(Junction) -> (Vec<u8>, Vec<u8>)) -> bool {
    let spans = |j: Junction| reads_i.iter().filter_map(|r| usage.get(r)).filter(|u| u.0 < j.0 && u.1 > j.1).count();
    let mut ji: BTreeMap<Junction, usize> = BTreeMap::new();
    for u in reads_i.iter().filter_map(|r| usage.get(r)) {
        for &j in &u.2 {
            *ji.entry(j).or_default() += 1;
        }
    }
    let n_foreign = ji
        .iter()
        .filter(|(&j, &k)| k >= MIN_READS && 2 * k >= spans(j))
        .filter(|(&j, _)| {
            let (d, a) = site(j);
            !CANONICAL.iter().any(|(cd, ca)| cd.as_bytes() == d.as_slice() && ca.as_bytes() == a.as_slice())
        })
        .count();
    n_foreign >= MIN_READS
}

/// >= half of haplotype i's reads with a record on the slot copy are soft-clipped there by >= CLIP_BP non-tail bases, and
/// >= MIN_SUPPORT of them at one boundary.
pub fn foreign_clipping(reads_i: &BTreeSet<ReadId>, recs: &Pile) -> bool {
    let tail = |q: &[u8]| {
        !q.is_empty() && {
            let a = q.iter().filter(|&&b| b == b'A').count();
            let t = q.iter().filter(|&&b| b == b'T').count();
            a.max(t) as f64 >= 0.7 * q.len() as f64
        }
    };
    let (mut seen, mut hit) = (0usize, 0usize);
    let mut bnd: Vec<((u8, i64), usize)> = Vec::new();
    for r in reads_i {
        let Some(rec) = recs.get(r) else { continue };
        if rec.cigar.is_empty() {
            continue;
        }
        seen += 1;
        let (q, c) = (&rec.seq[..], &rec.cigar);
        let mut ev = Vec::new();
        let (op0, n0) = c[0];
        if op0 == OP_S && n0 >= CLIP_BP && !tail(&q[..(n0 as usize).min(q.len())]) {
            ev.push((b'L', rec.start.div_euclid(CLIP_BIN)));
        }
        let (op1, n1) = c[c.len() - 1];
        if op1 == OP_S && n1 >= CLIP_BP && !tail(&q[q.len() - (n1 as usize).min(q.len())..]) {
            ev.push((b'R', rec.end.div_euclid(CLIP_BIN)));
        }
        hit += !ev.is_empty() as usize;
        for e in ev {
            match bnd.iter_mut().find(|(k, _)| *k == e) {
                Some(x) => x.1 += 1,
                None => bnd.push((e, 1)),
            }
        }
    }
    seen > 0 && 2 * hit >= seen && bnd.iter().map(|(_, n)| *n).max().is_some_and(|m| m >= MIN_SUPPORT)
}

/// >= half of haplotype i's reads hold an aligned block no block of the copy's other holders' reads overlaps.
pub fn foreign_coverage(reads_i: &BTreeSet<ReadId>, reads_co: &BTreeSet<ReadId>, blocks: &HashMap<ReadId, Vec<(i64, i64)>>) -> bool {
    let mut co: Vec<(i64, i64)> = reads_co.iter().filter_map(|r| blocks.get(r)).flatten().copied().collect();
    co.sort();
    if reads_co.iter().filter(|r| blocks.contains_key(r)).count() < MIN_READS || co.is_empty() {
        return false;
    }
    let free = |b: &(i64, i64)| !co.iter().any(|&(x, y)| x < b.1 && b.0 < y);
    let hit = reads_i.iter().filter_map(|r| blocks.get(r)).filter(|bs| bs.iter().any(free)).count();
    let n = reads_i.iter().filter(|r| blocks.contains_key(r)).count();
    2 * hit >= n && n > 0
}

/// One putative copy: (haplotype, its slot copy, its multi fraction, the other holders' multi fractions).
pub type Putative = (usize, String, f64, Vec<f64>);

pub fn putative_copies(
    demand: &[usize],
    matched: &Matching,
    reads_of: &BTreeMap<usize, BTreeSet<ReadId>>,
    n_piles: &HashMap<ReadId, usize>,
    foreign: Option<&dyn Fn(usize, &str, &[usize]) -> bool>,
) -> Vec<Putative> {
    let multi: BTreeMap<usize, f64> = demand
        .iter()
        .map(|&i| {
            let rs = &reads_of[&i];
            (i, rs.iter().filter(|r| n_piles.get(r).copied().unwrap_or(0) >= 2).count() as f64 / rs.len().max(1) as f64)
        })
        .collect();
    let mut out = Vec::new();
    for &i in demand {
        let x = matched.get(&i).cloned().flatten();
        let co: Vec<usize> = demand.iter().copied().filter(|&j| j != i && matched.get(&j).cloned().flatten() == x).collect();
        let Some(x) = x else { continue };
        if co.is_empty() {
            continue;
        }
        if (2.0 * multi[&i] >= 1.0 && co.iter().all(|j| 2.0 * multi[j] < 1.0)) || foreign.is_some_and(|f| f(i, &x, &co)) {
            out.push((i, x, multi[&i], co.iter().map(|j| multi[j]).collect()));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::copy_phasing::pile::Rec;
    use std::sync::Arc;

    fn usage4() -> Usage {
        (0..4).map(|i| (i, (100, 900, BTreeSet::from([(200, 300), (400, 600)])))).collect()
    }
    fn s2(a: &str, b: &str) -> (Vec<u8>, Vec<u8>) {
        (a.as_bytes().to_vec(), b.as_bytes().to_vec())
    }

    #[test]
    fn foreign_splicing_counts_non_canonical_majority_junctions() {
        let u = usage4();
        let all: BTreeSet<ReadId> = (0..4).collect();
        assert!(!foreign_splicing(&all, &u, |_| s2("GT", "AG")));
        assert!(!foreign_splicing(&all, &u, |j| if j == (200, 300) { s2("GT", "AG") } else { s2("AA", "TT") }));
        assert!(foreign_splicing(&all, &u, |_| s2("AA", "TT")));
        let mut minor = u.clone();
        for i in 10..18 {
            minor.insert(i, (100, 900, BTreeSet::from([(200, 300)]))); // (400, 600) now a minor isoform
        }
        let keys: BTreeSet<ReadId> = minor.keys().copied().collect();
        assert!(!foreign_splicing(&keys, &minor, |_| s2("AA", "TT")));
    }

    #[test]
    fn putative_needs_a_contrast_with_the_copys_own_allele() {
        // reads a1 a2 a3 = 0 1 2, e1 e2 e3 = 3 4 5, z1 z2 = 6 7
        let reads_of = BTreeMap::from([(0, BTreeSet::from([0, 1, 2])), (1, BTreeSet::from([3, 4, 5])), (2, BTreeSet::from([6, 7]))]);
        let n_piles = HashMap::from([(0, 1), (1, 1), (2, 2), (3, 2), (4, 2), (5, 1), (6, 2), (7, 2)]);
        let matched: Matching = BTreeMap::from([(0, Some("X".into())), (1, Some("X".into())), (2, Some("Y".into()))]);
        let got = putative_copies(&[0, 1, 2], &matched, &reads_of, &n_piles, None);
        assert_eq!(got.iter().map(|(i, c, _, _)| (*i, c.as_str())).collect::<Vec<_>>(), vec![(1, "X")]);
        assert_eq!(got[0].2, 2.0 / 3.0);
        assert_eq!(got[0].3, vec![1.0 / 3.0]);
    }

    #[test]
    fn putative_structure_flags_a_one_survivor_family() {
        let reads_of = BTreeMap::from([(0, BTreeSet::from([0, 1])), (1, BTreeSet::from([2, 3]))]);
        let n_piles = HashMap::from([(0, 1), (1, 1), (2, 1), (3, 1)]);
        let matched: Matching = BTreeMap::from([(0, Some("X".into())), (1, Some("X".into()))]);
        assert!(putative_copies(&[0, 1], &matched, &reads_of, &n_piles, None).is_empty());
        let f = |i: usize, _: &str, _: &[usize]| i == 1;
        let got = putative_copies(&[0, 1], &matched, &reads_of, &n_piles, Some(&f));
        assert_eq!(got.iter().map(|(i, c, _, _)| (*i, c.as_str())).collect::<Vec<_>>(), vec![(1, "X")]);
    }

    fn rec(read: ReadId, cigar: &[(u8, u32)], start: i64, end: i64, seq: &str) -> Rec {
        Rec {
            read,
            chrom: Arc::from("c"),
            start,
            end,
            as_: Some(0),
            nm: None,
            reverse: false,
            secondary: false,
            supplementary: false,
            cigar: cigar.to_vec(),
            seq: seq.as_bytes().to_vec(),
        }
    }

    #[test]
    fn foreign_clipping() {
        let g80c900 = format!("{}{}", "G".repeat(80), "C".repeat(900));
        let clip: Pile = (0..8).map(|k| (k, rec(k, &[(OP_S, 80), (OP_M, 900)], 1000 + k as i64, 1900, &g80c900))).collect();
        let clean: Pile = (8..16).map(|k| (k, rec(k, &[(OP_M, 980)], 1000, 1980, &"G".repeat(980)))).collect();
        let c900a80 = format!("{}{}", "C".repeat(900), "A".repeat(80));
        let tails: Pile = (16..24).map(|k| (k, rec(k, &[(OP_M, 900), (OP_S, 80)], 1000, 1900, &c900a80))).collect();
        let keys = |p: &Pile| p.keys().copied().collect::<BTreeSet<_>>();
        assert!(super::foreign_clipping(&keys(&clip), &clip));
        assert!(!super::foreign_clipping(&keys(&clean), &clean));
        assert!(!super::foreign_clipping(&keys(&tails), &tails));
        let mut mixed: Pile = clip.iter().take(3).map(|(k, v)| (*k, v.clone())).collect();
        mixed.extend(clean.iter().map(|(k, v)| (*k, v.clone())));
        assert!(!super::foreign_clipping(&keys(&mixed), &mixed));
        // a record without its sequence (a secondary) counts as not a tail
        let noseq: Pile = (30..38).map(|k| (k, rec(k, &[(OP_M, 900), (OP_S, 80)], 1000, 1900, ""))).collect();
        assert!(super::foreign_clipping(&keys(&noseq), &noseq));
    }

    #[test]
    fn foreign_coverage_needs_free_blocks_on_half_the_reads() {
        let blocks = HashMap::from([(0, vec![(0, 100)]), (1, vec![(50, 150)]), (2, vec![(500, 600)]), (3, vec![(10, 20)])]);
        let co = BTreeSet::from([0, 1]);
        assert!(super::foreign_coverage(&BTreeSet::from([2, 3]), &co, &blocks)); // 1 of 2 free: half
        assert!(!super::foreign_coverage(&BTreeSet::from([3]), &co, &blocks));
        assert!(!super::foreign_coverage(&BTreeSet::from([2]), &BTreeSet::from([0]), &blocks)); // < MIN_READS holders
    }

    #[test]
    fn junction_usage_reads_introns_off_the_cigar() {
        let p: Pile = [(0, rec(0, &[(OP_S, 5), (OP_M, 10), (OP_N, 100), (OP_M, 5), (OP_D, 2), (OP_N, 50), (OP_M, 3)], 1000, 1170, ""))].into_iter().collect();
        assert_eq!(junction_usage(&p)[&0], (1000, 1170, BTreeSet::from([(1010, 1110), (1117, 1167)])));
    }
}

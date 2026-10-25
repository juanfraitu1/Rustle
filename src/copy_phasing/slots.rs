//! Allele slots, consensus and explanation (o3port.py `landing_loci`, `b_matching`, `forced_edges`, `consensus`, `pile_consensus`,
//! `divergence`, `explanation`, `private_detail`, `private_alleles`).

use super::stats::{counter, most_common};
use super::{Hap, Obs, ReadId, MIN_READS};
use indexmap::IndexMap;
use std::collections::{BTreeMap, BTreeSet};

/// A haplotype after merging the piles it was found in (the Python's `merged` rows).
#[derive(Clone, Debug, PartialEq)]
pub struct Merged {
    pub reads: BTreeSet<ReadId>,
    pub piles: BTreeSet<String>,
    pub div: f64,
    /// {pile copy: the haplotype's alleles there}, in the order the piles joined.
    pub hap: IndexMap<String, Hap>,
}

/// {haplotype: its copy options}.
pub type Options = BTreeMap<usize, Vec<String>>;
/// {haplotype: its slot copy, or None}.
pub type Matching = BTreeMap<usize, Option<String>>;

/// A parent of the DADA2 rule: (name, {column: base, None where the Python stores None}).
pub type Parent = (String, BTreeMap<i64, Option<u8>>);

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Explain {
    Single(String),
    Join(String, String),
}

/// `landing_loci`: records `(chrom, start, end, read, near_best)` overlapping no copy span, merged into loci; those with
/// >= `min_reads` reads as (chrom, start, end, reads, near-best reads, weaker-only reads).
pub fn landing_loci(
    records: &[(String, i64, i64, ReadId, bool)],
    spans: &[(String, i64, i64)],
    min_reads: usize,
) -> Vec<(String, i64, i64, usize, usize, usize)> {
    let mut out: Vec<&(String, i64, i64, ReadId, bool)> =
        records.iter().filter(|r| !spans.iter().any(|(c, s, e)| r.0 == *c && r.1 < *e && *s < r.2)).collect();
    out.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut loci: Vec<(String, i64, i64, BTreeSet<ReadId>, BTreeSet<ReadId>)> = Vec::new();
    for (ch, s, e, n, near) in out {
        match loci.last_mut() {
            Some(l) if l.0 == *ch && *s < l.2 => l.2 = l.2.max(*e),
            _ => loci.push((ch.clone(), *s, *e, BTreeSet::new(), BTreeSet::new())),
        }
        let l = loci.last_mut().unwrap();
        l.3.insert(*n);
        if *near {
            l.4.insert(*n);
        }
    }
    loci.into_iter()
        .filter(|l| l.3.len() >= min_reads)
        .map(|(ch, s, e, r, nb)| {
            let weaker = r.difference(&nb).count();
            (ch, s, e, r.len(), nb.len(), weaker)
        })
        .collect()
}

/// Maximum b-matching by augmenting paths, haplotypes in `order` (default: fewest options first, then index).
pub fn b_matching(options: &Options, capacity: &BTreeMap<String, usize>, order: Option<&[usize]>) -> Matching {
    struct St<'a> {
        options: &'a Options,
        capacity: &'a BTreeMap<String, usize>,
        load: BTreeMap<String, Vec<usize>>,
        matched: Matching,
    }
    fn aug(st: &mut St, h: usize, avoid: Option<&str>, seen: &mut BTreeSet<String>) -> bool {
        let opts = st.options[&h].clone();
        for c in &opts {
            if Some(c.as_str()) == avoid || seen.contains(c) {
                continue;
            }
            seen.insert(c.clone());
            let load = st.load.entry(c.clone()).or_default();
            if load.len() < st.capacity[c] {
                load.push(h);
                st.matched.insert(h, Some(c.clone()));
                return true;
            }
            for other in load.clone() {
                if aug(st, other, Some(c), seen) {
                    let load = st.load.get_mut(c).unwrap();
                    let k = load.iter().position(|&x| x == other).unwrap();
                    load.remove(k);
                    load.push(h);
                    st.matched.insert(h, Some(c.clone()));
                    return true;
                }
            }
        }
        false
    }
    let default: Vec<usize>;
    let order = match order {
        Some(o) => o,
        None => {
            let mut hs: Vec<usize> = options.keys().copied().collect();
            hs.sort_by_key(|h| (options[h].len(), *h));
            default = hs;
            &default
        }
    };
    let mut st = St { options, capacity, load: BTreeMap::new(), matched: Matching::new() };
    for &h in order {
        if !aug(&mut st, h, None, &mut BTreeSet::new()) {
            st.matched.insert(h, None);
        }
    }
    st.matched
}

/// The matched edges present in every maximum b-matching.
pub fn forced_edges(options: &Options, capacity: &BTreeMap<String, usize>, matched: &Matching, order: &[usize]) -> BTreeMap<usize, String> {
    let size = matched.values().filter(|v| v.is_some()).count();
    let mut out = BTreeMap::new();
    for (&i, x) in matched {
        let Some(x) = x else { continue };
        let mut opt = options.clone();
        opt.insert(i, options[&i].iter().filter(|c| *c != x).cloned().collect());
        let alt = b_matching(&opt, capacity, Some(order));
        if alt.values().filter(|v| v.is_some()).count() < size {
            out.insert(i, x.clone());
        }
    }
    out
}

/// The majority consensus of one copy's alignments (`steps` of `observations`), positions in genome order.
pub fn consensus<'a>(alignments: impl IntoIterator<Item = &'a Vec<(i64, u8, Vec<u8>)>>) -> Vec<u8> {
    let mut at: BTreeMap<i64, Vec<(u8, &[u8])>> = BTreeMap::new();
    for steps in alignments {
        for (pos, base, ins) in steps {
            at.entry(*pos).or_default().push((*base, ins));
        }
    }
    let mut out = Vec::new();
    for xs in at.values() {
        if xs.len() < MIN_READS {
            continue;
        }
        let (b, _) = most_common(&counter(xs.iter().map(|x| x.0))).unwrap();
        if b != b'-' {
            out.push(b);
        }
        let mut ins: Vec<(&[u8], usize)> = Vec::new();
        for &(_, i) in xs.iter().filter(|x| !x.1.is_empty()) {
            match ins.iter_mut().find(|(s, _)| *s == i) {
                Some(e) => e.1 += 1,
                None => ins.push((i, 1)),
            }
        }
        if let Some((sq, n)) = most_common(&ins) {
            let total: usize = ins.iter().map(|(_, k)| k).sum();
            if 2 * total > xs.len() && n >= MIN_READS {
                out.extend_from_slice(sq);
            }
        }
    }
    out
}

/// `pile_consensus`: the most common base of every observed column.
pub fn pile_consensus(obs: &Obs) -> IndexMap<i64, u8> {
    obs.iter()
        .filter(|(_, xs)| !xs.is_empty())
        .map(|(&c, xs)| (c, most_common(&counter(xs.iter().map(|x| x.base))).unwrap().0))
        .collect()
}

/// Fraction of a haplotype's columns where its allele is not the reference base (columns with an N reference skipped).
pub fn divergence(h: &Hap, reference: &BTreeMap<i64, u8>) -> f64 {
    let cols: Vec<i64> = h.keys().copied().filter(|c| reference.get(c).copied().unwrap_or(b'N') != b'N').collect();
    if cols.is_empty() {
        return 0.0;
    }
    cols.iter().filter(|c| h[*c] != reference[*c]).count() as f64 / cols.len() as f64
}

/// DADA2's bimera rule on a haplotype's columns.
pub fn explanation(hp: &Hap, parents: &[Parent]) -> Option<Explain> {
    let mut cols: Vec<i64> = hp.keys().copied().collect();
    cols.sort();
    if cols.is_empty() {
        return None;
    }
    let (mut pre, mut suf): (Vec<(usize, &str)>, Vec<(usize, &str)>) = (Vec::new(), Vec::new());
    for (name, pb) in parents {
        let m: Vec<bool> = cols.iter().map(|c| pb.get(c).copied().flatten() == Some(hp[c])).collect();
        if m.iter().all(|&x| x) {
            return Some(Explain::Single(name.clone()));
        }
        let lp = m.iter().position(|&x| !x).unwrap_or(m.len());
        let ls = m.iter().rev().position(|&x| !x).unwrap_or(m.len());
        pre.push((lp, name));
        suf.push((ls, name));
    }
    pre.sort_by(|a, b| b.cmp(a));
    suf.sort_by(|a, b| b.cmp(a));
    for &(lp, a) in pre.iter().take(2) {
        for &(ls, b) in suf.iter().take(2) {
            if a != b && lp + ls >= cols.len() {
                return Some(Explain::Join(a.to_string(), b.to_string()));
            }
        }
    }
    None
}

/// (pile, [(column, allele)]): the alleles of unmatched haplotype i no slot-holding haplotype of its option copies carries.
pub fn private_detail(i: usize, merged: &[Merged], matched: &Matching, options: &Options) -> (Option<String>, Vec<(i64, u8)>) {
    let mut best: (Option<String>, Vec<(i64, u8)>) = (None, Vec::new());
    for (cid, hp) in &merged[i].hap {
        let holders: Vec<&Hap> = matched
            .iter()
            .filter(|(&j, m)| m.is_some() && j != i && merged[j].hap.contains_key(cid) && options[&i].contains(cid))
            .map(|(&j, _)| &merged[j].hap[cid])
            .collect();
        let pv: Vec<(i64, u8)> = hp.iter().filter(|(c, b)| holders.iter().all(|h| h.get(*c) != Some(*b))).map(|(&c, &b)| (c, b)).collect();
        if best.0.is_none() || pv.len() > best.1.len() {
            best = (Some(cid.clone()), pv);
        }
    }
    best
}

pub fn private_alleles(i: usize, merged: &[Merged], matched: &Matching, options: &Options) -> usize {
    private_detail(i, merged, matched, options).1.len()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn opts(xs: &[(usize, &[&str])]) -> Options {
        xs.iter().map(|(h, cs)| (*h, cs.iter().map(|c| c.to_string()).collect())).collect()
    }
    fn cap(xs: &[(&str, usize)]) -> BTreeMap<String, usize> {
        xs.iter().map(|(c, n)| (c.to_string(), *n)).collect()
    }
    fn hap(xs: &[(i64, u8)]) -> Hap {
        xs.iter().copied().collect()
    }

    #[test]
    fn b_matching_leaves_the_last_in_order_unmatched() {
        let m = b_matching(&opts(&[(0, &["c1"]), (1, &["c1"]), (2, &["c1"])]), &cap(&[("c1", 2)]), Some(&[0, 1, 2]));
        assert_eq!(m[&0].as_deref(), Some("c1"));
        assert_eq!(m[&1].as_deref(), Some("c1"));
        assert_eq!(m[&2], None);
        // augmenting: haplotype 0 can move to c2, so 2 fits in c1
        let m = b_matching(&opts(&[(0, &["c1", "c2"]), (1, &["c1"]), (2, &["c1"])]), &cap(&[("c1", 2), ("c2", 2)]), Some(&[0, 1, 2]));
        assert!(m.values().all(Option::is_some));
        // the default order: fewest options first
        let m = b_matching(&opts(&[(0, &["c1", "c2"]), (1, &["c1"])]), &cap(&[("c1", 1), ("c2", 1)]), None);
        assert_eq!((m[&0].as_deref(), m[&1].as_deref()), (Some("c2"), Some("c1")));
    }

    #[test]
    fn forced_edges() {
        let opt = opts(&[(0, &["X"]), (1, &["X", "Y"])]);
        let c = cap(&[("X", 1), ("Y", 1)]);
        let m = b_matching(&opt, &c, Some(&[0, 1]));
        assert_eq!(super::forced_edges(&opt, &c, &m, &[0, 1]), BTreeMap::from([(0, "X".to_string()), (1, "Y".to_string())]));
        let c = cap(&[("X", 2), ("Y", 2)]);
        let m = b_matching(&opt, &c, Some(&[0, 1]));
        assert_eq!(super::forced_edges(&opt, &c, &m, &[0, 1]), BTreeMap::from([(0, "X".to_string())]));
    }

    #[test]
    fn consensus_majority_with_deletion_and_insertion() {
        let al = vec![
            vec![(10, b'A', vec![]), (11, b'C', b"GG".to_vec()), (12, b'-', vec![])],
            vec![(10, b'A', vec![]), (11, b'C', b"GG".to_vec()), (12, b'-', vec![])],
            vec![(10, b'T', vec![]), (11, b'C', vec![]), (12, b'G', vec![])],
        ];
        assert_eq!(consensus(&al), b"ACGG");
    }

    #[test]
    fn divergence() {
        assert_eq!(super::divergence(&hap(&[(1, b'A'), (2, b'G')]), &BTreeMap::from([(1, b'A'), (2, b'T')])), 0.5);
        assert_eq!(super::divergence(&Hap::new(), &BTreeMap::new()), 0.0);
        assert_eq!(super::divergence(&hap(&[(1, b'A'), (2, b'G')]), &BTreeMap::from([(1, b'N'), (2, b'T')])), 1.0);
    }

    #[test]
    fn private_alleles_mosaic_versus_third_pattern() {
        let m = |h: Hap| Merged { reads: BTreeSet::new(), piles: BTreeSet::new(), div: 0.0, hap: IndexMap::from([("c".to_string(), h)]) };
        let merged = vec![
            m(hap(&[(1, b'A'), (2, b'C'), (3, b'G'), (4, b'T')])),
            m(hap(&[(1, b'T'), (2, b'T'), (3, b'C'), (4, b'A')])),
            m(hap(&[(1, b'A'), (2, b'C'), (3, b'C'), (4, b'A')])),
            m(hap(&[(1, b'G'), (2, b'C'), (3, b'G'), (4, b'G')])),
        ];
        let matched: Matching = BTreeMap::from([(0, Some("c".into())), (1, Some("c".into())), (2, None), (3, None)]);
        let options = opts(&[(0, &["c"]), (1, &["c"]), (2, &["c"]), (3, &["c"])]);
        assert_eq!(private_alleles(2, &merged, &matched, &options), 0);
        assert_eq!(private_alleles(3, &merged, &matched, &options), 2);
        assert_eq!(private_detail(3, &merged, &matched, &options), (Some("c".to_string()), vec![(1, b'G'), (4, b'G')]));
    }

    fn parent(name: &str, xs: &[(i64, u8)]) -> Parent {
        (name.to_string(), xs.iter().map(|&(c, b)| (c, Some(b))).collect())
    }

    #[test]
    fn explanation_dada2_rule() {
        let hp = hap(&[(1, b'A'), (2, b'C'), (3, b'G'), (4, b'T')]);
        let left = parent("p1", &[(1, b'A'), (2, b'C'), (3, b'A'), (4, b'A')]);
        let right = parent("p2", &[(1, b'G'), (2, b'G'), (3, b'G'), (4, b'T')]);
        assert_eq!(explanation(&hp, &[left.clone(), right]), Some(Explain::Join("p1".into(), "p2".into())));
        assert_eq!(explanation(&hp, &[parent("same", &[(1, b'A'), (2, b'C'), (3, b'G'), (4, b'T')])]), Some(Explain::Single("same".into())));
        assert_eq!(explanation(&hp, &[left]), None);
        assert_eq!(explanation(&hp, &[parent("x", &[(1, b'A'), (4, b'T')])]), None);
        // a None base (no allele anywhere) is a mismatch
        let mut p = parent("y", &[(1, b'A'), (2, b'C'), (3, b'G')]);
        p.1.insert(4, None);
        assert_eq!(explanation(&hp, &[p]), None);
    }

    #[test]
    fn landing_loci_merges_records_outside_the_copies() {
        let c = |s: &str| s.to_string();
        let recs = vec![
            (c("c1"), 100, 200, 1, true),
            (c("c1"), 150, 260, 2, false),
            (c("c1"), 150, 260, 3, false),
            (c("c1"), 1000, 1100, 1, false),
            (c("c2"), 0, 50, 4, true),
        ];
        assert_eq!(landing_loci(&recs, &[(c("c1"), 1000, 2000)], 3), vec![(c("c1"), 100, 260, 3, 1, 2)]);
    }

    #[test]
    fn pile_consensus_first_appearance_wins_a_tie() {
        use crate::copy_phasing::ObsX;
        let x = |read, base| ObsX { read, base, qual: 30, eps: None };
        let mut obs = Obs::new();
        obs.insert(7, vec![x(0, b'G'), x(1, b'A'), x(2, b'A'), x(3, b'G')]);
        obs.insert(3, vec![]);
        assert_eq!(pile_consensus(&obs), IndexMap::from([(7, b'G')]));
    }
}

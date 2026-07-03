//! O1 multi-copy family-definition predicate `distinct_loci`.
//!
//! Faithful Rust port of the Python reference
//! `bench/genome_family_def.py::distinct_loci` (LOCUS_OVERLAP = 0.50). This is
//! the first Python->Rust migration of the family-definition layer: it counts
//! the DISTINCT PHYSICAL LOCI among a family's member gene spans, collapsing
//! redundant/nested annotations at the same physical locus (e.g. MAGEA9 ==
//! a nested LOC copy) before the ">= 2 distinct loci" multi-copy test.
//!
//! The DNA multi-copy family predicate `|distinct_loci(C)| >= 2` is defined in
//! `bench/DNA_FAMILY_DEFINITION_FORMAL.md` (the >=2-loci CERTIFICATE); the
//! collapse step is the `LOCUS_OVERLAP` reciprocal-span rule documented there.
//!
//! Semantics preserved EXACTLY from the Python reference:
//!   (a) every member counts as its own locus unless collapsed (the Python
//!       `uf.find(i)` singleton guarantee) -- a member that overlaps no other
//!       is one locus;
//!   (b) two members collapse to one locus iff SAME contig AND overlap > 0 AND
//!       overlap >= 0.50 * len_a AND overlap >= 0.50 * len_b (RECIPROCAL);
//!   (c) collapse is transitive (union-find), so a chain A~B~C where A and C do
//!       not directly overlap enough still collapses to one locus;
//!   (d) members on different contigs never collapse.
//!
//! The count is invariant to member order and to union order, so the per-contig
//! `sort_by_key(start)` and the naive-union union-find below reproduce the
//! Python `len(uf.groups())` byte-for-byte. Verified against a golden fixture
//! generated from the shipped Python build (see `#[test] distinct_loci_parity`).

use std::collections::HashMap;

/// Two member spans on the SAME contig are the SAME physical locus iff their
/// spans reciprocally overlap by >= this fraction of BOTH lengths. Mirrors
/// `genome_family_def.LOCUS_OVERLAP`.
pub const LOCUS_OVERLAP: f64 = 0.50;

/// Number of distinct physical loci among `members`, each `(chrom, start, end)`
/// with 0-based half-open coords (matching the Python `genes` dict).
///
/// Members on the same contig whose spans reciprocally overlap by at least
/// `LOCUS_OVERLAP` of both lengths are collapsed (union-find) into one locus;
/// the return value is the count of distinct groups. See the module docs for
/// the exact semantics preserved from `genome_family_def.py::distinct_loci`.
pub fn distinct_loci(members: &[(String, i64, i64)]) -> usize {
    let mut uf = UnionFind::new(members.len());

    // Group member indices by contig (cross-contig pairs never collapse).
    let mut by_contig: HashMap<&str, Vec<usize>> = HashMap::new();
    for (i, (chrom, _, _)) in members.iter().enumerate() {
        by_contig.entry(chrom.as_str()).or_default().push(i);
    }

    for idxs in by_contig.values_mut() {
        // Start-sorted so the `sb >= ea` break is a valid overlap-prune, exactly
        // as in the Python `idxs.sort(key=lambda i: genes[i]["start"])`.
        idxs.sort_by_key(|&i| members[i].1);
        for a in 0..idxs.len() {
            let ia = idxs[a];
            let (sa, ea) = (members[ia].1, members[ia].2);
            for &ib in &idxs[a + 1..] {
                let (sb, eb) = (members[ib].1, members[ib].2);
                if sb >= ea {
                    break; // start-sorted: nothing further on this contig can overlap A
                }
                let ov = ea.min(eb) - sa.max(sb);
                // RECIPROCAL >= LOCUS_OVERLAP of both lengths (f64 to match Python).
                if ov > 0
                    && ov as f64 >= LOCUS_OVERLAP * (ea - sa) as f64
                    && ov as f64 >= LOCUS_OVERLAP * (eb - sb) as f64
                {
                    uf.union(ia, ib);
                }
            }
        }
    }
    uf.num_groups()
}

/// Minimal union-find (path-compression + naive union), a faithful port of the
/// Python `genome_family_def.UF`. Naive union (`parent[root_a] = root_b`) is
/// kept deliberately: the returned GROUP COUNT is union-order invariant, so the
/// exact merge direction is irrelevant to `distinct_loci`.
struct UnionFind {
    parent: Vec<usize>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n).collect(),
        }
    }

    fn find(&mut self, x: usize) -> usize {
        let mut root = x;
        while self.parent[root] != root {
            root = self.parent[root];
        }
        // Path compression (as in the Python `while self.p[x] != r` walk).
        let mut cur = x;
        while self.parent[cur] != root {
            let next = self.parent[cur];
            self.parent[cur] = root;
            cur = next;
        }
        root
    }

    fn union(&mut self, a: usize, b: usize) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra != rb {
            self.parent[ra] = rb;
        }
    }

    /// Number of distinct groups among all `n` members (each member is a node,
    /// matching the Python `uf.find(i)` singleton guarantee).
    fn num_groups(&mut self) -> usize {
        let n = self.parent.len();
        let mut roots = std::collections::HashSet::with_capacity(n);
        for i in 0..n {
            let r = self.find(i);
            roots.insert(r);
        }
        roots.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(chrom: &str, start: i64, end: i64) -> (String, i64, i64) {
        (chrom.to_string(), start, end)
    }

    #[test]
    fn empty_and_singleton() {
        assert_eq!(distinct_loci(&[]), 0);
        assert_eq!(distinct_loci(&[m("c1", 100, 200)]), 1);
    }

    #[test]
    fn reciprocal_and_boundary() {
        // exact 50% overlap of both -> collapse
        assert_eq!(distinct_loci(&[m("c1", 100, 200), m("c1", 150, 250)]), 1);
        // just under 50% -> distinct
        assert_eq!(distinct_loci(&[m("c1", 100, 200), m("c1", 151, 251)]), 2);
        // nested but reciprocal FAILS on the big span (10% of it) -> distinct
        assert_eq!(distinct_loci(&[m("c1", 100, 1100), m("c1", 100, 200)]), 2);
        // different contig, identical coords -> never collapse
        assert_eq!(distinct_loci(&[m("c1", 100, 200), m("c2", 100, 200)]), 2);
        // zero-length member never overlaps (ov > 0 fails) -> distinct
        assert_eq!(distinct_loci(&[m("c1", 100, 100), m("c1", 90, 110)]), 2);
        // transitive chain A~B~C (A,C don't directly overlap enough) -> one locus
        assert_eq!(
            distinct_loci(&[m("c1", 0, 100), m("c1", 30, 130), m("c1", 60, 160)]),
            1
        );
    }

    /// Byte-parity against the golden fixture emitted by
    /// `bench/gen_distinct_loci_fixture.py` from the shipped Python build. Every
    /// row's Rust `distinct_loci` MUST equal the Python golden count.
    #[test]
    fn distinct_loci_parity() {
        let fixture = include_str!("testdata/distinct_loci_fixture.tsv");
        let mut n_rows = 0usize;
        let mut n_collapse = 0usize;
        let mut n_xchrom = 0usize;
        for (i, line) in fixture.lines().enumerate() {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            assert_eq!(
                cols.len(),
                4,
                "fixture line {}: expected 4 tab-columns, got {}",
                i + 1,
                cols.len()
            );
            let label = cols[0];
            let n_members: usize = cols[1].parse().expect("n_members int");
            let golden: usize = cols[2].parse().expect("golden int");
            let members: Vec<(String, i64, i64)> = cols[3]
                .split(';')
                .map(|tok| {
                    let mut p = tok.split('|');
                    let chrom = p.next().expect("chrom").to_string();
                    let start: i64 = p.next().expect("start").parse().expect("start int");
                    let end: i64 = p.next().expect("end").parse().expect("end int");
                    assert!(p.next().is_none(), "member token has extra fields: {tok}");
                    (chrom, start, end)
                })
                .collect();
            assert_eq!(
                members.len(),
                n_members,
                "fixture line {} ({}): declared n_members {} != parsed {}",
                i + 1,
                label,
                n_members,
                members.len()
            );

            let got = distinct_loci(&members);
            assert_eq!(
                got, golden,
                "PARITY MISMATCH family '{}' (fixture line {}): rust distinct_loci={} \
                 != python golden={} (n_members={})",
                label,
                i + 1,
                got,
                golden,
                n_members
            );

            n_rows += 1;
            if golden < n_members {
                n_collapse += 1;
            }
            let chroms: std::collections::HashSet<&str> =
                members.iter().map(|(c, _, _)| c.as_str()).collect();
            if chroms.len() >= 2 {
                n_xchrom += 1;
            }
        }
        // Guard the fixture stays representative (diverse, exercises collapse).
        assert!(n_rows >= 200, "fixture too small: {n_rows} rows (< 200)");
        assert!(
            n_collapse >= 10,
            "fixture under-exercises same-locus collapse: only {n_collapse} collapse rows"
        );
        assert!(
            n_xchrom >= 10,
            "fixture under-exercises cross-chrom families: only {n_xchrom} rows"
        );
        eprintln!(
            "distinct_loci parity: {n_rows} rows OK ({n_collapse} collapse, {n_xchrom} cross-chrom)"
        );
    }
}

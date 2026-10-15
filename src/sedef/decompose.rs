//! Stage 4: decomposition of mosaic hubs.
//!
//! Refined blocks describe a whole block pair with a single identity, but
//! pericentromeric hubs are mosaic: a block pair is typically a chain of
//! constituent fragments copied from different sources. This stage splits
//! each block into its elementary collinear segments and then groups all
//! segment endpoints into duplication families (the core-duplicon analog):
//!
//! 1. **Multi-run chaining** -- repeatedly extract the best anchor run
//!    (in either strand space), emit it as an elementary segment pair,
//!    remove its anchors, and re-chain. Stops when the best run is shorter
//!    than `min_seg` bp (SDs are >= 1 kbp by definition).
//! 2. **Segment dedup** -- overlapping blocks rediscover the same cores;
//!    segment pairs are normalised (transposed duplicates merged) and
//!    near-identical pairs unioned.
//! 3. **Families** -- both endpoints of every segment become nodes; nodes
//!    on the same chromosome overlapping >= 80% are unioned (inverted
//!    copies included). Connected components are duplication families.
//!
//! Still no alignments: everything is k-mer anchoring + interval algebra.

use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::io::BufRead;

use crate::sedef::anchors::BlockAnchors;
use crate::sedef::refine::{best_run, extract_anchors};

pub struct DecomposeParams {
    pub k: usize,
    pub max_gap: u64,
    pub diag_bin: i64,
    /// Elementary segments shorter than this are not emitted.
    pub min_seg: u64,
    pub max_pairs: u64,
}

/// One refined block loaded from the stage-3 TSV (ok rows only).
pub struct InBlock {
    chr_a: String,
    a0: u64,
    a1: u64,
    chr_b: String,
    b0: u64,
    b1: u64,
}

#[derive(Debug, Clone)]
pub struct Segment {
    pub chr_a: String,
    pub a0: u64,
    pub a1: u64,
    pub chr_b: String,
    pub b0: u64,
    pub b1: u64,
    pub strand: char,
    pub anchors: usize,
    pub est_id: f64,
}

pub struct Family {
    pub id: usize,
    /// The family's segment endpoints, counted before duplicates are removed: two segments that share an endpoint count
    /// it twice. `copies.len()` is the number of distinct copies, and duplicons.tsv has one row per distinct copy. Families
    /// are ordered by `n_copies`, largest first.
    pub n_copies: usize,
    pub chrs: Vec<String>,
    pub span_bp: u64,
    /// Vendored change (2026-10-13): the family's distinct copies `(chrom, start, end)`, sorted.
    pub copies: Vec<(String, u64, u64)>,
}

pub fn load_ok_blocks(path: &std::path::Path) -> Result<Vec<InBlock>, Box<dyn Error>> {
    let f = std::fs::File::open(path)?;
    let reader = std::io::BufReader::new(f);
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let t: Vec<&str> = line.split('\t').collect();
        if t.len() < 11 || t[10] != "ok" {
            continue;
        }
        out.push(InBlock {
            chr_a: t[0].to_string(),
            a0: t[1].parse()?,
            a1: t[2].parse()?,
            chr_b: t[3].to_string(),
            b0: t[4].parse()?,
            b1: t[5].parse()?,
        });
    }
    Ok(out)
}

fn get_seq<'a>(seqs: &'a HashMap<String, Vec<u8>>, chr: &str) -> Option<&'a [u8]> {
    seqs.get(chr).map(|v| v.as_slice()).or_else(|| {
        let first = chr.split_whitespace().next()?;
        seqs.get(first).map(|v| v.as_slice())
    })
}

/// Normalise a segment pair so the two sides have a deterministic order
/// (transposed discoveries merge in dedup).
fn normalize(mut s: Segment) -> Segment {
    let a_key = (&s.chr_a, s.a0);
    let b_key = (&s.chr_b, s.b0);
    if b_key < a_key {
        std::mem::swap(&mut s.chr_a, &mut s.chr_b);
        std::mem::swap(&mut s.a0, &mut s.b0);
        std::mem::swap(&mut s.a1, &mut s.b1);
    }
    s
}

fn overlap_frac(x0: u64, x1: u64, y0: u64, y1: u64) -> f64 {
    let hi = x1.min(y1);
    let lo = x0.max(y0);
    if hi <= lo {
        return 0.0;
    }
    (hi - lo) as f64 / (x1 - x0).min(y1 - y0) as f64
}

struct Dsu {
    parent: Vec<usize>,
}

impl Dsu {
    fn new(n: usize) -> Self {
        Dsu {
            parent: (0..n).collect(),
        }
    }
    fn find(&mut self, x: usize) -> usize {
        let mut r = x;
        while self.parent[r] != r {
            r = self.parent[r];
        }
        let mut c = x;
        while self.parent[c] != r {
            let n = self.parent[c];
            self.parent[c] = r;
            c = n;
        }
        r
    }
    fn union(&mut self, a: usize, b: usize) {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra != rb {
            self.parent[ra] = rb;
        }
    }
}

/// Split one block into elementary segment pairs via repeated best-run
/// extraction in both strand spaces, using pre-extracted anchors.
fn block_segments(b: &BlockAnchors, p: &DecomposeParams) -> Vec<Segment> {
    let mut out = Vec::new();
    let (a0, a1) = (b.a0 as usize, b.a1 as usize);
    let (b0, b1) = (b.b0 as usize, b.b1 as usize);
    if a1 <= a0 || b1 <= b0 {
        return out;
    }
    let k = p.k;
    let sb_len = (b1 - b0) as i64;

    let mut remaining: Vec<(u32, u32)> = b.anchors.clone();
    loop {
        if remaining.len() < 10 {
            break;
        }
        let fwd = best_run(&remaining, k, p.max_gap, p.diag_bin, |pb| pb as i64);
        let rev = best_run(&remaining, k, p.max_gap, p.diag_bin, |pb| {
            sb_len - pb as i64 - k as i64
        });
        let pick = |r: &crate::sedef::refine::Run| r.cov_a.min(r.cov_b);
        let (strand, run) = match (fwd, rev) {
            (Some(f), Some(r)) => {
                if pick(&f) >= pick(&r) {
                    ('+', f)
                } else {
                    ('-', r)
                }
            }
            (Some(f), None) => ('+', f),
            (None, Some(r)) => ('-', r),
            (None, None) => break,
        };
        if run.pts.len() < 3 || run.span_a < p.min_seg || run.span_b < p.min_seg {
            break;
        }
        let mean_span = (run.span_a + run.span_b) as f64 / 2.0;
        let est_id =
            (100.0 * (run.cov_a + run.cov_b) as f64 / (2.0 * mean_span)).min(100.0);

        let pa_min = run.pts.first().unwrap().0 as u64;
        let pa_max = run.pts.last().unwrap().0 as u64 + k as u64;
        let pb_min = run.pts.iter().map(|&(_, pb)| pb as u64).min().unwrap();
        let pb_max = run.pts.iter().map(|&(_, pb)| pb as u64).max().unwrap() + k as u64;
        let (rb0, rb1) = if strand == '-' {
            (b0 as u64 + sb_len as u64 - pb_max, b0 as u64 + sb_len as u64 - pb_min)
        } else {
            (b0 as u64 + pb_min, b0 as u64 + pb_max)
        };
        out.push(Segment {
            chr_a: b.chr_a.clone(),
            a0: a0 as u64 + pa_min,
            a1: a0 as u64 + pa_max,
            chr_b: b.chr_b.clone(),
            b0: rb0,
            b1: rb1,
            strand,
            anchors: run.pts.len(),
            est_id,
        });

        // Remove the consumed anchors (map pb back to original coords).
        let consumed: rustc_hash::FxHashSet<(u32, u32)> = run
            .pts
            .iter()
            .map(|&(pa, pb)| {
                let pb_orig = if strand == '-' {
                    (sb_len - pb as i64 - k as i64) as u32
                } else {
                    pb
                };
                (pa, pb_orig)
            })
            .collect();
        remaining.retain(|x| !consumed.contains(x));
    }
    out
}

/// Merge near-identical segment pairs (same strand, both sides >= 80%
/// overlap). Returns representative segments. Group representatives are
/// chosen deterministically (max anchors, then lexicographically smallest
/// coordinates) so reruns are bit-identical.
fn dedup_segments(mut segs: Vec<Segment>) -> Vec<Segment> {
    let n = segs.len();
    if n <= 1 {
        return segs;
    }
    let mut dsu = Dsu::new(n);
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        (&segs[i].chr_a, segs[i].a0, &segs[i].chr_b, segs[i].b0)
            .cmp(&(&segs[j].chr_a, segs[j].a0, &segs[j].chr_b, segs[j].b0))
    });
    let mut active: Vec<usize> = Vec::new();
    for &i in &order {
        active.retain(|&j| segs[j].a1 > segs[i].a0);
        for &j in &active {
            if segs[i].strand == segs[j].strand
                && segs[i].chr_a == segs[j].chr_a
                && segs[i].chr_b == segs[j].chr_b
                && overlap_frac(segs[i].a0, segs[i].a1, segs[j].a0, segs[j].a1) >= 0.8
                && overlap_frac(segs[i].b0, segs[i].b1, segs[j].b0, segs[j].b1) >= 0.8
            {
                dsu.union(i, j);
            }
        }
        active.push(i);
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        groups.entry(dsu.find(i)).or_default().push(i);
    }
    let mut out: Vec<Segment> = groups
        .into_values()
        .map(|members| {
            // Deterministic representative: max anchors, then smallest coords.
            let mut best = members[0];
            for &m in &members[1..] {
                let (a, b) = (&segs[m], &segs[best]);
                let better = (a.anchors, 0u8) > (b.anchors, 0u8)
                    || (a.anchors == b.anchors
                        && (&a.chr_a, a.a0, &a.chr_b, a.b0)
                            < (&b.chr_a, b.a0, &b.chr_b, b.b0));
                if better {
                    best = m;
                }
            }
            let mut rep = segs[best].clone();
            for &m in &members {
                rep.anchors += segs[m].anchors;
                if segs[m].est_id > rep.est_id {
                    rep.est_id = segs[m].est_id;
                }
            }
            rep
        })
        .collect();
    out.sort_by(|x, y| {
        y.anchors.cmp(&x.anchors).then_with(|| {
            (&x.chr_a, x.a0, &x.chr_b, x.b0, x.a1, x.b1, x.strand)
                .cmp(&(&y.chr_a, y.a0, &y.chr_b, y.b0, y.a1, y.b1, y.strand))
        })
    });
    out
}

/// Group segment endpoints into duplication families: the two endpoints of a segment are one family, and nodes on
/// the same chromosome overlapping >= 80% (of the shorter; any strand) are one family. Connected components are
/// duplication families.
///
/// Vendored changes (2026-10-13): (1) the sandbox never joined a segment's two endpoints, so its families were
/// clusters of overlapping endpoints with every copy of a duplication in a family of its own; (2) `span_bp` sums the
/// covered bases per chromosome (the sandbox merged intervals of different chromosomes); (3) the distinct copies are
/// kept (`Family::copies`) and break the remaining ties of the order, so ids no longer follow HashMap order.
fn families(segs: &[Segment]) -> Vec<Family> {
    #[derive(Clone)]
    struct Node {
        chr: String,
        s: u64,
        e: u64,
    }
    let nodes: Vec<Node> = segs
        .iter()
        .flat_map(|s| {
            [
                Node { chr: s.chr_a.clone(), s: s.a0, e: s.a1 },
                Node { chr: s.chr_b.clone(), s: s.b0, e: s.b1 },
            ]
        })
        .collect();
    let n = nodes.len();
    let mut dsu = Dsu::new(n);
    for i in 0..segs.len() {
        dsu.union(2 * i, 2 * i + 1);
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| (&nodes[i].chr, nodes[i].s).cmp(&(&nodes[j].chr, nodes[j].s)));
    let mut active: Vec<usize> = Vec::new();
    for &i in &order {
        active.retain(|&j| nodes[j].e > nodes[i].s);
        for &j in &active {
            if nodes[i].chr == nodes[j].chr
                && overlap_frac(nodes[i].s, nodes[i].e, nodes[j].s, nodes[j].e) >= 0.8
            {
                dsu.union(i, j);
            }
        }
        active.push(i);
    }
    let mut agg: HashMap<usize, Vec<(String, u64, u64)>> = HashMap::new();
    for (i, nd) in nodes.iter().enumerate() {
        agg.entry(dsu.find(i)).or_default().push((nd.chr.clone(), nd.s, nd.e));
    }
    let mut out: Vec<Family> = agg
        .into_values()
        .map(|of_root| {
            let n_copies = of_root.len();
            let mut copies = of_root;
            copies.sort();
            copies.dedup();
            let mut chrs: Vec<String> = copies.iter().map(|c| c.0.clone()).collect();
            chrs.dedup();
            let span_bp = chrs
                .iter()
                .map(|c| {
                    let iv: Vec<(u64, u64)> = copies.iter().filter(|x| &x.0 == c).map(|x| (x.1, x.2)).collect();
                    crate::sedef::refine::covered_len(&iv)
                })
                .sum();
            Family { id: 0, n_copies, chrs, span_bp, copies }
        })
        .collect();
    out.sort_by(|x, y| {
        y.n_copies
            .cmp(&x.n_copies)
            .then_with(|| y.span_bp.cmp(&x.span_bp))
            .then_with(|| x.chrs.cmp(&y.chrs))
            .then_with(|| x.copies.cmp(&y.copies))
    });
    for (i, f) in out.iter_mut().enumerate() {
        f.id = i;
    }
    out
}

pub struct DecomposeStats {
    pub n_blocks: usize,
    pub n_segments: usize,
    pub n_deduped: usize,
    pub n_families: usize,
}

/// Build anchor records the expensive way (from sequences) -- used when no
/// anchors file is available.
pub fn anchor_from_fasta(
    blocks: Vec<InBlock>,
    seqs: &HashMap<String, Vec<u8>>,
    k: usize,
    max_pairs: u64,
) -> Vec<BlockAnchors> {
    let out: Vec<BlockAnchors> = blocks
        .par_iter()
        .map(|b| {
            let anchors = match (
                get_seq(seqs, &b.chr_a),
                get_seq(seqs, &b.chr_b),
            ) {
                (Some(sa), Some(sb)) => {
                    let (a0, a1) = ((b.a0 as usize).min(sa.len()), (b.a1 as usize).min(sa.len()));
                    let (b0, b1) = ((b.b0 as usize).min(sb.len()), (b.b1 as usize).min(sb.len()));
                    if a1 > a0 && b1 > b0 {
                        extract_anchors(&sa[a0..a1], &sb[b0..b1], k, max_pairs)
                    } else {
                        Vec::new()
                    }
                }
                _ => Vec::new(),
            };
            BlockAnchors {
                chr_a: b.chr_a.clone(),
                a0: b.a0,
                a1: b.a1,
                chr_b: b.chr_b.clone(),
                b0: b.b0,
                b1: b.b1,
                status: 1,
                anchors,
            }
        })
        .collect();
    out
}

pub fn decompose(
    blocks: Vec<BlockAnchors>,
    p: &DecomposeParams,
) -> (Vec<Segment>, Vec<Family>, DecomposeStats) {
    let n_blocks = blocks.len();
    let raw: Vec<Segment> = blocks
        .into_par_iter()
        .filter(|b| b.status == 1)
        .flat_map(|b| block_segments(&b, p))
        .map(normalize)
        .collect();
    let n_segments = raw.len();
    let deduped = dedup_segments(raw);
    let n_deduped = deduped.len();
    let fams = families(&deduped);
    let n_families = fams.len();
    (
        deduped,
        fams,
        DecomposeStats {
            n_blocks,
            n_segments,
            n_deduped,
            n_families,
        },
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seg(ca: &str, a0: u64, a1: u64, cb: &str, b0: u64, b1: u64, strand: char) -> Segment {
        Segment { chr_a: ca.into(), a0, a1, chr_b: cb.into(), b0, b1, strand, anchors: 100, est_id: 99.0 }
    }

    #[test]
    fn a_segment_joins_its_two_copies_and_overlap_joins_further_copies() {
        // S1 chrA:1000-3000 ~ chrB:5000-7000 (inverted); S2 chrB:5100-7000 ~ chrC:100-2000 (its chrB side overlaps
        // S1's by 1900 / 1900 >= 0.8); S3 chrA:50000-52000 ~ chrA:90000-92000 (a family of its own)
        let fams = families(&[
            seg("chrA", 1000, 3000, "chrB", 5000, 7000, '-'),
            seg("chrB", 5100, 7000, "chrC", 100, 2000, '+'),
            seg("chrA", 50000, 52000, "chrA", 90000, 92000, '+'),
        ]);
        assert_eq!(fams.len(), 2);
        assert_eq!((fams[0].id, fams[0].n_copies), (0, 4));
        assert_eq!(fams[0].chrs, vec!["chrA".to_string(), "chrB".to_string(), "chrC".to_string()]);
        let want: Vec<(String, u64, u64)> = vec![
            ("chrA".to_string(), 1000, 3000),
            ("chrB".to_string(), 5000, 7000),
            ("chrB".to_string(), 5100, 7000),
            ("chrC".to_string(), 100, 2000),
        ];
        assert_eq!(fams[0].copies, want);
        // span per chromosome: chrA 2000 + chrB 2000 (5000-7000 covers 5100-7000) + chrC 1900
        assert_eq!(fams[0].span_bp, 5900);
        assert_eq!((fams[1].id, fams[1].n_copies, fams[1].span_bp), (1, 2, 4000));
    }

    #[test]
    fn two_segments_sharing_an_endpoint_count_it_twice_but_keep_it_once() {
        // S1 chrA:1000-3000 ~ chrB:5000-7000; S2 chrB:5000-7000 ~ chrC:100-2100. The chrB interval is an endpoint of
        // both segments: two nodes (n_copies counts four endpoints), one distinct copy (three copies)
        let fams = families(&[
            seg("chrA", 1000, 3000, "chrB", 5000, 7000, '+'),
            seg("chrB", 5000, 7000, "chrC", 100, 2100, '+'),
        ]);
        assert_eq!(fams.len(), 1);
        assert_eq!((fams[0].n_copies, fams[0].copies.len()), (4, 3));
        assert_eq!(
            fams[0].copies,
            vec![("chrA".to_string(), 1000, 3000), ("chrB".to_string(), 5000, 7000), ("chrC".to_string(), 100, 2100)]
        );
        assert_eq!(fams[0].span_bp, 6000);
    }

    #[test]
    fn equal_families_are_ordered_by_their_copies() {
        // same copy count, span and chromosome: the family whose first copy comes first is F0, whatever the
        // HashMap's order (each call builds new maps with new random seeds)
        let segs = [seg("chrA", 9000, 11000, "chrA", 30000, 32000, '+'), seg("chrA", 1000, 3000, "chrA", 20000, 22000, '+')];
        for _ in 0..20 {
            let f = families(&segs);
            assert_eq!(f.len(), 2);
            assert_eq!(f[0].copies[0], ("chrA".to_string(), 1000, 3000));
            assert_eq!(f[1].copies[0], ("chrA".to_string(), 9000, 11000));
            assert_eq!(f.iter().map(|x| x.id).collect::<Vec<_>>(), vec![0, 1]);
        }
    }

    #[test]
    fn dedup_orders_tied_representatives_by_their_ends_and_strand() {
        // the same anchors and starts on opposite strands stay two segments (dedup joins one strand only); the
        // sandbox left their order to HashMap order
        for _ in 0..20 {
            let out = dedup_segments(vec![
                seg("chrA", 0, 2000, "chrB", 0, 2000, '-'),
                seg("chrA", 0, 2000, "chrB", 0, 2000, '+'),
            ]);
            assert_eq!(out.iter().map(|s| s.strand).collect::<Vec<_>>(), vec!['+', '-']);
        }
    }
}

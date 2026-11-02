//! Per-sample read counts on a final assembly (docs/superpowers/specs/2026-11-02-multi-sample-design.md). A primary mapped read
//! is compatible with a transcript when its intron chain is a contiguous block of the transcript's chain and its first / last
//! exon ends lie within END_TOL_BP outside the matching transcript exons (a single-exon read: inside one transcript exon
//! +/- END_TOL_BP). Strand is not a filter. A read compatible with k transcripts adds 1/k to each.

pub const END_TOL_BP: i64 = 10;

pub struct Tx {
    pub id: String,
    pub gene: String,
    pub chrom: String,
    /// 0-based half-open, sorted
    pub exons: Vec<(i64, i64)>,
}

pub fn compatible(read: &[(i64, i64)], tx: &Tx) -> bool {
    let (t, m) = (&tx.exons, read.len());
    if m == 0 || t.len() < m {
        return false;
    }
    let (r0, r1) = (read[0].0, read[m - 1].1);
    if m == 1 {
        return t.iter().any(|&(a, b)| r0 >= a - END_TOL_BP && r1 <= b + END_TOL_BP);
    }
    // the read's first intron starts at read[0].1; find the transcript exon ending there
    (0..=t.len() - m).any(|k| {
        (0..m - 1).all(|i| t[k + i].1 == read[i].1 && t[k + i + 1].0 == read[i + 1].0)
            && r0 >= t[k].0 - END_TOL_BP
            && r1 <= t[k + m - 1].1 + END_TOL_BP
    })
}

/// Transcripts spanning more than LONG_BP are scanned separately, so one giant model does not widen every read's window.
const LONG_BP: i64 = 100_000;

/// Per-transcript counts (same order as `txs`), reads seen, reads compatible. `reads`: (chrom, exons), streamed.
pub fn count(txs: &[Tx], reads: impl Iterator<Item = (String, Vec<(i64, i64)>)>) -> (Vec<f64>, u64, u64) {
    use std::collections::HashMap;
    // per chrom: short transcripts (start, end, index) sorted by start, and the long ones
    let mut by: HashMap<&str, (Vec<(i64, i64, usize)>, Vec<(i64, i64, usize)>)> = HashMap::new();
    for (i, t) in txs.iter().enumerate() {
        let (s, e) = (t.exons[0].0, t.exons[t.exons.len() - 1].1);
        let v = by.entry(t.chrom.as_str()).or_default();
        if e - s > LONG_BP { v.1.push((s, e, i)) } else { v.0.push((s, e, i)) }
    }
    for v in by.values_mut() {
        v.0.sort_unstable();
    }
    let mut out = vec![0.0; txs.len()];
    let (mut seen, mut comp) = (0u64, 0u64);
    let mut hits = Vec::new();
    for (chrom, ex) in reads {
        seen += 1;
        let Some((short, long)) = by.get(chrom.as_str()) else { continue };
        let (r0, r1) = (ex[0].0, ex[ex.len() - 1].1);
        hits.clear();
        let lo = short.partition_point(|x| x.0 < r0 - LONG_BP - END_TOL_BP);
        for &(s, e, i) in short[lo..].iter().take_while(|x| x.0 <= r1 + END_TOL_BP).chain(long.iter()) {
            if s <= r1 + END_TOL_BP && e >= r0 - END_TOL_BP && compatible(&ex, &txs[i]) {
                hits.push(i);
            }
        }
        if !hits.is_empty() {
            comp += 1;
            let w = 1.0 / hits.len() as f64;
            for &i in &hits {
                out[i] += w;
            }
        }
    }
    (out, seen, comp)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tx(e: &[(i64, i64)]) -> Tx {
        Tx { id: "t".into(), gene: "g".into(), chrom: "c".into(), exons: e.to_vec() }
    }

    #[test]
    fn chain_must_be_a_contiguous_block() {
        let r = [(100, 200), (300, 400)];
        assert!(compatible(&r, &tx(&[(50, 200), (300, 500)])));
        assert!(compatible(&r, &tx(&[(50, 200), (300, 500), (600, 700)])));
        assert!(!compatible(&r, &tx(&[(50, 200), (310, 500)])));
    }

    #[test]
    fn ends_within_ten_bp() {
        let r = [(100, 200), (300, 400)];
        assert!(compatible(&r, &tx(&[(95, 200), (300, 405)])));
        assert!(compatible(&r, &tx(&[(105, 200), (300, 400)])));
        assert!(!compatible(&r, &tx(&[(115, 200), (300, 400)])));
    }

    #[test]
    fn single_exon_read_inside_one_exon() {
        assert!(compatible(&[(120, 180)], &tx(&[(100, 200)])));
        assert!(!compatible(&[(120, 180)], &tx(&[(100, 150), (300, 400)])));
    }

    #[test]
    fn a_long_transcript_does_not_hide_short_ones_or_get_missed() {
        // review finding 4: a 2 Mb model on the contig; reads far from it still find their short transcript, and a read
        // inside the long model's last exon finds the long model
        let t = vec![tx(&[(0, 100), (1_999_000, 2_000_000)]), tx(&[(1_500_000, 1_500_200)])];
        let reads = vec![
            ("c".to_string(), vec![(1_500_050, 1_500_150)]),
            ("c".to_string(), vec![(1_999_100, 1_999_900)]),
        ];
        let (c, _, comp) = count(&t, reads.into_iter());
        assert_eq!((c, comp), (vec![1.0, 1.0], 2));
    }

    #[test]
    fn a_read_on_two_transcripts_adds_half_to_each() {
        let t = vec![tx(&[(50, 200), (300, 500)]), tx(&[(80, 200), (300, 450)])];
        let (c, seen, comp) = count(&t, vec![("c".to_string(), vec![(100, 200), (300, 400)])].into_iter());
        assert_eq!((c, seen, comp), (vec![0.5, 0.5], 1, 1));
    }
}

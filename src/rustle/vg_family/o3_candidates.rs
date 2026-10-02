//! O3 candidate copies (spec docs/superpowers/specs/2026-10-02-o3-candidates-design.md): the reference-absent-copy chain of
//! PREREG_rna_allele_haplotype_count Amendments 7-11 without IsoCon. Pure functions here; BAM passes and minimap2 calls in the binary.
//!
//! **STATUS:** TEST-ONLY — no caller yet (the `o3_candidates` binary and the driver stage land in later tasks of the plan; flips to OPT-IN then)  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use std::collections::HashMap;

pub const KMER_K: usize = 31;
pub const SKETCH_W: usize = 5;
pub const MIN_UNMAPPED_LEN: usize = 300;
const ATTRIB_MIN_FRAC: f64 = 0.30;
const ATTRIB_MIN_RATIO: f64 = 2.0;

fn code(b: u8) -> Option<u64> { match b { b'A' | b'a' => Some(0), b'C' | b'c' => Some(1), b'G' | b'g' => Some(2), b'T' | b't' => Some(3), _ => None } }

/// Canonical (min of forward / reverse-complement) 2-bit k-mers, k <= 32; windows with N are skipped.
pub fn canonical_kmers(seq: &[u8], k: usize) -> Vec<u64> {
    let mask: u64 = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let (mut fw, mut rv, mut valid) = (0u64, 0u64, 0usize);
    let mut out = Vec::with_capacity(seq.len().saturating_sub(k) + 1);
    for &b in seq {
        match code(b) {
            Some(c) => { fw = ((fw << 2) | c) & mask; rv = (rv >> 2) | ((3 - c) << (2 * (k - 1))); valid += 1; }
            None => { valid = 0; fw = 0; rv = 0; }
        }
        if valid >= k { out.push(fw.min(rv)); }
    }
    out
}

fn mix(x: u64) -> u64 { let mut z = x.wrapping_add(0x9E3779B97F4A7C15); z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9); z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB); z ^ (z >> 31) }

/// (k, w) minimizers over the canonical k-mers: the smallest hashed k-mer of every window of w consecutive k-mers, deduplicated
/// consecutively; sorted and deduplicated so sketches compare by merge.
pub fn minimizer_sketch(seq: &[u8], k: usize, w: usize) -> Vec<u64> {
    let km: Vec<u64> = canonical_kmers(seq, k).into_iter().map(mix).collect();
    let mut out = Vec::new();
    if km.len() < w { out.extend(km.iter().copied().min()); }
    for win in km.windows(w) { let m = *win.iter().min().unwrap(); if out.last() != Some(&m) { out.push(m); } }
    out.sort_unstable(); out.dedup(); out
}

/// Shared minimizers / size of the smaller sketch (both sorted, deduplicated).
pub fn sketch_share(a: &[u64], b: &[u64]) -> f64 {
    let (mut i, mut j, mut shared) = (0, 0, 0usize);
    while i < a.len() && j < b.len() { if a[i] == b[j] { shared += 1; i += 1; j += 1; } else if a[i] < b[j] { i += 1; } else { j += 1; } }
    shared as f64 / a.len().min(b.len()).max(1) as f64
}

/// Canonical 31-mers of every copy sequence -> the families carrying them (k-mers in > max_families families dropped as repeats).
pub struct FamilyKmerIndex { map: HashMap<u64, Vec<String>>, k: usize }
impl FamilyKmerIndex {
    pub fn build(copies: &[(String, Vec<u8>)], k: usize, max_families: usize) -> Self {
        let mut map: HashMap<u64, Vec<String>> = HashMap::new();
        for (fam, seq) in copies {
            for km in canonical_kmers(seq, k) { let v = map.entry(km).or_default(); if !v.contains(fam) { v.push(fam.clone()); } }
        }
        map.retain(|_, v| v.len() <= max_families);
        FamilyKmerIndex { map, k }
    }
    /// The family with the most k-mer hits when >= 30% of the read's k-mers hit it and it leads the runner-up >= 2x; reads < 300 bp: None.
    pub fn attribute(&self, read: &[u8]) -> Option<String> {
        if read.len() < MIN_UNMAPPED_LEN { return None; }
        let kms = canonical_kmers(read, self.k);
        let mut hits: HashMap<&str, usize> = HashMap::new();
        for km in &kms { if let Some(f) = self.map.get(km) { for fam in f { *hits.entry(fam).or_insert(0) += 1; } } }
        let mut v: Vec<(&str, usize)> = hits.into_iter().collect();
        v.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(b.0)));
        let (best, n) = *v.first()?;
        let second = v.get(1).map(|x| x.1).unwrap_or(0);
        if (n as f64) < ATTRIB_MIN_FRAC * kms.len() as f64 { return None; }
        if second > 0 && (n as f64) < ATTRIB_MIN_RATIO * second as f64 { return None; }
        Some(best.to_string())
    }
}

/// One minimap2 PAF alignment: columns 1-11 (mapq is not kept) plus the `de:f:` and `cs:Z:` tags. Coordinates are 0-based half-open as in
/// PAF; `strand` is `b'+'` or `b'-'`; `de` is minimap2's gap-compressed divergence (1.0 when the tag is absent); `cs` is the raw string.
#[derive(Clone, Debug, PartialEq)]
pub struct PafHit {
    pub q: String, pub qlen: usize, pub qs: usize, pub qe: usize, pub strand: u8,
    pub t: String, pub tlen: usize, pub ts: usize, pub te: usize,
    pub matches: usize, pub block: usize, pub de: f64, pub cs: Option<String>,
}

/// Parses minimap2 PAF text into one `PafHit` per alignment line. Lenient, like the repo's other PAF readers
/// (`shared_definition::parse_paf`): a line with fewer than 12 tab-separated columns, a non-numeric fixed column or a strand other than
/// `+`/`-` (blank lines, a truncated last line) is skipped. Tags are looked up after column 12: `de:f:` (absent or unreadable -> 1.0, the
/// worst divergence: such a hit never links and never merges) and `cs:Z:` (absent -> `None`).
pub fn parse_paf(text: &str) -> Vec<PafHit> { text.lines().filter_map(parse_paf_line).collect() }

fn parse_paf_line(line: &str) -> Option<PafHit> {
    let f: Vec<&str> = line.split('\t').collect();
    if f.len() < 12 { return None; }
    let tag = |prefix: &str| f[12..].iter().find_map(|x| x.strip_prefix(prefix));
    Some(PafHit {
        q: f[0].to_string(), qlen: f[1].parse().ok()?, qs: f[2].parse().ok()?, qe: f[3].parse().ok()?,
        strand: match f[4] { "+" => b'+', "-" => b'-', _ => return None },
        t: f[5].to_string(), tlen: f[6].parse().ok()?, ts: f[7].parse().ok()?, te: f[8].parse().ok()?,
        matches: f[9].parse().ok()?, block: f[10].parse().ok()?,
        de: tag("de:f:").and_then(|v| v.parse::<f64>().ok()).unwrap_or(1.0),
        cs: tag("cs:Z:").map(str::to_string),
    })
}

/// Per query name, the hit with the most matching bases; a tie keeps the first encountered (stable).
pub fn best_by_matches(hits: &[PafHit]) -> HashMap<String, PafHit> {
    let mut best: HashMap<String, PafHit> = HashMap::new();
    for h in hits { if best.get(&h.q).map_or(true, |b| h.matches > b.matches) { best.insert(h.q.clone(), h.clone()); } }
    best
}

// The three quantities below divide by lengths from the PAF; a zero length counts as 1 (no NaN / inf) and an end before its start as an
// empty span (no usize wrap-around).

/// Identity x coverage = `matches / block x (qe - qs) / qlen`, in the chain's order of evaluation (`bench/rna_allele/link_test.py:141`).
/// `>= 0.999` means the sequence is already in the reference; below it the sequence is flagged.
pub fn id_cov(h: &PafHit) -> f64 { h.matches as f64 / h.block.max(1) as f64 * h.qe.saturating_sub(h.qs) as f64 / h.qlen.max(1) as f64 }

/// Whole-length divergence `1 - matches / qlen` (`bench/rna_allele/yag_test.py:245`): the quantity the link rule compares with delta.
pub fn whole_length_d(h: &PafHit) -> f64 { 1.0 - h.matches as f64 / h.qlen.max(1) as f64 }

/// The aligned span on the SHORTER of query / target over that sequence's length; the query on a tie (`qlen <= tlen`, as `best_pairs` in
/// `bench/rna_allele/merge_test.py`). The merge and clustering rules join a pair when this is >= 0.5.
pub fn shorter_cov(h: &PafHit) -> f64 {
    if h.qlen <= h.tlen { h.qe.saturating_sub(h.qs) as f64 / h.qlen.max(1) as f64 } else { h.te.saturating_sub(h.ts) as f64 / h.tlen.max(1) as f64 }
}

/// One operation of a minimap2 `cs:Z:` string (the short form written by `--cs`; the long form of `--cs=long` is not supported).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CsOp {
    /// `:N` — N identical bases.
    Eq(usize),
    /// `*xy` — a substitution, reference base x -> query base y, bytes as written (lowercase from minimap2).
    Sub(u8, u8),
    /// `+seq` — bases in the query only (an insertion in the query).
    Ins(Vec<u8>),
    /// `-seq` — bases in the reference only (a deletion from the reference).
    Del(Vec<u8>),
    /// `~xxNNyy` — an intron of NN reference bases (the donor `xx` and acceptor `yy` dinucleotides are dropped).
    Intron(usize),
}

fn is_cs_base(b: u8) -> bool { matches!(b.to_ascii_lowercase(), b'a' | b'c' | b'g' | b't' | b'n') }

/// The decimal number at `b[i..]` and the index after it; `None` when there is no digit or the value overflows `usize`.
fn cs_number(b: &[u8], mut i: usize) -> Option<(usize, usize)> {
    let start = i;
    let mut n = 0usize;
    while i < b.len() && b[i].is_ascii_digit() { n = n.checked_mul(10)?.checked_add((b[i] - b'0') as usize)?; i += 1; }
    (i > start).then_some((n, i))
}

fn cs_err(b: &[u8], at: usize, what: &str) -> anyhow::Error {
    anyhow::anyhow!("malformed cs string: {what} at byte {at}, near `{}`", String::from_utf8_lossy(&b[at.saturating_sub(12)..(at + 12).min(b.len())]))
}

/// Parses a minimap2-produced `cs:Z:` short string into its operations (`:N` `*xy` `+seq` `-seq` `~xxNNyy`; bases are `acgtn`, either case).
/// Anything else is a hard error naming the byte offset and its neighbourhood, never a partial list and never a panic: an unknown operation
/// character, `:` or the intron length without digits, a number that overflows `usize`, `*` or `~` without its two bases, `+`/`-` without bases.
pub fn parse_cs(cs: &str) -> anyhow::Result<Vec<CsOp>> {
    let b = cs.as_bytes();
    let (mut ops, mut i) = (Vec::new(), 0usize);
    let two_bases = |i: usize| i + 1 < b.len() && is_cs_base(b[i]) && is_cs_base(b[i + 1]);
    while i < b.len() {
        let at = i;
        i += 1;
        match b[at] {
            b':' => {
                let (n, j) = cs_number(b, i).ok_or_else(|| cs_err(b, at, "`:` without a length"))?;
                ops.push(CsOp::Eq(n)); i = j;
            }
            b'*' => {
                if !two_bases(i) { return Err(cs_err(b, at, "`*` without two bases")); }
                ops.push(CsOp::Sub(b[i], b[i + 1])); i += 2;
            }
            op @ (b'+' | b'-') => {
                let j = i + b[i..].iter().take_while(|&&x| is_cs_base(x)).count();
                if j == i { return Err(cs_err(b, at, "`+`/`-` without bases")); }
                let seq = b[i..j].to_vec();
                ops.push(if op == b'+' { CsOp::Ins(seq) } else { CsOp::Del(seq) }); i = j;
            }
            b'~' => {
                if !two_bases(i) { return Err(cs_err(b, at, "`~` without a two-base donor")); }
                let (n, j) = cs_number(b, i + 2).ok_or_else(|| cs_err(b, at, "`~` without an intron length"))?;
                if !two_bases(j) { return Err(cs_err(b, at, "`~` without a two-base acceptor")); }
                ops.push(CsOp::Intron(n)); i = j + 2;
            }
            _ => return Err(cs_err(b, at, "unknown operation")),
        }
    }
    Ok(ops)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn rand_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed; (0..n).map(|_| { x ^= x << 13; x ^= x >> 7; x ^= x << 17; b"ACGT"[(x % 4) as usize] }).collect()
    }
    fn mutate(s: &[u8], rate_per_kb: usize, seed: u64) -> Vec<u8> {
        // substitutes every (1000 / rate_per_kb)-th base; the new base is always different from the original
        let mut v = s.to_vec(); let mut x = seed.max(1);
        for i in (0..v.len()).step_by(1000 / rate_per_kb.max(1)) {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            let orig = b"ACGT".iter().position(|&c| c == v[i]).expect("ACGT input");
            v[i] = b"ACGT"[(orig + 1 + (x % 3) as usize) % 4];   // index + 1..=3 (mod 4) never returns the original index
        }
        v
    }
    #[test]
    fn canonical_kmers_are_strand_free() {
        let s = rand_seq(500, 7);
        let rc = crate::vg_family::seq_utils::reverse_complement(&s);
        let mut a = canonical_kmers(&s, 31); let mut b = canonical_kmers(&rc, 31);
        a.sort(); b.sort();
        assert_eq!(a, b);
    }
    #[test]
    fn sketch_share_separates_copies_from_errors() {
        let copy_a = rand_seq(3000, 1);
        let same_copy_read = mutate(&copy_a, 2, 3);        // 0.2%: HiFi-like (every 500th base: 6 substitutions in 3 kb)
        let copy_b = mutate(&copy_a, 25, 5);               // 2.5% diverged paralog (every 40th base: 75 substitutions in 3 kb)
        let real_subs = |x: &[u8]| x.iter().zip(&copy_a).filter(|(p, q)| p != q).count();
        assert_eq!((real_subs(&same_copy_read), real_subs(&copy_b)), (6, 75), "every substitution must be real: exactly 0.2% and 2.5% of 3000 bp");
        let sa = minimizer_sketch(&copy_a, 31, 5);
        assert!(sketch_share(&sa, &minimizer_sketch(&same_copy_read, 31, 5)) > 0.8);
        assert!(sketch_share(&sa, &minimizer_sketch(&copy_b, 31, 5)) < 0.6);
    }
    #[test]
    fn attribution_needs_a_clear_winner_and_ignores_short_reads() {
        let fam1 = rand_seq(2000, 11); let fam2 = rand_seq(2000, 12);
        let idx = FamilyKmerIndex::build(&[("F1".into(), fam1.clone()), ("F2".into(), fam2.clone())], 31, 8);
        assert_eq!(idx.attribute(&mutate(&fam1, 2, 2)).as_deref(), Some("F1"));
        assert_eq!(idx.attribute(&rand_seq(1500, 99)), None);               // nothing hits
        let half: Vec<u8> = fam1[..1000].iter().chain(fam2[..1000].iter()).copied().collect();
        assert_eq!(idx.attribute(&half), None);                             // no 2x winner
        assert_eq!(idx.attribute(&fam1[..69]), None);                       // short_unmapped_reads_are_ignored
    }

    fn paf_hit(q: &str, t: &str, matches: usize) -> PafHit {
        PafHit { q: q.into(), qlen: 1000, qs: 0, qe: 1000, strand: b'+', t: t.into(), tlen: 1000, ts: 0, te: 1000, matches, block: 1000, de: 0.01, cs: None }
    }
    #[test]
    fn paf_quantities_match_the_chain() {
        let line = "r1\t1000\t10\t990\t+\tc1\t5000\t100\t1100\t960\t1000\t60\tNM:i:40\tde:f:0.0150\tcs:Z::500*ac:300+tt:179-g:0";
        let h = &parse_paf(line)[0];
        assert!((whole_length_d(h) - 0.04).abs() < 1e-9);          // 1 - 960/1000
        assert!((id_cov(h) - 0.96 * 0.98).abs() < 1e-9);
        assert!((h.de - 0.015).abs() < 1e-9);
        let ops = parse_cs(h.cs.as_deref().unwrap()).unwrap();
        assert!(matches!(ops[1], CsOp::Sub(b'a', b'c')));
        assert!(matches!(&ops[3], CsOp::Ins(v) if v == b"tt"));
        assert!(matches!(&ops[5], CsOp::Del(v) if v == b"g"));
    }
    #[test]
    fn parse_paf_reads_many_lines_and_defaults_absent_tags() {
        // a CRLF-terminated '-' hit with both tags, a blank line, a hit with another tag only, three unusable lines (truncated to 4
        // columns, strand '.', non-numeric qlen: all skipped), and a hit with exactly the 12 fixed columns and no tag at all
        let text = [
            "a\t120\t3\t117\t-\tb\t260\t50\t170\t105\t118\t60\tde:f:0.01\tcs:Z::100\r",
            "",
            "c\t300\t10\t290\t+\td\t350\t5\t285\t270\t281\t0\ttp:A:P",
            "e\t100\t0\t5",
            "x\t100\t0\t100\t.\ty\t100\t0\t100\t100\t100\t60\tde:f:0.0",
            "z\tNA\t0\t100\t+\ty\t100\t0\t100\t100\t100\t60\tde:f:0.0",
            "g\t50\t1\t49\t+\th\t60\t2\t52\t40\t48\t30",
        ].join("\n") + "\n";
        let hits = parse_paf(&text);
        assert_eq!(hits.len(), 3);
        let (a, c, g) = (&hits[0], &hits[1], &hits[2]);
        // every fixed column lands in its own field (all values distinct so a swapped column cannot pass)
        assert_eq!((a.q.as_str(), a.qlen, a.qs, a.qe, a.strand, a.t.as_str(), a.tlen, a.ts, a.te, a.matches, a.block), ("a", 120, 3, 117, b'-', "b", 260, 50, 170, 105, 118));
        assert!((a.de - 0.01).abs() < 1e-12);
        assert_eq!(a.cs.as_deref(), Some(":100"));              // the CR of CRLF is not part of the cs string
        assert_eq!((c.q.as_str(), c.qlen, c.qs, c.qe, c.strand, c.t.as_str(), c.tlen, c.ts, c.te, c.matches, c.block), ("c", 300, 10, 290, b'+', "d", 350, 5, 285, 270, 281));
        assert_eq!((c.de, c.cs.as_deref()), (1.0, None));       // absent de -> 1.0 (never clusters, never links), absent cs -> None
        assert_eq!((g.q.as_str(), g.t.as_str(), g.matches, g.de, g.cs.as_deref()), ("g", "h", 40, 1.0, None));
        assert!(parse_paf("").is_empty() && parse_paf("\n\n").is_empty());
    }
    #[test]
    fn parse_cs_walks_every_operation_and_rejects_malformed_strings() {
        let ops = parse_cs(":10*ag~gt1500ag+acg-t:5*nc").unwrap();
        assert_eq!(ops, vec![CsOp::Eq(10), CsOp::Sub(b'a', b'g'), CsOp::Intron(1500), CsOp::Ins(b"acg".to_vec()), CsOp::Del(b"t".to_vec()), CsOp::Eq(5), CsOp::Sub(b'n', b'c')]);
        assert_eq!(parse_cs("").unwrap(), vec![]);
        assert_eq!(parse_cs(":0").unwrap(), vec![CsOp::Eq(0)]);                       // a zero-length match is kept, not dropped
        assert_eq!(parse_cs("*AG+ACG").unwrap(), vec![CsOp::Sub(b'A', b'G'), CsOp::Ins(b"ACG".to_vec())]);   // either case, bytes kept as written
        // every malformed string is an error, including ones whose prefix is a valid operation list (never a partial result)
        for bad in [":", ":x", "*a", "*ac*", "*a-", "+", "-:3", "~", "~gt", "~gt15", "~gt15a", "~gt15ag7", "~g1ag", "~1215ag", "~gtag", ":10;", "10", " :3", "=ACGT", ":5\n", "*aé", ":99999999999999999999999"] {
            assert!(parse_cs(bad).is_err(), "{bad:?} must be rejected");
        }
        let msg = parse_cs(":10;").unwrap_err().to_string();
        assert!(msg.contains("byte 3") && msg.contains(":10;"), "the error names the offset and the neighbourhood: {msg}");
    }
    #[test]
    fn best_by_matches_keeps_the_hit_with_most_matches_per_query() {
        let hits = vec![
            paf_hit("r1", "t_low", 900), paf_hit("r2", "t_only", 500), paf_hit("r1", "t_high", 950), paf_hit("r1", "t_later_low", 940),
            paf_hit("r3", "t_first", 700), paf_hit("r3", "t_tie", 700),
        ];
        let best = best_by_matches(&hits);
        assert_eq!(best.len(), 3);
        assert_eq!(best["r1"].t, "t_high");                     // the higher-matches hit wins whatever its position
        assert_eq!(best["r2"].t, "t_only");
        assert_eq!(best["r3"].t, "t_first");                    // a tie keeps the first encountered
        assert!(best_by_matches(&[]).is_empty());
    }
    #[test]
    fn shorter_cov_uses_the_span_on_the_shorter_sequence() {
        let mut h = paf_hit("q", "t", 700);
        (h.qlen, h.qs, h.qe) = (3000, 100, 1000);               // query: 900 aligned bases of 3000
        (h.tlen, h.ts, h.te) = (1000, 150, 950);                // target: 800 aligned bases of 1000, the shorter sequence
        assert!((shorter_cov(&h) - 0.8).abs() < 1e-12);         // 800 / 1000, not 900 / 3000
        (h.qlen, h.qs, h.qe, h.tlen, h.ts, h.te) = (1000, 150, 950, 3000, 100, 1000);   // roles swapped: the query is the shorter
        assert!((shorter_cov(&h) - 0.8).abs() < 1e-12);
        (h.qlen, h.qs, h.qe, h.tlen, h.ts, h.te) = (1000, 0, 900, 1000, 0, 800);        // equal lengths: the query's span (the chain's qlen <= tlen)
        assert!((shorter_cov(&h) - 0.9).abs() < 1e-12);
    }
    #[test]
    fn quantities_divide_by_the_right_lengths() {
        // distinct qlen / tlen / block, so that a swapped denominator changes the value (the brief's hit has qlen == block)
        let mut h = paf_hit("q", "t", 950);
        (h.qlen, h.qs, h.qe, h.tlen, h.ts, h.te, h.block) = (1000, 20, 980, 4000, 100, 1100, 1100);
        assert!((id_cov(&h) - 950.0 / 1100.0 * 960.0 / 1000.0).abs() < 1e-12);
        assert!((whole_length_d(&h) - 0.05).abs() < 1e-12);               // 1 - 950 / qlen: neither / block nor / tlen
    }
    #[test]
    fn quantities_stay_finite_on_empty_or_inverted_hits() {
        let mut h = paf_hit("q", "t", 0);
        (h.qlen, h.qs, h.qe, h.tlen, h.ts, h.te, h.block) = (0, 0, 0, 0, 0, 0, 0);
        assert_eq!((id_cov(&h), whole_length_d(&h), shorter_cov(&h)), (0.0, 1.0, 0.0));
        h.qlen = 10;                                                                  // an empty target: the target-side division is guarded too
        assert_eq!(shorter_cov(&h), 0.0);
        (h.qlen, h.qs, h.qe, h.tlen, h.ts, h.te, h.block, h.matches) = (100, 60, 40, 100, 60, 40, 10, 5);   // end before start: zero span, not a wrapped one
        assert_eq!((id_cov(&h), shorter_cov(&h)), (0.0, 0.0));
    }
}

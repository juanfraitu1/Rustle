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

/// Indels of at least this many bases are STRUCTURE (isoforms, ruling R2), not sequencing errors; shorter ones follow the 50% majority.
const STRUCT_MIN_INDEL: usize = 20;
/// Members that must cover a column before the vote (or a < 20 bp indel majority) may change the template there.
const VOTE_MIN_COVER: usize = 3;
/// Members that must carry an insertion of >= 20 bp for it to be inserted (an exon the template lacks), whatever its share.
const STRUCT_MIN_SUPPORT: usize = 3;
/// Members that must cover a template column at either end of the consensus for it to survive the end trim (spec §5.4, ruling R5).
const END_MIN_COVER: usize = 2;

/// Clusters one family's reads (spec §5.3, with minimap2's all-vs-all hits in place of the greedy pass; ruling §9b): union-find over the reads,
/// joining a pair when its best hit covers >= 50% of the shorter read (`shorter_cov`) with `de <= delta`. The pair's best hit is the one with
/// the most matching bases, the first on a tie, whichever read is `q` (`best_pairs` of `bench/rna_allele/merge_test.py`, which this
/// reproduces). Self hits and hits naming a read that is not in `names` are ignored; the strand is not consulted (`de` and the spans do not
/// depend on orientation). `names` must be distinct. Returns the clusters as ascending index lists ordered by (size descending, first
/// index). The clusters are a function of the set of joined pairs, so neither the hash order of the pair table nor the order of `ava`
/// (beyond the first-wins tie rule) reaches the output.
pub fn cluster_reads(names: &[String], ava: &[PafHit], delta: f64) -> Vec<Vec<usize>> {
    let idx: HashMap<&str, usize> = names.iter().enumerate().map(|(i, n)| (n.as_str(), i)).collect();
    let mut best: HashMap<(usize, usize), &PafHit> = HashMap::new();
    for h in ava {
        let (Some(&a), Some(&b)) = (idx.get(h.q.as_str()), idx.get(h.t.as_str())) else { continue };
        if a == b { continue; }
        let key = (a.min(b), a.max(b));
        if best.get(&key).map_or(true, |o| h.matches > o.matches) { best.insert(key, h); }
    }
    fn find(p: &mut [usize], mut x: usize) -> usize { while p[x] != x { p[x] = p[p[x]]; x = p[x]; } x }
    let mut par: Vec<usize> = (0..names.len()).collect();
    for (&(a, b), h) in &best {
        if shorter_cov(h) >= 0.5 && h.de <= delta {
            let (ra, rb) = (find(&mut par, a), find(&mut par, b));
            if ra != rb { par[ra.max(rb)] = ra.min(rb); }
        }
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..names.len() { let r = find(&mut par, i); groups.entry(r).or_default().push(i); }   // pushed in index order: each group is ascending
    let mut out: Vec<Vec<usize>> = groups.into_values().collect();
    out.sort_by(|a, b| b.len().cmp(&a.len()).then(a[0].cmp(&b[0])));                                // first indices are distinct: a total order
    out
}

/// Template-and-vote consensus of one cluster (spec §5.4, rulings R2 and R5). `template` is the sequence the members were aligned to; each
/// member's hit carries its `cs` against it (`ts..te` on the template). The member sequences are not read (the `cs` holds every base the vote
/// needs) and the strand is not consulted (a `cs` always lies along the target's forward strand). A member without a `cs` abstains; a
/// malformed `cs` is an error naming the member. The template votes (and covers) only if it is itself among the members, as a self hit.
/// Per template column, over the members that COVER it: `Eq`, `Sub` and `< 20 bp` `Del` columns cover; an `Intron` that a member splices out
/// and a `>= 20 bp` `Del` (an isoform's absence) do not: such a member carries no base there, so it is neither a vote for the template base
/// nor a vote to delete (R5).
/// * a column covered by >= 3 members takes the substituted base (the query base, upper case) with the most votes when that count exceeds the
///   number of covering members that carry the template base (no substitution, no deletion there); a tie keeps the template; two variants
///   with equal support take the smaller base; fewer than 3 covering members keep the template base;
/// * an insertion < 20 bp before a column (after the last column when `t == n`) is inserted when its most frequent sequence (ties to the
///   smaller) is carried by >= 50% of >= 3 covering members; a deletion < 20 bp removes the column when >= 50% of >= 3 covering members
///   delete it (a plain majority of the covering members, also inside an exon that other members skip);
/// * indels >= 20 bp are STRUCTURE, not errors (R2): an insertion >= 20 bp whose most frequent sequence is carried by >= 3 members is inserted
///   whatever its share (an exon the template lacks); a deletion >= 20 bp is never applied (the template's exon stays), so the consensus is
///   the exon union of the cluster's isoforms with SNV-level majority voting.
///
/// End trim (spec §5.4, R5): after the vote, the leading and trailing template columns covered by fewer than 2 members are dropped. An
/// insertion goes with the column it precedes (the one after the last column, with the last column). When no column is covered by 2 members
/// nothing survives: the consensus is empty (as it is for an empty member list). The output is upper case. A `cs` that runs past the template
/// end is clamped (nothing beyond the end is voted or appended); the frame of the `cs` (`ts` and the template) is trusted, not checked.
pub fn consensus_from_template(template: &[u8], member_hits: &[(&[u8], &PafHit)]) -> anyhow::Result<Vec<u8>> {
    let n = template.len();
    let mut cover = vec![0usize; n];
    let mut subs: Vec<HashMap<u8, usize>> = vec![HashMap::new(); n];
    let mut dels = vec![0usize; n];
    let mut ins: Vec<HashMap<Vec<u8>, usize>> = vec![HashMap::new(); n + 1];
    for (_, h) in member_hits {
        let Some(cs) = h.cs.as_deref() else { continue };
        let ops = parse_cs(cs).map_err(|e| anyhow::anyhow!("member {} against {}: {e}", h.q, h.t))?;
        let mut t = h.ts;
        for op in ops {
            match op {
                CsOp::Eq(len) => { for p in t..t.saturating_add(len).min(n) { cover[p] += 1; } t = t.saturating_add(len); }
                CsOp::Sub(_, qb) => { if t < n { cover[t] += 1; *subs[t].entry(qb.to_ascii_uppercase()).or_insert(0) += 1; } t = t.saturating_add(1); }
                // a < 20 bp deletion is an error-sized gap: the member covers the columns and votes to delete them; a >= 20 bp one is an isoform's
                // absence (R5, as for an `Intron`): the member carries no base there, so it neither covers nor votes
                CsOp::Del(seq) => {
                    if seq.len() < STRUCT_MIN_INDEL { for p in t..t.saturating_add(seq.len()).min(n) { cover[p] += 1; dels[p] += 1; } }
                    t = t.saturating_add(seq.len());
                }
                CsOp::Ins(mut seq) => { if t <= n { seq.make_ascii_uppercase(); *ins[t].entry(seq).or_insert(0) += 1; } }
                CsOp::Intron(len) => { t = t.saturating_add(len); }
            }
        }
    }
    // spec §5.4 / R5 end trim: columns lo..=hi survive, the leading and trailing ones covered by fewer than END_MIN_COVER members do not; an
    // insertion goes with the column it precedes, so the one after the last column survives only with the last column. No covered column
    // (an empty template, no member) returns here, so n >= 1 below.
    let (Some(lo), Some(hi)) = (cover.iter().position(|&c| c >= END_MIN_COVER), cover.iter().rposition(|&c| c >= END_MIN_COVER)) else { return Ok(Vec::new()); };
    let last = if hi + 1 == n { n } else { hi };
    let mut out = Vec::with_capacity(hi - lo + 64);
    for p in lo..=last {
        // the most frequent insertion before column p: the larger count, then the smaller sequence (a total order: no hash order reaches the output)
        if let Some((seq, &cnt)) = ins[p].iter().max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0))) {
            let covering = cover[p.min(n - 1)];
            let structure = seq.len() >= STRUCT_MIN_INDEL && cnt >= STRUCT_MIN_SUPPORT;
            if structure || (seq.len() < STRUCT_MIN_INDEL && covering >= VOTE_MIN_COVER && 2 * cnt >= covering) { out.extend_from_slice(seq); }
        }
        if p == n { break; }
        if cover[p] >= VOTE_MIN_COVER && 2 * dels[p] >= cover[p] { continue; }
        let mut base = template[p].to_ascii_uppercase();
        if cover[p] >= VOTE_MIN_COVER {
            let same = cover[p].saturating_sub(subs[p].values().sum::<usize>() + dels[p]);
            if let Some((&b, &c)) = subs[p].iter().max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0))) { if c > same { base = b; } }
        }
        out.push(base);
    }
    Ok(out)
}

/// The number of substitutions in a pairwise `cs`: spec §5.5's k, the columns that distinguish two consensus sequences. Indels are not counted.
pub fn distinguishing_columns(ops: &[CsOp]) -> usize { ops.iter().filter(|o| matches!(o, CsOp::Sub(..))).count() }

/// P(X >= n_small) for X ~ Binomial(n_small + n_large, eps^k), in log space with a log-factorial table; the variant is real iff it is < alpha
/// (spec §5.5: n_small = the reads of the smaller cluster, n_large = the larger's, k = `distinguishing_columns`). No distinguishing column or
/// no reads in the small cluster is never real; an `eps^k` that underflows to 0 gives P = 0: real.
pub fn variant_is_real(n_small: usize, n_large: usize, k: usize, eps: f64, alpha: f64) -> bool {
    if k == 0 || n_small == 0 { return false; }
    let n = n_small + n_large;
    let p = eps.powi(i32::try_from(k).unwrap_or(i32::MAX));
    let (lp, lq) = (p.ln(), (1.0 - p).ln());
    let mut lf = vec![0f64; n + 1];
    for i in 1..=n { lf[i] = lf[i - 1] + (i as f64).ln(); }
    let tail: f64 = (n_small..=n).map(|x| (lf[n] - lf[x] - lf[n - x] + x as f64 * lp + (n - x) as f64 * lq).exp()).sum();
    tail < alpha
}

/// One refinement pass over a cluster (spec §5.4's one polishing pass, then the split). `members` pairs each member read with its hit against
/// the consensus; `template` is that consensus. Neither it nor the read sequences are consulted (the rule reads the hits only; the parameters
/// keep the signature parallel to `consensus_from_template`). Returns the ascending indices of the members that fit the consensus
/// (`de <= delta` and `shorter_cov >= 0.5`) and of the rest; the caller makes the rest a new cluster when it holds >= `--min-cluster` reads
/// and places a member that has no hit at all.
pub fn refine_cluster(_template: &[u8], members: &[(&[u8], &PafHit)], delta: f64) -> (Vec<usize>, Vec<usize>) {
    (0..members.len()).partition(|&i| { let h = members[i].1; h.de <= delta && shorter_cov(h) >= 0.5 })
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

    // ---- read clustering, template-and-vote consensus, the real-vs-error test, one refinement pass --------------------------------------

    /// A hit between two 3000 bp reads covering `cov` of the shorter one (both are 3000 bp), for the clustering tests.
    fn ava_hit(a: &str, b: &str, matches: usize, cov: f64, de: f64) -> PafHit {
        let span = (3000.0 * cov).round() as usize;
        PafHit { q: a.into(), qlen: 3000, qs: 0, qe: span, strand: b'+', t: b.into(), tlen: 3000, ts: 0, te: span, matches, block: 3000, de, cs: None }
    }
    fn names_of(s: &str) -> Vec<String> { s.split_whitespace().map(str::to_string).collect() }

    /// An edit of a synthetic member relative to the template, in template coordinates (ascending, non-overlapping).
    enum Edit<'a> { Sub(usize, u8), Ins(usize, &'a [u8]), Del(usize, usize) }
    impl Edit<'_> { fn at(&self) -> usize { match self { Edit::Sub(p, _) | Edit::Ins(p, _) | Edit::Del(p, _) => *p } } }
    /// A synthetic member: its sequence (the template with `edits` applied) and its full-length `+` hit against the template, with the short
    /// `cs` minimap2 would write (bases lower case). `matches`, `block` and `de` are fillers: the vote reads `ts` and `cs` only.
    fn member_of(template: &[u8], edits: &[Edit]) -> (Vec<u8>, PafHit) {
        let lower = |s: &[u8]| String::from_utf8(s.to_ascii_lowercase()).unwrap();
        let (mut seq, mut cs, mut t) = (Vec::new(), String::new(), 0usize);
        for e in edits {
            assert!(e.at() >= t, "edits must ascend and not overlap");
            if e.at() > t { cs += &format!(":{}", e.at() - t); seq.extend_from_slice(&template[t..e.at()]); }
            t = e.at();
            match e {
                Edit::Sub(p, b) => { cs += &format!("*{}{}", lower(&template[*p..*p + 1]), lower(&[*b])); seq.push(*b); t = p + 1; }
                Edit::Ins(_, s) => { cs += &format!("+{}", lower(s)); seq.extend_from_slice(s); }
                Edit::Del(p, len) => { cs += &format!("-{}", lower(&template[*p..p + len])); t = p + len; }
            }
        }
        if t < template.len() { cs += &format!(":{}", template.len() - t); seq.extend_from_slice(&template[t..]); }
        let hit = PafHit { q: "m".into(), qlen: seq.len(), qs: 0, qe: seq.len(), strand: b'+', t: "t".into(), tlen: template.len(), ts: 0, te: template.len(), matches: template.len(), block: template.len(), de: 0.0, cs: Some(cs) };
        (seq, hit)
    }
    /// A synthetic member that aligns to the window `template[ts..te]` only: `edits` are in the window's coordinates, the hit starts at `ts`.
    fn member_span(template: &[u8], ts: usize, te: usize, edits: &[Edit]) -> (Vec<u8>, PafHit) {
        let (seq, mut h) = member_of(&template[ts..te], edits);
        (h.ts, h.te, h.tlen) = (ts, te, template.len());
        (seq, h)
    }
    /// `total` members, the first `carriers` of them with `edits` and the rest equal to the template.
    fn cluster_of(template: &[u8], carriers: usize, total: usize, edits: &[Edit]) -> Vec<(Vec<u8>, PafHit)> {
        let none: &[Edit] = &[];
        (0..total).map(|i| member_of(template, if i < carriers { edits } else { none })).collect()
    }
    fn consensus_of(template: &[u8], members: &[(Vec<u8>, PafHit)]) -> Vec<u8> {
        let refs: Vec<(&[u8], &PafHit)> = members.iter().map(|(s, h)| (s.as_slice(), h)).collect();
        consensus_from_template(template, &refs).unwrap()
    }
    /// The k-th base different from `b`.
    fn other_base(b: u8, k: usize) -> u8 { *b"ACGT".iter().filter(|&&c| c != b.to_ascii_uppercase()).nth(k).unwrap() }

    #[test]
    fn clustering_splits_two_percent_and_joins_point_two_percent() {
        // 6 reads: r0..r2 from copy A (0.2% errors), r3..r5 from copy B (2% from A); hits written as the minimap2 ava would give them
        let names: Vec<String> = (0..6).map(|i| format!("r{i}")).collect();
        let hit = |a: usize, b: usize, de: f64| PafHit { q: format!("r{a}"), qlen: 3000, qs: 0, qe: 3000, strand: b'+', t: format!("r{b}"), tlen: 3000, ts: 0, te: 3000, matches: 2900, block: 3000, de, cs: None };
        let ava = vec![hit(0,1,0.002), hit(1,2,0.003), hit(0,2,0.002), hit(3,4,0.002), hit(4,5,0.002), hit(0,3,0.021), hit(2,5,0.019)];
        let cl = cluster_reads(&names, &ava, 0.00958);
        assert_eq!(cl, vec![vec![0,1,2], vec![3,4,5]]);
    }
    #[test]
    fn clustering_decides_on_the_pairs_best_hit_and_needs_half_of_the_shorter_read() {
        let names = names_of("a b c d e f g h i j k l m n o p q r");           // indices 0..=17
        let d = 0.00958;
        let ava = vec![
            ava_hit("a", "b", 2900, 1.0, 0.002), ava_hit("b", "a", 2950, 1.0, 0.05),    // the hit with the MOST MATCHES decides, and it is too divergent: not joined
            ava_hit("c", "d", 2900, 1.0, 0.002), ava_hit("d", "c", 2900, 1.0, 0.05),    // equal matches: the first hit is kept, here a good one: joined
            ava_hit("e", "f", 2900, 1.0, 0.05), ava_hit("f", "e", 2900, 1.0, 0.002),    // ... and here the first is the divergent one: not joined
            ava_hit("g", "h", 2900, 1.0, d),                                            // de == delta joins
            ava_hit("i", "j", 2900, 0.5, 0.001),                                        // exactly half of the shorter read joins
            ava_hit("k", "l", 2900, 0.49, 0.001),                                       // just under half does not
            PafHit { strand: b'-', ..ava_hit("n", "m", 2900, 1.0, 0.002) },             // an opposite-strand pair, named the other way round: joined
            ava_hit("o", "o", 3000, 1.0, 0.0), ava_hit("o", "zz", 3000, 1.0, 0.0),      // a self hit and a hit to an unknown read are ignored
            ava_hit("p", "q", 2900, 1.0, 0.002), ava_hit("q", "r", 2900, 1.0, 0.002),    // p-q-r is one cluster with no p-r hit at all
        ];
        let want: Vec<Vec<usize>> = vec![vec![15, 16, 17], vec![2, 3], vec![6, 7], vec![8, 9], vec![12, 13], vec![0], vec![1], vec![4], vec![5], vec![10], vec![11], vec![14]];
        assert_eq!(cluster_reads(&names, &ava, d), want);                               // by size descending, then first index
    }
    #[test]
    fn clustering_of_nothing_and_of_unlinked_reads() {
        assert!(cluster_reads(&[], &[], 0.01).is_empty());
        let names = names_of("x y z");
        assert_eq!(cluster_reads(&names, &[], 0.01), vec![vec![0], vec![1], vec![2]]);   // singletons, in index order
        assert_eq!(cluster_reads(&names, &[ava_hit("z", "x", 2900, 1.0, 0.9)], 0.01), vec![vec![0], vec![1], vec![2]]);
    }
    #[test]
    fn clustering_matches_a_brute_force_reference_on_random_hit_lists() {
        // reproducible xorshift; 40 reads, 50 random hits (either direction, repeated pairs, self hits) so that components stay small and varied
        let mut x = 0x2545F4914F6CDD1Du64;
        let mut next = move || { x ^= x << 13; x ^= x >> 7; x ^= x << 17; x };
        let (n, d) = (40usize, 0.00958);
        let names: Vec<String> = (0..n).map(|i| format!("r{i}")).collect();
        for round in 0..40 {
            let ava: Vec<PafHit> = (0..50).map(|_| {
                let (a, b) = ((next() % n as u64) as usize, (next() % n as u64) as usize);
                let (m, c, e) = (2800 + (next() % 200) as usize, [0.3, 0.45, 0.5, 0.75, 1.0][(next() % 5) as usize], [0.001, 0.005, d, 0.012, 0.03][(next() % 5) as usize]);
                ava_hit(&names[a], &names[b], m, c, e)
            }).collect();
            // reference: per unordered pair the first hit with the most matches; joined iff cov >= 0.5 && de <= delta; components by label propagation
            let mut joined = vec![];
            for a in 0..n { for b in a + 1..n {
                let mut best: Option<&PafHit> = None;
                for h in &ava { if (h.q == names[a] && h.t == names[b]) || (h.q == names[b] && h.t == names[a]) { if best.map_or(true, |o| h.matches > o.matches) { best = Some(h); } } }
                if let Some(h) = best { if shorter_cov(h) >= 0.5 && h.de <= d { joined.push((a, b)); } }
            } }
            let mut label: Vec<usize> = (0..n).collect();
            loop {
                let mut changed = false;
                for &(a, b) in &joined { let m = label[a].min(label[b]); if label[a] != m || label[b] != m { label[a] = m; label[b] = m; changed = true; } }
                if !changed { break; }
            }
            let mut groups: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
            for i in 0..n { groups.entry(label[i]).or_default().push(i); }
            let mut want: Vec<Vec<usize>> = groups.into_values().collect();
            want.sort_by(|p, q| q.len().cmp(&p.len()).then(p[0].cmp(&q[0])));
            assert_eq!(cluster_reads(&names, &ava, d), want, "round {round}");
        }
    }

    #[test]
    fn consensus_vote_fixes_errors_and_keeps_majority_indels() {
        let template = b"ACGTACGTACGTTTTTACGTACGT".to_vec();          // template carries a 1-base error at index 4 (A instead of G)
        let truth    = b"ACGTGCGTACGTTTTTACGTACGT".to_vec();
        // three members equal to truth: cs vs template = ":4*ag:19"  (4 eq, sub a->g, 19 eq)
        let h = PafHit { q: "m".into(), qlen: 24, qs: 0, qe: 24, strand: b'+', t: "t".into(), tlen: 24, ts: 0, te: 24, matches: 23, block: 24, de: 0.04, cs: Some(":4*ag:19".into()) };
        let members: Vec<(&[u8], &PafHit)> = vec![(&truth[..], &h), (&truth[..], &h), (&truth[..], &h)];
        assert_eq!(consensus_from_template(&template, &members).unwrap(), truth);
    }
    #[test]
    fn consensus_keeps_an_exon_that_two_of_five_members_skip_and_inserts_a_24bp_insertion_carried_by_three() {
        // ruling R2: the template carries an exon (40..64, 24 bp) that 2 of 5 members skip (a 24 bp `-` in their cs); 3 of 5 carry a 24 bp insertion at 100
        let template = rand_seq(120, 41);
        let ins: &[u8] = b"GATTACAGATTACACCGGTTAACC";
        assert_eq!(ins.len(), 24);
        let carrier = || member_of(&template, &[Edit::Ins(100, ins)]);
        let skipper = || member_of(&template, &[Edit::Del(40, 24)]);
        let members = vec![carrier(), carrier(), carrier(), skipper(), skipper()];
        let mut want = template.clone(); want.splice(100..100, ins.iter().copied());
        let got = consensus_of(&template, &members);
        assert_eq!(got.len(), 120 + 24);                                  // the exon stayed and the insertion went in
        assert_eq!(got, want);
    }
    #[test]
    fn consensus_treats_indels_of_20_bp_and_more_as_structure_whatever_their_share() {
        let template = rand_seq(120, 43);
        let ins24: &[u8] = b"GATTACAGATTACACCGGTTAACC";
        // deletions: 19 bp follows the 50% majority (4 of 5 -> applied) but 20 bp never does, even when every member carries it
        let mut minus19 = template.clone(); minus19.drain(40..59);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 4, 5, &[Edit::Del(40, 19)])), minus19);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 4, 5, &[Edit::Del(40, 20)])), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 5, 5, &[Edit::Del(40, 24)])), template);
        // insertions: 20 bp with 3 members is inserted although 3 of 8 is only 37.5%, with 2 members it is not ...
        let with = |ins: &[u8]| { let mut v = template.clone(); v.splice(100..100, ins.iter().copied()); v };
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 8, &[Edit::Ins(100, &ins24[..20])])), with(&ins24[..20]));
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 8, &[Edit::Ins(100, &ins24[..20])])), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 2, &[Edit::Ins(100, &ins24[..20])])), template);   // 100%, but only 2 members
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 3, &[Edit::Ins(100, &ins24[..20])])), template);   // 67% of 3: 20 bp needs 3 carriers
        // ... while 19 bp needs 50% of the covering members: 3 of 8 no, 4 of 8 yes
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 8, &[Edit::Ins(100, &ins24[..19])])), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 4, 8, &[Edit::Ins(100, &ins24[..19])])), with(&ins24[..19]));
    }
    #[test]
    fn consensus_small_indels_follow_the_fifty_percent_majority_of_three_or_more_covering_members() {
        let template = rand_seq(80, 81);
        let (mut minus2, mut plus3) = (template.clone(), template.clone());
        minus2.drain(30..32); plus3.splice(50..50, b"CCA".iter().copied());
        let (del2, ins3) = ([Edit::Del(30, 2)], [Edit::Ins(50, b"CCA")]);
        // 3 of 5 and exactly 2 of 4 (>= 50%) apply; 2 of 5 does not; two covering members are below the floor of 3
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 5, &del2)), minus2);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 4, &del2)), minus2);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 5, &del2)), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 2, &del2)), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 5, &ins3)), plus3);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 4, &ins3)), plus3);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 5, &ins3)), template);
        assert_eq!(consensus_of(&template, &cluster_of(&template, 2, 2, &ins3)), template);
        // an insertion before the first column and one after the last are placed there
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 4, &[Edit::Ins(0, b"GG")])), [&b"GG"[..], &template].concat());
        assert_eq!(consensus_of(&template, &cluster_of(&template, 3, 4, &[Edit::Ins(80, b"GG")])), [&template[..], &b"GG"[..]].concat());
        // the most frequent insertion sequence wins over a rarer one, whichever sorts first (3 TT + 1 AC) ...
        let mut skew = cluster_of(&template, 3, 3, &[Edit::Ins(50, b"TT")]);
        skew.extend(cluster_of(&template, 1, 1, &[Edit::Ins(50, b"AC")]));
        let mut plus_tt = template.clone(); plus_tt.splice(50..50, b"TT".iter().copied());
        assert_eq!(consensus_of(&template, &skew), plus_tt);
        // ... and equally frequent ones tie to the smaller sequence (3 + 3 of 6, each at 50%), whatever the hash order: every call builds its own
        // HashMaps with their own random seed, so repeating the call exposes a missing tie-break
        let mut mixed = cluster_of(&template, 3, 3, &[Edit::Ins(50, b"TT")]);
        mixed.extend(cluster_of(&template, 3, 3, &[Edit::Ins(50, b"AC")]));
        let mut plus_ac = template.clone(); plus_ac.splice(50..50, b"AC".iter().copied());
        for _ in 0..40 { assert_eq!(consensus_of(&template, &mixed), plus_ac); }
    }
    #[test]
    fn consensus_substitution_needs_three_covering_members_and_more_votes_than_the_template_base() {
        let template = rand_seq(50, 71);
        let alt = |k: usize| other_base(template[20], k);
        let variant = |k: usize| member_of(&template, &[Edit::Sub(20, alt(k))]);
        let plain = || member_of(&template, &[]);
        let mut with_alt = template.clone(); with_alt[20] = alt(0);
        assert_eq!(consensus_of(&template, &[variant(0), variant(0)]), template);                              // 2 covering members: below the floor
        assert_eq!(consensus_of(&template, &[variant(0), variant(0), variant(0)]), with_alt);                  // 3 of 3
        assert_eq!(consensus_of(&template, &[variant(0), variant(0), plain()]), with_alt);                     // 2 votes > 1 template base
        assert_eq!(consensus_of(&template, &[variant(0), variant(0), plain(), plain()]), template);            // a tie keeps the template
        assert_eq!(consensus_of(&template, &[variant(0), variant(0), variant(0), plain(), plain()]), with_alt);
        // two variants with the same support: the smaller base, provided it beats the template base (whatever the hash order: repeated calls)
        let mut with_min = template.clone(); with_min[20] = alt(0).min(alt(1));
        for _ in 0..40 {
            assert_eq!(consensus_of(&template, &[variant(0), variant(0), variant(1), variant(1), plain()]), with_min);
            assert_eq!(consensus_of(&template, &[variant(1), variant(1), variant(0), variant(0), plain()]), with_min);
        }
        // the better supported variant wins over a rarer one, also when it is the larger base (3 x alt(1) against 1 x alt(0))
        let mut with_big = template.clone(); with_big[20] = alt(1);
        assert!(alt(1) > alt(0));
        assert_eq!(consensus_of(&template, &[variant(1), variant(1), variant(1), variant(0)]), with_big);
        assert_eq!(consensus_of(&template, &[variant(1), variant(1), variant(1), variant(0), plain()]), with_big);
        // members that delete the column carry no base: they are neither template carriers nor votes (2 variants, 2 deleters, 1 template base)
        let del1 = || member_of(&template, &[Edit::Del(20, 1)]);
        assert_eq!(consensus_of(&template, &[variant(0), variant(0), del1(), del1(), plain()]), with_alt);
    }
    #[test]
    fn consensus_recovers_the_truth_from_ten_reads_with_private_errors_and_a_homopolymer_indel() {
        let mut truth = rand_seq(200, 17);
        truth[99] = b'A'; truth[100..106].fill(b'T'); truth[106] = b'C';                    // a 6-T homopolymer at 100..106
        let mut template = truth.clone(); template.insert(100, b'T');                       // the template has 7 T's (an insertion error) ...
        template[151] = if truth[150] == b'G' { b'A' } else { b'G' };                       // ... and a substitution error (truth 150 = template 151)
        let fix = truth[150];
        let mut members = vec![];
        for p in [10usize, 30, 60, 90, 120, 140, 170] {                                      // 7 members: their own error, the homopolymer fixed, the template error fixed
            let mut e = vec![Edit::Sub(p, other_base(template[p], 0)), Edit::Del(100, 1), Edit::Sub(151, fix)];
            e.sort_by_key(|x| x.at()); members.push(member_of(&template, &e));
        }
        members.push(member_of(&template, &[Edit::Del(100, 2), Edit::Sub(151, fix)]));       // one that also lost a second T (5 T's)
        members.push(member_of(&template, &[Edit::Sub(151, fix)]));                          // one that shares the template's 7 T's
        members.push(member_of(&template, &[Edit::Del(100, 1), Edit::Sub(151, fix), Edit::Sub(190, other_base(template[190], 1))]));
        assert_eq!(members.len(), 10);
        assert_eq!(consensus_of(&template, &members), truth);
    }
    #[test]
    fn consensus_intron_and_long_deletion_members_do_not_cover_but_a_short_deletion_member_does() {
        // ruling R5. 8 members; 4 lack columns 40..64; of the 4 that have them, 2 carry a 5 bp insertion inside (at 50)
        let template = rand_seq(100, 61);
        let mut plus5 = template.clone(); plus5.splice(50..50, b"ACGTA".iter().copied());
        let (carrier, plain) = (|| member_of(&template, &[Edit::Ins(50, b"ACGTA")]), || member_of(&template, &[]));
        let eight = |absent: &dyn Fn() -> (Vec<u8>, PafHit)| vec![carrier(), carrier(), plain(), plain(), absent(), absent(), absent(), absent()];
        // an intron, a 24 bp `-` and a 20 bp `-` all leave the 4 members out of the coverage: 2 of the 4 covering members carry the insertion (50%)
        let spliced = || { let (s, mut h) = member_of(&template, &[Edit::Del(40, 24)]); h.cs = Some(":40~gt24ag:36".into()); (s, h) };
        assert_eq!(consensus_of(&template, &eight(&spliced)), plus5);
        assert_eq!(consensus_of(&template, &eight(&|| member_of(&template, &[Edit::Del(40, 24)]))), plus5);
        assert_eq!(consensus_of(&template, &eight(&|| member_of(&template, &[Edit::Del(40, 20)]))), plus5);
        // a 19 bp deletion is error-sized: those 4 members cover the columns and delete them (4 of 8 = 50% removes 40..59), and the insertion
        // has only 2 of 8 (25%)
        let mut minus19 = template.clone(); minus19.drain(40..59);
        assert_eq!(consensus_of(&template, &eight(&|| member_of(&template, &[Edit::Del(40, 19)]))), minus19);
    }
    #[test]
    fn consensus_members_that_skip_an_exon_do_not_vote_inside_it_and_the_carriers_decide() {
        // ruling R5: a member with a >= 20 bp deletion has no base there (an isoform's absence): it neither covers the columns nor counts as
        // agreeing with the template, so the vote inside the exon is among the members that carry it and needs 3 of them
        let template = rand_seq(120, 53);                                  // exon 40..64, variant column 50
        let alt = other_base(template[50], 0);
        let mut with_alt = template.clone(); with_alt[50] = alt;
        let (skip, variant, plain) = (|| member_of(&template, &[Edit::Del(40, 24)]), || member_of(&template, &[Edit::Sub(50, alt)]), || member_of(&template, &[]));
        assert_eq!(consensus_of(&template, &[variant(), variant(), skip(), skip(), skip()]), template);                        // 3 of 5 skip: 2 carriers are below the floor of 3
        assert_eq!(consensus_of(&template, &[variant(), variant(), plain(), skip(), skip(), skip()]), with_alt);               // 3 of 6 skip: 2 against 1 among the carriers flips it
        assert_eq!(consensus_of(&template, &[variant(), variant(), variant(), skip(), skip(), skip()]), with_alt);             // 3 against 0
        assert_eq!(consensus_of(&template, &[variant(), variant(), plain(), plain(), skip(), skip(), skip()]), template);      // 2 against 2: the template stays
    }
    #[test]
    fn consensus_small_deletions_are_a_plain_majority_of_the_covering_members_even_inside_a_skipped_exon() {
        // ruling R5: no guard for columns that a member skips with a >= 20 bp deletion, and the skippers do not dilute the share
        let template = rand_seq(120, 53);
        let mut minus1 = template.clone(); minus1.remove(50);
        let (skip, del1, plain) = (|| member_of(&template, &[Edit::Del(40, 24)]), || member_of(&template, &[Edit::Del(50, 1)]), || member_of(&template, &[]));
        assert_eq!(consensus_of(&template, &[del1(), del1(), del1(), skip(), skip()]), minus1);          // 3 of the 3 carriers delete the column
        assert_eq!(consensus_of(&template, &[del1(), del1(), plain(), skip(), skip()]), minus1);         // 2 of 3
        assert_eq!(consensus_of(&template, &[del1(), plain(), plain(), skip(), skip()]), template);      // 1 of 3
        assert_eq!(consensus_of(&template, &[del1(), del1(), skip(), skip(), skip()]), template);        // 2 carriers: below the floor of 3
    }
    #[test]
    fn consensus_trims_template_ends_covered_by_fewer_than_two_members() {
        // spec §5.4 as ruled in R5: leading and trailing template columns with coverage < 2 are dropped (the template counts only as a member: here a self hit)
        let template = rand_seq(150, 33);
        let selfhit = || member_of(&template, &[]);
        let window = |ts: usize, te: usize| member_span(&template, ts, te, &[]);
        // the template is 30 bp longer at the 5' end than every other member (they start at 30 and at 45) and 20 bp longer at the 3' end:
        // the consensus starts where the second member starts and stops where the others stop
        assert_eq!(consensus_of(&template, &[selfhit(), window(30, 130), window(45, 130)]), template[30..130].to_vec());
        // coverage of exactly 2 survives, also at the ends ...
        assert_eq!(consensus_of(&template, &[selfhit(), selfhit()]), template);
        assert_eq!(consensus_of(&template, &[selfhit(), window(10, 140)]), template[10..140].to_vec());
        // ... one member covers nothing twice and no member covers nothing: an empty consensus
        assert!(consensus_of(&template, &[selfhit()]).is_empty());
        assert!(consensus_of(&template, &[]).is_empty());
        // only the ENDS are trimmed: a stretch in the middle that a single member covers stays
        assert_eq!(consensus_of(&template, &[selfhit(), window(0, 60), window(90, 150)]), template);
        // an insertion goes with the column it precedes: one before the first surviving column stays, one after the last surviving column (it
        // precedes a trimmed one) goes. (Alignments that begin or end with an insertion: minimap2 does not write them, a vote over `cs` must handle them.)
        let gg = |ts: usize, te: usize, at: usize| member_span(&template, ts, te, &[Edit::Ins(at, b"GG")]);
        let mut want = b"GG".to_vec(); want.extend_from_slice(&template[30..130]);
        assert_eq!(consensus_of(&template, &[selfhit(), gg(30, 130, 0), gg(30, 130, 0), gg(30, 130, 0)]), want);
        assert_eq!(consensus_of(&template, &[gg(20, 120, 100), gg(20, 120, 100), gg(20, 120, 100)]), template[20..120].to_vec());
    }
    #[test]
    fn consensus_places_a_member_that_starts_inside_the_template_at_its_ts() {
        let template = rand_seq(60, 5);
        let alt = other_base(template[30], 0);
        let mut with_alt = template.clone(); with_alt[30] = alt;
        let partial = || { let (s, mut h) = member_of(&template[25..], &[Edit::Sub(5, alt)]); h.ts = 25; h.te = 60; h.tlen = 60; (s, h) };   // cs `:5*xy:29` from ts = 25
        // two full-length members keep every column covered twice; the three carriers of the variant at 30 are one of them and the two partial members
        assert_eq!(consensus_of(&template, &[member_of(&template, &[Edit::Sub(30, alt)]), member_of(&template, &[]), partial(), partial()]), with_alt);
        // the two partial members alone cover 25..60 twice (below the vote's 3): the consensus is that stretch, unvoted, starting at their ts
        assert_eq!(consensus_of(&template, &[partial(), partial()]), template[25..].to_vec());
    }
    #[test]
    fn consensus_ignores_members_without_a_cs_and_rejects_a_malformed_one() {
        let template = rand_seq(40, 91);
        let alt = other_base(template[10], 0);
        let mut with_alt = template.clone(); with_alt[10] = alt;
        let (s, h) = member_of(&template, &[Edit::Sub(10, alt)]);
        let none = PafHit { cs: None, ..h.clone() };
        let read = &s[..];
        let mut three: Vec<(&[u8], &PafHit)> = vec![(read, &h), (read, &h), (read, &h)];
        three.extend([(read, &none), (read, &none), (read, &none), (read, &none), (read, &none)]);                  // five members with no alignment string abstain
        assert_eq!(consensus_from_template(&template, &three).unwrap(), with_alt);
        let two: Vec<(&[u8], &PafHit)> = vec![(read, &h), (read, &h), (read, &none), (read, &none), (read, &none)];
        assert_eq!(consensus_from_template(&template, &two).unwrap(), template);                                    // only 2 members count
        let bad = PafHit { q: "read7".into(), cs: Some(":10*a".into()), ..h.clone() };
        let err = consensus_from_template(&template, &[(read, &h), (read, &bad)]).unwrap_err().to_string();
        assert!(err.contains("read7") && err.contains("malformed cs"), "the error names the member and the cs problem: {err}");
    }
    #[test]
    fn consensus_survives_an_empty_template_and_a_cs_longer_than_the_template_and_returns_upper_case() {
        let hit = |cs: &str, ts: usize| PafHit { q: "m".into(), qlen: 0, qs: 0, qe: 0, strand: b'+', t: "t".into(), tlen: 0, ts, te: 0, matches: 0, block: 0, de: 0.0, cs: Some(cs.into()) };
        let read: &[u8] = b"";
        // an empty template: an insertion-only cs from three members must not panic (a bare cover[n - 1] would underflow)
        let ins = hit("+acgtacgtacgtacgtacgtacgt", 0);
        assert_eq!(consensus_from_template(b"", &[(read, &ins), (read, &ins), (read, &ins)]).unwrap(), b"".to_vec());
        assert_eq!(consensus_from_template(b"", &[]).unwrap(), b"".to_vec());
        // a cs that runs past the end of the template (garbage in) is clamped: nothing beyond the end is voted or appended
        let template = rand_seq(24, 5);
        let over = hit(":50*ac:10-acg+tt", 0);
        assert_eq!(consensus_from_template(&template, &[(read, &over), (read, &over), (read, &over)]).unwrap(), template);
        let after = hit(":5+aaaa", 30);                                                                      // members that start past the end cover nothing ...
        assert!(consensus_from_template(&template, &[(read, &after), (read, &after), (read, &after)]).unwrap().is_empty());   // ... and nothing survives the end trim
        let whole = hit(":9", 0);
        assert_eq!(consensus_from_template(b"acgtnacgt", &[(read, &whole), (read, &whole)]).unwrap(), b"ACGTNACGT".to_vec());   // lower-case template: upper-case consensus
    }

    #[test]
    fn significance_keeps_supported_variants_and_merges_singletons() {
        assert!(variant_is_real(5, 200, 3, 0.001, 0.01));     // 5 reads carrying 3 distinguishing bases: not error
        assert!(!variant_is_real(2, 200, 1, 0.001, 0.01));    // 2 reads, 1 column: P(X>=2) with p=0.001, n=202 ~ 0.018 >= 0.01 -> merge
        assert!(!variant_is_real(3, 10, 0, 0.001, 0.01));
    }
    #[test]
    fn significance_is_the_binomial_tail_of_the_small_cluster() {
        // exact P(X >= n_small), X ~ Binomial(n_small + n_large, eps^k), eps = 0.001, from rational arithmetic outside this file: the alpha pairs
        // bracket each value, so the p-value is pinned to about three digits (including n = 1000, where the log-factorial table matters)
        let real = |ns: usize, nl: usize, k: usize, alpha: f64| variant_is_real(ns, nl, k, 0.001, alpha);
        assert!(!real(2, 200, 1, 0.0177) && real(2, 200, 1, 0.0178));          // 0.0177860
        assert!(!real(1, 200, 1, 0.1821) && real(1, 200, 1, 0.1822));          // 0.1821698
        assert!(!real(3, 200, 1, 0.00118) && real(3, 200, 1, 0.00119));        // 0.00118318
        assert!(!real(2, 200, 2, 2.02e-8) && real(2, 200, 2, 2.04e-8));        // eps^2 = 1e-6: 2.02983e-8
        assert!(!real(30, 970, 1, 9.4e-34) && real(30, 970, 1, 9.6e-34));      // n = 1000: 9.5031e-34
        // no distinguishing column or no reads is never a real variant; a vanishing p (underflow of eps^k) is real; alpha = 0 admits nothing
        assert!(!real(5, 200, 0, 0.5) && !real(0, 200, 3, 0.5));
        assert!(real(2, 200, 400, 0.01) && !real(2, 200, 400, 0.0));
        // more supporting reads or more columns only make a variant more real
        assert!(!real(2, 200, 1, 0.01) && real(3, 200, 1, 0.01) && real(2, 200, 2, 0.01));
    }
    #[test]
    fn distinguishing_columns_count_substitutions_only() {
        assert_eq!(distinguishing_columns(&parse_cs(":10*ag:5+acgt:3-tt~gt100ag*ct:2").unwrap()), 2);
        assert_eq!(distinguishing_columns(&parse_cs(":50+aaaa-cc:20").unwrap()), 0);
        assert_eq!(distinguishing_columns(&[]), 0);
    }

    #[test]
    fn refine_splits_the_members_that_do_not_fit_the_consensus() {
        let consensus = rand_seq(3000, 3);
        let hit = |de: f64, span: usize| PafHit { q: "m".into(), qlen: 3000, qs: 0, qe: span, strand: b'+', t: "cons".into(), tlen: 3000, ts: 0, te: span, matches: span, block: span, de, cs: None };
        let d = 0.00958;
        // fits; de exactly delta fits; de above delta splits (also between delta and 2 delta); half of the shorter read fits; under half splits
        let hits = [hit(0.002, 3000), hit(d, 3000), hit(0.021, 3000), hit(0.001, 1500), hit(0.001, 1499), hit(0.015, 3000)];
        let read: &[u8] = b"";
        let members: Vec<(&[u8], &PafHit)> = hits.iter().map(|h| (read, h)).collect();
        assert_eq!(refine_cluster(&consensus, &members, d), (vec![0, 1, 3], vec![2, 4, 5]));
        assert_eq!(refine_cluster(&consensus, &members[..2], d), (vec![0, 1], vec![]));              // nothing to split
        assert_eq!(refine_cluster(&consensus, &[], d), (vec![], vec![]));
    }
}

//! Assignment of tied reads by distinguishing sites, and copy-to-copy homology (o3port.py `site_bases`, `column_bases`,
//! `semi_tied`, `anchor_bases`, `fixed_site`, `site_votes`, `homologous`, `one_site`, `disfavoured`, `rec_strand_base`, `homology`).

use super::pile::{aligned_pairs, cmp_base, oriented, Piles, Prim, Rec};
use super::{Copy, Params, ReadId, MIN_READS, MIN_SUPPORT, SITE_FLANK};
use crate::genome::IndexedFasta;
use anyhow::{Context, Result};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::Path;

/// {read index in the primary's orientation: (reference base, read base, reference position)}.
pub type SiteBases = HashMap<i64, (u8, u8, i64)>;
/// {(reference position, base): reads} on one copy.
pub type ColumnCounts = HashMap<(i64, u8), usize>;

/// The Python's `cache` dict: per-copy column counts of every pile read and of its anchors.
#[derive(Default)]
pub struct SiteCache {
    cols: HashMap<String, ColumnCounts>,
    anchors: HashMap<String, ColumnCounts>,
}

/// The read's bases at each aligned position of `rec`, both bases in the primary's orientation.
pub fn site_bases(rec: &Rec, pr: &Prim, genome: &IndexedFasta) -> SiteBases {
    let (seq, _, lead) = oriented(pr, rec);
    let flip = rec.reverse != pr.reverse;
    let rf = upper(genome.fetch(&rec.chrom, rec.start, rec.end).unwrap_or_default());
    let l = seq.len() as i64;
    let mut out = SiteBases::new();
    for (q, t) in aligned_pairs(rec) {
        let (mut i, mut rb, mut qb) = ((lead + q) as i64, rf[(t - rec.start) as usize], seq[lead + q]);
        if flip {
            i = l - 1 - i;
            rb = cmp_base(rb);
            qb = cmp_base(qb);
        }
        out.insert(i, (rb, qb, t));
    }
    out
}

/// `column_bases`: base counts at every position of copy c over its pile's reads, in each record's orientation.
pub fn column_bases<'c>(piles: &Piles, c: &str, prim: &HashMap<ReadId, Prim>, cache: &'c mut SiteCache) -> &'c ColumnCounts {
    if !cache.cols.contains_key(c) {
        cache.cols.insert(c.to_string(), counts(piles, c, prim, |_| true));
    }
    &cache.cols[c]
}

/// `anchor_bases`: as `column_bases`, over c's anchors (reads in no other copy's pile).
pub fn anchor_bases<'c>(piles: &Piles, c: &str, prim: &HashMap<ReadId, Prim>, cache: &'c mut SiteCache) -> &'c ColumnCounts {
    if !cache.anchors.contains_key(c) {
        let others: BTreeSet<ReadId> =
            if piles.len() > 1 { piles.iter().filter(|(d, _)| d.as_str() != c).flat_map(|(_, p)| p.keys().copied()).collect() } else { BTreeSet::new() };
        cache.anchors.insert(c.to_string(), counts(piles, c, prim, |n| !others.contains(&n)));
    }
    &cache.anchors[c]
}

/// The two alignments compare homologous reference bases at read position i (gap-free over i +- SITE_FLANK in both).
pub fn homologous(i: i64, a: &SiteBases, b: &SiteBases) -> bool {
    for x in [a, b] {
        let t = x[&i].2;
        for k in -SITE_FLANK..=SITE_FLANK {
            match x.get(&(i + k)) {
                Some(y) if (y.2 - t).abs() == k.abs() => {}
                _ => return false,
            }
        }
    }
    true
}

/// A base in the primary's orientation, back in this record's orientation.
pub fn rec_strand_base(b: u8, rec: &Rec, pr: &Prim) -> u8 {
    if rec.reverse != pr.reverse {
        cmp_base(b)
    } else {
        b
    }
}

/// --anchor: the site is fixed in the sample (>= MIN_READS anchors of each copy carry its own base, none the other's).
#[allow(clippy::too_many_arguments)]
pub fn fixed_site(ta: i64, ra: u8, rb: u8, tb: i64, piles: &Piles, ca: &str, cb: &str, prim: &HashMap<ReadId, Prim>, cache: &mut SiteCache) -> bool {
    let get = |m: &ColumnCounts, t: i64, b: u8| m.get(&(t, b)).copied().unwrap_or(0);
    let (a_ra, a_rb) = {
        let a = anchor_bases(piles, ca, prim, cache);
        (get(a, ta, ra), get(a, ta, rb))
    };
    let b = anchor_bases(piles, cb, prim, cache);
    a_ra >= MIN_READS && get(b, tb, rb) >= MIN_READS && a_rb == 0 && get(b, tb, ra) == 0
}

/// (sites voting copy ca, sites voting copy cb) for read r, each site counted only when >= MIN_SUPPORT reads of ca's pile carry
/// the read's base there.
#[allow(clippy::too_many_arguments)]
pub fn site_votes(r: ReadId, ca: &str, cb: &str, piles: &Piles, prim: &HashMap<ReadId, Prim>, genome: &IndexedFasta, cache: &mut SiteCache, fixed: bool) -> (i64, i64) {
    let a = site_bases(&piles[ca][&r], &prim[&r], genome);
    let b = site_bases(&piles[cb][&r], &prim[&r], genome);
    column_bases(piles, ca, prim, cache);
    let (mut va, mut vb) = (0i64, 0i64);
    for (i, &(ra, qa, ta)) in &a {
        let Some(&(rb, _, tb)) = b.get(i) else { continue };
        if ra != rb && b"ACGT".contains(&ra) && b"ACGT".contains(&rb) && homologous(*i, &a, &b) {
            if fixed && !fixed_site(ta, ra, rb, tb, piles, ca, cb, prim, cache) {
                continue;
            }
            // the read's base must be an allele of the sample: carried by >= MIN_SUPPORT reads of ca's pile at this column
            let fwd = rec_strand_base(qa, &piles[ca][&r], &prim[&r]);
            if cache.cols[ca].get(&(ta, fwd)).copied().unwrap_or(0) < MIN_SUPPORT {
                continue;
            }
            va += (qa == ra) as i64;
            vb += (qa == rb) as i64;
        }
    }
    (va, vb)
}

fn upper(mut v: Vec<u8>) -> Vec<u8> {
    v.make_ascii_uppercase();
    v
}

/// Base counts at each position of copy c over its pile's reads that `keep` admits (and that have a primary).
fn counts(piles: &Piles, c: &str, prim: &HashMap<ReadId, Prim>, keep: impl Fn(ReadId) -> bool) -> ColumnCounts {
    let mut col = ColumnCounts::new();
    for (&n, rec) in &piles[c] {
        let Some(pr) = prim.get(&n) else { continue };
        if !keep(n) {
            continue;
        }
        let (seq, _, lead) = oriented(pr, rec);
        for (q, t) in aligned_pairs(rec) {
            *col.entry((t, seq[lead + q])).or_default() += 1;
        }
    }
    col
}

/// The read's piles ordered as `sorted(((AS, -NM, copy) ...), reverse=True)`.
fn ranked(r: ReadId, piles: &Piles) -> Vec<(i32, i32, &str)> {
    let mut sc: Vec<(i32, i32, &str)> = piles
        .iter()
        .filter_map(|(c, p)| p.get(&r).map(|rec| (rec.as_.unwrap_or(0), -rec.nm.unwrap_or(0), c.as_str())))
        .collect();
    sc.sort_by(|a, b| b.cmp(a));
    sc
}

/// `semi_tied`: HOPELESS on an exact tie of the two best copies, the leader when it leads by >= MIN_READS supported sites, else
/// ABSTAIN (with --anchor a one-site lead at a fixed site suffices).
pub fn semi_tied(r: ReadId, piles: &Piles, prim: &HashMap<ReadId, Prim>, genome: &IndexedFasta, cache: &mut SiteCache, p: &Params) -> String {
    let sc = ranked(r, piles);
    if sc.len() < 2 || (sc[0].0, sc[0].1) == (sc[1].0, sc[1].1) {
        return "HOPELESS".into();
    }
    if !prim.contains_key(&r) {
        return "ABSTAIN".into();
    }
    let (ca, cb) = (sc[0].2, sc[1].2);
    let (va, vb) = site_votes(r, ca, cb, piles, prim, genome, cache, false);
    if va - vb >= MIN_READS as i64 {
        return ca.into();
    }
    if p.anchor && va - vb >= 1 {
        let (fa, fb) = site_votes(r, ca, cb, piles, prim, genome, cache, true);
        if fa >= 1 && fb == 0 {
            return ca.into();
        }
    }
    "ABSTAIN".into()
}

/// --one-site: a one-site lead suffices when the median lead over the read's haplotype's reads tied to the same copies is >= 1.
#[allow(clippy::too_many_arguments)]
pub fn one_site(
    r: ReadId,
    hap_reads: &BTreeSet<ReadId>,
    piles: &Piles,
    prim: &HashMap<ReadId, Prim>,
    genome: &IndexedFasta,
    cache: &mut SiteCache,
    hap_cache: &mut HashMap<(usize, String, String), bool>,
    h: usize,
) -> String {
    let sc = ranked(r, piles);
    if sc.len() < 2 || !prim.contains_key(&r) {
        return "ABSTAIN".into();
    }
    let (ca, cb) = (sc[0].2, sc[1].2);
    let (va, vb) = site_votes(r, ca, cb, piles, prim, genome, cache, false);
    if va - vb < 1 {
        return "ABSTAIN".into();
    }
    let key = (h, ca.to_string(), cb.to_string());
    if !hap_cache.contains_key(&key) {
        let mut d: Vec<i64> = Vec::new();
        for &q in hap_reads {
            if piles[ca].contains_key(&q) && piles[cb].contains_key(&q) && prim.contains_key(&q) {
                let (x, y) = site_votes(q, ca, cb, piles, prim, genome, cache, false);
                d.push(x - y);
            }
        }
        d.sort();
        let ok = d.len() >= MIN_SUPPORT && d[d.len() / 2] >= 1;
        hap_cache.insert(key.clone(), ok);
    }
    if hap_cache[&key] { ca.into() } else { "ABSTAIN".into() }
}

/// --prefer: the option copies a haplotype's reads clearly disfavour (median lead of another option >= MIN_READS).
pub fn disfavoured(reads: &BTreeSet<ReadId>, opts: &[String], piles: &Piles, prim: &HashMap<ReadId, Prim>, genome: &IndexedFasta, cache: &mut SiteCache) -> BTreeSet<String> {
    let mut out = BTreeSet::new();
    for c in opts {
        for z in opts {
            if c == z || out.contains(z) {
                continue;
            }
            let mut d: Vec<i64> = Vec::new();
            for &r in reads {
                if piles[c].contains_key(&r) && piles[z].contains_key(&r) && prim.contains_key(&r) {
                    let (x, y) = site_votes(r, c, z, piles, prim, genome, cache, false);
                    d.push(x - y);
                }
            }
            d.sort();
            if !d.is_empty() && d[d.len() / 2] >= MIN_READS as i64 {
                out.insert(z.clone());
            }
        }
    }
    out
}

/// {other copy: {home position: (other position, '+' strand)}} from `minimap2 -c --eqx -x asm20 -t 1` of each other copy's locus
/// against the home copy's (best-matching alignment first; a position keeps its first mapping).
pub fn homology(home: &Copy, others: &[&Copy], genome: &IndexedFasta, tmpdir: &Path) -> Result<HashMap<String, HashMap<i64, (i64, bool)>>> {
    use std::io::Write;
    let t = tmpdir.join("t.fa");
    let q = tmpdir.join("q.fa");
    let seq = |c: &Copy| genome.fetch(&c.chrom, c.start, c.end).unwrap_or_default();
    {
        let mut fh = std::fs::File::create(&t).with_context(|| format!("creating {}", t.display()))?;
        fh.write_all(format!(">{}\n", home.tid).as_bytes())?;
        fh.write_all(&seq(home))?;
        fh.write_all(b"\n")?;
        let mut fh = std::fs::File::create(&q).with_context(|| format!("creating {}", q.display()))?;
        for c in others {
            fh.write_all(format!(">{}\n", c.tid).as_bytes())?;
            fh.write_all(&seq(c))?;
            fh.write_all(b"\n")?;
        }
    }
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let out = std::process::Command::new(&mm2)
        .args(["-c", "--eqx", "-x", "asm20", "-t", "1"])
        .arg(&t)
        .arg(&q)
        .stderr(std::process::Stdio::null())
        .output()
        .with_context(|| format!("starting {mm2}"))?;
    let paf = String::from_utf8_lossy(&out.stdout);
    let start: HashMap<&str, i64> = others.iter().map(|c| (c.tid.as_str(), c.start)).collect();
    let mut lines: Vec<Vec<&str>> = paf.lines().map(|l| l.split('\t').collect()).collect();
    lines.sort_by_key(|p| -p[9].parse::<i64>().unwrap_or(0)); // stable: PAF order on a tie
    let mut maps: HashMap<String, HashMap<i64, (i64, bool)>> = HashMap::new();
    for p in &lines {
        let (cid, plus) = (p[0], p[4] == "+");
        let (qs, qe, ts): (i64, i64, i64) = (p[2].parse()?, p[3].parse()?, p[7].parse()?);
        let cg = p[12..].iter().find_map(|x| x.strip_prefix("cg:Z:")).context("a PAF line without cg:Z")?;
        let mp = maps.entry(cid.to_string()).or_default();
        let (mut tpos, mut qo) = (ts, 0i64);
        let mut num = 0i64;
        for ch in cg.bytes() {
            if ch.is_ascii_digit() {
                num = num * 10 + (ch - b'0') as i64;
                continue;
            }
            let n = std::mem::take(&mut num);
            match ch {
                b'=' | b'X' | b'M' => {
                    for k in 0..n {
                        let qq = if plus { qs + qo + k } else { qe - 1 - (qo + k) };
                        mp.entry(home.start + tpos + k).or_insert((start[cid] + qq, plus));
                    }
                    tpos += n;
                    qo += n;
                }
                b'I' => qo += n,
                b'D' => tpos += n,
                _ => {}
            }
        }
    }
    Ok(maps)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::copy_phasing::pile::{Pile, OP_M};
    use std::sync::Arc;

    /// A FASTA (+ .fai) in a fresh temp dir.
    pub(crate) fn fasta(tag: &str, seqs: &[(&str, &str)]) -> (std::path::PathBuf, IndexedFasta) {
        let dir = std::env::temp_dir().join(format!("copy_phasing_{tag}_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("g.fa");
        let (mut text, mut fai, mut off) = (String::new(), String::new(), 0usize);
        for (n, s) in seqs {
            let head = format!(">{n}\n");
            off += head.len();
            fai.push_str(&format!("{n}\t{}\t{off}\t{}\t{}\n", s.len(), s.len(), s.len() + 1));
            text.push_str(&head);
            text.push_str(s);
            text.push('\n');
            off += s.len() + 1;
        }
        std::fs::write(&fa, text).unwrap();
        std::fs::write(dir.join("g.fa.fai"), fai).unwrap();
        let g = IndexedFasta::open(fa.to_str().unwrap()).unwrap();
        (dir, g)
    }

    fn rec(read: ReadId, chrom: &str, len: u32) -> Rec {
        Rec {
            read,
            chrom: Arc::from(chrom),
            start: 0,
            end: len as i64,
            as_: Some(50),
            nm: Some(1),
            reverse: false,
            secondary: false,
            supplementary: false,
            cigar: vec![(OP_M, len)],
            seq: Vec::new(),
        }
    }

    /// Copies a / b (30 bp) differ at position 15 (C / T); read 0 carries b's T and aligns 30M to both; `carriers` more reads of a's
    /// pile carry T there.
    fn toy(tag: &str, carriers: u32) -> (std::path::PathBuf, IndexedFasta, Piles, HashMap<ReadId, Prim>) {
        let a = "ACGTACGTACGTACGCACGTACGTACGTAC";
        let b = "ACGTACGTACGTACGTACGTACGTACGTAC";
        let (dir, g) = fasta(&format!("{tag}{carriers}"), &[("a", a), ("b", b)]);
        let mut pa = Pile::new();
        let mut pb = Pile::new();
        let mut prim = HashMap::new();
        for k in 0..=carriers {
            pa.insert(k, rec(k, "a", 30));
            prim.insert(k, Prim::new(b.as_bytes().to_vec(), None, false, 50));
        }
        pb.insert(0, rec(0, "b", 30));
        let piles: Piles = [("a".to_string(), pa), ("b".to_string(), pb)].into_iter().collect();
        (dir, g, piles, prim)
    }

    #[test]
    fn site_votes_counts_supported_sites_only() {
        let (d, g, piles, prim) = toy("votes", 6); // 7 reads of a's pile carry T: supported
        assert_eq!(site_votes(0, "a", "b", &piles, &prim, &g, &mut SiteCache::default(), false), (0, 1));
        let _ = std::fs::remove_dir_all(d);
        let (d, g, piles, prim) = toy("votes", 4); // 5 carriers < MIN_SUPPORT
        assert_eq!(site_votes(0, "a", "b", &piles, &prim, &g, &mut SiteCache::default(), false), (0, 0));
        let _ = std::fs::remove_dir_all(d);
    }

    #[test]
    fn semi_tied_exact_tie_is_hopeless_and_a_one_site_lead_abstains() {
        let (d, g, mut piles, prim) = toy("semi", 6);
        assert_eq!(semi_tied(0, &piles, &prim, &g, &mut SiteCache::default(), &Params::default()), "HOPELESS");
        piles.get_mut("b").unwrap().get_mut(&0).unwrap().nm = Some(0); // b leads by edits, by one site only
        assert_eq!(semi_tied(0, &piles, &prim, &g, &mut SiteCache::default(), &Params::default()), "ABSTAIN");
        assert_eq!(semi_tied(9, &piles, &prim, &g, &mut SiteCache::default(), &Params::default()), "HOPELESS");
        let _ = std::fs::remove_dir_all(d);
    }

    fn line(shift_at: i64) -> SiteBases {
        (0..30).map(|i| (i, (b'A', b'A', if i >= shift_at { i + 1 } else { i }))).collect()
    }

    #[test]
    fn homologous_needs_flank() {
        let a = line(100);
        assert!(!homologous(15, &a, &line(18))); // an indel 3 bp from the site
        assert!(homologous(15, &a, &line(21))); // 6 bp: outside SITE_FLANK
        assert!(!homologous(2, &a, &a)); // too close to the read's start
    }

    #[test]
    fn homology_no_hit() {
        // a 2 kb pseudo-random sequence (an LCG), and a poly-T one
        let mut st = 12345u64;
        let x: String = (0..2000)
            .map(|_| {
                st = st.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(st >> 62) as usize] as char
            })
            .collect();
        let y = "T".repeat(2000);
        let (x, y) = (x.as_str(), y.as_str());
        let (d, g) = fasta("hom", &[("x", x), ("y", y)]);
        let hx = Copy { tid: "X".into(), chrom: "x".into(), start: 0, end: x.len() as i64 };
        let hy = Copy { tid: "Y".into(), chrom: "y".into(), start: 0, end: y.len() as i64 };
        let m = homology(&hx, &[&hy], &g, &d).unwrap();
        assert!(m.is_empty());
        // a copy against itself maps every position to itself on '+'
        let m = homology(&hx, &[&Copy { tid: "X2".into(), ..hx.clone() }], &g, &d).unwrap();
        assert_eq!(m["X2"][&10], (10, true));
        let _ = std::fs::remove_dir_all(d);
    }
}

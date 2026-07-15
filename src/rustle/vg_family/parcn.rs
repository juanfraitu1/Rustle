//! Assembly-based paralog-specific copy number (parCN). OPTIONAL assembly/DNA-side supplement:
//! projects catalog copy consensuses onto phased haplotype assemblies and counts per-copy genomic loci,
//! disambiguated by deterministic SUN witnesses. Consumes only copies.fa + the assemblies; never wires
//! into the RNA-exclusive core. See docs/superpowers/specs/2026-07-14-assembly-parcn-design.md.

use std::collections::BTreeMap;

use crate::vg_family::copy_assign_pipeline::banded_msa_pair;

#[derive(Clone, Debug)]
pub struct Copy {
    pub family_id: String,
    pub copy_id: String,
    pub seq: Vec<u8>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Tier { T1, T2, T3, NA }

#[derive(Clone, Debug)]
pub struct CopySun { pub copy_id: String, pub tier: Tier, pub private: Vec<(usize, u8)> }

/// For copy `b` vs sibling `s`: the SET of offsets in `b` (non-gap) whose aligned base in `s` differs
/// (substitution or gap), plus `(matches, aligned_cols)` for identity. Uses the banded 2-row MSA; if the
/// pair can't be aligned in-band, a difference from `s` cannot be CONFIRMED for any offset, so this
/// returns an EMPTY diff set. Callers intersect diff sets across siblings to build a copy's private-position
/// set (`private = ∩ diff_offsets(b, s_j)`); an empty set is absorbing under intersection (`∅ ∩ X = ∅`), so a
/// failed comparison correctly removes ALL candidates from privateness (conservative: nothing private).
/// Returning "every offset differs" instead would be a no-op under intersection (`U ∩ X = X`) and could
/// fabricate a spurious Tier-1 (SUN) call from zero evidence — the bug this fallback fixes.
fn diff_offsets(b: &[u8], s: &[u8], band: usize) -> (std::collections::HashSet<usize>, usize, usize) {
    let msa = match banded_msa_pair(b, s, band) {
        Some(m) => m,
        None => return (std::collections::HashSet::new(), 0, b.len().max(1)),
    };
    let (ab, asb) = (&msa[0], &msa[1]);
    let (mut boff, mut diff, mut matches, mut cols) = (0usize, std::collections::HashSet::new(), 0usize, 0usize);
    for k in 0..ab.len() {
        let (cb, cs) = (ab[k], asb[k]);
        if cb != b'-' {
            cols += 1;
            if cs == cb { matches += 1; } else { diff.insert(boff); }
            boff += 1;
        }
    }
    (diff, matches, cols.max(1))
}

/// Per-copy private positions (a SUN = an offset in copy B whose base differs from EVERY sibling) + tier.
/// Band scales with the family's copy-length spread. Threshold-free: T1 iff ≥1 private position; T3 iff a
/// sibling is ≥99.9% identical (indistinguishable); T2 otherwise; NA for a single-copy family.
pub fn sun_positions(copies: &[Copy], band: usize) -> Vec<CopySun> {
    let mut out = Vec::with_capacity(copies.len());
    for (i, b) in copies.iter().enumerate() {
        if copies.len() == 1 {
            out.push(CopySun { copy_id: b.copy_id.clone(), tier: Tier::NA, private: Vec::new() });
            continue;
        }
        // private = offsets differing from ALL siblings; max_id = closest sibling identity.
        let mut private_set: Option<std::collections::HashSet<usize>> = None;
        let mut max_id = 0.0f64;
        for (j, s) in copies.iter().enumerate() {
            if i == j { continue; }
            let (diff, matches, cols) = diff_offsets(&b.seq, &s.seq, band);
            max_id = max_id.max(matches as f64 / cols as f64);
            private_set = Some(match private_set { None => diff, Some(acc) => acc.intersection(&diff).copied().collect() });
        }
        let private_set = private_set.unwrap_or_default();
        let mut private: Vec<(usize, u8)> = private_set.iter().map(|&p| (p, b.seq[p])).collect();
        private.sort_unstable();
        let tier = if !private.is_empty() { Tier::T1 } else if max_id >= 0.999 { Tier::T3 } else { Tier::T2 };
        out.push(CopySun { copy_id: b.copy_id.clone(), tier, private });
    }
    out
}

/// Parse a `gw_family_catalog` `copies.fa` (`>{family}|{copy_idx}|{chrom}:{s}-{e}|{strand}|nexon={n}`) into
/// families → copies, preserving file order. Sequence lines are concatenated and upper-cased.
pub fn parse_copies_fa(path: &str) -> anyhow::Result<BTreeMap<String, Vec<Copy>>> {
    let text = std::fs::read_to_string(path)?;
    let mut fams: BTreeMap<String, Vec<Copy>> = BTreeMap::new();
    let (mut fam, mut cid, mut seq) = (String::new(), String::new(), Vec::<u8>::new());
    let mut have = false;
    let flush = |fams: &mut BTreeMap<String, Vec<Copy>>, fam: &str, cid: &str, seq: &mut Vec<u8>| {
        if !fam.is_empty() {
            fams.entry(fam.to_string()).or_default().push(Copy { family_id: fam.to_string(), copy_id: cid.to_string(), seq: std::mem::take(seq) });
        }
    };
    for line in text.lines() {
        if let Some(h) = line.strip_prefix('>') {
            if have { flush(&mut fams, &fam, &cid, &mut seq); }
            let mut it = h.split('|');
            fam = it.next().unwrap_or("").to_string();
            cid = it.next().unwrap_or("0").to_string();
            have = true;
        } else {
            seq.extend(line.trim().bytes().map(|b| b.to_ascii_uppercase()));
        }
    }
    if have { flush(&mut fams, &fam, &cid, &mut seq); }
    Ok(fams)
}

/// Walk a minimap2 `cs:Z:` short string and return, per requested QUERY offset, the aligned TARGET
/// (assembly) base. cs ops: `:N` = N matches (target base == query base); `*xy` = substitution, x = target
/// base, y = query base (advance query 1); `+seq` = insertion, query-only (target base None); `-seq` =
/// deletion, target-only (no query advance); `~` splice = target-only. Bases are returned upper-cased.
pub fn cs_bases_at(cs: &str, query: &[u8], positions: &[usize]) -> Vec<Option<u8>> {
    let want: std::collections::HashSet<usize> = positions.iter().copied().collect();
    let mut base_at: std::collections::HashMap<usize, Option<u8>> = std::collections::HashMap::new();
    let bytes = cs.as_bytes();
    let (mut k, mut qoff) = (0usize, 0usize);
    while k < bytes.len() {
        match bytes[k] {
            b':' => {
                let mut j = k + 1; let mut n = 0usize;
                while j < bytes.len() && bytes[j].is_ascii_digit() { n = n * 10 + (bytes[j] - b'0') as usize; j += 1; }
                for _ in 0..n {
                    if want.contains(&qoff) { base_at.insert(qoff, query.get(qoff).map(|b| b.to_ascii_uppercase())); }
                    qoff += 1;
                }
                k = j;
            }
            b'*' => {
                // *<target><query>
                let tgt = bytes.get(k + 1).map(|b| b.to_ascii_uppercase());
                if want.contains(&qoff) { base_at.insert(qoff, tgt); }
                qoff += 1;
                k += 3;
            }
            b'+' => {
                let mut j = k + 1;
                while j < bytes.len() && bytes[j].is_ascii_alphabetic() {
                    if want.contains(&qoff) { base_at.insert(qoff, None); }
                    qoff += 1; j += 1;
                }
                k = j;
            }
            b'-' | b'~' => { // target-only: skip the following letters/coords, no query advance
                let mut j = k + 1;
                while j < bytes.len() && bytes[j] != b':' && bytes[j] != b'*' && bytes[j] != b'+' && bytes[j] != b'-' && bytes[j] != b'~' { j += 1; }
                k = j;
            }
            _ => { k += 1; }
        }
    }
    positions.iter().map(|p| base_at.get(p).copied().flatten()).collect()
}

#[derive(Clone, Debug, PartialEq)]
pub enum Method { Sun, AlignFallback, Unresolved, SingleCopy }

#[derive(Clone, Debug)]
pub struct Locus {
    pub chrom: String, pub start: u64, pub end: u64,
    pub best_copy: String, pub identity: f64, pub runner_up_identity: f64, pub cs: String,
}

#[derive(Clone, Debug)]
pub struct Assignment { pub copy_id: Option<String>, pub method: Method }

const ALIGN_MARGIN: f64 = 0.002;

/// Hybrid assignment of a projected locus to its best copy. Tier-1: confirm the assembly carries the best
/// copy's private base (via cs) at ≥1 private position → deterministic SUN. Tier-2: assign to the best copy
/// iff its identity beats the runner-up by ≥ ALIGN_MARGIN (flagged fallback). Tier-3 / private-not-confirmed
/// / near-tie → UNRESOLVED. NA (single copy) → single_copy.
pub fn assign_locus(locus: &Locus, sun: &CopySun, best_copy_seq: &[u8]) -> Assignment {
    match sun.tier {
        Tier::NA => Assignment { copy_id: Some(locus.best_copy.clone()), method: Method::SingleCopy },
        Tier::T1 => {
            let positions: Vec<usize> = sun.private.iter().map(|&(p, _)| p).collect();
            let bases = cs_bases_at(&locus.cs, best_copy_seq, &positions);
            let confirmed = sun.private.iter().zip(bases.iter())
                .any(|(&(_, pb), got)| *got == Some(pb));
            if confirmed { Assignment { copy_id: Some(locus.best_copy.clone()), method: Method::Sun } }
            else { Assignment { copy_id: None, method: Method::Unresolved } }
        }
        Tier::T2 => {
            if locus.identity - locus.runner_up_identity >= ALIGN_MARGIN {
                Assignment { copy_id: Some(locus.best_copy.clone()), method: Method::AlignFallback }
            } else { Assignment { copy_id: None, method: Method::Unresolved } }
        }
        Tier::T3 => Assignment { copy_id: None, method: Method::Unresolved },
    }
}

#[derive(Clone, Debug)]
pub struct ParcnRow { pub family_id: String, pub copy_id: String, pub tier: Tier, pub loci_mat: usize, pub loci_pat: usize, pub method: Method }

fn tier_str(t: &Tier) -> &'static str { match t { Tier::T1 => "T1", Tier::T2 => "T2", Tier::T3 => "T3", Tier::NA => "NA" } }
fn method_str(m: &Method) -> &'static str { match m { Method::Sun => "SUN", Method::AlignFallback => "align_fallback", Method::Unresolved => "UNRESOLVED", Method::SingleCopy => "single_copy" } }

/// Count per-copy assigned loci (mat/pat), pick each copy's dominant assignment method, and total the
/// unresolved loci across both haplotypes. Copies with no assigned locus still get a row (parCN 0).
pub fn tabulate(family_id: &str, copies: &[Copy], suns: &[CopySun], mat: &[Assignment], pat: &[Assignment]) -> (Vec<ParcnRow>, usize) {
    use std::collections::HashMap;
    let tier_of: HashMap<&str, &Tier> = suns.iter().map(|s| (s.copy_id.as_str(), &s.tier)).collect();
    let mut mat_c: HashMap<String, usize> = HashMap::new();
    let mut pat_c: HashMap<String, usize> = HashMap::new();
    let mut method_of: HashMap<String, Method> = HashMap::new();
    let mut n_unres = 0usize;
    for (side, counts) in [(mat, &mut mat_c), (pat, &mut pat_c)] {
        for a in side {
            match &a.copy_id {
                Some(cp) => { *counts.entry(cp.clone()).or_insert(0) += 1; method_of.entry(cp.clone()).or_insert_with(|| a.method.clone()); }
                None => n_unres += 1,
            }
        }
    }
    let mut rows = Vec::with_capacity(copies.len());
    for c in copies {
        let tier = (*tier_of.get(c.copy_id.as_str()).unwrap_or(&&Tier::NA)).clone();
        let method = method_of.get(&c.copy_id).cloned().unwrap_or(Method::Unresolved);
        rows.push(ParcnRow { family_id: family_id.to_string(), copy_id: c.copy_id.clone(), tier,
            loci_mat: *mat_c.get(&c.copy_id).unwrap_or(&0), loci_pat: *pat_c.get(&c.copy_id).unwrap_or(&0), method });
    }
    (rows, n_unres)
}

pub fn format_parcn_row(r: &ParcnRow) -> String {
    format!("{}\t{}\t{}\t{}\t{}\t{}\t{}", r.family_id, r.copy_id, tier_str(&r.tier), r.loci_mat, r.loci_pat, r.loci_mat + r.loci_pat, method_str(&r.method))
}
pub fn format_family_row(family_id: &str, rows: &[ParcnRow], n_unresolved: usize) -> String {
    let famcn: usize = rows.iter().map(|r| r.loci_mat + r.loci_pat).sum();
    format!("{}\t{}\t{}\t{}", family_id, rows.len(), famcn, n_unresolved)
}

fn recip_overlap(a: &Locus, b: &Locus) -> f64 {
    if a.chrom != b.chrom { return 0.0; }
    let (lo, hi) = (a.start.max(b.start), a.end.min(b.end));
    if hi <= lo { return 0.0; }
    let ov = (hi - lo) as f64;
    let la = (a.end - a.start).max(1) as f64;
    let lb = (b.end - b.start).max(1) as f64;
    (ov / la).min(ov / lb)
}

/// Collapse reciprocal-overlap ≥ 0.50 loci into one, keeping the highest-identity member (its best_copy)
/// and recording the next-highest overlapping identity as runner_up_identity (for the Tier-2 margin gate).
pub fn dedup_loci(mut loci: Vec<Locus>) -> Vec<Locus> {
    loci.sort_by(|a, b| b.identity.partial_cmp(&a.identity).unwrap_or(std::cmp::Ordering::Equal));
    let mut kept: Vec<Locus> = Vec::new();
    for l in loci {
        if let Some(k) = kept.iter_mut().find(|k| recip_overlap(k, &l) >= 0.50) {
            // loci are sorted DESC by identity, so the kept member `k` already has the higher identity and
            // `l` is a runner-up candidate for that overlap group; record the highest runner-up seen.
            k.runner_up_identity = k.runner_up_identity.max(l.identity);
        } else {
            kept.push(l);
        }
    }
    kept
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_copies_fa_groups_by_family() {
        let dir = std::env::temp_dir();
        let p = dir.join(format!("parcn_copies_{}.fa", std::process::id()));
        std::fs::write(&p, ">RBMY|0|chrY:1-9|+|nexon=3\nACGTACGTA\n>RBMY|1|chrY:20-28|+|nexon=3\nACGTACGTT\n>DAZ|0|chrY:99-104|-|nexon=1\nGGGCCC\n").unwrap();
        let fams = parse_copies_fa(p.to_str().unwrap()).unwrap();
        std::fs::remove_file(&p).ok();
        assert_eq!(fams.len(), 2);
        assert_eq!(fams["RBMY"].len(), 2);
        assert_eq!(fams["RBMY"][0].copy_id, "0");
        assert_eq!(fams["RBMY"][1].seq, b"ACGTACGTT");
        assert_eq!(fams["DAZ"][0].seq, b"GGGCCC");
    }

    #[test]
    fn sun_positions_finds_private_snv_and_tiers() {
        // copy0 vs copy1 differ ONLY at offset 4 (A vs T) -> each has a private position there (Tier-1).
        // copy2 is identical to copy0 -> Tier-3 (indistinguishable), and offset 4 is no longer private to copy0.
        let copies = vec![
            Copy { family_id: "F".into(), copy_id: "0".into(), seq: b"ACGTAGGTCA".to_vec() },
            Copy { family_id: "F".into(), copy_id: "1".into(), seq: b"ACGTTGGTCA".to_vec() },
            Copy { family_id: "F".into(), copy_id: "2".into(), seq: b"ACGTAGGTCA".to_vec() },
        ];
        let suns = sun_positions(&copies, 8);
        let s0 = suns.iter().find(|s| s.copy_id == "0").unwrap();
        let s1 = suns.iter().find(|s| s.copy_id == "1").unwrap();
        let s2 = suns.iter().find(|s| s.copy_id == "2").unwrap();
        // copy1's 'T' at offset 4 is unique among the three -> private -> Tier-1.
        assert_eq!(s1.tier, Tier::T1);
        assert!(s1.private.iter().any(|&(p, b)| p == 4 && b == b'T'));
        // copy0 and copy2 are identical -> neither can have a private position -> Tier-3.
        assert_eq!(s0.tier, Tier::T3);
        assert_eq!(s2.tier, Tier::T3);
        assert!(s0.private.is_empty());
    }

    #[test]
    fn sun_positions_single_copy_is_na() {
        let copies = vec![Copy { family_id: "F".into(), copy_id: "0".into(), seq: b"ACGTACGT".to_vec() }];
        let suns = sun_positions(&copies, 8);
        assert_eq!(suns[0].tier, Tier::NA);
    }

    #[test]
    fn cs_bases_reads_match_and_substitution_and_insertion() {
        // query ACGTACGT (len 8). cs: 3 matches, sub (target g / query t) at q=3, 2 matches,
        // insertion of "AA" at q=6..8, tail is target-only del (does not advance query).
        // cs grammar: :N match run; *<tgt><qry> substitution; +<seq> insertion (query-only); -<seq> deletion (target-only).
        let cs = ":3*gt:2+aa-cc";
        let q = b"ACGTACGT";
        // q0 match -> assembly base 'A'(=query); q3 substitution -> assembly base 'G'(target, upper); q5 match 'C'; q6 insertion -> None.
        let got = cs_bases_at(cs, q, &[0, 3, 5, 6]);
        assert_eq!(got, vec![Some(b'A'), Some(b'G'), Some(b'C'), None]);
    }

    #[test]
    fn sun_positions_band_edge_yields_no_private() {
        // length difference (8 vs 12) exceeds a tiny band -> banded_msa_pair returns None ->
        // the conservative fallback must yield NO private positions (never a spurious Tier-1).
        let copies = vec![
            Copy { family_id: "F".into(), copy_id: "0".into(), seq: b"ACGTACGT".to_vec() },
            Copy { family_id: "F".into(), copy_id: "1".into(), seq: b"ACGTACGTACGT".to_vec() },
        ];
        let suns = sun_positions(&copies, 1); // band=1 << |8-12|
        for s in &suns {
            assert!(s.private.is_empty(), "band-edge failure must not fabricate private positions");
            assert_ne!(s.tier, super::Tier::T1, "a copy that could not be compared must not be Tier-1");
        }
    }

    #[test]
    fn assign_locus_hybrid_tiers() {
        let mk = |cs: &str, id: f64, ru: f64| Locus { chrom: "c".into(), start: 0, end: 9, best_copy: "0".into(), identity: id, runner_up_identity: ru, cs: cs.into() };
        let seq = b"ACGTAGGTCA"; // best copy's private base is 'A' at offset 4
        let sun_t1 = CopySun { copy_id: "0".into(), tier: Tier::T1, private: vec![(4, b'A')] };
        // cs shows a MATCH across offset 4 -> assembly carries 'A' -> SUN confirmed.
        assert_eq!(assign_locus(&mk(":10", 0.99, 0.90), &sun_t1, seq).method, Method::Sun);
        // cs shows a substitution at offset 4 (target g) -> assembly does NOT carry 'A' -> UNRESOLVED.
        assert_eq!(assign_locus(&mk(":4*ga:5", 0.99, 0.90), &sun_t1, seq).method, Method::Unresolved);
        // Tier-2 with a clear identity margin -> align_fallback.
        let sun_t2 = CopySun { copy_id: "0".into(), tier: Tier::T2, private: vec![] };
        assert_eq!(assign_locus(&mk(":10", 0.99, 0.90), &sun_t2, seq).method, Method::AlignFallback);
        // Tier-2 near-tie -> UNRESOLVED.
        assert_eq!(assign_locus(&mk(":10", 0.991, 0.990), &sun_t2, seq).method, Method::Unresolved);
        // Tier-3 -> UNRESOLVED. NA -> single_copy.
        let sun_t3 = CopySun { copy_id: "0".into(), tier: Tier::T3, private: vec![] };
        assert_eq!(assign_locus(&mk(":10", 0.99, 0.0), &sun_t3, seq).method, Method::Unresolved);
        let sun_na = CopySun { copy_id: "0".into(), tier: Tier::NA, private: vec![] };
        assert_eq!(assign_locus(&mk(":10", 0.99, 0.0), &sun_na, seq).method, Method::SingleCopy);
    }

    #[test]
    fn tabulate_counts_and_formats() {
        let copies = vec![
            Copy { family_id: "RBMY".into(), copy_id: "0".into(), seq: b"AAAA".to_vec() },
            Copy { family_id: "RBMY".into(), copy_id: "1".into(), seq: b"AAAT".to_vec() },
        ];
        let suns = vec![
            CopySun { copy_id: "0".into(), tier: Tier::T1, private: vec![(3, b'A')] },
            CopySun { copy_id: "1".into(), tier: Tier::T2, private: vec![] },
        ];
        let a = |cp: &str, m: Method| Assignment { copy_id: Some(cp.into()), method: m };
        let un = || Assignment { copy_id: None, method: Method::Unresolved };
        // mat: copy0 SUN once, one unresolved. pat: copy0 SUN once, copy1 fallback once.
        let mat = vec![a("0", Method::Sun), un()];
        let pat = vec![a("0", Method::Sun), a("1", Method::AlignFallback)];
        let (rows, n_unres) = tabulate("RBMY", &copies, &suns, &mat, &pat);
        let r0 = rows.iter().find(|r| r.copy_id == "0").unwrap();
        assert_eq!((r0.loci_mat, r0.loci_pat), (1, 1));       // parCN 2
        let r1 = rows.iter().find(|r| r.copy_id == "1").unwrap();
        assert_eq!((r1.loci_mat, r1.loci_pat), (0, 1));       // parCN 1
        assert_eq!(n_unres, 1);
        assert_eq!(format_parcn_row(r0), "RBMY\t0\tT1\t1\t1\t2\tSUN");
        assert_eq!(format_family_row("RBMY", &rows, n_unres), "RBMY\t2\t3\t1");
    }

    #[test]
    fn dedup_collapses_overlapping_keeps_best() {
        let mk = |s: u64, e: u64, cp: &str, id: f64| Locus { chrom: "c1".into(), start: s, end: e, best_copy: cp.into(), identity: id, runner_up_identity: 0.0, cs: ":1".into() };
        // Two heavily-overlapping hits (copy0 id .99, copy1 id .97) + one disjoint locus.
        let loci = vec![mk(1000, 2000, "0", 0.99), mk(1010, 1990, "1", 0.97), mk(50000, 51000, "3", 0.98)];
        let mut out = dedup_loci(loci);
        out.sort_by_key(|l| l.start);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].best_copy, "0");                 // highest identity wins
        assert!((out[0].runner_up_identity - 0.97).abs() < 1e-9); // runner-up recorded
        assert_eq!(out[1].best_copy, "3");
    }
}

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
/// pair can't be aligned in-band, treats EVERY b offset as differing (conservative: nothing private).
fn diff_offsets(b: &[u8], s: &[u8], band: usize) -> (std::collections::HashSet<usize>, usize, usize) {
    let msa = match banded_msa_pair(b, s, band) {
        Some(m) => m,
        None => return ((0..b.len()).collect(), 0, b.len().max(1)),
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
}

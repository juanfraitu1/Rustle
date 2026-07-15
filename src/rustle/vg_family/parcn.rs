//! Assembly-based paralog-specific copy number (parCN). OPTIONAL assembly/DNA-side supplement:
//! projects catalog copy consensuses onto phased haplotype assemblies and counts per-copy genomic loci,
//! disambiguated by deterministic SUN witnesses. Consumes only copies.fa + the assemblies; never wires
//! into the RNA-exclusive core. See docs/superpowers/specs/2026-07-14-assembly-parcn-design.md.

use std::collections::BTreeMap;

#[derive(Clone, Debug)]
pub struct Copy {
    pub family_id: String,
    pub copy_id: String,
    pub seq: Vec<u8>,
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
}

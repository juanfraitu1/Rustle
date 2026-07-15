//! Generalized projection (`--project-all-families`): project EVERY resolved copy's consensus (not one per
//! family) to localize members a single family-consensus projection misses. OPTIONAL, additive, DNA-localized
//! parCN leg — never changes the RNA-split catalog or the family definition. Pure extraction + dedup + row
//! formatting here; the minimap2 projection + read-support gate are wired in gw_family_catalog.

use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct CopyIn { pub seq: Vec<u8>, pub chrom: String, pub start: u64, pub end: u64 }

/// One `(family_id, consensus)` entry PER COPY, with the family_id repeated across its copies so
/// `project_families_batch` unions all copies' hits under that one family key.
pub fn all_copy_consensuses(fams: &[(String, Vec<CopyIn>)]) -> Vec<(String, Vec<u8>)> {
    fams.iter().flat_map(|(fid, copies)| copies.iter().map(move |c| (fid.clone(), c.seq.clone()))).collect()
}

/// Per-family copy spans, for the projection's `known` self-exclusion (a copy projecting back onto its own
/// catalogued locus is not a new localization).
pub fn known_from_fams(fams: &[(String, Vec<CopyIn>)]) -> HashMap<String, Vec<(String, u64, u64)>> {
    fams.iter().map(|(fid, copies)| (fid.clone(), copies.iter().map(|c| (c.chrom.clone(), c.start, c.end)).collect())).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extraction_repeats_fid_and_collects_spans() {
        let fams = vec![
            ("GWFAM0".to_string(), vec![
                CopyIn { seq: b"ACGT".to_vec(), chrom: "chr1".into(), start: 100, end: 200 },
                CopyIn { seq: b"ACGA".to_vec(), chrom: "chr1".into(), start: 500, end: 600 },
            ]),
            ("GWFAM1".to_string(), vec![
                CopyIn { seq: b"TTTT".to_vec(), chrom: "chr2".into(), start: 10, end: 20 },
            ]),
        ];
        let cons = all_copy_consensuses(&fams);
        // one entry per copy, family_id repeated per copy (so project_families_batch unions per family)
        assert_eq!(cons, vec![
            ("GWFAM0".to_string(), b"ACGT".to_vec()),
            ("GWFAM0".to_string(), b"ACGA".to_vec()),
            ("GWFAM1".to_string(), b"TTTT".to_vec()),
        ]);
        let known = known_from_fams(&fams);
        assert_eq!(known["GWFAM0"], vec![("chr1".to_string(),100,200), ("chr1".to_string(),500,600)]);
        assert_eq!(known["GWFAM1"], vec![("chr2".to_string(),10,20)]);
    }
}

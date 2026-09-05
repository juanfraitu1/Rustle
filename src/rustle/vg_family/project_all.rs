//! Generalized projection (`--project-all-families`): project EVERY resolved copy's consensus (not one per
//! family) to localize members a single family-consensus projection misses. OPTIONAL, additive, DNA-localized
//! parCN leg — never changes the RNA-split catalog or the family definition. Pure extraction + dedup + row
//! formatting here; the minimap2 projection + read-support gate are wired in gw_family_catalog.
//!
//! **STATUS:** OPT-IN — --project-all-families, `#[arg(long, default_value_t = false)]` at src/bin/gw_family_catalog.rs:195-196; equivalently env RUSTLE_PROJECT_ALL_FAMILIES=  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use std::collections::HashMap;

use crate::vg_family::genome_projection::CopyLocus;

fn recip_overlap(a: &CopyLocus, b: &CopyLocus) -> f64 {
    if a.chrom != b.chrom { return 0.0; }
    let (lo, hi) = (a.start.max(b.start), a.end.min(b.end));
    if hi <= lo { return 0.0; }
    let ov = (hi - lo) as f64;
    (ov / (a.end - a.start).max(1) as f64).min(ov / (b.end - b.start).max(1) as f64)
}

/// Collapse reciprocal-overlap ≥ 0.50 loci (from different sibling consensuses hitting one genomic locus)
/// into one, keeping the highest-identity survivor.
pub fn dedup_overlapping(mut loci: Vec<CopyLocus>) -> Vec<CopyLocus> {
    loci.sort_by(|a, b| b.identity.partial_cmp(&a.identity).unwrap_or(std::cmp::Ordering::Equal));
    let mut kept: Vec<CopyLocus> = Vec::new();
    for l in loci {
        if !kept.iter().any(|k| recip_overlap(k, &l) >= 0.50) { kept.push(l); }
    }
    kept
}

pub fn overlaps_any(chrom: &str, start: u64, end: u64, spans: &[(String, u64, u64)]) -> bool {
    spans.iter().any(|(c, s, e)| c == chrom && !(*s > end || *e < start))
}

pub fn format_allproj_row(family_id: &str, l: &CopyLocus, n_support: usize, overlaps_existing: bool) -> String {
    format!("{family_id}\t{}\t{}\t{}\t{:.3}\t{}\t{}", l.chrom, l.start, l.end, l.identity, n_support, overlaps_existing)
}

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
    fn dedup_overlap_and_row_format() {
        use crate::vg_family::genome_projection::CopyLocus;
        let mk = |s: u64, e: u64, id: f64| CopyLocus { chrom: "chr7".into(), start: s, end: e, identity: id, cov: 0.95 };
        // two overlapping hits (id .994 vs .982) + one disjoint -> 2 survivors, higher id kept
        let out = dedup_overlapping(vec![mk(75976253,75991692,0.994), mk(75976300,75991600,0.982), mk(76360590,76375995,0.990)]);
        assert_eq!(out.len(), 2);
        assert!(out.iter().any(|l| (l.identity - 0.994).abs() < 1e-9 && l.start == 75976253));
        assert!(overlaps_any("chr7", 75976300, 75991600, &[("chr7".into(), 75976253, 75991692)]));
        assert!(!overlaps_any("chr7", 76360590, 76375995, &[("chr7".into(), 75976253, 75991692)]));
        let row = format_allproj_row("GWFAM7", &mk(75976253,75991692,0.994), 41, false);
        assert_eq!(row, "GWFAM7\tchr7\t75976253\t75991692\t0.994\t41\tfalse");
    }

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

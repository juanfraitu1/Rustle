//! DEAD CODE -- StringTie-era assembler; slated for removal.
//! NOT part of the multi-copy-family thesis (O1 family-def / O2 copy-assign / O3 ASJ / O4 absent-copy).
//! See docs/RETIREMENT_AND_MIGRATION.md. Do not extend.
//! -G2 reports: vg_families.tsv (family merge audit) + vg_eval.tsv (head-to-head vs -G).
use crate::annotation_families::AnnotationFamily;

/// Render the family-merge audit. copy_chroms is read from each copy's own chrom
/// (NOT a family-graph node label, which collapses multi-chrom families to copies[0]).
pub fn render_vg_families(families: &[AnnotationFamily], threshold: f64) -> String {
    let mut out = String::from("family_id\tn_copies\tcopy_chroms\tmerge_threshold_T\tachieved_min_sim\tachieved_mean_sim\n");
    for f in families {
        let chroms: Vec<&str> = f.copies.iter().map(|c| c.chrom.as_str()).collect();
        out.push_str(&format!("{}\t{}\t{}\t{:.3}\t{:.3}\t{:.3}\n",
            f.family_id, f.copies.len(), chroms.join(";"),
            threshold, f.achieved_min_sim, f.achieved_mean_sim));
    }
    out
}

/// Per-copy verdict. in_vg = rustle-VG emitted a matching chain; in_st = -G (StringTie) had it.
pub fn verdict(in_vg: bool, in_st: bool) -> &'static str {
    match (in_vg, in_st) {
        (true, false) => "WIN",
        (true, true) => "tie",
        (false, false) => "miss",
        (false, true) => "regression",
    }
}

/// True if `chain` (ordered introns) matches any chain in `set` within +/-2 bp per coordinate.
pub fn chain_in_set(chain: &[(u64, u64)], set: &std::collections::HashSet<Vec<(u64, u64)>>) -> bool {
    for cand in set {
        if cand.len() != chain.len() { continue; }
        if chain.iter().zip(cand).all(|(a, b)| a.0.abs_diff(b.0) <= 2 && a.1.abs_diff(b.1) <= 2) {
            return true;
        }
    }
    false
}

/// Render the head-to-head scorecard. Each entry = (copy_id, family_id, in_st, in_vg).
pub fn render_vg_eval(entries: &[(String, String, bool, bool)]) -> String {
    let mut out = String::from("copy_id\tfamily_id\tin_annotation\tin_stringtie\tin_vg\tverdict\n");
    let (mut win, mut tie, mut miss, mut reg) = (0, 0, 0, 0);
    for (cid, fid, in_st, in_vg) in entries {
        let v = verdict(*in_vg, *in_st);
        match v {
            "WIN" => win += 1,
            "tie" => tie += 1,
            "miss" => miss += 1,
            _ => reg += 1,
        }
        let st = if *in_st { "yes" } else { "no" };
        let vg = if *in_vg { "yes" } else { "no" };
        // in_annotation is always "yes": every scorecard row is a -G2 annotation copy by construction.
        out.push_str(&format!("{}\t{}\tyes\t{}\t{}\t{}\n", cid, fid, st, vg, v));
    }
    out.push_str(&format!("SUMMARY\tmode=guided\tWIN={}\ttie={}\tmiss={}\tregression={}\n", win, tie, miss, reg));
    out
}

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

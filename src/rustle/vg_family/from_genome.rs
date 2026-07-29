//! DNA front-end: discover duplicated genomic loci by self-alignment and emit them as reps for the
//! shared homology-grouping core (`denovo_pipeline::homology_blocks`). Read-free and annotation-free —
//! the genome-only counterpart of the RNA read front-end. The reps differ from RNA reps in exactly one
//! way: genomic `seq` (introns included) and an EMPTY intron chain. That single difference is the
//! scientific claim (splicing discards the intron/flank sequence that separates near-identical copies).
use std::collections::{HashMap, HashSet};
use anyhow::Result;
use crate::vg_family::family_detect::DenovoTranscript;
use crate::genome::GenomeIndex;
use crate::vg_family::genome_projection::project_families_batch;

/// Parameters for the DNA-mode SD detector. `min_identity` is the LOCUS-DISCOVERY floor ("is this segment
/// duplicated enough to be a candidate locus"); it is distinct from the family-GROUPING identity floor
/// applied later inside `homology_blocks`. Env overrides let the SD floor be swept (e.g. 0.98 for SD98).
pub struct GenomeRepParams {
    pub min_identity: f64,
    pub min_block: u64,
    pub max_locus_span: u64,
    pub minimap2: String,
    pub threads: usize,
}
impl Default for GenomeRepParams {
    fn default() -> Self {
        Self { min_identity: 0.90, min_block: 1000, max_locus_span: 3_000_000,
                minimap2: "minimap2".into(), threads: 4 }
    }
}
impl GenomeRepParams {
    pub fn from_env() -> Self {
        let mut p = Self::default();
        if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_IDENTITY") { if let Ok(x) = v.parse() { p.min_identity = x; } }
        if let Ok(v) = std::env::var("RUSTLE_GENOME_MIN_BLOCK") { if let Ok(x) = v.parse() { p.min_block = x; } }
        if let Ok(v) = std::env::var("RUSTLE_MINIMAP2") { p.minimap2 = v; }
        p
    }
}

/// Discover duplicated genomic loci within `windows` and return them as reps (genomic sequence, empty
/// intron chain) for `homology_blocks`. Steps: (1) self-align each window's sequence against the genome
/// FASTA (`project_families_batch`, the same minimap2 primitive the famCN projection uses) to find every
/// locus it recurs at; (2) collect the hit loci across all windows, keep blocks in `[min_block,
/// max_locus_span]`, merge overlapping loci; (3) emit one rep per merged locus with its genomic sequence.
pub fn genome_reps(
    fasta_path: &str,
    windows: &[(String, u64, u64)],
    p: &GenomeRepParams,
) -> Result<Vec<DenovoTranscript>> {
    let contigs: HashSet<String> = windows.iter().map(|(c, _, _)| c.clone()).collect();
    let genome = GenomeIndex::from_fasta_contigs(fasta_path, &contigs)?;

    // (1) SD detector: each window's sequence is a query; project_families_batch returns every genome
    // locus it recurs at (identity >= min_identity). One batched minimap2 pass. `known` empty = keep all
    // hits (incl. the self locus — a window is itself a candidate locus; grouping decides families).
    let consensuses: Vec<(String, Vec<u8>)> = windows.iter().enumerate().filter_map(|(i, (c, s, e))| {
        genome.fetch_sequence(c, *s, *e).map(|seq| (format!("w{i}"), seq))
    }).collect();
    let known: HashMap<String, Vec<(String, u64, u64)>> = HashMap::new();
    let cov = 0.0_f64; // block length is gated by min_block below, not by fractional window coverage
    let hits = project_families_batch(&consensuses, fasta_path, &known, p.min_identity, cov,
                                      &p.minimap2, p.threads)?;

    // (2) rep construction. Each WINDOW is a locus to group — emit it directly so every member locus is a
    // node the quasi-clique can place (this is Soto's setup: given the loci, group them by shared sequence).
    // The self-alignment hits then ADD paralog loci OUTSIDE the windows, giving windowed singletons a sibling.
    // Windows are NOT merged (they are distinct member loci) and NOT min_block-filtered (a short member is
    // still a locus); only the extra discovered loci are size-gated and deduped, so dense paralog regions
    // can't fuse into giant blocks that miss a member's exact coordinate.
    let mut reps = Vec::new();
    let mut window_spans: Vec<(String, u64, u64)> = Vec::with_capacity(windows.len());
    for (c, s, e) in windows {
        if let Some(seq) = genome.fetch_sequence(c, *s, *e) {
            reps.push(DenovoTranscript {
                tid: format!("DN_{c}_{s}_1"),
                chrom: c.clone(), start: *s, end: *e, n_reads: 1, strand: '+',
                introns: vec![], seq, distinguishing_uniq: 0, tes: None,
            });
            window_spans.push((c.clone(), *s, *e));
        }
    }
    // discovered paralog loci that fall OUTSIDE every window (a paralog at a locus no member covers).
    let mut extra: Vec<(String, u64, u64)> = Vec::new();
    for hs in hits.into_values() {
        for h in hs {
            let len = h.end.saturating_sub(h.start);
            if len < p.min_block || len > p.max_locus_span { continue; }
            let inside_window = window_spans.iter().any(|(wc, ws, we)| *wc == h.chrom && h.start < *we && *ws < h.end);
            if inside_window { continue; }
            extra.push((h.chrom, h.start, h.end));
        }
    }
    extra.sort_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
    for (chrom, start, end) in merge_overlapping(&extra) {
        if let Some(seq) = genome.fetch_sequence(&chrom, start, end) {
            reps.push(DenovoTranscript {
                tid: format!("DN_{chrom}_{start}_1"),
                chrom, start, end, n_reads: 1, strand: '+',
                introns: vec![], seq, distinguishing_uniq: 0, tes: None,
            });
        }
    }
    Ok(reps)
}

/// Single-linkage merge of overlapping genomic intervals (input sorted by (chrom, start)).
fn merge_overlapping(loci: &[(String, u64, u64)]) -> Vec<(String, u64, u64)> {
    let mut out: Vec<(String, u64, u64)> = Vec::new();
    for (c, s, e) in loci.iter().cloned() {
        match out.last_mut() {
            Some((pc, _ps, pe)) if *pc == c && s <= *pe => { *pe = (*pe).max(e); }
            _ => out.push((c, s, e)),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_overlapping_joins_adjacent_and_keeps_disjoint() {
        let loci = vec![
            ("chr1".to_string(), 10, 100),
            ("chr1".to_string(), 90, 200),   // overlaps previous -> merge to 10..200
            ("chr1".to_string(), 500, 600),  // disjoint -> separate
            ("chr2".to_string(), 0, 50),     // different contig -> separate
        ];
        let m = merge_overlapping(&loci);
        assert_eq!(m, vec![
            ("chr1".to_string(), 10, 200),
            ("chr1".to_string(), 500, 600),
            ("chr2".to_string(), 0, 50),
        ]);
    }

    #[test]
    fn genome_reps_finds_family_copies_with_genomic_seq_and_no_introns() {
        // Real subset fixture: 3 near-identical NCF1 copies + 2 unrelated decoys, each as its own contig.
        // genome_reps must surface the duplicated NCF1 copies as reps (with genomic seq, empty introns).
        if std::process::Command::new("minimap2").arg("--version").output().is_err() { return; }
        let fa = "tests/fixtures/from_genome/subset.fa";
        if std::fs::metadata(fa).is_err() { eprintln!("fixture absent; skip"); return; }
        // windows = each contig full-length (as in windows.bed).
        let windows: Vec<(String, u64, u64)> = [
            ("NCF1", 15440u64), ("NCF1B", 15319), ("NCF1C", 15406), ("DECOY1", 10001), ("DECOY2", 10001),
        ].iter().map(|(c, l)| (c.to_string(), 0u64, *l)).collect();
        let p = GenomeRepParams { min_identity: 0.90, min_block: 400, ..Default::default() };
        let reps = genome_reps(fa, &windows, &p).unwrap();
        // the three NCF1 copies must all be discovered as duplicated loci.
        for want in ["NCF1", "NCF1B", "NCF1C"] {
            assert!(reps.iter().any(|r| r.chrom == want), "missing duplicated locus {want}; got {:?}",
                    reps.iter().map(|r| r.chrom.as_str()).collect::<Vec<_>>());
        }
        // every rep is genomic: empty intron chain and seq length == span.
        for r in &reps {
            assert!(r.introns.is_empty(), "DNA rep {} must have empty intron chain", r.chrom);
            assert_eq!(r.seq.len() as u64, r.end - r.start, "DNA rep {} seq must be the full genomic span", r.chrom);
        }
    }
}

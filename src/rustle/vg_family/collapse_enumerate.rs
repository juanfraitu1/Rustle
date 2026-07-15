//! K=0-collapsed family re-admission (behind `--collapse-enumerate`). A near-identical family that
//! collapses to <2 RNA-distinct loci is re-admitted as copy NUMBER iff it shows a LOCAL collapse:
//! a `hidden_copy` second-haplotype witness that is BALANCED (co-equal depth) AND projects to >=2 genomic loci.
use crate::vg_family::hidden_copy::HiddenCopyEvidence;
use crate::vg_family::genome_projection::CopyLocus;

use crate::vg_family::denovo_assemble::{BamRead, reads_in_region};
use crate::vg_family::hidden_copy::{ReadObs, HiddenCopyParams, detect_hidden_copy};
use crate::vg_family::genome_projection::project_family_copies;
use crate::genome::GenomeIndex;

/// Min balanced-alt fraction: a co-equal collapsed 2nd copy (~0.5), not a minor het/edit (~<=0.1).
pub const MIN_ALT_FRAC: f64 = 0.30;

/// famCN from a genome projection's loci count. `project_family_copies` is called with the candidate's
/// own span as the sole `known` entry, so its returned loci are the OTHER genomic copies -- the seed
/// locus itself is excluded by construction. Total genomic copy number = the seed copy + the projected
/// others, hence the `+ 1`.
pub fn famcn_from_projection(n_projection_loci: usize) -> usize {
    n_projection_loci + 1
}

/// The three-signal gate. ALL must hold (see spec / Global Constraints).
pub fn admit_collapse(ev: &HiddenCopyEvidence, n_projection_loci: usize) -> bool {
    ev.flagged && ev.alt_read_fraction >= MIN_ALT_FRAC && n_projection_loci >= 2
}

#[derive(Debug, Clone)]
pub struct CollapsedFamily {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub famcn: usize,          // total genomic copy number = seed locus + projected other loci
    pub n_alt_reads: usize,    // hidden 2nd-haplotype depth
    pub alt_read_fraction: f64,
    pub projection: Vec<CopyLocus>,
}

/// Per-read mismatch positions vs the reference window `[lo,hi)` — the `alts` a `detect_hidden_copy`
/// column analysis consumes. PRIMARY alignments only (project invariant, `-F 2308`): secondary and
/// supplementary records are additional placements of the SAME physical molecule, and `detect_hidden_copy`
/// counts every `ReadObs` as an independent molecule — including them here would double-count exactly the
/// multimapping reads this feature targets, corrupting `alt_read_fraction`. Only `M/=/X` CIGAR ops advance
/// both ref and query; `I/S` advance query, `D/N` advance ref. A position counts as an alt iff the read
/// base differs from the reference base. `refwin` can be shorter than `hi - lo` when `[lo,hi)` runs past
/// the contig end (subtelomeric candidates); `refwin.get(..)` skips positions past the fetched window
/// instead of panicking.
fn read_obs_from_bam_reads(reads: &[BamRead], chrom: &str, lo: u64, hi: u64, genome: &GenomeIndex) -> Vec<ReadObs> {
    let refwin = match genome.fetch_sequence(chrom, lo, hi) { Some(s) => s, None => return Vec::new() };
    let mut out = Vec::with_capacity(reads.len());
    for br in reads {
        if br.is_secondary || br.is_supplementary {
            continue;
        }
        let r = &br.read;
        let mut ref_pos = r.ref_start;
        let mut q = 0usize;
        let mut alts = Vec::new();
        for &(op, len) in &r.cigar {
            match op {
                'M' | '=' | 'X' => {
                    for k in 0..len {
                        let rp = ref_pos + k;
                        if rp >= lo && rp < hi {
                            if let Some(&rb) = refwin.get((rp - lo) as usize) {
                                if let Some(&qb) = r.seq.get(q + k as usize) {
                                    if qb.to_ascii_uppercase() != rb.to_ascii_uppercase() {
                                        alts.push(rp);
                                    }
                                }
                            }
                        }
                    }
                    ref_pos += len; q += len as usize;
                }
                'I' | 'S' => { q += len as usize; }
                'D' | 'N' => { ref_pos += len; }
                _ => {}
            }
        }
        out.push(ReadObs { start: r.ref_start, end: ref_pos, alts });
    }
    out
}

/// Re-admission driver for one dropped collapsed candidate locus. Three-signal gate:
/// local hidden-copy witness (balanced 2nd haplotype) + >=2 genome-projected loci.
/// Returns `Some(CollapsedFamily)` on admit, `None` otherwise (including any I/O failure — a
/// dropped candidate that cannot be evaluated stays dropped, exactly as today).
pub fn readmit_locus(
    bam_path: &str, chrom: &str, lo: u64, hi: u64, consensus: &[u8],
    genome: &GenomeIndex, fasta_path: &str, minimap2: &str, threads: usize,
) -> Option<CollapsedFamily> {
    let (_p, bam_reads) = reads_in_region(bam_path, chrom, lo, hi, threads).ok()?;
    let obs = read_obs_from_bam_reads(&bam_reads, chrom, lo, hi, genome);
    let ev = detect_hidden_copy(&obs, &HiddenCopyParams::default());
    if !ev.flagged { return None; }                       // short-circuit before the expensive projection
    let known = vec![(chrom.to_string(), lo, hi)];
    let loci = project_family_copies(consensus, fasta_path, &known, 0.98, 0.90, minimap2, threads).ok()?;
    if !admit_collapse(&ev, loci.len()) { return None; }
    Some(CollapsedFamily {
        chrom: chrom.to_string(), start: lo, end: hi,
        famcn: famcn_from_projection(loci.len()), n_alt_reads: ev.n_alt_reads, alt_read_fraction: ev.alt_read_fraction,
        projection: loci,
    })
}

/// One `<out>.collapsed.tsv` data row for a re-admitted K=0-collapsed family. Columns:
/// family_id, chrom, start, end, famCN, n_alt_reads, alt_frac(3dp), status, projection_loci
/// (`chrom:start-end@identity` joined by `;`). Copy-NUMBER only — these families never appear in copies.tsv.
pub fn format_collapsed_row(family_id: &str, f: &CollapsedFamily) -> String {
    let proj = f.projection.iter()
        .map(|c| format!("{}:{}-{}@{:.3}", c.chrom, c.start, c.end, c.identity))
        .collect::<Vec<_>>().join(";");
    format!("{family_id}\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{}\t{}",
        f.chrom, f.start, f.end, f.famcn, f.n_alt_reads, f.alt_read_fraction, "K0_COLLAPSED", proj)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn ev(flagged: bool, frac: f64) -> HiddenCopyEvidence {
        HiddenCopyEvidence { n_primary_reads: 300, n_alt_positions: 40, n_alt_reads: (300.0*frac) as usize, alt_read_fraction: frac, flagged }
    }
    #[test]
    fn admits_only_when_all_three_signals_hold() {
        assert!(admit_collapse(&ev(true, 0.50), 2), "flagged + balanced + >=2 loci -> admit");
        assert!(!admit_collapse(&ev(false, 0.50), 2), "not flagged -> reject");
        assert!(!admit_collapse(&ev(true, 0.10), 2), "minor 2nd haplotype (het/edit-like) -> reject");
        assert!(!admit_collapse(&ev(true, 0.50), 1), "single projection locus -> reject");
    }

    /// The gate is `n_projection_loci >= 2`: exactly `MIN_ALT_FRAC` and exactly 2 projection loci
    /// must ADMIT (not a strict `>` gate on either signal).
    #[test]
    fn admit_collapse_boundary_at_min_alt_frac_and_two_loci() {
        assert!(admit_collapse(&ev(true, MIN_ALT_FRAC), 2), "alt_read_fraction == MIN_ALT_FRAC exactly -> admit (>=)");
    }

    /// `famcn` is seed-inclusive: `project_family_copies` excludes the seed's own locus (it is the sole
    /// `known` entry), so the projection count is copies OTHER than the seed. Total famCN = seed + others.
    #[test]
    fn famcn_is_seed_inclusive() {
        assert_eq!(famcn_from_projection(0), 1, "no other projected loci -> just the seed copy");
        assert_eq!(famcn_from_projection(3), 4, "3 other projected loci + the seed copy");
    }

    #[test]
    fn readmit_decision_from_readobs_balanced_vs_het() {
        use crate::vg_family::hidden_copy::{ReadObs, HiddenCopyParams, detect_hidden_copy};
        // 20 candidate columns; ~half the reads carry every alt (a co-equal collapsed 2nd copy)
        let cols: Vec<u64> = (0..20).map(|i| 1000 + i * 10).collect();
        let mk = |carry: bool| ReadObs { start: 1000, end: 1200, alts: if carry { cols.clone() } else { vec![] } };
        let mut collapse: Vec<ReadObs> = (0..150).map(|_| mk(true)).collect();
        collapse.extend((0..150).map(|_| mk(false)));           // 0.50 balanced
        let ev = detect_hidden_copy(&collapse, &HiddenCopyParams::default());
        assert!(ev.flagged && ev.alt_read_fraction >= MIN_ALT_FRAC);
        assert!(admit_collapse(&ev, 2), "balanced collapse + 2 loci admits");
        // minor het: only 8% carry the alts
        let mut het: Vec<ReadObs> = (0..24).map(|_| mk(true)).collect();
        het.extend((0..276).map(|_| mk(false)));                // 0.08
        let ev2 = detect_hidden_copy(&het, &HiddenCopyParams::default());
        assert!(!admit_collapse(&ev2, 2), "minor het does not admit");
    }

    use crate::vg_family::copy_split::AlignedRead;

    /// Build a `BamRead` with a single mismatch (alt) baked into an otherwise-reference-matching
    /// `M`-only alignment at `mismatch_offset` (relative to `ref_start`), so `read_obs_from_bam_reads`
    /// always produces exactly one alt -- at `ref_start + mismatch_offset` -- per surviving read. A
    /// distinct offset per read gives each read an identifiable alt position, so tests can confirm
    /// WHICH reads survived a filter, not just how many.
    fn mk_bam_read(ref_start: u64, len: u64, mismatch_offset: u64, is_secondary: bool, is_supplementary: bool, name: &str) -> BamRead {
        let mut seq = vec![b'A'; len as usize];
        seq[mismatch_offset as usize] = b'C'; // mismatch vs an all-'A' reference at relative offset `mismatch_offset`
        BamRead {
            chrom: "c1".to_string(),
            read: AlignedRead { ref_start, cigar: vec![('M', len)], seq, qual: Vec::new() },
            mapq: 60,
            name: name.to_string(),
            as_score: 0,
            de: 0.0,
            is_supplementary,
            is_secondary,
        }
    }

    #[test]
    fn read_obs_skips_secondary_and_supplementary() {
        let seq = vec![b'A'; 200];
        let genome = GenomeIndex::from_seqs(&[("c1", &seq[..])]);
        // Each read carries a distinct, identifiable alt offset so the surviving ReadObs can be
        // matched back to the exact primary reads that produced them (not merely counted).
        let reads = vec![
            mk_bam_read(10, 50, 5, false, false, "primary1"),      // alt at 15
            mk_bam_read(10, 50, 10, true, false, "secondary"),     // alt at 20 (must NOT survive)
            mk_bam_read(10, 50, 15, false, true, "supplementary"), // alt at 25 (must NOT survive)
            mk_bam_read(10, 50, 20, false, false, "primary2"),     // alt at 30
        ];
        let obs = read_obs_from_bam_reads(&reads, "c1", 0, 200, &genome);
        assert_eq!(obs.len(), 2, "only the two primary (non-secondary, non-supplementary) reads survive");
        let alts: Vec<Vec<u64>> = obs.iter().map(|o| o.alts.clone()).collect();
        assert_eq!(
            alts,
            vec![vec![15], vec![30]],
            "survivors are exactly primary1 (alt@15) and primary2 (alt@30), in read order -- not the \
             secondary (alt@20) or supplementary (alt@25)"
        );
    }

    #[test]
    fn read_obs_window_past_contig_end_does_not_panic() {
        let seq = vec![b'A'; 100];
        let genome = GenomeIndex::from_seqs(&[("c1", &seq[..])]);
        // Read near the contig end; window requested [0, 100_000) runs far past the 100bp contig,
        // so the fetched reference window is truncated to 100bp while reads may extend past it.
        let reads = vec![mk_bam_read(90, 20, 5, false, false, "primary_near_end")];
        let obs = read_obs_from_bam_reads(&reads, "c1", 0, 100_000, &genome);
        assert_eq!(obs.len(), 1, "the read is kept; only its out-of-window positions are dropped");
        assert_eq!(
            obs[0].alts,
            vec![95],
            "the read's in-window alt (ref_start 90 + offset 5 = 95, within the truncated 100bp \
             reference window) is present -- only positions past the fetched window are dropped"
        );
    }

    #[test]
    fn collapsed_tsv_row_format() {
        use crate::vg_family::genome_projection::CopyLocus;
        let f = CollapsedFamily {
            chrom: "chr2".into(), start: 108994973, end: 109147842, famcn: 2, n_alt_reads: 600, alt_read_fraction: 0.49,
            projection: vec![
                CopyLocus { chrom: "chr2".into(), start: 108994973, end: 109147842, identity: 0.99,  cov: 0.95 },
                CopyLocus { chrom: "chr2".into(), start: 110869109, end: 110895544, identity: 0.993, cov: 0.92 },
            ],
        };
        let row = format_collapsed_row("GWFAMc0", &f);
        assert_eq!(row, "GWFAMc0\tchr2\t108994973\t109147842\t2\t600\t0.490\tK0_COLLAPSED\tchr2:108994973-109147842@0.990;chr2:110869109-110895544@0.993");
    }
}

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

/// The three-signal gate. ALL must hold (see spec / Global Constraints).
pub fn admit_collapse(ev: &HiddenCopyEvidence, n_projection_loci: usize) -> bool {
    ev.flagged && ev.alt_read_fraction >= MIN_ALT_FRAC && n_projection_loci >= 2
}

#[derive(Debug, Clone)]
pub struct CollapsedFamily {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub famcn: usize,          // genome-projected copy number
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
        famcn: loci.len(), n_alt_reads: ev.n_alt_reads, alt_read_fraction: ev.alt_read_fraction,
        projection: loci,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vg_family::hidden_copy::HiddenCopyEvidence;
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
}

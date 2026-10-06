//! K=0-collapsed family re-admission (behind `--collapse-enumerate`). A near-identical family that
//! collapses to <2 RNA-distinct loci is re-admitted as copy NUMBER iff it shows a LOCAL collapse:
//! a `hidden_copy` second-haplotype witness that is BALANCED (co-equal depth) AND projects to >=2 genomic loci.
//!
//! **STATUS:** OPT-IN — --collapse-enumerate (src/bin/gw_family_catalog.rs:177-178, default_value_t = false) or env RUSTLE_COLLAPSE_ENUMERATE=1 (denovo_pipeline.rs:180); sibl  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)
use std::collections::HashMap;
use crate::vg_family::collapse_enumerate::hidden_copy::HiddenCopyEvidence;
use crate::vg_family::genome_projection::CopyLocus;

use crate::vg_family::denovo_assemble::{BamRead, reads_in_region};
use crate::vg_family::collapse_enumerate::hidden_copy::{ReadObs, HiddenCopyParams, detect_hidden_copy};
use crate::vg_family::genome_projection::{project_family_copies, project_families_batch};
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

/// Minimum PRIMARY reads at a projected locus for it to count as an EXPRESSED copy (the EEF1A1 guard:
/// silent pseudogenes fall below this).
pub const MIN_LOCUS_READS: usize = 3;

/// Expressed-collapsed admission (PSV-free): a dropped candidate is a real multi-copy family iff it
/// projects to >= 2 genomic loci that are EACH read-supported. No hidden-copy witness required — these
/// families are exon-identical (0 PSVs) but their copies are all transcribed.
pub fn admit_expressed_collapse(read_supported_loci: usize) -> bool {
    read_supported_loci >= 2
}

#[derive(Debug, Clone)]
pub struct ExpressedCollapsedFamily {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub famcn: usize,           // seed-inclusive: read-supported projected loci + the seed
    pub min_locus_reads: usize, // weakest admitted locus's support (transparency)
    pub projection: Vec<CopyLocus>,
}

/// One `<out>.expressed_collapsed.tsv` row: family_id, chrom, start, end, famCN, min_locus_reads,
/// status (`K0_COLLAPSED_EXPRESSED`), projection_loci (`chrom:start-end@identity` joined by `;`).
pub fn format_expressed_collapsed_row(family_id: &str, f: &ExpressedCollapsedFamily) -> String {
    let proj = f.projection.iter()
        .map(|c| format!("{}:{}-{}@{:.3}", c.chrom, c.start, c.end, c.identity))
        .collect::<Vec<_>>().join(";");
    format!("{family_id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        f.chrom, f.start, f.end, f.famcn, f.min_locus_reads, "K0_COLLAPSED_EXPRESSED", proj)
}

/// One `<out>.dna_family.tsv` row for an RNA-orphan locus recovered by the DNA edge oracle (`--dna-family-
/// fallback`): same schema as `expressed_collapsed`, status `DNA_FAMILY_RNA_NONHOMOLOGOUS`. The locus's
/// EXPRESSED transcript is non-homologous to its paralogs (no RNA family forms), yet it projects to >= 2
/// DIVERGENT genomic copies. Copy NUMBER only; per-read resolution needs DNA parCN.
pub fn format_dna_family_row(family_id: &str, f: &ExpressedCollapsedFamily) -> String {
    let proj = f.projection.iter()
        .map(|c| format!("{}:{}-{}@{:.3}", c.chrom, c.start, c.end, c.identity))
        .collect::<Vec<_>>().join(";");
    format!("{family_id}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        f.chrom, f.start, f.end, f.famcn, f.min_locus_reads, "DNA_FAMILY_RNA_NONHOMOLOGOUS", proj)
}

/// Pure: assemble an `ExpressedCollapsedFamily` from already-read-supported projection loci (paired with
/// their support counts). Returns `None` unless `>= 2` loci are supported. Factored out so the admit/famCN
/// logic is unit-testable without minimap2/BAM I/O.
fn build_expressed_family(chrom: &str, lo: u64, hi: u64, loci: Vec<CopyLocus>, supports: &[usize]) -> Option<ExpressedCollapsedFamily> {
    if !admit_expressed_collapse(loci.len()) { return None; }
    let min_locus_reads = supports.iter().copied().min().unwrap_or(0);
    Some(ExpressedCollapsedFamily {
        chrom: chrom.to_string(), start: lo, end: hi,
        famcn: famcn_from_projection(loci.len()), min_locus_reads, projection: loci,
    })
}

/// Batched re-admission for ALL dropped expressed candidates in ONE minimap2 index load (vs re-indexing the
/// genome per candidate). `candidates` = `(unique_id, chrom, lo, hi, consensus)`. Projects every consensus
/// at `min_identity`/cov>=0.90 (each candidate's own span excluded via `known`), then per-candidate keeps
/// projection loci with >= MIN_LOCUS_READS primary reads (the EEF1A1 guard) and admits if >= 2 remain.
/// `min_identity` = 0.98 for the exon-identical collapse-expressed path; the DNA-family fallback passes a
/// looser floor (~0.90) to reach DIVERGENT genomic paralogs whose expressed transcript is non-homologous.
pub fn readmit_expressed_batch(
    candidates: &[(String, String, u64, u64, Vec<u8>)],
    bam_path: &str, fasta_path: &str, minimap2: &str, threads: usize, min_identity: f64,
) -> Vec<ExpressedCollapsedFamily> {
    if candidates.is_empty() { return Vec::new(); }
    let consensuses: Vec<(String, Vec<u8>)> = candidates.iter().map(|(id, _, _, _, seq)| (id.clone(), seq.clone())).collect();
    let known: HashMap<String, Vec<(String, u64, u64)>> =
        candidates.iter().map(|(id, ch, lo, hi, _)| (id.clone(), vec![(ch.clone(), *lo, *hi)])).collect();
    let proj = project_families_batch(&consensuses, fasta_path, &known, min_identity, 0.90, minimap2, threads).unwrap_or_default();
    let mut out = Vec::new();
    for (id, chrom, lo, hi, _seq) in candidates {
        let loci = match proj.get(id) { Some(l) => l.clone(), None => continue };
        let mut supported = Vec::new();
        let mut supports = Vec::new();
        for l in loci {
            let n = reads_in_region(bam_path, &l.chrom, l.start, l.end, threads).map(|(p, _)| p.len()).unwrap_or(0);
            if n >= MIN_LOCUS_READS { supports.push(n); supported.push(l); }
        }
        if let Some(f) = build_expressed_family(chrom, *lo, *hi, supported, &supports) { out.push(f); }
    }
    out
}

/// Soft-mask (RepeatMasker lowercase) fraction of a reference sequence: lowercase bases over total. The repeat
/// signal for the dna-family gate — Alu SINEs and low-complexity, which project to ~99%-identical dispersed
/// copies genome-wide and masquerade as SD paralogs, are soft-masked in the reference. An empty slice is 0.0.
pub fn softmask_frac(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let lc = seq.iter().filter(|b| b.is_ascii_lowercase()).count();
    lc as f64 / seq.len() as f64
}

/// DNA-family re-admission for RNA-orphan loci (the DNA edge oracle, `--dna-family-fallback`). Projects each
/// orphan consensus at `min_identity`/cov>=0.90; the seed's own span is excluded via `known`, so the returned
/// loci are the OTHER genomic paralog copies. Admits as a DNA-family iff >= 1 other genomic copy exists (total
/// famCN >= 2 = the >= 2-genomic-copy definition). UNLIKE `readmit_expressed_batch`, the projected loci need
/// NOT be read-supported: a DNA-family's paralogs are frequently SILENT — that is precisely WHY the orphan's
/// expressed transcript is non-homologous to them (DNA-family != RNA-family). Per-locus read counts are still
/// computed and reported (min over paralogs; typically 0 for the silent copies) but do NOT gate admission.
/// Copy-NUMBER only; per-read resolution of the silent copies requires DNA parCN.
///
/// REPEAT GATE (`max_softmask`): a candidate whose re-admitted locus is >= `max_softmask` soft-masked
/// (repeat/transposon-derived) is REJECTED — an adversarial audit showed the un-gated fallback projects orphan
/// loci onto Alus/low-complexity (11/12 off-benchmark, both controls) rather than real gene families. This
/// mirrors Soto's own RepeatMasker exclusion. Soft-mask is read verbatim via `IndexedFasta` (the loaded
/// `GenomeIndex` upper-cases and cannot feed it).
pub fn readmit_dna_family_batch(
    candidates: &[(String, String, u64, u64, Vec<u8>)],
    bam_path: &str, fasta_path: &str, minimap2: &str, threads: usize, min_identity: f64, max_softmask: f64,
) -> Vec<ExpressedCollapsedFamily> {
    if candidates.is_empty() { return Vec::new(); }
    let consensuses: Vec<(String, Vec<u8>)> = candidates.iter().map(|(id, _, _, _, seq)| (id.clone(), seq.clone())).collect();
    let known: HashMap<String, Vec<(String, u64, u64)>> =
        candidates.iter().map(|(id, ch, lo, hi, _)| (id.clone(), vec![(ch.clone(), *lo, *hi)])).collect();
    let proj = project_families_batch(&consensuses, fasta_path, &known, min_identity, 0.90, minimap2, threads).unwrap_or_default();
    let fa = crate::genome::IndexedFasta::open(fasta_path).ok();
    let mut out = Vec::new();
    for (id, chrom, lo, hi, _seq) in candidates {
        // repeat gate: reject re-admitted loci that are repeat-dominated (soft-mask >= max_softmask).
        if let Some(fa) = &fa {
            if let Some(bytes) = fa.fetch(chrom, *lo as i64, *hi as i64) {
                if softmask_frac(&bytes) >= max_softmask {
                    continue;
                }
            }
        }
        let loci = match proj.get(id) { Some(l) => l.clone(), None => continue };
        if loci.is_empty() { continue; } // >= 1 other genomic copy => famCN >= 2 => DNA-family
        let min_reads = loci.iter()
            .map(|l| reads_in_region(bam_path, &l.chrom, l.start, l.end, threads).map(|(p, _)| p.len()).unwrap_or(0))
            .min().unwrap_or(0);
        out.push(ExpressedCollapsedFamily {
            chrom: chrom.clone(), start: *lo, end: *hi,
            famcn: famcn_from_projection(loci.len()), min_locus_reads: min_reads, projection: loci,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    fn ev(flagged: bool, frac: f64) -> HiddenCopyEvidence {
        HiddenCopyEvidence { n_primary_reads: 300, n_alt_positions: 40, n_alt_reads: (300.0*frac) as usize, alt_read_fraction: frac, flagged }
    }
    #[test]
    fn softmask_frac_counts_lowercase_over_total() {
        assert_eq!(softmask_frac(b"ACGT"), 0.0);
        assert_eq!(softmask_frac(b"acgt"), 1.0);
        assert_eq!(softmask_frac(b"ACac"), 0.5);
        assert_eq!(softmask_frac(b""), 0.0);
        assert_eq!(softmask_frac(b"ACGTn"), 0.2, "lowercase n (masked) counts as soft-masked");
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
        use crate::vg_family::collapse_enumerate::hidden_copy::{ReadObs, HiddenCopyParams, detect_hidden_copy};
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
            is_secondary, reverse: false, ts: None }
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

    #[test]
    fn admit_expressed_needs_two_read_supported_loci() {
        assert!(!admit_expressed_collapse(0));
        assert!(!admit_expressed_collapse(1));
        assert!(admit_expressed_collapse(2));
        assert!(admit_expressed_collapse(5));
    }

    #[test]
    fn expressed_collapsed_row_format() {
        let f = ExpressedCollapsedFamily {
            chrom: "chr2".into(), start: 97950885, end: 98048181, famcn: 3, min_locus_reads: 32,
            projection: vec![
                CopyLocus { chrom: "chr2".into(), start: 97950885, end: 98048181, identity: 0.998, cov: 0.95 },
                CopyLocus { chrom: "chr2".into(), start: 99100000, end: 99198000, identity: 0.994, cov: 0.93 },
            ],
        };
        assert_eq!(format_expressed_collapsed_row("GWFAMe0", &f),
            "GWFAMe0\tchr2\t97950885\t98048181\t3\t32\tK0_COLLAPSED_EXPRESSED\tchr2:97950885-98048181@0.998;chr2:99100000-99198000@0.994");
    }

    #[test]
    fn expressed_driver_builds_family_from_supported_loci() {
        // exercises the post-projection assembly logic: given read-supported loci, famCN is seed-inclusive and
        // min_locus_reads is the weakest. (The minimap2/BAM I/O path is covered by the Soto/EEF1A1 live runs.)
        let loci = vec![
            CopyLocus { chrom: "c".into(), start: 0, end: 9, identity: 0.99, cov: 0.95 },
            CopyLocus { chrom: "c".into(), start: 50, end: 59, identity: 0.995, cov: 0.95 },
        ];
        let supports = vec![32usize, 4usize];
        let fam = build_expressed_family("c", 0, 9, loci.clone(), &supports);
        assert!(fam.is_some());
        let fam = fam.unwrap();
        assert_eq!(fam.famcn, 3);            // 2 supported loci + seed
        assert_eq!(fam.min_locus_reads, 4);  // weakest
        // fewer than 2 supported -> None
        assert!(build_expressed_family("c", 0, 9, loci[..1].to_vec(), &supports[..1]).is_none());
    }
}

// ---- merged 2026-10-05: was `vg_family/hidden_copy.rs`, now the inline module below (one component) ----
#[allow(clippy::all)]
pub mod hidden_copy {
//! Detect gene-family copies PRESENT in the reads but ABSENT from the reference genome
//! (collapsed segdup / assembly gap / CNV / private duplication).
//!
//! A hidden copy's reads have no correct home, so they mismap to the closest sibling reference
//! copy carrying their PRIVATE SNPs — a COHERENT second haplotype the reference (one copy at
//! that locus) cannot explain. This detector finds that second haplotype among the locus's
//! PRIMARY alignments and FLAGS the discrepancy. Per the DAZ3 discipline it DETECTS and reports
//! evidence ("the reads imply ≥2 copies; the reference models 1") and ABSTAINS from placing or
//! fabricating the missing copy — it never manufactures a copy.
//!
//! Design = synthesis of an independent design panel (statistical / algorithmic / honesty
//! lenses), which converged on: PRIMARY-alignments-only matrix (the paralog-bleed firewall, since an
//! in-reference paralog's reads are primary at THEIR locus), candidate columns where the
//! non-reference allele frequency sits in a balanced band (excludes 0.5% sequencing error and
//! fixed differences), and a co-segregation/block test (a hidden copy's alt columns co-occur on
//! ONE read subset — distinguishing it from scattered heterozygous SNPs by requiring many).
//!
//! **STATUS:** OPT-IN — `--collapse-enumerate` (src/bin/gw_family_catalog.rs:177-178, `#[arg(long, default_value_t = false)]`)  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

/// One primary alignment's observation at the locus: its covered span [start, end) and the positions
/// where it carries a non-reference allele (a mismatch).
#[derive(Debug, Clone)]
pub struct ReadObs {
    pub start: u64,
    pub end: u64,
    pub alts: Vec<u64>,
}

#[derive(Debug, Clone)]
pub struct HiddenCopyParams {
    pub balanced_lo: f64,        // min alt-allele fraction for a candidate column (≫ error rate)
    pub balanced_hi: f64,        // max alt-allele fraction (above = fixed diff / ref error)
    pub min_depth: usize,        // min coverage at a candidate column
    pub min_alt_positions: usize,// min candidate columns to call a hidden copy (≫ a few hets)
    pub min_alt_reads: usize,    // min reads in the alt haplotype (the hidden copy's depth)
    pub share_hi: f64,           // a read joins H if alt at ≥ this fraction of candidate cols it covers
}

impl Default for HiddenCopyParams {
    fn default() -> Self {
        HiddenCopyParams {
            balanced_lo: 0.20,
            balanced_hi: 0.60,
            min_depth: 8,
            min_alt_positions: 12, // het firewall: a diploid het is 1-2 columns; a copy is dozens
            min_alt_reads: 5,
            share_hi: 0.60,
        }
    }
}

impl HiddenCopyParams {
    pub fn from_env() -> Self {
        let mut p = HiddenCopyParams::default();
        let getf = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<f64>().ok());
        let getu = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<usize>().ok());
        if let Some(v) = getf("RUSTLE_VG_HIDDEN_ALT_LO") { p.balanced_lo = v; }
        if let Some(v) = getf("RUSTLE_VG_HIDDEN_ALT_HI") { p.balanced_hi = v; }
        if let Some(v) = getu("RUSTLE_VG_HIDDEN_MIN_DEPTH") { p.min_depth = v; }
        if let Some(v) = getu("RUSTLE_VG_HIDDEN_MIN_POSITIONS") { p.min_alt_positions = v; }
        if let Some(v) = getu("RUSTLE_VG_HIDDEN_MIN_READS") { p.min_alt_reads = v; }
        p
    }
}

/// Evidence for a copy not in the reference. DETECT + FLAG only — no placement, no sequence.
#[derive(Debug, Clone, PartialEq)]
pub struct HiddenCopyEvidence {
    pub n_primary_reads: usize,
    pub n_alt_positions: usize, // coherent second-haplotype columns
    pub n_alt_reads: usize,     // reads in the alt haplotype (the hidden copy's apparent depth)
    pub alt_read_fraction: f64, // n_alt_reads / n_primary_reads
    pub flagged: bool,          // evidence of an unmodeled copy at this locus
}

/// Pure detector over PRIMARY alignments at one reference-copy locus. The caller MUST pass primary
/// reads only (the paralog-bleed firewall). Deterministic; no I/O.
pub fn detect_hidden_copy(reads: &[ReadObs], p: &HiddenCopyParams) -> HiddenCopyEvidence {
    let n = reads.len();
    let none = HiddenCopyEvidence {
        n_primary_reads: n, n_alt_positions: 0, n_alt_reads: 0,
        alt_read_fraction: 0.0, flagged: false,
    };
    if n < p.min_depth {
        return none;
    }

    // Per-position alt-read count (only positions some read calls alt are candidates).
    let mut alt_count: crate::types::DetHashMap<u64, usize> = Default::default();
    for r in reads {
        for &pos in &r.alts {
            *alt_count.entry(pos).or_insert(0) += 1;
        }
    }

    // Candidate columns: balanced alt fraction over sufficient depth. Error (~0.5%) never
    // reaches balanced_lo; a fixed difference / reference error sits above balanced_hi.
    let mut candidates: Vec<u64> = Vec::new();
    for (&pos, &ac) in &alt_count {
        if ac < 2 {
            continue; // a singleton is error, not a haplotype
        }
        let cov = reads.iter().filter(|r| r.start <= pos && pos < r.end).count();
        if cov < p.min_depth {
            continue;
        }
        let frac = ac as f64 / cov as f64;
        if frac >= p.balanced_lo && frac <= p.balanced_hi {
            candidates.push(pos);
        }
    }
    candidates.sort_unstable();
    let n_alt_positions = candidates.len();

    // Co-segregation: a hidden COPY's candidate columns co-occur on ONE read subset (H). Partition
    // reads by their alt-share over the candidate columns they cover; H = the alt haplotype.
    let cand_set: crate::types::DetHashSet<u64> = candidates.iter().copied().collect();
    let n_alt_reads = reads.iter().filter(|r| {
        let covered = candidates.iter().filter(|&&pos| r.start <= pos && pos < r.end).count();
        if covered == 0 {
            return false;
        }
        let alt_at = r.alts.iter().filter(|pos| cand_set.contains(pos)).count();
        (alt_at as f64 / covered as f64) >= p.share_hi
    }).count();

    // Flag only with MANY co-segregating positions (≫ a het) AND a real alt-haplotype read group.
    let flagged = n_alt_positions >= p.min_alt_positions && n_alt_reads >= p.min_alt_reads;

    HiddenCopyEvidence {
        n_primary_reads: n,
        n_alt_positions,
        n_alt_reads,
        alt_read_fraction: if n > 0 { n_alt_reads as f64 / n as f64 } else { 0.0 },
        flagged,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p() -> HiddenCopyParams { HiddenCopyParams::default() }

    // n reads all spanning [0, span); `hap` reads carry alt at every position in `shared`,
    // the rest carry alt at `noise` random-but-distinct positions each (sequencing error).
    fn reads(n: usize, span: u64, hap: usize, shared: &[u64], noise_per_read: u64) -> Vec<ReadObs> {
        (0..n).map(|r| {
            let mut alts: Vec<u64> = if r < hap { shared.to_vec() } else { Vec::new() };
            // distinct error positions per read (no cross-read coherence)
            for k in 0..noise_per_read {
                alts.push(span - 1 - (r as u64 * 17 + k)); // deterministic, scattered, unique-ish
            }
            ReadObs { start: 0, end: span, alts }
        }).collect()
    }

    #[test]
    fn hidden_copy_is_flagged() {
        // 60 reads, 30 carry a 20-position shared alt haplotype → a hidden copy.
        let shared: Vec<u64> = (100..120).collect();
        let rs = reads(60, 2000, 30, &shared, 3);
        let e = detect_hidden_copy(&rs, &p());
        assert!(e.flagged);
        assert_eq!(e.n_alt_positions, 20);
        assert_eq!(e.n_alt_reads, 30);
        assert!((e.alt_read_fraction - 0.5).abs() < 1e-9);
    }

    #[test]
    fn sequencing_error_is_not_flagged() {
        // No shared haplotype, just ~9 random errors per read → no candidate columns.
        let rs = reads(60, 2000, 0, &[], 9);
        let e = detect_hidden_copy(&rs, &p());
        assert!(!e.flagged);
        assert_eq!(e.n_alt_positions, 0);
    }

    #[test]
    fn heterozygous_snp_is_not_flagged() {
        // 30 of 60 reads share just 2 alt positions (a diploid het) → below min_alt_positions.
        let rs = reads(60, 2000, 30, &[500, 900], 3);
        let e = detect_hidden_copy(&rs, &p());
        assert!(!e.flagged, "a het (2 positions) must not be called a hidden copy");
        assert!(e.n_alt_positions < p().min_alt_positions);
    }

    #[test]
    fn fixed_difference_above_band_is_not_a_candidate() {
        // A position alt in ~all reads (fixed diff / reference error) is above balanced_hi.
        let rs = reads(40, 2000, 40, &[700], 0); // all 40 alt at 700 → frac 1.0 > 0.60
        let e = detect_hidden_copy(&rs, &p());
        assert_eq!(e.n_alt_positions, 0);
        assert!(!e.flagged);
    }

    #[test]
    fn low_depth_abstains() {
        let shared: Vec<u64> = (100..120).collect();
        let rs = reads(4, 2000, 2, &shared, 0); // below min_depth
        assert!(!detect_hidden_copy(&rs, &p()).flagged);
    }
}
}

// ---- merged 2026-10-05: was `vg_family/collapse_gate.rs`, now the inline module below (one component) ----
#[allow(clippy::all)]
pub mod collapse_gate {
//! Collapse gate: is a single-rep locus actually several copies the aligner could not separate?
//!
//! SDA (Vollger et al., Nat Methods 2019) detects a collapse by read-depth excess and only THEN defines PSVs,
//! "requiring sequence coverages consistent with a single-copy locus in order to distinguish PSVs from allelic
//! variants". We ran that second stage alone, and a single-copy gene (TSPYL1) reported 12 collapsed copies
//! against DAZ's 3 — see `bench/COLLAPSED_COPY_GATE.md`. **Collapse first, haplotypes second.**
//!
//! In DNA a collapse shows as excess depth. In RNA depth is copy number × expression (Clair3-RNA: "the coverage
//! is uneven across genomic regions in RNA-seq"), and allele-specific expression destroys the allelic balance
//! that would otherwise separate a het allele from a PSV ("zygosity flipping can happen"). So we keep SDA's
//! structure and change its instrument: a collapse shows instead as **reads the aligner cannot place uniquely**.
//!
//! Measured on GGO Iso-Seq (primary records only): **0 MAPQ-0 primaries across 9449 reads** at five single-copy
//! loci (TSPYL1, DERPC, ATXN7L3B, GSPT2, EEF1A1), against 19/20 at DAZ2 and 30/34 at TSPY. The statistic is
//! expression-invariant — TSPYL1 has 2151 reads and no ambiguity; DAZ2 has 20 reads and 95%.
//!
//! ⚠⚠ **DEFAULT OFF. The instrument is not what this module's name claims, and a control proved it.**
//!
//! MAPQ 0 means "this read maps equally well somewhere else". It does NOT mean "this locus is collapsed". Run
//! genome-wide, the gate fires on **EEF1A1** — whose MAPQ-0 reads align to its processed pseudogenes on
//! NC_073224.2 and NC_073227.2, other chromosomes entirely — and reports `chi(H) = 7` for a locus with one copy.
//!
//! The logic actually inverts. If a copy were truly ABSENT from the reference, its reads would pile onto the
//! present copy at HIGH mapping quality, giving depth excess and *no ambiguity at all*. That is precisely why
//! SDA detects collapses by read depth and not by mapping quality. Ambiguity detects **unresolvable paralogy**
//! (which is what `read_conflict`'s E_c oracle already does), not collapse.
//!
//! On DAZ the gate emits `chi(H) = 2`, matching the annotation — but DAZ2 is present in the reference 20 kb
//! away, so even there the signal is paralogy, not collapse. What actually failed at DAZ is assembly: DAZ2 has
//! 20 primary reads and never becomes a rep.
//!
//! Kept, off by default, because the machinery and its tests are sound and the p-value is correct for the
//! question it truly answers. Do not enable it until `chi(H)` is shown to bound copies rather than haplotypes.
//! See `bench/COLLAPSE_GATE_VALIDATION.md`.
//!
//! **STATUS:** REFUTED  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use crate::vg_family::copy_assign::poisson_binomial_upper_tail;
use crate::vg_family::readonly_copy_number::chi_h;

/// Ambiguously-placed primary reads (`k`) out of primary reads (`n`) at a locus.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Ambiguity {
    pub n: usize,
    pub k: usize,
}

/// Background per-read ambiguity rate, under a Jeffreys prior so it is never exactly zero.
///
/// `bg` is pooled over the region's uniquely-mappable reps (those that are not gate candidates). The controls
/// observe `k = 0`, whose MLE is 0 — and a zero rate would make a single stray MAPQ-0 read infinitely
/// significant. `None` when there is no background to estimate from: the caller must then ABSTAIN, never fire.
pub fn estimate_eps_amb(bg: Ambiguity) -> Option<f64> {
    if bg.n == 0 {
        return None;
    }
    Some((bg.k as f64 + 0.5) / (bg.n as f64 + 1.0))
}

/// Upper-tail probability of seeing `obs.k` or more ambiguous reads among `obs.n`, if the locus were unique.
///
/// Under the null `k ~ Binomial(n, eps_amb)`. Reuses the shipped Poisson-binomial tail with a constant
/// probability vector rather than adding a second distribution implementation. That routine is O(n²) in the
/// number of trials, and here a trial is a READ rather than a distinguishing position — but `k == 0` returns
/// immediately, and a clean locus is exactly the `k == 0` case, so uniquely-mapping loci cost nothing.
pub fn collapse_pvalue(obs: Ambiguity, eps_amb: f64) -> f64 {
    if obs.k == 0 || obs.n == 0 {
        return 1.0;
    }
    let probs = vec![eps_amb; obs.n];
    poisson_binomial_upper_tail(obs.k, &probs)
}

/// What the gate decided about a locus that has only one assembled rep.
#[derive(Clone, Debug, PartialEq)]
pub enum CollapseVerdict {
    /// Collapsed, and its reads resolve into `chi_h >= min_copies` conflicting haplotypes.
    Fire { chi_h: usize, p_value: f64 },
    /// The locus places its reads unambiguously, or its haplotypes do not reach `min_copies`.
    NotCollapsed { p_value: f64 },
    /// Cannot decide — never fire on an unbounded statistic.
    Abstain(&'static str),
}

/// Two legs, in SDA's order. Leg 2 is NOT consulted unless leg 1 fires: a single-copy gene reports plenty of
/// haplotypes (a het allele *is* a haplotype), and gating on them alone made TSPYL1 report 12 copies.
///
/// `eps_amb` is the background per-read ambiguity rate. It must be a **genome-wide** quantity, not a
/// region-local one: SDA estimates its background from "unique regions" of the genome, and for good reason —
/// in the DAZ window the only reads that fall outside DAZ1's span are DAZ2's, every one of them ambiguous, so a
/// region-local background would be ~0.95 and the gate could never fire. Measured on `GGO_mm.bam`:
/// 5785 MAPQ-0 primaries in 4,404,440, i.e. `eps_amb = 0.0013`, itself conservative because it includes the
/// genuinely collapsed loci. `None` ⇒ ABSTAIN; never fire on an unbounded statistic.
///
/// `haplotypes` are the `allele_vector`s of the **identifiable** copies at the locus.
///
/// χ(H) is a **lower bound** on copy number, never an estimate: two copies × two alleles also yields four
/// haplotypes. That is why the caller emits a copy NUMBER with reads certified tied, not an assignment.
pub fn collapse_verdict(
    obs: Ambiguity,
    eps_amb: Option<f64>,
    haplotypes: &[Vec<Option<u8>>],
    alpha: f64,
    min_copies: usize,
) -> CollapseVerdict {
    // leg 1 — is this locus collapsed at all?
    let Some(eps_amb) = eps_amb else {
        return CollapseVerdict::Abstain("no background ambiguity rate (pass --eps-amb)");
    };
    let p_value = collapse_pvalue(obs, eps_amb);
    if p_value >= alpha {
        return CollapseVerdict::NotCollapsed { p_value };
    }
    // leg 2 — and how many copies collapsed?
    let chi = chi_h(haplotypes);
    if chi < min_copies {
        return CollapseVerdict::NotCollapsed { p_value };
    }
    CollapseVerdict::Fire { chi_h: chi, p_value }
}

/// Measured background on `GGO_mm.bam`: 5785 MAPQ-0 primaries out of 4,404,440. Conservative — it includes the
/// genuinely collapsed loci, so the true unique-region rate is lower and the gate is harder to fire, not easier.
/// Recompute per sample:
/// `echo $(( $(samtools view -c -F 2308 b.bam) - $(samtools view -c -F 2308 -q 1 b.bam) ))`
pub const GENOME_WIDE_EPS_AMB: f64 = 0.001313;

#[cfg(test)]
mod tests {
    use super::*;

    /// DAZ1 itself: 22 ambiguous of 200. Against the genome-wide background this is overwhelming.
    #[test]
    fn verdict_fires_on_daz1_against_the_genome_wide_background() {
        let v = collapse_verdict(
            Ambiguity { n: 200, k: 22 },
            Some(GENOME_WIDE_EPS_AMB),
            &haps(&[b"ACGT", b"ACGA"]),
            1e-3,
            2,
        );
        assert!(matches!(v, CollapseVerdict::Fire { chi_h: 2, .. }), "DAZ1 must fire, got {v:?}");
    }

    /// Three stray ambiguous reads in a well-covered locus must NOT fire at alpha = 1e-3.
    #[test]
    fn verdict_does_not_fire_on_a_few_stray_ambiguous_reads() {
        let v = collapse_verdict(
            Ambiguity { n: 500, k: 3 },
            Some(GENOME_WIDE_EPS_AMB),
            &haps(&[b"ACGT", b"ACGA"]),
            1e-3,
            2,
        );
        assert!(matches!(v, CollapseVerdict::NotCollapsed { .. }), "3 strays must not fire, got {v:?}");
    }

    #[test]
    fn eps_amb_is_never_zero_even_when_no_background_read_is_ambiguous() {
        // The five single-copy controls give 0 ambiguous reads in 9449. The MLE is 0, under which ONE stray
        // MAPQ-0 read would be infinitely significant. Jeffreys keeps it strictly positive.
        let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
        assert!(eps > 0.0, "eps_amb must be strictly positive, got {eps}");
        assert!((eps - 0.5 / 9450.0).abs() < 1e-12, "Jeffreys: (k + 1/2) / (n + 1)");
    }

    #[test]
    fn eps_amb_abstains_without_background_reads() {
        assert_eq!(estimate_eps_amb(Ambiguity { n: 0, k: 0 }), None, "no background => cannot estimate => abstain");
    }

    #[test]
    fn collapse_pvalue_is_significant_for_daz2_and_not_for_a_clean_locus() {
        let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
        let p_daz2 = collapse_pvalue(Ambiguity { n: 20, k: 19 }, eps); // DAZ2: 19 of 20 ambiguous
        assert!(p_daz2 < 1e-6, "DAZ2 must be overwhelmingly significant, got {p_daz2}");
        let p_clean = collapse_pvalue(Ambiguity { n: 2151, k: 0 }, eps); // TSPYL1
        assert!((p_clean - 1.0).abs() < 1e-12, "k = 0 => p = 1, got {p_clean}");
    }

    #[test]
    fn collapse_pvalue_of_a_single_stray_read_is_not_significant_at_alpha() {
        let eps = estimate_eps_amb(Ambiguity { n: 9449, k: 0 }).unwrap();
        let p = collapse_pvalue(Ambiguity { n: 500, k: 1 }, eps);
        assert!(p > 1e-3, "a single stray MAPQ-0 read must not fire the gate, got {p}");
    }

    /// Allele vectors: haplotypes differing at a shared column conflict, so `chi_h` counts them separately.
    fn haps(rows: &[&[u8]]) -> Vec<Vec<Option<u8>>> {
        rows.iter().map(|r| r.iter().map(|&b| if b == b'.' { None } else { Some(b) }).collect()).collect()
    }

    /// DAZ: 19/20 ambiguous, background clean, two distinguishable haplotypes.
    #[test]
    fn verdict_fires_on_a_collapsed_locus_with_two_haplotypes() {
        let v = collapse_verdict(
            Ambiguity { n: 20, k: 19 },
            Some(GENOME_WIDE_EPS_AMB),
            &haps(&[b"ACGT", b"ACGA"]),
            1e-3,
            2,
        );
        match v {
            CollapseVerdict::Fire { chi_h, .. } => assert_eq!(chi_h, 2),
            other => panic!("expected Fire, got {other:?}"),
        }
    }

    /// TSPYL1: a single-copy gene whose reads are ALL uniquely placed. Leg 2 would report many haplotypes
    /// (het alleles, editing, isoform noise) — leg 1 must stop it before leg 2 is ever consulted.
    #[test]
    fn verdict_rejects_a_unique_locus_however_many_haplotypes_it_reports() {
        let twelve: Vec<Vec<Option<u8>>> = (0..12u8).map(|i| vec![Some(b'A' + i), Some(b'C')]).collect();
        let v = collapse_verdict(Ambiguity { n: 2151, k: 0 }, Some(GENOME_WIDE_EPS_AMB), &twelve, 1e-3, 2);
        assert!(matches!(v, CollapseVerdict::NotCollapsed { .. }), "unique locus must never fire, got {v:?}");
    }

    #[test]
    fn verdict_abstains_without_a_background_estimate() {
        let v = collapse_verdict(Ambiguity { n: 20, k: 19 }, None, &haps(&[b"AC", b"AG"]), 1e-3, 2);
        assert!(matches!(v, CollapseVerdict::Abstain(_)), "no background => abstain, got {v:?}");
    }

    /// `min_copies` applies to χ(H), not to the rep count: one haplotype is not a family.
    #[test]
    fn verdict_rejects_a_collapse_that_resolves_to_one_haplotype() {
        let v = collapse_verdict(Ambiguity { n: 20, k: 19 }, Some(GENOME_WIDE_EPS_AMB), &haps(&[b"ACGT"]), 1e-3, 2);
        assert!(
            matches!(v, CollapseVerdict::NotCollapsed { .. }),
            "chi_h = 1 < min_copies => no family, got {v:?}"
        );
    }
}
}

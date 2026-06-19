//! De-novo assembly INPUT stages (integration) — the `bench/twopass_denovo_gw_pass1.py` Pass-1
//! read-coherence skeletons + the `bench/denovo_assemble_gate.py` general-purpose assembly gate.
//!
//! These produce the `DenovoTranscript` records that the family-detection / rescue / assignment cores
//! consume. Annotation-free, minimizer-free:
//!   - **Pass-1** groups PRIMARY reads by their exact intron chain → de-novo transcript skeletons.
//!   - **Gate** keeps skeletons with enough reads, a bounded genomic span, and ALL-canonical
//!     consistent-strand junctions, then builds the spliced sequence (reverse-complemented for `-`) and
//!     bounds its length.
//!
//! BAM/FASTA I/O is the outer edge: the integration driver parses primary alignments into `PrimaryRead`
//! (via `bam::exons_from_cigar`) and loads the `GenomeIndex`; these functions are the testable transforms.
//! The gate reuses `GenomeIndex::{is_canonical_junction, fetch_sequence}` and `vg::reverse_complement`.

use anyhow::Result;
use noodles_sam::alignment::RecordBuf;

use super::family_detect::DenovoTranscript;
use crate::genome::GenomeIndex;
use crate::vg::reverse_complement;

/// Pass-1 read support to keep a skeleton (`twopass_denovo_gw_pass1.py::MIN_READS`).
pub const PASS1_MIN_READS: u32 = 2;
/// Gate read support (`denovo_assemble_gate.py::MIN_READS`).
pub const GATE_MIN_READS: u32 = 3;
/// Reject skeletons spanning more than this genomic distance (chromosome-spanning artifacts).
pub const MAX_SPAN: u64 = 3_000_000;
/// Reject spliced products shorter than this.
pub const MIN_SPLICED: usize = 100;
/// Reject implausibly long spliced products.
pub const MAX_SPLICED: usize = 300_000;

/// A parsed PRIMARY alignment (input to Pass-1). Built from the BAM by the I/O edge via
/// `bam::exons_from_cigar` (introns = the gaps between consecutive exons).
#[derive(Clone, Debug)]
pub struct PrimaryRead {
    pub chrom: String,
    pub ref_start: u64,
    pub ref_end: u64,
    /// Intron `(donor, acceptor)` chain (empty for a single-exon read).
    pub introns: Vec<(u64, u64)>,
}

/// A de-novo read-coherence skeleton (Pass-1 output): a distinct `(chrom, intron-chain)` with its support.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Skeleton {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_reads: u32,
    pub introns: Vec<(u64, u64)>,
}

/// Pass-1: group PRIMARY reads by `(chrom, intron-chain)` → skeletons, keeping groups with `>= min_reads`.
/// `start`/`end` are the min `ref_start` / max `ref_end` over the group. Deterministic (sorted by
/// `(chrom, introns)`). Mirrors `twopass_denovo_gw_pass1.py`.
pub fn pass1_skeletons(reads: &[PrimaryRead], min_reads: u32) -> Vec<Skeleton> {
    use std::collections::BTreeMap;
    // key = (chrom, intron-chain); val = (n_reads, min_start, max_end).
    let mut groups: BTreeMap<(&str, Vec<(u64, u64)>), (u32, u64, u64)> = BTreeMap::new();
    for r in reads {
        let e = groups
            .entry((r.chrom.as_str(), r.introns.clone()))
            .or_insert((0, u64::MAX, 0));
        e.0 += 1;
        e.1 = e.1.min(r.ref_start);
        e.2 = e.2.max(r.ref_end);
    }
    groups
        .into_iter()
        .filter(|(_, (n, _, _))| *n >= min_reads)
        .map(|((chrom, introns), (n, start, end))| Skeleton {
            chrom: chrom.to_string(),
            start,
            end,
            n_reads: n,
            introns,
        })
        .collect()
}

/// Build a `PrimaryRead` from a mapped PRIMARY alignment record (the Pass-1 I/O edge, mirroring the
/// python `bam.fetch()` filter). Returns `None` for unmapped / secondary / supplementary records (the ones
/// Pass-1 skips) or a record with no exons. Coordinates are 0-based (matching `GenomeIndex` and pysam's
/// `reference_start`); introns are the gaps between consecutive `bam::exons_from_cigar` exons. `chrom` is
/// the record's reference-sequence name, resolved by the caller from the header.
pub fn primary_read_from_record(record: &RecordBuf, chrom: &str) -> Option<PrimaryRead> {
    let flags = record.flags();
    if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
        return None;
    }
    let ref_start = (record.alignment_start()?.get() as u64).saturating_sub(1); // 1-based -> 0-based
    let exons = crate::bam::exons_from_cigar(ref_start, record.cigar()).ok()?;
    if exons.is_empty() {
        return None;
    }
    // introns = gaps between consecutive exons (the shared exon-parser derivation).
    let introns: Vec<(u64, u64)> = exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
    Some(PrimaryRead {
        chrom: chrom.to_string(),
        ref_start,
        ref_end: exons.last().map(|e| e.1).unwrap_or(ref_start),
        introns,
    })
}

/// Scan every PRIMARY mapped alignment in a BAM into `PrimaryRead`s (the Pass-1 input). I/O driver:
/// opens the BAM, reads the header for reference names, and applies `primary_read_from_record`.
pub fn primary_reads_from_bam(bam_path: &str, threads: usize) -> Result<Vec<PrimaryRead>> {
    let mut reader = crate::bam::open_bam(bam_path, threads.max(1))?;
    let header = reader.read_header()?;
    let mut record = RecordBuf::default();
    let mut out = Vec::new();
    while reader.read_record_buf(&header, &mut record)? > 0 {
        let chrom = match record
            .reference_sequence_id()
            .and_then(|id| header.reference_sequences().get_index(id))
        {
            Some((name, _)) => format!("{name}"),
            None => continue,
        };
        if let Some(pr) = primary_read_from_record(&record, &chrom) {
            out.push(pr);
        }
    }
    Ok(out)
}

/// Gate parameters (defaults mirror `denovo_assemble_gate.py`).
#[derive(Clone, Copy, Debug)]
pub struct GateParams {
    pub min_reads: u32,
    pub max_span: u64,
    pub min_spliced: usize,
    pub max_spliced: usize,
}

impl Default for GateParams {
    fn default() -> Self {
        GateParams {
            min_reads: GATE_MIN_READS,
            max_span: MAX_SPAN,
            min_spliced: MIN_SPLICED,
            max_spliced: MAX_SPLICED,
        }
    }
}

/// Classify a junction's transcription strand from its donor/acceptor dinucleotides: `Some('+')` for a
/// plus-canonical motif (GT-AG/GC-AG/AT-AC), `Some('-')` for the reverse-strand motif, else `None`.
fn junction_strand(genome: &GenomeIndex, chrom: &str, donor: u64, acceptor: u64) -> Option<char> {
    if genome.is_canonical_junction(chrom, donor, acceptor, '+') {
        Some('+')
    } else if genome.is_canonical_junction(chrom, donor, acceptor, '-') {
        Some('-')
    } else {
        None
    }
}

/// Build a transcript's spliced sequence + strand from its exon structure `[start, end)` + intron chain:
/// require ALL-canonical consistent-strand junctions, fetch each exon, uppercase (so a later
/// `reverse_complement`, which maps lowercase → N, is safe), and reverse-complement for a `-` strand.
/// Returns `None` if a junction is non-canonical/inconsistent or a fetch fails. Length gating is the
/// caller's. Shared by the assemble gate and the rescue thin-locus scan.
pub fn build_spliced_seq(
    genome: &GenomeIndex,
    chrom: &str,
    start: u64,
    end: u64,
    introns: &[(u64, u64)],
) -> Option<(Vec<u8>, char)> {
    let mut strand: Option<char> = None;
    for &(d, a) in introns {
        match junction_strand(genome, chrom, d, a) {
            Some(st) => {
                if strand.is_some_and(|s0| s0 != st) {
                    return None;
                }
                strand = Some(st);
            }
            None => return None,
        }
    }
    let mut exons = Vec::with_capacity(introns.len() + 1);
    let mut prev = start;
    for &(d, a) in introns {
        exons.push((prev, d));
        prev = a;
    }
    exons.push((prev, end));
    let mut seq = Vec::new();
    for &(xs, xe) in &exons {
        let s = genome.fetch_sequence(chrom, xs, xe)?;
        seq.extend(s.iter().map(|b| b.to_ascii_uppercase()));
    }
    let strand = strand.unwrap_or('+');
    if strand == '-' {
        seq = reverse_complement(&seq);
    }
    Some((seq, strand))
}

/// Assemble gate: keep skeletons with `>= min_reads` reads, span `<= max_span`, and ALL-canonical
/// consistent-strand junctions; build the spliced sequence (reverse-complemented for a `-` strand) and
/// require its length in `[min_spliced, max_spliced]`. Returns the gated `DenovoTranscript`s in input
/// order. Mirrors `denovo_assemble_gate.py`.
pub fn assemble_gate(skeletons: &[Skeleton], genome: &GenomeIndex, p: &GateParams) -> Vec<DenovoTranscript> {
    let mut out = Vec::new();
    for sk in skeletons {
        if sk.n_reads < p.min_reads {
            continue;
        }
        if sk.end.saturating_sub(sk.start) > p.max_span {
            continue;
        }
        let (seq, strand) = match build_spliced_seq(genome, &sk.chrom, sk.start, sk.end, &sk.introns) {
            Some(v) => v,
            None => continue,
        };
        if seq.len() < p.min_spliced || seq.len() > p.max_spliced {
            continue;
        }
        let n_exon = sk.introns.len() + 1;
        out.push(DenovoTranscript {
            tid: format!("DN_{}_{}_{}", sk.chrom, sk.start, n_exon),
            chrom: sk.chrom.clone(),
            start: sk.start,
            end: sk.end,
            n_reads: sk.n_reads,
            strand,
            introns: sk.introns.clone(),
            seq,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pr(chrom: &str, start: u64, end: u64, introns: &[(u64, u64)]) -> PrimaryRead {
        PrimaryRead { chrom: chrom.into(), ref_start: start, ref_end: end, introns: introns.to_vec() }
    }

    use noodles_core::Position;
    use noodles_sam::alignment::record::cigar::{op::Kind, Op};
    use noodles_sam::alignment::record::Flags;
    use noodles_sam::alignment::record_buf::Cigar;

    /// Build a minimal mapped record with the given flags, 1-based alignment start, and CIGAR ops.
    fn rec(flags: Flags, start_1b: usize, ops: Vec<Op>) -> RecordBuf {
        let cigar: Cigar = ops.into_iter().collect();
        RecordBuf::builder()
            .set_flags(flags)
            .set_alignment_start(Position::try_from(start_1b).unwrap())
            .set_cigar(cigar)
            .build()
    }

    // ---- primary_read_from_record ----

    #[test]
    fn record_to_primary_read_0based_with_introns() {
        // 80M20N80M at 1-based 101 -> 0-based start 100, exon1 [100,180), intron [180,200), exon2 [200,280).
        let r = rec(
            Flags::default(),
            101,
            vec![Op::new(Kind::Match, 80), Op::new(Kind::Skip, 20), Op::new(Kind::Match, 80)],
        );
        let p = primary_read_from_record(&r, "c1").expect("primary mapped read");
        assert_eq!(p.chrom, "c1");
        assert_eq!(p.ref_start, 100);
        assert_eq!(p.ref_end, 280);
        assert_eq!(p.introns, vec![(180, 200)]);
    }

    #[test]
    fn record_to_primary_read_single_exon_no_introns() {
        let r = rec(Flags::default(), 101, vec![Op::new(Kind::Match, 100)]);
        let p = primary_read_from_record(&r, "c1").expect("primary mapped read");
        assert_eq!((p.ref_start, p.ref_end), (100, 200));
        assert!(p.introns.is_empty());
    }

    /// Real-data smoke check for the BAM adapter + Pass-1 (no FASTA needed). Ignored by default; run with:
    ///   RUSTLE_DENOVO_SMOKE_BAM=/path/to.bam cargo test --lib -- --ignored smoke_pass1
    #[test]
    #[ignore = "needs a real BAM via RUSTLE_DENOVO_SMOKE_BAM"]
    fn smoke_pass1_on_real_bam() {
        let bam = match std::env::var("RUSTLE_DENOVO_SMOKE_BAM") {
            Ok(p) => p,
            Err(_) => return,
        };
        let reads = primary_reads_from_bam(&bam, 4).expect("read the BAM");
        assert!(!reads.is_empty(), "a real BAM should yield primary reads");
        let sk = pass1_skeletons(&reads, PASS1_MIN_READS);
        let multi = sk.iter().filter(|s| !s.introns.is_empty()).count();
        assert!(!sk.is_empty(), "primary reads should form read-coherence skeletons");
        eprintln!(
            "smoke: {} primary reads -> {} skeletons (>= {} reads); {} multi-exon",
            reads.len(),
            sk.len(),
            PASS1_MIN_READS,
            multi
        );
    }

    #[test]
    fn record_to_primary_read_skips_secondary_supplementary_unmapped() {
        let ops = || vec![Op::new(Kind::Match, 80), Op::new(Kind::Skip, 20), Op::new(Kind::Match, 80)];
        assert!(primary_read_from_record(&rec(Flags::SECONDARY, 101, ops()), "c1").is_none());
        assert!(primary_read_from_record(&rec(Flags::SUPPLEMENTARY, 101, ops()), "c1").is_none());
        assert!(primary_read_from_record(&rec(Flags::UNMAPPED, 101, ops()), "c1").is_none());
    }
    fn skel(chrom: &str, start: u64, end: u64, n: u32, introns: &[(u64, u64)]) -> Skeleton {
        Skeleton { chrom: chrom.into(), start, end, n_reads: n, introns: introns.to_vec() }
    }
    /// Genome: exon1 [0,80), intron [80,100) with the given dinucleotides, exon2 [100,180). 160 bp spliced.
    fn genome_one_intron(donor: &[u8; 2], acc: &[u8; 2]) -> GenomeIndex {
        let mut s = vec![b'A'; 180];
        s[80] = donor[0];
        s[81] = donor[1];
        s[98] = acc[0];
        s[99] = acc[1];
        GenomeIndex::from_seqs(&[("c1", &s)])
    }

    // ---- pass1_skeletons ----

    #[test]
    fn pass1_groups_same_intron_chain() {
        let reads = [
            pr("c1", 100, 500, &[(200, 300)]),
            pr("c1", 110, 520, &[(200, 300)]),
        ];
        let sk = pass1_skeletons(&reads, 2);
        assert_eq!(sk, vec![skel("c1", 100, 520, 2, &[(200, 300)])]);
    }

    #[test]
    fn pass1_separates_distinct_chains() {
        let reads = [
            pr("c1", 100, 500, &[(200, 300)]),
            pr("c1", 100, 500, &[(250, 350)]),
        ];
        assert_eq!(pass1_skeletons(&reads, 1).len(), 2);
    }

    #[test]
    fn pass1_min_reads_filters_singletons() {
        let reads = [
            pr("c1", 100, 500, &[(200, 300)]),
            pr("c1", 100, 500, &[(200, 300)]),
            pr("c1", 100, 500, &[(250, 350)]), // only 1 read
        ];
        let sk = pass1_skeletons(&reads, 2);
        assert_eq!(sk.len(), 1, "the 1-read chain is dropped at min_reads=2");
        assert_eq!(sk[0].introns, vec![(200, 300)]);
    }

    #[test]
    fn pass1_single_exon_reads_group() {
        let reads = [pr("c1", 100, 500, &[]), pr("c1", 120, 480, &[])];
        let sk = pass1_skeletons(&reads, 2);
        assert_eq!(sk.len(), 1);
        assert!(sk[0].introns.is_empty());
        assert_eq!((sk[0].start, sk[0].end), (100, 500));
    }

    #[test]
    fn pass1_different_chrom_separate() {
        let reads = [pr("c1", 100, 500, &[(200, 300)]), pr("c2", 100, 500, &[(200, 300)])];
        assert_eq!(pass1_skeletons(&reads, 1).len(), 2);
    }

    // ---- assemble_gate ----

    #[test]
    fn gate_plus_strand_builds_transcript() {
        let g = genome_one_intron(b"GT", b"AG");
        let out = assemble_gate(&[skel("c1", 0, 180, 3, &[(80, 100)])], &g, &GateParams::default());
        assert_eq!(out.len(), 1);
        let t = &out[0];
        assert_eq!(t.tid, "DN_c1_0_2");
        assert_eq!(t.n_reads, 3);
        assert_eq!(t.introns, vec![(80, 100)]);
        assert_eq!(t.seq.len(), 160, "exon1 80 + exon2 80");
        assert_eq!(t.strand, '+');
        assert!(t.seq.iter().all(|&b| b == b'A'), "plus strand: not reverse-complemented");
    }

    #[test]
    fn gate_minus_strand_reverse_complements() {
        // CT-AC is the reverse-strand motif -> strand '-' -> spliced seq is reverse-complemented.
        let g = genome_one_intron(b"CT", b"AC");
        let out = assemble_gate(&[skel("c1", 0, 180, 3, &[(80, 100)])], &g, &GateParams::default());
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].strand, '-');
        assert!(out[0].seq.iter().all(|&b| b == b'T'), "revcomp of all-A is all-T");
    }

    #[test]
    fn gate_rejects_noncanonical_junction() {
        let g = genome_one_intron(b"AA", b"AG"); // donor AA is not a canonical motif
        assert!(assemble_gate(&[skel("c1", 0, 180, 3, &[(80, 100)])], &g, &GateParams::default()).is_empty());
    }

    #[test]
    fn gate_rejects_below_min_reads() {
        let g = genome_one_intron(b"GT", b"AG");
        assert!(assemble_gate(&[skel("c1", 0, 180, 2, &[(80, 100)])], &g, &GateParams::default()).is_empty());
    }

    #[test]
    fn gate_rejects_oversized_span() {
        let g = genome_one_intron(b"GT", b"AG");
        let p = GateParams { max_span: 100, ..GateParams::default() };
        assert!(assemble_gate(&[skel("c1", 0, 180, 3, &[(80, 100)])], &g, &p).is_empty());
    }

    #[test]
    fn gate_rejects_short_spliced() {
        // exon1 [0,10), intron [10,30), exon2 [30,40): 20 bp spliced < MIN_SPLICED.
        let mut s = vec![b'A'; 40];
        s[10] = b'G';
        s[11] = b'T';
        s[28] = b'A';
        s[29] = b'G';
        let g = GenomeIndex::from_seqs(&[("c1", &s)]);
        assert!(assemble_gate(&[skel("c1", 0, 40, 3, &[(10, 30)])], &g, &GateParams::default()).is_empty());
    }

    #[test]
    fn gate_single_exon_kept() {
        // no introns -> no junction check, strand '+', one exon [0,180).
        let g = GenomeIndex::from_seqs(&[("c1", &vec![b'A'; 180])]);
        let out = assemble_gate(&[skel("c1", 0, 180, 3, &[])], &g, &GateParams::default());
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].tid, "DN_c1_0_1");
        assert_eq!(out[0].seq.len(), 180);
    }

    #[test]
    fn gate_rejects_inconsistent_strand() {
        // intron1 [60,80) GT-AG ('+'), intron2 [140,160) CT-AC ('-') -> mixed strand -> rejected.
        let mut s = vec![b'A'; 220];
        s[60] = b'G';
        s[61] = b'T';
        s[78] = b'A';
        s[79] = b'G';
        s[140] = b'C';
        s[141] = b'T';
        s[158] = b'A';
        s[159] = b'C';
        let g = GenomeIndex::from_seqs(&[("c1", &s)]);
        let sk = skel("c1", 0, 220, 3, &[(60, 80), (140, 160)]);
        assert!(assemble_gate(&[sk], &g, &GateParams::default()).is_empty());
    }
}

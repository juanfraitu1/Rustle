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
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::RecordBuf;

use super::copy_split::AlignedRead;
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

/// AS-tie gate (the phantom-avoidance kernel for collapsed-copy recovery). Given a region's alignments as
/// `(read_name, is_secondary, alignment_score, PrimaryRead)`, emit the SECONDARY alignments whose AS is within
/// `as_ratio` of that read's BEST score in the region — i.e. genuine TIED co-located copies (minimap2 flagged
/// a sibling primary and this one secondary), NOT homology-shadow spillover (AS far below the primary). The
/// read's primary is already counted in the main set; these tied secondaries restore the STARVED sibling
/// copy's read support so the rescue can assemble it. A truly indistinguishable (no copy-specific feature)
/// phantom is still rejected downstream by the identifiability / PSV gate — this only feeds candidates in.
pub fn tied_secondary_reads(aln: &[(String, bool, i32, PrimaryRead)], as_ratio: f64) -> Vec<PrimaryRead> {
    use std::collections::HashMap;
    let mut best: HashMap<&str, i32> = HashMap::new();
    for (name, _, as_, _) in aln {
        let e = best.entry(name.as_str()).or_insert(i32::MIN);
        if *as_ > *e {
            *e = *as_;
        }
    }
    aln.iter()
        .filter(|(name, is_sec, as_, _)| *is_sec && (*as_ as f64) >= best[name.as_str()] as f64 * as_ratio)
        .map(|(_, _, _, pr)| pr.clone())
        .collect()
}

/// Read the `AS:i` alignment score from a record (the gate signal). `None` if absent.
fn record_as(record: &RecordBuf) -> Option<i32> {
    use noodles_sam::alignment::record::data::field::{Tag, Value};
    for entry in noodles_sam::alignment::Record::data(record).iter() {
        let (tag, value) = entry.ok()?;
        if tag == Tag::ALIGNMENT_SCORE {
            return match value {
                Value::Int8(v) => Some(v as i32),
                Value::UInt8(v) => Some(v as i32),
                Value::Int16(v) => Some(v as i32),
                Value::UInt16(v) => Some(v as i32),
                Value::Int32(v) => Some(v),
                Value::UInt32(v) => Some(v as i32),
                _ => None,
            };
        }
    }
    None
}

/// Read the `de:f` (gap-compressed per-base divergence) tag from a record — the conflict-criterion signal.
/// `None` if absent. `de` is a custom 2-char tag carrying a float.
fn record_de(record: &RecordBuf) -> Option<f32> {
    use noodles_sam::alignment::record::data::field::Tag;
    use noodles_sam::alignment::record_buf::data::field::Value;
    let de_tag = Tag::new(b'd', b'e');
    match record.data().get(&de_tag)? {
        Value::Float(v) => Some(*v),
        _ => None,
    }
}

/// Like `primary_read_from_record` but ALSO accepts SECONDARY alignments, returning `(PrimaryRead, read_name,
/// is_secondary, AS)` for the AS-tie gate. Unmapped / supplementary / no-AS / no-exon records are skipped.
pub fn any_read_from_record(record: &RecordBuf, chrom: &str) -> Option<(PrimaryRead, String, bool, i32)> {
    let flags = record.flags();
    if flags.is_unmapped() || flags.is_supplementary() {
        return None;
    }
    let ref_start = (record.alignment_start()?.get() as u64).saturating_sub(1);
    let exons = crate::bam::exons_from_cigar(ref_start, record.cigar()).ok()?;
    if exons.is_empty() {
        return None;
    }
    let introns: Vec<(u64, u64)> = exons.windows(2).map(|w| (w[0].1, w[1].0)).collect();
    let pr = PrimaryRead {
        chrom: chrom.to_string(),
        ref_start,
        ref_end: exons.last().map(|e| e.1).unwrap_or(ref_start),
        introns,
    };
    let name = record.name().map(|n| n.to_string()).unwrap_or_default();
    let as_ = record_as(record)?;
    Some((pr, name, flags.is_secondary(), as_))
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

fn cigar_kind_to_char(k: Kind) -> char {
    match k {
        Kind::Match => 'M',
        Kind::Insertion => 'I',
        Kind::Deletion => 'D',
        Kind::Skip => 'N',
        Kind::SoftClip => 'S',
        Kind::HardClip => 'H',
        Kind::Pad => 'P',
        Kind::SequenceMatch => '=',
        Kind::SequenceMismatch => 'X',
    }
}

/// A mapped alignment for copy assignment: its `AlignedRead`, reference name, mapping quality (the
/// silver-standard "unique mapper" signal: `mapq > 0`), read name (any ground-truth label), and the `AS:i`
/// alignment score (the conflict-graph tie criterion; 0 if absent from the BAM record).
#[derive(Clone, Debug)]
pub struct BamRead {
    pub chrom: String,
    pub read: AlignedRead,
    pub mapq: u8,
    pub name: String,
    pub as_score: i32,
    /// minimap2 `de:f` gap-compressed per-base divergence (0.0 if absent) — the conflict-tie signal.
    pub de: f32,
    /// chimeric/split alignment (SAM flag 0x800) — excluded from conflict placements.
    pub is_supplementary: bool,
}

/// Build an `AlignedRead` (ref_start 0-based, CIGAR ops as chars, read sequence) + mapping quality + read
/// NAME + `AS:i` alignment score from a mapped record — the per-read input copy ASSIGNMENT consumes
/// (`copy_assign_pipeline`). The sequence keeps soft-clipped bases (excludes hard-clips), matching
/// `copy_split::allele_at`'s CIGAR walk. `as_score` is 0 if the tag is absent. `None` if unmapped.
pub fn aligned_read_from_record(record: &RecordBuf) -> Option<(AlignedRead, u8, String, i32, f32, bool)> {
    if record.flags().is_unmapped() {
        return None;
    }
    let ref_start = (record.alignment_start()?.get() as u64).saturating_sub(1);
    let cigar: Vec<(char, u64)> = record
        .cigar()
        .as_ref()
        .iter()
        .map(|op| (cigar_kind_to_char(op.kind()), op.len() as u64))
        .collect();
    let seq: Vec<u8> = record.sequence().as_ref().to_vec();
    let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
    let name = record.name().map(|n| n.to_string()).unwrap_or_default();
    let as_score = record_as(record).unwrap_or(0);
    let de = record_de(record).unwrap_or(0.0);
    let is_supplementary = record.flags().is_supplementary();
    Some((AlignedRead { ref_start, cigar, seq }, mapq, name, as_score, de, is_supplementary))
}

/// Scan every mapped alignment in a BAM into `BamRead`s (the copy-assignment read input). I/O driver,
/// mirroring `primary_reads_from_bam`. Includes secondary/supplementary (multimappers are exactly what
/// assignment resolves); the chrom resolves from the header.
pub fn aligned_reads_from_bam(bam_path: &str, threads: usize) -> Result<Vec<BamRead>> {
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
        if let Some((read, mapq, name, as_score, de, is_supplementary)) = aligned_read_from_record(&record) {
            out.push(BamRead { chrom, read, mapq, name, as_score, de, is_supplementary });
        }
    }
    Ok(out)
}

/// Collect BOTH the primary reads (detection input) and ALL mapped reads (assignment input, incl.
/// secondary/supplementary multimappers) overlapping `[lo, hi)` on `chrom`. Uses the `.bai`-indexed region
/// query when available (fast — reads only the region, not the whole file), falling back to a full scan.
pub fn reads_in_region(
    bam_path: &str,
    chrom: &str,
    lo: u64,
    hi: u64,
    threads: usize,
) -> Result<(Vec<PrimaryRead>, Vec<BamRead>)> {
    match reads_in_region_indexed(bam_path, chrom, lo, hi) {
        Ok(r) => Ok(r),
        Err(_) => reads_in_region_scan(bam_path, chrom, lo, hi, threads),
    }
}

/// `.bai`-indexed region query (fast path). Errs if the index is missing/unreadable so the caller can scan.
fn reads_in_region_indexed(
    bam_path: &str,
    chrom: &str,
    lo: u64,
    hi: u64,
) -> Result<(Vec<PrimaryRead>, Vec<BamRead>)> {
    let bai_path = format!("{bam_path}.bai");
    anyhow::ensure!(std::path::Path::new(&bai_path).exists(), "no .bai index");
    let mut reader = noodles_bam::io::reader::Builder::default().build_from_path(bam_path)?;
    let header = reader.read_header()?;
    let index = noodles_bam::bai::read(&bai_path)?;
    let region: noodles_core::Region = format!("{chrom}:{}-{}", lo + 1, hi).parse()?;
    let query = reader.query(&header, &index, &region)?;
    let mut primary = Vec::new();
    let mut bam_reads = Vec::new();
    for result in query {
        let record = result?;
        let rb = RecordBuf::try_from_alignment_record(&header, &record)?;
        if let Some(pr) = primary_read_from_record(&rb, chrom) {
            primary.push(pr);
        }
        if let Some((read, mapq, name, as_score, de, is_supplementary)) = aligned_read_from_record(&rb) {
            bam_reads.push(BamRead { chrom: chrom.to_string(), read, mapq, name, as_score, de, is_supplementary });
        }
    }
    Ok((primary, bam_reads))
}

/// Collect AS-TIED SECONDARY reads overlapping `[lo, hi)` on `chrom`, as `PrimaryRead`s — the starved-copy
/// read support for collapsed-copy recovery (see [`tied_secondary_reads`]). Indexed query; errs without a
/// `.bai`. These augment the RESCUE input only (additive), never the primary-only main detection pass.
pub fn tied_secondary_reads_in_region(
    bam_path: &str,
    chrom: &str,
    lo: u64,
    hi: u64,
    as_ratio: f64,
) -> Result<Vec<PrimaryRead>> {
    let bai_path = format!("{bam_path}.bai");
    anyhow::ensure!(std::path::Path::new(&bai_path).exists(), "no .bai index");
    let mut reader = noodles_bam::io::reader::Builder::default().build_from_path(bam_path)?;
    let header = reader.read_header()?;
    let index = noodles_bam::bai::read(&bai_path)?;
    let region: noodles_core::Region = format!("{chrom}:{}-{}", lo + 1, hi).parse()?;
    let query = reader.query(&header, &index, &region)?;
    let mut aln = Vec::new();
    for result in query {
        let rb = RecordBuf::try_from_alignment_record(&header, &result?)?;
        if let Some((pr, name, is_sec, as_)) = any_read_from_record(&rb, chrom) {
            aln.push((name, is_sec, as_, pr));
        }
    }
    Ok(tied_secondary_reads(&aln, as_ratio))
}

/// PRIMARY mapped `AlignedRead`s (seq + CIGAR) overlapping `[lo, hi)` on `chrom` — the ASJ scan input
/// (secondary/supplementary excluded, matching the python). Indexed query; errs without a `.bai`.
pub fn primary_aligned_reads_in_region(
    bam_path: &str,
    chrom: &str,
    lo: u64,
    hi: u64,
) -> Result<Vec<crate::vg_family::copy_split::AlignedRead>> {
    let bai_path = format!("{bam_path}.bai");
    anyhow::ensure!(std::path::Path::new(&bai_path).exists(), "no .bai index");
    let mut reader = noodles_bam::io::reader::Builder::default().build_from_path(bam_path)?;
    let header = reader.read_header()?;
    let index = noodles_bam::bai::read(&bai_path)?;
    let region: noodles_core::Region = format!("{chrom}:{}-{}", lo + 1, hi).parse()?;
    let query = reader.query(&header, &index, &region)?;
    let mut out = Vec::new();
    for result in query {
        let rb = RecordBuf::try_from_alignment_record(&header, &result?)?;
        let f = rb.flags();
        if f.is_unmapped() || f.is_secondary() || f.is_supplementary() {
            continue;
        }
        if let Some((read, _, _, _, _, _)) = aligned_read_from_record(&rb) {
            out.push(read);
        }
    }
    Ok(out)
}

/// Fraction of PRIMARY reads covering `pos` (0-based) that are MAPQ-0 (multimapping) — the ASJ confound
/// control. A HIGH fraction at a het-ASJ anchor flags a collapsed-paralog masquerade (the two "alleles" are
/// paralog copies); a LOW fraction means genuine within-gene heterozygosity. (For copy-specific junctions the
/// reading inverts: multimapping is expected.) Capped at 600 reads, matching the python.
pub fn frac_mq0_at(bam_path: &str, chrom: &str, pos: u64) -> Result<f64> {
    let bai_path = format!("{bam_path}.bai");
    anyhow::ensure!(std::path::Path::new(&bai_path).exists(), "no .bai index");
    let mut reader = noodles_bam::io::reader::Builder::default().build_from_path(bam_path)?;
    let header = reader.read_header()?;
    let index = noodles_bam::bai::read(&bai_path)?;
    let region: noodles_core::Region = format!("{chrom}:{}-{}", pos + 1, pos + 1).parse()?;
    let query = reader.query(&header, &index, &region)?;
    let (mut n, mut z) = (0u32, 0u32);
    for result in query {
        let rb = RecordBuf::try_from_alignment_record(&header, &result?)?;
        let f = rb.flags();
        if f.is_unmapped() || f.is_secondary() || f.is_supplementary() {
            continue;
        }
        n += 1;
        if rb.mapping_quality().map(|q| q.get()).unwrap_or(0) == 0 {
            z += 1;
        }
        if n >= 600 {
            break;
        }
    }
    Ok(if n > 0 { z as f64 / n as f64 } else { 0.0 })
}

/// Full-scan fallback (no index): one pass, region-filtered.
fn reads_in_region_scan(
    bam_path: &str,
    chrom: &str,
    lo: u64,
    hi: u64,
    threads: usize,
) -> Result<(Vec<PrimaryRead>, Vec<BamRead>)> {
    let mut reader = crate::bam::open_bam(bam_path, threads.max(1))?;
    let header = reader.read_header()?;
    let mut record = RecordBuf::default();
    let mut primary = Vec::new();
    let mut bam_reads = Vec::new();
    while reader.read_record_buf(&header, &mut record)? > 0 {
        let rchrom = match record
            .reference_sequence_id()
            .and_then(|id| header.reference_sequences().get_index(id))
        {
            Some((name, _)) => format!("{name}"),
            None => continue,
        };
        if rchrom != chrom {
            continue;
        }
        let start = match record.alignment_start() {
            Some(p) => (p.get() as u64).saturating_sub(1),
            None => continue,
        };
        let end = start + record.cigar().alignment_span() as u64;
        if start >= hi || end <= lo {
            continue;
        }
        if let Some(pr) = primary_read_from_record(&record, chrom) {
            primary.push(pr);
        }
        if let Some((read, mapq, name, as_score, de, is_supplementary)) = aligned_read_from_record(&record) {
            bam_reads.push(BamRead { chrom: chrom.to_string(), read, mapq, name, as_score, de, is_supplementary });
        }
    }
    Ok((primary, bam_reads))
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

    #[test]
    fn tied_secondary_kept_spillover_dropped() {
        // r1: a TIED multimapper (primary @ copyA AS=500, secondary @ copyB AS=500) -> the secondary restores
        // copyB's support. r2: spillover (secondary AS 300 << 500) -> dropped. r3: primary-only -> never emitted.
        let aln = vec![
            ("r1".to_string(), false, 500, pr("c1", 0, 100, &[])),
            ("r1".to_string(), true, 500, pr("c1", 1000, 1100, &[])),
            ("r2".to_string(), false, 500, pr("c1", 0, 100, &[])),
            ("r2".to_string(), true, 300, pr("c1", 1000, 1100, &[])),
            ("r3".to_string(), false, 400, pr("c1", 0, 100, &[])),
        ];
        let out = tied_secondary_reads(&aln, 0.98);
        assert_eq!(out.len(), 1, "only the tied secondary (r1) survives");
        assert_eq!(out[0].ref_start, 1000, "and it sits at the STARVED copyB position");
    }

    #[test]
    fn tied_secondary_margin_boundary() {
        // secondary AS 495 vs best 500: kept at ratio 0.98 (>=490), dropped at 0.999 (<499.5).
        let aln = vec![
            ("r".to_string(), false, 500, pr("c1", 5, 105, &[])),
            ("r".to_string(), true, 495, pr("c1", 5, 105, &[])),
        ];
        assert_eq!(tied_secondary_reads(&aln, 0.98).len(), 1);
        assert_eq!(tied_secondary_reads(&aln, 0.999).len(), 0);
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
    fn aligned_read_from_record_extracts_cigar_and_seq() {
        use noodles_sam::alignment::record_buf::Sequence;
        let cigar: Cigar =
            vec![Op::new(Kind::Match, 5), Op::new(Kind::Skip, 10), Op::new(Kind::Match, 5)]
                .into_iter()
                .collect();
        let seq: Sequence = Sequence::from(b"AAAAACCCCC".to_vec());
        let r = RecordBuf::builder()
            .set_flags(Flags::default())
            .set_alignment_start(Position::try_from(101usize).unwrap())
            .set_cigar(cigar)
            .set_sequence(seq)
            .build();
        let (ar, _mapq, _name, _as_score, _de, _is_supp) = aligned_read_from_record(&r).expect("mapped read");
        assert_eq!(ar.ref_start, 100); // 1-based 101 -> 0-based
        assert_eq!(ar.cigar, vec![('M', 5), ('N', 10), ('M', 5)]);
        assert_eq!(ar.seq, b"AAAAACCCCC".to_vec());
    }

    #[test]
    fn aligned_read_from_record_skips_unmapped() {
        let r = RecordBuf::builder().set_flags(Flags::UNMAPPED).build();
        assert!(aligned_read_from_record(&r).is_none());
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

    // ---- record_de ----

    #[test]
    fn record_de_parses_float_tag() {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value as BufValue;
        let mut record = RecordBuf::builder().set_flags(Flags::default()).build();
        record.data_mut().insert(Tag::new(b'd', b'e'), BufValue::Float(0.0123));
        let de = record_de(&record).expect("de tag present");
        assert!((de - 0.0123).abs() < 1e-6);
    }

    #[test]
    fn record_de_absent_is_none() {
        let record = RecordBuf::builder().set_flags(Flags::default()).build();
        assert!(record_de(&record).is_none());
    }

    // ---- aligned_read_from_record — de + is_supplementary ----

    #[test]
    fn aligned_read_from_record_extracts_de_and_supplementary() {
        use noodles_sam::alignment::record::data::field::Tag;
        use noodles_sam::alignment::record_buf::data::field::Value as BufValue;
        use noodles_sam::alignment::record_buf::Sequence;
        let cigar: Cigar = vec![Op::new(Kind::Match, 4)].into_iter().collect();
        let mut record = RecordBuf::builder()
            .set_flags(Flags::SUPPLEMENTARY)
            .set_alignment_start(Position::try_from(1usize).unwrap())
            .set_cigar(cigar)
            .set_sequence(Sequence::from(b"ACGT".to_vec()))
            .build();
        record.data_mut().insert(Tag::new(b'd', b'e'), BufValue::Float(0.02));
        let (_ar, _mapq, _name, _as, de, is_supp) =
            aligned_read_from_record(&record).expect("mapped");
        assert!((de - 0.02).abs() < 1e-6);
        assert!(is_supp, "supplementary flag must be surfaced");
    }
}

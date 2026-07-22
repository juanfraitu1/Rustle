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
use crate::vg_family::seq_utils::reverse_complement;

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
    pass1_skeletons_robust(reads, min_reads, 1)
}

/// Like `pass1_skeletons`, but the transcript EXTENT must be reached by at least `min_terminal_support`
/// reads, so a single runaway read (a chimeric / intra-primed / mis-clipped terminal exon) cannot
/// artificially inflate the length. With `min_terminal_support == 1` this is exactly `pass1_skeletons`
/// (the union min-start / max-end). With `k > 1` the 5' boundary is the `k`-th smallest read start and the
/// 3' boundary is the `k`-th largest read end (the furthest position SUPPORTED by `k` reads), trimming the
/// outermost `k-1` outlier reads at each end. For FLNC IsoSeq (full-length reads), the true ends are reached
/// by many reads, so a small `k` trims only runaways and leaves real boundaries intact. If a group has fewer
/// than `k` reads the boundary falls back to the outermost available read.
pub fn pass1_skeletons_robust(reads: &[PrimaryRead], min_reads: u32, min_terminal_support: u32) -> Vec<Skeleton> {
    use std::collections::BTreeMap;
    let k = min_terminal_support.max(1) as usize;
    // key = (chrom, intron-chain); val = (n_reads, k-smallest starts asc, k-largest ends desc).
    let mut groups: BTreeMap<(&str, Vec<(u64, u64)>), (u32, Vec<u64>, Vec<u64>)> = BTreeMap::new();
    for r in reads {
        let e = groups
            .entry((r.chrom.as_str(), r.introns.clone()))
            .or_insert((0, Vec::new(), Vec::new()));
        e.0 += 1;
        // keep the k smallest starts (ascending)
        let pos = e.1.partition_point(|&x| x <= r.ref_start);
        if pos < k {
            e.1.insert(pos, r.ref_start);
            e.1.truncate(k);
        }
        // keep the k largest ends (descending)
        let pos = e.2.partition_point(|&x| x >= r.ref_end);
        if pos < k {
            e.2.insert(pos, r.ref_end);
            e.2.truncate(k);
        }
    }
    groups
        .into_iter()
        .filter(|(_, (n, _, _))| *n >= min_reads)
        .map(|((chrom, introns), (n, starts, ends))| {
            // robust boundary = the k-th supported value (or the outermost available if the group is smaller).
            let si = k.min(starts.len()).saturating_sub(1);
            let ei = k.min(ends.len()).saturating_sub(1);
            Skeleton {
                chrom: chrom.to_string(),
                start: starts[si],
                end: ends[ei],
                n_reads: n,
                introns,
            }
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
/// unique-mapper agreement's "unique mapper" signal: `mapq > 0`), read name (any ground-truth label), the `AS:i`
/// alignment score (the conflict-graph tie criterion; 0 if absent from the BAM record), the `de:f`
/// gap-compressed per-base divergence (the de-tie conflict criterion; 0.0 if absent), and
/// `is_supplementary` (chimeric/split flag — excluded from conflict placements).
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
    /// secondary alignment (SAM flag 0x100) — a multimapping read's non-primary placement; excluded,
    /// together with `is_supplementary`, before any per-read CIGAR statistic (project invariant: primary
    /// alignments only, `-F 2308`), so a single physical molecule is never double-counted as two witnesses.
    pub is_secondary: bool,
}

/// Build an `AlignedRead` (ref_start 0-based, CIGAR ops as chars, read sequence) + mapping quality + read
/// NAME + `AS:i` alignment score + `de:f` divergence + `is_supplementary` + `is_secondary` from a mapped
/// record — the per-read input copy ASSIGNMENT consumes (`copy_assign_pipeline`). The sequence keeps
/// soft-clipped bases (excludes hard-clips), matching `copy_split::allele_at`'s CIGAR walk. `as_score` is
/// 0 if the tag is absent. `None` if unmapped.
pub fn aligned_read_from_record(record: &RecordBuf) -> Option<(AlignedRead, u8, String, i32, f32, bool, bool)> {
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
    // Per-base Phred qualities, parallel to `seq` (both keep soft-clips). Empty if the BAM
    // carried no quality string -> the PSV likelihood falls back to the flat error rate.
    let qual: Vec<u8> = record.quality_scores().as_ref().to_vec();
    let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);
    let name = record.name().map(|n| n.to_string()).unwrap_or_default();
    let as_score = record_as(record).unwrap_or(0);
    // `de` is universally present on minimap2 output; absent → 0.0 is a fail-soft default.
    // Two genuinely-absent placements would both be 0.0 and would tie, but this is acceptable because
    // `de` is always present on the production BAM. A de-less BAM cannot use the de-tie criterion.
    let de = record_de(record).unwrap_or(0.0);
    let is_supplementary = record.flags().is_supplementary();
    let is_secondary = record.flags().is_secondary();
    Some((AlignedRead { ref_start, cigar, seq, qual }, mapq, name, as_score, de, is_supplementary, is_secondary))
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
        if let Some((read, mapq, name, as_score, de, is_supplementary, is_secondary)) = aligned_read_from_record(&record) {
            out.push(BamRead { chrom, read, mapq, name, as_score, de, is_supplementary, is_secondary });
        }
    }
    Ok(out)
}

/// Whether a record's flags mark it an UNMAPPED PRIMARY read — the filter behind
/// `unmapped_reads_from_bam`, i.e. the mirror image of `primary_read_from_record`'s mapped-primary
/// filter: unmapped is required, secondary/supplementary are excluded (a "primary" record here means
/// neither secondary nor supplementary, matching SAM's use of the term for unmapped records too).
fn is_unmapped_primary(flags: noodles_sam::alignment::record::Flags) -> bool {
    flags.is_unmapped() && !flags.is_secondary() && !flags.is_supplementary()
}

/// From a slice of records, collect `(read_name, sequence)` for every UNMAPPED PRIMARY record with a
/// non-empty sequence — the testable transform behind `unmapped_reads_from_bam` (this is the input to
/// the VG re-align supplement's unmapped-read minimizer routing stage). Unit-tested directly against
/// in-memory `RecordBuf`s, without needing a BAM fixture with an unmapped record.
fn collect_unmapped(records: &[RecordBuf]) -> Vec<(String, Vec<u8>)> {
    records
        .iter()
        .filter(|r| is_unmapped_primary(r.flags()))
        .filter_map(|r| {
            let seq: Vec<u8> = r.sequence().as_ref().to_vec();
            if seq.is_empty() {
                return None;
            }
            let name = r.name().map(|n| n.to_string()).unwrap_or_default();
            Some((name, seq))
        })
        .collect()
}

/// Scan every UNMAPPED PRIMARY record in a BAM into `(read_name, sequence)` pairs. I/O driver, mirroring
/// `primary_reads_from_bam`'s reader setup. Unmapped reads carry no alignment position, so there is no
/// `.bai` region query to use here — this is always a full-file scan. Applies the identical
/// `collect_unmapped` filter to each record while streaming, so the whole BAM is never buffered.
pub fn unmapped_reads_from_bam(bam_path: &str, threads: usize) -> Result<Vec<(String, Vec<u8>)>> {
    let mut reader = crate::bam::open_bam(bam_path, threads.max(1))?;
    let header = reader.read_header()?;
    let mut record = RecordBuf::default();
    let mut out = Vec::new();
    while reader.read_record_buf(&header, &mut record)? > 0 {
        out.extend(collect_unmapped(std::slice::from_ref(&record)));
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
        if let Some((read, mapq, name, as_score, de, is_supplementary, is_secondary)) = aligned_read_from_record(&rb) {
            bam_reads.push(BamRead { chrom: chrom.to_string(), read, mapq, name, as_score, de, is_supplementary, is_secondary });
        }
    }
    Ok((primary, bam_reads))
}

/// A reusable indexed-BAM handle that parses the `.bai` index and header ONCE, so a long sequence of region
/// queries (the per-region copy-assignment loop, thousands of regions at genome scale) does not re-parse the
/// multi-MB index for every region. Re-opening the file handle per query is ~free; parsing the index is not.
/// `reads_in_region` returns exactly what the free function [`reads_in_region`]'s indexed path returns.
pub struct BamIndexCache {
    header: noodles_sam::Header,
    index: noodles_bam::bai::Index,
}

impl BamIndexCache {
    /// Parse `{bam_path}.bai` and the header once. Errs (so the caller can fall back to per-region opens or a
    /// full scan) when the index is missing or unreadable.
    pub fn open(bam_path: &str) -> Result<Self> {
        let bai_path = format!("{bam_path}.bai");
        anyhow::ensure!(std::path::Path::new(&bai_path).exists(), "no .bai index");
        let mut reader = noodles_bam::io::reader::Builder::default().build_from_path(bam_path)?;
        let header = reader.read_header()?;
        let index = noodles_bam::bai::read(&bai_path)?;
        Ok(Self { header, index })
    }

    /// Query one region using the cached header+index (the reader seeks via the index, independent of the
    /// freshly-opened handle's position). Byte-identical to `reads_in_region_indexed`'s output.
    pub fn reads_in_region(
        &self,
        bam_path: &str,
        chrom: &str,
        lo: u64,
        hi: u64,
    ) -> Result<(Vec<PrimaryRead>, Vec<BamRead>)> {
        let mut reader = noodles_bam::io::reader::Builder::default().build_from_path(bam_path)?;
        let region: noodles_core::Region = format!("{chrom}:{}-{}", lo + 1, hi).parse()?;
        let query = reader.query(&self.header, &self.index, &region)?;
        let mut primary = Vec::new();
        let mut bam_reads = Vec::new();
        for result in query {
            let record = result?;
            let rb = RecordBuf::try_from_alignment_record(&self.header, &record)?;
            if let Some(pr) = primary_read_from_record(&rb, chrom) {
                primary.push(pr);
            }
            if let Some((read, mapq, name, as_score, de, is_supplementary, is_secondary)) = aligned_read_from_record(&rb) {
                bam_reads.push(BamRead { chrom: chrom.to_string(), read, mapq, name, as_score, de, is_supplementary, is_secondary });
            }
        }
        Ok((primary, bam_reads))
    }
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
        if let Some((read, _, _, _, _, _, _)) = aligned_read_from_record(&rb) {
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
        if let Some((read, mapq, name, as_score, de, is_supplementary, is_secondary)) = aligned_read_from_record(&record) {
            bam_reads.push(BamRead { chrom: chrom.to_string(), read, mapq, name, as_score, de, is_supplementary, is_secondary });
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
    /// Test `min_reads` against the LOCUS (junction-incidence component) rather than one isoform. Default on.
    /// `false` reproduces the pre-fix per-isoform gate exactly. See [`locus_support`].
    pub pool_locus_support: bool,
}

impl Default for GateParams {
    fn default() -> Self {
        GateParams {
            pool_locus_support: true,
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

/// Read support of the LOCUS each skeleton belongs to: connected components of the **junction-incidence
/// graph**, where two skeletons are adjacent iff they share an exact `(chrom, donor, acceptor)`. Returns a
/// vector parallel to `skeletons`, each entry the summed `n_reads` of that skeleton's component.
///
/// Why this exists. `assemble_gate` used to test `min_reads` against a single intron chain, i.e. against ONE
/// ISOFORM. A locus expressed as several minor isoforms is then shattered: at DAZ2 the 12 spliced primary
/// reads fragment into 9 distinct chains whose best support is 2, so every chain died at `GATE_MIN_READS = 3`
/// and the locus was never assembled — even though all 12 reads share the terminal junction
/// `(42939630, 42943604)` and none of them shares a junction with DAZ1. Pooling first and gating second keeps
/// the locus. This replaces a per-isoform threshold with a graph operation; the threshold that remains applies
/// to a locus, which is the object it was always meant to describe.
///
/// A single-exon skeleton has NO junctions, so it is adjacent to nothing and forms its own component: its
/// locus support is its own read count, exactly as before. The empty intron chain never pools — which matters,
/// because keying on it would union every unspliced read on a chromosome.
/// A skeleton is a CHIMERIC BRIDGE if it shares a junction with two skeletons whose genomic spans are
/// DISJOINT — it splices across two separate loci, so it belongs to neither.
///
/// Without this test, junction pooling merges paralogs. At GSTM a 2-read spliced transcript spanning
/// `129191743-129222260` carries junctions of BOTH GSTM5 (`129191742-129197751`) and GSTM1
/// (`129216297-129222748`). It bridged them in the incidence graph, inherited their combined support, cleared
/// the gate, and span-aware collapse then merged GSTM5 into it — destroying a real annotated copy.
///
/// Exact, not thresholded: two spans either intersect or they do not. Isoforms of ONE locus overlap each other,
/// so a genuine full-length isoform — even one that is the only chain carrying some junction pair — is never
/// flagged. DAZ2's nine chains all overlap within `42879944-42943604` and are unaffected.
fn is_chimeric_bridge(i: usize, skeletons: &[Skeleton], neighbours: &[Vec<usize>]) -> bool {
    let disjoint = |a: &Skeleton, b: &Skeleton| a.chrom != b.chrom || a.end <= b.start || b.end <= a.start;
    let nb = &neighbours[i];
    nb.iter().enumerate().any(|(x, &j)| nb[x + 1..].iter().any(|&k| disjoint(&skeletons[j], &skeletons[k])))
}

/// Read support of the LOCUS each skeleton belongs to (see the doc above). Chimeric bridges are excluded from
/// pooling and keep only their own read count, so a 2-read chimera can never inherit a real locus's support.
pub fn locus_support(skeletons: &[Skeleton]) -> Vec<u32> {
    // adjacency: two skeletons share an exact (chrom, donor, acceptor)
    let mut by_junction: std::collections::BTreeMap<(&str, u64, u64), Vec<usize>> = std::collections::BTreeMap::new();
    for (i, sk) in skeletons.iter().enumerate() {
        for &(d, a) in &sk.introns {
            by_junction.entry((sk.chrom.as_str(), d, a)).or_default().push(i);
        }
    }
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); skeletons.len()];
    for members in by_junction.values() {
        for (x, &i) in members.iter().enumerate() {
            for &j in &members[x + 1..] {
                neighbours[i].push(j);
                neighbours[j].push(i);
            }
        }
    }
    for nb in neighbours.iter_mut() {
        nb.sort_unstable();
        nb.dedup();
    }
    let chimeric: Vec<bool> = (0..skeletons.len()).map(|i| is_chimeric_bridge(i, skeletons, &neighbours)).collect();

    let mut parent: Vec<usize> = (0..skeletons.len()).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    for i in 0..skeletons.len() {
        if chimeric[i] {
            continue; // a chimera joins nothing; it is its own locus
        }
        for &j in &neighbours[i] {
            if chimeric[j] {
                continue;
            }
            let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
            if ri != rj {
                parent[ri] = rj;
            }
        }
    }
    let mut sum: std::collections::HashMap<usize, u32> = std::collections::HashMap::new();
    for (i, sk) in skeletons.iter().enumerate() {
        let r = find(&mut parent, i);
        *sum.entry(r).or_insert(0) += sk.n_reads;
    }
    (0..skeletons.len())
        .map(|i| if chimeric[i] { skeletons[i].n_reads } else { let r = find(&mut parent, i); sum[&r] })
        .collect()
}

/// Assemble gate: keep skeletons whose **locus** (junction-incidence component, see [`locus_support`]) has
/// `>= min_reads` reads, span `<= max_span`, and ALL-canonical consistent-strand junctions; build the spliced
/// sequence (reverse-complemented for a `-` strand) and require its length in `[min_spliced, max_spliced]`.
/// Returns the gated `DenovoTranscript`s in input order. Mirrors `denovo_assemble_gate.py`, except that the
/// read-count test is applied to the locus rather than to one isoform (`p.pool_locus_support`).
/// THE FORMAL TRANSCRIPT DEFINITION. A transcript is an exact-intron-chain cluster of primary reads whose LOCUS
/// (the junction-incidence component, see `locus_support`) carries >= `min_reads` (GATE_MIN_READS) reads, whose
/// junctions are all canonical + consistent-strand, and whose spliced length is in `[min_spliced, max_spliced]`.
/// A locus is a copy family of size >= 1; single-copy is the chi(H)=1 boundary case (see the `single_copy` module).
pub fn assemble_gate(skeletons: &[Skeleton], genome: &GenomeIndex, p: &GateParams) -> Vec<DenovoTranscript> {
    let support: Vec<u32> =
        if p.pool_locus_support { locus_support(skeletons) } else { skeletons.iter().map(|s| s.n_reads).collect() };
    let mut out = Vec::new();
    for (i, sk) in skeletons.iter().enumerate() {
        if support[i] < p.min_reads {
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
         ..Default::default() });
    }
    out
}

/// Position-aware seeding for UNSPLICED (empty-intron-chain) reads: single-linkage span-overlap
/// clustering per chromosome, so each pseudogene/retrocopy copy seeds its OWN skeleton instead of
/// every unspliced read on a chromosome pooling into one `(chrom, [])` group that overruns
/// MAX_SPLICED. Threshold-free (pure overlap), mirroring `thin_loci`'s single-linkage-by-span.
pub fn cluster_unspliced(reads: &[PrimaryRead], min_reads: u32, k: usize) -> Vec<Skeleton> {
    use std::collections::BTreeMap;
    let mut by_chrom: BTreeMap<&str, Vec<&PrimaryRead>> = BTreeMap::new();
    for r in reads {
        if r.introns.is_empty() {
            by_chrom.entry(r.chrom.as_str()).or_default().push(r);
        }
    }
    let mut out = Vec::new();
    for (chrom, mut rs) in by_chrom {
        rs.sort_by_key(|r| (r.ref_start, r.ref_end));
        let mut i = 0;
        while i < rs.len() {
            // single-linkage: extend the cluster while the next read starts before the running max end
            let mut cluster_end = rs[i].ref_end;
            let mut j = i + 1;
            while j < rs.len() && rs[j].ref_start < cluster_end {
                cluster_end = cluster_end.max(rs[j].ref_end);
                j += 1;
            }
            let cluster = &rs[i..j];
            let n = cluster.len() as u32;
            if n >= min_reads {
                // per-cluster robust boundaries: k-th smallest start, k-th largest end (matches
                // pass1_skeletons_robust's boundary rule).
                let mut starts: Vec<u64> = cluster.iter().map(|r| r.ref_start).collect();
                let mut ends: Vec<u64> = cluster.iter().map(|r| r.ref_end).collect();
                starts.sort_unstable();
                ends.sort_unstable_by(|a, b| b.cmp(a));
                let si = k.min(starts.len()).saturating_sub(1);
                let ei = k.min(ends.len()).saturating_sub(1);
                out.push(Skeleton {
                    chrom: chrom.to_string(),
                    start: starts[si],
                    end: ends[ei],
                    n_reads: n,
                    introns: Vec::new(),
                });
            }
            i = j;
        }
    }
    out
}

#[cfg(test)]
mod locus_support_tests {
    use super::*;

    fn sk(chrom: &str, start: u64, end: u64, n_reads: u32, introns: &[(u64, u64)]) -> Skeleton {
        Skeleton { chrom: chrom.into(), start, end, n_reads, introns: introns.to_vec() }
    }

    /// The real DAZ2 shape. Its 12 spliced primary reads fragment into 9 intron chains whose best support is
    /// 2, so every chain died at `GATE_MIN_READS = 3`. They all share the terminal junction
    /// (42939630, 42943604), so the junction-incidence graph pools them into one locus with support 12.
    #[test]
    fn locus_support_pools_isoform_fragments_sharing_a_junction() {
        const TERMINAL: (u64, u64) = (42_939_630, 42_943_604);
        let sks: Vec<Skeleton> = (0..3)
            .map(|i| {
                let d = 42_900_000 + i * 1_000;
                sk("NC_073248.2", 42_879_944, 42_943_604, 2, &[(d, d + 500), TERMINAL])
            })
            .collect();
        let sup = locus_support(&sks);
        assert_eq!(sup, vec![6, 6, 6], "three 2-read chains sharing a junction are ONE locus with 6 reads");
        assert!(sup.iter().all(|&s| s >= GATE_MIN_READS), "the pooled locus clears the gate that killed each chain");
    }

    /// DAZ1 and DAZ2 share ZERO junctions (measured). Pooling must not bridge them.
    #[test]
    fn locus_support_does_not_bridge_loci_with_no_shared_junction() {
        let daz1 = sk("NC_073248.2", 42_783_133, 42_859_657, 44, &[(42_800_000, 42_801_000)]);
        let daz2 = sk("NC_073248.2", 42_879_944, 42_943_604, 2, &[(42_939_630, 42_943_604)]);
        assert_eq!(locus_support(&[daz1, daz2]), vec![44, 2], "disjoint junction sets => disjoint loci");
    }

    /// THE safety property. A single-exon skeleton has no junctions, so it is adjacent to nothing and its
    /// locus support is its own read count. Keying on the empty intron chain would union every unspliced read
    /// on a chromosome (measured: 746 reads across 44.3 Mbp), which is exactly what must not happen.
    #[test]
    fn locus_support_never_pools_single_exon_skeletons() {
        let sks = vec![
            sk("c1", 100, 200, 1, &[]),
            sk("c1", 5_000_000, 5_000_100, 1, &[]),
            sk("c1", 900, 1_000, 1, &[]),
        ];
        assert_eq!(locus_support(&sks), vec![1, 1, 1], "the empty intron chain must never pool");
    }

    /// The real GSTM regression. A 2-read chain carrying junctions of BOTH GSTM5 and GSTM1 -- two genes whose
    /// spans are disjoint -- must not inherit their pooled support, or it clears the gate and span-aware
    /// collapse merges the two real copies into one locus.
    #[test]
    fn locus_support_refuses_to_pool_through_a_chimeric_bridge() {
        const J5: (u64, u64) = (129_193_000, 129_194_000);
        const J1: (u64, u64) = (129_218_000, 129_219_000);
        let gstm5 = sk("c1", 129_191_742, 129_197_751, 40, &[J5]);
        let gstm1 = sk("c1", 129_216_297, 129_222_748, 90, &[J1]);
        let bridge = sk("c1", 129_191_743, 129_222_260, 2, &[J5, J1]); // spliced GSTM5->GSTM1 readthrough
        let sup = locus_support(&[gstm5, gstm1, bridge]);
        assert_eq!(sup, vec![40, 90, 2], "the chimera keeps only its own 2 reads; the paralogs stay separate");
        assert!(sup[2] < GATE_MIN_READS, "and it therefore dies at the gate");
    }

    /// A genuine full-length isoform joins chains whose spans OVERLAP, so it is not a chimera and does pool.
    #[test]
    fn locus_support_pools_through_a_full_length_isoform_whose_neighbours_overlap() {
        let a = sk("c1", 1_000, 5_000, 2, &[(1_500, 2_000)]);
        let hub = sk("c1", 1_000, 6_000, 2, &[(1_500, 2_000), (4_000, 4_500)]);
        let b = sk("c1", 1_200, 6_000, 2, &[(4_000, 4_500)]); // overlaps `a`, so `hub` bridges nothing disjoint
        assert_eq!(locus_support(&[a, hub, b]), vec![6, 6, 6], "one locus, three isoforms");
    }

    #[test]
    fn locus_support_is_transitive_through_a_chain_of_shared_junctions() {
        // a—b share J1; b—c share J2; so a, b, c are one locus even though a and c share nothing.
        let a = sk("c1", 0, 100, 1, &[(10, 20)]);
        let b = sk("c1", 0, 200, 2, &[(10, 20), (30, 40)]);
        let c = sk("c1", 0, 300, 4, &[(30, 40)]);
        assert_eq!(locus_support(&[a, b, c]), vec![7, 7, 7]);
    }

    #[test]
    fn locus_support_keeps_same_junction_on_different_chroms_separate() {
        let a = sk("c1", 0, 100, 3, &[(10, 20)]);
        let b = sk("c2", 0, 100, 5, &[(10, 20)]);
        assert_eq!(locus_support(&[a, b]), vec![3, 5], "a junction is (chrom, donor, acceptor)");
    }
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
        let (ar, _mapq, _name, _as_score, _de, _is_supp, _is_secondary) = aligned_read_from_record(&r).expect("mapped read");
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
    fn pass1_robust_trims_a_runaway_terminal_read() {
        // four FLNC reads agree on the locus ~100..520; one runaway read extends far on both ends. With
        // k=2 the boundary is the 2nd-most-extreme position, so the single runaway cannot inflate the length.
        let reads = [
            pr("c1", 100, 500, &[(200, 300)]),
            pr("c1", 102, 515, &[(200, 300)]),
            pr("c1", 105, 520, &[(200, 300)]),
            pr("c1", 1, 9000, &[(200, 300)]), // runaway: starts way left, ends way right
        ];
        // k=1 (legacy union) lets the runaway dominate: 1..9000.
        assert_eq!(pass1_skeletons_robust(&reads, 2, 1), vec![skel("c1", 1, 9000, 4, &[(200, 300)])]);
        // k=2 trims it to the 2nd-smallest start (100) and 2nd-largest end (520).
        assert_eq!(pass1_skeletons_robust(&reads, 2, 2), vec![skel("c1", 100, 520, 4, &[(200, 300)])]);
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
        let (_ar, _mapq, _name, _as, de, is_supp, _is_secondary) =
            aligned_read_from_record(&record).expect("mapped");
        assert!((de - 0.02).abs() < 1e-6);
        assert!(is_supp, "supplementary flag must be surfaced");
    }

    // ---- collect_unmapped (the unmapped_reads_from_bam filter, unit-tested without a BAM fixture) ----

    /// Build a minimal record with the given flags, name, and sequence (no alignment position — mirrors
    /// an unmapped record straight off a BAM).
    fn unmapped_rec(flags: Flags, name: &str, seq: &[u8]) -> RecordBuf {
        use noodles_sam::alignment::record_buf::Sequence;
        RecordBuf::builder()
            .set_flags(flags)
            .set_name(name)
            .set_sequence(Sequence::from(seq.to_vec()))
            .build()
    }

    #[test]
    fn collect_unmapped_keeps_only_unmapped_primary() {
        let records = vec![
            // unmapped primary -> kept.
            unmapped_rec(Flags::UNMAPPED, "read_unmapped", b"ACGTACGT"),
            // mapped primary -> excluded (not unmapped).
            unmapped_rec(Flags::default(), "read_mapped", b"TTTTTTTT"),
            // unmapped but flagged secondary -> excluded (not a primary record).
            unmapped_rec(Flags::UNMAPPED | Flags::SECONDARY, "read_unmapped_secondary", b"GGGGGGGG"),
            // unmapped but flagged supplementary -> excluded.
            unmapped_rec(Flags::UNMAPPED | Flags::SUPPLEMENTARY, "read_unmapped_supp", b"CCCCCCCC"),
            // unmapped primary with an empty sequence -> excluded.
            unmapped_rec(Flags::UNMAPPED, "read_unmapped_noseq", b""),
        ];
        let out = collect_unmapped(&records);
        assert_eq!(
            out,
            vec![("read_unmapped".to_string(), b"ACGTACGT".to_vec())],
            "only the unmapped PRIMARY record with a non-empty sequence survives"
        );
    }

    #[test]
    fn is_unmapped_primary_flag_matrix() {
        assert!(is_unmapped_primary(Flags::UNMAPPED));
        assert!(!is_unmapped_primary(Flags::default()));
        assert!(!is_unmapped_primary(Flags::UNMAPPED | Flags::SECONDARY));
        assert!(!is_unmapped_primary(Flags::UNMAPPED | Flags::SUPPLEMENTARY));
    }

    // ---- cluster_unspliced ----

    #[test]
    fn cluster_unspliced_separates_distant_loci() {
        // two unspliced read groups 50kb apart -> TWO skeletons, not one chromosome-spanning giant
        let reads = vec![
            pr("c1", 1000, 3000, &[]), pr("c1", 1100, 3100, &[]), pr("c1", 1200, 3200, &[]),
            pr("c1", 51000, 53000, &[]), pr("c1", 51100, 53100, &[]), pr("c1", 51200, 53200, &[]),
        ];
        let sk = cluster_unspliced(&reads, 3, 1);
        assert_eq!(sk.len(), 2, "distant unspliced loci must seed as separate skeletons");
        assert!(sk.iter().all(|s| s.end - s.start < 300_000), "each skeleton spans one locus, not the chromosome");
        assert!(sk.iter().all(|s| s.introns.is_empty()));
    }

    #[test]
    fn cluster_unspliced_merges_overlapping_and_filters_min_reads() {
        // one overlapping pile -> ONE skeleton; a lone read below min_reads -> dropped
        let reads = vec![
            pr("c1", 100, 400, &[]), pr("c1", 150, 450, &[]), pr("c1", 200, 500, &[]),
            pr("c1", 90000, 90200, &[]), // singleton, below min_reads=2
        ];
        let sk = cluster_unspliced(&reads, 2, 1);
        assert_eq!(sk.len(), 1, "overlapping unspliced reads merge; the singleton is below min_reads");
        assert_eq!(sk[0].n_reads, 3);
    }
}

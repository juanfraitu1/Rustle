//! BAM reading: alignment blocks, exons, strand, weight.

use anyhow::Result;
use noodles_bam::io::Reader as BamReader;
use noodles_sam::alignment::{
    record::cigar::op::Kind,
    record::data::field::Tag,
    record::Cigar as _,
    RecordBuf,
};
use noodles_sam::alignment::record_buf::Cigar;
use std::io;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use crate::polya;
use crate::types::{BundleRead, Junction};

static NEXT_READ_UID: AtomicU64 = AtomicU64::new(1);

#[inline]
fn fnv1a64(s: &str) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in s.as_bytes() {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

#[derive(Debug, Default)]
struct CigarCompatParse {
    exons: Vec<(u64, u64)>,
    /// Per-intron deletion offsets: (left_del_before_N, right_del_after_N).
    junc_del_offsets: Vec<(u64, u64)>,
    clip_left: u32,
    clip_right: u32,
    insertion_sites: Vec<u64>,
}

/// Parse CIGAR with CIGAR coordinates following the original algorithm.
///
/// Coordinates are returned as 0-based half-open intervals.
fn parse_cigar_compat(ref_start: u64, cigar: &Cigar) -> CigarCompatParse {
    let mut out = CigarCompatParse::default();

    // Byte-for-byte port of StringTie's `GSamRecord::setupCoordinates`
    // (gclib/GSam.cpp:197-291). `juncsdel` entries are pushed by BOTH the
    // D-branch AND the M-branch whenever `intron` is true — for `N D M`
    // patterns this produces TWO entries per intron. We intentionally
    // replicate that shape so `juncsdel[i-1]` indexing (rlink.cpp:1226-1230)
    // matches StringTie's output byte-for-byte, including the cases where
    // the extra entry misaligns indexing for reads with `N ... D ... M ... N`
    // patterns — the goal is parity, not correctness-over-StringTie.
    let mut l = 0u64;
    let mut exstart = ref_start;
    let mut exon_started = false;
    let mut intron = false;
    let mut ins = false;
    let mut del = 0u64;
    let mut prevdel = 0u64;

    // StringTie's `GSeg(s, e)` constructor silently SWAPS `s` and `e` if `s > e`
    // so that `start <= end` always holds (gclib/GBase.h:381). That means the
    // `juncsdel` entries ST pushes are not semantically `(left_del, right_del)`
    // but rather `(min(a, b), max(a, b))`. Replicate exactly.
    let push_gseg = |out: &mut CigarCompatParse, a: u64, b: u64| {
        let (s, e) = if a > b { (b, a) } else { (a, b) };
        out.junc_del_offsets.push((s, e));
    };

    for result in cigar.iter() {
        let Ok(op) = result else {
            continue;
        };
        let kind = op.kind();
        let len = op.len() as u64;

        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                exon_started = true;
                l = l.saturating_add(len);
                if intron {
                    // GSam.cpp:223-226 — push GSeg(prevdel, 0) on first M after N.
                    push_gseg(&mut out, prevdel, 0);
                }
                intron = false;
                ins = false;
                del = 0;
            }
            Kind::Deletion => {
                del = len;
                l = l.saturating_add(del);
                if intron {
                    // GSam.cpp:234-237 — D immediately after N (possibly following I).
                    // `intron` stays true; the subsequent M will also push (prevdel, 0).
                    push_gseg(&mut out, prevdel, del);
                }
                ins = false;
            }
            Kind::Insertion => {
                // Same reference position tracking as mismatch_anchor CIGAR pass.
                out.insertion_sites.push(ref_start.saturating_add(l));
                ins = true;
            }
            Kind::Skip => {
                // guard: anomalous leading intron before any exon.
                if !exon_started {
                    // The original algorithm only skips this CIGAR op;
                    // it does not abort parsing the rest of the alignment.
                    continue;
                }
                // GSam.cpp:247 — `if(!ins || !intron)` closes the preceding exon.
                if !ins || !intron {
                    let exon_end = ref_start.saturating_add(l);
                    if exon_end > exstart {
                        out.exons.push((exstart, exon_end));
                    }
                }
                l = l.saturating_add(len);
                exstart = ref_start.saturating_add(l);
                prevdel = del;
                intron = true;
                del = 0;
            }
            Kind::SoftClip => {
                if l > 0 {
                    out.clip_right = len as u32;
                } else {
                    out.clip_left = len as u32;
                }
                ins = false;
                del = 0;
            }
            Kind::HardClip => {
                ins = false;
                del = 0;
            }
            Kind::Pad => {}
        }
    }

    // GSam.cpp:283-287 — close the final exon only if the read did not end in an intron.
    if !intron {
        let exon_end = ref_start.saturating_add(l);
        if exon_end > exstart {
            out.exons.push((exstart, exon_end));
        }
    }

    out
}

fn deletion_aware_junctions_from_parts(
    exons: &[(u64, u64)],
    junc_del_offsets: &[(u64, u64)],
) -> Vec<Junction> {
    let raw = junctions_from_exons(exons);
    let mut out = Vec::with_capacity(raw.len());
    for (i, j) in raw.into_iter().enumerate() {
        let (left_del, right_del) = junc_del_offsets.get(i).copied().unwrap_or((0, 0));
        out.push(Junction::new(
            j.donor.saturating_sub(left_del),
            j.acceptor.saturating_add(right_del),
        ));
    }
    out
}

fn parse_tag(tag_name: &str) -> Option<Tag> {
    let bytes = tag_name.as_bytes();
    if bytes.len() == 2 {
        Some(Tag::new(bytes[0], bytes[1]))
    } else {
        None
    }
}

fn exons_from_cigar_impl(ref_start: u64, cigar: &Cigar) -> io::Result<Vec<(u64, u64)>> {
    Ok(parse_cigar_compat(ref_start, cigar).exons)
}

/// Exons = alignment blocks (N = intron splits). 0-based [start, end).
pub fn exons_from_cigar(ref_start: u64, cigar: &Cigar) -> Result<Vec<(u64, u64)>> {
    Ok(exons_from_cigar_impl(ref_start, cigar)?)
}

/// Extract raw junctions from exons: (exon[i].1, exon[i+1].0).
pub fn junctions_from_exons(exons: &[(u64, u64)]) -> Vec<Junction> {
    let mut j = Vec::new();
    for i in 0..exons.len().saturating_sub(1) {
        j.push(Junction::new(exons[i].1, exons[i + 1].0));
    }
    j
}

/// NH tag (alignment hit count) for multi-mapper weight: weight = 1/NH ( 899; original algorithm 1/NH).
fn get_nh(record: &RecordBuf) -> u32 {
    let data = record.data();
    data.get(&Tag::ALIGNMENT_HIT_COUNT)
        .and_then(|v| v.as_int().and_then(|n| u32::try_from(n).ok()))
        .unwrap_or(1)
}

fn get_tag_int(record: &RecordBuf, tag_name: &str) -> Option<u32> {
    let tag = parse_tag(tag_name)?;
    let data = record.data();
    data.get(&tag)
        .and_then(|v| v.as_int().and_then(|n| u32::try_from(n).ok()))
}

/// HI tag (alignment hit index), used by the original algorithm to disambiguate mate matching keys.
fn get_hi(record: &RecordBuf) -> u32 {
    get_tag_int(record, "HI").unwrap_or(0)
}

/// NM/nM edit distance with adjustments.
///
/// behavior
/// - `nm = tag_int("NM")`
/// - if `nm == 0`, fallback to `tag_int("nM")`
/// - if fallback came from `nM` and read is paired, divide by 2
/// - add +1 if left clipped, +1 if right clipped
fn get_nm_adjusted(record: &RecordBuf, clip_left: u32, clip_right: u32) -> u32 {
    let mut nm = get_tag_int(record, "NM").unwrap_or(0);
    if nm == 0 {
        nm = get_tag_int(record, "nM").unwrap_or(0);
        if record.flags().is_segmented() {
            nm /= 2;
        }
    }
    if clip_left > 0 {
        nm = nm.saturating_add(1);
    }
    if clip_right > 0 {
        nm = nm.saturating_add(1);
    }
    nm
}

/// YC tag (alignment count): the original algorithm rdcount uses YC when present, then scales by NH.
fn get_yc(record: &RecordBuf) -> f64 {
    let Some(tag) = parse_tag("YC") else {
        return 1.0;
    };
    let data = record.data();
    let Some(v) = data.get(&tag) else {
        return 1.0;
    };
    match v {
        noodles_sam::alignment::record_buf::data::field::Value::Float(f) => {
            let x = *f as f64;
            if x > 0.0 {
                return x;
            }
        }
        _ => {
            if let Some(i) = v.as_int() {
                let x = i as f64;
                if x > 0.0 {
                    return x;
                }
            }
        }
    }
    1.0
}

/// YK tag: unitig coverage (unitig_cov=brec.tag_float("YK")).
fn get_yk(record: &RecordBuf) -> f64 {
    let Some(tag) = parse_tag("YK") else {
        return 0.0;
    };
    let data = record.data();
    let Some(v) = data.get(&tag) else {
        return 0.0;
    };
    match v {
        noodles_sam::alignment::record_buf::data::field::Value::Float(f) => {
            let x = *f as f64;
            if x > 0.0 {
                return x;
            }
        }
        _ => {
            if let Some(i) = v.as_int() {
                let x = i as f64;
                if x > 0.0 {
                    return x;
                }
            }
        }
    }
    0.0
}

/// TS tag (minimap2 transcript strand): ts:A:+ or ts:A:-. Combined with is_reverse gives transcript strand (reference script:713-724).
/// Returns None if absent; Some(true)=+ strand, Some(false)=- strand. For future use when assigning read to bundle by transcript strand.
#[allow(dead_code)]
fn get_ts_strand(record: &RecordBuf) -> Option<bool> {
    let tag = Tag::new(b't', b's');
    let data = record.data();
    let v = data.get(&tag)?;
    let ts = match v {
        noodles_sam::alignment::record_buf::data::field::Value::Character(c) => char::from(*c),
        noodles_sam::alignment::record_buf::data::field::Value::String(s) => {
            let bytes: &[u8] = s.as_ref();
            char::from(*bytes.first()?)
        }
        _ => return None,
    };
    let is_rev = record.flags().is_reverse_complemented();
    // Gate: RUSTLE_TS_DIRECT=1 treats ts tag as direct transcript strand
    // (minimap2's ts is already relative to reference; no XOR needed).
    // Default: match StringTie's GSam.cpp:spliceStrand() XOR with reverse flag.
    if std::env::var_os("RUSTLE_TS_DIRECT").is_some() {
        match ts {
            '+' => Some(true),
            '-' => Some(false),
            _ => None,
        }
    } else {
        match ts {
            '+' => Some(!is_rev),
            '-' => Some(is_rev),
            _ => None,
        }
    }
}

#[inline]
fn inferred_read_strand(record: &RecordBuf) -> char {
    if let Some(is_plus) = get_ts_strand(record) {
        if is_plus {
            '+'
        } else {
            '-'
        }
    } else {
        // can keep unknown strand (0 / '.') and infer it later from
        // junction/polyA evidence; do not force orientation fallback here.
        '.'
    }
}

#[inline]
fn seq_bytes_from_record(record: &RecordBuf) -> Vec<u8> {
    record.sequence()
        .as_ref()
        .iter()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' | b'C' | b'G' | b'T' | b'U' | b'N' => b.to_ascii_uppercase(),
            _ => b'N',
        })
        .collect()
}

const POLYA_MIN_CONSEC: u32 = 5;
const POLYA_WINDOW: usize = 20;
const POLYA_MIN_FRAC: f64 = 0.8;

// detect_polya_from_record removed — inlined at call site to share seq_bytes buffer.

/// check_last_exon_polyA/check_first_exon_polyT:
/// ratio of A on last exon, ratio of T on first exon (>=0.8 triggers trim candidate).
/// Exons are 0-based half-open, so exon length is end-start.
fn detect_terminal_exon_poly_flags_with_seq(
    seq: &[u8],
    exons: &[(u64, u64)],
    clip_left: u32,
    clip_right: u32,
) -> (bool, bool) {
    if exons.is_empty() || seq.is_empty() {
        return (false, false);
    }

    let seq_len = seq.len();
    let left_clip = (clip_left as usize).min(seq_len);
    let right_clip = (clip_right as usize).min(seq_len);

    let first_exon_len = exons
        .first()
        .map_or(0usize, |(s, e)| e.saturating_sub(*s) as usize);
    let first_start = left_clip.min(seq_len);
    let first_end = first_start
        .saturating_add(first_exon_len)
        .min(seq_len.saturating_sub(right_clip));
    let mut t_count = 0usize;
    let mut first_count = 0usize;
    for &b in seq.get(first_start..first_end).unwrap_or(&[]) {
        if b == b'T' {
            t_count += 1;
        }
        first_count += 1;
    }
    let first_exon_polyt = if first_count > 0 {
        (t_count as f64) / (first_count as f64) >= POLYA_MIN_FRAC
    } else {
        false
    };

    let last_exon_len = exons
        .last()
        .map_or(0usize, |(s, e)| e.saturating_sub(*s) as usize);
    let last_end = seq_len.saturating_sub(right_clip);
    let last_start = last_end.saturating_sub(last_exon_len);
    let mut a_count = 0usize;
    let mut last_count = 0usize;
    for &b in seq.get(last_start..last_end).unwrap_or(&[]) {
        if b == b'A' {
            a_count += 1;
        }
        last_count += 1;
    }
    let last_exon_polya = if last_count > 0 {
        (a_count as f64) / (last_count as f64) >= POLYA_MIN_FRAC
    } else {
        false
    };

    (last_exon_polya, first_exon_polyt)
}

/// For boundary detection only: ref span for any mapped alignment (incl. secondary/supplementary).
/// Returns (ref_start_0based_incl, ref_end_0based_excl). None if unmapped or invalid.
/// original algorithm.cpp:516 — mapped_len >= 10. SAM positions are 1-based so we convert.
pub fn record_ref_span(record: &RecordBuf) -> Option<(u64, u64)> {
    if record.flags().is_unmapped() {
        return None;
    }
    let start_1b = record.alignment_start().map(|p| p.get() as u64).unwrap_or(0);
    let end_1b = start_1b
        .saturating_add(record.cigar().alignment_span() as u64)
        .saturating_sub(1);
    if end_1b <= start_1b || end_1b.saturating_sub(start_1b) < 9 {
        return None;
    }
    let ref_start_0 = start_1b.saturating_sub(1);
    let ref_end_0_excl = end_1b;
    Some((ref_start_0, ref_end_0_excl))
}

/// Ref span capped to the portion before any N (intron) op larger than `max_intron`.
/// Used for bundle-boundary construction to prevent a single read with a huge intron
/// from bridging two separate loci into one mega-bundle.
/// Returns (ref_start_0, ref_end_0_excl) where ref_end is truncated to exclude
/// anything after the first huge intron.
pub fn record_ref_span_capped(record: &RecordBuf, max_intron: u64) -> Option<(u64, u64)> {
    use noodles_sam::alignment::record::cigar::op::Kind;
    if record.flags().is_unmapped() {
        return None;
    }
    let start_1b = record.alignment_start().map(|p| p.get() as u64).unwrap_or(0);
    if start_1b == 0 {
        return None;
    }
    let ref_start_0 = start_1b.saturating_sub(1);
    let mut cur = ref_start_0;
    let cigar = record.cigar();
    for op in cigar.as_ref().iter() {
        let len = op.len() as u64;
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch | Kind::Deletion => {
                cur = cur.saturating_add(len);
            }
            Kind::Skip => {
                if len > max_intron {
                    // Truncate here: the remaining portion of the read is treated as
                    // a separate locus contribution, not merged across the huge gap.
                    if cur <= ref_start_0 + 9 {
                        return None;
                    }
                    return Some((ref_start_0, cur));
                }
                cur = cur.saturating_add(len);
            }
            _ => {}
        }
    }
    if cur <= ref_start_0 + 9 {
        return None;
    }
    Some((ref_start_0, cur))
}

/// Extract per-base mismatches by comparing read sequence to reference FASTA
/// at each aligned position. Does not require the MD tag. Walks the CIGAR once
/// and emits `(ref_pos_0based, query_base_upper)` for each position where
/// read_base != ref_base.
///
/// Cheap: O(aligned_length). Genome lookups are O(1) on the in-memory FASTA.
pub fn extract_mismatches_vs_fasta(
    record: &RecordBuf,
    ref_start: u64,
    chrom: &str,
    genome: &crate::genome::GenomeIndex,
) -> Vec<(u64, u8)> {
    use noodles_sam::alignment::record::cigar::op::Kind;

    let seq = record.sequence();
    if seq.is_empty() {
        return Vec::new();
    }
    let seq_bytes = seq.as_ref();
    let mut mismatches = Vec::new();
    let mut ref_pos = ref_start;
    let mut query_offset = 0usize;
    let cigar = record.cigar();
    for result in cigar.iter() {
        let Ok(op) = result else {
            continue;
        };
        let kind = op.kind();
        let len = op.len() as u64;
        match kind {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                for i in 0..len {
                    let rp = ref_pos + i;
                    let qo = query_offset + i as usize;
                    if qo >= seq_bytes.len() {
                        break;
                    }
                    let q = seq_bytes[qo].to_ascii_uppercase();
                    // `seq_bytes[qo]` is raw from noodles (already 1 base per u8).
                    // Only record A/C/G/T mismatches; ignore N and non-canonical bases.
                    if !matches!(q, b'A' | b'C' | b'G' | b'T') {
                        continue;
                    }
                    if let Some(seq_ref) = genome_ref_base(genome, chrom, rp) {
                        if seq_ref != q {
                            mismatches.push((rp, q));
                        }
                    }
                }
                ref_pos += len;
                query_offset += len as usize;
            }
            Kind::Deletion | Kind::Skip => {
                ref_pos += len;
            }
            Kind::Insertion | Kind::SoftClip => {
                query_offset += len as usize;
            }
            Kind::HardClip | Kind::Pad => {}
        }
    }
    mismatches
}

/// Look up a single reference base at a 0-based position. Uses the public
/// fetch_sequence with a 1-bp window (cheap on a fully-resident FASTA).
#[inline]
fn genome_ref_base(
    genome: &crate::genome::GenomeIndex,
    chrom: &str,
    pos0: u64,
) -> Option<u8> {
    genome
        .fetch_sequence(chrom, pos0, pos0 + 1)
        .and_then(|v| v.first().copied())
}

/// Convert one BAM record to BundleRead (exons, weight, poly). 0-based coordinates.
/// Returns None for unmapped; callers should also skip secondary/supplementary for readlist
pub fn record_to_bundle_read(record: &RecordBuf) -> Option<BundleRead> {
    record_to_bundle_read_with_snp(record, None, None)
}

/// Same as [`record_to_bundle_read`] but additionally extracts per-base
/// mismatches vs the reference FASTA when both `chrom` and `genome` are
/// provided. Used in VG mode to enable copy-distinguishing SNP scoring
/// without depending on the MD tag.
pub fn record_to_bundle_read_with_snp(
    record: &RecordBuf,
    chrom: Option<&str>,
    genome: Option<&crate::genome::GenomeIndex>,
) -> Option<BundleRead> {
    if record.flags().is_unmapped() {
        return None;
    }
    let start_1b = record.alignment_start().map(|p| p.get() as u64).unwrap_or(0);
    let end_1b = start_1b
        .saturating_add(record.cigar().alignment_span() as u64)
        .saturating_sub(1);
    // original algorithm.cpp:516 — reject invalid/very short mappings globally.
    if end_1b <= start_1b || end_1b.saturating_sub(start_1b) < 9 {
        return None;
    }
    let ref_start = start_1b.saturating_sub(1);
    let parsed = parse_cigar_compat(ref_start, record.cigar());
    let exons = parsed.exons;
    if exons.is_empty() {
        return None;
    }
    let ref_end = exons.last().map(|e| e.1).unwrap_or(ref_start);
    let nh = get_nh(record);
    // MD tag: defer String allocation — only needed downstream for --vg-snp.
    // Store None here; vg.rs re-extracts from BAM when needed.
    let md: Option<String> = None;
    let yc = get_yc(record);
    let yk = get_yk(record);
    // rdcount = YC; if(YK) rdcount=YK; if(isunitig) unitig=true
    let mut weight = yc;
    if yk > 0.0 {
        weight = yk;
    }
    // nomulti: if(!nomulti) rdcount/=nh; - we always divide by nh
    weight = weight / (nh as f64).max(1.0);
    let is_reverse = record.flags().is_reverse_complemented();
    let strand = inferred_read_strand(record);
    let query_length = Some(record.sequence().len() as u64);
    let clip_left = parsed.clip_left;
    let clip_right = parsed.clip_right;
    // Extract sequence bytes once and reuse for both poly-A checks.
    let seq_bytes = seq_bytes_from_record(record);
    let (
        has_poly_start_aligned,
        has_poly_start_unaligned,
        has_poly_end_aligned,
        has_poly_end_unaligned,
    ) = if seq_bytes.is_empty() {
        (false, false, false, false)
    } else {
        polya::detect_polya_aligned_unaligned(
            &seq_bytes, clip_left, clip_right,
            POLYA_MIN_CONSEC, POLYA_WINDOW, POLYA_MIN_FRAC,
        )
    };
    let (has_last_exon_polya, has_first_exon_polyt) =
        detect_terminal_exon_poly_flags_with_seq(&seq_bytes, &exons, clip_left, clip_right);
    // the original algorithm CReadAln stores integer unaligned poly tail evidence counters.
    // One BAM alignment contributes one unit of support when tail evidence is present.
    let unaligned_poly_t: u16 = if has_poly_start_unaligned { 1 } else { 0 };
    let unaligned_poly_a: u16 = if has_poly_end_unaligned { 1 } else { 0 };
    let nm = get_nm_adjusted(record, clip_left, clip_right);
    let has_poly_start = has_poly_start_aligned || has_poly_start_unaligned;
    let has_poly_end = has_poly_end_aligned || has_poly_end_unaligned;
    let insertion_sites = parsed.insertion_sites;
    let junctions_raw = junctions_from_exons(&exons);
    let junctions_del = deletion_aware_junctions_from_parts(&exons, &parsed.junc_del_offsets);

    // RUSTLE_JUNCSDEL_DUMP=/path/file.tsv — emit per-read juncsdel array for
    // byte-level parity verification against StringTie's PARITY_JUNCSDEL_TSV.
    // Columns: read_name, ref_start, ref_end, n_exons, n_juncs, n_juncsdel,
    //          cigar_exons (as start-end,...), juncsdel (as start|end,...).
    if let Ok(path) = std::env::var("RUSTLE_JUNCSDEL_DUMP") {
        if !path.is_empty() {
            use std::io::Write;
            use std::sync::Mutex;
            use std::sync::OnceLock;
            static WR: OnceLock<Mutex<Option<std::fs::File>>> = OnceLock::new();
            static HDR: OnceLock<Mutex<bool>> = OnceLock::new();
            let wr = WR.get_or_init(|| {
                Mutex::new(
                    std::fs::OpenOptions::new()
                        .create(true)
                        .append(true)
                        .open(&path)
                        .ok(),
                )
            });
            let hdr = HDR.get_or_init(|| Mutex::new(false));
            if let (Ok(mut f_opt), Ok(mut hdr_w)) = (wr.lock(), hdr.lock()) {
                if let Some(f) = f_opt.as_mut() {
                    if !*hdr_w {
                        let _ = writeln!(
                            f,
                            "source\tread_name\tref_start\tref_end\tn_exons\tn_juncs\tn_juncsdel\texons\tjuncsdel"
                        );
                        *hdr_w = true;
                    }
                    let name = record
                        .name()
                        .map(|n| n.to_string())
                        .unwrap_or_default();
                    let ex_str = exons
                        .iter()
                        .map(|(s, e)| format!("{}-{}", s, e))
                        .collect::<Vec<_>>()
                        .join(",");
                    let jd_str = parsed
                        .junc_del_offsets
                        .iter()
                        .map(|(l, r)| format!("{}|{}", l, r))
                        .collect::<Vec<_>>()
                        .join(",");
                    let _ = writeln!(
                        f,
                        "rustle\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        name,
                        ref_start,
                        ref_end,
                        exons.len(),
                        junctions_raw.len(),
                        parsed.junc_del_offsets.len(),
                        ex_str,
                        jd_str
                    );
                }
            }
        }
    }

    if std::env::var_os("RUSTLE_READ_DEBUG_COORDS").is_some()
        && exons.iter().any(|(s, e)| *s == 17439907 || *e == 17433882)
    {
        eprintln!(
            "[READ_DEBUG] name={} start={} end={} weight={:.2} raw={:?} del={:?}",
            record.name().map(|n| n.to_string()).unwrap_or_default(),
            ref_start,
            ref_end,
            weight,
            junctions_raw,
            junctions_del
        );
    }

    // Use exon-derived (raw) junctions for the active pipeline to match
    // StringTie's simpler junction model (rd.juncs[i]->start/end align with
    // rd.segs[i].end/rd.segs[i+1].start, no deletion offsets).
    // Opt-out with RUSTLE_ACTIVE_JUNCTIONS_DEL=1 to fall back to the prior
    // deletion-aware active stream if a regression needs isolating.
    let use_raw = std::env::var_os("RUSTLE_ACTIVE_JUNCTIONS_DEL").is_none();
    let junctions = if use_raw {
        junctions_raw.clone()
    } else {
        junctions_del.clone()
    };
    // When the active pipeline uses raw junctions, also align `junctions_del` so
    // deljuncmatch-based dedup groups reads by the same exon-derived signature
    // StringTie effectively uses (it has no separate deletion-aware signature).
    let junctions_del = if use_raw { junctions_raw.clone() } else { junctions_del };
    let junction_valid = vec![true; junctions.len()];
    let mate_start = record
        .mate_alignment_start()
        .map(|p| (p.get() as u64).saturating_sub(1));
    let exon_bp: u64 = exons.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
    let nh_f = (nh as f64).max(1.0);
    // countFragment:
    //   frag_len += sum(exon_len)/NH
    //   num_fragments += 1/NH for unpaired OR read1 OR (read2 && mate_unmapped)
    let flags = record.flags();
    let countfrag_len = (exon_bp as f64) / nh_f;
    let countfrag_num = if !flags.is_segmented()
        || flags.is_first_segment()
        || (flags.is_last_segment() && flags.is_mate_unmapped())
    {
        1.0 / nh_f
    } else {
        0.0
    };

    let read_name: Arc<str> = record
        .name()
        .map(|n| Arc::<str>::from(n.to_string()))
        .unwrap_or_default();
    let read_name_hash = fnv1a64(&read_name);

    // VG fields: HP/PS phasing tags.
    let hp_tag = get_tag_int(record, "HP").and_then(|v| if v > 0 { Some(v as u8) } else { None });
    let ps_tag = get_tag_int(record, "PS").and_then(|v| if v > 0 { Some(v as u32) } else { None });

    // VG SNP: compute per-base mismatches vs FASTA when genome + chrom are
    // supplied. This replaces the old MD-tag-based path so we can work with
    // BAMs that lack MD tags. Non-VG runs pass `None` and pay nothing.
    let mismatches: Vec<(u64, u8)> = match (chrom, genome) {
        (Some(c), Some(g)) => extract_mismatches_vs_fasta(record, ref_start, c, g),
        _ => Vec::new(),
    };

    let is_primary_alignment =
        !record.flags().is_secondary() && !record.flags().is_supplementary();

    Some(BundleRead {
        read_uid: NEXT_READ_UID.fetch_add(1, Ordering::Relaxed),
        read_name,
        read_name_hash,
        ref_id: record.reference_sequence_id(),
        mate_ref_id: record.mate_reference_sequence_id(),
        mate_start,
        hi: get_hi(record),
        ref_start,
        ref_end,
        exons,
        junctions,
        junction_valid,
        junctions_raw,
        junctions_del,
        weight,
        is_reverse,
        strand,
        has_poly_start,
        has_poly_end,
        has_poly_start_aligned,
        has_poly_start_unaligned,
        has_poly_end_aligned,
        has_poly_end_unaligned,
        unaligned_poly_t,
        unaligned_poly_a,
        has_last_exon_polya,
        has_first_exon_polyt,
        query_length,
        clip_left,
        clip_right,
        nh,
        nm,
        md,
        insertion_sites,
        unitig: yk > 0.0,
        unitig_cov: yk,
        read_count_yc: yc,
        countfrag_len,
        countfrag_num,
        junc_mismatch_weight: 0.0,
        pair_idx: Vec::new(),
        pair_count: Vec::new(),
        mapq: 0, // Not used for filtering (minimap2 assigns MAPQ=0 to all multi-mappers)
        mismatches,
        hp_tag,
        ps_tag,
        is_primary_alignment,
    })
}

/// Open BAM from path. Plain .bam (uncompressed or bgzf) via noodles_bam.
pub fn open_bam<P: AsRef<Path>>(
    path: P,
    num_threads: usize,
) -> Result<BamReader<noodles_bgzf::MultithreadedReader<io::BufReader<std::fs::File>>>> {
    let path = path.as_ref();
    let file = std::fs::File::open(path)?;
    let buf = io::BufReader::with_capacity(1 << 20, file); // 1MB buffer for BAM I/O
    let worker_count = std::num::NonZeroUsize::new(num_threads).unwrap_or(std::num::NonZeroUsize::MIN);
    let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(worker_count, buf);
    let reader = BamReader::from(bgzf);
    Ok(reader)
}

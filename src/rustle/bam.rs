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
struct CigarParityParse {
    exons: Vec<(u64, u64)>,
    /// Per-intron deletion offsets: (left_del_before_N, right_del_after_N).
    junc_del_offsets: Vec<(u64, u64)>,
    clip_left: u32,
    clip_right: u32,
    insertion_sites: Vec<u64>,
}

/// Parse CIGAR with the reference assembler `GSamRecord::setupCoordinates` parity semantics.
///
/// Coordinates are returned as 0-based half-open intervals.
fn parse_cigar_parity(ref_start: u64, cigar: &Cigar) -> CigarParityParse {
    let mut out = CigarParityParse::default();

    // Mirrors GSam.cpp state machine variables.
    let mut l = 0u64;
    let mut exstart = ref_start;
    let mut exon_started = false;
    let mut intron = false;
    let mut ins = false;
    let mut del = 0u64;
    let mut prevdel = 0u64;
    // Deletion(s) immediately after an intron (N) before the first matched base of the next exon.
    // We accumulate these and apply them once to the acceptor coordinate for that intron.
    let mut post_n_del = 0u64;

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
                    // Close the intron boundary: one offset entry per N.
                    out.junc_del_offsets.push((prevdel, post_n_del));
                    post_n_del = 0;
                }
                intron = false;
                ins = false;
                del = 0;
            }
            Kind::Deletion => {
                del = len;
                l = l.saturating_add(del);
                if intron {
                    // Deletions between N and the first match belong to the acceptor boundary.
                    post_n_del = post_n_del.saturating_add(del);
                }
                ins = false;
            }
            Kind::Insertion => {
                // Same reference position tracking as C++ mismatch_anchor CIGAR pass.
                out.insertion_sites.push(ref_start.saturating_add(l));
                ins = true;
            }
            Kind::Skip => {
                // C++ guard: anomalous leading intron before any exon.
                if !exon_started {
                    // GSamRecord::setupCoordinates() only skips this CIGAR op;
                    // it does not abort parsing the rest of the alignment.
                    continue;
                }
                // C++: if (!ins || !intron) close previous exon.
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
                post_n_del = 0;
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
    Ok(parse_cigar_parity(ref_start, cigar).exons)
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

/// NH tag (alignment hit count) for multi-mapper weight: weight = 1/NH (C++ reference, 899; reference assembler 1/NH).
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

/// HI tag (alignment hit index), used by the reference assembler to disambiguate mate matching keys.
fn get_hi(record: &RecordBuf) -> u32 {
    get_tag_int(record, "HI").unwrap_or(0)
}

/// NM/nM edit distance with C++ parity adjustments.
///
/// C++ behavior (C++ reference):
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

/// MD mismatch/deletion string when present.
fn get_md(record: &RecordBuf) -> Option<String> {
    let tag = parse_tag("MD")?;
    let data = record.data();
    data.get(&tag)
        .and_then(|v| match v {
            noodles_sam::alignment::record_buf::data::field::Value::String(s)
            | noodles_sam::alignment::record_buf::data::field::Value::Hex(s) => {
                Some(s.to_string())
            }
            _ => None,
        })
}

/// Parse MD tag into per-base mismatches: Vec<(ref_position, query_base)>.
/// The MD tag format: "10A5^AC3T2" means 10 matches, mismatch at ref A, 5 matches,
/// 2bp deletion (^AC), 3 matches, mismatch at ref T, 2 matches.
/// We extract the mismatch positions (not deletions) and look up the query base
/// from the BAM record's sequence.
fn parse_md_mismatches(md: Option<&str>, ref_start: u64, record: &RecordBuf) -> Vec<(u64, u8)> {
    let Some(md_str) = md else {
        return Vec::new();
    };
    let seq = record.sequence();
    if seq.is_empty() {
        return Vec::new();
    }
    let bytes = md_str.as_bytes();
    let mut mismatches = Vec::new();
    let mut ref_pos = ref_start;
    let mut query_offset = 0usize; // Offset into the query sequence.
    let mut p = 0usize;

    while p < bytes.len() {
        let c = bytes[p];
        if c.is_ascii_digit() {
            // Parse number of matching bases.
            let mut n = 0u64;
            while p < bytes.len() && bytes[p].is_ascii_digit() {
                n = n * 10 + (bytes[p] - b'0') as u64;
                p += 1;
            }
            ref_pos += n;
            query_offset += n as usize;
        } else if c == b'^' {
            // Deletion from reference: skip ref bases, don't advance query.
            p += 1;
            while p < bytes.len() && bytes[p].is_ascii_alphabetic() {
                ref_pos += 1;
                p += 1;
            }
        } else if c.is_ascii_alphabetic() {
            // Mismatch: `c` is the reference base, query base is at query_offset.
            if query_offset < seq.len() {
                let query_base = seq.as_ref()[query_offset];
                mismatches.push((ref_pos, query_base));
            }
            ref_pos += 1;
            query_offset += 1;
            p += 1;
        } else {
            p += 1; // Skip unknown characters.
        }
    }

    mismatches
}

/// YC tag (alignment count): the reference assembler rdcount uses YC when present, then scales by NH.
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

/// YK tag: unitig coverage (C++ reference unitig_cov=brec.tag_float("YK")).
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
    match ts {
        '+' => Some(!is_rev),
        '-' => Some(is_rev),
        _ => None,
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
        // C++ can keep unknown strand (0 / '.') and infer it later from
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

/// Detect polyA/T from record: returns (aligned_start, unaligned_start, aligned_end, unaligned_end).
/// aligned_polyT/aligned_polyA = in-alignment body; unaligned_* = soft-clip tail.
fn detect_polya_from_record(
    record: &RecordBuf,
    clip_left: u32,
    clip_right: u32,
) -> (bool, bool, bool, bool) {
    let seq_bytes = seq_bytes_from_record(record);
    if seq_bytes.is_empty() {
        return (false, false, false, false);
    }
    polya::detect_polya_aligned_unaligned(
        &seq_bytes,
        clip_left,
        clip_right,
        POLYA_MIN_CONSEC,
        POLYA_WINDOW,
        POLYA_MIN_FRAC,
    )
}

/// C++ reference check_last_exon_polyA/check_first_exon_polyT parity:
/// ratio of A on last exon, ratio of T on first exon (>=0.8 triggers trim candidate).
/// Exons are 0-based half-open, so exon length is end-start.
fn detect_terminal_exon_poly_flags(
    record: &RecordBuf,
    exons: &[(u64, u64)],
    clip_left: u32,
    clip_right: u32,
) -> (bool, bool) {
    if exons.is_empty() {
        return (false, false);
    }
    let seq = seq_bytes_from_record(record);
    if seq.is_empty() {
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
/// reference assembler.cpp:516 — mapped_len >= 10. SAM positions are 1-based so we convert.
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

/// Convert one BAM record to BundleRead (exons, weight, poly). 0-based coordinates.
/// Returns None for unmapped; callers should also skip secondary/supplementary for readlist (C++ reference).
pub fn record_to_bundle_read(record: &RecordBuf) -> Option<BundleRead> {
    if record.flags().is_unmapped() {
        return None;
    }
    let start_1b = record.alignment_start().map(|p| p.get() as u64).unwrap_or(0);
    let end_1b = start_1b
        .saturating_add(record.cigar().alignment_span() as u64)
        .saturating_sub(1);
    // reference assembler.cpp:516 — reject invalid/very short mappings globally.
    if end_1b <= start_1b || end_1b.saturating_sub(start_1b) < 9 {
        return None;
    }
    let ref_start = start_1b.saturating_sub(1);
    let parsed = parse_cigar_parity(ref_start, record.cigar());
    let exons = parsed.exons;
    if exons.is_empty() {
        return None;
    }
    let ref_end = exons.last().map(|e| e.1).unwrap_or(ref_start);
    let nh = get_nh(record);
    let md = get_md(record);
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
    let (
        has_poly_start_aligned,
        has_poly_start_unaligned,
        has_poly_end_aligned,
        has_poly_end_unaligned,
    ) = detect_polya_from_record(record, clip_left, clip_right);
    let (has_last_exon_polya, has_first_exon_polyt) =
        detect_terminal_exon_poly_flags(record, &exons, clip_left, clip_right);
    // the reference assembler CReadAln stores integer unaligned poly tail evidence counters.
    // One BAM alignment contributes one unit of support when tail evidence is present.
    let unaligned_poly_t: u16 = if has_poly_start_unaligned { 1 } else { 0 };
    let unaligned_poly_a: u16 = if has_poly_end_unaligned { 1 } else { 0 };
    let nm = get_nm_adjusted(record, clip_left, clip_right);
    let has_poly_start = has_poly_start_aligned || has_poly_start_unaligned;
    let has_poly_end = has_poly_end_aligned || has_poly_end_unaligned;
    let insertion_sites = parsed.insertion_sites;
    let junctions_raw = junctions_from_exons(&exons);
    let junctions_del = deletion_aware_junctions_from_parts(&exons, &parsed.junc_del_offsets);

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

    // Keep active junction stream deletion-aware (existing pipeline behavior),
    // while retaining raw exon-boundary junctions for parity audits.
    let junctions = junctions_del.clone();
    let junction_valid = vec![true; junctions.len()];
    let mate_start = record
        .mate_alignment_start()
        .map(|p| (p.get() as u64).saturating_sub(1));
    let exon_bp: u64 = exons.iter().map(|(s, e)| e.saturating_sub(*s)).sum();
    let nh_f = (nh as f64).max(1.0);
    // C++ reference countFragment parity:
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

    // VG SNP: extract per-base mismatches from MD tag.
    // Only populated here — downstream code (vg.rs) uses it when --vg-snp is active.
    let mismatches = parse_md_mismatches(md.as_deref(), ref_start, record);

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
    })
}

/// Open BAM from path. Plain .bam (uncompressed or bgzf) via noodles_bam.
pub fn open_bam<P: AsRef<Path>>(
    path: P,
    num_threads: usize,
) -> Result<BamReader<noodles_bgzf::MultithreadedReader<io::BufReader<std::fs::File>>>> {
    let path = path.as_ref();
    let file = std::fs::File::open(path)?;
    let buf = io::BufReader::new(file);
    let worker_count = std::num::NonZeroUsize::new(num_threads).unwrap_or(std::num::NonZeroUsize::MIN);
    let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(worker_count, buf);
    let reader = BamReader::from(bgzf);
    Ok(reader)
}

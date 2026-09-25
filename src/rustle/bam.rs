//! BAM reading: open a BAM, turn an alignment's CIGAR into exon blocks, and read the aux tags the
//! pipeline uses (`AS:i`, `de:f`, `ts:A`) from either record representation.
//!
//! This is all that is left of the retired StringTie read-ingestion layer (bundle reads, poly-A
//! detection, NH/HI/YC/YK tags, SAM/CRAM transcoding — removed 2026-09-24, recover from tag
//! `notebook-2026-09-24`). The exon walk below is kept exactly as it was, because every assembled
//! intron chain goes through it.

use anyhow::Result;
use noodles_bam::io::Reader as BamReader;
use noodles_sam::alignment::record::cigar::op::Kind;
use noodles_sam::alignment::record::data::field::{Tag, Value};
use noodles_sam::alignment::record::Cigar as _;
use noodles_sam::alignment::record_buf::Cigar;
use std::io;
use std::path::Path;

/// Exon blocks of a CIGAR, following StringTie's `GSamRecord::setupCoordinates`
/// (gclib/GSam.cpp:197-291): `M`/`=`/`X`/`D` extend the current block, `N` closes it (unless an
/// insertion sits directly after the intron, `GSam.cpp:247`), a leading `N` before any aligned base is
/// skipped, and the final block is closed only when the read does not end in an intron.
///
/// Coordinates are 0-based half-open.
fn cigar_exons(ref_start: u64, cigar: &Cigar) -> Vec<(u64, u64)> {
    let mut exons = Vec::new();
    let mut l = 0u64;
    let mut exstart = ref_start;
    let mut exon_started = false;
    let mut intron = false;
    let mut ins = false;

    for result in cigar.iter() {
        let Ok(op) = result else {
            continue;
        };
        let len = op.len() as u64;

        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                exon_started = true;
                l = l.saturating_add(len);
                intron = false;
                ins = false;
            }
            Kind::Deletion => {
                l = l.saturating_add(len);
                ins = false;
            }
            Kind::Insertion => {
                ins = true;
            }
            Kind::Skip => {
                // guard: anomalous leading intron before any exon — skip the op, keep parsing.
                if !exon_started {
                    continue;
                }
                // GSam.cpp:247 — `if(!ins || !intron)` closes the preceding exon.
                if !ins || !intron {
                    let exon_end = ref_start.saturating_add(l);
                    if exon_end > exstart {
                        exons.push((exstart, exon_end));
                    }
                }
                l = l.saturating_add(len);
                exstart = ref_start.saturating_add(l);
                intron = true;
            }
            Kind::SoftClip | Kind::HardClip => {
                ins = false;
            }
            Kind::Pad => {}
        }
    }

    // GSam.cpp:283-287 — close the final exon only if the read did not end in an intron.
    if !intron {
        let exon_end = ref_start.saturating_add(l);
        if exon_end > exstart {
            exons.push((exstart, exon_end));
        }
    }

    exons
}

/// Exons = alignment blocks (N = intron splits). 0-based [start, end).
pub fn exons_from_cigar(ref_start: u64, cigar: &Cigar) -> Result<Vec<(u64, u64)>> {
    Ok(cigar_exons(ref_start, cigar))
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

// ---------------------------------------------------------------------------------------------------------
// Aux tags. One implementation for `RecordBuf` and for lazy `noodles_bam::Record`s: both implement
// `noodles_sam::alignment::Record`, whose `data().iter()` yields the fields in stored order. For a `RecordBuf`
// that is the order `Data::get` searches (first match), so these readers return exactly what `get` did; for
// a lazy record they return what `noodles_bam::record::Data::get` does, including its "a malformed field
// before the tag is an error" rule.
// ---------------------------------------------------------------------------------------------------------

/// minimap2's gap-compressed per-base divergence, `de:f`.
pub const DE_TAG: Tag = Tag::new(b'd', b'e');
/// Transcript strand, `ts:A` (`+`/`-`).
pub const TS_TAG: Tag = Tag::new(b't', b's');

/// `AS:i` at any integer width; `None` for a non-integer value.
fn value_as_i32(value: Value<'_>) -> Option<i32> {
    match value {
        Value::Int8(v) => Some(v as i32),
        Value::UInt8(v) => Some(v as i32),
        Value::Int16(v) => Some(v as i32),
        Value::UInt16(v) => Some(v as i32),
        Value::Int32(v) => Some(v),
        Value::UInt32(v) => Some(v as i32),
        _ => None,
    }
}

fn value_as_f32(value: Value<'_>) -> Option<f32> {
    match value {
        Value::Float(v) => Some(v),
        _ => None,
    }
}

fn value_as_char(value: Value<'_>) -> Option<char> {
    match value {
        Value::Character(c) => Some(c as char),
        _ => None,
    }
}

/// The FIRST `tag` field of `record`, converted by `conv` (`Ok(None)` if absent or of another type). Fields
/// are decoded in order and the scan stops at the first match, so `Err` means a field BEFORE it is malformed.
pub fn aux_first<R, T>(record: &R, tag: Tag, conv: fn(Value<'_>) -> Option<T>) -> io::Result<Option<T>>
where
    R: noodles_sam::alignment::Record + ?Sized,
{
    for entry in record.data().iter() {
        let (t, value) = entry?;
        if t == tag {
            return Ok(conv(value));
        }
    }
    Ok(None)
}

/// `AS:i` alignment score (the gate signal). `None` if absent, non-integer, or behind a malformed field.
pub fn record_as<R: noodles_sam::alignment::Record + ?Sized>(record: &R) -> Option<i32> {
    aux_first(record, Tag::ALIGNMENT_SCORE, value_as_i32).ok().flatten()
}

/// `de:f` divergence (the conflict-criterion signal). `None` if absent, non-float, or behind a malformed field.
pub fn record_de<R: noodles_sam::alignment::Record + ?Sized>(record: &R) -> Option<f32> {
    aux_first(record, DE_TAG, value_as_f32).ok().flatten()
}

/// `ts:A` transcript strand. `None` if absent, not a character, or behind a malformed field.
pub fn record_ts<R: noodles_sam::alignment::Record + ?Sized>(record: &R) -> Option<char> {
    aux_first(record, TS_TAG, value_as_char).ok().flatten()
}

/// `AS:i`, `de:f` and `ts:A` in ONE pass over EVERY field: the first occurrence of each tag wins (a later
/// duplicate is ignored even when the first could not be converted), and a malformed field ANYWHERE is an
/// error, as it is when the record is decoded into a `RecordBuf`.
pub fn record_as_de_ts<R>(record: &R) -> io::Result<(Option<i32>, Option<f32>, Option<char>)>
where
    R: noodles_sam::alignment::Record + ?Sized,
{
    let (mut as_score, mut de, mut ts) = (None, None, None);
    let (mut seen_as, mut seen_de, mut seen_ts) = (false, false, false);
    for entry in record.data().iter() {
        let (tag, value) = entry?;
        if tag == Tag::ALIGNMENT_SCORE && !seen_as {
            seen_as = true;
            as_score = value_as_i32(value);
        } else if tag == DE_TAG && !seen_de {
            seen_de = true;
            de = value_as_f32(value);
        } else if tag == TS_TAG && !seen_ts {
            seen_ts = true;
            ts = value_as_char(value);
        }
    }
    Ok((as_score, de, ts))
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_sam::alignment::record::cigar::op::{Kind, Op};

    fn cig(ops: &[(Kind, usize)]) -> Cigar {
        ops.iter().map(|&(k, n)| Op::new(k, n)).collect()
    }

    #[test]
    fn exons_split_on_introns_and_keep_deletions_inside_blocks() {
        let c = cig(&[(Kind::SoftClip, 5), (Kind::Match, 10), (Kind::Deletion, 2), (Kind::Match, 3),
                      (Kind::Skip, 100), (Kind::Match, 20), (Kind::HardClip, 4)]);
        assert_eq!(exons_from_cigar(1000, &c).unwrap(), vec![(1000, 1015), (1115, 1135)]);
    }

    #[test]
    fn leading_intron_is_skipped_and_trailing_intron_leaves_the_last_block_open() {
        let c = cig(&[(Kind::Skip, 50), (Kind::Match, 10), (Kind::Skip, 30)]);
        assert_eq!(exons_from_cigar(0, &c).unwrap(), vec![(0, 10)]);
    }

    #[test]
    fn insertion_right_after_an_intron_does_not_close_a_block() {
        // N I N: the second N sees `ins && intron` and does not push the (empty) block between them.
        let c = cig(&[(Kind::Match, 10), (Kind::Skip, 5), (Kind::Insertion, 2), (Kind::Skip, 5), (Kind::Match, 10)]);
        assert_eq!(exons_from_cigar(0, &c).unwrap(), vec![(0, 10), (20, 30)]);
    }
}

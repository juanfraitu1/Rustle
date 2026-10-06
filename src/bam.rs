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
    let worker_count =
        std::num::NonZeroUsize::new(num_threads).unwrap_or(std::num::NonZeroUsize::MIN);
    let bgzf = noodles_bgzf::MultithreadedReader::with_worker_count(worker_count, buf);
    let reader = BamReader::from(bgzf);
    Ok(reader)
}

/// The reader [`open_bam`] returns.
pub type BamFile = BamReader<noodles_bgzf::MultithreadedReader<io::BufReader<std::fs::File>>>;

/// One piece of a contig read on its own: records of `contig` whose alignment START (0-based) lies in `[lo, hi)`.
/// `lo = 0, hi = u64::MAX` is the whole contig. Used by the catalog's `--piecewise` mode.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ContigSpan {
    pub contig: String,
    pub lo: u64,
    pub hi: u64,
}

impl ContigSpan {
    pub fn whole(contig: &str) -> ContigSpan {
        ContigSpan {
            contig: contig.to_string(),
            lo: 0,
            hi: u64::MAX,
        }
    }
    pub fn is_whole(&self) -> bool {
        self.lo == 0 && self.hi == u64::MAX
    }
    /// `chr1` for a whole contig, `chr1:lo-hi` (0-based half-open, `end` for the open end) for a sub-range.
    pub fn label(&self) -> String {
        if self.is_whole() {
            self.contig.clone()
        } else if self.hi == u64::MAX {
            format!("{}:{}-end", self.contig, self.lo)
        } else {
            format!("{}:{}-{}", self.contig, self.lo, self.hi)
        }
    }
    /// Where a record of this contig falls relative to the span, from its 0-based alignment start (`None` = no
    /// position): `Less` = before `lo` (skip it), `Equal` = inside, `Greater` = at or past `hi` (stop: the BAM is
    /// sorted by start). A whole-contig span keeps every record, exactly as a whole-file scan does.
    pub fn place(&self, start0: Option<u64>) -> std::cmp::Ordering {
        use std::cmp::Ordering::*;
        if self.is_whole() {
            return Equal;
        }
        match start0 {
            None => Less, // unplaced within a placed contig cannot happen in a sorted BAM; never read as inside
            Some(s) if s < self.lo => Less,
            Some(s) if s >= self.hi => Greater,
            Some(_) => Equal,
        }
    }
}

/// Position a reader from [`open_bam`] (header already read) at the FIRST record of `contig`, through the `.bai`
/// next to the BAM. Returns the contig's reference-sequence id, or `None` when the index places no record on it.
///
/// The caller then reads records exactly as a whole-file scan does and STOPS at the first record whose reference
/// id is not this one. In a coordinate-sorted BAM (the only kind a `.bai` exists for) a contig's records are
/// contiguous, so this yields exactly the records a whole-file scan meets on `contig`, in the same order, with no
/// region-overlap filter in between (a region query drops a record whose alignment end it cannot compute; a
/// whole-file scan does not). The start is the index's per-reference metadata (samtools' pseudo-bin: the offset of
/// the contig's first record, placed-unmapped records included) when present, else the smallest chunk start of a
/// whole-contig query; the reference-id check makes either exact.
pub fn seek_to_contig(
    reader: &mut BamFile,
    header: &noodles_sam::Header,
    bam_path: &str,
    contig: &str,
) -> Result<Option<usize>> {
    seek_to_span(reader, header, bam_path, &ContigSpan::whole(contig))
}

/// [`seek_to_contig`] for a [`ContigSpan`]. A span starting past the contig's first base seeks to the smallest
/// chunk start of the index query for `[lo+1, ...]` (1-based), which precedes every record overlapping it; the
/// caller skips the records that start before `lo` ([`ContigSpan::place`]).
pub fn seek_to_span(
    reader: &mut BamFile,
    header: &noodles_sam::Header,
    bam_path: &str,
    span: &ContigSpan,
) -> Result<Option<usize>> {
    use noodles_bgzf::io::Seek as _;
    use noodles_csi::binning_index::{BinningIndex as _, ReferenceSequence as _};
    let contig = span.contig.as_str();
    let id = header
        .reference_sequences()
        .get_index_of(contig.as_bytes())
        .ok_or_else(|| anyhow::anyhow!("contig {contig} is not in the header of {bam_path}"))?;
    let bai = format!("{bam_path}.bai");
    let index = noodles_bam::bai::read(&bai).map_err(|e| anyhow::anyhow!("reading {bai}: {e}"))?;
    let meta = index
        .reference_sequences()
        .get(id)
        .and_then(|r| r.metadata())
        .map(|m| {
            (
                m.start_position(),
                m.mapped_record_count() + m.unmapped_record_count(),
            )
        });
    if let Some((_, 0)) = meta {
        return Ok(None);
    }
    let from_query = |lo1: usize| -> Result<Option<noodles_bgzf::VirtualPosition>> {
        let start = noodles_core::Position::try_from(lo1.max(1))?;
        Ok(index
            .query(id, noodles_core::region::Interval::from(start..))?
            .iter()
            .map(|c| c.start())
            .min())
    };
    let start = if span.lo == 0 {
        match meta {
            Some((pos, _)) => Some(pos),
            None => from_query(1)?,
        }
    } else {
        from_query(usize::try_from(span.lo + 1)?)?
    };
    let Some(start) = start else { return Ok(None) };
    reader.get_mut().seek_to_virtual_position(start)?;
    Ok(Some(id))
}

/// Cut `contig` into consecutive spans of at most about `max_records` records each, cutting ONLY at positions no
/// mapped record crosses (a record crosses `x` when its alignment starts before `x` and ends after it; every flag counts,
/// secondary and supplementary included, since the catalog reads them all). Every catalog step groups records that
/// overlap or share a junction, and neither can happen across such a cut, so the spans are independent: see
/// `denovo_pipeline::build_catalog_reps`. Each cut is made at the LAST read-free position before the span would
/// exceed `max_records`, so a span exceeds it only where no read-free position lies within reach (long genes and
/// giant introns link whole regions); that yields one larger span rather than an inexact cut. Returns
/// `(span, mapped records in it)`, covering the whole contig.
pub fn read_free_cuts(
    bam_path: &str,
    contig: &str,
    max_records: u64,
) -> Result<Vec<(ContigSpan, u64)>> {
    let mut reader = open_bam(bam_path, 2)?;
    let header = reader.read_header()?;
    let Some(id) = seek_to_contig(&mut reader, &header, bam_path, contig)? else {
        return Ok(vec![(ContigSpan::whole(contig), 0)]);
    };
    let mut out: Vec<(ContigSpan, u64)> = Vec::new();
    let (mut lo, mut n, mut max_end) = (0u64, 0u64, 0u64);
    // the latest read-free position of the current span, with the records before it
    let mut cand: Option<(u64, u64)> = None;
    let mut record = noodles_bam::Record::default();
    while reader.read_record(&mut record)? > 0 {
        match record.reference_sequence_id().transpose()? {
            Some(r) if r == id => {}
            Some(r) if r < id => continue,
            _ => break,
        }
        if record.flags().is_unmapped() {
            continue;
        }
        let Some(start) = record.alignment_start().transpose()? else {
            continue;
        };
        let s = (usize::from(start) as u64).saturating_sub(1);
        // `s >= max_end`: every earlier record ends at or before `s`, every later one starts at or after it
        if s >= max_end && s > lo && n > 0 {
            if n >= max_records {
                out.push((
                    ContigSpan {
                        contig: contig.to_string(),
                        lo,
                        hi: s,
                    },
                    n,
                ));
                lo = s;
                n = 0;
                cand = None;
            } else {
                cand = Some((s, n));
            }
        }
        let mut span = 0u64;
        for op in record.cigar().iter() {
            let op = op?;
            use noodles_sam::alignment::record::cigar::op::Kind::*;
            if matches!(
                op.kind(),
                Match | Deletion | Skip | SequenceMatch | SequenceMismatch
            ) {
                span += op.len() as u64;
            }
        }
        max_end = max_end.max(s + span);
        n += 1;
        if n > max_records {
            if let Some((x, m)) = cand.take() {
                out.push((
                    ContigSpan {
                        contig: contig.to_string(),
                        lo,
                        hi: x,
                    },
                    m,
                ));
                lo = x;
                n -= m;
            }
        }
    }
    out.push((
        ContigSpan {
            contig: contig.to_string(),
            lo,
            hi: u64::MAX,
        },
        n,
    ));
    if out.len() == 1 {
        out[0].0 = ContigSpan::whole(contig);
    }
    Ok(out)
}

/// Per-contig record counts from the `.bai` metadata (mapped + placed-unmapped), in HEADER order; `None` for a
/// contig whose index entry carries no metadata. Cheap: reads the header and the index, no record.
pub fn contig_record_counts(bam_path: &str) -> Result<Vec<(String, Option<u64>)>> {
    use noodles_csi::binning_index::ReferenceSequence as _;
    let mut reader = open_bam(bam_path, 1)?;
    let header = reader.read_header()?;
    let bai = format!("{bam_path}.bai");
    let index = noodles_bam::bai::read(&bai).map_err(|e| anyhow::anyhow!("reading {bai}: {e}"))?;
    Ok(header
        .reference_sequences()
        .keys()
        .enumerate()
        .map(|(i, name)| {
            // no metadata AND no bin = nothing indexed on this contig (samtools writes neither for an empty
            // reference); no metadata but bins = unknown count
            let n = index
                .reference_sequences()
                .get(i)
                .map_or(Some(0), |r| match r.metadata() {
                    Some(m) => Some(m.mapped_record_count() + m.unmapped_record_count()),
                    None if r.bins().is_empty() => Some(0),
                    None => None,
                });
            (format!("{name}"), n)
        })
        .collect())
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
pub fn aux_first<R, T>(
    record: &R,
    tag: Tag,
    conv: fn(Value<'_>) -> Option<T>,
) -> io::Result<Option<T>>
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
    aux_first(record, Tag::ALIGNMENT_SCORE, value_as_i32)
        .ok()
        .flatten()
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
        let c = cig(&[
            (Kind::SoftClip, 5),
            (Kind::Match, 10),
            (Kind::Deletion, 2),
            (Kind::Match, 3),
            (Kind::Skip, 100),
            (Kind::Match, 20),
            (Kind::HardClip, 4),
        ]);
        assert_eq!(
            exons_from_cigar(1000, &c).unwrap(),
            vec![(1000, 1015), (1115, 1135)]
        );
    }

    #[test]
    fn leading_intron_is_skipped_and_trailing_intron_leaves_the_last_block_open() {
        let c = cig(&[(Kind::Skip, 50), (Kind::Match, 10), (Kind::Skip, 30)]);
        assert_eq!(exons_from_cigar(0, &c).unwrap(), vec![(0, 10)]);
    }

    #[test]
    fn insertion_right_after_an_intron_does_not_close_a_block() {
        // N I N: the second N sees `ins && intron` and does not push the (empty) block between them.
        let c = cig(&[
            (Kind::Match, 10),
            (Kind::Skip, 5),
            (Kind::Insertion, 2),
            (Kind::Skip, 5),
            (Kind::Match, 10),
        ]);
        assert_eq!(exons_from_cigar(0, &c).unwrap(), vec![(0, 10), (20, 30)]);
    }
}

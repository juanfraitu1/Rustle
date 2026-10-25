//! The BAM side (o3port.py `read_primaries`, `oriented`, `pile_records`, `walk`, `observations`): records carry what the Python
//! reads off a pysam record; piles keep pysam's fetch order (Python dict insertion order: a better record replaces a read's
//! record in place).

use super::{Copy, Obs, ObsX, ReadId, AS_TIE};
use anyhow::{Context, Result};
use indexmap::IndexMap;
use std::collections::{BTreeSet, HashMap};
use std::sync::Arc;

// pysam CIGAR op codes
pub const OP_M: u8 = 0;
pub const OP_I: u8 = 1;
pub const OP_D: u8 = 2;
pub const OP_N: u8 = 3;
pub const OP_S: u8 = 4;
pub const OP_H: u8 = 5;
pub const OP_EQ: u8 = 7;
pub const OP_X: u8 = 8;

/// One alignment record of a read of the run (the fields the Python uses).
#[derive(Clone, Debug, PartialEq)]
pub struct Rec {
    pub read: ReadId,
    pub chrom: Arc<str>,
    pub start: i64,
    pub end: i64,
    pub as_: Option<i32>,
    pub nm: Option<i32>,
    pub reverse: bool,
    pub secondary: bool,
    pub supplementary: bool,
    pub cigar: Vec<(u8, u32)>,
    /// `query_sequence` ("" when the record has none, as a secondary often does).
    pub seq: Vec<u8>,
}

/// A read's primary record: Python's `(seq, quals, is_reverse, best AS)`, with both orientations precomputed.
#[derive(Clone, Debug)]
pub struct Prim {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub reverse: bool,
    pub best_as: i32,
    /// The primary record's (AS, intron count): the reference of the `--tie-structure` window (register 1404).
    pub primary: Option<(i32, usize)>,
    seq_rc: Vec<u8>,
    qual_rev: Vec<u8>,
}

impl Prim {
    /// Missing base qualities count as 40 each (the Python's `[40] * len(seq)`).
    pub fn new(seq: Vec<u8>, qual: Option<Vec<u8>>, reverse: bool, best_as: i32) -> Prim {
        let qual = qual.unwrap_or_else(|| vec![40; seq.len()]);
        let seq_rc = seq.iter().rev().map(|&b| comp(b)).collect();
        let qual_rev = qual.iter().rev().copied().collect();
        Prim { seq, qual, reverse, best_as, primary: None, seq_rc, qual_rev }
    }

    /// O2's tie window for one record of the read (AS, intron count): within AS_TIE of the read's best AS, or with `structure`
    /// (`--tie-structure`, register 1404) the primary's intron count and AS within AS_TIE of the primary's. Splice-mode AS charges each
    /// intron, so by best AS an unspliced processed-copy alignment can tie or outscore the read's own spliced one.
    pub fn near(&self, as_: i32, introns: usize, structure: bool) -> bool {
        if structure {
            self.primary.is_some_and(|(pa, pn)| introns == pn && (as_ - pa).abs() <= AS_TIE)
        } else {
            as_ >= self.best_as - AS_TIE
        }
    }
}

/// The record's intron count (CIGAR N operations).
pub fn introns(rec: &Rec) -> usize {
    rec.cigar.iter().filter(|&&(op, _)| op == OP_N).count()
}

/// `(chrom, start, end, AS, is secondary)` of every non-supplementary record with an AS tag.
pub type RecSummary = (Arc<str>, i64, i64, i32, bool, usize);
/// A pile: {read: its best record on the copy}, in fetch order.
pub type Pile = IndexMap<ReadId, Rec>;
/// {copy tid: pile}, in `--copies` order.
pub type Piles = IndexMap<String, Pile>;
/// `observations`' second product: {read: [(ref pos, base or b'-', inserted bases)]}.
pub type Steps = IndexMap<ReadId, Vec<(i64, u8, Vec<u8>)>>;

/// The read's bases and qualities in this record's orientation, and the offset of its first CIGAR query base (a leading hard clip).
pub fn oriented<'a>(p: &'a Prim, r: &Rec) -> (&'a [u8], &'a [u8], usize) {
    let lead = match r.cigar.first() {
        Some(&(OP_H, n)) => n as usize,
        _ => 0,
    };
    if r.reverse != p.reverse {
        (&p.seq_rc, &p.qual_rev, lead)
    } else {
        (&p.seq, &p.qual, lead)
    }
}

/// Python's `str.maketrans("ACGTNacgtn", "TGCANtgcan")`: other bytes unchanged.
pub fn comp(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        _ => b,
    }
}

/// Python's `CMP.get(b, "N")` (uppercase ACGT only).
pub fn cmp_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// `walk`'s items: M/=/X bases (ref pos, query index), D positions, I runs (the previous ref pos, query index, length).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Step {
    M(i64, usize),
    I(i64, usize, usize),
    D(i64),
}

pub fn walk(r: &Rec) -> Vec<Step> {
    let (mut qp, mut rp, mut out) = (0usize, r.start, Vec::new());
    for &(op, n) in &r.cigar {
        let (n, ni) = (n as usize, n as i64);
        match op {
            OP_M | OP_EQ | OP_X => {
                out.extend((0..n).map(|k| Step::M(rp + k as i64, qp + k)));
                qp += n;
                rp += ni;
            }
            OP_I => {
                out.push(Step::I(rp - 1, qp, n));
                qp += n;
            }
            OP_D => {
                out.extend((0..ni).map(|k| Step::D(rp + k)));
                rp += ni;
            }
            OP_N => rp += ni,
            OP_S => qp += n,
            _ => {}
        }
    }
    out
}

/// pysam `get_aligned_pairs(matches_only=True)`: (query index, ref pos) of every M/=/X base.
pub fn aligned_pairs(r: &Rec) -> impl Iterator<Item = (usize, i64)> + '_ {
    walk(r).into_iter().filter_map(|s| match s {
        Step::M(t, q) => Some((q, t)),
        _ => None,
    })
}

/// pysam `get_blocks()`: the aligned M/=/X blocks.
pub fn blocks(r: &Rec) -> Vec<(i64, i64)> {
    let (mut pos, mut out) = (r.start, Vec::new());
    for &(op, n) in &r.cigar {
        match op {
            OP_M | OP_EQ | OP_X => {
                out.push((pos, pos + n as i64));
                pos += n as i64;
            }
            OP_D | OP_N => pos += n as i64,
            _ => {}
        }
    }
    out
}

/// `observations`: ({column: [(read, base, qual, eps)]}, steps) for a pile, restricted to `[s0, e0)`.
pub fn observations<'a>(pile: impl IntoIterator<Item = &'a Rec>, prim: &HashMap<ReadId, Prim>, span: (i64, i64)) -> (Obs, Steps) {
    let (mut obs, mut steps) = (Obs::new(), Steps::new());
    let (s0, e0) = span;
    for rec in pile {
        let n = rec.read;
        let (seq, qual, lead) = oriented(&prim[&n], rec);
        let aligned: u64 = rec.cigar.iter().filter(|(o, _)| [OP_M, OP_I, OP_EQ, OP_X].contains(o)).map(|(_, k)| *k as u64).sum();
        let eps = match rec.nm {
            Some(nm) if aligned > 0 => Some(nm as f64 / aligned as f64),
            _ => None,
        };
        let mut st: Vec<(i64, u8, Vec<u8>)> = Vec::new();
        for step in walk(rec) {
            match step {
                Step::M(rp, q) if s0 <= rp && rp < e0 => {
                    let (b, bq) = (seq[q + lead], qual[q + lead]);
                    st.push((rp, b, Vec::new()));
                    obs.entry(rp).or_default().push(ObsX { read: n, base: b, qual: bq, eps });
                }
                Step::D(rp) if s0 <= rp && rp < e0 => st.push((rp, b'-', Vec::new())),
                Step::I(rp, qs, ln) if s0 <= rp && rp < e0 => {
                    if let Some(last) = st.last_mut().filter(|l| l.0 == rp) {
                        last.2.extend_from_slice(&seq[qs + lead..qs + lead + ln]);
                    }
                }
                _ => {}
            }
        }
        steps.insert(n, st);
    }
    (obs, steps)
}

/// The run's BAM: its path and header, opened per call (each worker thread reads on its own).
pub struct Bam {
    pub path: String,
    header: noodles_sam::Header,
    index: noodles_csi::binning_index::Index<noodles_csi::binning_index::index::reference_sequence::index::LinearIndex>,
    contigs: Vec<Arc<str>>,
}

impl Bam {
    pub fn open(path: &str) -> Result<Bam> {
        let mut reader = crate::bam::open_bam(path, 1).with_context(|| format!("opening --bam {path}"))?;
        let header = reader.read_header().with_context(|| format!("reading the header of {path}"))?;
        let bai = format!("{path}.bai");
        let index = noodles_bam::bai::read(&bai).with_context(|| format!("reading the BAM index {bai}"))?;
        let contigs = header.reference_sequences().keys().map(|k| Arc::<str>::from(k.to_string().as_str())).collect();
        Ok(Bam { path: path.to_string(), header, index, contigs })
    }

    /// One record as `Rec` (`read` = the id `names` gives it), or `None` when the read is not wanted or the record unmapped.
    fn rec(&self, r: &noodles_bam::Record, names: &HashMap<Vec<u8>, ReadId>, with_seq: bool) -> Result<Option<Rec>> {
        let flags = r.flags();
        if flags.is_unmapped() {
            return Ok(None);
        }
        let Some(&read) = r.name().and_then(|n| names.get(<&[u8]>::from(n))) else { return Ok(None) };
        let tid = r.reference_sequence_id().transpose()?.context("a mapped record without a reference")?;
        let start = r.alignment_start().transpose()?.context("a mapped record without a start")?.get() as i64 - 1;
        let mut cigar = Vec::new();
        let mut rl = 0i64;
        for op in r.cigar().iter() {
            use noodles_sam::alignment::record::cigar::op::Kind;
            let op = op?;
            let code = match op.kind() {
                Kind::Match => OP_M,
                Kind::Insertion => OP_I,
                Kind::Deletion => OP_D,
                Kind::Skip => OP_N,
                Kind::SoftClip => OP_S,
                Kind::HardClip => OP_H,
                Kind::Pad => 6,
                Kind::SequenceMatch => OP_EQ,
                Kind::SequenceMismatch => OP_X,
            };
            if [OP_M, OP_D, OP_N, OP_EQ, OP_X].contains(&code) {
                rl += op.len() as i64;
            }
            cigar.push((code, op.len() as u32));
        }
        let int = |tag: [u8; 2]| {
            crate::bam::aux_first(r, noodles_sam::alignment::record::data::field::Tag::new(tag[0], tag[1]), crate::bam::value_as_i32)
        };
        Ok(Some(Rec {
            read,
            chrom: self.contigs[tid].clone(),
            start,
            end: start + rl,
            as_: int(*b"AS")?,
            nm: int(*b"NM")?,
            reverse: flags.is_reverse_complemented(),
            secondary: flags.is_secondary(),
            supplementary: flags.is_supplementary(),
            cigar,
            seq: if with_seq { r.sequence().iter().collect() } else { Vec::new() },
        }))
    }

    /// `read_primaries` with `records`: one pass over the whole BAM. Returns ({read: primary}, {read: [record summary]}).
    pub fn read_primaries(&self, names: &HashMap<Vec<u8>, ReadId>) -> Result<(HashMap<ReadId, Prim>, HashMap<ReadId, Vec<RecSummary>>)> {
        let mut reader = crate::bam::open_bam(&self.path, 2)?;
        reader.read_header()?;
        let (mut prim, mut best, mut records): (HashMap<ReadId, (Vec<u8>, Option<Vec<u8>>, bool, Option<(i32, usize)>)>, HashMap<ReadId, i32>, HashMap<ReadId, Vec<RecSummary>>) =
            Default::default();
        let mut buf = noodles_bam::Record::default();
        while reader.read_record(&mut buf)? != 0 {
            let Some(rec) = self.rec(&buf, names, false)? else { continue };
            if let Some(a) = rec.as_ {
                if !rec.supplementary {
                    records.entry(rec.read).or_default().push((rec.chrom.clone(), rec.start, rec.end, a, rec.secondary, introns(&rec)));
                }
                let b = best.entry(rec.read).or_insert(-1_000_000_000);
                *b = (*b).max(a);
            }
            if !rec.secondary && !rec.supplementary && !buf.sequence().is_empty() {
                let seq: Vec<u8> = buf.sequence().iter().collect();
                let q = buf.quality_scores();
                let q: &[u8] = q.as_ref();
                let qual = if q.is_empty() || q[0] == 0xff { None } else { Some(q.to_vec()) };
                prim.insert(rec.read, (seq, qual, rec.reverse, rec.as_.map(|a| (a, introns(&rec)))));
            }
        }
        let prim = prim
            .into_iter()
            .map(|(n, (s, q, rev, pa))| {
                let mut p = Prim::new(s, q, rev, best.get(&n).copied().unwrap_or(-1_000_000_000));
                p.primary = pa;
                (n, p)
            })
            .collect();
        Ok((prim, records))
    }

    /// pysam `fetch(chrom, s, e)` (records overlapping `[s, e)`, file order) of the wanted reads, mapped only.
    pub fn fetch(&self, chrom: &str, s: i64, e: i64, names: &HashMap<Vec<u8>, ReadId>) -> Result<Vec<Rec>> {
        use noodles_core::{Position, Region};
        let mut reader = crate::bam::open_bam(&self.path, 1)?;
        reader.read_header()?;
        let region = Region::new(chrom, Position::try_from(s.max(0) as usize + 1)?..=Position::try_from(e.max(s + 1) as usize)?);
        let mut out = Vec::new();
        for r in reader.query(&self.header, &self.index, &region).with_context(|| format!("querying {} at {chrom}:{s}-{e}", self.path))? {
            if let Some(rec) = self.rec(&r?, names, true)? {
                if rec.start < e && rec.end > s {
                    out.push(rec);
                }
            }
        }
        Ok(out)
    }
}

/// `pile_records`: {copy: {read: record}} for the records in the tie window (`Prim::near`; `structure` = `--tie-structure`) that
/// fall on the copy's span.
pub fn pile_records(
    bam: &Bam,
    copies: &[Copy],
    reads: &BTreeSet<ReadId>,
    prim: &HashMap<ReadId, Prim>,
    names: &HashMap<Vec<u8>, ReadId>,
    structure: bool,
) -> Result<Piles> {
    let mut piles = Piles::new();
    for c in copies {
        let mut p = Pile::new();
        for rec in bam.fetch(&c.chrom, c.start, c.end, names)? {
            let n = rec.read;
            let (Some(pr), Some(a)) = (prim.get(&n), rec.as_) else { continue };
            if !reads.contains(&n) || !pr.near(a, introns(&rec), structure) {
                continue;
            }
            if p.get(&n).map_or(true, |old| a > old.as_.unwrap()) {
                p.insert(n, rec);
            }
        }
        piles.insert(c.tid.clone(), p);
    }
    Ok(piles)
}

#[cfg(test)]
mod tests {
    use super::*;

    pub(crate) fn rec(start: i64, cigar: &[(u8, u32)], reverse: bool) -> Rec {
        let rl: i64 = cigar.iter().filter(|(o, _)| [OP_M, OP_D, OP_N, OP_EQ, OP_X].contains(o)).map(|(_, n)| *n as i64).sum();
        Rec {
            read: 0,
            chrom: Arc::from("c"),
            start,
            end: start + rl,
            as_: Some(0),
            nm: None,
            reverse,
            secondary: false,
            supplementary: false,
            cigar: cigar.to_vec(),
            seq: Vec::new(),
        }
    }

    #[test]
    fn walk_matches_python() {
        // list(o3port.walk(r)) for reference_start 100, cigartuples [(5,5),(4,3),(0,4),(1,2),(0,3),(2,1),(3,2),(0,4)]
        let r = rec(100, &[(OP_H, 5), (OP_S, 3), (OP_M, 4), (OP_I, 2), (OP_M, 3), (OP_D, 1), (OP_N, 2), (OP_M, 4)], false);
        use Step::*;
        assert_eq!(
            walk(&r),
            vec![
                M(100, 3),
                M(101, 4),
                M(102, 5),
                M(103, 6),
                I(103, 7, 2),
                M(104, 9),
                M(105, 10),
                M(106, 11),
                D(107),
                M(110, 12),
                M(111, 13),
                M(112, 14),
                M(113, 15)
            ]
        );
    }

    #[test]
    fn walk_hard_clip_reverse() {
        let p = Prim::new(b"AACGTN".to_vec(), Some(vec![1, 2, 3, 4, 5, 6]), false, 0);
        let r = rec(0, &[(OP_H, 5), (OP_M, 1)], true);
        let (s, q, lead) = oriented(&p, &r);
        assert_eq!(s, b"NACGTT");
        assert_eq!(q, &[6, 5, 4, 3, 2, 1]);
        assert_eq!(lead, 5);
        let (s, q, lead) = oriented(&p, &rec(0, &[(OP_M, 6)], false));
        assert_eq!((s, q, lead), (&b"AACGTN"[..], &[1u8, 2, 3, 4, 5, 6][..], 0));
    }

    #[test]
    fn tie_window_by_best_as_or_by_the_primarys_structure() {
        // the P15 case: spliced primary (5 introns, AS 592) and an unspliced processed-copy secondary with a higher AS (626)
        let mut p = Prim::new(b"ACGT".to_vec(), None, false, 626);
        p.primary = Some((592, 5));
        // by best AS (the default): the secondary ties, the read's own spliced primary does not (592 < 626 - 10)
        assert!(p.near(626, 0, false));
        assert!(!p.near(592, 5, false));
        // by structure (--tie-structure, register 1404): same intron count as the primary and AS within 10 of it
        assert!(p.near(592, 5, true));
        assert!(!p.near(626, 0, true));
        assert!(p.near(600, 5, true) && p.near(582, 5, true));
        assert!(!p.near(603, 5, true) && !p.near(581, 5, true));
        assert!(!p.near(592, 4, true));
        // no primary AS: nothing ties by structure
        p.primary = None;
        assert!(!p.near(592, 5, true));
        assert_eq!(introns(&rec(0, &[(OP_M, 5), (OP_N, 100), (OP_M, 5), (OP_N, 7), (OP_M, 1)], false)), 2);
    }

    #[test]
    fn missing_quals_are_40() {
        let p = Prim::new(b"ACG".to_vec(), None, false, 0);
        assert_eq!(p.qual, vec![40, 40, 40]);
    }

    #[test]
    fn observations_insertions_deletions_and_span() {
        // read ACGTTA on 10M? : 2M 2I 1M 1D 1M from ref 10; span [10, 14)
        let mut r = rec(10, &[(OP_M, 2), (OP_I, 2), (OP_M, 1), (OP_D, 1), (OP_M, 1)], false);
        r.nm = Some(3);
        let prim = HashMap::from([(0, Prim::new(b"ACGTTA".to_vec(), Some(vec![30; 6]), false, 0))]);
        let (obs, steps) = observations([&r], &prim, (10, 14));
        assert_eq!(obs.keys().copied().collect::<Vec<_>>(), vec![10, 11, 12]);
        assert_eq!(obs[&12][0], ObsX { read: 0, base: b'T', qual: 30, eps: Some(3.0 / 6.0) });
        assert_eq!(steps[&0], vec![(10, b'A', Vec::new()), (11, b'C', b"GT".to_vec()), (12, b'T', Vec::new()), (13, b'-', Vec::new())]);
    }
}

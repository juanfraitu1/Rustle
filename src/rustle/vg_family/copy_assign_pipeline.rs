//! Per-read COPY ASSIGNMENT driver (integration stage 4a) — the orchestration of `bench/copy_assign.py`.
//!
//! Given a family's copies (`DenovoTranscript`: spliced sequence + exon structure + strand) and the reads
//! over the family region, assign each read — especially the hard multimappers minimap2 leaves at MAPQ 0 —
//! to a specific paralog COPY, via paralog-sequence-variant (PSV) bases + copy-specific junctions. This is
//! the "resolve the ambiguity, pick one mapping" step on top of the family layer.
//!
//! It builds the `CopyProfile`/`ReadFeatures` the already-ported `copy_assign::assign_read` consumes:
//!   1. `discover_psvs` — all-pairs alignment vs copy[0] → columns where copies differ (per-copy genomic
//!      position + transcription-strand base).
//!   2. spliced↔genomic exon map (`exon_map`/`gen2off`) + intron `copy_boundaries`.
//!   3. per read: read its base at each column's genomic position (reverse-complemented for a `-` copy) and
//!      its intron boundaries (mapped to spliced space) → a feature vector → `assign_read`.

use std::collections::{BTreeMap, BTreeSet};

use super::copy_assign::{assign_read, AssignParams, Assignment, CopyProfile, ReadFeatures};
use super::copy_split::{allele_at, intron_chain_of, AlignedRead};
use super::family_detect::DenovoTranscript;
use super::family_graph::poa_msa_with_costs;

fn rc_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        other => other,
    }
}
fn is_acgt(b: u8) -> bool {
    matches!(b, b'A' | b'C' | b'G' | b'T')
}

/// Exon coordinates `[start,donor0), [acc0,donor1), ..., [acc_last,end)` from the intron chain + extents.
fn exons_of(t: &DenovoTranscript) -> Vec<(u64, u64)> {
    let mut exons = Vec::with_capacity(t.introns.len() + 1);
    let mut prev = t.start;
    for &(d, a) in &t.introns {
        exons.push((prev, d));
        prev = a;
    }
    exons.push((prev, t.end));
    exons
}

/// Reference end (0-based exclusive) of an aligned read = ref_start + reference-consuming CIGAR.
fn read_ref_end(read: &AlignedRead) -> u64 {
    let mut end = read.ref_start;
    for &(op, len) in &read.cigar {
        if matches!(op, 'M' | '=' | 'X' | 'D' | 'N') {
            end += len;
        }
    }
    end
}

/// Map each spliced-sequence offset → its 0-based FORWARD-genome coordinate, in transcription order
/// (exons reversed, each descending, for a `-` copy). Mirrors `copy_idx2gen`; `len == spliced seq len`.
pub fn exon_map(t: &DenovoTranscript) -> Vec<u64> {
    let exons = exons_of(t);
    let mut g = Vec::new();
    if t.strand == '-' {
        for &(xs, xe) in exons.iter().rev() {
            g.extend((xs..xe).rev());
        }
    } else {
        for &(xs, xe) in &exons {
            g.extend(xs..xe);
        }
    }
    g
}

/// Inverse of `exon_map`: forward-genome coord → spliced offset.
pub fn gen2off(t: &DenovoTranscript) -> BTreeMap<u64, usize> {
    exon_map(t).into_iter().enumerate().map(|(i, gp)| (gp, i)).collect()
}

/// Intron-boundary offsets in transcription-spliced space (cumulative exon lengths, excluding the last).
/// Mirrors `copy_boundaries`.
pub fn copy_boundaries(t: &DenovoTranscript) -> Vec<i64> {
    let exons = exons_of(t);
    let mut lens: Vec<u64> = exons.iter().map(|&(xs, xe)| xe - xs).collect();
    if t.strand == '-' {
        lens.reverse();
    }
    let mut bnd = Vec::new();
    let mut cum = 0i64;
    for &l in &lens[..lens.len().saturating_sub(1)] {
        cum += l as i64;
        bnd.push(cum);
    }
    bnd
}

/// One PSV column: for each copy (by index) its `(forward-genome position, transcription-strand base)` at
/// this column, or `None` if the copy is gapped here.
pub type PsvColumn = Vec<Option<(u64, u8)>>;

/// Discover PSV columns by all-pairs alignment of every copy vs copy[0]: columns are ref offsets where
/// `>= 1` copy's aligned base differs. EVERY copy gets its aligned `(genomic, base)` at a column (a copy
/// matching the ref inherits that base, not `None` — the python bug guard). `exon_maps[c] = exon_map(copies[c])`.
///
/// Invariant (held by gate-produced transcripts): `exon_maps[c].len() == copies[c].seq.len()` (the spliced
/// sequence is built from the exon coordinates), so the per-column `exon_maps[c][off]` / `seq[off]` indexing
/// is in range — the analogue of python `assign_family`'s `len(s) == len(g)` filter.
pub fn discover_psvs(copies: &[&DenovoTranscript], exon_maps: &[Vec<u64>]) -> Vec<PsvColumn> {
    use poasta::aligner::scoring::GapAffine;
    let n = copies.len();
    if n < 2 {
        return Vec::new();
    }
    let ref_seq = &copies[0].seq;
    // for each non-ref copy: amap (ref_off -> other_off) over both-non-gap columns; collect differing offs.
    let mut amaps: Vec<BTreeMap<usize, usize>> = vec![BTreeMap::new(); n]; // index 0 unused
    let mut diff_off: BTreeSet<usize> = BTreeSet::new();
    for other in 1..n {
        // strong gap-open anchors the conserved core column-for-column (same config as contiguous_core_coverage).
        let msa = match poa_msa_with_costs(&[ref_seq.clone(), copies[other].seq.clone()], GapAffine::new(1, 1, 32)) {
            Ok(m) if m.len() == 2 => m,
            _ => continue,
        };
        let (r0, r1) = (&msa[0], &msa[1]);
        let (mut ro, mut oo) = (0usize, 0usize);
        for c in 0..r0.len().min(r1.len()) {
            let (ca, cb) = (r0[c], r1[c]);
            let (a_gap, b_gap) = (ca == b'-', cb == b'-');
            if !a_gap && !b_gap {
                amaps[other].insert(ro, oo);
                if is_acgt(ca) && is_acgt(cb) && ca != cb {
                    diff_off.insert(ro);
                }
            }
            if !a_gap {
                ro += 1;
            }
            if !b_gap {
                oo += 1;
            }
        }
    }
    diff_off
        .into_iter()
        .map(|ro| {
            let mut rec: PsvColumn = vec![None; n];
            rec[0] = Some((exon_maps[0][ro], ref_seq[ro]));
            for other in 1..n {
                if let Some(&oo) = amaps[other].get(&ro) {
                    rec[other] = Some((exon_maps[other][oo], copies[other].seq[oo]));
                }
            }
            rec
        })
        .collect()
}

/// Per-copy assignment profiles + the genomic PSV positions (to read a read in each copy's frame).
pub struct FamilyProfiles {
    pub profiles: Vec<CopyProfile>,
    /// `copy_gpos[c]` = `{column -> forward-genome PSV position}` for copy `c`.
    pub copy_gpos: Vec<BTreeMap<usize, u64>>,
    pub n_cols: usize,
}

/// Build the per-copy `CopyProfile`s (PSV alleles + intron boundaries) and genomic PSV positions.
pub fn build_family_profiles(copies: &[&DenovoTranscript]) -> FamilyProfiles {
    let exon_maps: Vec<Vec<u64>> = copies.iter().map(|c| exon_map(c)).collect();
    let cols = discover_psvs(copies, &exon_maps);
    let n = copies.len();
    let mut profiles = Vec::with_capacity(n);
    let mut copy_gpos = vec![BTreeMap::new(); n];
    for (ci, c) in copies.iter().enumerate() {
        let alleles: Vec<Option<u8>> = cols.iter().map(|col| col[ci].map(|(_, b)| b)).collect();
        for (j, col) in cols.iter().enumerate() {
            if let Some((g, _)) = col[ci] {
                copy_gpos[ci].insert(j, g);
            }
        }
        profiles.push(CopyProfile { copy_id: ci, alleles, junctions: copy_boundaries(c) });
    }
    FamilyProfiles { profiles, copy_gpos, n_cols: cols.len() }
}

/// Assign one read to a copy, reading its PSV bases in `mapped_copy`'s genomic frame (reverse-complemented
/// for a `-` copy) and its intron boundaries via that copy's `gen2off`.
pub fn assign_one_read(
    read: &AlignedRead,
    mapped_copy: usize,
    copies: &[&DenovoTranscript],
    fp: &FamilyProfiles,
    p: &AssignParams,
) -> Option<Assignment> {
    let mc = mapped_copy;
    let minus = copies[mc].strand == '-';
    // PSV observations: the read's base at each column's genomic position in the mapped copy's frame.
    let mut psv_obs = vec![None; fp.n_cols];
    for (&col, &g) in &fp.copy_gpos[mc] {
        if let Some(b) = allele_at(read, g) {
            // copies' alleles are transcription-strand; RC the forward-genome base for a '-' copy.
            psv_obs[col] = Some(if minus { rc_base(b) } else { b });
        }
    }
    // the read's intron boundaries in the mapped copy's spliced space, via gen2off. `checked_sub` (not
    // saturating) so a donor at genome coord 0 yields no key — exactly like python's `g2o.get(d0 - 1)` where
    // `d0 - 1 == -1` is never a key (a saturating 0 could spuriously match genome position 0).
    let g2o = gen2off(copies[mc]);
    let mut junctions = Vec::new();
    for (d0, a0) in intron_chain_of(read) {
        let o1 = d0.checked_sub(1).and_then(|k| g2o.get(&k));
        if let (Some(&o1), Some(&o2)) = (o1, g2o.get(&a0)) {
            junctions.push(o1.max(o2) as i64);
        }
    }
    assign_read(&ReadFeatures { psv_obs, junctions }, &fp.profiles, p)
}

/// Assign every read over a co-located family to a copy. Each read is mapped to the copy whose genomic span
/// it overlaps most (reads overlapping no copy are skipped). Returns `(read_index, Assignment)`. Mirrors
/// `assign_family`.
pub fn assign_family(
    copies: &[&DenovoTranscript],
    reads: &[AlignedRead],
    p: &AssignParams,
) -> Vec<(usize, Assignment)> {
    if copies.len() < 2 {
        return Vec::new();
    }
    let fp = build_family_profiles(copies);
    let mut out = Vec::new();
    for (ri, read) in reads.iter().enumerate() {
        let r_end = read_ref_end(read);
        // map the read to the copy whose genomic span it overlaps most.
        let mut best: Option<usize> = None;
        let mut best_ov = 0i64;
        for (ci, c) in copies.iter().enumerate() {
            let ov = (r_end.min(c.end) as i64) - (read.ref_start.max(c.start) as i64);
            if ov > best_ov {
                best_ov = ov;
                best = Some(ci);
            }
        }
        if let Some(mc) = best {
            if let Some(a) = assign_one_read(read, mc, copies, &fp, p) {
                out.push((ri, a));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::super::copy_assign::AssignStatus;
    use super::*;

    struct SplitMix64(u64);
    impl SplitMix64 {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }
    }
    fn rand_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut rng = SplitMix64(seed);
        const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
        (0..n).map(|_| B[(rng.next_u64() % 4) as usize]).collect()
    }
    fn copy_tx(tid: &str, s: u64, e: u64, strand: char, introns: &[(u64, u64)], seq: Vec<u8>) -> DenovoTranscript {
        DenovoTranscript {
            tid: tid.into(),
            chrom: "c1".into(),
            start: s,
            end: e,
            n_reads: 5,
            strand,
            introns: introns.to_vec(),
            seq,
        }
    }

    // ---- exon_map / gen2off / copy_boundaries ----

    #[test]
    fn exon_map_plus_minus_and_spliced() {
        let plus = copy_tx("p", 10, 15, '+', &[], vec![b'A'; 5]);
        assert_eq!(exon_map(&plus), vec![10, 11, 12, 13, 14]);
        let minus = copy_tx("m", 10, 15, '-', &[], vec![b'A'; 5]);
        assert_eq!(exon_map(&minus), vec![14, 13, 12, 11, 10]);
        let two = copy_tx("t", 10, 23, '+', &[(13, 20)], vec![b'A'; 6]);
        assert_eq!(exon_map(&two), vec![10, 11, 12, 20, 21, 22]);
    }

    #[test]
    fn gen2off_inverts_exon_map() {
        let t = copy_tx("t", 10, 23, '+', &[(13, 20)], vec![b'A'; 6]);
        let m = gen2off(&t);
        assert_eq!(m.get(&10), Some(&0));
        assert_eq!(m.get(&20), Some(&3)); // first base of the 2nd exon
        assert_eq!(m.get(&13), None); // inside the intron
    }

    #[test]
    fn copy_boundaries_are_cumulative_exon_lengths() {
        // exons [0,2)=2, [10,13)=3, [20,24)=4 -> boundaries 2, 5.
        let t = copy_tx("u", 0, 24, '+', &[(2, 10), (13, 20)], vec![b'A'; 9]);
        assert_eq!(copy_boundaries(&t), vec![2, 5]);
        // minus strand reverses the exon lengths: 4,3,2 -> boundaries 4, 7.
        let m = copy_tx("v", 0, 24, '-', &[(2, 10), (13, 20)], vec![b'A'; 9]);
        assert_eq!(copy_boundaries(&m), vec![4, 7]);
    }

    // ---- discover_psvs ----

    #[test]
    fn discover_psvs_finds_differing_columns() {
        let sa = rand_seq(200, 0xD1);
        let mut sb = sa.clone();
        sb[40] = if sa[40] == b'A' { b'C' } else { b'A' };
        sb[120] = if sa[120] == b'G' { b'T' } else { b'G' };
        let (a40, b40, a120, b120) = (sa[40], sb[40], sa[120], sb[120]);
        let ca = copy_tx("A", 0, 200, '+', &[], sa);
        let cb = copy_tx("B", 1000, 1200, '+', &[], sb);
        let copies = [&ca, &cb];
        let maps = [exon_map(&ca), exon_map(&cb)];
        let cols = discover_psvs(&copies, &maps);
        assert_eq!(cols.len(), 2, "two SNP columns");
        assert_eq!(cols[0], vec![Some((40, a40)), Some((1040, b40))]);
        assert_eq!(cols[1], vec![Some((120, a120)), Some((1120, b120))]);
    }

    // ---- assign_family (end to end) ----

    #[test]
    fn assign_resolves_multimapper_to_its_true_copy() {
        // two single-exon copies identical except at 3 PSV columns (A vs C). A read aligned to copyA's
        // region but carrying copyB's alleles must be assigned to copyB (the multimapper resolution).
        let base = rand_seq(300, 0xBA5E);
        let psv = [50usize, 150, 250];
        let mut sa = base.clone();
        let mut sb = base.clone();
        for &p in &psv {
            sa[p] = b'A';
            sb[p] = b'C';
        }
        let ca = copy_tx("A", 0, 300, '+', &[], sa);
        let cb = copy_tx("B", 1000, 1300, '+', &[], sb.clone());
        let copies = [&ca, &cb];
        let read = AlignedRead { ref_start: 0, cigar: vec![('M', 300)], seq: sb };
        let res = assign_family(&copies, &[read], &AssignParams::default());
        assert_eq!(res.len(), 1);
        let (ri, a) = &res[0];
        assert_eq!(*ri, 0);
        assert_eq!(a.best_copy, 1, "carries copyB alleles -> copyB, despite aligning to copyA's region");
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn assign_read_spanning_no_psv_is_tied() {
        // copies differ only at offsets >= 200; a read covering only [0,150) spans no PSV -> tied.
        let base = rand_seq(300, 0xB16);
        let mut sa = base.clone();
        let mut sb = base.clone();
        for &p in &[210usize, 250, 290] {
            sa[p] = b'A';
            sb[p] = b'C';
        }
        let ca = copy_tx("A", 0, 300, '+', &[], sa);
        let cb = copy_tx("B", 1000, 1300, '+', &[], sb);
        let copies = [&ca, &cb];
        // read covers [0,150): spans none of the PSV columns (all >= 210).
        let read = AlignedRead { ref_start: 0, cigar: vec![('M', 150)], seq: base[0..150].to_vec() };
        let res = assign_family(&copies, &[read], &AssignParams::default());
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].1.status, AssignStatus::Tied, "spans no decisive feature");
    }

    #[test]
    fn assign_minus_strand_copy_reverse_complements_read_base() {
        // copyB is on the '-' strand: its allele vector is in transcription space, so the read's
        // forward-genome base must be reverse-complemented before comparison. copyA '+', copyB '-'.
        // Build so copyB's transcription seq differs from copyA only at PSVs.
        let base = rand_seq(240, 0x515E);
        let psv = [60usize, 120, 180];
        let mut sa = base.clone(); // copyA transcription seq (= forward genome, + strand)
        let mut sb_txn = base.clone(); // copyB transcription seq (- strand)
        for &p in &psv {
            sa[p] = b'A';
            sb_txn[p] = b'C';
        }
        let ca = copy_tx("A", 0, 240, '+', &[], sa);
        let cb = copy_tx("B", 1000, 1240, '-', &[], sb_txn.clone());
        let copies = [&ca, &cb];
        // a true-copyB read aligned to copyB's region: forward-genome seq = revcomp(transcription seq).
        let fwd_b: Vec<u8> = sb_txn.iter().rev().map(|&b| rc_base(b)).collect();
        let read = AlignedRead { ref_start: 1000, cigar: vec![('M', 240)], seq: fwd_b };
        let res = assign_family(&copies, &[read], &AssignParams::default());
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].1.best_copy, 1, "minus-strand read RC'd to transcription space -> copyB");
        assert_eq!(res[0].1.status, AssignStatus::Assigned);
    }
}

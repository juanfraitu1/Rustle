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
use super::copy_split::{intron_chain_of, AlignedRead};
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
pub(crate) fn read_ref_end(read: &AlignedRead) -> u64 {
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

/// Per-copy assignment profiles + cached frames (so reading a read costs no per-read map rebuilds).
pub struct FamilyProfiles {
    pub profiles: Vec<CopyProfile>,
    /// `copy_gpos[c]` = `[(column, forward-genome PSV position)]` for copy `c` (sorted by genomic position
    /// for a single-pass CIGAR walk).
    pub copy_gpos: Vec<Vec<(usize, u64)>>,
    /// `gen2off[c]` = forward-genome coord → spliced offset for copy `c` (cached once, not per read).
    pub gen2off: Vec<BTreeMap<u64, usize>>,
    /// `strand[c]` for each copy.
    pub strand: Vec<char>,
    pub n_cols: usize,
}

/// Build the per-copy `CopyProfile`s (PSV alleles + intron boundaries), the genomic PSV positions, and the
/// per-copy `gen2off`/strand — everything a read needs, computed ONCE per family (not per read).
pub fn build_family_profiles(copies: &[&DenovoTranscript]) -> FamilyProfiles {
    let exon_maps: Vec<Vec<u64>> = copies.iter().map(|c| exon_map(c)).collect();
    let cols = discover_psvs(copies, &exon_maps);
    let n = copies.len();
    let mut profiles = Vec::with_capacity(n);
    let mut copy_gpos: Vec<Vec<(usize, u64)>> = vec![Vec::new(); n];
    for (ci, c) in copies.iter().enumerate() {
        let alleles: Vec<Option<u8>> = cols.iter().map(|col| col[ci].map(|(_, b)| b)).collect();
        for (j, col) in cols.iter().enumerate() {
            if let Some((g, _)) = col[ci] {
                copy_gpos[ci].push((j, g));
            }
        }
        copy_gpos[ci].sort_unstable_by_key(|&(_, g)| g); // genomic order for the single CIGAR sweep
        profiles.push(CopyProfile { copy_id: ci, alleles, junctions: copy_boundaries(c) });
    }
    FamilyProfiles {
        profiles,
        copy_gpos,
        gen2off: copies.iter().map(|c| gen2off(c)).collect(),
        strand: copies.iter().map(|c| c.strand).collect(),
        n_cols: cols.len(),
    }
}

/// Assign one read to a copy, reading its PSV bases in `mapped_copy`'s genomic frame (reverse-complemented
/// for a `-` copy) and its intron boundaries via that copy's `gen2off`.
pub fn assign_one_read(
    read: &AlignedRead,
    mapped_copy: usize,
    fp: &FamilyProfiles,
    p: &AssignParams,
) -> Option<Assignment> {
    assign_read(&read_features(read, mapped_copy, fp), &fp.profiles, p)
}

/// Fill `psv_obs[col] = Some(base)` for the read's base at each `(col, genomic-position)` (sorted by
/// position) in a SINGLE CIGAR walk — O(cigar + positions) instead of O(cigar * positions). `None` is left
/// for positions inside an intron/deletion or outside the read. `minus` reverse-complements the base.
fn fill_psv_obs(read: &AlignedRead, gpos: &[(usize, u64)], minus: bool, psv_obs: &mut [Option<u8>]) {
    let mut pi = 0;
    while pi < gpos.len() && gpos[pi].1 < read.ref_start {
        pi += 1; // before the read
    }
    let mut ref_cur = read.ref_start;
    let mut seq_cur = 0u64;
    for &(op, len) in &read.cigar {
        if pi >= gpos.len() {
            break;
        }
        match op {
            'M' | '=' | 'X' => {
                let block_end = ref_cur + len;
                while pi < gpos.len() && gpos[pi].1 < block_end {
                    let (col, g) = gpos[pi];
                    if let Some(&b) = read.seq.get((seq_cur + (g - ref_cur)) as usize) {
                        psv_obs[col] = Some(if minus { rc_base(b) } else { b });
                    }
                    pi += 1;
                }
                ref_cur = block_end;
                seq_cur += len;
            }
            'N' | 'D' => {
                let block_end = ref_cur + len;
                while pi < gpos.len() && gpos[pi].1 < block_end {
                    pi += 1; // inside intron/deletion -> unmatched
                }
                ref_cur = block_end;
            }
            'I' | 'S' => seq_cur += len,
            _ => {}
        }
    }
}

/// The read's feature vector in the mapped copy's frame: PSV bases at each column's genomic position
/// (reverse-complemented for a `-` copy, single CIGAR walk) + intron boundaries via the cached `gen2off`.
fn read_features(read: &AlignedRead, mc: usize, fp: &FamilyProfiles) -> ReadFeatures {
    let mut psv_obs = vec![None; fp.n_cols];
    fill_psv_obs(read, &fp.copy_gpos[mc], fp.strand[mc] == '-', &mut psv_obs);
    // intron boundaries in the mapped copy's spliced space, via the cached gen2off. `checked_sub` (not
    // saturating) so a donor at genome coord 0 yields no key — exactly like python's `g2o.get(d0 - 1)`.
    let g2o = &fp.gen2off[mc];
    let mut junctions = Vec::new();
    for (d0, a0) in intron_chain_of(read) {
        let o1 = d0.checked_sub(1).and_then(|k| g2o.get(&k));
        if let (Some(&o1), Some(&o2)) = (o1, g2o.get(&a0)) {
            junctions.push(o1.max(o2) as i64);
        }
    }
    ReadFeatures { psv_obs, junctions }
}

/// The copy whose genomic span the read overlaps most (`None` if it overlaps none).
fn best_overlap_copy(read: &AlignedRead, copies: &[&DenovoTranscript]) -> Option<usize> {
    let r_end = read_ref_end(read);
    let mut best = None;
    let mut best_ov = 0i64;
    for (ci, c) in copies.iter().enumerate() {
        let ov = (r_end.min(c.end) as i64) - (read.ref_start.max(c.start) as i64);
        if ov > best_ov {
            best_ov = ov;
            best = Some(ci);
        }
    }
    best
}

/// Two-pass per-read assignment: the mapped copy, the PSV-ONLY assignment, and the PSV+JUNCTION assignment
/// (the decision). Mirrors copy_assign.py::assign_family's two `assign_read` calls.
#[derive(Clone, Debug)]
pub struct ReadResult {
    pub read_index: usize,
    pub mapped_copy: usize,
    pub psv: Assignment,
    pub combined: Assignment,
    /// this read's base at each PSV column (None = uncovered) — the raw per-molecule evidence the assignment
    /// is built from (for the assignment-proof genotype visualization).
    pub psv_obs: Vec<Option<u8>>,
}

/// Two-pass assignment detail for a family: the per-read results + the PSV-column count + the SOFT per-copy
/// abundance (EM over per-read PSV likelihoods) with a normal-approx 95% CI half-width per copy.
#[derive(Clone, Debug)]
pub struct FamilyDetail {
    pub results: Vec<ReadResult>,
    pub n_cols: usize,
    /// `copy_abundance[c]` = estimated fraction of reads from copy `c` (sums to 1).
    pub copy_abundance: Vec<f64>,
    /// `copy_abundance_ci[c]` = 95% CI half-width for `copy_abundance[c]`.
    pub copy_abundance_ci: Vec<f64>,
    /// reads whose per-PSV copy pattern SWITCHES mid-molecule (gene-conversion / recombinant candidates).
    pub mosaic_reads: usize,
    /// family-confirmed gene conversions: the breakpoint RECURS across independent molecules (vs a one-off
    /// chimera). The enriched per-molecule signal the multimappers carry beyond presence/abundance.
    pub conversions: Vec<super::mosaic::ConversionEvent>,
    /// COPY-level historical gene conversions: a de-novo copy whose PSV-allele vector is itself a mosaic of two
    /// OTHER copies (the APOBEC3/RFPL signal — baked into the copy sequence, invisible to the read-level scan).
    pub copy_conversions: Vec<CopyConversion>,
    /// PSV column `j` → canonical forward-genome position (the column frame shared by reads + copies) — for the
    /// genotype-matrix visualization that proves the per-read assignment.
    pub psv_col_pos: Vec<Option<u64>>,
    /// `copy_psv_alleles[c][j]` = copy `c`'s allele at PSV column `j` (None = gapped) — the per-copy reference
    /// the reads are matched against in the genotype matrix.
    pub copy_psv_alleles: Vec<Vec<Option<u8>>>,
}

/// One de-novo COPY whose PSV-allele vector is a MOSAIC of two other copies — a HISTORICAL gene conversion
/// (the converted copy resembles donor A over its 5' tract, donor B over its 3' tract). Distinct from the
/// read-level recombinant-molecule scan: this is the copy-level signal the famous conversion families carry.
#[derive(Clone, Debug)]
pub struct CopyConversion {
    pub copy_c: usize,
    pub copy_a: usize,
    pub copy_b: usize,
    pub breakpoint: (u64, u64),
    pub n_decisive: usize,
}

/// For each copy in a (>=3-copy) family, ask whether its PSV-allele vector switches from matching one copy to
/// another at a breakpoint = a historical gene conversion. Runs the same `detect_mosaic` primitive on the COPY
/// consensus (self excluded from the match set, else it trivially matches itself). `col_gpos[col]` = the
/// column's canonical genomic position (for ordering + the breakpoint frame).
pub fn scan_copy_conversions(
    profiles: &[CopyProfile],
    col_gpos: &[Option<u64>],
    p: &super::mosaic::MosaicParams,
) -> Vec<CopyConversion> {
    use super::mosaic::{detect_mosaic, SiteObs};
    let n = profiles.len();
    if n < 3 {
        return Vec::new(); // a copy can only be a mosaic of two OTHERS if there are >= 3 copies
    }
    let mut out = Vec::new();
    for c in 0..n {
        let mut site_obs: Vec<SiteObs> = Vec::new();
        for col in 0..col_gpos.len() {
            let Some(rpos) = col_gpos[col] else { continue };
            let Some(ca) = profiles[c].alleles.get(col).copied().flatten() else { continue };
            let match_bits = (0..n).map(|o| o != c && profiles[o].alleles[col] == Some(ca)).collect();
            site_obs.push(SiteObs { ref_pos: rpos, match_bits });
        }
        site_obs.sort_by_key(|s| s.ref_pos);
        // a copy consensus is fully decisive at every PSV (no read error), so a divergent family yields
        // thousands of sites -> stride-sample to bound `detect_mosaic` (breakpoint detection needs spread, not
        // every column).
        const MAX_CONV_SITES: usize = 300;
        if site_obs.len() > MAX_CONV_SITES {
            let stride = site_obs.len() / MAX_CONV_SITES + 1;
            site_obs = site_obs.into_iter().step_by(stride).collect();
        }
        let call = detect_mosaic(&site_obs, n, 0.01, p);
        if let (true, Some(a), Some(b), Some(bp)) =
            (call.is_mosaic(), call.copy_a, call.copy_b, call.breakpoint_ref)
        {
            out.push(CopyConversion { copy_c: c, copy_a: a, copy_b: b, breakpoint: bp, n_decisive: call.n_decisive });
        }
    }
    out
}

/// Default per-base error for the soft-quant likelihood (HiFi-ish; match prob `1-e`, each mismatch `e/3`).
pub const QUANT_ERROR: f64 = 0.01;

/// SOFT per-copy abundance via EM over per-read PSV likelihoods — the estimator the quantification benchmark
/// showed beats hard assignment at sparse PSVs (it uses the partial evidence instead of discarding it).
/// `read_obs[r][j]` = read r's base at PSV column j (None = uncovered); `copy_alleles[c][j]` = copy c's allele
/// at column j (None = copy lacks it). Returns the per-copy abundance (sums to 1). At zero informative PSVs it
/// returns the uniform prior — the honest identifiability floor (no information => no resolution).
pub fn soft_quantify_em(
    read_obs: &[Vec<Option<u8>>],
    copy_alleles: &[Vec<Option<u8>>],
    error: f64,
    iters: usize,
) -> Vec<f64> {
    let k = copy_alleles.len();
    if k == 0 {
        return Vec::new();
    }
    if read_obs.is_empty() {
        return vec![1.0 / k as f64; k];
    }
    let (ln_match, ln_mis) = ((1.0 - error).ln(), (error / 3.0).ln());
    // per-read log P(read | copy c) over the covered PSV columns.
    let loglik: Vec<Vec<f64>> = read_obs
        .iter()
        .map(|obs| {
            copy_alleles
                .iter()
                .map(|ca| {
                    let mut ll = 0.0;
                    for (j, o) in obs.iter().enumerate() {
                        if let (Some(ob), Some(Some(al))) = (o, ca.get(j)) {
                            ll += if ob == al { ln_match } else { ln_mis };
                        }
                    }
                    ll
                })
                .collect()
        })
        .collect();
    let mut theta = vec![1.0 / k as f64; k];
    for _ in 0..iters {
        let mut sum = vec![0.0; k];
        for ll in &loglik {
            let logpost: Vec<f64> = (0..k).map(|c| theta[c].max(1e-12).ln() + ll[c]).collect();
            let m = logpost.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let z: f64 = logpost.iter().map(|&x| (x - m).exp()).sum();
            for c in 0..k {
                sum[c] += (logpost[c] - m).exp() / z;
            }
        }
        let n = read_obs.len() as f64;
        for c in 0..k {
            theta[c] = sum[c] / n;
        }
    }
    theta
}

/// Like `assign_family` but returns the TWO-PASS detail per read so callers can report how many reads a
/// copy-specific junction resolved that PSVs alone could not (`junction_only`), and — with read mapq — the
/// silver-standard unique-mapper agreement. Reads overlapping no copy are skipped.
pub fn assign_family_detailed(
    copies: &[&DenovoTranscript],
    reads: &[AlignedRead],
    p: &AssignParams,
) -> FamilyDetail {
    if copies.len() < 2 {
        return FamilyDetail {
            results: Vec::new(),
            n_cols: 0,
            copy_abundance: Vec::new(),
            copy_abundance_ci: Vec::new(),
            mosaic_reads: 0,
            conversions: Vec::new(),
            copy_conversions: Vec::new(),
            psv_col_pos: Vec::new(),
            copy_psv_alleles: Vec::new(),
        };
    }
    use super::mosaic::{aggregate_family, detect_mosaic, MosaicParams, SiteObs};
    const MOSAIC_EPS: f64 = 0.01; // HiFi per-base error for the mosaic likelihood
    const MAX_MOSAIC_SITES: usize = 250; // cap PSV sites per detect_mosaic (it is O(sites^2)); stride-sample
    let fp = build_family_profiles(copies);
    // canonical column -> genomic position (first copy that has the column) so every read's switch breakpoints
    // are in ONE frame and `aggregate_family` can cluster recurrences across molecules.
    let mut col_canon: Vec<Option<u64>> = vec![None; fp.n_cols];
    for gpos in &fp.copy_gpos {
        for &(col, g) in gpos {
            col_canon[col].get_or_insert(g);
        }
    }
    let mparams = MosaicParams::from_env();
    let mut results = Vec::new();
    let mut read_obs: Vec<Vec<Option<u8>>> = Vec::new(); // per-read PSV observations for the soft EM
    let mut mosaic_calls = Vec::new();
    let mut mosaic_reads = 0usize;
    for (ri, read) in reads.iter().enumerate() {
        let Some(mc) = best_overlap_copy(read, copies) else { continue };
        let feats = read_features(read, mc, &fp);
        // gene-conversion / mosaic: does this molecule's per-PSV copy match SWITCH mid-read?
        let mut site_obs: Vec<SiteObs> = Vec::new();
        for col in 0..fp.n_cols {
            if let (Some(ob), Some(rp)) = (feats.psv_obs[col], col_canon[col]) {
                let match_bits = fp.profiles.iter().map(|pr| pr.alleles[col] == Some(ob)).collect();
                site_obs.push(SiteObs { ref_pos: rp, match_bits });
            }
        }
        site_obs.sort_by_key(|s| s.ref_pos);
        if site_obs.len() > MAX_MOSAIC_SITES {
            let stride = site_obs.len() / MAX_MOSAIC_SITES + 1;
            site_obs = site_obs.into_iter().step_by(stride).collect();
        }
        let mcall = detect_mosaic(&site_obs, copies.len(), MOSAIC_EPS, &mparams);
        if mcall.is_mosaic() {
            mosaic_reads += 1;
        }
        mosaic_calls.push(mcall);
        let Some(combined) = assign_read(&feats, &fp.profiles, p) else { continue };
        let obs = feats.psv_obs.clone();
        read_obs.push(obs.clone());
        let psv_feats = ReadFeatures { psv_obs: feats.psv_obs, junctions: vec![] };
        let Some(psv) = assign_read(&psv_feats, &fp.profiles, p) else { continue };
        results.push(ReadResult { read_index: ri, mapped_copy: mc, psv, combined, psv_obs: obs });
    }
    let conversions = aggregate_family(&mosaic_calls, &mparams);
    // copy-level historical conversions: is any copy's PSV-allele vector a mosaic of two others?
    let copy_conversions = scan_copy_conversions(&fp.profiles, &col_canon, &mparams);
    // soft per-copy abundance (EM) + a normal-approx 95% CI half-width (theta(1-theta)/N).
    let copy_alleles: Vec<Vec<Option<u8>>> = fp.profiles.iter().map(|pr| pr.alleles.clone()).collect();
    let copy_abundance = soft_quantify_em(&read_obs, &copy_alleles, QUANT_ERROR, 100);
    let n = read_obs.len().max(1) as f64;
    let copy_abundance_ci = copy_abundance.iter().map(|&t| 1.96 * (t * (1.0 - t) / n).sqrt()).collect();
    FamilyDetail {
        results,
        n_cols: fp.n_cols,
        copy_abundance,
        copy_abundance_ci,
        mosaic_reads,
        conversions,
        copy_conversions,
        psv_col_pos: col_canon,
        copy_psv_alleles: copy_alleles,
    }
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
        if let Some(mc) = best_overlap_copy(read, copies) {
            if let Some(a) = assign_one_read(read, mc, &fp, p) {
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
    fn junction_resolves_read_when_psvs_cannot() {
        // two copies with IDENTICAL spliced sequence (so ZERO PSV columns) but DIFFERENT junction positions
        // (a copy-specific junction). A read whose intron boundary matches copyA is resolved by the junction
        // alone -- the two-pass `junction_only` case (the DSFAM43 10%->99% effect).
        let spliced = rand_seq(300, 0x5151);
        let ca = copy_tx("A", 0, 400, '+', &[(100, 200)], spliced.clone());
        let cb = copy_tx("B", 1000, 1400, '+', &[(1150, 1250)], spliced.clone());
        let copies = [&ca, &cb];
        // read aligned to copyA's region with copyA's intron structure (boundary at spliced offset 100).
        let read = AlignedRead { ref_start: 0, cigar: vec![('M', 100), ('N', 100), ('M', 200)], seq: spliced };
        let detail = assign_family_detailed(&copies, &[read], &AssignParams::default());
        assert_eq!(detail.n_cols, 0, "identical sequences -> no PSV columns");
        assert_eq!(detail.results.len(), 1);
        let r = &detail.results[0];
        assert_eq!(r.psv.n_decisive, 0, "PSVs alone cannot resolve");
        assert!(r.combined.n_decisive >= 1, "the copy-specific junction resolves it");
        assert_eq!(r.combined.best_copy, 0, "boundary matches copyA's junction");
        assert_eq!(r.combined.status, AssignStatus::Assigned);
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

    /// Ground-truth validation on the sim5x 5-copy synthetic dataset (the python `copy_assign.py sim5x`
    /// oracle): each read name encodes its TRUE copy, so we can report assignment accuracy. K = PSVs per
    /// copy (the identifiability ladder): K=0 has no PSVs → all tied; K>=2 → ~100% accurate. Ignored by
    /// default; run in release with RUSTLE_SIM5X_DIR set.
    #[test]
    #[ignore = "needs sim5x data via RUSTLE_SIM5X_DIR"]
    fn smoke_sim5x_ground_truth() {
        use crate::genome::GenomeIndex;
        use crate::vg_family::denovo_assemble::aligned_reads_from_bam;
        let dir = match std::env::var("RUSTLE_SIM5X_DIR") {
            Ok(d) => d,
            Err(_) => return,
        };
        fn true_copy(name: &str) -> Option<usize> {
            name.split("_c").nth(1)?.split('_').next()?.parse().ok()
        }
        eprintln!(
            "{:>2} {:>7} {:>8} {:>11} {:>13} {:>11} {:>6}",
            "K", "reads", "PSVcols", "resolvable%", "acc|assigned", "acc|argmax", "tied%"
        );
        for k in [0u32, 1, 2, 4] {
            let ref_fa = format!("{dir}/sim5x_K{k}.ref.fa");
            let bam = format!("{dir}/sim5x_K{k}.bam");
            if !std::path::Path::new(&bam).exists() {
                continue;
            }
            let genome = GenomeIndex::from_fasta(&ref_fa).expect("ref");
            let contig = format!("SIM5X_K{k}");
            let clen = genome.chrom_len(&contig);
            let len_g = (clen - 4 * 2000) / 5;
            let unit = len_g + 2000;
            let copies_owned: Vec<DenovoTranscript> = (0..5)
                .map(|c| {
                    let start = c as u64 * unit;
                    let end = start + len_g;
                    let seq = genome.fetch_sequence(&contig, start, end).unwrap();
                    DenovoTranscript {
                        tid: format!("copy{c}"),
                        chrom: contig.clone(),
                        start,
                        end,
                        n_reads: 0,
                        strand: '+',
                        introns: vec![],
                        seq,
                    }
                })
                .collect();
            let copies: Vec<&DenovoTranscript> = copies_owned.iter().collect();
            let fp = build_family_profiles(&copies);
            let reads = aligned_reads_from_bam(&bam, 4).expect("bam");
            let names: Vec<String> = reads.iter().map(|br| br.name.clone()).collect();
            let ars: Vec<AlignedRead> = reads.into_iter().map(|br| br.read).collect();
            let assigns = assign_family(&copies, &ars, &AssignParams::default());
            let (mut n, mut resolvable, mut assigned, mut tied) = (0usize, 0usize, 0usize, 0usize);
            let (mut corr_assigned, mut corr_argmax) = (0usize, 0usize);
            for (ri, a) in &assigns {
                let true_c = match true_copy(&names[*ri]) {
                    Some(t) => t,
                    None => continue,
                };
                n += 1;
                match a.status {
                    AssignStatus::Tied => tied += 1,
                    _ => {
                        resolvable += 1;
                        corr_argmax += (a.best_copy == true_c) as usize;
                        if a.status == AssignStatus::Assigned {
                            assigned += 1;
                            corr_assigned += (a.best_copy == true_c) as usize;
                        }
                    }
                }
            }
            let pct = |x: usize| if n > 0 { 100.0 * x as f64 / n as f64 } else { 0.0 };
            let acc = |c: usize, d: usize| if d > 0 { c as f64 / d as f64 } else { f64::NAN };
            eprintln!(
                "{k:>2} {n:>7} {:>8} {:>10.1}% {:>13.3} {:>11.3} {:>5.1}%",
                fp.n_cols,
                pct(resolvable),
                acc(corr_assigned, assigned),
                acc(corr_argmax, resolvable),
                pct(tied),
            );
        }
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

    #[test]
    fn soft_quantify_recovers_skewed_abundance() {
        // copy0 allele = A, copy1 allele = C at 3 PSV columns. 70 reads carry A, 30 carry C.
        // the EM must recover ~0.7 / 0.3 (the benchmark's skewed case the soft estimator nails).
        let copy_alleles = vec![vec![Some(b'A'); 3], vec![Some(b'C'); 3]];
        let mut reads: Vec<Vec<Option<u8>>> = Vec::new();
        for _ in 0..70 {
            reads.push(vec![Some(b'A'); 3]);
        }
        for _ in 0..30 {
            reads.push(vec![Some(b'C'); 3]);
        }
        let theta = soft_quantify_em(&reads, &copy_alleles, QUANT_ERROR, 100);
        assert!((theta[0] - 0.7).abs() < 0.02, "copy0 ~0.7, got {}", theta[0]);
        assert!((theta[1] - 0.3).abs() < 0.02, "copy1 ~0.3, got {}", theta[1]);
        assert!((theta.iter().sum::<f64>() - 1.0).abs() < 1e-9, "sums to 1");
    }

    #[test]
    fn mosaic_recombinant_read_detected_end_to_end() {
        // two copies differing at 8 PSVs; a RECOMBINANT read carries copy A's alleles at the first 4 PSVs and
        // copy B's at the last 4 -> a gene-conversion switch the wired detector must flag (mosaic_reads >= 1).
        let a = rand_seq(300, 0x0ABC_0001);
        let mut sa = a.clone();
        let mut sb = a.clone();
        let psv = [30usize, 60, 90, 120, 150, 180, 210, 240];
        for &p in &psv {
            sa[p] = b'A';
            sb[p] = b'C';
        }
        let ca = copy_tx("A", 0, 300, '+', &[], sa.clone());
        let cb = copy_tx("B", 1000, 1300, '+', &[], sb);
        let mut recomb = sa; // copy A everywhere...
        for &p in &psv[4..] {
            recomb[p] = b'C'; // ...then switch to copy B at the last 4 PSVs
        }
        let read = AlignedRead { ref_start: 0, cigar: vec![('M', 300)], seq: recomb };
        let detail = assign_family_detailed(&[&ca, &cb], &[read], &AssignParams::default());
        assert!(detail.mosaic_reads >= 1, "recombinant flagged as mosaic; got {}", detail.mosaic_reads);
    }

    #[test]
    fn copy_level_conversion_detected() {
        // 3 copies: A (allele X at all 8 PSVs), B (allele Y), and C = a MOSAIC (matches A over the first 4
        // PSVs, B over the last 4) -- a historical gene conversion baked into copy C's sequence.
        let base = rand_seq(300, 0x0C0C_0001);
        let psv = [30usize, 60, 90, 120, 150, 180, 210, 240];
        let (mut sa, mut sb, mut sc) = (base.clone(), base.clone(), base);
        for (i, &p) in psv.iter().enumerate() {
            sa[p] = b'A';
            sb[p] = b'C';
            sc[p] = if i < 4 { b'A' } else { b'C' };
        }
        let ca = copy_tx("A", 0, 300, '+', &[], sa);
        let cb = copy_tx("B", 1000, 1300, '+', &[], sb);
        let cc = copy_tx("C", 2000, 2300, '+', &[], sc);
        let detail = assign_family_detailed(&[&ca, &cb, &cc], &[], &AssignParams::default());
        assert_eq!(detail.copy_conversions.len(), 1, "exactly copy C is a mosaic of two others");
        assert_eq!(detail.copy_conversions[0].copy_c, 2, "copy C (index 2) is the converted copy");
    }

    #[test]
    fn soft_quantify_no_psv_returns_uniform_prior() {
        // no informative PSV columns (copies identical / reads uncovered) -> the honest identifiability floor:
        // the EM returns the uniform prior, not a false confident split.
        let copy_alleles = vec![vec![None; 2], vec![None; 2]];
        let reads = vec![vec![None; 2]; 40];
        let theta = soft_quantify_em(&reads, &copy_alleles, QUANT_ERROR, 100);
        assert!((theta[0] - 0.5).abs() < 1e-9 && (theta[1] - 0.5).abs() < 1e-9, "uniform: {theta:?}");
    }
}

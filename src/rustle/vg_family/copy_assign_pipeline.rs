//! Per-read COPY ASSIGNMENT driver (integration stage 4a) — the orchestration of `bench/copy_assign.py`.
//!
//! Given a family's copies (`DenovoTranscript`: spliced sequence + exon structure + strand) and the reads
//! over the family region, assign each read — especially the hard multimappers minimap2 leaves at MAPQ 0 —
//! to a specific paralog COPY. Equivalently (the Canzar flip): the copies are PATHS through one family
//! VARIATION GRAPH whose bubbles are the paralog-sequence-variant (PSV) columns, and each read is THREADED
//! through the graph — its PSV bases + copy-specific junctions select its maximum-likelihood copy-path, which
//! the significance gate then accepts or abstains on. "PSV votes ≡ path log-likelihood": the per-copy log-L
//! this driver sums IS the read's score along each copy-path, so votes-vs-threading is one computation, not two.
//!
//! It builds the `CopyProfile`/`ReadFeatures` the already-ported `copy_assign::assign_read` consumes:
//!   1. `discover_psvs` — all-pairs alignment vs copy[0] → columns where copies differ (per-copy genomic
//!      position + transcription-strand base).
//!   2. spliced↔genomic exon map (`exon_map`/`gen2off`) + intron `copy_boundaries`.
//!   3. per read: read its base at each column's genomic position (reverse-complemented for a `-` copy) and
//!      its intron boundaries (mapped to spliced space) → a feature vector → `assign_read`.
//!
//! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use std::collections::{BTreeMap, BTreeSet};

use rayon::prelude::*;

use super::copy_assign::{
    assign_read, assign_read_editing, boundary_present, copy_pair_significance, AssignParams, AssignStatus,
    Assignment, BubbleGraph, CopyProfile, ReadFeatures,
};
use super::copy_split::{intron_chain_of, AlignedRead};
use super::family_detect::DenovoTranscript;
use super::family_graph::poa_msa_with_costs;

/// erf via Abramowitz–Stegun 7.1.26 (|error| < 1.5e-7), enough for a tail test at alpha ~ 1e-3.
fn erf_approx(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * (-x * x).exp();
    sign * y
}

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

/// Pairwise alignment of `other` (query) to `ref_seq` (target) via **minimap2 asm20** (`-c --eqx`),
/// returned in the SAME 2-row gapped-MSA format as `poa_msa_with_costs` (`[gapped_ref, gapped_other]`,
/// `b'-'` = gap) so `discover_psvs` can consume it unchanged. minimap2 (heuristic seed-chain-align) is
/// far faster than poasta (exact Dijkstra DP) on the long (~10 kb) near-identical copy transcripts, and
/// yields the same substitution columns. Opt in with RUSTLE_PSV_MINIMAP2=1.
///
/// Unaligned ends (minimap2 may clip divergent termini) are rendered as separate ref-only / query-only
/// gap blocks so the offset walk in `discover_psvs` stays correct (ro spans `0..ref.len()`, oo spans
/// `0..other.len()`); only the aligned core has both-non-gap columns, so no spurious PSV is created in
/// a clipped end. A reverse-strand or absent alignment returns `Err` → the caller skips that pair (same
/// as a poasta failure).
fn minimap2_msa_pair(ref_seq: &[u8], other: &[u8]) -> anyhow::Result<Vec<Vec<u8>>> {
    use std::io::Write;
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let dir = std::env::temp_dir();
    let pid = std::process::id();
    // A PROCESS-UNIQUE nonce (atomic counter), not a length-derived one: two regions aligning equal-length
    // copies concurrently (region-parallel sweep) would otherwise collide on the same temp path. `pid` keeps
    // distinct processes disjoint. The filename does not affect alignment output, so this is byte-identical.
    use std::sync::atomic::{AtomicUsize, Ordering};
    static PSV_NONCE: AtomicUsize = AtomicUsize::new(0);
    let nonce = PSV_NONCE.fetch_add(1, Ordering::Relaxed);
    let tpath = dir.join(format!("rustle_psv_t_{pid}_{nonce}.fa"));
    let qpath = dir.join(format!("rustle_psv_q_{pid}_{nonce}.fa"));
    struct Cleanup(std::path::PathBuf, std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
            let _ = std::fs::remove_file(&self.1);
        }
    }
    let _cl = Cleanup(tpath.clone(), qpath.clone());
    {
        let mut t = std::fs::File::create(&tpath)?;
        t.write_all(b">t\n")?;
        t.write_all(ref_seq)?;
        t.write_all(b"\n")?;
        let mut q = std::fs::File::create(&qpath)?;
        q.write_all(b">q\n")?;
        q.write_all(other)?;
        q.write_all(b"\n")?;
    }
    let out = std::process::Command::new(&mm2)
        .args(["-c", "--eqx", "-x", "asm20", "-t", "1", "--secondary=no"])
        .arg(&tpath)
        .arg(&qpath)
        .output()?;
    if !out.status.success() {
        anyhow::bail!("minimap2 failed");
    }
    // pick the alignment line with the longest target span (te-ts); forward strand only.
    let text = String::from_utf8_lossy(&out.stdout);
    let mut best: Option<(usize, usize, u64, &str)> = None; // (ts, qs, span, cigar)
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        if f[4] != "+" {
            continue; // reverse-strand: copies are same orientation; skip (caller `continue`s)
        }
        let (qs, ts, te) = (
            f[2].parse::<usize>().unwrap_or(0),
            f[7].parse::<usize>().unwrap_or(0),
            f[8].parse::<usize>().unwrap_or(0),
        );
        let cigar = f.iter().find_map(|x| x.strip_prefix("cg:Z:")).unwrap_or("");
        if cigar.is_empty() {
            continue;
        }
        let span = (te - ts) as u64;
        if best.map_or(true, |b| span > b.2) {
            best = Some((ts, qs, span, cigar));
        }
    }
    let (ts, qs, _span, cigar) = best.ok_or_else(|| anyhow::anyhow!("no forward minimap2 alignment"))?;
    Ok(cigar_to_gapped_msa(ref_seq, other, ts, qs, cigar))
}

/// Reconstruct the 2-row gapped MSA (`[gapped_ref, gapped_other]`, `b'-'` = gap) from a minimap2 PAF
/// alignment: target-start `ts`, query-start `qs`, and the `cg:Z` CIGAR (`=`/`X`/`M` consume both; `I` =
/// query insertion → ref gap; `D` = query deletion → query gap). Unaligned ends become separate ref-only
/// and query-only gap blocks so the walk in `discover_psvs` yields the correct absolute offsets
/// (`ro ∈ 0..ref.len`, `oo ∈ 0..other.len`) with both-non-gap columns only in the aligned core. Pure +
/// binary-free — unit-tested directly.
fn cigar_to_gapped_msa(ref_seq: &[u8], other: &[u8], ts: usize, qs: usize, cigar: &str) -> Vec<Vec<u8>> {
    let mut rrow: Vec<u8> = Vec::with_capacity(ref_seq.len() + other.len());
    let mut orow: Vec<u8> = Vec::with_capacity(ref_seq.len() + other.len());
    // prefix: unaligned ref[0:ts] vs ref-gaps, then unaligned other[0:qs] vs query-gaps (no both-non-gap).
    rrow.extend_from_slice(&ref_seq[..ts]);
    orow.extend(std::iter::repeat(b'-').take(ts));
    rrow.extend(std::iter::repeat(b'-').take(qs));
    orow.extend_from_slice(&other[..qs]);
    let (mut ro, mut oo) = (ts, qs);
    let mut num = 0usize;
    for ch in cigar.bytes() {
        if ch.is_ascii_digit() {
            num = num * 10 + (ch - b'0') as usize;
            continue;
        }
        match ch {
            b'=' | b'X' | b'M' => {
                rrow.extend_from_slice(&ref_seq[ro..ro + num]);
                orow.extend_from_slice(&other[oo..oo + num]);
                ro += num;
                oo += num;
            }
            b'I' => {
                rrow.extend(std::iter::repeat(b'-').take(num));
                orow.extend_from_slice(&other[oo..oo + num]);
                oo += num;
            }
            b'D' => {
                rrow.extend_from_slice(&ref_seq[ro..ro + num]);
                orow.extend(std::iter::repeat(b'-').take(num));
                ro += num;
            }
            _ => {}
        }
        num = 0;
    }
    // suffix: unaligned ref[te:] vs gaps, then unaligned other[qe:] vs gaps.
    rrow.extend_from_slice(&ref_seq[ro..]);
    orow.extend(std::iter::repeat(b'-').take(ref_seq.len() - ro));
    rrow.extend(std::iter::repeat(b'-').take(other.len() - oo));
    orow.extend_from_slice(&other[oo..]);
    vec![rrow, orow]
}

/// Banded global affine-gap (Gotoh) alignment of two sequences, matching poasta's `GapAffine(mismatch=1,
/// gap_extend=1, gap_open=32)` cost (poasta's length-L gap = `open + L·extend`, same as here), restricted to the
/// diagonal band `|i - j| <= band`. Returns the same 2-row gapped MSA (`b'-'` = gap) that `discover_psvs`
/// consumes, or `None` when the band cannot contain the alignment (length difference > `band`, or the traceback
/// touches the band edge) — the caller then falls back to exact poasta.
///
/// **MEASURED (2026-07-12, opt-in `RUSTLE_PSV_BAND`).** A correct direct 3-matrix Gotoh (open-from-match,
/// STATEFUL traceback — the first version's stateless traceback was a bug that fragmented gaps and degraded
/// assignments; fixed) is an OPTIMAL affine-gap aligner and **~54× faster than poasta** (PCDHB 5-copy 27s → 0.5s):
/// the win is avoiding poasta's *graph/A\* machinery*, not the band (the band must be ≥ the copy length
/// difference, so indel-heavy pairs use a wide band and still win). It is **byte-identical to poasta on
/// unambiguous alignments** (verified on a synthetic deletion+SNV pair), and on real families **copy calls /
/// χ_H are UNCHANGED** (PCDHB 5cp both) with per-read assignments differing **< 1%** (PCDHB assigned 6641 vs 6707)
/// — the residual is **co-optimal** gap placement, where poasta's graph traceback and this DP each pick a
/// different equal-cost alignment (neither is "more correct"). Still OPT-IN pending a full known-family
/// re-validation before it could become the default engine.
pub(crate) fn banded_msa_pair(a: &[u8], b: &[u8], band: usize) -> Option<Vec<Vec<u8>>> {
    let (n, m) = (a.len(), b.len());
    if n == 0 || m == 0 || n.abs_diff(m) > band {
        return None; // a truncation/large-indel pair cannot stay in the band — fall back to exact
    }
    const GAP_OPEN: u32 = 32;
    const GAP_EXT: u32 = 1;
    const MISMATCH: u32 = 1;
    const INF: u32 = u32::MAX / 4;
    // Standard 3-matrix Gotoh affine gap (min-cost), banded to |i - j| <= band. mm = ending in a (mis)match at
    // (i,j); e = ending in a gap in `a` (horizontal, consumes b); f = ending in a gap in `b` (vertical, consumes
    // a). Gaps OPEN ONLY from the match state (mm). This matches poasta's gap cost.
    //
    // ⭐ §6em MEMORY: scores live in TWO rolling rows per matrix; the traceback is ONE BYTE per band cell
    // (bits 0-1: mm's predecessor state at (i-1,j-1); bit 2: e extended (1) or opened (0); bit 3: f extended
    // (1) or opened (0)), recorded with the SAME tie-breaks the value-comparison traceback used (mm ≻ f ≻ e for
    // mm; extend ≻ open for e/f), so the alignment is byte-identical to the three-full-matrix version
    // (`banded_msa_pair_full_matrix`, kept under cfg(test) as the oracle). A 25 kb copy against a 4 kb copy
    // (band 22 kb) needed 3 × 25k × 43k × 4 B = 13 GB of scores and took a 25 GB machine to its limit; it now
    // needs 1.1 GB of traceback bytes.
    let lo = |i: usize| i.saturating_sub(band);
    let hi = |i: usize| (i + band).min(m);
    let width = 2 * band + 1;
    let idx = |i: usize, j: usize| -> usize { j - lo(i) };
    let mut tb: Vec<u8> = vec![0; (n + 1) * width];
    let (mut mm_prev, mut e_prev, mut f_prev) = (vec![INF; width], vec![INF; width], vec![INF; width]);
    let (mut mm_cur, mut e_cur, mut f_cur) = (vec![INF; width], vec![INF; width], vec![INF; width]);
    mm_prev[idx(0, 0)] = 0;
    for j in 1..=hi(0) {
        e_prev[idx(0, j)] = GAP_OPEN + (j as u32) * GAP_EXT; // leading gap in a
        tb[j] |= 0b0100; // extended (matches the value traceback: e[0][j-1] + EXT == e[0][j])
    }
    for i in 1..=n {
        let (jl, jh) = (lo(i), hi(i));
        for v in mm_cur.iter_mut() {
            *v = INF;
        }
        for v in e_cur.iter_mut() {
            *v = INF;
        }
        for v in f_cur.iter_mut() {
            *v = INF;
        }
        let row = i * width;
        for j in jl..=jh {
            let mut t: u8 = 0;
            // f: gap in b (vertical) — open from mm[i-1][j], extend from f[i-1][j].
            if j >= lo(i - 1) && j <= hi(i - 1) {
                let o = mm_prev[idx(i - 1, j)].saturating_add(GAP_OPEN + GAP_EXT);
                let x = f_prev[idx(i - 1, j)].saturating_add(GAP_EXT);
                f_cur[idx(i, j)] = o.min(x);
                if x <= o {
                    t |= 0b1000;
                }
            }
            if j == 0 {
                f_cur[idx(i, 0)] = GAP_OPEN + (i as u32) * GAP_EXT; // leading gap in b
                t |= 0b1000;
                tb[row + idx(i, 0)] = t;
                continue;
            }
            // e: gap in a (horizontal) — open from mm[i][j-1], extend from e[i][j-1].
            if j - 1 >= jl {
                let o = mm_cur[idx(i, j - 1)].saturating_add(GAP_OPEN + GAP_EXT);
                let x = e_cur[idx(i, j - 1)].saturating_add(GAP_EXT);
                e_cur[idx(i, j)] = o.min(x);
                if x <= o {
                    t |= 0b0100;
                }
            }
            // mm: (mis)match — best of the three states at (i-1, j-1) plus the substitution cost.
            if j - 1 >= lo(i - 1) && j - 1 <= hi(i - 1) {
                let sub = if a[i - 1] == b[j - 1] { 0 } else { MISMATCH };
                let (pm, pe, pf) =
                    (mm_prev[idx(i - 1, j - 1)], e_prev[idx(i - 1, j - 1)], f_prev[idx(i - 1, j - 1)]);
                let prev = pm.min(pe).min(pf);
                mm_cur[idx(i, j)] = prev.saturating_add(sub);
                t |= if pm == prev {
                    0
                } else if pf == prev {
                    2
                } else {
                    1
                };
            }
            tb[row + idx(i, j)] = t;
        }
        std::mem::swap(&mut mm_prev, &mut mm_cur);
        std::mem::swap(&mut e_prev, &mut e_cur);
        std::mem::swap(&mut f_prev, &mut f_cur);
    }
    if m < lo(n) || m > hi(n) {
        return None;
    }
    let end_mm = mm_prev[idx(n, m)];
    let end_e = e_prev[idx(n, m)];
    let end_f = f_prev[idx(n, m)];
    let best = end_mm.min(end_e).min(end_f);
    if best >= INF {
        return None;
    }
    // STATEFUL traceback over the recorded predecessors (0=mm, 1=e, 2=f).
    let (mut i, mut j) = (n, m);
    let mut state: u8 = if best == end_mm {
        0
    } else if best == end_f {
        2
    } else {
        1
    };
    let (mut ra, mut rb): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
    let mut hit_edge = false;
    while i > 0 || j > 0 {
        if i.abs_diff(j) >= band {
            hit_edge = true; // path reached the band boundary — the true optimum may lie outside
        }
        let t = tb[i * width + idx(i, j)];
        match state {
            0 => {
                ra.push(a[i - 1]);
                rb.push(b[j - 1]);
                i -= 1;
                j -= 1;
                state = t & 0b11;
            }
            2 => {
                ra.push(a[i - 1]);
                rb.push(b'-');
                i -= 1;
                state = if t & 0b1000 != 0 { 2 } else { 0 };
            }
            _ => {
                ra.push(b'-');
                rb.push(b[j - 1]);
                j -= 1;
                state = if t & 0b0100 != 0 { 1 } else { 0 };
            }
        }
    }
    if hit_edge {
        return None; // fall back to exact — the band clipped the alignment
    }
    ra.reverse();
    rb.reverse();
    Some(vec![ra, rb])
}

#[cfg(test)]
pub(crate) fn banded_msa_pair_full_matrix(a: &[u8], b: &[u8], band: usize) -> Option<Vec<Vec<u8>>> {
    let (n, m) = (a.len(), b.len());
    if n == 0 || m == 0 || n.abs_diff(m) > band {
        return None; // a truncation/large-indel pair cannot stay in the band — fall back to exact
    }
    const GAP_OPEN: u32 = 32;
    const GAP_EXT: u32 = 1;
    const MISMATCH: u32 = 1;
    const INF: u32 = u32::MAX / 4;
    // Standard 3-matrix Gotoh affine gap (min-cost), banded to |i - j| <= band. mm = ending in a (mis)match at
    // (i,j); e = ending in a gap in `a` (horizontal, consumes b); f = ending in a gap in `b` (vertical, consumes
    // a). Gaps OPEN ONLY from the match state (mm) — the classic formulation that yields clean, non-fragmented
    // gaps (a length-L gap = open + L·extend, one open). This matches poasta's gap cost.
    let lo = |i: usize| i.saturating_sub(band);
    let hi = |i: usize| (i + band).min(m);
    let width = 2 * band + 1;
    let idx = |i: usize, j: usize| -> usize { j - lo(i) };
    let mut mm = vec![vec![INF; width]; n + 1];
    let mut e = vec![vec![INF; width]; n + 1];
    let mut f = vec![vec![INF; width]; n + 1];
    mm[0][idx(0, 0)] = 0;
    for j in 1..=hi(0) {
        e[0][idx(0, j)] = GAP_OPEN + (j as u32) * GAP_EXT; // leading gap in a
    }
    for i in 1..=n {
        let (jl, jh) = (lo(i), hi(i));
        for j in jl..=jh {
            // f: gap in b (vertical) — open from mm[i-1][j], extend from f[i-1][j].
            if j >= lo(i - 1) && j <= hi(i - 1) {
                let o = mm[i - 1][idx(i - 1, j)].saturating_add(GAP_OPEN + GAP_EXT);
                let x = f[i - 1][idx(i - 1, j)].saturating_add(GAP_EXT);
                f[i][idx(i, j)] = o.min(x);
            }
            if j == 0 {
                f[i][idx(i, 0)] = GAP_OPEN + (i as u32) * GAP_EXT; // leading gap in b
                continue;
            }
            // e: gap in a (horizontal) — open from mm[i][j-1], extend from e[i][j-1].
            if j - 1 >= jl {
                let o = mm[i][idx(i, j - 1)].saturating_add(GAP_OPEN + GAP_EXT);
                let x = e[i][idx(i, j - 1)].saturating_add(GAP_EXT);
                e[i][idx(i, j)] = o.min(x);
            }
            // mm: (mis)match — best of the three states at (i-1, j-1) plus the substitution cost.
            if j - 1 >= lo(i - 1) && j - 1 <= hi(i - 1) {
                let sub = if a[i - 1] == b[j - 1] { 0 } else { MISMATCH };
                let prev = mm[i - 1][idx(i - 1, j - 1)]
                    .min(e[i - 1][idx(i - 1, j - 1)])
                    .min(f[i - 1][idx(i - 1, j - 1)]);
                mm[i][idx(i, j)] = prev.saturating_add(sub);
            }
        }
    }
    if m < lo(n) || m > hi(n) {
        return None;
    }
    let end_mm = mm[n][idx(n, m)];
    let end_e = e[n][idx(n, m)];
    let end_f = f[n][idx(n, m)];
    let best = end_mm.min(end_e).min(end_f);
    if best >= INF {
        return None;
    }
    // STATEFUL traceback: carry the current matrix (0=mm, 1=e, 2=f) so a gap run continues as one gap.
    let (mut i, mut j) = (n, m);
    let mut state: u8 = if best == end_mm {
        0
    } else if best == end_f {
        2
    } else {
        1
    };
    let (mut ra, mut rb): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
    let mut hit_edge = false;
    let inband = |i: usize, j: usize| j >= lo(i) && j <= hi(i);
    while i > 0 || j > 0 {
        if i.abs_diff(j) >= band {
            hit_edge = true; // path reached the band boundary — the true optimum may lie outside
        }
        match state {
            0 => {
                // match/mismatch: consume both; decide predecessor state at (i-1, j-1).
                ra.push(a[i - 1]);
                rb.push(b[j - 1]);
                let sub = if a[i - 1] == b[j - 1] { 0 } else { MISMATCH };
                let target = mm[i][idx(i, j)].wrapping_sub(sub);
                i -= 1;
                j -= 1;
                state = if inband(i, j) && mm[i][idx(i, j)] == target {
                    0
                } else if inband(i, j) && f[i][idx(i, j)] == target {
                    2
                } else {
                    1
                };
            }
            2 => {
                // gap in b (vertical): consume a. continue f, or opened from mm[i-1][j].
                ra.push(a[i - 1]);
                rb.push(b'-');
                let cur = f[i][idx(i, j)];
                let from_ext = inband(i - 1, j) && f[i - 1][idx(i - 1, j)].saturating_add(GAP_EXT) == cur;
                i -= 1;
                state = if from_ext { 2 } else { 0 };
            }
            _ => {
                // gap in a (horizontal): consume b. continue e, or opened from mm[i][j-1].
                ra.push(b'-');
                rb.push(b[j - 1]);
                let cur = e[i][idx(i, j)];
                let from_ext = j >= 1 && inband(i, j - 1) && e[i][idx(i, j - 1)].saturating_add(GAP_EXT) == cur;
                j -= 1;
                state = if from_ext { 1 } else { 0 };
            }
        }
    }
    if hit_edge {
        return None; // fall back to exact — the band clipped the alignment
    }
    ra.reverse();
    rb.reverse();
    Some(vec![ra, rb])
}

/// ⭐ O2-9 / D3 — "THE READ IS THE STAR" (PREREG adj/d3; §6fa). One minimap2 batch per family: every
/// molecule's sequence (query) against every copy's spliced unit (target), all hits kept. For each read × copy
/// the hit's CIGAR (`--eqx`) gives the copy base aligned to each READ position, in the read's own
/// orientation (complemented for a `-` hit). Returns, per read, per copy, `None` (no hit) or a vector over
/// read positions of `Option<u8>` (the copy's aligned base, `None` where the read has an insertion or the
/// hit does not cover the position). The star-projected family columns (`discover_psvs`, copy[0] as the
/// star) misplace positions on copies whose spliced units differ in exon content (NPIP: 34 % of columns
/// disagree between two placements of the same read, §6fa); here every read carries its own columns.
pub fn read_star_alignments(copy_seqs: &[&[u8]], read_seqs: &[&[u8]]) -> Vec<Vec<Option<Vec<Option<u8>>>>> {
    // full (memory-hungry) form, kept for tests: per read x copy the aligned base at every read position
    let mut out: Vec<Vec<Option<Vec<Option<u8>>>>> = (0..read_seqs.len()).map(|_| vec![None; copy_seqs.len()]).collect();
    read_star_stream(copy_seqs, read_seqs, |ri, per_copy, _, _| out[ri] = per_copy.to_vec());
    out
}

/// One observation per molecule, streamed: the PAF is consumed query by query and each read is reduced to
/// `(candidate copies, read bases at its columns, one CopyProfile per candidate)` before the next read's
/// hits are parsed. Peak memory = one read's hits. `copy_id` in the profiles is the family copy index.
/// Per read: `(candidates, obs, profiles, unit_edits)` — `unit_edits[ci] = Some((edits, aligned))` of the hit
/// used for copy `ci` (edits = minimap2 `NM`, aligned = alignment block length), the input to the
/// genome-versus-unit certificate in the assigner.
/// `(candidates, obs, profiles, unit_edits, covered)` — `covered[k]` = the read positions candidate `k` has an
/// aligned base at, as sorted half-open intervals (a junction term is admitted only where BOTH candidates cover
/// the read position; an uncovered position is not "the copy lacks the junction").
pub type ReadStarObs = (Vec<usize>, Vec<Option<u8>>, Vec<CopyProfile>, Vec<Option<(u64, u64, u64, u64)>>, Vec<Vec<(u32, u32)>>);

/// `copy_bounds[c]` = the unit offsets where copy `c` starts a new exon (its splice boundaries in unit space);
/// each candidate's `CopyProfile::junctions` becomes the READ positions aligned to those offsets, so the read's
/// own junction positions (from its CIGAR) can be compared per candidate in one coordinate system (§6fc).
pub fn read_star_observations(copy_seqs: &[&[u8]], read_seqs: &[&[u8]], copy_bounds: &[Vec<u32>]) -> Vec<ReadStarObs> {
    let mut out: Vec<ReadStarObs> =
        (0..read_seqs.len()).map(|_| (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new())).collect();
    read_star_stream(copy_seqs, read_seqs, |ri, per_copy, per_pos, edits| {
        let (obs, profiles) = read_star_columns(read_seqs[ri], per_copy);
        let cand: Vec<usize> = per_copy.iter().enumerate().filter(|(_, h)| h.is_some()).map(|(i, _)| i).collect();
        // every column stays; a candidate with no base at a column is `None` there and the PAIRWISE certificate
        // (§6fc) scores each pair on the columns both carry — no intersection over all candidates (that rule made
        // one 13-read admitted copy turn 1,981 assignments into 233, and `-p 0` halve the assignments, row 705).
        let prof2: Vec<CopyProfile> = cand
            .iter()
            .map(|&ci| {
                // read positions aligned to this copy's exon boundaries
                let bounds: std::collections::BTreeSet<u32> = copy_bounds.get(ci).map(|v| v.iter().copied().collect()).unwrap_or_default();
                let junctions: Vec<i64> = per_pos[ci]
                    .as_ref()
                    .map(|pv| pv.iter().enumerate().filter_map(|(r, o)| o.filter(|o| bounds.contains(o)).map(|_| r as i64)).collect())
                    .unwrap_or_default();
                CopyProfile { copy_id: ci, alleles: profiles[ci].alleles.clone(), junctions }
            })
            .collect();
        let covered: Vec<Vec<(u32, u32)>> = cand
            .iter()
            .map(|&ci| {
                let mut iv: Vec<(u32, u32)> = Vec::new();
                if let Some(pv) = per_copy[ci].as_ref() {
                    for (r, b) in pv.iter().enumerate() {
                        if b.is_some() {
                            match iv.last_mut() {
                                Some(last) if last.1 == r as u32 => last.1 = r as u32 + 1,
                                _ => iv.push((r as u32, r as u32 + 1)),
                            }
                        }
                    }
                }
                iv
            })
            .collect();
        out[ri] = (cand, obs, prof2, edits.to_vec(), covered);
    });
    out
}

/// Shared driver: runs minimap2 (all copies as targets, all molecules as queries, every hit kept) and calls
/// `sink(read_index, per_copy)` once per read with that read's aligned copy bases per read position.
fn read_star_stream<F: FnMut(usize, &[Option<Vec<Option<u8>>>], &[Option<Vec<Option<u32>>>], &[Option<(u64, u64, u64, u64)>])>(copy_seqs: &[&[u8]], read_seqs: &[&[u8]], mut sink: F) {
    use std::io::{BufRead, Write};
    use std::sync::atomic::{AtomicUsize, Ordering};
    static NONCE: AtomicUsize = AtomicUsize::new(0);
    let (nc, nr) = (copy_seqs.len(), read_seqs.len());
    if nc == 0 || nr == 0 {
        return;
    }
    let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".to_string());
    let dir = std::env::temp_dir();
    let nonce = NONCE.fetch_add(1, Ordering::Relaxed);
    let refp = dir.join(format!("rustle_star_ref_{}_{nonce}.fa", std::process::id()));
    let qp = dir.join(format!("rustle_star_q_{}_{nonce}.fa", std::process::id()));
    struct Cleanup(std::path::PathBuf);
    impl Drop for Cleanup {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }
    let _c1 = Cleanup(refp.clone());
    let _c2 = Cleanup(qp.clone());
    let write_fa = |path: &std::path::Path, seqs: &[&[u8]]| -> Option<()> {
        let mut fh = std::io::BufWriter::new(std::fs::File::create(path).ok()?);
        for (i, s) in seqs.iter().enumerate() {
            writeln!(fh, ">{i}").ok()?;
            fh.write_all(s).ok()?;
            fh.write_all(b"\n").ok()?;
        }
        fh.flush().ok()
    };
    if write_fa(&refp, copy_seqs).is_none() || write_fa(&qp, read_seqs).is_none() {
        return;
    }
    let n_arg = nc.to_string();
    // reporting floor for a copy's hit relative to the read's best hit (minimap2 -p). 0.3 by default; a copy
    // below it is not a candidate. `RUSTLE_STAR_P=0` reports every chain (the test of that floor, §6fc).
    let p_arg = std::env::var("RUSTLE_STAR_P").unwrap_or_else(|_| "0.3".to_string());
    let mut child = match std::process::Command::new(&mm2)
        .args(["-c", "--eqx", "-x", "map-hifi", "--secondary=yes", "-p", &p_arg, "-N", &n_arg, "-t", "2"])
        .arg(&refp)
        .arg(&qp)
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::null())
        .spawn()
    {
        Ok(c) => c,
        Err(_) => return,
    };
    let stdout = child.stdout.take().expect("piped stdout");
    let reader = std::io::BufReader::new(stdout);
    // per-query accumulation: minimap2 emits a query's hits consecutively
    let mut cur: Option<usize> = None;
    let mut per_copy: Vec<Option<Vec<Option<u8>>>> = vec![None; nc];
    let mut per_pos: Vec<Option<Vec<Option<u32>>>> = vec![None; nc];
    let mut best_match: Vec<u64> = vec![0; nc];
    let mut edits: Vec<Option<(u64, u64, u64, u64)>> = vec![None; nc];
    let mut flush = |ri: Option<usize>, per_copy: &mut Vec<Option<Vec<Option<u8>>>>, per_pos: &mut Vec<Option<Vec<Option<u32>>>>, best_match: &mut Vec<u64>, edits: &mut Vec<Option<(u64, u64, u64, u64)>>| {
        if let Some(ri) = ri {
            sink(ri, per_copy, per_pos, edits);
        }
        for (((h, pp), b), e) in per_copy.iter_mut().zip(per_pos.iter_mut()).zip(best_match.iter_mut()).zip(edits.iter_mut()) {
            *h = None;
            *pp = None;
            *b = 0;
            *e = None;
        }
    };
    for line in reader.lines().map_while(Result::ok) {
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 12 {
            continue;
        }
        let (Ok(ri), Ok(qlen), Ok(qs), Ok(qe), Ok(ci), Ok(ts), Ok(nmatch), Ok(blk)) = (
            f[0].parse::<usize>(),
            f[1].parse::<usize>(),
            f[2].parse::<usize>(),
            f[3].parse::<usize>(),
            f[5].parse::<usize>(),
            f[7].parse::<usize>(),
            f[9].parse::<u64>(),
            f[10].parse::<u64>(),
        ) else {
            continue;
        };
        if ri >= nr || ci >= nc {
            continue;
        }
        if cur != Some(ri) {
            flush(cur, &mut per_copy, &mut per_pos, &mut best_match, &mut edits);
            cur = Some(ri);
        }
        if nmatch <= best_match[ci] {
            continue; // keep the hit with the most MATCHING bases per read x copy (not the longest block)
        }
        let Some(cg) = f[12..].iter().find_map(|t| t.strip_prefix("cg:Z:")) else { continue };
        let nm: u64 = f[12..].iter().find_map(|t| t.strip_prefix("NM:i:")).and_then(|v| v.parse().ok()).unwrap_or(blk - nmatch);
        let minus = f[4] == "-";
        let copy = copy_seqs[ci];
        let mut v: Vec<Option<u8>> = vec![None; qlen];
        let mut vp: Vec<Option<u32>> = vec![None; qlen];
        let (mut q, mut t) = (if minus { qlen - qe } else { qs }, ts);
        let mut num = 0usize;
        let (mut n_x, mut n_aligned) = (0u64, 0u64); // substitutions and aligned bases: the ORIGIN certificate's
                                                    // numerator/denominator (indels are isoform structure, not origin, §6fc)
        for ch in cg.bytes() {
            if ch.is_ascii_digit() {
                num = num * 10 + (ch - b'0') as usize;
                continue;
            }
            match ch {
                b'=' | b'X' | b'M' => {
                    if ch == b'X' {
                        n_x += num as u64;
                    }
                    n_aligned += num as u64;
                    for k in 0..num {
                        let (qq, tt) = (q + k, t + k);
                        if tt < copy.len() && qq < qlen {
                            let b = copy[tt].to_ascii_uppercase();
                            let (pos, base) = if minus { (qlen - 1 - qq, rc_base(b)) } else { (qq, b) };
                            v[pos] = Some(base);
                            vp[pos] = Some(tt as u32);
                        }
                    }
                    q += num;
                    t += num;
                }
                b'I' | b'S' => q += num,
                b'D' | b'N' => t += num,
                _ => {}
            }
            num = 0;
        }
        best_match[ci] = nmatch;
        edits[ci] = Some((n_x, n_aligned, nm, blk));
        per_copy[ci] = Some(v);
        per_pos[ci] = Some(vp);
    }
    flush(cur, &mut per_copy, &mut per_pos, &mut best_match, &mut edits);
    let _ = child.wait();
}

/// Per-read columns from `read_star_alignments`: read positions where ≥ 2 copies carry an aligned base and
/// not all of them agree. Returns the read's own bases at those positions (the observation, always `Some`)
/// and one `CopyProfile` per copy (its aligned base per column, `None` where it has no base).
pub fn read_star_columns(read: &[u8], per_copy: &[Option<Vec<Option<u8>>>]) -> (Vec<Option<u8>>, Vec<CopyProfile>) {
    let nc = per_copy.len();
    let mut cols: Vec<usize> = Vec::new();
    for pos in 0..read.len() {
        let mut seen: Option<u8> = None;
        let (mut n, mut differ) = (0usize, false);
        for pc in per_copy.iter().flatten() {
            if let Some(b) = pc.get(pos).copied().flatten() {
                n += 1;
                match seen {
                    None => seen = Some(b),
                    Some(x) if x != b => differ = true,
                    _ => {}
                }
            }
        }
        if n >= 2 && differ {
            cols.push(pos);
        }
    }
    let obs: Vec<Option<u8>> = cols.iter().map(|&p| Some(read[p].to_ascii_uppercase())).collect();
    let profiles = (0..nc)
        .map(|ci| CopyProfile {
            copy_id: ci,
            alleles: cols.iter().map(|&p| per_copy[ci].as_ref().and_then(|pc| pc.get(p).copied().flatten())).collect(),
            junctions: Vec::new(),
        })
        .collect();
    (obs, profiles)
}

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
    // PSV-column alignment engine. DEFAULT: our own direct banded Gotoh affine-gap DP (`banded_msa_pair`) — an
    // OPTIMAL pairwise aligner that is ~54x faster than poasta's graph/A* machinery on the 2-sequence case,
    // byte-identical on unambiguous alignments and differing only on co-optimal gap placement (<1% of reads, 0
    // change to copy calls). The band is ADAPTIVE per pair (|len_a - len_b| + PSV_BAND_MARGIN), so net-length
    // differences never overflow it; the rare large internal indel that exceeds the band falls back to EXACT
    // poasta, so correctness is preserved regardless of the margin. Escapes: RUSTLE_PSV_POASTA=1 forces exact
    // poasta (reproduces the pre-2026-07-12 engine); RUSTLE_PSV_BAND=<w> pins a FIXED absolute band;
    // RUSTLE_PSV_MINIMAP2=1 uses the minimap2 asm20 heuristic (fast but loses divergent-flank PSVs). All three
    // return the same 2-row gapped-MSA format.
    let use_mm2 = std::env::var_os("RUSTLE_PSV_MINIMAP2").is_some();
    let force_poasta = std::env::var_os("RUSTLE_PSV_POASTA").is_some();
    let fixed_band: Option<usize> =
        std::env::var("RUSTLE_PSV_BAND").ok().and_then(|v| v.parse::<usize>().ok()).filter(|&w| w > 0);
    const PSV_BAND_MARGIN: usize = 1024; // slack over |len_a-len_b| for internal (length-canceling) indels
    // Align each non-ref copy against the ref (the per-pair poasta DP is the dominant per-family cost — an
    // exact O(len^2) alignment of two long divergent transcripts). The (n-1) alignments are INDEPENDENT, so
    // run them concurrently and merge: each yields its own `amap` (ref_off -> other_off) plus the ref offsets
    // where the two bases DIFFER. The merge is order-independent (per-index amaps + a set-union of diffs), so
    // the result is byte-identical to the serial walk. poasta is a pure function (thread-safe); the opt-in
    // minimap2 path also stays serial HERE (kept simple), though its temp files are now process-unique
    // (atomic nonce) so it is safe under the region-parallel sweep.
    let align_one = |other: usize| -> (BTreeMap<usize, usize>, BTreeSet<usize>) {
        // strong gap-open anchors the conserved core column-for-column (same config as contiguous_core_coverage).
        let banded = if use_mm2 || force_poasta {
            None
        } else {
            let band =
                fixed_band.unwrap_or_else(|| ref_seq.len().abs_diff(copies[other].seq.len()) + PSV_BAND_MARGIN);
            banded_msa_pair(ref_seq, &copies[other].seq, band)
        };
        let aln = if use_mm2 {
            minimap2_msa_pair(ref_seq, &copies[other].seq)
        } else if let Some(band) = banded {
            Ok(band) // banded DP succeeded within the band
        } else {
            // forced poasta, or the band could not contain the alignment -> exact poasta (accuracy preserved)
            poa_msa_with_costs(&[ref_seq.clone(), copies[other].seq.clone()], GapAffine::new(1, 1, 32))
        };
        let mut amap = BTreeMap::new();
        let mut diffs = BTreeSet::new();
        if let Ok(msa) = aln {
            if msa.len() == 2 {
                let (r0, r1) = (&msa[0], &msa[1]);
                let (mut ro, mut oo) = (0usize, 0usize);
                for c in 0..r0.len().min(r1.len()) {
                    let (ca, cb) = (r0[c], r1[c]);
                    let (a_gap, b_gap) = (ca == b'-', cb == b'-');
                    if !a_gap && !b_gap {
                        amap.insert(ro, oo);
                        if is_acgt(ca) && is_acgt(cb) && ca != cb {
                            diffs.insert(ro);
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
        }
        (amap, diffs)
    };
    let per_copy: Vec<(BTreeMap<usize, usize>, BTreeSet<usize>)> = if use_mm2 {
        (1..n).map(&align_one).collect()
    } else {
        (1..n).into_par_iter().map(&align_one).collect()
    };
    let mut amaps: Vec<BTreeMap<usize, usize>> = vec![BTreeMap::new(); n]; // index 0 unused
    let mut diff_off: BTreeSet<usize> = BTreeSet::new();
    for (other, (amap, diffs)) in (1..n).zip(per_copy.into_iter()) {
        amaps[other] = amap;
        diff_off.extend(diffs);
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

/// INTRON-RETENTION lever (opt-in, `RUSTLE_INTRON_PSV=1`): discover PSV columns that fall in INTRONS, where
/// copies that are exon-identical may still differ. A read that RETAINS an intron carries M-aligned bases at
/// these genomic positions (a spliced read has none — it is `None`, the same as an unspanned column), so the
/// existing per-read CIGAR sweep fills them automatically and the significance gate uses them as extra
/// distinguishing columns. This is the last RNA-intrinsic lever for the K=0 (exon-identical) Tied reads.
///
/// Aligns each copy's FORWARD genomic span (exons + introns, fetched from `genome`) vs copy[0], keeps the
/// substitution columns whose copy[0] position is intronic, and emits `PsvColumn`s with forward-genome
/// positions and alleles in transcription orientation (complemented for a `-` copy, matching the read-base
/// reader's reverse-complement). Returns empty if the genome lookup fails or a span is too long for the aligner.
pub fn discover_intron_psvs(
    copies: &[&DenovoTranscript],
    genome: &crate::genome::GenomeIndex,
) -> Vec<PsvColumn> {
    discover_genomic_psvs(copies, genome, true)
}

/// ⭐ O2-8c (§6eo): PSV columns from the GENOMIC alignment of the copies' spans. With `intronic_only = false`
/// this REPLACES the spliced-sequence discovery: units whose exon composition differs (read-limited chains,
/// 89 bp to 19 kb) sent the spliced star projection to thousands of non-homologous "PSV" columns and min_p
/// 3e-270 on a wrong call (register 683); the genomic spans of SEDEF core hulls are collinear, so the
/// pairwise alignment is a thin band. If the forward alignment covers less than half of the shorter span the
/// reverse complement is tried (an inverted duplication), with positions mapped back accordingly.
/// ⭐ O2-8c′ (§6ep): per-copy SEDEF core hulls, registered by `catalog_input` from the optional `core_hull`
/// column and read by `discover_genomic_psvs`. A process-wide side table keyed by copy `tid`: adding a field
/// to `DenovoTranscript` would touch 81 constructors for an opt-in lever (same rationale as the env-var
/// levers). Only `--psv-genomic` consults it.
static CORE_HULLS: std::sync::OnceLock<std::sync::Mutex<std::collections::HashMap<String, (u64, u64)>>> =
    std::sync::OnceLock::new();
pub fn register_core_hull(tid: &str, hull: (u64, u64)) {
    CORE_HULLS.get_or_init(Default::default).lock().unwrap().insert(tid.to_string(), hull);
}
pub fn core_hull_of(tid: &str) -> Option<(u64, u64)> {
    CORE_HULLS.get().and_then(|m| m.lock().unwrap().get(tid).copied())
}

pub fn discover_genomic_psvs(
    copies: &[&DenovoTranscript],
    genome: &crate::genome::GenomeIndex,
    intronic_only: bool,
) -> Vec<PsvColumn> {
    let n = copies.len();
    if n < 2 {
        return Vec::new();
    }
    // Alignment span per copy: its registered SEDEF core hull (coextensive across the family by construction,
    // §6ei) when present, else its extent. §6eo showed extents manufacture columns from non-homologous flanks.
    let span_of: Vec<(u64, u64)> = copies
        .iter()
        .map(|c| core_hull_of(&c.tid).filter(|&(a, b)| b > a).unwrap_or((c.start, c.end)))
        .collect();
    let comp = |b: u8| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        x => x,
    };
    // Genomic spans (exons + introns) in TRANSCRIPTION orientation — the same allele space as the spliced
    // sequences and as the read observations (`fill_psv_obs` complements a read's base when its mapped copy is
    // '-'). A '-' copy's span is reverse-complemented, and its offset `o` maps back to genome `end - 1 - o`.
    let spans: Vec<Vec<u8>> = copies
        .iter()
        .zip(&span_of)
        .map(|(c, &(s0, e0))| {
            let fwd = genome
                .fetch_sequence(&c.chrom, s0, e0)
                .map(|s| s.to_ascii_uppercase())
                .unwrap_or_default();
            if c.strand == '-' { fwd.iter().rev().map(|&b| comp(b)).collect() } else { fwd }
        })
        .collect();
    if spans[0].is_empty() {
        return Vec::new();
    }
    let c0 = copies[0];
    let ref_span = &spans[0];
    let gpos = |ci: usize, off: usize| -> u64 {
        let (s0, e0) = span_of[ci];
        if copies[ci].strand == '-' { e0 - 1 - off as u64 } else { s0 + off as u64 }
    };
    let is_intronic = |gpos: u64| c0.introns.iter().any(|&(d, a)| gpos >= d && gpos < a);
    use poasta::aligner::scoring::GapAffine;
    // poasta is exact but O(n^2); above this, fall back to minimap2. Overridable via `RUSTLE_POA_CAP` (set from
    // `copy_assign --poa-cap`, default 20000). Read via env var rather than a fn parameter (see the rationale at
    // `RUSTLE_SKIP_POA_DIAGNOSTIC`).
    let poa_cap: usize = std::env::var("RUSTLE_POA_CAP").ok().and_then(|s| s.parse().ok()).unwrap_or(20_000);
    // The banded DP's traceback is one byte per band cell; refuse it above ~1.5 GB (a 225 kb span against a
    // 25 kb one asked for 77 GB, §6eo) and let poasta/minimap2 take the pair.
    const BAND_CELL_BUDGET: usize = 1_500_000_000;
    let mut amaps: Vec<BTreeMap<usize, usize>> = vec![BTreeMap::new(); n];
    let mut diff_off: BTreeSet<usize> = BTreeSet::new();
    for other in 1..n {
        if spans[other].is_empty() {
            continue;
        }
        let os = &spans[other];
        // §6ep: genomic hulls are NOT collinear end-to-end even between two copies of one family — a 16 kb
        // core at 96% identity sits beside flanking segments at 77% (repeat mosaics) or unshared. A GLOBAL
        // DP forces those into the alignment and manufactures thousands of "PSV" columns (register 684/685).
        // minimap2's chained alignment keeps the homologous block and renders the rest as gap-only ends, which
        // create no column; the exact aligners are the fallback for pairs minimap2 cannot seed (short spans).
        let aln = match minimap2_msa_pair(ref_span, os) {
            Ok(m) => Ok(m),
            Err(_) => {
                let band = ref_span.len().abs_diff(os.len()) + 1024;
                let banded = if (ref_span.len() + 1).saturating_mul(2 * band + 1) <= BAND_CELL_BUDGET {
                    banded_msa_pair(ref_span, os, band)
                } else {
                    None
                };
                if let Some(m) = banded {
                    Ok(m)
                } else if ref_span.len().max(os.len()) <= poa_cap {
                    poa_msa_with_costs(&[ref_span.clone(), os.clone()], GapAffine::new(1, 1, 32))
                } else {
                    Err(anyhow::anyhow!("no aligner could take this pair"))
                }
            }
        };
        if let Ok(msa) = aln {
            if msa.len() == 2 {
                let (r0, r1) = (&msa[0], &msa[1]);
                let (mut ro, mut oo) = (0usize, 0usize);
                for c in 0..r0.len().min(r1.len()) {
                    let (ca, cb) = (r0[c], r1[c]);
                    let (a_gap, b_gap) = (ca == b'-', cb == b'-');
                    if !a_gap && !b_gap {
                        amaps[other].insert(ro, oo);
                        if is_acgt(ca) && is_acgt(cb) && ca != cb && (!intronic_only || is_intronic(gpos(0, ro))) {
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
        }
    }
    diff_off
        .into_iter()
        .map(|ro| {
            let mut rec: PsvColumn = vec![None; n];
            rec[0] = Some((gpos(0, ro), ref_span[ro]));
            for other in 1..n {
                if let Some(&oo) = amaps[other].get(&ro) {
                    rec[other] = Some((gpos(other, oo), spans[other][oo]));
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
pub fn build_family_profiles(
    copies: &[&DenovoTranscript],
    genome: Option<&crate::genome::GenomeIndex>,
) -> FamilyProfiles {
    let exon_maps: Vec<Vec<u64>> = copies.iter().map(|c| exon_map(c)).collect();
    // ⭐ O2-8c: `RUSTLE_PSV_GENOMIC=1` (set by `copy_assign --psv-genomic`) discovers the columns on the GENOMIC
    // alignment of the copies' spans instead of their spliced sequences (register 683). Unset => unchanged.
    let genomic_mode = genome.is_some() && std::env::var_os("RUSTLE_PSV_GENOMIC").is_some();
    let mut cols = if let (true, Some(g)) = (genomic_mode, genome) {
        discover_genomic_psvs(copies, g, false)
    } else {
        discover_psvs(copies, &exon_maps)
    };
    // INTRON-RETENTION lever (opt-in): append intronic PSV columns so reads that retain an intron can use the
    // intronic sequence. OFF (env unset / no genome) => `cols` is unchanged => byte-identical to the exon-only gate.
    if let (false, Some(g)) = (genomic_mode, genome) {
        if std::env::var_os("RUSTLE_INTRON_PSV").is_some() {
            cols.extend(discover_intron_psvs(copies, g));
        }
    }
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

/// Minimum reads per allele for a column to count as a read-supported PSV (see `read_supported_columns`).
pub const PSV_MIN_ALLELE_READS: usize = 2;
/// Minimum coverage before a candidate column is JUDGED at all. A column with fewer reads than this is kept
/// (we cannot invalidate what we cannot see — e.g. copy-level analysis with no reads, or a rare copy); only
/// a WELL-covered monomorphic column is dropped as an artifact.
pub const PSV_MIN_JUDGE_COV: usize = 4;

/// Read-support PSV validation. A candidate column (from copy-vs-copy alignment in `discover_psvs`) is kept
/// unless the READS positively contradict it — i.e. it has `>= min_judge_cov` coverage yet fewer than 2
/// distinct alleles reach `>= min_reads` reads. This drops assembly-artifact columns where the copies'
/// reference sequences differ but every molecule agrees (a mis-assembled base fabricates a difference no read
/// supports), while never dropping a low-coverage column we cannot judge — so a PSV is validated by a pileup,
/// not only by copy-vs-copy alignment. `all_obs[r][col]` = read r's observed base at `col` (`None` = uncovered).
pub fn read_supported_columns(
    all_obs: &[Vec<Option<u8>>],
    n_cols: usize,
    min_reads: usize,
    min_judge_cov: usize,
) -> Vec<bool> {
    (0..n_cols)
        .map(|col| {
            let mut counts: BTreeMap<u8, usize> = BTreeMap::new();
            let mut cov = 0usize;
            for obs in all_obs {
                if let Some(&Some(b)) = obs.get(col) {
                    cov += 1;
                    *counts.entry(b).or_insert(0) += 1;
                }
            }
            // keep if we cannot judge (too little coverage) OR >= 2 read-supported alleles are present
            cov < min_judge_cov || counts.values().filter(|&&c| c >= min_reads).count() >= 2
        })
        .collect()
}

/// Restrict `FamilyProfiles` to the columns marked `true` in `keep` (the read-support PSV filter): filters
/// each copy's allele vector, re-indexes `copy_gpos` into the surviving-column space, and updates `n_cols`.
/// `gen2off`/`strand` are column-independent and copied through.
fn restrict_family_profiles(fp: &FamilyProfiles, keep: &[bool]) -> FamilyProfiles {
    let idx: Vec<usize> = (0..fp.n_cols).filter(|&j| keep[j]).collect();
    let mut new_of = vec![usize::MAX; fp.n_cols];
    for (nj, &oj) in idx.iter().enumerate() {
        new_of[oj] = nj;
    }
    let profiles = fp
        .profiles
        .iter()
        .map(|pr| CopyProfile {
            copy_id: pr.copy_id,
            alleles: idx.iter().map(|&j| pr.alleles[j]).collect(),
            junctions: pr.junctions.clone(),
        })
        .collect();
    let copy_gpos = fp
        .copy_gpos
        .iter()
        .map(|gp| {
            let mut v: Vec<(usize, u64)> =
                gp.iter().filter(|&&(c, _)| keep[c]).map(|&(c, g)| (new_of[c], g)).collect();
            v.sort_unstable_by_key(|&(_, g)| g);
            v
        })
        .collect();
    FamilyProfiles {
        profiles,
        copy_gpos,
        gen2off: fp.gen2off.clone(),
        strand: fp.strand.clone(),
        n_cols: idx.len(),
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
fn fill_psv_obs(
    read: &AlignedRead,
    gpos: &[(usize, u64)],
    minus: bool,
    psv_obs: &mut [Option<u8>],
    psv_qual: &mut [Option<u8>],
) {
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
                    let idx = (seq_cur + (g - ref_cur)) as usize;
                    if let Some(&b) = read.seq.get(idx) {
                        psv_obs[col] = Some(if minus { rc_base(b) } else { b });
                        // QV is per physical base, unaffected by the minus-strand RC of the base.
                        psv_qual[col] = read.qual.get(idx).copied();
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
    let mut psv_qual = vec![None; fp.n_cols];
    fill_psv_obs(read, &fp.copy_gpos[mc], fp.strand[mc] == '-', &mut psv_obs, &mut psv_qual);
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
    ReadFeatures { psv_obs, psv_qual, junctions }
}

/// The copy whose genomic span the read overlaps most (`None` if it overlaps none).
pub(crate) fn best_overlap_copy(read: &AlignedRead, copies: &[&DenovoTranscript]) -> Option<usize> {
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
    /// this read's intron-boundary offsets in `mapped_copy`'s spliced space (Task H2: the same per-read
    /// junction evidence `combined`'s junction term already used, exposed so callers can thread it into the
    /// EM engine alongside `psv_obs`).
    pub junctions: Vec<i64>,
}

/// Two-stage freeze for reference-ABSENT (collapsed) copy admission.
///
/// `stage1` is the per-read assignment over the `n_ref` reference copies only; `stage2` is the re-assignment
/// over the ref copies (indices `0..n_ref`, unchanged) PLUS the admitted absent copies (indices `>= n_ref`).
/// The result is `stage2`-based EXCEPT that a read which was **Stage-1 Assigned at a multi-ref-copy family**
/// keeps its Stage-1 `combined`/`psv` (frozen — its decision was already certified against the ref copies and
/// must not be perturbed by adding a copy). Reads that were Stage-1 `Tied`/`Ambiguous`, reads absent from
/// Stage-1, and ALL reads when the family had a single ref copy (`n_ref <= 1`, where there was nothing to
/// freeze) keep their Stage-2 result. Any surviving Stage-2 assignment to an absent copy (`best_copy >=
/// n_ref`) is flagged `discovery_coupled`. Matched by `read_index` (NOT zip position — adding a copy changes
/// which reads get a result, so the two stages are not position-aligned). Deterministic (`BTreeMap` index).
pub(crate) fn freeze_merge(stage1: &[ReadResult], stage2: Vec<ReadResult>, n_ref: usize) -> Vec<ReadResult> {
    let single_ref = n_ref <= 1;
    let s1: BTreeMap<usize, &ReadResult> = stage1.iter().map(|r| (r.read_index, r)).collect();
    stage2
        .into_iter()
        .map(|mut r2| {
            let frozen = !single_ref
                && s1
                    .get(&r2.read_index)
                    .is_some_and(|r1| r1.combined.status == AssignStatus::Assigned);
            if frozen {
                let r1 = s1[&r2.read_index];
                // Stage-1 only saw the ref copies (`best_copy < n_ref`), so the frozen decision is valid in
                // the copies2 frame; never `discovery_coupled` (it predates the absent copy). Restore the WHOLE
                // assignment side — combined, psv AND mapped_copy — to Stage-1's, so the uniq-agreement/collapsed
                // diagnostics (`best_copy == mapped_copy`, `by_copy[mapped_copy]`) stay Stage-1-consistent and
                // a frozen read's mapped_copy can't flip to a Stage-2 absent index.
                r2.combined = r1.combined.clone();
                r2.psv = r1.psv.clone();
                r2.mapped_copy = r1.mapped_copy;
                // NOTE: do NOT touch `r2.psv_obs` — it must stay the Stage-2 (copies2-frame) observations, which
                // the abundance EM in `assign_family_detailed` consumes over the full copies2 column frame.
            } else if r2.combined.status == AssignStatus::Assigned && r2.combined.best_copy >= n_ref {
                r2.combined.discovery_coupled = true;
            }
            r2
        })
        .collect()
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
    /// per-`conversions` unified verdict (gene conversion vs RT/template-switch artifact vs chimera/ambiguous)
    /// from the recurrence + microhomology + DNA legs. Empty unless the caller classifies (genome-dependent).
    pub conversion_class: Vec<super::mosaic::Classification>,
    /// COPY-level historical gene conversions: a de-novo copy whose PSV-allele vector is itself a mosaic of two
    /// OTHER copies (the APOBEC3/RFPL signal — baked into the copy sequence, invisible to the read-level scan).
    pub copy_conversions: Vec<CopyConversion>,
    /// PSV column `j` → canonical forward-genome position (the column frame shared by reads + copies) — for the
    /// genotype-matrix visualization that proves the per-read assignment.
    pub psv_col_pos: Vec<Option<u64>>,
    /// `copy_psv_alleles[c][j]` = copy `c`'s allele at PSV column `j` (None = gapped) — the per-copy reference
    /// the reads are matched against in the genotype matrix.
    pub copy_psv_alleles: Vec<Vec<Option<u8>>>,
    /// `copy_junctions[c]` = copy `c`'s copy-specific intron-boundary offsets (`CopyProfile.junctions`,
    /// i.e. `copy_boundaries(c)`), parallel to `copy_psv_alleles` (Task H2: threads O3 junction evidence
    /// into the same per-copy frame the EM engine consumes).
    pub copy_junctions: Vec<Vec<i64>>,
    /// Indices into the input `copies` slice that survived iterative pruning. When pruning is off this is
    /// simply `0..copies.len()`. Callers use this to align the output detail with the original copy roster.
    pub copy_indices: Vec<usize>,
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

/// A copy must show the MINOR A/G allele at >= this fraction AND read count among its reads for a column to be
/// called an A-to-I editing site — far above the ~3e-4 sequencing-error rate, so genuine A/G paralog SNVs
/// (each copy monomorphic) and sequencing error never trip it.
const EDIT_MIN_FRAC: f64 = 0.05;
const EDIT_MIN_READS: u32 = 2;

/// Flag the PSV columns that are A-to-I RNA-editing sites (Clair3-RNA-style), so `assign_read_editing`
/// downweights them in the significance certificate. PSV alleles are in transcription orientation, so editing
/// is uniformly A→G. A column `j` is flagged iff (1) it is an **A↔G** column — the copy-consensus alleles at
/// `j` are exactly `{A, G}` — and (2) some copy shows **within-copy A/G heterogeneity**: provisionally assign
/// each read to its argmax-match copy (real PSVs dominate the editing-column minority), then a copy with both
/// alleles present at `j` (minor `>= EDIT_MIN_READS` and `>= EDIT_MIN_FRAC`) marks editing. A real A/G paralog
/// SNV has each copy monomorphic at `j` → not flagged. `reads_obs[r][j]` = read `r`'s base at column `j`.
pub(crate) fn detect_editing_columns(reads_obs: &[Vec<Option<u8>>], copies: &[CopyProfile]) -> Vec<bool> {
    let n = copies.len();
    let n_cols = copies.iter().map(|c| c.alleles.len()).max().unwrap_or(0);
    let mut flag = vec![false; n_cols];
    if n == 0 || n_cols == 0 {
        return flag;
    }
    // per (copy, column) -> (n_A, n_G) among reads provisionally assigned to that copy.
    let mut counts = vec![vec![(0u32, 0u32); n_cols]; n];
    for obs in reads_obs {
        // Columns this read spans — both the argmax match-count and the A/G tally only depend on spanned
        // columns, so iterate those instead of re-walking every column for every copy (identical result).
        let spanned: Vec<usize> = (0..n_cols).filter(|&j| obs.get(j).copied().flatten().is_some()).collect();
        let (mut best, mut best_m) = (0usize, -1i64);
        for (ci, c) in copies.iter().enumerate() {
            let mut m = 0i64;
            for &j in &spanned {
                if let (Some(ob), Some(al)) = (obs.get(j).copied().flatten(), c.alleles.get(j).copied().flatten()) {
                    if ob == al {
                        m += 1;
                    }
                }
            }
            if m > best_m {
                best_m = m;
                best = ci;
            }
        }
        for &j in &spanned {
            match obs.get(j).copied().flatten() {
                Some(b'A') => counts[best][j].0 += 1,
                Some(b'G') => counts[best][j].1 += 1,
                _ => {}
            }
        }
    }
    for j in 0..n_cols {
        // (1) A<->G column: copy-consensus alleles at j are exactly {A, G} (both present, nothing else).
        let (mut has_a, mut has_g, mut other) = (false, false, false);
        for c in copies {
            match c.alleles.get(j).copied().flatten() {
                Some(b'A') => has_a = true,
                Some(b'G') => has_g = true,
                Some(_) => other = true,
                None => {}
            }
        }
        if other || !(has_a && has_g) {
            continue;
        }
        // (2) within-copy A/G heterogeneity in some copy
        for ci in 0..n {
            let (na, ng) = counts[ci][j];
            let tot = na + ng;
            if tot == 0 {
                continue;
            }
            let minor = na.min(ng);
            if minor >= EDIT_MIN_READS && (minor as f64 / tot as f64) >= EDIT_MIN_FRAC {
                flag[j] = true;
                break;
            }
        }
    }
    flag
}

/// Pairwise profile distance between two copies: number of PSV columns where both are defined and differ
/// plus the number of junction boundaries present in exactly one copy. Used for iterative pruning's
/// nearest-neighbor computation.
fn profile_distance(a: &CopyProfile, b: &CopyProfile) -> usize {
    let psv_diff = a
        .alleles
        .iter()
        .zip(b.alleles.iter())
        .filter(|(x, y)| matches!((x, y), (Some(xa), Some(ya)) if xa != ya))
        .count();
    let mut junc_diff = 0usize;
    for &jb in &a.junctions {
        if !boundary_present(jb, &b.junctions, 4) {
            junc_diff += 1;
        }
    }
    for &jb in &b.junctions {
        if !boundary_present(jb, &a.junctions, 4) {
            junc_diff += 1;
        }
    }
    psv_diff + junc_diff
}

/// For each copy, return the index of its nearest profile neighbor (smallest `profile_distance`).
/// Ties break to the lowest index. A copy with no defined neighbor (single-copy family) returns 0.
fn nearest_neighbor_profiles(profiles: &[CopyProfile]) -> Vec<usize> {
    let n = profiles.len();
    let mut out = vec![0usize; n];
    if n < 2 {
        return out;
    }
    for i in 0..n {
        let mut best_j = if i == 0 { 1 } else { 0 };
        let mut best_d = profile_distance(&profiles[i], &profiles[best_j]);
        for j in 0..n {
            if j == i {
                continue;
            }
            let d = profile_distance(&profiles[i], &profiles[j]);
            if d < best_d || (d == best_d && j < best_j) {
                best_j = j;
                best_d = d;
            }
        }
        out[i] = best_j;
    }
    out
}

/// Find copies that have no read with significant evidence distinguishing them from their nearest neighbor.
/// Returns a list of `(copy_index, merge_target_index)` pairs. The caller removes one (or more) of these
/// copies and re-runs assignment.
fn find_weak_copies(
    detail: &FamilyDetail,
    copies: &[&DenovoTranscript],
    reads: &[AlignedRead],
    p: &AssignParams,
    genome: Option<&crate::genome::GenomeIndex>,
) -> Vec<(usize, usize)> {
    let n = copies.len();
    if n < 2 {
        return Vec::new();
    }
    let fp = build_family_profiles(copies, genome);
    let nn = nearest_neighbor_profiles(&fp.profiles);
    let thr = p.alpha / (n.saturating_sub(1).max(1) as f64);

    let mut sig_count = vec![0usize; n];
    for r in &detail.results {
        if r.combined.best_copy >= n {
            continue;
        }
        let i = r.combined.best_copy;
        let j = nn[i];
        // Recompute read features in the current profile frame so the pairwise certificate uses the
        // surviving copy set and its PSV columns / junction boundaries.
        let feats = read_features(&reads[r.read_index], r.mapped_copy, &fp);
        let (p_ij, _) = copy_pair_significance(&feats, &fp.profiles[i], &fp.profiles[j], p, &[]);
        if p_ij < thr {
            sig_count[i] += 1;
        }
    }

    let mut out: Vec<(usize, usize)> = (0..n)
        .filter(|&i| sig_count[i] == 0)
        .map(|i| (i, nn[i]))
        .collect();
    // Deterministic removal order: merge the copy closest to its neighbor first (smallest profile distance),
    // then re-evaluate. This avoids arbitrary index ordering effects.
    out.sort_by_key(|&(i, j)| profile_distance(&fp.profiles[i], &fp.profiles[j]));
    out
}

/// Like `assign_family` but returns the TWO-PASS detail per read so callers can report how many reads a
/// copy-specific junction resolved that PSVs alone could not (`junction_only`), and — with read mapq — the
/// unique-mapper agreement. Reads overlapping no copy are skipped.
///
/// `mol_names` — the BAM query name of each entry of `reads`, when the caller can supply it. `reads` is a
/// list of alignment RECORDS, and a multimapping molecule contributes one record per placement (the shipped
/// BAM is aligned with `-Y`, so its secondary records carry the full SEQ and each is independently
/// scorable). With names supplied, the unit of an assignment becomes the MOLECULE: see
/// `assign_family_detailed_once`. `None` (unit tests, callers with no name vector) keeps the historical
/// one-row-per-record behaviour exactly.
pub fn assign_family_detailed(
    copies: &[&DenovoTranscript],
    reads: &[AlignedRead],
    p: &AssignParams,
    genome: Option<&crate::genome::GenomeIndex>,
    mol_names: Option<&[String]>,
) -> FamilyDetail {
    assign_family_detailed_once(copies, reads, p, genome, mol_names)
}

/// IsoCon-style iterative copy pruning: repeatedly assign reads, identify copies with no significant
/// read-backed evidence against their nearest neighbor, merge the weakest such copy into its neighbor,
/// and reassign until all surviving copies are defensible.
///
/// Returns a `FamilyDetail` whose `copy_indices` field lists the original indices (into the input `copies`
/// slice) that survived. The caller should use those indices to align downstream bookkeeping with the
/// reduced copy set.
pub fn assign_family_detailed_pruned(
    copies: &[&DenovoTranscript],
    reads: &[AlignedRead],
    p: &AssignParams,
    genome: Option<&crate::genome::GenomeIndex>,
    mol_names: Option<&[String]>,
) -> FamilyDetail {
    // Phase 1: iteratively decide which original copies survive.  `merge_target[i]` always points from
    // original index i to the original index that currently represents it.  Initially each copy represents
    // itself.
    let mut merge_target: Vec<usize> = (0..copies.len()).collect();
    let mut current_indices: Vec<usize> = (0..copies.len()).collect();
    let mut current_copies: Vec<&DenovoTranscript> = copies.to_vec();
    loop {
        let detail = assign_family_detailed_once(&current_copies, reads, p, genome, mol_names);
        let weak = find_weak_copies(&detail, &current_copies, reads, p, genome);
        if weak.is_empty() {
            break;
        }
        // `weak` is reported in *current* index space; convert back to original indices for the merge map.
        let (remove_current, target_current) = weak[0];
        let removed_orig = current_indices[remove_current];
        let target_orig = current_indices[target_current];
        // Anything currently represented by `removed_orig` now rolls up to `target_orig`.
        for t in merge_target.iter_mut() {
            if *t == removed_orig {
                *t = target_orig;
            }
        }
        current_indices.remove(remove_current);
        current_copies.remove(remove_current);
    }

    if current_indices.len() == copies.len() {
        // Nothing was pruned — use the detail from the last full-size iteration.
        let mut detail = assign_family_detailed_once(copies, reads, p, genome, mol_names);
        detail.copy_indices = current_indices;
        return detail;
    }

    // Phase 2: assign against the FULL original copy set so reads that only overlapped a removed copy
    // are not dropped.  Then remap every per-copy index to its surviving output position.
    let full_detail = assign_family_detailed_once(copies, reads, p, genome, mol_names);

    // Map original copy index -> output copy index (position in `current_indices`).
    let mut orig_to_out: Vec<Option<usize>> = vec![None; copies.len()];
    for (out_idx, &orig_idx) in current_indices.iter().enumerate() {
        orig_to_out[orig_idx] = Some(out_idx);
    }
    // Resolve an original copy index through the merge chain to a surviving output index.
    let resolve = |mut i: usize| -> Option<usize> {
        let mut seen = BTreeSet::new();
        while seen.insert(i) {
            if let Some(out) = orig_to_out[i] {
                return Some(out);
            }
            i = merge_target[i];
            if i >= orig_to_out.len() {
                return None;
            }
        }
        None
    };

    // Helper: remap an Assignment from original copy space to surviving output copy space.
    // If all copies with meaningful posterior mass collapse to the same output, the read is considered
    // assigned to that output even if the original call was Tied among them.
    let remap_assignment = |a: &Assignment| -> Option<Assignment> {
        let best_out = resolve(a.best_copy)?;
        // Sum posterior probabilities into output buckets.
        let mut post_out: Vec<f64> = vec![0.0; current_indices.len()];
        for (orig, &prob) in a.posterior.iter().enumerate() {
            if let Some(out) = resolve(orig) {
                post_out[out] += prob;
            }
        }
        // After merging, a read is resolvable to an output if the bulk of its posterior mass lands in that
        // single bucket.  Keep the original status when the mass is split across outputs.
        let conf = post_out[best_out];
        let (status, best_out) = if conf >= 0.9 {
            (AssignStatus::Assigned, best_out)
        } else {
            (a.status, best_out)
        };
        // Renormalize posterior so it still sums to 1.
        let sum: f64 = post_out.iter().sum();
        if sum > 0.0 {
            post_out.iter_mut().for_each(|v| *v /= sum);
        }
        Some(Assignment {
            best_copy: best_out,
            log_lr_margin: a.log_lr_margin,
            n_decisive: a.n_decisive,
            resolvable: a.resolvable,
            status,
            p_value: a.p_value,
            min_p_value: a.min_p_value,
            discovery_coupled: a.discovery_coupled,
            junction_conflict: a.junction_conflict,
            origin_rejected: a.origin_rejected,
            posterior: post_out,
        })
    };

    let mut results: Vec<ReadResult> = Vec::with_capacity(full_detail.results.len());
    let mut read_obs_for_em: Vec<Vec<Option<u8>>> = Vec::new();
    for r in full_detail.results {
        let mapped_out = resolve(r.mapped_copy);
        let psv_out = remap_assignment(&r.psv);
        let combined_out = remap_assignment(&r.combined);
        if mapped_out.is_none() || psv_out.is_none() || combined_out.is_none() {
            continue;
        }
        // Reads whose PSV observations contributed to the EM in the full assignment still contribute after
        // merging (the observations are in the shared canonical column frame).
        read_obs_for_em.push(r.psv_obs.clone());
        results.push(ReadResult {
            read_index: r.read_index,
            mapped_copy: mapped_out.unwrap(),
            psv: psv_out.unwrap(),
            combined: combined_out.unwrap(),
            psv_obs: r.psv_obs,
            junctions: r.junctions,
        });
    }

    // Recompute soft abundance on the surviving copy set.
    let surviving_alleles: Vec<Vec<Option<u8>>> = current_indices
        .iter()
        .map(|&orig| full_detail.copy_psv_alleles[orig].clone())
        .collect();
    let surviving_junctions: Vec<Vec<i64>> = current_indices
        .iter()
        .map(|&orig| full_detail.copy_junctions[orig].clone())
        .collect();
    let copy_abundance = soft_quantify_em(&read_obs_for_em, &surviving_alleles, QUANT_ERROR, 100);
    let n_eff = results.iter().filter(|r| r.combined.n_decisive >= 1).count();
    let copy_abundance_ci: Vec<f64> = if n_eff == 0 {
        vec![0.5; copy_abundance.len()]
    } else {
        let n = n_eff as f64;
        copy_abundance
            .iter()
            .map(|&t| (1.96 * (t * (1.0 - t) / n).sqrt()).min(0.5))
            .collect()
    };

    // Remap copy-level conversions to surviving output indices; drop any whose focal copy was removed.
    let mut copy_conversions: Vec<CopyConversion> = Vec::new();
    for cc in full_detail.copy_conversions {
        if let (Some(c_out), Some(a_out), Some(b_out)) =
            (resolve(cc.copy_c), resolve(cc.copy_a), resolve(cc.copy_b))
        {
            // Avoid self-conversions created by merging.
            if c_out != a_out && c_out != b_out && a_out != b_out {
                copy_conversions.push(CopyConversion {
                    copy_c: c_out,
                    copy_a: a_out,
                    copy_b: b_out,
                    breakpoint: cc.breakpoint,
                    n_decisive: cc.n_decisive,
                });
            }
        }
    }

    FamilyDetail {
        results,
        n_cols: full_detail.n_cols,
        copy_abundance,
        copy_abundance_ci,
        mosaic_reads: full_detail.mosaic_reads,
        conversions: full_detail.conversions,
        conversion_class: Vec::new(),
        copy_conversions,
        psv_col_pos: full_detail.psv_col_pos,
        copy_psv_alleles: surviving_alleles,
        copy_junctions: surviving_junctions,
        copy_indices: current_indices,
    }
}

/// # `mol_names`: RECORD in, MOLECULE out
///
/// `reads` is a list of alignment RECORDS. A multimapping molecule contributes one record per placement,
/// and the shipped BAMs are aligned with minimap2 `-Y`, which writes the FULL SEQ on secondary records —
/// so every one of a molecule's records is an independently scorable observation and each gets its own
/// verdict below. That is what the assignment stage WANTS as input (a molecule's records corroborate each
/// other; on the planted-truth `-Y` ladder 998/1000 molecules had ≥2 `Assigned` records and they AGREED),
/// but it is NOT a valid output unit: a molecule comes from exactly one copy, so emitting one row per
/// record makes one molecule several "reads", double-counts it in every per-copy read count, and — where
/// the records disagree — asserts that one molecule came from two copies at once.
///
/// With `mol_names` supplied, each molecule is reduced to ONE result, with no new constant:
///
/// * **contradiction ⇒ abstain.** If ≥2 of the molecule's records are `Assigned` and they name ≥2 DISTINCT
///   copies, the molecule is `Ambiguous`. This is O2's assign-or-abstain contract applied to the molecule,
///   and it is the only admissible arbitration: on the gorilla SHARP family every one of the 323 native
///   contradicting molecules (323/323) has its LARGEST margin on the record carrying the PRIMARY flag, so
///   any "keep the strongest record" rule is the aligner's primary flag in disguise — exactly the defect
///   that retired the `uniq_agree` metric.
/// * otherwise the representative is the record that OBSERVED the most (max `n_decisive`, then max margin,
///   then min `p_value`, then lowest record index). When the records agree, the decision is theirs
///   unanimously and this choice only picks which certificate is reported.
///
/// The reduction is applied before `results`, the EM observation vector and the mosaic calls are built, so
/// every downstream statistic (`.assignments.tsv` rows, `families.tsv` counters, `quant.tsv`
/// `n_reads_hard` and `abundance`, `collapsed_copies`, the iterative-prune decisions) counts molecules.
///
/// `None` keeps the historical one-row-per-record behaviour byte-for-byte.
fn assign_family_detailed_once(
    copies: &[&DenovoTranscript],
    reads: &[AlignedRead],
    p: &AssignParams,
    genome: Option<&crate::genome::GenomeIndex>,
    mol_names: Option<&[String]>,
) -> FamilyDetail {
    if copies.len() < 2 {
        return FamilyDetail {
            results: Vec::new(),
            n_cols: 0,
            copy_abundance: Vec::new(),
            copy_abundance_ci: Vec::new(),
            mosaic_reads: 0,
            conversions: Vec::new(),
            conversion_class: Vec::new(),
            copy_conversions: Vec::new(),
            psv_col_pos: Vec::new(),
            copy_psv_alleles: Vec::new(),
            copy_junctions: Vec::new(),
            copy_indices: Vec::new(),
        };
    }
    use super::mosaic::{aggregate_family, detect_mosaic, MosaicParams, SiteObs};
    const MOSAIC_EPS: f64 = 0.01; // HiFi per-base error for the mosaic likelihood
    const MAX_MOSAIC_SITES: usize = 250; // cap PSV sites per detect_mosaic (it is O(sites^2)); stride-sample
    let timing = std::env::var_os("RUSTLE_TIMING").is_some();
    let t_psv = std::time::Instant::now();
    let fp0 = build_family_profiles(copies, genome);
    if timing {
        eprintln!(
            "[timing]     build_family_profiles/discover_psvs ({} copies, {} cols): {:.1}s",
            copies.len(),
            fp0.n_cols,
            t_psv.elapsed().as_secs_f64()
        );
    }
    // READ-SUPPORT PSV validation. The candidate columns come from copy-vs-copy alignment (`discover_psvs`);
    // keep a column only if the READS observe >= 2 alleles there (>= PSV_MIN_ALLELE_READS each), dropping
    // assembly-artifact columns where the paralogs' reference sequences differ but every molecule agrees. So a
    // PSV is validated by a per-read pileup, not by alignment alone. This pass also builds `all_obs`, reused by
    // the editing pre-pass. Since 2026-09-05 the `copy_assign` binary sets RUSTLE_PSV_READFILTER=0 unless the user
    // passes `--psv-read-filter` or sets the variable (§6eu, register 689: the filter deletes every column of an
    // unexpressed paralogue). Library callers that leave the variable unset still get the filter ON.
    let mut all_obs: Vec<Vec<Option<u8>>> = Vec::with_capacity(reads.len());
    for read in reads {
        if let Some(mc) = best_overlap_copy(read, copies) {
            let mut psv_obs = vec![None; fp0.n_cols];
            let mut psv_qual = vec![None; fp0.n_cols];
            fill_psv_obs(read, &fp0.copy_gpos[mc], fp0.strand[mc] == '-', &mut psv_obs, &mut psv_qual);
            all_obs.push(psv_obs);
        }
    }
    let read_filter = std::env::var("RUSTLE_PSV_READFILTER").ok().as_deref() != Some("0");
    let keep: Vec<bool> = if read_filter {
        read_supported_columns(&all_obs, fp0.n_cols, PSV_MIN_ALLELE_READS, PSV_MIN_JUDGE_COV)
    } else {
        vec![true; fp0.n_cols]
    };
    let n_dropped = keep.iter().filter(|&&k| !k).count();
    let fp = if n_dropped > 0 {
        // reindex `all_obs` into the surviving-column space so the editing pre-pass stays aligned to `fp`
        let idx: Vec<usize> = (0..keep.len()).filter(|&j| keep[j]).collect();
        for obs in all_obs.iter_mut() {
            let filtered: Vec<Option<u8>> = idx.iter().map(|&j| obs[j]).collect();
            *obs = filtered;
        }
        restrict_family_profiles(&fp0, &keep)
    } else {
        fp0
    };
    if timing && n_dropped > 0 {
        eprintln!("[timing]     read-support PSV filter: dropped {}/{} candidate columns", n_dropped, keep.len());
    }
    // canonical column -> genomic position (first copy that has the column) so every read's switch breakpoints
    // are in ONE frame and `aggregate_family` can cluster recurrences across molecules.
    let mut col_canon: Vec<Option<u64>> = vec![None; fp.n_cols];
    for gpos in &fp.copy_gpos {
        for &(col, g) in gpos {
            col_canon[col].get_or_insert(g);
        }
    }
    let mparams = MosaicParams::from_env();
    // RNA-editing filter (Clair3-RNA): a pre-pass flags A↔G columns with within-copy heterogeneity so the
    // significance certificate downweights them. Reuses the per-read PSV observations built above.
    let editing_cols: Vec<bool> = if p.rna_editing_filter {
        detect_editing_columns(&all_obs, &fp.profiles)
    } else {
        Vec::new()
    };
    let t_assign = std::time::Instant::now();
    // Per-read assignment is independent across reads, so compute it in parallel and merge in read order
    // (rayon's indexed collect preserves order, keeping output identical to the serial loop). The merge
    // reproduces the serial sequencing exactly: a read with an overlap copy always contributes a mosaic
    // call; it adds an EM observation only if `combined` assigns; and a `ReadResult` only if `psv` also
    // assigns.
    struct PerRead {
        mcall: super::mosaic::MosaicCall,
        obs_for_em: Option<Vec<Option<u8>>>,
        result: Option<ReadResult>,
    }
    // Built once per family (from the family's copy roster) and shared by every read in the loop below.
    let graph = BubbleGraph::from_copies(&fp.profiles);
    // The per-record tail (mosaic call + PSV/junction assignment) shared by both observation modes.
    let finish = |ri: usize, mc: usize, feats: ReadFeatures| -> PerRead {
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
            let Some(combined) = assign_read_editing(&feats, &graph, &fp.profiles, p, &editing_cols) else {
                return PerRead { mcall, obs_for_em: None, result: None };
            };
            let obs = feats.psv_obs.clone();
            let junctions = feats.junctions.clone();
            let psv_feats = ReadFeatures { psv_obs: feats.psv_obs, psv_qual: feats.psv_qual, junctions: vec![] };
            let Some(psv) = assign_read_editing(&psv_feats, &graph, &fp.profiles, p, &editing_cols) else {
                return PerRead { mcall, obs_for_em: Some(obs), result: None };
            };
            PerRead {
                mcall,
                obs_for_em: Some(obs.clone()),
                result: Some(ReadResult {
                    read_index: ri,
                    mapped_copy: mc,
                    psv,
                    combined,
                    psv_obs: obs,
                    junctions,
                }),
            }
    };
    // ⭐ O2-9 / D3 (`p.molecule_pool`, "the read is the star"): one observation per MOLECULE, built from the
    // molecule's own sequence aligned to every copy's unit; the family's star-projected columns are not used.
    let pooled = p.molecule_pool && mol_names.is_some();
    let per_read: Vec<Option<PerRead>> = if pooled {
        let names = mol_names.expect("pooled requires mol_names");
        debug_assert_eq!(names.len(), reads.len());
        let mut groups: std::collections::HashMap<&str, Vec<usize>> = std::collections::HashMap::new();
        let mut order: Vec<&str> = Vec::new();
        for (i, n) in names.iter().enumerate().take(reads.len()) {
            let e = groups.entry(n.as_str()).or_insert_with(|| {
                order.push(n.as_str());
                Vec::new()
            });
            e.push(i);
        }
        // representative record per molecule: the one overlapping a copy, with the longest sequence
        let reps: Vec<usize> = order
            .iter()
            .filter_map(|name| {
                groups[name]
                    .iter()
                    .copied()
                    .filter(|&i| best_overlap_copy(&reads[i], copies).is_some())
                    .max_by(|&a, &b| reads[a].seq.len().cmp(&reads[b].seq.len()).then(b.cmp(&a)))
            })
            .collect();
        let copy_seqs: Vec<&[u8]> = copies.iter().map(|c| c.seq.as_slice()).collect();
        let read_seqs: Vec<&[u8]> = reps.iter().map(|&i| reads[i].seq.as_slice()).collect();
        // each copy's splice boundaries in unit offsets: where the spliced->genomic map jumps
        let copy_bounds: Vec<Vec<u32>> = copies
            .iter()
            .map(|c| {
                let m = exon_map(c);
                (1..m.len()).filter(|&o| m[o] != m[o - 1] + 1 && m[o] + 1 != m[o - 1]).map(|o| o as u32).collect()
            })
            .collect();
        let t_star = std::time::Instant::now();
        let alns = read_star_observations(&copy_seqs, &read_seqs, &copy_bounds);
        if timing {
            eprintln!("[timing]     read-star minimap2 ({} molecules x {} copies): {:.1}s", reps.len(), copies.len(), t_star.elapsed().as_secs_f64());
        }
        let done: Vec<(usize, PerRead)> = reps
            .par_iter()
            .zip(alns.par_iter())
            .map(|(&ri, (cand, obs, profiles, unit_edits, covered))| {
                let read = &reads[ri];
                let mc = best_overlap_copy(read, copies).expect("rep overlaps a copy");
                let (cand, obs, profiles) = (cand.clone(), obs.clone(), profiles.clone());
                // the molecule's best GENOME placement (fewest edits per aligned base over all its records;
                // needs `=`/`X` CIGARs, else None and the certificate below is skipped)
                let genome_rate: Option<f64> = groups[names[ri].as_str()]
                    .iter()
                    .filter_map(|&i| {
                        let (mut eq, mut ed) = (0u64, 0u64);
                        let mut has_eqx = false;
                        for &(op, n) in &reads[i].cigar {
                            match op {
                                '=' => { eq += n; has_eqx = true; }
                                'X' => { ed += n; has_eqx = true; }
                                'I' | 'D' => ed += n,
                                _ => {}
                            }
                        }
                        (has_eqx && eq + ed > 0).then(|| ed as f64 / (eq + ed) as f64)
                    })
                    .fold(None, |m: Option<f64>, r| Some(m.map_or(r, |x| x.min(r))));
                // a single candidate has no competitor to certify against: `assign_read` returns Tied for it
                // (margin ∞, nothing rejected), so the molecule stays in the output as a tie, never as a claim.
                if cand.is_empty() {
                    let mcall = detect_mosaic(&[], copies.len(), MOSAIC_EPS, &mparams);
                    return (ri, PerRead { mcall, obs_for_em: None, result: None });
                }
                if std::env::var_os("RUSTLE_STAR_DEBUG").is_some() {
                    eprintln!("[star] read {ri} len={} candidates={:?} cols={}", read.seq.len(), cand, obs.len());
                }
                // the read's own splice junctions as READ positions (the first base after each intron), from the
                // record whose sequence was aligned — the same coordinate system as the candidates' junctions
                let read_junctions: Vec<i64> = {
                    let mut out = Vec::new();
                    let mut qpos = 0i64;
                    for &(op, n) in &read.cigar {
                        match op {
                            'M' | '=' | 'X' | 'I' | 'S' => qpos += n as i64,
                            'N' => out.push(qpos),
                            _ => {}
                        }
                    }
                    out
                };
                let feats = ReadFeatures { psv_obs: obs, psv_qual: Vec::new(), junctions: if p.read_star_junctions { read_junctions } else { Vec::new() } };
                let mcall = detect_mosaic(&[], copies.len(), MOSAIC_EPS, &mparams);
                // ⭐ PAIRWISE (§6fc): best = the candidate with the most matching bases in its alignment; each other
                // candidate is certified against it on the columns BOTH carry (`copy_pair_significance`). Tied when
                // some competitor shares no distinguishing column (K = 0 for that pair); Ambiguous when some
                // competitor is not rejected at alpha/(n-1); Assigned when every competitor is.
                let matches = |k: usize| unit_edits.get(cand[k]).copied().flatten().map_or(0i64, |(_, _, nm, blk)| blk as i64 - nm as i64);
                let bk = (0..cand.len()).max_by_key(|&k| (matches(k), std::cmp::Reverse(k))).unwrap();
                let thr = p.alpha / (cand.len().saturating_sub(1).max(1) as f64);
                let (mut min_p, mut p_read, mut n_dec_min, mut margin_min, mut k0) = (1.0f64, 0.0f64, usize::MAX, f64::INFINITY, false);
                let mut posterior_ll = vec![0.0f64; cand.len()];
                let lr = ((1.0 - p.error_rate) / (p.error_rate / 3.0)).ln();
                for k in 0..cand.len() {
                    if k == bk {
                        continue;
                    }
                    let both: Vec<i64> = feats
                        .junctions
                        .iter()
                        .copied()
                        .filter(|&rj| {
                            let t = p.boundary_tol.max(0);
                            let cov = |iv: &[(u32, u32)]| iv.iter().any(|&(s, e)| (s as i64) <= rj + t && rj - t < e as i64);
                            cov(&covered[bk]) && cov(&covered[k])
                        })
                        .collect();
                    let feats_pair = ReadFeatures { psv_obs: feats.psv_obs.clone(), psv_qual: Vec::new(), junctions: both };
                    let (pk, _) = copy_pair_significance(&feats_pair, &profiles[bk], &profiles[k], p, &[]);
                    let (mut dec, mut llr) = (0usize, 0.0f64);
                    for j in 0..feats.psv_obs.len() {
                        if let (Some(o), Some(ba), Some(ca)) = (feats.psv_obs[j], profiles[bk].alleles[j], profiles[k].alleles[j]) {
                            if ba != ca {
                                dec += 1;
                                if o == ba { llr += lr } else if o == ca { llr -= lr }
                            }
                        }
                    }
                    // junction terms: a read junction present in one candidate and not the other is decisive
                    let lrj = ((1.0 - p.junction_err) / p.junction_err.max(1e-12)).ln();
                    let covers = |iv: &[(u32, u32)], r: i64| -> bool {
                        let t = p.boundary_tol.max(0);
                        iv.iter().any(|&(s, e)| (s as i64) <= r + t && r - t < e as i64)
                    };
                    for &rj in &feats.junctions {
                        // only where BOTH candidates have aligned bases around the junction position
                        if !(covers(&covered[bk], rj) && covers(&covered[k], rj)) {
                            continue;
                        }
                        let in_b = boundary_present(rj, &profiles[bk].junctions, p.boundary_tol);
                        let in_c = boundary_present(rj, &profiles[k].junctions, p.boundary_tol);
                        if in_b != in_c {
                            dec += 1;
                            llr += if in_b { lrj } else { -lrj };
                        }
                    }
                    if dec == 0 {
                        k0 = true;
                    }
                    min_p = min_p.min(pk);
                    p_read = p_read.max(pk);
                    n_dec_min = n_dec_min.min(dec);
                    margin_min = margin_min.min(llr);
                    posterior_ll[k] = -llr;
                }
                let single = cand.len() < 2;
                let status = if single || k0 {
                    AssignStatus::Tied
                } else if p_read < thr && margin_min > 0.0 {
                    AssignStatus::Assigned
                } else {
                    AssignStatus::Ambiguous
                };
                let m = posterior_ll.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let mut post: Vec<f64> = posterior_ll.iter().map(|&l| (l - m).exp()).collect();
                let z: f64 = post.iter().sum();
                for x in post.iter_mut() {
                    *x /= z;
                }
                let mut full = vec![0.0f64; copies.len()];
                for (k, &ci) in cand.iter().enumerate() {
                    full[ci] = post[k];
                }
                let mut a = Assignment {
                    best_copy: cand[bk],
                    log_lr_margin: if single { 0.0 } else { margin_min },
                    n_decisive: if single { 0 } else { n_dec_min },
                    resolvable: !single && !k0 && min_p < thr,
                    status,
                    p_value: if single { 1.0 } else { p_read },
                    min_p_value: if single { 1.0 } else { min_p },
                    discovery_coupled: false,
                    junction_conflict: false,
                    origin_rejected: false,
                    posterior: full,
                };
                // ⭐ origin certificate: under H0 "the best candidate's unit IS this read's origin", the read's
                // SUBSTITUTIONS against that unit are Binomial(aligned bases, error_rate) — indels are excluded on
                // purpose: a retained intron or an exon the majority chain lacks is isoform STRUCTURE, not evidence
                // about the copy of origin (with NM the certificate rejected minority isoforms of the right copy) — the same sequence up to sequencing error
                // (allelic variation, ~0.1–0.2 %, sits inside that tail). If the unit's edits are significantly
                // MORE (one-sided, alpha), no candidate explains the read: its origin's unit lacks the read's
                // content, or its origin is a copy the catalog does not hold (NPIP: anchor reads hit their own
                // copy's unit not at all and a paralogue's at 0.54; reads placed at MAPQ 60 with NM 333) ⟹ Ambiguous.
                // The read's genome placement is NOT the null: a placement at NM 333/4,615 is no origin either.
                let _ = genome_rate;
                if let Some(&Some((n_x, n_aligned, nm_all, blk_all))) = unit_edits.get(a.best_copy) {
                    let (nm, blk) = if p.origin_subst_only { (n_x, n_aligned) } else { (nm_all, blk_all) };
                    let r0 = p.error_rate;
                    let n = blk as f64;
                    let mean = n * r0;
                    let sd = (n * r0 * (1.0 - r0)).sqrt().max(1e-9);
                    let z = (nm as f64 - mean) / sd; // normal approximation of the binomial tail
                    let p_tail = 0.5 * (1.0 - erf_approx(z / std::f64::consts::SQRT_2));
                    if nm as f64 > mean && p_tail < p.alpha {
                        a.status = AssignStatus::Ambiguous;
                        a.resolvable = false;
                        a.origin_rejected = true;
                    }
                }
                if cand.len() < 2 {
                    // no competitor above the reporting floor: no certificate is possible, so no claim
                    a.status = AssignStatus::Tied;
                    a.resolvable = false;
                    a.log_lr_margin = 0.0;
                    a.n_decisive = 0;
                    a.p_value = 1.0;
                    a.min_p_value = 1.0;
                    a.posterior = vec![1.0 / copies.len() as f64; copies.len()];
                }
                if std::env::var_os("RUSTLE_STAR_DEBUG").is_some() && !feats.psv_obs.is_empty() {
                    let match_counts: Vec<usize> = profiles.iter().map(|pr| pr.alleles.iter().zip(feats.psv_obs.iter()).filter(|(x, o)| x.is_some() && *x == *o).count()).collect();
                    eprintln!("[star]   verdict best={} status={:?} n_dec={} margin={:.2} min_p={:.2e} p={:.2e} matches={:?}", a.best_copy, a.status, a.n_decisive, a.log_lr_margin, a.min_p_value, a.p_value, match_counts);
                }
                (
                    ri,
                    PerRead {
                        mcall,
                        obs_for_em: None, // the family EM works in the star-projected columns; not fed here
                        result: Some(ReadResult {
                            read_index: ri,
                            mapped_copy: mc,
                            psv: a.clone(),
                            combined: a,
                            psv_obs: vec![None; fp.n_cols],
                            junctions: Vec::new(),
                        }),
                    },
                )
            })
            .collect();
        let mut out: Vec<Option<PerRead>> = (0..reads.len()).map(|_| None).collect();
        for (ri, pr) in done {
            out[ri] = Some(pr);
        }
        out
    } else {
        reads
            .par_iter()
            .enumerate()
            .map(|(ri, read)| {
                let mc = best_overlap_copy(read, copies)?;
                Some(finish(ri, mc, read_features(read, mc, &fp)))
            })
            .collect()
    };
    if timing {
        eprintln!(
            "[timing]     parallel per-read assign ({} reads): {:.1}s",
            reads.len(),
            t_assign.elapsed().as_secs_f64()
        );
    }
    // RECORD -> MOLECULE reduction (see this function's doc comment). No-op when `mol_names` is `None`.
    let per_read: Vec<Option<PerRead>> = match mol_names {
        _ if pooled => per_read, // one result per molecule already; no representative, no contradiction rule
        None => per_read,
        Some(names) => {
            debug_assert_eq!(names.len(), reads.len());
            let mut groups: std::collections::HashMap<&str, Vec<usize>> =
                std::collections::HashMap::new();
            let mut order: Vec<&str> = Vec::new();
            for (i, n) in names.iter().enumerate().take(per_read.len()) {
                let e = groups.entry(n.as_str()).or_insert_with(|| {
                    order.push(n.as_str());
                    Vec::new()
                });
                e.push(i);
            }
            let mut keep: Vec<bool> = vec![false; per_read.len()];
            let mut downgrade: Vec<bool> = vec![false; per_read.len()];
            for name in &order {
                let idx = &groups[name];
                // present records only; a record the assigner produced nothing for cannot represent the
                // molecule unless no record did.
                let with_result: Vec<usize> =
                    idx.iter().copied().filter(|&i| per_read[i].as_ref().is_some_and(|pr| pr.result.is_some())).collect();
                let pool: &[usize] = if with_result.is_empty() { idx } else { &with_result };
                // "most-observed" representative: max n_decisive, then max margin, then min p, then index.
                let rep = pool
                    .iter()
                    .copied()
                    .max_by(|&a, &b| {
                        let key = |i: usize| {
                            per_read[i].as_ref().and_then(|pr| pr.result.as_ref()).map(|r| {
                                (r.combined.n_decisive, r.combined.log_lr_margin, -r.combined.p_value)
                            })
                        };
                        match (key(a), key(b)) {
                            (Some(ka), Some(kb)) => ka
                                .0
                                .cmp(&kb.0)
                                .then(ka.1.partial_cmp(&kb.1).unwrap_or(std::cmp::Ordering::Equal))
                                .then(ka.2.partial_cmp(&kb.2).unwrap_or(std::cmp::Ordering::Equal))
                                .then(b.cmp(&a)), // lower record index wins ties
                            (Some(_), None) => std::cmp::Ordering::Greater,
                            (None, Some(_)) => std::cmp::Ordering::Less,
                            (None, None) => b.cmp(&a),
                        }
                    })
                    .expect("group is non-empty");
                keep[rep] = true;
                // CONTRADICTION: ≥2 records `Assigned` to ≥2 distinct copies ⇒ the molecule abstains.
                let mut named: Vec<usize> = idx
                    .iter()
                    .filter_map(|&i| per_read[i].as_ref().and_then(|pr| pr.result.as_ref()))
                    .filter(|r| r.combined.status == AssignStatus::Assigned)
                    .map(|r| r.combined.best_copy)
                    .collect();
                named.sort_unstable();
                named.dedup();
                if named.len() >= 2 {
                    downgrade[rep] = true;
                }
            }
            per_read
                .into_iter()
                .enumerate()
                .map(|(i, pr)| {
                    if !keep[i] {
                        return None;
                    }
                    pr.map(|mut pr| {
                        if downgrade[i] {
                            if let Some(r) = pr.result.as_mut() {
                                r.combined.status = AssignStatus::Ambiguous;
                            }
                        }
                        pr
                    })
                })
                .collect()
        }
    };
    let mut results = Vec::new();
    let mut read_obs: Vec<Vec<Option<u8>>> = Vec::new(); // per-read PSV observations for the soft EM
    let mut mosaic_calls = Vec::new();
    let mut mosaic_reads = 0usize;
    for pr in per_read.into_iter().flatten() {
        if pr.mcall.is_mosaic() {
            mosaic_reads += 1;
        }
        mosaic_calls.push(pr.mcall);
        if let Some(o) = pr.obs_for_em {
            read_obs.push(o);
        }
        if let Some(r) = pr.result {
            results.push(r);
        }
    }
    // Breakpoints are in the shared `col_canon` frame (the co-located family's chromosome), so stamp
    // every call with the family chrom so each emitted event carries it (for the microhomology check).
    let fam_chrom = copies.first().map(|c| c.chrom.as_str()).unwrap_or("");
    let mosaic_chroms = vec![fam_chrom; mosaic_calls.len()];
    let conversions = aggregate_family(&mosaic_calls, &mosaic_chroms, &mparams);
    // copy-level historical conversions: is any copy's PSV-allele vector a mosaic of two others?
    let copy_conversions = scan_copy_conversions(&fp.profiles, &col_canon, &mparams);
    // soft per-copy abundance (EM) + a normal-approx 95% CI half-width.
    let copy_alleles: Vec<Vec<Option<u8>>> = fp.profiles.iter().map(|pr| pr.alleles.clone()).collect();
    let copy_abundance = soft_quantify_em(&read_obs, &copy_alleles, QUANT_ERROR, 100);
    // L8: the CI must track INFORMATIVE-PSV coverage (reads carrying >= 1 decisive feature), NOT the raw
    // read count. With raw N the half-width shrinks as 1/sqrt(N) even in the K=0 / non-identifiable regime
    // (all reads Tied, n_decisive=0) where the per-copy fractions are unidentifiable — false precision on
    // a default user-facing output. So: n_eff = #reads with a decisive (PSV-or-junction) feature; when
    // n_eff = 0 the abundance is unidentifiable and the CI is the full-simplex half-width (0.5); otherwise
    // clamp to 0.5 so it never claims more certainty than the [0,1] interval allows.
    let n_eff = results.iter().filter(|r| r.combined.n_decisive >= 1).count();
    let copy_abundance_ci: Vec<f64> = if n_eff == 0 {
        vec![0.5; copy_abundance.len()]
    } else {
        let n = n_eff as f64;
        copy_abundance
            .iter()
            .map(|&t| (1.96 * (t * (1.0 - t) / n).sqrt()).min(0.5))
            .collect()
    };
    FamilyDetail {
        results,
        n_cols: fp.n_cols,
        copy_abundance,
        copy_abundance_ci,
        mosaic_reads,
        conversions,
        conversion_class: Vec::new(), // populated by the caller via classify_conversions (genome-dependent)
        copy_conversions,
        psv_col_pos: col_canon,
        copy_psv_alleles: copy_alleles,
        copy_junctions: fp.profiles.iter().map(|pr| pr.junctions.clone()).collect(),
        copy_indices: (0..copies.len()).collect(),
    }
}

/// Classify each confirmed conversion event with the unified discriminator (recurrence already in the
/// event's `confirmed`; microhomology from the genome at the breakpoint bracket; DNA support optional).
/// Genome-dependent, so it runs at the pipeline call site (not inside `assign_family_detailed`). Each
/// event carries its own `chrom` (set by `aggregate_family`). The `dna_support` closure maps an event →
/// `Some(true/false)` if a DNA catalog was consulted, else `None`.
pub fn classify_conversions(
    detail: &super::copy_assign_pipeline::FamilyDetail,
    genome: &crate::genome::GenomeIndex,
    dna_support: impl Fn(&super::mosaic::ConversionEvent) -> Option<bool>,
) -> Vec<super::mosaic::Classification> {
    detail
        .conversions
        .iter()
        .map(|ev| super::mosaic::classify_event(ev, event_microhomology(genome, ev), dna_support(ev)))
        .collect()
}

/// Microhomology leg for one event: the template-switch direct-repeat signature at the breakpoint
/// bracket, on the event's own chromosome. `None` when the bracket is unusable (coord 0).
pub fn event_microhomology(
    genome: &crate::genome::GenomeIndex,
    ev: &super::mosaic::ConversionEvent,
) -> Option<bool> {
    const MH_KMIN: u64 = 6;
    const MH_KMAX: u64 = 12;
    let (lo, hi) = ev.breakpoint_ref;
    if lo > 0 && hi > 0 {
        Some(genome.breakpoint_microhomology(&ev.chrom, lo, hi, MH_KMIN, MH_KMAX))
    } else {
        None
    }
}

/// Assign every read over a co-located family to a copy. Each read is mapped to the copy whose genomic span
/// it overlaps most (reads overlapping no copy are skipped). Returns `(read_index, Assignment)`. Mirrors
/// `assign_family`.
pub fn assign_family(
    copies: &[&DenovoTranscript],
    reads: &[AlignedRead],
    p: &AssignParams,
    genome: Option<&crate::genome::GenomeIndex>,
) -> Vec<(usize, Assignment)> {
    if copies.len() < 2 {
        return Vec::new();
    }
    let fp = build_family_profiles(copies, genome);
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

    #[test]
    fn read_star_columns_find_the_differing_positions_and_assign_the_read() {
        // two copies: a 2-kb random core, copy B differs at 20 positions; the read is 1.2 kb of copy A
        // (positions 300..1500), plus one read that is the reverse complement of copy B's 300..1500.
        let mut seed = 0x9E37_79B9u64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            b"ACGT"[(seed % 4) as usize]
        };
        let a: Vec<u8> = (0..2000).map(|_| rnd()).collect();
        let mut b = a.clone();
        let psv: Vec<usize> = (0..20).map(|k| 350 + k * 55).collect();
        for &p in &psv {
            b[p] = match a[p] { b'A' => b'C', b'C' => b'G', b'G' => b'T', _ => b'A' };
        }
        let read_a = a[300..1500].to_vec();
        let read_b_rc: Vec<u8> = b[300..1500].iter().rev().map(|&x| rc_base(x)).collect();
        let alns = read_star_alignments(&[&a, &b], &[&read_a, &read_b_rc]);
        assert_eq!(alns.len(), 2);
        assert!(alns[0][0].is_some() && alns[0][1].is_some(), "read A must hit both copies: {:?}", alns[0].iter().map(|x| x.is_some()).collect::<Vec<_>>());
        let (obs, prof) = read_star_columns(&read_a, &alns[0]);
        assert_eq!(obs.len(), 20, "20 differing positions inside the read");
        let p = AssignParams::default();
        let asg = assign_read(&ReadFeatures { psv_obs: obs, psv_qual: vec![], junctions: vec![] }, &prof, &p).unwrap();
        assert_eq!((asg.best_copy, asg.status), (0, AssignStatus::Assigned), "{asg:?}");
        let (obs_b, prof_b) = read_star_columns(&read_b_rc, &alns[1]);
        assert_eq!(obs_b.len(), 20);
        let asg_b = assign_read(&ReadFeatures { psv_obs: obs_b, psv_qual: vec![], junctions: vec![] }, &prof_b, &p).unwrap();
        assert_eq!((asg_b.best_copy, asg_b.status), (1, AssignStatus::Assigned), "{asg_b:?}");
    }

    #[test]
    fn read_support_keeps_two_allele_columns_drops_unbacked() {
        // 4 reads (>= the judge-coverage floor, so every column IS judged). min_reads = 2, min_judge_cov = 3.
        // col0: reads split A/G — a real PSV (2 alleles, each >= 2 reads) → KEEP.
        // col1: every read shows C — the paralogs' reference "differs" but no molecule backs a 2nd allele
        //       (an assembly artifact) → DROP.
        // col2: only ONE read shows the minority base (T), below min_reads → DROP.
        let all_obs = vec![
            vec![Some(b'A'), Some(b'C'), Some(b'A')],
            vec![Some(b'A'), Some(b'C'), Some(b'A')],
            vec![Some(b'G'), Some(b'C'), Some(b'A')],
            vec![Some(b'G'), Some(b'C'), Some(b'T')],
        ];
        assert_eq!(read_supported_columns(&all_obs, 3, 2, 3), vec![true, false, false]);
    }

    #[test]
    fn read_support_keeps_low_coverage_columns_it_cannot_judge() {
        // Only 2 reads cover col0 (below min_judge_cov = 4), monomorphic — but we cannot judge it, so KEEP.
        // col1 is well-covered (4 reads) and monomorphic → an artifact → DROP.
        let all_obs = vec![
            vec![Some(b'A'), Some(b'C')],
            vec![Some(b'A'), Some(b'C')],
            vec![None, Some(b'C')],
            vec![None, Some(b'C')],
        ];
        assert_eq!(read_supported_columns(&all_obs, 2, 2, 4), vec![true, false]);
    }

    #[test]
    fn restrict_family_profiles_drops_and_reindexes() {
        let fp = FamilyProfiles {
            profiles: vec![
                CopyProfile { copy_id: 0, alleles: vec![Some(b'A'), Some(b'C'), Some(b'G')], junctions: vec![] },
                CopyProfile { copy_id: 1, alleles: vec![Some(b'T'), Some(b'C'), Some(b'A')], junctions: vec![] },
            ],
            copy_gpos: vec![vec![(0, 100), (1, 200), (2, 300)], vec![(0, 105), (1, 205), (2, 305)]],
            gen2off: vec![Default::default(), Default::default()],
            strand: vec!['+', '+'],
            n_cols: 3,
        };
        let out = restrict_family_profiles(&fp, &[true, false, true]); // drop col 1
        assert_eq!(out.n_cols, 2);
        assert_eq!(out.profiles[0].alleles, vec![Some(b'A'), Some(b'G')]);
        assert_eq!(out.profiles[1].alleles, vec![Some(b'T'), Some(b'A')]);
        assert_eq!(out.copy_gpos[0], vec![(0, 100), (1, 300)]); // col 2 re-indexed to 1
    }

    // exploration harness: compare our direct DP vs poasta on a controlled indel+SNV pair.
    fn aln_summary(tag: &str, msa: &[Vec<u8>]) -> (Vec<usize>, Vec<(usize, usize)>) {
        // returns (ref-offset diff columns, gap runs as (ref_off, len)). Prints a compact view.
        let (r0, r1) = (&msa[0], &msa[1]);
        let (mut ro, mut diffs, mut gaps) = (0usize, Vec::new(), Vec::new());
        let (mut gap_start, mut gap_len) = (0usize, 0usize);
        for c in 0..r0.len() {
            let (ca, cb) = (r0[c], r1[c]);
            if cb == b'-' {
                if gap_len == 0 {
                    gap_start = ro;
                }
                gap_len += 1;
            } else {
                if gap_len > 0 {
                    gaps.push((gap_start, gap_len));
                    gap_len = 0;
                }
                if ca != b'-' && ca != cb {
                    diffs.push(ro);
                }
            }
            if ca != b'-' {
                ro += 1;
            }
        }
        if gap_len > 0 {
            gaps.push((gap_start, gap_len));
        }
        eprintln!("[{tag}] cols={} diffs={} gap_runs={:?} first10diffs={:?}", r0.len(), diffs.len(), gaps, &diffs[..diffs.len().min(10)]);
        (diffs, gaps)
    }

    #[test]
    fn explore_directdp_vs_poasta_on_indel_pair() {
        use poasta::aligner::scoring::GapAffine;
        let a = rand_seq(400, 0xE0F1);
        // B = A with a 50 bp deletion at [150,200) and SNVs at ref offsets 40 and 300.
        let mut b: Vec<u8> = a[0..150].to_vec();
        b.extend_from_slice(&a[200..400]);
        let flip = |x: u8| if x == b'A' { b'C' } else { b'A' };
        // apply SNVs in B coords: offset 40 (shared prefix), offset 250 (= A offset 300 in the suffix)
        b[40] = flip(b[40]);
        b[250] = flip(b[250]);
        let poasta = poa_msa_with_costs(&[a.clone(), b.clone()], GapAffine::new(1, 1, 32)).expect("poasta");
        let banded = banded_msa_pair(&a, &b, 100).expect("banded fits");
        let (pd, pg) = aln_summary("poasta", &poasta);
        let (bd, bg) = aln_summary("direct", &banded);
        eprintln!("SAME diffs? {} | SAME gaps? {}", pd == bd, pg == bg);
        // the SNVs live at A offsets 40 and 300; a clean alignment reports exactly those 2 diffs + one 50bp gap.
        eprintln!("true SNVs at A-offsets [40, 300]; true deletion at A [150,200)");
    }

    // cost of a 2-row MSA under GapAffine(mismatch=1, gap_extend=1, gap_open=32): gap run = 32 + L, mismatch = 1.
    fn msa_pair_cost(msa: &[Vec<u8>]) -> u64 {
        let (a, b) = (&msa[0], &msa[1]);
        let (mut cost, mut in_a, mut in_b) = (0u64, false, false);
        for c in 0..a.len().min(b.len()) {
            let (ga, gb) = (a[c] == b'-', b[c] == b'-');
            if ga {
                cost += if in_a { 1 } else { 33 };
            }
            if gb {
                cost += if in_b { 1 } else { 33 };
            }
            if !ga && !gb && a[c] != b[c] {
                cost += 1;
            }
            in_a = ga;
            in_b = gb;
        }
        cost
    }

    #[test]
    fn banded_dp_never_costs_more_than_poasta_on_divergent_pairs() {
        // Load-bearing invariant: our exact banded DP is a provably optimal affine-gap aligner, so on any pair it
        // must cost NO MORE than poasta's graph aligner. This guards the DP's correctness across divergent,
        // multi-indel, truncated paralog-like pairs. (On the real GGO paralogs poasta is in fact measurably
        // suboptimal — GSTM banded 1181 < poasta 1331, PCDHB 3474 < 4152, both validity-checked — but that
        // reproduction is data-shape-dependent, so here we only assert the guarantee and report any gap.)
        // Both alignments are validated to reconstruct their inputs before their costs are compared.
        use poasta::aligner::scoring::GapAffine;
        let validate = |msa: &[Vec<u8>], a: &[u8], b: &[u8]| -> bool {
            msa.len() == 2 && msa[0].len() == msa[1].len() && {
                let ung = |r: &[u8]| -> Vec<u8> { r.iter().copied().filter(|&c| c != b'-').collect() };
                ung(&msa[0]) == a && ung(&msa[1]) == b
            }
        };
        let flip = |x: u8| match x {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            _ => b'A',
        };
        let mut poasta_suboptimal = 0usize;
        for seed in 0..4u64 {
            // REPETITIVE sequence is the key: real transcripts have low-complexity / repeat content, and that is
            // where poasta's graph search takes a garden path and bounds out. (On clean random sequence poasta is
            // optimal; on tiled-motif sequence — like here — it is not.) Tile 8 short motifs to ~900 bp.
            let motifs: Vec<Vec<u8>> = (0..8).map(|m| rand_seq(60, 0x1000 + seed * 97 + m)).collect();
            let mut r = SplitMix64(seed.wrapping_mul(0xABCD) | 1);
            let mut a: Vec<u8> = Vec::new();
            for _ in 0..15 {
                a.extend_from_slice(&motifs[(r.next_u64() % 8) as usize]);
            }
            // b: ~9% subs + a 90 bp internal deletion + a 240 bp terminal truncation (length disparity + repeats).
            let mut b: Vec<u8> = a[80..a.len() - 240].to_vec();
            b.drain(400..490);
            for pos in b.iter_mut() {
                if r.next_u64() % 100 < 9 {
                    *pos = flip(*pos);
                }
            }
            let band = a.len().abs_diff(b.len()) + 1024;
            let banded = banded_msa_pair(&a, &b, band).expect("banded fits");
            let poasta = poa_msa_with_costs(&[a.clone(), b.clone()], GapAffine::new(1, 1, 32)).expect("poasta");
            assert!(validate(&banded, &a, &b), "banded alignment must reconstruct both inputs");
            assert!(validate(&poasta, &a, &b), "poasta alignment must reconstruct both inputs");
            let (bc, pc) = (msa_pair_cost(&banded), msa_pair_cost(&poasta));
            // THE OPTIMALITY INVARIANT: our exact banded DP is never more expensive than poasta.
            assert!(bc <= pc, "seed {seed}: banded {bc} must be <= poasta {pc} (banded is exact-optimal)");
            if bc < pc {
                poasta_suboptimal += 1;
                eprintln!("seed {seed}: poasta SUBOPTIMAL — banded={bc} < poasta={pc} (Δ={})", pc - bc);
            }
        }
        // On repetitive (transcript-like) sequence poasta is measurably suboptimal — the same failure seen on the
        // real GGO paralogs (GSTM banded 1181 < poasta 1331; PCDHB 3474 < 4152, both validity-checked). The exact
        // banded DP recovers the true optimum every time.
        assert!(
            poasta_suboptimal > 0,
            "expected poasta to be suboptimal on >=1 repetitive divergent pair; got {poasta_suboptimal}/4"
        );
    }

    #[test]
    fn rolling_row_banded_dp_is_byte_identical_to_the_full_matrix_oracle() {
        // §6em: 200 random pairs with substitutions, insertions, deletions and length differences up to 60 bp,
        // bands from tight to loose: the rolling-row/traceback-byte DP must return EXACTLY the oracle's MSA
        // (same co-optimal gap placement) and the same None verdicts.
        let mut seed = 0x9E3779B97F4A7C15u64;
        let mut rnd = || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };
        let mut agree = 0usize;
        for t in 0..200 {
            let n = 60 + (rnd() % 400) as usize;
            let a: Vec<u8> = (0..n).map(|_| b"ACGT"[(rnd() % 4) as usize]).collect();
            let mut b = a.clone();
            let edits = 1 + (rnd() % 12) as usize;
            for _ in 0..edits {
                let pos = (rnd() % (b.len() as u64)) as usize;
                match rnd() % 3 {
                    0 => b[pos] = b"ACGT"[(rnd() % 4) as usize],
                    1 => {
                        let l = 1 + (rnd() % 20) as usize;
                        for _ in 0..l {
                            b.insert(pos, b"ACGT"[(rnd() % 4) as usize]);
                        }
                    }
                    _ => {
                        let l = (1 + (rnd() % 20) as usize).min(b.len() - pos);
                        b.drain(pos..pos + l);
                    }
                }
            }
            let band = a.len().abs_diff(b.len()) + [2usize, 8, 32, 128][t % 4];
            let x = banded_msa_pair(&a, &b, band);
            let y = banded_msa_pair_full_matrix(&a, &b, band);
            assert_eq!(x, y, "pair {t}: n={} m={} band={band}", a.len(), b.len());
            if x.is_some() {
                agree += 1;
            }
        }
        assert!(agree > 100, "most pairs must align within their band ({agree}/200)");
    }

    #[test]
    fn banded_msa_pair_aligns_unambiguous_pair_and_falls_back_on_truncation() {
        // Near-identical, same length, one mismatch -> unambiguous diagonal alignment, no gaps, the diff column
        // recovered exactly. (This is the regime where the banded DP is byte-identical to poasta.)
        let a = b"ACGTACGTACGT".to_vec();
        let mut b = a.clone();
        b[5] = b'A'; // single mismatch at offset 5 (C->A)
        let msa = banded_msa_pair(&a, &b, 4).expect("similar-length pair fits the band");
        assert_eq!(msa.len(), 2);
        assert_eq!(msa[0], a, "no gaps in the aligned ref");
        assert_eq!(msa[1], b, "no gaps in the aligned other");
        // the only differing column is offset 5.
        let diffs: Vec<usize> = (0..msa[0].len()).filter(|&c| msa[0][c] != msa[1][c]).collect();
        assert_eq!(diffs, vec![5]);
        // a truncation (length difference > band) cannot stay in the band -> None -> caller uses exact poasta.
        assert!(banded_msa_pair(&a, &a[..3], 4).is_none(), "length diff 9 > band 4 -> fall back");
    }

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
         ..Default::default() }
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

    #[test]
    fn detect_editing_columns_flags_editing_not_real_psv() {
        // 2 copies, 2 A↔G columns. col0 = a REAL PSV (copy0=A monomorphic, copy1=G monomorphic).
        // col1 = an EDITING site: copy0-assigned reads carry A but ~30% show an edited G.
        let copies = vec![
            CopyProfile { copy_id: 0, alleles: vec![Some(b'A'), Some(b'A')], junctions: vec![] },
            CopyProfile { copy_id: 1, alleles: vec![Some(b'G'), Some(b'G')], junctions: vec![] },
        ];
        let mut reads_obs: Vec<Vec<Option<u8>>> = Vec::new();
        for _ in 0..7 {
            reads_obs.push(vec![Some(b'A'), Some(b'A')]); // copy0, unedited
        }
        for _ in 0..3 {
            reads_obs.push(vec![Some(b'A'), Some(b'G')]); // copy0, edited at col1
        }
        for _ in 0..10 {
            reads_obs.push(vec![Some(b'G'), Some(b'G')]); // copy1
        }
        let flag = detect_editing_columns(&reads_obs, &copies);
        assert_eq!(flag, vec![false, true], "col0 real PSV not flagged; col1 editing flagged");
    }

    #[test]
    fn detect_editing_columns_ignores_sequencing_error() {
        // a real PSV where a couple of error G's appear in copy0's reads must NOT be flagged (frac < 0.05).
        let copies = vec![
            CopyProfile { copy_id: 0, alleles: vec![Some(b'A')], junctions: vec![] },
            CopyProfile { copy_id: 1, alleles: vec![Some(b'G')], junctions: vec![] },
        ];
        let mut reads_obs: Vec<Vec<Option<u8>>> = Vec::new();
        for _ in 0..98 {
            reads_obs.push(vec![Some(b'A')]); // copy0
        }
        for _ in 0..2 {
            reads_obs.push(vec![Some(b'G')]); // 2 error G's among 100 copy0 reads = 2% < EDIT_MIN_FRAC (5%)
        }
        for _ in 0..50 {
            reads_obs.push(vec![Some(b'G')]); // copy1
        }
        // the 2 error-G reads argmax to copy1 (they match G), so copy0's reads stay monomorphic A -> not flagged.
        assert_eq!(detect_editing_columns(&reads_obs, &copies), vec![false]);
    }

    #[test]
    fn fill_psv_obs_carries_per_base_quality() {
        // read covers ref 10..20 with per-base QVs; a PSV column at genome pos 13 -> col 0.
        let read = AlignedRead {
            ref_start: 10,
            cigar: vec![('M', 10)],
            seq: b"ACGTACGTAC".to_vec(),
            qual: vec![20, 21, 22, 40, 24, 25, 26, 27, 28, 29], // QV 40 at offset 3 (ref 13)
        };
        let gpos = vec![(0usize, 13u64)];
        let mut obs = vec![None];
        let mut q = vec![None];
        fill_psv_obs(&read, &gpos, false, &mut obs, &mut q);
        assert_eq!(obs[0], Some(b'T'), "base at ref 13 is seq[3]='T'");
        assert_eq!(q[0], Some(40), "QV travels with the base (qual[3]=40)");
        // empty qual -> None (flat-rate fallback), base still resolved
        let read2 = AlignedRead { qual: vec![], ..read.clone() };
        let mut obs2 = vec![None];
        let mut q2 = vec![None];
        fill_psv_obs(&read2, &gpos, false, &mut obs2, &mut q2);
        assert_eq!(obs2[0], Some(b'T'));
        assert_eq!(q2[0], None);
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

    #[test]
    fn cigar_to_gapped_msa_reconstructs_alignment() {
        // CIGAR-walk reconstruction (no minimap2 binary). ref has a 3bp insert vs query and one mismatch.
        // ref:  A C G T A C G T   (8)        cigar (target-relative): 4= 3D 1X
        // query:A C G T  - - - C   -> wait: build a clear case below.
        let r = b"ACGTACGT".to_vec(); // 8
        let q = b"ACGTC".to_vec(); // 5: matches ref[0:4], deletes ref[4:7], mismatch at ref[7] (T->C)
        // ts=0 qs=0 cigar: 4= (ACGT) 3D (ref ACG, query gap) 1X (ref T vs query C)
        let msa = cigar_to_gapped_msa(&r, &q, 0, 0, "4=3D1X");
        assert_eq!(msa.len(), 2);
        // both rows equal length; ref row has no gaps (global over ref), query row has the 3 deletion gaps.
        assert_eq!(msa[0].len(), msa[1].len());
        assert_eq!(msa[0], b"ACGTACGT", "ref row = full ref, no gaps");
        assert_eq!(msa[1], b"ACGT---C", "query row: matched core, 3 gaps for the deletion, mismatch C at end");
        // walking both-non-gap columns recovers exactly the mismatch at ref offset 7.
        let diff: Vec<usize> = (0..msa[0].len())
            .filter(|&c| msa[0][c] != b'-' && msa[1][c] != b'-' && msa[0][c] != msa[1][c])
            .collect();
        // column 7 in the alignment == ref offset 7 (no insertions before it)
        assert_eq!(diff, vec![7], "the lone PSV is at ref offset 7");
    }

    #[test]
    fn cigar_to_gapped_msa_handles_clipped_ends() {
        // ts=2, qs=1: ref[0:2] and query[0:1] are unaligned ends -> rendered as gap blocks, offsets stay correct.
        let r = b"GGACGT".to_vec(); // ref, aligned core starts at ts=2 (ACGT)
        let q = b"TACGT".to_vec(); // query, aligned core starts at qs=1 (ACGT)
        let msa = cigar_to_gapped_msa(&r, &q, 2, 1, "4=");
        assert_eq!(msa[0].len(), msa[1].len());
        // ref row contains all of ref (no ref base lost); query row contains all of query.
        assert_eq!(msa[0].iter().filter(|&&b| b != b'-').count(), r.len());
        assert_eq!(msa[1].iter().filter(|&&b| b != b'-').count(), q.len());
        // no spurious PSV: the only both-non-gap columns are the aligned ACGT core (all matches).
        let diff = (0..msa[0].len())
            .filter(|&c| msa[0][c] != b'-' && msa[1][c] != b'-' && msa[0][c] != msa[1][c])
            .count();
        assert_eq!(diff, 0, "clipped divergent ends produce NO spurious PSV");
    }

    /// minimap2 PSV-discovery must find the SAME columns as poasta. Needs the minimap2 binary + long
    /// enough seqs to seed (asm20), so it is ignored by default; run with `--ignored` where minimap2 is on
    /// PATH (or RUSTLE_MINIMAP2 set). Validates the gapped-MSA reconstruction from the minimap2 CIGAR.
    #[test]
    #[ignore = "needs the minimap2 binary"]
    fn minimap2_psv_discovery_matches_poasta() {
        // ~3 kb near-identical copies with 4 planted substitutions (asm20 seeds on this length).
        let sa = rand_seq(3000, 0x9E1);
        let mut sb = sa.clone();
        let flips = [300usize, 1100, 1900, 2600];
        for &p in &flips {
            sb[p] = match sa[p] { b'A' => b'C', b'C' => b'A', b'G' => b'T', _ => b'G' };
        }
        let ca = copy_tx("A", 0, 3000, '+', &[], sa);
        let cb = copy_tx("B", 10000, 13000, '+', &[], sb);
        let copies = [&ca, &cb];
        let maps = [exon_map(&ca), exon_map(&cb)];
        let poasta = discover_psvs(&copies, &maps);
        std::env::set_var("RUSTLE_PSV_MINIMAP2", "1");
        let mm2 = discover_psvs(&copies, &maps);
        std::env::remove_var("RUSTLE_PSV_MINIMAP2");
        let off = |cols: &[PsvColumn]| -> Vec<u64> {
            cols.iter().filter_map(|c| c[0].map(|(g, _)| g)).collect()
        };
        assert_eq!(off(&poasta), off(&mm2), "minimap2 finds the SAME PSV ref-genome positions as poasta");
        assert_eq!(mm2.len(), flips.len(), "all 4 planted substitutions found");
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
        let read = AlignedRead { ref_start: 0, cigar: vec![('M', 300)], seq: sb, qual: vec![] };
        let res = assign_family(&copies, &[read], &AssignParams::default(), None);
        assert_eq!(res.len(), 1);
        let (ri, a) = &res[0];
        assert_eq!(*ri, 0);
        assert_eq!(a.best_copy, 1, "carries copyB alleles -> copyB, despite aligning to copyA's region");
        assert_eq!(a.status, AssignStatus::Assigned);
    }

    #[test]
    fn synthetic_ground_truth_assigns_every_read_to_its_true_copy() {
        // Self-contained (NO external fixture) ground-truth panel for the O2 accuracy claim — the
        // always-on, data-free complement to `smoke_sim5x_ground_truth` (which no-ops without
        // RUSTLE_SIM5X_DIR). Three single-exon copies share a base sequence and differ at 5 planted PSV
        // columns, a distinct private allele per copy, so every pair is separated by K=5 >> 2
        // distinguishing sites. Each copy emits reads carrying ITS OWN alleles at ITS OWN locus, tagged
        // with the true copy index. The IsoCon certificate must assign every read to its true copy with
        // ZERO misassignment (assign-or-abstain, never 1/k).
        const N_COPIES: usize = 3;
        const READS_PER_COPY: usize = 8;
        let psv = [50usize, 100, 150, 200, 250];
        let base = rand_seq(300, 0x6D0D);
        let alleles = [b'A', b'C', b'G']; // one private allele per copy at every PSV column
        let seqs: Vec<Vec<u8>> = (0..N_COPIES)
            .map(|c| {
                let mut s = base.clone();
                for &p in &psv {
                    s[p] = alleles[c];
                }
                s
            })
            .collect();
        let copies_owned: Vec<DenovoTranscript> = (0..N_COPIES)
            .map(|c| {
                let start = c as u64 * 1000;
                copy_tx(&format!("copy{c}"), start, start + 300, '+', &[], seqs[c].clone())
            })
            .collect();
        let copies: Vec<&DenovoTranscript> = copies_owned.iter().collect();
        // Each copy's reads are aligned to THAT copy's locus and carry THAT copy's sequence.
        let mut reads = Vec::new();
        let mut truth = Vec::new();
        for c in 0..N_COPIES {
            for _ in 0..READS_PER_COPY {
                reads.push(AlignedRead {
                    ref_start: c as u64 * 1000,
                    cigar: vec![('M', 300)],
                    seq: seqs[c].clone(),
                    qual: vec![],
                });
                truth.push(c);
            }
        }
        let res = assign_family(&copies, &reads, &AssignParams::default(), None);
        assert_eq!(res.len(), reads.len());
        let (mut assigned, mut misassigned) = (0usize, 0usize);
        for (ri, a) in &res {
            if a.status == AssignStatus::Assigned {
                assigned += 1;
                misassigned += (a.best_copy != truth[*ri]) as usize;
            }
        }
        assert_eq!(misassigned, 0, "no read may be assigned to the WRONG copy (ground truth)");
        assert_eq!(assigned, reads.len(), "every K=5 read spans 5 private PSVs -> all assignable");
    }

    #[test]
    fn synthetic_ground_truth_identical_copies_all_tie() {
        // The identifiability floor (sim5x K=0, data-free): two EXONICALLY IDENTICAL copies carry no PSV,
        // so no read can be resolved even in principle -- every read must be certified Tied, never guessed.
        let seq = rand_seq(300, 0xF10A7);
        let ca = copy_tx("A", 0, 300, '+', &[], seq.clone());
        let cb = copy_tx("B", 1000, 1300, '+', &[], seq.clone());
        let copies = [&ca, &cb];
        let reads: Vec<AlignedRead> = (0..12)
            .map(|_| AlignedRead { ref_start: 0, cigar: vec![('M', 300)], seq: seq.clone(), qual: vec![] })
            .collect();
        let res = assign_family(&copies, &reads, &AssignParams::default(), None);
        assert_eq!(res.len(), reads.len());
        assert!(
            res.iter().all(|(_, a)| a.status == AssignStatus::Tied),
            "identical copies (0 PSVs) => every read Tied, none assigned"
        );
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
        let read = AlignedRead { ref_start: 0, cigar: vec![('M', 100), ('N', 100), ('M', 200)], seq: spliced, qual: vec![] };
        // junction_err=1e-4 < alpha=1e-3: a single copy-specific junction resolves under the
        // significance gate (p_read=1e-4 < 1e-3, margin >> 0 -> Assigned).
        let detail = assign_family_detailed(&copies, &[read], &AssignParams::default(), None, None);
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
        let read = AlignedRead { ref_start: 0, cigar: vec![('M', 150)], seq: base[0..150].to_vec(), qual: vec![] };
        let res = assign_family(&copies, &[read], &AssignParams::default(), None);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].1.status, AssignStatus::Tied, "spans no decisive feature");
    }

    #[test]
    fn abundance_ci_is_full_simplex_when_unidentifiable() {
        // L8: when NO read carries a decisive feature (K=0 regime), the per-copy abundance is
        // unidentifiable; the CI must be the full-simplex half-width (0.5), NOT a sqrt(N)-shrunk false
        // precision. Two copies differing only at offsets >= 200; all reads cover only [0,150) -> tied.
        let base = rand_seq(300, 0xC10);
        let mut sa = base.clone();
        let mut sb = base.clone();
        for &p in &[210usize, 250, 290] {
            sa[p] = b'A';
            sb[p] = b'C';
        }
        let ca = copy_tx("A", 0, 300, '+', &[], sa);
        let cb = copy_tx("B", 1000, 1300, '+', &[], sb);
        let copies = [&ca, &cb];
        // 30 reads, all spanning no decisive feature -> n_eff = 0 despite N=30.
        let reads: Vec<AlignedRead> = (0..30)
            .map(|_| AlignedRead { ref_start: 0, cigar: vec![('M', 150)], seq: base[0..150].to_vec(), qual: vec![] })
            .collect();
        let detail = assign_family_detailed(&copies, &reads, &AssignParams::default(), None, None);
        assert!(
            detail.copy_abundance_ci.iter().all(|&ci| (ci - 0.5).abs() < 1e-9),
            "unidentifiable (no decisive reads) -> CI = 0.5 full-simplex, got {:?}",
            detail.copy_abundance_ci
        );
    }

    /// Ground-truth validation on the sim5x 5-copy synthetic dataset (the python `copy_assign.py sim5x`
    /// oracle): each read name encodes its TRUE copy, so we can report assignment accuracy. K = PSVs per
    /// copy (the identifiability ladder): K=0 has no PSVs → all tied; K>=2 → ~100% accurate. Ignored by
    /// default; run in release with RUSTLE_SIM5X_DIR set.
    // PINNED in CI (P3): runs in the normal suite. Early-returns (passes) when RUSTLE_SIM5X_DIR is unset, so
    // it is a no-op without the data; when the sim5x dir IS present it ASSERTS the identifiability ladder
    // invariants (K=0 -> 100% tied; K>=2 -> acc|assigned >= 0.99). This is the non-circular accuracy point.
    #[test]
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
                     ..Default::default() }
                })
                .collect();
            let copies: Vec<&DenovoTranscript> = copies_owned.iter().collect();
            let fp = build_family_profiles(&copies, None);
            let reads = aligned_reads_from_bam(&bam, 4).expect("bam");
            let names: Vec<String> = reads.iter().map(|br| br.name.clone()).collect();
            let ars: Vec<AlignedRead> = reads.into_iter().map(|br| br.read).collect();
            let assigns = assign_family(&copies, &ars, &AssignParams::default(), None);
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
            // PINNED INVARIANTS (the identifiability ladder; assert when sim5x data is present):
            //   K=0  -> exonically identical, no PSVs => 100% Tied, 0 assigned.
            //   K>=2 -> reads span >=2 distinguishing PSVs => near-perfect accuracy on the ASSIGNED reads.
            if n > 0 {
                if k == 0 {
                    assert_eq!(resolvable, 0, "sim5x K=0 must be 100% tied (no PSVs); got {resolvable} resolvable / {n}");
                    assert_eq!(assigned, 0, "sim5x K=0 must assign 0 reads; got {assigned}");
                } else if k >= 2 && assigned > 0 {
                    let acc_assigned = corr_assigned as f64 / assigned as f64;
                    assert!(
                        acc_assigned >= 0.99,
                        "sim5x K={k}: acc|assigned {acc_assigned:.4} < 0.99 ({corr_assigned}/{assigned} correct)"
                    );
                }
            }
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
        let read = AlignedRead { ref_start: 1000, cigar: vec![('M', 240)], seq: fwd_b, qual: vec![] };
        let res = assign_family(&copies, &[read], &AssignParams::default(), None);
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
        let read = AlignedRead { ref_start: 0, cigar: vec![('M', 300)], seq: recomb, qual: vec![] };
        let detail = assign_family_detailed(&[&ca, &cb], &[read], &AssignParams::default(), None, None);
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
        let detail = assign_family_detailed(&[&ca, &cb, &cc], &[], &AssignParams::default(), None, None);
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

    // ---- freeze_merge (two-stage freeze) ----

    fn asg(best: usize, status: AssignStatus) -> Assignment {
        Assignment {
            best_copy: best,
            log_lr_margin: 0.0,
            n_decisive: 1,
            resolvable: status != AssignStatus::Tied,
            status,
            p_value: 0.0,
            min_p_value: 0.0,
            discovery_coupled: false,
            junction_conflict: false,
                        origin_rejected: false,
            posterior: vec![],
        }
    }
    fn rr(read_index: usize, mapped_copy: usize, combined: Assignment) -> ReadResult {
        ReadResult { read_index, mapped_copy, psv: combined.clone(), combined, psv_obs: vec![], junctions: vec![] }
    }
    /// fetch the merged result for a given read_index (matching is by read_index, not position).
    fn by_idx(v: &[ReadResult], idx: usize) -> &ReadResult {
        v.iter().find(|r| r.read_index == idx).expect("read_index present")
    }

    #[test]
    fn freeze_keeps_stage1_assigned_at_multi_ref() {
        // n_ref=2. read 0 was Stage-1 Assigned to copy 0; Stage-2 (with an absent copy) would move it to 1.
        // Multi-ref => Stage-1 wins (frozen): the WHOLE assignment side (best_copy AND mapped_copy) is restored
        // to Stage-1's, NOT discovery_coupled. read 1 was Stage-1 Tied => non-frozen, keeps Stage-2's
        // mapped_copy. Stage-1/Stage-2 carry DIFFERENT mapped_copy values to prove the freeze.
        let stage1 = vec![
            rr(0, 0, asg(0, AssignStatus::Assigned)),
            rr(1, 0, asg(0, AssignStatus::Tied)),
        ];
        let stage2 = vec![
            rr(0, 1, asg(1, AssignStatus::Assigned)),
            rr(1, 2, asg(2, AssignStatus::Assigned)),
        ];
        let merged = freeze_merge(&stage1, stage2, 2);
        let r0 = by_idx(&merged, 0);
        assert_eq!(r0.combined.best_copy, 0, "Stage-1 assignment frozen");
        assert_eq!(r0.mapped_copy, 0, "frozen read also restores Stage-1 mapped_copy");
        assert!(!r0.combined.discovery_coupled);
        let r1 = by_idx(&merged, 1);
        assert_eq!(r1.mapped_copy, 2, "non-frozen read keeps Stage-2 mapped_copy");
        assert!(r1.combined.discovery_coupled, "non-frozen read assigned to absent copy 2 is coupled");
    }

    #[test]
    fn freeze_tied_read_takes_stage2_absent_and_flags_coupled() {
        // n_ref=2. read 1 was Stage-1 Tied; Stage-2 resolves it to absent copy index 2 (>= n_ref).
        // => take Stage-2, flag discovery_coupled.
        let stage1 = vec![rr(1, 0, asg(0, AssignStatus::Tied))];
        let stage2 = vec![rr(1, 2, asg(2, AssignStatus::Assigned))];
        let merged = freeze_merge(&stage1, stage2, 2);
        let r = by_idx(&merged, 1);
        assert_eq!(r.combined.best_copy, 2, "Stage-2 absent assignment kept");
        assert!(r.combined.discovery_coupled, "assigned to absent copy => coupled");
    }

    #[test]
    fn freeze_stage2_to_ref_copy_is_not_coupled() {
        // Stage-1 Ambiguous read resolves under Stage-2 to a REF copy (best_copy 1 < n_ref=2) => not coupled.
        let stage1 = vec![rr(3, 0, asg(0, AssignStatus::Ambiguous))];
        let stage2 = vec![rr(3, 1, asg(1, AssignStatus::Assigned))];
        let merged = freeze_merge(&stage1, stage2, 2);
        let r = by_idx(&merged, 3);
        assert_eq!(r.combined.best_copy, 1);
        assert!(!r.combined.discovery_coupled, "ref-copy assignment is never discovery_coupled");
    }

    #[test]
    fn freeze_single_ref_takes_all_stage2() {
        // n_ref=1 (single ref copy): nothing to freeze, even a Stage-1-Assigned read takes Stage-2.
        let stage1 = vec![rr(0, 0, asg(0, AssignStatus::Assigned))];
        let stage2 = vec![rr(0, 1, asg(1, AssignStatus::Assigned))];
        let merged = freeze_merge(&stage1, stage2, 1);
        let r = by_idx(&merged, 0);
        assert_eq!(r.combined.best_copy, 1, "single_ref => Stage-2 wins");
        assert!(r.combined.discovery_coupled, "best_copy 1 >= n_ref 1 => coupled");
    }

    #[test]
    fn freeze_matches_by_read_index_not_position() {
        // Stage-1 and Stage-2 disagree on membership AND order, to prove matching is by read_index:
        //   stage1 = [read 0 Assigned->0, read 3 Assigned->1]   (read 3 only exists in Stage-1)
        //   stage2 = [read 2 ->0, read 0 ->1, read 1 Tied]      (read 0 sits at position 1)
        // read 0 must freeze to Stage-1's best_copy 0 (by index), NOT to position-1's Stage-1 read 3 (best 1).
        let stage1 = vec![
            rr(0, 0, asg(0, AssignStatus::Assigned)),
            rr(3, 1, asg(1, AssignStatus::Assigned)),
        ];
        let stage2 = vec![
            rr(2, 0, asg(0, AssignStatus::Assigned)),
            rr(0, 1, asg(1, AssignStatus::Assigned)),
            rr(1, 0, asg(0, AssignStatus::Tied)),
        ];
        let merged = freeze_merge(&stage1, stage2, 2);
        // Stage-2 order/membership preserved (read 3 from Stage-1 is dropped; Stage-2 is the base).
        assert_eq!(merged.iter().map(|r| r.read_index).collect::<Vec<_>>(), vec![2, 0, 1]);
        // read 0 frozen to Stage-1 best_copy 0 (index match), not position-1's read-3 value (1).
        assert_eq!(by_idx(&merged, 0).combined.best_copy, 0);
        // read 2 has no Stage-1 entry => keeps Stage-2 (best 0, ref copy, not coupled).
        let r2 = by_idx(&merged, 2);
        assert_eq!(r2.combined.best_copy, 0);
        assert!(!r2.combined.discovery_coupled);
        // read 1 has no Stage-1 entry and Stage-2 left it Tied => unchanged Tied.
        assert_eq!(by_idx(&merged, 1).combined.status, AssignStatus::Tied);
    }

    #[test]
    fn intron_psv_finds_a_divergent_intron_column_when_exons_are_identical() {
        // The intron-retention lever: two copies with IDENTICAL exons but a DIVERGENT intron. discover_psvs
        // (exon-only) finds nothing; discover_intron_psvs finds the intronic distinguishing column.
        use crate::genome::GenomeIndex;
        let mm2 = std::env::var("RUSTLE_MINIMAP2").unwrap_or_else(|_| "minimap2".into());
        let ok = std::process::Command::new(&mm2).arg("--version").output().map(|o| o.status.success());
        if !matches!(ok, Ok(true)) {
            return; // minimap2 not runnable -> skip (discover_intron_psvs aligns the genomic spans with it)
        }
        // deterministic pseudo-random ACGT (non-repetitive, so minimap2 seeds the 200 bp spans cleanly)
        let rseq = |n: usize, seed: u64| -> Vec<u8> {
            let mut x = seed;
            (0..n)
                .map(|_| {
                    x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    b"ACGT"[((x >> 33) & 3) as usize]
                })
                .collect()
        };
        let exon1 = rseq(60, 1);
        let exon2 = rseq(40, 2);
        let mut intron_a = rseq(100, 3);
        let mut intron_b = intron_a.clone();
        intron_a[50] = b'A';
        intron_b[50] = b'C'; // the single divergent intronic position (span offset 60+50 = 110)
        // contig = A-span(200) + filler(20) + B-span(200)
        let mut contig = Vec::new();
        contig.extend_from_slice(&exon1);
        contig.extend_from_slice(&intron_a);
        contig.extend_from_slice(&exon2); // A span [0,200)
        contig.extend_from_slice(&rseq(20, 9));
        let b_start = contig.len() as u64; // 220
        contig.extend_from_slice(&exon1);
        contig.extend_from_slice(&intron_b);
        contig.extend_from_slice(&exon2); // B span [220,420)
        let fa = std::env::temp_dir().join(format!("rustle_intronpsv_{}.fa", std::process::id()));
        std::fs::write(&fa, format!(">c1\n{}\n", String::from_utf8(contig).unwrap())).unwrap();
        let genome = GenomeIndex::from_fasta(fa.to_str().unwrap()).unwrap();
        let spliced: Vec<u8> = exon1.iter().chain(exon2.iter()).cloned().collect(); // identical for both copies
        let mk = |start: u64| DenovoTranscript {
            tid: format!("c{start}"),
            chrom: "c1".into(),
            start,
            end: start + 200,
            n_reads: 0,
            strand: '+',
            introns: vec![(start + 60, start + 160)],
            seq: spliced.clone(),
         ..Default::default() };
        let a = mk(0);
        let b = mk(b_start);
        let copies = [&a, &b];
        // exon-only PSVs: none (exons identical)
        let exon_maps: Vec<Vec<u64>> = copies.iter().map(|c| exon_map(c)).collect();
        assert!(discover_psvs(&copies, &exon_maps).is_empty(), "exons are identical -> no exonic PSV");
        // intronic PSV: exactly one, at A:110 / B:330, with differing bases
        let icols = discover_intron_psvs(&copies, &genome);
        let _ = std::fs::remove_file(&fa);
        assert_eq!(icols.len(), 1, "one divergent intron column expected");
        let col = &icols[0];
        assert_eq!(col[0].unwrap().0, 110, "copy A intronic PSV genomic pos");
        assert_eq!(col[1].unwrap().0, 330, "copy B intronic PSV genomic pos");
        assert_ne!(col[0].unwrap().1, col[1].unwrap().1, "the intronic bases differ (A vs G)");
    }

    #[test]
    fn freeze_overwrites_psv_too_when_frozen() {
        // freezing replaces BOTH combined and psv with Stage-1's.
        let mut s1_combined = asg(0, AssignStatus::Assigned);
        s1_combined.log_lr_margin = 7.5; // a Stage-1 signature to detect
        let mut s1 = rr(5, 0, s1_combined);
        s1.psv = asg(0, AssignStatus::Assigned);
        s1.psv.log_lr_margin = 3.3;
        let stage2 = vec![rr(5, 1, asg(1, AssignStatus::Assigned))];
        let merged = freeze_merge(&[s1], stage2, 2);
        let r = by_idx(&merged, 5);
        assert_eq!(r.combined.log_lr_margin, 7.5, "combined frozen from Stage-1");
        assert_eq!(r.psv.log_lr_margin, 3.3, "psv frozen from Stage-1");
    }

    // ---- iterative copy pruning ----

    fn seq_with(base: &[u8], pos_allele: &[(usize, u8)]) -> Vec<u8> {
        let mut s = base.to_vec();
        for &(p, a) in pos_allele {
            s[p] = a;
        }
        s
    }

    fn q40_read(ref_start: u64, cigar: Vec<(char, u64)>, seq: Vec<u8>) -> AlignedRead {
        let qual = vec![40; seq.len()];
        AlignedRead { ref_start, cigar, seq, qual }
    }

    #[test]
    fn iterative_prune_merges_duplicate_copy() {
        // Three single-exon copies. c0 differs from c1/c2 at offset 50; c1 and c2 are identical.
        // Q40 reads carrying each allele are aligned to each copy's region. c2 has no evidence that
        // significantly distinguishes it from c1, so it should be merged into c1.
        let base = rand_seq(100, 0xC0D0);
        let c0_seq = seq_with(&base, &[(50, b'A')]);
        let c1_seq = seq_with(&base, &[(50, b'C')]);
        let c2_seq = c1_seq.clone();
        let c0 = copy_tx("c0", 0, 100, '+', &[], c0_seq);
        let c1 = copy_tx("c1", 1000, 1100, '+', &[], c1_seq);
        let c2 = copy_tx("c2", 2000, 2100, '+', &[], c2_seq);
        let copies = [&c0, &c1, &c2];
        let r0 = q40_read(0, vec![('M', 100)], seq_with(&base, &[(50, b'A')]));
        let r1 = q40_read(1000, vec![('M', 100)], seq_with(&base, &[(50, b'C')]));
        let r2 = q40_read(2000, vec![('M', 100)], seq_with(&base, &[(50, b'C')]));
        let reads = [r0, r1, r2];
        let mut p = AssignParams::default();
        p.iterative_prune = true;
        let detail = assign_family_detailed_pruned(&copies, &reads, &p, None, None);
        assert_eq!(detail.copy_psv_alleles.len(), 2, "duplicate c2 merged -> 2 surviving copies");
        // The surviving copies should be c0 and c1 (c2 removed, not c0).
        assert!(detail.results.iter().all(|r| r.combined.best_copy < 2));
    }

    #[test]
    fn iterative_prune_keeps_distinct_supported_copies() {
        // Three copies with a symmetric 3-column PSV design: every pair differs at two columns,
        // and each copy's read matches its own allele at all three columns. No copy should be merged.
        let base = rand_seq(200, 0xD150);
        let c0_seq = seq_with(&base, &[(50, b'A'), (100, b'C'), (150, b'C')]);
        let c1_seq = seq_with(&base, &[(50, b'C'), (100, b'A'), (150, b'C')]);
        let c2_seq = seq_with(&base, &[(50, b'C'), (100, b'C'), (150, b'A')]);
        let c0 = copy_tx("c0", 0, 200, '+', &[], c0_seq);
        let c1 = copy_tx("c1", 1000, 1200, '+', &[], c1_seq);
        let c2 = copy_tx("c2", 2000, 2200, '+', &[], c2_seq);
        let copies = [&c0, &c1, &c2];
        let r0 = q40_read(0, vec![('M', 200)], seq_with(&base, &[(50, b'A'), (100, b'C'), (150, b'C')]));
        let r1 = q40_read(1000, vec![('M', 200)], seq_with(&base, &[(50, b'C'), (100, b'A'), (150, b'C')]));
        let r2 = q40_read(2000, vec![('M', 200)], seq_with(&base, &[(50, b'C'), (100, b'C'), (150, b'A')]));
        let reads = [r0, r1, r2];
        let mut p = AssignParams::default();
        p.iterative_prune = true;
        let detail = assign_family_detailed_pruned(&copies, &reads, &p, None, None);
        assert_eq!(detail.copy_psv_alleles.len(), 3, "all distinct supported copies kept");
    }

    #[test]
    fn iterative_prune_default_off_is_unchanged() {
        // Same setup as iterative_prune_merges_duplicate_copy, but with default params (prune off).
        let base = rand_seq(100, 0xC0D0);
        let c0_seq = seq_with(&base, &[(50, b'A')]);
        let c1_seq = seq_with(&base, &[(50, b'C')]);
        let c2_seq = c1_seq.clone();
        let c0 = copy_tx("c0", 0, 100, '+', &[], c0_seq);
        let c1 = copy_tx("c1", 1000, 1100, '+', &[], c1_seq);
        let c2 = copy_tx("c2", 2000, 2100, '+', &[], c2_seq);
        let copies = [&c0, &c1, &c2];
        let r0 = q40_read(0, vec![('M', 100)], seq_with(&base, &[(50, b'A')]));
        let r1 = q40_read(1000, vec![('M', 100)], seq_with(&base, &[(50, b'C')]));
        let r2 = q40_read(2000, vec![('M', 100)], seq_with(&base, &[(50, b'C')]));
        let reads = [r0, r1, r2];
        let detail = assign_family_detailed(&copies, &reads, &AssignParams::default(), None, None);
        assert_eq!(detail.copy_psv_alleles.len(), 3, "default (prune off) keeps all 3 copies");
    }
}

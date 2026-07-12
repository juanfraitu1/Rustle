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

use rayon::prelude::*;

use super::copy_assign::{
    assign_read, assign_read_editing, boundary_present, copy_pair_significance, AssignParams, AssignStatus,
    Assignment, CopyProfile, ReadFeatures,
};
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
fn banded_msa_pair(a: &[u8], b: &[u8], band: usize) -> Option<Vec<Vec<u8>>> {
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
        if std::env::var_os("RUSTLE_PSV_BAND_DEBUG").is_some() {
            eprintln!(
                "[psv-band] copy{other} len {} vs {}: {}",
                ref_seq.len(),
                copies[other].seq.len(),
                if banded.is_some() { "BANDED" } else { "fell back to exact" }
            );
        }
        // Diagnostic (RUSTLE_PSV_COST_AUDIT=1): per-pair alignment cost + structure for the banded DP vs poasta,
        // each validated to reconstruct both inputs. This is how poasta was shown to be SUBOPTIMAL on divergent
        // GGO paralogs (GSTM copy1 banded 1181 < poasta 1331; PCDHB copy2 3474 < 4152) — the exact banded DP finds
        // a strictly cheaper VALID alignment. Cost model = GapAffine(mismatch=1, gap_extend=1, gap_open=32): a gap
        // run of length L costs 32 + L, a mismatch 1.
        if std::env::var_os("RUSTLE_PSV_COST_AUDIT").is_some() {
            // returns (total_cost, mismatches, gap_runs, gap_bases, lead_gap, trail_gap)
            let msa_cost = |msa: &[Vec<u8>]| -> (u64, u64, u64, u64, u64, u64) {
                if msa.len() != 2 {
                    return (u64::MAX, 0, 0, 0, 0, 0);
                }
                let (a, b) = (&msa[0], &msa[1]);
                let (mut cost, mut mism, mut runs, mut glen) = (0u64, 0u64, 0u64, 0u64);
                let (mut in_gap_a, mut in_gap_b) = (false, false);
                let cols = a.len().min(b.len());
                for c in 0..cols {
                    let (ga, gb) = (a[c] == b'-', b[c] == b'-');
                    if ga {
                        cost += if in_gap_a { 1 } else { 33 };
                        glen += 1;
                        if !in_gap_a {
                            runs += 1;
                        }
                    }
                    if gb {
                        cost += if in_gap_b { 1 } else { 33 };
                        glen += 1;
                        if !in_gap_b {
                            runs += 1;
                        }
                    }
                    if !ga && !gb && a[c] != b[c] {
                        cost += 1;
                        mism += 1;
                    }
                    in_gap_a = ga;
                    in_gap_b = gb;
                }
                // leading / trailing gap length (either row) — flags free-end-gap behavior
                let lead = (0..cols).take_while(|&c| a[c] == b'-' || b[c] == b'-').count() as u64;
                let trail = (0..cols).rev().take_while(|&c| a[c] == b'-' || b[c] == b'-').count() as u64;
                (cost, mism, runs, glen, lead, trail)
            };
            // valid iff each row, with gaps removed, exactly reconstructs its input sequence.
            let validate = |msa: &[Vec<u8>]| -> bool {
                if msa.len() != 2 {
                    return false;
                }
                let ungap = |r: &[u8]| -> Vec<u8> { r.iter().copied().filter(|&c| c != b'-').collect() };
                msa[0].len() == msa[1].len()
                    && ungap(&msa[0]) == *ref_seq
                    && ungap(&msa[1]) == copies[other].seq
            };
            let banded_m = banded_msa_pair(
                ref_seq,
                &copies[other].seq,
                ref_seq.len().abs_diff(copies[other].seq.len()) + PSV_BAND_MARGIN,
            );
            let poasta_m = poa_msa_with_costs(
                &[ref_seq.clone(), copies[other].seq.clone()],
                GapAffine::new(1, 1, 32),
            )
            .ok();
            let fmt = |m: &Vec<Vec<u8>>| {
                let (c, mm, r, g, ld, tr) = msa_cost(m);
                format!(
                    "cost={c} mism={mm} gap_runs={r} gap_bases={g} lead={ld} trail={tr} valid={} cols={}",
                    validate(m),
                    m.first().map(|x| x.len()).unwrap_or(0)
                )
            };
            eprintln!("[psv-cost] copy{other}:");
            if let Some(m) = &banded_m {
                eprintln!("    banded: {}", fmt(m));
            }
            if let Some(m) = &poasta_m {
                eprintln!("    poasta: {}", fmt(m));
            }
        }
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
    let n = copies.len();
    if n < 2 {
        return Vec::new();
    }
    // forward genomic spans (exons + introns), uppercased
    let spans: Vec<Vec<u8>> = copies
        .iter()
        .map(|c| {
            genome
                .fetch_sequence(&c.chrom, c.start, c.end)
                .map(|s| s.to_ascii_uppercase())
                .unwrap_or_default()
        })
        .collect();
    if spans[0].is_empty() {
        return Vec::new();
    }
    let c0 = copies[0];
    let ref_span = &spans[0];
    let is_intronic = |gpos: u64| c0.introns.iter().any(|&(d, a)| gpos >= d && gpos < a);
    let comp = |b: u8| match b {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        x => x,
    };
    use poasta::aligner::scoring::GapAffine;
    const POA_CAP: usize = 20_000; // poasta is exact but O(n^2); above this, fall back to minimap2 (which needs
                                   // >~200 bp to seed, so poasta covers the short cases minimap2 cannot).
    let mut amaps: Vec<BTreeMap<usize, usize>> = vec![BTreeMap::new(); n];
    let mut diff_off: BTreeSet<usize> = BTreeSet::new();
    for other in 1..n {
        if spans[other].is_empty() {
            continue;
        }
        let aln = if ref_span.len().max(spans[other].len()) <= POA_CAP {
            poa_msa_with_costs(&[ref_span.clone(), spans[other].clone()], GapAffine::new(1, 1, 32))
        } else {
            minimap2_msa_pair(ref_span, &spans[other])
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
                        if is_acgt(ca) && is_acgt(cb) && ca != cb && is_intronic(c0.start + ro as u64) {
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
            let g0 = c0.start + ro as u64;
            let b0 = if c0.strand == '-' { comp(ref_span[ro]) } else { ref_span[ro] };
            rec[0] = Some((g0, b0));
            for other in 1..n {
                if let Some(&oo) = amaps[other].get(&ro) {
                    let gc = copies[other].start + oo as u64;
                    let bc = if copies[other].strand == '-' {
                        comp(spans[other][oo])
                    } else {
                        spans[other][oo]
                    };
                    rec[other] = Some((gc, bc));
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
    let mut cols = discover_psvs(copies, &exon_maps);
    // INTRON-RETENTION lever (opt-in): append intronic PSV columns so reads that retain an intron can use the
    // intronic sequence. OFF (env unset / no genome) => `cols` is unchanged => byte-identical to the exon-only gate.
    if let Some(g) = genome {
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
pub fn assign_family_detailed(
    copies: &[&DenovoTranscript],
    reads: &[AlignedRead],
    p: &AssignParams,
    genome: Option<&crate::genome::GenomeIndex>,
) -> FamilyDetail {
    assign_family_detailed_once(copies, reads, p, genome)
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
) -> FamilyDetail {
    // Phase 1: iteratively decide which original copies survive.  `merge_target[i]` always points from
    // original index i to the original index that currently represents it.  Initially each copy represents
    // itself.
    let mut merge_target: Vec<usize> = (0..copies.len()).collect();
    let mut current_indices: Vec<usize> = (0..copies.len()).collect();
    let mut current_copies: Vec<&DenovoTranscript> = copies.to_vec();
    loop {
        let detail = assign_family_detailed_once(&current_copies, reads, p, genome);
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
        let mut detail = assign_family_detailed_once(copies, reads, p, genome);
        detail.copy_indices = current_indices;
        return detail;
    }

    // Phase 2: assign against the FULL original copy set so reads that only overlapped a removed copy
    // are not dropped.  Then remap every per-copy index to its surviving output position.
    let full_detail = assign_family_detailed_once(copies, reads, p, genome);

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

fn assign_family_detailed_once(
    copies: &[&DenovoTranscript],
    reads: &[AlignedRead],
    p: &AssignParams,
    genome: Option<&crate::genome::GenomeIndex>,
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
    let fp = build_family_profiles(copies, genome);
    if timing {
        eprintln!(
            "[timing]     build_family_profiles/discover_psvs ({} copies, {} cols): {:.1}s",
            copies.len(),
            fp.n_cols,
            t_psv.elapsed().as_secs_f64()
        );
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
    // significance certificate downweights them. Built once from all reads' PSV observations.
    let t_edit = std::time::Instant::now(); // its own timer — reusing t_psv double-counted this as ~4.5s
    let editing_cols: Vec<bool> = if p.rna_editing_filter {
        let mut all_obs: Vec<Vec<Option<u8>>> = Vec::with_capacity(reads.len());
        for read in reads {
            if let Some(mc) = best_overlap_copy(read, copies) {
                let mut psv_obs = vec![None; fp.n_cols];
                let mut psv_qual = vec![None; fp.n_cols];
                fill_psv_obs(read, &fp.copy_gpos[mc], fp.strand[mc] == '-', &mut psv_obs, &mut psv_qual);
                all_obs.push(psv_obs);
            }
        }
        detect_editing_columns(&all_obs, &fp.profiles)
    } else {
        Vec::new()
    };
    if timing {
        eprintln!("[timing]     editing-filter pre-pass: {:.1}s", t_edit.elapsed().as_secs_f64());
    }
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
    let per_read: Vec<Option<PerRead>> = reads
        .par_iter()
        .enumerate()
        .map(|(ri, read)| {
            let mc = best_overlap_copy(read, copies)?;
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
            let Some(combined) = assign_read_editing(&feats, &fp.profiles, p, &editing_cols) else {
                return Some(PerRead { mcall, obs_for_em: None, result: None });
            };
            let obs = feats.psv_obs.clone();
            let junctions = feats.junctions.clone();
            let psv_feats = ReadFeatures { psv_obs: feats.psv_obs, psv_qual: feats.psv_qual, junctions: vec![] };
            let Some(psv) = assign_read_editing(&psv_feats, &fp.profiles, p, &editing_cols) else {
                return Some(PerRead { mcall, obs_for_em: Some(obs), result: None });
            };
            Some(PerRead {
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
            })
        })
        .collect();
    if timing {
        eprintln!(
            "[timing]     parallel per-read assign ({} reads): {:.1}s",
            reads.len(),
            t_assign.elapsed().as_secs_f64()
        );
    }
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
        let detail = assign_family_detailed(&copies, &[read], &AssignParams::default(), None);
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
        let detail = assign_family_detailed(&copies, &reads, &AssignParams::default(), None);
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
                    }
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
        let detail = assign_family_detailed(&[&ca, &cb], &[read], &AssignParams::default(), None);
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
        let detail = assign_family_detailed(&[&ca, &cb, &cc], &[], &AssignParams::default(), None);
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
        };
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
        let detail = assign_family_detailed_pruned(&copies, &reads, &p, None);
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
        let detail = assign_family_detailed_pruned(&copies, &reads, &p, None);
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
        let detail = assign_family_detailed(&copies, &reads, &AssignParams::default(), None);
        assert_eq!(detail.copy_psv_alleles.len(), 3, "default (prune off) keeps all 3 copies");
    }
}

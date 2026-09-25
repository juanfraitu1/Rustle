//! Contiguous-core homology kernel shared by family detection, edge confirmation and rescue.
//!
//! `poa_msa_with_costs` (poasta POA MSA with explicit affine costs), the contiguous-core coverage
//! criterion built on it (`contiguous_core_coverage*`, with the bounded/LCS fallback and the
//! process-wide memo), `longest_common_substring`, and `core_coverage_reaches`.
//!
//! The per-exon `FamilyGraph` object this module was named after (exon-equivalence classes, junction
//! edges, minimizer-Jaccard merge, `RUSTLE_VG_FAMILY_MERGE_*` / `RUSTLE_VG_FAMILY_MIN_CORE_COVERAGE`)
//! had no caller in any binary and was removed 2026-09-24; recover it from tag `notebook-2026-09-24`.
//!
//! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use anyhow::Result;

/// Build a POA graph by progressively aligning each sequence under the given poasta `AlignmentConfig`.
/// Generic over the config so `poa_msa_with_costs` can pick exact Dijkstra or exact-A* (`AffineMinGapCost`)
/// without duplicating the loop — both compute the same OPTIMAL-cost alignment (A* only prunes the search).
fn build_poa_graph<C>(seqs: &[Vec<u8>], config: C) -> Result<poasta::graphs::poa::POAGraph<u32>>
where
    C: poasta::aligner::config::AlignmentConfig,
{
    use anyhow::anyhow;
    use poasta::aligner::scoring::AlignmentType;
    use poasta::aligner::PoastaAligner;
    use poasta::graphs::poa::{POAGraph, POANodeIndex};

    let aligner = PoastaAligner::new(config, AlignmentType::Global);
    let mut graph: POAGraph<u32> = POAGraph::new();
    let unit_weights: Vec<usize> = vec![1; seqs.iter().map(|s| s.len()).max().unwrap_or(0)];
    for (i, seq) in seqs.iter().enumerate() {
        let w = &unit_weights[..seq.len()];
        if graph.is_empty() {
            // First sequence: add unaligned (no existing graph to align to).
            graph
                .add_alignment_with_weights(&format!("seq{i}"), seq, None, w)
                .map_err(|e| anyhow!("poasta add_alignment_with_weights failed: {e}"))?;
        } else {
            // Subsequent sequences: align against the current graph.
            let aln_result = aligner.align::<u32, _>(&graph, seq);
            let alignment = if aln_result.alignment.is_empty() {
                None
            } else {
                Some(&aln_result.alignment as &poasta::aligner::Alignment<POANodeIndex<u32>>)
            };
            graph
                .add_alignment_with_weights(&format!("seq{i}"), seq, alignment, w)
                .map_err(|e| anyhow!("poasta add_alignment_with_weights failed: {e}"))?;
        }
    }
    Ok(graph)
}

/// Multiple sequence alignment via poasta POA graph traversal with EXPLICIT
/// affine-gap costs: one row per input sequence, in input order, padded with
/// b'-' for gaps; errors if fewer than 2 sequences are given. The caller picks
/// the scoring, so a specialised alignment (e.g. the contiguous-core gate) can
/// anchor a conserved core against divergent flanks.
///
/// `GapAffine::new(cost_mismatch, cost_gap_extend, cost_gap_open)`.
pub fn poa_msa_with_costs(
    seqs: &[Vec<u8>],
    gap_costs: poasta::aligner::scoring::GapAffine,
) -> Result<Vec<Vec<u8>>> {
    use anyhow::anyhow;

    if seqs.len() < 2 {
        return Err(anyhow!("poa_msa requires at least 2 sequences"));
    }

    // Exact gap-affine global alignment. Default HERE = AffineDijkstra (uniform-cost search).
    // RUSTLE_POA_ASTAR=1 swaps in AffineMinGapCost = the SAME optimal-cost search with an admissible
    // min-gap-cost A* heuristic.
    // VERIFIED (2026-07-12) on real families GSTM(745 cols)/DAZ/RBMY(1542)/PCDHB(3293): A* is BYTE-IDENTICAL to
    // Dijkstra (families/quant/assignments diff = 0) — so the co-optimal-traceback concern does not bite here —
    // but gives NO measurable speedup UNDER THE PSV GAP COSTS (0.0/0.5/12.7s, PCDHB marginally slower):
    // poasta's min-gap heuristic does not prune those paralog alignments. So Dijkstra stays the default on
    // THIS path.
    // ⚠ THE SCOPE OF THAT "no speedup" FINDING WAS TOO WIDE. Re-measured 2026-08-09 on the contiguous-core
    // path, whose gap_open is 32 rather than the PSV costs: A* is 32-44% FASTER there (paired interleaved
    // n=5, 15/15 wins; TFRC 0.563x, TBP 0.593x, HERC2 0.678x) with all 75 control-panel output files
    // byte-identical. That path therefore selects A* explicitly via `contiguous_core_coverage`, which calls
    // `poa_msa_with_costs_cfg` — it does NOT go through this env switch.
    // (minimap2 via RUSTLE_PSV_MINIMAP2 is fast but LOSES PSVs — 109 vs 3293 on PCDHB, 0 vs 745 on GSTM — as it
    // clips divergent flanks; not a safe default. The exact-DP cost is inherent to accurate PSV discovery.)
    // ⚠ WAS `var_os(..).is_some()`, i.e. `RUSTLE_POA_ASTAR=0` TURNED A* ON here. That is the tree's
    // opt-in idiom (`set and not "0"`) inverted, and it silently invalidates any A/B that disables the
    // toggle: measured 2026-08-09, a copy_assign GSTM run with `RUSTLE_POA_ASTAR=0` was ~8% FASTER than
    // the unset default, because the "off" arm had switched THIS path to A*. Parsed as an ordinary
    // opt-in flag now. The DEFAULT (unset) is unchanged, so no shipped configuration moves.
    let astar = matches!(std::env::var("RUSTLE_POA_ASTAR"), Ok(ref v) if v != "0" && !v.is_empty());
    poa_msa_with_costs_cfg_inner(seqs, gap_costs, astar)
}

/// [`poa_msa_with_costs`] with the aligner variant passed EXPLICITLY instead of read from the environment.
///
/// Exists because the two poasta consumers want different defaults and must not share a global switch:
/// the PSV/MSA path measured A* as no faster under ITS gap costs (see the note above), while the
/// contiguous-core path (gap_open=32) measured 32-44% faster with byte-identical output. A caller that
/// picks the variant can be memoized, because the variant is then an explicit part of the cache key
/// rather than ambient state.
pub fn poa_msa_with_costs_cfg(
    seqs: &[Vec<u8>],
    gap_costs: poasta::aligner::scoring::GapAffine,
    astar: bool,
) -> Result<Vec<Vec<u8>>> {
    use anyhow::anyhow;
    if seqs.len() < 2 {
        return Err(anyhow!("poa_msa requires at least 2 sequences"));
    }
    poa_msa_with_costs_cfg_inner(seqs, gap_costs, astar)
}

fn poa_msa_with_costs_cfg_inner(
    seqs: &[Vec<u8>],
    gap_costs: poasta::aligner::scoring::GapAffine,
    astar: bool,
) -> Result<Vec<Vec<u8>>> {
    use anyhow::anyhow;
    use poasta::aligner::config::{AffineDijkstra, AffineMinGapCost};
    use poasta::graphs::poa::POAGraph;
    use poasta::io::fasta::poa_graph_to_fasta;

    let graph: POAGraph<u32> = if astar {
        build_poa_graph(seqs, AffineMinGapCost(gap_costs))?
    } else {
        build_poa_graph(seqs, AffineDijkstra(gap_costs))?
    };

    // Extract the MSA rows by using poasta's public poa_graph_to_fasta function,
    // which walks the graph in topological order and writes aligned FASTA records.
    // We write into a Vec<u8> buffer and parse the resulting FASTA lines.
    let mut fasta_buf: Vec<u8> = Vec::new();
    poa_graph_to_fasta(&graph, &mut fasta_buf)
        .map_err(|e| anyhow!("poa_graph_to_fasta failed: {e}"))?;

    // Parse the FASTA output: each sequence record is a set of lines starting with '>'.
    // We collect the sequences in order and convert them to Vec<u8> rows.
    let mut msa: Vec<Vec<u8>> = Vec::with_capacity(seqs.len());
    let mut current_seq: Option<Vec<u8>> = None;
    for line in fasta_buf.split(|&b| b == b'\n') {
        if line.is_empty() {
            continue;
        }
        if line[0] == b'>' {
            if let Some(seq) = current_seq.take() {
                msa.push(seq);
            }
            current_seq = Some(Vec::new());
        } else if let Some(ref mut seq) = current_seq {
            seq.extend_from_slice(line);
        }
    }
    if let Some(seq) = current_seq {
        msa.push(seq);
    }

    if msa.len() != seqs.len() {
        return Err(anyhow!(
            "poa_msa: expected {} rows from graph but got {}",
            seqs.len(),
            msa.len()
        ));
    }

    // poasta's fasta_aln_for_seq has an off-by-one in the trailing-gap fill for
    // sequences whose path ends before max_col, so rows can occasionally come
    // back one column short. Pad with trailing gaps to the max row length —
    // this is content-neutral for an MSA (gaps are not part of the ungapped
    // sequence) and the round-trip ungap(row) ≡ original input is preserved.
    let n_cols = msa.iter().map(|r| r.len()).max().unwrap_or(0);
    for row in msa.iter_mut() {
        if row.len() < n_cols {
            row.extend(std::iter::repeat(b'-').take(n_cols - row.len()));
        }
    }
    Ok(msa)
}

/// POA-derived CONTIGUOUS-CORE coverage of two sequences.
///
/// Aligns `a` and `b` via the 2-sequence POA instance (`poa_msa`, which for a
/// pair reduces to a single global alignment), then finds the LONGEST run of
/// consecutive alignment columns where BOTH rows are non-gap AND carry the same
/// base, and returns that run length divided by `min(a.len(), b.len())`.
///
/// This is the validated criterion from `bench/poa_family_definition.py`: a true
/// paralog copy shares ONE long contiguous homologous core (HIGH coverage, the
/// prototype's `>= 0.13` band), whereas a domain-sharer co-aligns over only a
/// short block and then diverges (LOW coverage, indistinguishable from random
/// cross-family controls). Unlike all-column reciprocal coverage it is NOT
/// inflated by a global aligner's scattered chance-match filler, because it
/// requires the matching columns to be CONTIGUOUS.
///
/// Purely alignment-column-derived (POA-only): no DNA/protein domain annotation,
/// no BLAST, no k-mer/minimizer is used. Deterministic (POA is deterministic).
///
/// ROBUSTNESS: this uses a DEDICATED alignment config — a STRONG gap-open
/// (`GapAffine::new(1, 1, 32)`: mismatch=1, gap_extend=1, gap_open=32) — instead
/// of the main `poa_msa` scoring (gap_open=2). With the weak default gap-open,
/// poasta is content-dependently unstable on two copies that share a long core
/// but have DIVERGENT 5' AND 3' flanks: a few CHEAP single-base gaps in the
/// divergent flanks frame-shift the otherwise-on-diagonal alignment THROUGH the
/// conserved core, so the core never re-anchors and the longest equal run
/// collapses to a chance-match ~3 bp (observed core coverage 0.71 -> ~0.01 on
/// true copies). A strong gap-open makes those scattered frame-shift gaps
/// uneconomical, so the alignment stays on one diagonal and the identical core
/// aligns column-for-column (coverage restored to ~ the core fraction).
///
/// This mirrors the validated python prototype (`bench/poa_family_definition.py`,
/// BioPython NW gap_open=-5) which avoided the same instability by gapping the
/// divergent flank and anchoring the core. The strong gap-open ONLY fixes the
/// false-collapse of true copies; it leaves domain-sharers / disjoint pairs LOW
/// (measured: short-block 0.05, reordered-chunk 0.01, disjoint 0.01 — unchanged
/// from the weak-gap default), so the separation is preserved, not blurred. The
/// gap-open is the documented robustness lever: gap_open>=16 anchors all tested
/// divergent-flank cases; 32 is a comfortable margin (poasta's internal Score
/// arithmetic overflows for gap_open<2 and RAISING MISMATCH was counterproductive
/// — both empirically falsified, so the strong gap-open is the chosen knob).
/// poasta 0.1.0's `EndsFree`/semi-global mode is `todo!()` (panics), so an
/// ends-free alignment was not available; the strong-gap-open config is the
/// POA-only fix.
///
/// Edge cases: returns 0.0 if either sequence is empty or the alignment fails;
/// identical sequences return ~1.0 (the whole sequence is one matched run).
pub fn contiguous_core_coverage(a: &[u8], b: &[u8]) -> f64 {
    contiguous_core_coverage_with(a, b, core_astar())
}

/// [`contiguous_core_coverage`] with the aligner variant chosen by the caller rather than by
/// [`core_astar`]. Both variants are exact optimal-cost searches over the same graph and scoring; which
/// one is faster depends on the SHAPE of the pair, which is a property of the CALL SITE (see
/// [`contiguous_core_coverage_bounded_with`]).
pub fn contiguous_core_coverage_with(a: &[u8], b: &[u8], astar: bool) -> f64 {
    use poasta::aligner::scoring::GapAffine;
    let minlen = a.len().min(b.len());
    if minlen == 0 {
        return 0.0;
    }
    // Dedicated strong-gap-open scoring (mismatch=1, gap_extend=1, gap_open=32).
    // See the doc comment: anchors the conserved core against divergent flanks.
    let core_gap_costs = GapAffine::new(1, 1, 32);
    let msa = match poa_msa_with_costs_cfg(&[a.to_vec(), b.to_vec()], core_gap_costs, astar) {
        Ok(m) => m,
        Err(_) => return 0.0,
    };
    if msa.len() != 2 {
        return 0.0;
    }
    let (row_a, row_b) = (&msa[0], &msa[1]);
    let n_cols = row_a.len().min(row_b.len());
    let mut longest = 0usize;
    let mut run = 0usize;
    for c in 0..n_cols {
        let (ca, cb) = (row_a[c], row_b[c]);
        if ca != b'-' && cb != b'-' && ca == cb {
            run += 1;
            if run > longest {
                longest = run;
            }
        } else {
            run = 0;
        }
    }
    longest as f64 / minlen as f64
}

/// Longest common SUBSTRING (contiguous, EXACT) of `a` and `b`, in O(min(|a|,|b|)) memory via a suffix
/// automaton built on the SHORTER sequence and scanned by the longer.
///
/// This is the memory-bounded equivalent of [`contiguous_core_coverage`]'s "longest ungapped equal run": a
/// run there is a maximal stretch of alignment columns that are BOTH non-gap AND equal — i.e. a substring
/// shared by both sequences (a mismatch resets the run, so there are no internal mismatches). poasta finds it
/// via a memory-hungry graph alignment that OOMs on a long (e.g. 228 kb read-through) operand; the suffix
/// automaton finds the GLOBAL longest common substring directly in LINEAR memory, so it never OOMs.
pub fn longest_common_substring(a: &[u8], b: &[u8]) -> usize {
    use std::collections::BTreeMap;
    let (s, t) = if a.len() <= b.len() { (a, b) } else { (b, a) };
    if s.is_empty() || t.is_empty() {
        return 0;
    }
    // suffix automaton over the shorter string `s`. State = (longest len in its endpos class, suffix link,
    // outgoing transitions). Deterministic (BTreeMap transitions; the result never depends on map order).
    struct St {
        len: i32,
        link: i32,
        next: BTreeMap<u8, i32>,
    }
    let mut st: Vec<St> = Vec::with_capacity(2 * s.len());
    st.push(St { len: 0, link: -1, next: BTreeMap::new() });
    let mut last = 0i32;
    for &c in s {
        let cur = st.len() as i32;
        let cur_len = st[last as usize].len + 1;
        st.push(St { len: cur_len, link: -1, next: BTreeMap::new() });
        let mut p = last;
        while p != -1 && !st[p as usize].next.contains_key(&c) {
            st[p as usize].next.insert(c, cur);
            p = st[p as usize].link;
        }
        if p == -1 {
            st[cur as usize].link = 0;
        } else {
            let q = st[p as usize].next[&c];
            if st[p as usize].len + 1 == st[q as usize].len {
                st[cur as usize].link = q;
            } else {
                let clone = st.len() as i32;
                let clone_len = st[p as usize].len + 1;
                let (qlink, qnext) = (st[q as usize].link, st[q as usize].next.clone());
                st.push(St { len: clone_len, link: qlink, next: qnext });
                while p != -1 && st[p as usize].next.get(&c) == Some(&q) {
                    st[p as usize].next.insert(c, clone);
                    p = st[p as usize].link;
                }
                st[q as usize].link = clone;
                st[cur as usize].link = clone;
            }
        }
        last = cur;
    }
    // scan `t`, tracking the longest match ending at each position (canonical SAM-LCS walk).
    let (mut v, mut l, mut best) = (0i32, 0i32, 0i32);
    for &c in t {
        while v != 0 && !st[v as usize].next.contains_key(&c) {
            v = st[v as usize].link;
            l = st[v as usize].len;
        }
        match st[v as usize].next.get(&c) {
            Some(&nx) => {
                v = nx;
                l += 1;
            }
            None => l = 0, // dead-ended at the root
        }
        if l > best {
            best = l;
        }
    }
    best as usize
}

/// Which exact aligner the CONTIGUOUS-CORE path uses. A* (`AffineMinGapCost`) is the default here — see the
/// re-measurement note on [`poa_msa_with_costs`]. `RUSTLE_POA_ASTAR=0` forces the old Dijkstra search back on
/// this path so the change can be A/B'd without a rebuild; any other value (set or unset) means A*.
///
/// Both are EXACT optimal-cost searches over the same graph and scoring, so this selects only how the
/// optimum is found, not what it is. It is nonetheless part of the memo key below, because co-optimal
/// traceback is not proven unique and a cache must never serve a value the current setting did not produce.
pub fn core_astar() -> bool {
    !matches!(std::env::var("RUSTLE_POA_ASTAR"), Ok(v) if v == "0")
}

/// ASCII-uppercase VIEW of `s` that allocates only when `s` actually contains a lowercase ASCII letter.
///
/// The POA collapse loop uppercased both operands of every candidate pair on every sweep, i.e. two full
/// sequence copies per pair per iteration — while `GenomeIndex` already uppercases every base at load, so
/// the copies were almost always identical to their input. The uppercase is still APPLIED (soft-masked
/// input reaching these paths from elsewhere must behave exactly as before); only the allocation is
/// conditional, so the bytes handed downstream are unchanged.
pub fn upper_cow(s: &[u8]) -> std::borrow::Cow<'_, [u8]> {
    if s.iter().any(|c| c.is_ascii_lowercase()) {
        std::borrow::Cow::Owned(s.to_ascii_uppercase())
    } else {
        std::borrow::Cow::Borrowed(s)
    }
}

// ---------------------------------------------------------------------------------------------------
// Memo for the contiguous-core kernel.
//
// WHY: `collapse_parent` re-runs its full O(r^2) candidate sweep until an entire sweep merges nothing, so
// the LAST sweep — by definition the one that changes no state — re-aligns every surviving pair and throws
// every result away, and each earlier sweep re-aligns every pair it failed to merge. The kernel is a pure
// function of its arguments, so those repeats are recomputation with no other effect.
//
// THE KEY IS EVERY ARGUMENT THAT CAN CHANGE THE VALUE, and nothing else:
//   (interned id of `a`'s exact bytes, interned id of `b`'s exact bytes, `poasta_cap`, `core_astar()`)
// - The two sequences are keyed SEPARATELY AND IN ORDER. `poa_msa_with_costs` builds the graph from the
//   FIRST sequence, so f(a,b) == f(b,a) is UNPROVEN; normalising the pair order into one key would be a
//   result-changing edit disguised as a cache. (`memo_order_is_part_of_the_key` pins this.)
// - `poasta_cap` selects the poasta-vs-suffix-automaton BRANCH, so it changes the value, not just the cost.
// - The aligner variant is ambient (env), so it is resolved to a bool and stored explicitly.
// - The gap costs are NOT in the key because they are a hard-coded constant of `contiguous_core_coverage`
//   (1,1,32). If they ever become a parameter they MUST be added here.
// - The merge THRESHOLDS (`collapse_span_core`, COLLAPSE_CONTAIN_FRAC, `t_core`) are deliberately ABSENT:
//   they are applied by the callers to this function's return value and cannot change it. That separation
//   is what makes a threshold sweep reuse the alignments.
//
// Sequences are interned so the table stores each distinct sequence once (O(r) bytes) rather than once per
// pair (O(r^2) bytes), and the key is the exact bytes — not a hash — so there is no collision path to a
// silently wrong value.
// ---------------------------------------------------------------------------------------------------

/// Interned-byte budget. Past this the memo stops INSERTING (it still serves hits), so a genome-wide run
/// degrades to today's recompute-everything behaviour instead of growing without bound.
const CORE_MEMO_MAX_BYTES: usize = 512 * 1024 * 1024;
/// Entry-count budget, same degradation rule. 4M entries of a 4-word key + f64 is ~200 MB.
const CORE_MEMO_MAX_ENTRIES: usize = 4_000_000;

/// FxHash, not SipHash. The intern map is probed with a FULL SEQUENCE (up to `len_cap` = 20 kb) on EVERY
/// call including misses, so the hash pass is the memo's entire per-call overhead — free next to a POA
/// alignment, not free next to a memo HIT or the cheap suffix-automaton branch.
///
/// MEASURED, not assumed (medians of 3, warm cache, interleaved builds): on HERC2 — the panel region
/// where the memo actually fires — SipHash 13.423 s [13.377-13.474] vs FxHash 9.867 s [9.803-10.091],
/// a 26% difference with both spreads under 0.3 s. Neutral on the regions where the memo rarely hits
/// (GSTM 0.976, MAGEA 0.987, RABL2 1.001, SDHA 0.997, TBP 0.996).
///
/// Exactness is unaffected: `HashMap` still compares keys with `Eq`, so a hash collision costs a byte
/// comparison, never a wrong value. (Also deterministic — FxHash has no random seed — though nothing
/// here depends on iteration order.)
#[derive(Default)]
struct CoreMemo {
    ids: crate::types::DetHashMap<Vec<u8>, u32>,
    interned_bytes: usize,
    vals: crate::types::DetHashMap<(u32, u32, usize, bool), f64>,
    hits: u64,
    misses: u64,
}

fn core_memo() -> &'static std::sync::Mutex<CoreMemo> {
    static MEMO: std::sync::OnceLock<std::sync::Mutex<CoreMemo>> = std::sync::OnceLock::new();
    MEMO.get_or_init(|| std::sync::Mutex::new(CoreMemo::default()))
}

/// A poisoned memo is still a valid cache (the data is plain values, and a panic mid-update cannot leave a
/// wrong entry — insertion is the last step), so recover the guard rather than cascading the panic.
fn core_memo_lock() -> std::sync::MutexGuard<'static, CoreMemo> {
    core_memo().lock().unwrap_or_else(|e| e.into_inner())
}

/// `(hits, misses)` on the contiguous-core memo since process start. Test/diagnostic hook — this is how
/// `memo_*` tests prove that a changed key component MISSES rather than silently reusing a stale value.
pub fn core_memo_stats() -> (u64, u64) {
    let m = core_memo_lock();
    (m.hits, m.misses)
}

/// Is `(a, b, poasta_cap, astar)` — the FULL key — currently cached? Diagnostic, and the hook the memo
/// key tests use: asserting on presence per key is immune to sibling tests sharing the process-global
/// table, which a hit/miss counter is not.
pub fn core_memo_contains(a: &[u8], b: &[u8], poasta_cap: usize, astar: bool) -> bool {
    let m = core_memo_lock();
    match (m.ids.get(a), m.ids.get(b)) {
        (Some(&ia), Some(&ib)) => m.vals.contains_key(&(ia, ib, poasta_cap, astar)),
        _ => false,
    }
}

/// Plant a value under an exact key WITHOUT running the kernel. Test-only: a planted value the kernel
/// would never return is what makes "this call was served from the cache under THIS key" a deterministic
/// assertion rather than an inference from a counter.
#[cfg(test)]
fn core_memo_plant(a: &[u8], b: &[u8], poasta_cap: usize, astar: bool, v: f64) {
    let mut m = core_memo_lock();
    let ia = match m.ids.get(a).copied() {
        Some(i) => i,
        None => {
            let i = m.ids.len() as u32;
            m.ids.insert(a.to_vec(), i);
            m.interned_bytes += a.len();
            i
        }
    };
    let ib = match m.ids.get(b).copied() {
        Some(i) => i,
        None => {
            let i = m.ids.len() as u32;
            m.ids.insert(b.to_vec(), i);
            m.interned_bytes += b.len();
            i
        }
    };
    m.vals.insert((ia, ib, poasta_cap, astar), v);
}

/// Drop every cached value and interned sequence. Only for tests that need a known starting point.
pub fn core_memo_clear() {
    let mut m = core_memo_lock();
    m.ids.clear();
    m.vals.clear();
    m.interned_bytes = 0;
    m.hits = 0;
    m.misses = 0;
}

/// [`contiguous_core_coverage`] with a MEMORY GUARD: when the larger sequence exceeds `poasta_cap`, poasta's
/// graph aligner would OOM, so use the linear-memory longest-common-substring metric instead (faithful — a
/// poasta ungapped-equal run IS a common substring). At or below the cap, the exact poasta path is used, so
/// the validated small-transcript behaviour is byte-identical.
///
/// MEMOIZED — see the key documentation above. The memo is transparent: it never changes what this function
/// returns, only how often the kernel underneath it runs.
pub fn contiguous_core_coverage_bounded(a: &[u8], b: &[u8], poasta_cap: usize) -> f64 {
    contiguous_core_coverage_bounded_with(a, b, poasta_cap, core_astar())
}

/// `contiguous_core_coverage_bounded(a, b, cap) >= threshold`, with an EXACT early exit.
///
/// The core is the longest run of equal, non-gap alignment columns, which is a common SUBSTRING of `a` and
/// `b`; so `core <= longest_common_substring(a, b) / minlen`, and both sides use the same division, so the
/// float comparison is monotone. When that bound is already below `threshold` the verdict is `false`
/// without aligning: on divergent 4-6 kb pairs the poasta call costs seconds to tens of seconds (gorilla
/// NC_073244.2 locus collapse: 331 of 331 POA CPU-s were such pairs before the bound, 89 after), while the
/// suffix-automaton bound costs milliseconds. Byte-identical verdicts by construction.
pub fn core_coverage_reaches(a: &[u8], b: &[u8], poasta_cap: usize, threshold: f64) -> bool {
    let minlen = a.len().min(b.len());
    // above the cap the kernel IS the LCS ratio, so the bound would only compute it twice
    if minlen > 0
        && a.len().max(b.len()) <= poasta_cap
        && (longest_common_substring(a, b) as f64 / minlen as f64) < threshold
    {
        return false;
    }
    contiguous_core_coverage_bounded(a, b, poasta_cap) >= threshold
}

/// A* IS NOT UNIFORMLY BETTER — WHICH VARIANT WINS IS A PROPERTY OF THE CALL SITE, SO THE CALL SITE PICKS.
///
/// Measured 2026-08-09, both directions on real data:
/// - LOCUS COLLAPSE (`collapse_parent` / `distinct_locus_reps`), which aligns long OVERLAPPING transcript
///   models up to `len_cap` = 20 kb: A* is a large win. 25-region control panel 93.0 s -> 40.8 s overall;
///   isolating the aligner alone on the POA-bound regions, HERC2 34.9 -> 17.6 s (0.504x).
/// - EDGE CONFIRMATION (`confirm_edge`, rayon-parallel over candidate rep pairs) and the rescue path: A*
///   LOSES. `copy_assign` on GSTM, where `RUSTLE_LOCUS_JUNCTION_ONLY=1` shows the collapse costs nothing,
///   is 30.20 s with A* off vs 31.68 s with it on (medians of 3, interleaved, rotated order) — the
///   min-gap heuristic's per-node overhead is not repaid on those pairs.
///
/// So `EDGE_CONFIRM_ASTAR = false` is not a tuning constant; it is the pre-existing behaviour of that path,
/// left in place because the measurement that justified changing the collapse does not extend to it.
pub fn contiguous_core_coverage_bounded_with(
    a: &[u8],
    b: &[u8],
    poasta_cap: usize,
    astar: bool,
) -> f64 {
    // `RUSTLE_POA_MEMO=0` bypasses the memo entirely. Kept because the ONLY way to be sure a cache is
    // transparent on real data is to be able to re-run the same job without it, and because it is what
    // measures the memo's contribution separately from the other speed changes.
    if matches!(std::env::var("RUSTLE_POA_MEMO"), Ok(ref v) if v == "0") {
        return contiguous_core_coverage_bounded_uncached_with(a, b, poasta_cap, astar);
    }
    // Look up first; the lock is held only across two hash probes, never across an alignment.
    let ids = {
        let mut m = core_memo_lock();
        match (m.ids.get(a).copied(), m.ids.get(b).copied()) {
            (Some(ia), Some(ib)) => {
                if let Some(&v) = m.vals.get(&(ia, ib, poasta_cap, astar)) {
                    m.hits += 1;
                    return v;
                }
                m.misses += 1;
                Some((ia, ib))
            }
            _ => {
                m.misses += 1;
                None
            }
        }
    };

    let v = contiguous_core_coverage_bounded_uncached_with(a, b, poasta_cap, astar);

    let mut m = core_memo_lock();
    if m.vals.len() >= CORE_MEMO_MAX_ENTRIES {
        return v;
    }
    let (ia, ib) = match ids {
        Some(p) => p,
        None => {
            if m.interned_bytes.saturating_add(a.len() + b.len()) > CORE_MEMO_MAX_BYTES {
                return v;
            }
            let ia = match m.ids.get(a).copied() {
                Some(i) => i,
                None => {
                    let i = m.ids.len() as u32;
                    m.ids.insert(a.to_vec(), i);
                    m.interned_bytes += a.len();
                    i
                }
            };
            let ib = match m.ids.get(b).copied() {
                Some(i) => i,
                None => {
                    let i = m.ids.len() as u32;
                    m.ids.insert(b.to_vec(), i);
                    m.interned_bytes += b.len();
                    i
                }
            };
            (ia, ib)
        }
    };
    m.vals.insert((ia, ib, poasta_cap, astar), v);
    v
}

/// The kernel itself, with no memo in front of it. Kept separate so the memo can be proven transparent by
/// comparing the two on the same inputs.
pub fn contiguous_core_coverage_bounded_uncached(a: &[u8], b: &[u8], poasta_cap: usize) -> f64 {
    contiguous_core_coverage_bounded_uncached_with(a, b, poasta_cap, core_astar())
}

/// [`contiguous_core_coverage_bounded_uncached`] with the aligner variant chosen by the caller.
pub fn contiguous_core_coverage_bounded_uncached_with(
    a: &[u8],
    b: &[u8],
    poasta_cap: usize,
    astar: bool,
) -> f64 {
    if a.len().max(b.len()) > poasta_cap {
        let minlen = a.len().min(b.len());
        if minlen == 0 {
            return 0.0;
        }
        longest_common_substring(a, b) as f64 / minlen as f64
    } else {
        contiguous_core_coverage_with(a, b, astar)
    }
}

/// The aligner variant used by EDGE CONFIRMATION and the rescue path (`confirm_edge`,
/// `family_rescue`): plain Dijkstra, i.e. exactly what those paths did before the collapse path moved to
/// A*. See [`contiguous_core_coverage_bounded_with`] for the measurement in both directions.
pub const EDGE_CONFIRM_ASTAR: bool = false;

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Mutex;

    /// SplitMix64 PRNG for deterministic test sequences (no test-time RNG dependency).
    struct SplitMix64(u64);
    impl SplitMix64 {
        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
            z ^ (z >> 31)
        }
    }

    // ----------------------------------------------------------------------
    // POA-derived CONTIGUOUS-CORE coverage criterion (see
    // bench/poa_family_definition.py). contiguous_core_coverage(a,b) aligns the
    // two sequences via the 2-sequence POA instance (poa_msa) and reads off the
    // LONGEST run of consecutive columns where BOTH rows are non-gap AND equal,
    // divided by the shorter sequence. A true copy shares ONE long homologous
    // core (HIGH); a domain-sharer shares only a short block then diverges (LOW).
    // ----------------------------------------------------------------------

    /// Deterministic random DNA of length `n` (SplitMix64, no test-time RNG dep).
    fn core_rand_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut rng = SplitMix64(seed);
        const B: [u8; 4] = [b'A', b'C', b'G', b'T'];
        (0..n).map(|_| B[(rng.next_u64() % 4) as usize]).collect()
    }

    // ----------------------------------------------------------------------
    // MEMO KEY TESTS for `contiguous_core_coverage_bounded`.
    //
    // The rule these enforce: the cache key contains EVERY argument that can change the returned value —
    // both sequences (in order), the poasta cap, and the resolved aligner variant — and a change to any
    // of them MISSES. Each test plants a sentinel value (0.4242, which the kernel cannot produce for
    // these inputs) under one exact key, then shows that only that key is served and every neighbouring
    // key recomputes.
    //
    // Every test uses sequences unique to itself, so the process-global table cannot be perturbed by a
    // sibling test running concurrently; presence is asserted per key rather than via a global counter.
    // ----------------------------------------------------------------------

    /// RUSTLE_POA_ASTAR is process-global. Serialize the tests that mutate it.
    static ASTAR_ENV_LOCK: Mutex<()> = Mutex::new(());
    const MEMO_SENTINEL: f64 = 0.4242;

    #[test]
    fn memo_serves_its_key_and_a_changed_cap_misses() {
        let a = core_rand_seq(300, 0x0BAD_1001);
        let b = core_rand_seq(300, 0x0BAD_1002);
        let astar = core_astar();

        // Plant a value the kernel would never produce, under cap = 10_000 only.
        core_memo_plant(&a, &b, 10_000, astar, MEMO_SENTINEL);
        assert!(core_memo_contains(&a, &b, 10_000, astar));
        assert_eq!(
            contiguous_core_coverage_bounded(&a, &b, 10_000),
            MEMO_SENTINEL,
            "a repeat call under the SAME key must be served from the memo"
        );

        // A DIFFERENT cap is a different key. It must miss -- and it must miss for a real reason: the cap
        // selects the poasta-vs-suffix-automaton branch, so it genuinely changes the value.
        assert!(!core_memo_contains(&a, &b, 100, astar));
        let v_small_cap = contiguous_core_coverage_bounded(&a, &b, 100);
        assert_ne!(
            v_small_cap, MEMO_SENTINEL,
            "changing poasta_cap must MISS the cache, not reuse the value cached for another cap"
        );
        assert_eq!(
            v_small_cap,
            contiguous_core_coverage_bounded_uncached(&a, &b, 100),
            "the value computed after the miss must equal the uncached kernel"
        );
        // ...and the original key is untouched by the miss.
        assert_eq!(contiguous_core_coverage_bounded(&a, &b, 10_000), MEMO_SENTINEL);
    }

    #[test]
    fn memo_order_is_part_of_the_key() {
        // poasta builds its graph from the FIRST sequence, so f(a,b) == f(b,a) is UNPROVEN. The key must
        // therefore be ordered: swapping the operands must recompute, never reuse.
        let a = core_rand_seq(280, 0x0BAD_2001);
        let b = core_rand_seq(280, 0x0BAD_2002);
        let astar = core_astar();
        core_memo_plant(&a, &b, 10_000, astar, MEMO_SENTINEL);

        assert!(!core_memo_contains(&b, &a, 10_000, astar));
        let swapped = contiguous_core_coverage_bounded(&b, &a, 10_000);
        assert_ne!(
            swapped, MEMO_SENTINEL,
            "(b,a) must MISS the entry cached for (a,b) -- pair order is part of the key"
        );
        assert_eq!(swapped, contiguous_core_coverage_bounded_uncached(&b, &a, 10_000));
    }

    /// The aligner variant is part of the key: the two variants have SEPARATE entries, and the lookup
    /// uses whichever variant `core_astar()` currently resolves to.
    ///
    /// Deliberately does NOT mutate RUSTLE_POA_ASTAR while an alignment is in flight — that would change
    /// the aligner under any concurrently-running test. The env plumbing is pinned separately, by
    /// `core_astar_defaults_to_true_and_only_zero_disables`, which touches no alignment.
    #[test]
    fn memo_astar_variant_is_part_of_the_key() {
        let a = core_rand_seq(260, 0x0BAD_3001);
        let b = core_rand_seq(260, 0x0BAD_3002);
        let current = core_astar();

        // Plant DIFFERENT sentinels under the two variants of the same (a, b, cap).
        core_memo_plant(&a, &b, 10_000, current, MEMO_SENTINEL);
        core_memo_plant(&a, &b, 10_000, !current, MEMO_SENTINEL + 0.1);
        assert!(core_memo_contains(&a, &b, 10_000, true));
        assert!(core_memo_contains(&a, &b, 10_000, false));

        assert_eq!(
            contiguous_core_coverage_bounded(&a, &b, 10_000),
            MEMO_SENTINEL,
            "the lookup must use the CURRENTLY resolved aligner variant, not the other entry"
        );

        // And with only the OTHER variant cached, the current one must miss and recompute.
        let c = core_rand_seq(260, 0x0BAD_3003);
        core_memo_plant(&a, &c, 10_000, !current, MEMO_SENTINEL);
        assert!(!core_memo_contains(&a, &c, 10_000, current));
        let v = contiguous_core_coverage_bounded(&a, &c, 10_000);
        assert_ne!(
            v, MEMO_SENTINEL,
            "a value cached for the other aligner variant must NOT be served"
        );
        assert_eq!(v, contiguous_core_coverage_bounded_uncached(&a, &c, 10_000));
    }

    /// `core_astar()` is the env plumbing only — no alignment runs inside the locked window, so mutating
    /// the process-global variable here cannot perturb a concurrent test's aligner.
    #[test]
    fn core_astar_defaults_to_true_and_only_zero_disables() {
        let _guard = ASTAR_ENV_LOCK.lock().unwrap_or_else(|p| p.into_inner());
        let prior = std::env::var("RUSTLE_POA_ASTAR").ok();
        std::env::remove_var("RUSTLE_POA_ASTAR");
        let unset = core_astar();
        std::env::set_var("RUSTLE_POA_ASTAR", "0");
        let zero = core_astar();
        std::env::set_var("RUSTLE_POA_ASTAR", "1");
        let one = core_astar();
        match prior {
            Some(v) => std::env::set_var("RUSTLE_POA_ASTAR", v),
            None => std::env::remove_var("RUSTLE_POA_ASTAR"),
        }
        assert!(unset, "A* is the DEFAULT on the contiguous-core path (measured 32-44% faster there)");
        assert!(!zero, "RUSTLE_POA_ASTAR=0 must restore the Dijkstra search on this path");
        assert!(one);
    }

    #[test]
    fn memo_is_transparent_versus_the_uncached_kernel() {
        // The memo must never change what the function returns -- only how often the kernel runs. Cover
        // both branches of the cap guard, an identical pair, an empty operand, and a repeat call.
        let core = core_rand_seq(200, 0x0BAD_4000);
        let mut a = core_rand_seq(60, 0x0BAD_4001);
        a.extend_from_slice(&core);
        let mut b = core_rand_seq(60, 0x0BAD_4002);
        b.extend_from_slice(&core);
        let empty: Vec<u8> = Vec::new();
        let cases: [(&[u8], &[u8], usize); 6] = [
            (&a, &b, 10_000), // poasta branch
            (&a, &b, 10),     // suffix-automaton branch (cap exceeded)
            (&a, &a, 10_000), // identical operands
            (&a, &b, 10_000), // repeat -> served from cache, must still agree
            (&empty, &b, 10_000),
            (&a, &empty, 10),
        ];
        for (x, y, cap) in cases {
            assert_eq!(
                contiguous_core_coverage_bounded(x, y, cap),
                contiguous_core_coverage_bounded_uncached(x, y, cap),
                "memoized result must equal the uncached kernel (cap={cap})"
            );
        }
    }

    #[test]
    fn upper_cow_borrows_when_already_uppercase_and_still_uppercases() {
        use std::borrow::Cow;
        let up = b"ACGTNACGT".to_vec();
        assert!(matches!(upper_cow(&up), Cow::Borrowed(_)), "no lowercase -> no allocation");
        let mixed = b"acgtNACgt".to_vec();
        let got = upper_cow(&mixed);
        assert!(matches!(got, Cow::Owned(_)), "lowercase present -> must allocate");
        assert_eq!(
            got.as_ref(),
            mixed.to_ascii_uppercase().as_slice(),
            "the bytes handed downstream must be identical to to_ascii_uppercase()"
        );
    }

    #[test]
    fn contiguous_core_coverage_long_shared_core_is_high() {
        // Two "true copy" sequences: a long identical 400 bp middle core,
        // divergent (independent random) flanks on each side.
        let core = core_rand_seq(400, 0xC0FE_0001);
        let mut a = core_rand_seq(80, 0xAAAA_0001);
        a.extend_from_slice(&core);
        a.extend(core_rand_seq(80, 0xAAAA_0002));
        let mut b = core_rand_seq(80, 0xBBBB_0001);
        b.extend_from_slice(&core);
        b.extend(core_rand_seq(80, 0xBBBB_0002));
        let cov = contiguous_core_coverage(&a, &b);
        assert!(cov >= 0.13,
            "long shared core should give HIGH contiguous-core coverage (got {cov:.3})");
    }

    /// ROBUSTNESS REGRESSION GUARD (the known defect). Two TRUE copies that share
    /// a long internal core (400 bp) but have LONG DIVERGENT flanks on BOTH the 5'
    /// AND 3' ends (120 bp each, independent random per copy) must STILL score a
    /// HIGH contiguous-core coverage (~ the core fraction 400/640 ≈ 0.62).
    ///
    /// Before the fix, poasta's default affine gap costs threaded the divergent
    /// flanks diagonally rather than gapping them, so the aligner never re-anchored
    /// the conserved core and coverage collapsed (observed 0.71 -> 0.01). A
    /// content-dependent collapse like that would wrongly SPLIT real paralog copies
    /// at the family-merge gate. We require coverage >= 0.5 here (well above the
    /// 0.13 merge bar, and close to the 0.62 core fraction).
    #[test]
    fn contiguous_core_coverage_divergent_flanks_still_high() {
        // Construction chosen (seed=4) to RELIABLY trigger the threading collapse
        // under poasta's default affine gap costs: with the old `poa_msa` scoring
        // this exact input scored ~0.005. The fix must restore it to ~0.62.
        const FLANK: usize = 120;
        const S: u64 = 4;
        let core = core_rand_seq(400, 0xC0FE_0000 ^ (S * 0x1001));
        // 120 bp DIVERGENT (independent random) flanks on BOTH ends of each copy.
        let mut a = core_rand_seq(FLANK, 0xAAAA_0000 ^ (S * 0x2002));
        a.extend_from_slice(&core);
        a.extend(core_rand_seq(FLANK, 0xAAAA_1000 ^ (S * 0x3003)));
        let mut b = core_rand_seq(FLANK, 0xBBBB_0000 ^ (S * 0x4004));
        b.extend_from_slice(&core);
        b.extend(core_rand_seq(FLANK, 0xBBBB_1000 ^ (S * 0x5005)));
        let cov = contiguous_core_coverage(&a, &b);
        assert!(cov >= 0.5,
            "two true copies sharing a 400 bp core with divergent 5' AND 3' flanks \
             must score HIGH contiguous-core coverage (~0.62); got {cov:.3} — the \
             aligner threaded the flanks and lost the core");
    }

    #[test]
    fn contiguous_core_coverage_short_domain_sharer_is_low() {
        // A "domain-sharer": shares only a short ~30 bp block, then both
        // sequences are otherwise independent random (long, so 30/min is small).
        let domain = core_rand_seq(30, 0xD0D0_0001);
        let mut a = core_rand_seq(300, 0xAAAA_1001);
        a.extend_from_slice(&domain);
        a.extend(core_rand_seq(300, 0xAAAA_1002));
        let mut b = core_rand_seq(300, 0xBBBB_1001);
        b.extend_from_slice(&domain);
        b.extend(core_rand_seq(300, 0xBBBB_1002));
        let cov = contiguous_core_coverage(&a, &b);
        assert!(cov < 0.13,
            "short shared block should give LOW contiguous-core coverage (got {cov:.3})");
    }

    #[test]
    fn contiguous_core_coverage_identical_is_one() {
        let s = core_rand_seq(200, 0x1234_5678);
        let cov = contiguous_core_coverage(&s, &s);
        assert!(cov >= 0.99,
            "identical sequences should give contiguous-core coverage ~1.0 (got {cov:.3})");
    }

    #[test]
    fn contiguous_core_coverage_disjoint_is_near_zero() {
        // Two independent random sequences: no long shared run, only short
        // chance-match runs -> near zero.
        let a = core_rand_seq(300, 0x1111_0001);
        let b = core_rand_seq(300, 0x2222_0001);
        let cov = contiguous_core_coverage(&a, &b);
        assert!(cov < 0.13,
            "disjoint random sequences should give near-zero contiguous-core coverage (got {cov:.3})");
    }

    #[test]
    fn longest_common_substring_basic() {
        assert_eq!(longest_common_substring(b"ABCDE", b"XBCDY"), 3, "BCD");
        assert_eq!(longest_common_substring(b"AAAA", b"AAAA"), 4, "identical");
        assert_eq!(longest_common_substring(b"ABC", b"XYZ"), 0, "disjoint");
        assert_eq!(longest_common_substring(b"", b"ABC"), 0, "empty");
        // asymmetric: the short string is a contiguous substring of the long one (the read-through case).
        assert_eq!(longest_common_substring(b"GATTACA", b"TTTGATTACAGGG"), 7, "short ⊂ long");
        assert_eq!(longest_common_substring(b"TTTGATTACAGGG", b"GATTACA"), 7, "order-independent");
    }

    #[test]
    fn longest_common_substring_matches_naive_on_random() {
        // cross-check the suffix-automaton LCS against an O(n^2) naive on small random strings.
        fn naive(a: &[u8], b: &[u8]) -> usize {
            let mut best = 0;
            for i in 0..a.len() {
                for j in 0..b.len() {
                    let mut k = 0;
                    while i + k < a.len() && j + k < b.len() && a[i + k] == b[j + k] {
                        k += 1;
                    }
                    if k > best {
                        best = k;
                    }
                }
            }
            best
        }
        for seed in 0..8u64 {
            let a = core_rand_seq(120, 0x5151_0000 ^ seed);
            let b = core_rand_seq(140, 0x6262_0000 ^ (seed * 7));
            assert_eq!(longest_common_substring(&a, &b), naive(&a, &b), "seed {seed}");
        }
    }

    #[test]
    fn bounded_core_coverage_below_cap_is_exact_poasta() {
        // a true-copy pair (400 bp exact core + divergent flanks); below the cap the bounded form must equal
        // the exact poasta path bit-for-bit (the validated small-transcript behaviour is untouched).
        let core = core_rand_seq(400, 0xC0FE_4242);
        let mut a = core_rand_seq(80, 0xAAAA_4242);
        a.extend_from_slice(&core);
        a.extend(core_rand_seq(80, 0xAAAA_4243));
        let mut b = core_rand_seq(80, 0xBBBB_4242);
        b.extend_from_slice(&core);
        b.extend(core_rand_seq(80, 0xBBBB_4243));
        let poasta = contiguous_core_coverage(&a, &b);
        assert_eq!(
            contiguous_core_coverage_bounded(&a, &b, 100_000),
            poasta,
            "below cap = exact poasta"
        );
        // above the cap the LCS fallback recovers the SAME 400 bp core fraction (no OOM, no loss).
        let lcs = contiguous_core_coverage_bounded(&a, &b, 10);
        assert!(lcs >= 0.6, "fallback recovers the 400 bp core fraction (got {lcs:.3})");
        assert!(
            (lcs - poasta).abs() < 0.1,
            "fallback ~ poasta on a clean exact core (lcs {lcs:.3} vs poasta {poasta:.3})"
        );
    }

    #[test]
    fn bounded_core_coverage_fallback_separates_copies_from_disjoint() {
        // the fallback must preserve the family/non-family separation poasta gives. Two LONG disjoint random
        // sequences (forcing the fallback) share only chance substrings -> well under the 0.13 family bar.
        let a = core_rand_seq(5000, 0x1111_9999);
        let b = core_rand_seq(5000, 0x2222_9999);
        assert!(
            contiguous_core_coverage_bounded(&a, &b, 100) < 0.13,
            "disjoint long sequences stay below the family bar under the fallback"
        );
        // a long read-through that CONTAINS a 2 kb copy as a substring -> high coverage of the copy (the
        // DSFAM45 hub case: the copy joins the family via the bounded fallback instead of OOMing).
        let copy = core_rand_seq(2000, 0x7777_0001);
        let mut hub = core_rand_seq(3000, 0x8888_0001);
        hub.extend_from_slice(&copy);
        hub.extend(core_rand_seq(3000, 0x8888_0002));
        assert!(
            contiguous_core_coverage_bounded(&copy, &hub, 100) >= 0.9,
            "a copy embedded in a long read-through hub is confirmed via the fallback"
        );
    }
}

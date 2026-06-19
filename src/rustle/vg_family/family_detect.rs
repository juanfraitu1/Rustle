//! Strand-aware de-novo multi-copy FAMILY DETECTION — the `bench/denovo_families.py` core.
//!
//! Annotation-free, minimizer-free. Given the de-novo assembled transcripts (built by the integration
//! layer from BAM/FASTA), find which gene loci form multi-copy families:
//!
//!   1. **Collapse isoforms → gene loci** by shared intron junctions (union-find on identical
//!      `(chrom, donor, acceptor)`). NOT raw span overlap — dense genes overlap transitively and a
//!      span-merge chains a whole chromosome into one bogus locus.
//!   2. **Canonical exact-k-mer ownership pre-filter** (the "counting-bloom done exactly"): a k-mer owned
//!      by `[cnt_min, cnt_max]` distinct reps is *family-informative*; a rep with `>= k_share` informative
//!      k-mers is a candidate. Single-copy genes (unique k-mers) are rejected here and never reach POA.
//!   3. **Contiguous-span pre-filter**: keep a candidate pair only if its shared informative k-mers SPAN
//!      `>= t_core * min(len)` (a cheap, recall-safe proxy for POA's contiguous core — only pre-rejects
//!      pairs POA would reject).
//!   4. **POA contiguous-core grading is THE criterion** (`contiguous_core_coverage >= t_core`), strand-aware
//!      via a reverse-complement fallback (copies assembled on opposite strands). Confirmed pairs are edges.
//!
//! The k-mer encoder is shared with `family_rescue` (canonical base-4 KMER=18); the POA primitive is the
//! already-ported `family_graph::contiguous_core_coverage`. The BAM/FASTA orchestration that builds the
//! `DenovoTranscript` records and parallelises the POA is the integration layer.

use std::collections::{BTreeMap, BTreeSet};

use super::family_graph::contiguous_core_coverage;
use super::family_rescue::window_canon_code;
use crate::vg::reverse_complement;

/// Canonical k-mer length (matches `family_rescue::KMER` and `denovo_families.py::KMER`).
pub const KMER: usize = 18;
/// POA contiguous-core coverage threshold to confirm a family edge.
pub const T_CORE: f64 = 0.13;
/// Skip POA on pairs whose shorter sequence exceeds this (cost guard).
pub const LEN_CAP: usize = 20_000;
/// A k-mer must be owned by `>= CNT_MIN` distinct reps to be family-informative.
pub const CNT_MIN: usize = 2;
/// Drop pervasive k-mers owned by `> CNT_MAX` reps (mobile-element/repeat, not family-specific).
pub const CNT_MAX: usize = 40;
/// Skip pair generation from an inverted-index bucket owning `> PAIR_CAP` reps (cost guard).
pub const PAIR_CAP: usize = 40;
/// Propose a candidate pair only if the two reps co-own `>= K_SHARE` informative k-mers.
pub const K_SHARE: usize = 6;
/// Hard guard on the distinct-pair set (prevents OOM at genome scale).
pub const MAX_PAIRS: usize = 8_000_000;

/// A de-novo assembled transcript (one isoform). Built by the integration layer from the BAM/FASTA;
/// this module operates purely on these in-memory records.
#[derive(Clone, Debug)]
pub struct DenovoTranscript {
    pub tid: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_reads: u32,
    /// Intron `(donor, acceptor)` genomic coordinates.
    pub introns: Vec<(u64, u64)>,
    /// Spliced sequence in transcription orientation.
    pub seq: Vec<u8>,
}

/// Tunable detection parameters (defaults mirror `denovo_families.py`).
#[derive(Clone, Copy, Debug)]
pub struct DetectParams {
    pub cnt_min: usize,
    pub cnt_max: usize,
    pub pair_cap: usize,
    pub k_share: usize,
    pub t_core: f64,
    pub len_cap: usize,
    pub max_pairs: usize,
}

impl Default for DetectParams {
    fn default() -> Self {
        DetectParams {
            cnt_min: CNT_MIN,
            cnt_max: CNT_MAX,
            pair_cap: PAIR_CAP,
            k_share: K_SHARE,
            t_core: T_CORE,
            len_cap: LEN_CAP,
            max_pairs: MAX_PAIRS,
        }
    }
}

/// Canonical k-mer codes of `seq` with the FIRST-occurrence position of each in N-FILTERED window space.
/// Mirrors `np.unique(kmer_hashes(seq), return_index=True)`: python's `kmer_hashes` DROPS N-touching
/// windows BEFORE indexing, so positions are compacted over surviving windows (a dropped window leaves no
/// gap in the counter). For all-ACGT input this equals the absolute window index. `BTreeMap` keeps the
/// keys sorted for deterministic iteration.
pub fn canonical_kmer_first_pos(seq: &[u8]) -> BTreeMap<u64, u32> {
    let mut map = BTreeMap::new();
    if seq.len() < KMER {
        return map;
    }
    let mut pos = 0u32; // compacted index: only advances on SURVIVING (non-N) windows
    for w in seq.windows(KMER) {
        if let Some(code) = window_canon_code(w) {
            map.entry(code).or_insert(pos); // keep the FIRST occurrence
            pos += 1;
        }
    }
    map
}

/// Iterative path-halving union-find (matches `annotation_families::find`).
fn uf_find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}
fn uf_union(parent: &mut [usize], a: usize, b: usize) {
    let ra = uf_find(parent, a);
    let rb = uf_find(parent, b);
    if ra != rb {
        parent[ra] = rb;
    }
}

/// (1) Collapse isoforms to GENE LOCI by shared intron junctions (union-find on identical
/// `(chrom, donor, acceptor)`). Returns the rep INDEX (into `transcripts`) for each locus, sorted
/// ascending. Rep = most reads, tie-break longest span, then earliest index. Deterministic.
pub fn collapse_loci(transcripts: &[DenovoTranscript]) -> Vec<usize> {
    let n = transcripts.len();
    let mut parent: Vec<usize> = (0..n).collect();
    // (chrom, donor, acceptor) -> first transcript index that used it.
    let mut junc_owner: BTreeMap<(&str, u64, u64), usize> = BTreeMap::new();
    for (i, t) in transcripts.iter().enumerate() {
        for &(d, a) in &t.introns {
            match junc_owner.get(&(t.chrom.as_str(), d, a)) {
                Some(&owner) => uf_union(&mut parent, i, owner),
                None => {
                    junc_owner.insert((t.chrom.as_str(), d, a), i);
                }
            }
        }
    }
    // group indices by component root (members in ascending index order).
    let mut comp: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = uf_find(&mut parent, i);
        comp.entry(r).or_default().push(i);
    }
    // rep = max by (n_reads, span); on a tie the EARLIEST index (python `max` keeps the first maximal).
    let mut reps: Vec<usize> = comp
        .into_values()
        .map(|members| {
            *members
                .iter()
                .max_by(|&&a, &&b| {
                    let ka = (transcripts[a].n_reads, transcripts[a].end - transcripts[a].start);
                    let kb = (transcripts[b].n_reads, transcripts[b].end - transcripts[b].start);
                    // break key ties by preferring the smaller index (so it is the "maximum").
                    ka.cmp(&kb).then_with(|| b.cmp(&a))
                })
                .unwrap()
        })
        .collect();
    reps.sort_unstable();
    reps
}

/// (2) Candidate homologous rep pairs. Exact canonical-k-mer ownership pre-filter (family-informative =
/// owned by `[cnt_min, cnt_max]` reps; a rep needs `>= k_share` informative k-mers) + contiguous-span
/// filter (shared informative-k-mer position span `>= t_core * min(len)`). Returns `(i, j)` index pairs
/// into `reps` with `i < j`, sorted.
pub fn candidate_pairs(reps: &[DenovoTranscript], p: &DetectParams) -> Vec<(usize, usize)> {
    let n = reps.len();
    // per-rep canonical k-mers with first-occurrence positions.
    let rep_kmers: Vec<BTreeMap<u64, u32>> =
        reps.iter().map(|r| canonical_kmer_first_pos(&r.seq)).collect();

    // exact ownership: code -> number of distinct reps owning it.
    let mut owner_count: BTreeMap<u64, usize> = BTreeMap::new();
    for km in &rep_kmers {
        for code in km.keys() {
            *owner_count.entry(*code).or_insert(0) += 1;
        }
    }
    // family-informative k-mers: owned by [cnt_min, cnt_max] reps.
    let info: BTreeSet<u64> = owner_count
        .iter()
        .filter(|(_, &c)| c >= p.cnt_min && c <= p.cnt_max)
        .map(|(&code, _)| code)
        .collect();

    // per-rep informative signature (code -> pos); a rep is a candidate iff it has >= k_share of them.
    let mut sig: Vec<BTreeMap<u64, u32>> = Vec::with_capacity(n);
    let mut is_candidate = vec![false; n];
    for i in 0..n {
        let s: BTreeMap<u64, u32> = rep_kmers[i]
            .iter()
            .filter(|(code, _)| info.contains(*code))
            .map(|(&c, &pos)| (c, pos))
            .collect();
        is_candidate[i] = s.len() >= p.k_share;
        sig.push(s);
    }

    // inverted index over informative k-mers of CANDIDATE reps (bounded -> no OOM).
    let mut inv: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        if !is_candidate[i] {
            continue;
        }
        for code in sig[i].keys() {
            inv.entry(*code).or_default().push(i);
        }
    }
    // distinct co-occurring candidate pairs; skip pervasive buckets (> pair_cap); MAX_PAIRS guard.
    let mut pair_set: BTreeSet<(usize, usize)> = BTreeSet::new();
    'outer: for lst in inv.values() {
        if lst.len() < 2 || lst.len() > p.pair_cap {
            continue;
        }
        for a in 0..lst.len() {
            for b in (a + 1)..lst.len() {
                let (ri, rj) = (lst[a], lst[b]);
                pair_set.insert((ri.min(rj), ri.max(rj)));
                if pair_set.len() > p.max_pairs {
                    break 'outer;
                }
            }
        }
    }

    // contiguous-span filter: shared informative-k-mer positions must SPAN >= t_core * min(len) in BOTH
    // reps (the shorter span). A real copy shares a contiguous block; a domain-sharer a short one.
    let mut out: Vec<(usize, usize)> = Vec::new();
    for &(a, b) in &pair_set {
        let (sa, sb) = (&sig[a], &sig[b]);
        let (amin, amax, bmin, bmax, common) = {
            let (mut amin, mut amax, mut bmin, mut bmax, mut common) =
                (u32::MAX, 0u32, u32::MAX, 0u32, 0usize);
            for (code, &pa) in sa {
                if let Some(&pb) = sb.get(code) {
                    common += 1;
                    amin = amin.min(pa);
                    amax = amax.max(pa);
                    bmin = bmin.min(pb);
                    bmax = bmax.max(pb);
                }
            }
            (amin, amax, bmin, bmax, common)
        };
        if common < p.k_share {
            continue;
        }
        let core = (amax - amin).min(bmax - bmin) as usize;
        let minlen = reps[a].seq.len().min(reps[b].seq.len());
        if core as f64 >= p.t_core * minlen as f64 {
            out.push((a, b));
        }
    }
    out.sort_unstable();
    out
}

/// (3) POA edge confirmation: `contiguous_core_coverage` in both orientations (reverse-complement
/// fallback for opposite-strand copies), `Some(core_recip)` iff `>= t_core`. `len_cap` guard. Uppercases
/// the POA operands (the soft-mask lesson from `family_rescue`: `reverse_complement` maps lowercase → N).
pub fn confirm_edge(a: &[u8], b: &[u8], p: &DetectParams) -> Option<f64> {
    if a.len().min(b.len()) > p.len_cap {
        return None;
    }
    let au = a.to_ascii_uppercase();
    let bu = b.to_ascii_uppercase();
    let mut cr = contiguous_core_coverage(&au, &bu);
    if cr < p.t_core {
        // opposite orientation (copies assembled on different strands).
        let rc = contiguous_core_coverage(&au, &reverse_complement(&bu));
        if rc > cr {
            cr = rc;
        }
    }
    if cr >= p.t_core {
        Some(cr)
    } else {
        None
    }
}

/// Convenience driver: `candidate_pairs` then `confirm_edge` over `reps`. Returns confirmed
/// `(i, j, core_recip)` edges (`i < j`).
///
/// The POA confirmation (the expensive step) runs in PARALLEL via rayon — `confirm_edge` is independent
/// per pair and `contiguous_core_coverage` is pure/deterministic, so this matches the python's
/// process-pool POA. Order is preserved (rayon's indexed `collect`), so the edge list is deterministic
/// and identical to the serial version (candidate-pair order).
pub fn detect_edges(reps: &[DenovoTranscript], p: &DetectParams) -> Vec<(usize, usize, f64)> {
    use rayon::prelude::*;
    candidate_pairs(reps, p)
        .into_par_iter()
        .filter_map(|(a, b)| confirm_edge(&reps[a].seq, &reps[b].seq, p).map(|cr| (a, b, cr)))
        .collect()
}

#[cfg(test)]
mod tests {
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
    fn cat(parts: &[&[u8]]) -> Vec<u8> {
        parts.iter().flat_map(|p| p.iter().copied()).collect()
    }
    fn tx(
        tid: &str,
        chrom: &str,
        start: u64,
        end: u64,
        n_reads: u32,
        introns: &[(u64, u64)],
        seq: Vec<u8>,
    ) -> DenovoTranscript {
        DenovoTranscript {
            tid: tid.into(),
            chrom: chrom.into(),
            start,
            end,
            n_reads,
            introns: introns.to_vec(),
            seq,
        }
    }
    fn homolog_tx_flank(tid: &str, fs1: u64, core: &[u8], fs2: u64, flank: usize, n: u32) -> DenovoTranscript {
        let seq = cat(&[&rand_seq(flank, fs1), core, &rand_seq(flank, fs2)]);
        let len = seq.len() as u64;
        tx(tid, "c1", 0, len, n, &[(10, 20)], seq)
    }
    fn homolog_tx(tid: &str, fs1: u64, core: &[u8], fs2: u64, n: u32) -> DenovoTranscript {
        homolog_tx_flank(tid, fs1, core, fs2, 80, n)
    }

    // ---- canonical_kmer_first_pos ----

    #[test]
    fn first_pos_are_window_indices() {
        let s = rand_seq(40, 0xABCD);
        let m = canonical_kmer_first_pos(&s);
        let nw = (40 - KMER + 1) as u32;
        assert_eq!(m.len(), nw as usize, "distinct random 18-mers");
        let mut positions: Vec<u32> = m.values().copied().collect();
        positions.sort_unstable();
        assert_eq!(positions, (0..nw).collect::<Vec<_>>());
    }

    #[test]
    fn first_pos_keeps_earliest_for_repeat() {
        // tandem repeat of a 20 bp unit -> 18-mers recur each period; first occurrence kept.
        let unit = rand_seq(20, 0xBEEF);
        let s = cat(&[&unit, &unit, &unit]);
        let m = canonical_kmer_first_pos(&s);
        let code0 = window_canon_code(&s[0..KMER]).unwrap();
        assert_eq!(m.get(&code0), Some(&0));
    }

    #[test]
    fn first_pos_compacts_over_dropped_n_windows() {
        // Positions are indices in N-FILTERED window space (python np.unique over kmer_hashes, which DROPS
        // N-touching windows BEFORE indexing). An internal N must NOT leave a gap in the position counter.
        // seq: 5 clean bases, an N at index 5, then a clean tail -> windows starting 0..=5 touch the N and
        // are dropped; the first SURVIVING window starts at absolute index 6 but must get COMPACTED pos 0.
        let mut s = rand_seq(5, 0x9999);
        s.push(b'N');
        s.extend_from_slice(&rand_seq(40, 0x1234));
        let m = canonical_kmer_first_pos(&s);
        let first_clean = window_canon_code(&s[6..6 + KMER]).unwrap();
        assert_eq!(m.get(&first_clean), Some(&0), "first surviving window -> compacted position 0");
    }

    #[test]
    fn first_pos_compacted_equals_absolute_without_n() {
        // no N -> compacted index == absolute window index (existing all-ACGT behaviour preserved).
        let s = rand_seq(35, 0x4321);
        let m = canonical_kmer_first_pos(&s);
        let mut ps: Vec<u32> = m.values().copied().collect();
        ps.sort_unstable();
        assert_eq!(ps, (0..(35 - KMER + 1) as u32).collect::<Vec<_>>());
    }

    // ---- collapse_loci ----

    #[test]
    fn collapse_groups_shared_junction() {
        let txs = [
            tx("A", "c1", 100, 500, 5, &[(200, 300)], rand_seq(400, 1)),
            tx("B", "c1", 100, 520, 9, &[(200, 300)], rand_seq(420, 2)),
        ];
        let reps = collapse_loci(&txs);
        assert_eq!(reps.len(), 1);
        assert_eq!(reps[0], 1, "rep is B (more reads)");
    }

    #[test]
    fn collapse_separates_distinct_junctions() {
        let txs = [
            tx("A", "c1", 100, 500, 5, &[(200, 300)], rand_seq(400, 1)),
            tx("B", "c1", 1000, 1500, 5, &[(1200, 1300)], rand_seq(400, 2)),
        ];
        assert_eq!(collapse_loci(&txs).len(), 2);
    }

    #[test]
    fn collapse_single_exon_are_singletons() {
        let txs = [
            tx("A", "c1", 100, 500, 5, &[], rand_seq(400, 1)),
            tx("B", "c1", 100, 500, 5, &[], rand_seq(400, 2)),
        ];
        assert_eq!(collapse_loci(&txs).len(), 2, "no introns -> never merged");
    }

    #[test]
    fn collapse_transitive_chain() {
        let txs = [
            tx("A", "c1", 100, 900, 3, &[(200, 300)], rand_seq(400, 1)),
            tx("B", "c1", 100, 900, 7, &[(200, 300), (600, 700)], rand_seq(400, 2)),
            tx("C", "c1", 100, 900, 5, &[(600, 700)], rand_seq(400, 3)),
        ];
        let reps = collapse_loci(&txs);
        assert_eq!(reps.len(), 1, "A-B and B-C share junctions -> one locus");
        assert_eq!(reps[0], 1, "rep is B (most reads)");
    }

    #[test]
    fn collapse_junction_requires_same_chrom() {
        let txs = [
            tx("A", "c1", 100, 500, 5, &[(200, 300)], rand_seq(400, 1)),
            tx("B", "c2", 100, 500, 5, &[(200, 300)], rand_seq(400, 2)),
        ];
        assert_eq!(collapse_loci(&txs).len(), 2, "same coords, different chrom -> not merged");
    }

    #[test]
    fn collapse_rep_tiebreak_longest_span() {
        let txs = [
            tx("A", "c1", 100, 500, 5, &[(200, 300)], rand_seq(400, 1)),
            tx("B", "c1", 100, 900, 5, &[(200, 300)], rand_seq(800, 2)),
        ];
        assert_eq!(collapse_loci(&txs), vec![1], "equal reads -> longest span wins");
    }

    #[test]
    fn collapse_rep_tiebreak_earliest_on_full_tie() {
        // identical n_reads AND identical span -> the EARLIEST index wins (python max keeps first maximal).
        let txs = [
            tx("A", "c1", 100, 500, 5, &[(200, 300)], rand_seq(400, 1)),
            tx("B", "c1", 100, 500, 5, &[(200, 300)], rand_seq(400, 2)),
        ];
        assert_eq!(collapse_loci(&txs), vec![0], "full tie -> smallest index");
    }

    // ---- candidate_pairs ----

    #[test]
    fn candidate_pairs_finds_homologous_reps() {
        let core = rand_seq(400, 0x0C0FE);
        let reps = [
            homolog_tx("r0", 0xA1, &core, 0xA2, 5),
            homolog_tx("r1", 0xB1, &core, 0xB2, 5),
            tx("r2", "c1", 0, 560, 5, &[], rand_seq(560, 0xDEAD)),
        ];
        assert_eq!(candidate_pairs(&reps, &DetectParams::default()), vec![(0, 1)]);
    }

    #[test]
    fn candidate_pairs_rejects_all_single_copy() {
        let reps = [
            tx("r0", "c1", 0, 560, 5, &[], rand_seq(560, 0x01)),
            tx("r1", "c1", 0, 560, 5, &[], rand_seq(560, 0x02)),
            tx("r2", "c1", 0, 560, 5, &[], rand_seq(560, 0x03)),
        ];
        assert!(candidate_pairs(&reps, &DetectParams::default()).is_empty());
    }

    #[test]
    fn candidate_pairs_span_filter_rejects_short_domain_sharer() {
        // 40 bp shared block -> 23 shared 18-mers (>= k_share, pre-filter PASSES), but the span (~22) is
        // < t_core*min(len) (~109 for 840 bp transcripts) -> contiguous-span filter rejects.
        let block = rand_seq(40, 0xD0D0);
        let reps = [
            homolog_tx_flank("r0", 0xA1, &block, 0xA2, 400, 5),
            homolog_tx_flank("r1", 0xB1, &block, 0xB2, 400, 5),
        ];
        assert!(candidate_pairs(&reps, &DetectParams::default()).is_empty());
    }

    #[test]
    fn candidate_pairs_cnt_max_drops_pervasive_kmers() {
        let core = rand_seq(400, 0xCAFE);
        let reps = [
            homolog_tx("r0", 0xA1, &core, 0xA2, 5),
            homolog_tx("r1", 0xB1, &core, 0xB2, 5),
            homolog_tx("r2", 0xC1, &core, 0xC2, 5),
        ];
        // default cnt_max=40: core owned by 3 reps -> informative -> all three pairs found.
        assert_eq!(
            candidate_pairs(&reps, &DetectParams::default()),
            vec![(0, 1), (0, 2), (1, 2)]
        );
        // cnt_max=2: core k-mers owned by 3 > 2 -> NOT informative -> no candidates.
        let p = DetectParams { cnt_max: 2, ..DetectParams::default() };
        assert!(candidate_pairs(&reps, &p).is_empty());
    }

    #[test]
    fn candidate_pairs_cnt_max_upper_bound_isolated() {
        // exactly two reps share a 400 bp core -> the core k-mers are owned by exactly 2 reps. cnt_min=1
        // makes the flanks (owned by 1) informative too, so this toggles ONLY the cnt_max upper bound:
        // cnt_max=1 excludes the core (owned 2 > 1) -> only singleton flank buckets -> no pair; cnt_max=2
        // includes it -> the pair is proposed.
        let core = rand_seq(400, 0x5A5A);
        let reps = [
            homolog_tx("r0", 0xA1, &core, 0xA2, 5),
            homolog_tx("r1", 0xB1, &core, 0xB2, 5),
        ];
        let excl = DetectParams { cnt_min: 1, cnt_max: 1, ..DetectParams::default() };
        assert!(candidate_pairs(&reps, &excl).is_empty(), "core owned by 2 > cnt_max=1 -> excluded");
        let incl = DetectParams { cnt_min: 1, cnt_max: 2, ..DetectParams::default() };
        assert_eq!(candidate_pairs(&reps, &incl), vec![(0, 1)]);
    }

    #[test]
    fn candidate_pairs_strand_symmetric_via_canonical_kmers() {
        // r1 is the reverse-complement of r0's homolog: the canonical (min(fwd,rc)) encoder still matches,
        // so the opposite-strand copy survives the k-mer pre-filter and the pair is proposed (this is what
        // feeds confirm_edge's RC fallback in the real pipeline).
        let core = rand_seq(400, 0x57A7);
        let r0 = homolog_tx("r0", 0xA1, &core, 0xA2, 5);
        let fwd = cat(&[&rand_seq(80, 0xB1), &core, &rand_seq(80, 0xB2)]);
        let rc = reverse_complement(&fwd);
        let r1 = tx("r1", "c1", 0, rc.len() as u64, 5, &[(10, 20)], rc);
        assert_eq!(candidate_pairs(&[r0, r1], &DetectParams::default()), vec![(0, 1)]);
    }

    // ---- confirm_edge ----

    #[test]
    fn confirm_edge_confirms_homologous() {
        let core = rand_seq(400, 0xC0FE_9001);
        let a = cat(&[&rand_seq(80, 0xA1), &core, &rand_seq(80, 0xA2)]);
        let b = cat(&[&rand_seq(80, 0xB1), &core, &rand_seq(80, 0xB2)]);
        let cr = confirm_edge(&a, &b, &DetectParams::default()).expect("homologous pair confirms");
        assert!(cr >= T_CORE, "core_recip {cr} >= {T_CORE}");
    }

    #[test]
    fn confirm_edge_rejects_disjoint() {
        let a = rand_seq(560, 0x111);
        let b = rand_seq(560, 0x222);
        assert!(confirm_edge(&a, &b, &DetectParams::default()).is_none());
    }

    #[test]
    fn confirm_edge_reverse_complement() {
        let core = rand_seq(400, 0xC0FE_9301);
        let a = cat(&[&rand_seq(80, 0xA1), &core, &rand_seq(80, 0xA2)]);
        let bfwd = cat(&[&rand_seq(80, 0xB1), &core, &rand_seq(80, 0xB2)]);
        let b = reverse_complement(&bfwd);
        let cr = confirm_edge(&a, &b, &DetectParams::default()).expect("opposite-strand copy confirms");
        assert!(cr >= T_CORE);
    }

    #[test]
    fn confirm_edge_uppercases_softmasked_input() {
        let core = rand_seq(400, 0xC0FE_9401);
        let a = cat(&[&rand_seq(80, 0xA1), &core, &rand_seq(80, 0xA2)]).to_ascii_lowercase();
        let bfwd = cat(&[&rand_seq(80, 0xB1), &core, &rand_seq(80, 0xB2)]);
        let b = reverse_complement(&bfwd).to_ascii_lowercase();
        let cr = confirm_edge(&a, &b, &DetectParams::default())
            .expect("lowercase opposite-strand copy must still confirm");
        assert!(cr >= T_CORE);
    }

    #[test]
    fn confirm_edge_len_cap_guard() {
        let p = DetectParams { len_cap: 100, ..DetectParams::default() };
        let a = rand_seq(560, 0x1);
        let b = rand_seq(560, 0x2);
        assert!(confirm_edge(&a, &b, &p).is_none(), "min(len) 560 > len_cap 100 -> skipped");
    }

    // ---- detect_edges ----

    #[test]
    fn detect_edges_end_to_end() {
        let core = rand_seq(400, 0xEE01);
        let reps = [
            homolog_tx("r0", 0xA1, &core, 0xA2, 5),
            homolog_tx("r1", 0xB1, &core, 0xB2, 5),
            tx("r2", "c1", 0, 560, 5, &[], rand_seq(560, 0xDEAD)),
        ];
        let edges = detect_edges(&reps, &DetectParams::default());
        assert_eq!(edges.len(), 1);
        assert_eq!((edges[0].0, edges[0].1), (0, 1));
        assert!(edges[0].2 >= T_CORE);
    }

    #[test]
    fn detect_edges_is_deterministic_under_parallelism() {
        // three paralogs sharing one core -> three candidate pairs -> three edges; the PARALLEL POA must
        // return them in candidate-pair order, identically across runs.
        let core = rand_seq(400, 0xEE77);
        let reps = [
            homolog_tx("r0", 0xA1, &core, 0xA2, 5),
            homolog_tx("r1", 0xB1, &core, 0xB2, 5),
            homolog_tx("r2", 0xC1, &core, 0xC2, 5),
        ];
        let e1 = detect_edges(&reps, &DetectParams::default());
        let e2 = detect_edges(&reps, &DetectParams::default());
        assert_eq!(
            e1.iter().map(|&(a, b, _)| (a, b)).collect::<Vec<_>>(),
            vec![(0, 1), (0, 2), (1, 2)],
            "edges in sorted candidate-pair order"
        );
        assert_eq!(e1, e2, "parallel POA must be order-deterministic");
    }
}

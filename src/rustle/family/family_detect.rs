//! Strand-aware de-novo multi-copy FAMILY DETECTION — the `bench/denovo_families.py` core.
//!
//! **Terminology.** This module uses "locus" in the **gene-locus** sense: a set of
//! isoforms collapsed by shared splice junctions. It is NOT the physical `(chrom,
//! start, end)` span used for the ≥2-distinct-loci certificate (that is
//! `family_definition::distinct_loci`). See `docs/REFERENCE.md` for the
//! canonical vocabulary.
//!
//! Annotation-free, minimizer-free. Given the de-novo assembled transcripts (built by the integration
//! layer from BAM/FASTA), find which gene loci form multi-copy families:
//!
//!   1. **Collapse isoforms → gene loci** by shared intron junctions (union-find on identical
//!      `(chrom, donor, acceptor)`). NOT raw span overlap — dense genes overlap transitively and a
//!      span-merge chains a whole chromosome into one bogus locus. A span-aware recovery pass can
//!      additionally merge isoforms with disjoint junctions but strong span homology/containment
//!      (`collapse_loci_span_aware`).
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
//!
//! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

use crate::types::{DetHashMap, DetHashSet};
use std::collections::{BTreeMap, BTreeSet};

use super::family_rescue::window_canon_code;
use crate::family::seq_utils::reverse_complement;
use family_graph::{core_coverage_reaches, upper_cow};

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
/// For span-aware locus collapse: a shorter transcript must be contained in the longer by at least this
/// fraction of its own length to be merged as a same-gene isoform without running POA.
pub const COLLAPSE_CONTAIN_FRAC: f64 = 0.5;
/// For span-aware locus collapse: minimum POA contiguous-core coverage for two span-overlapping,
/// same-strand transcripts with disjoint junctions to be merged as same-gene isoforms. Conservative
/// (well above the family-detection `T_CORE` of 0.13) to avoid merging adjacent paralogs.
pub const COLLAPSE_SPAN_CORE: f64 = 0.50;

/// A de-novo assembled transcript (one isoform). Built by the integration layer from the BAM/FASTA;
/// this module operates purely on these in-memory records.
#[derive(Clone, Debug)]
pub struct DenovoTranscript {
    pub tid: String,
    pub chrom: String,
    pub start: u64,
    pub end: u64,
    pub n_reads: u32,
    /// Transcription strand (`'+'`/`'-'`), from the gate's canonical-junction classification. The `seq` is
    /// in this orientation; copy assignment needs it for the spliced↔genomic map and read-base orientation.
    pub strand: char,
    /// Intron `(donor, acceptor)` genomic coordinates.
    pub introns: Vec<(u64, u64)>,
    /// Spliced sequence in transcription orientation.
    pub seq: Vec<u8>,
    /// Per-copy unique-mapper support: the number of reads that place at this copy's locus UNIQUELY
    /// (`mapq > 0` — the aligner found no competing tied placement). This is the χ(H) read-distinguishability
    /// signal `distinct_locus_reps`'s same-strand merge guard consumes via `read_conflict::reads_distinguish`
    /// (Task 3, identifiability-merge): a co-located same-strand pair collapses only when NEITHER copy clears
    /// the `min_reads` floor here. Populated where BAM read placements already exist (`detect_and_assign`'s
    /// `placements`, and the genome-wide homology catalog's own placement pass before it drops its reads);
    /// 0 (the safe/collapsing default) elsewhere — synthetic/test copies, and admitted/rescued copies built
    /// after the merge decision has already run.
    pub distinguishing_uniq: usize,

    /// Bases of this copy's span carried by at least `RUSTLE_ER_CORE_DEPTH` reads — the READ-SUPPORTED
    /// CORE, as opposed to the called span. 0 means "not measured" and callers must fall back to the span.
    ///
    /// Why it exists: the E_r coverage test divides the aligned length by the SHORTER COPY'S SPAN, so the
    /// criterion that decides membership is measured against a quantity the pipeline itself produces, and
    /// both boundary errors corrupt it. An over-extended locus inflates its own denominator and LOSES real
    /// edges (measured: 14 of 521 true co-family pairs fail the 0.50 floor purely for this, NOTCH2NLB at
    /// 0.436 and AMY2A at 0.232 where true sizes give 0.939 and 1.000). A truncated locus shrinks it and
    /// passes trivially, which is why coverage never detected the 0.55x truncation problem at all.
    ///
    /// A readthrough tail is low-depth, so it never enters the core and cannot inflate the denominator.
    /// Measured edge-level on chr1+chr15: `aligned / min(core at depth >= 10) >= 0.80` gives precision
    /// 0.830 and recall 0.883 against the span rule's 0.822 / 0.845 — better on BOTH axes, where every
    /// denominator-free alternative (absolute aligned length, exonic coverage) was worse on at least one.
    pub core_bp: u64,

    /// This copy is a STUB: its representative is single-exon, but reads at its locus DO assert a splice
    /// junction, so the gene is spliced and the representative is a fragment of it. Distinct from a
    /// genuinely intronless copy (histone, retrocopy), which is also single-exon but has no spliced
    /// evidence at all -- measured on chr1+chr15, 324 of 432 single-exon copies are stubs and 108 are not.
    ///
    /// Why the distinction is worth carrying: a stub is short, so it matches a FRAGMENT of a real gene
    /// almost perfectly and can form an E_r edge on that basis. In the NPIP family every false positive was
    /// a single-exon stub joining through a short high-identity match -- and identity and coverage cannot
    /// reject them, because those pairs score HIGHER than the true NPIP pairs (coverage 0.979 vs 0.912,
    /// 13,992 aligned bp vs 2,923). The real members carried 4, 8 and 21 exons.
    ///
    /// Deleting stubs outright is worse than keeping them (F1 0.704 -> 0.608): 70% of copies are
    /// single-exon, so removal dissolves 65 of 95 families through the >=2-loci gate, and the families that
    /// die are disproportionately the small pure ones. The defensible use is to deny a stub the right to
    /// CREATE family membership while keeping it as a locus.
    pub stub: bool,

    /// Observed 3' terminus (TES) of this copy, when the reads' 3'-end distribution is SHARP; `None` when it
    /// is broad (differential coverage rather than a real polyadenylation site) or when no read evidence was
    /// available. Strand-aware: the genomic coordinate of the transcript's LAST base, so `>= end` on `+` and
    /// `<= start` on `-`.
    ///
    /// Exists because the 3' end carries ~5x more copy-discriminating signal than the 5' end
    /// (`bench/soto/bam_tie_signals.md` §9: 14/42 sibling pairs have distinct TES vs 3/40 for TSS, median
    /// shift 4.9 kb, up to 58 kb for truncated paralogs like NOTCH2 vs NOTCH2NLB). `refine_copy_seq` extends
    /// the terminal exon of the exon-sum to here, which is exactly the sequence that distinguishes a
    /// truncated duplicate from its parent and that the quantile boundary discards.
    pub tes: Option<u64>,
}

impl Default for DenovoTranscript {
    /// All-zero/empty — NOT a real transcript, only a base for struct-update syntax (`..Default::default()`)
    /// in call sites that specify every semantically-meaningful field explicitly.
    fn default() -> Self {
        DenovoTranscript {
            tid: String::new(),
            chrom: String::new(),
            start: 0,
            end: 0,
            n_reads: 0,
            strand: '+',
            introns: Vec::new(),
            seq: Vec::new(),
            distinguishing_uniq: 0,
            core_bp: 0,
            stub: false,
            tes: None,
        }
    }
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
    /// If true, `collapse_loci` additionally merges isoforms that share no junction but are
    /// span-overlapping and either strongly contained or highly homologous. Default true.
    pub collapse_span_aware: bool,
    /// POA contiguous-core threshold used by span-aware collapse for disjoint-junction isoforms.
    /// Conservative to avoid merging adjacent paralogs.
    pub collapse_span_core: f64,
    /// Which contiguous core `confirm_edge` measures. `Lcs` (default since docs/o1_ledger.md §6jd) = the longest
    /// exact common substring over `min(len)`, forward then reverse complement: linear time, no alignment.
    /// `Poa` = the global poasta alignment core (the pre-§6jd behaviour; escape hatch `RUSTLE_EDGE_CORE=poa`).
    /// Evidence, all pre-registered (`docs/PREREG_core_definition_2026-09-12.md`): human pairs LCS F1 0.861 /
    /// AUC 0.948 vs POA 0.688 / 0.796 (§6ja); human families ARI 0.681 vs 0.525 (§6jb); held-out gorilla
    /// precision 1.000 for both with recall 0.143 vs 0.063 (§6jb) and 0.140 vs 0.060 on fresh pairs (§6jc).
    pub edge_core: EdgeCore,
}

/// Contiguous-core definition used by `confirm_edge` (see [`DetectParams::edge_core`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EdgeCore {
    Poa,
    Lcs,
}

/// `RUSTLE_EDGE_CORE` value -> [`EdgeCore`]: `poa` (case-insensitive) selects the old poasta core, anything else
/// (including unset) is the LCS default.
pub fn edge_core_from_env_value(v: Option<&str>) -> EdgeCore {
    match v {
        Some(s) if s.eq_ignore_ascii_case("poa") => EdgeCore::Poa,
        _ => EdgeCore::Lcs,
    }
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
            collapse_span_aware: true,
            collapse_span_core: COLLAPSE_SPAN_CORE,
            edge_core: EdgeCore::Lcs,
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
                    let ka = (
                        transcripts[a].n_reads,
                        transcripts[a].end - transcripts[a].start,
                    );
                    let kb = (
                        transcripts[b].n_reads,
                        transcripts[b].end - transcripts[b].start,
                    );
                    // break key ties by preferring the smaller index (so it is the "maximum").
                    ka.cmp(&kb).then_with(|| b.cmp(&a))
                })
                .unwrap()
        })
        .collect();
    reps.sort_unstable();
    reps
}

/// Like `collapse_loci`, but returns the GENE rep index for EVERY transcript (transcript i belongs to the
/// gene whose representative is `groups[i]`) — so a FLAIR-style emitter can group isoforms under one
/// `gene_id`. Same union-find on shared `(chrom, donor, acceptor)` junctions and the same rep tie-break as
/// `collapse_loci`, so the chosen reps are identical; this just also reports the membership.
pub fn collapse_loci_groups(transcripts: &[DenovoTranscript]) -> Vec<usize> {
    let n = transcripts.len();
    let mut parent: Vec<usize> = (0..n).collect();
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
    let mut comp: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = uf_find(&mut parent, i);
        comp.entry(r).or_default().push(i);
    }
    let mut group = vec![0usize; n];
    for members in comp.values() {
        let rep = *members
            .iter()
            .max_by(|&&a, &&b| {
                let ka = (
                    transcripts[a].n_reads,
                    transcripts[a].end - transcripts[a].start,
                );
                let kb = (
                    transcripts[b].n_reads,
                    transcripts[b].end - transcripts[b].start,
                );
                ka.cmp(&kb).then_with(|| b.cmp(&a))
            })
            .unwrap();
        for &m in members {
            group[m] = rep;
        }
    }
    group
}

/// Compute the representative index for each union-find component in `parent`, using the same tie-break
/// as `collapse_loci`: most reads, then longest span, then earliest index. Returns rep indices sorted
/// ascending.
/// Locus GROUPS: the transcript indices collapsed onto each locus, in the same order `locus_reps` emits its
/// representatives. `locus_reps` discards this grouping by returning one index per locus; the exon-union
/// substrate needs it, because the union is taken over the whole group.
pub(crate) fn locus_groups(transcripts: &[DenovoTranscript], p: &DetectParams) -> Vec<Vec<usize>> {
    let parent = collapse_parent(transcripts, p);
    let mut find = |mut x: usize| {
        while parent[x] != x {
            x = parent[x];
        }
        x
    };
    let mut by_root: std::collections::BTreeMap<usize, Vec<usize>> =
        std::collections::BTreeMap::new();
    for i in 0..transcripts.len() {
        by_root.entry(find(i)).or_default().push(i);
    }
    by_root.into_values().collect()
}

/// Union the EXON intervals of a locus group into one geometry `(start, end, introns)`.
///
/// `pick_locus_rep` returns the single highest-support intron chain. In segmental duplications reads shatter
/// into many partial chains — NOTCH2: 490 reads over 240 distinct chains, the largest holding 5% — so that
/// one chain describes a FRAGMENT of the locus (measured at median 0.54x the truth span, 104 truncated vs 12
/// over-extended). Fragments of *different* genes then align to each other, which is how a single component
/// came to fuse 40 of 83 Soto families.
///
/// Only members with `>= min_chain_reads` contribute. That floor is the mechanism, not a tuning knob:
/// unioning EVERY read instead rebuilt the blob (39 families fused vs 2 with the floor).
///
/// The returned `introns` are the gaps BETWEEN merged exons, so they are not necessarily any single
/// transcript's junctions and may be chimeric — the caller must therefore build the sequence directly and
/// must NOT re-apply the canonical-junction gate, whose job was already done per-chain upstream.
pub(crate) fn union_locus_geometry(
    transcripts: &[DenovoTranscript],
    members: &[usize],
    min_chain_reads: u32,
) -> Option<(u64, u64, Vec<(u64, u64)>)> {
    union_locus_geometry_corroborated(transcripts, members, min_chain_reads, 0)
}

/// ⭐ §6y9/r1040: the union rep, restricted to CORROBORATED exons.
///
/// `RUSTLE_LOCUS_EXON_UNION` was refuted (register 303: Soto recall 65.5% -> 44.8%) because "widened reps
/// inflate their own coverage denominator and fall below the 0.50 floor". Measured on chr20
/// (`a119b_polished.gtf`, 505 multi-transcript loci), the union it builds is mostly noise: of the 3,483
/// exons a non-rep transcript contributes, **63.5% appear in exactly ONE transcript** while only **36.5%
/// appear in >= 2**. The existing `min_chain_reads` floor (default 3) does not separate these -- it admits
/// **74.8%** of them, because a single chain can carry plenty of reads.
///
/// `min_exon_tx` filters per EXON instead of per TRANSCRIPT: a merged interval survives only if at least
/// that many member transcripts place an exon on it. The REP's own exons are always kept, so the result is
/// never smaller than today's single-chain rep. `0` = off = the historical behaviour, byte-identical.
pub(crate) fn union_locus_geometry_corroborated(
    transcripts: &[DenovoTranscript],
    members: &[usize],
    min_chain_reads: u32,
    min_exon_tx: u32,
) -> Option<(u64, u64, Vec<(u64, u64)>)> {
    let mut ex: Vec<(u64, u64)> = Vec::new();
    for &i in members {
        if transcripts[i].n_reads >= min_chain_reads {
            ex.extend(exons_of(&transcripts[i]));
        }
    }
    if ex.is_empty() {
        return None;
    }
    ex.sort_unstable();
    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(ex.len());
    for (a, b) in ex {
        match merged.last_mut() {
            Some(last) if a <= last.1 => last.1 = last.1.max(b),
            _ => merged.push((a, b)),
        }
    }
    if min_exon_tx > 1 && !merged.is_empty() {
        // the locus rep, chosen exactly as `pick_locus_rep` does, so its exons are never dropped
        let rep_i = *members
            .iter()
            .max_by_key(|&&i| {
                (
                    transcripts[i].n_reads,
                    transcripts[i].end - transcripts[i].start,
                )
            })
            .expect("locus group is never empty");
        let rep_ex = exons_of(&transcripts[rep_i]);
        let overlaps =
            |iv: &(u64, u64), xs: &[(u64, u64)]| xs.iter().any(|e| e.0 < iv.1 && iv.0 < e.1);
        let kept: Vec<(u64, u64)> = merged
            .iter()
            .filter(|iv| {
                if overlaps(iv, &rep_ex) {
                    return true;
                }
                let n = members
                    .iter()
                    .filter(|&&i| transcripts[i].n_reads >= min_chain_reads)
                    .filter(|&&i| overlaps(iv, &exons_of(&transcripts[i])))
                    .count() as u32;
                n >= min_exon_tx
            })
            .copied()
            .collect();
        if !kept.is_empty() {
            merged = kept;
        }
    }
    let start = merged[0].0;
    let end = merged[merged.len() - 1].1;
    let introns: Vec<(u64, u64)> = merged.windows(2).map(|w| (w[0].1, w[1].0)).collect();
    Some((start, end, introns))
}

/// Build a locus geometry as the maximum-weight path through a READ-WITNESSED splice graph.
///
/// WHY NOT `pick_locus_rep`. It returns the single best-supported intron chain, and inside segmental
/// duplications there is no good chain to return: NOTCH2NLC's 156 spliced reads occupy 117 distinct chains.
/// `RUSTLE_SPLICED_REP` fixed the spliced-vs-unspliced comparison and was still net-marginal, because the
/// winning class must still nominate one OBSERVED chain. Shattering is combinatorial in chains and not in
/// junctions -- a read must match EVERY junction to join a chain, while each junction is supported
/// independently, and at that same locus 55 junctions carry >= 3 reads and 32 carry >= 10.
///
/// WHY NOT `union_locus_geometry`. Unioning every member's exons lost 20 recall points, because it
/// lengthens a representative and coverage is aligned/min(qlen,tlen), so the rep inflates its own
/// denominator. It also produces gaps between merged exons that no read asserts.
///
/// THE CONSTRUCTION. Nodes are the locus's junctions with >= `min_reads` support. An edge A -> B exists
/// only when some member transcript carries A immediately followed by B, so every step is witnessed by
/// reads; a transcript groups reads by EXACT intron chain, so its consecutive introns are co-observed by
/// construction. The representative is the maximum-weight path, computed exactly by one pass over the
/// junctions in coordinate order (edges only ever run forward, so that order is a topological one).
///
/// The co-observation constraint is what makes a low floor safe. Selecting a maximum-weight set of merely
/// NON-OVERLAPPING junctions has no locality: junctions from a neighbouring gene do not overlap either, so
/// the chain hops into it. Measured on the 52 stub-represented Soto members at a floor of 3 reads, the
/// unconstrained set put 35% of loci above 2x their true size (DNM1P50 reached 22x) while the co-observed
/// path put 27% there and raised correctly-sized loci from 22 to 26. The two arms converge at a floor of
/// 20, which is the proof that raising the floor never fixed the chaining -- it removed the loci where bad
/// chaining could happen.
///
/// Terminal exons come from evidence: `start`/`end` are the extreme bounds of the transcripts that
/// contributed a junction to the chosen path, so they cannot run past what reads at this locus support.
pub(crate) fn cothread_locus_geometry(
    transcripts: &[DenovoTranscript],
    members: &[usize],
    min_reads: u32,
) -> Option<(u64, u64, Vec<(u64, u64)>)> {
    let mut weight: BTreeMap<(u64, u64), u64> = BTreeMap::new();
    let mut succ_reads: BTreeMap<(u64, u64), BTreeMap<(u64, u64), u64>> = BTreeMap::new();
    for &i in members {
        let t = &transcripts[i];
        for &j in &t.introns {
            *weight.entry(j).or_insert(0) += t.n_reads as u64;
        }
        for w in t.introns.windows(2) {
            *succ_reads.entry(w[0]).or_default().entry(w[1]).or_insert(0) += t.n_reads as u64;
        }
    }
    // The SAME floor applies to edges as to nodes. Requiring only that SOME read witnesses A -> B is not
    // enough: a readthrough read genuinely witnesses successions ACROSS genes, so co-observation alone let
    // the path walk out of a locus into its neighbour. Measured at NOTCH2NLC (annotated span 81 kb, 3
    // introns): the path ran to 161 kb and 88 introns, NONE of them inside the gene, threading small
    // junctions in the downstream NOTCH2NL/NBPF repeat region. Neither existing guard catches that -- the
    // readthrough filter only targets SINGLE-EXON transcripts, and the mis-chain filter only targets GIANT
    // introns, while these were a 709 bp median. A readthrough succession is carried by few reads, so the
    // node floor applied to edges removes it and leaves real successions alone.
    let succ: BTreeMap<(u64, u64), BTreeSet<(u64, u64)>> = succ_reads
        .into_iter()
        .map(|(a, nxt)| {
            (
                a,
                nxt.into_iter()
                    .filter(|&(_, r)| r >= min_reads as u64)
                    .map(|(b, _)| b)
                    .collect(),
            )
        })
        .collect();
    let nodes: Vec<(u64, u64)> = weight
        .iter()
        .filter(|&(_, &w)| w >= min_reads as u64)
        .map(|(&j, _)| j)
        .collect();
    if nodes.is_empty() {
        return None;
    }
    let index: BTreeMap<(u64, u64), usize> =
        nodes.iter().enumerate().map(|(i, &j)| (j, i)).collect();
    let mut best: Vec<u64> = nodes.iter().map(|j| weight[j]).collect();
    let mut prev: Vec<Option<usize>> = vec![None; nodes.len()];
    for i in 0..nodes.len() {
        let (_, acceptor) = nodes[i];
        let Some(nexts) = succ.get(&nodes[i]) else {
            continue;
        };
        for nxt in nexts {
            let Some(&k) = index.get(nxt) else { continue };
            // Forward only, and the introns must not overlap, or the chain is not a valid transcript.
            if k <= i || nxt.0 < acceptor {
                continue;
            }
            let cand = best[i] + weight[nxt];
            if cand > best[k] {
                best[k] = cand;
                prev[k] = Some(i);
            }
        }
    }
    let mut end_i = 0usize;
    for i in 1..nodes.len() {
        if best[i] > best[end_i] {
            end_i = i;
        }
    }
    let mut chain: Vec<(u64, u64)> = Vec::new();
    let mut cur = Some(end_i);
    while let Some(i) = cur {
        chain.push(nodes[i]);
        cur = prev[i];
    }
    chain.reverse();
    let on_path: BTreeSet<(u64, u64)> = chain.iter().copied().collect();
    // Terminal exons come only from transcripts whose WHOLE chain is a sub-chain of the path. Accepting any
    // transcript that merely SHARES a junction lets a rejected readthrough donate its span: at NOTCH2NLC the
    // model kept a single 80 kb first exon that way, even though the cross-gene succession itself had
    // already been refused. A contributor must be consistent with the model, not merely touch it.
    let (mut start, mut stop) = (u64::MAX, 0u64);
    for &i in members {
        let t = &transcripts[i];
        if !t.introns.is_empty() && t.introns.iter().all(|j| on_path.contains(j)) {
            start = start.min(t.start);
            stop = stop.max(t.end);
        }
    }
    if start == u64::MAX || stop <= start {
        return None;
    }
    // The path's own extent must fit inside the contributing transcripts' bounds.
    if chain[0].0 < start || chain[chain.len() - 1].1 > stop {
        return None;
    }
    Some((start, stop, chain))
}

fn locus_reps(transcripts: &[DenovoTranscript], parent: &[usize]) -> Vec<usize> {
    let n = transcripts.len();
    let mut comp: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    let mut ptmp = parent.to_vec();
    for i in 0..n {
        let r = uf_find(&mut ptmp, i);
        comp.entry(r).or_default().push(i);
    }
    let mut reps: Vec<usize> = comp
        .into_values()
        .map(|members| pick_locus_rep(transcripts, &members))
        .collect();
    reps.sort_unstable();
    reps
}

/// Choose one representative transcript for a collapsed locus.
///
/// DEFAULT: the maximum `(n_reads, span)`. That is wrong in a specific, systematic way inside segmental
/// duplications. `pass1_skeletons_robust` groups spliced reads by EXACT intron chain, while
/// `cluster_unspliced` pools unspliced reads by span overlap with no chain constraint. So at a mis-chaining
/// -prone locus the spliced reads shatter — NOTCH2's 490 spliced reads occupy 240 distinct chains, the
/// largest holding 5% of them — while every unspliced read lands in ONE cluster. Comparing each small
/// spliced chain individually against that pooled cluster means the unspliced side wins by construction:
/// SRGAP2C emits a 1-exon/96-read copy though its biggest spliced chain has 46 reads. The locus is then
/// represented by what is frequently INTRONIC pre-mRNA (§14), and the family fails any isoform test.
///
/// (`RUSTLE_SPLICED_REP`, which summed spliced vs unspliced evidence before choosing, was removed on 2026-09-24:
/// register M row "chr7 F1 0.570 -> 0.411, chr16 0.910 -> 0.761".)
fn pick_locus_rep(transcripts: &[DenovoTranscript], members: &[usize]) -> usize {
    let best = |cands: &[usize]| -> usize {
        *cands
            .iter()
            .max_by(|&&a, &&b| {
                let ka = (
                    transcripts[a].n_reads,
                    transcripts[a].end - transcripts[a].start,
                );
                let kb = (
                    transcripts[b].n_reads,
                    transcripts[b].end - transcripts[b].start,
                );
                // break key ties by preferring the smaller index (so it is the "maximum").
                ka.cmp(&kb).then_with(|| b.cmp(&a))
            })
            .unwrap()
    };
    best(members)
}

/// Exon intervals implied by a transcript's `(start, end, introns)`.
pub(crate) fn exons_of(t: &DenovoTranscript) -> Vec<(u64, u64)> {
    let mut out = Vec::with_capacity(t.introns.len() + 1);
    let mut prev = t.start;
    for &(d, a) in &t.introns {
        if d > prev {
            out.push((prev, d));
        }
        prev = a;
    }
    if t.end > prev {
        out.push((prev, t.end));
    }
    out
}

/// Overlap between two transcripts' EXONS, plus each one's exonic length.
///
/// The span-based alternative counts INTRONIC space as shared, which is how a mis-chained model whose giant
/// intron happens to cross a locus absorbs that locus during collapse. At NPIPB12 the surviving rep spans
/// 152 kb while carrying only 2,540 bp of exon — 98% of its "overlap" with the true 23 kb locus is intron,
/// i.e. sequence the model itself asserts is spliced OUT. Measuring on exons makes containment mean what it
/// is supposed to: the two models describe the same transcribed sequence.
fn exonic_overlap(a: &DenovoTranscript, b: &DenovoTranscript) -> (u64, u64, u64) {
    let (ea, eb) = (exons_of(a), exons_of(b));
    let la: u64 = ea.iter().map(|(s, e)| e - s).sum();
    let lb: u64 = eb.iter().map(|(s, e)| e - s).sum();
    let mut ov = 0u64;
    for &(s1, e1) in &ea {
        for &(s2, e2) in &eb {
            let (lo, hi) = (s1.max(s2), e1.min(e2));
            if hi > lo {
                ov += hi - lo;
            }
        }
    }
    (ov, la, lb)
}

/// Span-aware isoform-to-gene-locus collapse.
///
/// First performs the standard junction-based collapse (`collapse_loci`). Then, if
/// `p.collapse_span_aware` is true, iteratively merges locus representatives that share no junction but
/// are on the same chromosome and strand, have overlapping spans, and satisfy either:
///   - strong containment: the overlap covers at least `COLLAPSE_CONTAIN_FRAC` of the shorter transcript, or
///   - high sequence homology: POA contiguous-core coverage >= `p.collapse_span_core`.
///
/// This recovers genuine alternative-splice isoforms whose intron sets are disjoint (e.g., alternative
/// first/last exons) without chaining adjacent paralogs, which typically fail both the containment and
/// the conservative core-coverage bars.
/// Union-find grouping of transcripts into loci (junction-share + span-aware merge). The `parent` array is the
/// grouping the reps are derived from; exposed so a caller can sum reads over a whole LOCUS, not just its
/// representative isoform (see `collapse_loci_span_aware_with_totals`).
fn collapse_parent(transcripts: &[DenovoTranscript], p: &DetectParams) -> Vec<usize> {
    let n = transcripts.len();
    let mut parent: Vec<usize> = (0..n).collect();
    if n == 0 {
        return parent;
    }

    // Phase 1: junction-based collapse (identical to collapse_loci).
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

    // `RUSTLE_LOCUS_JUNCTION_ONLY=1` stops here, leaving a locus defined ONLY by phase 1: the connected
    // components of "two transcripts share a read-witnessed junction".
    //
    // WHY THIS IS A DIFFERENT KIND OF DEFINITION, not just a stricter one. Phase 1 is a relation on
    // EVIDENCE and needs no representative. Phase 2 below is a fixed-point iteration -- it calls
    // `locus_reps` to pick representatives, merges them by SPAN, then recomputes -- so membership depends
    // on the representative and the representative depends on membership. Two consequences follow from
    // that, and only from that:
    //   - span overlap is satisfied VACUOUSLY by a giant intron, so a transcript that splices straight OVER
    //     a gene is admitted as a member of it (NPIPB9: the selected rep has 24 aligned blocks, NONE inside
    //     the gene, joined by one 104,410 bp intron spanning the whole of it -- membership by absence);
    //   - the merge is order-dependent through the rep it happens to pick at each round.
    // Phase 1 has neither property. It cannot admit a transcript that shares no junction with the locus, so
    // it cannot be satisfied by an intron.
    //
    // The cost is the thing phase 2 was written for: genuine alternative-first/last-exon isoforms whose
    // intron sets are DISJOINT stay separate loci, splitting one gene into several. That trade is what this
    // knob exists to measure. Default unset = today's behaviour, byte-identical.
    let junction_only =
        matches!(std::env::var("RUSTLE_LOCUS_JUNCTION_ONLY"), Ok(v) if v != "0" && !v.is_empty());
    if !p.collapse_span_aware || junction_only {
        return parent;
    }

    // Phase 2: iterative span-aware merge of locus representatives.
    //
    // Under `RUSTLE_COLLAPSE_UNSTRANDED`, an INTRONLESS transcript is treated as having UNKNOWN strand rather than
    // `+`. Strand is derived from canonical junction motifs, so a cluster with no junctions has none to
    // derive from and its `+` is a placeholder. The merge condition `a.strand == b.strand` nevertheless
    // treats that placeholder as evidence, which makes an unspliced cluster unmergeable with any MINUS-strand
    // gene's spliced skeletons. Measured on the Soto catalog: all 740 single-exon copies are `+`, while
    // spliced copies split 149 `+` / 218 `-`. So at a minus-strand gene -- SRGAP2 among them -- the intronic
    // unspliced cluster never even meets the spliced evidence, survives as an independent locus, and is
    // emitted as its own single-exon copy (§14, §20).
    //
    // History: this clause used to be welded to `RUSTLE_SPLICED_REP`, which also changed WHICH transcript
    // represents a locus and was refuted end to end (chr7 F1 0.570 -> 0.411, chr16 0.910 -> 0.761; removed
    // 2026-09-24). `RUSTLE_COLLAPSE_UNSTRANDED` turns on the STRAND clause alone. Threshold-free in the
    // `is_chimeric_bridge` sense: it asks whether the strand field was MEASURED, not how large a number is.
    // Unset, the `else` branch below is reached unchanged, so the OFF arm is byte-identical by construction.
    let unstranded_unspliced =
        matches!(std::env::var("RUSTLE_COLLAPSE_UNSTRANDED"), Ok(v) if v != "0" && !v.is_empty());
    let exonic_collapse =
        matches!(std::env::var("RUSTLE_COLLAPSE_EXONIC"), Ok(v) if v != "0" && !v.is_empty());
    // SPEED (2026-09-24, byte-identical): the merge decisions of one round depend only on that round's
    // representatives, never on unions made earlier in the same round, so the round's union SET — and hence
    // the partition `locus_reps` reads — is independent of pair order. Two consequences are exploited:
    //   - only SPAN-OVERLAPPING same-contig pairs are visited (a sweep over reps sorted by (chrom, start)):
    //     a pair with no span overlap has `ov == 0` under both the span and the exonic measure, so it can
    //     never merge; the old all-pairs loop visited n^2/2 pairs per round (chr16: ~25k reps);
    //   - the POA core-coverage calls (the expensive step: 257 of 292 s on a gorilla contig, ~3,300 of
    //     3,495 s on human chr16) run in parallel via rayon, exactly as `confirm_edge` already does, and
    //     their verdicts are applied afterwards. Each pair keeps its original (lower rep index, higher)
    //     orientation, so every call sees the same arguments as before.
    use rayon::prelude::*;
    let t_collapse = std::time::Instant::now();
    let (mut n_rounds, mut n_visited, mut n_poa) = (0usize, 0usize, 0usize);
    let mut slowest: (f64, usize, usize) = (0.0, 0, 0); // (seconds, len a, len b) of the slowest POA call
    loop {
        n_rounds += 1;
        let reps = locus_reps(transcripts, &parent);
        let mut order: Vec<usize> = (0..reps.len()).collect();
        order.sort_by(|&x, &y| {
            let (a, b) = (&transcripts[reps[x]], &transcripts[reps[y]]);
            (a.chrom.as_str(), a.start, x).cmp(&(b.chrom.as_str(), b.start, y))
        });
        let mut poa_pairs: Vec<(usize, usize)> = Vec::new();
        let mut merges: Vec<(usize, usize)> = Vec::new();
        for (oi, &x) in order.iter().enumerate() {
            let (sx_chrom, sx_end) = (&transcripts[reps[x]].chrom, transcripts[reps[x]].end);
            for &y in &order[oi + 1..] {
                let ty = &transcripts[reps[y]];
                if &ty.chrom != sx_chrom || ty.start >= sx_end {
                    break; // sorted by (chrom, start): no later rep overlaps x's span
                }
                let (i, j) = if x < y { (x, y) } else { (y, x) };
                n_visited += 1;
                let a = &transcripts[reps[i]];
                let b = &transcripts[reps[j]];
                // Strand blocks a merge only when BOTH sides actually have a junction-derived strand.
                let strand_conflict = if unstranded_unspliced {
                    !a.introns.is_empty() && !b.introns.is_empty() && a.strand != b.strand
                } else {
                    a.strand != b.strand
                };
                if a.chrom != b.chrom || strand_conflict {
                    continue;
                }
                // `RUSTLE_COLLAPSE_EXONIC`: measure containment on EXONS rather than the genomic span, so a
                // model whose giant intron crosses a locus cannot absorb it. This must accompany
                // RUSTLE_JUNCTION_MAJORITY: relaxing the junction gate alone admits longer models that then
                // BRIDGE separate loci through exactly this span-based rule (chr16 66 -> 34 copies).
                // ⚠ RUSTLE_JUNCTION_MAJORITY's default flipped to ON 2026-09-21 WITHOUT also flipping this
                // flag: the validated full chr16 arm (`bench/CHR16_JUNCTION_MAJORITY_ARM.md`, §6n4) ran
                // `gw_family_catalog` with junction-majority alone (this flag at its own default, off,
                // in both arms) and the feared fusion did not occur there (families 282->290, max size
                // unchanged at 71). That test used `gw_family_catalog`, not this exact code path's every
                // caller, so the concern this comment states remains logically live beyond chr16 --
                // re-test before trusting this pairing is unnecessary on a new substrate.
                let (ov, minlen) = if exonic_collapse {
                    let (eov, la, lb) = exonic_overlap(a, b);
                    (eov, la.min(lb).max(1))
                } else {
                    let ov = a.end.min(b.end).saturating_sub(a.start.max(b.start));
                    (ov, (a.end - a.start).min(b.end - b.start).max(1))
                };
                if ov == 0 {
                    continue;
                }
                let containment = ov as f64 / minlen as f64;
                if containment >= COLLAPSE_CONTAIN_FRAC {
                    merges.push((i, j));
                    continue;
                }
                // Disjoint-junction isoforms with similar length: require strong POA core coverage
                // (evaluated below, in parallel).
                poa_pairs.push((i, j));
            }
        }
        n_poa += poa_pairs.len();
        let poa_merge: Vec<(bool, f64)> = poa_pairs
            .par_iter()
            .map(|&(i, j)| {
                let t = std::time::Instant::now();
                let au = upper_cow(&transcripts[reps[i]].seq);
                let bu = upper_cow(&transcripts[reps[j]].seq);
                // exact substring bound first (skips the hopeless alignments), then the POA core
                let m = core_coverage_reaches(&au, &bu, p.len_cap, p.collapse_span_core);
                (m, t.elapsed().as_secs_f64())
            })
            .collect();
        for (&(i, j), &(_, secs)) in poa_pairs.iter().zip(&poa_merge) {
            if secs > slowest.0 {
                slowest = (
                    secs,
                    transcripts[reps[i]].seq.len(),
                    transcripts[reps[j]].seq.len(),
                );
            }
        }
        let poa_secs: f64 = poa_merge.iter().map(|x| x.1).sum();
        if collapse_stats_enabled() {
            eprintln!(
                "[collapse] round {n_rounds}: {} reps, {} POA pairs, {poa_secs:.1} CPU-s in POA",
                reps.len(),
                poa_pairs.len()
            );
        }
        for (&(i, j), &(m, _)) in poa_pairs.iter().zip(&poa_merge) {
            if m {
                merges.push((i, j));
            }
        }
        // Replay the round's unions in (i, j) lexicographic order — exactly the order the all-pairs loop
        // applied them — so `parent` (and every component ROOT, which `locus_groups` orders its output by
        // under RUSTLE_LOCUS_EXON_UNION / RUSTLE_COTHREAD_REP) is identical, not just the partition.
        merges.sort_unstable();
        for &(i, j) in &merges {
            uf_union(&mut parent, reps[i], reps[j]);
        }
        let merged = !merges.is_empty();
        if !merged {
            break;
        }
    }
    if collapse_stats_enabled() || t_collapse.elapsed().as_secs() >= 30 {
        eprintln!(
            "[collapse] {} transcripts: {n_rounds} round(s), {n_visited} span-overlapping pairs, {n_poa} POA calls, {:.1} s wall; slowest POA {:.1} s ({} x {} bp)",
            transcripts.len(), t_collapse.elapsed().as_secs_f64(), slowest.0, slowest.1, slowest.2
        );
    }

    parent
}

/// `RUSTLE_COLLAPSE_STATS=1`: per-round POA counts of the locus collapse on stderr (the summary line is
/// printed anyway whenever the collapse takes >= 30 s).
fn collapse_stats_enabled() -> bool {
    matches!(std::env::var("RUSTLE_COLLAPSE_STATS"), Ok(v) if v != "0" && !v.is_empty())
}

/// One representative transcript index per collapse locus. Behavior-preserving thin wrapper over
/// `collapse_parent` + `locus_reps` (unchanged public signature).
pub fn collapse_loci_span_aware(transcripts: &[DenovoTranscript], p: &DetectParams) -> Vec<usize> {
    if transcripts.is_empty() {
        return Vec::new();
    }
    let parent = collapse_parent(transcripts, p);
    locus_reps(transcripts, &parent)
}

/// Like `collapse_loci_span_aware`, but also returns, PARALLEL to the reps, the MEMBER INDICES of each
/// rep's locus -- every isoform that collapsed into it, including the rep itself.
///
/// Needed because the pipeline otherwise discards the non-representative isoforms, and they carry exon
/// sequence the representative may lack: 46% of the representatives covering a known family member are
/// single-exon stubs, yet 53% of those loci have a spliced chain supported by >= 3 reads. Uses the SAME
/// `collapse_parent` union-find as `collapse_loci_span_aware`, so reps and members cannot disagree.
pub fn collapse_loci_span_aware_with_members(
    transcripts: &[DenovoTranscript],
    p: &DetectParams,
) -> (Vec<usize>, Vec<Vec<usize>>) {
    if transcripts.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let parent = collapse_parent(transcripts, p);
    let reps = locus_reps(transcripts, &parent);
    let mut ptmp = parent.clone();
    let mut by_root: DetHashMap<usize, Vec<usize>> = DetHashMap::default();
    for i in 0..transcripts.len() {
        let r = uf_find(&mut ptmp, i);
        by_root.entry(r).or_default().push(i);
    }
    let members: Vec<Vec<usize>> = reps
        .iter()
        .map(|&rep| {
            let r = uf_find(&mut ptmp, rep);
            by_root.get(&r).cloned().unwrap_or_else(|| vec![rep])
        })
        .collect();
    (reps, members)
}

/// Like `collapse_loci_span_aware`, but also returns, PARALLEL to the reps, the LOCUS TOTAL read count = the
/// summed `n_reads` over ALL isoforms collapsed into each rep's locus (not just the representative isoform).
/// This is the correct single-copy expression basis for λ_global: a gene's total expression, not its dominant
/// isoform's count. `totals[k]` corresponds to `reps[k]`.
pub fn collapse_loci_span_aware_with_totals(
    transcripts: &[DenovoTranscript],
    p: &DetectParams,
) -> (Vec<usize>, Vec<u32>) {
    if transcripts.is_empty() {
        return (Vec::new(), Vec::new());
    }
    let parent = collapse_parent(transcripts, p);
    let reps = locus_reps(transcripts, &parent);
    let mut ptmp = parent.clone();
    let mut sum_by_root: DetHashMap<usize, u32> = DetHashMap::default();
    for (i, t) in transcripts.iter().enumerate() {
        let r = uf_find(&mut ptmp, i);
        *sum_by_root.entry(r).or_insert(0) += t.n_reads;
    }
    let totals: Vec<u32> = reps
        .iter()
        .map(|&rep| {
            let r = uf_find(&mut ptmp, rep);
            sum_by_root[&r]
        })
        .collect();
    (reps, totals)
}

/// (2) Candidate homologous rep pairs. Exact canonical-k-mer ownership pre-filter (family-informative =
/// owned by `[cnt_min, cnt_max]` reps; a rep needs `>= k_share` informative k-mers) + contiguous-span
/// filter (shared informative-k-mer position span `>= t_core * min(len)`). Returns `(i, j)` index pairs
/// into `reps` with `i < j`, sorted.
pub fn candidate_pairs(reps: &[DenovoTranscript], p: &DetectParams) -> Vec<(usize, usize)> {
    let n = reps.len();
    // per-rep canonical k-mers with first-occurrence positions.
    let rep_kmers: Vec<BTreeMap<u64, u32>> = reps
        .iter()
        .map(|r| canonical_kmer_first_pos(&r.seq))
        .collect();

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

/// (3) POA edge confirmation: contiguous-core coverage in both orientations (reverse-complement fallback for
/// opposite-strand copies), `Some(core_recip)` iff `>= t_core`. Uppercases the operands (the soft-mask lesson
/// from `family_rescue`: `reverse_complement` maps lowercase → N).
///
/// `len_cap` is the poasta MEMORY THRESHOLD, not a hard skip: when the LARGER operand exceeds it, poasta's
/// graph aligner would OOM (it allocates with the longer sequence — a (5 kb, 228 kb) pair blows up on the
/// 228 kb side even though `min = 5 kb`), so confirmation falls back to a linear-memory longest-common-
/// substring metric (`contiguous_core_coverage_bounded`). The fallback is faithful — a poasta ungapped-equal
/// run IS a common substring — so a large read-through "hub" that homologously contains a copy is still
/// confirmed (the DSFAM45 case) instead of being lost. Below the cap the exact poasta path is unchanged.
pub fn confirm_edge(a: &[u8], b: &[u8], p: &DetectParams) -> Option<f64> {
    use family_graph::{
        contiguous_core_coverage_bounded_with, longest_common_substring, EDGE_CONFIRM_ASTAR,
    };
    let au = upper_cow(a);
    let bu = upper_cow(b);
    if p.edge_core == EdgeCore::Lcs {
        let minlen = au.len().min(bu.len());
        if minlen == 0 {
            return None;
        }
        let mut cr = longest_common_substring(&au, &bu) as f64 / minlen as f64;
        if cr < p.t_core {
            let rc = longest_common_substring(&au, &reverse_complement(&bu)) as f64 / minlen as f64;
            if rc > cr {
                cr = rc;
            }
        }
        return (cr >= p.t_core).then_some(cr);
    }
    let mut cr = contiguous_core_coverage_bounded_with(&au, &bu, p.len_cap, EDGE_CONFIRM_ASTAR);
    if cr < p.t_core {
        // opposite orientation (copies assembled on different strands).
        let rc = contiguous_core_coverage_bounded_with(
            &au,
            &reverse_complement(&bu),
            p.len_cap,
            EDGE_CONFIRM_ASTAR,
        );
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
    detect_edges_reporting(reps, p).0
}

/// POA-CORE COMPLETION — extend read-conflict families with loosely-related paralogs at NEW loci.
///
/// The read-conflict graph links only copies that reads CONFUSE (empirically down to ~87% identity); a
/// divergent paralog at another locus (reads resolve it, so it raises no conflict edge) is missed. Seeded by
/// the conflict families ("when assignment is determined needed"), this attaches any genome rep that shares a
/// contiguous POA core `>= p.t_core` with a family member — reaching paralogs that retain ONE conserved
/// exon-core even as their flanks diverge (the loosely-related case the conflict graph cannot reach).
///
/// BOUNDED ("only when needed"): `candidate_pairs` (the minimizer-LSH prefilter) restricts the expensive POA
/// `confirm_edge` to homologous candidate pairs, and we grade ONLY pairs with exactly one endpoint already in
/// a conflict family (the other a FREE rep) — so POA never runs all-pairs, and free-vs-free pairs (which would
/// be NEW families, not seeded by read-conflict) are skipped. A free rep matching several families is attached
/// to the one with the strongest core.
///
/// `families[f]` = the rep indices of conflict family `f`. Returns `adds[f]` = the extra rep indices to append
/// to family `f` (the homology-extension copies, distinct from the confusable core).
pub fn poa_core_completion_adds(
    reps: &[DenovoTranscript],
    families: &[Vec<usize>],
    p: &DetectParams,
) -> Vec<Vec<usize>> {
    use DetHashMap;
    let mut in_fam: DetHashMap<usize, usize> = DetHashMap::default();
    for (f, members) in families.iter().enumerate() {
        for &m in members {
            in_fam.insert(m, f);
        }
    }
    // FAST contiguous-core confirm: the longest-common-substring / min-len — the LINEAR-MEMORY faithful
    // equivalent of `contiguous_core_coverage`'s "longest ungapped equal run" (see `longest_common_substring`
    // docs). Using poasta (`confirm_edge`) here is too slow genome-wide (the documented POA bottleneck); LCS
    // is O(n) per pair, so the family-adjacent grading stays cheap. Checks both orientations.
    let core_cov = |a: &[u8], b: &[u8]| -> f64 {
        let minlen = a.len().min(b.len());
        if minlen == 0 {
            return 0.0;
        }
        let au = a.to_ascii_uppercase();
        let bu = b.to_ascii_uppercase();
        let fwd = crate::family::family_detect::family_graph::longest_common_substring(&au, &bu);
        let rev = crate::family::family_detect::family_graph::longest_common_substring(
            &au,
            &reverse_complement(&bu),
        );
        fwd.max(rev) as f64 / minlen as f64
    };
    // best (family, core) per FREE rep over all family-adjacent candidate pairs that confirm a POA core.
    let mut best: DetHashMap<usize, (usize, f64)> = DetHashMap::default();
    for (i, j) in candidate_pairs(reps, p) {
        let (fam, free) = match (in_fam.get(&i), in_fam.get(&j)) {
            (Some(&f), None) => (f, j),
            (None, Some(&f)) => (f, i),
            _ => continue, // both in families (no merge here) or both free (not seeded by read-conflict)
        };
        let member = if free == j { i } else { j };
        let core = core_cov(&reps[member].seq, &reps[free].seq);
        if core >= p.t_core {
            let e = best.entry(free).or_insert((fam, 0.0));
            if core > e.1 {
                *e = (fam, core);
            }
        }
    }
    let mut adds: Vec<Vec<usize>> = vec![Vec::new(); families.len()];
    for (free, (fam, _)) in best {
        adds[fam].push(free);
    }
    for a in &mut adds {
        a.sort_unstable();
    }
    adds
}

/// Like [`detect_edges`] but ALSO returns the candidate pairs whose larger transcript exceeded the poasta
/// memory threshold (`p.len_cap`) and were therefore confirmed via the memory-bounded longest-common-substring
/// FALLBACK rather than the exact poasta path. These edges are still produced (the fallback confirms them — no
/// OOM, no loss); reporting them lets the caller AUDIT which family edges rest on the approximate large-
/// sequence metric. `partition`/candidate order is preserved, so both lists are deterministic.
pub fn detect_edges_reporting(
    reps: &[DenovoTranscript],
    p: &DetectParams,
) -> (Vec<(usize, usize, f64)>, Vec<(usize, usize)>) {
    use rayon::prelude::*;
    let pairs = candidate_pairs(reps, p);
    let fallback: Vec<(usize, usize)> = pairs
        .iter()
        .copied()
        .filter(|&(a, b)| reps[a].seq.len().max(reps[b].seq.len()) > p.len_cap)
        .collect();
    let edges = pairs
        .into_par_iter()
        .filter_map(|(a, b)| confirm_edge(&reps[a].seq, &reps[b].seq, p).map(|cr| (a, b, cr)))
        .collect();
    (edges, fallback)
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
            strand: '+',
            introns: introns.to_vec(),
            seq,
            ..Default::default()
        }
    }
    #[test]
    fn cothread_chain_spans_more_junctions_than_any_single_transcript() {
        // The shattering case: no observed chain carries the whole locus, but consecutive junctions are
        // co-observed pairwise, so the graph path recovers all four.
        let t = vec![
            tx(
                "a",
                "c",
                100,
                900,
                5,
                &[(200, 300), (400, 500)],
                vec![b'A'; 10],
            ),
            tx(
                "b",
                "c",
                100,
                900,
                4,
                &[(400, 500), (600, 700)],
                vec![b'A'; 10],
            ),
            tx(
                "c",
                "c",
                100,
                900,
                4,
                &[(600, 700), (750, 800)],
                vec![b'A'; 10],
            ),
        ];
        let (s, e, chain) = cothread_locus_geometry(&t, &[0, 1, 2], 3).unwrap();
        assert_eq!((s, e), (100, 900));
        assert_eq!(
            chain,
            vec![(200, 300), (400, 500), (600, 700), (750, 800)],
            "path threads all four junctions though no transcript carries more than two"
        );
    }

    #[test]
    fn cothread_edges_obey_the_same_read_floor_as_nodes() {
        // A readthrough carrying 2 reads witnesses the succession into the neighbouring gene. Both
        // junctions clear the NODE floor on their own support, so only an EDGE floor can stop the walk.
        let t = vec![
            tx("here", "c", 100, 900, 40, &[(200, 300)], vec![b'A'; 10]),
            tx("nbr", "c", 5000, 6000, 40, &[(5200, 5300)], vec![b'A'; 10]),
            tx(
                "readthrough",
                "c",
                100,
                6000,
                2,
                &[(200, 300), (5200, 5300)],
                vec![b'A'; 10],
            ),
        ];
        let (s, e, chain) = cothread_locus_geometry(&t, &[0, 1, 2], 3).unwrap();
        assert_eq!(
            chain.len(),
            1,
            "a 2-read readthrough must not license the cross-gene step"
        );
        assert!(
            e - s < 5000,
            "extent must not run into the neighbour, got {s}-{e}"
        );
        // With enough reads behind it the succession is real and IS taken.
        let t2 = vec![
            tx("here", "c", 100, 900, 40, &[(200, 300)], vec![b'A'; 10]),
            tx(
                "both",
                "c",
                100,
                6000,
                9,
                &[(200, 300), (5200, 5300)],
                vec![b'A'; 10],
            ),
        ];
        assert_eq!(cothread_locus_geometry(&t2, &[0, 1], 3).unwrap().2.len(), 2);
    }

    #[test]
    fn cothread_refuses_to_chain_junctions_no_read_ever_co_observed() {
        // THE POINT OF THE CONSTRAINT. (200,300) and (5000,5100) do not overlap, so a max-weight
        // NON-OVERLAPPING set would happily take both and produce a locus spanning the neighbour. No
        // transcript carries them together, so the co-observed path must not.
        let t = vec![
            tx("here", "c", 100, 400, 9, &[(200, 300)], vec![b'A'; 10]),
            tx("nbr", "c", 4900, 5200, 9, &[(5000, 5100)], vec![b'A'; 10]),
        ];
        let (s, e, chain) = cothread_locus_geometry(&t, &[0, 1], 3).unwrap();
        assert_eq!(
            chain.len(),
            1,
            "no read witnesses the succession, so no chaining across genes"
        );
        assert!(
            e - s <= 300,
            "extent stays at the contributing transcript, got {s}-{e}"
        );
    }

    #[test]
    fn cothread_applies_the_junction_floor_and_abstains_when_nothing_clears_it() {
        let t = vec![
            tx(
                "weak",
                "c",
                100,
                900,
                2,
                &[(200, 300), (400, 500)],
                vec![b'A'; 10],
            ),
            tx("stub", "c", 100, 900, 50, &[], vec![b'A'; 10]),
        ];
        assert!(
            cothread_locus_geometry(&t, &[0, 1], 3).is_none(),
            "2 reads/junction is below the floor; caller must fall back to today's representative"
        );
        // The same junctions, now carried by enough reads, do produce a chain.
        let t2 = vec![tx(
            "ok",
            "c",
            100,
            900,
            3,
            &[(200, 300), (400, 500)],
            vec![b'A'; 10],
        )];
        assert_eq!(
            cothread_locus_geometry(&t2, &[0], 3).unwrap().2,
            vec![(200, 300), (400, 500)]
        );
    }

    #[test]
    fn cothread_prefers_the_heavier_path_when_two_compete() {
        // Both successions are read-witnessed from (200,300); the heavier one wins.
        let t = vec![
            tx(
                "light",
                "c",
                100,
                900,
                3,
                &[(200, 300), (400, 500)],
                vec![b'A'; 10],
            ),
            tx(
                "heavy",
                "c",
                100,
                900,
                30,
                &[(200, 300), (600, 700)],
                vec![b'A'; 10],
            ),
        ];
        let chain = cothread_locus_geometry(&t, &[0, 1], 3).unwrap().2;
        assert_eq!(chain, vec![(200, 300), (600, 700)]);
    }

    #[test]
    fn union_geometry_merges_partial_chains_into_the_whole_locus() {
        // Two chains covering opposite halves of one locus -- the shattering case. Either alone is a
        // fragment; the union is the locus.
        let t = vec![
            tx("t100_5", "c", 100, 400, 5, &[(200, 300)], vec![b'A'; 10]), // exons 100-200, 300-400
            tx("t500_4", "c", 500, 900, 4, &[(600, 700)], vec![b'A'; 10]), // exons 500-600, 700-900
        ];
        let (s, e, introns) = union_locus_geometry(&t, &[0, 1], 3).unwrap();
        assert_eq!((s, e), (100, 900));
        // gaps between merged exons: 200-300, 400-500, 600-700
        assert_eq!(introns, vec![(200, 300), (400, 500), (600, 700)]);
    }

    #[test]
    fn union_geometry_drops_chains_below_the_read_floor() {
        // The floor is the mechanism, not a knob: without it the union readmits noise and rebuilds the
        // 40-family blob. A 2-read chain must contribute nothing at floor 3.
        let t = vec![
            tx("t100_5", "c", 100, 200, 5, &[], vec![b'A'; 10]),
            tx("t900_2", "c", 900, 1000, 2, &[], vec![b'A'; 10]),
        ];
        let (s, e, introns) = union_locus_geometry(&t, &[0, 1], 3).unwrap();
        assert_eq!(
            (s, e),
            (100, 200),
            "the 2-read chain must not extend the locus"
        );
        assert!(introns.is_empty());
    }

    #[test]
    fn union_geometry_returns_none_when_no_chain_clears_the_floor() {
        // Whole locus is lost rather than silently rebuilt from sub-threshold evidence. This is the
        // measured cost of the substrate (Soto 299 -> 263 loci).
        let t = vec![
            tx("t100_1", "c", 100, 200, 1, &[], vec![b'A'; 10]),
            tx("t150_2", "c", 150, 250, 2, &[], vec![b'A'; 10]),
        ];
        assert!(union_locus_geometry(&t, &[0, 1], 3).is_none());
    }

    #[test]
    fn union_geometry_of_a_single_chain_reproduces_that_chain() {
        // With one chain above the floor the union must be a no-op, so enabling the substrate cannot
        // perturb loci that never shattered.
        let t = vec![tx(
            "t100_9",
            "c",
            100,
            400,
            9,
            &[(200, 300)],
            vec![b'A'; 10],
        )];
        let (s, e, introns) = union_locus_geometry(&t, &[0], 3).unwrap();
        assert_eq!((s, e), (t[0].start, t[0].end));
        assert_eq!(introns, t[0].introns);
    }

    #[test]
    fn union_geometry_absorbs_a_retained_intron_isoform() {
        // One chain splices 200-300, another retains it. The union treats retained sequence as exonic,
        // which is the point: include everything the reads support.
        let t = vec![
            tx("t100_6", "c", 100, 400, 6, &[(200, 300)], vec![b'A'; 10]),
            tx("t100_4", "c", 100, 400, 4, &[], vec![b'A'; 10]),
        ];
        let (s, e, introns) = union_locus_geometry(&t, &[0, 1], 3).unwrap();
        assert_eq!((s, e), (100, 400));
        assert!(introns.is_empty(), "retained-intron support fills the gap");
    }

    #[test]
    fn exonic_containment_refuses_to_merge_a_locus_swallowed_by_a_giant_intron() {
        // The NPIPB12 shape: a 10 kb model whose single intron spans a separate 1 kb locus. By SPAN the small
        // one is 100% contained; by EXONS they share nothing, because the "overlap" is sequence the long
        // model asserts is spliced out.
        let long = DenovoTranscript {
            chrom: "c1".into(),
            start: 0,
            end: 10_000,
            strand: '+',
            introns: vec![(500, 9_500)],
            seq: b"AC".to_vec(),
            ..Default::default()
        };
        let inner = DenovoTranscript {
            chrom: "c1".into(),
            start: 4_000,
            end: 5_000,
            strand: '+',
            introns: vec![],
            seq: b"AC".to_vec(),
            ..Default::default()
        };
        // span containment: the inner locus is fully inside the long model's span
        let span_ov = long
            .end
            .min(inner.end)
            .saturating_sub(long.start.max(inner.start));
        assert_eq!(span_ov, 1_000);
        let span_minlen = (long.end - long.start).min(inner.end - inner.start);
        assert_eq!(
            span_ov as f64 / span_minlen as f64,
            1.0,
            "span rule says fully contained -> would merge"
        );

        // exonic containment: zero shared exonic sequence
        let (eov, la, lb) = exonic_overlap(&long, &inner);
        assert_eq!(eov, 0, "no shared EXONIC sequence");
        assert_eq!(la, 1_000, "long model's exons are 0-500 and 9500-10000");
        assert_eq!(lb, 1_000);

        // and it still merges genuine isoforms that share exonic sequence
        let iso = DenovoTranscript {
            chrom: "c1".into(),
            start: 0,
            end: 400,
            strand: '+',
            introns: vec![],
            seq: b"AC".to_vec(),
            ..Default::default()
        };
        let (eov2, _, lb2) = exonic_overlap(&long, &iso);
        assert_eq!(eov2, 400);
        assert!(
            eov2 as f64 / lb2 as f64 >= COLLAPSE_CONTAIN_FRAC,
            "real isoform still merges"
        );
    }

    #[test]
    fn collapse_loci_groups_maps_isoforms_to_their_gene_rep() {
        // two isoforms of gene A share the junction (100,200); a third transcript at a disjoint locus is its
        // own gene. groups[i] must be the rep index of i's gene (rep = most reads, here A_iso1 with 10).
        let txs = vec![
            tx("A_iso1", "c1", 0, 300, 10, &[(100, 200)], vec![b'A'; 50]),
            tx(
                "A_iso2",
                "c1",
                0,
                400,
                4,
                &[(100, 200), (250, 320)],
                vec![b'A'; 60],
            ),
            tx(
                "B_iso1",
                "c1",
                1000,
                1300,
                8,
                &[(1100, 1200)],
                vec![b'C'; 50],
            ),
        ];
        let g = collapse_loci_groups(&txs);
        assert_eq!(
            g[0], g[1],
            "isoforms sharing a junction collapse to one gene"
        );
        assert_ne!(g[0], g[2], "a disjoint locus is its own gene");
        assert_eq!(g[0], 0, "gene rep = the higher-read isoform (A_iso1)");
        assert_eq!(g[2], 2, "B is its own rep");
    }
    fn homolog_tx_flank(
        tid: &str,
        fs1: u64,
        core: &[u8],
        fs2: u64,
        flank: usize,
        n: u32,
    ) -> DenovoTranscript {
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
        assert_eq!(
            m.get(&first_clean),
            Some(&0),
            "first surviving window -> compacted position 0"
        );
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
            tx(
                "B",
                "c1",
                100,
                900,
                7,
                &[(200, 300), (600, 700)],
                rand_seq(400, 2),
            ),
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
        assert_eq!(
            collapse_loci(&txs).len(),
            2,
            "same coords, different chrom -> not merged"
        );
    }

    #[test]
    fn collapse_rep_tiebreak_longest_span() {
        let txs = [
            tx("A", "c1", 100, 500, 5, &[(200, 300)], rand_seq(400, 1)),
            tx("B", "c1", 100, 900, 5, &[(200, 300)], rand_seq(800, 2)),
        ];
        assert_eq!(
            collapse_loci(&txs),
            vec![1],
            "equal reads -> longest span wins"
        );
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

    // ---- collapse_loci_span_aware (disjoint-junction isoform recovery) ----

    #[test]
    fn span_aware_collapse_merges_disjoint_junction_isoforms() {
        // Two isoforms of the same gene with NO shared junction, but span-overlapping and sharing a
        // 300 bp core (75% contiguous-core coverage). The span-aware pass should collapse them.
        let core = rand_seq(300, 0xBEEF);
        let iso1_seq = cat(&[&rand_seq(50, 0xA1), &core, &rand_seq(50, 0xA2)]);
        let iso2_seq = cat(&[&rand_seq(50, 0xB1), &core, &rand_seq(50, 0xB2)]);
        let txs = [
            tx("iso1", "c1", 100, 600, 5, &[(200, 300)], iso1_seq),
            tx("iso2", "c1", 250, 750, 5, &[(450, 550)], iso2_seq),
        ];
        let p = DetectParams::default();
        assert_eq!(
            collapse_loci_span_aware(&txs, &p).len(),
            1,
            "disjoint-junction isoforms with high span homology collapse to one locus"
        );
        // With span-aware disabled, they stay separate.
        let off = DetectParams {
            collapse_span_aware: false,
            ..DetectParams::default()
        };
        assert_eq!(collapse_loci_span_aware(&txs, &off).len(), 2);
    }

    #[test]
    fn span_aware_collapse_merges_by_containment() {
        // A retained-intron / shorter isoform contained within a longer one, with disjoint junctions.
        // Containment alone (no POA needed) should merge them.
        let long_seq = rand_seq(500, 0xCAFE);
        let short_seq = long_seq[50..350].to_vec(); // 300 bp contained within long
        let txs = [
            tx(
                "long",
                "c1",
                100,
                800,
                8,
                &[(200, 300), (500, 600)],
                long_seq,
            ),
            tx("short", "c1", 200, 500, 5, &[(250, 350)], short_seq), // span 200..500 is 75% inside long
        ];
        assert_eq!(
            collapse_loci_span_aware(&txs, &DetectParams::default()).len(),
            1,
            "contained isoform merges into the longer one"
        );
    }

    #[test]
    fn span_aware_collapse_preserves_adjacent_paralogs() {
        // Two adjacent paralog copies with disjoint junctions and low sequence homology.
        // The conservative POA core threshold should NOT merge them.
        let txs = [
            tx(
                "para1",
                "c1",
                1000,
                11000,
                10,
                &[(2000, 3000), (4000, 5000)],
                rand_seq(1000, 0xC1),
            ),
            tx(
                "para2",
                "c1",
                7000,
                17000,
                9,
                &[(8000, 9000), (10000, 11000)],
                rand_seq(1000, 0xC2),
            ),
        ];
        assert_eq!(
            collapse_loci_span_aware(&txs, &DetectParams::default()).len(),
            2,
            "adjacent paralogs with low homology stay as two loci"
        );
    }

    // ---- RUSTLE_COLLAPSE_UNSTRANDED (the placeholder-strand un-weld) ----

    /// Serializes the two tests below against each other. `collapse_parent` reads its switches from the
    /// PROCESS environment, which cargo's harness shares across test threads. Every other collapse test in
    /// this module builds strand `'+'` on both sides, where the two branches of `strand_conflict` are
    /// equivalent, so none of them can be perturbed by what these two set.
    static COLLAPSE_ENV_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

    /// Save-and-restore for the three variables `collapse_parent` reads, so an ambient value in the
    /// caller's shell cannot decide the result and neither test leaks into the other.
    fn with_collapse_env<T>(set: &[(&str, &str)], f: impl FnOnce() -> T) -> T {
        // Held across the whole set/run/restore window. (The lock was declared but never taken, so the two
        // tests below raced each other under the multi-threaded harness; fixed 2026-09-24.)
        let _g = COLLAPSE_ENV_LOCK.lock().unwrap_or_else(|p| p.into_inner());
        const VARS: [&str; 2] = ["RUSTLE_COLLAPSE_UNSTRANDED", "RUSTLE_COLLAPSE_EXONIC"];
        let prior: Vec<(&str, Option<String>)> =
            VARS.iter().map(|k| (*k, std::env::var(k).ok())).collect();
        for k in VARS {
            std::env::remove_var(k);
        }
        for (k, v) in set {
            std::env::set_var(k, v);
        }
        let out = f();
        for (k, v) in prior {
            match v {
                Some(v) => std::env::set_var(k, v),
                None => std::env::remove_var(k),
            }
        }
        out
    }

    /// Build a transcript with an EXPLICIT strand. `tx` above hard-codes `'+'`, which is exactly the
    /// placeholder this clause is about, so these tests cannot use it.
    fn tx_strand(
        tid: &str,
        chrom: &str,
        start: u64,
        end: u64,
        n_reads: u32,
        strand: char,
        introns: &[(u64, u64)],
        seq: Vec<u8>,
    ) -> DenovoTranscript {
        DenovoTranscript {
            strand,
            ..tx(tid, chrom, start, end, n_reads, introns, seq)
        }
    }

    #[test]
    fn collapse_unstranded_clause_merges_intronless_stub_into_minus_strand_gene() {
        // THE CLAUSE. An INTRONLESS model has no canonical junction motif to derive a strand from, so its
        // `'+'` is a placeholder, not a measurement. Strictly enclosed inside a spliced `'-'` gene
        // (containment 1.0, same chrom), the ONLY thing separating the two is that placeholder.
        let gene = tx_strand(
            "gene",
            "c1",
            1_000,
            9_000,
            8,
            '-',
            &[(2_000, 3_000), (5_000, 6_000)],
            rand_seq(900, 0xD1),
        );
        let stub = tx_strand("stub", "c1", 4_000, 4_600, 4, '+', &[], rand_seq(600, 0xD2));
        let txs = [gene, stub];
        let p = DetectParams::default();

        let off = with_collapse_env(&[], || collapse_loci_span_aware(&txs, &p).len());
        assert_eq!(
            off, 2,
            "default: the placeholder '+' is treated as evidence and blocks the merge"
        );

        let on = with_collapse_env(&[("RUSTLE_COLLAPSE_UNSTRANDED", "1")], || {
            collapse_loci_span_aware(&txs, &p).len()
        });
        assert_eq!(
            on, 1,
            "RUSTLE_COLLAPSE_UNSTRANDED: an UNMEASURED strand cannot block the merge"
        );
    }

    #[test]
    fn collapse_unstranded_keeps_genuine_antisense_spliced_pairs_apart_in_both_arms() {
        // THE VALUE, pinned separately from the clause. When BOTH sides carry a junction-derived strand
        // the field WAS measured, so opposite strands are real antisense and must stay two loci whether or
        // not the variable is set. This is the mirror of the genuine-antisense cases the un-weld must not
        // touch.
        let minus = tx_strand(
            "minus",
            "c1",
            1_000,
            9_000,
            8,
            '-',
            &[(2_000, 3_000), (5_000, 6_000)],
            rand_seq(900, 0xE1),
        );
        let plus = tx_strand(
            "plus",
            "c1",
            3_500,
            4_500,
            4,
            '+',
            &[(3_800, 4_000)],
            rand_seq(700, 0xE2),
        );
        let anti = [minus.clone(), plus.clone()];
        let p = DetectParams::default();

        let off = with_collapse_env(&[], || collapse_loci_span_aware(&anti, &p).len());
        let on = with_collapse_env(&[("RUSTLE_COLLAPSE_UNSTRANDED", "1")], || {
            collapse_loci_span_aware(&anti, &p).len()
        });
        assert_eq!(
            (off, on),
            (2, 2),
            "two MEASURED opposite strands stay two loci in BOTH arms"
        );

        // Control: the strand is the SOLE blocker above. Flip the enclosed model to '-' -- nothing else
        // changes -- and the same pair merges in both arms, so the assertion above is not passing because
        // containment or chrom happened to fail.
        let same = [
            minus,
            DenovoTranscript {
                strand: '-',
                ..plus
            },
        ];
        let ctrl_off = with_collapse_env(&[], || collapse_loci_span_aware(&same, &p).len());
        let ctrl_on = with_collapse_env(&[("RUSTLE_COLLAPSE_UNSTRANDED", "1")], || {
            collapse_loci_span_aware(&same, &p).len()
        });
        assert_eq!(
            (ctrl_off, ctrl_on),
            (1, 1),
            "same-strand control: containment and chrom do pass"
        );
    }

    #[test]
    fn collapse_with_totals_sums_reads_over_the_whole_locus() {
        // Two isoforms of one gene share junction (100,200) -> one locus; a third transcript is a disjoint
        // locus. The locus total is the SUM of every isoform's reads (32+20=52), not the rep isoform's 32.
        let txs = vec![
            tx("A_iso1", "c1", 0, 300, 32, &[(100, 200)], vec![b'A'; 50]),
            tx(
                "A_iso2",
                "c1",
                0,
                400,
                20,
                &[(100, 200), (250, 320)],
                vec![b'A'; 60],
            ),
            tx(
                "B_iso1",
                "c1",
                1000,
                1300,
                10,
                &[(1100, 1200)],
                vec![b'C'; 50],
            ),
        ];
        let (reps, totals) = collapse_loci_span_aware_with_totals(&txs, &DetectParams::default());
        // reps align with the behavior-preserving wrapper.
        assert_eq!(
            reps,
            collapse_loci_span_aware(&txs, &DetectParams::default())
        );
        assert_eq!(reps.len(), 2, "two loci: gene A and gene B");
        // rep of gene A is the max-reads isoform (A_iso1, idx 0); its locus total is 32+20=52.
        let a_pos = reps.iter().position(|&r| r == 0).expect("A_iso1 is a rep");
        assert_eq!(
            totals[a_pos], 52,
            "locus total = sum over all isoforms of the gene"
        );
        let b_pos = reps.iter().position(|&r| r == 2).expect("B_iso1 is a rep");
        assert_eq!(
            totals[b_pos], 10,
            "single-isoform locus total = its own reads"
        );
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
        assert_eq!(
            candidate_pairs(&reps, &DetectParams::default()),
            vec![(0, 1)]
        );
    }

    #[test]
    fn poa_core_completion_attaches_a_divergent_paralog_at_a_new_locus() {
        // read-conflict family {r0,r1} (share the conserved core). r2 is a FREE rep that ALSO shares the
        // conserved core but at a new locus (divergent flanks) — the loosely-related paralog the read-conflict
        // graph misses; the completion attaches it. r3 (unrelated) is not added; free-vs-free is not seeded.
        let core = rand_seq(400, 0x0C0FE);
        let reps = [
            homolog_tx("r0", 0xA1, &core, 0xA2, 5),
            homolog_tx("r1", 0xB1, &core, 0xB2, 5),
            homolog_tx("r2", 0xC1, &core, 0xC2, 5),
            tx("r3", "c1", 0, 560, 5, &[], rand_seq(560, 0xDEAD)),
        ];
        let families = vec![vec![0usize, 1usize]];
        let adds = poa_core_completion_adds(&reps, &families, &DetectParams::default());
        assert_eq!(
            adds,
            vec![vec![2usize]],
            "r2 attaches via the conserved core; r3 does not"
        );
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
        let p = DetectParams {
            cnt_max: 2,
            ..DetectParams::default()
        };
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
        let excl = DetectParams {
            cnt_min: 1,
            cnt_max: 1,
            ..DetectParams::default()
        };
        assert!(
            candidate_pairs(&reps, &excl).is_empty(),
            "core owned by 2 > cnt_max=1 -> excluded"
        );
        let incl = DetectParams {
            cnt_min: 1,
            cnt_max: 2,
            ..DetectParams::default()
        };
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
        assert_eq!(
            candidate_pairs(&[r0, r1], &DetectParams::default()),
            vec![(0, 1)]
        );
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
        let cr =
            confirm_edge(&a, &b, &DetectParams::default()).expect("opposite-strand copy confirms");
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
    fn confirm_edge_over_cap_uses_fallback_not_skip() {
        // over the poasta cap, confirm_edge switches to the linear-memory LCS fallback instead of skipping.
        // DISJOINT random sequences still fail the t_core bar under the fallback -> None (low coverage), NOT
        // an OOM and NOT a blind drop.
        // The length cap only exists on the POA path (LCS is linear), so pin the POA core explicitly.
        let p = DetectParams {
            len_cap: 100,
            edge_core: EdgeCore::Poa,
            ..DetectParams::default()
        };
        let a = rand_seq(560, 0x1);
        let b = rand_seq(560, 0x2);
        assert!(
            confirm_edge(&a, &b, &p).is_none(),
            "disjoint over-cap pair -> low fallback coverage -> None"
        );
        // a true homologous over-cap pair (shared 400 bp core) must now be CONFIRMED via the fallback (the
        // recall the old hard-skip lost): the embedded core gives high LCS coverage.
        let core = rand_seq(400, 0x5);
        let mut h0 = rand_seq(120, 0x6);
        h0.extend_from_slice(&core);
        h0.extend(rand_seq(120, 0x7));
        let mut h1 = rand_seq(120, 0x8);
        h1.extend_from_slice(&core);
        h1.extend(rand_seq(120, 0x9));
        assert!(
            confirm_edge(&h0, &h1, &p).is_some(),
            "over-cap homologous pair confirmed by the bounded fallback (no OOM, no loss)"
        );
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

    #[test]
    fn detect_edges_reporting_uses_fallback_for_oversized_pairs() {
        // a homologous pair that candidate_pairs finds. Under a NORMAL cap it is POA-confirmed and no fallback
        // is used; under a tiny cap the SAME pair is over the poasta threshold -> it is confirmed by the
        // bounded fallback (NOT skipped) AND reported, so the family edge survives and is auditable.
        let core = rand_seq(400, 0xEE88);
        let reps = [
            homolog_tx("r0", 0xA1, &core, 0xA2, 5),
            homolog_tx("r1", 0xB1, &core, 0xB2, 5),
        ];
        let poa = DetectParams {
            edge_core: EdgeCore::Poa,
            ..DetectParams::default()
        };
        let (edges, fb) = detect_edges_reporting(&reps, &poa);
        assert_eq!(edges.len(), 1, "homologous pair confirmed under normal cap");
        assert!(fb.is_empty(), "no fallback used under normal cap");
        // the .0 list still matches plain detect_edges (faithfulness of the delegation).
        assert_eq!(edges, detect_edges(&reps, &poa));

        let p = DetectParams {
            len_cap: 5,
            edge_core: EdgeCore::Poa,
            ..DetectParams::default()
        };
        let (edges2, fb2) = detect_edges_reporting(&reps, &p);
        assert_eq!(
            edges2.len(),
            1,
            "fallback still confirms the homologous pair (no OOM, no loss)"
        );
        assert_eq!((edges2[0].0, edges2[0].1), (0, 1));
        assert_eq!(
            fb2,
            vec![(0, 1)],
            "the fallback-confirmed pair is reported for audit"
        );
    }

    /// AD-HOC MEASUREMENT (not CI, needs external fixture files) -- o1_ledger.md §6iz/§6j0: dump real
    /// `confirm_edge` (T_CORE=0.13 POA core-coverage) output over the 31 real curated-NPIP genomic loci
    /// (`§5e` rung-1 oracle nodes) and the 5 real spliced exon-sum reps (`§5g` rung-1b), to check the
    /// derived `d_max(L) = 1 - exp(-ln(L)/(0.13*L))` formula against the ACTUAL production mechanism
    /// (poasta POA / its bounded fallback) instead of the minimap2 proxy the original 20/20 validation used.
    /// Run with: `RUSTLE_ORACLE_DIR=/mnt/linuxdisk/home/juanfraitu/o1_oracle cargo test --release
    /// --lib family::family_detect::tests::dump_oracle_npip_core_recip -- --ignored --nocapture
    /// > /tmp/oracle_core_recip.tsv`
    #[test]
    #[ignore = "ad-hoc measurement against external NPIP oracle fixtures, not a CI assertion"]
    fn dump_oracle_npip_core_recip() {
        fn read_fasta(path: &str) -> Vec<(String, Vec<u8>)> {
            let content =
                std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
            let mut out = Vec::new();
            let mut name = String::new();
            let mut seq = Vec::new();
            for line in content.lines() {
                if let Some(rest) = line.strip_prefix('>') {
                    if !name.is_empty() {
                        out.push((name.clone(), seq.clone()));
                    }
                    name = rest.split_whitespace().next().unwrap_or("").to_string();
                    seq.clear();
                } else {
                    seq.extend(line.trim().bytes());
                }
            }
            if !name.is_empty() {
                out.push((name, seq));
            }
            out
        }

        let dir = std::env::var("RUSTLE_ORACLE_DIR")
            .unwrap_or_else(|_| "/mnt/linuxdisk/home/juanfraitu/o1_oracle".to_string());
        let p = DetectParams {
            edge_core: EdgeCore::Poa,
            ..DetectParams::default()
        };
        println!("dataset\tname_a\tname_b\tlen_a\tlen_b\tlen_min\tcore_recip\tt_core_pass");
        for (dataset, file) in [
            ("genomic31", "oracle_nodes.fa"),
            ("exonsum5", "oracle_exonsum.fa"),
        ] {
            let path = format!("{dir}/{file}");
            let seqs = read_fasta(&path);
            for i in 0..seqs.len() {
                for j in (i + 1)..seqs.len() {
                    let (na, a) = &seqs[i];
                    let (nb, b) = &seqs[j];
                    let core = confirm_edge(a, b, &p);
                    let len_min = a.len().min(b.len());
                    let (core_str, pass) = match core {
                        Some(v) => (format!("{v:.6}"), v >= p.t_core),
                        None => ("NA".to_string(), false),
                    };
                    println!(
                        "{dataset}\t{na}\t{nb}\t{}\t{}\t{len_min}\t{core_str}\t{pass}",
                        a.len(),
                        b.len()
                    );
                }
            }
        }
    }

    /// AD-HOC MEASUREMENT (not CI) -- o1_ledger.md §6iz/§6j0, part 2: the genomic-scale NPIP oracle above
    /// (25.7 kb mean) is far larger than what `confirm_edge` actually sees in production (spliced
    /// `DenovoTranscript` reps, real median ~3.1 kb per a real gorilla `gw_family_catalog` run) and some
    /// pairs there are pathologically slow for the bounded core-coverage fallback. This measurement instead
    /// uses a MANIFEST of real rep pairs at the REAL production length scale: same-family pairs from a real
    /// catalog run (`ggo_reps.copies.fa`/`.tsv`, size-2..8 families) as a same-family signal, plus an
    /// equal-sized random cross-family sample as a contrast population.
    /// Run with: `RUSTLE_PAIRS_MANIFEST=<tsv with key_a,key_b,ground_truth_same_family>
    /// RUSTLE_PAIRS_FASTA=<fasta with '>family_id|copy_idx' headers> cargo test --release --lib
    /// family::family_detect::tests::dump_real_catalog_core_recip -- --ignored --nocapture`
    #[test]
    #[ignore = "ad-hoc measurement against an external real-catalog pairs manifest, not a CI assertion"]
    fn dump_real_catalog_core_recip() {
        use DetHashMap;

        fn read_fasta_map(path: &str) -> DetHashMap<String, Vec<u8>> {
            let content =
                std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
            let mut out = DetHashMap::default();
            let mut name = String::new();
            let mut seq = Vec::new();
            for line in content.lines() {
                if let Some(rest) = line.strip_prefix('>') {
                    if !name.is_empty() {
                        out.insert(name.clone(), seq.clone());
                    }
                    name = rest.split_whitespace().next().unwrap_or("").to_string();
                    seq.clear();
                } else {
                    seq.extend(line.trim().bytes());
                }
            }
            if !name.is_empty() {
                out.insert(name, seq);
            }
            out
        }

        let manifest_path = std::env::var("RUSTLE_PAIRS_MANIFEST")
            .unwrap_or_else(|_| "/tmp/pairs_manifest.tsv".to_string());
        let fasta_path = std::env::var("RUSTLE_PAIRS_FASTA")
            .unwrap_or_else(|_| "/tmp/reps_subset.fa".to_string());
        let seqs = read_fasta_map(&fasta_path);
        let p = DetectParams::default();

        let manifest = std::fs::read_to_string(&manifest_path)
            .unwrap_or_else(|e| panic!("read {manifest_path}: {e}"));
        println!("key_a\tkey_b\tground_truth_same_family\tlen_a\tlen_b\tlen_min\tcore_recip\tt_core_pass");
        for line in manifest.lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() < 3 {
                continue;
            }
            let (ka, kb, gt) = (f[0].to_string(), f[1].to_string(), f[2].to_string());
            let a = seqs
                .get(&ka)
                .unwrap_or_else(|| panic!("missing seq for {ka}"))
                .clone();
            let b = seqs
                .get(&kb)
                .unwrap_or_else(|| panic!("missing seq for {kb}"))
                .clone();
            let len_a = a.len();
            let len_b = b.len();
            let len_min = len_a.min(len_b);
            // Default edge core (LCS) is linear-time, so no per-pair guard is needed. The old thread +
            // recv_timeout guard leaked uncancellable poasta threads (docs/o1_ledger.md §6j9); for exact POA
            // values use `from_genome`'s `bridge` phase, one killable process per pair.
            let core_str_pass = match confirm_edge(&a, &b, &p) {
                Some(v) => (format!("{v:.6}"), v >= p.t_core),
                None => ("NA".to_string(), false),
            };
            println!(
                "{ka}\t{kb}\t{gt}\t{len_a}\t{len_b}\t{len_min}\t{}\t{}",
                core_str_pass.0, core_str_pass.1
            );
        }
    }

    #[test]
    fn poa_edge_core_escape_hatch_confirms_homologous_and_rc_and_rejects_disjoint() {
        let p = DetectParams {
            edge_core: EdgeCore::Poa,
            ..DetectParams::default()
        };
        let core = rand_seq(400, 0xC0FE_7001);
        let a = cat(&[&rand_seq(80, 0xA1), &core, &rand_seq(80, 0xA2)]);
        let b = cat(&[&rand_seq(80, 0xB1), &core, &rand_seq(80, 0xB2)]);
        assert!(confirm_edge(&a, &b, &p).expect("POA: homologous pair confirms") >= T_CORE);
        assert!(
            confirm_edge(&a, &reverse_complement(&b), &p).is_some(),
            "POA: opposite-strand copy confirms"
        );
        assert!(
            confirm_edge(&rand_seq(560, 0x111), &rand_seq(560, 0x222), &p).is_none(),
            "POA: disjoint rejected"
        );
    }

    #[test]
    fn edge_core_defaults_to_lcs() {
        assert_eq!(DetectParams::default().edge_core, EdgeCore::Lcs);
    }

    #[test]
    fn edge_core_env_value_parsing() {
        assert_eq!(edge_core_from_env_value(None), EdgeCore::Lcs);
        assert_eq!(edge_core_from_env_value(Some("lcs")), EdgeCore::Lcs);
        assert_eq!(edge_core_from_env_value(Some("poa")), EdgeCore::Poa);
        assert_eq!(edge_core_from_env_value(Some("POA")), EdgeCore::Poa);
        assert_eq!(edge_core_from_env_value(Some("")), EdgeCore::Lcs);
    }

    #[test]
    fn lcs_edge_core_admits_a_shared_core_at_different_offsets_and_lengths() {
        let core = rand_seq(600, 11);
        let a = cat(&[&rand_seq(400, 12), &core, &rand_seq(300, 13)]);
        let b = cat(&[&rand_seq(2500, 14), &core, &rand_seq(900, 15)]);
        let p = DetectParams {
            edge_core: EdgeCore::Lcs,
            ..DetectParams::default()
        };
        let v = confirm_edge(&a, &b, &p).expect("shared 600 bp core must confirm under LCS");
        assert!(
            (v - 600.0 / a.len() as f64).abs() < 0.01,
            "core fraction {v}"
        );
    }

    #[test]
    fn lcs_edge_core_uses_the_reverse_complement_orientation() {
        let core = rand_seq(500, 21);
        let a = cat(&[&rand_seq(300, 22), &core, &rand_seq(300, 23)]);
        let b = cat(&[
            &rand_seq(700, 24),
            &reverse_complement(&core),
            &rand_seq(200, 25),
        ]);
        let p = DetectParams {
            edge_core: EdgeCore::Lcs,
            ..DetectParams::default()
        };
        let v = confirm_edge(&a, &b, &p).expect("reverse-complement core must confirm under LCS");
        assert!(v >= 500.0 / a.len() as f64 - 1e-9);
    }

    #[test]
    fn lcs_edge_core_rejects_unrelated_sequences() {
        let a = rand_seq(1500, 31);
        let b = rand_seq(2200, 32);
        let p = DetectParams {
            edge_core: EdgeCore::Lcs,
            ..DetectParams::default()
        };
        assert_eq!(confirm_edge(&a, &b, &p), None);
    }

    #[test]
    fn lcs_edge_core_is_case_insensitive() {
        let core = rand_seq(400, 41);
        let a = cat(&[&core, &rand_seq(600, 42)]);
        let b_lower: Vec<u8> = cat(&[&rand_seq(300, 43), &core]).to_ascii_lowercase();
        let p = DetectParams {
            edge_core: EdgeCore::Lcs,
            ..DetectParams::default()
        };
        assert!(confirm_edge(&a, &b_lower, &p).is_some());
    }
}

// ---- merged 2026-10-05: was `vg_family/family_graph.rs`, now the inline module below (one component) ----
#[allow(clippy::all)]
pub mod family_graph {
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
    use crate::types::{DetHashMap, DetHashSet};

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
        let astar =
            matches!(std::env::var("RUSTLE_POA_ASTAR"), Ok(ref v) if v != "0" && !v.is_empty());
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
        st.push(St {
            len: 0,
            link: -1,
            next: BTreeMap::new(),
        });
        let mut last = 0i32;
        for &c in s {
            let cur = st.len() as i32;
            let cur_len = st[last as usize].len + 1;
            st.push(St {
                len: cur_len,
                link: -1,
                next: BTreeMap::new(),
            });
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
                    st.push(St {
                        len: clone_len,
                        link: qlink,
                        next: qnext,
                    });
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
    /// Exactness is unaffected: `DetHashMap` still compares keys with `Eq`, so a hash collision costs a byte
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
            assert_eq!(
                contiguous_core_coverage_bounded(&a, &b, 10_000),
                MEMO_SENTINEL
            );
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
            assert_eq!(
                swapped,
                contiguous_core_coverage_bounded_uncached(&b, &a, 10_000)
            );
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
            assert!(
                unset,
                "A* is the DEFAULT on the contiguous-core path (measured 32-44% faster there)"
            );
            assert!(
                !zero,
                "RUSTLE_POA_ASTAR=0 must restore the Dijkstra search on this path"
            );
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
            assert!(
                matches!(upper_cow(&up), Cow::Borrowed(_)),
                "no lowercase -> no allocation"
            );
            let mixed = b"acgtNACgt".to_vec();
            let got = upper_cow(&mixed);
            assert!(
                matches!(got, Cow::Owned(_)),
                "lowercase present -> must allocate"
            );
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
            assert!(
                cov >= 0.13,
                "long shared core should give HIGH contiguous-core coverage (got {cov:.3})"
            );
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
            assert!(
                cov >= 0.5,
                "two true copies sharing a 400 bp core with divergent 5' AND 3' flanks \
             must score HIGH contiguous-core coverage (~0.62); got {cov:.3} — the \
             aligner threaded the flanks and lost the core"
            );
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
            assert!(
                cov < 0.13,
                "short shared block should give LOW contiguous-core coverage (got {cov:.3})"
            );
        }

        #[test]
        fn contiguous_core_coverage_identical_is_one() {
            let s = core_rand_seq(200, 0x1234_5678);
            let cov = contiguous_core_coverage(&s, &s);
            assert!(
                cov >= 0.99,
                "identical sequences should give contiguous-core coverage ~1.0 (got {cov:.3})"
            );
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
            assert_eq!(
                longest_common_substring(b"GATTACA", b"TTTGATTACAGGG"),
                7,
                "short ⊂ long"
            );
            assert_eq!(
                longest_common_substring(b"TTTGATTACAGGG", b"GATTACA"),
                7,
                "order-independent"
            );
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
                assert_eq!(
                    longest_common_substring(&a, &b),
                    naive(&a, &b),
                    "seed {seed}"
                );
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
            assert!(
                lcs >= 0.6,
                "fallback recovers the 400 bp core fraction (got {lcs:.3})"
            );
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
}

// ---- merged 2026-10-05: was `vg_family/family_split.rs`, now the inline module below (one component) ----
#[allow(clippy::all)]
pub mod family_split {
    //! De-novo family DECOMPOSITION — the `bench/denovo_family_split.py` core.
    //!
    //! Connected components of the POA-homology edge graph transitively close a SUPERFAMILY (e.g. KRAB-ZNF)
    //! into one blob, because genes that merely share a tandem DOMAIN pairwise-confirm. A real recent-duplicate
    //! family is a DENSE, mutually-homologous subgraph; an over-merge is a SPARSE chain bridged by lower
    //! `core_recip` domain edges. We keep the same POA edges but split large components by **weighted
    //! modularity** (uses the `core_recip` edge weight; no arbitrary similarity cutoff), then FLAG sparse
    //! homology webs (`size >= web_min_size & density < web_max_density`) rather than claim them as families.
    //!
    //! Community detection is a **hand-rolled, deterministic** weighted-modularity Louvain (the Python uses
    //! `networkx.community.louvain_communities`, seeded). Byte-matching networkx is infeasible (its greedy
    //! order + RNG); the Python pipeline established the partition is robust (resolution 1→8 holds 96–96.7%
    //! concordance), so the faithful target is the same modularity DEFINITION with deterministic tie-breaks
    //! (sorted node/community order, strict-improvement moves), not a byte-identical partition.
    //!
    //! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)
    use crate::types::{DetHashMap, DetHashSet};

    use std::collections::{BTreeMap, BTreeSet};

    /// Only decompose connected components this size or larger (small families left intact).
    pub const MIN_DECOMP: usize = 6;
    /// Modularity resolution (`gamma`); 1.0 = standard.
    pub const RESOLUTION: f64 = 1.0;
    /// A final community of `>= WEB_MIN_SIZE` nodes ...
    pub const WEB_MIN_SIZE: usize = 10;
    /// ... and `< WEB_MAX_DENSITY` density is flagged a homology web (a sparse domain-sharing over-merge),
    /// not a family — and Web families are EXCLUDED from copy-assignment (denovo_pipeline.rs).
    ///
    /// Set to 0.30 to ALIGN with the validated DNA cDNA-homology manifest bar (`make_dna_family_manifest.py`
    /// drops `overmerge_sparse` at `density < 0.30`). Measured on the real de-novo split (698 families): the
    /// family-density distribution is BIMODAL — median 1.000 (cliques) with NOTHING legitimate in the
    /// [0.15, 0.30) band; the only 4 large (n>=10) families there are multi-chromosome domain-sharing
    /// over-merges (DSFAM0 = a 164-member ZNF spanning 19 chromosomes, etc.) that the old 0.15 bar let
    /// through as "families." `WEB_MIN_SIZE` stays 10 (not the manifest's 4) ON PURPOSE: a SMALL sparse
    /// group can be a real divergent family (e.g. a 7-copy single-chrom MAGEB at density 0.24), so only
    /// LARGE-and-sparse is called a web. See loose end L2 in bench/LOOSE_ENDS_AUDIT.md.
    pub const WEB_MAX_DENSITY: f64 = 0.30;

    /// Tunable decomposition parameters (defaults mirror `denovo_family_split.py`).
    #[derive(Clone, Copy, Debug)]
    pub struct SplitParams {
        pub min_decomp: usize,
        pub resolution: f64,
        pub web_min_size: usize,
        pub web_max_density: f64,
    }

    impl Default for SplitParams {
        fn default() -> Self {
            SplitParams {
                min_decomp: MIN_DECOMP,
                resolution: RESOLUTION,
                web_min_size: WEB_MIN_SIZE,
                web_max_density: WEB_MAX_DENSITY,
            }
        }
    }

    /// A discrete recent-duplicate family vs a non-discretizable homology web.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub enum FamilyClass {
        Family,
        Web,
    }

    /// Structural diagnostics of a community in the edge graph.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct CommunityStats {
        pub n: usize,
        pub n_edges: usize,
        pub density: f64,
        pub avg_core_recip: f64,
        pub n_articulation: usize,
        /// Global edge connectivity λ of the induced subgraph — the minimum number of edges whose removal
        /// disconnects this community. `λ >= 2` certifies that no single alignment record's loss can split
        /// it. **Reported, never used to decide membership** (a 2-node family has λ = 1 necessarily); see
        /// `edge_connectivity`.
        pub lambda: usize,
    }

    /// A final decomposed family with its structural diagnostics and class.
    #[derive(Clone, Debug)]
    pub struct SplitFamily {
        pub members: Vec<usize>,
        pub stats: CommunityStats,
        pub class: FamilyClass,
    }

    /// Connected components (each a sorted node list) of the edge graph, keeping components of `>= min_size`.
    /// Nodes are the ids appearing in `edges`; isolated nodes (no edge) are not represented.
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

    pub fn connected_components(edges: &[(usize, usize, f64)], min_size: usize) -> Vec<Vec<usize>> {
        let max_id = match edges.iter().map(|&(a, b, _)| a.max(b)).max() {
            Some(m) => m,
            None => return Vec::new(),
        };
        let mut parent: Vec<usize> = (0..=max_id).collect();
        let mut present = vec![false; max_id + 1];
        for &(a, b, _) in edges {
            present[a] = true;
            present[b] = true;
            uf_union(&mut parent, a, b);
        }
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for node in 0..=max_id {
            if present[node] {
                let r = uf_find(&mut parent, node);
                groups.entry(r).or_default().push(node);
            }
        }
        let mut comps: Vec<Vec<usize>> = groups
            .into_values()
            .filter(|c| c.len() >= min_size)
            .collect();
        for c in &mut comps {
            c.sort_unstable();
        }
        // deterministic: size desc, then smallest member.
        comps.sort_by(|a, b| b.len().cmp(&a.len()).then_with(|| a[0].cmp(&b[0])));
        comps
    }

    /// Hand-rolled deterministic weighted-modularity Louvain over a graph with local node ids `0..n`.
    /// Returns the communities (each a sorted Vec of local node ids). `resolution` is the modularity gamma.
    pub fn louvain_communities(
        n: usize,
        edges: &[(usize, usize, f64)],
        resolution: f64,
    ) -> Vec<Vec<usize>> {
        if n == 0 {
            return Vec::new();
        }
        let labels = louvain_labels(n, edges, resolution);
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (node, &lab) in labels.iter().enumerate() {
            groups.entry(lab).or_default().push(node);
        }
        let mut comms: Vec<Vec<usize>> = groups.into_values().collect();
        for c in &mut comms {
            c.sort_unstable();
        }
        comms.sort_by(|a, b| a[0].cmp(&b[0])); // deterministic: by smallest member
        comms
    }

    /// Multi-level Louvain: run local-moving to convergence, aggregate communities into super-nodes, repeat
    /// until a level makes no moves. Returns a community label (renumbered `0..k`) for each original node.
    fn louvain_labels(n: usize, edges: &[(usize, usize, f64)], resolution: f64) -> Vec<usize> {
        let mut node2super: Vec<usize> = (0..n).collect();
        let mut cur_n = n;
        let mut cur_edges: Vec<(usize, usize, f64)> = edges.to_vec();
        loop {
            let (labels, moved) = louvain_one_level(cur_n, &cur_edges, resolution);
            for s in node2super.iter_mut() {
                *s = labels[*s];
            }
            let k = labels.iter().copied().max().map(|m| m + 1).unwrap_or(0);
            if !moved || k >= cur_n {
                break;
            }
            // aggregate: each community becomes a super-node; internal edges become self-loops.
            let mut agg: BTreeMap<(usize, usize), f64> = BTreeMap::new();
            for &(i, j, w) in &cur_edges {
                let (ci, cj) = (labels[i], labels[j]);
                *agg.entry((ci.min(cj), ci.max(cj))).or_insert(0.0) += w;
            }
            cur_edges = agg.into_iter().map(|((a, b), w)| (a, b, w)).collect();
            cur_n = k;
        }
        // renumber the final super-node ids to a compact 0..K.
        let mut remap: BTreeMap<usize, usize> = BTreeMap::new();
        let mut next = 0usize;
        node2super
            .iter()
            .map(|&s| {
                *remap.entry(s).or_insert_with(|| {
                    let v = next;
                    next += 1;
                    v
                })
            })
            .collect()
    }

    /// One Louvain level: greedy local moving to a local modularity optimum. Returns `(labels 0..k, moved)`.
    /// Weighted-modularity gain (python-louvain formulation): moving node `i` into community `C` gains
    /// `remove_cost + w(i,C) - γ·Σ_tot(C)·k_i/(2m)`, where `remove_cost = -w(i,C_old) + γ·(Σ_tot(C_old)-k_i)·k_i/(2m)`.
    /// Determinism: nodes scanned in id order; neighbour communities in sorted (BTreeMap) order; only STRICT
    /// improvements move, so the lowest-id community wins ties and the node stays put on a non-positive best.
    fn louvain_one_level(
        n: usize,
        edges: &[(usize, usize, f64)],
        resolution: f64,
    ) -> (Vec<usize>, bool) {
        let mut deg = vec![0.0f64; n];
        let mut neigh: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
        let mut m = 0.0f64;
        for &(i, j, w) in edges {
            m += w;
            if i == j {
                deg[i] += 2.0 * w; // a self-loop counts twice toward degree
            } else {
                deg[i] += w;
                deg[j] += w;
                neigh[i].push((j, w));
                neigh[j].push((i, w));
            }
        }
        if m <= 0.0 {
            return ((0..n).collect(), false);
        }
        let two_m = 2.0 * m;
        let mut com: Vec<usize> = (0..n).collect();
        let mut com_deg: Vec<f64> = deg.clone(); // Σ_tot of each community (community id ∈ 0..n)
        let mut any_moved = false;
        let mut improved = true;
        while improved {
            improved = false;
            for node in 0..n {
                let c_old = com[node];
                let dctw = deg[node] / two_m;
                let mut ncw: BTreeMap<usize, f64> = BTreeMap::new();
                for &(nb, w) in &neigh[node] {
                    *ncw.entry(com[nb]).or_insert(0.0) += w;
                }
                let dnc_old = *ncw.get(&c_old).unwrap_or(&0.0);
                let remove_cost = -dnc_old + resolution * (com_deg[c_old] - deg[node]) * dctw;
                com_deg[c_old] -= deg[node]; // isolate the node
                let mut best_com = c_old;
                let mut best_inc = 0.0f64;
                for (&c, &dnc) in &ncw {
                    let inc = remove_cost + dnc - resolution * com_deg[c] * dctw;
                    if inc > best_inc {
                        best_inc = inc;
                        best_com = c;
                    }
                }
                com_deg[best_com] += deg[node];
                com[node] = best_com;
                if best_com != c_old {
                    improved = true;
                    any_moved = true;
                }
            }
        }
        let mut remap: BTreeMap<usize, usize> = BTreeMap::new();
        let mut next = 0usize;
        let labels: Vec<usize> = com
            .iter()
            .map(|&c| {
                *remap.entry(c).or_insert_with(|| {
                    let v = next;
                    next += 1;
                    v
                })
            })
            .collect();
        (labels, any_moved)
    }

    /// Articulation (cut) points of an undirected graph with local node ids `0..n`. Returns sorted ids.
    pub fn articulation_points(n: usize, edges: &[(usize, usize)]) -> Vec<usize> {
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for &(a, b) in edges {
            if a != b {
                adj[a].push(b);
                adj[b].push(a);
            }
        }
        let mut visited = vec![false; n];
        let mut disc = vec![0usize; n];
        let mut low = vec![0usize; n];
        let mut is_ap = vec![false; n];
        let mut timer = 0usize;
        // iterative Tarjan (avoids recursion-depth blowups on large families).
        for start in 0..n {
            if visited[start] {
                continue;
            }
            let mut stack: Vec<(usize, isize, usize)> = vec![(start, -1, 0)]; // (node, parent, child cursor)
            let mut root_children = 0usize;
            while let Some(&(u, parent, ci)) = stack.last() {
                if ci == 0 {
                    visited[u] = true;
                    timer += 1;
                    disc[u] = timer;
                    low[u] = timer;
                }
                if ci < adj[u].len() {
                    stack.last_mut().unwrap().2 += 1;
                    let v = adj[u][ci];
                    if v as isize == parent {
                        continue;
                    }
                    if !visited[v] {
                        if parent == -1 {
                            root_children += 1;
                        }
                        stack.push((v, u as isize, 0));
                    } else {
                        low[u] = low[u].min(disc[v]);
                    }
                } else {
                    stack.pop();
                    if let Some(&(p, pparent, _)) = stack.last() {
                        low[p] = low[p].min(low[u]);
                        if pparent != -1 && low[u] >= disc[p] {
                            is_ap[p] = true; // non-root cut vertex
                        }
                    }
                }
            }
            if root_children > 1 {
                is_ap[start] = true;
            }
        }
        (0..n).filter(|&i| is_ap[i]).collect()
    }

    /// Global EDGE CONNECTIVITY λ of an undirected graph on local ids `0..n` (Stoer–Wagner). Every entry of
    /// `edges` contributes weight 1, so PARALLEL EDGES ADD — de-duplicate before calling if that is not wanted
    /// (`community_stats` passes a de-duplicated induced set).
    ///
    /// λ = the minimum number of edges whose removal disconnects the graph. Returns **0** when `n < 2` or when
    /// the graph is already disconnected: in both cases no edge has to be paid to separate two nodes.
    ///
    /// WHY THIS EXISTS — and why it is NOT part of the family definition. `λ >= 2` is a per-family
    /// CERTIFICATE: it states that **no single alignment record's loss can split this family**. It is
    /// deliberately not a membership criterion, because a 2-copy family has `λ = 1` NECESSARILY (one edge is
    /// all a 2-node graph can have), so gating membership on `λ >= 2` would delete every 2-copy family —
    /// the most common family there is. It reports confidence; it never decides membership.
    /// See `docs/seeded_family_definition.md` §1★.5.
    ///
    /// Deterministic: the maximum-adjacency search breaks ties toward the smallest node id.
    pub fn edge_connectivity(n: usize, edges: &[(usize, usize)]) -> usize {
        if n < 2 {
            return 0;
        }
        let mut w = vec![vec![0usize; n]; n];
        for &(a, b) in edges {
            if a != b && a < n && b < n {
                w[a][b] += 1;
                w[b][a] += 1;
            }
        }
        let mut active: Vec<usize> = (0..n).collect();
        let mut best = usize::MAX;
        while active.len() > 1 {
            // maximum-adjacency search over the surviving supernodes
            let mut in_a = vec![false; n];
            let mut wsum = vec![0usize; n];
            let (mut prev, mut last) = (usize::MAX, usize::MAX);
            for _ in 0..active.len() {
                let mut sel = usize::MAX;
                for &v in &active {
                    // strict `>` with `active` ascending => ties go to the smallest id (deterministic)
                    if !in_a[v] && (sel == usize::MAX || wsum[v] > wsum[sel]) {
                        sel = v;
                    }
                }
                in_a[sel] = true;
                prev = last;
                last = sel;
                for &v in &active {
                    if !in_a[v] {
                        wsum[v] += w[sel][v];
                    }
                }
            }
            // cut-of-the-phase = weight from the last-added node to everything else (all of which are in A)
            let cut: usize = active
                .iter()
                .filter(|&&v| v != last)
                .map(|&v| w[last][v])
                .sum();
            best = best.min(cut);
            if best == 0 {
                return 0; // already disconnected; no smaller cut exists
            }
            if prev == usize::MAX {
                break;
            }
            // merge `last` into `prev`
            for &v in &active {
                if v != last && v != prev {
                    w[prev][v] += w[last][v];
                    w[v][prev] = w[prev][v];
                }
            }
            active.retain(|&v| v != last);
        }
        if best == usize::MAX {
            0
        } else {
            best
        }
    }

    /// Structural diagnostics of the subgraph induced by `members` (global node ids) over `edges`.
    pub(crate) fn community_stats(
        members: &[usize],
        edges: &[(usize, usize, f64)],
    ) -> CommunityStats {
        let set: BTreeSet<usize> = members.iter().copied().collect();
        let n = set.len();
        let mut n_edges = 0usize;
        let mut wsum = 0.0f64;
        let mut internal: Vec<(usize, usize)> = Vec::new();
        for &(a, b, w) in edges {
            if a != b && set.contains(&a) && set.contains(&b) {
                n_edges += 1;
                wsum += w;
                internal.push((a, b));
            }
        }
        let density = if n > 1 {
            2.0 * n_edges as f64 / (n as f64 * (n as f64 - 1.0))
        } else {
            1.0
        };
        let avg_core_recip = if n_edges > 0 {
            wsum / n_edges as f64
        } else {
            0.0
        };
        // Local relabelling `member -> 0..n`, shared by the articulation and λ passes.
        let idx: BTreeMap<usize, usize> =
            set.iter().enumerate().map(|(i, &node)| (node, i)).collect();
        let local: Vec<(usize, usize)> =
            internal.iter().map(|&(a, b)| (idx[&a], idx[&b])).collect();
        // articulation only for n > 2 (python `arts = ... if n > 2 else 0`).
        let n_articulation = if n > 2 {
            articulation_points(n, &local).len()
        } else {
            0
        };
        // λ over the DE-DUPLICATED induced edge set: `edges` may list a pair more than once, and
        // `edge_connectivity` counts every entry as weight 1, which would inflate the cut.
        let dedup: BTreeSet<(usize, usize)> =
            local.iter().map(|&(a, b)| (a.min(b), a.max(b))).collect();
        let lambda = edge_connectivity(n, &dedup.iter().copied().collect::<Vec<_>>());
        CommunityStats {
            n,
            n_edges,
            density,
            avg_core_recip,
            n_articulation,
            lambda,
        }
    }

    /// Classify a community by size + density (web iff `n >= web_min_size && density < web_max_density`).
    pub fn classify(n: usize, density: f64, p: &SplitParams) -> FamilyClass {
        if n >= p.web_min_size && density < p.web_max_density {
            FamilyClass::Web
        } else {
            FamilyClass::Family
        }
    }

    /// Decompose the POA-homology edge graph into final families: connected components, with components
    /// `>= min_decomp` split by weighted modularity (communities of `>= 2` kept), each annotated with
    /// structural diagnostics and a family/web class. Sorted by size desc, then smallest member.
    pub fn decompose_families(edges: &[(usize, usize, f64)], p: &SplitParams) -> Vec<SplitFamily> {
        let comps = connected_components(edges, 2);
        let mut final_members: Vec<Vec<usize>> = Vec::new();
        for comp in comps {
            if comp.len() < p.min_decomp {
                final_members.push(comp);
                continue;
            }
            // decompose: relabel the component to local ids 0..k, Louvain, map communities (>= 2) back.
            let set: BTreeSet<usize> = comp.iter().copied().collect();
            let idx: BTreeMap<usize, usize> =
                comp.iter().enumerate().map(|(i, &nd)| (nd, i)).collect();
            let local_edges: Vec<(usize, usize, f64)> = edges
                .iter()
                .filter(|&&(a, b, _)| a != b && set.contains(&a) && set.contains(&b))
                .map(|&(a, b, w)| (idx[&a], idx[&b], w))
                .collect();
            for community in louvain_communities(comp.len(), &local_edges, p.resolution) {
                if community.len() >= 2 {
                    final_members.push(community.iter().map(|&li| comp[li]).collect());
                }
            }
        }
        let mut out: Vec<SplitFamily> = final_members
            .into_iter()
            .map(|mut members| {
                members.sort_unstable();
                let stats = community_stats(&members, edges);
                let class = classify(stats.n, stats.density, p);
                SplitFamily {
                    members,
                    stats,
                    class,
                }
            })
            .collect();
        // deterministic output order: size desc, then smallest member.
        out.sort_by(|a, b| {
            b.members
                .len()
                .cmp(&a.members.len())
                .then_with(|| a.members[0].cmp(&b.members[0]))
        });
        out
    }

    /// Connected components of the graph on nodes `0..n` over `edges`, INCLUDING size-1 (degree-0) singletons.
    /// Unlike `conflict_families` (which drops <2-node components), this guarantees every node in `0..n`
    /// appears in exactly one returned component, so callers can partition ALL of `0..n`. Components and their
    /// members are returned in ascending id order (deterministic).
    fn all_components(n: usize, edges: &[(usize, usize, f64)]) -> Vec<Vec<usize>> {
        let mut parent: Vec<usize> = (0..n).collect();
        for &(a, b, _) in edges {
            if a < n && b < n {
                uf_union(&mut parent, a, b);
            }
        }
        let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for node in 0..n {
            let r = uf_find(&mut parent, node);
            groups.entry(r).or_default().push(node);
        }
        groups.into_values().collect()
    }

    /// Internal edge density of a node block over the induced subgraph: 2|E|/(|C|(|C|-1)); <=1 node = 1.0.
    fn induced_density(members: &[usize], edges: &[(usize, usize, f64)]) -> f64 {
        let set: DetHashSet<usize> = members.iter().copied().collect();
        let n = members.len();
        if n <= 1 {
            return 1.0;
        }
        let m = edges
            .iter()
            .filter(|(a, b, _)| set.contains(a) && set.contains(b))
            .count();
        2.0 * m as f64 / (n as f64 * (n as f64 - 1.0))
    }

    /// Restrict edges to those with both endpoints in `members`, remapped to local indices 0..members.len().
    fn induced_edges(
        members: &[usize],
        edges: &[(usize, usize, f64)],
    ) -> (usize, Vec<(usize, usize, f64)>) {
        let idx: DetHashMap<usize, usize> =
            members.iter().enumerate().map(|(i, &g)| (g, i)).collect();
        let local = edges
            .iter()
            .filter_map(|&(a, b, w)| match (idx.get(&a), idx.get(&b)) {
                (Some(&la), Some(&lb)) => Some((la, lb, w)),
                _ => None,
            })
            .collect();
        (members.len(), local)
    }

    /// Guaranteed-progress splitter: adaptive-resolution Louvain, else component split, else deterministic halving.
    /// Returns global-index blocks; never a single block equal to the input when |members| > 2.
    fn split_once(members: &[usize], edges: &[(usize, usize, f64)]) -> Vec<Vec<usize>> {
        let (n, local) = induced_edges(members, edges);
        for res in [1.0, 2.0, 4.0, 8.0] {
            let parts = louvain_communities(n, &local, res);
            if parts.len() >= 2 {
                return parts
                    .into_iter()
                    .map(|p| p.into_iter().map(|l| members[l]).collect())
                    .collect();
            }
        }
        // connected-component fallback (also catches a disconnected block). Use `all_components` (NOT
        // conflict_families) so an isolated member is kept as its own component, never dropped.
        let comps = all_components(n, &local);
        if comps.len() >= 2 {
            return comps
                .into_iter()
                .map(|c| c.into_iter().map(|l| members[l]).collect())
                .collect();
        }
        // deterministic halving.
        let h = members.len() / 2;
        vec![members[..h].to_vec(), members[h..].to_vec()]
    }

    /// gamma-quasi-clique partition: keep a block whole iff it is already a gamma-quasi-clique (or <=2 nodes),
    /// else split (guaranteed-progress) and recurse. Blocks partition 0..n. Deterministic (Louvain is
    /// deterministic here).
    pub fn gamma_quasi_clique_partition(
        n: usize,
        edges: &[(usize, usize, f64)],
        gamma: f64,
    ) -> Vec<Vec<usize>> {
        // start from raw connected components (INCLUDING singletons, so the output partitions ALL of 0..n),
        // then refine each.
        let comps = all_components(n, edges);
        let mut out = Vec::new();
        let mut stack: Vec<Vec<usize>> = comps;
        while let Some(block) = stack.pop() {
            if block.len() <= 2 || induced_density(&block, edges) >= gamma {
                out.push(block);
                continue;
            }
            let parts = split_once(&block, edges);
            if parts.len() == 1 && parts[0].len() == block.len() {
                // provably unreachable: split_once always returns >=2 strictly-smaller parts for len>=3.
                out.push(block);
            } else {
                stack.extend(parts);
            }
        }
        out
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        fn we(a: usize, b: usize) -> (usize, usize, f64) {
            (a, b, 1.0)
        }

        // ---- connected_components ----

        #[test]
        fn cc_two_separate_edges_two_components() {
            let e = [we(0, 1), we(2, 3)];
            let mut c = connected_components(&e, 2);
            c.sort();
            assert_eq!(c, vec![vec![0, 1], vec![2, 3]]);
        }

        #[test]
        fn cc_transitive_one_component() {
            let e = [we(0, 1), we(1, 2)];
            assert_eq!(connected_components(&e, 2), vec![vec![0, 1, 2]]);
        }

        // ---- louvain_communities ----

        #[test]
        fn louvain_single_edge_one_community() {
            assert_eq!(louvain_communities(2, &[we(0, 1)], 1.0), vec![vec![0, 1]]);
        }

        #[test]
        fn louvain_single_clique_one_community() {
            let e = [we(0, 1), we(1, 2), we(0, 2)];
            assert_eq!(louvain_communities(3, &e, 1.0), vec![vec![0, 1, 2]]);
        }

        #[test]
        fn louvain_two_cliques_bridge_splits() {
            // two triangles {0,1,2} and {3,4,5}, strong internal weight, bridged by ONE weak edge.
            let e = [
                we(0, 1),
                we(1, 2),
                we(0, 2),
                we(3, 4),
                we(4, 5),
                we(3, 5),
                (0, 3, 0.05),
            ];
            let mut comms = louvain_communities(6, &e, 1.0);
            comms.sort();
            assert_eq!(comms, vec![vec![0, 1, 2], vec![3, 4, 5]]);
        }

        // ---- articulation_points ----

        #[test]
        fn artic_path_middle_is_cut() {
            assert_eq!(articulation_points(3, &[(0, 1), (1, 2)]), vec![1]);
        }

        #[test]
        fn artic_triangle_none() {
            assert_eq!(
                articulation_points(3, &[(0, 1), (1, 2), (0, 2)]),
                Vec::<usize>::new()
            );
        }

        #[test]
        fn artic_path_of_four_two_cuts() {
            assert_eq!(
                articulation_points(4, &[(0, 1), (1, 2), (2, 3)]),
                vec![1, 2]
            );
        }

        // ---- community_stats ----

        #[test]
        fn stats_clique_density_one() {
            let e = [(0, 1, 0.5), (1, 2, 0.7), (0, 2, 0.9)];
            let s = community_stats(&[0, 1, 2], &e);
            assert_eq!(s.n, 3);
            assert_eq!(s.n_edges, 3);
            assert!((s.density - 1.0).abs() < 1e-9);
            assert!((s.avg_core_recip - (0.5 + 0.7 + 0.9) / 3.0).abs() < 1e-9);
            assert_eq!(s.n_articulation, 0);
            assert_eq!(s.lambda, 2, "K3 is 2-edge-connected");
        }

        // ── λ (EDGE CONNECTIVITY) — the per-family certificate ────────────────────────────────────────
        //
        // λ is REPORTED, never used to decide membership. The `lambda_two_node_is_one` case is the reason:
        // a 2-copy family cannot have λ >= 2, so gating membership on λ would delete every 2-copy family.

        #[test]
        fn lambda_complete_graph_is_n_minus_one() {
            for n in 2..=6usize {
                let edges: Vec<(usize, usize)> = (0..n)
                    .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
                    .collect();
                assert_eq!(
                    edge_connectivity(n, &edges),
                    n - 1,
                    "K{n} must have lambda = {}",
                    n - 1
                );
            }
        }

        #[test]
        fn lambda_two_node_is_one() {
            // THE REASON λ IS NOT A MEMBERSHIP CRITERION: one edge is all a 2-node graph can hold.
            assert_eq!(edge_connectivity(2, &[(0, 1)]), 1);
        }

        #[test]
        fn lambda_path_and_cycle() {
            // path 0-1-2-3 hangs on any single edge; the 4-cycle needs two cuts.
            assert_eq!(edge_connectivity(4, &[(0, 1), (1, 2), (2, 3)]), 1);
            assert_eq!(edge_connectivity(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]), 2);
        }

        #[test]
        fn lambda_disconnected_and_degenerate_are_zero() {
            assert_eq!(
                edge_connectivity(4, &[(0, 1), (2, 3)]),
                0,
                "already disconnected"
            );
            assert_eq!(edge_connectivity(3, &[(0, 1)]), 0, "isolated node 2");
            assert_eq!(edge_connectivity(1, &[]), 0, "single node: no cut to pay");
            assert_eq!(edge_connectivity(0, &[]), 0);
        }

        #[test]
        fn lambda_two_cliques_joined_by_a_bridge_is_one() {
            // K3 {0,1,2} -- bridge 2-3 -- K3 {3,4,5}: dense, but one record's loss splits it.
            let edges = [(0, 1), (0, 2), (1, 2), (2, 3), (3, 4), (3, 5), (4, 5)];
            assert_eq!(edge_connectivity(6, &edges), 1);
            let s = community_stats(
                &[0, 1, 2, 3, 4, 5],
                &edges.iter().map(|&(a, b)| (a, b, 1.0)).collect::<Vec<_>>(),
            );
            assert_eq!(s.lambda, 1);
            assert_eq!(
                s.n_articulation, 2,
                "the bridge's two endpoints are cut vertices"
            );
        }

        #[test]
        fn lambda_ignores_duplicate_records_for_the_same_pair() {
            // `community_stats` de-duplicates before λ: three alignment records for one pair is still ONE
            // edge whose loss splits the family, so λ must stay 1 and not be inflated to 3.
            let e = [(0, 1, 1.0), (0, 1, 1.0), (0, 1, 1.0)];
            assert_eq!(community_stats(&[0, 1], &e).lambda, 1);
        }

        #[test]
        fn stats_path_density_and_articulation() {
            // path 0-1-2-3: 3 edges, density = 2*3/(4*3) = 0.5, articulation points {1,2}.
            let e = [(0, 1, 0.2), (1, 2, 0.2), (2, 3, 0.2)];
            let s = community_stats(&[0, 1, 2, 3], &e);
            assert_eq!(s.n_edges, 3);
            assert!((s.density - 0.5).abs() < 1e-9);
            assert_eq!(s.n_articulation, 2);
        }

        #[test]
        fn stats_only_counts_internal_edges() {
            // members {0,1} but an edge 1-2 leaves the set -> not counted.
            let e = [(0, 1, 0.4), (1, 2, 0.9)];
            let s = community_stats(&[0, 1], &e);
            assert_eq!(s.n_edges, 1);
            assert!((s.avg_core_recip - 0.4).abs() < 1e-9);
        }

        // ---- classify ----

        #[test]
        fn classify_web_vs_family() {
            let p = SplitParams::default(); // web_max_density = 0.30, web_min_size = 10
            assert_eq!(classify(10, 0.10, &p), FamilyClass::Web); // size>=10 & sparse
            assert_eq!(classify(12, 0.167, &p), FamilyClass::Web); // the DSFAM0-class over-merge (164 ZNF/19 chr): was Family at the old 0.15 bar
            assert_eq!(classify(10, 0.29, &p), FamilyClass::Web); // just under the aligned 0.30 bar
            assert_eq!(classify(10, 0.40, &p), FamilyClass::Family); // dense -> family
            assert_eq!(classify(10, 1.00, &p), FamilyClass::Family); // clique
            assert_eq!(classify(5, 0.10, &p), FamilyClass::Family); // too SMALL -> family even if sparse (small divergent fam, e.g. MAGEB)
            assert_eq!(classify(7, 0.24, &p), FamilyClass::Family); // n<10 stays family (the MAGEB-class protection)
        }

        // ---- decompose_families ----

        #[test]
        fn decompose_small_component_kept_intact() {
            let e = [we(0, 1), we(1, 2), we(0, 2)]; // 3 < MIN_DECOMP
            let fams = decompose_families(&e, &SplitParams::default());
            assert_eq!(fams.len(), 1);
            assert_eq!(fams[0].members, vec![0, 1, 2]);
            assert_eq!(fams[0].class, FamilyClass::Family);
        }

        #[test]
        fn decompose_large_component_splits() {
            // two K4 cliques bridged by one weak edge; 8 nodes >= MIN_DECOMP -> split into two size-4 families.
            let mut e = vec![];
            for &(a, b) in &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)] {
                e.push(we(a, b));
            }
            for &(a, b) in &[(4, 5), (4, 6), (4, 7), (5, 6), (5, 7), (6, 7)] {
                e.push(we(a, b));
            }
            e.push((0, 4, 0.05));
            let fams = decompose_families(&e, &SplitParams::default());
            assert_eq!(fams.len(), 2);
            let mut sizes: Vec<usize> = fams.iter().map(|f| f.members.len()).collect();
            sizes.sort();
            assert_eq!(sizes, vec![4, 4]);
            assert!(fams.iter().all(|f| f.class == FamilyClass::Family));
        }

        #[test]
        fn decompose_flags_sparse_web() {
            // a star: center 0 with 14 leaves. Louvain keeps it one community; n=15, density ~0.133 -> web.
            let e: Vec<(usize, usize, f64)> = (1..15).map(|leaf| we(0, leaf)).collect();
            let fams = decompose_families(&e, &SplitParams::default());
            assert_eq!(fams.len(), 1);
            assert_eq!(fams[0].members.len(), 15);
            assert!(fams[0].stats.density < 0.15);
            assert_eq!(fams[0].class, FamilyClass::Web);
        }

        #[test]
        fn louvain_resolution_controls_granularity() {
            // K6 clique. Low resolution keeps it whole; high resolution shatters it into singletons -- this
            // directly pins the `gamma` term in the modularity gain.
            let mut e = vec![];
            for a in 0..6 {
                for b in (a + 1)..6 {
                    e.push(we(a, b));
                }
            }
            assert_eq!(
                louvain_communities(6, &e, 0.5).len(),
                1,
                "low resolution keeps the clique whole"
            );
            assert_eq!(
                louvain_communities(6, &e, 2.0).len(),
                6,
                "high resolution shatters K6 into singletons"
            );
        }

        #[test]
        fn louvain_multilevel_aggregation_merges_supernodes() {
            // Hierarchical: two DISCONNECTED groups, each = 4 triangles whose anchor nodes form a strong K4.
            // Level-1 local-moving forms the per-triangle communities; only the AGGREGATION level coalesces the
            // 4 triangle super-nodes of a group into one community. The correct partition is 2 groups of 12.
            let mut e: Vec<(usize, usize, f64)> = vec![];
            for g in 0..2 {
                let base = g * 12;
                let anchors = [base, base + 3, base + 6, base + 9];
                for t in 0..4 {
                    let (x, y, z) = (base + 3 * t, base + 3 * t + 1, base + 3 * t + 2);
                    e.push(we(x, y));
                    e.push(we(y, z));
                    e.push(we(x, z)); // triangle (weight 1.0)
                }
                for i in 0..anchors.len() {
                    for j in (i + 1)..anchors.len() {
                        e.push((anchors[i], anchors[j], 5.0)); // strong inter-triangle K4 among anchors
                    }
                }
            }
            let mut comms = louvain_communities(24, &e, 1.0);
            comms.sort();
            assert_eq!(comms.len(), 2, "two groups");
            assert!(
                comms.iter().all(|c| c.len() == 12),
                "each group is its 12 nodes: {comms:?}"
            );
            assert_eq!(comms[0], (0..12).collect::<Vec<_>>());
            assert_eq!(comms[1], (12..24).collect::<Vec<_>>());
        }

        #[test]
        fn louvain_beats_weakest_edge_cut() {
            // The single WEAKEST edge (0-2, weight 0.1) is INSIDE community A; the correct split is across the
            // bridge (2-3, weight 0.3). A "cut the weakest edge" heuristic would wrongly split A; only a real
            // modularity optimiser recovers {0,1,2},{3,4,5}.
            let e = [
                (0, 1, 1.0),
                (1, 2, 1.0),
                (0, 2, 0.1), // community A: one weak intra edge
                (3, 4, 1.0),
                (4, 5, 1.0),
                (3, 5, 1.0), // community B: strong triangle
                (2, 3, 0.3), // bridge (stronger than the weakest intra edge)
            ];
            let mut comms = louvain_communities(6, &e, 1.0);
            comms.sort();
            assert_eq!(comms, vec![vec![0, 1, 2], vec![3, 4, 5]]);
        }

        #[test]
        fn decompose_drops_singleton_communities() {
            // K7 (strong, weight 10) plus a pendant node 7 weakly attached (0-7, weight 0.05). At resolution
            // 1.1 the clique survives but the pendant detaches into its OWN singleton community, which
            // decompose drops (`len >= 2`), so node 7 vanishes from all families.
            let mut e = vec![];
            for a in 0..7 {
                for b in (a + 1)..7 {
                    e.push((a, b, 10.0));
                }
            }
            e.push((0, 7, 0.05));
            let p = SplitParams {
                resolution: 1.1,
                ..SplitParams::default()
            };
            let fams = decompose_families(&e, &p);
            assert_eq!(fams.len(), 1, "only the clique family survives");
            assert_eq!(fams[0].members, (0..7).collect::<Vec<_>>());
            assert!(
                !fams[0].members.contains(&7),
                "the dropped singleton node is gone"
            );
        }

        // ---- gamma_quasi_clique_partition ----

        #[test]
        fn gamma_quasi_clique_keeps_array_whole_splits_repeat_chain() {
            // A 5-node dense clique (a tandem array) stays ONE block.
            let clique: Vec<(usize, usize, f64)> = {
                let mut e = Vec::new();
                for i in 0..5 {
                    for j in (i + 1)..5 {
                        e.push((i, j, 1.0));
                    }
                }
                e
            };
            let blocks = gamma_quasi_clique_partition(5, &clique, 0.20);
            assert_eq!(blocks.len(), 1, "a dense array is one gamma-quasi-clique");
            assert_eq!(blocks[0].len(), 5);

            // A LONG sparse bridge chain (12-node path: density 2*11/(12*11)=0.167 < gamma=0.20) must split.
            let chain: Vec<(usize, usize, f64)> = (0..11).map(|i| (i, i + 1, 1.0)).collect();
            let blocks = gamma_quasi_clique_partition(12, &chain, 0.20);
            assert!(
                blocks.len() >= 2,
                "a sparse repeat-bridge chain is split, got {:?}",
                blocks
            );
        }

        #[test]
        fn gamma_quasi_clique_partition_preserves_isolated_nodes() {
            // node 2 has no edge (degree 0). The partition must still COVER all of 0..3.
            let blocks = gamma_quasi_clique_partition(3, &[(0, 1, 1.0)], 0.2);
            let mut all: Vec<usize> = blocks.iter().flatten().copied().collect();
            all.sort_unstable();
            assert_eq!(
                all,
                vec![0, 1, 2],
                "partition must cover all of 0..n incl isolated node 2"
            );
            assert!(
                blocks.iter().any(|b| b == &vec![2]),
                "isolated node 2 is its own block, got {blocks:?}"
            );
        }

        #[test]
        fn community_stats_relabels_noncontiguous_ids() {
            // members are NON-contiguous global ids forming a path 10-3-27-5; exercises the global->local
            // relabel the integration layer relies on. Articulation points are the two middle nodes.
            let e = [(10, 3, 0.2), (3, 27, 0.2), (27, 5, 0.2)];
            let s = community_stats(&[10, 3, 27, 5], &e);
            assert_eq!(s.n, 4);
            assert_eq!(s.n_edges, 3);
            assert!((s.density - 0.5).abs() < 1e-9);
            assert_eq!(s.n_articulation, 2);
        }
    }
}

pub mod mosaic {
    //! Gene-conversion mosaic-read detection (audit theme: "VG finds unusual exon combinations").
    //!
    //! A paralog family has near-identical copies that differ at a set of *diagnostic sites*.
    //! For a read, at each diagnostic site it covers we know which copies' expected base it
    //! matches. A GENE-CONVERSION recombinant read's per-site copy pattern SWITCHES from one
    //! copy to another at a contiguous breakpoint (copy A for a run of sites, then copy B) —
    //! an "unusual combination" directly observed in ONE read, not enumerated. This module
    //! detects that switch and rejects the look-alikes: sequencing-error flips, non-identifiable
    //! reads, low-power reads, and (at the family layer) one-off chimeras.
    //!
    //! Design = synthesis of an independent design panel (statistical / algorithmic / biological
    //! lenses). The per-read core `detect_mosaic` is a PURE function over the per-site match
    //! structure; `aggregate_family` confirms a conversion only when the breakpoint RECURS across
    //! independent molecules (the conversion-vs-chimera discriminator). Default-OFF in the
    //! pipeline (RUSTLE_VG_MOSAIC_ON); additive — the existing per-copy EM scoring is untouched.
    //!
    //! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    /// Per-site observation handed to the detector: which copies the read matches at one
    /// diagnostic site. `match_bits[c]` = (read base == copy c's expected base). Sites are in
    /// genomic order.
    #[derive(Debug, Clone, PartialEq)]
    pub struct SiteObs {
        pub ref_pos: u64,
        pub match_bits: Vec<bool>,
    }

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum MosaicStatus {
        LowPower,        // too few decisive sites to call anything
        NonIdentifiable, // copies tie at too many sites — no per-site signal; abstain
        SingleCopy,      // all decisive sites point at one copy
        NoSwitch,        // a switch was scored but failed the gates (error scatter / weak)
        Mosaic,          // a confident contiguous copy switch
    }

    /// Per-read verdict. Abstain/SingleCopy/NoSwitch are all `is_mosaic()==false`; the EM's
    /// existing per-copy scoring is independent of this (additive metadata).
    #[derive(Debug, Clone, PartialEq)]
    pub struct MosaicCall {
        pub status: MosaicStatus,
        pub n_sites: usize,
        pub n_decisive: usize,
        pub copy_a: Option<usize>,              // 5'-proximal tract copy
        pub copy_b: Option<usize>,              // 3'-proximal tract copy
        pub breakpoint_ref: Option<(u64, u64)>, // (last A-site, first B-site) — honest bracket
        pub tract_a_sites: usize,
        pub tract_b_sites: usize,
        pub left_purity: f64,
        pub right_purity: f64,
        pub margin: i32,    // decisive sites the split explains beyond best single
        pub lr_switch: f64, // 2*(two-segment loglik − best single-copy loglik)
        pub threshold_used: f64,
        pub score: f64, // lr_switch − threshold_used (>0 iff Mosaic)
    }

    impl MosaicCall {
        pub fn is_mosaic(&self) -> bool {
            self.status == MosaicStatus::Mosaic
        }
    }

    #[derive(Debug, Clone)]
    pub struct MosaicParams {
        pub min_decisive_sites: usize, // D below this → LowPower
        pub min_tract_sites: usize,    // each flank must have ≥ this many agreeing decisive sites
        pub min_improvement: i32,      // split must explain ≥ this many more decisive sites
        pub min_seg_purity: f64,       // per-tract agreement fraction on both sides
        pub max_ambig_frac: f64,       // > this fraction Ambig → NonIdentifiable
        pub alpha_target: f64,         // target per-read FPR (Bonferroni-corrected internally)
        pub bic_penalty_coeff: f64,    // coeff on ln(S) in the model-complexity penalty
        pub eps_floor: f64,            // lower clamp on per-read error rate
        pub eps_cap: f64,              // upper clamp (and fail-safe when read.de is absent)
        // family aggregation
        pub family_min_supporting_reads: usize,
        pub breakpoint_tol: u64,
        pub max_breakpoint_dispersion: u64,
    }

    impl Default for MosaicParams {
        fn default() -> Self {
            MosaicParams {
                min_decisive_sites: 6,
                min_tract_sites: 3,
                min_improvement: 3,
                min_seg_purity: 0.85,
                max_ambig_frac: 0.5,
                alpha_target: 0.01,
                bic_penalty_coeff: 2.0,
                eps_floor: 0.005,
                eps_cap: 0.05,
                family_min_supporting_reads: 3,
                breakpoint_tol: 50,
                max_breakpoint_dispersion: 50,
            }
        }
    }

    impl MosaicParams {
        /// Build from RUSTLE_VG_MOSAIC_* env overrides, falling back to defaults.
        pub fn from_env() -> Self {
            let mut p = MosaicParams::default();
            let getf = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<f64>().ok());
            let getu = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<usize>().ok());
            if let Some(v) = getu("RUSTLE_VG_MOSAIC_MIN_DECISIVE") {
                p.min_decisive_sites = v;
            }
            if let Some(v) = getu("RUSTLE_VG_MOSAIC_MIN_TRACT") {
                p.min_tract_sites = v;
            }
            if let Some(v) = getu("RUSTLE_VG_MOSAIC_MIN_IMPROVEMENT") {
                p.min_improvement = v as i32;
            }
            if let Some(v) = getf("RUSTLE_VG_MOSAIC_MIN_PURITY") {
                p.min_seg_purity = v;
            }
            if let Some(v) = getf("RUSTLE_VG_MOSAIC_ALPHA") {
                p.alpha_target = v;
            }
            if let Some(v) = getu("RUSTLE_VG_MOSAIC_MIN_READS") {
                p.family_min_supporting_reads = v;
            }
            p
        }
    }

    /// Pure per-read detector. `obs` = ordered per-site match structure; `eps` = per-site error
    /// rate (clamped by the caller). Deterministic, no I/O. With decisive-site agreements the
    /// likelihood-ratio reduces to `2·margin·ln((1−eps)/eps)`, and the χ²(2df) null quantile is
    /// exactly `−2·ln(α)` — so the threshold is closed-form (the spec's bootstrap_reps=0 path).
    pub fn detect_mosaic(
        obs: &[SiteObs],
        n_copies: usize,
        eps: f64,
        p: &MosaicParams,
    ) -> MosaicCall {
        let s = obs.len();
        let mk = |status: MosaicStatus, d: usize| MosaicCall {
            status,
            n_sites: s,
            n_decisive: d,
            copy_a: None,
            copy_b: None,
            breakpoint_ref: None,
            tract_a_sites: 0,
            tract_b_sites: 0,
            left_purity: 0.0,
            right_purity: 0.0,
            margin: 0,
            lr_switch: 0.0,
            threshold_used: 0.0,
            score: 0.0,
        };

        // Tokenize: a DECISIVE site has exactly one matching copy; ≥2 = Ambig (non-identifiable
        // here); 0 = Novel (read base matches no modeled copy) — a wildcard, ignored.
        let mut decisive_copy: Vec<usize> = Vec::with_capacity(s);
        let mut decisive_pos: Vec<u64> = Vec::with_capacity(s);
        let mut n_ambig = 0usize;
        for o in obs {
            let mut matched = usize::MAX;
            let mut n_match = 0usize;
            for c in 0..n_copies {
                if o.match_bits.get(c).copied().unwrap_or(false) {
                    n_match += 1;
                    matched = c;
                }
            }
            if n_match == 1 {
                decisive_copy.push(matched);
                decisive_pos.push(o.ref_pos);
            } else if n_match >= 2 {
                n_ambig += 1;
            }
        }
        let d = decisive_copy.len();

        // Gates (abstain before any switch test). NonIdentifiable takes precedence over
        // LowPower: a read swamped by ties has no per-site signal even if it covers many sites
        // (and an all-ambiguous read has d=0, which would otherwise read as LowPower).
        if s > 0 && (n_ambig as f64) > p.max_ambig_frac * s as f64 {
            return mk(MosaicStatus::NonIdentifiable, d);
        }
        if d < p.min_decisive_sites {
            return mk(MosaicStatus::LowPower, d);
        }
        // Distinct copies among decisive sites, ASCENDING (a dense Vec, not a BTreeSet: avoids a
        // per-read heap allocation + tree walk in the hot scan below; ascending order preserves the
        // old tie-break — `max_by_key` keeps the highest copy id among ties).
        let mut distinct: Vec<usize> = decisive_copy.clone();
        distinct.sort_unstable();
        distinct.dedup();
        if distinct.len() <= 1 {
            return mk(MosaicStatus::SingleCopy, d);
        }

        // Prefix counts make a range agreement count O(1): `prefix[c][i]` = #copy-c decisive sites
        // in `decisive_copy[0..i]`. Replaces the per-(k, copy) linear re-scan
        // (O(d²·copies) → O(d·copies)); one allocation reused across all k.
        let ncap = distinct.last().copied().unwrap_or(0) + 1;
        let mut prefix: Vec<Vec<u32>> = vec![vec![0u32; d + 1]; ncap];
        for i in 0..d {
            for c in 0..ncap {
                prefix[c][i + 1] = prefix[c][i];
            }
            prefix[decisive_copy[i]][i + 1] += 1;
        }
        let agree_for =
            |c: usize, lo: usize, hi: usize| -> usize { (prefix[c][hi] - prefix[c][lo]) as usize };
        let best_single = distinct
            .iter()
            .map(|&c| agree_for(c, 0, d))
            .max()
            .unwrap_or(0);

        // Exhaustive single-changepoint scan over decisive-site split indices, maximizing the
        // number of decisive sites explained by (copy a left of k, copy b right of k), a≠b.
        let mut best_split = best_single;
        let mut best_kab: Option<(usize, usize, usize)> = None;
        for k in 1..d {
            let a = *distinct
                .iter()
                .max_by_key(|&&c| agree_for(c, 0, k))
                .unwrap();
            let b = *distinct
                .iter()
                .max_by_key(|&&c| agree_for(c, k, d))
                .unwrap();
            if a == b {
                continue;
            }
            let explained = agree_for(a, 0, k) + agree_for(b, k, d);
            // First k to reach a new best (or to match the single-copy best) wins the tie.
            if explained > best_split || (explained == best_split && best_kab.is_none()) {
                best_split = explained;
                best_kab = Some((k, a, b));
            }
        }

        let (k, a, b) = match best_kab {
            Some(x) => x,
            None => return mk(MosaicStatus::SingleCopy, d), // no a≠b split improves on single copy
        };

        let tract_a_sites = agree_for(a, 0, k);
        let tract_b_sites = agree_for(b, k, d);
        let left_purity = tract_a_sites as f64 / k as f64;
        let right_purity = tract_b_sites as f64 / (d - k) as f64;
        let margin = best_split as i32 - best_single as i32;

        // Likelihood-ratio over decisive sites: 2·margin·ln((1−eps)/eps).
        let lr = 2.0 * margin as f64 * ((1.0 - eps) / eps).ln();
        // Threshold = BIC complexity penalty + χ²(2df) quantile at the Bonferroni-corrected α.
        let n_pairs = (n_copies * n_copies.saturating_sub(1) / 2).max(1) as f64;
        let alpha_corr = (p.alpha_target / (s as f64 * n_pairs)).max(1e-12);
        let threshold = p.bic_penalty_coeff * (s as f64).ln() + (-2.0 * alpha_corr.ln());

        let mut call = mk(MosaicStatus::NoSwitch, d);
        call.copy_a = Some(a);
        call.copy_b = Some(b);
        call.breakpoint_ref = Some((decisive_pos[k - 1], decisive_pos[k]));
        call.tract_a_sites = tract_a_sites;
        call.tract_b_sites = tract_b_sites;
        call.left_purity = left_purity;
        call.right_purity = right_purity;
        call.margin = margin;
        call.lr_switch = lr;
        call.threshold_used = threshold;
        call.score = lr - threshold;

        // Dual gate: calibrated LR AND hard integer/run-length backstops (the backstops hold
        // even when the i.i.d. error null is violated by correlated/homopolymer error bursts).
        let pass = a != b
            && tract_a_sites >= p.min_tract_sites
            && tract_b_sites >= p.min_tract_sites
            && margin >= p.min_improvement
            && left_purity >= p.min_seg_purity
            && right_purity >= p.min_seg_purity
            && lr > threshold;
        if pass {
            call.status = MosaicStatus::Mosaic;
        }
        call
    }

    /// A family-level confirmed (or suspected) gene-conversion event.
    #[derive(Debug, Clone, PartialEq)]
    pub struct ConversionEvent {
        pub copy_a: usize,
        pub copy_b: usize,
        pub chrom: String, // chromosome the breakpoint coordinates live on (for the microhomology check)
        pub breakpoint_ref: (u64, u64), // consensus bracket (min last-A, max first-B)
        pub n_supporting_reads: usize,
        pub breakpoint_dispersion: u64, // spread of per-read breakpoint midpoints
        pub confirmed: bool,            // false = ChimeraSuspect
    }

    /// Family aggregation: a genuine conversion RECURS at a fixed breakpoint across independent
    /// molecules; a one-off chimera does not. Inputs are per-read Mosaic calls (caller dedupes
    /// to distinct molecules first). Clusters by oriented (a→b) pair and breakpoint midpoint
    /// within `breakpoint_tol`; confirms a cluster with ≥ `family_min_supporting_reads` molecules
    /// and tight dispersion.
    pub fn aggregate_family(
        calls: &[MosaicCall],
        chroms: &[&str],
        p: &MosaicParams,
    ) -> Vec<ConversionEvent> {
        // `chroms[i]` is the chromosome `calls[i]`'s breakpoint coordinates live on (parallel array). The
        // chrom is part of the cluster KEY: breakpoints at coincidentally-similar positions on DIFFERENT
        // chromosomes (multi-chrom paralog families, e.g. RABL2A/RABL2B) must not cluster together, and the
        // emitted event needs its chrom for the downstream microhomology check.
        let mut mosaics: Vec<(usize, usize, &str, u64, (u64, u64))> = calls
            .iter()
            .enumerate()
            .filter(|(_, c)| c.is_mosaic())
            .filter_map(|(i, c)| match (c.copy_a, c.copy_b, c.breakpoint_ref) {
                (Some(a), Some(b), Some(br)) => Some((
                    a,
                    b,
                    chroms.get(i).copied().unwrap_or(""),
                    (br.0 + br.1) / 2,
                    br,
                )),
                _ => None,
            })
            .collect();
        // Stable order by (a, b, chrom, midpoint).
        mosaics.sort_by(|x, y| (x.0, x.1, x.2, x.3).cmp(&(y.0, y.1, y.2, y.3)));

        let mut events: Vec<ConversionEvent> = Vec::new();
        let mut i = 0usize;
        while i < mosaics.len() {
            let (a, b, chrom, _, _) = mosaics[i];
            // Greedily grow a cluster of same oriented pair ON THE SAME CHROM within breakpoint_tol.
            let mut j = i;
            let mut mids: Vec<u64> = Vec::new();
            let mut br_lo = u64::MAX;
            let mut br_hi = 0u64;
            while j < mosaics.len() {
                let (aj, bj, chromj, midj, brj) = mosaics[j];
                if aj != a || bj != b || chromj != chrom {
                    break;
                }
                if let Some(&last) = mids.last() {
                    if midj.saturating_sub(last) > p.breakpoint_tol {
                        break;
                    }
                }
                mids.push(midj);
                br_lo = br_lo.min(brj.0);
                br_hi = br_hi.max(brj.1);
                j += 1;
            }
            let n = mids.len();
            let dispersion = mids.last().copied().unwrap_or(0).saturating_sub(mids[0]);
            let confirmed =
                n >= p.family_min_supporting_reads && dispersion <= p.max_breakpoint_dispersion;
            events.push(ConversionEvent {
                copy_a: a,
                copy_b: b,
                chrom: chrom.to_string(),
                breakpoint_ref: (br_lo, br_hi),
                n_supporting_reads: n,
                breakpoint_dispersion: dispersion,
                confirmed,
            });
            i = j;
        }
        events
    }

    /// Unified gene-conversion-vs-artifact verdict for a mosaic family event. `aggregate_family`'s
    /// `confirmed` flag captures only ONE leg (recurrence across molecules); but recurrence alone is
    /// insufficient — a sequence-driven template-switch hotspot (microhomology at the same point)
    /// produces RT-switch chimeras that ALSO recur, so they pass the recurrence gate. The discriminator
    /// therefore needs two ORTHOGONAL legs in addition to recurrence:
    ///   * **microhomology** at the breakpoint = the RT/template-switch signature (a direct repeat
    ///     flanking the switch point — `genome::is_rt_switch` applied to the breakpoint bracket);
    ///   * **DNA support** = heritability: a real (historical) gene conversion is in the genome and so
    ///     recurs in matched DNA reads; an RT/template switch is an RNA-library artifact, absent from DNA.
    /// Both legs are passed in as `Option<bool>` so "no evidence available" (`None`) is distinct from
    /// "negative evidence" (`Some(false)`). The DNA leg can act as a **veto** (`Some(false)` → `Ambiguous`)
    /// when a RELIABLE absence source exists; this lets the two cheap legs (recurrence + microhomology)
    /// ship without the DNA catalog wired, while a reliable DNA source strengthens or vetoes the call.
    ///
    /// MEASURED (bench/mosaic_discriminator/dna_support.py, T2T DNA PSV catalog): the catalog signal is
    /// real (42% of multi-copy families show a heritable-conversion DNA mosaic) but SPARSE and ref0-centric
    /// (only ~2.9% of the genome is in a ref0 interval; localized mosaics return "absent" almost everywhere
    /// even in families that HAVE one). So catalog "absent" is UNRELIABLE negative evidence — a catalog-
    /// backed DNA closure must return `Some(true)` / `None` only (positive corroboration), NEVER
    /// `Some(false)`, or it would wrongly downgrade real conversions. The production paths pass `None`.
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum Classification {
        /// recurrent + present in DNA + no template-switch signature → real biological gene conversion.
        GeneConversion,
        /// direct-repeat/microhomology at the breakpoint and NOT DNA-confirmed → RT/template-switch artifact.
        RtSwitchArtifact,
        /// sporadic (did not recur) with no template-switch signature → one-off chimera, lean artifact.
        ChimeraSuspect,
        /// conflicting or insufficient evidence (e.g. microhomology AND DNA support; or recurrent but DNA unknown).
        Ambiguous,
    }

    /// Classify one family event from the three orthogonal legs (recurrence via `ev.confirmed`,
    /// `microhomology`, `dna_supported`). Pure: the caller supplies the two genome/DNA-derived signals.
    /// Evaluation order matters — the microhomology-artifact rule fires BEFORE the gene-conversion rule,
    /// so a recurrent-but-RT-signature event is correctly called an artifact rather than a conversion.
    pub fn classify_event(
        ev: &ConversionEvent,
        microhomology: Option<bool>,
        dna_supported: Option<bool>,
    ) -> Classification {
        let mh = microhomology == Some(true);
        let dna_present = dna_supported == Some(true);
        let dna_absent = dna_supported == Some(false);
        if mh && !dna_present {
            // template-switch signature, not rescued by positive DNA support → artifact (even if it recurs).
            Classification::RtSwitchArtifact
        } else if ev.confirmed && !mh && !dna_absent {
            // recurrent + no template signature + DNA not contradicting (present or unchecked) → conversion.
            Classification::GeneConversion
        } else if !ev.confirmed && !mh {
            // sporadic, no signature → one-off chimera.
            Classification::ChimeraSuspect
        } else {
            // microhomology∧DNA-present conflict, or recurrent-but-DNA-absent (heritability contradicted).
            Classification::Ambiguous
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        // Build per-site obs from a label string where each char is the copy index that the read
        // uniquely matches; '*' = Ambig (all copies match), '.' = Novel (no copy matches).
        fn obs_from(labels: &str, n_copies: usize) -> Vec<SiteObs> {
            labels
                .chars()
                .enumerate()
                .map(|(i, ch)| {
                    let mut match_bits = vec![false; n_copies];
                    match ch {
                        '*' => match_bits.iter_mut().for_each(|m| *m = true),
                        '.' => {}
                        c => {
                            let idx = c.to_digit(10).unwrap() as usize;
                            match_bits[idx] = true;
                        }
                    }
                    SiteObs {
                        ref_pos: 1000 + i as u64 * 10,
                        match_bits,
                    }
                })
                .collect()
        }

        fn p() -> MosaicParams {
            MosaicParams::default()
        }

        #[test]
        fn clean_switch_is_mosaic() {
            let c = detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p());
            assert_eq!(c.status, MosaicStatus::Mosaic);
            assert_eq!((c.copy_a, c.copy_b), (Some(0), Some(1)));
            assert_eq!((c.tract_a_sites, c.tract_b_sites), (3, 3));
            assert_eq!(c.margin, 3);
            // breakpoint bracket between the last 0-site and the first 1-site.
            assert_eq!(c.breakpoint_ref, Some((1000 + 2 * 10, 1000 + 3 * 10)));
        }

        #[test]
        fn reversed_switch_records_direction() {
            let c = detect_mosaic(&obs_from("111000", 2), 2, 0.005, &p());
            assert!(c.is_mosaic());
            assert_eq!((c.copy_a, c.copy_b), (Some(1), Some(0)));
        }

        #[test]
        fn isolated_error_flip_is_not_mosaic() {
            // D=8, one interior flip; no a≠b contiguous split improves → margin 0.
            let c = detect_mosaic(&obs_from("00000100", 2), 2, 0.005, &p());
            assert!(!c.is_mosaic());
            assert!(c.margin < p().min_improvement);
        }

        #[test]
        fn pure_copy_is_single_copy() {
            let c = detect_mosaic(&obs_from("000000", 2), 2, 0.005, &p());
            assert_eq!(c.status, MosaicStatus::SingleCopy);
            assert!(!c.is_mosaic());
        }

        #[test]
        fn all_ambiguous_abstains_nonidentifiable() {
            let c = detect_mosaic(&obs_from("******", 2), 2, 0.005, &p());
            assert_eq!(c.status, MosaicStatus::NonIdentifiable);
        }

        #[test]
        fn too_few_sites_abstains_lowpower() {
            let c = detect_mosaic(&obs_from("0011", 2), 2, 0.005, &p());
            assert_eq!(c.status, MosaicStatus::LowPower);
        }

        #[test]
        fn short_right_tract_fails_min_tract() {
            // D=7, right tract only 2 sites (< min_tract_sites=3) → not Mosaic.
            let c = detect_mosaic(&obs_from("0000011", 2), 2, 0.005, &p());
            assert!(!c.is_mosaic());
            assert!(c.tract_b_sites < p().min_tract_sites);
        }

        #[test]
        fn novel_wildcard_does_not_break_tracts() {
            // '.' (Novel) between clean 3+3 tracts is neutral.
            let c = detect_mosaic(&obs_from("000.111", 2), 2, 0.005, &p());
            assert!(c.is_mosaic());
            assert_eq!((c.tract_a_sites, c.tract_b_sites), (3, 3));
        }

        #[test]
        fn ambiguous_sites_are_neutral_in_tracts() {
            let c = detect_mosaic(&obs_from("00*00111*1", 2), 2, 0.005, &p());
            assert!(c.is_mosaic());
            assert_eq!((c.copy_a, c.copy_b), (Some(0), Some(1)));
        }

        #[test]
        fn three_copy_family_picks_discriminating_pair() {
            // copies 0 and 2 are the switching pair in a 3-copy family.
            let c = detect_mosaic(&obs_from("000222", 3), 3, 0.005, &p());
            assert!(c.is_mosaic());
            assert_eq!((c.copy_a, c.copy_b), (Some(0), Some(2)));
        }

        #[test]
        fn de_none_uses_wider_eps_is_fail_safe() {
            // A borderline call must clear a HARDER bar at the larger (fail-safe) eps.
            let lo = detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p());
            let hi = detect_mosaic(&obs_from("000111", 2), 2, 0.05, &p());
            assert!(lo.score > hi.score); // larger eps → smaller LR → smaller margin-over-threshold
        }

        #[test]
        fn family_confirms_reproducible_breakpoint() {
            // 3 independent reads with the SAME switch/breakpoint → Confirmed.
            let calls: Vec<MosaicCall> = (0..3)
                .map(|_| detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p()))
                .collect();
            let events = aggregate_family(&calls, &vec!["c1"; calls.len()], &p());
            assert_eq!(events.len(), 1);
            assert!(events[0].confirmed);
            assert_eq!(events[0].n_supporting_reads, 3);
        }

        #[test]
        fn family_event_carries_its_chrom() {
            let calls: Vec<MosaicCall> = (0..3)
                .map(|_| detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p()))
                .collect();
            let events = aggregate_family(&calls, &vec!["chrX"; calls.len()], &p());
            assert_eq!(events.len(), 1);
            assert_eq!(events[0].chrom, "chrX");
        }

        #[test]
        fn same_breakpoint_on_different_chroms_does_not_cluster() {
            // 3 reads on chrom A + 3 reads on chrom B, all with the IDENTICAL switch/midpoint. They must
            // form TWO separate events (one per chrom), not one merged cluster — the multi-chrom paralog
            // (e.g. RABL2A/RABL2B) case. Chrom is part of the cluster key.
            let calls: Vec<MosaicCall> = (0..6)
                .map(|_| detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p()))
                .collect();
            let chroms = ["cA", "cA", "cA", "cB", "cB", "cB"];
            let events = aggregate_family(&calls, &chroms, &p());
            assert_eq!(
                events.len(),
                2,
                "different chroms must not cluster together"
            );
            let mut cs: Vec<&str> = events.iter().map(|e| e.chrom.as_str()).collect();
            cs.sort();
            assert_eq!(cs, vec!["cA", "cB"]);
            assert!(events
                .iter()
                .all(|e| e.confirmed && e.n_supporting_reads == 3));
        }

        #[test]
        fn family_rejects_singleton_as_chimera_suspect() {
            let calls = vec![detect_mosaic(&obs_from("000111", 2), 2, 0.005, &p())];
            let events = aggregate_family(&calls, &vec!["c1"; calls.len()], &p());
            assert_eq!(events.len(), 1);
            assert!(!events[0].confirmed); // 1 molecule < family_min_supporting_reads
        }

        // ----- classify_event (unified gene-conversion-vs-artifact discriminator) -----

        fn ev(confirmed: bool) -> ConversionEvent {
            ConversionEvent {
                copy_a: 0,
                copy_b: 1,
                chrom: "c1".to_string(),
                breakpoint_ref: (1000, 1010),
                n_supporting_reads: if confirmed { 5 } else { 1 },
                breakpoint_dispersion: 0,
                confirmed,
            }
        }

        #[test]
        fn classify_recurrent_dna_no_microhomology_is_gene_conversion() {
            assert_eq!(
                classify_event(&ev(true), Some(false), Some(true)),
                Classification::GeneConversion
            );
        }

        #[test]
        fn classify_microhomology_without_dna_is_rt_switch_even_if_recurrent() {
            // the load-bearing case: recurrence ALONE would have called this confirmed, but the
            // template-switch signature (microhomology) + no DNA support overrides it to artifact.
            assert_eq!(
                classify_event(&ev(true), Some(true), Some(false)),
                Classification::RtSwitchArtifact
            );
            assert_eq!(
                classify_event(&ev(true), Some(true), None),
                Classification::RtSwitchArtifact
            );
        }

        #[test]
        fn classify_sporadic_no_signature_is_chimera_suspect() {
            assert_eq!(
                classify_event(&ev(false), Some(false), None),
                Classification::ChimeraSuspect
            );
        }

        #[test]
        fn classify_microhomology_and_dna_conflict_is_ambiguous() {
            // direct repeat AND present in DNA: could be a real conversion at a repeat-prone site — abstain.
            assert_eq!(
                classify_event(&ev(true), Some(true), Some(true)),
                Classification::Ambiguous
            );
        }

        #[test]
        fn classify_recurrent_no_microhomology_dna_unchecked_is_gene_conversion() {
            // DNA is a VETO, not a requirement: unchecked (None) does not block the two cheap legs.
            assert_eq!(
                classify_event(&ev(true), Some(false), None),
                Classification::GeneConversion
            );
        }

        #[test]
        fn classify_recurrent_but_dna_absent_is_ambiguous() {
            // DNA was CHECKED and the breakpoint is ABSENT from the genome → contradicts heritability.
            assert_eq!(
                classify_event(&ev(true), Some(false), Some(false)),
                Classification::Ambiguous
            );
        }

        // ----- ground-truth confusion matrix: full real path (detect -> aggregate -> genome
        //       microhomology -> classify), BOTH directions, over a constructed genome -----

        /// Build N recurrent recombinant reads that switch copy 0 -> copy 1 at site index `bp`, with
        /// their breakpoint bracket placed at genome positions `(left, right)` (so the genome's sequence
        /// at those coords decides microhomology). Returns the aggregated (single) confirmed event.
        fn recurrent_event_at(bp_left: u64, bp_right: u64, n_reads: usize) -> ConversionEvent {
            // 6 decisive sites in ASCENDING genomic order: three copy-0 sites ending exactly at `bp_left`,
            // then three copy-1 sites starting exactly at `bp_right` = a clean "000111" 0->1 switch whose
            // breakpoint bracket is (bp_left, bp_right).
            let positions = [
                (bp_left - 20, 0usize),
                (bp_left - 10, 0),
                (bp_left, 0),
                (bp_right, 1),
                (bp_right + 10, 1),
                (bp_right + 20, 1),
            ];
            let mut calls = Vec::new();
            for _ in 0..n_reads {
                let obs: Vec<SiteObs> = positions
                    .iter()
                    .map(|&(ref_pos, copy)| {
                        let mut mb = vec![false; 2];
                        mb[copy] = true;
                        SiteObs {
                            ref_pos,
                            match_bits: mb,
                        }
                    })
                    .collect();
                calls.push(detect_mosaic(&obs, 2, 0.005, &p()));
            }
            let mut events = aggregate_family(&calls, &vec!["c1"; calls.len()], &p());
            assert_eq!(events.len(), 1, "one oriented breakpoint cluster");
            let ev = events.pop().unwrap();
            assert!(
                ev.confirmed,
                "recurrent across {} molecules -> confirmed",
                n_reads
            );
            ev
        }

        #[test]
        fn ground_truth_conversion_vs_rt_switch_confusion_matrix() {
            use crate::genome::GenomeIndex;
            // One genome, two breakpoint loci:
            //  * CONVERSION locus: breakpoint flanks DIFFER -> no microhomology -> GeneConversion.
            //  * RT-SWITCH locus: breakpoint sits at an exact direct repeat -> microhomology -> artifact.
            let mut seq = vec![b'A'; 400];
            // RT-switch direct repeat: 8 bp ending at 200 == 8 bp ending at 240.
            seq[192..200].copy_from_slice(b"CGTACGTA");
            seq[232..240].copy_from_slice(b"CGTACGTA");
            // Conversion locus: flanks ending at 300 vs 340 deliberately DIFFER.
            seq[292..300].copy_from_slice(b"CGTACGTA");
            seq[332..340].copy_from_slice(b"TTGGAACC");
            let g = GenomeIndex::from_seqs(&[("c1", &seq[..])]);

            let mh = |left: u64, right: u64| g.breakpoint_microhomology("c1", left, right, 6, 12);

            // RT-switch breakpoint bracket (200, 240): direct repeat present.
            let rt_ev = recurrent_event_at(200, 240, 5);
            assert!(mh(200, 240), "direct repeat at the RT-switch breakpoint");
            assert_eq!(
                classify_event(&rt_ev, Some(mh(200, 240)), None),
                Classification::RtSwitchArtifact,
                "recurrent + microhomology + DNA-unchecked -> artifact (recurrence alone would have mis-called it)"
            );

            // Conversion breakpoint bracket (300, 340): no direct repeat.
            let gc_ev = recurrent_event_at(300, 340, 5);
            assert!(
                !mh(300, 340),
                "no direct repeat at the conversion breakpoint"
            );
            assert_eq!(
                classify_event(&gc_ev, Some(mh(300, 340)), None),
                Classification::GeneConversion,
                "recurrent + no microhomology + DNA-not-contradicting -> gene conversion"
            );

            // DNA leg as a veto: same clean conversion, but ABSENT from DNA -> downgraded to Ambiguous.
            assert_eq!(
                classify_event(&gc_ev, Some(mh(300, 340)), Some(false)),
                Classification::Ambiguous,
            );
        }
    }
}

pub mod read_conflict {
    //! Read-conflict (mutual-mappability) graph — the operational SCOPE for copy-resolution and the SEED for the
    //! family catalog. NOT a rival definition of a family (`bench/family_def_readconflict.md`).
    //!
    //! DEFINITION vs SCOPE. The family DEFINITION is homology: a γ-quasi-clique component of the transcribed-
    //! homology graph over ≥ 2 loci (`family_definition.rs`, the E_r oracle; `docs/RETIREMENT_AND_MIGRATION.md`,
    //! `DEFINITIONS_FORMAL.md`). This module builds a DIFFERENT graph — over the SAME loci but with a different
    //! edge — that answers the question copy-assignment actually cares about: **do reads cross-map between these
    //! loci?** Two loci are linked iff some read has a placement in BOTH with TIED alignment scores (a genuine
    //! alternative placement — the multimapping conflict). Homology says which loci ARE one family; read-conflict
    //! says how many COPIES they hide and how the assignment decomposes. So it is the SCOPE (which loci must be
    //! co-resolved) and the SEED for the catalog, not the definition of membership.
    //!
    //! Why this is the right SCOPE: (1) no tuned similarity threshold — the boundary is the alignment-score tie
    //! (with `RUSTLE_CONFLICT_SIG`, the SAME significance level α the assignment gate uses), a property of the
    //! data; (2) reads never cross-map outside their component, so the assignment problem decomposes EXACTLY
    //! across families with no information lost; (3) it never groups domain-sharers (a read over a shared exon
    //! maps to one locus — no alternative placement), and it picks out exactly the families where assignment is
    //! needed (validated on Compara labels: 0 conflict on 7/7 domain-sharers, fires on RABL2/APOBEC3, silent on
    //! the resolvable RFPL).
    //!
    //! This is the portable KERNEL: `conflict_edges` (read placements → weighted edges) + `conflict_families`
    //! (edges → connected-component families). The remaining integration is plumbing per-locus secondary
    //! placements (`secondary_index` / `tied_secondary_reads_in_region`) into the detection stage.
    //!
    //! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    /// One read's placement on a candidate locus: the locus index plus the signals the conflict criterion uses —
    /// `de` (gap-compressed divergence, the tie discriminant), `mapq` (both-0 = genuine-multimapper corroboration),
    /// and `as_score` (kept only to log the `de-tie ⊆ AS-tie` invariant).
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct Placement {
        pub locus: usize,
        pub de: f32,
        pub mapq: u8,
        pub as_score: i32,
        /// Aligned-block length (number of M/=/X columns) of this placement. Used only by the SIGNIFICANCE
        /// de-tie criterion to convert the divergence rate `de` into a mismatch COUNT (`m = de * aln_len`);
        /// the legacy delta-based `de_tied` ignores it, so existing behaviour is unaffected when sig is off.
        pub aln_len: u32,
    }

    /// One read's placements over the family's candidate loci. Built from the BAM by the adapter.
    pub type ReadPlacements = Vec<Placement>;

    /// Tunables for the read-conflict criterion (de-tie). Defaults are the bake-off operating point
    /// (`bench/family_criterion_bakeoff.md`); env-overridable via `RUSTLE_CONFLICT_DE_DELTA/DE_MAX/MIN_READS`.
    #[derive(Clone, Copy, Debug)]
    pub struct ConflictParams {
        /// Two placements conflict iff their divergences are within `delta` AND both `<= de_max` (both fit).
        pub delta: f64,
        pub de_max: f64,
        /// Minimum conflicting-read count for an edge (guards the noise floor).
        pub min_reads: usize,
        /// SIGNIFICANCE de-tie (the UNIFICATION with the copy-assignment gate): when `Some((eps, alpha))`, a read
        /// de-ties between two loci iff it CANNOT significantly distinguish them under the SAME IsoCon real-vs-error
        /// test the assignment gate uses — `eps^delta >= alpha`, where `delta = |m_a - m_b|` is the excess mismatch
        /// count (the per-read distinguishing-column proxy, mirroring Theorem 4's `min_p = eps^delta`). The
        /// `de_max` quality floor still applies. `None` (default) = the legacy fixed-`delta` `de_tied`, so OFF is
        /// byte-identical. Replaces the arbitrary `delta=0.005` with the error-model-derived tie threshold.
        pub sig: Option<(f64, f64)>,
    }

    impl Default for ConflictParams {
        fn default() -> Self {
            ConflictParams {
                delta: 0.005,
                de_max: 0.05,
                min_reads: 3,
                sig: None,
            }
        }
    }

    impl ConflictParams {
        /// Read overrides from `RUSTLE_CONFLICT_DE_DELTA`, `RUSTLE_CONFLICT_DE_MAX`, `RUSTLE_CONFLICT_MIN_READS`.
        pub fn from_env() -> Self {
            let d = Self::default();
            let f = |k: &str, v: f64| {
                std::env::var(k)
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(v)
            };
            let u = |k: &str, v: usize| {
                std::env::var(k)
                    .ok()
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(v)
            };
            // Significance de-tie is DEFAULT ON: the conflict edge uses the SAME IsoCon real-vs-error test (and the
            // SAME level `alpha`) as the assignment gate, so no hand-set score-gap `delta` decides the conflict
            // scope. eps = per-distinguishing-column error proxy (e/3 ~ 0.001 HiFi), alpha = significance.
            // RUSTLE_CONFLICT_SIG=0 reverts to the legacy `delta` tie-width (A/B comparison / legacy reproduction).
            let sig = if std::env::var("RUSTLE_CONFLICT_SIG").ok().as_deref() == Some("0") {
                None
            } else {
                Some((
                    f("RUSTLE_CONFLICT_EPS", 0.001),
                    f("RUSTLE_CONFLICT_ALPHA", 1e-3),
                ))
            };
            ConflictParams {
                delta: f("RUSTLE_CONFLICT_DE_DELTA", d.delta),
                de_max: f("RUSTLE_CONFLICT_DE_MAX", d.de_max),
                min_reads: u("RUSTLE_CONFLICT_MIN_READS", d.min_reads),
                sig,
            }
        }
    }

    /// de-tie: `|de_a − de_b| <= delta` AND `max(de_a, de_b) <= de_max` (the read fits both copies, comparably).
    fn de_tied(a: &Placement, b: &Placement, p: &ConflictParams) -> bool {
        let (da, db) = (a.de as f64, b.de as f64);
        (da - db).abs() <= p.delta && da.max(db) <= p.de_max
    }

    /// SIGNIFICANCE de-tie — the unification with the copy-assignment gate. A read counts as conflict evidence
    /// between two loci iff it CANNOT significantly distinguish them: with `m_x = de_x * aln_len_x` the mismatch
    /// count to locus `x`, the excess `delta = |m_a - m_b|` is the per-read distinguishing-column proxy, and the
    /// read is tied iff `eps^delta >= alpha` — exactly the assignment gate's `min_p >= alpha` (Theorem 4). The
    /// `de_max` quality floor still applies (both alignments must genuinely fit). No arbitrary `delta` constant.
    fn sig_tied(a: &Placement, b: &Placement, de_max: f64, eps: f64, alpha: f64) -> bool {
        let (da, db) = (a.de as f64, b.de as f64);
        if da.max(db) > de_max {
            return false;
        }
        let ma = da * a.aln_len as f64;
        let mb = db * b.aln_len as f64;
        let delta_cols = (ma - mb).abs();
        eps.powf(delta_cols) >= alpha
    }

    /// Whether a read's two placements conflict (de-tie), under either criterion (significance if `p.sig` is set).
    fn tied(a: &Placement, b: &Placement, p: &ConflictParams) -> bool {
        match p.sig {
            Some((eps, alpha)) => sig_tied(a, b, p.de_max, eps, alpha),
            None => de_tied(a, b, p),
        }
    }

    /// A read's alignment-score evidence: the best and runner-up `AS:i` over its placements, both raw and
    /// normalized by aligned length.
    ///
    /// AS is REPORTED, never decisive — `de` decides. Raw AS is an absolute score that grows with aligned
    /// length, so a genuine multimapper whose second placement is a PARTIAL alignment scores far lower there
    /// even at identical per-base quality. Measured on GGO Iso-Seq: median runner-up/best AS ratio 0.713 while
    /// the aligned-length ratio is 0.897. That length confound is why `de` (a rate) replaced AS (a total), and
    /// why `de-tie ⊆ AS-tie` does NOT hold on real data. `per_base` divides out the confound and is the value
    /// to compare across placements.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub struct AsEvidence {
        pub best: i32,
        pub second: Option<i32>,
        pub best_per_base: f32,
        pub second_per_base: Option<f32>,
    }

    impl AsEvidence {
        /// `best - second`; `None` when the read has a single placement (nothing to be ambiguous between).
        pub fn margin(&self) -> Option<i32> {
            self.second.map(|s| self.best - s)
        }
    }

    /// Best and runner-up alignment score over one read's placements, as `(as_score, aligned_len)` pairs.
    /// Ranks on RAW `AS` (that is the quantity people quote and Eichler's `AS >= 10` rule uses); the per-base
    /// values are carried alongside for the length-fair comparison. `None` if the read has no placements.
    /// Zero-length placements get a per-base score of 0.0 rather than a division by zero.
    pub fn as_evidence(placements: &[(i32, u32)]) -> Option<AsEvidence> {
        let per_base = |(s, l): (i32, u32)| if l == 0 { 0.0 } else { s as f32 / l as f32 };
        let mut sorted: Vec<(i32, u32)> = placements.to_vec();
        sorted.sort_by(|a, b| b.0.cmp(&a.0));
        let &first = sorted.first()?;
        let second = sorted.get(1).copied();
        Some(AsEvidence {
            best: first.0,
            second: second.map(|s| s.0),
            best_per_base: per_base(first),
            second_per_base: second.map(per_base),
        })
    }

    /// `min(a,b) >= as_tie * max(a,b)` — the legacy AS-tie predicate, kept only for the audit edge-set.
    fn as_tied(a: i32, b: i32, as_tie: f64) -> bool {
        let (hi, lo) = (a.max(b), a.min(b));
        hi > 0 && (lo as f64) >= as_tie * (hi as f64)
    }

    /// Build the read-conflict edges over `n_loci` under the **de-tie** criterion. For each read, every pair of its
    /// placements that de-ties contributes one conflict observation to that locus pair. Returns `(i, j, weight)`
    /// with `i < j` for pairs whose count `>= p.min_reads`, sorted. Self-pairs (same locus) ignored.
    pub fn conflict_edges(
        n_loci: usize,
        reads: &[ReadPlacements],
        p: &ConflictParams,
    ) -> Vec<(usize, usize, usize)> {
        use std::collections::BTreeMap;
        let mut weight: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        for placements in reads {
            for a in 0..placements.len() {
                for b in (a + 1)..placements.len() {
                    let (pa, pb) = (&placements[a], &placements[b]);
                    if pa.locus == pb.locus || pa.locus >= n_loci || pb.locus >= n_loci {
                        continue;
                    }
                    if tied(pa, pb, p) {
                        let key = (pa.locus.min(pb.locus), pa.locus.max(pb.locus));
                        *weight.entry(key).or_insert(0) += 1;
                    }
                }
            }
        }
        weight
            .into_iter()
            .filter(|&(_, w)| w >= p.min_reads)
            .map(|((i, j), w)| (i, j, w))
            .collect()
    }

    /// AS-tie edge node-pairs over `n_loci` — the legacy criterion, kept only as a logged comparison edge-set.
    ///
    /// ⚠ `de-tie ⊆ AS-tie` holds only when both placements are FULL-LENGTH (the unit tests and the planted sims).
    /// On real Iso-Seq it is FALSE: secondary placements are partial, raw AS scales with aligned length, so a
    /// de-tied pair can fail `as_tied` outright. Every real GGO region logs `de⊆AS=false`. Do not treat the
    /// audit line as a regression invariant on real data — it is a diagnostic.
    pub fn as_tie_edges(
        n_loci: usize,
        reads: &[ReadPlacements],
        as_tie: f64,
        min_reads: usize,
    ) -> std::collections::BTreeSet<(usize, usize)> {
        use std::collections::BTreeMap;
        let mut weight: BTreeMap<(usize, usize), usize> = BTreeMap::new();
        for placements in reads {
            for a in 0..placements.len() {
                for b in (a + 1)..placements.len() {
                    let (pa, pb) = (&placements[a], &placements[b]);
                    if pa.locus == pb.locus || pa.locus >= n_loci || pb.locus >= n_loci {
                        continue;
                    }
                    if as_tied(pa.as_score, pb.as_score, as_tie) {
                        *weight
                            .entry((pa.locus.min(pb.locus), pa.locus.max(pb.locus)))
                            .or_insert(0) += 1;
                    }
                }
            }
        }
        weight
            .into_iter()
            .filter(|&(_, w)| w >= min_reads)
            .map(|((i, j), _)| (i, j))
            .collect()
    }

    /// Count reads supporting a conflict family and how many of those reads have BOTH placements at mapq==0
    /// (the genuine-multimapper corroboration). A read is counted when it contributes at least one de-tied pair
    /// whose two loci are BOTH in `family`; of those, a `both_mapq0` read has mapq==0 on BOTH placements in
    /// that pair. Returns `(supporting_reads, both_mapq0_reads)`. Log-only: does NOT gate any edge.
    ///
    /// For a read with placements on >=3 family loci the mapq0 check is applied to the FIRST qualifying pair by
    /// iteration order, so the `both_mapq0` count is CONSERVATIVE for multi-copy families (it can undercount but
    /// never overstate multimapper evidence) — exact for the common 2-locus family.
    pub fn family_mapq0_support(
        reads: &[ReadPlacements],
        family: &[usize],
        p: &ConflictParams,
    ) -> (usize, usize) {
        let fset: std::collections::BTreeSet<usize> = family.iter().copied().collect();
        let mut support = 0usize;
        let mut both_mapq0 = 0usize;
        'read: for placements in reads {
            // Scan every pair within this read; stop at the first pair that fires (count once per read).
            for ai in 0..placements.len() {
                for bi in (ai + 1)..placements.len() {
                    let (pa, pb) = (&placements[ai], &placements[bi]);
                    if !fset.contains(&pa.locus) || !fset.contains(&pb.locus) {
                        continue;
                    }
                    if pa.locus == pb.locus {
                        continue;
                    }
                    if tied(pa, pb, p) {
                        support += 1;
                        if pa.mapq == 0 && pb.mapq == 0 {
                            both_mapq0 += 1;
                        }
                        continue 'read;
                    }
                }
            }
        }
        (support, both_mapq0)
    }

    /// True ⟹ the two co-located copies are DISTINGUISHABLE by reads and must be kept separate;
    /// false ⟹ no read separates them (true K=0) and they may collapse. This is the χ(H) edge
    /// predicate restricted to a co-located pair. No new threshold: `min_reads` is the conflict
    /// floor, and the PSV/junction flag is gated upstream at `PSV_MIN_ALLELE_READS`.
    pub fn reads_distinguish(
        uniq_i: usize,
        uniq_j: usize,
        shared_psv_or_junction: bool,
        min_reads: usize,
    ) -> bool {
        uniq_i >= min_reads || uniq_j >= min_reads || shared_psv_or_junction
    }

    /// Per-locus UNIQUE-mapper count from placements: for each locus, how many reads place there with
    /// `mapq > 0` — the aligner's own uniqueness verdict (it found no competing placement to tie against). This
    /// is the raw per-copy signal carried on `DenovoTranscript::distinguishing_uniq` and consumed by
    /// `distinct_locus_reps`'s same-strand merge guard via `reads_distinguish`. Independent of `ConflictParams`:
    /// MAPQ, not the de-tie criterion, decides uniqueness here. Placements at an out-of-range locus are ignored
    /// (defensive; callers always size `n_loci` to the rep set the placements were built against).
    pub fn locus_unique_mapper_counts(reads: &[ReadPlacements], n_loci: usize) -> Vec<usize> {
        let mut counts = vec![0usize; n_loci];
        for placements in reads {
            for p in placements {
                if p.locus < n_loci && p.mapq > 0 {
                    counts[p.locus] += 1;
                }
            }
        }
        counts
    }

    /// Connected-component families over the conflict edges (union-find). Returns components of size `>= 2`
    /// (a locus with no conflict needs no resolution — it is not a family), each sorted ascending, the list
    /// sorted by first member (deterministic).
    pub fn conflict_families(n_loci: usize, edges: &[(usize, usize, usize)]) -> Vec<Vec<usize>> {
        let mut parent: Vec<usize> = (0..n_loci).collect();
        fn find(parent: &mut [usize], x: usize) -> usize {
            let mut r = x;
            while parent[r] != r {
                r = parent[r];
            }
            let mut c = x;
            while parent[c] != r {
                let next = parent[c];
                parent[c] = r;
                c = next;
            }
            r
        }
        for &(a, b, _) in edges {
            if a < n_loci && b < n_loci {
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                if ra != rb {
                    parent[ra.max(rb)] = ra.min(rb);
                }
            }
        }
        let mut groups: std::collections::BTreeMap<usize, Vec<usize>> =
            std::collections::BTreeMap::new();
        for x in 0..n_loci {
            let r = find(&mut parent, x);
            groups.entry(r).or_default().push(x);
        }
        let mut out: Vec<Vec<usize>> = groups.into_values().filter(|g| g.len() >= 2).collect();
        for g in &mut out {
            g.sort_unstable();
        }
        out.sort_by_key(|g| g[0]);
        out
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        fn p(locus: usize, de: f32) -> Placement {
            Placement {
                locus,
                de,
                mapq: 0,
                as_score: 100,
                aln_len: 2000,
            }
        }

        #[test]
        fn de_tied_placements_make_an_edge_and_a_family() {
            // both copies fit comparably (de 0.010 vs 0.012, both < 0.05) -> tie -> edge -> family.
            let reads = vec![vec![p(0, 0.010), p(1, 0.012)]];
            let edges = conflict_edges(
                2,
                &reads,
                &ConflictParams {
                    delta: 0.005,
                    de_max: 0.05,
                    min_reads: 1,
                    sig: None,
                },
            );
            assert_eq!(edges, vec![(0, 1, 1)]);
            assert_eq!(conflict_families(2, &edges), vec![vec![0, 1]]);
        }

        #[test]
        fn sig_criterion_ties_ambiguous_resolves_distinguishing() {
            // The UNIFICATION: significance edge mirrors the assignment gate's min_p>=alpha (eps^delta>=alpha).
            let sig = ConflictParams {
                delta: 0.005,
                de_max: 0.05,
                min_reads: 1,
                sig: Some((1e-3, 1e-3)),
            };
            // tau = floor(ln(alpha)/ln(eps)) = floor(ln(1e-3)/ln(1e-3)) = 1: tied iff excess-mismatches <= 1.
            // AMBIGUOUS read: equal divergence on both copies (delta_cols = 0) -> tied -> edge.
            let amb = vec![vec![
                Placement {
                    locus: 0,
                    de: 0.0050,
                    mapq: 0,
                    as_score: 100,
                    aln_len: 2000,
                },
                Placement {
                    locus: 1,
                    de: 0.0050,
                    mapq: 0,
                    as_score: 100,
                    aln_len: 2000,
                },
            ]];
            assert_eq!(conflict_edges(2, &amb, &sig), vec![(0, 1, 1)]);
            // DISTINGUISHING read: m_a=0.0005*2000=1, m_b=0.0050*2000=10 -> delta_cols=9 >> tau -> NO edge
            // (the boundary is now eps^delta>=alpha, not a hand-set 0.005).
            let dist = vec![vec![
                Placement {
                    locus: 0,
                    de: 0.0005,
                    mapq: 0,
                    as_score: 100,
                    aln_len: 2000,
                },
                Placement {
                    locus: 1,
                    de: 0.0050,
                    mapq: 0,
                    as_score: 100,
                    aln_len: 2000,
                },
            ]];
            assert!(conflict_edges(2, &dist, &sig).is_empty());
        }

        #[test]
        fn sig_edge_is_a_refinement_of_de_tied_equal_length() {
            // Exhaustive over a divergence grid (equal aligned length, default eps/alpha): every SIG-tie is also a
            // de-tie, so the significance edge set is a SUBSET of the de-tie edge set -> SIG can only shrink/split
            // families, never invent them (the rigorous refinement claim behind the 81->71 catalog narrowing).
            let de_p = ConflictParams {
                delta: 0.005,
                de_max: 0.05,
                min_reads: 1,
                sig: None,
            };
            let (eps, alpha, l) = (1e-3f64, 1e-3f64, 2000u32);
            for ia in 0..=60u32 {
                for ib in 0..=60u32 {
                    let (da, db) = (ia as f32 * 0.001, ib as f32 * 0.001);
                    let a = Placement {
                        locus: 0,
                        de: da,
                        mapq: 0,
                        as_score: 0,
                        aln_len: l,
                    };
                    let b = Placement {
                        locus: 1,
                        de: db,
                        mapq: 0,
                        as_score: 0,
                        aln_len: l,
                    };
                    if sig_tied(&a, &b, de_p.de_max, eps, alpha) {
                        assert!(
                            de_tied(&a, &b, &de_p),
                            "sig-tie not a de-tie at de=({da},{db})"
                        );
                    }
                }
            }
        }

        #[test]
        fn sig_off_default_is_byte_identical_to_de_tied() {
            // With sig: None (default), `tied` == `de_tied` exactly: same edges as the legacy criterion.
            let reads = vec![
                vec![p(0, 0.010), p(1, 0.012)], // tied under delta=0.005
                vec![p(0, 0.001), p(1, 0.020)], // resolved under delta=0.005
            ];
            let legacy = ConflictParams {
                delta: 0.005,
                de_max: 0.05,
                min_reads: 1,
                sig: None,
            };
            assert_eq!(conflict_edges(2, &reads, &legacy), vec![(0, 1, 1)]); // only the first read ties
            assert_eq!(ConflictParams::default().sig, None); // default ships OFF
        }

        #[test]
        fn divergence_gap_beyond_delta_is_not_a_conflict() {
            // read fits copy 0 (de 0.001) far better than copy 1 (de 0.020): |Δ|=0.019 > 0.005 -> resolvable.
            let reads = vec![vec![p(0, 0.001), p(1, 0.020)]];
            let edges = conflict_edges(
                2,
                &reads,
                &ConflictParams {
                    delta: 0.005,
                    de_max: 0.05,
                    min_reads: 1,
                    sig: None,
                },
            );
            assert!(edges.is_empty());
        }

        #[test]
        fn both_high_divergence_blocked_by_ceiling() {
            // de_a 0.06 ~ de_b 0.061 (tied within delta) but both exceed de_max 0.05 -> read fits neither.
            let reads = vec![vec![p(0, 0.060), p(1, 0.061)]];
            let edges = conflict_edges(
                2,
                &reads,
                &ConflictParams {
                    delta: 0.005,
                    de_max: 0.05,
                    min_reads: 1,
                    sig: None,
                },
            );
            assert!(edges.is_empty());
        }

        #[test]
        fn single_placement_read_is_a_singleton_not_a_family() {
            let reads = vec![vec![p(0, 0.01)], vec![p(1, 0.01)]];
            let edges = conflict_edges(
                2,
                &reads,
                &ConflictParams {
                    delta: 0.005,
                    de_max: 0.05,
                    min_reads: 1,
                    sig: None,
                },
            );
            assert!(edges.is_empty());
            assert!(conflict_families(2, &edges).is_empty());
        }

        #[test]
        fn min_reads_threshold_drops_thin_conflicts() {
            let one = vec![vec![p(0, 0.01), p(1, 0.012)]];
            let pr = ConflictParams {
                delta: 0.005,
                de_max: 0.05,
                min_reads: 3,
                sig: None,
            };
            assert!(conflict_edges(2, &one, &pr).is_empty());
            let three = vec![one[0].clone(), one[0].clone(), one[0].clone()];
            assert_eq!(conflict_edges(2, &three, &pr), vec![(0, 1, 3)]);
        }

        #[test]
        fn transitive_conflict_closes_into_one_family() {
            let reads = vec![
                vec![p(0, 0.010), p(1, 0.012)],
                vec![p(1, 0.010), p(2, 0.013)],
            ];
            let edges = conflict_edges(
                3,
                &reads,
                &ConflictParams {
                    delta: 0.005,
                    de_max: 0.05,
                    min_reads: 1,
                    sig: None,
                },
            );
            assert_eq!(conflict_families(3, &edges), vec![vec![0, 1, 2]]);
        }

        #[test]
        fn disjoint_conflicts_form_separate_families() {
            let reads = vec![
                vec![p(0, 0.01), p(1, 0.012)],
                vec![p(2, 0.01), p(3, 0.012)],
                vec![p(4, 0.01)],
            ];
            let edges = conflict_edges(
                5,
                &reads,
                &ConflictParams {
                    delta: 0.005,
                    de_max: 0.05,
                    min_reads: 1,
                    sig: None,
                },
            );
            assert_eq!(conflict_families(5, &edges), vec![vec![0, 1], vec![2, 3]]);
        }

        #[test]
        fn default_params_are_the_operating_point() {
            let d = ConflictParams::default();
            assert!((d.delta - 0.005).abs() < 1e-9);
            assert!((d.de_max - 0.05).abs() < 1e-9);
            assert_eq!(d.min_reads, 3);
        }

        #[test]
        fn deterministic_under_placement_order() {
            // Shuffling placement order within a read and across reads must not change the edge/family output.
            let pr = ConflictParams {
                delta: 0.005,
                de_max: 0.05,
                min_reads: 1,
                sig: None,
            };
            // two reads, each with placements on loci {0,1,2} in different orders.
            let forward = vec![
                vec![p(0, 0.010), p(1, 0.012), p(2, 0.030)],
                vec![p(1, 0.011), p(0, 0.013)],
            ];
            let reversed = vec![
                vec![p(2, 0.030), p(1, 0.012), p(0, 0.010)],
                vec![p(0, 0.013), p(1, 0.011)],
            ];
            let reads_swapped = vec![forward[1].clone(), forward[0].clone()];
            let e_fwd = conflict_edges(3, &forward, &pr);
            let e_rev = conflict_edges(3, &reversed, &pr);
            let e_swp = conflict_edges(3, &reads_swapped, &pr);
            assert_eq!(
                e_fwd, e_rev,
                "placement order within reads must not change edges"
            );
            assert_eq!(e_fwd, e_swp, "read order must not change edges");
            assert_eq!(conflict_families(3, &e_fwd), conflict_families(3, &e_rev));
            assert_eq!(conflict_families(3, &e_fwd), conflict_families(3, &e_swp));
        }

        #[test]
        fn de_max_boundary_exactly_at_threshold_fires_just_over_is_blocked() {
            // Safely inside: de 0.049 vs 0.049 — max(de)=0.049 < 0.05, |Δ|=0 ≤ 0.005 → should fire.
            let pr = ConflictParams {
                delta: 0.005,
                de_max: 0.05,
                min_reads: 1,
                sig: None,
            };
            let inside = vec![vec![p(0, 0.049), p(1, 0.049)]];
            let e_inside = conflict_edges(2, &inside, &pr);
            assert_eq!(
                e_inside,
                vec![(0, 1, 1)],
                "de=0.049 <= de_max=0.05 must fire"
            );
            // NOTE: exactly 0.05f32 widens above 0.05f64 after f32→f64 cast → on the precision boundary;
            // we test the clearly-over case (0.051) which is blocked regardless.
            let over = vec![vec![p(0, 0.051), p(1, 0.051)]];
            let e_over = conflict_edges(2, &over, &pr);
            assert!(e_over.is_empty(), "de=0.051 > de_max=0.05 must be blocked");
        }

        #[test]
        fn family_mapq0_support_counts_correctly() {
            let pr = ConflictParams {
                delta: 0.005,
                de_max: 0.05,
                min_reads: 1,
                sig: None,
            };
            let family = vec![0usize, 1];
            // read A: both loci de-tied, both mapq==0 → counts in support AND both_mapq0.
            let read_a: ReadPlacements = vec![
                Placement {
                    locus: 0,
                    de: 0.010,
                    mapq: 0,
                    as_score: 100,
                    aln_len: 2000,
                },
                Placement {
                    locus: 1,
                    de: 0.012,
                    mapq: 0,
                    as_score: 100,
                    aln_len: 2000,
                },
            ];
            // read B: both loci de-tied, but mapq>0 on one → counts in support but NOT both_mapq0.
            let read_b: ReadPlacements = vec![
                Placement {
                    locus: 0,
                    de: 0.010,
                    mapq: 60,
                    as_score: 100,
                    aln_len: 2000,
                },
                Placement {
                    locus: 1,
                    de: 0.012,
                    mapq: 0,
                    as_score: 100,
                    aln_len: 2000,
                },
            ];
            // read C: de NOT tied (gap too large) → not counted at all.
            let read_c: ReadPlacements = vec![
                Placement {
                    locus: 0,
                    de: 0.001,
                    mapq: 0,
                    as_score: 100,
                    aln_len: 2000,
                },
                Placement {
                    locus: 1,
                    de: 0.020,
                    mapq: 0,
                    as_score: 100,
                    aln_len: 2000,
                },
            ];
            let reads = vec![read_a, read_b, read_c];
            let (support, mapq0) = family_mapq0_support(&reads, &family, &pr);
            assert_eq!(
                support, 2,
                "read_a and read_b both contribute a de-tied pair in the family"
            );
            assert_eq!(mapq0, 1, "only read_a has both placements mapq==0");
        }

        #[test]
        fn as_tie_edges_superset_of_de_edges() {
            // AS ties two placements that de SPLITS (de 0.001 vs 0.020): AS-edge exists, de-edge does not.
            let reads = vec![vec![
                Placement {
                    locus: 0,
                    de: 0.001,
                    mapq: 0,
                    as_score: 500,
                    aln_len: 2000,
                },
                Placement {
                    locus: 1,
                    de: 0.020,
                    mapq: 0,
                    as_score: 498,
                    aln_len: 2000,
                },
            ]];
            let de_edges = conflict_edges(
                2,
                &reads,
                &ConflictParams {
                    delta: 0.005,
                    de_max: 0.05,
                    min_reads: 1,
                    sig: None,
                },
            );
            let as_edges = as_tie_edges(2, &reads, 0.9, 1);
            assert!(de_edges.is_empty());
            assert_eq!(as_edges, std::collections::BTreeSet::from([(0, 1)]));
        }

        #[test]
        fn as_evidence_none_without_placements() {
            assert_eq!(as_evidence(&[]), None);
        }

        #[test]
        fn as_evidence_single_placement_has_no_second_or_margin() {
            let e = as_evidence(&[(1842, 916)]).unwrap();
            assert_eq!(e.best, 1842);
            assert_eq!(e.second, None);
            assert_eq!(e.second_per_base, None);
            assert_eq!(
                e.margin(),
                None,
                "one placement => nothing to be ambiguous between"
            );
        }

        #[test]
        fn as_evidence_ranks_best_and_runner_up_over_three_placements() {
            let e = as_evidence(&[(1310, 800), (1842, 916), (900, 500)]).unwrap();
            assert_eq!(e.best, 1842);
            assert_eq!(e.second, Some(1310));
            assert_eq!(e.margin(), Some(532));
            assert!((e.best_per_base - 1842.0 / 916.0).abs() < 1e-6);
            assert!((e.second_per_base.unwrap() - 1310.0 / 800.0).abs() < 1e-6);
        }

        #[test]
        fn as_evidence_zero_length_placement_does_not_divide_by_zero() {
            let e = as_evidence(&[(50, 0)]).unwrap();
            assert_eq!(e.best_per_base, 0.0);
        }

        /// The measured real-data confound, pinned: two placements of EQUAL per-base quality where the runner-up
        /// is a partial alignment. Raw AS says "not a tie" (ratio 0.60 < 0.90) while per-base AS says they are
        /// indistinguishable. This is why AS is reported and `de` decides.
        #[test]
        fn raw_as_misjudges_a_partial_placement_that_per_base_as_calls_equal() {
            let (full, partial) = ((2000, 1000), (1200, 600)); // both exactly 2.0 AS per aligned base
            let e = as_evidence(&[full, partial]).unwrap();
            assert_eq!(e.margin(), Some(800), "raw AS shows a large margin...");
            assert!(
                !as_tied(full.0, partial.0, 0.9),
                "...so the legacy AS-tie predicate rejects the pair"
            );
            assert_eq!(
                e.best_per_base,
                e.second_per_base.unwrap(),
                "...yet per-aligned-base they are identical: the margin is pure length confound"
            );
        }

        #[test]
        fn reads_distinguish_keeps_separate_when_unique_mappers_present() {
            // one copy has 40 unique reads (the ID_26 case): distinguishable -> keep separate
            assert!(reads_distinguish(40, 0, false, 3));
            // both sides above the floor: distinguishable
            assert!(reads_distinguish(11, 8, false, 3));
            // a read-supported PSV/junction separates them even with no unique mappers
            assert!(reads_distinguish(0, 0, true, 3));
        }

        #[test]
        fn locus_unique_mapper_counts_counts_mapq_positive_placements_per_locus() {
            let reads = vec![
                vec![Placement {
                    locus: 0,
                    de: 0.01,
                    mapq: 60,
                    as_score: 100,
                    aln_len: 2000,
                }], // uniq @0
                vec![Placement {
                    locus: 0,
                    de: 0.01,
                    mapq: 0,
                    as_score: 100,
                    aln_len: 2000,
                }], // ambiguous, not counted
                vec![Placement {
                    locus: 1,
                    de: 0.01,
                    mapq: 40,
                    as_score: 100,
                    aln_len: 2000,
                }], // uniq @1
                vec![Placement {
                    locus: 1,
                    de: 0.01,
                    mapq: 40,
                    as_score: 100,
                    aln_len: 2000,
                }], // uniq @1
            ];
            assert_eq!(locus_unique_mapper_counts(&reads, 2), vec![1, 2]);
        }

        #[test]
        fn reads_distinguish_merges_true_k0() {
            // no unique mappers either side, no distinguishing PSV/junction -> K=0 -> merge
            assert!(!reads_distinguish(0, 0, false, 3));
            // unique support below the floor is noise, not a distinction -> merge
            assert!(!reads_distinguish(2, 1, false, 3));
        }
    }
}

pub mod family_rescue {
    //! Family-aware copy RESCUE (borrow strength / partial pooling): recover an
    //! under-ASSEMBLED copy that falls below the general-assembly read gate but is
    //! sequence-homologous to a CONFIRMED de-novo family.
    //!
    //! Rationale (mirrors `bench/family_rescue.py`): a single confident
    //! canonical-junction read forming a multi-exon chain that POA-confirms against
    //! an existing family is strong evidence of a real copy -- the family PRIOR
    //! overcomes thin read support. This recovers expression-limited copies (e.g.
    //! RFPL2, a single read) WITHOUT lowering the genome-wide `>= 3`-read gate
    //! (fuzzy / low-gate experiments only added noise; see the pipeline record).
    //!
    //! This module is the testable DECISION CORE. The BAM-neighbourhood scan that
    //! builds the thin candidate loci (intron-chain collapse, member-span exclusion,
    //! canonical-junction strand check, spliced-sequence construction) is the
    //! integration layer (`gen2off`-style boundary mapping + BAM/FASTA), kept out of
    //! here so the rescue logic is unit-testable against `contiguous_core_coverage`.
    //!
    //! The decision, per candidate thin locus:
    //!   1. A CANONICAL exact base-4 k-mer pre-filter (`canonical_kmer_set`,
    //!      `KMER = 18`, strand-symmetric) selects the SINGLE best family member by
    //!      shared-k-mer overlap (`>= K_RESCUE`). This bounds POA to true candidates
    //!      and never decides membership -- POA does.
    //!   2. POA-confirm the thin locus against that one member via the validated
    //!      `contiguous_core_coverage` primitive, trying BOTH orientations (forward,
    //!      then a reverse-complement fallback for copies assembled on the opposite
    //!      strand). Rescue iff the contiguous-core coverage `>= T_CORE = 0.13`.
    //!
    //! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

    use crate::family::seq_utils::reverse_complement;
    use crate::types::DetHashSet;

    /// Exact k-mer length for the canonical pre-filter. `KMER = 18` gives a
    /// `4^18 ~ 6.9e10` space (negligible coincidental sharing) and fits in `u64`
    /// (36 bits). Matches `bench/denovo_families.py::KMER`.
    pub const KMER: usize = 18;

    /// POA contiguous-core coverage threshold to confirm a rescued copy. A true
    /// recent-duplicate copy shares one long homologous core (`>= 0.13` of the
    /// shorter sequence); a domain-sharer co-aligns over only a short block.
    /// Matches `bench/family_rescue.py::T_CORE`.
    pub const T_CORE: f64 = 0.13;

    /// Minimum number of shared canonical k-mers a thin locus must have with a
    /// family member to even reach POA (real copies share many; common-domain
    /// sharers share few). Matches `bench/family_rescue.py::K_RESCUE`.
    pub const K_RESCUE: usize = 20;

    /// Length cap (bp) on the POA inputs (cost guard; POA is O(L^2)). Matches
    /// `bench/family_rescue.py::LEN_CAP`.
    pub const LEN_CAP: usize = 9000;

    /// 2-bit base code A=0, C=1, G=2, T=3 (lowercase accepted), matching the python
    /// `_CODE` table. Returns `None` for any non-ACGT base (so a window touching it
    /// is dropped, as in the python `badwin` mask).
    fn base_code(b: u8) -> Option<u8> {
        match b {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }

    /// CANONICAL base-4 code of a single window = `min(forward, reverse-complement)`
    /// so a window and its reverse-complement collapse to the SAME code (strand
    /// symmetry). Returns `None` if the window touches any non-ACGT base.
    ///
    /// Faithful to `bench/denovo_families.py::kmer_hashes`: forward Horner code
    /// `Σ base_t · 4^(k-1-t)`; reverse-complement code substitutes `3 - base` and
    /// reverses the position order.
    ///
    /// `pub(crate)` so the family-detection layer can reuse the exact same canonical
    /// encoder for its position-aware k-mer signatures.
    pub(crate) fn window_canon_code(window: &[u8]) -> Option<u64> {
        let mut fwd: u64 = 0;
        let mut rc: u64 = 0;
        for (i, &b) in window.iter().enumerate() {
            let c = base_code(b)?; // None drops a window touching a non-ACGT base
            fwd = (fwd << 2) | c as u64;
            // reverse-complement code: complement (3 - c), placed in reverse order.
            rc |= ((3 - c) as u64) << (2 * i);
        }
        Some(fwd.min(rc))
    }

    /// CANONICAL exact base-4 k-mer set (`KMER`-mers, `min(fwd, rc)` per window;
    /// windows touching a non-ACGT base dropped). A sequence and its reverse
    /// complement yield the SAME set, so rescue is strand-symmetric. Mirrors
    /// `set(kmer_hashes(seq))` in the python pipeline.
    pub fn canonical_kmer_set(seq: &[u8]) -> DetHashSet<u64> {
        let mut set = DetHashSet::default();
        if seq.len() < KMER {
            return set;
        }
        for w in seq.windows(KMER) {
            if let Some(code) = window_canon_code(w) {
                set.insert(code);
            }
        }
        set
    }

    /// Number of shared canonical k-mers between two sets (the pre-filter overlap).
    pub fn kmer_overlap(a: &DetHashSet<u64>, b: &DetHashSet<u64>) -> usize {
        // Iterate the smaller set for speed; the count is order-independent.
        let (small, large) = if a.len() <= b.len() { (a, b) } else { (b, a) };
        small.iter().filter(|k| large.contains(k)).count()
    }

    /// A confirmed de-novo family member the rescue compares thin loci against. Its
    /// canonical k-mer set is precomputed once (the python pipeline caches
    /// `memkmers[t]`).
    #[derive(Clone, Debug)]
    pub struct FamilyMember {
        /// Member transcript id (e.g. a de-novo transcript id).
        pub tid: String,
        /// The family this member belongs to.
        pub family_id: String,
        /// The member's spliced sequence in transcription orientation.
        pub seq: Vec<u8>,
        /// Precomputed canonical k-mer set of `seq`.
        pub kmers: DetHashSet<u64>,
    }

    impl FamilyMember {
        /// Build a member, computing its canonical k-mer set from `seq`.
        pub fn new(tid: String, family_id: String, seq: Vec<u8>) -> Self {
            let kmers = canonical_kmer_set(&seq);
            FamilyMember {
                tid,
                family_id,
                seq,
                kmers,
            }
        }
    }

    /// Which orientation of the member the contiguous core was found in.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub enum Orientation {
        /// The member's stored orientation aligned the core directly.
        Forward,
        /// The reverse-complement of the member aligned the core (opposite-strand
        /// assembly).
        ReverseComplement,
    }

    /// Tunable rescue parameters (defaults mirror `bench/family_rescue.py`).
    #[derive(Clone, Copy, Debug)]
    pub struct RescueParams {
        /// POA contiguous-core coverage required to confirm a copy.
        pub t_core: f64,
        /// Minimum shared canonical k-mers with the best member to reach POA.
        pub k_rescue: usize,
        /// Length cap on POA inputs.
        pub len_cap: usize,
    }

    impl Default for RescueParams {
        fn default() -> Self {
            RescueParams {
                t_core: T_CORE,
                k_rescue: K_RESCUE,
                len_cap: LEN_CAP,
            }
        }
    }

    /// A rescued copy: which family it joins, the member that confirmed it, the POA
    /// contiguous-core coverage, and the orientation the core aligned in.
    #[derive(Clone, Debug)]
    pub struct RescueOutcome {
        pub family_id: String,
        pub best_member: String,
        pub core_recip: f64,
        pub orientation: Orientation,
    }

    /// Decide whether a thin-locus spliced sequence is a rescued copy of one of the
    /// confirmed families represented by `members`.
    ///
    /// `members` is the candidate neighbourhood the integration layer already
    /// windowed (the python `WIN = 1 Mb` member set near the locus). The function
    /// pre-filters to the SINGLE best member by canonical-k-mer overlap
    /// (`>= p.k_rescue`), then POA-confirms via `contiguous_core_coverage` in the
    /// forward orientation, falling back to the reverse complement if forward is
    /// below `p.t_core`. Returns the rescue iff the best coverage reaches
    /// `p.t_core`. Deterministic: members are scanned in slice order and the best
    /// overlap ties to the earliest member (matching the python `if ov > best_ov`).
    pub fn rescue_thin_locus(
        thin_seq: &[u8],
        members: &[FamilyMember],
        p: &RescueParams,
    ) -> Option<RescueOutcome> {
        let kset = canonical_kmer_set(thin_seq);
        if kset.is_empty() {
            return None;
        }

        // Pre-filter: the SINGLE best member by canonical-k-mer overlap. `best_ov`
        // starts at `k_rescue - 1` so a member must clear `>= k_rescue` to be picked
        // (python `best_ov = K_RESCUE - 1; if ov > best_ov`). Strict `>` ties to the
        // earliest member in slice order (determinism).
        let mut best: Option<&FamilyMember> = None;
        let mut best_ov = p.k_rescue.saturating_sub(1);
        for m in members {
            // POA cost guard (python `_poa_rescue`: skip if min(len) > LEN_CAP).
            if thin_seq.len().min(m.seq.len()) > p.len_cap {
                continue;
            }
            let ov = kmer_overlap(&kset, &m.kmers);
            if ov > best_ov {
                best_ov = ov;
                best = Some(m);
            }
        }
        let m = best?;

        // POA-confirm against that one member. The k-mer pre-filter is
        // case-insensitive (like python's `_CODE` table), so uppercase the POA
        // operands here to match it: python's `poa_pair_stats` uppercases its inputs,
        // and `reverse_complement` maps lowercase -> N, which would otherwise silently
        // void the RC fallback on soft-masked input. Forward first; reverse-complement
        // fallback only if forward is below threshold, keeping the max (mirrors the
        // python `_poa_rescue` RC retry).
        use crate::family::family_detect::family_graph::{
            contiguous_core_coverage_bounded_with, EDGE_CONFIRM_ASTAR,
        };
        let thin_up = crate::family::family_detect::family_graph::upper_cow(thin_seq);
        let mem_up = crate::family::family_detect::family_graph::upper_cow(&m.seq);
        let mut core_recip =
            contiguous_core_coverage_bounded_with(&thin_up, &mem_up, p.len_cap, EDGE_CONFIRM_ASTAR);
        let mut orientation = Orientation::Forward;
        if core_recip < p.t_core {
            let rc = contiguous_core_coverage_bounded_with(
                &thin_up,
                &reverse_complement(&mem_up),
                p.len_cap,
                EDGE_CONFIRM_ASTAR,
            );
            if rc > core_recip {
                core_recip = rc;
                orientation = Orientation::ReverseComplement;
            }
        }

        if core_recip >= p.t_core {
            Some(RescueOutcome {
                family_id: m.family_id.clone(),
                best_member: m.tid.clone(),
                core_recip,
                orientation,
            })
        } else {
            None
        }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        // Deterministic DNA (SplitMix64), mirroring family_graph's core-coverage tests.
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
        fn member(tid: &str, fid: &str, seq: Vec<u8>) -> FamilyMember {
            FamilyMember::new(tid.to_string(), fid.to_string(), seq)
        }

        // ---- canonical k-mer encoding (faithful to python base-4 codes) ----

        #[test]
        fn window_canon_code_is_base4_and_strand_canonical() {
            // A=0,C=1,G=2,T=3 ; fwd("AC")=0*4+1=1 ; RC("AC")="GT"=2*4+3=11 ; min=1.
            assert_eq!(window_canon_code(b"AC"), Some(1));
            // "GT" is the reverse-complement of "AC" -> SAME canonical code.
            assert_eq!(window_canon_code(b"GT"), Some(1));
            // a window containing a non-ACGT base is dropped.
            assert_eq!(window_canon_code(b"AN"), None);
        }

        #[test]
        fn canonical_kmer_set_is_strand_symmetric() {
            let s = rand_seq(200, 0x5EED_0001);
            let rc = reverse_complement(&s);
            assert_eq!(canonical_kmer_set(&s), canonical_kmer_set(&rc));
            assert!(!canonical_kmer_set(&s).is_empty());
        }

        #[test]
        fn canonical_kmer_set_drops_only_windows_touching_n() {
            // 60 bp random -> 60-18+1 = 43 windows; random 18-mers are distinct.
            let s = rand_seq(60, 0x5EED_0002);
            let full = canonical_kmer_set(&s);
            assert_eq!(full.len(), 60 - KMER + 1, "distinct random 18-mers");
            let mut withn = s.clone();
            withn[30] = b'N';
            let dropped = canonical_kmer_set(&withn);
            // windows with start in [13, 30] cover index 30 -> 18 windows removed.
            assert_eq!(full.len() - dropped.len(), 18);
        }

        #[test]
        fn kmer_overlap_counts_shared_canonical_kmers() {
            let s = rand_seq(60, 0x7001);
            let a = canonical_kmer_set(&s);
            let b = canonical_kmer_set(&reverse_complement(&s));
            assert_eq!(kmer_overlap(&a, &b), a.len()); // identical canonical sets
            let c = canonical_kmer_set(&rand_seq(60, 0x7002));
            assert!(kmer_overlap(&a, &c) < a.len()); // independent seqs share ~none
        }

        // ---- rescue decision ----

        #[test]
        fn rescue_confirms_homologous_copy() {
            // thin locus and a member share a 400 bp identical core (>> K_RESCUE
            // shared 18-mers, contiguous-core coverage >= 0.13), divergent flanks.
            let core = rand_seq(400, 0xC0FE_2001);
            let thin = cat(&[
                &rand_seq(80, 0xAAAA_2001),
                &core,
                &rand_seq(80, 0xAAAA_2002),
            ]);
            let mseq = cat(&[
                &rand_seq(80, 0xBBBB_2001),
                &core,
                &rand_seq(80, 0xBBBB_2002),
            ]);
            let members = [member("M1", "FAM7", mseq)];
            let out = rescue_thin_locus(&thin, &members, &RescueParams::default())
                .expect("homologous copy should be rescued");
            assert_eq!(out.family_id, "FAM7");
            assert_eq!(out.best_member, "M1");
            assert_eq!(out.orientation, Orientation::Forward);
            assert!(
                out.core_recip >= T_CORE,
                "core_recip {} >= {}",
                out.core_recip,
                T_CORE
            );
        }

        #[test]
        fn rescue_rejects_domain_sharer_below_core_threshold() {
            // shares a 40 bp block (40-18+1 = 23 shared 18-mers >= K_RESCUE so the
            // PRE-FILTER passes) but in long otherwise-independent sequences, so the
            // POA contiguous-core coverage is < T_CORE -> NOT rescued (POA decides).
            let block = rand_seq(40, 0xD0D0_3001);
            let thin = cat(&[
                &rand_seq(350, 0xAAAA_3001),
                &block,
                &rand_seq(350, 0xAAAA_3002),
            ]);
            let mseq = cat(&[
                &rand_seq(350, 0xBBBB_3001),
                &block,
                &rand_seq(350, 0xBBBB_3002),
            ]);
            let members = [member("M1", "FAM3", mseq)];
            // precondition: the pre-filter DOES pass (>= K_RESCUE shared k-mers).
            let thin_k = canonical_kmer_set(&thin);
            assert!(
                kmer_overlap(&thin_k, &members[0].kmers) >= K_RESCUE,
                "pre-filter precondition: shared k-mers >= K_RESCUE"
            );
            assert!(rescue_thin_locus(&thin, &members, &RescueParams::default()).is_none());
        }

        #[test]
        fn rescue_pre_filter_rejects_too_few_shared_kmers() {
            // shares only a 25 bp block -> 25-18+1 = 8 shared 18-mers < K_RESCUE
            // -> NO member clears the pre-filter -> no POA, no rescue.
            let block = rand_seq(25, 0xD0D0_4001);
            let thin = cat(&[
                &rand_seq(300, 0xAAAA_4001),
                &block,
                &rand_seq(300, 0xAAAA_4002),
            ]);
            let mseq = cat(&[
                &rand_seq(300, 0xBBBB_4001),
                &block,
                &rand_seq(300, 0xBBBB_4002),
            ]);
            let members = [member("M1", "FAM4", mseq)];
            assert!(rescue_thin_locus(&thin, &members, &RescueParams::default()).is_none());
        }

        #[test]
        fn rescue_detects_reverse_complement_copy() {
            // A member stored on the OPPOSITE strand: canonical k-mers still match
            // (strand-symmetric), forward POA is LOW, RC fallback POA is HIGH.
            let core = rand_seq(400, 0xC0FE_5001);
            let thin = cat(&[
                &rand_seq(80, 0xAAAA_5001),
                &core,
                &rand_seq(80, 0xAAAA_5002),
            ]);
            let mfwd = cat(&[
                &rand_seq(80, 0xBBBB_5001),
                &core,
                &rand_seq(80, 0xBBBB_5002),
            ]);
            let mseq = reverse_complement(&mfwd); // member assembled on the other strand
            let members = [member("M1", "FAM5", mseq)];
            let out = rescue_thin_locus(&thin, &members, &RescueParams::default())
                .expect("RC homologous copy should be rescued via the RC fallback");
            assert_eq!(out.orientation, Orientation::ReverseComplement);
            assert!(
                out.core_recip >= T_CORE,
                "core_recip {} >= {}",
                out.core_recip,
                T_CORE
            );
        }

        #[test]
        fn rescue_picks_best_kmer_matching_member() {
            // thin shares its full 400 bp core with M2 (many k-mers) and nothing with
            // M1 (unrelated); POA must confirm against M2 and report it.
            let core = rand_seq(400, 0xC0FE_6001);
            let thin = cat(&[
                &rand_seq(80, 0xAAAA_6001),
                &core,
                &rand_seq(80, 0xAAAA_6002),
            ]);
            let m1 = rand_seq(660, 0x1111_6001); // fully unrelated
            let m2 = cat(&[
                &rand_seq(80, 0xBBBB_6001),
                &core,
                &rand_seq(80, 0xBBBB_6002),
            ]);
            let members = [member("M1", "FAMx", m1), member("M2", "FAM6", m2)];
            let out = rescue_thin_locus(&thin, &members, &RescueParams::default())
                .expect("should rescue via M2");
            assert_eq!(out.best_member, "M2");
            assert_eq!(out.family_id, "FAM6");
        }

        #[test]
        fn rescue_empty_members_returns_none() {
            let thin = rand_seq(300, 0x9001);
            assert!(rescue_thin_locus(&thin, &[], &RescueParams::default()).is_none());
        }

        #[test]
        fn rescue_thin_seq_shorter_than_kmer_returns_none() {
            let thin = rand_seq(10, 0x9002); // < KMER -> empty k-mer set -> None
            let m = member("M1", "F", rand_seq(400, 0x9003));
            assert!(rescue_thin_locus(&thin, &[m], &RescueParams::default()).is_none());
        }

        #[test]
        fn rescue_detects_reverse_complement_copy_lowercase() {
            // Soft-masked (lowercase) inputs must behave IDENTICALLY to uppercase. The
            // k-mer pre-filter is case-insensitive, so a lowercase opposite-strand copy
            // reaches POA; the RC fallback must not be silently killed by a
            // case-sensitive reverse_complement (which maps lowercase -> N). Same
            // geometry as rescue_detects_reverse_complement_copy, just lowercased.
            let core = rand_seq(400, 0xC0FE_5001);
            let thin = cat(&[
                &rand_seq(80, 0xAAAA_5001),
                &core,
                &rand_seq(80, 0xAAAA_5002),
            ]);
            let mfwd = cat(&[
                &rand_seq(80, 0xBBBB_5001),
                &core,
                &rand_seq(80, 0xBBBB_5002),
            ]);
            let thin_lc = thin.to_ascii_lowercase();
            let mseq_lc = reverse_complement(&mfwd).to_ascii_lowercase();
            let members = [member("M1", "FAM5", mseq_lc)];
            let out = rescue_thin_locus(&thin_lc, &members, &RescueParams::default())
                .expect("lowercase RC homologous copy should still be rescued via the RC fallback");
            assert_eq!(out.orientation, Orientation::ReverseComplement);
            assert!(
                out.core_recip >= T_CORE,
                "core_recip {} >= {}",
                out.core_recip,
                T_CORE
            );
        }

        #[test]
        fn rescue_tie_break_picks_earliest_member() {
            // Two members with the IDENTICAL sequence (so provably EQUAL k-mer overlap
            // with the thin locus); the strict `>` update must keep the FIRST in slice
            // order. (Distinct random flanks do NOT guarantee equal overlap — a flank
            // k-mer can coincidentally collide — so identical seqs pin the tie exactly.)
            let core = rand_seq(400, 0xC0FE_7001);
            let thin = cat(&[
                &rand_seq(80, 0xAAAA_7001),
                &core,
                &rand_seq(80, 0xAAAA_7002),
            ]);
            let mseq = cat(&[
                &rand_seq(80, 0x1111_7001),
                &core,
                &rand_seq(80, 0x1111_7002),
            ]);
            let members = [
                member("FIRST", "FAM7", mseq.clone()),
                member("SECOND", "FAM7", mseq),
            ];
            // precondition: the two overlaps are genuinely equal (the tie under test).
            let tk = canonical_kmer_set(&thin);
            assert_eq!(
                kmer_overlap(&tk, &members[0].kmers),
                kmer_overlap(&tk, &members[1].kmers),
                "members must have equal overlap for this to test the tie-break"
            );
            let out = rescue_thin_locus(&thin, &members, &RescueParams::default())
                .expect("should rescue");
            assert_eq!(
                out.best_member, "FIRST",
                "equal overlap must tie to the earliest member"
            );
        }

        #[test]
        fn rescue_threshold_boundary_is_inclusive() {
            // Pin the `core_recip >= t_core` semantics: at t_core == cr the copy is
            // ACCEPTED (inclusive), just above cr it is REJECTED. `cr` is read off the
            // SAME deterministic primitive the function uses (forward orientation).
            use crate::family::family_detect::family_graph::contiguous_core_coverage;
            let core = rand_seq(400, 0xC0FE_8001);
            let thin = cat(&[
                &rand_seq(80, 0xAAAA_8001),
                &core,
                &rand_seq(80, 0xAAAA_8002),
            ]);
            let mseq = cat(&[
                &rand_seq(80, 0xBBBB_8001),
                &core,
                &rand_seq(80, 0xBBBB_8002),
            ]);
            let cr = contiguous_core_coverage(&thin, &mseq);
            let members = [member("M1", "FAM8", mseq)];
            let at = RescueParams {
                t_core: cr,
                ..RescueParams::default()
            };
            assert!(
                rescue_thin_locus(&thin, &members, &at).is_some(),
                "t_core == cr ({cr}) must be ACCEPTED (inclusive >=)"
            );
            let above = RescueParams {
                t_core: cr + 1e-6,
                ..RescueParams::default()
            };
            assert!(
                rescue_thin_locus(&thin, &members, &above).is_none(),
                "t_core just above cr ({cr}) must be REJECTED"
            );
        }
    }

    // ---- merged 2026-10-05: was `vg_family/rescue_pipeline.rs`, now the inline module below (one component) ----
    #[allow(clippy::all)]
    pub mod rescue_pipeline {
        //! Family-aware RESCUE thin-locus scan (integration stage 4b) — the BAM-neighbourhood orchestration of
        //! `bench/family_rescue.py`.
        //!
        //! Recover an under-ASSEMBLED copy that falls below the general-assembly `>= 3`-read gate but is
        //! sequence-homologous to a CONFIRMED de-novo family (borrow strength / partial pooling). The decision core
        //! (`family_rescue::rescue_thin_locus`: canonical-k-mer pre-filter + both-orientation POA) is already
        //! ported; this builds the thin candidate loci from reads in the family neighbourhood:
        //!   1. group reads by exact intron chain, collapse OVERLAPPING chains into loci (best-supported chain wins);
        //!   2. drop loci overlapping an already-assembled family member span;
        //!   3. build each thin locus's spliced sequence (canonical-junction gate, reverse-complement for `-`);
        //!   4. POA-confirm against the family members; dedup by locus, keep the best core coverage.
        //!
        //! **STATUS:** SHIPPED-DEFAULT  (docs/MODULE_STATUS.md; assigned by reachability, not by this header)

        use std::collections::{BTreeMap, BTreeSet};

        use super::{rescue_thin_locus, FamilyMember, RescueOutcome, RescueParams};
        use crate::family::denovo_assemble::{build_spliced_seq, PrimaryRead};
        use crate::genome::GenomeIndex;

        /// Minimum read support for a thin candidate locus (`family_rescue.py::MIN_SUPPORT`).
        pub const RESCUE_MIN_SUPPORT: u32 = 1;
        /// Minimum spliced length of a thin locus to attempt rescue.
        pub const RESCUE_MIN_LEN: usize = 200;

        /// A thin candidate locus: a (collapsed) intron chain with read support, below the general-assembly gate.
        #[derive(Clone, Debug, PartialEq, Eq)]
        pub struct ThinLocus {
            pub chrom: String,
            pub start: u64,
            pub end: u64,
            pub support: u32,
            pub introns: Vec<(u64, u64)>,
        }

        /// An already-assembled family member's genomic span (used to EXCLUDE overlapping thin loci).
        #[derive(Clone, Debug)]
        pub struct MemberSpan {
            pub chrom: String,
            pub start: u64,
            pub end: u64,
        }

        /// A rescued copy: the thin locus, the family/member that confirmed it (+ core coverage + orientation), the
        /// transcription strand, and the spliced sequence (so the rescued copy is assignable).
        #[derive(Clone, Debug)]
        pub struct RescuedCopy {
            pub locus: ThinLocus,
            pub outcome: RescueOutcome,
            pub strand: char,
            pub seq: Vec<u8>,
        }

        /// Group reads by exact intron chain, then collapse OVERLAPPING chains (by span) into loci, keeping the
        /// best-supported chain (by support, then span) as each locus's representative. Multi-exon only. Mirrors
        /// `family_rescue.py`'s per-interval locus collapse. Deterministic (sorted by start, ties by chain).
        pub fn thin_loci(reads: &[PrimaryRead], min_support: u32) -> Vec<ThinLocus> {
            // group reads by (chrom, exact intron chain) -> (support, min_start, max_end); multi-exon only.
            let mut chains: BTreeMap<(&str, Vec<(u64, u64)>), (u32, u64, u64)> = BTreeMap::new();
            for r in reads {
                if r.introns.is_empty() {
                    continue;
                }
                let e = chains
                    .entry((r.chrom.as_str(), r.introns.clone()))
                    .or_insert((0, u64::MAX, 0));
                e.0 += 1;
                e.1 = e.1.min(r.ref_start);
                e.2 = e.2.max(r.ref_end);
            }
            // per chrom: collapse overlapping chains (single-linkage by span); the rep is the best (support, span).
            let mut by_chrom: BTreeMap<&str, Vec<(Vec<(u64, u64)>, u32, u64, u64)>> =
                BTreeMap::new();
            for ((chrom, chain), (sup, s, e)) in chains {
                by_chrom.entry(chrom).or_default().push((chain, sup, s, e));
            }
            struct Loc {
                s: u64,
                e: u64,
                sup: u32,
                s2: u64,
                e2: u64,
                intr: Vec<(u64, u64)>,
            }
            let mut out = Vec::new();
            for (chrom, mut group) in by_chrom {
                group.sort_by_key(|&(_, _, s, _)| s); // stable: ties keep (chrom,chain) order
                let mut loci: Vec<Loc> = Vec::new();
                for (intr, sup, s, e) in group {
                    if sup < min_support {
                        continue;
                    }
                    let mut merged = false;
                    for l in loci.iter_mut() {
                        if s <= l.e && e >= l.s {
                            l.s = l.s.min(s);
                            l.e = l.e.max(e);
                            if (sup, e - s) > (l.sup, l.e2 - l.s2) {
                                l.sup = sup;
                                l.s2 = s;
                                l.e2 = e;
                                l.intr = intr.clone();
                            }
                            merged = true;
                            break;
                        }
                    }
                    if !merged {
                        loci.push(Loc {
                            s,
                            e,
                            sup,
                            s2: s,
                            e2: e,
                            intr,
                        });
                    }
                }
                for l in loci {
                    out.push(ThinLocus {
                        chrom: chrom.to_string(),
                        start: l.s2,
                        end: l.e2,
                        support: l.sup,
                        introns: l.intr,
                    });
                }
            }
            out
        }

        /// Rescue under-assembled copies: for each thin locus NOT overlapping a member span, build its spliced
        /// sequence and POA-confirm it against the family `members` (`rescue_thin_locus`). Dedup by locus, keeping
        /// the best `core_recip`. Mirrors `family_rescue.py`'s scan + dedup.
        pub fn rescue_thin_loci(
            loci: &[ThinLocus],
            members: &[FamilyMember],
            member_spans: &[MemberSpan],
            genome: &GenomeIndex,
            p: &RescueParams,
        ) -> Vec<RescuedCopy> {
            let mut by_key: BTreeMap<(String, u64, u64), RescuedCopy> = BTreeMap::new();
            for locus in loci {
                // exclude a thin locus that overlaps an already-assembled family member span.
                if member_spans
                    .iter()
                    .any(|m| m.chrom == locus.chrom && locus.start < m.end && locus.end > m.start)
                {
                    continue;
                }
                let (seq, strand) = match build_spliced_seq(
                    genome,
                    &locus.chrom,
                    locus.start,
                    locus.end,
                    &locus.introns,
                    None,
                ) {
                    Some(v) => v,
                    None => continue,
                };
                if seq.len() < RESCUE_MIN_LEN || seq.len() > p.len_cap {
                    continue;
                }
                if let Some(outcome) = rescue_thin_locus(&seq, members, p) {
                    let key = (locus.chrom.clone(), locus.start, locus.end);
                    let better = by_key
                        .get(&key)
                        .map_or(true, |prev| outcome.core_recip > prev.outcome.core_recip);
                    if better {
                        by_key.insert(
                            key,
                            RescuedCopy {
                                locus: locus.clone(),
                                outcome,
                                strand,
                                seq,
                            },
                        );
                    }
                }
            }
            by_key.into_values().collect()
        }

        /// ITERATIVE family-aware rescue (borrow strength across passes). A rescued copy is itself a new family
        /// member that can bridge to OTHER under-assembled copies the first pass couldn't reach (homologous to the
        /// rescued copy but not the original family). So: rescue, fold the rescued copies in as members + spans,
        /// rescue the remaining loci again, until a pass recovers nothing new (or `max_iters`). Deterministic.
        pub fn rescue_thin_loci_iterative(
            loci: &[ThinLocus],
            members: &[FamilyMember],
            member_spans: &[MemberSpan],
            genome: &GenomeIndex,
            p: &RescueParams,
            max_iters: usize,
        ) -> Vec<RescuedCopy> {
            let mut all_members: Vec<FamilyMember> = members.to_vec();
            let mut all_spans: Vec<MemberSpan> = member_spans.to_vec();
            let mut remaining: Vec<ThinLocus> = loci.to_vec();
            let mut rescued_all: Vec<RescuedCopy> = Vec::new();
            for _ in 0..max_iters {
                let rescued = rescue_thin_loci(&remaining, &all_members, &all_spans, genome, p);
                if rescued.is_empty() {
                    break;
                }
                // each rescued copy becomes a member (so it can bridge) and a span (so it is not re-rescued).
                let done: BTreeSet<(String, u64, u64)> = rescued
                    .iter()
                    .map(|rc| (rc.locus.chrom.clone(), rc.locus.start, rc.locus.end))
                    .collect();
                for rc in &rescued {
                    all_members.push(FamilyMember::new(
                        format!("RC_{}_{}", rc.locus.chrom, rc.locus.start),
                        rc.outcome.family_id.clone(),
                        rc.seq.clone(),
                    ));
                    all_spans.push(MemberSpan {
                        chrom: rc.locus.chrom.clone(),
                        start: rc.locus.start,
                        end: rc.locus.end,
                    });
                }
                remaining.retain(|l| !done.contains(&(l.chrom.clone(), l.start, l.end)));
                rescued_all.extend(rescued);
                if remaining.is_empty() {
                    break;
                }
            }
            rescued_all
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
            fn read(chrom: &str, s: u64, e: u64, introns: &[(u64, u64)]) -> PrimaryRead {
                PrimaryRead {
                    chrom: chrom.into(),
                    ref_start: s,
                    ref_end: e,
                    introns: introns.to_vec(),
                    reverse: false,
                }
            }

            /// Genome with one thin gene: exon1 [0,200), canonical intron [200,220) GT..AG, exon2 [220,420);
            /// spliced = flank(50) + `core`(300) + flank(50). Returns the genome.
            fn thin_gene_genome(core: &[u8]) -> GenomeIndex {
                let mut g = vec![b'A'; 500];
                g[0..50].copy_from_slice(&rand_seq(50, 0x71));
                g[50..200].copy_from_slice(&core[0..150]);
                g[200] = b'G';
                g[201] = b'T';
                g[218] = b'A';
                g[219] = b'G';
                g[220..370].copy_from_slice(&core[150..300]);
                g[370..420].copy_from_slice(&rand_seq(50, 0x72));
                GenomeIndex::from_seqs(&[("c1", &g)])
            }

            // ---- thin_loci ----

            #[test]
            fn thin_loci_groups_and_keeps_multi_exon() {
                let reads = [
                    read("c1", 0, 400, &[(100, 200)]),
                    read("c1", 5, 410, &[(100, 200)]),
                    read("c1", 0, 400, &[]), // single-exon -> ignored
                ];
                let loci = thin_loci(&reads, 1);
                assert_eq!(loci.len(), 1);
                assert_eq!(loci[0].support, 2);
                assert_eq!(loci[0].introns, vec![(100, 200)]);
            }

            #[test]
            fn thin_loci_collapses_overlapping_chains_best_supported() {
                // chain A [0,400) support 2, chain B [100,500) support 1 overlap -> one locus, rep = chain A.
                let reads = [
                    read("c1", 0, 400, &[(100, 200)]),
                    read("c1", 0, 400, &[(100, 200)]),
                    read("c1", 100, 500, &[(150, 250)]),
                ];
                let loci = thin_loci(&reads, 1);
                assert_eq!(loci.len(), 1, "overlapping chains collapse to one locus");
                assert_eq!(loci[0].support, 2, "best-supported chain wins");
                assert_eq!(loci[0].introns, vec![(100, 200)]);
                assert_eq!((loci[0].start, loci[0].end), (0, 400));
            }

            #[test]
            fn thin_loci_separates_nonoverlapping() {
                let reads = [
                    read("c1", 0, 300, &[(100, 200)]),
                    read("c1", 1000, 1300, &[(1100, 1200)]),
                ];
                assert_eq!(thin_loci(&reads, 1).len(), 2);
            }

            // ---- rescue_thin_loci ----

            #[test]
            fn rescue_recovers_thin_locus_homologous_to_family() {
                let core = rand_seq(300, 0xC0FE_F00D);
                let genome = thin_gene_genome(&core);
                // a family member sharing the same core (distinct flanks)
                let mseq = cat(&[&rand_seq(50, 0x81), &core, &rand_seq(50, 0x82)]);
                let members = [FamilyMember::new("M1".into(), "FAM1".into(), mseq)];
                // a SINGLE read (support 1, below the >=3 gate) forming the thin locus
                let reads = [read("c1", 0, 420, &[(200, 220)])];
                let loci = thin_loci(&reads, RESCUE_MIN_SUPPORT);
                assert_eq!(loci.len(), 1);
                let rescued =
                    rescue_thin_loci(&loci, &members, &[], &genome, &RescueParams::default());
                assert_eq!(rescued.len(), 1, "the thin copy is rescued into the family");
                assert_eq!(rescued[0].outcome.family_id, "FAM1");
                assert!(rescued[0].outcome.core_recip >= 0.13);
                assert_eq!(rescued[0].seq.len(), 400);
            }

            #[test]
            fn rescue_excludes_loci_overlapping_a_member_span() {
                let core = rand_seq(300, 0xC0FE_F00D);
                let genome = thin_gene_genome(&core);
                let mseq = cat(&[&rand_seq(50, 0x81), &core, &rand_seq(50, 0x82)]);
                let members = [FamilyMember::new("M1".into(), "FAM1".into(), mseq)];
                let reads = [read("c1", 0, 420, &[(200, 220)])];
                let loci = thin_loci(&reads, RESCUE_MIN_SUPPORT);
                // a member already assembled across the locus span -> the thin locus is NOT a new copy.
                let spans = [MemberSpan {
                    chrom: "c1".into(),
                    start: 100,
                    end: 300,
                }];
                let rescued =
                    rescue_thin_loci(&loci, &members, &spans, &genome, &RescueParams::default());
                assert!(
                    rescued.is_empty(),
                    "locus overlapping an assembled member is excluded"
                );
            }

            #[test]
            fn rescue_rejects_non_homologous_thin_locus() {
                let genome = thin_gene_genome(&rand_seq(300, 0xC0FE_F00D));
                // a family member with a DIFFERENT, unrelated core -> no rescue.
                let mseq = cat(&[
                    &rand_seq(50, 0x81),
                    &rand_seq(300, 0xDEAD_BEEF),
                    &rand_seq(50, 0x82),
                ]);
                let members = [FamilyMember::new("M1".into(), "FAM1".into(), mseq)];
                let reads = [read("c1", 0, 420, &[(200, 220)])];
                let loci = thin_loci(&reads, RESCUE_MIN_SUPPORT);
                let rescued =
                    rescue_thin_loci(&loci, &members, &[], &genome, &RescueParams::default());
                assert!(
                    rescued.is_empty(),
                    "a non-homologous thin locus is not rescued"
                );
            }

            #[test]
            fn iterative_rescue_recovers_a_bridged_copy() {
                // M1 = flank + core1. L1 spliced = core1 + core2 (rescued pass 1 via core1, shared with M1).
                // L2 spliced = core2 + flank (homologous to L1's core2 but NOT to M1) -> rescued ONLY in pass 2,
                // once L1 is a member. Single-pass recovers 1; iterative recovers 2.
                let core1 = rand_seq(200, 0xC0DE_0001);
                let core2 = rand_seq(200, 0xC0DE_0002);
                let mut g = vec![b'A'; 1600];
                // L1 = core1 | core2 ; L2 = flankL | core2 ; M1 = core1 | flankM. Shared cores sit at the SAME
                // relative position (no offset) so POA anchors them cleanly.
                g[0..200].copy_from_slice(&core1);
                g[200] = b'G';
                g[201] = b'T';
                g[218] = b'A';
                g[219] = b'G';
                g[220..420].copy_from_slice(&core2);
                g[1000..1200].copy_from_slice(&rand_seq(200, 0xF00D)); // flankL
                g[1200] = b'G';
                g[1201] = b'T';
                g[1218] = b'A';
                g[1219] = b'G';
                g[1220..1420].copy_from_slice(&core2);
                let genome = GenomeIndex::from_seqs(&[("c1", &g)]);
                let m1 = FamilyMember::new(
                    "M1".into(),
                    "FAM1".into(),
                    cat(&[&core1, &rand_seq(200, 0xAA)]),
                );
                let loci = [
                    ThinLocus {
                        chrom: "c1".into(),
                        start: 0,
                        end: 420,
                        support: 1,
                        introns: vec![(200, 220)],
                    },
                    ThinLocus {
                        chrom: "c1".into(),
                        start: 1000,
                        end: 1420,
                        support: 1,
                        introns: vec![(1200, 1220)],
                    },
                ];
                let single = rescue_thin_loci(
                    &loci,
                    std::slice::from_ref(&m1),
                    &[],
                    &genome,
                    &RescueParams::default(),
                );
                assert_eq!(
                    single.len(),
                    1,
                    "single-pass rescues only the directly-homologous locus"
                );
                let iter = rescue_thin_loci_iterative(
                    &loci,
                    &[m1],
                    &[],
                    &genome,
                    &RescueParams::default(),
                    5,
                );
                assert_eq!(
                    iter.len(),
                    2,
                    "iterative recovers the bridged locus via the first rescued copy"
                );
            }
        }
    }
}

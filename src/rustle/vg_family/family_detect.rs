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

use std::collections::{BTreeMap, BTreeSet};

use super::family_graph::{contiguous_core_coverage_bounded, upper_cow};
use super::family_rescue::window_canon_code;
use crate::vg_family::seq_utils::reverse_complement;

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
                let ka = (transcripts[a].n_reads, transcripts[a].end - transcripts[a].start);
                let kb = (transcripts[b].n_reads, transcripts[b].end - transcripts[b].start);
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
    let mut by_root: std::collections::BTreeMap<usize, Vec<usize>> = std::collections::BTreeMap::new();
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
            (a, nxt.into_iter().filter(|&(_, r)| r >= min_reads as u64).map(|(b, _)| b).collect())
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
    let index: BTreeMap<(u64, u64), usize> = nodes.iter().enumerate().map(|(i, &j)| (j, i)).collect();
    let mut best: Vec<u64> = nodes.iter().map(|j| weight[j]).collect();
    let mut prev: Vec<Option<usize>> = vec![None; nodes.len()];
    for i in 0..nodes.len() {
        let (_, acceptor) = nodes[i];
        let Some(nexts) = succ.get(&nodes[i]) else { continue };
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
    // Opt-in `RUSTLE_SPLICED_REP`: compare SPLICED evidence to UNSPLICED evidence in aggregate before
    // choosing, instead of letting one pooled unspliced cluster outvote each spliced isoform individually.
    // See `bench/soto/merge_quality_analysis.md` §20. Default off -> byte-identical.
    let spliced_rep = matches!(std::env::var("RUSTLE_SPLICED_REP"), Ok(v) if v != "0" && !v.is_empty());
    let mut reps: Vec<usize> = comp
        .into_values()
        .map(|members| pick_locus_rep(transcripts, &members, spliced_rep))
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
/// With `spliced_rep`, spliced chains at a locus are treated as what they are — ISOFORMS of one
/// transcription unit, not competitors — and their reads are summed before the comparison. The winning
/// CLASS then supplies the representative. A genuinely intronless locus has no spliced members, so its
/// unspliced cluster still wins; nothing forces a spliced representative where none exists.
fn pick_locus_rep(transcripts: &[DenovoTranscript], members: &[usize], spliced_rep: bool) -> usize {
    let best = |cands: &[usize]| -> usize {
        *cands
            .iter()
            .max_by(|&&a, &&b| {
                let ka = (transcripts[a].n_reads, transcripts[a].end - transcripts[a].start);
                let kb = (transcripts[b].n_reads, transcripts[b].end - transcripts[b].start);
                // break key ties by preferring the smaller index (so it is the "maximum").
                ka.cmp(&kb).then_with(|| b.cmp(&a))
            })
            .unwrap()
    };
    if spliced_rep {
        let spliced: Vec<usize> = members
            .iter()
            .copied()
            .filter(|&i| !transcripts[i].introns.is_empty())
            .collect();
        if !spliced.is_empty() {
            let spl: u64 = spliced.iter().map(|&i| transcripts[i].n_reads as u64).sum();
            let unspl: u64 = members
                .iter()
                .filter(|&&i| transcripts[i].introns.is_empty())
                .map(|&i| transcripts[i].n_reads as u64)
                .sum();
            if spl >= unspl {
                return best(&spliced);
            }
        }
    }
    best(members)
}

/// Exon intervals implied by a transcript's `(start, end, introns)`.
pub(crate) fn exons_of(t: &DenovoTranscript) -> Vec<(u64, u64)> {
    let mut out = Vec::with_capacity(t.introns.len() + 1);
    let mut prev = t.start;
    for &(d, a) in &t.introns {
        if d > prev { out.push((prev, d)); }
        prev = a;
    }
    if t.end > prev { out.push((prev, t.end)); }
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
            if hi > lo { ov += hi - lo; }
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
    let junction_only = matches!(std::env::var("RUSTLE_LOCUS_JUNCTION_ONLY"), Ok(v) if v != "0" && !v.is_empty());
    if !p.collapse_span_aware || junction_only {
        return parent;
    }

    // Phase 2: iterative span-aware merge of locus representatives.
    //
    // Under `RUSTLE_SPLICED_REP`, an INTRONLESS transcript is treated as having UNKNOWN strand rather than
    // `+`. Strand is derived from canonical junction motifs, so a cluster with no junctions has none to
    // derive from and its `+` is a placeholder. The merge condition `a.strand == b.strand` nevertheless
    // treats that placeholder as evidence, which makes an unspliced cluster unmergeable with any MINUS-strand
    // gene's spliced skeletons. Measured on the Soto catalog: all 740 single-exon copies are `+`, while
    // spliced copies split 149 `+` / 218 `-`. So at a minus-strand gene -- SRGAP2 among them -- the intronic
    // unspliced cluster never even meets the spliced evidence, survives as an independent locus, and is
    // emitted as its own single-exon copy (§14, §20).
    //
    // `RUSTLE_COLLAPSE_UNSTRANDED` UN-WELDS that clause from `RUSTLE_SPLICED_REP`. `RUSTLE_SPLICED_REP`
    // also changes WHICH transcript represents a locus (`pick_locus_rep` above), and it is that leg the
    // end-to-end regression was attributed to (chr7 F1 0.570 -> 0.411, chr16 0.910 -> 0.761). The two
    // effects cannot be separated while one variable gates both, so this one turns on the STRAND clause
    // ALONE, leaving rep selection at its default. Threshold-free in the `is_chimeric_bridge` sense: it
    // asks whether the strand field was MEASURED, not how large a number is. With NEITHER variable set the
    // `else` branch below is reached unchanged, so the OFF arm is byte-identical by construction.
    let unstranded_unspliced =
        matches!(std::env::var("RUSTLE_SPLICED_REP"), Ok(v) if v != "0" && !v.is_empty())
            || matches!(std::env::var("RUSTLE_COLLAPSE_UNSTRANDED"), Ok(v) if v != "0" && !v.is_empty());
    let exonic_collapse =
        matches!(std::env::var("RUSTLE_COLLAPSE_EXONIC"), Ok(v) if v != "0" && !v.is_empty());
    loop {
        let reps = locus_reps(transcripts, &parent);
        let mut merged = false;
        for i in 0..reps.len() {
            let a = &transcripts[reps[i]];
            for j in (i + 1)..reps.len() {
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
                    uf_union(&mut parent, reps[i], reps[j]);
                    merged = true;
                    continue;
                }
                // Disjoint-junction isoforms with similar length: require strong POA core coverage.
                let au = upper_cow(&a.seq);
                let bu = upper_cow(&b.seq);
                let core = contiguous_core_coverage_bounded(&au, &bu, p.len_cap);
                if core >= p.collapse_span_core {
                    uf_union(&mut parent, reps[i], reps[j]);
                    merged = true;
                }
            }
        }
        if !merged {
            break;
        }
    }

    parent
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
    let mut by_root: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
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
    let mut sum_by_root: std::collections::HashMap<usize, u32> = std::collections::HashMap::new();
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
    use super::family_graph::{contiguous_core_coverage_bounded_with, EDGE_CONFIRM_ASTAR};
    let au = upper_cow(a);
    let bu = upper_cow(b);
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
    use std::collections::HashMap;
    let mut in_fam: HashMap<usize, usize> = HashMap::new();
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
        let fwd = crate::vg_family::family_graph::longest_common_substring(&au, &bu);
        let rev = crate::vg_family::family_graph::longest_common_substring(&au, &reverse_complement(&bu));
        fwd.max(rev) as f64 / minlen as f64
    };
    // best (family, core) per FREE rep over all family-adjacent candidate pairs that confirm a POA core.
    let mut best: HashMap<usize, (usize, f64)> = HashMap::new();
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
         ..Default::default() }
    }
    #[test]
    fn cothread_chain_spans_more_junctions_than_any_single_transcript() {
        // The shattering case: no observed chain carries the whole locus, but consecutive junctions are
        // co-observed pairwise, so the graph path recovers all four.
        let t = vec![
            tx("a", "c", 100, 900, 5, &[(200, 300), (400, 500)], vec![b'A'; 10]),
            tx("b", "c", 100, 900, 4, &[(400, 500), (600, 700)], vec![b'A'; 10]),
            tx("c", "c", 100, 900, 4, &[(600, 700), (750, 800)], vec![b'A'; 10]),
        ];
        let (s, e, chain) = cothread_locus_geometry(&t, &[0, 1, 2], 3).unwrap();
        assert_eq!((s, e), (100, 900));
        assert_eq!(chain, vec![(200, 300), (400, 500), (600, 700), (750, 800)],
                   "path threads all four junctions though no transcript carries more than two");
    }

    #[test]
    fn cothread_edges_obey_the_same_read_floor_as_nodes() {
        // A readthrough carrying 2 reads witnesses the succession into the neighbouring gene. Both
        // junctions clear the NODE floor on their own support, so only an EDGE floor can stop the walk.
        let t = vec![
            tx("here", "c", 100, 900, 40, &[(200, 300)], vec![b'A'; 10]),
            tx("nbr", "c", 5000, 6000, 40, &[(5200, 5300)], vec![b'A'; 10]),
            tx("readthrough", "c", 100, 6000, 2, &[(200, 300), (5200, 5300)], vec![b'A'; 10]),
        ];
        let (s, e, chain) = cothread_locus_geometry(&t, &[0, 1, 2], 3).unwrap();
        assert_eq!(chain.len(), 1, "a 2-read readthrough must not license the cross-gene step");
        assert!(e - s < 5000, "extent must not run into the neighbour, got {s}-{e}");
        // With enough reads behind it the succession is real and IS taken.
        let t2 = vec![
            tx("here", "c", 100, 900, 40, &[(200, 300)], vec![b'A'; 10]),
            tx("both", "c", 100, 6000, 9, &[(200, 300), (5200, 5300)], vec![b'A'; 10]),
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
        assert_eq!(chain.len(), 1, "no read witnesses the succession, so no chaining across genes");
        assert!(e - s <= 300, "extent stays at the contributing transcript, got {s}-{e}");
    }

    #[test]
    fn cothread_applies_the_junction_floor_and_abstains_when_nothing_clears_it() {
        let t = vec![
            tx("weak", "c", 100, 900, 2, &[(200, 300), (400, 500)], vec![b'A'; 10]),
            tx("stub", "c", 100, 900, 50, &[], vec![b'A'; 10]),
        ];
        assert!(cothread_locus_geometry(&t, &[0, 1], 3).is_none(),
                "2 reads/junction is below the floor; caller must fall back to today's representative");
        // The same junctions, now carried by enough reads, do produce a chain.
        let t2 = vec![tx("ok", "c", 100, 900, 3, &[(200, 300), (400, 500)], vec![b'A'; 10])];
        assert_eq!(cothread_locus_geometry(&t2, &[0], 3).unwrap().2,
                   vec![(200, 300), (400, 500)]);
    }

    #[test]
    fn cothread_prefers_the_heavier_path_when_two_compete() {
        // Both successions are read-witnessed from (200,300); the heavier one wins.
        let t = vec![
            tx("light", "c", 100, 900, 3, &[(200, 300), (400, 500)], vec![b'A'; 10]),
            tx("heavy", "c", 100, 900, 30, &[(200, 300), (600, 700)], vec![b'A'; 10]),
        ];
        let chain = cothread_locus_geometry(&t, &[0, 1], 3).unwrap().2;
        assert_eq!(chain, vec![(200, 300), (600, 700)]);
    }

    #[test]
    fn union_geometry_merges_partial_chains_into_the_whole_locus() {
        // Two chains covering opposite halves of one locus -- the shattering case. Either alone is a
        // fragment; the union is the locus.
        let t = vec![
            tx("t100_5", "c", 100, 400, 5, &[(200, 300)], vec![b'A'; 10]),   // exons 100-200, 300-400
            tx("t500_4", "c", 500, 900, 4, &[(600, 700)], vec![b'A'; 10]),   // exons 500-600, 700-900
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
        let t = vec![tx("t100_5", "c", 100, 200, 5, &[], vec![b'A'; 10]), tx("t900_2", "c", 900, 1000, 2, &[], vec![b'A'; 10])];
        let (s, e, introns) = union_locus_geometry(&t, &[0, 1], 3).unwrap();
        assert_eq!((s, e), (100, 200), "the 2-read chain must not extend the locus");
        assert!(introns.is_empty());
    }

    #[test]
    fn union_geometry_returns_none_when_no_chain_clears_the_floor() {
        // Whole locus is lost rather than silently rebuilt from sub-threshold evidence. This is the
        // measured cost of the substrate (Soto 299 -> 263 loci).
        let t = vec![tx("t100_1", "c", 100, 200, 1, &[], vec![b'A'; 10]), tx("t150_2", "c", 150, 250, 2, &[], vec![b'A'; 10])];
        assert!(union_locus_geometry(&t, &[0, 1], 3).is_none());
    }

    #[test]
    fn union_geometry_of_a_single_chain_reproduces_that_chain() {
        // With one chain above the floor the union must be a no-op, so enabling the substrate cannot
        // perturb loci that never shattered.
        let t = vec![tx("t100_9", "c", 100, 400, 9, &[(200, 300)], vec![b'A'; 10])];
        let (s, e, introns) = union_locus_geometry(&t, &[0], 3).unwrap();
        assert_eq!((s, e), (t[0].start, t[0].end));
        assert_eq!(introns, t[0].introns);
    }

    #[test]
    fn union_geometry_absorbs_a_retained_intron_isoform() {
        // One chain splices 200-300, another retains it. The union treats retained sequence as exonic,
        // which is the point: include everything the reads support.
        let t = vec![tx("t100_6", "c", 100, 400, 6, &[(200, 300)], vec![b'A'; 10]), tx("t100_4", "c", 100, 400, 4, &[], vec![b'A'; 10])];
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
            chrom: "c1".into(), start: 0, end: 10_000, strand: '+',
            introns: vec![(500, 9_500)], seq: b"AC".to_vec(), ..Default::default()
        };
        let inner = DenovoTranscript {
            chrom: "c1".into(), start: 4_000, end: 5_000, strand: '+',
            introns: vec![], seq: b"AC".to_vec(), ..Default::default()
        };
        // span containment: the inner locus is fully inside the long model's span
        let span_ov = long.end.min(inner.end).saturating_sub(long.start.max(inner.start));
        assert_eq!(span_ov, 1_000);
        let span_minlen = (long.end - long.start).min(inner.end - inner.start);
        assert_eq!(span_ov as f64 / span_minlen as f64, 1.0, "span rule says fully contained -> would merge");

        // exonic containment: zero shared exonic sequence
        let (eov, la, lb) = exonic_overlap(&long, &inner);
        assert_eq!(eov, 0, "no shared EXONIC sequence");
        assert_eq!(la, 1_000, "long model's exons are 0-500 and 9500-10000");
        assert_eq!(lb, 1_000);

        // and it still merges genuine isoforms that share exonic sequence
        let iso = DenovoTranscript {
            chrom: "c1".into(), start: 0, end: 400, strand: '+',
            introns: vec![], seq: b"AC".to_vec(), ..Default::default()
        };
        let (eov2, _, lb2) = exonic_overlap(&long, &iso);
        assert_eq!(eov2, 400);
        assert!(eov2 as f64 / lb2 as f64 >= COLLAPSE_CONTAIN_FRAC, "real isoform still merges");
    }

    #[test]
    fn spliced_rep_beats_a_read_richer_unspliced_pool_but_only_when_spliced_evidence_wins() {
        // The SRGAP2C shape: three spliced isoform-chains of 20/16/10 reads (46 total) at one locus, plus a
        // single pooled unspliced cluster of 40. The default rule takes the 40-read unspliced cluster
        // because it beats every individual chain; the aggregate rule sees 46 > 40 and takes the best chain.
        let mk = |i: usize, reads: u32, introns: Vec<(u64, u64)>| DenovoTranscript {
            tid: format!("t{i}"), chrom: "c1".into(), start: 100, end: 900,
            n_reads: reads, strand: '+', introns, seq: b"ACGT".to_vec(), ..Default::default()
        };
        let tx = vec![
            mk(0, 20, vec![(200, 300)]),
            mk(1, 16, vec![(200, 400)]),
            mk(2, 10, vec![(250, 300)]),
            mk(3, 40, vec![]),            // pooled unspliced cluster — read-richest single member
        ];
        let members: Vec<usize> = (0..4).collect();
        assert_eq!(pick_locus_rep(&tx, &members, false), 3, "default rule takes the unspliced pool");
        assert_eq!(pick_locus_rep(&tx, &members, true), 0, "aggregate rule takes the best spliced chain");

        // A genuinely INTRONLESS locus has no spliced member, so nothing is forced: the unspliced cluster
        // still represents it. This is what keeps histones and processed pseudogenes correct.
        let only_unspliced = vec![mk(0, 5, vec![]), mk(1, 30, vec![])];
        assert_eq!(pick_locus_rep(&only_unspliced, &[0, 1], true), 1);

        // When unspliced evidence genuinely dominates (40 vs 12), the spliced side does NOT win -- the rule
        // is a fair aggregate comparison, not an unconditional preference for splicing.
        let weak_spliced = vec![mk(0, 12, vec![(200, 300)]), mk(1, 40, vec![])];
        assert_eq!(pick_locus_rep(&weak_spliced, &[0, 1], true), 1);
    }

    #[test]
    fn collapse_loci_groups_maps_isoforms_to_their_gene_rep() {
        // two isoforms of gene A share the junction (100,200); a third transcript at a disjoint locus is its
        // own gene. groups[i] must be the rep index of i's gene (rep = most reads, here A_iso1 with 10).
        let txs = vec![
            tx("A_iso1", "c1", 0, 300, 10, &[(100, 200)], vec![b'A'; 50]),
            tx("A_iso2", "c1", 0, 400, 4, &[(100, 200), (250, 320)], vec![b'A'; 60]),
            tx("B_iso1", "c1", 1000, 1300, 8, &[(1100, 1200)], vec![b'C'; 50]),
        ];
        let g = collapse_loci_groups(&txs);
        assert_eq!(g[0], g[1], "isoforms sharing a junction collapse to one gene");
        assert_ne!(g[0], g[2], "a disjoint locus is its own gene");
        assert_eq!(g[0], 0, "gene rep = the higher-read isoform (A_iso1)");
        assert_eq!(g[2], 2, "B is its own rep");
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
            tx("long", "c1", 100, 800, 8, &[(200, 300), (500, 600)], long_seq),
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
            tx("para1", "c1", 1000, 11000, 10, &[(2000, 3000), (4000, 5000)], rand_seq(1000, 0xC1)),
            tx("para2", "c1", 7000, 17000, 9, &[(8000, 9000), (10000, 11000)], rand_seq(1000, 0xC2)),
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
        const VARS: [&str; 3] =
            ["RUSTLE_SPLICED_REP", "RUSTLE_COLLAPSE_UNSTRANDED", "RUSTLE_COLLAPSE_EXONIC"];
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
        DenovoTranscript { strand, ..tx(tid, chrom, start, end, n_reads, introns, seq) }
    }

    #[test]
    fn collapse_unstranded_clause_merges_intronless_stub_into_minus_strand_gene() {
        // THE CLAUSE. An INTRONLESS model has no canonical junction motif to derive a strand from, so its
        // `'+'` is a placeholder, not a measurement. Strictly enclosed inside a spliced `'-'` gene
        // (containment 1.0, same chrom), the ONLY thing separating the two is that placeholder.
        let gene = tx_strand("gene", "c1", 1_000, 9_000, 8, '-',
            &[(2_000, 3_000), (5_000, 6_000)], rand_seq(900, 0xD1));
        let stub = tx_strand("stub", "c1", 4_000, 4_600, 4, '+', &[], rand_seq(600, 0xD2));
        let txs = [gene, stub];
        let p = DetectParams::default();

        let off = with_collapse_env(&[], || collapse_loci_span_aware(&txs, &p).len());
        assert_eq!(off, 2, "default: the placeholder '+' is treated as evidence and blocks the merge");

        let on = with_collapse_env(&[("RUSTLE_COLLAPSE_UNSTRANDED", "1")], || {
            collapse_loci_span_aware(&txs, &p).len()
        });
        assert_eq!(on, 1, "RUSTLE_COLLAPSE_UNSTRANDED: an UNMEASURED strand cannot block the merge");
    }

    #[test]
    fn collapse_unstranded_keeps_genuine_antisense_spliced_pairs_apart_in_both_arms() {
        // THE VALUE, pinned separately from the clause. When BOTH sides carry a junction-derived strand
        // the field WAS measured, so opposite strands are real antisense and must stay two loci whether or
        // not the variable is set. This is the mirror of the genuine-antisense cases the un-weld must not
        // touch.
        let minus = tx_strand("minus", "c1", 1_000, 9_000, 8, '-',
            &[(2_000, 3_000), (5_000, 6_000)], rand_seq(900, 0xE1));
        let plus = tx_strand("plus", "c1", 3_500, 4_500, 4, '+',
            &[(3_800, 4_000)], rand_seq(700, 0xE2));
        let anti = [minus.clone(), plus.clone()];
        let p = DetectParams::default();

        let off = with_collapse_env(&[], || collapse_loci_span_aware(&anti, &p).len());
        let on = with_collapse_env(&[("RUSTLE_COLLAPSE_UNSTRANDED", "1")], || {
            collapse_loci_span_aware(&anti, &p).len()
        });
        assert_eq!((off, on), (2, 2), "two MEASURED opposite strands stay two loci in BOTH arms");

        // Control: the strand is the SOLE blocker above. Flip the enclosed model to '-' -- nothing else
        // changes -- and the same pair merges in both arms, so the assertion above is not passing because
        // containment or chrom happened to fail.
        let same = [minus, DenovoTranscript { strand: '-', ..plus }];
        let ctrl_off = with_collapse_env(&[], || collapse_loci_span_aware(&same, &p).len());
        let ctrl_on = with_collapse_env(&[("RUSTLE_COLLAPSE_UNSTRANDED", "1")], || {
            collapse_loci_span_aware(&same, &p).len()
        });
        assert_eq!((ctrl_off, ctrl_on), (1, 1), "same-strand control: containment and chrom do pass");
    }

    #[test]
    fn collapse_with_totals_sums_reads_over_the_whole_locus() {
        // Two isoforms of one gene share junction (100,200) -> one locus; a third transcript is a disjoint
        // locus. The locus total is the SUM of every isoform's reads (32+20=52), not the rep isoform's 32.
        let txs = vec![
            tx("A_iso1", "c1", 0, 300, 32, &[(100, 200)], vec![b'A'; 50]),
            tx("A_iso2", "c1", 0, 400, 20, &[(100, 200), (250, 320)], vec![b'A'; 60]),
            tx("B_iso1", "c1", 1000, 1300, 10, &[(1100, 1200)], vec![b'C'; 50]),
        ];
        let (reps, totals) = collapse_loci_span_aware_with_totals(&txs, &DetectParams::default());
        // reps align with the behavior-preserving wrapper.
        assert_eq!(reps, collapse_loci_span_aware(&txs, &DetectParams::default()));
        assert_eq!(reps.len(), 2, "two loci: gene A and gene B");
        // rep of gene A is the max-reads isoform (A_iso1, idx 0); its locus total is 32+20=52.
        let a_pos = reps.iter().position(|&r| r == 0).expect("A_iso1 is a rep");
        assert_eq!(totals[a_pos], 52, "locus total = sum over all isoforms of the gene");
        let b_pos = reps.iter().position(|&r| r == 2).expect("B_iso1 is a rep");
        assert_eq!(totals[b_pos], 10, "single-isoform locus total = its own reads");
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
        assert_eq!(adds, vec![vec![2usize]], "r2 attaches via the conserved core; r3 does not");
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
    fn confirm_edge_over_cap_uses_fallback_not_skip() {
        // over the poasta cap, confirm_edge switches to the linear-memory LCS fallback instead of skipping.
        // DISJOINT random sequences still fail the t_core bar under the fallback -> None (low coverage), NOT
        // an OOM and NOT a blind drop.
        let p = DetectParams { len_cap: 100, ..DetectParams::default() };
        let a = rand_seq(560, 0x1);
        let b = rand_seq(560, 0x2);
        assert!(confirm_edge(&a, &b, &p).is_none(), "disjoint over-cap pair -> low fallback coverage -> None");
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
        let (edges, fb) = detect_edges_reporting(&reps, &DetectParams::default());
        assert_eq!(edges.len(), 1, "homologous pair confirmed under normal cap");
        assert!(fb.is_empty(), "no fallback used under normal cap");
        // the .0 list still matches plain detect_edges (faithfulness of the delegation).
        assert_eq!(edges, detect_edges(&reps, &DetectParams::default()));

        let p = DetectParams { len_cap: 5, ..DetectParams::default() };
        let (edges2, fb2) = detect_edges_reporting(&reps, &p);
        assert_eq!(edges2.len(), 1, "fallback still confirms the homologous pair (no OOM, no loss)");
        assert_eq!((edges2[0].0, edges2[0].1), (0, 1));
        assert_eq!(fb2, vec![(0, 1)], "the fallback-confirmed pair is reported for audit");
    }
}

#!/usr/bin/env python3
"""Full read-coherence way, genome-wide: from the de-novo assembled transcripts (annotation-free), build
the multi-copy FAMILY layer. Minimizer-FREE end to end.

Pipeline:
  1. collapse isoforms -> one representative per GENE LOCUS by SHARED INTRON JUNCTIONS (union-find on
     identical (donor,acceptor) coords). NOT raw span overlap -- dense GGO genes overlap transitively and
     a span-merge chains a whole chromosome into one bogus "locus" (the old bug: 93k tx -> 27 reps).
  2. COUNTING-BLOOM pre-filter ("only POA what is likely a family member"): one linear pass deposits every
     rep's k-mers into a counting Bloom filter; a k-mer owned by >=2 (and <= CNT_MAX) distinct reps is
     "family-informative". A rep whose informative-k-mer fraction >= TAU is a candidate family member;
     single-copy genes (unique k-mers) are rejected here and never reach POA. This also fixes the prior OOM
     (the old all-k-mer inverted index exploded on the millions of UNIQUE single-copy k-mers; here the
     index is built ONLY over family-informative k-mers, a small bounded set).
  3. POA contiguous-core grading is THE criterion (core_recip >= T_CORE). The Bloom is a recall-safe
     pre-screen only -- it can over-include (costs an extra POA) but never decides membership.

Frames the result as a GENERAL assembler (all transcripts) with a specialized family/copy layer on the
multi-copy subset. Run with /home/juanfra/miniforge3/bin/python (Bio.Align + numpy). Parallel POA.

⚠ PROVENANCE / SCOPE (read before citing the catalog). This genome-wide catalog (denovo_families.tsv)
is produced by the SIMILARITY criterion core_recip >= T_CORE (= 0.13), an ARBITRARY threshold. It is
NOT the principled de-tie READ-CONFLICT-GRAPH family definition (the thesis's interest-I criterion,
family = clique / connected-component of the read-conflict graph, no similarity threshold). That
conflict definition exists and is TDD'd, but currently runs only PER-REGION
(src/rustle/vg_family/denovo_pipeline.rs detect_and_assign over local windows; family_def_readconflict.tsv
= a 12-pair labelled panel only) and has NEVER been run genome-wide. Consequence: this 0.13 catalog still
contains the over-merges the conflict criterion is meant to remove (e.g. DNFAM0 = 728 members spanning
chr1..chrY). When reporting family numbers, state which criterion produced them; the de-tie genome-wide
catalog is an open TODO.
"""
import argparse
import multiprocessing as mp
import os
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import poa_family_definition as P

FA = "/home/juanfra/winloci_scratch/denovo_transcripts.fa"
META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
SKEL = "/home/juanfra/winloci_scratch/denovo_skeletons.tsv"
OUT = os.path.join(os.path.dirname(__file__), "denovo_families.tsv")
EDGES = os.path.join(os.path.dirname(__file__), "denovo_family_edges.tsv")
GENES = os.path.join(os.path.dirname(__file__), "gene_meta_strand.tsv")

# --- POA decision (THE criterion; principled, validated) ---
T_CORE = 0.13      # contiguous-core reciprocal-coverage threshold
LEN_CAP = 20000    # skip POA on pairs where the shorter seq exceeds this (cost guard)

# --- exact k-mer counting pre-filter (the "bloom filter" idea, done exactly: cheap pre-screen so only
#     likely family members reach POA; recall-safe -- it NEVER decides membership, POA does) ---
# KMER=18 gives a 4^18 ~ 6.9e10 k-mer space vs ~47M k-mers here (~0.07% load), so coincidental sharing is
# negligible -- two genes sharing >= K_SHARE exact 18-mers is essentially always real homology, and a
# single-copy gene has < K_SHARE informative k-mers so it is truly rejected (KMER=13/16 were too saturated
# -> every gene looked "shared"). Counts are EXACT (np.unique distinct-rep ownership), no Bloom collisions.
KMER = 18
CNT_MIN = 2        # a k-mer must be owned by >= this many distinct gene-reps to be "family-informative"
CNT_MAX = 40       # drop pervasive k-mers owned by > this many reps (mobile-element/repeat, not family-spec)
PAIR_CAP = 40      # skip PAIR generation from a k-mer owned by > this many reps (cost guard; consistent
                   #   with CNT_MAX -- big families (> this) are rare and logged)
K_SHARE = 6        # propose a candidate pair only if the two reps co-own >= this many informative k-mers
K_SHARE_RESCUE = 2 # relaxed k-mer sharing threshold for pairs involving an annotation-rescued rep
MAX_PAIRS = 8_000_000  # hard guard on the distinct-pair set (prevents any OOM; logged if hit)

# For annotation-aware split: a transcript is assigned to its dominant RefSeq gene only if the gene covers
# at least this fraction of the transcript's spliced exon span.  This avoids mis-assigning a transcript to
# a gene it merely grazes, while still splitting readthroughs that span two distinct genes.
GENE_ASSIGN_FRAC = 0.5

# For graph-aware readthrough split: an exon is assigned to a RefSeq gene only if the gene covers at least
# this fraction of the exon's length.  Exons with no confident assignment are dropped from the split.
EXON_GENE_ASSIGN_FRAC = 0.5
# Additional span gate: only split a transcript if >=2 distinct genes each cover at least this fraction of
# the transcript's total spliced exon span.  This prevents cutting transcripts that merely graze a neighbor.
GENE_SPAN_SPLIT_FRAC = 0.25

_CODE = np.full(256, -1, dtype=np.int64)
for _i, _b in enumerate(b"ACGT"):
    _CODE[_b] = _i
    _CODE[_b + 32] = _i  # lowercase


def _horner(codes, n):
    h = np.zeros(n, dtype=np.int64)
    for t in range(KMER):                       # Horner: h = h*4 + base, vectorized over all windows
        h = (h << 2) + codes[t:t + n]
    return h


def kmer_hashes(seq: str):
    """CANONICAL base-4 codes of every KMER-mer of seq (min of the k-mer and its reverse-complement, so a
    sequence and its RC yield the SAME codes -> family detection is strand-symmetric; paralogs assembled on
    opposite strands still share k-mers). Windows containing N are dropped. KMER<=31 fits int64. Minimizer-
    free. Vectorized Horner rolling hash -- NOT an int64 matmul (numpy integer matmul is a slow non-BLAS loop)."""
    codes = _CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    L = codes.size
    if L < KMER:
        return np.empty(0, dtype=np.int64)
    n = L - KMER + 1
    bad = codes < 0
    safe = np.where(bad, 0, codes)
    fwd = _horner(safe, n)
    rc = _horner((3 - safe)[::-1], n)[::-1]      # rc[i] = hash of the reverse-complement of window i
    canon = np.minimum(fwd, rc)
    cb = np.concatenate(([0], np.cumsum(bad.astype(np.int64))))
    badwin = (cb[KMER:] - cb[:n]) > 0            # drop any window touching an N
    return canon[~badwin]


# ---------- POA grading worker (parallel) ----------
_S = None
_A = None


def _init():
    global _S, _A
    _S = P.load_fasta(FA)
    _A = P.make_aligner()


_RC = str.maketrans("ACGTN", "TGCAN")


def _poa(t):
    a, b = t
    sa, sb = _S.get(a), _S.get(b)
    if not sa or not sb or min(len(sa), len(sb)) > LEN_CAP:
        return None
    cr = P.poa_pair_stats(_A, sa, sb)["core_recip"]
    if cr < T_CORE:                              # try opposite orientation (copies assembled on diff strands)
        cr = max(cr, P.poa_pair_stats(_A, sa, sb.translate(_RC)[::-1])["core_recip"])
    return (a, b, round(cr, 4)) if cr >= T_CORE else None


def load_genes(path):
    """Load RefSeq gene coordinates from gene_meta_strand.tsv."""
    genes = []
    if not os.path.exists(path):
        return genes
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            chrom, s, e, strand, name, biotype = parts[:6]
            genes.append((chrom, int(s), int(e), strand, name, biotype))
    return genes


def _exon_intervals(tid, meta, skel_introns):
    c, s, e, _strand, _ne, _nr = meta[tid]
    introns = skel_introns.get((c, s, e), [])
    exons = []
    cur = s
    for d, a in introns:
        exons.append((cur, d))      # half-open [cur, d)
        cur = a
    exons.append((cur, e + 1))      # last exon half-open [cur, e+1)
    return exons


def _exon_len(exons):
    return sum(ee - es for es, ee in exons)


def dominant_gene(tid, meta, skel_introns, genes, min_frac=GENE_ASSIGN_FRAC):
    """Return the RefSeq gene that covers >= min_frac of tid's spliced exon span, if any."""
    c, _s, _e, _strand, _ne, _nr = meta[tid]
    exons = _exon_intervals(tid, meta, skel_introns)
    total = _exon_len(exons)
    if total == 0:
        return None
    best_name = None
    best_ov = 0
    for gc, gs, ge, gstrand, name, _bio in genes:
        if gc != c:
            continue
        ov = 0
        for es, ee in exons:
            ov += max(0, min(ee, ge) - max(es, gs))
        if ov > best_ov:
            best_ov = ov
            best_name = name
    if best_name is not None and best_ov >= min_frac * total:
        return best_name
    return None


def _gene_of_exon(exon, genes, min_frac=EXON_GENE_ASSIGN_FRAC):
    """Return the RefSeq gene that best overlaps a single exon (s,e), if it covers >= min_frac of the exon."""
    es, ee = exon
    best_name = None
    best_ov = 0
    for gc, gs, ge, _gstrand, name, _bio in genes:
        if ee <= gs or ge <= es:
            continue
        ov = min(ee, ge) - max(es, gs)
        if ov > best_ov:
            best_ov = ov
            best_name = name
    if best_name is not None and best_ov >= min_frac * (ee - es):
        return best_name
    return None


def _split_transcript_by_genes(tid, meta, skel_introns, seqs, genes):
    """Split a readthrough transcript into gene-pure sub-transcripts.

    Uses the exon path of the transcript and the RefSeq gene annotations.  Consecutive exons assigned to
    the same gene become one sub-transcript; unassigned exons are skipped.  Returns a list of
    (new_tid, new_meta_tuple, new_seq, new_introns) for every gene-pure segment.
    """
    c, s, e, strand, ne, nr = meta[tid]
    introns = skel_introns.get((c, s, e), [])
    exons = _exon_intervals(tid, meta, skel_introns)
    seq = seqs[tid]
    total_exon_len = _exon_len(exons)
    if total_exon_len == 0:
        return []

    # Span gate: only split if >=2 distinct genes each cover >= GENE_SPAN_SPLIT_FRAC of spliced exon span.
    gene_span = defaultdict(int)
    for gc, gs, ge, _gstrand, name, _bio in genes:
        if gc != c:
            continue
        for es, ee in exons:
            gene_span[name] += max(0, min(ee, ge) - max(es, gs))
    strong_genes = [g for g, ov in gene_span.items() if ov >= GENE_SPAN_SPLIT_FRAC * total_exon_len]
    if len(strong_genes) < 2:
        return []

    # per-exon gene assignment (only among the strong genes to avoid noise from weak overlaps)
    strong_set = set(strong_genes)

    def gene_of_exon_restricted(exon):
        es, ee = exon
        best_name = None
        best_ov = 0
        for gc, gs, ge, _gstrand, name, _bio in genes:
            if name not in strong_set or gc != c or ee <= gs or ge <= es:
                continue
            ov = min(ee, ge) - max(es, gs)
            if ov > best_ov:
                best_ov = ov
                best_name = name
        if best_name is not None and best_ov >= EXON_GENE_ASSIGN_FRAC * (ee - es):
            return best_name
        return None

    exon_gene = [gene_of_exon_restricted(ex) for ex in exons]

    # group consecutive exons with the same assigned gene
    segments = []
    i = 0
    while i < len(exons):
        if exon_gene[i] is None:
            i += 1
            continue
        g = exon_gene[i]
        j = i
        while j < len(exons) and exon_gene[j] == g:
            j += 1
        segments.append((g, i, j - 1))
        i = j

    if len(segments) < 2:
        return []  # not a readthrough by gene annotation

    out = []
    for seg_idx, (g, ei, ej) in enumerate(segments):
        seg_exons = exons[ei:ej + 1]
        seg_start = seg_exons[0][0]
        seg_end = seg_exons[-1][1]
        seg_introns = introns[ei:ej]
        seg_ne = len(seg_exons)
        # extract sub-sequence in transcription orientation
        exon_lens = [ee - es for es, ee in exons]
        seg_exon_lens = [ee - es for es, ee in seg_exons]
        if strand == "+":
            pre = sum(exon_lens[:ei])
            sub_seq = seq[pre:pre + sum(seg_exon_lens)]
        else:
            # minus: transcription order is reverse genomic order
            pre = sum(exon_lens[ej + 1:])
            sub_seq = seq[pre:pre + sum(seg_exon_lens)]
        new_tid = f"{tid}_seg{seg_idx}_{g}"
        new_meta = (c, seg_start, seg_end - 1, strand, seg_ne, nr)
        out.append((new_tid, new_meta, sub_seq, seg_introns))
    return out


def graph_split_transcripts(seqs, meta, skel_introns, genes):
    """Graph-aware readthrough split: decompose transcripts whose exon path visits >=2 distinct genes.

    Returns new seqs, meta, skel_introns dicts with readthrough transcripts replaced by their gene-pure
    sub-transcripts.  Transcripts that do not traverse multiple genes are kept unchanged.
    """
    new_seqs = {}
    new_meta = {}
    new_skel = {}
    n_split = 0
    n_sub = 0
    for tid in seqs:
        sub = _split_transcript_by_genes(tid, meta, skel_introns, seqs, genes)
        if not sub:
            new_seqs[tid] = seqs[tid]
            new_meta[tid] = meta[tid]
            continue
        n_split += 1
        n_sub += len(sub)
        for new_tid, new_meta_tuple, sub_seq, seg_introns in sub:
            new_seqs[new_tid] = sub_seq
            new_meta[new_tid] = new_meta_tuple
            c, s, e, _strand, _ne, _nr = new_meta_tuple
            new_skel[(c, s, e)] = seg_introns
    # copy over skeletons for unsplit transcripts
    for key, intr in skel_introns.items():
        if key not in new_skel:
            new_skel[key] = intr
    return new_seqs, new_meta, new_skel, n_split, n_sub


def annotation_split_components(comp, meta, skel_introns, genes):
    """Split collapsed components whose members are assigned to distinct RefSeq genes.

    Each returned component is a list of tids.  Transcripts with no confident gene assignment are kept
    together in an 'unassigned' sub-component per original component.

    Returns (components, rescued) where rescued is the set of tids that landed in a component created by
    the split (these are candidates for the downstream rescue pass).
    """
    out = []
    rescued = set()
    for members in comp.values():
        by_gene = defaultdict(list)
        for tid in members:
            g = dominant_gene(tid, meta, skel_introns, genes)
            by_gene[g].append(tid)
        # Only split if there are at least two distinct assigned genes.  A single assigned gene + unassigned
        # transcripts is kept as one component (do not create spurious splits for one stray unassigned tx).
        assigned_genes = [g for g in by_gene if g is not None]
        if len(assigned_genes) >= 2:
            for g, group in by_gene.items():
                if group:
                    out.append(group)
                    rescued.update(group)
        else:
            out.append(members)
    return out, rescued


def main(args=None):
    if args is None:
        args = parse_args()
    import time
    t0 = time.perf_counter()

    def tick(msg):
        print(f"[{time.perf_counter() - t0:6.1f}s] {msg}", flush=True)

    seqs = P.load_fasta(args.fa)
    meta = {}
    for line in open(args.meta):
        if line.startswith("id\t"):
            continue
        tid, c, s, e, strand, ne, nr = line.rstrip("\n").split("\t")
        meta[tid] = (c, int(s), int(e), strand, int(ne), int(nr))
    n_assembly = len(seqs)

    # introns per skeleton, keyed by (chrom, start, end) so distinct chains at the same start don't collide
    skel_introns = {}
    for line in open(args.skel):
        if line.startswith("chrom\t"):
            continue
        c, s, e, ne, nr, intr = line.rstrip("\n").split("\t")
        ivs = [tuple(map(int, x.split("-"))) for x in intr.split(";")] if intr else []
        skel_introns[(c, int(s), int(e))] = ivs

    # Optional graph-aware readthrough split: cut transcripts whose exon path visits >=2 distinct genes.
    if args.graph_split:
        genes = load_genes(args.genes)
        seqs, meta, skel_introns, n_split, n_sub = graph_split_transcripts(seqs, meta, skel_introns, genes)
        n_assembly = len(seqs)
        tick(f"graph-split: {n_split:,} readthrough transcripts -> {n_sub:,} gene-pure sub-transcripts")

    # ---- (1) collapse isoforms -> GENE loci by SHARED INTRON JUNCTIONS (union-find) ----
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        r = x
        while parent[r] != r:
            r = parent[r]
        while parent[x] != r:
            parent[x], x = r, parent[x]
        return r

    def union(a, b):
        parent[find(a)] = find(b)

    junc_owner = {}  # (chrom, donor, acceptor) -> a tid that uses it
    for tid, (c, s, e, _strand, _ne, nr) in meta.items():
        find(tid)  # register every tid (single-exon genes stay singletons)
        for d, a in skel_introns.get((c, s, e), []):
            key = (c, d, a)
            if key in junc_owner:
                union(tid, junc_owner[key])
            else:
                junc_owner[key] = tid

    comp = defaultdict(list)
    for tid in meta:
        comp[find(tid)].append(tid)

    # Optional annotation-aware split of readthrough components that span multiple RefSeq genes.
    rescued_tids = set()
    if args.annotation_split:
        genes = load_genes(args.genes)
        comp_groups, rescued_tids = annotation_split_components(comp, meta, skel_introns, genes)
        n_before = len(comp)
        comp = {i: group for i, group in enumerate(comp_groups)}
        tick(f"annotation-split: {n_before:,} components -> {len(comp):,} after splitting by RefSeq gene-of; "
             f"{len(rescued_tids):,} rescued tids")
    else:
        comp = dict(comp)

    reps = []
    for members in comp.values():
        # rep = most reads, tie-break longest span
        reps.append(max(members, key=lambda t: (meta[t][5], meta[t][2] - meta[t][1])))
    reps = [r for r in reps if r in seqs]
    n_loci = len(reps)
    tick(f"collapsed {n_assembly:,} transcripts -> {n_loci:,} gene loci (intron-junction union-find)")

    # ---- (2) exact k-mer counting pre-filter: family-informative k-mers, then candidate reps ----
    # per-rep UNIQUE k-mers (dedup within rep) -> np.unique(concatenation, return_counts) gives EXACT
    # distinct-rep ownership of every observed k-mer (no Bloom-array collisions; KMER=16 space is sparse).
    rep_kmers = {}   # r -> (sorted unique k-mer values, their first-occurrence positions along the transcript)
    chunks = []
    for r in reps:
        vals, pos = np.unique(kmer_hashes(seqs[r]), return_index=True)
        rep_kmers[r] = (vals, pos)
        if vals.size:
            chunks.append(vals)
    allk = np.concatenate(chunks) if chunks else np.empty(0, dtype=np.int64)
    uniq, cnt = np.unique(allk, return_counts=True)        # cnt = #distinct gene-reps owning each k-mer
    info_kmers = uniq[(cnt >= CNT_MIN) & (cnt <= CNT_MAX)]  # family-informative k-mers (sorted)
    del allk, chunks, uniq, cnt
    tick(f"exact k-mer counts built; {len(info_kmers):,} family-informative k-mers "
         f"(owned by {CNT_MIN}..{CNT_MAX} reps)")

    # per-rep informative signature = its k-mers that are family-informative (with positions); a rep can
    # only be in a >= K_SHARE pair if it has >= K_SHARE informative k-mers -> cheap, threshold-free gate.
    def members(vals):                          # boolean membership of vals in the SORTED info_kmers
        idx = np.searchsorted(info_kmers, vals)
        ok = idx < info_kmers.size
        res = np.zeros(vals.size, dtype=bool)
        res[ok] = info_kmers[idx[ok]] == vals[ok]
        return res
    informative_sig = {}   # r -> (sig values sorted, sig positions)
    rescued_candidates = set()
    for r in reps:
        vals, pos = rep_kmers[r]
        if vals.size == 0:
            continue
        m = members(vals)
        sig_size = int(m.sum())
        if sig_size >= K_SHARE or (r in rescued_tids and sig_size >= K_SHARE_RESCUE):
            informative_sig[r] = (vals[m], pos[m])
            if r in rescued_tids and sig_size < K_SHARE:
                rescued_candidates.add(r)
    candidates = list(informative_sig)
    del rep_kmers
    seqlen = {r: len(seqs[r]) for r in candidates}
    tick(f"candidate family-member reps: {len(candidates):,} / {n_loci:,} loci "
         f"({len(rescued_candidates):,} annotation-rescued); "
         f"single-copy genes have < {K_SHARE} informative k-mers -> rejected, never reach POA")

    # inverted index ONLY over family-informative k-mers of CANDIDATE reps (bounded -> no OOM)
    inv = defaultdict(list)
    for r in candidates:
        for h in informative_sig[r][0]:
            inv[int(h)].append(r)
    # dedup co-occurring pairs into a SET (a paralog pair recurs across many buckets but is added ONCE ->
    # bounded; counting per-k-mer is what blew up before)
    pair_set = set()
    truncated = False
    for h, lst in inv.items():
        if len(lst) < 2 or len(lst) > PAIR_CAP:    # skip common-domain buckets (cost guard)
            continue
        for i in range(len(lst)):
            ri = lst[i]
            for j in range(i + 1, len(lst)):
                rj = lst[j]
                pair_set.add((ri, rj) if ri < rj else (rj, ri))
        if len(pair_set) > MAX_PAIRS:
            truncated = True
            break
    tick(f"distinct co-occurring candidate pairs: {len(pair_set):,}{' (TRUNCATED at guard)' if truncated else ''}")

    # CONTIGUOUS-SPAN pre-filter: a real paralog pair shares >= K_SHARE k-mers that SPAN a contiguous
    # block; domain-sharers share k-mers confined to one short domain. The span of shared-k-mer positions
    # in the shorter transcript / its length is a cheap proxy for POA's core_recip -- use the SAME T_CORE
    # threshold POA decides on, so this only pre-rejects pairs POA would reject (recall-safe), at ~0 cost.
    # This is what lets POA run on the real-family subset instead of every domain-sharer (237k -> far less).
    cands = []
    n_kshare = 0
    for a, b in pair_set:
        va, pa = informative_sig[a]
        vb, pb = informative_sig[b]
        common, ia, ib = np.intersect1d(va, vb, assume_unique=True, return_indices=True)
        rescue_pair = (a in rescued_tids or b in rescued_tids)
        min_share = K_SHARE_RESCUE if rescue_pair else K_SHARE
        if common.size < min_share:
            continue
        n_kshare += 1
        sa, sb = pa[ia], pb[ib]
        core = min(int(sa.max() - sa.min()), int(sb.max() - sb.min()))   # contiguous-block proxy (bp)
        if core >= T_CORE * min(seqlen[a], seqlen[b]):
            cands.append((a, b))
    print(f"general assembly transcripts: {n_assembly:,}; gene loci (intron-junction collapse): {n_loci:,}; "
          f"candidate reps: {len(candidates):,} ({len(rescued_candidates):,} annotation-rescued); "
          f"pairs >= {K_SHARE} shared k-mers: {n_kshare:,}; "
          f"after contiguous-span filter: {len(cands):,} POA pairs")
    tick("starting POA grading")

    # ---- (3) POA contiguous-core grading = THE family criterion ----
    with mp.Pool(5, initializer=_init) as pool:
        conf = [r for r in pool.map(_poa, cands, chunksize=64) if r]
    tick(f"POA confirmed {len(conf):,} / {len(cands):,} pairs")

    # persist the confirmed-edge graph (a, b, core_recip) so the family decomposition can be iterated
    # WITHOUT re-running the expensive POA (community/density split reads this file).
    with open(args.edges, "w") as eh:
        eh.write("a\tb\tcore_recip\n")
        for a, b, cr in conf:
            eh.write(f"{a}\t{b}\t{cr}\n")
    tick(f"wrote {len(conf):,} confirmed edges -> {args.edges}")

    # families = connected components of POA-confirmed pairs
    fparent = {}

    def ffind(x):
        fparent.setdefault(x, x)
        r = x
        while fparent[r] != r:
            r = fparent[r]
        while fparent[x] != r:
            fparent[x], x = r, fparent[x]
        return r
    for a, b, cr in conf:
        fparent[ffind(a)] = ffind(b)
    fcomp = defaultdict(set)
    for a, b, cr in conf:
        fcomp[ffind(a)].add(a)
        fcomp[ffind(a)].add(b)
    fams = [c for c in fcomp.values() if len(c) >= 2]
    in_fam = sum(len(c) for c in fams)

    with open(args.out, "w") as fh:
        fh.write("family_id\tn_copies\tmembers\n")
        for i, c in enumerate(sorted(fams, key=lambda c: -len(c))):
            fh.write(f"DNFAM{i}\t{len(c)}\t{','.join(sorted(c))}\n")

    # ---- validation: do de-novo families recover KNOWN ones? (map known loci -> reps by overlap) ----
    known = {"RABL2": [("NC_073235.2", 15131675, 15147533), ("NC_086018.1", 48818440, 48831690)],
             "APOBEC3": [("NC_086018.1", 37023493, 37058621)],
             "RFPL": [("NC_086018.1", 27445867, 30376055)]}
    repmeta = [(meta[r][0], meta[r][1], meta[r][2], r) for r in reps]
    g2fam = {}
    for i, c in enumerate(fams):
        for m in c:
            g2fam[m] = i

    def reps_in(chrom, s, e):
        return [r for (c, rs, re, r) in repmeta if c == chrom and rs < e and re > s]
    print("\nde-novo families:", len(fams), "; transcripts in a multi-copy family:", in_fam)
    print("validation (known families recovered as de-novo families):")
    for name, regions in known.items():
        comps = defaultdict(int)
        for c, s, e in regions:
            for r in reps_in(c, s, e):
                if r in g2fam:
                    comps[g2fam[r]] += 1
        best = max(comps.values()) if comps else 0
        print(f"  {name}: {'RECOVERED' if best >= 2 else ('partial/' + str(best))} (reps co-grouped={best})")
    print(f"\n[wrote {args.out}]")
    print(f"FRAMING: general read-coherence assembler = {n_assembly:,} transcripts; the family/copy layer "
          f"applies to {in_fam:,} transcripts in {len(fams):,} multi-copy families; the rest are single-copy "
          f"(assembled normally).")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Genome-wide de-novo multi-copy family catalog from read-coherence transcripts.")
    parser.add_argument("--fa", default=FA, help="de-novo transcript FASTA")
    parser.add_argument("--meta", default=META, help="de-novo transcript metadata TSV")
    parser.add_argument("--skel", default=SKEL, help="skeleton intron-chain TSV")
    parser.add_argument("--out", default=OUT, help="output families TSV")
    parser.add_argument("--edges", default=EDGES, help="output confirmed edges TSV")
    parser.add_argument("--genes", default=GENES,
                        help="RefSeq gene coordinate TSV (chrom start end strand name biotype)")
    parser.add_argument("--annotation-split", action="store_true",
                        help="split collapsed locus components whose members map to distinct RefSeq genes")
    parser.add_argument("--graph-split", action="store_true",
                        help="split individual readthrough transcripts whose exon path visits >=2 distinct "
                             "RefSeq genes into gene-pure sub-transcripts before locus collapse")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    main(args)

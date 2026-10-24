#!/usr/bin/env python3
"""O3 port (2026-10-23): Strainy's phasing with IsoCon's significance test, and an allele-slot certificate.
Design: docs/superpowers/specs/2026-10-23-o3-strainy-isocon-design.md; decision rules: docs/archive/2026-10/merged/PREREG_o3_port_2026-10-23.md.

  o3port.py run --bam B --copies COPIES.tsv --genome G.fa --net net_reads.tsv --out PREFIX [--families F1,F2]

Per family (the reads of its net with used = 1, as `candidate_copies` built it):
  1. piles: a read joins the pile of every copy holding one of its alignments within AS_TIE of its best (O2's tie rule);
  2. informative columns (Strainy's T-informative SNPs, IsoCon's test): >= 2 alleles each carried by >= MIN_READS reads at
     base quality >= MIN_Q, a minor allele improbable as sequencing error (Poisson tail of the carriers' error probabilities,
     Bonferroni over the tested columns, ALPHA); a base that is no supported allele counts as missing;
  3. read graph (clusteringDivergence 0): two reads are linked when they share >= 1 informative column and disagree at none;
  4. clusters: deterministic label propagation (each distinct allele pattern votes with its read count);
  5. refinement: each cluster is re-clustered on the columns informative inside it, until nothing splits;
  6. certificate: per pile, the largest set of pairwise-conflicting clusters (distinct haplotypes present together); over the
     family, a maximum b-matching of those haplotypes into allele slots (2 per copy) of the copies whose pile holds >= half of
     their reads; an unmatched haplotype is a certified extra copy;
  7. its contig: the consensus of its reads (majority base, deletion or insertion at every position >= MIN_READS of them cover)
     along their alignments on the copy whose pile holds most of them.
Products (the `candidate_copies` contract the acceptance harness reads): PREFIX.candidates.tsv (one row per certified haplotype,
flagged = 1), PREFIX.contigs.fa, PREFIX.clusters.tsv, PREFIX.reads.tsv, PREFIX.families.tsv. Beside them, putative copies (not
certified; 2026-10-23): PREFIX.putative.tsv / .fa (slot holders whose reads tie to other copies while the slot copy's other
allele's do not, or that the slot copy cannot have spliced, `putative_copies`) and PREFIX.loci.tsv (where >= MIN_SUPPORT of the family's reads land outside its copies,
`landing_loci`: an unannotated similar locus, expressed or silent).
"""
import argparse
import collections
import csv
import json
import math
import sys

MIN_READS = 2      # the project's ">= 2 independent observations"
MIN_Q = 10         # base quality below which a base counts as missing
AS_TIE = 10        # O2's tie rule (Eichler): alignments within 10 AS of the read's best
ALPHA = 0.01       # IsoCon-style test, Bonferroni over the tested columns
MAX_ROUNDS = 100   # label propagation stops earlier when no label changes
CLIPPING = True    # milestone 2, third signal: reads the slot copy's DNA cannot hold end to end (--no-clipping turns it off)
CLIP_BP = 50       # a soft clip this long that is not a poly-A / poly-T tail
CLIP_BIN = 20      # clip boundaries in the same CLIP_BIN-bp bin on the same side count as one boundary
ANCHOR = False     # --anchor: a one-site lead suffices at a site the copies' anchor reads show fixed (development, 2026-10-24)
ONE_SITE = False   # --one-site: a one-site lead suffices when the read's haplotype agrees (development, 2026-10-24)
GUARD = False      # --guard: O2 abstains on copies where O3 flags a putative copy (development, 2026-10-24)
PREFER = False     # --prefer: a haplotype cannot take a slot of a copy its reads clearly disfavour (development, 2026-10-24)
RECONSTRUCT = True # O2 uses O3's certified copies (--no-reconstruct: O2 as if O3 had found none; the O3 -> O2 comparison)
SEMI = False       # O2 --semi (opt-in, register 1388): semi-tied reads go to their best copy, exact ties are HOPELESS
COVERAGE = False   # milestone 2 retry: foreign coverage as a putative signal (--coverage)
CORRECT = False    # milestone 2b: IsoCon-style read correction before clustering (--correct)
MIN_SUPPORT = 6    # reads a haplotype needs to take a slot or be certified: the shipped stage's own flag floor
                   # (candidate_copies --min-support 6, ruling R1), so the port is judged on the stage's terms

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


# ---------------------------------------------------------------- statistics and columns (pure; tested)

def poisson_tail(k, lam):
    """P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    term = math.exp(-lam)
    cdf = term
    for j in range(1, k):
        term *= lam / j
        cdf += term
    return max(0.0, 1.0 - cdf) if cdf < 1.0 - 1e-12 else _small_tail(k, lam)


def _small_tail(k, lam):
    """The tail by its leading terms when 1 - cdf underflows (lam small); log space, so a large k cannot overflow."""
    lt = -lam + k * math.log(lam) - math.lgamma(k + 1)
    s, j = 0.0, k
    while lt > -700 and j < k + 50:
        s += math.exp(lt)
        j += 1
        lt += math.log(lam) - math.log(j)
    return s


def err(x):
    """A read's per-base error probability at a column: its empirical rate on the copy when known (IsoCon estimates error
    rates from the reads' own alignments; PacBio base qualities are capped at Q93 and overconfident), else its base quality."""
    return x[3] if len(x) > 3 else 10 ** (-x[2] / 10)


def informative(obs, n_tests=None):
    """obs: {column: [(read, base, qual[, error rate]), ...]}. Returns {column: set of supported alleles} for the columns where
    >= 2 alleles have >= MIN_READS carriers at quality >= MIN_Q and a minor allele passes the error test (Bonferroni over the
    columns tested here, i.e. those with >= 2 such alleles, unless n_tests is given)."""
    cand = {}
    for col, xs in obs.items():
        good = [x for x in xs if x[2] >= MIN_Q and x[1] in "ACGT"]
        cnt = collections.Counter(x[1] for x in good)
        alleles = [b for b, n in cnt.items() if n >= MIN_READS]
        if len(alleles) >= 2:
            cand[col] = (good, cnt, alleles)
    m = n_tests if n_tests is not None else max(1, len(cand))
    out = {}
    for col, (good, cnt, alleles) in cand.items():
        lam = sum(err(x) / 3 for x in good)
        major = max(alleles, key=lambda b: (cnt[b], b))
        real = {major} | {b for b in alleles if b != major and poisson_tail(cnt[b], lam) < ALPHA / m}
        if len(real) >= 2:
            out[col] = real
    return out


def patterns(reads, obs, cols):
    """{read: {col: allele}} restricted to informative columns, a base that is no supported allele counted as missing."""
    pat = {r: {} for r in reads}
    for col, real in cols.items():
        for x in obs.get(col, ()):
            r, b, q = x[0], x[1], x[2]
            if r in pat and q >= MIN_Q and b in real:
                pat[r][col] = b
    return pat


def compatible(a, b):
    shared = a.keys() & b.keys()
    return bool(shared) and all(a[c] == b[c] for c in shared)


def label_propagation(pat):
    """Deterministic label propagation on the zero-disagreement read graph. Reads with the same allele pattern have the same
    neighbourhood, so the graph is built on the distinct patterns, each voting with its read count (a read also hears the other
    reads of its own pattern). Returns a list of read sets (reads with no informative column are left out)."""
    groups = collections.defaultdict(list)
    for r, p in pat.items():
        if p:
            groups[tuple(sorted(p.items()))].append(r)
    keys = sorted(groups, key=lambda k: (-len(groups[k]), k))
    idx = {k: i for i, k in enumerate(keys)}
    nbr = {i: [] for i in range(len(keys))}
    dicts = [dict(k) for k in keys]
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if compatible(dicts[i], dicts[j]):
                nbr[i].append(j)
                nbr[j].append(i)
    label = list(range(len(keys)))
    w = [len(groups[k]) for k in keys]
    for _ in range(MAX_ROUNDS):
        changed = False
        for i in range(len(keys)):
            votes = collections.Counter({label[i]: w[i] - 1})
            for j in nbr[i]:
                votes[label[j]] += w[j]
            top = max(votes.values())
            if votes[label[i]] == top:
                continue
            new = min(l for l, v in votes.items() if v == top)
            label[i] = new
            changed = True
        if not changed:
            break
    clusters = collections.defaultdict(set)
    for k, i in idx.items():
        clusters[label[i]].update(groups[k])
    return sorted(clusters.values(), key=lambda s: (-len(s), sorted(s)))


def correlation_clustering(pat):
    """Correlation clustering of reads by their allele patterns (2026-10-23; the advisor's suggestion; Strainy itself uses label
    propagation). Reads with the same pattern are one node weighted by their count. Two nodes that share >= 1 column are similar
    (+) when they differ at < MIN_READS of them and dissimilar (-) when they differ at >= MIN_READS (the project's two-observation
    rule: one difference can be one error); nodes that share no column are neutral. The partition minimizes the read pairs that
    disagree inside a cluster plus the read pairs that agree across clusters. Deterministic: a pivot start (KwikCluster; the
    heaviest unclustered node takes its similar unclustered neighbours), then each node in weight order moves to the cluster of
    highest net agreement (or alone, when every cluster is net disagreement) until no move lowers the cost. Returns read sets,
    as label_propagation (reads with no informative column are left out)."""
    groups = collections.defaultdict(list)
    for r, p in pat.items():
        if p:
            groups[tuple(sorted(p.items()))].append(r)
    keys = sorted(groups, key=lambda k: (-len(groups[k]), k))
    w = [len(groups[k]) for k in keys]
    dicts = [dict(k) for k in keys]
    n = len(keys)
    nbr = [{} for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            shared = dicts[i].keys() & dicts[j].keys()
            if shared:
                sg = w[i] * w[j] * (1 if sum(dicts[i][c] != dicts[j][c] for c in shared) < MIN_READS else -1)
                nbr[i][j] = nbr[j][i] = sg
    label = [-1] * n
    for i in range(n):
        if label[i] < 0:
            label[i] = i
            for j, sg in nbr[i].items():
                if sg > 0 and label[j] < 0:
                    label[j] = i
    size = collections.Counter(label)
    for _ in range(MAX_ROUNDS):
        moved = False
        for i in range(n):
            score = collections.Counter()
            for j, sg in nbr[i].items():
                score[label[j]] += sg
            cur = score[label[i]]                            # net agreement with the rest of its cluster (0 when alone)
            best, best_s = label[i], cur
            for lab, sc in sorted(score.items()):
                if lab != label[i] and sc > best_s:
                    best, best_s = lab, sc
            if best_s < 0:                                   # every cluster, its own included, is net disagreement:
                best, best_s = n + i, 0                      # alone, under a fresh label of its own
            if best != label[i] and best_s > cur:
                size[label[i]] -= 1
                size[best] += 1
                label[i] = best
                moved = True
        if not moved:
            break
    clusters = collections.defaultdict(set)
    for k, i in zip(keys, range(n)):
        clusters[label[i]].update(groups[k])
    return sorted(clusters.values(), key=lambda s: (-len(s), sorted(s)))


CLUSTERER = {"lp": label_propagation, "cc": correlation_clustering}


def correct_patterns(pat):
    """IsoCon-style correction (milestone 2b, 2026-10-23): a read's allele at a column is replaced when >= MIN_READS of its
    neighbours (reads that share >= 1 other column with it and agree at every other shared column) carry one other allele and
    none carries its own. One pass on the observed patterns (computed per distinct pattern, weighted by its reads); corrected
    patterns feed the clustering only."""
    groups = collections.defaultdict(list)
    for r, p in pat.items():
        groups[tuple(sorted(p.items()))].append(r)
    keys = list(groups)
    dicts = [dict(k) for k in keys]
    w = [len(groups[k]) for k in keys]
    fixed = {}
    for i, p in enumerate(dicts):
        q = dict(p)
        for c, b in p.items():
            votes = collections.Counter()
            for j, o in enumerate(dicts):
                if c not in o:
                    continue
                other = (p.keys() & o.keys()) - {c}
                if other and all(p[x] == o[x] for x in other):
                    votes[o[c]] += w[j] - (1 if j == i else 0)
            votes = +votes
            if votes and votes[b] == 0 and len(votes) == 1:
                alt, n = votes.most_common(1)[0]
                if n >= MIN_READS:
                    q[c] = alt
        fixed[i] = q
    return {r: fixed[i] for i, k in enumerate(keys) for r in groups[k]}


def cluster_pile(reads, obs, method="cc"):
    """Strainy's clustering with refinement: label propagation on the pile's informative columns, then each cluster again on
    the columns informative inside it, until nothing splits. Returns (clusters of >= MIN_READS reads, every column used, the
    pile-level columns). Only the pile-level columns may prove two haplotypes distinct (the certificate): a column that passes
    the test only inside one cluster (smaller Bonferroni family, lower error mass) can come from a recurrent error or chimera."""
    cols = informative({c: [x for x in xs if x[0] in reads] for c, xs in obs.items()})
    used = dict(cols)
    pat = patterns(reads, obs, cols)
    if CORRECT:
        pat = correct_patterns(pat)
    todo = [c for c in CLUSTERER[method](pat) if len(c) >= MIN_READS]
    done = []
    while todo:
        c = todo.pop()
        inner = informative({col: [x for x in xs if x[0] in c] for col, xs in obs.items()})
        if not inner:
            done.append(c)
            continue
        parts = [p for p in CLUSTERER[method](patterns(c, obs, inner)) if len(p) >= MIN_READS]
        if len(parts) >= 2:
            for col, real in inner.items():
                used.setdefault(col, set()).update(real)
            todo.extend(parts)
        else:
            done.append(c)
    return sorted(done, key=lambda s: (-len(s), sorted(s))), used, cols


def haplotype(cluster, obs, cols):
    """A cluster's allele at each column: the strict majority of its members' supported alleles, carried by >= MIN_READS."""
    h = {}
    for col, real in cols.items():
        cnt = collections.Counter(x[1] for x in obs.get(col, ()) if x[0] in cluster and x[2] >= MIN_Q and x[1] in real)
        if cnt:
            b, n = cnt.most_common(1)[0]
            if n >= MIN_READS and 2 * n > sum(cnt.values()):
                h[col] = b
    return h


def conflicting_set(haps):
    """The largest set of pairwise-conflicting haplotypes: distinct haplotypes present together. Two haplotypes conflict when
    they differ at >= MIN_READS shared columns (two independent observations of a difference, the same >= 2 rule as for reads;
    a single differing column is one observation, e.g. a recurrent error splitting off a few reads). Exact by branch and
    bound; piles hold few clusters."""
    n = len(haps)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            shared = haps[i].keys() & haps[j].keys()
            adj[i][j] = adj[j][i] = sum(haps[i][c] != haps[j][c] for c in shared) >= MIN_READS
    best = []

    def grow(cur, rest):
        nonlocal best
        if len(cur) > len(best):
            best = list(cur)
        for k, v in enumerate(rest):
            if len(cur) + len(rest) - k <= len(best):
                return
            grow(cur + [v], [u for u in rest[k + 1:] if adj[v][u]])

    grow([], list(range(n)))
    return sorted(best)


def absorb(hs, keep):
    """v5.1 (design step 5: a haplotype split where the reads show no heterozygosity is one haplotype): {cluster outside the
    conflicting set: the kept cluster it joins}. It joins the one kept haplotype it shares >= 1 column with and is not in
    conflict with (differs at < MIN_READS shared columns, the conflict rule); with none or several such haplotypes it stays out."""
    out = {}
    for j in range(len(hs)):
        if j in keep:
            continue
        fit = [k for k in keep if hs[j].keys() & hs[k].keys()
               and sum(hs[j][c] != hs[k][c] for c in hs[j].keys() & hs[k].keys()) < MIN_READS]
        if len(fit) == 1:
            out[j] = fit[0]
    return out


def junction_usage(pile):
    """{read: (alignment start, alignment end, frozenset of introns)} from each read's record on a pile (CIGAR N runs)."""
    out = {}
    for n, rec in pile.items():
        pos, js = rec.reference_start, []
        for op, k in rec.cigartuples:
            if op == 3:
                js.append((pos, pos + k))
            if op in (0, 2, 3, 7, 8):
                pos += k
        out[n] = (rec.reference_start, rec.reference_end, frozenset(js))
    return out


CANONICAL = {("GT", "AG"), ("GC", "AG"), ("AT", "AC"), ("CT", "AC"), ("CT", "GC"), ("GT", "AT")}   # both strands


def foreign_splicing(reads_i, usage, site):
    """Structure as evidence of a putative copy (2026-10-23): >= MIN_READS of haplotype i's own junctions on the slot copy are non-canonical
    on that copy's DNA (donor / acceptor dinucleotides not GT-AG, GC-AG or AT-AC on either strand): two independent observations,
    the project's rule (a first version asked >= half of them; P5's GWFAM101 had 7 of 17, so this threshold was set after looking
    at a held-out miss and P5 is no evidence for it). i's own junctions = those used
    by >= MIN_READS of its reads and by >= half of its reads that span them (its majority structure, not a minor isoform). Every
    isoform of a copy splices at that copy's own sites, so this cannot be raised by alternative splicing, only by reads the copy
    cannot have produced; it flags, it never certifies. site(j) -> (donor, acceptor) dinucleotides of junction j on the copy."""
    spans = lambda j: sum(1 for r in reads_i if r in usage and usage[r][0] < j[0] and usage[r][1] > j[1])
    ji = collections.Counter(j for r in reads_i if r in usage for j in usage[r][2])
    J = [j for j, k in ji.items() if k >= MIN_READS and 2 * k >= spans(j)]
    return sum(1 for j in J if site(j) not in CANONICAL) >= MIN_READS


def foreign_clipping(reads_i, recs):
    """Clipping as evidence of a putative copy (2026-10-23, milestone 2, third signal): >= half of haplotype i's reads with a record
    on the slot copy are soft-clipped there by >= CLIP_BP bases that are not a poly-A / poly-T tail (>= 70% A or T; a record
    without its sequence, a secondary, counts as not a tail), and >= MIN_SUPPORT of them clip at one boundary (same side, same
    CLIP_BIN-bp bin). Where the deleted copy's UTR or an exon diverges from the slot copy, its reads stop matching at the same
    place; a copy's own alleles align end to end. Like the splicing signal it compares the reads with the slot copy's DNA,
    so a survivor allele with another isoform mix cannot raise it. recs: {read: record on the copy}. Development (after P6's
    verdict, P5 / P6 spent): P6 free-slot recall 6 / 13 with 0 flags on either control arm."""
    tail = lambda q: bool(q) and max(q.count("A"), q.count("T")) >= 0.7 * len(q)
    seen, hit, bnd = 0, 0, collections.Counter()
    for r in reads_i:
        rec = recs.get(r)
        if rec is None or not rec.cigartuples:
            continue
        seen += 1
        q, c, ev = rec.query_sequence or "", rec.cigartuples, []
        if c[0][0] == 4 and c[0][1] >= CLIP_BP and not tail(q[:c[0][1]]):
            ev.append(("L", rec.reference_start // CLIP_BIN))
        if c[-1][0] == 4 and c[-1][1] >= CLIP_BP and not tail(q[-c[-1][1]:]):
            ev.append(("R", rec.reference_end // CLIP_BIN))
        hit += bool(ev)
        bnd.update(ev)
    return seen > 0 and 2 * hit >= seen and bool(bnd) and bnd.most_common(1)[0][1] >= MIN_SUPPORT


def foreign_coverage(reads_i, reads_co, blocks):
    """Coverage as evidence of a putative copy (2026-10-23, milestone 2 retry): >= half of haplotype i's reads each hold an aligned
    block that overlaps no aligned block of the slot copy's other holders' reads (>= MIN_READS of them): i's reads cover sequence
    the copy's own allele is never seen to cover. blocks: {read: [(start, end)]} on the slot copy."""
    co = sorted(b for r in reads_co if r in blocks for b in blocks[r])
    if sum(1 for r in reads_co if r in blocks) < MIN_READS or not co:
        return False
    def free(b):
        return not any(x < b[1] and b[0] < y for x, y in co)
    hit = sum(1 for r in reads_i if r in blocks and any(free(b) for b in blocks[r]))
    return 2 * hit >= sum(1 for r in reads_i if r in blocks) > 0


def putative_copies(demand, match, reads_of, n_piles, foreign=None):
    """The haplotypes RNA cannot prove extra (they hold an allele slot, Theorem A) but whose reads behave like another copy's:
    >= half of them tie (within AS_TIE of their best) to another copy as well as to the slot's, while < half of the reads of
    every other haplotype holding a slot of the same copy do. It compares the sample with itself (tie patterns), never with the
    reference sequence, so a reference from another individual does not move it. A slot copy with no other holder gives no
    contrast: no flag. reads_of: {haplotype: reads}; n_piles: {read: piles holding it}. With `foreign` (a function (haplotype,
    copy, other holders) -> bool, foreign_splicing), a holder whose reads the copy cannot have spliced is flagged too, the signal a
    one-survivor family still offers. Returns [(haplotype, copy, its multi fraction, the other
    holders' multi fractions)]."""
    multi = {i: sum(1 for r in reads_of[i] if n_piles.get(r, 0) >= 2) / max(1, len(reads_of[i])) for i in demand}
    out = []
    for i in demand:
        X = match.get(i)
        co = [j for j in demand if j != i and match.get(j) == X]
        if X is None or not co:
            continue
        if (2 * multi[i] >= 1 and all(2 * multi[j] < 1 for j in co)) or (foreign is not None and foreign(i, X, co)):
            out.append((i, X, multi[i], [multi[j] for j in co]))
    return out


def site_bases(rec, pr, genome):
    """{read index in the primary's orientation: (reference base, read base)} over the record's aligned pairs, both bases in the
    primary's orientation (so two records of one read on different strands compare position by position)."""
    seq, _, lead = oriented(pr, rec)
    flip = rec.is_reverse != pr[2]
    ref = genome.fetch(rec.reference_name, rec.reference_start, rec.reference_end).upper()
    out, L = {}, len(seq)
    for q, t in rec.get_aligned_pairs(matches_only=True):
        i, rb, qb = lead + q, ref[t - rec.reference_start], seq[lead + q]
        if flip:
            i, rb, qb = L - 1 - i, CMP.get(rb, "N"), CMP.get(qb, "N")
        out[i] = (rb, qb, t)
    return out


def column_bases(piles, c, prim, genome, cache):
    """{reference position on copy c: Counter of read bases} over every read in c's pile (each in its own record's orientation, the
    copy's strand), computed once per copy."""
    if c not in cache:
        col = collections.defaultdict(collections.Counter)
        for n, rec in piles[c].items():
            if n not in prim:
                continue
            seq, _, lead = oriented(prim[n], rec)
            for q, t in rec.get_aligned_pairs(matches_only=True):
                col[t][seq[lead + q]] += 1
        cache[c] = col
    return cache[c]


def semi_tied(r, piles, prim, genome, cache):
    """O2 on a tied read nothing forces (user rule, 2026-10-23: "the totally tied ones are hopeless"). Its two best copies by
    (AS, fewest edits): HOPELESS when they give the same AS and edit distance (an exact tie). Otherwise count the distinguishing
    sites the read covers on both: read positions where the two copies' reference bases differ; each votes for the copy whose base
    the read carries. The read goes to the leader only when it leads by >= MIN_READS (2) sites: a conversion in the sample flips one
    site, and an edit-distance lead also counts read errors (P11, register 1390), so neither can place a read alone."""
    sc = sorted(((rec.get_tag("AS"), -rec.get_tag("NM"), c) for c, p in piles.items() if r in p for rec in [p[r]]), reverse=True)
    if len(sc) < 2 or sc[0][:2] == sc[1][:2]:
        return "HOPELESS"
    if r not in prim:
        return "ABSTAIN"
    va, vb = site_votes(r, sc[0][2], sc[1][2], piles, prim, genome, cache)
    if va - vb >= MIN_READS:
        return sc[0][2]
    if ANCHOR and va - vb >= 1:
        fa, fb = site_votes(r, sc[0][2], sc[1][2], piles, prim, genome, cache, fixed=True)
        if fa >= 1 and fb == 0:
            return sc[0][2]
    return "ABSTAIN"


def anchor_bases(piles, c, prim, cache):
    """{reference position on copy c: Counter of bases} over c's ANCHORS: reads in c's pile and in no other copy's pile (no rival
    within AS_TIE, so minimap2's placement is the read's own evidence), in the copy's strand."""
    key = ("anchor", c)
    if key not in cache:
        others = set().union(*(set(p) for d, p in piles.items() if d != c)) if len(piles) > 1 else set()
        col = collections.defaultdict(collections.Counter)
        for n, rec in piles[c].items():
            if n in others or n not in prim:
                continue
            seq, _, lead = oriented(prim[n], rec)
            for q, t in rec.get_aligned_pairs(matches_only=True):
                col[t][seq[lead + q]] += 1
        cache[key] = col
    return cache[key]


def fixed_site(ta, ra, rb, tb, piles, ca, cb, prim, cache):
    """--anchor: the site is fixed in the sample when >= MIN_READS anchors of each copy cover it, each copy's anchors carry its own
    reference base, and none carries the other copy's (a converted allele seen in one anchor read unfixes it)."""
    A, B = anchor_bases(piles, ca, prim, cache)[ta], anchor_bases(piles, cb, prim, cache)[tb]
    return A[ra] >= MIN_READS and B[rb] >= MIN_READS and A[rb] == 0 and B[ra] == 0


def site_votes(r, ca, cb, piles, prim, genome, cache, fixed=False):
    """(sites voting copy ca, sites voting copy cb) for read r: the distinguishing sites (the copies' reference bases differ) it
    covers on both, each counted only when the read's base is carried by >= MIN_SUPPORT reads of ca's pile at that column."""
    a, b = (site_bases(piles[c][r], prim[r], genome) for c in (ca, cb))
    col = column_bases(piles, ca, prim, genome, cache)
    va = vb = 0
    for i in a.keys() & b.keys():
        (ra, qa, ta), (rb, _, tb) = a[i], b[i]
        if ra != rb and ra in "ACGT" and rb in "ACGT" and homologous(i, a, b):
            if fixed and not fixed_site(ta, ra, rb, tb, piles, ca, cb, prim, cache):
                continue
            # the read's base must be an allele of the sample: carried by >= MIN_SUPPORT reads of the family at this column (a
            # read error at a distinguishing site is not; P11's GWFAM118 / dev GWFAM248 r058: conversion + error looked like two sites)
            fwd = rec_strand_base(qa, piles[ca][r], prim[r])
            if col[ta][fwd] < MIN_SUPPORT:
                continue
            va += qa == ra
            vb += qa == rb
    return va, vb


SITE_FLANK = 5     # a distinguishing site needs this many gap-free aligned bases on each side in both alignments


def homologous(i, a, b):
    """Read position i compares homologous reference bases of the two copies: both alignments are gap-free over i +- SITE_FLANK
    (an indel beside a site shifts one alignment by a base, and the shifted positions then differ without being a PSV; P13,
    register 1395)."""
    for x in (a, b):
        t = x[i][2]
        for k in range(-SITE_FLANK, SITE_FLANK + 1):
            if i + k not in x or abs(x[i + k][2] - t) != abs(k):
                return False
    return True


def one_site(r, hap_reads, piles, prim, genome, cache, hap_cache, h):
    """--one-site (IsoCon's rule, decided per haplotype, not per read): a semi-tied read with a lead of >= 1 supported site goes to
    the leader when its haplotype (>= MIN_SUPPORT reads) agrees: the median lead over the haplotype's reads tied to the same two
    copies is >= 1 for that copy. One site seen on many reads of one haplotype is an allele, not a read error."""
    sc = sorted(((rec.get_tag("AS"), -rec.get_tag("NM"), c) for c, p in piles.items() if r in p for rec in [p[r]]), reverse=True)
    if len(sc) < 2 or r not in prim:
        return "ABSTAIN"
    ca, cb = sc[0][2], sc[1][2]
    va, vb = site_votes(r, ca, cb, piles, prim, genome, cache)
    if va - vb < 1:
        return "ABSTAIN"
    key = (h, ca, cb)
    if key not in hap_cache:
        d = sorted(x - y for q in hap_reads if q in piles[ca] and q in piles[cb] and q in prim
                   for x, y in [site_votes(q, ca, cb, piles, prim, genome, cache)])
        hap_cache[key] = len(d) >= MIN_SUPPORT and d[len(d) // 2] >= 1
    return ca if hap_cache[key] else "ABSTAIN"


def disfavoured(reads, opts, piles, prim, genome, cache):
    """--prefer: the option copies a haplotype's reads clearly disfavour: z, when for some other option c the median over its reads
    on both of (sites voting c - sites voting z) is >= MIN_READS. Such a haplotype cannot hold a slot of z (a sample copy the
    reference lacks, X*, otherwise fills a near copy's free slot: Theorem A's blind spot)."""
    out = set()
    for c in opts:
        for z in opts:
            if c == z or z in out:
                continue
            d = sorted(x - y for r in reads if r in piles[c] and r in piles[z] and r in prim
                       for x, y in [site_votes(r, c, z, piles, prim, genome, cache)])
            if d and d[len(d) // 2] >= MIN_READS:
                out.add(z)
    return out


def rec_strand_base(b, rec, pr):
    """A base in the primary's orientation, back in this record's orientation (the copy's strand)."""
    return CMP.get(b, "N") if rec.is_reverse != pr[2] else b


def landing_loci(records, spans, min_reads):
    """Where a family's reads land outside its copies: records = [(chrom, start, end, read, near_best)] of the family's reads (every
    primary and secondary record; near_best = within AS_TIE of the read's best AS); spans = [(chrom, start, end)] of the copies.
    Records overlapping no copy are merged into loci (overlapping spans); a locus with >= min_reads distinct reads is reported as
    (chrom, start, end, reads, reads with a near-best record there = an expressed copy's own reads, reads with only weaker records
    there = paralogs' reads reaching it, as they reach a silent pseudogene)."""
    out = [r for r in records if not any(r[0] == c and r[1] < e and s < r[2] for c, s, e in spans)]
    loci = []
    for ch, s, e, n, near in sorted(out):
        if loci and loci[-1][0] == ch and s < loci[-1][2]:
            loci[-1][2] = max(loci[-1][2], e)
        else:
            loci.append([ch, s, e, set(), set()])
        L = loci[-1]
        L[3].add(n)
        if near:
            L[4].add(n)
    return [(ch, s, e, len(r), len(nb), len(r - nb)) for ch, s, e, r, nb in loci if len(r) >= min_reads]


def b_matching(options, capacity, order=None):
    """options: {haplotype: [copies]}; capacity: {copy: slots}. Maximum b-matching by augmenting paths, haplotypes taken in
    `order` (default: fewest options first). An augmenting path never unmatches a matched haplotype, so the ones left
    unmatched are the last in the order. Returns {haplotype: copy or None}."""
    load = collections.defaultdict(list)
    match = {}

    def augment(h, seen):
        for c in options[h]:
            if c in seen:
                continue
            seen.add(c)
            if len(load[c]) < capacity[c]:
                load[c].append(h)
                match[h] = c
                return True
            for other in list(load[c]):
                if augment_from(other, c, seen):
                    load[c].remove(other)
                    load[c].append(h)
                    match[h] = c
                    return True
        return False

    def augment_from(h, avoid, seen):
        for c in options[h]:
            if c == avoid or c in seen:
                continue
            seen.add(c)
            if len(load[c]) < capacity[c]:
                load[c].append(h)
                match[h] = c
                return True
            for other in list(load[c]):
                if augment_from(other, c, seen):
                    load[c].remove(other)
                    load[c].append(h)
                    match[h] = c
                    return True
        return False

    for h in (order if order is not None else sorted(options, key=lambda h: (len(options[h]), h))):
        if not augment(h, set()):
            match[h] = None
    return match


def forced_edges(options, capacity, match, order):
    """O2 (milestone 3, 2026-10-23): {haplotype: copy} for the matched edges present in every maximum b-matching. An edge (i, X) is
    forced when the maximum matching without it is smaller: no reshuffling of the slots can place i elsewhere and keep every
    other haplotype placed."""
    size = sum(1 for v in match.values() if v is not None)
    out = {}
    for i, X in match.items():
        if X is None:
            continue
        opt = dict(options)
        opt[i] = [c for c in options[i] if c != X]
        alt = b_matching(opt, capacity, order)
        if sum(1 for v in alt.values() if v is not None) < size:
            out[i] = X
    return out


def consensus(alignments):
    """alignments: {read: [(ref_pos, base or '-', inserted string)]} on one copy. The majority at every position >= MIN_READS
    of them cover (a base, or '-' for a deletion), followed by the majority insertion when most of the covering reads carry
    one. Positions in genome order."""
    at = collections.defaultdict(list)
    for read, steps in alignments.items():
        for pos, base, ins in steps:
            at[pos].append((base, ins))
    out = []
    for pos in sorted(at):
        xs = at[pos]
        if len(xs) < MIN_READS:
            continue
        b, _ = collections.Counter(b for b, _ in xs).most_common(1)[0]
        if b != "-":
            out.append(b)
        ins = collections.Counter(i for _, i in xs if i)
        if ins:
            s, n = ins.most_common(1)[0]
            if 2 * sum(ins.values()) > len(xs) and n >= MIN_READS:
                out.append(s)
    return "".join(out)


# ---------------------------------------------------------------- BAM side

def load_copies(path):
    by_fam = collections.defaultdict(list)
    for r in csv.DictReader(open(path), delimiter="\t"):
        blocks = sorted(tuple(map(int, b.split("-"))) for b in r["exons"].split(",") if b)
        by_fam[r["family_id"]].append((r["tid"], r["chrom"], blocks[0][0], blocks[-1][1]))
    return by_fam


def load_net(path, families):
    net = collections.defaultdict(set)
    n_net = collections.Counter()
    for r in csv.DictReader(open(path), delimiter="\t"):
        if families and r["family"] not in families:
            continue
        n_net[r["family"]] += 1
        if r["used"] == "1":
            net[r["family"]].add(r["read"])
    return net, n_net


def read_primaries(bam, wanted, records=None):
    """{read: (seq, quals, is_reverse, best AS over every record of the read)} in one pass over the BAM; with `records` (a dict),
    also every non-supplementary record of the wanted reads as (chrom, start, end, AS, is secondary) per read."""
    prim, best = {}, collections.defaultdict(lambda: -10 ** 9)
    for rec in bam.fetch(until_eof=True):
        n = rec.query_name
        if n not in wanted or rec.is_unmapped:
            continue
        if records is not None and not rec.is_supplementary and rec.has_tag("AS"):
            records.setdefault(n, []).append((rec.reference_name, rec.reference_start, rec.reference_end, rec.get_tag("AS"),
                                              rec.is_secondary))
        if rec.has_tag("AS"):
            best[n] = max(best[n], rec.get_tag("AS"))
        if not rec.is_secondary and not rec.is_supplementary and rec.query_sequence:
            q = rec.query_qualities
            prim[n] = (rec.query_sequence, list(q) if q is not None else [40] * len(rec.query_sequence), rec.is_reverse)
    return {n: v + (best[n],) for n, v in prim.items()}


def oriented(prim, rec):
    """The read's bases and qualities in this record's orientation, and the offset of its first CIGAR query base."""
    seq, qual, rev = prim[:3]
    if rec.is_reverse != rev:
        seq, qual = seq.translate(COMP)[::-1], qual[::-1]
    lead = rec.cigartuples[0][1] if rec.cigartuples and rec.cigartuples[0][0] == 5 else 0
    return seq, qual, lead


def pile_records(bam, copies, reads, prim):
    """{copy: {read: record}} for the records within AS_TIE of the read's best that fall on the copy's span."""
    piles = {}
    for cid, chrom, s, e in copies:
        p = {}
        for rec in bam.fetch(chrom, s, e):
            n = rec.query_name
            if n not in reads or n not in prim or rec.is_unmapped or not rec.has_tag("AS"):
                continue
            if rec.get_tag("AS") < prim[n][3] - AS_TIE:
                continue
            if n not in p or rec.get_tag("AS") > p[n].get_tag("AS"):
                p[n] = rec
        piles[cid] = p
    return piles


def walk(rec):
    """(kind, ref position, record query index) along the CIGAR: M/=/X bases, D positions, I runs (attached to the previous
    reference position, with the run's length); N (introns), S and H advance without emitting."""
    qp, rp = 0, rec.reference_start
    for op, n in rec.cigartuples:
        if op in (0, 7, 8):
            for k in range(n):
                yield "M", rp + k, qp + k
            qp += n
            rp += n
        elif op == 1:
            yield "I", rp - 1, (qp, n)
            qp += n
        elif op == 2:
            for k in range(n):
                yield "D", rp + k, None
            rp += n
        elif op == 3:
            rp += n
        elif op == 4:
            qp += n


def observations(pile, prim, span):
    """{column: [(read, base, qual)]} and {read: [(pos, base or '-', insertion)]} for a pile, restricted to the copy's span."""
    obs = collections.defaultdict(list)
    steps = {}
    s0, e0 = span
    for n, rec in pile.items():
        seq, qual, lead = oriented(prim[n], rec)
        aligned = sum(k for op, k in rec.cigartuples if op in (0, 1, 7, 8))
        eps = rec.get_tag("NM") / aligned if rec.has_tag("NM") and aligned else None
        st = []
        for kind, rp, q in walk(rec):
            if not (s0 <= rp < e0):
                continue
            if kind == "M":
                b, bq = seq[q + lead], qual[q + lead]
                st.append((rp, b, ""))
                obs[rp].append((n, b, bq) if eps is None else (n, b, bq, eps))
            elif kind == "D":
                st.append((rp, "-", ""))
            elif st and st[-1][0] == rp:
                qs, ln = q
                pos, b, ins = st[-1]
                st[-1] = (pos, b, ins + seq[qs + lead:qs + lead + ln])
        steps[n] = st
    return obs, steps


def explanation(hp, parents):
    """DADA2's bimera rule on a haplotype's columns. parents: [(name, {col: base})], each more abundant than the haplotype or a
    reference copy. Returns ("single", name) when one parent carries every allele, ("join", (left, right)) when a left part
    matches one parent and the rest another (a chimeric molecule), else None. A column a parent does not cover counts as a
    mismatch for that parent."""
    cols = sorted(hp)
    if not cols:
        return None
    pre, suf = [], []
    for name, pb in parents:
        m = [pb.get(c) == hp[c] for c in cols]
        if all(m):
            return "single", name
        lp = next((k for k, x in enumerate(m) if not x), len(m))
        ls = next((k for k, x in enumerate(reversed(m)) if not x), len(m))
        pre.append((lp, name))
        suf.append((ls, name))
    pre.sort(reverse=True)
    suf.sort(reverse=True)
    for lp, a in pre[:2]:
        for ls, b in suf[:2]:
            if a != b and lp + ls >= len(cols):
                return "join", (a, b)
    return None


def homology(home, others, genome, tmpdir):
    """{other copy: {home position: (other position, strand)}} from minimap2 -c --eqx -x asm20 of each other copy's locus
    against the home copy's locus (best-matching alignment first; a position keeps its first mapping)."""
    import os
    import re
    import subprocess
    hid, hch, hs, he = home
    t = os.path.join(tmpdir, "t.fa")
    q = os.path.join(tmpdir, "q.fa")
    with open(t, "w") as fh:
        fh.write(f">{hid}\n{genome.fetch(hch, hs, he)}\n")
    with open(q, "w") as fh:
        for cid, ch, s, e in others:
            fh.write(f">{cid}\n{genome.fetch(ch, s, e)}\n")
    paf = subprocess.run(["minimap2", "-c", "--eqx", "-x", "asm20", "-t", "1", t, q], capture_output=True, text=True).stdout
    start = {cid: s for cid, _, s, _ in others}
    maps = {}
    lines = [ln.split("\t") for ln in paf.splitlines()]
    for p in sorted(lines, key=lambda p: -int(p[9])):
        cid, strand, qs, qe, ts = p[0], p[4], int(p[2]), int(p[3]), int(p[7])
        cg = next(x[5:] for x in p[12:] if x.startswith("cg:Z:"))
        mp = maps.setdefault(cid, {})
        tpos, qo = ts, 0
        for n, op in re.findall(r"(\d+)([=XIDM])", cg):
            n = int(n)
            if op in "=XM":
                for k in range(n):
                    qq = qs + qo + k if strand == "+" else qe - 1 - (qo + k)
                    mp.setdefault(hs + tpos + k, (start[cid] + qq, strand))
                tpos += n
                qo += n
            elif op == "I":
                qo += n
            elif op == "D":
                tpos += n
    return maps


def pile_consensus(obs):
    return {c: collections.Counter(x[1] for x in xs).most_common(1)[0][0] for c, xs in obs.items() if xs}


CMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def private_detail(i, merged, match, options):
    """(pile, [(column, allele)]): the alleles of unmatched haplotype i that no slot-holding haplotype of its option copies
    carries, on the pile where i has the most of them (the first such pile on a tie)."""
    best = (None, [])
    for cid, hp in merged[i]["hap"].items():
        holders = [merged[j]["hap"][cid] for j in match if match[j] is not None and j != i
                   and cid in merged[j]["hap"] and cid in options[i]]
        priv = [(col, b) for col, b in hp.items() if all(h.get(col) != b for h in holders)]
        if best[0] is None or len(priv) > len(best[1]):
            best = (cid, priv)
    return best


def private_alleles(i, merged, match, options):
    """Alleles of unmatched haplotype i that no slot-holding haplotype of its option copies carries, counted per pile where
    i has a pattern (the largest count). A haplotype whose every allele some slot holder carries is a mosaic of them (a
    chimeric molecule), not sequence the reference copies cannot hold: an extra copy needs >= MIN_READS private alleles."""
    return len(private_detail(i, merged, match, options)[1])


def divergence(h, ref):
    """Fraction of a haplotype's columns where its allele is not the reference base."""
    cols = [c for c in h if ref.get(c, "N") != "N"]
    return sum(h[c] != ref[c] for c in cols) / len(cols) if cols else 0.0


def run_family(fam, copies, reads, bam, prim, genome, slot_order="divergence", tmpdir=".", dump=None, absorb_fragments=False,
               floor_first=True, clusterer="cc", structure=True):
    """dump: a dict filled with the family's internal states (piles, clusters, haplotype fates, private alleles) when given;
    the products do not depend on it."""
    piles = pile_records(bam, copies, reads, prim)
    span = {cid: (s, e) for cid, _, s, e in copies}
    chrom = {cid: ch for cid, ch, _, _ in copies}
    haps, per_copy, allhaps, cons, refs = [], {}, [], {}, {}
    if dump is not None:
        dump.update(family=fam, copies={c: {"pile": sorted(p)} for c, p in piles.items()}, clusters=[], merged=[])
    for cid, pile in piles.items():
        if len(pile) < 2 * MIN_READS:
            continue
        obs, steps = observations(pile, prim, span[cid])
        clusters, cols, top = cluster_pile(set(pile), obs, clusterer)
        hs = [haplotype(c, obs, top) for c in clusters]
        if floor_first:   # v6: only clusters at the support floor can be counted as distinct haplotypes; with tiny clusters in
            big = [i for i, c in enumerate(clusters) if len(c) >= MIN_SUPPORT]   # the count, three 2-3-read artifacts beat a
            keep = [big[k] for k in conflicting_set([hs[i] for i in big])] if big else []   # 93-read allele (largest set wins)
        else:
            keep = conflicting_set(hs) if hs else []
        joined = absorb(hs, keep) if absorb_fragments else {}
        for j, k in joined.items():
            clusters[k] = clusters[k] | clusters[j]
        for k in set(joined.values()):
            hs[k] = haplotype(clusters[k], obs, top)
        per_copy[cid] = (len(pile), len(clusters), len(keep))
        ref = {col: genome.fetch(chrom[cid], col, col + 1).upper() for col in top}
        cons[cid] = pile_consensus(obs)
        for i, c in enumerate(clusters):
            allhaps.append((cid, len(c), hs[i]))
        for i in keep:
            haps.append((cid, clusters[i], steps, divergence(hs[i], ref), hs[i]))
        if dump is not None:
            dump["copies"][cid]["columns"] = len(top)
            dump["clusters"] += [{"copy": cid, "reads": sorted(c), "hap": {str(k): v for k, v in hs[i].items()},
                                  "ref": {str(k): ref[k] for k in hs[i]}, "kept": i in keep, "joined": i in joined}
                                 for i, c in enumerate(clusters)]
    # one haplotype found in several piles (tied reads) is one haplotype: merge by read overlap (>= half of the smaller)
    merged = []
    for cid, reads_h, steps, dv, hp in sorted(haps, key=lambda x: (-len(x[1]), x[0])):
        for m in merged:
            if len(m["reads"] & reads_h) * 2 >= min(len(m["reads"]), len(reads_h)):
                if absorb_fragments:            # v5.1: the merged haplotype holds the reads of every pile it was found in
                    m["reads"] |= reads_h
                m["piles"].add(cid)
                m["div"] = min(m["div"], dv)
                m["hap"].setdefault(cid, hp)
                break
        else:
            merged.append({"reads": set(reads_h), "piles": {cid}, "div": dv, "hap": {cid: hp}})
    options = {}
    for i, m in enumerate(merged):
        options[i] = sorted(c for c, p in piles.items() if len(m["reads"] & p.keys()) * 2 >= len(m["reads"]))
    # DADA2's rule before the slots: a duplicate of a more abundant haplotype shares its slot, a chimera (left part one parent,
    # right part another) needs none; a haplotype equal to copy Z's reference sequence is Z's allele (Z becomes an option)
    maps_of, why, fate = {}, collections.Counter(), {}
    demand = []
    for i, m in enumerate(merged):
        if len(m["reads"]) < MIN_SUPPORT:
            why["below the floor"] += 1
            fate[i] = ("below the floor", None)
            continue
        home = max(m["hap"], key=lambda c: (len(m["reads"] & piles[c].keys()), c))
        hp = m["hap"][home]
        if home not in maps_of:
            maps_of[home] = homology((home, chrom[home]) + span[home], [c for c in copies if c[0] != home], genome, tmpdir)
        mp = maps_of[home]
        parents = [(f"ref:{home}", {c: genome.fetch(chrom[home], c, c + 1).upper() for c in hp})]
        for cid, ch, _, _ in copies:
            if cid == home or cid not in mp:
                continue
            pb = {}
            for c in hp:
                if c in mp[cid]:
                    z, st = mp[cid][c]
                    b = genome.fetch(ch, z, z + 1).upper()
                    pb[c] = b if st == "+" else CMP.get(b, "N")
            parents.append((f"ref:{cid}", pb))
        for k, (pid, n, g) in enumerate(allhaps):
            if n <= len(m["reads"]):
                continue
            pb = {}
            for c in hp:
                if pid == home:
                    pb[c] = g.get(c, cons[pid].get(c))
                elif pid in mp and c in mp[pid]:
                    z, st = mp[pid][c]
                    b = g.get(z, cons[pid].get(z, genome.fetch(chrom[pid], z, z + 1).upper()))
                    pb[c] = b if st == "+" else CMP.get(b, "N")
            parents.append((f"hap{k}", pb))
        ex = explanation(hp, parents)
        if ex is None:
            demand.append(i)
            why["unexplained"] += 1
            fate[i] = ("unexplained", None)
        elif ex[0] == "single" and ex[1].startswith("ref:"):
            z = ex[1][4:]
            if z not in options[i] and z in piles:
                options[i] = sorted(options[i] + [z])
            demand.append(i)
            why["equals a reference copy"] += 1
            fate[i] = ("equals a reference copy", z)
        else:
            why["duplicate" if ex[0] == "single" else "chimera"] += 1
            fate[i] = ("duplicate" if ex[0] == "single" else "chimera", ex[1])
    if PREFER:
        pc = {}
        for i in demand:
            drop = disfavoured(merged[i]["reads"], options[i], piles, prim, genome, pc)
            if drop and len(drop) < len(options[i]):
                options[i] = [c for c in options[i] if c not in drop]
    capacity = {c: 2 for c in piles}
    if slot_order == "support":
        order = sorted(demand, key=lambda i: (-len(merged[i]["reads"]), merged[i]["div"], i))
    else:
        order = sorted(demand, key=lambda i: (merged[i]["div"], -len(merged[i]["reads"]), i))
    match = b_matching({i: options[i] for i in demand}, capacity, order)
    forced = forced_edges({i: options[i] for i in demand}, capacity, match, order)
    certified = []
    n_piles = collections.Counter(r for p in piles.values() for r in p)
    usage = {c: junction_usage(piles[c]) for c in per_copy}
    site = lambda X: (lambda j: (genome.fetch(chrom[X], j[0], j[0] + 2).upper(), genome.fetch(chrom[X], j[1] - 2, j[1]).upper()))
    blk = {c: {n: rec.get_blocks() for n, rec in piles[c].items()} for c in per_copy}
    foreign = (lambda i, X, co: X in usage and (foreign_splicing(merged[i]["reads"], usage[X], site(X)) or
               (CLIPPING and foreign_clipping(merged[i]["reads"], piles[X])) or
               (COVERAGE and foreign_coverage(merged[i]["reads"], set().union(*(merged[j]["reads"] for j in co)), blk[X]))))\
        if structure else None
    looks = []
    for i, X, mi, mc in putative_copies(demand, match, {i: merged[i]["reads"] for i in demand}, n_piles, foreign):
        if len(merged[i]["reads"]) < MIN_SUPPORT:
            continue
        _, steps = observations({r: piles[X][r] for r in sorted(merged[i]["reads"]) if r in piles[X]}, prim, span[X])
        looks.append((X, merged[i]["reads"], consensus(steps), mi, mc))
    for i in demand:
        m = merged[i]
        if match[i] is None and private_alleles(i, merged, match, options) >= MIN_READS:
            home = max(m["piles"], key=lambda c: (len(m["reads"] & piles[c].keys()), c))
            # sorted: the consensus breaks a tie by read order, which must not follow the process's string-hash order
            _, steps = observations({r: piles[home][r] for r in sorted(m["reads"]) if r in piles[home]}, prim, span[home])
            certified.append((home, m["reads"], consensus(steps)))
    if dump is not None:
        for i, m in enumerate(merged):
            row = {"reads": sorted(m["reads"]), "piles": sorted(m["piles"]), "div": m["div"], "options": options[i],
                   "fate": fate[i][0], "by": fate[i][1], "slot": match.get(i) if i in demand else None,
                   "certified": any(rs is m["reads"] for _, rs, _ in certified),
                   "putative": any(rs is m["reads"] for _, rs, _, _, _ in looks)}
            if i in demand and match[i] is None:
                cid, priv = private_detail(i, merged, match, options)
                row["private"] = [(cid, col, b, genome.fetch(chrom[cid], col, col + 1).upper()) for col, b in priv]
            dump["merged"].append(row)
    piled = set().union(*(set(piles[c]) for c in per_copy)) if per_copy else set()
    # O2: each AS-tied read (in >= 2 piles) goes to its haplotype's forced copy, or to the certified copy it belongs to; else abstain
    # (v2, 2026-10-23: only when the haplotype was weighed for every copy the read ties to; a haplotype seen on one pile is
    # trivially forced there, though its tied reads fit the other copy too)
    n_in = collections.Counter(r for p in piles.values() for r in p)
    home_of = {}
    for i in demand:
        for r in merged[i]["reads"]:
            home_of.setdefault(r, i)
    cert_reads = {r: k for k, (_, rs, _) in enumerate(certified) for r in rs} if RECONSTRUCT else {}
    # --guard: where O3 flags a putative copy the copy set is not reconstructed: abstain on every tied read touching that copy or
    # the copies its haplotype could hold a slot on (a free-slot X* makes forced edges and site leads there unsafe, register 1393)
    unsure = {X for X, *_ in looks} | {c for _, rs, *_ in looks for i in demand if merged[i]["reads"] is rs for c in options[i]}
    o2, semi_cache, hap_cache = [], {}, {}
    for r, k in n_in.items():
        if k < 2:
            continue
        if r in cert_reads:
            o2.append((r, f"new:{cert_reads[r]}", k))
        elif GUARD and {c for c, p in piles.items() if r in p} & unsure:
            o2.append((r, "ABSTAIN", k))
        elif r in home_of and home_of[r] in forced and {c for c, p in piles.items() if r in p} <= set(options[home_of[r]]):
            o2.append((r, forced[home_of[r]], k))
        elif SEMI:
            a = semi_tied(r, piles, prim, genome, semi_cache)
            if a == "ABSTAIN" and ONE_SITE and r in home_of:
                a = one_site(r, merged[home_of[r]]["reads"], piles, prim, genome, semi_cache, hap_cache, home_of[r])
            o2.append((r, a, k))
        else:
            o2.append((r, "ABSTAIN", k))
    return certified, per_copy, len(merged), why, piled, looks, o2


def isocon_side(prefix, net_path):
    """The IsoCon side (2026-10-23, design step 'reads that align to no copy'): the reference-free candidates of the shipped stage
    (`candidate_copies`, in-house IsoCon: clusters of the whole net -> consensus -> significance merge -> flag) built mostly from
    reads the Strainy side never sees (`isocon_keep`). Returns {family: [(stage candidate, members, contig, members via B)]}
    for every flagged candidate."""
    via = {(r["family"], r["read"]): r["via"] for r in csv.DictReader(open(net_path), delimiter="\t")}
    cand = {r["candidate"]: r for r in csv.DictReader(open(f"{prefix}.candidates.tsv"), delimiter="\t")}
    owner = {r["cluster"]: r["candidate"] for r in csv.DictReader(open(f"{prefix}.clusters.tsv"), delimiter="\t")
             if r["candidate"] in cand}
    members = collections.defaultdict(set)
    for r in csv.DictReader(open(f"{prefix}.reads.tsv"), delimiter="\t"):
        if r["cluster"] in owner:
            members[owner[r["cluster"]]].add((r["family"], r["read"]))
    seqs, name = {}, None
    for ln in open(f"{prefix}.contigs.fa"):
        ln = ln.strip()
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(ln)
    out = collections.defaultdict(list)
    for c, r in sorted(cand.items()):
        m = members[c]
        if r["flagged"] == "1" and m and c in seqs:
            out[r["family"]].append((c, {x[1] for x in m}, "".join(seqs[c]), sum(via.get(x) == "B" for x in m)))
    return out


def isocon_keep(members, n_via_b, piled, rule):
    """v5 (rule no-pile): at least half of the members are in no pile the Strainy side phased (no record within AS_TIE of their
    best on a copy with >= 2 * MIN_READS pile reads: reads that align better elsewhere, or not at all). v4 (rule via-b): at
    least half joined the net by attribution (no record on a copy at all)."""
    if rule == "via-b":
        return 2 * n_via_b >= len(members)
    return 2 * len(members - piled) >= len(members)


def cmd_run(a):
    import pysam
    fams = set(a.families.split(",")) if a.families else None
    by_fam = load_copies(a.copies)
    net, n_net = load_net(a.net, fams)
    bam = pysam.AlignmentFile(a.bam)
    wanted = set().union(*net.values()) if net else set()
    records = {}
    prim = read_primaries(bam, wanted, records)
    genome = pysam.FastaFile(a.genome)
    iso = isocon_side(a.isocon, a.net) if a.isocon else {}
    import os
    import tempfile
    tmpdir = tempfile.mkdtemp(dir=os.path.dirname(os.path.abspath(a.out)))
    dump_fh = open(a.dump, "w") if a.dump else None
    with open(f"{a.out}.candidates.tsv", "w") as ct, open(f"{a.out}.contigs.fa", "w") as cf, \
            open(f"{a.out}.clusters.tsv", "w") as cl, open(f"{a.out}.reads.tsv", "w") as rd, \
            open(f"{a.out}.families.tsv", "w") as ft:
        ct.write("family\tcandidate\tn_clusters\tn_reads\tflagged\tunion_len\tnearest_locus\td\tn_net\tn_used\n")
        cl.write("family\tcluster\tcandidate\tn_reads\tconsensus_len\tfate\tlinked_to\td\n")
        rd.write("read\tfamily\tcluster\n")
        ft.write("family\tn_net\tn_used\tpiles\tpile_reads\tclusters\tconflicting\thaplotypes\tcertified\n")
        lh = open(f"{a.out}.putative.tsv", "w")
        lf = open(f"{a.out}.putative.fa", "w")
        lo = open(f"{a.out}.loci.tsv", "w")
        lh.write("family\tcandidate\tcopy\tn_reads\tmulti\tother_holders_multi\tconsensus_len\n")
        lo.write("family\tlocus\tn_reads\tnear_best\tweaker_only\n")
        o2f = open(f"{a.out}.o2.tsv", "w")
        o2f.write("family\tread\tassigned\tn_piles\n")
        for fam in sorted(net):
            reads = net[fam]
            dump = {} if dump_fh else None
            certified, per_copy, n_h, why, piled, looks, o2 = run_family(fam, by_fam[fam], reads, bam, prim, genome, a.slot_order,
                                                                     tmpdir, dump, a.absorb, not a.v5_clique, a.clusterer,
                                                                     not a.no_structure)
            for r, c, k in sorted(o2):
                o2f.write(f"{fam}\t{r}\t{c}\t{k}\n")
            for k, (X, rs, seq, mi, mc) in enumerate(looks):
                lh.write(f"{fam}\tputative_{fam}_{k}\t{X}\t{len(rs)}\t{mi:.3f}\t{','.join(f'{x:.3f}' for x in mc)}\t{len(seq)}\n")
                lf.write(f">putative_{fam}_{k}\n{seq}\n")
            recs = []
            for n in reads:
                b = max((x[3] for x in records.get(n, [])), default=None)
                recs += [(ch, s0, e0, n, asc >= b - AS_TIE) for ch, s0, e0, asc, _ in records.get(n, [])]
            for ch, s0, e0, nr, nb, so in landing_loci(recs, [(c[1], c[2], c[3]) for c in by_fam[fam]], MIN_SUPPORT):
                lo.write(f"{fam}\t{ch}:{s0}-{e0}\t{nr}\t{nb}\t{so}\n")
            side = [(c, rs, seq) for c, rs, seq, vb in iso.get(fam, []) if isocon_keep(rs, vb, piled, a.isocon_rule)]
            if dump_fh:
                dump["isocon"] = [{"candidate": c, "reads": sorted(rs)} for c, rs, _ in side]
                dump_fh.write(json.dumps(dump) + "\n")
            certified = certified + [(f"isocon:{c}", rs, seq) for c, rs, seq in side]
            ft.write(f"{fam}\t{n_net[fam]}\t{len(reads)}\t{len(per_copy)}\t{sum(v[0] for v in per_copy.values())}\t"
                     f"{sum(v[1] for v in per_copy.values())}\t{sum(v[2] for v in per_copy.values())}\t{n_h}\t{len(certified)}\n")
            for k, (home, rs, seq) in enumerate(certified):
                cand = f"cand_{fam}_{k}"
                ct.write(f"{fam}\t{cand}\t1\t{len(rs)}\t1\t{len(seq)}\t{home}\tNA\t{n_net[fam]}\t{len(reads)}\n")
                cf.write(f">{cand}\n{seq}\n")
                cl.write(f"{fam}\t{fam}:h{k}\t{cand}\t{len(rs)}\t{len(seq)}\tnew_copy\t-\tNA\n")
                for r in sorted(rs):
                    rd.write(f"{r}\t{fam}\t{fam}:h{k}\n")
            print(f"[o3port] {fam}: reads {len(reads)}, piles {len(per_copy)}, haplotypes {n_h} {dict(why)}, "
                  f"certified {len(certified)}, putative {len(looks)}", file=sys.stderr, flush=True)
        for fh in (lh, lf, lo, o2f):
            fh.close()


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--bam", required=True)
    r.add_argument("--copies", required=True)
    r.add_argument("--genome", required=True)
    r.add_argument("--net", required=True)
    r.add_argument("--out", required=True)
    r.add_argument("--families", default="")
    r.add_argument("--slot-order", choices=["divergence", "support"], default="divergence")
    r.add_argument("--isocon", default="", help="the stage's product prefix (candidate_copies) for the IsoCon side")
    r.add_argument("--clusterer", choices=["cc", "lp"], default="cc", help="read clustering: correlation clustering of agreeing "
                   "/ disagreeing read pairs (default since 2026-10-23, after the P5 held-out), or label propagation on the "
                   "zero-disagreement graph (Strainy's; v6 as pre-registered)")
    r.add_argument("--coverage", action="store_true", help="putative copies also from coverage the slot copy's own allele never has")
    r.add_argument("--semi", action="store_true", help="O2: semi-tied reads to their best AS / fewest-edit copy, exact ties "
                   "HOPELESS (unsafe where the sample carries the other copy's base at the deciding site, register 1388)")
    r.add_argument("--anchor", action="store_true", help="O2: a one-site lead suffices at a site fixed by both copies' anchor reads")
    r.add_argument("--one-site", action="store_true", help="O2: a one-site lead suffices when the read's haplotype agrees")
    r.add_argument("--guard", action="store_true", help="O2 abstains on tied reads touching a copy with a putative flag")
    r.add_argument("--prefer", action="store_true", help="slots only on copies the haplotype's reads do not disfavour by >= 2 sites")
    r.add_argument("--no-reconstruct", action="store_true", help="O2 ignores O3's certified copies (O3 -> O2 baseline)")
    r.add_argument("--no-clipping", action="store_true", help="putative copies without the soft-clip signal (P6's port)")
    r.add_argument("--correct", action="store_true", help="milestone 2b: IsoCon-style neighbour correction of read alleles")
    r.add_argument("--no-structure", action="store_true", help="putative copies from tie patterns only (no non-canonical splicing signal)")
    r.add_argument("--v5-clique", action="store_true", help="v5: clusters below the support floor also enter the conflicting "
                   "set (they can push a real allele out of it)")
    r.add_argument("--absorb", action="store_true", help="v5.1 (rejected on development data): join clusters outside the "
                   "conflicting set to their one compatible haplotype, and a haplotype merged across piles keeps every pile's reads")
    r.add_argument("--isocon-rule", choices=["no-pile", "via-b"], default="no-pile",
                   help="which stage candidates the IsoCon side keeps: mostly reads in no phased pile (v5), or mostly via-B reads (v4)")
    r.add_argument("--dump", default="", help="JSON lines of each family's internal states (piles, clusters, haplotype "
                   "fates, private alleles); the products do not change")
    a = ap.parse_args(argv)
    global CORRECT, COVERAGE, CLIPPING, SEMI, RECONSTRUCT, PREFER, GUARD, ONE_SITE, ANCHOR
    CORRECT, COVERAGE, CLIPPING, SEMI, RECONSTRUCT, PREFER, GUARD, ONE_SITE, ANCHOR = a.correct, a.coverage, not a.no_clipping, a.semi, not a.no_reconstruct, a.prefer, a.guard, a.one_site, a.anchor
    cmd_run(a)


if __name__ == "__main__":
    main(sys.argv[1:])

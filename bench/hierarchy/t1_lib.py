#!/usr/bin/env python3
"""T1 (human, DNA level) scorer library of docs/PREREG_hierarchy_duplicon_expansion_2026-10-07.md (KEY=dupfamboundary).

Parts: 1 labels and unit strings (DupMasker colors BED -> exon labels -> rungs U0 U1 U2), 2 statistics (m1, m1c, m1t, m2, path
defect, audit, exceptions), 3 strata and the homology-preserving null, 4 verdict classes, T1-e and the universe loaders.
Plain Python plus numpy for T1-e. Every tie is sorted, so the output does not depend on PYTHONHASHSEED. Run the tests with
    python3 -B -m unittest bench/hierarchy/test_t1_lib.py
"""
import bisect
import csv
import hashlib
import json
import os
import random
from collections import Counter, defaultdict

SEED_T1 = 20260930
B_DRAWS = 10000
N_MIN = 8


# ---------------------------------------------------------------------------------------------------------------- utilities
def merge_intervals(iv):
    """Merge touching or overlapping [a, b] intervals (a <= end of the previous one merges, as the frozen code does)."""
    out = []
    for a, b in sorted((int(x[0]), int(x[1])) for x in iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def sha256_file(path, n=16):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n] if n else h.hexdigest()


def count_lines(path):
    n = 0
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            n += chunk.count(b"\n")
    return n


# ----------------------------------------------------------------------------------------------- part 1: labels, unit strings
class Colors:
    """DupMasker colors segments per chromosome, sorted; `maxlen` bounds how far upstream a segment can start and still overlap."""

    def __init__(self, segs):
        self.segs = {c: sorted(v) for c, v in segs.items()}
        self.starts = {c: [s for s, _, _ in v] for c, v in self.segs.items()}
        self.maxlen = {c: max(e - s for s, e, _ in v) for c, v in self.segs.items()}

    @classmethod
    def from_rows(cls, rows):
        segs = defaultdict(list)
        for c, s, e, d in rows:
            segs[c].append((int(s), int(e), d))
        return cls(segs)

    @classmethod
    def load(cls, path):
        segs = defaultdict(list)
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                segs[f[0]].append((int(f[1]), int(f[2]), f[3]))
        return cls(segs)


def composition(colors, chrom, a, b, union=False):
    """Label -> bases of the exon [a, b) covered. Every segment counts, so same-label overlaps add (the frozen comp()). With
    union=True the segments of one label are merged first, so a base counts once per label."""
    out = Counter()
    segs = colors.segs.get(chrom)
    if not segs:
        return out
    k = bisect.bisect_left(colors.starts[chrom], a - colors.maxlen[chrom])
    pieces = defaultdict(list)
    while k < len(segs) and segs[k][0] < b:
        s, e, d = segs[k]
        lo, hi = max(a, s), min(b, e)
        if hi > lo:
            if union:
                pieces[d].append([lo, hi])
            else:
                out[d] += hi - lo
        k += 1
    for d, iv in pieces.items():
        out[d] = sum(y - x for x, y in merge_intervals(iv))
    return out


def _num_key(label):
    digits = "".join(ch for ch in label if ch.isdigit())
    return (int(digits) if digits else -1, label)


def pick_dominant(comp, tie="lex"):
    """Dominant label(s) of an exon composition. lex: most bases, ties to the lexicographically smallest label string (the
    key of the frozen code). num: ties to the smallest numeric ID. set: every tied label. band: every label with bases >= 0.95 x the top."""
    if not comp:
        return frozenset()
    top = max(comp.values())
    if tie == "lex":
        return frozenset({min(comp.items(), key=lambda kv: (-kv[1], kv[0]))[0]})
    if tie == "num":
        return frozenset({min(comp.items(), key=lambda kv: (-kv[1], _num_key(kv[0])))[0]})
    if tie == "set":
        return frozenset(d for d, v in comp.items() if v == top)
    if tie == "band":
        return frozenset(d for d, v in comp.items() if 100 * v >= 95 * top)
    raise ValueError(f"unknown tie mode {tie!r}")


class GeneStrings:
    __slots__ = ("pooled", "exon_dom", "u0", "u1", "ordered", "any_overlap", "n_exons")

    def __init__(self, pooled, exon_dom, u0, u1, ordered, any_overlap, n_exons):
        self.pooled, self.exon_dom, self.u0, self.u1 = pooled, exon_dom, u0, u1
        self.ordered, self.any_overlap, self.n_exons = ordered, any_overlap, n_exons


def gene_strings(exons, chrom, colors, tie="lex", union=False):
    """Unit strings of one gene from its exons: U1 = the set of dominant labels over the merged exons, U0 = the label with
    most bases pooled over the gene (ties by the same key), `ordered` = dom(e) in genomic order with consecutive repeats
    collapsed and unlabelled exons skipped. A tie mode other than lex leaves U0 on the lex key when it is a set mode."""
    merged = merge_intervals(exons)
    pooled = Counter()
    exon_dom = []
    for a, b in merged:
        comp = composition(colors, chrom, a, b, union=union)
        pooled.update(comp)
        exon_dom.append(pick_dominant(comp, tie))
    u1 = frozenset(d for s in exon_dom for d in s)
    pooled_tie = tie if tie in ("lex", "num") else "lex"
    u0 = next(iter(pick_dominant(pooled, pooled_tie))) if pooled else None
    ordered = []
    for s in exon_dom:
        if not s:
            continue
        label = next(iter(s)) if len(s) == 1 else "+".join(sorted(s))
        if not ordered or ordered[-1] != label:
            ordered.append(label)
    return GeneStrings(pooled, exon_dom, u0, u1, tuple(ordered), frozenset(pooled), len(merged))


def cn_classes(values):
    """Single-linkage copy-number classes: sorted values, a new class where the step to the previous value is >= 2 (the
    comparator of soto_m2_unit_cover_cn.py is `m - previous < 2` for the same class)."""
    out = []
    for v in sorted(values):
        if out and v - out[-1][-1] < 2:
            out[-1].append(v)
        else:
            out.append([v])
    return out


def u2_index(u1_by_gene, cn_by_gene):
    """Copy-number classes per label: {d: {copy number: class index}}. The classes of d are built from every gene of `u1_by_gene` that has a copy
    number and d in its U1 set, by single linkage (cn_classes)."""
    vals = defaultdict(list)
    for g, s in u1_by_gene.items():
        c = cn_by_gene.get(g)
        if c is None:
            continue
        for d in s:
            vals[d].append(c)
    index = {}
    for d, v in vals.items():
        m = {}
        for k, cls in enumerate(cn_classes(v)):
            for x in cls:
                m[x] = k
        index[d] = m
    return index


def u2_strings(u1_by_gene, cn_by_gene):
    """Rung U2: for each gene with a copy number, the set {(d, k)} with k the index of the copy-number class of label d that contains the gene's own
    value. Genes without a copy number are outside U2 (absent from the result)."""
    index = u2_index(u1_by_gene, cn_by_gene)
    out = {}
    for g, s in u1_by_gene.items():
        c = cn_by_gene.get(g)
        if c is None or not s:
            continue
        out[g] = frozenset((d, index[d][c]) for d in s)
    return out


def u2_ordered(ordered_by_gene, cn_by_gene, index):
    """The ordered strings at rung U2: each unit of the ordered string of a gene with a copy number becomes (label, class of the gene's copy number)."""
    out = {}
    for g, seq in ordered_by_gene.items():
        c = cn_by_gene.get(g)
        if c is None:
            continue
        out[g] = tuple((d, index[d][c]) for d in seq)
    return out


# ---------------------------------------------------------------------------------------------------- part 2: statistics
def _skey(s):
    """A deterministic sort key of a unit string (a frozenset of labels or of (label, class) pairs)."""
    return tuple(sorted(s))


def _universe(fam, strings, genes=None):
    pool = set(fam) & set(strings) if genes is None else genes
    return sorted(g for g in pool if g in fam and g in strings)


def _family_counts(fam, strings, genes):
    fc = defaultdict(Counter)
    for g in genes:
        fc[fam[g]][strings[g]] += 1
    return fc


def _string_counts(fam, strings, genes):
    sc = defaultdict(Counter)
    for g in genes:
        sc[strings[g]][fam[g]] += 1
    return sc


def chain_weight(cnt):
    """W: the largest total count over chains of distinct strings under proper-subset inclusion (dynamic programme over the poset;
    its value is unique and does not depend on any tie)."""
    items = sorted(cnt, key=lambda s: (len(s), _skey(s)))
    best = {}
    top = 0
    for i, s in enumerate(items):
        b = cnt[s] + max([best[t] for t in items[:i] if t < s] + [0])
        best[s] = b
        top = max(top, b)
    return top


def m1(fam, strings, genes=None):
    """Strict D1: sum over families with at least two genes of (n_f - the count of the most frequent string)."""
    g = _universe(fam, strings, genes)
    return sum(sum(c.values()) - max(c.values()) for c in _family_counts(fam, strings, g).values() if sum(c.values()) >= 2)


def m1c(fam, strings, genes=None):
    """Chain-tolerant D1: sum over families with at least two genes of (n_f - W_f)."""
    g = _universe(fam, strings, genes)
    return sum(sum(c.values()) - chain_weight(c) for c in _family_counts(fam, strings, g).values() if sum(c.values()) >= 2)


def m2(fam, strings, genes=None):
    """Strict D2: sum over strings with at least two genes of (n_s - the largest count in one family)."""
    g = _universe(fam, strings, genes)
    return sum(sum(c.values()) - max(c.values()) for c in _string_counts(fam, strings, g).values() if sum(c.values()) >= 2)


def macro_d1(fam, strings, genes=None):
    g = _universe(fam, strings, genes)
    xs = [(sum(c.values()) - chain_weight(c)) / sum(c.values()) for c in _family_counts(fam, strings, g).values() if sum(c.values()) >= 2]
    return sum(xs) / len(xs) if xs else 0.0


def macro_d2(fam, strings, genes=None):
    g = _universe(fam, strings, genes)
    xs = [(sum(c.values()) - max(c.values())) / sum(c.values()) for c in _string_counts(fam, strings, g).values() if sum(c.values()) >= 2]
    return sum(xs) / len(xs) if xs else 0.0


def _comparable(a, b):
    """Ordered strings are comparable iff equal or one is a prefix or a suffix of the other."""
    if len(a) > len(b):
        a, b = b, a
    return b[: len(a)] == a or (len(a) > 0 and b[len(b) - len(a):] == a)


def _max_weight_clique(weight, adj):
    nodes = sorted(weight, key=lambda v: (-weight[v], v))
    best = 0

    def rec(cand, cur):
        nonlocal best
        best = max(best, cur)
        cand = list(cand)
        while cand:
            if cur + sum(weight[v] for v in cand) <= best:
                return
            v = cand.pop()
            rec([u for u in cand if u in adj[v]], cur + weight[v])

    rec(nodes, 0)
    return best


def m1t(fam, ordered, genes=None):
    """Gene-level D1 with nesting restricted to end truncation: W_f^t is the largest number of genes of the family whose ordered
    strings are pairwise comparable (equal, or one a prefix or suffix of the other); m1t = sum over families of (n_f - W_f^t)."""
    g = _universe(fam, ordered, genes)
    tot = 0
    for c in _family_counts(fam, ordered, g).values():
        n = sum(c.values())
        if n < 2:
            continue
        weight = dict(c)
        adj = {a: {b for b in weight if b != a and _comparable(a, b)} for a in weight}
        tot += n - _max_weight_clique(weight, adj)
    return tot


def path_defect(fam, ordered, strings, genes=None):
    """(Delta, Mult). H_f: simple undirected graph on the units of the members' ordered strings, edges between consecutive units.
    delta(f) = (c - 1) + sum_v max(deg(v) - 2, 0) + (|E| - |V| + c); Delta sums delta over families with at least two genes.
    Mult = sum over units u with mu(u) >= 1 of (mu(u) - 1), mu(u) = the number of scored families holding a member whose string contains u."""
    g = _universe(fam, ordered, genes)
    members = defaultdict(list)
    for x in g:
        members[fam[x]].append(x)
    delta = 0
    mu = Counter()
    for f, gs in members.items():
        if len(gs) < 2:
            continue
        verts, edges = set(), set()
        for x in gs:
            seq = ordered[x]
            verts.update(seq)
            for a, b in zip(seq, seq[1:]):
                if a != b:
                    edges.add(frozenset((a, b)))
        if verts:
            parent = {v: v for v in verts}

            def find(v):
                while parent[v] != v:
                    parent[v] = parent[parent[v]]
                    v = parent[v]
                return v

            deg = Counter()
            for e in edges:
                a, b = tuple(e)
                deg[a] += 1
                deg[b] += 1
                parent[find(a)] = find(b)
            c = len({find(v) for v in verts})
            delta += (c - 1) + sum(max(d - 2, 0) for d in deg.values()) + (len(edges) - len(verts) + c)
        held = set()
        for x in gs:
            held.update(strings.get(x, ()))
        for u in held:
            mu[u] += 1
    return delta, sum(m - 1 for m in mu.values() if m >= 1)


TIED_CAP = 200


def _best_chains(cnt, cap=None):
    """(W, chains, n_tied): the maximum chain weight, the maximum-weight chains of distinct strings listed from the smallest key (strings as sorted
    tuples, smallest string first, compared as lists) and at most `cap` of them, and the exact number of maximum-weight chains. The search runs upward
    from the smallest string of a chain, so the kept chain is found without listing every tie."""
    items = sorted(cnt, key=lambda s: (len(s), _skey(s)))
    n = len(items)
    up = [0] * n                                              # best weight of a chain whose smallest string is items[i]
    succ = [[] for _ in range(n)]                             # successors that keep the optimum
    ways = [0] * n
    for i in range(n - 1, -1, -1):
        s = items[i]
        above = [j for j in range(i + 1, n) if s < items[j]]
        up[i] = cnt[s] + max([up[j] for j in above] + [0])
        succ[i] = sorted((j for j in above if up[i] == cnt[s] + up[j]), key=lambda j: _skey(items[j]))
        ways[i] = sum(ways[j] for j in succ[i]) if succ[i] else 1
    top = max(up)
    starts = sorted((i for i in range(n) if up[i] == top), key=lambda i: _skey(items[i]))
    chains = []

    def walk(i, prefix):
        if cap is not None and len(chains) >= cap:
            return
        prefix = prefix + [items[i]]
        if not succ[i]:
            chains.append(prefix)
            return
        for j in succ[i]:
            walk(j, prefix)

    for i in starts:
        walk(i, [])
    return top, chains, sum(ways[i] for i in starts)


def _label_gene(g, nex_info, kept_genes):
    """short: at most two merged exons. ends-only: at least three merged exons and the set of dominant labels over all merged exons
    except the first and the last equals that set of a kept gene with at least three exons. unexplained: the rest."""
    inf = nex_info[g]
    if inf["nex"] <= 2:
        return "short"
    if any(nex_info[h]["nex"] >= 3 and nex_info[h]["interior"] == inf["interior"] for h in kept_genes):
        return "ends-only"
    return "unexplained"


def audit_d1(fam, strings, info, genes=None, cap=TIED_CAP):
    """Per family with at least two genes: W, the kept chain (maximum weight; ties to the lexicographically smallest key, the key being
    the strings of the chain as sorted tuples from the smallest to the largest compared as lists), every tied chain with its deficit
    genes, the deficit genes of the kept chain and their audit labels. The audit never changes a statistic."""
    g = _universe(fam, strings, genes)
    members = defaultdict(list)
    for x in g:
        members[fam[x]].append(x)
    out = {}
    for f in sorted(members):
        gs = members[f]
        if len(gs) < 2:
            continue
        cnt = Counter(strings[x] for x in gs)
        top, chains, n_tied = _best_chains(cnt, cap=cap)
        tied = []
        for ch in chains:
            inside = set(ch)
            tied.append(dict(chain=[_skey(s) for s in ch], deficit_genes=sorted(x for x in gs if strings[x] not in inside)))
        kept = set(chains[0])
        kept_genes = sorted(x for x in gs if strings[x] in kept)
        deficit = sorted(x for x in gs if strings[x] not in kept)
        out[f] = dict(n=len(gs), W=top, kept_chain=[_skey(s) for s in chains[0]], tied_chains=tied, n_tied_chains=n_tied, tied_truncated=n_tied > len(chains),
                      deficit_genes=deficit,
                      labels={x: _label_gene(x, info, kept_genes) for x in deficit})
    return out


def audit_d2(fam, strings, info, genes=None):
    """Per string with at least two genes: the best family (most genes; ties to the smaller family id), the deficit genes (the genes of the
    string outside it) and their labels (short or unexplained; ends-only is a D1 label)."""
    g = _universe(fam, strings, genes)
    out = {}
    sc = _string_counts(fam, strings, g)
    for s in sorted(sc, key=_skey):
        cnt = sc[s]
        n = sum(cnt.values())
        if n < 2:
            continue
        best = min(cnt.items(), key=lambda kv: (-kv[1], str(kv[0])))[0]
        tied = sorted((f for f, v in cnt.items() if v == cnt[best]), key=str)
        deficit = sorted(x for x in g if strings[x] == s and fam[x] != best)
        out[_skey(s)] = dict(n=n, best_family=best, tied_best_families=tied, deficit_genes=deficit,
                             labels={x: ("short" if info[x]["nex"] <= 2 else "unexplained") for x in deficit})
    return out


def exceptions_d1(fam, strings, info, genes=None):
    """The deficit genes of the kept chains of all families, sorted by (family, gene): the exception list of T1-a for D1 (length = m1c)."""
    rows = []
    for f, rec in audit_d1(fam, strings, info, genes).items():
        for x in rec["deficit_genes"]:
            rows.append(dict(gene=x, family=f, string=_skey(strings[x]), label=rec["labels"][x], kept_chain=rec["kept_chain"]))
    return sorted(rows, key=lambda r: (str(r["family"]), r["gene"]))


def exceptions_d2(fam, strings, info, genes=None):
    """The deficit genes of all strings with at least two genes, sorted by (string, gene): the exception list of T1-a for D2 (length = m2)."""
    rows = []
    for s, rec in audit_d2(fam, strings, info, genes).items():
        for x in rec["deficit_genes"]:
            rows.append(dict(gene=x, family=fam[x], string=s, best_family=rec["best_family"], label=rec["labels"][x]))
    return sorted(rows, key=lambda r: (r["string"], r["gene"]))


# --------------------------------------------------------------------------------- part 3: strata, the null, verdict classes
def graph_components(edges):
    """Connected components of an edge list [(u, v, w)]; a component is labelled by its smallest node string."""
    parent = {}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for e in edges:
        u, v = e[0], e[1]
        parent.setdefault(u, u)
        parent.setdefault(v, v)
        parent[find(u)] = find(v)
    groups = defaultdict(list)
    for n in parent:
        groups[find(n)].append(n)
    comp = {}
    for nodes in groups.values():
        rid = min(nodes)
        for n in nodes:
            comp[n] = rid
    return comp


def build_strata(vf_genes, gene_node, comp):
    """Strata of the decisive universe: V_F genes grouped by the component of their node. A stratum is named after its smallest gene id
    and the list is ordered by that id (string order). Genes without a node are outside every stratum."""
    by = defaultdict(list)
    for g in vf_genes:
        n = gene_node.get(g)
        if n is None or n not in comp:
            continue
        by[comp[n]].append(g)
    strata = [(min(gs), sorted(gs)) for gs in by.values()]
    strata.sort()
    return strata


def variation_strata(fam, strata):
    """The strata that hold at least two distinct family labels: the only ones where the null can move."""
    return [(sid, genes) for sid, genes in strata if len({fam[g] for g in genes}) >= 2]


def p_lower(obs, draws):
    """p = (1 + number of draws at least as extreme as observed, that is not larger than it) / (B + 1)."""
    return (1 + sum(1 for d in draws if d <= obs)) / (len(draws) + 1)


CLASS_ORDER = ["NOT-DISTINGUISHED", "PARTIAL", "CUT-ALIGNED"]


def content_class(p):
    if p <= 0.01:
        return "CUT-ALIGNED"
    if p <= 0.05:
        return "PARTIAL"
    return "NOT-DISTINGUISHED"


def lower_class(a, b):
    return CLASS_ORDER[min(CLASS_ORDER.index(a), CLASS_ORDER.index(b))]


def classify_cut(p_all, p_removed, n_var, n_live, valid=True, n_min=N_MIN):
    """Verdict algebra of T1-b and T1-c in the registered order. INVALID; UNDERPOWERED when fewer than n_min strata are variation-bearing;
    UNINFORMATIVE when there are at least n_min of them and none is live; UNDERPOWERED when fewer than n_min strata are live; else the
    content class, the lower of the class of p_all and the class of p_removed (the class of the p of the remaining strata alone; without a
    removal, the class of p_all). Dead or live is counted by the caller from the draws."""
    ca = content_class(p_all)
    cr = content_class(p_removed) if p_removed is not None else ca
    if not valid:
        cls = "INVALID"
    elif n_var < n_min:
        cls = "UNDERPOWERED"
    elif n_live == 0:
        cls = "UNINFORMATIVE"
    elif n_live < n_min:
        cls = "UNDERPOWERED"
    else:
        cls = lower_class(ca, cr)
    return dict(**{"class": cls}, class_all=ca, class_removed=cr)


def pick_removed(z):
    """The stratum with the most negative z (ties: the smaller stratum id); None when no z is negative (the removal step is skipped)."""
    neg = {k: v for k, v in z.items() if v < 0}
    if not neg:
        return None
    return min(neg.items(), key=lambda kv: (kv[1], kv[0]))[0]


class NullEngine:
    """The homology-preserving null of section 4: within each variation-bearing stratum the cluster labels of its V_F genes are permuted
    (the list of observed labels in gene-id order is shuffled and assigned back in that order), strata in the order of their smallest
    gene id, one shuffle per stratum per draw, random.Random(seed) consumed in that order. The statistics m1c and m2 are recomputed on
    the full label vector; families outside the variation-bearing strata keep their labels, so their part is constant.

    Draw b of a block that starts at `start` is the same draw b of the whole run: the first `start` draws are replayed without scoring."""

    def __init__(self, fam, strings, strata, seed=SEED_T1):
        self.fam, self.strings, self.seed = fam, strings, seed
        self.strata = [(sid, list(genes)) for sid, genes in strata]
        stratum_of_family = {}
        for sid, genes in self.strata:
            for g in genes:
                if stratum_of_family.setdefault(fam[g], sid) != sid:
                    raise ValueError(f"family {fam[g]!r} spans strata {stratum_of_family[fam[g]]!r} and {sid!r}")
        self.var = variation_strata(fam, self.strata)
        var_sids = {sid for sid, _ in self.var}
        self.sets = []
        self._id = {}

        def sid_of(s):
            if s not in self._id:
                self._id[s] = len(self.sets)
                self.sets.append(s)
            return self._id[s]

        self.vlab0 = [[fam[g] for g in genes] for _, genes in self.var]
        self.vstr = [[sid_of(strings[g]) for g in genes] for _, genes in self.var]
        for _, genes in self.strata:
            for g in genes:
                sid_of(strings[g])
        self._prepare_chains()
        # fixed part: families in strata where the null cannot move
        self.n_str = Counter()
        self.fix_max = {}
        self.fix_m1c = 0
        for sid, genes in self.strata:
            for g in genes:
                self.n_str[sid_of(strings[g])] += 1
            if sid in var_sids:
                continue
            fc = defaultdict(Counter)
            for g in genes:
                fc[fam[g]][sid_of(strings[g])] += 1
            for c in fc.values():
                n = sum(c.values())
                if n >= 2:
                    self.fix_m1c += n - self._w(c)
                for k, v in c.items():
                    if v > self.fix_max.get(k, 0):
                        self.fix_max[k] = v
        self.var_strings = {k for row in self.vstr for k in row}
        self.fix_m2_const = sum(self.n_str[k] - self.fix_max.get(k, 0) for k in self.n_str if k not in self.var_strings and self.n_str[k] >= 2)

    def _prepare_chains(self):
        """Proper-subset relations among the strings that occur, as lists of smaller string ids, and a rank that orders strings by size."""
        n = len(self.sets)
        order = sorted(range(n), key=lambda k: (len(self.sets[k]), _skey(self.sets[k])))
        self.rank = [0] * n
        for r, k in enumerate(order):
            self.rank[k] = r
        self.below = [[] for _ in range(n)]
        for i, a in enumerate(order):
            sa = self.sets[a]
            for b in order[:i]:
                if self.sets[b] < sa:
                    self.below[a].append(b)

    def _w(self, c):
        """W for a family given as {string id: count}: the longest-weight chain under proper-subset inclusion."""
        if len(c) == 1:
            return next(iter(c.values()))
        best = {}
        top = 0
        for k in sorted(c, key=self.rank.__getitem__):
            m = 0
            for j in self.below[k]:
                b = best.get(j)
                if b is not None and b > m:
                    m = b
            v = c[k] + m
            best[k] = v
            if v > top:
                top = v
        return top

    def _score(self, label_lists):
        """(total m1c, total m2, per-stratum m1c, per-stratum m2) for one assignment of labels to the variation-bearing strata."""
        tot_m1c = self.fix_m1c
        varmax = {}
        st_m1c, st_m2 = [], []
        for labs, strs in zip(label_lists, self.vstr):
            fcs = {}
            for lab, k in zip(labs, strs):
                c = fcs.get(lab)
                if c is None:
                    c = fcs[lab] = {}
                c[k] = c.get(k, 0) + 1
            a = 0
            s_tot, s_max = {}, {}
            for c in fcs.values():
                n = sum(c.values())
                if n >= 2:
                    a += n - self._w(c)
                for k, v in c.items():
                    s_tot[k] = s_tot.get(k, 0) + v
                    if v > s_max.get(k, 0):
                        s_max[k] = v
                    if v > varmax.get(k, 0):
                        varmax[k] = v
            st_m1c.append(a)
            st_m2.append(sum(t - s_max[k] for k, t in s_tot.items() if t >= 2))
            tot_m1c += a
        tot_m2 = self.fix_m2_const
        for k in self.var_strings:
            if self.n_str[k] >= 2:
                tot_m2 += self.n_str[k] - max(self.fix_max.get(k, 0), varmax.get(k, 0))
        return tot_m1c, tot_m2, st_m1c, st_m2

    def observed(self):
        a, b, sa, sb = self._score(self.vlab0)
        return dict(m1c=a, m2=b, st_m1c=sa, st_m2=sb)

    def _permuted(self, rng):
        out = []
        for lab0 in self.vlab0:
            perm = lab0[:]
            rng.shuffle(perm)
            out.append(perm)
        return out

    def draws(self, B, start=0, end=None):
        """Scores of draws start..end-1 of a run of B draws."""
        end = B if end is None else end
        rng = random.Random(self.seed)
        for _ in range(start):
            self._permuted(rng)
        res = dict(tot_m1c=[], tot_m2=[], st_m1c=[], st_m2=[])
        for _ in range(start, end):
            a, b, sa, sb = self._score(self._permuted(rng))
            res["tot_m1c"].append(a)
            res["tot_m2"].append(b)
            res["st_m1c"].append(sa)
            res["st_m2"].append(sb)
        return res

    def naive_draws(self, B):
        """The full label dictionary of each draw (reference for the tests)."""
        rng = random.Random(self.seed)
        for _ in range(B):
            fam_b = dict(self.fam)
            for (sid, genes), perm in zip(self.var, self._permuted(rng)):
                fam_b.update(zip(genes, perm))
            yield fam_b

    def draws_without(self, sid, B, start=0, end=None):
        """m1c and m2 of each draw on the genes outside stratum `sid` (same random stream)."""
        end = B if end is None else end
        drop = {g for s, genes in self.strata if s == sid for g in genes}
        keep = [(s, genes) for s, genes in self.strata if s != sid]
        sub = NullEngine({g: f for g, f in self.fam.items() if g not in drop}, self.strings, keep, self.seed)
        pos = [j for j, (s, _) in enumerate(self.var) if s != sid]
        rng = random.Random(self.seed)
        for _ in range(start):
            self._permuted(rng)
        out = dict(m1c=[], m2=[])
        for _ in range(start, end):
            perm = self._permuted(rng)
            a, b, _, _ = sub._score([perm[j] for j in pos])
            out["m1c"].append(a)
            out["m2"].append(b)
        out["observed"] = sub.observed()
        return out


def _cell(eng, obs, dr, stat, B, ranks, keep_draws=False):
    """One decisive cell from the observed statistics and the draws of an engine."""
    tot = dr["tot_" + stat]
    st_obs = obs["st_" + stat]
    n_var = len(eng.var)
    cols = [[row[j] for row in dr["st_" + stat]] for j in range(n_var)]
    z, rows, live = {}, [], []
    for j, (sid, genes) in enumerate(eng.var):
        mean = sum(cols[j]) / len(cols[j])
        sd = (sum((x - mean) ** 2 for x in cols[j]) / len(cols[j])) ** 0.5
        z[sid] = 0.0 if sd == 0 else (st_obs[j] - mean) / sd
        is_live = len(set(cols[j])) >= 2
        if is_live:
            live.append(sid)
        rows.append(dict(stratum=sid, rank=ranks.get(sid), n_genes=len(genes), n_clusters=len({eng.fam[g] for g in genes}),
                         obs=st_obs[j], mean=mean, sd=sd, z=z[sid], live=is_live))
    p_all = p_lower(obs[stat], tot)
    removed = pick_removed(z)
    p_removed = rem_obs = None
    if removed is not None:
        j = [sid for sid, _ in eng.var].index(removed)
        if stat == "m1c":
            rem_draws = [t - r[j] for t, r in zip(tot, dr["st_m1c"])]
            rem_obs = obs["m1c"] - st_obs[j]
        else:
            w = eng.draws_without(removed, B)
            rem_draws, rem_obs = w["m2"], w["observed"]["m2"]
        p_removed = p_lower(rem_obs, rem_draws)
    res = classify_cut(p_all, p_removed, n_var, len(live))
    res.update(stat=stat, obs=obs[stat], p=p_all, p_upper=(1 + sum(1 for d in tot if d >= obs[stat])) / (len(tot) + 1),
               p_removed=p_removed, obs_removed=rem_obs, removed_stratum=removed, removed_rank=ranks.get(removed) if removed else None,
               z=z, ranks={sid: ranks.get(sid) for sid, _ in eng.var}, z_table=rows, n_var=n_var, n_live=len(live), live=live, B=B,
               draw_mean=sum(tot) / len(tot), draw_min=min(tot), draw_max=max(tot))
    if keep_draws:
        res["draws"] = tot
    return res


def run_cells(fam, strings, strata, B=B_DRAWS, seed=SEED_T1, stats=("m1c", "m2"), keep_draws=False):
    """The decisive cells of a universe from one set of draws (cells that share a universe share draws): observed statistic, lower-tail p
    over B draws, the z table over the variation-bearing strata, the number of live strata, the removal of the stratum with the most
    negative z (skipped when none is negative) with p recomputed on the same draws, and the issued class, for each of `stats`."""
    eng = NullEngine(fam, strings, strata, seed)
    obs = eng.observed()
    dr = eng.draws(B)
    ranks = {sid: i + 1 for i, (sid, _) in enumerate(strata)}
    return {stat: _cell(eng, obs, dr, stat, B, ranks, keep_draws) for stat in stats}


def run_cell(fam, strings, strata, stat, B=B_DRAWS, seed=SEED_T1, keep_draws=False):
    """One decisive cell (T1-b: stat='m1c', T1-c: stat='m2')."""
    return run_cells(fam, strings, strata, B=B, seed=seed, stats=(stat,), keep_draws=keep_draws)[stat]


def synthetic_world(seed, n_strata=12, mode="independent", sizes=None, nclu=(2, 4), vocab="ABCDEFGH"):
    """A made-up universe for the controls of Gate 4: strata of genes, clusters inside each stratum and unit strings, none of it real.
    Modes: independent (strings unrelated to the cut), follow / identity (one string per cluster), merge (neighbouring clusters share a
    string), split (each cluster carries two random sub-strings), one_string_per_stratum, one_cluster_per_stratum."""
    rng = random.Random(seed)
    fam, strings, strata, info = {}, {}, [], {}
    for i in range(n_strata):
        n = sizes[i] if sizes else rng.randint(6, 20)
        k = 1 if mode == "one_cluster_per_stratum" else min(rng.randint(*nclu), n)
        base = [2] * k if n >= 2 * k else [1] * k
        for _ in range(n - sum(base)):
            base[rng.randrange(k)] += 1
        labels = [f"MCL{i:03d}_{c}" for c, sz in enumerate(base) for _ in range(sz)]
        rng.shuffle(labels)
        genes = [f"g{i:03d}_{j:03d}" for j in range(n)]
        for g, lab in zip(genes, labels):
            fam[g] = lab
            c = int(lab.rsplit("_", 1)[1])
            if mode in ("follow", "identity"):
                s = frozenset({f"S{i:03d}_{c}"})
            elif mode == "merge":
                s = frozenset({f"S{i:03d}_{c // 2}"})
            elif mode == "split":
                s = frozenset({f"S{i:03d}_{c}_{rng.randint(0, 1)}"})
            elif mode == "one_string_per_stratum":
                s = frozenset({f"T{i:03d}"})
            else:
                s = frozenset(rng.sample(vocab, rng.randint(1, 3)))
            strings[g] = s
            info[g] = dict(nex=rng.randint(1, 8), interior=frozenset())
        strata.append((genes[0], genes))
    return dict(fam=fam, strings=strings, strata=strata, info=info)


class ChainTable:
    """Proper-subset relations among the unit strings that occur, and the chain weight W of a family given as {string id: count}."""

    def __init__(self, strings):
        self.sets, self._id = [], {}
        for s in strings:
            self.sid(s)
        order = sorted(range(len(self.sets)), key=lambda k: (len(self.sets[k]), _skey(self.sets[k])))
        self.rank = [0] * len(self.sets)
        for r, k in enumerate(order):
            self.rank[k] = r
        self.below = [[] for _ in self.sets]
        for i, a in enumerate(order):
            for b in order[:i]:
                if self.sets[b] < self.sets[a]:
                    self.below[a].append(b)

    def sid(self, s):
        if s not in self._id:
            self._id[s] = len(self.sets)
            self.sets.append(s)
        return self._id[s]

    def w(self, c):
        if len(c) == 1:
            return next(iter(c.values()))
        best, top = {}, 0
        for k in sorted(c, key=self.rank.__getitem__):
            m = 0
            for j in self.below[k]:
                b = best.get(j)
                if b is not None and b > m:
                    m = b
            v = c[k] + m
            best[k] = v
            if v > top:
                top = v
        return top


class GeneralNull:
    """The permutation null for strata that families may span (chromosome, SD98 region, F_Q group): the same draw procedure as NullEngine (strata in the
    order given, which is the order of their smallest gene id; one shuffle of the observed label list in gene-id order per stratum that holds at
    least two distinct labels; random.Random(seed) consumed in that order), but m1c and m2 are recomputed on the whole label vector at each draw.
    Slower than NullEngine and used only for the sensitivity nulls, which carry no class."""

    def __init__(self, fam, strings, strata, seed=SEED_T1):
        self.fam, self.strings, self.seed = fam, strings, seed
        self.strata = [(sid, list(genes)) for sid, genes in strata]
        genes = [g for _, gs in self.strata for g in gs]
        if len(set(genes)) != len(genes):
            raise ValueError("a gene is in two strata")
        self.genes = sorted(genes)
        index = {g: i for i, g in enumerate(self.genes)}
        self.chains = ChainTable([strings[g] for g in self.genes])
        self.sid = [self.chains.sid(strings[g]) for g in self.genes]
        self.n_str = Counter(self.sid)
        self.var = [(sid, gs) for sid, gs in self.strata if len({fam[g] for g in gs}) >= 2]
        self.lab0 = [[fam[g] for g in gs] for _, gs in self.var]
        self.pos = [[index[g] for g in gs] for _, gs in self.var]
        self.base = [fam[g] for g in self.genes]

    def _score(self, labels):
        fams = {}
        for lab, k in zip(labels, self.sid):
            c = fams.get(lab)
            if c is None:
                c = fams[lab] = {}
            c[k] = c.get(k, 0) + 1
        m1c, best = 0, {}
        for c in fams.values():
            n = sum(c.values())
            if n >= 2:
                m1c += n - self.chains.w(c)
            for k, v in c.items():
                if v > best.get(k, 0):
                    best[k] = v
        return m1c, sum(self.n_str[k] - best[k] for k in self.n_str if self.n_str[k] >= 2)

    def observed(self):
        a, b = self._score(self.base)
        return dict(m1c=a, m2=b)

    def _labels(self, rng):
        labels = self.base[:]
        for lab0, pos in zip(self.lab0, self.pos):
            perm = lab0[:]
            rng.shuffle(perm)
            for i, lab in zip(pos, perm):
                labels[i] = lab
        return labels

    def draws(self, B, start=0, end=None):
        end = B if end is None else end
        rng = random.Random(self.seed)
        for _ in range(start):
            self._labels(rng)
        res = dict(tot_m1c=[], tot_m2=[])
        for _ in range(start, end):
            a, b = self._score(self._labels(rng))
            res["tot_m1c"].append(a)
            res["tot_m2"].append(b)
        return res

    def naive_draws(self, B):
        """The full label dictionary of each draw (reference for the tests)."""
        rng = random.Random(self.seed)
        for _ in range(B):
            fam_b = dict(self.fam)
            for (sid, genes), lab0 in zip(self.var, self.lab0):
                perm = lab0[:]
                rng.shuffle(perm)
                fam_b.update(zip(genes, perm))
            yield fam_b


def run_general_cells(fam, strings, strata, B=B_DRAWS, seed=SEED_T1):
    """The sensitivity null: p of m1c and of m2 (lower tail, upper tail printed) under GeneralNull. No class, no robustness clause."""
    eng = GeneralNull(fam, strings, strata, seed)
    obs = eng.observed()
    dr = eng.draws(B)
    out = {}
    for stat in ("m1c", "m2"):
        tot = dr["tot_" + stat]
        out[stat] = dict(stat=stat, obs=obs[stat], p=p_lower(obs[stat], tot), p_upper=(1 + sum(1 for d in tot if d >= obs[stat])) / (len(tot) + 1), B=B,
                         n_strata=len(eng.strata), n_var=len(eng.var), draw_mean=sum(tot) / len(tot), draw_min=min(tot), draw_max=max(tot))
    return out


# --------------------------------------------------------------------------- part 3b: helpers of the gates and the descriptive arms
def d_partition(blocks, cell_of, ignore_ungrouped=False):
    """D[X -> Y] = sum over blocks x of X with at least two genes of (|x| - max over cells y of |x & y|). cell_of(g) is the cell of Y that holds g, or
    None for a gene that Y leaves ungrouped: such genes are singleton cells of their own, or leave the block first when ignore_ungrouped
    (the at-least-two test is then applied to the grouped genes). Returns the genes to move, the blocks split over more than one cell and the
    number of blocks scored."""
    moved = split = scored = 0
    for b in blocks:
        genes = [g for g in blocks[b] if not (ignore_ungrouped and cell_of(g) is None)]
        if len(genes) < 2:
            continue
        cnt = Counter(cell_of(g) if cell_of(g) is not None else ("single", g) for g in genes)
        moved += len(genes) - max(cnt.values())
        split += len(cnt) > 1
        scored += 1
    return dict(moved=moved, split=split, blocks=scored)


def exon_label_rows(colors, chrom, exons, union=False):
    """For each merged exon: whether it has a label, whether its two labels with most bases tie exactly, and the dominant label under the lex and
    the numeric-ID tie rules."""
    out = []
    for a, b in exons:
        comp = composition(colors, chrom, a, b, union=union)
        if not comp:
            out.append(dict(labelled=False, tie=False, lex=None, num=None))
            continue
        vals = sorted(comp.values(), reverse=True)
        out.append(dict(labelled=True, tie=len(vals) >= 2 and vals[0] == vals[1], lex=next(iter(pick_dominant(comp, "lex"))),
                        num=next(iter(pick_dominant(comp, "num")))))
    return out


def cluster_sizes(fam, strata):
    """The cluster sizes (descending) inside each variation-bearing stratum, in stratum order. Family-only information."""
    return [sorted(Counter(fam[g] for g in genes).values(), reverse=True) for _, genes in variation_strata(fam, strata)]


def world_from_sizes(sizes, mode="independent", seed=0, vocab="ABCDEFGH"):
    """A made-up universe on a given shape: sizes[i] lists the cluster sizes of stratum i (sizes only, nothing real). Modes: independent
    (strings unrelated to the cut), identity (one string per cluster), pairs (adjacent clusters share a string, only in strata with at least
    three clusters, identity elsewhere), split (a cluster of at least three genes carries two alternating types), largest_follows (the stratum
    with most genes follows the strings, ties to the smaller id, the rest is independent)."""
    rng = random.Random(seed)
    totals = [sum(cl) for cl in sizes]
    largest = totals.index(max(totals)) if totals else None
    fam, strings, strata, info = {}, {}, [], {}
    for i, cl in enumerate(sizes):
        genes, j = [], 0
        for c, n in enumerate(sorted(cl, reverse=True)):
            for pos in range(n):
                g = f"g{i:03d}_{j:03d}"
                j += 1
                genes.append(g)
                fam[g] = f"MCL{i:03d}_{c}"
                m = mode
                if mode == "largest_follows":
                    m = "identity" if i == largest else "independent"
                if m == "independent":
                    s = frozenset(rng.sample(vocab, rng.randint(1, 3)))
                elif m == "identity":
                    s = frozenset({f"S{i:03d}_{c}"})
                elif m == "pairs":
                    s = frozenset({f"S{i:03d}_{c // 2}"}) if len(cl) >= 3 else frozenset({f"S{i:03d}_{c}"})
                elif m == "split":
                    s = frozenset({f"S{i:03d}_{c}_{pos % 2}"}) if n >= 3 else frozenset({f"S{i:03d}_{c}"})
                else:
                    raise ValueError(f"unknown mode {mode!r}")
                strings[g] = s
                info[g] = dict(nex=rng.randint(1, 8), interior=frozenset())
        strata.append((genes[0], genes))
    return dict(fam=fam, strings=strings, strata=strata, info=info)


def same_locus_pairs(info, genes):
    """Same-strand pairs of the listed genes whose spans overlap (a start before the other's end), and the gene deleted from each pair: the one
    with the smaller exon union, ties to the larger gene id. The pairs are independent of each other. Returns (sorted pairs, set of deleted genes)."""
    by = defaultdict(list)
    for g in genes:
        by[(info[g]["chrom"], info[g]["strand"])].append(g)
    size = {g: sum(b - a for a, b in merge_intervals(info[g]["exons"])) for g in genes}
    pairs, losers = [], set()
    for key in sorted(by):
        order = sorted(by[key], key=lambda g: (info[g]["start"], info[g]["end"], g))
        for i, a in enumerate(order):
            for b in order[i + 1:]:
                if info[b]["start"] >= info[a]["end"]:
                    break
                pairs.append(tuple(sorted((a, b))))
                x, y = pairs[-1]
                losers.add(x if size[x] < size[y] else y)
    return sorted(pairs), losers


# the port of the arms of soto_m2_unit_cover_cn.py (Gate 1): the same votes, units and owners; the permutation p-values are not ported
def _fam_num(f):
    return int(f.split("_")[1])


def _jacc(s, p):
    u = s | p
    return len(s & p) / len(u) if u else 0.0


def _argmax_family(cnt):
    cnt = +cnt
    return min(cnt.items(), key=lambda kv: (-kv[1], _fam_num(kv[0])))[0] if cnt else None


def _cn_classes_members(members):
    """members: [(cn, family, bp)] -> classes (lists of members), cut where consecutive copy numbers differ by 2 or more."""
    out = []
    for m in sorted(members, key=lambda m: m[0]):
        if out and m[0] - out[-1][-1][0] < 2:
            out[-1].append(m)
        else:
            out.append([m])
    return out


def unit_cover(genes, arm, cnsrc):
    """Arm A (any-overlap usage), B (a duplicon is used only if dominant on an exon), C (A with duplicon x copy-number class units), D (B with them).
    genes: dicts with id, S (non-Unassigned Soto families), comp (Counter of bases per label), dom (set of dominant labels), cn ({'s1c', 'ours'}),
    half. Returns the mean Jaccard per half, exact / recall / precision over the multi-family genes, the clean-gene leave-one-out shares, the
    prediction of every multi-family gene and the number of units."""
    clean = [g for g in genes if len(g["S"]) == 1]
    multi = [g for g in genes if len(g["S"]) >= 2]
    sliver_free, use_cn = arm in ("B", "D"), arm in ("C", "D")
    used = lambda g: [d for d in g["comp"] if (d in g["dom"] or not sliver_free)]
    votes = defaultdict(list)
    for g in clean:
        F = next(iter(g["S"]))
        for d in used(g):
            votes[d].append((g, F, g["comp"][d], g["cn"][cnsrc]))

    def units_of_D(d, minus=None):
        vs = [v for v in votes.get(d, []) if v[0] is not minus]
        if not vs:
            return []
        withcn = [(v[3], v[1], v[2]) for v in vs if v[3] is not None]
        if not use_cn or not withcn:
            cnt = Counter()
            for v in vs:
                cnt[v[1]] += v[2]
            o = _argmax_family(cnt)
            return [((d, "all"), o, None)] if o else []
        out = []
        for k, cls in enumerate(_cn_classes_members(withcn)):
            cnt = Counter()
            for cn, F, bp in cls:
                cnt[F] += bp
            o = _argmax_family(cnt)
            if o:
                out.append(((d, k), o, [m[0] for m in cls]))
        return out

    def gene_units(g, minus=None):
        out = []
        c = g["cn"][cnsrc]
        for d in used(g):
            for key, o, cns in units_of_D(d, minus):
                if cns is None or c is None or any(abs(c - x) < 2 for x in cns):
                    out.append((key, o))
        return out

    owner = {}
    for d in votes:
        for key, o, _ in units_of_D(d):
            owner[key] = o
    pred = {g["id"]: {owner[k] for k, _ in gene_units(g) if k in owner} for g in multi}
    obs = {}
    for h in ("dev", "heldout"):
        xs = [_jacc(g["S"], pred[g["id"]]) for g in multi if g["half"] == h]
        obs[h] = (sum(xs) / len(xs) if xs else 0.0, len(xs))
    prs = [len(g["S"] & pred[g["id"]]) / len(pred[g["id"]]) for g in multi if pred[g["id"]]]
    ok = mul = emp = 0
    for g in clean:
        P = {o for _, o in gene_units(g, minus=g)}
        ok += P == {next(iter(g["S"]))}
        mul += len(P) >= 2
        emp += not P
    n = len(multi)
    return dict(obs=obs, exact=sum(g["S"] == pred[g["id"]] for g in multi) / n if n else 0.0,
                recall=sum(len(g["S"] & pred[g["id"]]) / len(g["S"]) for g in multi) / n if n else 0.0,
                precision=sum(prs) / len(prs) if prs else 0.0, clean_ok=ok / len(clean), clean_multi=mul / len(clean),
                clean_empty=emp / len(clean), pred=pred, n_units=len(owner))


# ------------------------------------------------------------------------------ part 4: T1-e (descriptive, independent evidence)
class Regions:
    """SD98 regions (chrom, start, end) per chromosome. The region id is the tuple itself, so the order of ids is the tuple order."""

    def __init__(self, rows):
        by = defaultdict(list)
        for c, a, b in rows:
            by[c].append((int(a), int(b)))
        self.by = {c: sorted(v) for c, v in by.items()}
        self.starts = {c: [a for a, _ in v] for c, v in self.by.items()}
        self.maxlen = {c: max(b - a for a, b in v) for c, v in self.by.items()}

    @classmethod
    def load(cls, path):
        rows = []
        with open(path) as fh:
            for line in fh:
                f = line.split()
                if len(f) >= 3:
                    rows.append((f[0], int(f[1]), int(f[2])))
        return cls(rows)


def best_region(regions, chrom, exons):
    """The region with the largest total overlap with the gene's exons; ties go to the smaller region id (chrom, start, end)."""
    v = regions.by.get(chrom)
    if not v:
        return None
    ov = Counter()
    for x, y in exons:
        k = bisect.bisect_left(regions.starts[chrom], x - regions.maxlen[chrom])
        while k < len(v) and v[k][0] < y:
            a, b = v[k]
            o = min(y, b) - max(x, a)
            if o > 0:
                ov[(chrom, a, b)] += o
            k += 1
    return min(ov.items(), key=lambda kv: (-kv[1], kv[0]))[0] if ov else None


def average_ranks(values):
    """1-based ranks, tied values share the average of their ranks."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        for k in range(i, j + 1):
            out[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return out


def _stratum_blocks(rows, y, xs, strata=None):
    """Per-stratum demeaned regressors and outcome. `strata` lists the strata to resample (default: those that hold rows); a listed stratum
    without rows is an empty block."""
    import numpy as np
    by = defaultdict(list)
    for r in rows:
        by[r["stratum"]].append(r)
    ids = sorted(by) if strata is None else list(strata)
    blocks = []
    for sid in ids:
        rs = by.get(sid, [])
        if not rs:
            blocks.append((sid, np.zeros((0, len(xs))), np.zeros(0)))
            continue
        X = np.array([[r[k] for k in xs] for r in rs], dtype=float)
        yv = np.array([r[y] for r in rs], dtype=float)
        blocks.append((sid, X - X.mean(axis=0), yv - yv.mean()))
    return blocks


def fe_ols(rows, y, xs):
    """Ordinary least squares of y on xs after subtracting the stratum mean of every variable (stratum fixed effects, no intercept).
    Collinear regressors give the minimum-norm solution and rank_deficient = True."""
    import numpy as np
    blocks = _stratum_blocks(rows, y, xs)
    X = np.vstack([b[1] for b in blocks])
    yv = np.concatenate([b[2] for b in blocks])
    beta, _, rank, _ = np.linalg.lstsq(X, yv, rcond=None)
    return dict(coef=dict(zip(xs, (float(b) for b in beta))), n=len(yv), strata=len(blocks), rank=int(rank), rank_deficient=bool(rank < len(xs)))


def percentile_interval(vals, B):
    """95% percentile interval by the nearest rank: the round(0.025 B)-th and round(0.975 B)-th smallest replicates."""
    v = sorted(vals)
    lo = max(int(round(0.025 * B)) - 1, 0)
    hi = min(int(round(0.975 * B)) - 1, len(v) - 1)
    return v[lo], v[hi]


def bootstrap_fe(rows, y, xs, B=2000, seed=SEED_T1, strata=None):
    """Bootstrap over strata: the listed strata (default: those with rows, sorted by id; the registration resamples all 30 variation-bearing strata,
    those without a pair included) are resampled with replacement, the model is refitted on the normal equations (a resampled stratum counts as
    many times as it is drawn). Returns the point estimate, the replicates and the percentile intervals."""
    import numpy as np
    blocks = _stratum_blocks(rows, y, xs, strata)
    K = len(blocks)
    XtX = [b[1].T @ b[1] for b in blocks]
    Xty = [b[1].T @ b[2] for b in blocks]
    point = np.linalg.lstsq(sum(XtX), sum(Xty), rcond=None)[0]
    rng = random.Random(seed)
    reps = []
    for _ in range(B):
        idx = [rng.randrange(K) for _ in range(K)]
        a = sum(XtX[i] for i in idx)
        b = sum(Xty[i] for i in idx)
        reps.append(np.linalg.lstsq(a, b, rcond=None)[0])
    out = dict(point={k: float(point[j]) for j, k in enumerate(xs)}, B=B, replicates={}, interval={})
    for j, k in enumerate(xs):
        col = [float(r[j]) for r in reps]
        out["replicates"][k] = col
        out["interval"][k] = percentile_interval(col, B)
    return out


def t1e_pairs(var_strata, fam, u1, co, gene_node, edges, region_of):
    """The pairs of T1-e: unordered pairs of genes of one variation-bearing stratum, both with a copy number, a U1 label and a node, whose
    nodes are joined by a direct edge of the dumped graph (two genes on one node are not joined). `edges` maps a node pair, in either order, to
    its weight (the maximum if listed twice)."""
    w_of = {}
    for (a, b), w in edges.items():
        key = tuple(sorted((a, b)))
        w_of[key] = max(w, w_of.get(key, w))
    rows = []
    for sid, genes in var_strata:
        ok = [g for g in sorted(genes) if co.get(g) is not None and u1.get(g) and gene_node.get(g) is not None]
        for i, g in enumerate(ok):
            for h in ok[i + 1:]:
                na, nb = gene_node[g], gene_node[h]
                if na == nb:
                    continue
                w = w_of.get(tuple(sorted((na, nb))))
                if w is None:
                    continue
                inter, union = len(u1[g] & u1[h]), len(u1[g] | u1[h])
                rg, rh = region_of.get(g), region_of.get(h)
                rows.append(dict(stratum=sid, g=g, h=h, w=float(w), x_clu=int(fam[g] == fam[h]), x_str=int(u1[g] == u1[h]),
                                 x_jac=inter / union if union else 0.0, y_cn_raw=abs(co[g] - co[h]),
                                 y_loc=int(rg is not None and rg == rh)))
    return rows


def add_cn_rank(rows):
    """y_cn = the rank of |co_g - co_h| over all pairs of all strata (average ranks for ties)."""
    for r, k in zip(rows, average_ranks([r["y_cn_raw"] for r in rows])):
        r["y_cn"] = k
    return rows


def t1e_summary(rows, B=2000, seed=SEED_T1, strata=None):
    """The registered T1-e report: for each outcome, the coefficients of x_str and x_clu with 95% percentile intervals (x_jac in place of
    x_str as the second model), and the descriptive shares. No p-value and no class."""
    import numpy as np
    add_cn_rank(rows)
    out = dict(n_pairs=len(rows), strata=len({r["stratum"] for r in rows}))
    if not rows:
        return out
    out["share_x_str"] = sum(r["x_str"] for r in rows) / len(rows)
    out["share_x_clu"] = sum(r["x_clu"] for r in rows) / len(rows)
    a, b = np.array([r["x_str"] for r in rows], float), np.array([r["x_clu"] for r in rows], float)
    out["corr_x_str_x_clu"] = float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else None
    for outcome in ("y_cn", "y_loc"):
        for name, xs in (("x_str", ["w", "x_clu", "x_str"]), ("x_jac", ["w", "x_clu", "x_jac"])):
            bs = bootstrap_fe(rows, outcome, xs, B=B, seed=seed, strata=strata)
            fit = fe_ols(rows, outcome, xs)
            out[f"{outcome}|{name}"] = dict(coef=fit["coef"], interval=bs["interval"], rank_deficient=fit["rank_deficient"])
    return out


# ---------------------------------------------------------------------------------- part 5: registered inputs and the universe
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
HOME = "/mnt/linuxdisk/home/juanfraitu"
RUN = "/mnt/linuxdisk/tmp/hier_run/t1"
INPUTS = RUN + "/inputs"
M2 = HOME + "/rustle_m2/m2data/out"
WL = HOME + "/winloci_data"

# the F_O files are read from the hashed copy (INPUTS) and were copied from here
FO_SOURCES = {
    "names.tsv": "/mnt/linuxdisk/tmp/ourhom/names.tsv",
    "loci.gff3": "/mnt/linuxdisk/tmp/ourhom/loci.gff3",
    "loci.paf": "/mnt/linuxdisk/tmp/ourhom/loci.paf",
    "ours.clusters.tsv": "/mnt/linuxdisk/tmp/ourhom/ours.clusters.tsv",
    "ours.loci.tsv": "/mnt/linuxdisk/tmp/ourhom/ours.loci.tsv",
    "ours.pairs.tsv": "/mnt/linuxdisk/tmp/ourhom/ours.pairs.tsv",
    "ours.params.tsv": "/mnt/linuxdisk/tmp/ourhom/ours.params.tsv",
    "ours.graph.tsv": "/mnt/linuxdisk/tmp/prereg_revision/fo_graph/ours.graph.tsv",
}


def default_paths():
    return dict(
        colors=WL + "/duplicons/chm13.draft_v1.0_plus38Y_dupmasker_colors.bed",
        duplicons_bed=WL + "/duplicons/chm13.draft_v1.0_plus38Y.duplicons.bed",
        s1c=REPO + "/bench/soto/soto_famCN_S1C.tsv",
        split=REPO + "/bench/soto/soto_split_2026-09-29.tsv",
        mapback=REPO + "/bench/soto/shared_exons_5154_exon_mapback.tsv",
        sd98_exons=WL + "/soto_replication/sd98_gene_exons.tsv",
        sd98_regions=WL + "/soto_replication/sd98_v1.bed",
        famcn_ours=WL + "/soto_replication/famcn_ours_allwssd.tsv",
        gene_v2=WL + "/soto_replication/gene_v2_coords.bed",
        families_cn=M2 + "/families_cn.json",
        elig=M2 + "/elig.tsv",
        full=M2 + "/full.tsv",
        rt_genes="/mnt/linuxdisk/tmp/rustle_figures_dev/soto_losses/frozen/rt_genes.tsv",
        fo_dir=INPUTS,
        cat_genes=WL + "/gencode_chm13/chm13v2.0_CAT_Liftoff.genes.tsv",
        cat_refseq=WL + "/gencode_chm13/chm13v2.0_CAT_Liftoff.refseq_map.tsv",
        cat_v4=WL + "/soto_replication/cat_v4.bed",
        t2t_v1=WL + "/soto_replication/t2t-chm13-v1.0.fa.gz",
        mcl_binary="/mnt/linuxdisk/tmp/ourhom/bin/mcl_families",
    )


def load_s1c(path):
    """S1C by Gene ID (never by name): name, biotype, eligibility ('In Table S1 (SD98 gene set)' == Yes on any row) and the set of
    non-Unassigned Soto family labels over its rows."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            g = r["Gene ID"]
            if not g:
                continue
            o = out.setdefault(g, dict(name=r["Gene Name"], biotype=r["Biotype"], eligible=False, fams=set()))
            if r["In Table S1 (SD98 gene set)"] == "Yes":
                o["eligible"] = True
            f = r["Family ID"]
            if f and not f.startswith("Unassigned"):
                o["fams"].add(f)
    for o in out.values():
        o["fams"] = frozenset(o["fams"])
        o["cls"] = "clean" if len(o["fams"]) == 1 else ("multi" if len(o["fams"]) >= 2 else "famless")
    return out


def load_split(path):
    with open(path) as fh:
        return {r["gene_id"]: r["half"] for r in csv.DictReader(fh, delimiter="\t")}


def load_gene_exons(path, wanted):
    ex = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene_id"] in wanted:
                ex[r["gene_id"]].append((r["chrom"], int(r["start"]), int(r["end"])))
    return ex


class FO:
    """The F_O layer: the loci of the 2,334 genes, the cluster rows, the fold of annotation records into representatives and the
    dumped pre-MCL graph."""

    def __init__(self, d):
        self.names, self.loc2gene = {}, {}
        with open(d + "/names.tsv") as fh:
            for line in fh:
                g, loc = line.rstrip("\n").split("\t")
                self.names[g] = loc
                self.loc2gene[loc] = g
        self.cluster_of, self.rows = {}, Counter()
        with open(d + "/ours.clusters.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                self.cluster_of[f"{r['chrom']}:{r['start']}-{r['end']}"] = r["cluster_id"]
                self.rows[r["cluster_id"]] += 1
        self.fold = {}
        with open(d + "/ours.loci.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                self.fold[r["annotation"]] = r["representative"]
        self.edges = []
        with open(d + "/ours.graph.tsv") as fh:
            for line in fh:
                u, v, w = line.rstrip("\n").split("\t")
                self.edges.append((u, v, float(w)))
        self.nodes = {n for u, v, _ in self.edges for n in (u, v)}
        self.comp = graph_components(self.edges)
        self.weight = {}
        for u, v, w in self.edges:
            self.weight[(u, v)] = max(w, self.weight.get((u, v), w))

    def family(self, g):
        """(cluster id, rule): rule 1 the gene's own locus is a cluster row, rule 2 the representative of the locus is, else (None, None)."""
        loc = self.names.get(g)
        if loc in self.cluster_of:
            return self.cluster_of[loc], 1
        rep = self.fold.get(loc)
        if rep in self.cluster_of:
            return self.cluster_of[rep], 2
        return None, None

    def node(self, g):
        loc = self.names.get(g)
        if loc in self.nodes:
            return loc
        rep = self.fold.get(loc)
        return rep if rep in self.nodes else None


class Universe:
    """Everything the gates and the cells read, built once from the registered inputs."""


def build_universe(paths, tie="lex", union=False):
    U = Universe()
    s1c = load_s1c(paths["s1c"])
    with open(paths["families_cn"]) as fh:
        fc = json.load(fh)["genes"]
    U.genes = [x["id"] for x in fc]                                    # the order of families_cn.json (the frozen scripts' order)
    if set(U.genes) != set(s1c):
        raise ValueError("families_cn.json and S1C do not list the same genes")
    ex_rows = load_gene_exons(paths["sd98_exons"], set(s1c))
    missing = sorted(set(s1c) - set(ex_rows))
    if missing:
        raise ValueError(f"{len(missing)} S1C genes have no row in the exon table, for example {missing[:3]}")
    U.colors = Colors.load(paths["colors"])
    half = load_split(paths["split"])
    with open(paths["rt_genes"]) as fh:
        U.readthrough = {line.split("\t")[0] for line in list(fh)[1:] if line.strip()}
    U.regions = Regions.load(paths["sd98_regions"])
    U.info, U.strings = {}, {}
    U.exons_match_families_cn = True
    for x in fc:
        g = x["id"]
        rows = ex_rows[g]
        chrom = rows[0][0]
        merged = merge_intervals([[a, b] for _, a, b in rows])
        if merged != merge_intervals(x["x"]) or chrom != x["c"]:
            U.exons_match_families_cn = False
        st = gene_strings(merged, chrom, U.colors, tie=tie, union=union)
        U.strings[g] = st
        interior = frozenset(d for s in st.exon_dom[1:-1] for d in s)
        U.info[g] = dict(s1c[g], chrom=chrom, start=x["s"], end=x["e"], strand=x["st"], exons=merged, nex=len(merged), cs=x.get("cs"),
                         co=x.get("co"), half=half.get(g), interior=interior, sf=tuple(x["sf"]))
    U.V = sorted(g for g in U.genes if U.info[g]["cls"] != "multi" and U.strings[g].pooled)
    U.V_clean = [g for g in U.V if U.info[g]["cls"] == "clean"]
    U.V_famless = [g for g in U.V if U.info[g]["cls"] == "famless"]
    U.region = {g: best_region(U.regions, U.info[g]["chrom"], U.info[g]["exons"]) for g in U.genes}
    U.fo = FO(paths["fo_dir"])
    U.fo_label, U.fo_rule, U.gene_node = {}, {}, {}
    for g in U.genes:
        U.fo_label[g], U.fo_rule[g] = U.fo.family(g)
        U.gene_node[g] = U.fo.node(g)
    U.VF = [g for g in U.V if U.fo_label[g] is not None]
    U.strata = build_strata(U.VF, U.gene_node, U.fo.comp)
    return U

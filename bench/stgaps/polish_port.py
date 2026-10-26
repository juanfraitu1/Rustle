#!/usr/bin/env python3
"""Python port of `polish_drop_set` (src/copy_assign_bin.rs), fuzzy = 0, no ISM escape / 3' anchoring.

Used to measure StringTie-style variants of the isoform-fraction step (rlink.cpp:20058-20120) on an
existing `--assembly-polish none` GTF without rebuilding the binary. `polish()` returns
{tid: step} for every dropped transcript; the default arguments reproduce the shipped polish.

frac_mode:
  gene   shipped rule: r < F * max(reads) over the transcript's gene_id
  local  StringTie-like: r < F * max(reads) over kept transcripts on the same contig+strand that share
         >= 1 exonic bp with it
"""
import re
import sys
from collections import defaultdict

ATTR = {k: re.compile(r'(?:^|; )%s "([^"]*)"' % k) for k in ("transcript_id", "gene_id", "reads")}


def attr(s, k):
    m = ATTR[k].search(s)
    return m.group(1) if m else None


def load(path):
    exons, key, reads, gene = defaultdict(list), {}, {}, {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        t = attr(f[8], "transcript_id")
        if t is None:
            continue
        r = attr(f[8], "reads")
        if r is not None and r.isdigit():
            reads[t] = max(reads.get(t, 0), int(r))
        g = attr(f[8], "gene_id")
        if g is not None:  # undo bridge regrouping (it runs after the polish): .rgN / .fusN suffixes
            gene.setdefault(t, re.sub(r"\.(rg|fus)\d+$", "", g))
        if f[2] != "exon":
            continue
        exons[t].append((int(f[3]) - 1, int(f[4])))
        key.setdefault(t, (f[0], f[6]))
    chain, span = {}, {}
    for t, ex in exons.items():
        ex.sort()
        chain[t] = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
        span[t] = (ex[0][0], ex[-1][1])
    return dict(exons), key, reads, gene, chain, span


def quantile(vals, q):
    vals = sorted(vals)
    if not vals or q <= 0:
        return 0
    return vals[min(int(q * len(vals)), len(vals) - 1)]


def overlap_groups(ts, exons, key):
    """For each transcript, the set of transcripts (same contig+strand) sharing >= 1 exonic bp."""
    by = defaultdict(list)
    for t in ts:
        for a, b in exons[t]:
            by[key[t]].append((a, b, t))
    nb = defaultdict(set)
    for k, iv in by.items():
        iv.sort()
        active = []  # (end, tid)
        for a, b, t in iv:
            active = [(e, u) for e, u in active if e > a]
            for _, u in active:
                if u != t:
                    nb[t].add(u)
                    nb[u].add(t)
            active.append((b, t))
    return nb


def polish(data, mono_quantile=0.82, frac=0.02, shadow=True, ism_ratio=0.7, retained=10.0,
           frac_mode="gene", credit_ism=False):
    exons, key, reads, gene, chain, span = data
    R = lambda t: reads.get(t, 0)
    sites = {}
    drop = set()
    credit = defaultdict(float)  # reads of ISM-dropped fragments credited to their container(s)

    def mark(step):
        for t in drop:
            sites.setdefault(t, step)

    groups = defaultdict(list)
    for t, k in key.items():
        groups[k].append(t)
    for _, tids in sorted(groups.items()):
        multi = sorted((t for t in tids if chain[t]), key=lambda t: (-len(chain[t]), t))
        by_first = defaultdict(list)
        for y in multi:
            by_first[chain[y][0]].append(y)
        for x in multi:
            if x in drop:
                continue
            cx = chain[x]
            for j in cx:
                for y in by_first.get(j, ()):
                    if x == y or y in drop:
                        continue
                    cy = chain[y]
                    if len(cy) >= len(cx):
                        continue
                    sub = any(cx[k:k + len(cy)] == cy for k in range(len(cx) - len(cy) + 1))
                    rc = R(x)
                    if sub and not (rc > 0 and R(y) >= ism_ratio * rc):
                        drop.add(y)
                        credit[x] += R(y)
        alive = [m for m in multi if m not in drop]
        for t in sorted(t for t in tids if not chain[t] and t not in drop):
            s = span[t]
            for h in alive:  # first surviving multi in `multi` order whose span contains the mono
                if span[h][0] <= s[0] and span[h][1] >= s[1]:
                    rc = R(h)
                    if not (rc > 0 and R(t) >= ism_ratio * rc):
                        drop.add(t)
                    break
    for t in drop:
        sites.setdefault(t, "ism" if chain[t] else "ism_mono")

    floor = quantile([R(t) for t, c in chain.items() if c and t not in drop], mono_quantile)
    if floor > 0:
        for t, c in chain.items():
            if not c and t not in drop and R(t) < floor:
                drop.add(t)
    mark("floor")
    if shadow:
        ex_any, sp_same = defaultdict(list), defaultdict(list)
        for t, c in chain.items():
            if c and t not in drop:
                sp_same[key[t]].append(span[t])
                ex_any[key[t][0]].extend(exons[t])

        def merge(v):
            out = []
            for a, b in sorted(v):
                if out and a <= out[-1][1]:
                    out[-1][1] = max(out[-1][1], b)
                else:
                    out.append([a, b])
            return out

        ex_any = {k: merge(v) for k, v in ex_any.items()}
        sp_same = {k: merge(v) for k, v in sp_same.items()}
        import bisect

        def hits(v, s, e):
            hi = bisect.bisect_right([a for a, _ in v], e)
            return hi > 0 and v[hi - 1][1] >= s

        for t in sorted(t for t, c in chain.items() if not c and t not in drop):
            s, e = span[t]
            if (key[t][0] in ex_any and hits(ex_any[key[t][0]], s, e)) or \
               (key[t] in sp_same and hits(sp_same[key[t]], s, e)):
                drop.add(t)
    mark("shadow")

    if frac > 0:
        S = (lambda t: R(t) + credit[t]) if credit_ism else R
        kept = [t for t in gene if t not in drop]
        if frac_mode == "gene":
            best = defaultdict(float)
            for t in kept:
                best[gene[t]] = max(best[gene[t]], S(t))
            ref = {t: best[gene[t]] for t in kept}
        else:
            nb = overlap_groups([t for t in kept if t in exons], exons, key)
            ref = {t: max([S(t)] + [S(u) for u in nb.get(t, ())]) for t in kept}
        for t in sorted(kept):
            b = ref[t]
            if b > 0 and S(t) < b and S(t) < frac * b:
                drop.add(t)
    mark("frac")

    if retained > 0:
        support, gj = defaultdict(int), defaultdict(set)
        for t, c in chain.items():
            if t in drop:
                continue
            for j in c:
                support[(key[t], j)] += R(t)
            if t in gene:
                gj[(gene[t], key[t][1])].update(c)
        newly = set()
        for t in sorted(t for t in chain if t not in drop):
            if t not in gene:
                continue
            own = set(chain[t])
            for j in sorted(gj.get((gene[t], key[t][1]), ())):
                if j in own:
                    continue
                if any(a < j[0] and j[1] < b for a, b in exons[t]) and support[(key[t], j)] >= retained * R(t):
                    newly.add(t)
                    break
        drop |= newly
    mark("ret")
    return sites


if __name__ == "__main__":
    data = load(sys.argv[1])
    s = polish(data)
    from collections import Counter
    c = Counter(s.values())
    print(len(data[0]), "transcripts;", dict(c), "kept", len(data[0]) - len(s))

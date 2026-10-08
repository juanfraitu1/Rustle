#!/usr/bin/env python3
"""Attribution of a cluster consensus (or a single read) to a family (docs/PREREG_unmapped_rescue_2026-10-08.md section 1, step 3)."""
import collections
import subprocess


def family_of(target):
    """GWFAM12:3 -> GWFAM12"""
    return target.split(":")[0]


def attribute(rows, max_evalue=1e-5):
    """rows: [(query, target, bit score, e-value)] of a translated search. -> {query: (family | None, best bits, runner-up family's best bits)}:
    the family of the best-scoring target iff its bit score is STRICTLY above every other family's, else None (abstain: a tie).
    Hits above max_evalue are ignored; a query without a hit is absent."""
    best = collections.defaultdict(dict)
    for q, t, bits, ev in rows:
        if ev > max_evalue:
            continue
        f = family_of(t)
        if bits > best[q].get(f, -1.0):
            best[q][f] = bits
    out = {}
    for q, fam in best.items():
        ranked = sorted(fam.items(), key=lambda kv: -kv[1])
        top, second = ranked[0], (ranked[1][1] if len(ranked) > 1 else 0.0)
        out[q] = (top[0] if top[1] > second else None, top[1], second)
    return out


def read_m8(path):
    """mmseqs convertalis --format-output query,target,bits,evalue"""
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        yield (f[0], f[1], float(f[2]), float(f[3]))


def mmseqs_translated(query_fa, target_fa, out_m8, tmp, threads=4, evalue=1e-3, sens=7.5):
    """mmseqs2 translated-vs-translated search (tblastx-like, --search-type 2); the e-value cut of the registered rule is applied in attribute()"""
    cmd = (f"mmseqs easy-search {query_fa} {target_fa} {out_m8} {tmp} --search-type 2 -e {evalue} -s {sens} --threads {threads} "
           f"--format-output query,target,bits,evalue -v 1")
    subprocess.run(cmd, shell=True, check=True)


# ---- the frozen rule (prereg Amendment 1): nucleotide cover of the whole consensus, with a relative margin ----

def union_length(intervals):
    """number of bases covered by the union of inclusive (start, end) intervals"""
    tot, cur = 0, None
    for a, b in sorted(intervals):
        if cur is None or a > cur[1] + 0:
            if cur is not None:
                tot += cur[1] - cur[0] + 1
            cur = [a, b]
        else:
            cur[1] = max(cur[1], b)
    return tot + (cur[1] - cur[0] + 1 if cur else 0)


def cover_scores(hsps):
    """hsps: [(query, target, qstart, qend)] (1-based, inclusive, either order). -> {query: {family: covered query bases of the family's best target}}"""
    by = collections.defaultdict(list)
    for q, t, a, b in hsps:
        by[(q, t)].append((min(a, b), max(a, b)))
    out = collections.defaultdict(dict)
    for (q, t), iv in by.items():
        f = family_of(t)
        out[q][f] = max(out[q].get(f, 0), union_length(iv))
    return dict(out)


def attribute_cover(scores, margin=1.10):
    """{query: {family: score}} -> {query: (family | None, top, runner-up)}: the best family iff top >= margin x runner-up and top > 0, else None"""
    out = {}
    for q, fam in scores.items():
        ranked = sorted(fam.items(), key=lambda kv: -kv[1])
        if not ranked:
            out[q] = (None, 0, 0)
            continue
        top, second = ranked[0][1], (ranked[1][1] if len(ranked) > 1 else 0)
        out[q] = (ranked[0][0] if top > 0 and top + 1e-9 >= margin * second and top > second else None, top, second)
    return out


def read_blastn(path):
    """blastn -outfmt '6 qseqid sseqid bitscore evalue length pident qstart qend' -> [(query, target, qstart, qend)]; the consensus header '|n=' suffix is dropped"""
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        yield (f[0].split("|")[0], f[1], int(f[6]), int(f[7]))


def specific_cover_scores(hsps):
    """Amendment 2 (exploratory): repeat-robust family score. hsps: [(query, target, qstart, qend)] (1-based inclusive, either order).
    A family covers a consensus base if ANY of its copies has an HSP over it; a base covered by n families is credited 1/n to each of them, so bases
    that many families share (repeats) count little and bases only one family covers count in full. -> {query: {family: score}}"""
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for q, t, a, b in hsps:
        by[q][family_of(t)].append((min(a, b), max(a, b)))
    out = {}
    for q, fams in by.items():
        merged = {}
        for f, iv in fams.items():
            m = []
            for a, b in sorted(iv):
                if m and a <= m[-1][1]:
                    m[-1][1] = max(m[-1][1], b)
                else:
                    m.append([a, b])
            merged[f] = m
        events = sorted({p for m in merged.values() for a, b in m for p in (a, b + 1)})
        score = {f: 0.0 for f in merged}
        for lo, hi in zip(events, events[1:]):
            cov = [f for f, m in merged.items() if any(a <= lo and hi - 1 <= b for a, b in m)]
            for f in cov:
                score[f] += (hi - lo) / len(cov)
        out[q] = score
    return out


def capped_queries(hsps, cap):
    """queries whose BLAST output reached `cap` distinct targets: -max_target_seqs truncates in an order that is not by score, so a query at the cap
    can lose its own family's copy. Review of 2026-10-08 found cap 500 truncating 87 of 188 dev consensus sequences."""
    seen = collections.defaultdict(set)
    for q, t, _a, _b in hsps:
        seen[q].add(t)
    return sorted(q for q, ts in seen.items() if len(ts) >= cap)

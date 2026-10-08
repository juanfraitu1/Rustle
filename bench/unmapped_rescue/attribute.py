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

#!/usr/bin/env python3
"""Amendment 45: IsoCon-style read correction. A read is corrected by majority vote of its neighbours' alignments to it (long cs strings, so secondary
alignments count too). Pure: pileup_cs, majority, correct. minimap2-backed: correct_set."""
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polish as P  # noqa: E402

TOK = re.compile(r"(=[ACGTN]+|\*[acgtn][acgtn]|\+[acgtn]+|-[acgtn]+|~[acgtn]{2}\d+[acgtn]{2})")
IDX = {b: i for i, b in enumerate("ACGT")}


def pileup_cs(L, recs):
    """recs = [(target start, long cs)] of neighbours aligned to a target of length L -> (cnt, ins, ngap) as polish.pileup"""
    cnt = [[0] * 5 for _ in range(L)]
    ins = [dict() for _ in range(L + 1)]
    span = [0] * (L + 2)
    for ts, cs in recs:
        r = ts
        for t in TOK.findall(cs):
            c = t[0]
            if c == "=":
                for b in t[1:]:
                    if r < L and b in IDX:
                        cnt[r][IDX[b]] += 1
                    r += 1
            elif c == "*":
                b = t[2].upper()
                if r < L and b in IDX:
                    cnt[r][IDX[b]] += 1
                r += 1
            elif c == "+":
                s = t[1:].upper()
                if r <= L:
                    ins[r][s] = ins[r].get(s, 0) + 1
            elif c == "-":
                for _ in range(len(t) - 1):
                    if r < L:
                        cnt[r][4] += 1
                    r += 1
            else:
                r += int(re.findall(r"\d+", t)[0])
        span[min(L, ts + 1)] += 1
        span[max(min(L, r), ts + 1)] -= 1
    ngap, run = [0] * (L + 1), 0
    for g in range(L + 1):
        run += span[g]
        ngap[g] = run
    return cnt, ins, ngap


def majority(read, pile, min_cov=2):
    """corrections ('sub', i, base) | ('del', i, None) | ('ins', gap, string): another allele held by strictly more than half of >= min_cov covering neighbours"""
    cnt, ins, ngap = pile
    out = []
    for i, c in enumerate(cnt):
        n = sum(c)
        if n < min_cov:
            continue
        own = IDX.get(read[i])
        for a in range(5):
            if a != own and 2 * c[a] > n:
                out.append(("del", i, None) if a == 4 else ("sub", i, "ACGT"[a]))
                break
    for g in range(1, len(read)):
        n = ngap[g]
        if n < min_cov or not ins[g]:
            continue
        s, k = max(ins[g].items(), key=lambda kv: (kv[1], kv[0]))
        if 2 * k > n:
            out.append(("ins", g, s))
    return out


def correct(read, recs, min_cov=2):
    return P.apply(read, majority(read, pileup_cs(len(read), recs), min_cov))


def correct_set(reads, work, threads=4):
    """every read of the set corrected by the others (one all-vs-all, map-hifi, -N 50 -p 0.1, long cs)"""
    os.makedirs(work, exist_ok=True)
    with open(f"{work}/s.fa", "w") as o:
        for n, s in reads.items():
            o.write(f">{n}\n{s}\n")
    lines = subprocess.run(f"minimap2 -c --cs=long -x map-hifi -N 50 -p 0.1 -t {threads} {work}/s.fa {work}/s.fa", shell=True, stdout=subprocess.PIPE,
                           stderr=subprocess.DEVNULL, text=True).stdout.splitlines()
    by = {}
    for ln in lines:
        f = ln.split("\t")
        if f[0] == f[5]:
            continue
        cs = next((t[5:] for t in f[12:] if t.startswith("cs:Z:")), None)
        if cs:
            by.setdefault(f[5], []).append((int(f[7]), cs.rstrip("\n")))
    return {n: correct(s, by.get(n, [])) for n, s in reads.items()}

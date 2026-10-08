#!/usr/bin/env python3
"""IsoCon-style correction of a cluster consensus (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 9). Pure functions; pipeline in run_polish.py.

pileup(sam_lines, cons) -> (cnt, ins, ngap): cnt[i] = [A, C, G, T, deleted] counts of the primary read alignments at consensus position i; ins[g] = {string: reads}
inserted at gap g (between positions g-1 and g); ngap[g] = reads spanning that gap."""
import math
import re

BASES = "ACGT"
IDX = {b: i for i, b in enumerate(BASES)}
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
CS = re.compile(r":[0-9]+|\*[a-z][a-z]|[+\-][A-Za-z]+|~[a-z]{2}[0-9]+[a-z]{2}")
ALPHA = 0.05
MAX_ROUNDS = 5


def binom_tail(k, n, p):
    """P(Binomial(n, p) >= k)"""
    if k <= 0:
        return 1.0
    if k > n or p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    lp, lq = math.log(p), math.log1p(-p)
    terms = [math.lgamma(n + 1) - math.lgamma(j + 1) - math.lgamma(n - j + 1) + j * lp + (n - j) * lq for j in range(k, n + 1)]
    m = max(terms)
    return min(1.0, math.exp(m) * math.fsum(math.exp(t - m) for t in terms))


def pileup(sam_lines, cons):
    L = len(cons)
    cnt = [[0] * 5 for _ in range(L)]
    ins = [dict() for _ in range(L + 1)]
    span = [0] * (L + 2)
    for ln in sam_lines:
        if ln[0] == "@":
            continue
        f = ln.rstrip("\n").split("\t")
        if int(f[1]) & 2308 or f[2] == "*":
            continue
        r, q, seq = int(f[3]) - 1, 0, f[9]
        start = r
        for n, op in CIG.findall(f[5]):
            n = int(n)
            if op in "=XM":
                for _ in range(n):
                    b = cons[r] if op == "=" else seq[q]
                    if op == "M":
                        b = seq[q]
                    if b in IDX:
                        cnt[r][IDX[b]] += 1
                    r += 1
                    q += 1
            elif op == "I":
                s = seq[q:q + n]
                ins[r][s] = ins[r].get(s, 0) + 1
                q += n
            elif op == "D":
                for _ in range(n):
                    cnt[r][4] += 1
                    r += 1
            elif op == "N":
                r += n
            elif op == "S":
                q += n
        span[start + 1] += 1      # gaps start+1 .. r-1 are spanned by this read
        span[max(r, start + 1)] -= 1
    ngap, run = [0] * (L + 1), 0
    for g in range(L + 1):
        run += span[g]
        ngap[g] = run
    return cnt, ins, ngap


def corrections(cons, pile, alpha=ALPHA):
    """-> (applied, minority, e). e = non-consensus events (substitutions, deletions, insertion events) / aligned columns. An alternative allele held by k of n
    reads is significant iff P(Binomial(n, e) >= k) < alpha / (3 L); a significant allele held by more reads than the consensus allele is applied
    ('sub', i, base) | ('del', i, None) | ('ins', gap, string); a significant allele held by fewer is recorded as a minority variant (kind, pos, base, k, n)."""
    cnt, ins, ngap = pile
    L = len(cons)
    obs = alt = 0
    for i in range(L):
        n = sum(cnt[i])
        obs += n
        alt += n - (cnt[i][IDX[cons[i]]] if cons[i] in IDX else 0)
    for g in range(1, L):
        alt += sum(ins[g].values())
    e = alt / obs if obs else 0.0
    thr = alpha / (3 * L)
    applied, minority = [], []
    for i in range(L):
        c, n = cnt[i], sum(cnt[i])
        if n == 0:
            continue
        ci = IDX.get(cons[i])
        ck = c[ci] if ci is not None else 0
        best = None
        for a in range(5):
            k = c[a]
            if a == ci or k == 0 or binom_tail(k, n, e) >= thr:
                continue
            kind = ("del", i, None) if a == 4 else ("sub", i, BASES[a])
            if k > ck:
                if best is None or k > best[0]:
                    best = (k, kind)
            else:
                minority.append(kind + (k, n))
        if best:
            applied.append(best[1])
    for g in range(1, L):
        n = ngap[g]
        if n == 0 or not ins[g]:
            continue
        none = max(0, n - sum(ins[g].values()))
        best = None
        for s, k in ins[g].items():
            if binom_tail(k, n, e) >= thr:
                continue
            if k > none:
                if best is None or k > best[0]:
                    best = (k, ("ins", g, s))
            else:
                minority.append(("ins", g, s, k, n))
        if best:
            applied.append(best[1])
    return applied, minority, e


def apply(cons, applied):
    subs = {i: b for kind, i, b in applied if kind == "sub"}
    dels = {i for kind, i, b in applied if kind == "del"}
    inss = {g: s for kind, g, s in applied if kind == "ins"}
    out = []
    for i, b in enumerate(cons):
        if i in inss:
            out.append(inss[i])
        if i in dels:
            continue
        out.append(subs.get(i, b))
    return "".join(out)


def polish(cons, align_fn, max_rounds=MAX_ROUNDS):
    """-> (consensus, rounds run, all corrections applied, minority variants of the last round). align_fn(consensus) -> SAM lines of the cluster's reads"""
    total, rounds, minority = [], 0, []
    for _ in range(max_rounds):
        rounds += 1
        applied, minority, _e = corrections(cons, pileup(align_fn(cons), cons))
        if not applied:
            break
        cons = apply(cons, applied)
        total += applied
    return cons, rounds, total, minority


def cs_edit_columns(cs, qstart, qend, qlen, strand):
    """consensus columns (0-based, forward consensus coordinates) of the edits in a minimap2 short `cs` string: mismatches, inserted bases, and the column
    after a deletion. On the reverse strand the alignment runs from the end of the aligned query interval."""
    q, rel = 0, set()
    for tok in CS.findall(cs):
        t = tok[0]
        if t == ":":
            q += int(tok[1:])
        elif t == "*":
            rel.add(q)
            q += 1
        elif t == "+":
            rel.update(range(q, q + len(tok) - 1))
            q += len(tok) - 1
        elif t == "-":
            rel.add(q)
    return {qstart + x for x in rel} if strand == "+" else {qend - 1 - x for x in rel}


def allele_like_share(edit_cols, minority, window=1):
    """share of the edit columns that lie within `window` of a significant minority-variant position; None if there are no edits"""
    if not edit_cols:
        return None
    pos = [m[1] for m in minority]
    return sum(any(abs(c - p) <= window for p in pos) for c in edit_cols) / len(edit_cols)

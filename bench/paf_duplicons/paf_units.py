#!/usr/bin/env python3
"""Our own duplicon units from the locus self-alignment (docs/PREREG_paf_duplicon_units_2026-10-10.md).

    paf_units.py <loci.paf> <loci.gff3> <out.bed>      writes the unit BED (chrom start end PU<n>) and <out.bed>.stats.json

Alignments between different loci, not the same DNA, identity >= 0.90 and block >= 1 kb (the SD definition), each side carrying an exonic base of its locus;
breakpoints = alignment ends, projected through the alignments' CIGARs until none is added (snapped within 1 kb, at most 20 rounds); atoms = covered
stretches cut at the breakpoints, >= 1 kb; units = components of atoms joined by reciprocal-majority alignment. Pure functions are tested in test_paf_units.py."""
import bisect
import collections
import json
import re
import sys

MIN_IDENT = 0.90
MIN_BLOCK = 1000
TOL = 1000
MIN_ATOM = 1000
MAX_ROUNDS = 20


def parse_region(name):
    """'chr1:112072-113065' (1-based inclusive) -> ('chr1', 112071, 113065), 0-based half-open"""
    c, rng = name.rsplit(":", 1)
    s, e = rng.split("-")
    return c, int(s) - 1, int(e)


def cigar_ops(cg):
    return [(int(n), "M" if op in "=X" else op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cg)]


class Aln:
    """One alignment in genomic coordinates; the query runs reversed against the CIGAR on the minus strand (PAF convention)."""

    def __init__(self, qc, qs, qe, tc, ts, te, strand, ops, ident, block, qn=None, tn=None):
        self.qc, self.qs, self.qe, self.tc, self.ts, self.te = qc, qs, qe, tc, ts, te
        self.strand, self.ops, self.ident, self.block = strand, ops, ident, block
        self.qn, self.tn = qn or qc, tn or tc
        qo = to = 0
        self._q, self._t = ([], [], [], []), ([], [], [], [])     # (own starts, other starts, op, length) of the ops consuming each side
        for n, op in ops:
            if op in "MI":
                for lst, v in zip(self._q, (qo, to, op, n)):
                    lst.append(v)
            if op in "MD":
                for lst, v in zip(self._t, (to, qo, op, n)):
                    lst.append(v)
            qo += n if op in "MI" else 0
            to += n if op in "MD" else 0

    @staticmethod
    def _map(idx, x, gap):
        own, other, kind, ln = idx
        i = max(0, bisect.bisect_right(own, x) - 1)
        return other[i] if kind[i] == gap else other[i] + min(max(x - own[i], 0), ln[i])


def aln_from_fields(f):
    qc, q0, _ = parse_region(f[0])
    tc, t0, _ = parse_region(f[5])
    cg = next(t[5:] for t in f[12:] if t.startswith("cg:Z:"))
    return Aln(qc, q0 + int(f[2]), q0 + int(f[3]), tc, t0 + int(f[7]), t0 + int(f[8]), f[4], cigar_ops(cg), int(f[9]) / int(f[10]), int(f[10]), f[0], f[5])


def project(a, side, pos):
    """a breakpoint (boundary before base `pos`) on one side -> (chrom, position) on the other"""
    if side == "q":
        x = pos - a.qs if a.strand == "+" else a.qe - pos
        return a.tc, a.ts + Aln._map(a._q, x, "I")
    x = Aln._map(a._t, pos - a.ts, "D")
    return a.qc, (a.qs + x) if a.strand == "+" else (a.qe - x)


def _has_exon(iv, s, e):
    return any(es < e and s < ee for es, ee in iv)


def keep(a, exons):
    if a.qn == a.tn or (a.qc == a.tc and a.qs < a.te and a.ts < a.qe):
        return False
    if a.ident < MIN_IDENT or a.block < MIN_BLOCK:
        return False
    return _has_exon(exons.get(a.qn, ()), a.qs, a.qe) and _has_exon(exons.get(a.tn, ()), a.ts, a.te)


def add_bp(bps, chrom, pos, tol):
    lst = bps.setdefault(chrom, [])
    i = bisect.bisect_left(lst, pos)
    if (i > 0 and pos - lst[i - 1] < tol) or (i < len(lst) and lst[i] - pos < tol):
        return False
    lst.insert(i, pos)
    return True


def _sides(a):
    return (("q", a.qc, a.qs, a.qe), ("t", a.tc, a.ts, a.te))


def initial_breakpoints(alns, tol):
    bps = {}
    for c, p in sorted({(c, p) for a in alns for _, c, s, e in _sides(a) for p in (s, e)}):
        add_bp(bps, c, p, tol)
    return bps


def propagate(alns, bps, tol, max_rounds):
    """-> (rounds run, fixpoint reached)"""
    for r in range(1, max_rounds + 1):
        added = 0
        for a in alns:
            for side, c, s, e in _sides(a):
                lst = bps.get(c, [])
                for p in lst[bisect.bisect_right(lst, s):bisect.bisect_left(lst, e)]:
                    added += add_bp(bps, *project(a, side, p), tol)
        if not added:
            return r, True
    return max_rounds, False


def covered(alns):
    iv = collections.defaultdict(list)
    for a in alns:
        for _, c, s, e in _sides(a):
            iv[c].append((s, e))
    out = {}
    for c, v in iv.items():
        v.sort()
        m = [list(v[0])]
        for s, e in v[1:]:
            if s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        out[c] = [tuple(x) for x in m]
    return out


def atoms(cov, bps, min_len):
    out = []
    for c in sorted(cov):
        lst = bps.get(c, [])
        for s, e in cov[c]:
            cuts = [s] + lst[bisect.bisect_right(lst, s):bisect.bisect_left(lst, e)] + [e]
            out += [(c, x, y) for x, y in zip(cuts, cuts[1:]) if y - x >= min_len]
    return out


def _find(par, i):
    while par[i] != i:
        par[i] = par[par[i]]
        i = par[i]
    return i


def join_units(alns, at):
    """-> per atom (input order) the id of its unit (a union-find root); atoms joined by reciprocal-majority alignment"""
    by = collections.defaultdict(list)
    for i, (c, s, e) in enumerate(at):
        by[c].append((s, e, i))
    ends = {}
    for c in by:
        by[c].sort()
        ends[c] = [e for _, e, _ in by[c]]

    def over(c, lo, hi):
        lst = by.get(c, [])
        j = bisect.bisect_right(ends.get(c, []), lo)
        while j < len(lst) and lst[j][0] < hi:
            yield lst[j]
            j += 1

    par = list(range(len(at)))
    for a in alns:
        for s, e, i in list(over(a.qc, a.qs, a.qe)):
            cs, ce = max(s, a.qs), min(e, a.qe)
            if ce <= cs:
                continue
            p1, p2 = project(a, "q", cs)[1], project(a, "q", ce)[1]
            lo, hi = min(p1, p2), max(p1, p2)
            for bs, be, j in over(a.tc, lo, hi):
                ts_, te_ = max(lo, bs), min(hi, be)
                if te_ <= ts_ or 2 * (te_ - ts_) < be - bs:
                    continue
                b1, b2 = project(a, "t", ts_)[1], project(a, "t", te_)[1]
                if 2 * max(0, min(max(b1, b2), e) - max(min(b1, b2), s)) >= e - s:
                    ri, rj = _find(par, i), _find(par, j)
                    if ri != rj:
                        par[max(ri, rj)] = min(ri, rj)
    return [_find(par, i) for i in range(len(at))]


def bed_lines(at, units):
    name, out = {}, []
    for (c, s, e), u in sorted(zip(at, units)):
        name.setdefault(u, f"PU{len(name) + 1}")
        out.append(f"{c}\t{s}\t{e}\t{name[u]}")
    return out


def load_exons(gff):
    """locus name 'chrom:start-end' (the gene span) -> its exons, genomic 0-based half-open"""
    region, exons = {}, collections.defaultdict(list)
    rows = [ln.rstrip("\n").split("\t") for ln in open(gff) if not ln.startswith("#") and ln.strip()]
    for f in rows:
        if f[2] == "gene":
            gid = re.search(r"ID=([^;]+)", f[8]).group(1)
            region[gid] = f"{f[0]}:{f[3]}-{f[4]}"
    for f in rows:
        if f[2] == "exon":
            m = re.search(r"Parent=([^;]+)", f[8])
            if m and m.group(1) in region:
                exons[region[m.group(1)]].append((int(f[3]) - 1, int(f[4])))
    return exons


def main(paf, gff, out):
    exons = load_exons(gff)
    n, alns = 0, []
    for ln in open(paf):
        n += 1
        a = aln_from_fields(ln.rstrip("\n").split("\t"))
        if keep(a, exons):
            alns.append(a)
    bps = initial_breakpoints(alns, TOL)
    rounds, fix = propagate(alns, bps, TOL, MAX_ROUNDS)
    at = atoms(covered(alns), bps, MIN_ATOM)
    units = join_units(alns, at)
    lines = bed_lines(at, units)
    open(out, "w").write("".join(x + "\n" for x in lines))
    copies = collections.Counter(collections.Counter(x.split("\t")[3] for x in lines).values())
    st = dict(paf_records=n, kept=len(alns), rounds=rounds, fixpoint=fix, breakpoints=sum(len(v) for v in bps.values()), atoms=len(at),
              atom_bp=sum(e - s for _, s, e in at), units=len(set(units)), units_by_copies=dict(sorted(copies.items())))
    json.dump(st, open(out + ".stats.json", "w"), indent=1)
    print(json.dumps(st))


if __name__ == "__main__":
    main(*sys.argv[1:4])

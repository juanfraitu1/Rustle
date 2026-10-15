#!/usr/bin/env python3
"""DEVELOPMENT prototype (seeds 11-12 only) of the structure projection, before the rule is frozen in the pre-registration.

    proto.py WORLD VARIANT        # project every flagged union's exon structure, score it against the truth, print per family

Projection (route G, the genome): the union is aligned to the reference genome with the stage's spliced preset; every primary or
supplementary record contributes exon blocks in union coordinates (an `N` CIGAR op = an intron = a boundary); union stretches no record
covers, and insertions >= 20 bp, are segments the reference lacks (their edges are boundaries). Truth: the true copy's exon lengths along
its transcript; the union is mapped onto the true transcript and each projected boundary is carried through that alignment.
"""
import collections
import csv
import json
import re
import subprocess
import sys

MIN_SEG = 20
TOL = (0, 2, 5)


def rc(x):
    return x.translate(str.maketrans("ACGTacgtN", "TGCAtgcaN"))[::-1]


def read_fa(path):
    out, name = {}, None
    for ln in open(path):
        if ln.startswith(">"):
            name = ln[1:].split()[0]; out[name] = []
        else:
            out[name].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


def mm2(args, target, query):
    return subprocess.run(f"minimap2 {args} {target} {query} 2>/dev/null", shell=True, check=True, capture_output=True, text=True,
                          executable="/bin/bash").stdout.splitlines()


def paf(lines):
    for ln in lines:
        f = ln.split("\t")
        tags = dict((t[:2], t[5:]) for t in f[12:])
        yield dict(q=f[0], ql=int(f[1]), qs=int(f[2]), qe=int(f[3]), st=f[4], t=f[5], ts=int(f[7]), te=int(f[8]), m=int(f[9]), bl=int(f[10]),
                   tp=tags.get("tp", "P"), cg=tags.get("cg", ""))


def cigar(cg):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cg)]


SEG = 50        # an insertion or uncovered stretch of at least this length is re-aligned on its own (segment pass)
MATCH, MISMATCH = 1, -3


def walk(r):
    """one + strand record -> list of ops in query order: ('M', qs, qe, ts, te, longest '=' run) | ('I', qs, qe, t) |
    ('N', q, ts, te) | ('D', q, ts, te)"""
    assert r["st"] == "+", "prototype: + strand records only (every sim union is +)"
    out, q, t, run = [], r["qs"], r["ts"], 0
    for n, op in cigar(r["cg"]):
        if op in "M=X":
            run = run + n if op != "X" else 0
            if out and out[-1][0] == "M":
                out[-1] = ("M", out[-1][1], q + n, out[-1][3], t + n, max(out[-1][5], run))
            else:
                out.append(("M", q, q + n, t, t + n, run))
            q += n; t += n
        else:
            run = 0
            if op == "I":
                out.append(("I", q, q + n, t)); q += n
            elif op == "N":
                out.append(("N", q, t, t + n)); t += n
            elif op == "D":
                out.append(("D", q, t, t + n)); t += n
    return out


def extend(a, b):
    """ungapped extension score (match +1, mismatch -3): the length k <= min(len) maximising the score of a[:k] vs b[:k]"""
    best, k, sc = 0, 0, 0
    for i in range(min(len(a), len(b))):
        sc += MATCH if a[i] == b[i] else MISMATCH
        if sc > best:
            best, k = sc, i + 1
    return k


ANCHOR = 8       # a matched run shorter than this between gaps is chance agreement inside an insertion, not an anchor
EDGE_TOL = 5     # a hop between copies is kept when its far side lands this close to an annotated exon edge of a family copy
SLACK = 2000     # an end block may lie this far outside its home copy's span before it counts as a hop onto another copy


def outer_boundaries(r, U, G, home, edges=(set(), set())):
    """boundaries of one outer record and its residual inserted segments >= SEG (for the segment pass).
    Solid blocks = matched runs >= ANCHOR. Between two consecutive solid blocks, a gap holding an intron is a junction: the query gap
    (inserted bases, however the aligner shredded them) is normalised against the genome gap by head / tail extension; what is left
    is a novel segment (two boundaries) or nothing (one boundary). A gap without an intron is an indel inside an exon. End blocks that
    lie outside the home copy (+- SLACK) or leave < SEG union bases beyond the junction are dropped (overhang)."""
    ops = walk(r)
    solid = [o for o in ops if o[0] == "M" and o[5] >= ANCHOR]
    if not solid:
        return set(), [], (r["qs"], r["qs"])
    has_n = []                               # has_n[i]: an intron lies between solid[i] and solid[i + 1]
    for a, b in zip(solid, solid[1:]):
        has_n.append(any(o[0] == "N" and a[2] <= o[1] <= b[1] and a[4] <= o[2] < b[3] for o in ops))
    # exons = runs of solid blocks without an intron between them: [first, last] index pairs
    exons, i0 = [], 0
    for i, n in enumerate(has_n):
        if n:
            exons.append((i0, i)); i0 = i + 1
    exons.append((i0, len(solid) - 1))
    inside = lambda e: home is None or (solid[e[1]][4] > home[0] - SLACK and solid[e[0]][3] < home[1] + SLACK)
    span = lambda e: solid[e[1]][2] - solid[e[0]][1]
    near = lambda x, pool: any(abs(x - y) <= EDGE_TOL for y in pool)
    # an end exon off the home copy is a hop: kept only when it lands on an annotated exon edge (its far side: the leading exon's
    # genomic end is an exon end, the trailing exon's genomic start an exon start)
    lead_ok = lambda e: inside(e) or near(solid[e[1]][4], edges[1])
    trail_ok = lambda e: inside(e) or near(solid[e[0]][3], edges[0])
    while len(exons) > 1 and (not lead_ok(exons[0]) or span(exons[0]) < SEG):
        exons.pop(0)
    while len(exons) > 1 and (not trail_ok(exons[-1]) or span(exons[-1]) < SEG):
        exons.pop()
    cuts, segs = set(), []
    for x, y in zip(exons, exons[1:]):
        a, b = solid[x[1]], solid[y[0]]
        qL, qR, tL, tR = a[2], b[1], a[4], b[3]
        ins, gap = U[qL:qR], G[tL:tR]
        k1 = extend(ins, gap)
        k2 = extend(ins[k1:][::-1], gap[::-1])
        qs2, qe2 = qL + k1, qR - k2
        if qe2 - qs2 >= MIN_SEG:
            cuts.update((qs2, qe2))
            if qe2 - qs2 >= SEG:
                segs.append((qs2, qe2))
        else:
            cuts.add(qs2)
    return cuts, segs, (solid[exons[0][0]][1], solid[exons[-1][1]][2])


def project(hits, U, G, realign, copies=()):
    """union -> (sorted internal boundaries, note). hits: the union's genome records; G: the best hit's contig sequence;
    realign(list of (qs, qe)) -> {(qs, qe): [records]} aligns segments of the union on their own (same preset)."""
    qlen = len(U)
    if not hits:
        return None, "no_hit"
    prim = [h for h in hits if h["tp"] == "P"] or hits    # the structure follows minimap2's primary placement
    best = max(prim, key=lambda h: (h["m"] / h["bl"]) * (h["qe"] - h["qs"]) / h["ql"] if h["bl"] else 0)
    if best["st"] == "-":
        # work on the reverse complement of the union, where the primary is a + record; boundaries are mapped back at the end
        flip = lambda h: dict(h, qs=qlen - h["qe"], qe=qlen - h["qs"], st="+" if h["st"] == "-" else "-")
        cuts, note = project([flip(h) for h in hits], rc(U), G, realign, copies)
        return (sorted(qlen - c for c in cuts) if cuts is not None else None), note + " (- strand)"
    keep = sorted((h for h in hits if h["tp"] == "P" and h["t"] == best["t"] and h["st"] == "+"), key=lambda h: h["qs"]) or [best]
    ov = lambda c: max(0, min(c[2], best["te"]) - max(c[1], best["ts"]))
    home = max((c for c in copies if c[0] == best["t"]), key=ov, default=None)
    home = (home[1], home[2]) if home and ov(home) > 0 else None
    edges = ({x for c in copies if c[0] == best["t"] for x, _ in c[3]}, {y for c in copies if c[0] == best["t"] for _, y in c[3]})
    cuts, segs = set(), []
    span = []
    for h in keep:
        c, s, sp = outer_boundaries(h, U, G, home, edges)
        cuts |= c; segs += s; span.append(sp)
    # between records: a breakpoint (records abut) or an uncovered stretch; ends: terminal overhangs
    for a, b in zip(keep, keep[1:]):
        cuts.update((a["qe"], b["qs"]) if b["qs"] - a["qe"] >= MIN_SEG else (b["qs"],))
        if b["qs"] - a["qe"] >= SEG:
            segs.append((a["qe"], b["qs"]))
    first, last = span[0][0], span[-1][1]          # the union stretch the kept exons anchor; the rest is terminal overhang
    term = [(0, first)] if first >= SEG else []
    term += [(last, qlen)] if qlen - last >= SEG else []
    res = realign(segs + term, U) if segs or term else {}
    for (qs, qe), recs in res.items():
        recs = sorted((r for r in recs if r["tp"] == "P" and r["t"] == best["t"] and r["st"] == "+"), key=lambda r: r["qs"])
        if (qs, qe) in term and home is not None:      # a terminal overhang counts only where it lands on its home copy (no hop)
            recs = [r for r in recs if r["te"] > home[0] - SLACK and r["ts"] < home[1] + SLACK]
        if not recs:
            continue
        for r in recs:                                  # internal introns of the segment's own alignment
            for o in walk(r):
                if o[0] == "N":
                    cuts.add(qs + o[1])
        for a, b in zip(recs, recs[1:]):                # breakpoints between the segment's records
            cuts.add(qs + b["qs"])
        if (qs, qe) in term:                            # a terminal overhang that aligns elsewhere: its inner edge is a boundary
            cuts.add(qe if qs == 0 else qs)
    return sorted(c for c in cuts if 0 < c < qlen), f"records {len(keep)} segments {len(segs)} terminal {len(term)}"


def qmap(r):
    """query position -> target position along one alignment record (+ strand assumed)"""
    m, q, t = {}, r["qs"], r["ts"]
    for n, op in cigar(r["cg"]):
        if op in "M=X":
            for i in range(n):
                m[q + i] = t + i
            q += n; t += n
        elif op == "I":
            for i in range(n):
                m[q + i] = t
            q += n
        elif op in "DN":
            t += n
    m[q] = t
    return m


def truth_bounds(world, fam, copy):
    ex, strand = [], "+"
    for ln in open(f"{world}/fam/{fam}/truth.gtf"):
        r = ln.split("\t")
        if len(r) > 8 and r[2] == "exon" and f'transcript_id "{copy}.iso0"' in r[8]:
            ex.append((int(r[3]), int(r[4]))); strand = r[6]
    ex.sort()
    L = [e - s + 1 for s, e in ex]
    if strand == "-":
        L = L[::-1]                                     # transcript order on the - strand
    return [sum(L[:i + 1]) for i in range(len(L) - 1)], sum(L)


def main():
    world, v = sys.argv[1], sys.argv[2]
    P = f"{world}/{v}/P.cand"
    scen = json.load(open(f"{world}/world.json"))
    kinds = json.load(open(f"{world}/kinds.json"))
    cands = [r for r in csv.DictReader(open(f"{P}.candidates.tsv"), delimiter="\t") if r["flagged"] == "1"]
    U = read_fa(f"{P}.contigs.fa")
    T = read_fa(f"{world}/transcripts.fa")
    genome = read_fa(f"{world}/genome.ref.fa")
    fam_copies = collections.defaultdict(list)
    for row in csv.DictReader(open(f"{world}/P.fam.copies.tsv"), delimiter="\t"):
        fam_copies[row["family_id"]].append((row["chrom"], int(row["start"]), int(row["end"]),
                                             [tuple(map(int, x.split("-"))) for x in row["exons"].split(",")]))
    ghits = collections.defaultdict(list)
    for h in paf(mm2("-c --eqx -x splice:hq -uf -N 20", f"{world}/genome.ref.splice.mmi", f"{P}.contigs.fa")):
        ghits[h["q"]].append(h)
    thits = {}
    for h in paf(mm2("-c --eqx -x map-hifi -N 5", f"{world}/transcripts.fa", f"{P}.contigs.fa")):
        if h["tp"] == "P" and (h["q"] not in thits or h["m"] > thits[h["q"]]["m"]):
            thits[h["q"]] = h
    tot = collections.Counter()
    seen = set()
    for r in cands:
        c, fam = r["candidate"], r["family"]
        seen.add(fam)
        copy = fam + "E" if fam + "E" in kinds else None
        hs = ghits.get(c, [])
        G = genome[max(hs, key=lambda h: h["m"])["t"]] if hs else ""

        def realign(segs, seq):
            tmp = f"{P}.seg.fa"
            with open(tmp, "w") as fh:
                for qs, qe in segs:
                    fh.write(f">{qs}_{qe}\n{seq[qs:qe]}\n")
            out = collections.defaultdict(list)
            for h in paf(mm2("-c --eqx -x splice:hq -uf -N 20", f"{world}/genome.ref.splice.mmi", tmp)):
                a, b = map(int, h["q"].split("_"))
                out[(a, b)].append(h)
            return {sg: out.get(sg, []) for sg in segs}
        proj, note = project(hs, U[c], G, realign, fam_copies[fam])
        th = thits.get(c)
        if copy is None or th is None or th["t"] != copy or th["st"] != "+":
            print(f"{fam}\t{scen[fam]}\t{c}\tUNSCORED (truth hit {th['t'] if th else '-'})\t{note}")
            continue
        tb, tlen = truth_bounds(world, fam, copy)
        span = (th["ts"], th["te"])
        inside = [b for b in tb if span[0] + MIN_SEG <= b <= span[1] - MIN_SEG]
        mp = qmap(th)
        pt = [mp.get(b) for b in (proj or [])]
        res = []
        for tol in TOL:
            rec = sum(any(p is not None and abs(p - b) <= tol for p in pt) for b in inside)
            prec = sum(any(abs(p - b) <= tol for b in tb) for p in pt if p is not None)
            res.append((rec, prec))
        exact = res[0][0] == len(inside) and res[0][1] == len(pt)
        if not exact and "-v" in sys.argv:
            print(f"   DEBUG {c}: true {inside} projected->truth {pt} raw {proj}")
            for h in ghits.get(c, []):
                print(f"     rec {h['tp']} {h['t']}:{h['ts']}-{h['te']} q {h['qs']}-{h['qe']} {h['cg'][:300]}")
        tot["scored"] += 1; tot["exact"] += exact; tot["tol5"] += res[2][0] == len(inside) and res[2][1] == len(pt)
        print(f"{fam}\t{scen[fam]}\t{c}\ttrue_in_span {len(inside)}\tprojected {len(pt)}\t" +
              "\t".join(f"±{t}: rec {a}/{len(inside)} prec {b}/{len(pt)}" for t, (a, b) in zip(TOL, res)) +
              f"\t{'EXACT' if exact else ''}\t{note}\tspan {span[0]}-{span[1]}/{tlen}")
    missing = [f"{f}({scen[f]})" for f in sorted(scen) if f + "E" in kinds and f not in seen]
    print(f"TOTAL scored {tot['scored']} exact {tot['exact']} within±5 {tot['tol5']}; E families not flagged: {missing}")


if __name__ == "__main__":
    main()

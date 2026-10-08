#!/usr/bin/env python3
"""Prereg Addenda J2/K0: copy-level DNA nodes (atoms from SD alignment boundaries) and edges from the SD alignments.

usage: dna_sd_atoms.py <sd.bed> <human|gorilla> <contig,contig,...|@contigs.txt> <out_prefix> [linear|cigar]
  linear (J2): side A projected onto side B by linear interpolation (gap-blind; the only option without a CIGAR)
  cigar  (K0): exact aligned blocks from the SEDEF CIGAR (gorilla column 33); a row without a CIGAR aborts
writes <out>.nodes.tsv (idx chrom start end), <out>.edges.tsv (i j identity coverage; coverage = min over both atoms,
the pre-registered form) and <out>.edges_shorter.tsv (same, coverage of the SHORTER atom; disclosed secondary).
"""
import bisect
import collections
import heapq
import re
import sys

MIN_BLOCK = 1000
sd_path, fmt, out = sys.argv[1], sys.argv[2], sys.argv[4]
contigs = set(open(sys.argv[3][1:]).read().split()) if sys.argv[3].startswith("@") else set(sys.argv[3].split(","))
MODE = sys.argv[5] if len(sys.argv) > 5 else "linear"
# column indices (0-based): matches, mismatches, fracMatch
COLS = {"human": (18, 19, 22), "gorilla": (16, 17, 20)}[fmt]
CIGAR_COL = {"gorilla": 32}.get(fmt)
if MODE == "cigar" and CIGAR_COL is None:
    sys.exit("ABORT: cigar mode needs a SEDEF table with a CIGAR column")

pairs = []
bad = 0
for line in open(sd_path):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) <= max(COLS) or f[0] not in contigs or f[3] not in contigs:
        continue
    try:
        a1, a2, b1, b2 = int(f[1]), int(f[2]), int(f[4]), int(f[5])
        m, mm, frac = float(f[COLS[0]]), float(f[COLS[1]]), float(f[COLS[2]])
    except ValueError:
        continue
    ident = m / (m + mm) if m + mm > 0 else 0.0
    if abs(ident - frac) > 1e-4:
        bad += 1
    if a2 <= a1 or b2 <= b1:
        continue
    cig = None
    if MODE == "cigar":
        if len(f) <= CIGAR_COL or not re.fullmatch(r"(\d+[MID])+", f[CIGAR_COL]):
            sys.exit(f"ABORT: row without a CIGAR: {line[:120]}")
        cig = f[CIGAR_COL]
    pairs.append((f[0], a1, a2, f[8], f[3], b1, b2, f[9], ident, cig))
if bad:
    sys.exit(f"ABORT: {bad} rows where matches/(matches+mismatches) != fracMatch (column mapping wrong)")
print(f"[atoms] {len(pairs)} SD pairs on {len(contigs)} contigs; identity column check passed; edge mode {MODE}", file=sys.stderr)

# ---------- atoms ----------
sides = collections.defaultdict(list)  # chrom -> (start, end, sig_id)
for k, (ca, a1, a2, _, cb, b1, b2, _, _, _) in enumerate(pairs):
    sides[ca].append((a1, a2, (k, 0)))
    sides[cb].append((b1, b2, (k, 1)))

atoms = []  # (chrom, start, end)
for chrom in sorted(sides):
    ev = collections.defaultdict(lambda: ([], []))
    for s, e, sid in sides[chrom]:
        ev[s][0].append(sid)
        ev[e][1].append(sid)
    bps = sorted(ev)
    cover = set()
    segs = []  # [start, end, frozenset]
    for x, y in zip(bps, bps[1:]):
        for sid in ev[x][1]:
            cover.discard(sid)
        for sid in ev[x][0]:
            cover.add(sid)
        if not cover:
            continue
        sig = frozenset(cover)
        if segs and segs[-1][1] == x and segs[-1][2] == sig:
            segs[-1][1] = y
        else:
            segs.append([x, y, sig])
    # sliver absorption over a doubly linked list, shortest first
    n = len(segs)
    prev = list(range(-1, n - 1))
    nxt = list(range(1, n + 1))
    nxt[-1] = -1 if n else -1
    alive = [True] * n
    heap = [(s[1] - s[0], i) for i, s in enumerate(segs) if s[1] - s[0] < MIN_BLOCK]
    heapq.heapify(heap)

    def jac(p, q):
        u = len(p | q)
        return len(p & q) / u if u else 0.0

    def unlink(i):
        alive[i] = False
        p, q = prev[i], nxt[i]
        if p != -1:
            nxt[p] = q
        if q != -1:
            prev[q] = p

    def merge_identical(i):
        # merge i with touching neighbours carrying the identical signature; returns surviving index
        p = prev[i]
        if p != -1 and segs[p][1] == segs[i][0] and segs[p][2] == segs[i][2]:
            segs[p][1] = segs[i][1]
            unlink(i)
            i = p
        q = nxt[i]
        if q != -1 and segs[i][1] == segs[q][0] and segs[i][2] == segs[q][2]:
            segs[i][1] = segs[q][1]
            unlink(q)
        return i

    while heap:
        ln, i = heapq.heappop(heap)
        if not alive[i] or segs[i][1] - segs[i][0] != ln or ln >= MIN_BLOCK:
            continue
        p, q = prev[i], nxt[i]
        left = p if p != -1 and segs[p][1] == segs[i][0] else -1
        right = q if q != -1 and segs[i][1] == segs[q][0] else -1
        if left == -1 and right == -1:
            unlink(i)
            continue
        if left != -1 and (right == -1 or jac(segs[i][2], segs[left][2]) >= jac(segs[i][2], segs[right][2])):
            segs[left][1] = segs[i][1]
            unlink(i)
            j = merge_identical(left)
        else:
            segs[right][0] = segs[i][0]
            unlink(i)
            j = merge_identical(right)
        L = segs[j][1] - segs[j][0]
        if L < MIN_BLOCK:
            heapq.heappush(heap, (L, j))
    for i, s in enumerate(segs):
        if alive[i] and s[1] - s[0] >= MIN_BLOCK:
            atoms.append((chrom, s[0], s[1]))

by_chrom = collections.defaultdict(list)
for idx, (c, s, e) in enumerate(atoms):
    by_chrom[c].append((s, e, idx))
starts = {c: [x[0] for x in v] for c, v in by_chrom.items()}
maxlen = {c: max(e - s for s, e, _ in v) for c, v in by_chrom.items()}
L = [e - s for _, s, e in atoms]
print(f"[atoms] {len(atoms)} atoms; length median {sorted(L)[len(L) // 2] if L else 0} bp, max {max(L) if L else 0} bp",
      file=sys.stderr)


def overlapping(c, s, e):
    if c not in by_chrom:
        return []
    lo = bisect.bisect_left(starts[c], s - maxlen[c])
    hi = bisect.bisect_left(starts[c], e)
    return [x for x in by_chrom[c][lo:hi] if x[1] > s]


# ---------- edges ----------
cov_iv = collections.defaultdict(lambda: collections.defaultdict(list))  # (i,j) -> atom -> intervals
wsum = collections.Counter()
wid = collections.Counter()
def add(xi, x1, x2, yi, y1, y2, ident):
    key = (min(xi, yi), max(xi, yi))
    cov_iv[key][xi].append((x1, x2))
    cov_iv[key][yi].append((y1, y2))
    wsum[key] += y2 - y1
    wid[key] += (y2 - y1) * ident


def cigar_blocks(a1, b1, b2, sb, cig):
    """Exact gapless blocks (A start, B genomic start, length, reversed) of one SEDEF alignment."""
    pa = pb = 0
    for n, op in re.findall(r"(\d+)([MID])", cig):
        n = int(n)
        if op == "M":
            yield a1 + pa, (b1 + pb if sb == "+" else b2 - pb - n), n
            pa += n
            pb += n
        elif op == "D":
            pa += n
        else:
            pb += n


if MODE == "cigar":
    for ca, a1, a2, sa, cb, b1, b2, sb, ident, cig in pairs:
        rev = sb == "-"
        for As, Bs, n in cigar_blocks(a1, b1, b2, sb, cig):
            for xs, xe, xi in overlapping(ca, As, As + n):
                u1, u2 = max(xs, As), min(xe, As + n)
                if u2 <= u1:
                    continue
                o1, o2 = u1 - As, u2 - As  # offsets into the block
                m1, m2 = (Bs + o1, Bs + o2) if not rev else (Bs + n - o2, Bs + n - o1)
                for ys, ye, yi in overlapping(cb, m1, m2):
                    if yi == xi:
                        continue
                    v1, v2 = max(ys, m1), min(ye, m2)
                    if v2 <= v1:
                        continue
                    x1, x2 = (As + (v1 - Bs), As + (v2 - Bs)) if not rev else (As + (Bs + n - v2), As + (Bs + n - v1))
                    add(xi, x1, x2, yi, v1, v2, ident)
for ca, a1, a2, sa, cb, b1, b2, sb, ident, _ in (pairs if MODE == "linear" else []):
    for (c1, p1, p2, s1), (c2, q1, q2, s2) in (((ca, a1, a2, sa), (cb, b1, b2, sb)), ((cb, b1, b2, sb), (ca, a1, a2, sa))):
        LA, LB = p2 - p1, q2 - q1
        same = s1 == s2
        for xs, xe, xi in overlapping(c1, p1, p2):
            u1, u2 = max(xs, p1), min(xe, p2)
            r1, r2 = (u1 - p1) / LA, (u2 - p1) / LA
            m1, m2 = (q1 + r1 * LB, q1 + r2 * LB) if same else (q2 - r2 * LB, q2 - r1 * LB)
            for ys, ye, yi in overlapping(c2, int(m1), int(m2) + 1):
                if yi == xi:
                    continue
                o1, o2 = max(ys, m1), min(ye, m2)
                if o2 <= o1:
                    continue
                if same:
                    x1, x2 = p1 + (o1 - q1) / LB * LA, p1 + (o2 - q1) / LB * LA
                else:
                    x1, x2 = p1 + (q2 - o2) / LB * LA, p1 + (q2 - o1) / LB * LA
                add(xi, x1, x2, yi, o1, o2, ident)


def union_len(iv):
    tot, cur_s, cur_e = 0.0, None, None
    for s, e in sorted(iv):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                tot += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        tot += cur_e - cur_s
    return tot


with open(f"{out}.nodes.tsv", "w") as fh:
    fh.write("idx\tchrom\tstart\tend\n")
    for idx, (c, s, e) in enumerate(atoms):
        fh.write(f"{idx}\t{c}\t{s}\t{e}\n")
n_edges = 0
with open(f"{out}.edges.tsv", "w") as fa, open(f"{out}.edges_shorter.tsv", "w") as fb:
    fa.write("i\tj\tidentity\tcoverage\n")
    fb.write("i\tj\tidentity\tcoverage_shorter\n")
    for (i, j), d in sorted(cov_iv.items()):
        ci = min(1.0, union_len(d[i]) / L[i])
        cj = min(1.0, union_len(d[j]) / L[j])
        ident = wid[(i, j)] / wsum[(i, j)]
        short = ci if L[i] <= L[j] else cj
        fa.write(f"{i}\t{j}\t{ident:.6f}\t{min(ci, cj):.6f}\n")
        fb.write(f"{i}\t{j}\t{ident:.6f}\t{short:.6f}\n")
        n_edges += 1
print(f"[atoms] {n_edges} atom pairs with SD support -> {out}.{{nodes,edges,edges_shorter}}.tsv", file=sys.stderr)

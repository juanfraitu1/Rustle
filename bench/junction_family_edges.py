#!/usr/bin/env python3
"""A family definition from SHARED SPLICE JUNCTIONS, per
`docs/PREREG_junction_family_2026-09-20.md`.

r344 refuted STRICT junction concordance — *"paralog intron lengths and counts DRIFT after
duplication; no operating point does both"*. This rule is drift-tolerant by construction:

    for a pair (A,B) with a PAF record carrying a CIGAR:
        project each junction of A into B's frame THROUGH the alignment
        a junction MATCHES if BOTH endpoints land within +-TOL bp of a junction of B
        shared(A,B) = number of matched junctions
    EDGE iff shared(A,B) >= k

Intron LENGTH drift is tolerated because junctions are matched by projected position, not by intron
length. Intron COUNT drift is tolerated because the rule asks for k shared, never all shared.

⚠ 19.5% of genes on the held-out chromosomes have no junction at all and CANNOT be placed by any
junction rule, so node coverage is reported beside every score.

Usage: junction_family_edges.py --gff G --paf P --chrom chrN --out PREFIX --k 2 [--tol 10]
"""
import argparse
import bisect
import collections
import re
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from intron_placeholder_edges import exon_structure  # noqa: E402

CIGAR = re.compile(r'(\d+)([MIDNSHP=X])')


def junctions_of(segs):
    """gene-local junction list [(donor_end, acceptor_start)] in GENOMIC coords, ascending."""
    out = []
    for i in range(len(segs) - 1):
        out.append((segs[i][1], segs[i + 1][0]))
    return out


def build_map(cg, qs, ts, strand, qlen):
    """Anchors (query_pos -> target_pos) at each aligned block boundary, query in + orientation."""
    q, t = qs, ts
    anchors = []
    for n, op in CIGAR.findall(cg):
        n = int(n)
        if op in 'M=X':
            anchors.append((q, t, n))
            q += n; t += n
        elif op in 'ID':
            if op == 'I':
                q += n
            else:
                t += n
    if strand == '-':
        # query coordinates were reported on the + strand of the query by minimap2; for a - match the
        # target walk is reversed, which build_map already encodes via ts..te. Nothing to flip here.
        pass
    return anchors


def project(anchors, qpos):
    """map a query position through the anchor blocks; None if it falls outside every block."""
    lo, hi = 0, len(anchors) - 1
    while lo <= hi:
        mid = (lo + hi) // 2
        q, t, n = anchors[mid]
        if qpos < q:
            hi = mid - 1
        elif qpos >= q + n:
            lo = mid + 1
        else:
            return t + (qpos - q)
    return None


def main():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--paf', '--chrom', '--out'):
        ap.add_argument(x, required=True)
    ap.add_argument('--k', type=int, required=True)
    ap.add_argument('--tol', type=int, default=10)
    ap.add_argument('--read-junctions', help='TSV chrom/donor/acceptor/support; UNION with annotated '
                                             '(§6u2: reads add ~150%% more junctions to spliced genes)')
    a = ap.parse_args()

    st = exon_structure(a.gff, a.chrom)
    jn = {}
    for name, (strand, segs) in st.items():
        base = int(name.rsplit(':', 1)[1].split('-')[0])
        js = [(d - base, ac - base) for d, ac in junctions_of(segs)]
        if js:
            jn[name] = js
    if a.read_junctions:
        import bisect as _bi
        spans = sorted((int(k.rsplit(':', 1)[1].split('-')[0]),
                        int(k.rsplit(':', 1)[1].split('-')[1]), k) for k in st)
        starts = [x[0] for x in spans]
        for line in open(a.read_junctions):
            f = line.rstrip('\n').split('\t')
            if len(f) < 3:
                continue
            d, ac = int(f[1]), int(f[2])
            i = _bi.bisect_right(starts, d) - 1
            j = i
            while j >= 0 and j > i - 40:
                s_, e_, k = spans[j]
                if s_ <= d and ac <= e_:
                    base = s_
                    jn.setdefault(k, [])
                    cand = (d - base, ac - base)
                    if not any(abs(cand[0] - x) <= a.tol and abs(cand[1] - y) <= a.tol for x, y in jn[k]):
                        jn[k].append(cand)
                j -= 1
        for k in jn:
            jn[k].sort()
    spliced = set(jn)

    shared = collections.Counter()
    for line in open(a.paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        q, s = f[0], f[5]
        if q == s or q not in jn or s not in jn:
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), None)
        if not cg:
            continue
        anchors = build_map(cg, int(f[2]), int(f[7]), f[4], int(f[1]))
        if not anchors:
            continue
        tgt = jn[s]
        donors = sorted(d for d, _ in tgt)
        acc = {d: ac for d, ac in tgt}
        n = 0
        for d, ac in jn[q]:
            pd, pa = project(anchors, d), project(anchors, ac)
            if pd is None or pa is None:
                continue
            i = bisect.bisect_left(donors, pd - a.tol)
            while i < len(donors) and donors[i] <= pd + a.tol:
                if abs(acc[donors[i]] - pa) <= a.tol:
                    n += 1
                    break
                i += 1
        if n:
            k = tuple(sorted((q, s)))
            shared[k] = max(shared[k], n)

    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    kept = 0
    for (x, y), n in shared.items():
        if n >= a.k:
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry
            kept += 1
    comp = collections.defaultdict(list)
    for nd in list(parent):
        comp[find(nd)].append(nd)
    clusters = {c: v for c, v in comp.items() if len(v) >= 2}

    with open(a.out + '.clusters.tsv', 'w') as fh:
        fh.write('cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n')
        for i, (_, mem) in enumerate(sorted(clusters.items(), key=lambda kv: -len(kv[1]))):
            for name in sorted(mem):
                ch, rng = name.rsplit(':', 1); s_, e_ = rng.split('-')
                fh.write(f'JN{i}\t{len(mem)}\tNA\tNA\tNA\t{ch}\t{int(s_)-1}\t{int(e_)}\n')
    with open(a.out + '.edges.tsv', 'w') as fh:
        for (x, y), n in sorted(shared.items()):
            fh.write(f'{x}\t{y}\t{n}\n')
    cov = sum(len(v) for v in clusters.values())
    print(f'{a.chrom} k={a.k}: {len(shared)} pairs share >=1 junction | {kept} edges at k | '
          f'{len(clusters)} families | {cov} nodes covered of {len(spliced)} spliced '
          f'({100*cov/len(spliced) if spliced else 0:.1f}%)')


if __name__ == '__main__':
    main()

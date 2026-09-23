#!/usr/bin/env python3
"""Does a readthrough-fused de novo locus show a DISCONTINUITY in its PAF alignment-partner set
along its own coordinate axis?

Every node-split trigger refuted so far (registers 845/846/937-940/948-949/967/968) used READ-level
evidence inside the locus, and all failed the same way: a readthrough molecule is a genuine, abundant,
full-length transcript, so the "bridge" IS the population, not an outlier. Graph-structural splits were
tried on the graph's TOPOLOGY (r300 bridges 0.3%, r301 connectivity inverts, r827 lambda>=2) but never
on the POSITION of a partner's alignment along the locus's own axis, which is what this measures.

    For a cut c:  L = partners aligning left of c, R = partners aligning right of c
    score(c) = |L & R| / |L | R|          (Jaccard; low = the two halves have different relatives)
    discontinuity(locus) = min over c of score(c)

⚠ The naive form is DEGENERATE: with a weak floor the minimum is always found at an extreme cut where
one side has 2 partners. Both sides must carry real, EXCLUSIVE evidence, and the comparison must be
LENGTH-MATCHED -- a long locus has more partners and more chances to find a low-J cut, so an unmatched
fused-vs-clean comparison measures length, not fusion.

Usage:
  partner_discontinuity.py --paf dn16.paf --gff chr16.genes.gff [--locus X]...
"""
import argparse
import collections
import re
import statistics

MIN_EXCL = 5          # each side needs this many partners the other side does NOT have
CUT_LO, CUT_HI = 0.15, 0.85   # ignore cuts near the ends
NBINS = 100
MIN_FRAC = 0.05       # a partner counts for a side if it covers this much of that side


def load_paf(paf):
    hits = collections.defaultdict(list)
    qlen = {}
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 11:
            continue
        q, t = f[0], f[5]
        qlen[q] = int(f[1]); qlen[t] = int(f[6])
        if q == t:
            continue
        hits[q].append((int(f[2]), int(f[3]), t))
        hits[t].append((int(f[7]), int(f[8]), q))
    return hits, qlen


def genes_per_locus(gff, loci):
    """locus -> [gene names whose OWN span is >=50% inside the locus]."""
    spans = []
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('gene', 'pseudogene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            spans.append((int(f[3]), int(f[4]), m.group(1)))
    spans.sort()
    out = {}
    for name in loci:
        m = re.match(r'^\S+:(\d+)-(\d+)$', name)
        if not m:
            continue
        s, e = int(m.group(1)), int(m.group(2))
        inside = []
        for (gs, ge, g) in spans:
            if ge < s:
                continue
            if gs > e:
                break
            ov = min(e, ge) - max(s, gs)
            if ov > 0 and ov >= 0.5 * (ge - gs + 1):
                inside.append(g)
        out[name] = inside
    return out


def discontinuity(records, n):
    if not records or n <= 0:
        return None
    if len({p for (_, _, p) in records}) < 2 * MIN_EXCL:
        return None
    best = None
    for b in range(1, NBINS):
        frac = b / NBINS
        if frac < CUT_LO or frac > CUT_HI:
            continue
        c = n * frac
        lenL, lenR = c, n - c
        L, R = set(), set()
        for (qs, qe, p) in records:
            if min(qe, c) - qs >= MIN_FRAC * lenL:
                L.add(p)
            if qe - max(qs, c) >= MIN_FRAC * lenR:
                R.add(p)
        if len(L - R) < MIN_EXCL or len(R - L) < MIN_EXCL:
            continue
        j = len(L & R) / len(L | R)
        if best is None or j < best[0]:
            best = (j, frac, len(L - R), len(R - L), len(L & R))
    return best


def band(n):
    for hi, lbl in ((20000, '<20kb'), (50000, '20-50kb'), (100000, '50-100kb')):
        if n < hi:
            return lbl
    return '>=100kb'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--paf', required=True)
    ap.add_argument('--gff', required=True)
    ap.add_argument('--locus', action='append', default=[])
    a = ap.parse_args()

    hits, qlen = load_paf(a.paf)
    gpl = genes_per_locus(a.gff, list(qlen))

    rows = []
    for name, n in qlen.items():
        d = discontinuity(hits.get(name, []), n)
        if d is None:
            continue
        rows.append((name, n, len(gpl.get(name, [])), d))

    print(f"scoreable loci (>= {2*MIN_EXCL} partners, both sides >= {MIN_EXCL} exclusive): "
          f"{len(rows)} of {len(qlen)}\n")

    print("=== named loci of interest ===")
    print(f"{'locus':32} {'len':>8} {'genes':>6} {'minJ':>7} {'cut@':>6} {'exL':>5} {'exR':>5} {'shared':>7}")
    for name in a.locus:
        hit = [r for r in rows if r[0] == name]
        if not hit:
            print(f"{name:32} {qlen.get(name,0):>8} {'-':>6} {'not scoreable':>7}")
            continue
        _, n, ng, (j, frac, exL, exR, sh) = hit[0]
        print(f"{name:32} {n:>8} {ng:>6} {j:>7.3f} {frac:>6.2f} {exL:>5} {exR:>5} {sh:>7}")

    print("\n=== minJ by LENGTH BAND x GENES-INSIDE (the confound-controlled view) ===")
    print(f"{'band':>10} {'genes':>7} {'n':>5} {'median minJ':>12} {'mean':>7}")
    grp = collections.defaultdict(list)
    for (_, n, ng, d) in rows:
        cls = '1 gene' if ng == 1 else ('>=2 genes' if ng >= 2 else '0 genes')
        grp[(band(n), cls)].append(d[0])
    for b in ('<20kb', '20-50kb', '50-100kb', '>=100kb'):
        for cls in ('1 gene', '>=2 genes'):
            v = grp.get((b, cls), [])
            if len(v) >= 3:
                print(f"{b:>10} {cls:>7} {len(v):>5} {statistics.median(v):>12.3f} "
                      f"{statistics.mean(v):>7.3f}")


if __name__ == '__main__':
    main()

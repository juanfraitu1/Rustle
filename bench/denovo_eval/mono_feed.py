#!/usr/bin/env python3
"""Row 1665 prototype: give `families` single-exon loci it has no model for, additively.

  mono_feed.py FAMILIES.gtf FINAL.gtf most|span OUT.gtf     (OUT = FAMILIES.gtf + the added models; prints counts)
  mono_feed.py --selftest
Pile = overlapping same-strand single-exon cluster models of FINAL.gtf (attribute mono_model "cluster"). A pile is added
only when no FAMILIES.gtf exon overlaps it on the same strand ('.' matches both). One model per pile:
  most  the cluster with the most reads (ties: longest, then leftmost)
  span  the pile's merged span, unless that is > 2x its longest cluster (array blobs, row 1518): then the longest cluster
reads = sum over the pile. Nothing in FAMILIES.gtf is changed or removed.
"""
import sys, re, collections

def piles(models):
    """models [(start, end, reads)] sorted -> lists of overlapping ones."""
    out, end = [], -1
    for m in sorted(models):
        if out and m[0] < end: out[-1].append(m); end = max(end, m[1])
        else: out.append([m]); end = m[1]
    return out

def pick(pile, arm):
    longest = max(pile, key=lambda m: (m[1] - m[0], -m[0]))
    if arm == 'most': s, e, _ = max(pile, key=lambda m: (m[2], m[1] - m[0], -m[0]))
    else:
        s, e = min(m[0] for m in pile), max(m[1] for m in pile)
        if e - s > 2 * (longest[1] - longest[0]): s, e = longest[:2]
    return s, e, sum(m[2] for m in pile)

def free(s, e, exons):
    """no exon of the sorted non-overlapping-ish list overlaps [s, e)"""
    return not any(a < e and s < b for a, b in exons)

def main(fam, final, arm, out):
    ex = collections.defaultdict(list); mono = collections.defaultdict(list)
    for l in open(fam):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[(f[0], f[6])].append((int(f[3]) - 1, int(f[4])))
    for l in open(final):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'transcript' and 'mono_model "cluster"' in f[8]:
            mono[(f[0], f[6])].append((int(f[3]) - 1, int(f[4]), int(re.search(r'reads "(\d+)"', f[8]).group(1))))
    import bisect
    n = k = 0
    with open(out, 'w') as w:
        w.write(open(fam).read())
        for (c, st), ms in sorted(mono.items()):
            near = sorted(ex[(c, st)] + ex[(c, '.')] + (ex[(c, '+')] + ex[(c, '-')] if st == '.' else []))
            starts = [a for a, _ in near]; mx = max((b - a for a, b in near), default=0)
            for p in piles(ms):
                n += 1; ps, pe = min(m[0] for m in p), max(m[1] for m in p)
                lo = bisect.bisect_left(starts, ps - mx); hi = bisect.bisect_left(starts, pe)
                if not free(ps, pe, near[lo:hi]): continue
                s, e, r = pick(p, arm); k += 1; i = f'MONOFEED_{c}_{st}_{s}_{e}'
                a = f'gene_id "{i}"; transcript_id "{i}.1"; reads "{r}"; mono_feed "{arm}";'
                w.write(f'{c}\trustle\ttranscript\t{s + 1}\t{e}\t.\t{st}\t.\t{a}\n{c}\trustle\texon\t{s + 1}\t{e}\t.\t{st}\t.\t{a}\n')
    print(f'piles {n} added {k}')

def selftest():
    p = piles([(0, 100, 3), (50, 300, 9), (290, 320, 1), (500, 600, 2)])
    assert [len(x) for x in p] == [3, 1]
    assert pick(p[0], 'most') == (50, 300, 13) and pick(p[0], 'span') == (0, 320, 13)
    assert pick([(0, 100, 3), (90, 200, 3), (190, 300, 3)], 'span') == (90, 200, 9)     # chained blob: longest cluster
    assert free(10, 20, [(20, 30)]) and not free(10, 21, [(20, 30)])
    print('selftest ok')

if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main(*sys.argv[1:5])

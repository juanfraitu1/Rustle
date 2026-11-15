#!/usr/bin/env python3
"""PREREG row 1725 step 0: emitted chains that are the weaker copy of a stronger emitted chain with one exon end moved.

  shift_class.py REF_GTF DEFAULT_GTF FASTA [EXTRA_GFF3 | -] [OUT_TSV]   |   --selftest
S1 = one intron differs and shares one boundary (S1t when the moved boundary belongs to a terminal exon);
S2 = two adjacent introns differ, the exon between them moved at both ends and still overlaps."""
import sys, re, gzip, collections
from alt_events import load
op = lambda p: gzip.open(p, 'rt') if p.endswith('.gz') else open(p)

def partners(c, one, two):
    """-> [(class, differing introns of c, partner chain, shift bp)] over chains in the index."""
    out = []
    for i in range(len(c)):
        for x in one.get((c[:i], c[i + 1:]), ()):
            if (x[0] == c[i][0]) != (x[1] == c[i][1]):
                donor_side = x[1] == c[i][1]  # the intron start moved = end of the exon before it
                term = (i == 0 and donor_side) or (i == len(c) - 1 and not donor_side)
                out.append(('S1t' if term else 'S1', (c[i],), c[:i] + (x,) + c[i + 1:], abs(x[0] - c[i][0]) + abs(x[1] - c[i][1])))
        if i + 1 < len(c):
            for x, y in two.get((c[:i], c[i + 2:]), ()):
                a, b = c[i], c[i + 1]
                if x[0] == a[0] and y[1] == b[1] and x[1] != a[1] and y[0] != b[0] and x[1] < b[0] and a[1] < y[0]:
                    out.append(('S2', (a, b), c[:i] + (x, y) + c[i + 2:], max(abs(x[1] - a[1]), abs(y[0] - b[0]))))
    return out

def index(chains):
    one = collections.defaultdict(set); two = collections.defaultdict(set)
    for d in chains:
        for i in range(len(d)):
            one[d[:i], d[i + 1:]].add(d[i])
            if i + 1 < len(d): two[d[:i], d[i + 2:]].add((d[i], d[i + 1]))
    return one, two

def motif(seq, strand, i):
    m = (seq[i[0]:i[0] + 2] + seq[i[1] - 2:i[1]]).upper()
    if strand == '-': m = m[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    return {'GTAG': 'GT-AG', 'GCAG': 'GC-AG', 'ATAC': 'AT-AC'}.get(m, 'other')

def gff3_introns(p):
    ex = collections.defaultdict(list)
    for l in op(p):
        f = l.split('\t')
        if len(f) > 8 and f[2] == 'exon': ex[f[0], f[6], re.search(r'Parent=([^;\s]+)', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
    out = collections.defaultdict(set)
    for k, e in ex.items():
        e.sort(); out[k[:2]].update((a[1], b[0]) for a, b in zip(e, e[1:]))
    return out

def selftest():
    A = ((10, 20), (30, 40), (50, 60)); one, two = index([A])
    assert [p[:2] + p[3:] for p in partners(((10, 20), (32, 40), (50, 60)), one, two)] == [('S1', ((32, 40),), 2)]
    assert partners(((12, 20), (30, 40), (50, 60)), one, two)[0][0] == 'S1t' and partners(((10, 22), (30, 40), (50, 60)), one, two)[0][0] == 'S1'
    assert partners(((10, 20), (30, 40), (50, 63)), one, two)[0][0] == 'S1t'
    assert [p[0] for p in partners(((10, 22), (33, 40), (50, 60)), one, two)] == ['S2'] and partners(A, one, two) == []
    assert partners(((10, 20), (30, 60)), one, two) == [] and motif('xxCTxxxxACxx', '-', (2, 10)) == 'GT-AG'
    print('selftest ok')

if __name__ == '__main__':
    if sys.argv[1] == '--selftest': selftest(); sys.exit()
    import pysam
    ref = load(sys.argv[1]); D = load(sys.argv[2], 'reads'); fa = pysam.FastaFile(sys.argv[3])
    extra = gff3_introns(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] != '-' else {}
    out = open(sys.argv[5], 'w') if len(sys.argv) > 5 else None
    if out: print('contig\tstrand\tclass\tshift\treads\tpartner_reads\tin_ref\tintron_annotated\tmotif\tintrons\tpartner_in_ref\tchain\tpartner_introns', file=out)  # every emitted chain; class - = not in a class
    T = collections.Counter(); allc = []; rows = []
    for k in sorted(D):
        d = {c: max(t[2] or 0 for t in v) for c, v in D[k].items()}; r = ref.get(k, {}); ri = {i for ch in r for i in ch} | extra.get(k, set())
        one, two = index(d); seq = None
        for c, n in d.items():
            ps = [p for p in partners(c, one, two) if d[p[2]] > n or (d[p[2]] == n and p[2] < c)]  # weaker side; ties: the later chain
            allc.append((n, c in r, bool(ps) and any(p[0] != 'S1t' for p in ps)))
            if not ps:
                if out: print(k[0], k[1], '-', 0, n, 0, int(c in r), 0, '-', '-', 0, c, '-', sep='\t', file=out)
                continue
            cls, new, par, sh = min(ps, key=lambda p: (p[0] == 'S1t', -d[p[2]]))  # internal classes first, then the strongest partner
            if seq is None: seq = fa.fetch(k[0]) if k[0] in fa.references else ''
            ia = all(i in ri for i in new); mo = '/'.join(sorted({motif(seq, k[1], i) for i in new})) if seq else 'na'; q = n / d[par]
            if out: print(k[0], k[1], cls, sh, n, d[par], int(c in r), int(ia), mo, new, int(par in r), c, tuple(sorted(set(par) - set(c))), sep='\t', file=out)
            rows.append((cls, sh, q, mo, c in r, ia))
    N = len(allc); ok = sum(a[1] for a in allc)
    def tab(name, sel):
        v = [x for x in rows if sel(x)]; n = len(v) or 1
        print(f'{name}\t{len(v)}\t{sum(x[4] for x in v)}\t{sum(x[4] for x in v) / n:.3f}\t{sum(x[5] for x in v)}\t{sum(x[5] for x in v) / n:.3f}')
    print(f'emitted chains {N}, in RefSeq {ok} (precision {ok / N:.3f})\nstratum\tchains\tin_RefSeq\tshare\tintron_annotated\tshare')
    for cls in ('S1', 'S2', 'S1t'):
        tab(cls, lambda x: x[0] == cls)
        for nm, f in (('<=10bp', lambda x: x[1] <= 10), ('11-50bp', lambda x: 10 < x[1] <= 50), ('>50bp', lambda x: x[1] > 50), ('multiple of 3', lambda x: x[1] % 3 == 0), ('not multiple of 3', lambda x: x[1] % 3 != 0),
                      ('ratio <0.1', lambda x: x[2] < 0.1), ('ratio 0.1-0.5', lambda x: 0.1 <= x[2] < 0.5), ('ratio >=0.5', lambda x: x[2] >= 0.5),
                      ('<=10bp, ratio <0.1', lambda x: x[1] <= 10 and x[2] < 0.1), ('<=10bp, ratio >=0.5', lambda x: x[1] <= 10 and x[2] >= 0.5)) + tuple((m, lambda x, m=m: x[3] == m) for m in sorted({x[3] for x in rows if x[0] == cls})):
            tab(f'  {cls} {nm}', lambda x: x[0] == cls and f(x))
    drop = sum(a[2] and not a[1] for a in allc); una = sum(1 for x in rows if x[0] != 'S1t' and not x[5])
    allc.sort(key=lambda a: a[0]); curve = sum(a[1] for a in allc[drop:]) / (N - drop)  # stable sort = file order within a read count
    print(f'S1 + S2 chains {sum(a[2] for a in allc)}; with an unannotated differing intron {una} ({una / N:.4f} of emitted; bar >= 0.01)')
    print(f'ceiling: drop the {drop} S1 + S2 chains not in RefSeq -> precision {ok / (N - drop):.4f}; read-count curve at {drop} dropped {curve:.4f}; difference {ok / (N - drop) - curve:+.4f} (bar +0.02)')

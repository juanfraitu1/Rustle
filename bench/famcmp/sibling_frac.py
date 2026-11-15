#!/usr/bin/env python3
"""PREREG row 1719: chains the isoform-fraction rule drops that differ from an emitted chain by one event.

  sibling_frac.py REF_GTF DEFAULT_GTF FRACTION_OFF_GTF FASTA [PAIRS_TSV of alt_events.py] [SKIP_CONTIG,...]   |   --selftest"""
import sys, csv, ast, collections
from alt_events import load
MIN = 5

def index(chains):
    one = collections.defaultdict(set); two = collections.defaultdict(set)
    for d in chains:
        for i in range(len(d)):
            one[d[:i], d[i + 1:]].add(d[i])
            if i + 1 < len(d): two[d[:i], d[i + 2:]].add((d[i][0], d[i + 1][1]))
    return one, two

def sibling(c, one, two):
    """-> introns of c absent from a sibling chain that differs by one event (shifted boundary / cassette in or out), or None."""
    for i in range(len(c)):
        if any((x[0] == c[i][0]) != (x[1] == c[i][1]) for x in one.get((c[:i], c[i + 1:]), ())): return [c[i]]
        if c[i] in two.get((c[:i], c[i + 1:]), ()): return [c[i]]                                   # cassette exon skipped
        if i + 1 < len(c) and (c[i][0], c[i + 1][1]) in one.get((c[:i], c[i + 2:]), ()): return [c[i], c[i + 1]]   # cassette exon included
    return None

def canonical(seq, strand, introns):
    return all((seq[a:a + 2].upper(), seq[b - 2:b].upper()) == (('GT', 'AG') if strand == '+' else ('CT', 'AC')) for a, b in introns)

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        A = ((10, 20), (30, 40), (50, 60)); one, two = index([A])
        assert sibling(((10, 20), (32, 40), (50, 60)), one, two) == [(32, 40)] and sibling(((10, 20), (30, 60)), one, two) == [(30, 60)]
        assert sibling(((10, 20), (30, 34), (36, 40), (50, 60)), one, two) == [(30, 34), (36, 40)] and sibling(((10, 20), (31, 41), (50, 60)), one, two) is None
        assert sibling(A, one, two) is None and canonical('xxGTxxxxAGxx', '+', [(2, 10)]) and canonical('xxCTxxxxACxx', '-', [(2, 10)]) and not canonical('xxGTxxxxAGxx', '-', [(2, 10)])
        print('selftest ok'); sys.exit()
    import pysam
    ref_p, def_p, off_p, fa_p = sys.argv[1:5]; pairs = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] != '-' else None
    skip = set(sys.argv[6].split(',')) if len(sys.argv) > 6 else set()
    ref = load(ref_p); D = load(def_p); N = load(off_p, 'reads'); fa = pysam.FastaFile(fa_p); T = collections.Counter(); kept = set()
    for k in sorted(set(D) | set(N)):
        if k[0] in skip: continue
        d = D.get(k, {}); n = N.get(k, {}); r = ref.get(k, {}); ri = {i for ch in r for i in ch}
        T['default'] += len(d); T['default annotated'] += sum(ch in r for ch in d); T['default chains absent from the fraction-off run'] += sum(ch not in n for ch in d)
        cand = [ch for ch, v in n.items() if ch not in d and max(t[2] or 0 for t in v) >= MIN]
        if not cand: continue
        one, two = index(d); seq = fa.fetch(k[0]) if k[0] in fa.references else None
        for ch in cand:
            new = sibling(ch, one, two)
            g = 'not a sibling' if new is None else 'sibling, GT-AG' if seq and canonical(seq, k[1], new) else 'sibling, not GT-AG'
            T[g] += 1; T[g, 'ann'] += ch in r; T[g, 'introns'] += all(i in ri for i in ch)
            if g == 'sibling, GT-AG': kept.add((k[0], k[1], ch))
    print(f"default chains {T['default']}, annotated {T['default annotated']} (precision {T['default annotated'] / T['default']:.3f}); absent from the fraction-off run {T['default chains absent from the fraction-off run']}")
    for g in ('sibling, GT-AG', 'sibling, not GT-AG', 'not a sibling'):
        n = T[g] or 1; print(f"{g}: {T[g]} candidates; annotated {T[g, 'ann']} ({T[g, 'ann'] / n:.3f}); every intron annotated {T[g, 'introns']} ({T[g, 'introns'] / n:.3f})")
    if pairs:
        csv.field_size_limit(1 << 30); P = collections.Counter()
        for r in csv.DictReader(open(pairs), delimiter='\t'):
            a, b = int(r['reads_a']), int(r['reads_b'])
            if min(a, b) < MIN or r['default'] != 'major' or r['contig'] in skip: continue
            g = 'SITE' if r['class'].startswith('SITE') else r['class']; m = ast.literal_eval(r['chain_b'] if b <= a else r['chain_a'])
            P[g] += 1; P[g, 'kept'] += (r['contig'], r['strand'], m) in kept
        print('row 1714 lost pairs (minor >= 5) restored: ' + '; '.join(f"{g} {P[g, 'kept']} of {P[g]}" for g in ('SITE', 'CASSETTE', 'AFE', 'ALE', 'NEST5', 'NEST3')))

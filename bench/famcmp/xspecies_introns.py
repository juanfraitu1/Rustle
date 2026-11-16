#!/usr/bin/env python3
"""Addendum 2 of PREREG_ptr_ppa_frozen_rules_2026-11-15: introns of one sample's default chains, lifted to another
species' assembly, that the other sample's reads also splice.

  xspecies_introns.py Q.gtf Q.refseq.gtf Q.cat.gff3 MAP.paf LIFT.sh T.bam T.refseq.gtf T.cat.gff3 WORKDIR   |   --selftest"""
import sys, re, collections, subprocess, pysam

def models(p, key, reads=False):
    ex = collections.defaultdict(list); rd = {}; sd = {}
    for l in open(p):
        f = l.split('\t')
        if len(f) < 9: continue
        k = (f[0], re.search(key, f[8]).group(1))
        if f[2] == 'exon': ex[k].append((int(f[3]), int(f[4]))); sd[k] = f[6]
        elif reads and f[2] == 'transcript': rd[k] = int(re.search(r'reads "(\d+)"', f[8]).group(1))
    out = {}
    for k, e in ex.items():
        if len(e) > 1:
            e.sort(); c = (k[0], sd[k], tuple((a[1], b[0] - 1) for a, b in zip(e, e[1:])))  # introns 0-based half-open
            out[c] = max(out.get(c, 0), rd.get(k, 1))
    return out

def classify(ch, c, s, A, I, D, Ac):
    """-> {intron: class} for one chain: Y exact, N1 the single new intron sharing a site, N0 no annotated intron, N2 2+ new."""
    if (c, s, ch) in A: return {i: 'Y' for i in ch}
    new = [i for i in ch if (c, s, i) not in I]
    if new and len(new) == len(ch): return {i: 'N0' for i in ch}
    if len(new) == 1: return {new[0]: 'N1'} if ((c, s, new[0][0]) in D or (c, s, new[0][1]) in Ac) else {}
    return {i: 'N2' for i in new} if len(new) > 1 else {}

def annot(*paths_keys):
    A = set()
    for p, k in paths_keys: A |= set(models(p, k))
    I = {(c, s, i) for c, s, ch in A for i in ch}
    return A, I, {(c, s, i[0]) for c, s, i in I}, {(c, s, i[1]) for c, s, i in I}

def main(qg, qr, qc, paf, lift, tbam, tr, tc, wd):
    T = r'transcript_id "([^"]+)"'; P = r'Parent=(\S+)'
    A, I, D, Ac = annot((qr, T), (qc, P)); cls = {}
    for (c, s, ch), rd in models(qg, T, True).items():
        if rd < 6: continue
        st = '21+' if rd >= 21 else '6-20'
        for i, k in classify(ch, c, s, A, I, D, Ac).items():
            o = cls.get((c, i))
            if o is None or (st, k) < o: cls[(c, i)] = (st, k)   # an intron keeps its best stratum, then Y before N
    with open(f'{wd}/q.bed', 'w') as o:
        for (c, (a, b)) in sorted(cls): o.write(f'{c}\t{a}\t{b}\n')
    lifted = subprocess.run(['bash', lift, paf, f'{wd}/q.bed'], capture_output=True, text=True, check=True).stdout
    sup = collections.Counter()
    for r in pysam.AlignmentFile(tbam, threads=8):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        pos = r.reference_start
        for op, n in r.cigartuples:
            if op == 3: sup[(r.reference_name, pos, pos + n)] += 1
            if op in (0, 2, 3, 7, 8): pos += n
    _, TI, _, _ = annot((tr, T), (tc, P)); TI = {(c, i) for c, _, i in TI}
    n = collections.Counter(cls.values()); L = collections.Counter(); S = collections.Counter(); AN = collections.Counter()
    for l in lifted.splitlines():
        f = l.split('\t'); c, a, b = f[3].rsplit('_', 2); k = cls[(c, (int(a), int(b)))]; t = (f[0], int(f[1]), int(f[2]))
        L[k] += 1; S[k] += sup[t] >= 2; AN[k] += (t[0], (t[1], t[2])) in TI
    print('stratum\tclass\tintrons\tlifted\tshared (>= 2 target reads)\trate\tx yardstick\tannotated in the target')
    for st in ('21+', '6-20'):
        y = S[(st, 'Y')] / max(1, L[(st, 'Y')])
        for k in ('Y', 'N1', 'N0', 'N2'):
            r = S[(st, k)] / max(1, L[(st, k)]); print(f"{st}\t{k}\t{n[(st, k)]}\t{L[(st, k)]}\t{S[(st, k)]}\t{r:.3f}\t{r / y if y else 0:.2f}\t{AN[(st, k)] / max(1, L[(st, k)]):.3f}")

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        A = {('c', '+', ((10, 20), (30, 40)))}; I = {('c', '+', (10, 20)), ('c', '+', (30, 40))}; D = {('c', '+', 10), ('c', '+', 30)}; Ac = {('c', '+', 20), ('c', '+', 40)}
        assert classify(((10, 20), (30, 40)), 'c', '+', A, I, D, Ac) == {(10, 20): 'Y', (30, 40): 'Y'}
        assert classify(((10, 20), (30, 45)), 'c', '+', A, I, D, Ac) == {(30, 45): 'N1'}
        assert classify(((11, 21), (31, 41)), 'c', '+', A, I, D, Ac) == {(11, 21): 'N0', (31, 41): 'N0'}
        assert classify(((10, 20), (31, 41)), 'c', '+', A, I, D, Ac) == {}          # one new intron, no shared site
        assert classify(((10, 20), (30, 45), (50, 60)), 'c', '+', A, I, D, Ac) == {(30, 45): 'N2', (50, 60): 'N2'}
        print('selftest ok')
    else: main(*sys.argv[1:10])

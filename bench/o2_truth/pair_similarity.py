#!/usr/bin/env python3
"""Are two annotated genes similar over their non-repeat exonic sequence? (register 1957)

    pair_similarity.py PAIRS.tsv RELATIONS.tsv RMSK.bed GENOME.fa OUT_DIR [THREADS [BAND]]   (BAND = the band's lower bound as place_classes wrote it, default 0.9)

PAIRS.tsv = place_classes.py's PLACE_CLASSES_PAIRS (band, class, place, source gene, genes under, repeat share).
Each gene: exons merged and spliced, repeat bases as N; all against all with minimap2; a pair is similar when the aligned
non-N bases cover >= 30% of the shorter sequence's non-N bases and >= 200 bases. Needs DENOVO_CAT_GFF / DENOVO_CAT_TSV."""
import sys, os, random, subprocess, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import place_classes as PC

def union_len(iv): return sum(b - a for a, b in PC.merged(iv))

def similar(cov, la, lb): return cov >= 200 and cov >= 0.3 * min(la, lb)

def main(pairs, relp, rmsk, fa, out, thr=16, band='0.9'):
    import pysam, famtruth as F
    os.makedirs(out, exist_ok=True); CAT = F.load_cat(); rows = [l.rstrip('\n').split('\t') for l in open(pairs)]
    rows = [r for r in rows if r[0] == band and r[3] != '.']; genes = sorted({r[3] for r in rows} | {g for r in rows for g in r[4].split(',') if g != '.'})
    rp = collections.defaultdict(list)
    for l in open(rmsk):
        f = l.split('\t', 8)
        if f[6] in PC.REP: rp[f[0]].append((int(f[1]), int(f[2])))
    rp = {c: PC.merged(v) for c, v in rp.items()}; rs = {c: [a for a, _ in v] for c, v in rp.items()}
    FA = pysam.FastaFile(fa); L = {}
    with open(out + '/genes.fa', 'w') as w:
        for g in genes:
            d = CAT[g]; c = d['chrom']; seq = []
            for a, b in PC.merged(d['exons']):
                s = list(FA.fetch(c, a, b).upper()); k = max(0, __import__('bisect').bisect_right(rs.get(c, []), a) - 1); v = rp.get(c, [])
                while k < len(v) and v[k][0] < b:
                    for i in range(max(a, v[k][0]) - a, min(b, v[k][1]) - a): s[i] = 'N'
                    k += 1
                seq.append(''.join(s))
            seq = ''.join(seq); L[g] = len(seq) - seq.count('N'); w.write(f'>{g}\n{seq}\n')
    subprocess.run(f'minimap2 -c -X -x asm20 -k15 -w10 -N 200 -p 0.05 -t {thr} {out}/genes.fa {out}/genes.fa > {out}/all.paf 2> {out}/mm.err', shell=True, check=True)
    iv = collections.defaultdict(list)
    for l in open(out + '/all.paf'):
        f = l.split('\t'); a, b = sorted((f[0], f[5])); iv[(a, b)].append((int(f[2]), int(f[3])) if f[0] == a else (int(f[7]), int(f[8])))
    cov = {k: union_len(v) for k, v in iv.items()}   # shortcut: aligned span on one side, N bases inside it are few because minimap2 does not seed in N
    sim = lambda a, b: a != b and similar(cov.get(tuple(sorted((a, b))), 0), L[a], L[b])
    rel = collections.defaultdict(set)
    for l in list(open(relp))[1:]:
        f = l.rstrip('\n').split('\t'); rel[f[0]].add(f[2]); rel[f[2]].add(f[0])
    print(f'genes {len(genes)}; with fewer than 200 non-repeat exonic bases {sum(v < 200 for v in L.values())}; aligned pairs {len(cov)}')
    print(f'rows of band {band}\trows\tsimilar\tshare\tdistinct places\tplaces similar')
    T = collections.defaultdict(lambda: [0, 0, set(), set()])
    for band, k, place, s, under, share in rows:
        key = k if float(share) < 0.5 or k[0] == '4' else k + ' (repeat-dominated)'
        ok = any(sim(s, g) for g in under.split(',') if g != '.'); t = T[key]; t[0] += 1; t[1] += ok; t[2].add(place); ok and t[3].add(place)
    for k in sorted(T): t = T[k]; print(f'{k}\t{t[0]}\t{t[1]}\t{t[1] / t[0]:.3f}\t{len(t[2])}\t{len(t[3])}')
    keep = {k: t for k, t in T.items() if 'repeat' not in k}; n = sum(t[0] for t in keep.values()); r = sum(t[0] if k[0] in '123' else t[1] if k[0] == '5' else 0 for k, t in keep.items())
    print(f'rows with a repeat share under 0.5: {n}; relatives (classes 1-3, or class 5 and similar): {r} ({r / max(1, n):.3f})')
    random.seed(1957); n = hit = 0
    while n < 5000:
        a, b = random.sample(genes, 2)
        if b in rel[a]: continue
        n += 1; hit += sim(a, b)
    print(f'random pairs of the same genes that are not listed relatives: {hit} of {n} similar ({hit / n:.4f})')

def selftest():
    assert union_len([(0, 10), (5, 20)]) == 20 and similar(300, 1000, 5000) and not similar(299, 1000, 5000) and not similar(150, 300, 300)
    print('pair_similarity selftest ok')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']: selftest()
    else: main(*sys.argv[1:])

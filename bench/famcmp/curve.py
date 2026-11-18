#!/usr/bin/env python3
"""curve.py REF OURS ST FLAIR: exact intron-chain precision/recall of each tool, and ours cut by reads / locus fraction."""
import sys, re, gzip, collections
def load(path, key=None):
    ex = collections.defaultdict(list); val = {}
    for l in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        f = l.split('\t')
        if len(f) < 9: continue
        if f[2] == 'exon': ex[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))].append((int(f[3]), int(f[4])))
        elif f[2] == 'transcript' and key:
            m = re.search(r'\b' + key + r' "([^"]+)"', f[8]); g = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
            val[(f[0], f[6], re.search(r'transcript_id "([^"]+)"', f[8]).group(1))] = (float(m.group(1)) if m else 0.0, g)
    if key and not val: sys.exit(f'{path}: no transcript line, cannot cut on {key}')
    out = {}
    for k, e in ex.items():
        e.sort()
        if len(e) > 1:
            ch = (k[0], k[1], tuple((e[i][1], e[i + 1][0] - 1) for i in range(len(e) - 1)))
            v = val.get(k, (0.0, '')); out[ch] = max(out.get(ch, v), v)
    return out
if __name__ == '__main__':
    ref = set(load(sys.argv[1])); ours = load(sys.argv[2], 'reads'); st = load(sys.argv[3], 'cov'); fl = load(sys.argv[4])
    rep = lambda n, s: print(f'{n}\t{len(s)}\t{len(s & ref)}\t{len(s & ref) / max(1, len(s)):.3f}')
    print('set\tchains\tcorrect\tprecision'); rep('ours', set(ours)); rep('stringtie', set(st)); rep('flair', set(fl))
    rep('ours&st', set(ours) & set(st)); rep('ours&flair', set(ours) & set(fl)); rep('st&flair', set(st) & set(fl)); rep('all3', set(ours) & set(st) & set(fl))
    rep('ours only', set(ours) - set(st) - set(fl)); rep('st only', set(st) - set(ours) - set(fl)); rep('flair only', set(fl) - set(ours) - set(st))
    for name, d in (('ours reads', ours), ('stringtie cov', st)):
        for t in (0, 2, 3, 4, 6, 11, 21):
            rep(f'{name} >= {t}', {c for c, v in d.items() if v[0] >= t})
    gmax = collections.Counter()
    for c, (r, g) in ours.items(): gmax[g] = max(gmax[g], r)
    for fr in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0):
        rep(f'ours frac of locus top >= {fr}', {c for c, (r, g) in ours.items() if r >= fr * gmax[g]})

#!/usr/bin/env python3
"""One-read loci by class against an annotation (register 1987, third addendum).   one_read_score.py ONELOCI.gtf TRUTH.gtf"""
import sys, os, re, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread

def main(gtf, truth):
    by = collections.defaultdict(list); ann = set()
    for c, s, ex, _ in oneread.models(truth).values():
        ch = oneread.chain(ex)
        for i in ch: by[(c, i)].append(ch); ann.add((c, i))
    G = collections.defaultdict(lambda: [0, 0, 0, 0])
    for c, s, ex, lines in oneread.models(gtf).values():
        if len(ex) < 2: continue
        ch = oneread.chain(ex); m = re.search(r'one_read_locus "([^"]+)"', lines[0]); r = re.search(r' reads "([0-9.]+)"', lines[0])
        k = ('one read, supported' if m.group(1) == 'supported' else 'one read, other') if m else 'no count' if not r else 'main, 2 reads' if float(r.group(1)) < 3 else 'main, 3 reads' if float(r.group(1)) < 4 else 'main, > 3 reads'
        cand = by.get((c, ch[0]), ()); g = G[k]; g[0] += 1; e = any(a == ch for a in cand); g[1] += e; g[2] += e or any(oneread.run(ch, a) for a in cand); g[3] += all((c, i) in ann for i in ch)
    print('class\tchains\texact\texact or contiguous part\tall introns annotated')
    for k in sorted(G): n, e, p, a = G[k]; print(k, n, f'{e / n:.3f}', f'{p / n:.3f}', f'{a / n:.3f}', sep='\t')

if __name__ == '__main__': main(*sys.argv[1:])

#!/usr/bin/env python3
"""PREREG r1629 (docs/PREREG_sibling_filter_2026-11-11.md): default 2-read chains at informed family members, split by
whether a sibling model projects to the same intron chain.

  sibling_filter.py REFSEQ_GTF[.gz] DEFAULT_GTF PROJECT_OUTDIR"""
import sys, re, gzip, collections

def models(path):
    ex = collections.defaultdict(list); reads = {}
    op = gzip.open if path.endswith('.gz') else open
    for l in op(path, 'rt'):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('exon', 'transcript'): continue
        k = (f[0], f[6], re.search(r'transcript_id "([^"]*)"', f[8]).group(1))
        if f[2] == 'exon': ex[k].append((int(f[3]) - 1, int(f[4])))
        else:
            m = re.search(r'reads "(\d+)"', f[8]); reads[k] = int(m.group(1)) if m else 0
    out = {}
    for k, e in ex.items():
        e.sort(); c = tuple((a[1], b[0]) for a, b in zip(e, e[1:]))
        if c:
            key = (k[0], k[1], c); out[key] = (max(out.get(key, (0,))[0], reads.get(k, 0)), e[0][0], e[-1][1])
    return out

def main():
    ref, default, d = sys.argv[1:4]
    truth = set(models(ref)); dflt = models(default)
    conc = set()
    for l in open(f'{d}/concordant.tsv'):
        c, st, ints = l.rstrip('\n').split('\t')
        conc.add((c, st, tuple(tuple(map(int, x.split('-'))) for x in ints.split(','))))
    inf = collections.defaultdict(list)
    for l in open(f'{d}/informed.tsv'):
        c, s, e, st = l.split(); inf[(c, st)].append((int(s), int(e)))
    at = {k: v for k, v in dflt.items() if any(v[1] < e and s < v[2] for s, e in inf.get((k[0], k[1]), ()))}
    two = {k for k, v in at.items() if v[0] == 2}
    D = {k for k in two if k not in conc}; K = two - D
    ok = lambda S: sum(k in truth for k in S)
    sh = lambda S: ok(S) / max(len(S), 1)
    print('set\tchains\tcorrect\tshare')
    for n, S in (('all default at informed members', set(at)), ('2-read', two), ('D dropped (no sibling)', D), ('K kept (sibling agrees)', K),
                 ('3+ reads, sibling agrees', {k for k in at if at[k][0] >= 3 and k in conc}), ('3+ reads, no sibling', {k for k in at if at[k][0] >= 3 and k not in conc})):
        print(f'{n}\t{len(S)}\t{ok(S)}\t{sh(S):.3f}')
    after = set(at) - D
    print(f"B1 share D {sh(D):.3f} <= 0.5 x 2-read {0.5 * sh(two):.3f}: {'PASS' if sh(D) <= 0.5 * sh(two) else 'FAIL'}")
    print(f"B2 lost {ok(D)} <= 2% of {ok(at)} = {0.02 * ok(at):.1f}: {'PASS' if ok(D) <= 0.02 * ok(at) else 'FAIL'}")
    print(f"B3 precision {sh(at):.4f} -> {sh(after):.4f} (+{sh(after) - sh(at):.4f}, bar +0.01): {'PASS' if sh(after) - sh(at) >= 0.01 else 'FAIL'}")

if __name__ == '__main__': main()

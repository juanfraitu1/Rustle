#!/usr/bin/env python3
"""PREREG r1791 (docs/PREREG_tie_split_2026-11-18.md): split the tied near-tie reads into SEVERAL COPIES / ONE PLACE.
Port of the workspace script of r1791 (r1796, docs/PREREG_famwise_origin_2026-11-18.md): it reads a truth table made with
or without NT_REPAIR (set the same value here), and a tie its own locus test cannot reproduce is written as UNSPLIT.
  tie_split.py READS.tsv READS.fa PLACEMENTS.fa MAT.mmi PAT.mmi READ_FAMILY.tsv OUT.tsv"""
import sys, os, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import neutral_truth as NT
rt, rfa, pfa, mat, pat, famtsv, out = sys.argv[1:8]
STRATA = os.environ.get('TS_STRATA', 'S2').split(',')  # r1829: S1,S2 (default unchanged: S2)
ties = {l.split('\t')[0] for l in open(rt) if l.split('\t')[1] in STRATA and l.split('\t')[3] == 'tie'}
tmp = os.environ.get('TMPDIR', '.')
def subset(src, dst):
    keep = False
    with open(dst, 'w') as w:
        for l in open(src):
            if l[0] == '>': keep = l[1:].rstrip('\n').split('|')[0] in ties
            if keep: w.write(l)
subset(rfa, f'{tmp}/r.fa'); subset(pfa, f'{tmp}/p.fa')
H = collections.defaultdict(list)
for hap, mmi in (('mat', mat), ('pat', pat)):
    for src in (f'{tmp}/r.fa', f'{tmp}/p.fa'):
        for k, v in NT.hap_hits(src, mmi, hap).items(): H[k] += v
plc = collections.defaultdict(list)
for k in H:
    if '|' in k: plc[k.rsplit('|', 1)[0]].append(k)
near_perfect = lambda h: h[4] <= max(2, 0.003 * h[5])
fam = dict(l.rstrip('\n').split('\t') for l in open(famtsv))
C = collections.defaultdict(collections.Counter)
with open(out, 'w') as w:
    w.write('read\tfamily\tclass\tbest_hits_mat\tbest_hits_pat\tmatched_mat\tmatched_pat\torigin_loci\n')
    for q in sorted(ties):
        f = fam.get(q, 'NA'); hs = H.get(q, []); best = min(hs, key=lambda h: (h[4], -h[5]), default=None)
        if best is None or not near_perfect(best): C[f]['UNSPLIT'] += 1; C['ALL']['UNSPLIT'] += 1; w.write(f"{q}\t{f}\tUNSPLIT\t0\t0\t0\t0\t\n"); continue
        org = [h for h in hs if h[4] == best[4] and h[5] >= 0.98 * best[5]]; origin = set(); matched = set()
        for k in plc[q]:  # the scorer's locus test, unchanged
            ch = H[k]; bh = {}
            for h in ch: bh[h[0]] = min(bh.get(h[0], (10 ** 9, 0)), (h[4], h[5]))
            loc = [h for h in ch if h[4] <= bh[h[0]][0] + 2 and bh[h[0]][0] <= max(2, 0.003 * bh[h[0]][1])]
            m = {j for j, b in enumerate(org) if any(NT.same_place(a, b) for a in loc)}
            if m: origin.add(k.rsplit('|', 1)[1]); matched |= m
        if len(origin) < 2:  # the truth's status is not reproduced by this test: not split
            C[f]['UNSPLIT'] += 1; C['ALL']['UNSPLIT'] += 1; w.write(f"{q}\t{f}\tUNSPLIT\t0\t0\t0\t0\t{','.join(sorted(origin))}\n"); continue
        D = {hp: sum(1 for j in matched if org[j][0] == hp) for hp in ('mat', 'pat')}
        N = {hp: sum(1 for h in org if h[0] == hp) for hp in ('mat', 'pat')}
        cls = 'SEVERAL COPIES' if max(D.values()) >= 2 else 'ONE PLACE'
        C[f][cls] += 1; C['ALL'][cls] += 1
        if cls == 'ONE PLACE' and any(N[hp] > max(D[hp], 1) for hp in N): C[f]['one place, further hits'] += 1; C['ALL']['one place, further hits'] += 1
        w.write(f"{q}\t{f}\t{cls}\t{N['mat']}\t{N['pat']}\t{D['mat']}\t{D['pat']}\t{','.join(sorted(origin))}\n")
print('family\ttrue ties\tSEVERAL COPIES\tONE PLACE\tshare several\tone place with further unreached best hits')
for f in sorted(C, key=lambda f: -sum(C[f][c] for c in ('SEVERAL COPIES', 'ONE PLACE'))):
    s, o = C[f]['SEVERAL COPIES'], C[f]['ONE PLACE']
    if s + o >= 10: print(f"{f}\t{s + o}\t{s}\t{o}\t{s / (s + o):.1%}\t{C[f]['one place, further hits']}")
print(f"UNSPLIT (not reproduced by the locus test): {C['ALL']['UNSPLIT']}")

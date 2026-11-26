#!/usr/bin/env python3
"""The members stage on a grid of its two settings: precision against reach (register 1967).

    members_grid.py ROWS.tsv RELATIONS.tsv MEMBERS.tsv

ROWS.tsv = pair_similarity.py's rows.tsv run with BAND all over place_classes.py's bands 0.3 0.5 0.7 0.9 (band, class,
place, source gene, genes under, repeat share, similar). Per threshold x repeat filter: rows with an annotated source, the
share at a relative (classes 1-3, or similar), the reach of the listed relatives of the sources that have any row, by
kind, and the size of the table (MEMBERS.tsv gives the distinct places and those on no assembled locus)."""
import sys, collections
KINDS = ('paralog_recent', 'paralog_older', 'pseudogene_of', 'same_source')

def grid(rows, rel):
    """rows [(band, class, place, src, under set, repeat share, similar)]; rel {gene: {gene: kind}} -> list of result dicts"""
    srcs = {r[3] for r in rows}; den = collections.Counter(); pairs = set()
    for s in srcs:
        for g, k in rel.get(s, {}).items(): pairs.add((s, g, k)); den[k] += 1
    out = []
    for thr in (0.3, 0.5, 0.7, 0.9):
        for filt in (False, True):
            keep = [r for r in rows if r[0] >= thr and (not filt or r[5] < 0.5)]; good = sum(1 for r in keep if r[1][0] in '123' or r[6])
            hit = {(r[3], g) for r in keep for g in r[4]}; reach = {k: sum(1 for s, g, kk in pairs if kk == k and (s, g) in hit) for k in KINDS}
            out.append(dict(thr=thr, filt=filt, rows=len(keep), places={r[2] for r in keep}, prec=good / max(1, len(keep)), reach=reach, den=den))
    return out

def main(rowsp, relp, mem):
    rel = collections.defaultdict(dict)
    for l in list(open(relp))[1:]:
        f = l.rstrip('\n').split('\t'); rel[f[0]][f[2]] = f[4]; rel[f[2]].setdefault(f[0], f[4])
    rows = []
    for l in open(rowsp):
        b, k, p, s, u, sh, sim = l.rstrip('\n').split('\t'); rows.append((float(b), k, p, s, set(u.split(',')) - {'.'}, float(sh), sim == '1'))
    novel = {l.split('\t')[0] for l in list(open(mem))[1:] if l.split('\t')[10] == '.'}
    R = grid(rows, rel); d = R[0]['den']
    print('listed relatives of the sources that have a row: ' + ', '.join(f'{k} {d[k]}' for k in KINDS))
    print('threshold\trepeat filter\trows (annotated source)\tdistinct places\ton no assembled locus\tat a relative\t' + '\t'.join(f'reach {k}' for k in KINDS))
    for r in R:
        print(r['thr'], 'on' if r['filt'] else 'off', r['rows'], len(r['places']), len(r['places'] & novel), f"{r['prec']:.3f}", *[f"{r['reach'][k]} ({r['reach'][k] / max(1, d[k]):.3f})" for k in KINDS], sep='\t')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        rel = {'A': {'B': 'paralog_recent', 'P': 'pseudogene_of'}}
        rows = [(0.9, '1 listed relative', 'p1', 'A', {'B'}, 0.0, True), (0.3, '1 listed relative', 'p2', 'A', {'P'}, 0.0, False), (0.9, '4 repeat', 'p3', 'A', {'X'}, 0.9, False)]
        R = grid(rows, rel); g = {(r['thr'], r['filt']): r for r in R}
        assert g[(0.3, False)]['rows'] == 3 and abs(g[(0.3, False)]['prec'] - 2 / 3) < 1e-9 and g[(0.3, False)]['reach']['pseudogene_of'] == 1
        assert g[(0.9, True)]['rows'] == 1 and g[(0.9, True)]['prec'] == 1.0 and g[(0.9, True)]['reach'] == {'paralog_recent': 1, 'paralog_older': 0, 'pseudogene_of': 0, 'same_source': 0}
        print('members_grid selftest ok')
    else: main(*sys.argv[1:])

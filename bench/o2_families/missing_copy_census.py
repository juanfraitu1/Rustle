#!/usr/bin/env python3
"""Where the reads of hidden copies land on a reference that lacks them (register 1982).

    missing_copy_census.py PLANT_DIR SCENARIO HIDDEN.tsv BAM OUT.tsv [RECURRENCE.tsv]

HIDDEN.tsv from hide_copy.py; BAM = the planted reads on the reduced reference. A read's landing = the planted copy
whose span holds most of its primary's aligned bases, else `elsewhere`, else `unmapped`. Differences per kb = NM over
the aligned reference bases. With RECURRENCE.tsv (recurrence_test.py on that run) a hidden copy is `seen` when a called
sequence has >= 5 reads and >= 80% of them are its own. OUT.tsv: one row per hidden copy; stdout: by edits class."""
import sys, collections, statistics, csv

def klass(edits, twin): return '0 twin' if twin or edits == 0 else '1-5' if edits <= 5 else '6-19' if edits < 20 else '>= 20' if edits >= 20 else 'no alignment'
def klass_of(edits, twin): return 'no alignment' if edits < 0 else klass(edits, twin)
def med(v): return statistics.median(v) if v else float('nan')

def seen(rec_rows, origin, min_reads=5, purity=0.8):
    """rec_rows [(reads count, [names])]; origin {read: copy} -> {copy: best own-read count of a pure sequence}"""
    out = collections.Counter()
    for n, names in rec_rows:
        c = collections.Counter(origin.get(x) for x in names); g, k = c.most_common(1)[0]
        if g and len(names) >= min_reads and k >= purity * len(names): out[g] = max(out[g], k)
    return out

def main(P, sc, hid, bam, out, rec=None):
    import pysam
    U = [l.rstrip('\n').split('\t') for l in open(f'{P}/units.tsv')][1:]; H = {r['copy']: r for r in csv.DictReader(open(hid), delimiter='\t')}
    spans = collections.defaultdict(list)
    for f in U:
        if f[1] not in H: spans[f[5]].append((int(f[6]) - 1, int(f[7]), f[1]))
    fam = {f[1]: f[0] for f in U}; name = {f[1]: f[2] for f in U}
    origin = {l.split('\t')[0]: l.split('\t')[1] for l in open(f'{P}/{sc}.truth.tsv')}
    land = collections.defaultdict(collections.Counter); d = collections.defaultdict(lambda: collections.defaultdict(list)); clip = collections.defaultdict(list); q0 = collections.Counter()
    for r in pysam.AlignmentFile(bam):
        if r.is_secondary or r.is_supplementary: continue
        g = origin[r.query_name]
        if r.is_unmapped: land[g]['unmapped'] += 1; continue
        ov = collections.Counter()
        for a, b in r.get_blocks():
            for s, e, c in spans.get(r.reference_name, ()):
                if s < b and a < e: ov[c] += min(b, e) - max(a, s)
        where = max(ov, key=lambda k: (ov[k], k)) if ov else 'elsewhere'; land[g][where] += 1
        d[where][g].append(1000 * r.get_tag('NM') / max(1, r.reference_length if not r.cigartuples else sum(n for o, n in r.cigartuples if o in (0, 2, 7, 8))))
        if g in H: clip[g].append(sum(n for o, n in r.cigartuples if o == 4) / max(1, r.query_length)); q0[g] += r.mapping_quality == 0
    S = collections.Counter()
    if rec:
        csv.field_size_limit(1 << 30)
        S = seen([(int(r['reads']), r['read_names'].split(',')) for r in csv.DictReader(open(rec), delimiter='\t')], origin)
    rows = []
    with open(out, 'w') as w:
        w.write('family\tcopy\tname\tclass\tedits\treads\tunmapped\telsewhere\ton_family\ttop_absorber\ttop_share\ttop_is_nearest\tdiff_per_kb\tabsorber_own_diff_per_kb\tratio\tclipped\tmapq0\tseen_reads\n')
        for g, h in sorted(H.items()):
            L = land[g]; n = sum(L.values())
            if not n: continue
            famc = {k: v for k, v in L.items() if fam.get(k) == h['family']}; placed = n - L['unmapped']
            top = max(famc, key=lambda k: (famc[k], k)) if famc else '.'
            mine = med(d[top][g]) if famc else float('nan'); own = med(d[top][top]) if famc else float('nan')
            r = dict(k=klass_of(int(h['edits']), h['twin'] == '1'), n=n, un=L['unmapped'], el=n - L['unmapped'] - sum(famc.values()), fam=sum(famc.values()), top=top, ts=famc.get(top, 0) / max(1, placed),
                     near=int(top == h['nearest']), mine=mine, own=own, ratio=mine / own if own and own == own else float('nan'), clip=med(clip[g]), q0=q0[g] / max(1, placed), seen=S[g], placed=placed)
            rows.append(r)
            w.write(f"{h['family']}\t{g}\t{h['name']}\t{r['k']}\t{h['edits']}\t{n}\t{r['un']}\t{r['el']}\t{r['fam']}\t{name.get(top, top)}\t{r['ts']:.3f}\t{r['near']}\t{mine:.2f}\t{own:.2f}\t{r['ratio']:.2f}\t{r['clip']:.3f}\t{r['q0']:.3f}\t{r['seen']}\n")
    print('class\thidden copies\treads\ton a family copy\tunmapped\telsewhere\ttop absorber share (median)\ttop = nearest\tdiff ratio (median)\tclipped (median)\tMAPQ 0 (median)\tseen (of copies with >= 20 placed reads)')
    for k in ('0 twin', '1-5', '6-19', '>= 20', 'no alignment', 'ALL'):
        v = [r for r in rows if k in (r['k'], 'ALL')]
        if not v: continue
        n = sum(r['n'] for r in v); big = [r for r in v if r['placed'] >= 20]
        print(k, len(v), n, f"{sum(r['fam'] for r in v) / n:.3f}", f"{sum(r['un'] for r in v) / n:.3f}", f"{sum(r['el'] for r in v) / n:.3f}", f"{med([r['ts'] for r in v]):.3f}", f"{sum(r['near'] for r in v)} of {len(v)}",
              f"{med([r['ratio'] for r in v if r['ratio'] == r['ratio']]):.2f}", f"{med([r['clip'] for r in v if r['clip'] == r['clip']]):.3f}", f"{med([r['q0'] for r in v]):.3f}", f"{sum(1 for r in big if r['seen'])} of {len(big)}" if rec else 'not run', sep='\t')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        assert klass_of(0, True) == '0 twin' and klass_of(3, False) == '1-5' and klass_of(19, False) == '6-19' and klass_of(20, False) == '>= 20' and klass_of(-1, False) == 'no alignment'
        o = {'a1': 'A', 'a2': 'A', 'a3': 'A', 'a4': 'A', 'a5': 'A', 'b1': 'B'}
        assert seen([(5, ['a1', 'a2', 'a3', 'a4', 'b1']), (4, ['a1', 'a2', 'a3', 'a4'])], o) == {'A': 4} and not seen([(5, ['a1', 'a2', 'a3', 'b1', 'b1'])], o)
        print('missing_copy_census selftest ok')
    else: main(*sys.argv[1:])

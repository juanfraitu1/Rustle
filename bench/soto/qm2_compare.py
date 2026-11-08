#!/usr/bin/env python3
"""Register row 1588 (docs/PREREG_soto_quickmer2_rerun_2026-11-08.md): compare our QuicK-mer2 run of Soto's released
pipeline (4 1KGP genomes) with Table S1E (2,504 genomes).

Per gene and sample: unweighted mean CN of the est windows (chrom start end CN) overlapping the gene's S1E `SD98_v1.0`
interval(s) (column 12, chr:start-end, comma-separated if several). Genes with no overlapping window are reported.
  Bar A  among parCN class = Fixed genes with a value, >= 90% of gene x sample values within 0.5 of S1E Median parCN
  Bar B  Spearman(our 4-sample median, S1E Median parCN) over all genes with a value >= 0.80
usage: qm2_compare.py S1E.tsv SAMPLE.CN.bed [SAMPLE.CN.bed ...] [--out per_gene.tsv]
       qm2_compare.py --selftest
"""
import bisect, collections, statistics, sys


def parse_ivs(s):
    out = []
    for part in s.split(','):
        part = part.strip()
        if ':' not in part:
            continue
        c, r = part.split(':')
        a, b = r.replace(',', '').split('-')
        out.append((c, int(a), int(b)))
    return out


def load_windows(path):
    w = collections.defaultdict(list)
    for line in open(path):
        f = line.split('\t')
        if len(f) >= 4 and not line.startswith(('#', 'track')):
            w[f[0]].append((int(f[1]), int(f[2]), float(f[3])))
    for c in w:
        w[c].sort()
    return w


def gene_cn(win, ivs):
    """unweighted mean CN of windows overlapping any interval; None if no window."""
    vals = []
    for c, a, b in ivs:
        ws = win.get(c, [])
        i = bisect.bisect_left(ws, (a, -1, 0.0))
        i = max(0, i - 1)                      # the window starting before a may overlap it
        while i < len(ws) and ws[i][0] < b:
            if ws[i][1] > a:
                vals.append(ws[i][2])
            i += 1
    return sum(vals) / len(vals) if vals else None


def ranks(x):
    order = sorted(range(len(x)), key=lambda i: x[i])
    r = [0.0] * len(x)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and x[order[j + 1]] == x[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.fmean(rx), statistics.fmean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float('nan')


def load_s1e(path):
    genes = []
    with open(path) as fh:
        cols = fh.readline().rstrip('\n').split('\t')
        ix = {c: i for i, c in reversed(list(enumerate(cols)))}   # first occurrence wins (S1E repeats some names)
        for line in fh:
            f = line.rstrip('\n').split('\t')
            try:
                med = float(f[ix['Median parCN']])
            except ValueError:
                continue
            genes.append((f[ix['Gene Name']], f[ix['Gene ID']], parse_ivs(f[ix['SD98_v1.0']]), med, f[ix['parCN class']]))
    return genes


def compare(genes, samples):
    """samples: [(name, windows)]. Returns (rows, summary)."""
    rows = []
    for name, gid, ivs, med, cls in genes:
        vals = [gene_cn(w, ivs) for _, w in samples]
        rows.append((name, gid, cls, med, vals))
    with_val = [r for r in rows if any(v is not None for v in r[4])]
    fixed = [(r[3], v) for r in with_val if r[2] == 'Fixed' for v in r[4] if v is not None]
    a_frac = sum(abs(v - m) <= 0.5 for m, v in fixed) / len(fixed) if fixed else float('nan')
    med_ours = [statistics.median([v for v in r[4] if v is not None]) for r in with_val]
    rho = spearman(med_ours, [r[3] for r in with_val]) if len(with_val) > 2 else float('nan')
    by_cls = {}
    for cls in ('Fixed', 'Nearly-Fixed', 'Polymorphic'):
        pv = [(r[3], v) for r in with_val if r[2] == cls for v in r[4] if v is not None]
        by_cls[cls] = (sum(abs(v - m) <= 0.5 for m, v in pv) / len(pv) if pv else float('nan'), len(pv))
    return rows, dict(genes=len(rows), with_value=len(with_val), no_window=len(rows) - len(with_val),
                      fixed_values=len(fixed), barA=a_frac, spearman=rho, by_class=by_cls)


def main(argv):
    out = None
    if '--out' in argv:
        i = argv.index('--out'); out = argv[i + 1]; argv = argv[:i] + argv[i + 2:]
    genes = load_s1e(argv[0])
    samples = [(p.rsplit('/', 1)[-1].split('.')[0], load_windows(p)) for p in argv[1:]]
    rows, s = compare(genes, samples)
    if out:
        with open(out, 'w') as fh:
            fh.write('gene\tgene_id\tclass\ts1e_median_parCN\t' + '\t'.join(n for n, _ in samples) + '\n')
            for name, gid, cls, med, vals in rows:
                fh.write(f"{name}\t{gid}\t{cls}\t{med}\t" + '\t'.join('NA' if v is None else f'{v:.3f}' for v in vals) + '\n')
    print(f"samples {[n for n, _ in samples]}")
    print(f"genes {s['genes']}  with a window {s['with_value']}  no window {s['no_window']}")
    for cls, (fr, n) in s['by_class'].items():
        print(f"  {cls:13s} within 0.5 of S1E median: {fr:.3f} of {n} gene x sample values")
    print(f"Bar A (Fixed >= 0.90): {s['barA']:.3f} -> {'PASS' if s['barA'] >= 0.90 else 'FAIL'}")
    print(f"Bar B (Spearman >= 0.80): {s['spearman']:.3f} -> {'PASS' if s['spearman'] >= 0.80 else 'FAIL'}")
    print('ROW 1588:', 'PASS' if s['barA'] >= 0.90 and s['spearman'] >= 0.80 else 'FAIL')


def selftest():
    assert parse_ivs('chr1:100-200') == [('chr1', 100, 200)]
    w = {'chr1': [(0, 100, 1.0), (100, 200, 2.0), (200, 300, 4.0)]}
    assert gene_cn(w, [('chr1', 150, 250)]) == 3.0           # windows 2 and 3
    assert gene_cn(w, [('chr1', 99, 101)]) == 1.5
    assert gene_cn(w, [('chr2', 0, 10)]) is None
    assert abs(spearman([1, 2, 3, 4], [10, 20, 30, 40]) - 1) < 1e-9
    assert abs(spearman([1, 2, 3], [3, 2, 1]) + 1) < 1e-9
    genes = [('G1', 'i1', [('chr1', 0, 100)], 1.0, 'Fixed'), ('G2', 'i2', [('chr1', 100, 200)], 2.0, 'Fixed'),
             ('G3', 'i3', [('chr1', 200, 300)], 4.0, 'Polymorphic'), ('G4', 'i4', [('chr9', 0, 1)], 2.0, 'Fixed')]
    _, s = compare(genes, [('S', w)])
    assert (s['with_value'], s['no_window'], s['fixed_values'], s['barA']) == (3, 1, 2, 1.0), s
    assert abs(s['spearman'] - 1) < 1e-9
    print('selftest ok')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(sys.argv[1:])

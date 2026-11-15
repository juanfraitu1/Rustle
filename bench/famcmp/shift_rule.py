#!/usr/bin/env python3
"""PREREG row 1726: graded drop rule for weaker shifted-site chains (class table of shift_class.py).

  shift_rule.py CLASS_TSV                    all 24 cells
  shift_rule.py CLASS_TSV FRAME Q SHIFT      one cell (FRAME breaking|any, SHIFT bp or any)   |   --selftest
  shift_rule.py CLASS_TSV FRAME Q SHIFT IN_GTF OUT_GTF    also write IN_GTF without the dropped chains' transcripts"""
import sys, csv, re, collections

def cell(rows, frame, q, shift):
    """-> indices dropped: class S1 / S2, reads ratio < q, shift <= shift (None = any), frame-breaking only if asked."""
    return [i for i, r in enumerate(rows) if r[0] in ('S1', 'S2') and r[2] < q * r[3] and (shift is None or r[1] <= shift) and (frame == 'any' or r[1] % 3)]

def score(rows, drop):
    n = len(drop); N = len(rows); ok = sum(r[4] for r in rows); d_ok = sum(rows[i][4] for i in drop)
    by = sorted(range(N), key=lambda i: rows[i][2]); t = rows[by[n - 1]][2] if n else 0; lo = [i for i in by if rows[i][2] < t]; tie = [i for i in by if rows[i][2] == t]
    c_ok = sum(rows[i][4] for i in lo) + (n - len(lo)) * sum(rows[i][4] for i in tie) / max(1, len(tie))  # tied read count: its genome-wide share, not table order
    pairs = [i for i, r in enumerate(rows) if r[0] in ('S1', 'S2') and r[4] and r[6]]; ds = set(drop)
    return dict(dropped=n, in_ref=d_ok / max(1, n), cut_in_ref=c_ok / max(1, n), pairs_lost=sum(i in ds for i in pairs) / max(1, len(pairs)), correct_lost=d_ok / ok,
                intron_ann=sum(rows[i][5] for i in drop) / max(1, n), prec=ok / N, prec_rule=(ok - d_ok) / (N - n), prec_cut=(ok - c_ok) / (N - n))

def read(p):
    return [(r['class'], int(r['shift']), int(r['reads']), int(r['partner_reads']), int(r['in_ref']), int(r['intron_annotated']), int(r['partner_in_ref'])) for r in csv.DictReader(open(p), delimiter='\t')]

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        rows = [('S1', 4, 1, 100, 0, 0, 1), ('S1', 3, 1, 100, 1, 1, 1), ('S1t', 4, 1, 100, 0, 0, 1), ('-', 0, 2, 0, 1, 0, 0), ('S2', 40, 50, 100, 1, 1, 1)]
        assert cell(rows, 'breaking', 0.1, 10) == [0] and cell(rows, 'any', 0.1, None) == [0, 1] and cell(rows, 'any', 0.6, None) == [0, 1, 4]
        s = score(rows, [0, 1]); assert s['dropped'] == 2 and s['in_ref'] == 0.5 and abs(s['cut_in_ref'] - 1 / 3) < 1e-9 and s['pairs_lost'] == 0.5 and abs(s['correct_lost'] - 1 / 3) < 1e-9
        print('selftest ok'); sys.exit()
    rows = read(sys.argv[1]); cells = [(f, q, s) for f in ('breaking', 'any') for q in (0.02, 0.05, 0.1, 0.2) for s in (10, 50, None)]
    if len(sys.argv) > 2: cells = [(sys.argv[2], float(sys.argv[3]), None if sys.argv[4] == 'any' else int(sys.argv[4]))]
    print('frame\tQ\tshift\tdropped\tin_RefSeq\tread_cut_in_RefSeq\tratio\tpairs_lost\tcorrect_lost\tintron_annotated\tprecision\tafter_rule\tafter_read_cut')
    for f, q, s in cells:
        r = score(rows, cell(rows, f, q, s))
        print(f, q, s or 'any', r['dropped'], *(f'{r[k]:.3f}' for k in ('in_ref', 'cut_in_ref')), f"{r['in_ref'] / max(1e-9, r['cut_in_ref']):.2f}", *(f'{r[k]:.4f}' for k in ('pairs_lost', 'correct_lost')), f"{r['intron_ann']:.3f}", *(f'{r[k]:.4f}' for k in ('prec', 'prec_rule', 'prec_cut')), sep='\t')
    if len(sys.argv) > 6:
        csv.field_size_limit(1 << 30); meta = [(r['contig'], r['strand'], r['chain']) for r in csv.DictReader(open(sys.argv[1]), delimiter='\t')]
        gone = {meta[i] for i in cell(rows, *cells[0])}; ex = collections.defaultdict(list); T = re.compile(r'transcript_id "([^"]+)"'); L = open(sys.argv[5]).readlines()
        for l in L:
            f = l.split('\t')
            if len(f) > 8 and f[2] == 'exon': ex[f[0], f[6], T.search(f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
        bad = {k[2] for k, e in ex.items() if (e.sort() or True) and (k[0], k[1], str(tuple((a[1], b[0]) for a, b in zip(e, e[1:])))) in gone}
        with open(sys.argv[6], 'w') as w: w.writelines(l for l in L if not ((m := T.search(l)) and m.group(1) in bad))
        print(f'wrote {sys.argv[6]}: {len(bad)} transcripts of {len(gone)} chains removed', file=sys.stderr)

#!/usr/bin/env python3
"""Row 1662: the complete copy table (tsv + fa + regions) for families with extra locus representatives attached.

  reps_table.py A_PREFIX C_PREFIX CM_COPIES.tsv ATTACHED.tsv OUT_PREFIX
  reps_table.py --selftest
A_PREFIX / C_PREFIX = the two families runs' `X.fam.copies` (reads .tsv and .fa); CM_COPIES = C after piece_merge.py;
ATTACHED = piece_attach.py's table. Writes OUT_PREFIX.tsv/.fa/.regions as mcl_families does: a family the attach step
changed is renumbered 0..n-1 (its kept rows in their order, then the attached representatives by position) and its FASTA
records come from the run that made each row; an unchanged family is copied as it is. Regions = per family and contig, the
hull of the copies' start..end, 5 kb each side (checked equal to mcl_families' file on the shipped A119b table).
"""
import sys, collections

PAD = 5000


def read_fa(path):
    out, k = {}, None
    for l in open(path):
        if l[0] == '>':
            k = tuple(l[1:].split('|')[:2]); out[k] = []
        else:
            out[k].append(l)
    return {k: ''.join(v) for k, v in out.items()}


def build(hdr, a_lines, cm_lines, att_lines):
    """lines: copy-table rows (no header). -> [(out_row_fields, source ('A' | 'C', family, idx))] in output order."""
    ix = {c: i for i, c in enumerate(hdr)}
    f = lambda l: l.rstrip('\n').split('\t')
    a_set = set(a_lines)
    c_of = {f(l)[ix['tid']]: (f(l)[0], f(l)[1]) for l in cm_lines}
    order = list(dict.fromkeys(f(l)[0] for l in a_lines))
    kept, added = collections.defaultdict(list), collections.defaultdict(list)
    for l in att_lines:
        r = f(l)
        (kept if l in a_set else added)[r[0]].append(r)
    n_a = collections.Counter(f(l)[0] for l in a_lines)
    have = {r[ix['tid']] for v in kept.values() for r in v}
    added = collections.defaultdict(list, {k: [r for r in v if r[ix['tid']] not in have] for k, v in added.items()})
    out = []
    for fam in order:
        changed = bool(added[fam]) or len(kept[fam]) != n_a[fam]
        rows = [(r, ('A', fam, r[1])) for r in kept[fam]]
        rows += [(r, ('C',) + c_of[r[ix['tid']]]) for r in sorted(added[fam], key=lambda r: (r[ix['chrom']], int(r[ix['start']]), r[ix['tid']]))]
        for i, (r, src) in enumerate(rows):
            out.append((r[:1] + [str(i)] + r[2:] if changed else r, src))
    return out


def regions(hdr, rows):
    ix = {c: i for i, c in enumerate(hdr)}; hull = {}
    for r in rows:
        k = (r[0], r[ix['chrom']]); s, e = int(r[ix['start']]), int(r[ix['end']])
        hull[k] = (min(s, hull[k][0]), max(e, hull[k][1])) if k in hull else (s, e)
    return [f'{fam}\t{c}:{max(0, s - PAD)}-{e + PAD}\n' for (fam, c), (s, e) in hull.items()]


def main(a_prefix, c_prefix, cm, att, out):
    lines = lambda p: open(p).read().splitlines(keepends=True)
    a, attl = lines(a_prefix + '.tsv'), lines(att)
    hdr = a[0].rstrip('\n').split('\t'); ix = {c: i for i, c in enumerate(hdr)}
    rows = build(hdr, a[1:], lines(cm)[1:], attl[1:])
    assert len(rows) <= len(attl) - 1 and len({r[ix['tid']] for r, _ in rows}) == len(rows), 'a transcript is a copy of two families'
    fa = {'A': read_fa(a_prefix + '.fa'), 'C': read_fa(c_prefix + '.fa')}
    with open(out + '.tsv', 'w') as t, open(out + '.fa', 'w') as q:
        t.write(a[0])
        for r, (run, fam, idx) in rows:
            t.write('\t'.join(r) + '\n')
            q.write(f">{r[0]}|{r[1]}|{r[ix['chrom']]}:{r[ix['start']]}-{r[ix['end']]}|{r[ix['strand']]}|nexon={r[ix['n_exon']]}\n" + fa[run][(fam, idx)])
    open(out + '.regions', 'w').writelines(regions(hdr, [r for r, _ in rows]))
    print(f'[reps_table] {len(a) - 1} -> {len(rows)} copies, {sum(1 for _, s in rows if s[0] == "C")} attached representatives -> {out}.tsv/.fa/.regions')


def selftest():
    hdr = ['family_id', 'copy_idx', 'tid', 'chrom', 'start', 'end', 'n_exon', 'strand']
    mk = lambda *v: '\t'.join(map(str, v)) + '\n'
    a = [mk('F1', 0, 'L', 'c', 100, 200, 1, '+'), mk('F1', 1, 'M', 'c', 9000, 9100, 1, '+'), mk('F2', 0, 'Q', 'c', 50, 60, 1, '+'), mk('F2', 1, 'R', 'd', 5, 6, 1, '+')]
    cm = [mk('C7', 3, 'Lp', 'c', 300, 400, 2, '-'), mk('C7', 4, 'M', 'c', 9000, 9100, 1, '+')]
    att = [a[1], a[2], a[3], mk('F1', 3, 'Lp', 'c', 300, 400, 2, '-')]          # L replaced by its piece Lp; F2 untouched
    rows = build(hdr, a, cm, att)
    assert [(r[0], r[1], r[2], s) for r, s in rows] == [('F1', '0', 'M', ('A', 'F1', '1')), ('F1', '1', 'Lp', ('C', 'C7', '3')),
                                                        ('F2', '0', 'Q', ('A', 'F2', '0')), ('F2', '1', 'R', ('A', 'F2', '1'))], rows
    assert regions(hdr, [r for r, _ in rows]) == ['F1\tc:0-14100\n', 'F2\tc:0-5060\n', 'F2\td:0-5006\n']
    # a representative whose transcript is already a kept copy of ANOTHER family is not attached (one copy, one family)
    att2 = att + [mk('F2', 5, 'M', 'c', 9000, 9100, 1, '+')]
    assert [(r[0], r[2]) for r, _ in build(hdr, a, cm, att2)] == [('F1', 'M'), ('F1', 'Lp'), ('F2', 'Q'), ('F2', 'R')]
    print('selftest ok')


if __name__ == '__main__':
    selftest() if '--selftest' in sys.argv else main(*sys.argv[1:6])

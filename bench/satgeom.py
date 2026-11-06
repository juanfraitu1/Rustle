#!/usr/bin/env python3
"""Register row 1561: flag satellite-like families by copy geometry alone, score against a truth list.
Rule (docs/PREREG_satgeom_2026-11-06.md, fixed): >= 6 copies, one chromosome, span <= 5 Mb, median copy length <= 2 kb.
usage: satgeom.py COPIES_TSV TRUTH_TXT   (truth = tools/annotated_families.py COPIES censat.bed --majority-in)"""
import collections, statistics, sys

MIN_COPIES, MAX_SPAN, MAX_MEDIAN_LEN = 6, 5_000_000, 2_000


def failed_clauses(copies):
    """copies = [(chrom, start, end)]; returns the rule clauses the family fails (empty = satellite-like)."""
    fail = []
    if len(copies) < MIN_COPIES:
        fail.append('copies<6')
    if len({c for c, _, _ in copies}) > 1:
        fail.append('multi_chrom')
    elif max(e for _, _, e in copies) - min(s for _, s, _ in copies) > MAX_SPAN:
        fail.append('span>5Mb')
    if statistics.median(e - s for _, s, e in copies) > MAX_MEDIAN_LEN:
        fail.append('median_len>2kb')
    return fail


def main(copies_tsv, truth_txt):
    fams, genes = collections.defaultdict(list), collections.defaultdict(collections.Counter)
    with open(copies_tsv) as fh:
        col = {k: i for i, k in enumerate(fh.readline().rstrip('\n').split('\t'))}
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < len(col):
                continue
            fid = f[col['family_id']]
            fams[fid].append((f[col['chrom']], int(f[col['start']]), int(f[col['end']])))
            genes[fid][f[col['gene_id']]] += 1
    truth = {l.strip() for l in open(truth_txt) if l.strip()}
    assert truth <= fams.keys(), 'truth lists families absent from the copies table'
    fail = {fid: failed_clauses(c) for fid, c in fams.items()}
    flagged = {fid for fid, x in fail.items() if not x}
    tp, fp, fn = flagged & truth, flagged - truth, truth - flagged
    p, r = len(tp) / max(len(flagged), 1), len(tp) / len(truth)
    print(f"families\t{len(fams)}\ntruth\t{len(truth)}\nflagged\t{len(flagged)}\nTP\t{len(tp)}\nFP\t{len(fp)}\nFN\t{len(fn)}")
    print(f"precision\t{p:.3f}\t{'PASS' if p >= 0.9 else 'FAIL'} (bar 0.9)\nrecall\t{r:.3f}\t{'PASS' if r >= 0.8 else 'FAIL'} (bar 0.8)")
    print('FN by failing clause (a family may fail several):')
    for k, v in collections.Counter(x for fid in fn for x in fail[fid]).most_common():
        print(f"  {k}\t{v}")
    print('FP families (family, copies, chrom, span, median_len, top gene_ids):')
    for fid in sorted(fp, key=lambda x: -len(fams[x])):
        c = fams[fid]
        print(f"  {fid}\t{len(c)}\t{c[0][0]}\t{max(e for _, _, e in c) - min(s for _, s, _ in c)}\t"
              f"{statistics.median(e - s for _, s, e in c):.0f}\t{','.join(g for g, _ in genes[fid].most_common(3))}")


# self-check of the rule on toy families
assert failed_clauses([('chr1', i * 1000, i * 1000 + 500) for i in range(6)]) == []
assert failed_clauses([('chr1', 0, 500)] * 5) == ['copies<6']
assert failed_clauses([('chr1', 0, 500)] * 5 + [('chr2', 0, 500)]) == ['multi_chrom']
assert failed_clauses([('chr1', 0, 3000)] * 5 + [('chr1', 6_000_000, 6_003_000)]) == ['span>5Mb', 'median_len>2kb']

if __name__ == '__main__':
    main(*sys.argv[1:3])

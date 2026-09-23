#!/usr/bin/env python3
"""Does block colinearity separate real containments from repeat/domain ones?
Per `docs/PREREG_colinearity_conjunct_2026-09-21.md`.

r919 swept containment/identity from 0.30 to 0.99 on the 76 pairs the shipped rule rejects at
containment >= 0.90 (19 TRUE / 57 FALSE, held-out chr2/chr8/chr10) and found a FLAT precision curve
(~0.25) -- no pairwise scalar separates them. Every earlier pairwise script keeps only the BEST PAF
record per gene pair; 95.6% of chr2 gene pairs have more than one record (minimap2 run with -P), so
this is genuinely unused information.

    colinearity(pair) = Kendall-tau concordance of block order between A's and B's coordinates,
                         over every surviving PAF record for that pair (not just the best one)

Usage: colinearity_conjunct.py --gff G --soto S1C --pafs chr2=P1,chr8=P2,chr10=P3 --graphs DIR
"""
import argparse
import collections
import csv
import gzip
import itertools
import re
import statistics

MIN_BLOCK_NMATCH = 30


def _open(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def gene_names(gff, chrom):
    n = {}
    for line in _open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        m = re.search(r'Name=([^;]+)', f[8])
        if m:
            n[f'{f[0]}:{f[3]}-{f[4]}'] = m.group(1)
    return n


def soto(s1c):
    ids = collections.defaultdict(set)
    for r in csv.DictReader(open(s1c), delimiter='\t'):
        fid = (r.get('Family ID') or '').strip(); nm = (r.get('Gene Name') or '').strip()
        if fid and fid != 'N/A' and nm:
            ids[nm].add(fid)
    return {nm: next(iter(v)) for nm, v in ids.items() if len(v) == 1}


def load_all_blocks(paf):
    """(a,b) sorted tuple -> list of (a_start,a_end,b_start,b_end,strand) in a's/b's own coords,
    plus the best (m, la, lb) per pair for the containment gate (matches vg_multiplicity_containment.py)."""
    blocks = collections.defaultdict(list)
    best = {}
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        a, b = f[0], f[5]
        if a == b:
            continue
        la, lb, m = int(f[1]), int(f[6]), int(f[9])
        qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
        strand = f[4]
        k = tuple(sorted((a, b)))
        if m >= MIN_BLOCK_NMATCH:
            if k[0] == a:
                blocks[k].append((qs, qe, ts, te, strand))
            else:
                blocks[k].append((ts, te, qs, qe, strand))
        rec = (m, la, lb, qs, qe, ts, te) if k[0] == a else (m, lb, la, ts, te, qs, qe)
        if k not in best or m > best[k][0]:
            best[k] = rec
    return best, blocks


def colinearity(pairblocks):
    """Kendall-tau concordance of block order; None if <2 blocks or mixed strand."""
    if len(pairblocks) < 2:
        return None, len(pairblocks)
    strands = {s for (_, _, _, _, s) in pairblocks}
    if len(strands) > 1:
        return 0.0, len(pairblocks)
    strand = strands.pop()
    ordered = sorted(pairblocks, key=lambda x: x[0])
    bstarts = [x[2] for x in ordered]
    pairs = list(itertools.combinations(bstarts, 2))
    if not pairs:
        return None, len(pairblocks)
    if strand == '+':
        concordant = sum(1 for x, y in pairs if x <= y)
    else:
        concordant = sum(1 for x, y in pairs if x >= y)
    return concordant / len(pairs), len(pairblocks)


def main():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--soto', '--pafs', '--graphs'):
        ap.add_argument(x, required=True)
    a = ap.parse_args()
    fam = soto(a.soto)

    rows = []
    for spec in a.pafs.split(','):
        chrom, paf = spec.split('=')
        names = gene_names(a.gff, chrom)
        shipped = set()
        for line in open(f'{a.graphs}/{chrom}.graph.tsv'):
            f = line.rstrip('\n').split('\t')
            if len(f) == 3 and f[0] != f[1]:
                shipped.add(tuple(sorted((f[0], f[1]))))
        best, blocks = load_all_blocks(paf)
        for k, (m, la, lb, qs, qe, ts, te) in best.items():
            if k in shipped:
                continue
            ga, gb = names.get(k[0]), names.get(k[1])
            if not ga or not gb:
                continue
            fa, fb = fam.get(ga), fam.get(gb)
            if not fa or not fb:
                continue
            lo = min(la, lb)
            if not lo or m / lo < 0.90:
                continue
            score, nblocks = colinearity(blocks.get(k, []))
            rows.append((fa == fb, score, nblocks, ga, gb))

    T = [r for r in rows if r[0]]
    F = [r for r in rows if not r[0]]
    print(f"population: {len(rows)} rejected containment>=0.90 pairs — {len(T)} TRUE, {len(F)} FALSE\n")

    Tb = [r[2] for r in T]; Fb = [r[2] for r in F]
    print(f"  block count   TRUE median {statistics.median(Tb) if Tb else 0:.1f}"
          f"   FALSE median {statistics.median(Fb) if Fb else 0:.1f}")
    Tunder = sum(1 for r in T if r[2] < 2); Funder = sum(1 for r in F if r[2] < 2)
    print(f"  <2 blocks (colinearity undefined): TRUE {Tunder}/{len(T)} ({100*Tunder/len(T) if T else 0:.1f}%)"
          f"   FALSE {Funder}/{len(F)} ({100*Funder/len(F) if F else 0:.1f}%)\n")

    Ts = [r[1] for r in T if r[1] is not None]
    Fs = [r[1] for r in F if r[1] is not None]
    print(f"  colinearity (n={len(Ts)+len(Fs)} scoreable)  TRUE median "
          f"{statistics.median(Ts) if Ts else float('nan'):.3f}   FALSE median "
          f"{statistics.median(Fs) if Fs else float('nan'):.3f}\n")

    print(f"  {'score >= t':>10} {'TRUE kept':>10} {'FALSE kept':>11} {'precision':>10} {'recall':>8}")
    best_p = 0.0
    for t in (0.50, 0.70, 0.80, 0.90, 0.95, 0.99, 1.00):
        t_kept = sum(1 for r in T if r[1] is not None and r[1] >= t)
        f_kept = sum(1 for r in F if r[1] is not None and r[1] >= t)
        p = t_kept / (t_kept + f_kept) if (t_kept + f_kept) else 0.0
        if t_kept >= 10:
            best_p = max(best_p, p)
        print(f"  {t:>10.2f} {t_kept:>10} {f_kept:>11} {p:>10.3f} {t_kept/len(T) if T else 0:>8.3f}")

    pairs = list(itertools.product(Ts, Fs))
    wins = sum(1 for x, y in pairs if x > y) + 0.5 * sum(1 for x, y in pairs if x == y)
    auc = wins / len(pairs) if pairs else float('nan')
    print(f"\n  AUC (colinearity separating TRUE from FALSE, scoreable only) = {auc:.3f}"
          f"   [r906's multiplicity AUC was 0.681]")
    print(f"  baseline precision (containment alone) = 0.250")
    print(f"  best precision at >=10 TRUE retained (scoreable-only threshold sweep) = {best_p:.3f}")


if __name__ == '__main__':
    main()

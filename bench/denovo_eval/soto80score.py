#!/usr/bin/env python3
"""Score family catalogs against a Soto 2025 family panel BED (chrom start end GENE|ID_n ...; e.g. 80_fams.chr.bed, the
eligible members of 83 Soto S1C families, CHM13 v2.0). Descriptive, advisor-facing (register 1580 amendment 1).

A copy (PREFIX.fam.copies.tsv row) is placed on the panel gene whose interval holds >= half of its exon bases (BED has no
strand, so strand is ignored; ties -> the larger overlap, then the gene name). Truth pairs = panel genes sharing a Soto ID.
Predicted pairs = panel genes with copies in one of our families. Reported per catalog:
  precision  hit / predicted (only pairs of two panel genes are judged)
  recall     hit / all truth pairs;  recall_reach = hit / truth pairs whose two genes both carry a copy in that catalog
  families_complete  Soto families (>= 2 genes in the panel) whose reachable genes all share one of our families

usage: soto80score.py PANEL.bed COPIES.tsv [COPIES.tsv ...]
       soto80score.py --selftest
"""
import collections, itertools, sys


def load_panel(path):
    genes = []
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 4 or not f[3]:
            continue
        name, fam = f[3].rsplit('|', 1)
        genes.append((f[0], int(f[1]), int(f[2]), name, fam))
    return genes


def place(genes, chrom, exons):
    """exons: [(s, e)] 0-based half-open. Returns the panel gene holding >= half of the exon bases, or None."""
    tot = sum(e - s for s, e in exons)
    best = None
    for c, gs, ge, name, _ in genes:
        if c != chrom:
            continue
        ov = sum(max(0, min(e, ge) - max(s, gs)) for s, e in exons)
        if tot and 2 * ov >= tot and (best is None or (ov, name) > best):
            best = (ov, name)
    return best[1] if best else None


def score(genes, rows):
    """rows: [(family_id, chrom, exons)]. Returns a dict of counts."""
    fam_of = {name: fam for *_, name, fam in genes}
    by_chrom = collections.defaultdict(list)
    for g in genes:
        by_chrom[g[0]].append(g)
    ours = collections.defaultdict(set)          # our family -> panel genes
    reach = set()
    for fid, chrom, exons in rows:
        g = place(by_chrom.get(chrom, []), chrom, exons)
        if g:
            ours[fid].add(g)
            reach.add(g)
    truth = {frozenset(p) for p in itertools.combinations(sorted(fam_of), 2) if fam_of[p[0]] == fam_of[p[1]]}
    pred = {frozenset(p) for gs in ours.values() for p in itertools.combinations(sorted(gs), 2)}
    hit = pred & truth
    truth_reach = {p for p in truth if p <= reach}
    members = collections.defaultdict(set)
    for name, fam in fam_of.items():
        members[fam].add(name)
    complete = judged = 0
    for fam, ms in members.items():
        r = ms & reach
        if len(ms) < 2 or len(r) < 2:
            continue
        judged += 1
        complete += any(r <= gs for gs in ours.values())
    return dict(pred=len(pred), hit=len(hit), truth=len(truth), truth_reach=len(truth_reach), reach=len(reach),
                families_judged=judged, families_complete=complete)


def read_copies(path):
    rows = []
    with open(path) as fh:
        cols = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            r = dict(zip(cols, line.rstrip('\n').split('\t')))
            exons = [(int(a), int(b)) for a, b in (x.split('-') for x in r['exons'].split(','))]
            rows.append((r['family_id'], r['chrom'], exons))
    return rows


def main(panel, paths):
    genes = load_panel(panel)
    print(f"# panel {panel}: {len(genes)} genes, {len({g[4] for g in genes})} Soto families")
    print("catalog\tpred\thit\tprecision\trecall\trecall_reach\tgenes_reached\tfamilies_complete")
    for p in paths:
        s = score(genes, read_copies(p))
        print(f"{p}\t{s['pred']}\t{s['hit']}\t{s['hit'] / max(1, s['pred']):.3f}\t{s['hit'] / max(1, s['truth']):.3f}\t"
              f"{s['hit'] / max(1, s['truth_reach']):.3f}\t{s['reach']}\t{s['families_complete']}/{s['families_judged']}")


def selftest():
    genes = [('chr1', 0, 100, 'A1', 'ID_1'), ('chr1', 200, 300, 'A2', 'ID_1'), ('chr1', 400, 500, 'A3', 'ID_1'),
             ('chr2', 0, 100, 'B1', 'ID_2'), ('chr2', 200, 300, 'B2', 'ID_2')]
    assert place(genes, 'chr1', [(10, 20), (150, 160)]) == 'A1'          # half the exon bases inside A1
    assert place(genes, 'chr1', [(90, 100), (150, 170)]) is None          # only a third inside
    # our F1 joins A1+A2 (A3 never seen), F2 joins B1 with A3's neighbour nothing, F3 = B2 alone
    rows = [('F1', 'chr1', [(10, 50)]), ('F1', 'chr1', [(210, 250)]), ('F2', 'chr2', [(10, 50)]), ('F3', 'chr2', [(210, 250)])]
    s = score(genes, rows)
    assert (s['truth'], s['pred'], s['hit'], s['truth_reach'], s['reach']) == (4, 1, 1, 2, 4), s
    assert (s['families_judged'], s['families_complete']) == (2, 1), s     # ID_1 reachable A1+A2 together; ID_2 split
    print('selftest ok')


if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        selftest()
    else:
        main(sys.argv[1], sys.argv[2:])

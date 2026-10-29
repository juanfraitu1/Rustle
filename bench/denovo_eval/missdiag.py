"""Why did a Compara true pair miss? Per missed pair, classify the gated all-vs-all records between the two genes' loci.

  missdiag.py PREFIX COPIES_TSV OUT_TSV [--levels recent]

PREFIX = mcl_families --out prefix (PREFIX.loci.gff3, PREFIX.loci.gated.paf). Reasons, first that applies:
  no_locus    a gene has no locus of its own (merged into another gene's locus)
  no_records  no PAF record joins any locus of A to any locus of B
  short       every record < 300 bp (min_bp)
  low_id      every >= 300 bp record < 0.70 identity
  low_shared  best single record's exon-to-exon bases < 0.60 of the smaller locus's exonic length
  mcl_split   (with --graph) an edge joins the two genes' loci: lost in MCL
  low_cov_or_mcl  passes the exon gate on the best record, yet no edge (coverage) or, without --graph, an MCL split
Also reports, per pair, the UNION of exon-to-exon bases over all >= 0.70 identity records of any length (the retrocopy case).
"""
import argparse
import collections
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F  # noqa: E402


def merged(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def exonic(blocks, s, e):
    return sum(max(0, min(b, e) - max(a, s)) for a, b in blocks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('prefix')
    ap.add_argument('copies')
    ap.add_argument('out')
    ap.add_argument('--levels', default='recent')
    ap.add_argument('--graph', help='mcl_families --dump-graph TSV: a pair with an edge between its loci is mcl_split')
    ap.add_argument('--species', default='human', choices=['human', 'gorilla'], help='gorilla: RefSeq genes by symbol (famtruth.load_ggo)')
    a = ap.parse_args()
    rep = {}
    for line in open(a.prefix + '.loci.tsv'):
        x, y = line.rstrip('\n').split('\t')[:2]
        rep[x] = y
    edges = set()
    if a.graph:
        for line in open(a.graph):
            x, y, _ = line.split('\t')
            edges.add((x, y))
            edges.add((y, x))
    G = F.load_cat()
    if a.species == 'gorilla':
        G = F.load_ggo(G)
    idx = F.exon_index(G)
    loci = {}  # name -> (chrom, start0, strand, local exon blocks)
    cur = None
    for line in open(a.prefix + '.loci.gff3'):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        if f[2] == 'gene':
            cur = f'{f[0]}:{f[3]}-{f[4]}'
            loci[cur] = [f[0], int(f[3]) - 1, f[6] if f[6] in '+-' else '+', []]
        elif f[2] == 'exon':
            loci[cur][3].append((int(f[3]) - 1, int(f[4])))
    gene_loci = collections.defaultdict(list)
    for n, (c, s0, st, ex) in loci.items():
        g = F.map_copy(idx, c, st, sorted(ex))
        if g:
            gene_loci[g].append(n)
        loci[n][3] = [(x - s0, y - s0) for x, y in sorted(ex)]
    fams, _, _ = F.load_catalog(a.copies, G, idx)
    fam_of = {g: k for k, v in fams.items() for g in v}
    in_fam = set(fam_of)
    pairs, universe = F.truth_pairs(G, F.PARALOGS, F.LEVELS[a.levels])
    missed = [p for p in pairs if p[0] in in_fam | set(gene_loci) and p[1] in in_fam | set(gene_loci)
              and p[0] in universe and p[1] in universe and (fam_of.get(p[0]) is None or fam_of.get(p[0]) != fam_of.get(p[1]))]
    want = {}
    for x, y in missed:
        for la in gene_loci.get(x, []):
            for lb in gene_loci.get(y, []):
                want[(la, lb)] = (x, y)
                want[(lb, la)] = (x, y)
    recs = collections.defaultdict(list)
    for line in open(a.prefix + '.loci.gated.paf'):
        f = line.split('\t', 12)
        k = (f[0], f[5])
        if k in want:
            q, t = (f[0], f[5]) if f[0] < f[5] else (f[5], f[0])
            qs, qe, ts, te = (int(f[2]), int(f[3]), int(f[7]), int(f[8])) if f[0] < f[5] else (int(f[7]), int(f[8]), int(f[2]), int(f[3]))
            recs[want[k]].append((q, t, qs, qe, ts, te, int(f[9]), int(f[10])))
    rows = collections.Counter()
    with open(a.out, 'w') as w:
        w.write('gene_a\tgene_b\tname_a\tname_b\treason\tn_rec\tbest_bl\tbest_id\tbest_ee_frac\tunion_ee_frac\n')
        for x, y in sorted(missed):
            rr = recs.get((x, y), []) + recs.get((y, x), [])
            nm = (G[x]['name'], G[y]['name'])
            if not gene_loci.get(x) or not gene_loci.get(y):
                reason, stats = 'no_locus', ('', '', '', '')
            elif any((rep.get(la, la), rep.get(lb, lb)) in edges for la in gene_loci[x] for lb in gene_loci[y]):
                reason, stats = 'mcl_split', ('', '', '', '')
            elif not rr:
                reason, stats = 'no_records', ('', '', '', '')
            else:
                best_ee = 0.0
                uq, ut = collections.defaultdict(list), collections.defaultdict(list)
                ok_len = [r for r in rr if r[7] >= 300]
                ok_id = [r for r in ok_len if r[6] / r[7] >= 0.70]
                for q, t, qs, qe, ts, te, m, bl in rr:
                    if m / bl < 0.70:
                        continue
                    qx, tx = exonic(loci[q][3], qs, qe), exonic(loci[t][3], ts, te)
                    small = max(1, min(sum(b - a for a, b in loci[q][3]), sum(b - a for a, b in loci[t][3])))
                    if bl >= 300:
                        best_ee = max(best_ee, min(qx, tx) / small)
                    if qx and tx:
                        uq[(q, t)].append((qs, qe))
                        ut[(q, t)].append((ts, te))
                uee = 0.0
                for k in uq:
                    q, t = k
                    small = max(1, min(sum(b - a for a, b in loci[q][3]), sum(b - a for a, b in loci[t][3])))
                    eq = sum(exonic(loci[q][3], s, e) for s, e in merged(uq[k]))
                    et = sum(exonic(loci[t][3], s, e) for s, e in merged(ut[k]))
                    uee = max(uee, min(eq, et) / small)
                if not ok_len:
                    reason = 'short'
                elif not ok_id:
                    reason = 'low_id'
                elif best_ee < 0.60:
                    reason = 'low_shared'
                else:
                    reason = 'low_cov_or_mcl'
                b = max(rr, key=lambda r: r[7])
                stats = (b[7], f'{b[6] / b[7]:.3f}', f'{best_ee:.3f}', f'{uee:.3f}')
            rows[reason] += 1
            w.write('\t'.join(map(str, (x, y, *nm, reason, len(rr), *stats))) + '\n')
    print(f'missed true pairs ({a.levels}): {len(missed)}')
    for k, v in rows.most_common():
        print(f'  {k} {v}')


if __name__ == '__main__':
    main()

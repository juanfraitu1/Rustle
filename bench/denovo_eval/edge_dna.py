#!/usr/bin/env python3
"""Lift edges against the DNA-level evidence of the gated all-vs-all (vault idea 2026-10-29: Bailey 2002 / SEDEF segmental duplication =
>= 1 kb at >= 90% identity is 'edge evidence of paralogy for YOUNG copies'; exons of old paralogues stay alignable under purifying
selection, introns and flanks do not). Per edge: records in the PAF, max identity, aligned columns split into exon-exon, mixed and
intron-intron (neither base exonic), the longest single-record intron-intron stretch, and the placement identity / coverage of the edge.

  edge_dna.py LOCI_GFF3 GATED_PAF EDGES_TSV [--species gorilla]   (edges keyed chrom:start-end; cols 3-4 = weight, identity)
"""
import collections, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from family_container import aligned_runs, mask_ranges, exon_columns  # noqa: E402
from exon_rule_judge import Judge  # noqa: E402

gff, paf, edges = sys.argv[1:4]
species = sys.argv[sys.argv.index('--species') + 1] if '--species' in sys.argv else 'gorilla'
blocks = collections.defaultdict(list); gene_of_span = {}
for line in open(gff):
    if line[0] == '#':
        continue
    r = line.rstrip('\n').split('\t')
    if r[2] == 'exon':
        blocks[re.search(r'gene=([^;]+)', r[8]).group(1)].append((int(r[3]) - 1, int(r[4])))
    elif r[2] == 'gene':
        gene_of_span[f'{r[0]}:{r[3]}-{r[4]}'] = re.search(r'ID=gene-([^;]+)', r[8]).group(1)
L = {}


def locus(name):
    if name in L:
        return L[name]
    c, s, e = re.match(r'(.+):(\d+)-(\d+)$', name).groups()
    m = []
    for s0, e0 in sorted(blocks[gene_of_span[name]]):
        if m and s0 <= m[-1][1]:
            m[-1][1] = max(m[-1][1], e0)
        else:
            m.append([s0, e0])
    L[name] = (int(s) - 1, [x[0] for x in m], [x[1] for x in m], sum(b - a for a, b in m), int(e) - int(s) + 1)
    return L[name]


want = {}
for l in open(edges):
    a, b, w, ident = l.rstrip('\n').split('\t')[:4]
    want[(min(a, b), max(a, b))] = (float(w), float(ident))
acc = {k: dict(n=0, maxid=0.0, ee=0, mixed=0, ii=0, best_ii=0, best_ii_id=0.0) for k in want}
for line in open(paf):
    f = line.rstrip('\n').split('\t')
    q, t = f[0], f[5]
    k = (min(q, t), max(q, t))
    if q == t or k not in acc:
        continue
    a = acc[k]
    ident = int(f[9]) / max(1, int(f[10]))
    qs, qe, ts, te = int(f[2]), int(f[3]), int(f[7]), int(f[8])
    cg = next((x[5:] for x in f[12:] if x.startswith('cg:Z:')), '')
    qoff, qst, qen, _, _ = locus(q); toff, tst, ten, _, _ = locus(t)
    a['n'] += 1; a['maxid'] = max(a['maxid'], ident)
    ee = sum(h - l for l, h, *_ in exon_columns(cg, f[4], qs, qe, ts, te, qoff, toff, (qst, qen), (tst, ten)))
    tot = 0; qex = 0; tex = 0; rec_ii = 0
    for t0, q0, n in aligned_runs(cg, f[4], qs, qe, ts, te):
        tot += n
        tg, qg = toff + t0, qoff + q0
        tm = mask_ranges(tst, ten, tg, tg + n); qm = mask_ranges(qst, qen, qg, qg + n)
        tex += sum(h - l for l, h, _ in tm); qex += sum(h - l for l, h, _ in qm)
        # intron-intron columns of this run: columns with no exon base on either side
        te_cols = set()
        for l, h, _ in tm:
            te_cols.update(range(l - tg, h - tg))
        qe_cols = set()
        for l, h, _ in qm:
            if f[4] == '+':
                qe_cols.update(range(l - qg, h - qg))
            else:
                qe_cols.update(range(qg + n - h, qg + n - l))
        rec_ii += n - len(te_cols | qe_cols)
    a['ee'] += ee; a['ii'] += rec_ii; a['mixed'] += tot - ee - rec_ii
    if rec_ii > a['best_ii']:
        a['best_ii'], a['best_ii_id'] = rec_ii, ident
J = Judge(species, gff)
print('a\tb\tverdict\tnames\tplace_id\tplace_cov\tn_rec\tmax_id\tee_cols\tmixed_cols\tii_cols\tbest_ii\tbest_ii_id\tsmaller_intronic_bp')
for (a_, b_), (w, pid) in want.items():
    x = acc[(a_, b_)]
    A, B = locus(a_), locus(b_)
    intr = min(A[4] - A[3], B[4] - B[3])
    print(f"{a_}\t{b_}\t{J.verdict(a_, b_)}\t{J.name(a_)}/{J.name(b_)}\t{pid:.3f}\t{min(1.0, w / pid if pid else 0):.3f}\t{x['n']}\t{x['maxid']:.3f}\t"
          f"{x['ee']}\t{x['mixed']}\t{x['ii']}\t{x['best_ii']}\t{x['best_ii_id']:.3f}\t{intr}")

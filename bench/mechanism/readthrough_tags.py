#!/usr/bin/env python3
"""Per-read BAM tags: readthrough (RT) reads vs annotated-junction reads. Needs <tag>.junctions.tsv from
readthrough_rules.py. usage: readthrough_tags.py TAG BAM CHROM OUTDIR"""
import sys, csv, random, bisect, statistics, collections
import pysam
tag, bam, chrom, out = sys.argv[1:5]
random.seed(1)
J = {(int(r['start']), int(r['end'])): r['class'] for r in csv.DictReader(open(f'{out}/{tag}.junctions.tsv'), delimiter='\t')}
feats = ('de', 'NM_per_bp', 'AS_per_bp', 's2_over_s1', 'cm', 'mapq', 'softclip', 'min_anchor', 'local_err')
reads = []                      # (label, junction, feature dict)
for rd in pysam.AlignmentFile(bam).fetch(chrom):
    if rd.flag & 2308:
        continue
    pos, qp = rd.reference_start, 0
    blocks, introns, ops = [], [], rd.cigartuples
    cur = None; errs = []       # (refpos, 1) for X / I / D events
    for op, L in ops:
        if op in (0, 7, 8):
            if cur is None: cur = pos
            if op == 8: errs.extend(range(pos, pos + L))
            pos += L
        elif op == 2:
            errs.append(pos); pos += L
        elif op == 1:
            errs.append(pos)
        elif op == 3:
            blocks.append((cur, pos)); cur = None
            introns.append((pos + 1, pos + L, len(blocks) - 1)); pos += L
    if cur is not None: blocks.append((cur, pos))
    if not introns:
        continue
    cls = [J.get((s, e)) for s, e, _ in introns]
    if 'RT' in cls:
        lab = 'RT'; k = cls.index('RT')
    elif all(c in ('ANN_lt50k', 'ANN_ge50k') for c in cls) and random.random() < 0.10:
        lab = 'CTRL'; k = random.randrange(len(introns))
    else:
        continue
    s, e, bi = introns[k]
    b1, b2 = blocks[bi], blocks[bi + 1]
    errs.sort()
    loc = (bisect.bisect_left(errs, b1[1]) - bisect.bisect_left(errs, b1[1] - 20)) + \
          (bisect.bisect_left(errs, b2[0] + 20) - bisect.bisect_left(errs, b2[0]))
    alen = rd.query_alignment_length or 1
    g = lambda t, d=None: rd.get_tag(t) if rd.has_tag(t) else d
    s1 = g('s1'); s2 = g('s2', 0)
    sc = sum(L for op, L in ops if op == 4)
    reads.append((lab, (s, e), dict(de=g('de'), NM_per_bp=g('NM', 0) / alen, AS_per_bp=g('AS', 0) / alen,
                  s2_over_s1=(s2 / s1) if s1 else None, cm=g('cm'), mapq=rd.mapping_quality, softclip=sc,
                  min_anchor=min(b1[1] - b1[0], b2[1] - b2[0]), local_err=loc)))


def byj0():
    d = collections.defaultdict(list)
    for lab, j, f in reads:
        d[(lab, j)].append(f)
    return d


def auc(p, n):
    ns = sorted(n)
    a = sum(bisect.bisect_left(ns, x) + 0.5 * (bisect.bisect_right(ns, x) - bisect.bisect_left(ns, x)) for x in p) / (len(p) * len(n))
    return a


with open(f'{out}/{tag}.junction_tags.tsv', 'w') as fo:
    fo.write('label\tstart\tend\tn\t' + '\t'.join(feats) + '\n')
    for (lab, (s0, e0)), fs in byj0().items():
        med = [statistics.median([f[ft] for f in fs if f[ft] is not None] or [float('nan')]) for ft in feats]
        fo.write(f'{lab}\t{s0}\t{e0}\t{len(fs)}\t' + '\t'.join(f'{m:.6g}' for m in med) + '\n')
print(f'[{tag}] RT reads {sum(r[0]=="RT" for r in reads)}, control reads {sum(r[0]=="CTRL" for r in reads)}')
byj = collections.defaultdict(list)
for lab, j, f in reads:
    byj[(lab, j)].append(f)
print(f'{"feature":12s} {"RT median":>10s} {"ctrl median":>11s} {"read AUC":>9s} {"junction AUC":>13s}')
for ft in feats:
    p = [f[ft] for lab, _, f in reads if lab == 'RT' and f[ft] is not None]
    n = [f[ft] for lab, _, f in reads if lab == 'CTRL' and f[ft] is not None]
    jp = [statistics.median(x) for x in ([f[ft] for f in fs if f[ft] is not None] for (lab, _), fs in byj.items() if lab == 'RT') if x]
    jn = [statistics.median(x) for x in ([f[ft] for f in fs if f[ft] is not None] for (lab, _), fs in byj.items() if lab == 'CTRL') if x]
    if not p or not n:
        print(f'{ft:12s} (tag absent)'); continue
    a, aj = auc(p, n), auc(jp, jn)
    print(f'{ft:12s} {statistics.median(p):10.4g} {statistics.median(n):11.4g} {max(a, 1 - a):9.3f} {max(aj, 1 - aj):13.3f}')

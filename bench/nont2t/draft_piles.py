#!/usr/bin/env python3
"""PREREG r1539: blind divergence piles on a draft genome (docs/PREREG_absent_on_draft_2026-11-05.md).
  draft_piles.py NT.bam read_origin.tsv OUT [TAU=0.005]   -> OUT.piles.tsv OUT.copies.tsv, summary on stdout
Blind: every primary record clustered by landing (gaps <= 10 kb); pile = >= 10 reads at de > TAU and >= 25% of the cluster.
Scoring with origin labels: copy FOUND when >= 50% of its reads are high-de in a pile or unmapped; pile TRUE when >= 50% of its
labelled high-de reads come from missing copies."""
import collections
import sys

import pysam

nt, origin_tsv, out = sys.argv[1:4]
TAU = float(sys.argv[4]) if len(sys.argv) > 4 else 0.005
org = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1:] for l in open(origin_tsv) if not l.startswith('read\t')}
recs, unmapped = [], set()
for a in pysam.AlignmentFile(nt):
    if a.is_secondary or a.is_supplementary:
        continue
    if a.is_unmapped:
        unmapped.add(a.query_name)
        continue
    recs.append((a.reference_name, a.reference_start, a.reference_end, a.query_name, a.get_tag('de')))
recs.sort()
clusters, cur = [], None
for c, s, e, q, de in recs:
    if cur and cur[0] == c and s <= cur[2] + 10000:
        cur[2] = max(cur[2], e)
        cur[3].append((q, de))
    else:
        cur = [c, s, e, [(q, de)]]
        clusters.append(cur)
piles, in_pile = [], {}
for c, s, e, rs in clusters:
    hi = [q for q, de in rs if de > TAU]
    if len(hi) >= 10 and len(hi) >= 0.25 * len(rs):
        lab = [org[q][3] for q in hi if q in org]
        miss = sum(x == 'missing' for x in lab)
        piles.append((c, s, e, len(rs), len(hi), len(lab), miss, 'TRUE' if lab and miss >= 0.5 * len(lab) else 'FALSE',
                      collections.Counter(org[q][1] for q in hi if q in org).most_common(3)))
        for q in hi:
            in_pile[q] = f'{c}:{s}-{e}'
with open(out + '.piles.tsv', 'w') as w:
    w.write('contig\tstart\tend\treads\thigh_de\tlabelled\tfrom_missing\tverdict\ttop_origins\n')
    for p in piles:
        w.write('\t'.join(map(str, p)) + '\n')
by_copy = collections.defaultdict(list)
for q, o in org.items():
    by_copy[o[0]].append(q)
found = collections.Counter()
with open(out + '.copies.tsv', 'w') as w:
    w.write('gene\tname\tstatus\tn\tin_pile\tunmapped\tfound\n')
    for g, qs in sorted(by_copy.items(), key=lambda x: (x[1] and org[x[1][0]][3], -len(x[1]))):
        name, st = org[qs[0]][1], org[qs[0]][3]
        ip, um = sum(q in in_pile for q in qs), sum(q in unmapped for q in qs)
        f = (ip + um) >= 0.5 * len(qs)
        w.write(f'{g}\t{name}\t{st}\t{len(qs)}\t{ip}\t{um}\t{int(f)}\n')
        if len(qs) >= 10:
            found[(st, f)] += 1
sm, sp = found[('missing', True)], found[('missing', True)] + found[('missing', False)]
pt = sum(p[7] == 'TRUE' for p in piles)
print(f'TAU {TAU}: piles {len(piles)} true {pt} false {len(piles) - pt} (precision {pt / len(piles) if piles else float("nan"):.2f}); '
      f'unmapped reads {len(unmapped)}; scorable missing found {sm}/{sp} (recall {sm / sp if sp else float("nan"):.2f}); '
      f'present found {found[("present", True)]}/{found[("present", True)] + found[("present", False)]}')

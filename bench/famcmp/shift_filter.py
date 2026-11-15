#!/usr/bin/env python3
"""PREREG row 1732: drop weaker shifted-site chains (frame cell of shift_rule.py and / or indel cell of shift_reads.py).

  shift_filter.py CLASS_TSV READS_TSV IN_GTF OUT_GTF union|indel   (READS_TSV rows = the S1 rows of CLASS_TSV, in order)"""
import sys, csv, re, collections
from shift_rule import score
csv.field_size_limit(1 << 30); cls_p, rd_p, gin, gout, mode = sys.argv[1:6]
C = list(csv.DictReader(open(cls_p), delimiter='\t')); RD = iter(csv.DictReader(open(rd_p), delimiter='\t')); drop = []
for i, r in enumerate(C):
    rd = next(RD) if r['class'] == 'S1' else None; sh, n, pn = int(r['shift']), int(r['reads']), int(r['partner_reads'])
    frame = mode == 'union' and r['class'] in ('S1', 'S2') and sh <= 10 and sh % 3 and n < 0.05 * pn
    indel = rd is not None and int(rd['nK']) > 0 and int(rd['indelK']) / int(rd['nK']) >= 0.5
    if frame or indel: drop.append(i)
assert next(RD, None) is None, 'READS_TSV does not match CLASS_TSV'
rows = [(r['class'], int(r['shift']), int(r['reads']), int(r['partner_reads']), int(r['in_ref']), int(r['intron_annotated']), int(r['partner_in_ref'])) for r in C]
s = score(rows, drop); gone = {(C[i]['contig'], C[i]['strand'], C[i]['chain']) for i in drop}
ex = collections.defaultdict(list); T = re.compile(r'transcript_id "([^"]+)"'); L = open(gin).readlines()
for l in L:
    f = l.split('\t')
    if len(f) > 8 and f[2] == 'exon': ex[f[0], f[6], T.search(f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
bad = {k[2] for k, e in ex.items() if (e.sort() or True) and (k[0], k[1], str(tuple((a[1], b[0]) for a, b in zip(e, e[1:])))) in gone}
with open(gout, 'w') as w: w.writelines(l for l in L if not ((m := T.search(l)) and m.group(1) in bad))
print(f"{mode}: dropped {s['dropped']} chains ({len(bad)} transcripts), in RefSeq {s['in_ref']:.3f}; precision {s['prec']:.4f} -> {s['prec_rule']:.4f}, read cut {s['prec_cut']:.4f}; correct lost {s['correct_lost']:.4f}; annotated pairs lost (GTF count) {s['pairs_lost']:.4f}")

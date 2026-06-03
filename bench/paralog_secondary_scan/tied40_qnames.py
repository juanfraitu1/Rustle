#!/usr/bin/env python3
"""Reproduce the EXACT 40-read tied set using the same definition as
tied_reads_daz.py: a read with a placement in BOTH DAZ1 and DAZ3 whose
best-AS placements satisfy |AS_DAZ1 - AS_DAZ3| <= 5.  Print qnames.
"""
import re, collections

DAZ1 = (42783133, 42859657)   # - strand
DAZ3 = (42879918, 42945552)   # + strand (inverted paralog)

def copy_of(pos, endp):
    mid = (pos + endp) // 2
    if DAZ1[0] <= mid <= DAZ1[1]: return 'DAZ1'
    if DAZ3[0] <= mid <= DAZ3[1]: return 'DAZ3'
    return None

def ref_len(cig):
    return sum(int(n) for n, op in re.findall(r'(\d+)([MIDNSH=X])', cig) if op in 'MDN=X')

def tags(fields):
    d = {}
    for t in fields[11:]:
        p = t.split(':')
        d[p[0]] = p[-1]
    return d

reads = collections.defaultdict(dict)   # qname -> copy -> best record (max AS)
for line in open('/tmp/daz_aln.sam'):
    f = line.rstrip('\n').split('\t')
    qn, flag, pos, cig = f[0], int(f[1]), int(f[3]), f[5]
    if cig == '*': continue
    cp = copy_of(pos, pos + ref_len(cig))
    if cp is None: continue
    tg = tags(f)
    rec = dict(AS=int(tg.get('AS', -1)), NM=int(tg.get('NM', -1)),
               de=float(tg.get('de', 'nan')), pos=pos, cig=cig,
               strand='-' if flag & 0x10 else '+')
    cur = reads[qn].get(cp)
    if cur is None or rec['AS'] > cur['AS']:
        reads[qn][cp] = rec

both = {q: c for q, c in reads.items() if 'DAZ1' in c and 'DAZ3' in c}
tied = sorted(q for q, c in both.items() if abs(c['DAZ1']['AS'] - c['DAZ3']['AS']) <= 5)
print(f"both={len(both)}  tied(|dAS|<=5)={len(tied)}")
for q in tied:
    print(q)

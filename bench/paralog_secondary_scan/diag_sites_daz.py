#!/usr/bin/env python3
"""Do the AS-tied DAZ reads even COVER copy-distinguishing positions?
1. parse the DAZ1-vs-DAZ3rc cs alignment -> diagnostic (substitution) positions
2. for each read placed in both copies, count diagnostic positions inside its
   DAZ1 aligned blocks, split by whether AS is tied vs has a score gap.
Inputs: /tmp/daz1_vs_daz3.sam (cs=long), /tmp/daz_aln.sam
"""
import re, collections

DAZ1_OFF = 42783133            # 1-based genomic start of the daz1 contig
DAZ1 = (42783133, 42859657); DAZ3 = (42879918, 42945552)

# ---- 1. diagnostic positions from the copy-copy cs tag ----
rec = open('/tmp/daz1_vs_daz3.sam').readline().rstrip('\n').split('\t')
cs = next(t[5:] for t in rec[11:] if t.startswith('cs:Z:'))
refpos = int(rec[3])                      # 1-based in contig
diag = []                                 # genomic positions of substitutions
indel = 0
for op, val in re.findall(r'([=:*+\-~])([A-Za-z0-9]+)', cs):
    if op == '=':                          # run of matches (long form: bases)
        refpos += len(val)
    elif op == ':':                        # run of matches (short form: count)
        refpos += int(val)
    elif op == '*':                        # substitution ref->alt : 2 letters
        diag.append(DAZ1_OFF + refpos - 1)
        refpos += 1
    elif op == '-':                        # deletion from query (ref bases consumed)
        indel += 1; refpos += len(val)
    elif op == '+':                        # insertion (no ref consumed)
        indel += 1
    elif op == '~':                        # intron/splice ~gt###ag : digits = len
        refpos += int(re.search(r'\d+', val).group())
diag.sort()
print(f"diagnostic substitution positions between DAZ1 and DAZ3: {len(diag)}  "
      f"(+{indel} indel events)  over ~{DAZ1[1]-DAZ1[0]} bp")
print(f"  => ~1 diagnostic SNP every {(DAZ1[1]-DAZ1[0])//max(1,len(diag))} bp")

import bisect
def count_diag(blocks):
    n = 0
    for s, e in blocks:
        lo = bisect.bisect_left(diag, s); hi = bisect.bisect_right(diag, e)
        n += hi - lo
    return n

# ---- 2. per-read: AS gap + diagnostic positions covered in DAZ1 blocks ----
def ref_len(cig):
    return sum(int(n) for n,o in re.findall(r'(\d+)([MIDNSH=X])',cig) if o in 'MDN=X')
def blocks(pos, cig):                       # genomic aligned (M/=/X) intervals
    out=[]; cur=pos
    for n,o in re.findall(r'(\d+)([MIDNSH=X])',cig):
        n=int(n)
        if o in 'M=X': out.append((cur,cur+n)); cur+=n
        elif o in 'DN': cur+=n
    return out
def copy_of(pos,endp):
    mid=(pos+endp)//2
    if DAZ1[0]<=mid<=DAZ1[1]: return 'DAZ1'
    if DAZ3[0]<=mid<=DAZ3[1]: return 'DAZ3'
    return None

reads=collections.defaultdict(dict)
for line in open('/tmp/daz_aln.sam'):
    f=line.rstrip('\n').split('\t'); cig=f[5]
    if cig=='*': continue
    pos=int(f[3]); cp=copy_of(pos,pos+ref_len(cig))
    if cp is None: continue
    AS=next((int(t[5:]) for t in f[11:] if t.startswith('AS:i:')),-1)
    cur=reads[f[0]].get(cp)
    if cur is None or AS>cur['AS']:
        reads[f[0]][cp]=dict(AS=AS,pos=pos,cig=cig)

both={q:c for q,c in reads.items() if 'DAZ1' in c and 'DAZ3' in c}
buck=collections.defaultdict(list)
for q,c in both.items():
    dAS=c['DAZ1']['AS']-c['DAZ3']['AS']
    nd=count_diag(blocks(c['DAZ1']['pos'],c['DAZ1']['cig']))
    key='TIED |dAS|<=5' if abs(dAS)<=5 else 'gap 6-20' if abs(dAS)<=20 else \
        'gap 21-100' if abs(dAS)<=100 else 'gap >100'
    buck[key].append(nd)

print("\ndiagnostic positions COVERED by the read, grouped by AS gap:")
print(f"   {'group':<16} {'n_reads':>7} {'mean':>6} {'median':>7} {'cover 0':>8} {'cover>=3':>9}")
for k in ['TIED |dAS|<=5','gap 6-20','gap 21-100','gap >100']:
    v=buck.get(k,[])
    if not v: continue
    v.sort(); med=v[len(v)//2]; mean=sum(v)/len(v)
    z=sum(1 for x in v if x==0); ge3=sum(1 for x in v if x>=3)
    print(f"   {k:<16} {len(v):>7} {mean:>6.1f} {med:>7} {z:>8} {ge3:>9}")

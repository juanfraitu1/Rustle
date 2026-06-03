#!/usr/bin/env python3
"""Deep-dive the two anomalies:
 A) read 50660416: 291 bp inside the copy-specific block yet 'tied'. Show its
    full DAZ1 exon chain near the block, its DAZ3 exon chain, and whether AS/NM
    really tie. If it genuinely covers copy-specific bp, that's a SIGNAL.
 B) the 9 colinear-FALSE & not-direct-matched junction incidences: list them,
    locate relative to the 10.8 kb block start, and decide artifact vs signal.
"""
import re,json,collections,bisect
DAZ1=(42783133,42859657); DAZ3=(42879918,42945552)
cm=json.load(open('/tmp/copymap.json')); DAZ1_OFF=cm['daz1_off']
struct=[tuple(x) for x in cm['struct']]; shared_segs=cm['shared_segs']
DAZ3_LEN=DAZ3[1]-DAZ3[0]; S0,S1=struct[0]
seg_starts=[s for s,_,_ in shared_segs]
def delta(local):
    i=bisect.bisect_right(seg_starts,local)-1
    if i<0: return None
    s,e,d=shared_segs[i]; return d if s<=local<e else None
def copy_of(p,e):
    m=(p+e)//2
    return 'DAZ1' if DAZ1[0]<=m<=DAZ1[1] else 'DAZ3' if DAZ3[0]<=m<=DAZ3[1] else None
def ref_len(c): return sum(int(n) for n,o in re.findall(r'(\d+)([MIDNSH=X])',c) if o in 'MDN=X')
def exon_blocks(p,c):
    out=[];cur=p;ss=p
    for n,o in re.findall(r'(\d+)([MIDNSH=X])',c):
        n=int(n)
        if o in 'M=XD': cur+=n
        elif o=='N': out.append((ss,cur));cur+=n;ss=cur
    out.append((ss,cur)); return out
def introns(p,c):
    out=[];cur=p
    for n,o in re.findall(r'(\d+)([MIDNSH=X])',c):
        n=int(n)
        if o in 'M=XD': cur+=n
        elif o=='N': out.append((cur,cur+n));cur+=n
    return out
def d1_to_d3(g):
    local=g-DAZ1_OFF+1; d=delta(local)
    if d is None: return None
    q=local+d; span=DAZ3_LEN-q+1; return DAZ3[0]+span-1

# load ALL records (incl tags) for the two anomalous reads
TIED=set(l.strip() for l in open('/tmp/tied40.txt'))
recs=collections.defaultdict(list)
for line in open('/tmp/daz_aln.sam'):
    f=line.rstrip('\n').split('\t');cig=f[5]
    if cig=='*'or f[0] not in TIED: continue
    p=int(f[3]);cp=copy_of(p,p+ref_len(cig))
    if cp is None: continue
    tg={t.split(':')[0]:t.split(':')[-1] for t in f[11:]}
    recs[f[0]].append(dict(cp=cp,pos=p,cig=cig,AS=int(tg.get('AS',-1)),
        NM=int(tg.get('NM',-1)),de=tg.get('de','?'),sec=bool(int(f[1])&0x100),
        strand='-' if int(f[1])&0x10 else '+'))

def best(qn,cp):
    cand=[r for r in recs[qn] if r['cp']==cp]
    return max(cand,key=lambda r:r['AS'])

print("===== A) read 50660416: does it cover copy-specific sequence? =====")
qn=next(q for q in TIED if q.endswith('50660416/ccs'))
d1=best(qn,'DAZ1'); d3=best(qn,'DAZ3')
print(f"  DAZ1: pos={d1['pos']} strand={d1['strand']} AS={d1['AS']} NM={d1['NM']} de={d1['de']}")
print(f"  DAZ3: pos={d3['pos']} strand={d3['strand']} AS={d3['AS']} NM={d3['NM']} de={d3['de']}")
print(f"  AS tie: {d1['AS']}=={d3['AS']} -> {d1['AS']==d3['AS']}; NM tie: {d1['NM']==d3['NM']}")
bl=exon_blocks(d1['pos'],d1['cig'])
print(f"  DAZ1 exon chain ({len(bl)} exons); block start S0={S0}:")
for i,(s,e) in enumerate(bl):
    lo=max(s,S0);hi=min(e,S1); inb=max(0,hi-lo)
    tag='  <-- in COPY BLOCK' if inb>0 else ''
    print(f"    e{i+1:>2}: {s}-{e} len={e-s} in_block={inb}{tag}")
bl3=exon_blocks(d3['pos'],d3['cig'])
print(f"  DAZ3 exon chain ({len(bl3)} exons) (plus strand):")
for i,(s,e) in enumerate(bl3):
    print(f"    e{i+1:>2}: {s}-{e} len={e-s}")
# Are the in-block exons mapped to real DAZ3 exons (so they're NOT copy-unique)?
print("  in-block DAZ1 exons mapped to DAZ3 via copy map:")
for s,e in bl:
    if max(s,S0)<min(e,S1):
        ms=d1_to_d3(s); me=d1_to_d3(e-1)
        print(f"    DAZ1 {s}-{e}  -> DAZ3 {ms},{me} (None means no colinear backbone = TRULY copy-specific)")

print("\n===== B) the 9 colinear-FALSE & no-direct-match junctions =====")
print(f"{'read':<20}{'DAZ1 junc':>26}{'dist_donor-to-S0':>18}{'dist_accept-to-S0':>18}")
n=0
for qn in TIED:
    d1=best(qn,'DAZ1'); d3=best(qn,'DAZ3')
    j3s=[tuple(sorted(x)) for x in introns(d3['pos'],d3['cig'])]
    for gd,ga in introns(d1['pos'],d1['cig']):
        ld=delta(gd-DAZ1_OFF); la=delta(ga-DAZ1_OFF+1)
        colinear=(ld is not None and la is not None and ld==la)
        if colinear: continue
        a=d1_to_d3(gd);b=d1_to_d3(ga); direct=False
        if a is not None and b is not None:
            lo,hi=sorted((a,b)); direct=any(abs(x-lo)<=12 and abs(y-hi)<=12 for x,y in j3s)
        if direct: continue
        n+=1
        print(f"{qn[-18:]:<20}{f'{gd}-{ga}':>26}{gd-S0:>18}{ga-S0:>18}")
print(f"  total = {n}")
print(f"  (block spans S0={S0}..S1={S1}; a junction whose donor/acceptor straddle")
print(f"   the block START is one whose intron is the copy-specific deletion edge)")

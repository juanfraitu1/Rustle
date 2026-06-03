#!/usr/bin/env python3
"""ADVERSARIAL probes:
 (P1) which tied reads put ANY bp in the 10.8 kb copy-specific block, how much,
      and is it a real internal exon or just an alignment-edge overhang at the
      block boundary?  If a read genuinely covers copy-specific sequence it
      would NOT be tied (it would carry copy signal) -> resolve the paradox.
 (P2) the 'shared_in_both=False' junctions in channel B: are they genuinely
      copy-specific, or just artifacts of the colinear-delta test failing at
      SNV/indel-induced segment boundaries?  Cross-check with the DIRECT
      same-read DAZ3-junction match.
 (P3) Do the tied reads' DAZ1 and DAZ3 placements have the SAME number of
      junctions and same intron-size signature (structure recoverable either way)?
"""
import re, json, collections, bisect

DAZ1=(42783133,42859657); DAZ3=(42879918,42945552)
cm=json.load(open('/tmp/copymap.json')); DAZ1_OFF=cm['daz1_off']
struct=[tuple(x) for x in cm['struct']]; shared_segs=cm['shared_segs']
DAZ3_LEN=DAZ3[1]-DAZ3[0]
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

TIED=set(l.strip() for l in open('/tmp/tied40.txt'))
reads=collections.defaultdict(dict)
for line in open('/tmp/daz_aln.sam'):
    f=line.rstrip('\n').split('\t');cig=f[5]
    if cig=='*'or f[0] not in TIED: continue
    p=int(f[3]);cp=copy_of(p,p+ref_len(cig))
    if cp is None: continue
    AS=next((int(t[5:]) for t in f[11:] if t.startswith('AS:i:')),-1)
    cu=reads[f[0]].get(cp)
    if cu is None or AS>cu['AS']: reads[f[0]][cp]=dict(AS=AS,pos=p,cig=cig)

S0,S1=struct[0]
def bp_in_blk(bl):
    t=0
    for s,e in bl:
        lo=max(s,S0);hi=min(e,S1)
        if lo<hi:t+=hi-lo
    return t

print("=== P1: tied reads overlapping the 10.8 kb copy-specific block "
      f"({S0}-{S1}) in their DAZ1 placement ===")
print(f"{'read':<24}{'bp_in_blk':>9}{'overlap_exon(s)':>40}{'dist_to_block_edge':>20}")
tot=0
for qn,c in reads.items():
    bl=exon_blocks(c['DAZ1']['pos'],c['DAZ1']['cig'])
    b=bp_in_blk(bl);
    if b==0: continue
    tot+=b
    ovx=[(s,e) for s,e in bl if max(s,S0)<min(e,S1)]
    # distance from the overlapping exon's far boundary into the block
    desc=[]
    for s,e in ovx:
        lo=max(s,S0);hi=min(e,S1)
        # how far does it protrude past block start S0?
        desc.append(f"{s}-{e}(into+{hi-S0})")
    print(f"{qn[-20:]:<24}{b:>9}{(' '.join(desc)):>40}{hi-S0:>20}")
print(f"  TOTAL bp_in_block over tied reads = {tot}")
print("  Interpretation check: all overlaps are <=? bp protrusions at the block")
print("  START edge (42828911). A read truly INSIDE the 10.8kb block would be")
print("  DAZ1-unique and could NOT be tied.\n")

print("=== P2: per-read DAZ1 vs DAZ3 junction COUNT + size signature ===")
print(f"{'read':<22}{'nJ_DAZ1':>8}{'nJ_DAZ3':>8}{'same_count':>11}{'direct_match/J1':>16}")
def d1_to_d3(g):
    local=g-DAZ1_OFF+1; d=delta(local)
    if d is None: return None
    q=local+d; span=DAZ3_LEN-q+1; return DAZ3[0]+span-1
TOL=12
agg_inc=agg_match=0; same_count=0
for qn,c in sorted(reads.items()):
    j1=introns(c['DAZ1']['pos'],c['DAZ1']['cig'])
    j3=introns(c['DAZ3']['pos'],c['DAZ3']['cig'])
    j3s=[tuple(sorted(x)) for x in j3]
    m=0
    for gd,ga in j1:
        a=d1_to_d3(gd);b=d1_to_d3(ga)
        if a is None or b is None: continue
        lo,hi=sorted((a,b))
        if any(abs(x-lo)<=TOL and abs(y-hi)<=TOL for x,y in j3s): m+=1
    agg_inc+=len(j1);agg_match+=m
    sc=len(j1)==len(j3); same_count+=sc
    print(f"{qn[-18:]:<22}{len(j1):>8}{len(j3):>8}{str(sc):>11}{f'{m}/{len(j1)}':>16}")
print(f"\n  reads with identical junction COUNT in both copies: {same_count}/{len(reads)}")
print(f"  DIRECT junction match (mapped DAZ1 junc found in same read's DAZ3): "
      f"{agg_match}/{agg_inc} ({100*agg_match/max(1,agg_inc):.1f}%)")

print("\n=== P3: of the 9 colinear-FALSE junction incidences, how many are still")
print("    DIRECT-matched (i.e. the test, not the biology, was at fault)? ===")
cb_false_but_direct=0; cb_false_and_nodirect=0
for qn,c in reads.items():
    j1=introns(c['DAZ1']['pos'],c['DAZ1']['cig'])
    j3s=[tuple(sorted(x)) for x in introns(c['DAZ3']['pos'],c['DAZ3']['cig'])]
    for gd,ga in j1:
        ld=delta(gd-DAZ1_OFF); la=delta(ga-DAZ1_OFF+1)
        colinear = (ld is not None and la is not None and ld==la)
        if colinear: continue
        a=d1_to_d3(gd);b=d1_to_d3(ga)
        direct=False
        if a is not None and b is not None:
            lo,hi=sorted((a,b))
            direct=any(abs(x-lo)<=TOL and abs(y-hi)<=TOL for x,y in j3s)
        if direct: cb_false_but_direct+=1
        else: cb_false_and_nodirect+=1
print(f"  colinear-FALSE but DIRECT-matched (test artifact): {cb_false_but_direct}")
print(f"  colinear-FALSE and NOT direct-matched (genuinely copy-edge): "
      f"{cb_false_and_nodirect}")

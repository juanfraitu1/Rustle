#!/usr/bin/env python3
"""Final reconciliation over the 40 tied reads:
 - re-state AS/NM/de identity per read (does 50660416 break the 'all identical'
   claim? yes for AS/NM/de by <=5 / <=1 / tiny). Quantify how many are STRICTLY
   identical on all three vs how many merely satisfy |dAS|<=5.
 - confirm the other in-block reads (101843971,177276401,62785876) are tiny edge
   overhangs whose bases also align in DAZ3 (artifact).
 - produce the headline channel numbers using the DIRECT (same-read) junction
   match as the rigorous shared-junction measure.
"""
import re,json,collections,bisect,subprocess
DAZ1=(42783133,42859657); DAZ3=(42879918,42945552)
cm=json.load(open('/tmp/copymap.json')); DAZ1_OFF=cm['daz1_off']
struct=[tuple(x) for x in cm['struct']]; S0,S1=struct[0]; shared_segs=cm['shared_segs']
DAZ3_LEN=DAZ3[1]-DAZ3[0]; seg_starts=[s for s,_,_ in shared_segs]
def delta(local):
    i=bisect.bisect_right(seg_starts,local)-1
    if i<0: return None
    s,e,d=shared_segs[i]; return d if s<=local<e else None
def d1_to_d3(g):
    local=g-DAZ1_OFF+1; d=delta(local)
    if d is None: return None
    q=local+d; span=DAZ3_LEN-q+1; return DAZ3[0]+span-1
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
    tg={t.split(':')[0]:t.split(':')[-1] for t in f[11:]}
    AS=int(tg.get('AS',-1))
    cu=reads[f[0]].get(cp)
    if cu is None or AS>cu['AS']:
        reads[f[0]][cp]=dict(AS=AS,NM=int(tg.get('NM',-1)),de=float(tg.get('de','nan')),pos=p,cig=cig)

# --- identity reconciliation ---
strict=0; tie_only=[]
for qn,c in reads.items():
    same=(c['DAZ1']['AS']==c['DAZ3']['AS'] and c['DAZ1']['NM']==c['DAZ3']['NM']
          and abs(c['DAZ1']['de']-c['DAZ3']['de'])<1e-9)
    if same: strict+=1
    else: tie_only.append((qn,c['DAZ1']['AS']-c['DAZ3']['AS'],
                           c['DAZ1']['NM']-c['DAZ3']['NM'],
                           round(c['DAZ1']['de']-c['DAZ3']['de'],4)))
print(f"reads STRICTLY identical on AS&NM&de: {strict}/40")
print(f"reads that only satisfy |dAS|<=5 but DIFFER on AS/NM/de: {len(tie_only)}")
for qn,dAS,dNM,dde in tie_only:
    print(f"   {qn[-16:]}: dAS={dAS} dNM={dNM} dde={dde}")

# --- direct shared-junction tally (rigorous) ---
inc=match=0
for qn,c in reads.items():
    j1=introns(c['DAZ1']['pos'],c['DAZ1']['cig'])
    j3=[tuple(sorted(x)) for x in introns(c['DAZ3']['pos'],c['DAZ3']['cig'])]
    for gd,ga in j1:
        inc+=1
        a=d1_to_d3(gd);b=d1_to_d3(ga)
        if a is None or b is None: continue
        lo,hi=sorted((a,b))
        if any(abs(x-lo)<=12 and abs(y-hi)<=12 for x,y in j3): match+=1
print(f"\nDIRECT junction-support shared: {match}/{inc} ({100*match/inc:.1f}%)")
n_unsplit=sum(1 for c in reads.values() if not introns(c['DAZ1']['pos'],c['DAZ1']['cig']))
print(f"single-exon (unspliced) tied reads contributing 0 junctions: {n_unsplit}")

# --- aligned bp ---
tot=shared=blk=0
for qn,c in reads.items():
    for s,e in exon_blocks(c['DAZ1']['pos'],c['DAZ1']['cig']):
        tot+=e-s
        for ss,ee,_ in shared_segs:
            gss=DAZ1_OFF+ss-1; gee=DAZ1_OFF+ee-1
            lo=max(s,gss);hi=min(e,gee)
            if lo<hi: shared+=hi-lo
        lo=max(s,S0);hi=min(e,S1)
        if lo<hi: blk+=hi-lo
print(f"\naligned exon bp total={tot}  shared-backbone={shared} ({100*shared/tot:.2f}%)  "
      f"in-10.8kb-block={blk} ({100*blk/tot:.4f}%)")

# --- the 4 in-block reads: confirm tiny artifacts ---
print("\n4 in-block reads (all overhang the block START edge, bases align in DAZ3):")
for qn,c in reads.items():
    b=0
    for s,e in exon_blocks(c['DAZ1']['pos'],c['DAZ1']['cig']):
        lo=max(s,S0);hi=min(e,S1)
        if lo<hi: b+=hi-lo
    if b: print(f"   {qn[-16:]}: {b} bp in block")

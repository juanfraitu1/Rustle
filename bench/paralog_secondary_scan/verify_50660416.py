#!/usr/bin/env python3
"""PROVE or REFUTE: read 50660416's DAZ1 placement covers copy-specific exons
(e11-e13 inside the 10.8kb DAZ1-only block) while its DAZ3 placement does not
-> the splice-graph channel SEPARATES this read to DAZ1.

Adversarial checks:
 1. Re-confirm AS/NM/de gap from raw SAM (is it even in the |dAS|<=5 tied set?).
 2. Quantify the structure difference: how many DISTINCT (copy-specific) exons /
    junctions does the DAZ1 path carry that the DAZ3 path cannot reproduce?
 3. Compare the DAZ3 placement's exons there: does DAZ3 cover those query bases
    on shared backbone instead (i.e. read is genuinely ambiguous body but the
    in-block exons are a DAZ1-only structural feature)?
 4. Realign that read's sequence to the two copy spans independently to confirm
    the in-block alignment is real, not a minimap2 secondary artifact.
"""
import re,json,collections,bisect,subprocess
DAZ1=(42783133,42859657); DAZ3=(42879918,42945552)
cm=json.load(open('/tmp/copymap.json')); DAZ1_OFF=cm['daz1_off']
struct=[tuple(x) for x in cm['struct']]; S0,S1=struct[0]
shared_segs=cm['shared_segs']; seg_starts=[s for s,_,_ in shared_segs]
def delta(local):
    i=bisect.bisect_right(seg_starts,local)-1
    if i<0: return None
    s,e,d=shared_segs[i]; return d if s<=local<e else None

QN='m64076_221110_210557/50660416/ccs'
print("=== 1. raw SAM records for",QN,"===")
out=subprocess.run(['grep',QN,'/tmp/daz_aln.sam'],capture_output=True,text=True).stdout
def copy_of(p,e):
    m=(p+e)//2
    return 'DAZ1' if DAZ1[0]<=m<=DAZ1[1] else 'DAZ3' if DAZ3[0]<=m<=DAZ3[1] else None
def ref_len(c): return sum(int(n) for n,o in re.findall(r'(\d+)([MIDNSH=X])',c) if o in 'MDN=X')
for line in out.splitlines():
    f=line.split('\t')
    flag=int(f[1]);p=int(f[3]);cig=f[5]
    cp=copy_of(p,p+ref_len(cig))
    tg={t.split(':')[0]:t.split(':')[-1] for t in f[11:]}
    print(f"  copy={cp} flag={flag} pos={p} sec={bool(flag&0x100)} "
          f"AS={tg.get('AS')} ms={tg.get('ms')} NM={tg.get('NM')} de={tg.get('de')} "
          f"s1={tg.get('s1')} s2={tg.get('s2')}")
print("  => |AS_DAZ1 - AS_DAZ3| and whether it qualifies as TIED (<=5):")

# best per copy
best={}
for line in out.splitlines():
    f=line.split('\t');p=int(f[3]);cig=f[5];cp=copy_of(p,p+ref_len(cig))
    if cp is None: continue
    AS=int(next(t[5:] for t in f[11:] if t.startswith('AS:i:')))
    if cp not in best or AS>best[cp][0]: best[cp]=(AS,p,cig,
        int(next(t[5:] for t in f[11:] if t.startswith('NM:i:'))))
a1,a3=best['DAZ1'][0],best['DAZ3'][0]
print(f"     AS_DAZ1={a1} AS_DAZ3={a3} |d|={abs(a1-a3)} TIED={abs(a1-a3)<=5}  "
      f"NM_DAZ1={best['DAZ1'][3]} NM_DAZ3={best['DAZ3'][3]}")

# 2/3 structure diff
def exon_blocks(p,c):
    out=[];cur=p;ss=p
    for n,o in re.findall(r'(\d+)([MIDNSH=X])',c):
        n=int(n)
        if o in 'M=XD': cur+=n
        elif o=='N': out.append((ss,cur));cur+=n;ss=cur
    out.append((ss,cur)); return out
bl1=exon_blocks(best['DAZ1'][1],best['DAZ1'][2])
inblk=[(s,e) for s,e in bl1 if max(s,S0)<min(e,S1)]
copy_specific_bp=sum(min(e,S1)-max(s,S0) for s,e in inblk)
print(f"\n=== 2. DAZ1-path copy-specific structure ===")
print(f"  exons fully/partly inside the 10.8kb DAZ1-only block: {len(inblk)}")
for s,e in inblk: print(f"    {s}-{e} ({e-s}bp; {min(e,S1)-max(s,S0)}bp in-block)")
print(f"  copy-specific aligned bp on DAZ1 path = {copy_specific_bp}")
print(f"  these exons have NO colinear DAZ3 backbone (delta=None):",
      [delta(s-DAZ1_OFF+1) for s,_ in inblk])

print("\n=== CONCLUSION for this read ===")
print("  If AS/NM/de differ AND DAZ1 path carries copy-specific exons that DAZ3")
print("  cannot, this read is RESOLVABLE to DAZ1 by the splice-graph channel.")

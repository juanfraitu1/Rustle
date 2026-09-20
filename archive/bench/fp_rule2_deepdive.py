#!/usr/bin/env python
"""Deep-dive Rule 2 (genomic reciprocal overlap) + 5th candidate (length asymmetry).
Print coords/strand/E_p-family for every overlap hit + collateral, to judge whether
'collateral' pairs are real copies or adjacent/nested-gene truth-axis artifacts."""
import os
from collections import defaultdict
from itertools import combinations
BENCH="/mnt/c/Users/jfris/Desktop/Rustle/bench"; SCRATCH="/home/juanfra/winloci_scratch"
exec(open("/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/03f0156d-63b6-4d0f-bf09-862bddee491c/scratchpad/_shared.py").read())

def geom(k):
    a,b=tuple(k); ma,mb=gmeta.get(a),gmeta.get(b)
    if not ma or not mb: return None
    if ma[0]!=mb[0]: return (False,None,0.0,ma,mb)
    opp=ma[3]!=mb[3]; ov=min(ma[2],mb[2])-max(ma[1],mb[1])
    la,lb=ma[2]-ma[1],mb[2]-mb[1]; rof=ov/max(1,min(la,lb)) if ov>0 else 0.0
    return (True,opp,rof,ma,mb)

def epfam(g): return g2f.get(g,"-")

print("========== Rule2: ALL same-contig reciprocal-overlap>0 within-family pairs (FP + REAL) ==========")
rows=[]
for k in list(fp_pairs)+list(hom_pairs):
    g=geom(k)
    if not g or not g[0] or g[2]<=0: continue
    a,b=tuple(k); ma,mb=g[3],g[4]
    is_fp = k in fp_set
    rows.append((g[2], g[1], is_fp, a,b, ma,mb, fp_type.get(k,"REAL"), epfam(a),epfam(b)))
rows.sort(reverse=True)
print(f"{'rof':>5} {'opp':>4} {'FP?':>4} {'typeOrREAL':<22} {'Aep':>7}{'Bep':>7}  pair (coords)")
for rof,opp,is_fp,a,b,ma,mb,typ,fa,fb in rows:
    same_ep = "SAME_Ep" if (fa==fb and fa!="-") else ""
    print(f"{rof:5.2f} {str(opp):>4} {str(is_fp):>4} {typ:<22} {fa[:6]:>7}{fb[:6]:>7}  {a}[{ma[3]} {ma[1]}-{ma[2]}] + {b}[{mb[3]} {mb[1]}-{mb[2]}] {same_ep}")

print("\n========== 5th candidate: gene-span length-ratio as an FP rule ==========")
def lenratio(k):
    a,b=tuple(k); ma,mb=gmeta.get(a),gmeta.get(b)
    if not ma or not mb: return None
    la,lb=ma[2]-ma[1],mb[2]-mb[1]; return min(la,lb)/max(la,lb)
for thr in (0.1,0.2,0.3):
    fp_hit=[k for k in genuine_fp if (lenratio(k) is not None and lenratio(k)<thr)]
    real_hit=[k for k in hom_pairs if (lenratio(k) is not None and lenratio(k)<thr)]
    orc=[k for k in fp_hit+real_hit if all(g in oracle_multi for g in k)]
    from collections import Counter
    print(f" lenratio<{thr}: genuineFP={len(fp_hit)} REALcollat={len(real_hit)} oracle={len(orc)} bytype={dict(Counter(fp_type[k] for k in fp_hit))}")

#!/usr/bin/env python
import os
from collections import Counter, defaultdict
exec(open("/tmp/claude-1000/-mnt-c-Users-jfris-Desktop/03f0156d-63b6-4d0f-bf09-862bddee491c/scratchpad/_shared.py").read())

def geom(k):
    a,b=tuple(k); ma,mb=gmeta.get(a),gmeta.get(b)
    if not ma or not mb: return None
    if ma[0]!=mb[0]: return (False,None,0.0)
    opp=ma[3]!=mb[3]; ov=min(ma[2],mb[2])-max(ma[1],mb[1])
    la,lb=ma[2]-ma[1],mb[2]-mb[1]; rof=ov/max(1,min(la,lb)) if ov>0 else 0.0
    return (True,opp,rof)
def report(name,pred,note=""):
    fp=[k for k in fp_pairs if pred(k)]; real=[k for k in hom_pairs if pred(k)]
    gen=[k for k in fp if fp_type[k] not in TRUTH_ARTIFACT]
    orc=[k for k in fp+real if all(g in oracle_multi for g in k)]
    sameep=[k for k in real if g2f.get(tuple(k)[0]) and g2f.get(tuple(k)[0])==g2f.get(tuple(k)[1]) and g2f.get(tuple(k)[0]) not in mega]
    print(f"\n=== {name} === {note}")
    print(f"  FP={len(fp)} genuine={len(gen)} | REAL(E_p/DNA)collat={len(real)} (SAME_Ep-confirmed-paralog={len(sameep)}) | oracle-pair={len(orc)}")
    print(f"  FP-by-type={dict(Counter(fp_type[k] for k in fp))}")
    if real: print("  REAL collat: "+" ; ".join(f"{tuple(k)[0]}+{tuple(k)[1]}"+("[SAME_Ep]" if k in sameep else "") for k in real))
    return fp,real,gen

# Refined Rule 2: opposite strand AND high reciprocal span overlap
for X in (0.5,0.6,0.7):
    report(f"R2ref opp-strand & recip-overlap>={X}", lambda k,X=X:(geom(k) and geom(k)[0] and geom(k)[1] and geom(k)[2]>=X),
           "genomic-architecture axis (coords+strand only, NO similarity)")

# Narrow biotype rule: protein x SMALL structural RNA (exclude lncRNA + pseudogene)
SMALL={"snRNA","snoRNA","tRNA","rRNA","miRNA","misc_RNA","ncRNA"}
def r1narrow(k):
    a,b=tuple(k); ba,bb=biotype.get(a),biotype.get(b)
    return (ba=="protein_coding" and bb in SMALL) or (bb=="protein_coding" and ba in SMALL)
report("R1narrow protein x SMALL-structural-RNA (no lncRNA/pseudo)", r1narrow, "targets tRNA/snRNA-insertion FP")

# lncRNA-only mix (where collateral lived)
def r1lnc(k):
    a,b=tuple(k); ba,bb=biotype.get(a),biotype.get(b)
    return (ba=="protein_coding" and bb=="lncRNA") or (bb=="protein_coding" and ba=="lncRNA")
report("R1lnc protein x lncRNA", r1lnc, "for contrast (expect real collateral = lncRNA paralogs)")

# Composition of NONCODING_ELEMENT class
print("\n\n##### NONCODING_ELEMENT class composition (genuine FP) #####")
nc=[k for k in genuine_fp if fp_type[k]=="NONCODING_ELEMENT"]
comp=Counter()
for k in nc:
    a,b=tuple(k); bs=tuple(sorted((biotype.get(a,'NA'),biotype.get(b,'NA'))))
    comp[bs]+=1
for bs,c in comp.most_common(): print(f"  {bs}: {c}")

# Combined coverage: how much of genuine-FP do the CLEAN rules (R2ref.5 + R1narrow) cover, vs shipped-miss
r2set=set(k for k in genuine_fp if (geom(k) and geom(k)[0] and geom(k)[1] and geom(k)[2]>=0.5))
r1nset=set(k for k in genuine_fp if r1narrow(k))
print(f"\n##### CLEAN-rule coverage of genuine FP ({len(genuine_fp)}): R2ref.5={len(r2set)} R1narrow={len(r1nset)} union={len(r2set|r1nset)} #####")
print("  R2ref.5 catches:", " ; ".join(f"{tuple(k)[0]}+{tuple(k)[1]}[{fp_type[k]}]" for k in r2set))
print("  R1narrow catches:", " ; ".join(f"{tuple(k)[0]}+{tuple(k)[1]}[{fp_type[k]}]" for k in r1nset))

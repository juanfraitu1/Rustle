#!/usr/bin/env python3
"""How far down in identity can the homology edge actually reach, and is asm20 the binding constraint?

The E_r edge uses `minimap2 -x asm20` plus a sensitive tier (`-k 11 -w 5`, identity >= 0.70). Families that
failed the known-gene benchmark (S100A, SERPINA, CCL, MMP, SIGLEC) showed ~0 member pairs at >= 0.80 identity
under asm20 -- but that measurement CONFLATES two very different failures:
   (a) the paralogs really are too divergent, or
   (b) they are alignable, and asm20's k=19/w=10 seeding simply cannot find the anchors.
Only (b) is fixable by changing the preset. This measures pairs FOUND and their identity per preset.
"""
import subprocess, sys
from collections import defaultdict
O="/home/juanfra/winloci_scratch/knownfam"; MM2="/home/juanfra/miniforge3/bin/minimap2"
PRESETS=[("asm20",["-x","asm20"]),
         ("sensitive k11w5",["-k","11","-w","5"]),
         ("k9w3",["-k","9","-w","3"]),
         ("sr",["-x","sr"]),
         ("map-ont",["-x","map-ont"])]
seqs=defaultdict(dict); name=None
for ln in open(f"{O}/queries.fa"):
    if ln.startswith(">"):
        root,sym=ln[1:].strip().split("|"); name=(root,sym); seqs[root][sym]=[]
    else: seqs[name[0]][name[1]].append(ln.strip())
FAMS=["S100A","SERPINA","CCL","MMP","SIGLEC","H2BC","KRTAP","GSTM","PCDHB","TUBA"]
print(f"{'family':9s} {'n':>3} {'pairs':>6} " + " ".join(f"{p[0]:>17}" for p in PRESETS))
print(f"{'':9s} {'':>3} {'':>6} " + " ".join(f"{'aligned/>=0.8':>17}" for _ in PRESETS))
for root in FAMS:
    d={k:"".join(v) for k,v in seqs.get(root,{}).items()}
    if len(d)<2: continue
    fa=f"{O}/_ps_{root}.fa"
    with open(fa,"w") as fh:
        for k,v in d.items(): fh.write(f">{k}\n{v}\n")
    npair=len(d)*(len(d)-1)//2
    cells=[]
    for label,args in PRESETS:
        paf=subprocess.run([MM2,"-c","--eqx","-X",*args,"-N","50","-p","0.01",fa,fa],
                           capture_output=True,text=True).stdout
        best=defaultdict(float)
        for ln in paf.splitlines():
            f=ln.split("\t")
            if len(f)<12: continue
            a,b=sorted([f[0],f[5]]); al=int(f[10])
            if al: best[(a,b)]=max(best[(a,b)],int(f[9])/al)
        hi=sum(1 for v in best.values() if v>=0.80)
        cells.append(f"{len(best):>7}/{hi:<9}")
    print(f"{root:9s} {len(d):>3} {npair:>6} " + " ".join(cells))
print("\n  cell = pairs with ANY alignment / of those, pairs at >= 0.80 identity")

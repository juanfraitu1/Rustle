#!/usr/bin/env python3
"""IN-SCOPE (T2T genome + long RNA reads only) test for "this unspliced cluster is a real EXON".

Section 15 tested only NEGATIVE evidence -- reads splicing OVER a cluster -- and reached 0.53 precision,
because in an SD a read skipping the interval may belong to a paralog. This tests POSITIVE evidence instead,
which is not symmetric: a spliced junction whose acceptor lands ON the cluster's start, or whose donor lands
ON its end, is direct observed proof that this interval is retained in a mature transcript. Pre-mRNA in an
intron has no reason to produce junctions at its own boundaries.

Signals, all from the BAM + CHM13 FASTA:
  ACCEPTOR  reads whose intron ENDS   within TOL of the cluster start
  DONOR     reads whose intron STARTS within TOL of the cluster end
  GT/AG     the genome's dinucleotides at those junctions (canonical splice motif)

Label (annotation, grading only): intronic = 0 annotated exonic bases in the cluster.
"""
import subprocess, gzip, re, bisect, sys
from collections import defaultdict
import numpy as np

SAM="/home/juanfra/miniforge3/bin/samtools"
BAM="/home/juanfra/winloci_scratch/soto_cache/soto_regions.bam"
FA ="/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa"
GFF="/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz"
CAT="/home/juanfra/winloci_scratch/soto_cache/definitive.copies.tsv"
BED="bench/soto/80_fams.gene_preferred.bed"
TOL=25

members=[]
for ln in open(BED):
    c,s,e,name,*_=ln.rstrip("\n").split("\t"); members.append((c,int(s),int(e)))
exons=defaultdict(list)
with gzip.open(GFF,"rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f=ln.rstrip("\n").split("\t")
        if len(f)<9 or f[2]!="exon": continue
        exons[f[0]].append((int(f[3])-1,int(f[4])))
for c in exons:
    iv=sorted(exons[c]); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    exons[c]=o
def exbp(c,s,e):
    iv=exons.get(c,[]); i=bisect.bisect_left([x[1] for x in iv],s); t=0
    while i<len(iv) and iv[i][0]<e:
        t+=max(0,min(e,iv[i][1])-max(s,iv[i][0])); i+=1
    return t

clusters=[]
for i,ln in enumerate(open(CAT)):
    if i==0: continue
    f=ln.rstrip("\n").split("\t")
    if len(f)<9 or int(f[6])!=1: continue
    c,s,e=f[3],int(f[4]),int(f[5])
    if any(tc==c and ts<e and s<te for tc,ts,te in members): clusters.append((c,s,e))

def junctions(c,s,e):
    """introns observed in reads overlapping a window around the cluster"""
    lo,hi=max(0,s-2000),e+2000
    out=subprocess.run([SAM,"view","-F","2308",BAM,f"{c}:{lo+1}-{hi}"],capture_output=True,text=True).stdout
    js=defaultdict(int)
    for ln in out.splitlines():
        f=ln.split("\t")
        if len(f)<6 or "N" not in f[5]: continue
        ref=int(f[3])-1
        for num,op in re.findall(r"(\d+)([MIDNSHP=X])",f[5]):
            num=int(num)
            if op=="N": js[(ref,ref+num)]+=1; ref+=num
            elif op in "M=XD": ref+=num
    return js

rows=[]
for c,s,e in clusters:
    js=junctions(c,s,e)
    acc=sum(n for (a,b),n in js.items() if abs(b-s)<=TOL)     # intron ends at cluster start
    don=sum(n for (a,b),n in js.items() if abs(a-e)<=TOL)     # intron starts at cluster end
    rows.append((exbp(c,s,e)==0, acc, don, c, s, e))

P=sum(1 for r in rows if r[0]); n=len(rows)
print(f"unspliced clusters inside a Soto member: {n}")
print(f"  label INTRONIC {P}  EXONIC {n-P}   base rate(exonic) {(n-P)/n:.3f}")
print(f"\nPOSITIVE rule: junction abuts the cluster boundary (TOL={TOL}bp) => predict EXONIC")
print(f"  {'rule':<34} {'TP':>5} {'FP':>5} {'FN':>5} {'precision':>10} {'recall':>7} {'F1':>6}")
def ev(name,fn):
    tp=sum(1 for l,a,d,*_ in rows if not l and fn(a,d))
    fp=sum(1 for l,a,d,*_ in rows if     l and fn(a,d))
    fn_=sum(1 for l,a,d,*_ in rows if not l and not fn(a,d))
    pr=tp/max(tp+fp,1); rc=tp/max(tp+fn_,1)
    print(f"  {name:<34} {tp:>5} {fp:>5} {fn_:>5} {pr:>10.3f} {rc:>7.3f} {2*pr*rc/max(pr+rc,1e-9):>6.3f}")
ev("acceptor>=1 OR donor>=1", lambda a,d: a>=1 or d>=1)
ev("acceptor>=1 AND donor>=1", lambda a,d: a>=1 and d>=1)
ev("acceptor>=3 OR donor>=3", lambda a,d: a>=3 or d>=3)
ev("acceptor>=3 AND donor>=3", lambda a,d: a>=3 and d>=3)
ev("acceptor+donor >= 10", lambda a,d: a+d>=10)

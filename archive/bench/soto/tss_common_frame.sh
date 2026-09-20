#!/usr/bin/env bash
# COMMON-FRAME TSS TEST. Question: do sibling copies of a family use DIFFERENT transcription start sites,
# or the same homologous TSS? Answering needs a shared coordinate frame — raw offsets within each copy are
# not comparable because the copies differ in length/strand.
#
# Method: for each family with >=2 well-covered member copies, align copy B's genomic span onto copy A's
# (minimap2 -a), build a B->A position map from the CIGAR, project every read's 5' end into A's frame, then
# compare the TSS distributions. Between-copy median shift is judged against the WITHIN-copy spread (MAD),
# so "different TSS" means the shift exceeds the noise the k=2 boundary quantile is designed to absorb.
set -uo pipefail
SAM=/home/juanfra/miniforge3/bin/samtools
MM2=/home/juanfra/miniforge3/bin/minimap2
FA=/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa
BAM=/home/juanfra/winloci_scratch/soto_cache/soto_regions.bam
BED=/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/80_fams.chr.bed
O=/home/juanfra/winloci_scratch/tss_frame; rm -rf "$O"; mkdir -p "$O"

python3 - "$SAM" "$MM2" "$FA" "$BAM" "$BED" "$O" <<'PY'
import sys, subprocess, re, statistics
from collections import defaultdict, Counter
SAM,MM2,FA,BAM,BED,O = sys.argv[1:7]

# --- member loci per family
fams=defaultdict(list)
for ln in open(BED):
    c,s,e,name,*_=ln.rstrip("\n").split("\t")
    fams[name.split("|")[1]].append((name.split("|")[0],c,int(s),int(e)))

def five_prime_ends(c,s,e):
    """reads' 5' ends (genomic coord) for primary reads inside the locus, with strand handled."""
    out=subprocess.run([SAM,"view","-F","2308",BAM,f"{c}:{s}-{e}"],capture_output=True,text=True).stdout
    ends=[]
    for ln in out.splitlines():
        f=ln.split("\t")
        flag=int(f[1]); pos=int(f[3]); cig=f[5]
        span=sum(int(n) for n,op in re.findall(r"(\d+)([MIDNSHP=X])",cig) if op in "MDN=X")
        five = pos+span if (flag & 16) else pos
        if s<=five<=e: ends.append(five)
    return ends

def fetch(c,s,e,path):
    subprocess.run(f"{SAM} faidx {FA} {c}:{s+1}-{e} > {path}", shell=True, capture_output=True)

def bmap(a,b,tag):
    """align B onto A; return function mapping a B genomic coord -> A genomic coord (or None)."""
    (an,ac,as_,ae)=a; (bn,bc,bs,be)=b
    fa=f"{O}/{tag}_A.fa"; fb=f"{O}/{tag}_B.fa"
    fetch(ac,as_,ae,fa); fetch(bc,bs,be,fb)
    sam=subprocess.run([MM2,"-a","-x","asm20","--eqx",fa,fb],capture_output=True,text=True).stdout
    best=None
    for ln in sam.splitlines():
        if ln.startswith("@"): continue
        f=ln.split("\t")
        if len(f)<6 or f[2]=="*" or (int(f[1]) & 0x900): continue
        nm=sum(int(n) for n,op in re.findall(r"(\d+)([MIDN=X])",f[5]) if op in "M=")
        if best is None or nm>best[0]: best=(nm,f)
    if best is None: return None
    f=best[1]; rev=bool(int(f[1]) & 16); apos=int(f[3])-1; bpos=0
    pairs=[]                       # (b_local, a_local) anchors at match blocks
    for n,op in re.findall(r"(\d+)([MIDNSHP=X])",f[5]):
        n=int(n)
        if op in "M=X":
            for k in range(0,n,50): pairs.append((bpos+k, apos+k))
            bpos+=n; apos+=n
        elif op in "IS": bpos+=n
        elif op in "DN": apos+=n
    if not pairs: return None
    pairs.sort()
    blocal=[p[0] for p in pairs]; alocal=[p[1] for p in pairs]
    blen=be-bs
    def f_map(bgen):
        bl = (blen-1-(bgen-bs)) if rev else (bgen-bs)   # B local, flipped if the alignment is reverse
        import bisect
        i=bisect.bisect_left(blocal,bl)
        if i>=len(blocal): i=len(blocal)-1
        return as_+alocal[i]
    return f_map

def mad(xs):
    if len(xs)<2: return 0.0
    m=statistics.median(xs); return statistics.median(abs(x-m) for x in xs)

print(f"{'family':8s} {'copyA':12s} {'copyB':12s} {'nA':>5} {'nB':>5} {'shift(bp)':>10} {'withinMAD':>10} {'ratio':>7}  verdict")
tested=0; distinct=0; same=0
for fid,mem in sorted(fams.items()):
    cov=[(nm,c,s,e,five_prime_ends(c,s,e)) for (nm,c,s,e) in mem]
    cov=[x for x in cov if len(x[4])>=40]                 # need real depth on both copies
    if len(cov)<2: continue
    cov.sort(key=lambda x:-len(x[4]))
    (an,ac,as_,ae,aends)=cov[0]; (bn,bc,bs,be,bends)=cov[1]
    fm=bmap((an,ac,as_,ae),(bn,bc,bs,be),fid)
    if fm is None: continue
    bproj=[fm(x) for x in bends]
    shift=abs(statistics.median(aends)-statistics.median(bproj))
    within=max(mad(aends), mad(bproj), 1.0)
    ratio=shift/within
    tested+=1
    verdict = "DISTINCT TSS" if ratio>=3 else ("same TSS" if ratio<=1 else "ambiguous")
    if ratio>=3: distinct+=1
    elif ratio<=1: same+=1
    print(f"{fid:8s} {an:12s} {bn:12s} {len(aends):>5} {len(bproj):>5} {shift:>10.0f} {within:>10.0f} {ratio:>7.1f}  {verdict}")
print(f"\ntested {tested} families with 2 well-covered copies: DISTINCT {distinct}, same {same}, ambiguous {tested-distinct-same}")
print("(shift = |median 5' end of A  -  median 5' end of B projected into A's frame|;")
print(" withinMAD = larger of the two copies' own 5'-end median-absolute-deviations)")
PY

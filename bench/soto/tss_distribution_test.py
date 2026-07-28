#!/usr/bin/env python3
"""Distributional TSS test: is the exon-sum enough, or should it be EXTENDED to encode TSS/TES?

The earlier median-shift test (bam_tie_signals.md §7) is weak: a BIMODAL TSS used by both copies, or copies
using the SAME peaks in DIFFERENT PROPORTIONS, produces no median shift yet is genuinely discriminating.
This compares the full 5'-end DISTRIBUTIONS in a common coordinate frame, against an EMPIRICAL NULL.

Method per family (>=2 copies with >=40 primary reads):
  1. Align copy B's genomic span onto copy A's (minimap2 -a -x asm20); build a reverse-aware B->A map.
  2. REJECT pairs whose alignment covers <50% of B — those are non-equivalent copies (different genes
     sharing a duplicated block), which produced 5/7 spurious "distinct" calls in the median test.
  3. Project B's read 5' ends into A's frame.
  4. observed = Wasserstein-1(A_5ends, Bproj_5ends).
  5. NULL = Wasserstein-1 between two random halves of the SAME copy, resampled `NPERM` times (captures
     sampling noise AND the copy's own 5'-end heterogeneity/degradation). Pool A's and B's nulls.
  6. Call DISTINCT only if observed > 99th percentile of the null AND exceeds MIN_BP (biological floor —
     a statistically clean 9 bp shift is not worth representing).

Usage: tss_distribution_test.py [--nperm 200]
"""
import subprocess, re, sys, json
import numpy as np
from collections import defaultdict
from scipy.stats import wasserstein_distance

SAMTOOLS = "/home/juanfra/miniforge3/bin/samtools"
MINIMAP2 = "/home/juanfra/miniforge3/bin/minimap2"
FASTA    = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa"
BAM      = "/home/juanfra/winloci_scratch/soto_cache/soto_regions.bam"
BED      = "bench/soto/80_fams.chr.bed"
TMP      = "/home/juanfra/winloci_scratch/tss_dist"
MIN_READS   = 40      # per copy
MIN_ALN_FRAC= 0.50    # B must align this fraction to A, else the pair is non-equivalent
MIN_BP      = 50      # biological floor on a meaningful TSS difference
MIN_SHARP   = 0.30    # both copies must have >=30% of 5' ends in one 400bp window to count as a real TSS
NPERM       = int(sys.argv[sys.argv.index("--nperm")+1]) if "--nperm" in sys.argv else 200

subprocess.run(["mkdir","-p",TMP], check=False)
rng = np.random.default_rng(0)

def five_prime_ends(chrom, start, end):
    """5' genomic coordinate of each primary read within the locus (strand-aware)."""
    out = subprocess.run([SAMTOOLS,"view","-F","2308",BAM,f"{chrom}:{start}-{end}"],
                         capture_output=True, text=True).stdout
    ends=[]
    for ln in out.splitlines():
        f=ln.split("\t")
        flag=int(f[1]); pos=int(f[3])
        span=sum(int(n) for n,op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]) if op in "MDN=X")
        five = pos+span if (flag & 16) else pos
        if start <= five <= end: ends.append(five)
    return np.array(ends, dtype=float)

def fetch(chrom, start, end, path):
    subprocess.run(f"{SAMTOOLS} faidx {FASTA} {chrom}:{start+1}-{end} > {path}",
                   shell=True, capture_output=True)

def build_map(a, b, tag):
    """Align B onto A. Returns (map_fn, aligned_fraction_of_B) or (None, 0)."""
    (_an, ac, as_, ae) = a; (_bn, bc, bs, be) = b
    fa, fb = f"{TMP}/{tag}_A.fa", f"{TMP}/{tag}_B.fa"
    fetch(ac, as_, ae, fa); fetch(bc, bs, be, fb)
    sam = subprocess.run([MINIMAP2,"-a","-x","asm20","--eqx",fa,fb],
                         capture_output=True, text=True).stdout
    best=None
    for ln in sam.splitlines():
        if ln.startswith("@"): continue
        f=ln.split("\t")
        if len(f)<6 or f[2]=="*" or (int(f[1]) & 0x900): continue
        matched=sum(int(n) for n,op in re.findall(r"(\d+)([MIDN=X])", f[5]) if op in "M=")
        if best is None or matched>best[0]: best=(matched,f)
    if best is None: return None, 0.0
    matched, f = best
    blen = be-bs
    aln_frac = matched/max(blen,1)
    rev = bool(int(f[1]) & 16); apos=int(f[3])-1; bpos=0
    b_anchor=[]; a_anchor=[]
    for n,op in re.findall(r"(\d+)([MIDNSHP=X])", f[5]):
        n=int(n)
        if op in "M=X":
            for k in range(0,n,20): b_anchor.append(bpos+k); a_anchor.append(apos+k)
            bpos+=n; apos+=n
        elif op in "IS": bpos+=n
        elif op in "DN": apos+=n
    if not b_anchor: return None, aln_frac
    b_anchor=np.array(b_anchor); a_anchor=np.array(a_anchor)
    order=np.argsort(b_anchor); b_anchor=b_anchor[order]; a_anchor=a_anchor[order]
    def fmap(bgen):
        bl = (blen-1-(bgen-bs)) if rev else (bgen-bs)
        i = np.clip(np.searchsorted(b_anchor, bl), 0, len(b_anchor)-1)
        return as_ + a_anchor[i]
    return fmap, aln_frac

def sharpness(v):
    """Fraction of reads whose 5' end falls in the densest 400 bp window. A real TSS gives a sharp peak;
    broad scatter means the '5' end' is just wherever coverage fell — differential partial COVERAGE, not a
    promoter. Without this control the EMD test calls 28/40 pairs DISTINCT, of which 25 are coverage artifacts."""
    if len(v) < 10: return 0.0
    v = np.sort(v); best = 0
    for x in v[::max(1, len(v)//200)]:
        best = max(best, int(np.sum((v >= x-200) & (v <= x+200))))
    return best/len(v)

def null_emd(x, nperm):
    """Within-copy null: EMD between two random halves of the same sample."""
    if len(x) < 8: return None
    vals=[]
    for _ in range(nperm):
        p = rng.permutation(len(x)); h = len(x)//2
        vals.append(wasserstein_distance(x[p[:h]], x[p[h:2*h]]))
    return np.array(vals)

fams=defaultdict(list)
for ln in open(BED):
    c,s,e,name,*_ = ln.rstrip("\n").split("\t")
    fams[name.split("|")[1]].append((name.split("|")[0], c, int(s), int(e)))

print(f"{'family':8s} {'copyA':12s} {'copyB':12s} {'nA':>5} {'nB':>5} {'alnB':>5} "
      f"{'EMD_obs':>9} {'null_p99':>9} {'ratio':>7} {'shpA':>6} {'shpB':>6}  verdict")
res=[]
for fid, mem in sorted(fams.items()):
    cov=[(nm,c,s,e,five_prime_ends(c,s,e)) for (nm,c,s,e) in mem]
    cov=[x for x in cov if len(x[4])>=MIN_READS]
    if len(cov)<2: continue
    cov.sort(key=lambda x: -len(x[4]))
    (an,ac,as_,ae,A)=cov[0]; (bn,bc,bs,be,B)=cov[1]
    fmap, aln = build_map((an,ac,as_,ae),(bn,bc,bs,be), fid)
    if fmap is None:
        continue
    if aln < MIN_ALN_FRAC:
        print(f"{fid:8s} {an:12s} {bn:12s} {len(A):>5} {len(B):>5} {aln:>5.2f} "
              f"{'-':>9} {'-':>9} {'-':>7}  SKIP non-equivalent pair")
        res.append((fid,an,bn,None,None,None,"skip")); continue
    Bp = np.array([fmap(x) for x in B], dtype=float)
    obs = wasserstein_distance(A, Bp)
    nulls=[n for n in (null_emd(A,NPERM), null_emd(Bp,NPERM)) if n is not None]
    if not nulls:
        continue
    null = np.concatenate(nulls); p99 = float(np.percentile(null, 99))
    ratio = obs/max(p99, 1e-9)
    sa, sb = sharpness(A), sharpness(Bp)
    sharp = (sa >= MIN_SHARP) and (sb >= MIN_SHARP)
    distinct = (obs > p99) and (obs >= MIN_BP) and sharp
    if obs <= p99:                       verdict = "same TSS"
    elif obs < MIN_BP:                   verdict = f"significant but <{MIN_BP}bp"
    elif not sharp:                      verdict = "BROAD (coverage, not TSS)"
    else:                                verdict = "DISTINCT TSS"
    print(f"{fid:8s} {an:12s} {bn:12s} {len(A):>5} {len(B):>5} {aln:>5.2f} "
          f"{obs:>9.0f} {p99:>9.0f} {ratio:>7.1f} {sa:>6.2f} {sb:>6.2f}  {verdict}")
    res.append((fid,an,bn,obs,p99,ratio,verdict))

ok=[r for r in res if r[6]!="skip"]
d=sum(1 for r in ok if r[6]=="DISTINCT TSS")
s_=sum(1 for r in ok if r[6]=="same TSS")
t=sum(1 for r in ok if r[6].startswith("significant"))
sk=sum(1 for r in res if r[6]=="skip")
print(f"\ntested {len(ok)} homologous copy pairs ({sk} skipped as non-equivalent, alnB<{MIN_ALN_FRAC}):")
print(f"  DISTINCT TSS (EMD > null p99 AND >= {MIN_BP} bp): {d}")
print(f"  same TSS (EMD within null):                      {s_}")
print(f"  significant but sub-{MIN_BP}bp (not worth encoding): {t}")
print(f"  BROAD 5' scatter (differential coverage, NOT a promoter): {b_}")
print(f"\nnull = EMD between random halves of the same copy ({NPERM} perms) — captures sampling noise")
print("and each copy's own 5'-end heterogeneity (degradation), which is what the k=2 boundary quantile absorbs.")
json.dump([{"fam":r[0],"A":r[1],"B":r[2],"emd":r[3],"p99":r[4],"ratio":r[5],"verdict":r[6]} for r in res],
          open(f"{TMP}/result.json","w"), indent=1)
print(f"\nwrote {TMP}/result.json")

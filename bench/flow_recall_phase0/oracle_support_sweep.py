"""Read-support sweep for the 1-edit near-miss types. For each edit, measure the junction
read support (from the BAM) of the edit junction vs the transcript body, to see whether a rule
can fire the edit WITHOUT mis-firing on the matched refs.

Reuses the gffcompare tmaps already in /tmp/oracle_nm/<chrom>/ (run oracle_nearmiss.py first).

Per edit type:
  k_trim   : ratio = extra_junc_reads / median_body_reads. Low => droppable. Precision check:
             how often does a MATCHED ref carry a real junction with the same low ratio (would be
             wrongly dropped)?
  c_extend : the MISSING junction's read support. If >=K reads use it, it's recoverable (rustle
             skipped a supported junction); if ~0, genuinely absent.
  altsplice: ref_junc_reads vs rustle_junc_reads at the shared boundary. Rule = prefer higher
             support; report fraction where ref > rustle (rule correct).
"""
import os, sys, glob, re, subprocess
from collections import Counter, defaultdict
from statistics import median
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib

CHROMS = ["NC_073227.2", "NC_073231.2", "NC_073235.2", "NC_073239.2"]
_CIG = re.compile(r'(\d+)([MIDNSHP=X])')

def junction_counts(chrom):
    """(intron_first, intron_last) 1-based inclusive -> # reads spanning that junction."""
    jc = Counter()
    bam = f"{lib.PERCHROM}/{chrom}/c.bam"
    p = subprocess.Popen(["samtools", "view", "-F", "0x100", bam, chrom],
                         stdout=subprocess.PIPE, text=True)
    for line in p.stdout:
        f = line.split("\t", 6)
        if len(f) < 6:
            continue
        pos = int(f[3]); cig = f[5]
        if "N" not in cig:
            continue
        ref = pos
        for n, op in _CIG.findall(cig):
            n = int(n)
            if op in "MD=X":
                ref += n
            elif op == "N":
                jc[(ref, ref + n - 1)] += 1   # intron 1-based inclusive
                ref += n
    p.wait()
    return jc

def chains(path):
    tx = {}
    for line in open(path):
        f = line.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if m:
            tx.setdefault(m.group(1), []).append((int(f[3]), int(f[4])))
    out = {}
    for k, ex in tx.items():
        ex.sort()
        out[k] = tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
    return out

def main():
    ktrim_ratio = []; cext_support = []; alt = []  # alt: (ref_reads, ru_reads)
    matched_min_ratio = []  # for precision check: min junc/body ratio among matched refs
    for c in CHROMS:
        wd = f"/tmp/oracle_nm/{c}"
        tmaps = glob.glob(f"{wd}/g.*.tmap")
        if not tmaps:
            continue
        ref_ch = chains(f"{lib.PERCHROM}/{c}/ref.gtf")
        ru_ch = chains(f"/tmp/strand_safety/all/{c}/off.gtf")
        jc = junction_counts(c)
        def cnt(j):
            return jc.get(j, 0)
        ref_hits = defaultdict(list); matched = set()
        for line in open(tmaps[0]):
            f = line.rstrip("\n").split("\t")
            if f[0] == "ref_gene_id":
                continue
            if f[2] == "=":
                matched.add(f[1])
            elif f[2] in ("j", "k", "c", "m", "n"):
                ref_hits[f[1]].append(f[4])
        # precision baseline: min junction/body ratio among matched refs (real junctions)
        for ref_id in matched:
            R = ref_ch.get(ref_id)
            if not R or len(R) < 2:
                continue
            counts = [cnt(j) for j in R]
            med = median(counts) or 1
            matched_min_ratio.append(min(counts) / med)
        # edit cases
        for ref_id, qrys in ref_hits.items():
            if ref_id in matched or ref_id not in ref_ch:
                continue
            R = set(ref_ch[ref_id])
            best = None; bs = -1
            for q in qrys:
                Q = ru_ch.get(q)
                if not Q:
                    continue
                sh = len(R & set(Q))
                if sh > bs:
                    bs = sh; best = set(Q)
            if not best or bs <= 0:
                continue
            re_, qe = R - best, best - R
            body = [cnt(j) for j in (R & best)]
            med = median(body) if body else 1
            med = med or 1
            if not qe and len(re_) == 1:                       # c_extend (1 missing)
                cext_support.append(cnt(next(iter(re_))))
            elif not re_ and len(qe) == 1:                     # k_trim (1 extra)
                ktrim_ratio.append(cnt(next(iter(qe))) / med)
            elif len(re_) == 1 and len(qe) == 1:               # altsplice
                ra, qa = next(iter(re_)), next(iter(qe))
                if ra[0] == qa[0] or ra[1] == qa[1]:
                    alt.append((cnt(ra), cnt(qa)))

    print(f"=== read-support sweep (baseline, {len(CHROMS)} chroms) ===\n")
    if ktrim_ratio:
        ks = sorted(ktrim_ratio)
        lo = sum(1 for r in ks if r < 0.1)
        print(f"k_trim (n={len(ks)}): extra-junc reads / median-body reads")
        print(f"  median ratio={ks[len(ks)//2]:.3f}; <0.1: {lo} ({100*lo/len(ks):.0f}%)  <- droppable-looking")
        mm = sorted(matched_min_ratio)
        mlo = sum(1 for r in mm if r < 0.1)
        print(f"  PRECISION CHECK: matched refs whose WEAKEST real junction has ratio <0.1: "
              f"{mlo}/{len(mm)} ({100*mlo/max(1,len(mm)):.0f}%)  <- a <0.1 drop-rule would break these")
    if cext_support:
        cs = sorted(cext_support)
        sup = sum(1 for x in cs if x >= 2)
        print(f"\nc_extend (n={len(cs)}): read support of the MISSING junction rustle skipped")
        print(f"  median={cs[len(cs)//2]}; with >=2 reads (supported, recoverable): {sup} ({100*sup/len(cs):.0f}%)")
        print(f"  with 0 reads (genuinely absent): {sum(1 for x in cs if x==0)} ({100*sum(1 for x in cs if x==0)/len(cs):.0f}%)")
    if alt:
        ref_higher = sum(1 for r, q in alt if r > q)
        tie = sum(1 for r, q in alt if r == q)
        print(f"\naltsplice (n={len(alt)}): ref junction reads vs rustle's shifted-pick reads")
        print(f"  ref HIGHER support (rule 'prefer ref' correct): {ref_higher} ({100*ref_higher/len(alt):.0f}%)")
        print(f"  rustle higher (rule WRONG): {len(alt)-ref_higher-tie} ; tie: {tie}")

if __name__ == "__main__":
    main()

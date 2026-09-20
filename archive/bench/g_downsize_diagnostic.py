"""Level-1 -G diagnostic: from the EXISTING 200k BAM, count primary alignments whose CIGAR
carries an intron (N op) > 50k / > 100k, and classify each as a REAL large intron (both flanks
in ONE annotated gene) vs an inter-copy CHIMERA candidate (flanks in different genes = a
read-through a smaller -G would break). No re-alignment."""
import pysam, bisect, collections, json, os
BAM="/home/juanfra/winloci_scratch/GGO_mm.bam"
GENES="/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
T1,T2=50000,100000
# gene intervals per chrom
starts=collections.defaultdict(list); ivs=collections.defaultdict(list)
for ln in open(GENES):
    p=ln.rstrip("\n").split("\t")
    if len(p)<4: continue
    ivs[p[0]].append((int(p[1]),int(p[2]),p[3]))
for c in ivs:
    ivs[c].sort(); starts[c]=[x[0] for x in ivs[c]]
def genes_at(c,pos):
    """annotated gene names overlapping pos on chrom c"""
    out=set(); L=ivs.get(c);
    if not L: return out
    i=bisect.bisect_right(starts[c],pos)
    for j in range(max(0,i-40),min(len(L),i+1)):
        s,e,n=L[j]
        if s<=pos<=e: out.add(n)
    return out
N_OP=3
cnt=collections.Counter()
examples=[]
bam=pysam.AlignmentFile(BAM,"rb")
nread=0
for a in bam.fetch(until_eof=True):
    if a.is_secondary or a.is_supplementary or a.is_unmapped: continue
    nread+=1
    ref=a.reference_start; c=a.reference_name; big=[]
    for op,ln in a.cigartuples or []:
        if op in (0,2,7,8): ref+=ln          # M/D/=/X consume ref
        elif op==N_OP:
            if ln>T1: big.append((ref,ref+ln,ln))
            ref+=ln
        # I/S/H don't consume ref
    for (istart,iend,ilen) in big:
        thr = ">100k" if ilen>T2 else "50-100k"
        gup=genes_at(c,istart-1); gdn=genes_at(c,iend+1)
        if gup and gdn and (gup & gdn):      cls="real_large_intron"
        elif gup and gdn and not (gup&gdn):  cls="chimera_diff_genes"
        else:                                cls="chimera_or_unannot"
        cnt[(thr,cls)]+=1
        if len(examples)<25 and ilen>T2 and cls=="chimera_diff_genes":
            examples.append(dict(read=a.query_name,chrom=c,istart=istart,iend=iend,ilen=ilen,up=list(gup)[:2],dn=list(gdn)[:2]))
bam.close()
res=dict(n_primary=nread, counts={f"{t}|{c}":n for (t,c),n in sorted(cnt.items())}, examples=examples)
json.dump(res,open("g_downsize_diagnostic.json","w"),indent=1)
print(f"primary alignments scanned: {nread}")
for (t,c),n in sorted(cnt.items()): print(f"  intron {t:8} {c:20} : {n}")
tot100=sum(n for (t,c),n in cnt.items() if t==">100k")
chim100=sum(n for (t,c),n in cnt.items() if t==">100k" and c.startswith("chimera"))
tot50=sum(n for (t,c),n in cnt.items() if t=="50-100k")
print(f"\n  >100k introns: {tot100} total, {chim100} chimera-candidate ({100*chim100/max(1,tot100):.0f}%) <- what -G 100k would break")
print(f"  50-100k introns: {tot50} total <- additionally broken by -G 50k")

#!/usr/bin/env python3
"""Size the PRIMARY-read PSV-phasing opportunity: at each paralog-family copy locus on the
2 chroms, do the primary reads split into >=2 PSV-linked haplotype groups that DIFFER in
junction usage? That = allele/copy-specific splicing the merged (StringTie/primary-flow)
assembly collapses, recoverable by PSV-phasing of primary reads (no secondaries)."""
import pysam, csv, collections
BAM="/home/juanfra/winloci_scratch/GGO.bam"
EX="/tmp/cre_guided/exons.tsv"
CHROMS={"NC_073235.2","NC_086018.1"}
MINDP=12; ALT_FRAC=0.25; MIN_HET=2; MIN_GRP=5

# candidate loci = annotated copy spans on the 2 chroms (dedup overlapping)
spans=[]
for r in csv.DictReader(open(EX),delimiter="\t"):
    if r["chrom"] not in CHROMS: continue
    s=[int(x) for x in r["exon_starts"].split(",") if x]; e=[int(x) for x in r["exon_ends"].split(",") if x]
    if not s: continue
    spans.append((r["chrom"], s[0], e[-1]))
spans.sort()
loci=[]
for ch,s,e in spans:
    if loci and loci[-1][0]==ch and s<=loci[-1][2]:
        loci[-1]=(ch, loci[-1][1], max(loci[-1][2],e))
    else: loci.append((ch,s,e))

bam=pysam.AlignmentFile(BAM,"rb")
n_loci=n_phasable=n_struct=0
examples=[]
for ch,s,e in loci:
    if e-s>200000: continue
    # primary reads only, collect per-read bases (for het) + introns (for structure)
    rbases=collections.defaultdict(dict); rintr=collections.defaultdict(set); cov=collections.Counter(); alt=collections.defaultdict(collections.Counter)
    nreads=0
    for read in bam.fetch(ch,s,e):
        if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
        nreads+=1
        bl=read.get_blocks()
        for i in range(len(bl)-1): rintr[read.query_name].add((bl[i][1],bl[i+1][0]))
        for qpos,rpos in read.get_aligned_pairs(matches_only=True):
            if rpos is None or qpos is None or not (s<=rpos<e): continue
            b=read.query_sequence[qpos]
            rbases[read.query_name][rpos]=b; cov[rpos]+=1; alt[rpos][b]+=1
    if nreads<MINDP: continue
    n_loci+=1
    # het/PSV sites: 2nd allele >= ALT_FRAC and depth>=MINDP (filters seq errors which are < ALT_FRAC)
    het=[]
    for pos,c in cov.items():
        if c<MINDP: continue
        top=alt[pos].most_common(2)
        if len(top)==2 and top[1][1]>=ALT_FRAC*c: het.append(pos)
    if len(het)<MIN_HET: continue
    # phase reads into 2 groups by their alleles at het sites (majority haplotype split)
    # build the two dominant haplotype patterns
    pat=collections.Counter()
    for rn,bm in rbases.items():
        key=tuple(bm.get(p,'-') for p in het)
        if key.count('-')<=len(het)//2: pat[key]+=1
    top=[k for k,c in pat.most_common(2) if c>=MIN_GRP]
    if len(top)<2: continue
    n_phasable+=1
    # structural divergence: do the 2 haplotype groups differ in junction usage?
    g={0:set(),1:set()}; gi={0:collections.Counter(),1:collections.Counter()}; gn={0:0,1:0}
    for rn,bm in rbases.items():
        key=tuple(bm.get(p,'-') for p in het)
        # assign read to whichever top pattern it matches better
        d0=sum(1 for a,b in zip(key,top[0]) if a!='-' and a==b)
        d1=sum(1 for a,b in zip(key,top[1]) if a!='-' and a==b)
        if d0==d1: continue
        grp=0 if d0>d1 else 1; gn[grp]+=1
        for j in rintr[rn]: gi[grp][j]+=1
    if gn[0]<MIN_GRP or gn[1]<MIN_GRP: continue
    diverged=False
    alljn=set(gi[0])|set(gi[1])
    for j in alljn:
        f0=gi[0][j]/max(1,gn[0]); f1=gi[1][j]/max(1,gn[1])
        if (f0>=0.6 and f1<=0.2) or (f1>=0.6 and f0<=0.2):
            diverged=True; break
    if diverged:
        n_struct+=1
        if len(examples)<12: examples.append((ch,s,e,len(het),gn[0],gn[1]))

print(f"copy loci analyzed (>= {MINDP} primary reads): {n_loci}")
print(f"  with >=2 PSV-linked haplotype groups (co-mingled forms): {n_phasable}")
print(f"  of those, groups DIFFER in junction usage (PSV-linked alt-splicing): {n_struct}")
print(f"  -> primary-read PSV-phasing structural headroom on these 2 chroms: {n_struct}")
print("\nexamples (chrom, start, end, #het, grpA_reads, grpB_reads):")
for x in examples: print(f"  {x[0]}:{x[1]}-{x[2]}  het={x[3]} A={x[4]} B={x[5]}")

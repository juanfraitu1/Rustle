#!/usr/bin/env python3
"""Junction-QC the 69 novel (missed-by-StringTie, not-in-RefSeq) divergent junctions:
canonical splice site (GT-AG / CT-AC), RT-switch direct-repeat signature, and read support.
Credible real novel = canonical + supported + not RT-switch."""
import pysam, csv, collections, re, os
BAM="/home/juanfra/winloci_scratch/GGO.bam"; FA="/home/juanfra/winloci_scratch/GGO.fasta"
EX="/tmp/cre_guided/exons.tsv"; GW="/tmp/gw"
CHROMS={"NC_073235.2","NC_086018.1"}; MINDP=12; ALT_FRAC=0.25; MIN_HET=2; MIN_GRP=5; TOL=5; RPT=8

def introns(path, gff=False):
    tx=collections.defaultdict(list)
    if not os.path.exists(path): return []
    for l in open(path):
        if l.startswith("#"):continue
        f=l.rstrip().split("\t")
        if len(f)<9 or f[2]!="exon":continue
        m=(re.search(r'Parent=([^;]+)',f[8]) if gff else re.search(r'transcript_id "([^"]+)"',f[8]))
        if not m:continue
        tx[m.group(1)].append((int(f[3]),int(f[4])))
    out=[]
    for t,e in tx.items():
        e.sort()
        for i in range(len(e)-1): out.append((e[i][1], e[i+1][0]-1))
    return out
def match(j,lst): return any(abs(j[0]-k[0])<=TOL and abs(j[1]-k[1])<=TOL for k in lst)
ST={ch:introns(f"{GW}/st_{ch}.gtf") for ch in CHROMS}
REF={ch:introns(f"{GW}/ref_{ch}.gff3",gff=True) for ch in CHROMS}

spans=[]
for r in csv.DictReader(open(EX),delimiter="\t"):
    if r["chrom"] not in CHROMS: continue
    s=[int(x) for x in r["exon_starts"].split(",") if x]; e=[int(x) for x in r["exon_ends"].split(",") if x]
    if s: spans.append((r["chrom"],s[0],e[-1]))
spans.sort(); loci=[]
for ch,s,e in spans:
    if loci and loci[-1][0]==ch and s<=loci[-1][2]: loci[-1]=(ch,loci[-1][1],max(loci[-1][2],e))
    else: loci.append((ch,s,e))

bam=pysam.AlignmentFile(BAM,"rb"); fa=pysam.FastaFile(FA)
def dinuc(ch,a,b): return fa.fetch(ch,a,b).upper()
def canon(ch,d,a):  # intron [d,a) 0-based; donor=fa[d:d+2], acceptor=fa[a-2:a]
    don=dinuc(ch,d,d+2); acc=dinuc(ch,a-2,a)
    if (don,acc) in (("GT","AG"),("CT","AC")): return "canonical"
    if (don,acc) in (("GC","AG"),("CT","GC"),("AT","AC"),("GT","AT")): return "semi"
    return "noncanonical"
def rtswitch(ch,d,a):  # direct repeat flanking the junction (RTS signature)
    try:
        up=dinuc(ch,d-RPT,d); down=dinuc(ch,a-RPT,a)
        return up==down and "N" not in up
    except Exception: return False

juncs=[]  # (chrom, d, a, support)
for ch,s,e in loci:
    if e-s>200000: continue
    rbases=collections.defaultdict(dict); rintr=collections.defaultdict(set); cov=collections.Counter(); alt=collections.defaultdict(collections.Counter); nreads=0
    for read in bam.fetch(ch,s,e):
        if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
        nreads+=1; bl=read.get_blocks()
        for i in range(len(bl)-1): rintr[read.query_name].add((bl[i][1],bl[i+1][0]))
        for qpos,rpos in read.get_aligned_pairs(matches_only=True):
            if rpos is None or not (s<=rpos<e): continue
            b=read.query_sequence[qpos]; rbases[read.query_name][rpos]=b; cov[rpos]+=1; alt[rpos][b]+=1
    if nreads<MINDP: continue
    het=[p for p,c in cov.items() if c>=MINDP and len(alt[p].most_common(2))==2 and alt[p].most_common(2)[1][1]>=ALT_FRAC*c]
    if len(het)<MIN_HET: continue
    pat=collections.Counter()
    for rn,bm in rbases.items():
        key=tuple(bm.get(p,'-') for p in het)
        if key.count('-')<=len(het)//2: pat[key]+=1
    top=[k for k,c in pat.most_common(2) if c>=MIN_GRP]
    if len(top)<2: continue
    gi={0:collections.Counter(),1:collections.Counter()}; gn={0:0,1:0}
    for rn,bm in rbases.items():
        key=tuple(bm.get(p,'-') for p in het)
        d0=sum(1 for a,b in zip(key,top[0]) if a!='-' and a==b); d1=sum(1 for a,b in zip(key,top[1]) if a!='-' and a==b)
        if d0==d1: continue
        g=0 if d0>d1 else 1; gn[g]+=1
        for j in rintr[rn]: gi[g][j]+=1
    if gn[0]<MIN_GRP or gn[1]<MIN_GRP: continue
    for j in set(gi[0])|set(gi[1]):
        f0=gi[0][j]/max(1,gn[0]); f1=gi[1][j]/max(1,gn[1])
        if (f0>=0.6 and f1<=0.2) or (f1>=0.6 and f0<=0.2):
            if not match(j,ST[ch]) and not match(j,REF[ch]):
                sup=max(gi[0][j],gi[1][j])
                juncs.append((ch,j[0],j[1],sup))

cat=collections.Counter(); rts=0; sup_ok=0; credible=0; supdist=[]
for ch,d,a,sup in juncs:
    c=canon(ch,d,a); cat[c]+=1; supdist.append(sup)
    r=rtswitch(ch,d,a); rts+= r
    if sup>=5: sup_ok+=1
    if c=="canonical" and sup>=5 and not r: credible+=1
import statistics
print(f"novel divergent junctions QC'd: {len(juncs)}")
print(f"  splice-site: canonical(GT-AG/CT-AC)={cat['canonical']}  semi(GC-AG/AT-AC)={cat['semi']}  noncanonical={cat['noncanonical']}")
print(f"  RT-switch direct-repeat suspect: {rts}")
print(f"  read support >=5: {sup_ok}   (median support={statistics.median(supdist) if supdist else 0})")
print(f"\n  >>> CREDIBLE real novel (canonical + support>=5 + not RT-switch): {credible} / {len(juncs)}")

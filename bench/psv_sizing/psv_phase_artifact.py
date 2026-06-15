#!/usr/bin/env python3
"""Are the 68 non-corroborated divergent junctions cross-mapping artifacts or real novel?
For each, take the reads carrying the divergent junction (primary at THIS locus) and ask:
do they align at least as well to a SISTER family copy (AS_sister >= AS_here)? If most do,
the 'second haplotype' is sister-copy reads cross-mapping in (artifact). If most are
locus-specific (AS_here strictly best), it's a genuine novel form."""
import pysam, csv, collections, re, os
BAM="/home/juanfra/winloci_scratch/GGO.bam"; EX="/tmp/cre_guided/exons.tsv"
UNI="/tmp/cre_guided/universe.tsv"; GW="/tmp/gw"
CHROMS={"NC_073235.2","NC_086018.1"}; MINDP=12; ALT_FRAC=0.25; MIN_HET=2; MIN_GRP=5; TOL=5

def introns_from_gtf(path, gff=False):
    tx=collections.defaultdict(list); st={}
    if not os.path.exists(path): return []
    for l in open(path):
        if l.startswith("#"):continue
        f=l.rstrip().split("\t")
        if len(f)<9 or f[2]!="exon":continue
        m=(re.search(r'Parent=([^;]+)',f[8]) if gff else re.search(r'transcript_id "([^"]+)"',f[8]))
        if not m:continue
        tx[m.group(1)].append((int(f[3]),int(f[4]))); st[m.group(1)]=f[0]
    out=[]
    for t,e in tx.items():
        e.sort()
        for i in range(len(e)-1): out.append((e[i][1], e[i+1][0]-1))
    return out
def match(j,lst): return any(abs(j[0]-k[0])<=TOL and abs(j[1]-k[1])<=TOL for k in lst)
ST={ch:introns_from_gtf(f"{GW}/st_{ch}.gtf") for ch in CHROMS}
REF={ch:introns_from_gtf(f"{GW}/ref_{ch}.gff3",gff=True) for ch in CHROMS}

# exons.tsv: transcript -> (chrom, span); universe: transcript -> family
exspan={}
for r in csv.DictReader(open(EX),delimiter="\t"):
    s=[int(x) for x in r["exon_starts"].split(",") if x]; e=[int(x) for x in r["exon_ends"].split(",") if x]
    if s: exspan[r["transcript_id"]]=(r["chrom"], s[0], e[-1])
fam_of={}; fam_loci=collections.defaultdict(list)
for r in csv.DictReader(open(UNI),delimiter="\t"):
    t=r["transcript_id"]
    if t in exspan: fam_of[t]=r["family_id"]; fam_loci[r["family_id"]].append(exspan[t])
# locus -> family (by overlap with a family copy)
def sisters_of(ch,s,e):
    for fam,locs in fam_loci.items():
        if any(c==ch and ss< e and s< ee for (c,ss,ee) in locs):
            return [(c,ss,ee) for (c,ss,ee) in locs if not (c==ch and ss< e and s< ee)]
    return []

spans=sorted(set(exspan.values()))
loci=[]
for ch,s,e in spans:
    if ch not in CHROMS: continue
    if loci and loci[-1][0]==ch and s<=loci[-1][2]: loci[-1]=(ch,loci[-1][1],max(loci[-1][2],e))
    else: loci.append((ch,s,e))

bam=pysam.AlignmentFile(BAM,"rb")
def best_as_at(reg):
    d={}
    try: it=bam.fetch(reg[0],reg[1],reg[2])
    except (ValueError,KeyError): return d
    for r in it:
        if r.is_unmapped: continue
        try: a=r.get_tag("AS")
        except KeyError: continue
        if r.query_name not in d or a>d[r.query_name]: d[r.query_name]=a
    return d

n68=0; artifact=0; novel=0; examples=[]
for ch,s,e in loci:
    if e-s>200000: continue
    rbases=collections.defaultdict(dict); rintr=collections.defaultdict(set); rAS={}; cov=collections.Counter(); alt=collections.defaultdict(collections.Counter); nreads=0
    for read in bam.fetch(ch,s,e):
        if read.is_unmapped or read.is_secondary or read.is_supplementary: continue
        nreads+=1; bl=read.get_blocks()
        for i in range(len(bl)-1): rintr[read.query_name].add((bl[i][1],bl[i+1][0]))
        try: rAS[read.query_name]=read.get_tag("AS")
        except KeyError: pass
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
    grp={}; gi={0:collections.Counter(),1:collections.Counter()}; gn={0:0,1:0}
    for rn,bm in rbases.items():
        key=tuple(bm.get(p,'-') for p in het)
        d0=sum(1 for a,b in zip(key,top[0]) if a!='-' and a==b); d1=sum(1 for a,b in zip(key,top[1]) if a!='-' and a==b)
        if d0==d1: continue
        g=0 if d0>d1 else 1; grp[rn]=g; gn[g]+=1
        for j in rintr[rn]: gi[g][j]+=1
    if gn[0]<MIN_GRP or gn[1]<MIN_GRP: continue
    # divergent junctions missed-by-ST and NOT-in-ref (the 68 population)
    div=[]
    for j in set(gi[0])|set(gi[1]):
        f0=gi[0][j]/max(1,gn[0]); f1=gi[1][j]/max(1,gn[1])
        if (f0>=0.6 and f1<=0.2) or (f1>=0.6 and f0<=0.2):
            hi=0 if f0>=f1 else 1
            div.append((j,hi))
    div=[(j,hi) for (j,hi) in div if not match(j,ST[ch]) and not match(j,REF[ch])]
    if not div: continue
    n68+=1
    sis=sisters_of(ch,s,e)
    if not sis: novel+=1; continue   # no sister to cross-map from
    sisAS={}
    for sreg in sis:
        for k,v in best_as_at(sreg).items():
            if k not in sisAS or v>sisAS[k]: sisAS[k]=v
    # reads carrying a divergent junction (in its high group)
    djreads=set()
    for (j,hi) in div:
        for rn,g in grp.items():
            if g==hi and j in rintr[rn]: djreads.add(rn)
    if not djreads: novel+=1; continue
    nx=sum(1 for rn in djreads if rn in sisAS and sisAS[rn] >= rAS.get(rn,-1e9))
    frac=nx/len(djreads)
    if frac>=0.5: artifact+=1; tag="ARTIFACT"
    else: novel+=1; tag="novel"
    if len(examples)<14: examples.append((ch,s,e,len(djreads),round(frac,2),tag))

print(f"non-corroborated divergent-junction loci examined (the 68): {n68}")
print(f"  CROSS-MAPPING ARTIFACT (>=50% of div-junction reads align >= as well to a sister): {artifact}")
print(f"  locus-specific / novel (div-junction reads best at this locus):                    {novel}")
print("\nexamples (chrom:start-end, #div-junction-reads, sister-AS-fraction, call):")
for x in examples: print(f"  {x[0]}:{x[1]}-{x[2]}  reads={x[3]} sister_frac={x[4]} {x[5]}")

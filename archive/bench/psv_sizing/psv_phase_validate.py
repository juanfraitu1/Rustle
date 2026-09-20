#!/usr/bin/env python3
"""Validate the primary-read PSV-phasing headroom: for each locus where primary-read
PSV-groups splice differently, record the divergent junction and check it is (a) MISSED by
StringTie (real headroom) and (b) CORROBORATED by RefSeq (a real isoform, not a cross-map
artifact). Defensible count = divergent junctions missed-by-StringTie AND in-RefSeq."""
import pysam, csv, collections, glob, re, os, sys
# Defaults are the paths this instrument was written against; each is overridable so it can be pointed
# at a real corpus. A missing one is a HARD ABORT below, never a silent fallback (ledger §6am).
BAM=os.environ.get("PSV_BAM","/home/juanfra/winloci_scratch/GGO.bam")
EX=os.environ.get("PSV_EXONS","/tmp/cre_guided/exons.tsv"); GW=os.environ.get("PSV_GW","/tmp/gw")
CHROMS={"NC_073235.2","NC_086018.1"}
MINDP=12; ALT_FRAC=0.25; MIN_HET=2; MIN_GRP=5; TOL=5

def die(msg):
    # §6am: this instrument's verdict is a RATE. A missing or empty input turns it into 0-of-0 printed
    # as a finding, which is indistinguishable from a real result. Abort loudly instead of scoring it.
    sys.stderr.write("PSV-PHASE ABORT: %s\n" % msg); sys.exit(2)

# §6am: prevents an absent BAM / exon table from yielding an empty locus list that still prints a verdict.
for _lbl, _p, _env in (("BAM", BAM, "PSV_BAM"), ("exon table", EX, "PSV_EXONS")):
    if not os.path.exists(_p): die("%s not found: %s (override with %s=...)" % (_lbl, _p, _env))
    if os.path.getsize(_p) == 0: die("%s is present but EMPTY: %s" % (_lbl, _p))

def gtf_introns(path):
    by=collections.defaultdict(list); tx=collections.defaultdict(list); st={}
    if not os.path.exists(path): return by
    for l in open(path):
        if l.startswith("#"):continue
        f=l.rstrip().split("\t")
        if len(f)<9 or f[2]!="exon":continue
        m=re.search(r'transcript_id "([^"]+)"',f[8]) or re.search(r'Name=([^;]+)',f[8])
        if not m:continue
        tx[m.group(1)].append((int(f[3]),int(f[4]))); st[m.group(1)]=f[0]
    for t,e in tx.items():
        e.sort()
        for i in range(len(e)-1): by[st[t]].append((e[i][1], e[i+1][0]-1))
    return by
def gff_introns(path):
    by=collections.defaultdict(list); tx=collections.defaultdict(list); st={}
    if not os.path.exists(path): return by
    for l in open(path):
        f=l.split("\t")
        if len(f)<9 or f[2]!="exon":continue
        m=re.search(r'Parent=([^;]+)',f[8])
        if not m:continue
        tx[m.group(1)].append((int(f[3]),int(f[4]))); st[m.group(1)]=f[0]
    for t,e in tx.items():
        e.sort()
        for i in range(len(e)-1): by[st[t]].append((e[i][1], e[i+1][0]-1))
    return by
def match(j, lst): return any(abs(j[0]-k[0])<=TOL and abs(j[1]-k[1])<=TOL for k in lst)

ST={}; REF={}
for ch in CHROMS:
    st_p=f"{GW}/st_{ch}.gtf"; ref_p=f"{GW}/ref_{ch}.gff3"
    # §6am false-PASS fix: the old `gtf_introns(...) if os.path.exists(...) else []` made BOTH comparators
    # empty when $GW was absent, so match() was False for everything and the novelty/agreement verdict
    # below was scored against a comparator of size ZERO. Missing comparator annotation now aborts.
    for _lbl, _p in (("StringTie GTF", st_p), ("RefSeq GFF3", ref_p)):
        if not os.path.exists(_p): die("%s for %s not found: %s (set the annotation dir with PSV_GW=...)" % (_lbl, ch, _p))
        if os.path.getsize(_p) == 0: die("%s for %s is present but EMPTY: %s" % (_lbl, ch, _p))
    ST[ch]=gtf_introns(st_p)[ch]; REF[ch]=gff_introns(ref_p)[ch]
    # §6am guards the EVIDENCE, not just the file: a parsable-but-intron-free comparator (or one holding
    # no records for this contig) scores EVERY junction as "missed by StringTie" / "not in RefSeq".
    if not ST[ch]: die("StringTie comparator for %s holds 0 introns (%s) — every junction would score as MISSED" % (ch, st_p))
    if not REF[ch]: die("RefSeq comparator for %s holds 0 introns (%s) — no junction could be CORROBORATED" % (ch, ref_p))

spans=[]
for r in csv.DictReader(open(EX),delimiter="\t"):
    if r["chrom"] not in CHROMS: continue
    s=[int(x) for x in r["exon_starts"].split(",") if x]; e=[int(x) for x in r["exon_ends"].split(",") if x]
    if s: spans.append((r["chrom"], s[0], e[-1]))
spans.sort(); loci=[]
for ch,s,e in spans:
    if loci and loci[-1][0]==ch and s<=loci[-1][2]: loci[-1]=(ch,loci[-1][1],max(loci[-1][2],e))
    else: loci.append((ch,s,e))
# §6am: prevents an exon table with no rows on these contigs from printing "0 ... 0 headroom" as a
# finding — that is 0-of-0, an empty computation, not a measured negative.
if not loci: die("no loci on %s in %s — nothing to score" % (sorted(CHROMS), EX))

bam=pysam.AlignmentFile(BAM,"rb")
n_struct=0; missed_st=0; missed_and_ref=0; missed_not_ref=0
ex_real=[]
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
    divj=[]
    for j in set(gi[0])|set(gi[1]):
        f0=gi[0][j]/max(1,gn[0]); f1=gi[1][j]/max(1,gn[1])
        if (f0>=0.6 and f1<=0.2) or (f1>=0.6 and f0<=0.2): divj.append(j)
    if not divj: continue
    n_struct+=1
    # is ANY divergent junction missed by StringTie? corroborated by RefSeq?
    miss=[j for j in divj if not match(j, ST[ch])]
    if miss:
        missed_st+=1
        real=[j for j in miss if match(j, REF[ch])]
        if real:
            missed_and_ref+=1
            if len(ex_real)<14: ex_real.append((ch,s,e,len(het),gn[0],gn[1],len(real)))
        else: missed_not_ref+=1

# §6am: the comparator size is part of the verdict — a rate scored against a small or absent annotation
# is not readable without it. Print it next to the numbers it conditions.
print(f"comparator sizes (introns loaded from {GW}; loci scanned: {len(loci)}):")
for ch in sorted(CHROMS): print(f"  {ch}: StringTie={len(ST[ch])}  RefSeq={len(REF[ch])}")
print(f"loci with PSV-linked structural divergence: {n_struct}")
print(f"  divergent junction MISSED by StringTie:                 {missed_st}")
print(f"    of those, CORROBORATED by RefSeq (real, not artifact): {missed_and_ref}  <-- DEFENSIBLE headroom")
print(f"    missed but NOT in RefSeq (novel / possible artifact):  {missed_not_ref}")
print(f"  (the rest: divergent junction already in StringTie = no headroom)")
print("\ndefensible examples (chrom:start-end, #het, grpA, grpB, #real-missed-junctions):")
for x in ex_real: print(f"  {x[0]}:{x[1]}-{x[2]}  het={x[3]} A={x[4]} B={x[5]} real_missed={x[6]}")

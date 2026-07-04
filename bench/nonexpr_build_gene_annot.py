import sys, re, bisect, pickle
import pysam

GFF="/home/juanfra/winloci_scratch/GGO_genomic.gff"
FAM_DNA="/mnt/c/Users/jfris/Desktop/Rustle/bench/genome_families_refined.tsv"
FAM_RNA="/mnt/c/Users/jfris/Desktop/Rustle/bench/family_rna_refine.tsv"
SEGDUP="/mnt/c/Users/jfris/Desktop/final.bed"
FASTA="/home/juanfra/winloci_scratch/GGO.fasta"

# ---------- 1. genes ----------
genes=[]  # dict per gene
name_re=re.compile(r'(?:^|;)Name=([^;]+)')
bt_re=re.compile(r'gene_biotype=([^;]+)')
with open(GFF) as fh:
    for line in fh:
        if line.startswith('#'): continue
        f=line.rstrip('\n').split('\t')
        if len(f)<9 or f[2]!='gene': continue
        chrom=f[0]; start=int(f[3]); end=int(f[4]); strand=f[6]
        m=name_re.search(f[8]); nm=m.group(1) if m else None
        b=bt_re.search(f[8]); bt=b.group(1) if b else 'NA'
        genes.append(dict(name=nm, chrom=chrom, start=start, end=end,
                          strand=strand, biotype=bt))
print("genes(feature=gene):", len(genes), file=sys.stderr)
names=[g['name'] for g in genes]
print("unique gene Names:", len(set(names)), "missing name:", sum(1 for n in names if n is None), file=sys.stderr)

# ---------- 2. family gene-name set ----------
famset=set()
with open(FAM_DNA) as fh:
    hdr=fh.readline().rstrip('\n').split('\t')
    gi=hdr.index('genes')
    for line in fh:
        f=line.rstrip('\n').split('\t')
        if len(f)<=gi: continue
        for gn in f[gi].split(','):
            gn=gn.strip()
            if gn: famset.add(gn)
n_dna=len(famset)
with open(FAM_RNA) as fh:
    hdr=fh.readline().rstrip('\n').split('\t')
    mi=hdr.index('member_gene'); di=hdr.index('dominant_gene')
    for line in fh:
        f=line.rstrip('\n').split('\t')
        if len(f)<=max(mi,di): continue
        if f[mi].strip(): famset.add(f[mi].strip())
        if f[di].strip(): famset.add(f[di].strip())
print("family gene names (DNA union RNA):", len(famset), "(DNA-only had", n_dna, ")", file=sys.stderr)

# ---------- 3. segdup merged intervals per chrom ----------
seg=dict()  # chrom -> list of (start,end) 0-based
with open(SEGDUP) as fh:
    for line in fh:
        f=line.rstrip('\n').split('\t')
        if len(f)<6: continue
        for (c,s,e) in ((f[0],f[1],f[2]),(f[3],f[4],f[5])):
            try:
                s=int(s); e=int(e)
            except ValueError:
                continue
            if e<=s: continue
            seg.setdefault(c,[]).append((s,e))
# merge
seg_merged=dict()
for c,ivs in seg.items():
    ivs.sort()
    m=[]
    cs,ce=ivs[0]
    for s,e in ivs[1:]:
        if s<=ce:
            if e>ce: ce=e
        else:
            m.append((cs,ce)); cs,ce=s,e
    m.append((cs,ce))
    seg_merged[c]=m
print("segdup chroms:", len(seg_merged), "e.g. chr1 intervals:", len(seg_merged.get('NC_073224.2',[])), file=sys.stderr)

def seg_cover_frac(chrom, s0, e0):
    # s0,e0 0-based half-open
    m=seg_merged.get(chrom)
    if not m: return 0.0
    starts=[x[0] for x in m]
    i=bisect.bisect_right(starts, e0)  # intervals with start< e0
    covered=0
    j=i-1
    # walk back while overlap possible
    while j>=0:
        s,e=m[j]
        if e<=s0:
            # could still be earlier overlapping (merged non-overlapping sorted, once e<=s0 and going back ends increase? not monotone). Safety: continue a bit
            # merged intervals are disjoint & sorted; if e<=s0, earlier ones have smaller e -> no overlap. break
            break
        ov=min(e,e0)-max(s,s0)
        if ov>0: covered+=ov
        j-=1
    L=e0-s0
    return covered/L if L>0 else 0.0

# ---------- 4. repeat (softmask lowercase) fraction ----------
fa=pysam.FastaFile(FASTA)
faset=set(fa.references)

for gi_,g in enumerate(genes):
    s0=g['start']-1; e0=g['end']
    g['seg_frac']=seg_cover_frac(g['chrom'], s0, e0)
    if g['chrom'] in faset:
        seq=fa.fetch(g['chrom'], s0, e0)
        low=sum(1 for ch in seq if ch=='a' or ch=='c' or ch=='g' or ch=='t' or ch=='n')
        # n is softmasked-ambiguous; count lowercase incl n. Uppercase N would be gap.
        lowmask=sum(1 for ch in seq if ch.islower())
        g['rep_frac']=lowmask/len(seq) if len(seq)>0 else 0.0
    else:
        g['rep_frac']=0.0
    g['in_famname']= (g['name'] in famset) if g['name'] else False
    g['multi']= bool(g['in_famname'] or g['seg_frac']>=0.5)
    if (gi_+1)%5000==0:
        print("annot", gi_+1, file=sys.stderr)

# save
with open("/home/juanfra/winloci_scratch/gene_annot.pkl","wb") as out:
    pickle.dump(genes, out)

# quick summary
import statistics
multi=sum(1 for g in genes if g['multi'])
byname=sum(1 for g in genes if g['in_famname'])
byseg=sum(1 for g in genes if g['seg_frac']>=0.5)
rep=sum(1 for g in genes if g['rep_frac']>=0.10)
print("=== ANNOT SUMMARY ===", file=sys.stderr)
print("total genes:", len(genes), file=sys.stderr)
print("multi-copy:", multi, " unique:", len(genes)-multi, file=sys.stderr)
print("  by family-name:", byname, " by segdup>=0.5:", byseg, file=sys.stderr)
print("repeat-overlap(rep_frac>=0.10):", rep, file=sys.stderr)

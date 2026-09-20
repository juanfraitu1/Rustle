"""Do any large-intron read-through chimeras bridge two PARALOG COPIES of one family
(a real fusion artifact in the copy structure), vs unrelated neighbor genes?
Targeted: fetch only over paralog-family gene regions; for each read with a >50k intron,
test whether the two flanking genes are a KNOWN paralog pair."""
import pysam, bisect, collections, json
BAM="/home/juanfra/winloci_scratch/GGO_mm.bam"
GENES="/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
EDGES="family_def_dna_pr_edges.tsv"   # cols: gene_a gene_b in_rna in_dna_loose ...
T=50000; N_OP=3

# paralog pairs (DNA-loose truth OR E_r) + the set of paralog genes
para=set(); pgenes=set()
import csv
rd=csv.DictReader(open(EDGES),delimiter='\t')
for r in rd:
    a,b=r['gene_a'],r['gene_b']
    is_par = r.get('in_dna_loose','0')=='1' or r.get('in_rna','0')=='1'
    if is_par:
        para.add(frozenset((a,b))); pgenes.add(a); pgenes.add(b)
print(f"paralog pairs: {len(para)}  paralog genes: {len(pgenes)}")

# gene spans (all genes, for flank lookup) + paralog gene regions to fetch
ivs=collections.defaultdict(list); starts=collections.defaultdict(list)
pregions=[]
for ln in open(GENES):
    p=ln.rstrip("\n").split("\t")
    if len(p)<4: continue
    c,s,e,n=p[0],int(p[1]),int(p[2]),p[3]
    ivs[c].append((s,e,n))
    if n in pgenes: pregions.append((c,s,e,n))
for c in ivs:
    ivs[c].sort(); starts[c]=[x[0] for x in ivs[c]]
def genes_at(c,pos):
    out=set(); L=ivs.get(c)
    if not L: return out
    i=bisect.bisect_right(starts[c],pos)
    for j in range(max(0,i-40),min(len(L),i+1)):
        s,e,n=L[j]
        if s<=pos<=e: out.add(n)
    return out

bam=pysam.AlignmentFile(BAM,"rb")
seen=set(); para_fusion=[]; unrelated=0; n_bigintron_reads=0
for (c,s,e,gname) in pregions:
    for a in bam.fetch(c,max(0,s-1000),e+1000):
        if a.is_secondary or a.is_supplementary or a.is_unmapped: continue
        key=(a.query_name,a.reference_start)
        if key in seen: continue
        ref=a.reference_start; big=[]
        for op,ln in a.cigartuples or []:
            if op in (0,2,7,8): ref+=ln
            elif op==N_OP:
                if ln>T: big.append((ref,ref+ln))
                ref+=ln
        if not big: continue
        seen.add(key); n_bigintron_reads+=1
        for (istart,iend) in big:
            gu=genes_at(a.reference_name,istart-1); gd=genes_at(a.reference_name,iend+1)
            hit=None
            for u in gu:
                for d in gd:
                    if u!=d and frozenset((u,d)) in para: hit=(u,d)
            if hit:
                para_fusion.append(dict(read=a.query_name,chrom=a.reference_name,
                    istart=istart,iend=iend,ilen=iend-istart,gene_up=hit[0],gene_dn=hit[1]))
            elif gu and gd and not (gu&gd):
                unrelated+=1
bam.close()

# summarise paralog fusions by gene-pair
byp=collections.Counter(frozenset((x['gene_up'],x['gene_dn'])) for x in para_fusion)
res=dict(n_bigintron_reads_in_para_regions=n_bigintron_reads,
         n_paralog_fusion_reads=len(para_fusion),
         n_unrelated_neighbor_reads=unrelated,
         paralog_fusion_pairs={"|".join(sorted(k)):v for k,v in byp.most_common()})
json.dump(res,open("g_chimera_family_xref.json","w"),indent=1)
print(f"\nreads w/ >50k intron in paralog-gene regions: {n_bigintron_reads}")
print(f"  -> PARALOG-COPY fusion (both flanks = known paralog pair): {len(para_fusion)} reads, {len(byp)} distinct pairs")
print(f"  -> unrelated-neighbor read-through: {unrelated} reads")
print("\ntop paralog-fusion gene-pairs (reads):")
for k,v in byp.most_common(15): print(f"  {v:4}  {' | '.join(sorted(k))}")

import pickle, sys, time
import pysam
import numpy as np

BAM="/home/juanfra/winloci_scratch/GGO_mm.bam"
genes=pickle.load(open("/home/juanfra/winloci_scratch/gene_annot.pkl","rb"))
ng=len(genes)
B=100000

# per-gene arrays
gstart=np.array([g['start']-1 for g in genes], dtype=np.int64)  # 0-based
gend  =np.array([g['end']     for g in genes], dtype=np.int64)
gchrom=[g['chrom'] for g in genes]

# bins per chrom
bins={}   # chrom -> dict(bin_id -> list of gidx)
for i in range(ng):
    c=gchrom[i]; s=gstart[i]; e=gend[i]
    d=bins.setdefault(c,{})
    b0=s//B; b1=(e-1)//B
    for bb in range(b0,b1+1):
        d.setdefault(bb,[]).append(i)

expr=np.zeros(ng,dtype=np.int64)  # primary mapq>0
p0  =np.zeros(ng,dtype=np.int64)  # primary mapq==0
sec =np.zeros(ng,dtype=np.int64)  # secondary
supp=np.zeros(ng,dtype=np.int64)  # supplementary

bam=pysam.AlignmentFile(BAM,'rb')
t0=time.time()
n=0
gs_l=gstart.tolist(); ge_l=gend.tolist()
for r in bam:
    if r.is_unmapped: continue
    n+=1
    if n % 1000000==0:
        print("reads", n, "%.0fs"%(time.time()-t0), file=sys.stderr)
    c=r.reference_name
    d=bins.get(c)
    if d is None: continue
    blocks=r.get_blocks()
    if not blocks:
        blocks=[(r.reference_start, r.reference_end)]
    # candidate gene idx
    cand=set()
    for bs,be in blocks:
        b0=bs//B; b1=(be-1)//B
        for bb in range(b0,b1+1):
            lst=d.get(bb)
            if lst: cand.update(lst)
    if not cand: continue
    # category
    if r.is_secondary:
        arr=sec
    elif r.is_supplementary:
        arr=supp
    else:
        arr= expr if r.mapping_quality>0 else p0
    for gi in cand:
        gs=gs_l[gi]; ge=ge_l[gi]
        # any block overlaps [gs,ge)
        for bs,be in blocks:
            if bs<ge and be>gs:
                arr[gi]+=1
                break
bam.close()
print("done reads", n, "%.0fs"%(time.time()-t0), file=sys.stderr)

np.savez("/home/juanfra/winloci_scratch/bam_counts.npz",
         expr=expr, p0=p0, sec=sec, supp=supp)
print("saved counts. totals: expr", int(expr.sum()), "p0", int(p0.sum()),
      "sec", int(sec.sum()), "supp", int(supp.sum()), file=sys.stderr)

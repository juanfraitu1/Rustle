import pickle, numpy as np
g=pickle.load(open('/home/juanfra/winloci_scratch/gene_annot.pkl','rb'))
d=np.load('/home/juanfra/winloci_scratch/bam_counts.npz')
expr=d['expr']; p0=d['p0']; sec=d['sec']; supp=d['supp']
mm=p0+sec+supp
multi=np.array([x['multi'] for x in g]); rep=np.array([x['rep_frac'] for x in g])
ng=len(g)
L=[]
out=[]
out.append("NON-EXPRESSED x COPY-STATUS x REPEAT : RNA multimap-spillover per gene (GGO_mm.bam, gorilla)")
out.append("Definitions: expression proxy = # primary(non-sec,non-supp) alignments MAPQ>0 covering gene (splice-aware block overlap).")
out.append("multimap/spillover = # alignments covering gene that are secondary(0x100) OR supplementary(0x800) OR primary-MAPQ0.")
out.append("copy-status multi = gene in genome_families_refined OR family_rna_refine (by Name) OR >=50% covered by SEDEF final.bed segdup; else unique.")
out.append("repeat_frac = lowercase(soft-masked) fraction of gene span in GGO.fasta; repeat-free = repeat_frac<0.10.")
out.append(f"BAM totals (alignment-gene records): expr(primMAPQ>0)={int(expr.sum())} primMAPQ0={int(p0.sum())} secondary={int(sec.sum())} supplementary={int(supp.sum())}")
out.append(f"N genes(feature=gene) = {ng}  | multi={int(multi.sum())} unique={int((~multi).sum())}")
out.append("")
def cell(mask):
    n=int(mask.sum()); s=mm[mask]
    return (n, float(np.median(s)), float(np.mean(s)),
            float((s>=1).mean()), float((s>=5).mean()), float((s>=50).mean()))
for T in [0,1,2]:
    ne=expr<=T; rf=rep<0.10
    out.append(f"=== NON-EXPRESSED = primary-MAPQ>0 reads <= {T}  (n={int(ne.sum())}) | repeat-free=repeat_frac<0.10 ===")
    out.append(f"{'cell':30s} {'n':>6} {'mm_med':>7} {'mm_mean':>9} {'f>=1':>6} {'f>=5':>6} {'f>=50':>6}")
    for cs,cm in [('UNIQUE',~multi),('MULTI',multi)]:
        for rl,rm in [('repeat-free',rf),('repeat-ovl',~rf)]:
            n,med,mean,f1,f5,f50=cell(ne&cm&rm)
            out.append(f"{cs+' x '+rl:30s} {n:6d} {med:7.1f} {mean:9.1f} {f1:6.3f} {f5:6.3f} {f50:6.3f}")
    out.append("")
out.append("HEADLINE (T<=1): non-expr UNIQUE repeat-free median multimap=0, only 9.7% have>=1 read (top1% genes=78% of mass, all uncatalogued paralogs eg ferritin-like).")
out.append("Dissociation held-repeat-free: f>=1  UNIQUE=0.097 vs MULTI=0.652 (paralog spillover, 6.7x).  Held-unique: repeat-free 0.097 vs repeat-ovl 0.286 (repeat spillover, 2.9x).")
txt="\n".join(out)
open('/home/juanfra/winloci_scratch/nonexpr_multimap_summary.txt','w').write(txt+"\n")
print(txt)

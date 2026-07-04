import pickle, sys
import numpy as np

genes=pickle.load(open("/home/juanfra/winloci_scratch/gene_annot.pkl","rb"))
d=np.load("/home/juanfra/winloci_scratch/bam_counts.npz")
expr=d['expr']; p0=d['p0']; sec=d['sec']; supp=d['supp']
ng=len(genes)
assert len(expr)==ng

length=np.array([g['end']-g['start']+1 for g in genes],dtype=np.float64)
multimap=p0+sec+supp
mm_per_kb=multimap/(length/1000.0)
multi=np.array([g['multi'] for g in genes])
repfrac=np.array([g['rep_frac'] for g in genes])
segfrac=np.array([g['seg_frac'] for g in genes])
byname=np.array([g['in_famname'] for g in genes])

# ---- write per-gene TSV ----
OUT="/home/juanfra/winloci_scratch/nonexpr_multimap_pergene.tsv"
with open(OUT,"w") as o:
    o.write("gene\tchrom\tstart\tend\tlength\tstrand\tbiotype\t"
            "copy_status\tin_family_name\tsegdup_frac\trepeat_frac\t"
            "expr_primary_mapq_gt0\tprim_mapq0\tsecondary\tsupplementary\t"
            "multimap\tmultimap_per_kb\n")
    for i,g in enumerate(genes):
        cs="multi" if multi[i] else "unique"
        o.write(f"{g['name']}\t{g['chrom']}\t{g['start']}\t{g['end']}\t"
                f"{int(length[i])}\t{g['strand']}\t{g['biotype']}\t{cs}\t"
                f"{int(byname[i])}\t{segfrac[i]:.3f}\t{repfrac[i]:.3f}\t"
                f"{int(expr[i])}\t{int(p0[i])}\t{int(sec[i])}\t{int(supp[i])}\t"
                f"{int(multimap[i])}\t{mm_per_kb[i]:.2f}\n")
print("wrote", OUT, file=sys.stderr)

def stats(mask):
    n=int(mask.sum())
    if n==0: return dict(n=0)
    mm=multimap[mask]
    pk=mm_per_kb[mask]
    return dict(n=n,
        med=float(np.median(mm)), mean=float(np.mean(mm)),
        med_pk=float(np.median(pk)), mean_pk=float(np.mean(pk)),
        f_ge1=float((mm>=1).mean()), f_ge5=float((mm>=5).mean()),
        f_ge10=float((mm>=10).mean()), f_ge50=float((mm>=50).mean()))

def show(label,mask):
    s=stats(mask)
    if s['n']==0:
        print(f"{label:52s} n=0"); return
    print(f"{label:52s} n={s['n']:6d} | mm med={s['med']:7.1f} mean={s['mean']:9.1f} "
          f"| /kb med={s['med_pk']:7.2f} mean={s['mean_pk']:8.2f} "
          f"| f>=1={s['f_ge1']:.3f} f>=5={s['f_ge5']:.3f} f>=10={s['f_ge10']:.3f} f>=50={s['f_ge50']:.3f}")

names=np.array([g['name'] for g in genes])
biot=np.array([g['biotype'] for g in genes])

for EXPR_T in [0,1,2]:
    nonexpr = expr<=EXPR_T
    print("\n"+"="*140)
    print(f"### NON-EXPRESSED = primary(MAPQ>0) reads <= {EXPR_T}   (n_nonexpr={int(nonexpr.sum())} of {ng})")
    for REP_T in [0.10]:
        repfree = repfrac<REP_T
        print(f"--- repeat-free = repeat_frac < {REP_T} ; repeat-overlap = repeat_frac >= {REP_T} ---")
        for cs,cmask in [("UNIQUE",~multi),("MULTI",multi)]:
            for rl,rmask in [("repeat-free",repfree),("repeat-ovl",~repfree)]:
                m = nonexpr & cmask & rmask
                show(f"  nonexpr x {cs:6s} x {rl}", m)
    # example genes for the decisive cells at EXPR_T for repeat-free 0.10
    if EXPR_T==1:
        repfree=repfrac<0.10
        print("  -- example genes (highest & typical) per decisive cell --")
        for cs,cmask,rl,rmask in [
            ("UNIQUE",~multi,"repeat-free",repfree),
            ("UNIQUE",~multi,"repeat-ovl",~repfree),
            ("MULTI", multi, "repeat-free",repfree),
            ("MULTI", multi, "repeat-ovl",~repfree)]:
            m=nonexpr&cmask&rmask
            idx=np.where(m)[0]
            if len(idx)==0:
                print(f"    {cs} {rl}: none"); continue
            order=idx[np.argsort(-multimap[idx])]
            top=order[:4]
            ex="; ".join(f"{names[j]}({biot[j]},mm={int(multimap[j])},rep={repfrac[j]:.2f},L={int(length[j])})" for j in top)
            print(f"    {cs} {rl} (n={len(idx)}): TOP {ex}")

# gradient of multimap vs repeat, within non-expressed UNIQUE
print("\n"+"="*140)
print("### GRADIENT: non-expressed(<=1) UNIQUE genes, multimap by repeat_frac bin")
nonexpr=expr<=1
uni=~multi
bins=[(0,0.05),(0.05,0.10),(0.10,0.20),(0.20,0.30),(0.30,0.50),(0.50,1.01)]
for lo,hi in bins:
    m=nonexpr&uni&(repfrac>=lo)&(repfrac<hi)
    show(f"  rep[{lo:.2f},{hi:.2f})", m)
print("### same gradient for non-expressed(<=1) MULTI genes")
for lo,hi in bins:
    m=nonexpr&multi&(repfrac>=lo)&(repfrac<hi)
    show(f"  rep[{lo:.2f},{hi:.2f})", m)

# overall class sizes 2x2x2 (all genes, expr_t=1)
print("\n"+"="*140)
print("### 2x2x2 CELL SIZES (all genes) expr_thr<=1, repfree<0.10")
nonexpr=expr<=1; repfree=repfrac<0.10
for et in [False,True]:
    emask = nonexpr if et else ~nonexpr
    lab= "NON-EXPR" if et else "EXPRESSED"
    for cs,cmask in [("unique",~multi),("multi",multi)]:
        for rl,rmask in [("rep-free",repfree),("rep-ovl",~repfree)]:
            print(f"  {lab:9s} {cs:6s} {rl:8s}: {int((emask&cmask&rmask).sum())}")

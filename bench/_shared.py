import os
from collections import Counter, defaultdict
from itertools import combinations
BENCH="/mnt/c/Users/jfris/Desktop/Rustle/bench"; SCRATCH="/home/juanfra/winloci_scratch"
ANNOT=os.path.join(SCRATCH,"annot_intervals.tsv"); GMETA=os.path.join(SCRATCH,"gene_meta_strand.tsv")
ALN=os.path.join(SCRATCH,"aln.m8"); RNA_CAT=os.path.join(BENCH,"family_rna_refine.tsv")
PRFAM=os.path.join(BENCH,"protein_families_refined.tsv"); TRUTH=os.path.join(BENCH,"family_def_dna_pr_edges.tsv")
COLIN=os.path.join(BENCH,"colinear_multiexon_gate.tsv"); ORACLE=os.path.join(BENCH,"diploid_cn_oracle.tsv")
ALPHA_P,C_P,MEGA_MIN=1e-5,0.50,500
NONCODING={"lncRNA","snRNA","snoRNA","tRNA","rRNA","miRNA","misc_RNA","ncRNA","V_segment","C_region"}
PSEUDO={"pseudogene","transcribed_pseudogene"}
biotype={}
with open(ANNOT) as fh:
    next(fh)
    for ln in fh:
        f=ln.rstrip("\n").split("\t"); g=f[4]
        if g not in biotype or f[3]=="protein_coding": biotype[g]=f[3]
gmeta={}
with open(GMETA) as fh:
    for ln in fh:
        c,s,e,st,n,b=ln.rstrip("\n").split("\t"); gmeta[n]=(c,int(s),int(e),st,b)
rna=defaultdict(set)
with open(RNA_CAT) as fh:
    h=fh.readline().rstrip("\n").split("\t"); fi,gi=h.index("family_id"),h.index("member_gene")
    for ln in fh:
        f=ln.rstrip("\n").split("\t"); g=f[gi]
        if g and not g.startswith("DN_"): rna[f[fi]].add(g)
rna=dict(rna)
g2f,sizes={},{}
with open(PRFAM) as fh:
    next(fh)
    for ln in fh:
        f=ln.rstrip("\n").split("\t"); fid=f[0]; genes=f[11].split(",")
        sizes[fid]=len(genes)
        for g in genes: g2f[g]=fid
mega={fid for fid,n in sizes.items() if n>=MEGA_MIN}
ep_all={}
with open(ALN) as fh:
    for ln in fh:
        f=ln.rstrip("\n").split("\t")
        if len(f)<6: continue
        q,t=f[0],f[1]
        if q==t: continue
        try: ev,pid,qc,tc=float(f[2]),float(f[3]),float(f[4]),float(f[5])
        except ValueError: continue
        ep_all[(q,t)]=(ev,qc,tc)
def ep_edge(a,b):
    for x,y in ((a,b),(b,a)):
        e=ep_all.get((x,y))
        if not e or e[0]>ALPHA_P: return False
    e1,e2=ep_all[(a,b)],ep_all[(b,a)]
    return min(e1[1],e1[2])>=C_P and min(e2[1],e2[2])>=C_P
def ep_hom(a,b):
    if ep_edge(a,b): return True
    fa,fb=g2f.get(a),g2f.get(b); return fa is not None and fa==fb and fa not in mega
dna_loose=set()
with open(TRUTH) as fh:
    next(fh)
    for ln in fh:
        f=ln.rstrip("\n").split("\t")
        if int(f[3])==1: dna_loose.add(frozenset((f[0],f[1])))
colin,gstm2f,magef={},{},{}
with open(COLIN) as fh:
    h=fh.readline().rstrip("\n").split("\t"); ai,bi=h.index("gene_a"),h.index("gene_b")
    coi,gi,mgi=h.index("colinear"),h.index("gstm2"),h.index("mage")
    for ln in fh:
        x=ln.rstrip("\n").split("\t"); k=frozenset((x[ai],x[bi]))
        try: colin[k]=int(x[coi])
        except ValueError: pass
        gstm2f[k]=(x[gi]=="1"); magef[k]=(x[mgi]=="1")
oracle_multi=set()
with open(ORACLE) as fh:
    h=fh.readline().rstrip("\n").split("\t")
    for ln in fh:
        f=ln.rstrip("\n").split("\t"); d=dict(zip(h,f))
        if d["class"]=="multi_copy": oracle_multi.add(d["gene"])
pair_fams=defaultdict(set)
for fid,gs in rna.items():
    if len(gs)<2: continue
    for a,b in combinations(sorted(gs),2): pair_fams[frozenset((a,b))].add(fid)
def classify_type(k):
    a,b=tuple(k); ba,bb=biotype.get(a),biotype.get(b); col=colin.get(k,0)
    if gstm2f.get(k) or magef.get(k): return "CARDINALITY_ARRAY"
    if col>=2: return "EP_OVERSPLIT_realdup"
    if ba in PSEUDO or bb in PSEUDO: return "RETROCOPY_PSEUDOGENE"
    if ba in NONCODING or bb in NONCODING: return "NONCODING_ELEMENT"
    if ba=="protein_coding" and bb=="protein_coding": return "SINGLE_DOMAIN_SHARE"
    return "UNANNOTATED_OTHER"
fp_pairs,hom_pairs=[],[]
for k in pair_fams:
    a,b=tuple(k)
    if ep_hom(a,b) or k in dna_loose: hom_pairs.append(k)
    else: fp_pairs.append(k)
fp_set=set(fp_pairs)
fp_type={k:classify_type(k) for k in fp_pairs}
TRUTH_ARTIFACT={"CARDINALITY_ARRAY","EP_OVERSPLIT_realdup"}
genuine_fp=[k for k in fp_pairs if fp_type[k] not in TRUTH_ARTIFACT]

#!/usr/bin/env python
"""Test candidate NEW FP-exclusion rules beyond the shipped 3 gates.
Operates on the opt-out (floor-off) RNA catalog family_rna_refine.tsv (dca64cbd, 573 fams),
post-shipped-gates (so anything caught here is NEW-beyond-shipped).
Independent truth = E_p (protein homology) + DNA-loose (cDNA id>=.90) + diploid CN oracle.
"""
import os
from collections import Counter, defaultdict
from itertools import combinations

BENCH = "/mnt/c/Users/jfris/Desktop/Rustle/bench"
SCRATCH = "/home/juanfra/winloci_scratch"
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
GMETA = os.path.join(SCRATCH, "gene_meta_strand.tsv")
ALN = os.path.join(SCRATCH, "aln.m8")
AVA = os.path.join(SCRATCH, "rep_ava.tsv")
RNA_CAT = os.path.join(BENCH, "family_rna_refine.tsv")
PRFAM = os.path.join(BENCH, "protein_families_refined.tsv")
TRUTH = os.path.join(BENCH, "family_def_dna_pr_edges.tsv")
FEAT = os.path.join(BENCH, "rna_only_edge_features.tsv")
COLIN = os.path.join(BENCH, "colinear_multiexon_gate.tsv")
ORACLE = os.path.join(BENCH, "diploid_cn_oracle.tsv")

ALPHA_P, C_P, MEGA_MIN = 1e-5, 0.50, 500
NONCODING = {"lncRNA","snRNA","snoRNA","tRNA","rRNA","miRNA","misc_RNA","ncRNA","V_segment","C_region"}
PSEUDO = {"pseudogene","transcribed_pseudogene"}

# ---------- loaders (mirror fp_by_type) ----------
biotype, gchrom = {}, {}
with open(ANNOT) as fh:
    next(fh)
    for ln in fh:
        f = ln.rstrip("\n").split("\t"); g = f[4]
        if g not in biotype or f[3] == "protein_coding":
            biotype[g] = f[3]; gchrom[g] = f[0]

# GFF gene meta: strand + coords (authoritative genome annotation)
gmeta = {}  # gene -> (chrom,start,end,strand,biotype)
with open(GMETA) as fh:
    for ln in fh:
        c,s,e,st,n,b = ln.rstrip("\n").split("\t")
        gmeta[n] = (c,int(s),int(e),st,b)

def load_rna():
    fams = defaultdict(set); fchrom = defaultdict(set); fmembers = defaultdict(list)
    with open(RNA_CAT) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        fi, gi, ci, si, ei = (h.index("family_id"), h.index("member_gene"), h.index("chrom"),
                              h.index("start"), h.index("end"))
        for ln in fh:
            f = ln.rstrip("\n").split("\t"); g = f[gi]
            if g and not g.startswith("DN_"):
                fams[f[fi]].add(g)
            fchrom[f[fi]].add(f[ci])
            fmembers[f[fi]].append((g, f[ci], int(f[si]), int(f[ei])))
    return dict(fams), fchrom, fmembers
rna, fchrom, fmembers = load_rna()

g2f, sizes = {}, {}
with open(PRFAM) as fh:
    next(fh)
    for ln in fh:
        f = ln.rstrip("\n").split("\t"); fid = f[0]; genes = f[11].split(",")
        sizes[fid] = len(genes)
        for g in genes: g2f[g] = fid
mega = {fid for fid, n in sizes.items() if n >= MEGA_MIN}

def load_ep():
    fwd = {}
    with open(ALN) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6: continue
            q, t = f[0], f[1]
            if q == t: continue
            try: ev, pid, qc, tc = float(f[2]), float(f[3]), float(f[4]), float(f[5])
            except ValueError: continue
            fwd[(q, t)] = (ev, qc, tc)
    return fwd
ep_all = load_ep()

def ep_edge(a, b):
    for x, y in ((a, b), (b, a)):
        e = ep_all.get((x, y))
        if not e or e[0] > ALPHA_P: return False
    e1, e2 = ep_all[(a, b)], ep_all[(b, a)]
    return min(e1[1], e1[2]) >= C_P and min(e2[1], e2[2]) >= C_P
def ep_hom(a, b):
    if ep_edge(a, b): return True
    fa, fb = g2f.get(a), g2f.get(b)
    return fa is not None and fa == fb and fa not in mega

dna_loose = set()
with open(TRUTH) as fh:
    next(fh)
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if int(f[3]) == 1: dna_loose.add(frozenset((f[0], f[1])))

# colinear per gene pair
colin, gstm2f, magef = {}, {}, {}
nmatched, nexon_min, exon_frac, isodom = {}, {}, {}, {}
with open(COLIN) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ai, bi = h.index("gene_a"), h.index("gene_b")
    coi, gi, mgi = h.index("colinear"), h.index("gstm2"), h.index("mage")
    nmi, nemi, efi, idi = (h.index("n_matched"), h.index("n_exon_min"),
                           h.index("exon_match_frac"), h.index("isolated_domain"))
    for ln in fh:
        x = ln.rstrip("\n").split("\t"); k = frozenset((x[ai], x[bi]))
        try: colin[k] = int(x[coi])
        except ValueError: pass
        gstm2f[k] = (x[gi] == "1"); magef[k] = (x[mgi] == "1")
        try:
            nmatched[k] = int(x[nmi]); nexon_min[k] = int(x[nemi])
            exon_frac[k] = float(x[efi]); isodom[k] = int(x[idi])
        except ValueError: pass

# AVA inverted-alignment strand per edge (cDNA sense-vs-antisense)
ID_MIN, COV_MIN = 0.90, 0.30
ava_strand = {}
with open(AVA) as fh:
    best = {}
    for ln in fh:
        q,t,ql,tl,m,bl,qs,qe,strand,mapq,de = ln.rstrip("\n").split("\t")
        if q==t: continue
        ql=int(ql); m=int(m); bl=int(bl); qs=int(qs); qe=int(qe)
        if bl==0 or ql==0: continue
        ident=m/bl; cov=(qe-qs)/ql
        if ident<ID_MIN or cov<COV_MIN: continue
        e=frozenset((q,t))
        if e not in best or ident>best[e][0]: best[e]=(ident,strand)
    ava_strand = {e:s for e,(i,s) in best.items()}

# diploid oracle: 57 multi-copy genes (gold DNA)
oracle_multi = set()
with open(ORACLE) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        d = dict(zip(h,f))
        if d["class"] == "multi_copy":
            oracle_multi.add(d["gene"])

# ---------- reconstruct pairs ----------
pair_fams = defaultdict(set)
for fid, gs in rna.items():
    if len(gs) < 2: continue
    for a,b in combinations(sorted(gs),2):
        pair_fams[frozenset((a,b))].add(fid)

def classify_type(k):
    a,b = tuple(k); ba,bb = biotype.get(a), biotype.get(b); col = colin.get(k,0)
    if gstm2f.get(k) or magef.get(k): return "CARDINALITY_ARRAY"
    if col >= 2: return "EP_OVERSPLIT_realdup"
    if ba in PSEUDO or bb in PSEUDO: return "RETROCOPY_PSEUDOGENE"
    if ba in NONCODING or bb in NONCODING: return "NONCODING_ELEMENT"
    if ba=="protein_coding" and bb=="protein_coding": return "SINGLE_DOMAIN_SHARE"
    return "UNANNOTATED_OTHER"

fp_pairs, hom_pairs = [], []
for k in pair_fams:
    a,b = tuple(k)
    if ep_hom(a,b) or k in dna_loose: hom_pairs.append(k)
    else: fp_pairs.append(k)
fp_type = {k: classify_type(k) for k in fp_pairs}
print(f"pairs={len(pair_fams)} HOM(real)={len(hom_pairs)} FP={len(fp_pairs)}")
TRUTH_ARTIFACT = {"CARDINALITY_ARRAY","EP_OVERSPLIT_realdup"}
genuine_fp = [k for k in fp_pairs if fp_type[k] not in TRUTH_ARTIFACT]
print(f"genuine FP (excl truth-artifact) = {len(genuine_fp)}")

# ---------- geometry helpers ----------
def geom(k):
    """return (same_chrom, opp_strand, recip_overlap_frac) using GFF gene coords, or None."""
    a,b = tuple(k)
    ma,mb = gmeta.get(a), gmeta.get(b)
    if not ma or not mb: return None
    if ma[0]!=mb[0]: return (False,None,0.0)
    opp = ma[3]!=mb[3]
    ov = min(ma[2],mb[2]) - max(ma[1],mb[1])
    la,lb = ma[2]-ma[1], mb[2]-mb[1]
    rof = ov/max(1,min(la,lb)) if ov>0 else 0.0
    return (True, opp, rof)

def report_rule(name, predicate, note=""):
    fp_hit = [k for k in fp_pairs if predicate(k)]
    real_hit = [k for k in hom_pairs if predicate(k)]
    by_type = Counter(fp_type[k] for k in fp_hit)
    genuine_hit = [k for k in fp_hit if fp_type[k] not in TRUTH_ARTIFACT]
    # oracle collateral: pairs among oracle multi-copy genes that are flagged (would cut a real dup edge)
    oracle_hit = [k for k in fp_hit+real_hit if all(g in oracle_multi for g in k)]
    print(f"\n=== RULE: {name} ===  {note}")
    print(f"  FP flagged: {len(fp_hit)} (genuine {len(genuine_hit)}) | REAL(E_p/DNA-loose) collateral: {len(real_hit)} | oracle-gene-pair collateral: {len(oracle_hit)}")
    print(f"  FP-by-type: {dict(by_type)}")
    if real_hit[:8]:
        print("  REAL collateral examples: " + " ; ".join(f"{tuple(k)[0]}({biotype.get(tuple(k)[0],'?')[:5]})+{tuple(k)[1]}({biotype.get(tuple(k)[1],'?')[:5]})" for k in real_hit[:8]))
    if genuine_hit[:6]:
        print("  genuine-FP caught examples: " + " ; ".join(f"{tuple(k)[0]}+{tuple(k)[1]}[{fp_type[k]}]" for k in genuine_hit[:6]))
    return fp_hit, real_hit, genuine_hit

# ---------- RULE 1: biotype-consistency (protein x noncoding) ----------
def r1_strict(k):
    a,b=tuple(k); ba,bb=biotype.get(a),biotype.get(b)
    return (ba=="protein_coding" and bb in NONCODING) or (bb=="protein_coding" and ba in NONCODING)
report_rule("1a biotype-mix protein x noncoding(lnc/sn/t/r...)", r1_strict,
            "targets NONCODING_ELEMENT")

def r1_incl_pseudo(k):
    a,b=tuple(k); ba,bb=biotype.get(a),biotype.get(b)
    pc = {"protein_coding"}
    def cls(x): return "pc" if x in pc else ("nc" if x in NONCODING else ("ps" if x in PSEUDO else "other"))
    ca,cb = cls(ba),cls(bb)
    return {ca,cb}=={"pc","nc"} or {ca,cb}=={"pc","ps"}
report_rule("1b biotype-mix protein x (noncoding OR pseudogene)", r1_incl_pseudo,
            "adds RETROCOPY_PSEUDOGENE (RISK: retrocopy is a real copy)")

# ---------- RULE 2: genomic antisense + reciprocal overlap ----------
for X in (0.10, 0.30, 0.50):
    def r2(k, X=X):
        g = geom(k)
        return bool(g) and g[0] and g[1] and g[2] >= X
    report_rule(f"2 same-contig OPP-strand recip-overlap>={X}", r2, "targets antisense/overlap")

# also: same-contig ANY-strand reciprocal overlap (the plain 'overlap' feature)
def r2_anystrand(k):
    g = geom(k); return bool(g) and g[0] and g[2] >= 0.30
report_rule("2' same-contig ANY-strand recip-overlap>=0.30", r2_anystrand, "(no strand filter, for contrast)")

# ---------- RULE 3: inverted cDNA alignment (AVA strand '-') ----------
def r3(k):
    return ava_strand.get(k) == "-"
report_rule("3 inverted cDNA alignment (AVA strand '-')", r3, "targets inverted/palindromic self-aln")

# ---------- RULE 4: terminal/single shared exon ----------
def r4_single(k):
    return nmatched.get(k, 99) == 1 and nexon_min.get(k, 0) >= 2
report_rule("4a single shared exon (n_matched==1, both multi-exon)", r4_single,
            "close to falsified colinear gate")
def r4_isodom(k):
    return isodom.get(k, 0) == 1
report_rule("4b isolated_domain flag==1", r4_isodom, "the shipped isolated-domain signal")

# ---------- 5th-candidate SCAN: enumerate genuine-FP regularities ----------
print("\n\n########## 5TH-CANDIDATE SCAN over genuine FP ##########")
# a) same-gene-name-family antisense breakdown
g_same_chrom = sum(1 for k in genuine_fp if (geom(k) or (0,))[0])
g_opp = sum(1 for k in genuine_fp if (geom(k) and geom(k)[0] and geom(k)[1]))
g_ov = sum(1 for k in genuine_fp if (geom(k) and geom(k)[0] and geom(k)[2]>0))
print(f"genuine FP: same_chrom={g_same_chrom}/{len(genuine_fp)} opp_strand(same-chrom)={g_opp} coord-overlap>0={g_ov}")
# b) how many genuine FP have NO cDNA AVA alignment at all (id<.90/cov<.30) => joined by non-cdna signal
no_ava = [k for k in genuine_fp if k not in ava_strand]
print(f"genuine FP with NO qualifying cDNA AVA aln (id>=.90&cov>=.30): {len(no_ava)}/{len(genuine_fp)}")
# c) exon-count asymmetry: single-exon locus (n_exon_min<=1) joined to multi-exon
mono = [k for k in genuine_fp if nexon_min.get(k, 99) <= 1]
print(f"genuine FP touching a MONO-exon locus (n_exon_min<=1): {len(mono)}")
# d) distribution of n_matched among genuine FP
nm_dist = Counter(nmatched.get(k, "NA") for k in genuine_fp)
print(f"genuine FP n_matched exon distribution: {dict(sorted((str(x),c) for x,c in nm_dist.items()))}")
# e) coverage asymmetry / length ratio via gene span
def lenratio(k):
    a,b=tuple(k); ma,mb=gmeta.get(a),gmeta.get(b)
    if not ma or not mb: return None
    la,lb=ma[2]-ma[1],mb[2]-mb[1]
    return min(la,lb)/max(la,lb)
lr = [lenratio(k) for k in genuine_fp if lenratio(k) is not None]
import statistics
print(f"genuine FP gene-span length-ratio median={statistics.median(lr):.2f} (real HOM median={statistics.median([lenratio(k) for k in hom_pairs if lenratio(k) is not None]):.2f})")
# f) tiny-gene noncoding: very short loci (< 500 bp) = tRNA/snRNA scale
tiny = [k for k in genuine_fp if any((gmeta.get(g) and gmeta[g][2]-gmeta[g][1] < 500) for g in k)]
print(f"genuine FP touching a <500bp locus: {len(tiny)}")

# overlap of rule outputs with each other (are rules redundant?)
print("\n########## RULE OVERLAP among genuine-FP ##########")
sets = {
 "R1a": set(k for k in genuine_fp if r1_strict(k)),
 "R2opp10": set(k for k in genuine_fp if (geom(k) and geom(k)[0] and geom(k)[1] and geom(k)[2]>=0.10)),
 "R3inv": set(k for k in genuine_fp if r3(k)),
 "R4single": set(k for k in genuine_fp if r4_single(k)),
}
union = set().union(*sets.values())
print(f"genuine FP total={len(genuine_fp)}; union caught by R1a|R2|R3|R4 = {len(union)}")
for n,s in sets.items():
    print(f"  {n}: {len(s)}")

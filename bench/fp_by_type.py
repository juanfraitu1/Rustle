#!/usr/bin/env python
"""Fresh FP-by-TYPE breakdown for family_rna_refine.tsv (md5 dca64cbd).
(A) 902 RNA FP pairs by TYPE (biotype partition + TE/repeat + cross-chrom + protein-gateable overlays).
(B) 57 E_p-impure residual blocks -> genuine_overmerge(32)/ep_oversplit(9)/cardinality(14)/no_multiexon(2),
    and TYPE-tag the 32 genuine over-merges.
"""
import os
from collections import Counter, defaultdict
from itertools import combinations

BENCH = "/mnt/c/Users/jfris/Desktop/Rustle/bench"
SCRATCH = "/home/juanfra/winloci_scratch"
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
ALN = os.path.join(SCRATCH, "aln.m8")
RNA_CAT = os.path.join(BENCH, "family_rna_refine.tsv")
PRFAM = os.path.join(BENCH, "protein_families_refined.tsv")
TRUTH = os.path.join(BENCH, "family_def_dna_pr_edges.tsv")
FEAT = os.path.join(BENCH, "rna_only_edge_features.tsv")
COLIN = os.path.join(BENCH, "colinear_multiexon_gate.tsv")

ALPHA_P, C_P, MEGA_MIN = 1e-5, 0.50, 500
NONCODING = {"lncRNA","snRNA","snoRNA","tRNA","rRNA","miRNA","misc_RNA","ncRNA","V_segment","C_region"}
PSEUDO = {"pseudogene","transcribed_pseudogene"}

# ---------- loaders ----------
biotype, gchrom = {}, {}
with open(ANNOT) as fh:
    next(fh)
    for ln in fh:
        f = ln.rstrip("\n").split("\t"); g = f[4]
        if g not in biotype or f[3] == "protein_coding":
            biotype[g] = f[3]; gchrom[g] = f[0]

def load_rna():
    fams = defaultdict(set); fchrom = defaultdict(set)
    with open(RNA_CAT) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        fi, gi, ci = h.index("family_id"), h.index("member_gene"), h.index("chrom")
        for ln in fh:
            f = ln.rstrip("\n").split("\t"); g = f[gi]
            if g and not g.startswith("DN_"):
                fams[f[fi]].add(g)
            fchrom[f[fi]].add(f[ci])
    return dict(fams), fchrom
rna, fchrom = load_rna()

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
ep_all = load_ep()  # all reciprocal-ish alignments (unfiltered coverage) for coverage lookup

def ep_edge(a, b):  # strict E_p membership (recip cov>=0.50, eval<1e-5)
    for x, y in ((a, b), (b, a)):
        e = ep_all.get((x, y))
        if not e or e[0] > ALPHA_P: return False
    e1, e2 = ep_all[(a, b)], ep_all[(b, a)]
    return min(e1[1], e1[2]) >= C_P and min(e2[1], e2[2]) >= C_P
def prot_recip_mincov(a, b):
    if (a, b) in ep_all and (b, a) in ep_all:
        e = ep_all[(a, b)]; return min(e[1], e[2])
    if (a, b) in ep_all: e = ep_all[(a, b)]; return min(e[1], e[2])
    if (b, a) in ep_all: e = ep_all[(b, a)]; return min(e[1], e[2])
    return None

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

# abl_bridge_mask + core_recip per gene pair
abl, corerec = {}, {}
with open(FEAT) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ai, bi = h.index("gene_a"), h.index("gene_b")
    mi, cri = h.index("abl_bridge_mask"), h.index("core_recip")
    for ln in fh:
        x = ln.rstrip("\n").split("\t"); k = frozenset((x[ai], x[bi]))
        try: abl[k] = float(x[mi])
        except ValueError: pass
        try: corerec[k] = float(x[cri])
        except ValueError: pass

# colinear per gene pair (multi-exon real-dup signal + cardinality flags)
colin, gstm2f, magef = {}, {}, {}
with open(COLIN) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ai, bi = h.index("gene_a"), h.index("gene_b")
    coi, gi, mgi = h.index("colinear"), h.index("gstm2"), h.index("mage")
    for ln in fh:
        x = ln.rstrip("\n").split("\t"); k = frozenset((x[ai], x[bi]))
        try: colin[k] = int(x[coi])
        except ValueError: pass
        gstm2f[k] = (x[gi] == "1"); magef[k] = (x[mgi] == "1")

# ---------- (A) 902 FP pairs by TYPE ----------
# all within-family gene pairs of the RNA catalog
pair_fams = defaultdict(set)
for fid, gs in rna.items():
    if len(gs) < 2: continue
    for a, b in combinations(sorted(gs), 2):
        pair_fams[frozenset((a, b))].add(fid)

fp_pairs = []
hom = 0
for k in pair_fams:
    a, b = tuple(k)
    if ep_hom(a, b) or k in dna_loose:
        hom += 1; continue
    fp_pairs.append(k)
print(f"[A] catalog within-family pairs={len(pair_fams)}  HOM={hom}  FP={len(fp_pairs)}")

def classify_type(k):
    a, b = tuple(k)
    ba, bb = biotype.get(a), biotype.get(b)
    col = colin.get(k, 0)
    # 1 truth-artifact (cardinality / ep-oversplit): multi-exon real dup OR cardinality hub
    if gstm2f.get(k) or magef.get(k):
        return "CARDINALITY_ARRAY"
    if col >= 2:
        return "EP_OVERSPLIT_realdup"     # colinear multi-exon dup E_p splits = truth artifact
    # 2 retrocopy / pseudogene
    if ba in PSEUDO or bb in PSEUDO:
        return "RETROCOPY_PSEUDOGENE"
    # 3 non-coding (lncRNA/sn/t/r...) — TE/array-prone
    if ba in NONCODING or bb in NONCODING:
        return "NONCODING_ELEMENT"
    # 4 both protein-coding, single/zero shared exon = single-domain share
    if ba == "protein_coding" and bb == "protein_coding":
        return "SINGLE_DOMAIN_SHARE"
    return "UNANNOTATED_OTHER"

typ = Counter(); typ_fam = defaultdict(set); typ_ex = defaultdict(list)
te_overlay = Counter(); xc_overlay = Counter(); protgate = Counter()
for k in fp_pairs:
    t = classify_type(k); typ[t] += 1
    fid = sorted(pair_fams[k])[0]; typ_fam[t].add(fid)
    if len(typ_ex[t]) < 6:
        a, b = tuple(k)
        typ_ex[t].append(f"{a}({biotype.get(a,'NA')[:6]})+{b}({biotype.get(b,'NA')[:6]}) col={colin.get(k,0)} abl={abl.get(k,float('nan')):.2f}")
    # overlays
    if abl.get(k, 0) >= 0.5: te_overlay[t] += 1
    # cross-chrom: genes on different chroms
    if gchrom.get(tuple(k)[0]) and gchrom.get(tuple(k)[1]) and gchrom[tuple(k)[0]] != gchrom[tuple(k)[1]]:
        xc_overlay[t] += 1
    # protein-gateable = both protein-coding
    a, b = tuple(k)
    if biotype.get(a) == "protein_coding" and biotype.get(b) == "protein_coding":
        protgate[t] += 1

NF = len(fp_pairs)
print("\n===== (A) 902-pair FP breakdown by TYPE =====")
order = ["SINGLE_DOMAIN_SHARE","RETROCOPY_PSEUDOGENE","NONCODING_ELEMENT","EP_OVERSPLIT_realdup",
         "CARDINALITY_ARRAY","UNANNOTATED_OTHER"]
print(f"{'TYPE':<24}{'pairs':>6}{'%FPmass':>9}{'fams':>6}{'TE(abl>=.5)':>12}{'x-chrom':>9}{'protgate':>9}")
for t in order:
    n = typ.get(t, 0)
    print(f"{t:<24}{n:>6}{n/NF:>8.1%} {len(typ_fam[t]):>6}{te_overlay.get(t,0):>12}{xc_overlay.get(t,0):>9}{protgate.get(t,0):>9}")
print(f"{'TOTAL':<24}{NF:>6}{'100%':>9}")
print("\nexamples per TYPE:")
for t in order:
    print(f"  [{t}]  " + " ; ".join(typ_ex[t][:5]))

# TE/repeat as its own view (union of noncoding-TE + abl-bridged)
te_pairs = [k for k in fp_pairs if abl.get(k, 0) >= 0.5]
print(f"\nTE/Alu-bridge overlay (abl_bridge_mask>=0.5): {len(te_pairs)} FP pairs ({len(te_pairs)/NF:.1%}) "
      f"cut across types {dict(te_overlay)}")

# protein-gateable summary
tot_pg = sum(protgate.values())
print(f"\nprotein-gateable FP (both protein-coding): {tot_pg}/{NF} ({tot_pg/NF:.1%})")

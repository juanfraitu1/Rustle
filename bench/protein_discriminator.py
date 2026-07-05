#!/usr/bin/env python
"""Test whether PROTEIN reciprocal coverage separates single-domain-shares (genuine E_r over-merge)
from divergent real paralogs (truthbar) where the NUCLEOTIDE coverage axis (core_recip) cannot.

Substrate: rna_only_edge_features.tsv (E_r-admitted candidate edges, cls in {TP,truthbar,genuine}).
Protein axis: aln.m8 reciprocal mmseqs (min(qcov,tcov)).  Nucleotide axis: core_recip.
Independent truth: in_dna_loose (cDNA 90% bar) = TP.  Semi-circular: in_ep defines truthbar/genuine.
"""
import os, statistics as st
from collections import Counter, defaultdict

BENCH = "/mnt/c/Users/jfris/Desktop/Rustle/bench"
SCRATCH = "/home/juanfra/winloci_scratch"
ALN = os.path.join(SCRATCH, "aln.m8")
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
FEAT = os.path.join(BENCH, "rna_only_edge_features.tsv")

# ---- biotype ----
biotype = {}
with open(ANNOT) as fh:
    next(fh)
    for ln in fh:
        f = ln.rstrip("\n").split("\t"); g = f[4]
        if g not in biotype or f[3] == "protein_coding":
            biotype[g] = f[3]

# ---- edge features ----
edges = []  # (a,b,cls,in_dna,in_ep,core_recip)
genes_needed = set()
with open(FEAT) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ai, bi = h.index("gene_a"), h.index("gene_b")
    ci, di, ei = h.index("cls"), h.index("in_dna_loose"), h.index("in_ep")
    cri = h.index("core_recip")
    for ln in fh:
        x = ln.rstrip("\n").split("\t")
        a, b = x[ai], x[bi]
        try: cr = float(x[cri])
        except ValueError: cr = float("nan")
        edges.append((a, b, x[ci], x[di], x[ei], cr))
        genes_needed.add(a); genes_needed.add(b)

# ---- aln.m8 reciprocal protein coverage (only for needed genes) ----
fwd = {}  # (q,t) -> (pident,qcov,tcov)
with open(ALN) as fh:
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if len(f) < 6: continue
        q, t = f[0], f[1]
        if q == t or q not in genes_needed or t not in genes_needed: continue
        try: fwd[(q, t)] = (float(f[3]), float(f[4]), float(f[5]))
        except ValueError: pass

def prot_cov(a, b):
    """reciprocal min/max coverage; None if not a reciprocal protein alignment at all."""
    if (a, b) in fwd and (b, a) in fwd:
        _, qc1, tc1 = fwd[(a, b)]
        # min over both directions' min-coverage; use best-symmetric
        mn = min(qc1, tc1)
        return mn, max(qc1, tc1)
    # one-directional alignment still informative (partial)
    if (a, b) in fwd:
        _, qc, tc = fwd[(a, b)]; return min(qc, tc), max(qc, tc)
    if (b, a) in fwd:
        _, qc, tc = fwd[(b, a)]; return min(qc, tc), max(qc, tc)
    return None

def pc(g):
    return biotype.get(g) == "protein_coding"

# ============================================================ distributions
def summ(vals):
    vals = [v for v in vals if v == v]
    if not vals: return "n=0"
    return f"n={len(vals):>4} med={st.median(vals):.3f} mean={sum(vals)/len(vals):.3f}"

print("############### A. PROTEIN min-cov distribution by class (ALL pairs) ###############")
by_cls_prot = defaultdict(list); by_cls_prot_none = Counter(); by_cls_nuc = defaultdict(list)
for a, b, cls, dna, ep, cr in edges:
    pcov = prot_cov(a, b)
    if pcov is None:
        by_cls_prot_none[cls] += 1
    else:
        by_cls_prot[cls].append(pcov[0])
    by_cls_nuc[cls].append(cr)
for cls in ["TP", "truthbar", "genuine"]:
    total = sum(1 for e in edges if e[2] == cls)
    print(f"  {cls:<9} PROT-min-cov {summ(by_cls_prot[cls])}  | no-protein-edge={by_cls_prot_none[cls]}/{total}"
          f"  || NUCLEO core_recip {summ(by_cls_nuc[cls])}")

print("\n############### B. protein-coding-only subset (gate applicable) ###############")
# restrict to BOTH protein-coding
bc = [(a, b, cls, dna, ep, cr) for a, b, cls, dna, ep, cr in edges if pc(a) and pc(b)]
print(f"  both-protein-coding edges: {len(bc)} of {len(edges)}")
for cls in ["TP", "truthbar", "genuine"]:
    sub = [e for e in bc if e[2] == cls]
    pv = [prot_cov(a, b)[0] for a, b, *_ in sub if prot_cov(a, b) is not None]
    none = sum(1 for a, b, *_ in sub if prot_cov(a, b) is None)
    nv = [e[5] for e in sub]
    print(f"  {cls:<9} PROT-min-cov {summ(pv)} no-edge={none}/{len(sub)} || NUC {summ(nv)}")

print("\n############### C. THE ORDERING FLIP (domain-share vs divergent-real) ###############")
# genuine(protein-coding) = domain-share FP ; truthbar(protein-coding) = divergent real
gen_pc = [e for e in bc if e[2] == "genuine"]
tb_pc = [e for e in bc if e[2] == "truthbar"]
def cov_or_zero(a, b):
    p = prot_cov(a, b); return p[0] if p else 0.0
gen_prot = [cov_or_zero(a, b) for a, b, *_ in gen_pc]
tb_prot = [cov_or_zero(a, b) for a, b, *_ in tb_pc]
gen_nuc = [e[5] for e in gen_pc if e[5] == e[5]]
tb_nuc = [e[5] for e in tb_pc if e[5] == e[5]]
print(f"  NUCLEOTIDE core_recip:  domain-share(genuine) med={st.median(gen_nuc):.3f}  "
      f"divergent-real(truthbar) med={st.median(tb_nuc):.3f}  "
      f"-> {'INVERTED (FP>=real, cannot gate)' if st.median(gen_nuc)>=st.median(tb_nuc) else 'ok'}")
print(f"  PROTEIN min-cov(0 if no edge): domain-share med={st.median(gen_prot):.3f}  "
      f"divergent-real med={st.median(tb_prot):.3f}  "
      f"-> {'CORRECT (FP<real)' if st.median(gen_prot)<st.median(tb_prot) else 'inverted'}")

# AUC-style separation on protein axis
def auc(pos, neg):  # prob(random pos < random neg) ; pos=FP(genuine) neg=real(truthbar)
    import bisect
    s = sorted(neg); n = len(s); tot = 0
    for p in pos:
        lo = bisect.bisect_left(s, p); hi = bisect.bisect_right(s, p)
        tot += lo + (hi - lo) * 0.5
    return tot / (len(pos) * n)
print(f"  separation P(FP protein-cov < real protein-cov) = {auc(gen_prot, tb_prot):.3f}  (1.0=perfect)")
print(f"  separation P(FP nucleo    < real nucleo)        = {auc(gen_nuc, tb_nuc):.3f}  (0.5=none / <0.5 inverted)")

print("\n############### D. NON-CIRCULAR collateral: cDNA-confirmed reals a protein gate would CUT ###############")
# TP = in_dna_loose=1 (cDNA, independent of aln.m8). how many have NO/low protein cov?
tp = [e for e in edges if e[2] == "TP"]
tp_noprot = [(a, b) for a, b, *_ in tp if prot_cov(a, b) is None]
tp_lowprot = [(a, b) for a, b, *_ in tp if prot_cov(a, b) is not None and prot_cov(a, b)[0] < 0.50]
tp_pc = [(a, b) for a, b, cls, dna, ep, cr in tp if pc(a) and pc(b)]
tp_pc_noprot = [(a, b) for a, b in tp_pc if prot_cov(a, b) is None]
tp_pc_low = [(a, b) for a, b in tp_pc if prot_cov(a, b) is not None and prot_cov(a, b)[0] < 0.50]
print(f"  TP (cDNA-confirmed real) edges: {len(tp)}")
print(f"    no protein edge at all      : {len(tp_noprot)}  ({len(tp_noprot)/len(tp):.1%})")
print(f"    protein min-cov < 0.50      : {len(tp_lowprot)}")
print(f"    -> would be CUT by prot-gate: {len(tp_noprot)+len(tp_lowprot)}  ({(len(tp_noprot)+len(tp_lowprot))/len(tp):.1%})")
print(f"  of those, BOTH-protein-coding (true collateral, not out-of-scope):")
print(f"    protein-coding TP           : {len(tp_pc)}")
print(f"    ...no-edge {len(tp_pc_noprot)} + low {len(tp_pc_low)} = CUT {len(tp_pc_noprot)+len(tp_pc_low)} "
      f"({(len(tp_pc_noprot)+len(tp_pc_low))/max(1,len(tp_pc)):.1%} of protein-coding cDNA-reals)")
# biotype of the non-coding TP collateral
nc_bt = Counter()
for a, b in tp_noprot:
    if not (pc(a) and pc(b)):
        nc_bt[tuple(sorted((biotype.get(a,'NA'), biotype.get(b,'NA'))))] += 1
print(f"  biotypes of non-protein-coding TP-no-protein-edge (out of gate scope): {dict(nc_bt.most_common(8))}")
print(f"  examples protein-coding cDNA-real CUT by prot-gate: {(tp_pc_noprot+tp_pc_low)[:10]}")

print("\n############### E. threshold sweep on protein min-cov (FP excluded vs real cut) ###############")
# population = both-protein-coding edges. FP=genuine, real=TP+truthbar
real_bc = [(a, b) for a, b, cls, *_ in bc if cls in ("TP", "truthbar")]
fp_bc = [(a, b) for a, b, cls, *_ in bc if cls == "genuine"]
print(f"  both-protein-coding: real(TP+truthbar)={len(real_bc)}  FP(genuine)={len(fp_bc)}")
for thr in [0.0001, 0.10, 0.25, 0.40, 0.50, 0.60, 0.75]:
    fp_excl = sum(1 for a, b in fp_bc if (prot_cov(a, b) is None) or prot_cov(a, b)[0] < thr)
    real_cut = sum(1 for a, b in real_bc if (prot_cov(a, b) is None) or prot_cov(a, b)[0] < thr)
    print(f"    gate min-cov>={thr:<6}: FP excluded {fp_excl:>4}/{len(fp_bc)} ({fp_excl/len(fp_bc):5.1%})  "
          f"real cut {real_cut:>4}/{len(real_bc)} ({real_cut/len(real_bc):5.1%})")

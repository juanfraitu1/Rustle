#!/usr/bin/env python
"""MECHANISM of the DNA-vs-RNA specificity gap.

PART 1  DNA-admits: classify every within-family gene-pair of genome_families_refined into
        HOMOLOGOUS (Ep or cDNA truth) vs FP-CLASS (genuine over-merge). Priority partition of FP:
          (1) repeat-bridge / mega-blob chaining   (family repeat_mosaic=1 OR spans >=5 non-mega Ep fams)
          (2) pseudogene / retrocopy link          (>=1 gene biotype pseudogene/transcribed_pseudogene)
          (3) non-coding-element segdup            (>=1 gene non-protein-coding: lncRNA/sn/sno/t/r/miRNA...)
          (4) single-domain share (protein-coding) (both protein_coding, distinct genes, no whole-mol homology)
          (5) unannotated / other
        Same partition applied to the RNA catalog FP pairs -> which classes are DNA-EXCLUSIVE vs shared.

PART 2  RNA-misses: DNA multi-copy families with NO co-grouping RNA family. Decompose missed genes into
          SILENT (not in expressed universe) vs DIVERGENT-BELOW-FLOOR (expressed but E_r floor rejected).
        Supporting direct count from rna_only_edge_features.tsv: in_ep=1 pairs E_r rejected.

Truths (non-circular, same as head-to-head):  Ep = protein_families_refined co-membership (non-mega) +
aln.m8 reciprocal edges;  cDNA = family_def_dna_pr_edges.tsv in_dna_loose.
Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/family_fp_mechanism.py
"""
import json, os, re
from collections import defaultdict, Counter
from itertools import combinations

BENCH = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"
ALN = os.path.join(SCRATCH, "aln.m8")
ANNOT = os.path.join(SCRATCH, "annot_intervals.tsv")
DNA_CAT = os.path.join(BENCH, "genome_families_refined.tsv")
RNA_CAT = os.path.join(BENCH, "family_rna_refine.tsv")
TRUTH = os.path.join(BENCH, "family_def_dna_pr_edges.tsv")
PRFAM = os.path.join(BENCH, "protein_families_refined.tsv")
RNAFEAT = os.path.join(BENCH, "rna_only_edge_features.tsv")

ALPHA_P, C_P, MEGA_MIN = 1e-5, 0.50, 500
NONCODING = {"lncRNA","snRNA","snoRNA","tRNA","rRNA","miRNA","misc_RNA","ncRNA","V_segment","C_region"}
PSEUDO = {"pseudogene","transcribed_pseudogene"}
CHAIN_SPAN = 5   # a family spanning >=5 distinct non-mega Ep families = transitive chaining


# ------------------------------------------------------------------ loaders
def load_dna_catalog():
    fams, mosaic = {}, {}
    with open(DNA_CAT) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        gi, fi, mi = hdr.index("genes"), hdr.index("family_id"), hdr.index("repeat_mosaic")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            fams[f[fi]] = {g for g in f[gi].split(",") if g}
            mosaic[f[fi]] = (f[mi] == "1")
    return fams, mosaic

def load_rna_catalog():
    fams = defaultdict(set)
    with open(RNA_CAT) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        fi, mgi = hdr.index("family_id"), hdr.index("member_gene")
        for ln in fh:
            f = ln.rstrip("\n").split("\t"); g = f[mgi]
            if g and not g.startswith("DN_"):
                fams[f[fi]].add(g)
    return dict(fams)

def load_prfam():
    g2f, sizes = {}, {}
    with open(PRFAM) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t"); fid = f[0]; genes = f[11].split(",")
            sizes[fid] = len(genes)
            for g in genes: g2f[g] = fid
    mega = {fid for fid, n in sizes.items() if n >= MEGA_MIN}
    return g2f, mega

def load_ep_direct():
    fwd = {}
    with open(ALN) as fh:
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if len(f) < 6: continue
            q, t = f[0], f[1]
            if q == t: continue
            try: ev, pid, qc, tc = float(f[2]), float(f[3]), float(f[4]), float(f[5])
            except ValueError: continue
            if ev > ALPHA_P or qc < C_P or tc < C_P: continue
            fwd[(q, t)] = pid
    return {frozenset((q, t)) for (q, t) in fwd if (t, q) in fwd}

def load_truth():
    dna_loose = set()
    with open(TRUTH) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if int(f[3]) == 1: dna_loose.add(frozenset((f[0], f[1])))
    return dna_loose

def load_biotype():
    bt = {}
    with open(ANNOT) as fh:
        next(fh)
        for ln in fh:
            f = ln.rstrip("\n").split("\t"); g = f[4]
            # prefer protein_coding if any interval says so; else first seen
            if g not in bt or f[3] == "protein_coding":
                bt[g] = f[3]
    return bt

def load_expressed_universe():
    """rna_only_edge_features.tsv = E_r-ADMITTED edges only (in_er_refined==1). Its gene set = genes
    E_r placed into a refined family = 'present in RNA transcript output' proxy. cls semantics:
      TP       = E_r edge that is cDNA-confirmed (in_dna_loose)
      truthbar = E_r edge, in_ep=1 (real protein paralog) but below cDNA-90% bar = DIVERGENT paralog
                 E_r CORRECTLY grouped (evidence E_r's floor reaches divergent paralogs when expressed)
      genuine  = E_r edge, not homologous = E_r's own over-merge FP
    So this file measures E_r's ADMITTED edges, NOT its misses."""
    genes = set(); ep_divergent_grouped = []; ep_tp = 0
    with open(RNAFEAT) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ai, bi = hdr.index("gene_a"), hdr.index("gene_b")
        epi, cls = hdr.index("in_ep"), hdr.index("cls")
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            genes.add(f[ai]); genes.add(f[bi])
            if f[epi] == "1":
                if f[cls] == "TP": ep_tp += 1
                else: ep_divergent_grouped.append((f[ai], f[bi], f[cls]))
    return genes, ep_divergent_grouped, ep_tp


# ------------------------------------------------------------------ helpers
def ep_hom(key, ep_direct, g2f, mega):
    if key in ep_direct: return True
    a, b = tuple(key); fa, fb = g2f.get(a), g2f.get(b)
    return fa is not None and fa == fb and fa not in mega

def fam_ep_span(gs, g2f, mega):
    return len({g2f[g] for g in gs if g in g2f and g2f[g] not in mega})

def classify_fp(key, fam_is_blob, biotype):
    """priority partition of an FP (over-merge) pair."""
    if fam_is_blob:                     return "1_repeat_bridge_blob"
    a, b = tuple(key); ba, bb = biotype.get(a), biotype.get(b)
    if ba in PSEUDO or bb in PSEUDO:    return "2_pseudogene_retrocopy"
    if ba in NONCODING or bb in NONCODING: return "3_noncoding_element"
    if ba == "protein_coding" and bb == "protein_coding": return "4_domain_share_pc"
    return "5_unannotated_other"


def analyze_catalog(name, fams, mosaic, ep_direct, g2f, mega, dna_loose, biotype):
    # family-level: which families are "blob" (repeat_mosaic OR ep-span>=CHAIN_SPAN)
    fam_blob = {}
    for fid, gs in fams.items():
        fam_blob[fid] = (mosaic.get(fid, False) or fam_ep_span(gs, g2f, mega) >= CHAIN_SPAN)
    # all within-family pairs, tagged with family blob status
    pair_fam = {}
    for fid, gs in fams.items():
        if len(gs) < 2: continue
        for a, b in combinations(sorted(gs), 2):
            k = frozenset((a, b))
            # a pair may sit in multiple families (RNA); keep blob if any blob
            pair_fam[k] = pair_fam.get(k, False) or fam_blob[fid]
    total_pairs = len(pair_fam)
    hom = 0; fp_class = Counter(); fp_class_fams = defaultdict(set); fp_examples = defaultdict(list)
    # map pair->one family id for example reporting / family counting
    pair_one_fam = {}
    for fid, gs in fams.items():
        if len(gs) < 2: continue
        for a, b in combinations(sorted(gs), 2):
            pair_one_fam.setdefault(frozenset((a, b)), fid)
    for k, isblob in pair_fam.items():
        if ep_hom(k, ep_direct, g2f, mega) or k in dna_loose:
            hom += 1; continue
        c = classify_fp(k, isblob, biotype)
        fp_class[c] += 1
        fid = pair_one_fam[k]; fp_class_fams[c].add(fid)
        if len(fp_examples[c]) < 8:
            a, b = tuple(k)
            fp_examples[c].append((fid, a, biotype.get(a, "NA"), b, biotype.get(b, "NA")))
    fp_total = sum(fp_class.values())
    return dict(name=name, total_pairs=total_pairs, hom_pairs=hom, fp_pairs=fp_total,
                fp_class={c: fp_class[c] for c in sorted(fp_class)},
                fp_class_families={c: len(fp_class_fams[c]) for c in sorted(fp_class_fams)},
                fp_examples={c: fp_examples[c] for c in sorted(fp_examples)},
                n_blob_families=sum(1 for v in fam_blob.values() if v))


def main():
    dna, mosaic_d = load_dna_catalog()
    rna = load_rna_catalog()
    mosaic_r = {}  # RNA catalog has no repeat_mosaic flag
    g2f, mega = load_prfam()
    ep_direct = load_ep_direct()
    dna_loose = load_truth()
    biotype = load_biotype()
    expr_genes, ep_divergent_grouped, ep_tp = load_expressed_universe()
    # broaden expressed proxy with RNA-catalog member genes (grouped multi-copy set)
    for gs in rna.values(): expr_genes |= gs
    print(f"[load] DNA_fam={len(dna)} RNA_fam={len(rna)} Ep_genes={len(g2f)} mega={sorted(mega)} "
          f"ep_edges={len(ep_direct)} dna_loose={len(dna_loose)} biotype_genes={len(biotype)} "
          f"expressed_universe={len(expr_genes)}")

    # ============ PART 1 ============
    D = analyze_catalog("DNA", dna, mosaic_d, ep_direct, g2f, mega, dna_loose, biotype)
    R = analyze_catalog("RNA", rna, mosaic_r, ep_direct, g2f, mega, dna_loose, biotype)

    CLASS_LABEL = {"1_repeat_bridge_blob":"repeat-bridge / mega-blob chaining",
                   "2_pseudogene_retrocopy":"pseudogene / retrocopy link",
                   "3_noncoding_element":"non-coding-element segdup",
                   "4_domain_share_pc":"single-domain share (protein-coding)",
                   "5_unannotated_other":"unannotated / other"}
    def show_part1(res):
        print(f"\n== {res['name']}  total within-family pairs={res['total_pairs']}  "
              f"HOMOLOGOUS(Ep|cDNA)={res['hom_pairs']}  FP(over-merge)={res['fp_pairs']}  "
              f"(edge precision={res['hom_pairs']/max(1,res['total_pairs']):.3f})")
        for c in ["1_repeat_bridge_blob","2_pseudogene_retrocopy","3_noncoding_element",
                  "4_domain_share_pc","5_unannotated_other"]:
            n = res['fp_class'].get(c,0); nf = res['fp_class_families'].get(c,0)
            frac_fp = n/max(1,res['fp_pairs']); frac_all = n/max(1,res['total_pairs'])
            print(f"   {CLASS_LABEL[c]:<40} pairs={n:>6} ({frac_fp:5.1%} of FP, {frac_all:5.1%} of all)  families={nf}")
    print("\n############### PART 1  FP CLASSES DNA-admits (specificity leak) ###############")
    show_part1(D); show_part1(R)
    print("\n-- concrete DNA examples per class (fam, geneA[biotype], geneB[biotype]) --")
    for c in ["4_domain_share_pc","2_pseudogene_retrocopy","3_noncoding_element","1_repeat_bridge_blob"]:
        print(f"  [{CLASS_LABEL[c]}]")
        for ex in D['fp_examples'].get(c, [])[:5]:
            print(f"     {ex[0]:>8}  {ex[1]}({ex[2]}) -- {ex[3]}({ex[4]})")

    # ============ PART 2 ============
    print("\n############### PART 2  REAL COPIES RNA MISSES (recall/completeness) ###############")
    # RNA gene co-grouping: gene->set(rna family), and set of all RNA genes
    rna_genes = set()
    for gs in rna.values(): rna_genes |= gs
    # RNA within-family pairs (for co-group test)
    rna_pairs = set()
    for gs in rna.values():
        for a, b in combinations(sorted(gs), 2): rna_pairs.add(frozenset((a, b)))

    dna_mc = {fid: gs for fid, gs in dna.items() if len(gs) >= 2}
    miss_no_overlap = []      # DNA fam sharing 0 genes with any RNA gene
    miss_no_cogroup = []      # DNA fam where RNA groups <2 of its genes together
    recovered = []            # DNA fam with >=1 within-family pair also an RNA pair
    for fid, gs in dna_mc.items():
        shared = gs & rna_genes
        cogrouped = any(frozenset((a, b)) in rna_pairs for a, b in combinations(sorted(gs), 2))
        if cogrouped: recovered.append(fid)
        else:
            miss_no_cogroup.append(fid)
            if not shared: miss_no_overlap.append(fid)

    def ep_real_family(gs):   # >=2 genes in same non-mega Ep family = genuine paralog family
        c = Counter(g2f[g] for g in gs if g in g2f and g2f[g] not in mega)
        return any(v >= 2 for v in c.values())
    miss_real = [fid for fid in miss_no_cogroup if ep_real_family(dna_mc[fid])]
    miss_real_noov = [fid for fid in miss_no_overlap if ep_real_family(dna_mc[fid])]

    print(f"DNA multi-copy families                     : {len(dna_mc)}")
    print(f"  recovered by RNA (>=2 genes co-grouped)   : {len(recovered)}  ({len(recovered)/len(dna_mc):.1%})")
    print(f"  NOT co-grouped by RNA                     : {len(miss_no_cogroup)}  ({len(miss_no_cogroup)/len(dna_mc):.1%})")
    print(f"    of which share 0 genes with any RNA fam : {len(miss_no_overlap)}  (entirely invisible to RNA)")
    print(f"  MISSED families that are Ep-REAL paralogs : {len(miss_real)}  (genuine families DNA sees, RNA misses)")
    print(f"    Ep-real AND 0 RNA-gene-overlap          : {len(miss_real_noov)}")

    # decompose missed genes: silent vs divergent-below-floor
    missed_genes = set()
    for fid in miss_no_cogroup: missed_genes |= dna_mc[fid]
    silent = {g for g in missed_genes if g not in expr_genes}
    expressed_ungrouped = {g for g in missed_genes if g in expr_genes}
    print(f"\n  missed genes (in non-co-grouped DNA fams) : {len(missed_genes)}")
    print(f"    NOT in RNA output (silent/unassembled)  : {len(silent)}  ({len(silent)/max(1,len(missed_genes)):.1%})")
    print(f"    present in RNA output but DNA-fam not rec: {len(expressed_ungrouped)}  "
          f"({len(expressed_ungrouped)/max(1,len(missed_genes)):.1%})")

    print(f"\n  E_r-ADMITTED edge audit (rna_only_edge_features, in_er_refined==1):")
    print(f"    E_r edges that are cDNA-confirmed (TP)          : {ep_tp}")
    print(f"    E_r edges = real Ep-paralog below cDNA-bar      : {len(ep_divergent_grouped)}  "
          f"(E_r CORRECTLY grouped divergent paralogs; NOT misses)")
    print(f"      examples: {[(a,b) for a,b,_ in ep_divergent_grouped[:8]]}")

    # examples of Ep-real missed families
    print(f"\n  examples of Ep-REAL DNA families RNA misses entirely (fid: genes):")
    shown = 0
    for fid in miss_real:
        gs = dna_mc[fid]
        if 2 <= len(gs) <= 8:
            print(f"     {fid}: {sorted(gs)}  silent={sum(1 for g in gs if g not in expr_genes and g not in rna_genes)}/{len(gs)}")
            shown += 1
        if shown >= 10: break

    # ============ PART 3 ============
    print("\n############### PART 3  WHICH DIRECTION DOMINATES ###############")
    dna_fp = D['fp_pairs']; rna_fp = R['fp_pairs']
    print(f"SPECIFICITY cost (DNA extra over-merge): DNA FP pairs={dna_fp} vs RNA FP pairs={rna_fp}")
    print(f"  DNA FP class shares: " + ", ".join(
        f"{CLASS_LABEL[c].split(' /')[0].split(' (')[0]}={D['fp_class'].get(c,0)/max(1,dna_fp):.0%}"
        for c in ["1_repeat_bridge_blob","4_domain_share_pc","3_noncoding_element","2_pseudogene_retrocopy","5_unannotated_other"]))
    print(f"RECALL cost (families RNA misses): {len(miss_no_cogroup)}/{len(dna_mc)} DNA families "
          f"({len(miss_real)} Ep-real); missed genes {len(silent)} silent vs {len(expressed_ungrouped)} expressed-ungrouped")

    out = dict(part1_dna=D, part1_rna=R,
               part2=dict(dna_mc=len(dna_mc), recovered=len(recovered),
                          miss_no_cogroup=len(miss_no_cogroup), miss_no_overlap=len(miss_no_overlap),
                          miss_real=len(miss_real), miss_real_noov=len(miss_real_noov),
                          missed_genes=len(missed_genes), silent=len(silent),
                          expressed_ungrouped=len(expressed_ungrouped),
                          ep_tp=ep_tp, ep_divergent_grouped=len(ep_divergent_grouped)))
    # drop example tuples that aren't json-friendly-nested (keep as lists)
    with open(os.path.join(BENCH, "family_fp_mechanism.json"), "w") as fh:
        json.dump(out, fh, indent=2, default=list)
    print("\nwrote bench/family_fp_mechanism.json")


if __name__ == "__main__":
    main()

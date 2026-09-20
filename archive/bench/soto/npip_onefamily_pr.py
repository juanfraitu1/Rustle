#!/usr/bin/env python3
"""NPIP-as-ONE-family pairwise co-membership precision/recall across catalog arms.

⚠ DENOMINATORS ARE STATED, NOT IMPLIED (the project has retracted 7 numbers to a denominator
conditioned on the prediction). A predicted co-member pair is scored only when BOTH endpoints carry a
Soto gene label; pairs with an UNLABELLED endpoint are reported in their own bucket and are neither
numerator nor denominator. Recall's denominator is the truth pair set, which does not depend on any
prediction.

Two independent locus->gene assignment rules are reported. If they disagree the number is fragile and
must not be quoted as a single figure.
  A) 'genefrac' : best fraction-of-GENE covered (the rule ledger section 6co used)
  B) 'recip50'  : reciprocal overlap >= 0.50 of BOTH gene and copy (guards the large-gene attractor
                  that produced the retracted 'titin' instrument in section 6br)
"""
import csv, sys, itertools, collections

BENCH = "/mnt/c/Users/jfris/Desktop/Rustle/bench/soto"
ARMS = "/mnt/linuxdisk/home/juanfraitu/soto_adj"

def load_truth(npip_bed):
    genes = []   # (chrom, start, end, gene, family)
    for line in open(f"{BENCH}/{npip_bed}"):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        genes.append((f[0], int(f[1]), int(f[2]), f[3].split("|")[0], "ID_154"))
    for line in open(f"{BENCH}/80_fams.chr.bed"):
        f = line.rstrip("\n").split("\t")
        g, fid = f[3].split("|")
        if fid != "ID_154":
            genes.append((f[0], int(f[1]), int(f[2]), g, fid))
    return genes

def assign(chrom, start, end, genes, rule):
    best = None
    for gc, gs, ge, g, fid in genes:
        if gc != chrom: continue
        ov = min(ge, end) - max(gs, start)
        if ov <= 0: continue
        fg = ov / max(1, ge - gs)
        fc = ov / max(1, end - start)
        if rule == "recip50":
            if fg < 0.50 or fc < 0.50: continue
            score = min(fg, fc)
        else:
            score = fg
        if best is None or score > best[0]: best = (score, g, fid)
    return (best[1], best[2]) if best else (None, None)

def score(arm, genes, rule, npip_bed):
    path = f"{ARMS}/{arm}/cat.copies.tsv"
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    fam_of, lab = {}, {}
    byfam = collections.defaultdict(list)
    for i, r in enumerate(rows):
        ch = r.get("chrom") or r.get("chr")
        s, e = int(r["start"]), int(r["end"])
        g, fid = assign(ch, s, e, genes, rule)
        lab[i] = (g, fid)
        byfam[r["family_id"]].append(i)
    TP = FP = UNK = 0
    for fid_pred, members in byfam.items():
        for a, b in itertools.combinations(members, 2):
            (ga, fa), (gb, fb) = lab[a], lab[b]
            if ga is None or gb is None:
                if fa == "ID_154" or fb == "ID_154": UNK += 1
                continue
            if ga == gb: continue                      # same gene: not a cross-locus claim
            if fa == "ID_154" and fb == "ID_154": TP += 1
            elif fa == "ID_154" or fb == "ID_154":     FP += 1
    # recall denominator: NPIP truth genes DETECTED at all, paired
    npip_genes = {g for _, _, _, g, f in genes if f == "ID_154"}
    found = {lab[i][0] for i in lab if lab[i][1] == "ID_154" and lab[i][0]}
    truth_pairs = len(npip_genes) * (len(npip_genes) - 1) // 2
    # recovered truth gene-pairs actually co-predicted
    rec = set()
    for fid_pred, members in byfam.items():
        gs = sorted({lab[i][0] for i in members if lab[i][1] == "ID_154" and lab[i][0]})
        for a, b in itertools.combinations(gs, 2): rec.add((a, b))
    P = TP / (TP + FP) if (TP + FP) else float("nan")
    R = len(rec) / truth_pairs if truth_pairs else float("nan")
    return dict(arm=arm, rule=rule, TP=TP, FP=FP, UNK=UNK, P=P,
                rec_pairs=len(rec), truth_pairs=truth_pairs, R=R,
                genes_found=len(found), genes_total=len(npip_genes), copies=len(rows))

if __name__ == "__main__":
    npip_bed = sys.argv[1] if len(sys.argv) > 1 else "npip_ID154.adjudicated.bed"
    genes = load_truth(npip_bed)
    print(f"TRUTH = {npip_bed}; NPIP treated as ONE family (ID_154)\n")
    for rule in ("genefrac", "recip50"):
        print(f"--- assignment rule: {rule} ---")
        print(f"{'arm':<14}{'copies':>7}{'TP':>6}{'FP':>5}{'unlab':>7}{'precision':>11}{'recall':>9}   (P denom = TP+FP; R denom = truth gene-pairs)")
        for arm in ("arm_off", "arm_cl30", "arm_nostub", "arm_span", "arm_spancl30", "arm_repmask"):
            try: d = score(arm, genes, rule, npip_bed)
            except Exception as ex:
                print(f"{arm:<14}  ERROR {ex}"); continue
            print(f"{d['arm']:<14}{d['copies']:>7}{d['TP']:>6}{d['FP']:>5}{d['UNK']:>7}"
                  f"{d['P']:>11.4f}{d['R']:>9.4f}   [{d['TP']}/{d['TP']+d['FP']}] "
                  f"[{d['rec_pairs']}/{d['truth_pairs']}] genes {d['genes_found']}/{d['genes_total']}")
        print()

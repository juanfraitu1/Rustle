#!/usr/bin/env python3
"""Prereg Addendum M: guided-mode leave-out on a truth table — candidate loci from seed-sequence hits, scored for
family breadth (sensitivity, precision, under/over-merge, pairwise, bipartite) and locus width (Jaccard, boundary
offsets, truncated / overextended).

usage: guided_leaveout.py <truth.tsv> <units.paf> <genes.gff.gz> <out_prefix> [min_identity=0.80] [min_qcov=0.50]
"""
import collections
import csv
import gzip
import math
import random
import re
import statistics
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

truth_path, paf_path, gff_path, out = sys.argv[1:5]
MIN_ID = float(sys.argv[5]) if len(sys.argv) > 5 else 0.80
MIN_QCOV = float(sys.argv[6]) if len(sys.argv) > 6 else 0.50
REPS = 5

truth = list(csv.DictReader(open(truth_path), delimiter="\t"))
for t in truth:
    t["start0"], t["end"] = int(t["start0"]), int(t["end"])
rec = {t["name"]: t for t in truth}
families = sorted({t["family"] for t in truth})

hits = []
for line in open(paf_path):
    f = line.rstrip("\n").split("\t")
    qlen, qs, qe, nm, bl = int(f[1]), int(f[2]), int(f[3]), int(f[9]), int(f[10])
    ident, qcov = nm / bl, (qe - qs) / qlen
    if ident >= MIN_ID and qcov >= MIN_QCOV:
        hits.append({"q": f[0], "chrom": f[5], "s": int(f[7]), "e": int(f[8]), "nm": nm, "ident": ident, "qcov": qcov})

genes = collections.defaultdict(list)  # chrom -> (start0, end, name)
hit_chroms = {h["chrom"] for h in hits}
with gzip.open(gff_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene") or f[0] not in hit_chroms:
            continue
        m = re.search(r"(?:^|;)Name=([^;]+)", f[8])
        genes[f[0]].append((int(f[3]) - 1, int(f[4]), m.group(1) if m else "?"))


def ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def pairwise(pred, true):
    tp = fp = fn = 0
    for i in range(len(pred)):
        for j in range(i + 1, len(pred)):
            sp, st = pred[i] == pred[j], true[i] == true[j]
            tp += sp and st
            fp += sp and not st
            fn += st and not sp
    return (tp / (tp + fn) if tp + fn else float("nan")), (tp / (tp + fp) if tp + fp else float("nan"))


def bipartite(pred, true):
    P, T = sorted(set(pred), key=str), sorted(set(true), key=str)
    M = np.zeros((len(T), len(P)), dtype=int)
    for p, t in zip(pred, true):
        M[T.index(t), P.index(p)] += 1
    r, c = linear_sum_assignment(-M)
    matched = sum(M[i, j] for i, j in zip(r, c))
    size_p = M.sum(axis=0)
    msize = sum(size_p[j] for i, j in zip(r, c) if M[i, j] > 0)
    exact = sum(1 for i, j in zip(r, c) if M[i, j] > 0 and M[i, j] == M[i].sum() == size_p[j])
    return matched / len(pred), (matched / msize if msize else float("nan")), exact, len(T)


def run(level, rep):
    seeds, hidden = set(), set()
    for fi, fam in enumerate(families):
        names = sorted(t["name"] for t in truth if t["family"] == fam)
        random.Random(1000 * rep + fi).shuffle(names)
        k = math.ceil(len(names) / 2) if level == "half" else 1
        seeds.update(names[:k])
        hidden.update(names[k:])
    cand_hits = [h for h in hits if h["q"] in seeds and not any(
        rec[s]["chrom"] == h["chrom"] and ov(rec[s]["start0"], rec[s]["end"], h["s"], h["e"]) > 0 for s in seeds)]
    # single-linkage clusters of overlapping hit spans, per chromosome
    cands = []
    by_chrom = collections.defaultdict(list)
    for h in cand_hits:
        by_chrom[h["chrom"]].append(h)
    for chrom, hs in by_chrom.items():
        hs.sort(key=lambda h: h["s"])
        cluster, cur_end = [], -1
        for h in hs + [None]:
            if h is None or (cluster and h["s"] >= cur_end):
                best = max(cluster, key=lambda x: (x["nm"], -x["s"]))
                cands.append({"chrom": chrom, "s": best["s"], "e": best["e"], "family": rec[best["q"]]["family"],
                              "seed": best["q"]})
                cluster, cur_end = [], -1
            if h is not None:
                cluster.append(h)
                cur_end = max(cur_end, h["e"])
    for c in cands:
        hov = [(ov(rec[n]["start0"], rec[n]["end"], c["s"], c["e"]), n) for n in hidden if rec[n]["chrom"] == c["chrom"]]
        hov = [x for x in hov if x[0] > 0]
        if hov:
            c["match"], c["class"] = max(hov)[1], "hidden"
        else:
            g = [nm for s, e, nm in genes.get(c["chrom"], []) if ov(s, e, c["s"], c["e"]) > 0]
            c["match"], c["class"] = None, ("other_gene" if g else "unannotated")
            c["genes"] = sorted(set(g))
    # hidden record -> best overlapping candidate
    rec_pred = {}
    for n in hidden:
        t = rec[n]
        o = [(ov(t["start0"], t["end"], c["s"], c["e"]), i) for i, c in enumerate(cands) if c["chrom"] == t["chrom"]]
        o = [x for x in o if x[0] > 0]
        rec_pred[n] = cands[max(o)[1]] if o else None
    rows = []
    for fam in families + ["ALL"]:
        H = [n for n in hidden if fam in ("ALL", rec[n]["family"])]
        C = [c for c in cands if fam in ("ALL", c["family"])]
        tp = sum(1 for n in H if rec_pred[n] and rec_pred[n]["family"] == rec[n]["family"])
        wrong = sum(1 for n in H if rec_pred[n] and rec_pred[n]["family"] != rec[n]["family"])
        missed = sum(1 for n in H if not rec_pred[n])
        good_c = sum(1 for c in C if c["class"] == "hidden" and rec[c["match"]]["family"] == c["family"])
        other_fam = sum(1 for c in C if c["class"] == "hidden" and rec[c["match"]]["family"] != c["family"])
        other_gene = [c for c in C if c["class"] == "other_gene"]
        unann = sum(1 for c in C if c["class"] == "unannotated")
        items_pred, items_true = [], []
        for n in H:
            p = rec_pred[n]
            items_pred.append(p["family"] if p else f"missed:{n}")
            items_true.append(rec[n]["family"])
        for i, c in enumerate(C):
            if c["class"] != "hidden":
                items_pred.append(c["family"])
                items_true.append(f"nonmember:{i}")
        ps, pp = pairwise(items_pred, items_true) if items_pred else (float("nan"),) * 2
        bmr, bmp, bex, bnt = bipartite(items_pred, items_true) if items_pred else (float("nan"),) * 4
        widths = []
        for n in H:
            p = rec_pred[n]
            if not p or p["family"] != rec[n]["family"]:
                continue
            t = rec[n]
            L = t["end"] - t["start0"]
            inter = ov(t["start0"], t["end"], p["s"], p["e"])
            union = max(t["end"], p["e"]) - min(t["start0"], p["s"])
            left, right = t["start0"] - p["s"], p["e"] - t["end"]  # positive = beyond truth
            off5, off3 = (left, right) if t["strand"] == "+" else (right, left)
            widths.append((inter / union, off5, off3, inter / L < 0.90, (max(0, left) + max(0, right)) / L > 0.10))
        rows.append({
            "level": level, "rep": rep, "family": fam, "hidden": len(H), "candidates": len(C),
            "sens": tp / len(H) if H else float("nan"), "under_missed": missed, "under_wrong_family": wrong,
            "prec": good_c / len(C) if C else float("nan"), "over_other_family": other_fam,
            "over_other_gene": len(other_gene), "over_unannotated": unann,
            "other_gene_names": ",".join(sorted({g for c in other_gene for g in c["genes"]})),
            "pair_sens": ps, "pair_prec": pp, "bip_micro_R": bmr, "bip_micro_P": bmp, "bip_exact": f"{bex}/{bnt}",
            "n_width": len(widths),
            "jaccard_median": statistics.median([w[0] for w in widths]) if widths else float("nan"),
            "off5_median": statistics.median([w[1] for w in widths]) if widths else float("nan"),
            "off3_median": statistics.median([w[2] for w in widths]) if widths else float("nan"),
            "truncated": sum(w[3] for w in widths), "overextended": sum(w[4] for w in widths),
        })
    return rows, cands, rec_pred, seeds, hidden


all_rows = []
with open(f"{out}.candidates.tsv", "w") as fc:
    fc.write("level\trep\tchrom\tstart\tend\tfamily\tseed\tclass\tmatch\tgenes\n")
    for level in ("half", "keep1"):
        for rep in range(REPS):
            rows, cands, rec_pred, seeds, hidden = run(level, rep)
            all_rows += rows
            for c in cands:
                fc.write(f"{level}\t{rep}\t{c['chrom']}\t{c['s']}\t{c['e']}\t{c['family']}\t{c['seed']}\t{c['class']}\t"
                         f"{c['match'] or ''}\t{','.join(c.get('genes', []))}\n")
cols = list(all_rows[0].keys())
with open(f"{out}.per_rep.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in all_rows:
        fh.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k]) for k in cols) + "\n")


def ms(vals):
    v = [x for x in vals if not (isinstance(x, float) and math.isnan(x))]
    if not v:
        return "nan"
    return f"{statistics.mean(v):.3f}±{statistics.stdev(v):.3f}" if len(v) > 1 else f"{v[0]:.3f}"


print(f"# rule: identity >= {MIN_ID}, query coverage >= {MIN_QCOV}; {len(hits)} hits pass; {REPS} replicates")
for level in ("half", "keep1"):
    print(f"\n## level {level}")
    for fam in families + ["ALL"]:
        R = [r for r in all_rows if r["level"] == level and r["family"] == fam]
        g = lambda k: ms([r[k] for r in R])
        print(f"{fam:7s} hidden {g('hidden')} candidates {g('candidates')} | BREADTH sens {g('sens')} prec {g('prec')} "
              f"| under: missed {g('under_missed')} wrong-family {g('under_wrong_family')} "
              f"| over: other-family {g('over_other_family')} other-gene {g('over_other_gene')} unannotated {g('over_unannotated')} "
              f"| pairwise sens {g('pair_sens')} prec {g('pair_prec')} | bipartite micro R {g('bip_micro_R')} P {g('bip_micro_P')}")
        print(f"{'':7s} WIDTH (n={g('n_width')}): Jaccard median {g('jaccard_median')} | 5' offset median {g('off5_median')} bp "
              f"| 3' offset median {g('off3_median')} bp | truncated {g('truncated')} | overextended {g('overextended')}")
        names = sorted({x for r in R for x in r["other_gene_names"].split(",") if x})
        if names:
            print(f"{'':7s} other genes hit: {', '.join(names[:25])}{' ...' if len(names) > 25 else ''}")

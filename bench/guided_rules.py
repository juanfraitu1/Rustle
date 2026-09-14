#!/usr/bin/env python3
"""Prereg Addendum N: guided-mode rules on a leave-out — M0 (baseline), W1 gene-span projection, W2 isoform union,
I iterative expansion, I+W2. Same truth, leave-out sets, floors and scores as Addendum M (bench/guided_leaveout.py).

usage: guided_rules.py <workdir with truth.tsv, units.paf, isoforms.paf, genespan.paf> <genome.fa> <genome.mmi> <genes.gff.gz>
"""
import collections
import csv
import gzip
import math
import os
import random
import re
import statistics
import subprocess
import sys

import numpy as np
import pysam
from scipy.optimize import linear_sum_assignment

W, FASTA, MMI, GFF = sys.argv[1:5]
MIN_ID, MIN_QCOV, REPS, MAX_ROUNDS = 0.80, 0.50, 5, 10
ARMS = ("M0", "W1", "W2", "I", "I+W2")
FAMNAME = {"NPIP": re.compile(r"nuclear pore complex[- %2C]*interacting protein|NPIP", re.I),
           "TBC1D3": re.compile(r"TBC1 domain family member 3|TBC1D3", re.I)}

truth = list(csv.DictReader(open(f"{W}/truth.tsv"), delimiter="\t"))
for t in truth:
    t["start0"], t["end"] = int(t["start0"]), int(t["end"])
rec = {t["name"]: t for t in truth}
families = sorted({t["family"] for t in truth})
genome = pysam.FastaFile(FASTA)
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def parse_paf(path, gene_of=lambda q: q):
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        qlen, qs, qe, nm, bl = int(f[1]), int(f[2]), int(f[3]), int(f[9]), int(f[10])
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
        h = {"q": f[0], "gene": gene_of(f[0]), "strand": f[4], "chrom": f[5], "s": int(f[7]), "e": int(f[8]),
             "nm": nm, "ident": nm / bl, "qcov": (qe - qs) / qlen, "cg": cg}
        if h["ident"] >= MIN_ID and h["qcov"] >= MIN_QCOV:
            out.append(h)
    return out


U = parse_paf(f"{W}/units.paf")
ISO = parse_paf(f"{W}/isoforms.paf", lambda q: q.split("|")[0])
SPAN = parse_paf(f"{W}/genespan.paf")

genes = collections.defaultdict(list)
with gzip.open(GFF, "rt") as fh:
    for line in fh:
        f = line.split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
            continue
        n = re.search(r"(?:^|;)Name=([^;]+)", f[8])
        d = re.search(r"description=([^;]+)", f[8])
        genes[f[0]].append((int(f[3]) - 1, int(f[4]), n.group(1) if n else "?", d.group(1) if d else ""))


def spliced_target(h):
    """Target sequence of a hit's aligned exon blocks, in the hit's transcript orientation."""
    blocks, pos, cur = [], h["s"], h["s"]
    for n, op in re.findall(r"(\d+)([MIDNSHP=X])", h["cg"]):
        n = int(n)
        if op in "M=XD":
            pos += n
        elif op == "N":
            if pos > cur:
                blocks.append((cur, pos))
            pos += n
            cur = pos
    if pos > cur:
        blocks.append((cur, pos))
    seq = "".join(genome.fetch(h["chrom"], s, e) for s, e in blocks).upper()
    return seq.translate(COMP)[::-1] if h["strand"] == "-" else seq


def make_candidates(hits, blocked, width_mode, round_no, src_family):
    """Cluster passing hits that overlap nothing blocked; one candidate per cluster (best hit decides)."""
    free = [h for h in hits if not any(c == h["chrom"] and ov(s, e, h["s"], h["e"]) > 0 for c, s, e in blocked)]
    by = collections.defaultdict(list)
    for h in free:
        by[h["chrom"]].append(h)
    cands = []
    for chrom, hs in by.items():
        hs.sort(key=lambda h: h["s"])
        cluster, end = [], -1
        for h in hs + [None]:
            if h is None or (cluster and h["s"] >= end):
                best = max(cluster, key=lambda x: (x["nm"], -x["s"]))
                s, e = best["s"], best["e"]
                if width_mode == "W2":
                    same = [x for x in cluster if x["gene"] == best["gene"]]
                    s, e = min(x["s"] for x in same), max(x["e"] for x in same)
                elif width_mode == "W1":
                    for x in SPAN:
                        if x["gene"] == best["gene"] and x["chrom"] == chrom and ov(x["s"], x["e"], s, e) > 0:
                            s, e = min(s, x["s"]), max(e, x["e"])
                cands.append({"chrom": chrom, "s": s, "e": e, "family": src_family(best), "gene": best["gene"],
                              "best": best, "round": round_no})
                cluster, end = [], -1
            if h is not None:
                cluster.append(h)
                end = max(end, h["e"])
    return cands


# ---------------- leave-out states ----------------
states = {}
for level in ("half", "keep1"):
    for rep in range(REPS):
        seeds, hidden = set(), set()
        for fi, fam in enumerate(families):
            names = sorted(t["name"] for t in truth if t["family"] == fam)
            random.Random(1000 * rep + fi).shuffle(names)
            k = math.ceil(len(names) / 2) if level == "half" else 1
            seeds.update(names[:k])
            hidden.update(names[k:])
        blocked0 = [(rec[s]["chrom"], rec[s]["start0"], rec[s]["end"]) for s in seeds]
        for arm in ARMS:
            pool = ISO if arm in ("W2", "I+W2") else U
            wm = "W2" if arm in ("W2", "I+W2") else ("W1" if arm == "W1" else "W0")
            hits = [h for h in pool if h["gene"] in seeds]
            cands = make_candidates(hits, blocked0, wm, 0, lambda b: rec[b["gene"]]["family"])
            states[(level, rep, arm)] = {"seeds": seeds, "hidden": hidden, "blocked0": blocked0, "cands": cands,
                                         "new": list(range(len(cands))), "rounds": [len(cands)]}

# ---------------- iteration ----------------
os.makedirs(f"{W}/iter", exist_ok=True)
for rnd in range(1, MAX_ROUNDS + 1):
    queries = {}
    for key, st in states.items():
        if key[2] not in ("I", "I+W2") or not st["new"]:
            continue
        for ci in st["new"]:
            seq = spliced_target(st["cands"][ci]["best"])
            if len(seq) >= 50:
                queries[f"{key[0]}|{key[1]}|{key[2]}|{rnd}|{ci}"] = seq
        st["new"] = []
    if not queries:
        break
    fa, paf = f"{W}/iter/round{rnd}.fa", f"{W}/iter/round{rnd}.paf"
    with open(fa, "w") as fh:
        for q, s in queries.items():
            fh.write(f">{q}\n{s}\n")
    subprocess.run(["minimap2", "-c", "-x", "splice", "-N", "100", "-p", "0.1", "-t", "4", "-o", paf, MMI, fa],
                   check=True, stderr=subprocess.DEVNULL)
    grouped = collections.defaultdict(list)
    for h in parse_paf(paf):
        level, rep, arm, _, ci = h["q"].split("|")
        grouped[(level, int(rep), arm)].append((int(ci), h))
    added = 0
    for key, lst in grouped.items():
        st = states[key]
        blocked = st["blocked0"] + [(c["chrom"], c["s"], c["e"]) for c in st["cands"]]
        hits = []
        for ci, h in lst:
            h = dict(h)
            h["src_family"] = st["cands"][ci]["family"]
            h["gene"] = st["cands"][ci]["gene"]  # root seed gene, for provenance only
            hits.append(h)
        new = make_candidates(hits, blocked, "W0", rnd, lambda b: b["src_family"])
        base = len(st["cands"])
        st["cands"] += new
        st["new"] = list(range(base, base + len(new)))
        added += len(new)
    for key, st in states.items():
        if key[2] in ("I", "I+W2"):
            st["rounds"].append(len(st["new"]))
    print(f"[iter] round {rnd}: {len(queries)} queries -> {added} new candidates", file=sys.stderr)


# ---------------- scoring ----------------
def pairwise(pred, true):
    tp = fp = fn = 0
    for i in range(len(pred)):
        for j in range(i + 1, len(pred)):
            sp, st_ = pred[i] == pred[j], true[i] == true[j]
            tp += sp and st_
            fp += sp and not st_
            fn += st_ and not sp
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
    return matched / len(pred), (matched / msize if msize else float("nan"))


rows = []
overmerges = []
for (level, rep, arm), st in states.items():
    hidden, cands = st["hidden"], st["cands"]
    for c in cands:
        hov = [(ov(rec[n]["start0"], rec[n]["end"], c["s"], c["e"]), n) for n in hidden if rec[n]["chrom"] == c["chrom"]]
        hov = [x for x in hov if x[0] > 0]
        if hov:
            c["match"], c["class"] = max(hov)[1], "hidden"
        else:
            g = [x for x in genes.get(c["chrom"], []) if ov(x[0], x[1], c["s"], c["e"]) > 0]
            c["match"], c["class"] = None, ("other_gene" if g else "unannotated")
            c["named"] = any(FAMNAME[c["family"]].search(x[2] + " " + x[3]) for x in g)
            c["genes"] = sorted({x[2] for x in g})
    rec_pred = {}
    for n in hidden:
        t = rec[n]
        o = [(ov(t["start0"], t["end"], c["s"], c["e"]), i) for i, c in enumerate(cands) if c["chrom"] == t["chrom"]]
        o = [x for x in o if x[0] > 0]
        rec_pred[n] = cands[max(o)[1]] if o else None
    for fam in families + ["ALL"]:
        H = [n for n in hidden if fam in ("ALL", rec[n]["family"])]
        C = [c for c in cands if fam in ("ALL", c["family"])]
        tp = sum(1 for n in H if rec_pred[n] and rec_pred[n]["family"] == rec[n]["family"])
        wrong = sum(1 for n in H if rec_pred[n] and rec_pred[n]["family"] != rec[n]["family"])
        good = sum(1 for c in C if c["class"] == "hidden" and rec[c["match"]]["family"] == c["family"])
        other_fam = sum(1 for c in C if c["class"] == "hidden" and rec[c["match"]]["family"] != c["family"])
        named = sum(1 for c in C if c["class"] == "other_gene" and c["named"])
        bad_other = [c for c in C if c["class"] == "other_gene" and not c["named"]]
        unann = [c for c in C if c["class"] == "unannotated"]
        if fam != "ALL":
            for c in bad_other + unann + [c for c in C if c["class"] == "hidden" and rec[c["match"]]["family"] != c["family"]]:
                overmerges.append((level, rep, arm, fam, c["chrom"], c["s"], c["e"], c["class"], c.get("match") or ",".join(c.get("genes", []))))
        ip, it = [], []
        for n in H:
            p = rec_pred[n]
            ip.append(p["family"] if p else f"missed:{n}")
            it.append(rec[n]["family"])
        for i, c in enumerate(C):
            if c["class"] != "hidden":
                ip.append(c["family"])
                it.append(f"nonmember:{i}")
        ps, pp = pairwise(ip, it)
        br, bp = bipartite(ip, it)
        wd = []
        for n in H:
            p = rec_pred[n]
            if not p or p["family"] != rec[n]["family"]:
                continue
            t = rec[n]
            L = t["end"] - t["start0"]
            inter = ov(t["start0"], t["end"], p["s"], p["e"])
            union = max(t["end"], p["e"]) - min(t["start0"], p["s"])
            left, right = t["start0"] - p["s"], p["e"] - t["end"]
            o5, o3 = (left, right) if t["strand"] == "+" else (right, left)
            wd.append((inter / union, o5, o3, inter / L < 0.90, (max(0, left) + max(0, right)) / L > 0.10))
        sub = collections.Counter(rec[n]["level1"] for n in H if rec_pred[n] and rec_pred[n]["family"] == rec[n]["family"])
        rows.append({"level": level, "rep": rep, "arm": arm, "family": fam, "hidden": len(H), "candidates": len(C),
                     "sens": tp / len(H) if H else float("nan"), "missed": sum(1 for n in H if not rec_pred[n]),
                     "wrong_family": wrong, "prec": good / len(C) if C else float("nan"),
                     "prec_named": (good + named) / len(C) if C else float("nan"), "over_other_family": other_fam,
                     "over_other_gene_unnamed": len(bad_other), "over_unannotated": len(unann),
                     "pair_sens": ps, "pair_prec": pp, "bip_R": br, "bip_P": bp, "n_width": len(wd),
                     "jaccard": statistics.median([w[0] for w in wd]) if wd else float("nan"),
                     "off5": statistics.median([w[1] for w in wd]) if wd else float("nan"),
                     "off3": statistics.median([w[2] for w in wd]) if wd else float("nan"),
                     "truncated": sum(w[3] for w in wd), "overextended": sum(w[4] for w in wd),
                     "rounds": ",".join(map(str, st["rounds"])),
                     "subfamilies_recovered": ";".join(f"{k}:{v}" for k, v in sorted(sub.items()))})

cols = list(rows[0].keys())
with open(f"{W}/n.per_rep.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k]) for k in cols) + "\n")
with open(f"{W}/n.overmerges.tsv", "w") as fh:
    fh.write("level\trep\tarm\tfamily\tchrom\tstart\tend\tclass\tmatch_or_genes\n")
    for o in overmerges:
        fh.write("\t".join(map(str, o)) + "\n")


def ms(v):
    v = [x for x in v if not (isinstance(x, float) and math.isnan(x))]
    return "nan" if not v else (f"{statistics.mean(v):.3f}±{statistics.stdev(v):.3f}" if len(v) > 1 else f"{v[0]:.3f}")


# width ceilings: each record's self projection under each width rule
def self_cov(pool, name, mode):
    t = rec[name]
    hs = [h for h in pool if h["gene"] == name and h["chrom"] == t["chrom"] and ov(h["s"], h["e"], t["start0"], t["end"]) > 0]
    if not hs:
        return 0.0
    if mode == "W2":
        s, e = min(h["s"] for h in hs), max(h["e"] for h in hs)
    else:
        b = max(hs, key=lambda h: ov(h["s"], h["e"], t["start0"], t["end"]))
        s, e = b["s"], b["e"]
        if mode == "W1":
            for x in SPAN:
                if x["gene"] == name and x["chrom"] == t["chrom"] and ov(x["s"], x["e"], s, e) > 0:
                    s, e = min(s, x["s"]), max(e, x["e"])
    return ov(s, e, t["start0"], t["end"]) / (t["end"] - t["start0"])


print(f"# Addendum N: floors identity >= {MIN_ID}, query coverage >= {MIN_QCOV}; {REPS} replicates")
for mode, pool in (("W0", U), ("W1", U), ("W2", ISO)):
    cov = [self_cov(pool, t["name"], mode) for t in truth]
    print(f"width ceiling {mode}: records with self-projection covering >= 0.90 of own gene span: "
          f"{sum(c >= 0.90 for c in cov)}/{len(cov)}; median {statistics.median(cov):.3f}")
for level in ("half", "keep1"):
    print(f"\n## level {level}")
    for fam in families + ["ALL"]:
        for arm in ARMS:
            R = [r for r in rows if r["level"] == level and r["family"] == fam and r["arm"] == arm]
            g = lambda k: ms([r[k] for r in R])
            print(f"{fam:6s} {arm:5s} hid {g('hidden')} cand {g('candidates')} | sens {g('sens')} prec {g('prec')} "
                  f"prec_named {g('prec_named')} | missed {g('missed')} wrongfam {g('wrong_family')} | over: otherfam "
                  f"{g('over_other_family')} unnamed-gene {g('over_other_gene_unnamed')} unannot {g('over_unannotated')} "
                  f"| pair {g('pair_sens')}/{g('pair_prec')} bip {g('bip_R')}/{g('bip_P')} | W n={g('n_width')} "
                  f"J {g('jaccard')} 5' {g('off5')} 3' {g('off3')} trunc {g('truncated')} overext {g('overextended')}")
        if fam == "NPIP":
            for arm in ("I", "I+W2"):
                R = [r for r in rows if r["level"] == level and r["family"] == fam and r["arm"] == arm]
                print(f"        {arm} rounds (new candidates per round) {[r['rounds'] for r in R]}; subfamilies "
                      f"{[r['subfamilies_recovered'] for r in R]}")
print(f"\nover-merge candidates listed: {len(overmerges)} -> {W}/n.overmerges.tsv")

#!/usr/bin/env python3
"""Prereg Addendum O: guided mode at the literature's two levels.
G1 — family members found by the seed GENE BODY (chained asm20 alignments of the gene span), width clip / extrap.
G2 — subfamilies inside a G1 family from exon-masked gene-body identity, partitioned by the identity-gap rule of
bench/identity_gap.py (largest interior gap vs gauss/beta/smooth nulls, worst p governs; numpy draws).

usage: guided_genebody.py <workdir: truth.tsv units.paf genespan.paf> <refseq gff (uncompressed, has isoforms)>
       <genes.gff.gz (full, for classification)> <genome.fa>
"""
import collections
import csv
import gzip
import math
import random
import re
import statistics
import subprocess
import sys

import numpy as np
import pysam
from scipy.optimize import linear_sum_assignment

W, GFF_ISO, GFF_FULL, FASTA = sys.argv[1:5]
MIN_ID, MIN_COV, REPS, NULL_DRAWS = 0.80, 0.50, 5, 10000
FAMNAME = {"NPIP": re.compile(r"nuclear pore complex[- %2C]*interacting protein|NPIP", re.I),
           "TBC1D3": re.compile(r"TBC1 domain family member 3|TBC1D3", re.I)}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
genome = pysam.FastaFile(FASTA)

truth = list(csv.DictReader(open(f"{W}/truth.tsv"), delimiter="\t"))
for t in truth:
    t["start0"], t["end"] = int(t["start0"]), int(t["end"])
rec = {t["name"]: t for t in truth}
families = sorted({t["family"] for t in truth})


def ov(a0, a1, b0, b1):
    return max(0, min(a1, b1) - max(a0, b0))


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


# ---------------- union exons per truth gene (all isoforms) ----------------
gene_id, tx_parent, exons_by_parent = {}, {}, collections.defaultdict(list)
for line in open(GFF_ISO):
    if line.startswith("#"):
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 9:
        continue
    a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
    if f[2] in ("gene", "pseudogene") and a.get("Name") in rec:
        gene_id[a["ID"]] = a["Name"]
    elif f[2] in ("mRNA", "transcript", "lnc_RNA", "ncRNA", "primary_transcript"):
        tx_parent[a["ID"]] = a.get("Parent", "")
    elif f[2] == "exon":
        exons_by_parent[a.get("Parent", "")].append((int(f[3]) - 1, int(f[4])))
union_exons = {}
for gid, name in gene_id.items():
    ex = [e for tx, p in tx_parent.items() if p == gid for e in exons_by_parent.get(tx, [])]
    union_exons[name] = merge(ex or exons_by_parent.get(gid, []))
for t in truth:
    union_exons.setdefault(t["name"], [])


def exons_in_query(name):
    t = rec[name]
    out = []
    for s, e in union_exons[name]:
        s, e = max(s, t["start0"]), min(e, t["end"])
        if e > s:
            out.append((s - t["start0"], e - t["start0"]) if t["strand"] == "+" else (t["end"] - e, t["end"] - s))
    return merge(out)


# ---------------- genes for classification ----------------
genes = collections.defaultdict(list)
with gzip.open(GFF_FULL, "rt") as fh:
    for line in fh:
        f = line.split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
            continue
        n = re.search(r"(?:^|;)Name=([^;]+)", f[8])
        d = re.search(r"description=([^;]+)", f[8])
        genes[f[0]].append((int(f[3]) - 1, int(f[4]), n.group(1) if n else "?", d.group(1) if d else ""))


# ---------------- M0 hits (transcript units) ----------------
def parse_paf(path):
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        out.append({"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4], "chrom": f[5],
                    "clen": int(f[6]), "ts": int(f[7]), "te": int(f[8]), "nm": int(f[9]), "bl": int(f[10]),
                    "cg": next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")})
    return out


M0_hits = [dict(h, s=h["ts"], e=h["te"]) for h in parse_paf(f"{W}/units.paf")
           if h["nm"] / h["bl"] >= MIN_ID and (h["qe"] - h["qs"]) / h["qlen"] >= MIN_COV]

# ---------------- G1 chains ----------------
chains = []
by_key = collections.defaultdict(list)
for r in parse_paf(f"{W}/genespan.paf"):
    by_key[(r["q"], r["chrom"], r["strand"])].append(r)
for (q, chrom, strand), rs in by_key.items():
    rs.sort(key=lambda r: (r["ts"], r["te"]))
    Lq = rs[0]["qlen"]
    cur = []
    for r in rs + [None]:
        ok = False
        if r is not None and cur:
            last = cur[-1]
            gap = r["ts"] - max(x["te"] for x in cur)
            span = r["te"] - min(x["ts"] for x in cur)
            order = r["qs"] >= last["qs"] if strand == "+" else r["qs"] <= last["qs"]
            ok = gap <= Lq and span <= 2 * Lq and order
        if r is None or (cur and not ok):
            nm, bl = sum(x["nm"] for x in cur), sum(x["bl"] for x in cur)
            qiv = merge([(x["qs"], x["qe"]) for x in cur])
            aligned = sum(e - s for s, e in qiv)
            qmin, qmax = qiv[0][0], qiv[-1][1]
            ts, te = min(x["ts"] for x in cur), max(x["te"] for x in cur)
            if strand == "+":
                xs, xe = ts - qmin, te + (Lq - qmax)
            else:
                xs, xe = ts - (Lq - qmax), te + qmin
            clen = cur[0]["clen"]
            xs, xe = max(0, xs), min(clen, xe)
            Lt = xe - xs
            chains.append({"q": q, "chrom": chrom, "strand": strand, "recs": cur, "nm": nm, "ident": nm / bl,
                           "aligned": aligned, "Lq": Lq, "Lt": Lt, "s": ts, "e": te, "xs": xs, "xe": xe,
                           "pass": nm / bl >= MIN_ID and aligned >= MIN_COV * min(Lq, Lt)})
            cur = []
        if r is not None:
            cur.append(r)
G1_hits = [c for c in chains if c["pass"]]
print(f"[G1] {len(chains)} chains from {sum(len(v) for v in by_key.values())} records; {len(G1_hits)} pass", file=sys.stderr)


def project_exons(chain):
    """Target intervals of the seed's union exons, mapped base-to-base through the chain's CIGARs."""
    qex = exons_in_query(chain["q"])
    out = []
    for r in chain["recs"]:
        t, q = r["ts"], (r["qs"] if r["strand"] == "+" else r["qe"])
        for n, op in re.findall(r"(\d+)([MID])", r["cg"]):
            n = int(n)
            if op == "M":
                qa, qb = (q, q + n) if r["strand"] == "+" else (q - n, q)
                for es, ee in qex:
                    lo, hi = max(es, qa), min(ee, qb)
                    if hi > lo:
                        if r["strand"] == "+":
                            out.append((t + (lo - qa), t + (hi - qa)))
                        else:
                            out.append((t + (qb - hi), t + (qb - lo)))
                t += n
                q = q + n if r["strand"] == "+" else q - n
            elif op == "I":
                q = q + n if r["strand"] == "+" else q - n
            else:
                t += n
    return merge(out)


# ---------------- candidates ----------------
def candidates(hits, seeds, nm_key, span_key):
    blocked = [(rec[s]["chrom"], rec[s]["start0"], rec[s]["end"]) for s in seeds]
    free = [h for h in hits if h["q"] in seeds and not any(c == h["chrom"] and ov(s, e, h["s"], h["e"]) > 0 for c, s, e in blocked)]
    by = collections.defaultdict(list)
    for h in free:
        by[h["chrom"]].append(h)
    out = []
    for chrom, hs in by.items():
        hs.sort(key=lambda h: h["s"])
        cl, end = [], -1
        for h in hs + [None]:
            if h is None or (cl and h["s"] >= end):
                b = max(cl, key=lambda x: (x[nm_key], -x["s"]))
                s, e = span_key(b)
                out.append({"chrom": chrom, "s": s, "e": e, "clip": (b["s"], b["e"]), "family": rec[b["q"]]["family"],
                            "seed": b["q"], "best": b})
                cl, end = [], -1
            if h is not None:
                cl.append(h)
                end = max(end, h["e"])
    return out


def classify(cands, hidden):
    for c in cands:
        hov = [(ov(rec[n]["start0"], rec[n]["end"], c["s"], c["e"]), n) for n in hidden if rec[n]["chrom"] == c["chrom"]]
        hov = [x for x in hov if x[0] > 0]
        if hov:
            c["match"], c["class"], c["named"] = max(hov)[1], "hidden", False
        else:
            g = [x for x in genes.get(c["chrom"], []) if ov(x[0], x[1], c["s"], c["e"]) > 0]
            c["match"], c["class"] = None, ("other_gene" if g else "unannotated")
            c["named"] = any(FAMNAME[c["family"]].search(x[2] + " " + x[3]) for x in g)
            c["genes"] = sorted({x[2] for x in g})
    rp = {}
    for n in hidden:
        t = rec[n]
        o = [(ov(t["start0"], t["end"], c["s"], c["e"]), i) for i, c in enumerate(cands) if c["chrom"] == t["chrom"]]
        o = [x for x in o if x[0] > 0]
        rp[n] = cands[max(o)[1]] if o else None
    return rp


def pairwise(pred, true):
    tp = fp = fn = 0
    for i in range(len(pred)):
        for j in range(i + 1, len(pred)):
            a, b = pred[i] == pred[j], true[i] == true[j]
            tp += a and b
            fp += a and not b
            fn += b and not a
    return (tp / (tp + fn) if tp + fn else float("nan")), (tp / (tp + fp) if tp + fp else float("nan"))


def bipartite(pred, true):
    P, T = sorted(set(pred), key=str), sorted(set(true), key=str)
    M = np.zeros((len(T), len(P)), dtype=int)
    for p, t in zip(pred, true):
        M[T.index(t), P.index(p)] += 1
    r, c = linear_sum_assignment(-M)
    matched = sum(M[i, j] for i, j in zip(r, c))
    sp = M.sum(axis=0)
    msize = sum(sp[j] for i, j in zip(r, c) if M[i, j] > 0)
    exact = sum(1 for i, j in zip(r, c) if M[i, j] > 0 and M[i, j] == M[i].sum() == sp[j])
    return matched / len(pred), (matched / msize if msize else float("nan")), exact, len(T)


def width_stats(t, s, e):
    L = t["end"] - t["start0"]
    inter = ov(t["start0"], t["end"], s, e)
    left, right = t["start0"] - s, e - t["end"]
    o5, o3 = (left, right) if t["strand"] == "+" else (right, left)
    return (inter / (max(t["end"], e) - min(t["start0"], s)), o5, o3, inter / L < 0.90,
            (max(0, left) + max(0, right)) / L > 0.10)


# ---------------- G2 ----------------
def seq_minus(chrom, s, e, remove, strand):
    pieces, cur = [], s
    for a, b in merge([(max(s, a), min(e, b)) for a, b in remove if min(e, b) > max(s, a)]):
        if a > cur:
            pieces.append((cur, a))
        cur = max(cur, b)
    if e > cur:
        pieces.append((cur, e))
    seq = "".join(genome.fetch(chrom, a, b) for a, b in pieces).upper()
    return seq.translate(COMP)[::-1] if strand == "-" else seq


def largest_gap(xs):
    lo, hi = max(1, int(0.1 * len(xs))), min(len(xs) - 1, int(0.9 * len(xs)))
    best, at = 0.0, None
    for i in range(lo, hi):
        d = xs[i] - xs[i - 1]
        if d > best:
            best, at = d, (xs[i - 1] + xs[i]) / 2
    return best, at


def gap_matrix(M):
    n = M.shape[1]
    lo, hi = max(1, int(0.1 * n)), min(n - 1, int(0.9 * n))
    if hi <= lo:
        return np.zeros(M.shape[0])
    S = np.sort(M, axis=1)
    return np.diff(S, axis=1)[:, lo - 1:hi - 1].max(axis=1)


def identity_gap_partition(pairs, seed=0):
    v = sorted(p[2] for p in pairs)
    n = len(v)
    members = sorted({x for a, b, _ in pairs for x in (a, b)})
    if n < 6:
        return None, None, [members], "too few pairs"
    obs, cut = largest_gap(v)
    rng = np.random.default_rng(seed)
    mu, sd = statistics.mean(v), statistics.pstdev(v)
    pv = {"gauss": float(np.mean(gap_matrix(np.clip(rng.normal(mu, sd, (NULL_DRAWS, n)), 0, 1)) >= obs))}
    if sd > 0 and 0 < mu < 1:
        tt = mu * (1 - mu) / (sd * sd) - 1
        if tt > 0:
            pv["beta"] = float(np.mean(gap_matrix(rng.beta(mu * tt, (1 - mu) * tt, (NULL_DRAWS, n))) >= obs))
        q = statistics.quantiles(v, n=4)
        h = 0.9 * min(sd, (q[2] - q[0]) / 1.34) * n ** -0.2
        if h > 0:
            Mv = np.array(v)[rng.integers(0, n, (NULL_DRAWS, n))] + rng.normal(0, h, (NULL_DRAWS, n))
            pv["smooth"] = float(np.mean(gap_matrix(np.clip(Mv, 0, 1)) >= obs))
    p = max(pv.values())
    if p >= 0.05:
        return p, cut, [members], f"no split (gap {obs:.4f} at {cut:.4f})"
    parent = {m: m for m in members}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b, i in pairs:
        if i >= cut:
            parent[find(a)] = find(b)
    comp = collections.defaultdict(list)
    for m in members:
        comp[find(m)].append(m)
    return p, cut, sorted(comp.values(), key=len, reverse=True), f"SPLIT (gap {obs:.4f} at {cut:.4f})"


def g2_run(tag, member_seqs):
    fa = f"{W}/g2/{tag}.fa"
    with open(fa, "w") as fh:
        for k, s in member_seqs.items():
            fh.write(f">{k}\n{s}\n")
    paf = subprocess.run(["minimap2", "-c", "-x", "asm20", "-X", "-N", "50", "-p", "0.1", "-t", "4", fa, fa],
                         capture_output=True, text=True, check=True).stdout
    nm, bl = collections.Counter(), collections.Counter()
    for line in paf.splitlines():
        f = line.split("\t")
        if f[0] == f[5]:
            continue
        key = tuple(sorted((f[0], f[5])))
        nm[key] += int(f[9])
        bl[key] += int(f[10])
    pairs = [(a, b, nm[(a, b)] / bl[(a, b)]) for (a, b) in bl]
    unaligned = len(member_seqs) * (len(member_seqs) - 1) // 2 - len(pairs)
    p, cut, groups, verdict = identity_gap_partition(pairs)
    for m in member_seqs:
        if not any(m in g for g in groups):
            groups.append([m])
    return p, cut, groups, verdict, len(pairs), unaligned


import os
os.makedirs(f"{W}/g2", exist_ok=True)


def score_partition(groups, label_of):
    """label_of: member -> truth record name (or None for non-truth members)."""
    lab = {m: gi for gi, g in enumerate(groups) for m in g}
    items = [(m, label_of[m]) for m in lab if label_of.get(m)]
    out = {}
    for lvl in ("level1", "level2"):
        pred = [lab[m] for m, _ in items]
        true = [rec[r][lvl] for _, r in items]
        ps, pp = pairwise(pred, true)
        br, bp, ex, nt = bipartite(pred, true) if items else (float("nan"),) * 4
        out[lvl] = (ps, pp, br, bp, ex, nt)
    return out, len(items)


# reference: all truth records of each family, annotation exons removed
print("## G2 REFERENCE (no leave-out; all truth gene bodies, annotated exons removed)")
ref_results = {}
for fam in families:
    seqs = {t["name"]: seq_minus(t["chrom"], t["start0"], t["end"], union_exons[t["name"]], t["strand"])
            for t in truth if t["family"] == fam}
    p, cut, groups, verdict, npairs, unal = g2_run(f"ref_{fam}", seqs)
    sc, n = score_partition(groups, {k: k for k in seqs})
    ref_results[fam] = (p, verdict, groups, sc)
    print(f"{fam}: {n} members, {npairs} aligned pairs ({unal} unaligned); p = {p}; {verdict}")
    for g in groups:
        print("   group:", " ".join(f"{m}[{rec[m]['level1']}/{rec[m]['level2']}]" for m in sorted(g)))
    for lvl in ("level1", "level2"):
        ps, pp, br, bp, ex, nt = sc[lvl]
        print(f"   vs {lvl}: pairwise sens {ps:.3f} prec {pp:.3f} | bipartite micro R {br:.3f} P {bp:.3f} exact {ex}/{nt}")

# ---------------- leave-out ----------------
rows, g2rows = [], []
for level in ("half", "keep1"):
    for rep in range(REPS):
        seeds, hidden = set(), set()
        for fi, fam in enumerate(families):
            names = sorted(t["name"] for t in truth if t["family"] == fam)
            random.Random(1000 * rep + fi).shuffle(names)
            k = math.ceil(len(names) / 2) if level == "half" else 1
            seeds.update(names[:k])
            hidden.update(names[k:])
        arms = {
            "M0": candidates(M0_hits, seeds, "nm", lambda b: (b["s"], b["e"])),
            "G1": candidates(G1_hits, seeds, "nm", lambda b: (b["xs"], b["xe"])),
        }
        for arm, cands in arms.items():
            rp = classify(cands, hidden)
            for fam in families + ["ALL"]:
                H = [n for n in hidden if fam in ("ALL", rec[n]["family"])]
                C = [c for c in cands if fam in ("ALL", c["family"])]
                tp = sum(1 for n in H if rp[n] and rp[n]["family"] == rec[n]["family"])
                good = sum(1 for c in C if c["class"] == "hidden" and rec[c["match"]]["family"] == c["family"])
                named = sum(1 for c in C if c["class"] == "other_gene" and c["named"])
                crossfam = sum(1 for c in C if c["class"] == "hidden" and rec[c["match"]]["family"] != c["family"])
                unnamed = [",".join(c["genes"]) for c in C if c["class"] == "other_gene" and not c["named"]]
                unann = sum(1 for c in C if c["class"] == "unannotated")
                ip = [rp[n]["family"] if rp[n] else f"missed:{n}" for n in H] + [c["family"] for c in C if c["class"] != "hidden"]
                it = [rec[n]["family"] for n in H] + [f"nonmember:{i}" for i, c in enumerate(C) if c["class"] != "hidden"]
                ps, pp = pairwise(ip, it)
                br, bp, _, _ = bipartite(ip, it) if ip else (float("nan"),) * 4
                widths = {}
                for wname in (("W0",) if arm == "M0" else ("clip", "extrap")):
                    wd = []
                    for n in H:
                        c = rp[n]
                        if not c or c["family"] != rec[n]["family"]:
                            continue
                        s, e = (c["s"], c["e"]) if wname in ("W0", "extrap") else c["clip"]
                        wd.append(width_stats(rec[n], s, e))
                    widths[wname] = wd
                sub = collections.Counter(rec[n]["level1"] for n in H if rp[n] and rp[n]["family"] == rec[n]["family"])
                for wname, wd in widths.items():
                    rows.append({"level": level, "rep": rep, "arm": arm if arm == "M0" else f"G1-{wname}", "family": fam,
                                 "hidden": len(H), "candidates": len(C), "sens": tp / len(H) if H else float("nan"),
                                 "missed": sum(1 for n in H if not rp[n]), "prec": good / len(C) if C else float("nan"),
                                 "prec_named": (good + named) / len(C) if C else float("nan"), "cross_family": crossfam,
                                 "unnamed_other_gene": len(unnamed), "unannotated": unann, "unnamed_list": ";".join(unnamed),
                                 "pair_sens": ps, "pair_prec": pp, "bip_R": br, "bip_P": bp, "n_width": len(wd),
                                 "jaccard": statistics.median([w[0] for w in wd]) if wd else float("nan"),
                                 "off5": statistics.median([w[1] for w in wd]) if wd else float("nan"),
                                 "off3": statistics.median([w[2] for w in wd]) if wd else float("nan"),
                                 "truncated": sum(w[3] for w in wd), "overextended": sum(w[4] for w in wd),
                                 "subfamilies": ";".join(f"{k}:{v}" for k, v in sorted(sub.items()))})
            if arm != "G1":
                continue
            # G2 inside each G1 family
            for fam in families:
                members, label_of = {}, {}
                for s in seeds:
                    if rec[s]["family"] == fam:
                        t = rec[s]
                        members[s] = seq_minus(t["chrom"], t["start0"], t["end"], union_exons[s], t["strand"])
                        label_of[s] = s
                for i, c in enumerate(cands):
                    if c["family"] != fam:
                        continue
                    key = f"cand{i}"
                    members[key] = seq_minus(c["chrom"], c["s"], c["e"], project_exons(c["best"]), c["best"]["strand"])
                    if c["class"] == "hidden" and rec[c["match"]]["family"] == fam:
                        label_of[key] = c["match"]
                if len(members) < 3:
                    g2rows.append({"level": level, "rep": rep, "family": fam, "members": len(members), "verdict": "too few members"})
                    continue
                p, cut, groups, verdict, npairs, unal = g2_run(f"{level}_{rep}_{fam}", members)
                sc, n_truth = score_partition(groups, label_of)
                named_groups = [sorted(label_of.get(m, m) for m in g) for g in groups]
                g2rows.append({"level": level, "rep": rep, "family": fam, "members": len(members), "truth_members": n_truth,
                               "p": p, "verdict": verdict, "sc": sc, "groups": named_groups})

# ---------------- report ----------------
def ms(v):
    v = [x for x in v if not (isinstance(x, float) and math.isnan(x))]
    return "nan" if not v else (f"{statistics.mean(v):.3f}±{statistics.stdev(v):.3f}" if len(v) > 1 else f"{v[0]:.3f}")


print(f"\n## G1 breadth and width (floors identity >= {MIN_ID}, coverage >= {MIN_COV} of min(L_q, L_t); {REPS} replicates)")
for level in ("half", "keep1"):
    print(f"\n### level {level}")
    for fam in families + ["ALL"]:
        for arm in ("M0", "G1-clip", "G1-extrap"):
            R = [r for r in rows if r["level"] == level and r["family"] == fam and r["arm"] == arm]
            g = lambda k: ms([r[k] for r in R])
            print(f"{fam:6s} {arm:9s} hid {g('hidden')} cand {g('candidates')} | sens {g('sens')} prec {g('prec')} named {g('prec_named')} "
                  f"| missed {g('missed')} crossfam {g('cross_family')} unnamed {g('unnamed_other_gene')} unannot {g('unannotated')} "
                  f"| pair {g('pair_sens')}/{g('pair_prec')} bip {g('bip_R')}/{g('bip_P')} | W n={g('n_width')} J {g('jaccard')} "
                  f"5' {g('off5')} 3' {g('off3')} trunc {g('truncated')} overext {g('overextended')}")
            if arm == "G1-extrap" and fam == "NPIP":
                print(f"          subfamilies per rep: {[r['subfamilies'] for r in R]}")
            lists = sorted({x for r in R for x in r["unnamed_list"].split(";") if x})
            if lists and arm != "G1-clip":
                print(f"          unnamed other-gene candidates: {lists[:10]}")

print("\n## G2 inside G1 families (leave-out)")
for r in g2rows:
    if "sc" not in r:
        print(f"{r['level']} rep{r['rep']} {r['family']}: {r['members']} members -> {r['verdict']}")
        continue
    s1, s2 = r["sc"]["level1"], r["sc"]["level2"]
    recovered = r["p"] is not None and r["p"] < 0.05 and s1[4] == 2 and s1[5] == 2
    print(f"{r['level']} rep{r['rep']} {r['family']}: {r['members']} members ({r['truth_members']} truth), p={r['p']}, {r['verdict']} "
          f"| L1 pair {s1[0]:.3f}/{s1[1]:.3f} bip {s1[2]:.3f}/{s1[3]:.3f} exact {s1[4]}/{s1[5]} "
          f"| L2 pair {s2[0]:.3f}/{s2[1]:.3f} bip {s2[2]:.3f}/{s2[3]:.3f} exact {s2[4]}/{s2[5]} | L1 RECOVERED={recovered}")
    print("      groups:", " | ".join(" ".join(g) for g in r["groups"]))
for fam in families:
    R = [r for r in g2rows if r["family"] == fam and "sc" in r]
    rec_n = sum(1 for r in R if r["p"] is not None and r["p"] < 0.05 and r["sc"]["level1"][4] == 2 and r["sc"]["level1"][5] == 2)
    split_n = sum(1 for r in R if r["p"] is not None and r["p"] < 0.05)
    print(f"G2 {fam}: level-1 recovered in {rec_n}/{len(R)} leave-out runs; split called in {split_n}/{len(R)}")

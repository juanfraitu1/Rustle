#!/usr/bin/env python3
"""Guided O1 pipeline with the Addendum T fixes, on a leave-out of an annotated truth table.

From an annotation (the seeds) it finds new candidate loci with two finders — the seed TRANSCRIPT (spliced alignment)
and the seed GENE BODY (F2: CDS envelope, asm20 chains) — builds candidates from both without fusing neighbouring
copies (F1: reciprocal-overlap leader clustering), takes each candidate's width from its transcript hit when it has
one and from its gene body otherwise, and calls subfamilies as supported clades of exon and intron trees built from a
reference-projected alignment (F3: the reference member is the one aligned to the most members) with IQ-TREE.

usage: guided_pipeline.py --workdir DIR --gff refseq.gff --genes-gff-gz full.gff.gz --genome genome.fa
                          --mmi genome.mmi --iqtree iqtree3 [--expected-units units.fa] [--reps 5] [--threads 4]
DIR/truth.tsv columns: family name chrom start0 end strand biotype level1 level2 (width = start0..end).
Outputs in DIR: units.fa/.paf, envelope.fa/.paf, t.candidates.tsv, t.out (report), tree_t/ (alignments, trees).
Evaluated in docs/o1_ledger.md §6js (prereg Addendum T).
"""
import argparse
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

MIN_ID, MIN_COV, RECIP, CONTAIN = 0.80, 0.50, 0.50, 0.90
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
FAMNAME = {"NPIP": re.compile(r"nuclear pore complex[- %2C]*interacting protein|NPIP", re.I),
           "TBC1D3": re.compile(r"TBC1 domain family member 3|TBC1D3", re.I),
           "AMY": re.compile(r"(?<!gluco)amylase|\bAMY", re.I)}


# ---------------------------------------------------------------- small helpers
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


def recip(a0, a1, b0, b1):
    o = ov(a0, a1, b0, b1)
    return o >= RECIP * (a1 - a0) and o >= RECIP * (b1 - b0)


def rc(s):
    return s.translate(COMP)[::-1]


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
    return matched / len(pred), (matched / msize if msize else float("nan"))


def ms(v):
    v = [x for x in v if not (isinstance(x, float) and math.isnan(x))]
    return "nan" if not v else (f"{statistics.mean(v):.3f}±{statistics.stdev(v):.3f}" if len(v) > 1 else f"{v[0]:.3f}")


# ---------------------------------------------------------------- annotation
class Annotation:
    def __init__(self, gff, truth, genome):
        self.truth, self.genome = truth, genome
        self.rec = {t["name"]: t for t in truth}
        gene_id, tx_parent, tx_name = {}, {}, {}
        exons, cds = collections.defaultdict(list), collections.defaultdict(list)
        for line in open(gff):
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] in ("gene", "pseudogene") and a.get("Name") in self.rec:
                gene_id[a["ID"]] = a["Name"]
            elif f[2] in ("mRNA", "transcript", "lnc_RNA", "ncRNA", "primary_transcript"):
                tx_parent[a["ID"]] = a.get("Parent", "")
                tx_name[a["ID"]] = a.get("Name", a["ID"])
            elif f[2] == "exon":
                exons[a.get("Parent", "")].append((int(f[3]) - 1, int(f[4])))
            elif f[2] == "CDS":
                cds[a.get("Parent", "")].append((int(f[3]) - 1, int(f[4])))
        self.unit_seq, self.unit_model, self.union_exons, self.envelope = {}, {}, {}, {}
        for name, t in self.rec.items():
            gids = [g for g, n in gene_id.items() if n == name]
            txs = [tx for tx, p in tx_parent.items() if p in gids and exons.get(tx)]
            span = lambda tx: sum(e - s for s, e in merge(exons[tx]))
            pool = [tx for tx in txs if tx_name[tx].startswith(("NM_", "NR_"))] or txs
            if pool:
                tx = max(pool, key=lambda x: (span(x), x))
                blocks, model = merge(exons[tx]), tx_name[tx]
            else:
                blocks, model = merge([e for g in gids for e in exons.get(g, [])]), "gene-exons"
            self.union_exons[name] = merge([e for tx in txs for e in exons[tx]] or blocks)
            if blocks:
                seq = "".join(genome.fetch(t["chrom"], s, e) for s, e in blocks).upper()
            else:
                seq, model = genome.fetch(t["chrom"], t["start0"], t["end"]).upper(), "gene span"
            self.unit_seq[name] = rc(seq) if t["strand"] == "-" else seq
            self.unit_model[name] = model
            cd = [c for tx in txs for c in cds.get(tx, [])]
            if cd:
                env, kind = (min(s for s, _ in cd), max(e for _, e in cd)), "CDS envelope"
            elif blocks:
                env, kind = (blocks[0][0], blocks[-1][1]), "exon envelope"
            else:
                env, kind = (t["start0"], t["end"]), "gene span"
            self.envelope[name] = (env[0], env[1], kind)

    def envelope_seq(self, name):
        t = self.rec[name]
        s, e, _ = self.envelope[name]
        seq = self.genome.fetch(t["chrom"], s, e).upper()
        return rc(seq) if t["strand"] == "-" else seq

    def exons_in_envelope(self, name):
        t = self.rec[name]
        s0, e0, _ = self.envelope[name]
        out = []
        for s, e in self.union_exons[name]:
            s, e = max(s, s0), min(e, e0)
            if e > s:
                out.append((s - s0, e - s0) if t["strand"] == "+" else (e0 - e, e0 - s))
        return merge(out)

    def seq_minus(self, chrom, s, e, remove, strand):
        pieces, cur = [], s
        for a, b in merge([(max(s, a), min(e, b)) for a, b in remove if min(e, b) > max(s, a)]):
            if a > cur:
                pieces.append((cur, a))
            cur = max(cur, b)
        if e > cur:
            pieces.append((cur, e))
        seq = "".join(self.genome.fetch(chrom, a, b) for a, b in pieces).upper()
        return rc(seq) if strand == "-" else seq


# ---------------------------------------------------------------- alignments
def parse_paf(path):
    out = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        out.append({"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4], "chrom": f[5],
                    "clen": int(f[6]), "ts": int(f[7]), "te": int(f[8]), "nm": int(f[9]), "bl": int(f[10]),
                    "cg": next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")})
    return out


def transcript_hits(paf):
    return [dict(h, s=h["ts"], e=h["te"]) for h in parse_paf(paf)
            if h["nm"] / h["bl"] >= MIN_ID and (h["qe"] - h["qs"]) / h["qlen"] >= MIN_COV]


def gene_body_chains(paf):
    by_key = collections.defaultdict(list)
    for r in parse_paf(paf):
        by_key[(r["q"], r["chrom"], r["strand"])].append(r)
    chains = []
    for (q, chrom, strand), rs in by_key.items():
        rs.sort(key=lambda r: (r["ts"], r["te"]))
        Lq, cur = rs[0]["qlen"], []
        for r in rs + [None]:
            ok = False
            if r is not None and cur:
                gap = r["ts"] - max(x["te"] for x in cur)
                span = r["te"] - min(x["ts"] for x in cur)
                order = r["qs"] >= cur[-1]["qs"] if strand == "+" else r["qs"] <= cur[-1]["qs"]
                ok = gap <= Lq and span <= 2 * Lq and order
            if r is None or (cur and not ok):
                nm, bl = sum(x["nm"] for x in cur), sum(x["bl"] for x in cur)
                qiv = merge([(x["qs"], x["qe"]) for x in cur])
                aligned = sum(e - s for s, e in qiv)
                ts, te = min(x["ts"] for x in cur), max(x["te"] for x in cur)
                xs, xe = (ts - qiv[0][0], te + (Lq - qiv[-1][1])) if strand == "+" else (ts - (Lq - qiv[-1][1]), te + qiv[0][0])
                xs, xe = max(0, xs), min(cur[0]["clen"], xe)
                if nm / bl >= MIN_ID and aligned >= MIN_COV * min(Lq, xe - xs):
                    chains.append({"q": q, "chrom": chrom, "strand": strand, "recs": cur, "nm": nm, "s": ts, "e": te,
                                   "xs": xs, "xe": xe})
                cur = []
            if r is not None:
                cur.append(r)
    return chains


def project_exons(chain, qex):
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
                        out.append((t + (lo - qa), t + (hi - qa)) if r["strand"] == "+" else (t + (qb - hi), t + (qb - lo)))
                t += n
                q = q + n if r["strand"] == "+" else q - n
            elif op == "I":
                q = q + n if r["strand"] == "+" else q - n
            else:
                t += n
    return merge(out)


def tx_exon_blocks(h):
    blocks, pos, cur = [], h["ts"], h["ts"]
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
    return blocks


# ---------------------------------------------------------------- candidate construction
def blocked_by(seeds, rec, h):
    return any(rec[s]["chrom"] == h["chrom"] and ov(rec[s]["start0"], rec[s]["end"], h["s"], h["e"]) > 0 for s in seeds)


def single_linkage(hits, seeds, rec):
    """§6jm-§6jr construction (kept as the M0 / G1 references)."""
    free = [h for h in hits if h["q"] in seeds and not blocked_by(seeds, rec, h)]
    by = collections.defaultdict(list)
    for h in free:
        by[h["chrom"]].append(h)
    out = []
    for chrom, hs in by.items():
        hs.sort(key=lambda h: h["s"])
        cl, end = [], -1
        for h in hs + [None]:
            if h is None or (cl and h["s"] >= end):
                b = max(cl, key=lambda x: (x["nm"], -x["s"]))
                out.append({"chrom": chrom, "s": b["s"], "e": b["e"], "family": rec[b["q"]]["family"], "tx": None, "chain": None})
                cl, end = [], -1
            if h is not None:
                cl.append(h)
                end = max(end, h["e"])
    return out


def leaders(hits):
    """Leader clustering by reciprocal overlap >= RECIP, in decreasing nmatch."""
    lead = []
    for h in sorted(hits, key=lambda h: (-h["nm"], h["chrom"], h["s"])):
        if not any(L["chrom"] == h["chrom"] and recip(L["s"], L["e"], h["s"], h["e"]) for L in lead):
            lead.append(h)
    return lead


def union_fixed(tx_hits, chains, seeds, rec):
    """F1: transcript leaders and chain leaders, joined only as the same locus."""
    T = leaders([h for h in tx_hits if h["q"] in seeds and not blocked_by(seeds, rec, h)])
    C = leaders([c for c in chains if c["q"] in seeds and not blocked_by(seeds, rec, c)])
    parent = list(range(len(T) + len(C)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for j, c in enumerate(C):
        inside = [i for i, t in enumerate(T) if t["chrom"] == c["chrom"]
                  and ov(t["s"], t["e"], c["xs"], c["xe"]) >= CONTAIN * (t["e"] - t["s"])]
        for i, t in enumerate(T):
            if t["chrom"] != c["chrom"]:
                continue
            if recip(t["s"], t["e"], c["s"], c["e"]) or (i in inside and len(inside) == 1):
                parent[find(i)] = find(len(T) + j)
    comp = collections.defaultdict(list)
    for k in range(len(T) + len(C)):
        comp[find(k)].append(k)
    out = []
    for ks in comp.values():
        ts = [T[k] for k in ks if k < len(T)]
        cs = [C[k - len(T)] for k in ks if k >= len(T)]
        bt = max(ts, key=lambda h: (h["nm"], -h["s"])) if ts else None
        bc = max(cs, key=lambda h: (h["nm"], -h["s"])) if cs else None
        lead = bt or bc
        out.append({"chrom": lead["chrom"], "s": lead["s"], "e": lead["e"], "family": rec[lead["q"]]["family"],
                    "tx": bt, "chain": bc})
    return out


def classify(cands, hidden, rec, genes):
    for c in cands:
        hov = [(ov(rec[n]["start0"], rec[n]["end"], c["s"], c["e"]), n) for n in hidden if rec[n]["chrom"] == c["chrom"]]
        hov = [x for x in hov if x[0] > 0]
        if hov:
            c["match"], c["class"], c["named"], c["genes"] = max(hov)[1], "hidden", False, []
        else:
            g = [x for x in genes.get(c["chrom"], []) if ov(x[0], x[1], c["s"], c["e"]) > 0]
            c["match"], c["class"] = None, ("other_gene" if g else "unannotated")
            c["named"] = bool(FAMNAME.get(c["family"])) and any(FAMNAME[c["family"]].search(x[2] + " " + x[3]) for x in g)
            c["genes"] = sorted({x[2] for x in g})
    rp = {}
    for n in hidden:
        t = rec[n]
        o = [(ov(t["start0"], t["end"], c["s"], c["e"]), i) for i, c in enumerate(cands) if c["chrom"] == t["chrom"]]
        o = [x for x in o if x[0] > 0]
        rp[n] = cands[max(o)[1]] if o else None
    return rp


def width_stats(t, s, e):
    L = t["end"] - t["start0"]
    inter = ov(t["start0"], t["end"], s, e)
    left, right = t["start0"] - s, e - t["end"]
    return (inter / (max(t["end"], e) - min(t["start0"], s)), inter / L < 0.90, (max(0, left) + max(0, right)) / L > 0.10)


# ---------------------------------------------------------------- trees
def parse_newick(s):
    s = s.strip().rstrip(";")
    children, parent, label = {}, {}, {}
    stack, nid, i = [], [0], 0

    def new(p):
        n = nid[0]
        nid[0] += 1
        children[n], parent[n], label[n] = [], p, ""
        if p is not None:
            children[p].append(n)
        return n
    root = cur = new(None)
    while i < len(s):
        ch = s[i]
        if ch == "(":
            stack.append(cur)
            cur = new(cur)
            i += 1
        elif ch == ",":
            cur = new(stack[-1])
            i += 1
        elif ch == ")":
            cur = stack.pop()
            i += 1
        elif ch == ":":
            j = i + 1
            while j < len(s) and s[j] not in ",()":
                j += 1
            i = j
        else:
            j = i
            while j < len(s) and s[j] not in ",():":
                j += 1
            label[cur] = s[i:j]
            i = j
    memo = {}

    def leaves(n):
        if n not in memo:
            memo[n] = {label[n]} if not children[n] else set().union(*(leaves(c) for c in children[n]))
        return memo[n]
    allleaves = leaves(root)
    splits = []
    for n in children:
        if children[n] and parent[n] is not None:
            side = leaves(n)
            if 1 < len(side) < len(allleaves) - 1:
                sup = label[n].split("/")
                sh = float(sup[0]) if sup[0] else float("nan")
                splits.append((frozenset(side), sh, float(sup[1]) if len(sup) > 1 else float("nan")))
    return splits, allleaves


def projected_tree(tag, seqs, outdir, iqtree, threads):
    names = sorted(seqs)
    fa = f"{outdir}/{tag}.fa"
    with open(fa, "w") as fh:
        for n in names:
            fh.write(f">{n}\n{seqs[n]}\n")
    paf = subprocess.run(["minimap2", "-c", "-x", "asm20", "-X", "-N", "50", "-p", "0.1", "-t", str(threads), fa, fa],
                         capture_output=True, text=True, check=True).stdout
    partners, aligned = collections.defaultdict(set), collections.Counter()
    for line in paf.splitlines():
        f = line.split("\t")
        if f[0] == f[5]:
            continue
        partners[f[0]].add(f[5])
        partners[f[5]].add(f[0])
        aligned[f[0]] += int(f[3]) - int(f[2])
        aligned[f[5]] += int(f[8]) - int(f[7])
    ref = min(names, key=lambda n: (-len(partners[n]), -aligned[n], n))  # F3
    with open(f"{outdir}/{tag}.ref.fa", "w") as fh:
        fh.write(f">{ref}\n{seqs[ref]}\n")
    with open(f"{outdir}/{tag}.others.fa", "w") as fh:
        for n in names:
            if n != ref:
                fh.write(f">{n}\n{seqs[n]}\n")
    paf = subprocess.run(["minimap2", "-c", "-x", "asm20", "-N", "50", "-p", "0.1", "-t", str(threads),
                          f"{outdir}/{tag}.ref.fa", f"{outdir}/{tag}.others.fa"], capture_output=True, text=True, check=True).stdout
    by_q = collections.defaultdict(list)
    for line in paf.splitlines():
        f = line.split("\t")
        AS = next((int(x[5:]) for x in f[12:] if x.startswith("AS:i:")), 0)
        by_q[f[0]].append((AS, int(f[1]), int(f[2]), int(f[3]), f[4], int(f[7]),
                           next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")))
    Lr = len(seqs[ref])
    rows = {ref: list(seqs[ref])}
    for n in names:
        if n == ref:
            continue
        row, filled = ["-"] * Lr, bytearray(Lr)
        for AS, qlen, qs, qe, strand, ts, cg in sorted(by_q.get(n, []), key=lambda r: -r[0]):
            s = seqs[n] if strand == "+" else rc(seqs[n])
            q, t = (qs if strand == "+" else qlen - qe), ts
            for k, op in re.findall(r"(\d+)([MID])", cg):
                k = int(k)
                if op == "M":
                    for j in range(k):
                        if not filled[t + j]:
                            row[t + j], filled[t + j] = s[q + j], 1
                    t += k
                    q += k
                elif op == "D":
                    t += k
                else:
                    q += k
        rows[n] = row
    keep = [j for j in range(Lr) if sum(1 for n in names if rows[n][j] != "-") > 0.5 * len(names)]
    kept_rows = {n: "".join(rows[n][j] for j in keep) for n in names}
    dropped = [n for n, r in kept_rows.items() if not r.strip("-")]
    members = [n for n in names if n not in dropped]
    if len(members) < 4 or not keep:
        return None, ref, len(keep), dropped
    aln = f"{outdir}/{tag}.proj.fa"
    with open(aln, "w") as fh:
        for n in members:
            fh.write(f">{n}\n{kept_rows[n]}\n")
    subprocess.run([iqtree, "-s", aln, "-m", "MFP", "-B", "1000", "-alrt", "1000", "-T", str(threads), "--seed", "1",
                    "--prefix", f"{outdir}/{tag}", "-redo", "-quiet"], check=True, stdout=subprocess.DEVNULL,
                   stderr=subprocess.DEVNULL, timeout=600)
    return parse_newick(open(f"{outdir}/{tag}.treefile").read()), ref, len(keep), dropped


def literature_groups(truth):
    groups = {}
    for fam in sorted({t["family"] for t in truth}):
        T = [t for t in truth if t["family"] == fam]
        if fam == "NPIP":
            g = {"L1 NPIPA|NPIPB": {t["name"] for t in T if t["level1"] == "NPIPA"}}
            for v in ("A6-9", "B3-5", "B6-9", "B12/13"):
                g[f"L2 {v}"] = {t["name"] for t in T if t["level2"] == v}
            g["named B3,B4,B5,B11,B12,B13"] = {"NPIPB3", "NPIPB4", "NPIPB5", "NPIPB11", "NPIPB12", "NPIPB13"}
        elif fam == "TBC1D3":
            g = {v: {t["name"] for t in T if t["level2"] == v} for v in ("AE", "CDKL")}
        else:
            g = {}
            for lvl in ("level1", "level2"):
                for v in sorted({t[lvl] for t in T}):
                    m = {t["name"] for t in T if t[lvl] == v}
                    if len(m) >= 2 and f"L1 {v}" not in g:
                        g[f"{'L1' if lvl == 'level1' else 'L2'} {v}"] = m
        groups[fam] = g
    return groups


def clade_calls(splits, allleaves, label_of, groups, positional):
    truth_leaves = {m for m in allleaves if label_of.get(m)}
    present_all = {label_of[m] for m in truth_leaves}
    res = {}
    for gname, G in groups.items():
        present = present_all & G
        if len(present) < 2:
            res[gname] = "absent"
            continue
        call = "no split"
        for side, sh, _ in splits:
            ts = {label_of[m] for m in side if m in truth_leaves}
            tc = present_all - ts
            if ts == present or tc == present:
                call = "RECOVERED" if sh > 75 else "unsupported"
        res[gname] = call
    pos = None
    if positional:
        P = positional & present_all
        pos = any(sh > 75 and P and ({label_of[m] for m in side if m in truth_leaves} in (P, present_all - P))
                  for side, sh, _ in splits)
    return res, pos


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genes-gff-gz", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--mmi", required=True)
    ap.add_argument("--iqtree", required=True)
    ap.add_argument("--expected-units")
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    W = a.workdir
    genome = pysam.FastaFile(a.genome)
    truth = list(csv.DictReader(open(f"{W}/truth.tsv"), delimiter="\t"))
    for t in truth:
        t["start0"], t["end"] = int(t["start0"]), int(t["end"])
    rec = {t["name"]: t for t in truth}
    families = sorted({t["family"] for t in truth})
    ann = Annotation(a.gff, truth, genome)
    log = open(f"{W}/t.out", "w")

    def say(*x):
        print(*x, file=log, flush=True)

    if a.expected_units:
        exp, cur = {}, None
        for line in open(a.expected_units):
            if line.startswith(">"):
                cur = line[1:].strip()
                exp[cur] = ""
            else:
                exp[cur] += line.strip()
        bad = [n for n in rec if exp.get(n) != ann.unit_seq[n]]
        if bad:
            sys.exit(f"ABORT: seed units differ from {a.expected_units}: {bad}")
        say(f"# seed units identical to {a.expected_units} ({len(rec)} records)")
    with open(f"{W}/units.fa", "w") as fu, open(f"{W}/envelope.fa", "w") as fe:
        for n in sorted(rec):
            fu.write(f">{n}\n{ann.unit_seq[n]}\n")
            fe.write(f">{n}\n{ann.envelope_seq(n)}\n")
    for t in truth:
        s, e, kind = ann.envelope[t["name"]]
        say(f"# unit {t['name']}: {ann.unit_model[t['name']]}; gene-body query {kind} {e - s} bp (gene span {t['end'] - t['start0']} bp)")
    for fa, preset in (("units", "splice"), ("envelope", "asm20")):
        with open(f"{W}/{fa}.paf", "w") as fh:
            subprocess.run(["minimap2", "-c", "-x", preset, "-N", "100", "-p", "0.1", "-t", str(a.threads), a.mmi, f"{W}/{fa}.fa"],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
    tx_hits = transcript_hits(f"{W}/units.paf")
    chains = gene_body_chains(f"{W}/envelope.paf")
    say(f"# transcript hits passing {len(tx_hits)}; gene-body chains passing {len(chains)}")

    genes = collections.defaultdict(list)
    with gzip.open(a.genes_gff_gz, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
                continue
            n = re.search(r"(?:^|;)Name=([^;]+)", f[8])
            d = re.search(r"description=([^;]+)", f[8])
            genes[f[0]].append((int(f[3]) - 1, int(f[4]), n.group(1) if n else "?", d.group(1) if d else ""))

    groups = literature_groups(truth)
    positional = {"TBC1D3": {t["name"] for t in truth if t["family"] == "TBC1D3" and t["level1"] == "cluster1"}}
    tdir = f"{W}/tree_t"
    os.makedirs(tdir, exist_ok=True)
    rows, runs, cand_rows = [], [], []

    def exon_intron_record(n):
        t = rec[n]
        s0, e0, _ = ann.envelope[n]
        ex = [(max(s, s0), min(e, e0)) for s, e in ann.union_exons[n] if min(e, e0) > max(s, s0)]
        exon = "".join(genome.fetch(t["chrom"], s, e) for s, e in merge(ex)).upper()
        return (rc(exon) if t["strand"] == "-" else exon), ann.seq_minus(t["chrom"], s0, e0, ann.union_exons[n], t["strand"])

    def exon_intron_candidate(c):
        if c["tx"] is not None:
            b = tx_exon_blocks(c["tx"])
            exon = "".join(genome.fetch(c["chrom"], s, e) for s, e in b).upper()
            exon = rc(exon) if c["tx"]["strand"] == "-" else exon
        else:
            b = project_exons(c["chain"], ann.exons_in_envelope(c["chain"]["q"]))
            exon = "".join(genome.fetch(c["chrom"], s, e) for s, e in b).upper()
            exon = rc(exon) if c["chain"]["strand"] == "-" else exon
        if c["chain"] is not None:
            ch = c["chain"]
            intron = ann.seq_minus(c["chrom"], ch["xs"], ch["xe"], project_exons(ch, ann.exons_in_envelope(ch["q"])), ch["strand"])
        else:
            intron = ann.seq_minus(c["chrom"], c["tx"]["s"], c["tx"]["e"], tx_exon_blocks(c["tx"]), c["tx"]["strand"])
        return exon, intron

    for fam in families:
        ex, it = {}, {}
        for t in truth:
            if t["family"] == fam:
                ex[t["name"]], it[t["name"]] = exon_intron_record(t["name"])
        runs.append(("ref_" + fam, fam, {k: k for k in ex}, ex, it))

    for level in ("half", "keep1"):
        for rep in range(a.reps):
            seeds, hidden = set(), set()
            for fi, fam in enumerate(families):
                names = sorted(t["name"] for t in truth if t["family"] == fam)
                random.Random(1000 * rep + fi).shuffle(names)
                k = math.ceil(len(names) / 2) if level == "half" else 1
                seeds.update(names[:k])
                hidden.update(names[k:])
            arms = {"M0": single_linkage(tx_hits, seeds, rec),
                    "G1": single_linkage([dict(c, s=c["s"], e=c["e"]) for c in chains], seeds, rec),
                    "U": union_fixed(tx_hits, chains, seeds, rec)}
            for arm, cands in arms.items():
                rp = classify(cands, hidden, rec, genes)
                for fam in families + ["ALL"]:
                    H = [n for n in hidden if fam in ("ALL", rec[n]["family"])]
                    C = [c for c in cands if fam in ("ALL", c["family"])]
                    tp = sum(1 for n in H if rp[n] and rp[n]["family"] == rec[n]["family"])
                    good = sum(1 for c in C if c["class"] == "hidden" and rec[c["match"]]["family"] == c["family"])
                    named = sum(1 for c in C if c["class"] == "other_gene" and c["named"])
                    cross = sum(1 for c in C if c["class"] == "hidden" and rec[c["match"]]["family"] != c["family"])
                    unnamed = sorted({g for c in C if c["class"] == "other_gene" and not c["named"] for g in c["genes"]})
                    ip = [rp[n]["family"] if rp[n] else f"missed:{n}" for n in H] + [c["family"] for c in C if c["class"] != "hidden"]
                    itr = [rec[n]["family"] for n in H] + [f"nonmember:{i}" for i, c in enumerate(C) if c["class"] != "hidden"]
                    ps, pp = pairwise(ip, itr)
                    br, bp = bipartite(ip, itr) if ip else (float("nan"), float("nan"))
                    wd = [width_stats(rec[n], rp[n]["s"], rp[n]["e"]) for n in H if rp[n] and rp[n]["family"] == rec[n]["family"]]
                    rows.append({"level": level, "rep": rep, "arm": arm, "family": fam, "hidden": len(H), "cands": len(C),
                                 "sens": tp / len(H) if H else float("nan"), "prec": good / len(C) if C else float("nan"),
                                 "named": (good + named) / len(C) if C else float("nan"), "cross": cross,
                                 "unnamed": len(unnamed), "unnamed_genes": ",".join(unnamed),
                                 "pair_s": ps, "pair_p": pp, "bip_R": br, "bip_P": bp, "n_w": len(wd),
                                 "J": statistics.median([w[0] for w in wd]) if wd else float("nan"),
                                 "trunc": sum(w[1] for w in wd), "overext": sum(w[2] for w in wd)})
                if arm != "U":
                    continue
                for c in cands:
                    cand_rows.append((level, rep, c["chrom"], c["s"], c["e"], c["family"], c["class"], c["match"] or ",".join(c["genes"]),
                                      "tx" if c["tx"] else "chain"))
                for fam in families:
                    ex, it, lab = {}, {}, {}
                    for s in sorted(seeds):
                        if rec[s]["family"] == fam:
                            ex[s], it[s] = exon_intron_record(s)
                            lab[s] = s
                    for i, c in enumerate(cands):
                        if c["family"] == fam:
                            k = f"cand{i}"
                            ex[k], it[k] = exon_intron_candidate(c)
                            if c["class"] == "hidden" and rec[c["match"]]["family"] == fam:
                                lab[k] = c["match"]
                    runs.append((f"{level}_{rep}_{fam}", fam, lab, ex, it))

    with open(f"{W}/t.candidates.tsv", "w") as fh:
        fh.write("level\trep\tchrom\tstart\tend\tfamily\tclass\tmatch_or_genes\twidth_from\n")
        for r in cand_rows:
            fh.write("\t".join(map(str, r)) + "\n")

    say("\n## breadth and width (M0 transcript and G1 gene body: single-linkage as §6jr; U: F1 construction)")
    for level in ("half", "keep1"):
        for fam in families + ["ALL"]:
            for arm in ("M0", "G1", "U"):
                R = [r for r in rows if r["level"] == level and r["family"] == fam and r["arm"] == arm]
                g = lambda k: ms([r[k] for r in R])
                say(f"{level:5s} {fam:6s} {arm:2s} hid {g('hidden')} cand {g('cands')} | sens {g('sens')} prec {g('prec')} "
                    f"named {g('named')} | cross {g('cross')} unnamed {g('unnamed')} | pair {g('pair_s')}/{g('pair_p')} "
                    f"bip {g('bip_R')}/{g('bip_P')} | W n={g('n_w')} J {g('J')} trunc {g('trunc')} overext {g('overext')}")
            un = sorted({x for r in rows if r["level"] == level and r["family"] == fam and r["arm"] == "U" for x in r["unnamed_genes"].split(",") if x})
            if un:
                say(f"            U unnamed other genes: {un}")
    say("\n## B1 (U sensitivity >= max(M0, G1) - 0.02 at both levels, named precision >= 0.95, 0 cross-family)")
    for fam in families:
        verdict = []
        for level in ("half", "keep1"):
            m = {arm: statistics.mean([r["sens"] for r in rows if r["level"] == level and r["family"] == fam and r["arm"] == arm])
                 for arm in ("M0", "G1", "U")}
            pn = statistics.mean([r["named"] for r in rows if r["level"] == level and r["family"] == fam and r["arm"] == "U"])
            cr = sum(r["cross"] for r in rows if r["level"] == level and r["family"] == fam and r["arm"] == "U")
            ok = m["U"] >= max(m["M0"], m["G1"]) - 0.02 and pn >= 0.95 and cr == 0
            verdict.append(ok)
            say(f"  {fam} {level}: U {m['U']:.3f} vs M0 {m['M0']:.3f} / G1 {m['G1']:.3f}; named {pn:.3f}; cross {cr} -> {'ok' if ok else 'FAIL'}")
        say(f"  {fam}: B1 {'PASS' if all(verdict) else 'FAIL'}")

    say("\n## trees (reference-projected alignment, F3 reference; SH-aLRT > 75)")
    summary = collections.defaultdict(collections.Counter)
    for tag, fam, lab, ex, it in runs:
        calls_by_cls = {}
        for cls, seqs in (("exon", ex), ("intron", it)):
            seqs = {k: v for k, v in seqs.items() if len(v) >= 100}
            if len(seqs) < 4:
                say(f"{tag} {cls}: {len(seqs)} members -> not treed")
                continue
            tree, ref, kept, dropped = projected_tree(f"{tag}_{cls}", seqs, tdir, a.iqtree, a.threads)
            if tree is None:
                say(f"{tag} {cls}: reference {lab.get(ref, ref)}, {kept} columns, dropped {dropped} -> not treed")
                continue
            res, pos = clade_calls(*tree, {k: v for k, v in lab.items() if k in seqs and k not in dropped}, groups[fam], positional.get(fam))
            calls_by_cls[cls] = (res, pos)
            say(f"{tag} {cls}: {len(seqs) - len(dropped)} leaves, reference {lab.get(ref, ref)}, {kept} columns"
                f"{', dropped ' + str(dropped) if dropped else ''} | " + "; ".join(f"{g}: {v}" for g, v in res.items())
                + ("" if pos is None else f" | positional supported: {pos}"))
        kind = "reference" if tag.startswith("ref") else "leave-out"
        for g in groups[fam]:
            calls = {cls: v[0][g] for cls, v in calls_by_cls.items()}
            if not calls or all(c == "absent" for c in calls.values()):
                continue
            for cls, c in calls.items():
                summary[(fam, kind, g, cls)][c] += 1
            summary[(fam, kind, g, "either")]["RECOVERED" if "RECOVERED" in calls.values() else "not"] += 1
        for cls, v in calls_by_cls.items():
            if v[1] is not None:
                summary[(fam, kind, "positional split", cls)]["WRONG" if v[1] else "CORRECT"] += 1
    say("\n## tree summary (counts over runs)")
    for k in sorted(summary):
        say(f"{k[0]:6s} {k[1]:9s} {k[2]:28s} {k[3]:6s} {dict(summary[k])}")


if __name__ == "__main__":
    main()

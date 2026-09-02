#!/usr/bin/env python3
"""Score two arms of the frozen 150-window HUMAN negative panel and diff them.

Same two levels as `o1neg/score.py`, kept separate because they are different assertions:

  FAMILY level: every co-membership assertion surviving into `copies.tsv`. Each window holds
                exactly one locus, so ANY family emitted is a false merge. FALSE-by-coordinates
                if two copies of one family overlap genomically.
  EDGE   level: every E_r edge emitted. SELF-IDENTITY CERTIFIED if the two node spans overlap
                and the aligned block IS that genomic intersection.

Neither truth uses a gene group, an aligner, or an annotation label.

Usage: o1neg_score_arms.py <panel_dir> <armA> <armB>
       e.g. o1neg_score_arms.py /home/juanfra/winloci_scratch/o1neg off cl30
"""
import csv
import glob
import os
import sys
import json
import math
import itertools
import collections


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def score(panel_dir, arm):
    """Return the family- and edge-level tallies for one arm."""
    A = os.path.join(panel_dir, f"arm_{arm}")
    panel = json.load(open(os.path.join(panel_dir, "panel.json")))
    exits = os.path.join(A, "run_exit.tsv")
    done = [l.split("\t")[0] for l in open(exits)] if os.path.exists(exits) else []
    done = [d for d in done if d.startswith("W")]

    emitting, bad = [], []
    nfam = ncop = nass = 0
    for i, p in enumerate(panel):
        wid = f"W{i:03d}"
        f = os.path.join(A, "out", f"{wid}.copies.tsv")
        if not os.path.exists(f):
            continue
        by = collections.defaultdict(list)
        for r in csv.DictReader(open(f), delimiter="\t"):
            by[r["family_id"]].append((r["chrom"], int(r["start"]), int(r["end"])))
        if by:
            emitting.append((wid, p["sym"], {k: len(v) for k, v in by.items()}))
        for fid, v in by.items():
            nfam += 1
            ncop += len(v)
            for (c1, s1, e1), (c2, s2, e2) in itertools.combinations(sorted(v), 2):
                nass += 1
                if c1 == c2 and s1 < e2 and s2 < e1:
                    bad.append((wid, p["sym"], fid, f"{c1}:{s1}-{e1}", f"{c2}:{s2}-{e2}",
                                min(e1, e2) - max(s1, s2)))

    edges, ov_edges, certified = 0, 0, 0
    edge_rows, cert_rows = [], []
    for ef in sorted(glob.glob(os.path.join(A, "dump", "*.edges.tsv"))):
        pid = os.path.basename(ef).split(".")[0]
        nf = os.path.join(A, "dump", f"{pid}.nodes.tsv")
        replen = {}
        if os.path.exists(nf):
            for r in csv.DictReader(open(nf), delimiter="\t"):
                replen[r["idx"]] = int(r["exon_sum_len"])
        for r in csv.DictReader(open(ef), delimiter="\t"):
            edges += 1
            edge_rows.append((pid, r["rep_i"], r["rep_j"], float(r["identity"]),
                              float(r["coverage"]), r.get("cov_longer", "")))
            if r["chrom_i"] != r["chrom_j"]:
                continue
            s1, e1 = int(r["start_i"]), int(r["end_i"])
            s2, e2 = int(r["start_j"]), int(r["end_j"])
            ov = min(e1, e2) - max(s1, s2)
            if ov <= 0:
                continue
            ov_edges += 1
            li, lj = replen.get(r["rep_i"]), replen.get(r["rep_j"])
            if not li or not lj:
                continue
            if abs(round(float(r["coverage"]) * min(li, lj)) - ov) <= 2:
                certified += 1
                cert_rows.append((pid, r["node_key_i"], r["node_key_j"], ov))

    return dict(arm=arm, done=len(done), emitting=emitting, nfam=nfam, ncop=ncop, nass=nass,
                bad=bad, edges=edges, edge_rows=edge_rows, ov_edges=ov_edges,
                certified=certified, cert_rows=cert_rows)


def report(s):
    n = max(1, s["done"])
    lo, hi = wilson(len(s["emitting"]), n)
    print(f"--- arm {s['arm']}: {s['done']} windows completed")
    print(f"    windows emitting >=1 family: {len(s['emitting'])}/{s['done']} = "
          f"{len(s['emitting'])/n:.4f}  Wilson95 [{lo:.4f},{hi:.4f}]")
    for e in s["emitting"]:
        print(f"       {e}")
    print(f"    families={s['nfam']} copies={s['ncop']} co-membership assertions={s['nass']}"
          f"  self-overlap(false)={len(s['bad'])}")
    for b in s["bad"]:
        print(f"       {b}")
    print(f"    E_r edges={s['edges']}  overlapping-span={s['ov_edges']}  "
          f"self-identity CERTIFIED={s['certified']}")
    for r in s["edge_rows"]:
        print(f"       edge {r[0]} {r[1]}-{r[2]} id={r[3]:.6f} cov={r[4]:.6f} cov_longer={r[5]}")


def main():
    panel_dir, a, b = sys.argv[1], sys.argv[2], sys.argv[3]
    A, B = score(panel_dir, a), score(panel_dir, b)
    for s in (A, B):
        report(s)
        print()

    ea = {(r[0], r[1], r[2]) for r in A["edge_rows"]}
    eb = {(r[0], r[1], r[2]) for r in B["edge_rows"]}
    wa = {e[0] for e in A["emitting"]}
    wb = {e[0] for e in B["emitting"]}
    print("=== DIFF ===")
    print(f"edges {A['edges']} -> {B['edges']}   removed={sorted(ea - eb)}   ADDED={sorted(eb - ea)}")
    print(f"emitting windows {sorted(wa)} -> {sorted(wb)}")
    print(f"  lost={sorted(wa - wb)}   GAINED={sorted(wb - wa)}")
    json.dump(dict(a=A, b=B), open(os.path.join(panel_dir, f"score_{a}_vs_{b}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

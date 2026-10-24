#!/usr/bin/env python3
"""Score `candidate_copies` on a few-read world (docs/archive/2026-10/merged/PREREG_o3_fewread_flags_2026-10-14.md, Phase 1).

    score.py WORLD VARIANT [--flags OUT.tsv]

Per flagged candidate: its reads' true source copies (read names `<copy>|<iso>|<i>`), top_frac = the share of the most common
source, the union sequence's best true-transcript hit (minimap2 -c -x splice:hq -uf; identity x coverage as the stage computes it),
and the class (true / allele / M1 mixed / M2 pure-present). Prints one summary line per family type, the ERASED recall by E's
expression, and every M1/M2 flag. A flag is few-read iff n_reads <= 20.
"""
import argparse
import collections
import csv
import json
import subprocess
import sys

FEW = 20
TOP = 0.9


def classify(scenario, top_copy, top_frac, kinds):
    if top_frac < TOP:
        return "M1_mixed"
    k = kinds.get(top_copy, "?")
    if k == "present":
        return "M2_pure_present"
    if k == "copy" and scenario == "ERASED":
        return "true"
    if k == "allele" and scenario.startswith("ALLELE"):
        return "allele"
    return "other"


def best_hits(query_fa, target_fa):
    """best true-transcript hit per query by identity x coverage (matches/block x (qe-qs)/qlen)"""
    out = subprocess.run(f"minimap2 -c -x splice:hq -uf -N 20 {target_fa} {query_fa} 2>/dev/null", shell=True, check=True,
                         capture_output=True, text=True, executable="/bin/bash").stdout
    best = {}
    for ln in out.splitlines():
        r = ln.split("\t")
        q, ql, qs, qe, t, m, bl = r[0], int(r[1]), int(r[2]), int(r[3]), r[5], int(r[9]), int(r[10])
        ic = (m / bl) * (qe - qs) / ql if bl and ql else 0.0
        if ic > best.get(q, ("", -1.0))[1]:
            best[q] = (t, ic)
    return best


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("world"); ap.add_argument("variant"); ap.add_argument("--flags", help="write one row per flagged candidate")
    a = ap.parse_args(argv)
    w, v = a.world, a.variant
    P = f"{w}/{v}/P.cand"
    scen = json.load(open(f"{w}/world.json"))
    kinds = json.load(open(f"{w}/kinds.json"))
    cl2cand = {r["cluster"]: r["candidate"] for r in csv.DictReader(open(f"{P}.clusters.tsv"), delimiter="\t")}
    src = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(f"{P}.reads.tsv"), delimiter="\t"):
        c = cl2cand.get(r["cluster"])
        if c and c != "-":
            src[c][r["read"].split("|")[0]] += 1
    cands = [r for r in csv.DictReader(open(f"{P}.candidates.tsv"), delimiter="\t") if r["flagged"] == "1"]
    hits = best_hits(f"{P}.contigs.fa", f"{w}/transcripts.fa") if cands else {}

    rows = []
    for r in cands:
        fam, c = r["family"], r["candidate"]
        cnt = src[c]
        tot = sum(cnt.values())
        top_copy, top_n = cnt.most_common(1)[0] if cnt else ("-", 0)
        top_frac = top_n / tot if tot else 0.0
        hit = hits.get(c, ("-", 0.0))
        rows.append({"family": fam, "scenario": scen[fam], "candidate": c, "n_reads": int(r["n_reads"]), "n_clusters": int(r["n_clusters"]),
                     "d": r["d"], "union_len": r["union_len"], "top_copy": top_copy, "top_frac": round(top_frac, 3),
                     "sources": ",".join(f"{k}:{n}" for k, n in cnt.most_common()), "best_hit": hit[0], "best_id_cov": round(hit[1], 4),
                     "class": classify(scen[fam], top_copy, top_frac, kinds)})
    if a.flags:
        with open(a.flags, "w") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["family"], delimiter="\t")
            wr.writeheader(); wr.writerows(rows)

    types = ("NULL3", "NULLclose", "NULLlow", "ALLELE", "ERASED")
    fams = collections.defaultdict(set)
    for f, s in scen.items():
        fams[s.split("@")[0]].add(f)
    by = collections.defaultdict(lambda: collections.defaultdict(set))
    for x in rows:
        by[x["scenario"].split("@")[0]][x["class"]].add(x["family"])
        if x["class"] in ("M1_mixed", "M2_pure_present") and x["n_reads"] <= FEW:
            by[x["scenario"].split("@")[0]]["few_M1M2"].add(x["family"])
    for t in types:
        cls = by[t]
        print(f"{w.rstrip('/').split('/')[-1]}\t{v}\t{t}\tfamilies {len(fams[t])}\t" +
              "\t".join(f"{k} {len(cls[k])}" for k in ("true", "allele", "M1_mixed", "M2_pure_present", "other", "few_M1M2")))
    # ERASED recall by E's expression (from the world's reads)
    nE = collections.Counter()
    for ln in list(open(f"{w}/reads.truth.tsv"))[1:]:
        copy = ln.split("\t")[1]
        if kinds.get(copy) == "copy":
            nE[copy] += 1
    found = {x["top_copy"] for x in rows if x["class"] == "true"}
    lv = collections.defaultdict(lambda: [0, 0])
    for copy, n in nE.items():
        lv[n][0] += 1; lv[n][1] += copy in found
    print("  ERASED recall by E reads: " + "  ".join(f"{n}:{lv[n][1]}/{lv[n][0]}" for n in sorted(lv)))
    print("  true flags' n_reads: " + ",".join(str(x["n_reads"]) for x in sorted(rows, key=lambda x: x["n_reads"]) if x["class"] == "true"))
    for x in rows:
        if x["class"] not in ("true", "allele"):
            print("  FLAG\t" + "\t".join(f"{k}={x[k]}" for k in ("family", "scenario", "candidate", "class", "n_reads", "n_clusters", "d",
                                                               "top_frac", "sources", "best_hit", "best_id_cov")))


if __name__ == "__main__":
    sys.exit(main())

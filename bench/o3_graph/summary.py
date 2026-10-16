#!/usr/bin/env python3
"""Summaries the graph-O2/O3 verdicts are read from (round 4 review I4: these were scratch scripts before).

    summary.py o3 OUT.json WORLD ...      per scenario and arm: detections (REMOVED: of the DETECTABLE removed copies, median
                                          primary `de` on the world's genome > DELTA) and false-call families (score.py's rule)
    summary.py o2bars PREFIX WORLD ...    round 4's A1-A4: padded union O2 (o2.unionp.m2.tsv) against the current O2 and locus O2
                                          (o2.locus.m2.tsv), per world; writes PREFIX_<arm>.json (o2graph.py score)
    summary.py o2bars-unionf PREFIX WORLD ...   the same for round 5's graph O2 with the locus fallback (o2.unionf.m2.tsv)
"""
import collections
import csv
import importlib.util
import json
import os
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DELTA = 0.00958
ARMS = ("CUR", "S23", "S3", "S123m", "G", "H", "Hn", "K", "Kn")
NONCAT_WORLDS = ("OUTSIDE", "REMOVED", "BRANCH")


def detectable(des, delta=DELTA):
    """O3's definition: nothing in the genome explains the reads within delta (their median primary `de` > delta)."""
    return bool(des) and statistics.median(des) > delta


def o2_bars(worlds):
    """worlds: {world name: counters (cur_correct, graph_correct, graph_assigned, graph_wrong_family, locus_correct,
    cur_assigned_noncat, graph_assigned_noncat)} -> A1 (lowest correct ratio), A2 (lowest precision), A3 (fewer in REMOVED and
    BRANCH; pooled ratio over OUTSIDE, REMOVED, BRANCH), A4 (padded union / locus correct, pooled)."""
    a1 = min(v["graph_correct"] / max(1, v["cur_correct"]) for v in worlds.values())
    a2 = min(v["graph_correct"] / max(1, v["graph_assigned"] + v.get("graph_wrong_family", 0)) for v in worlds.values())
    nc = {n: (v.get("graph_assigned_noncat", 0), v.get("cur_assigned_noncat", 0)) for n, v in worlds.items()
          if n.split("_")[0] in NONCAT_WORLDS}
    fewer = all(g < c for n, (g, c) in nc.items() if n.split("_")[0] in ("REMOVED", "BRANCH"))
    pooled = sum(g for g, _ in nc.values()) / max(1, sum(c for _, c in nc.values()))
    a4 = sum(v["graph_correct"] for v in worlds.values()) / max(1, sum(v["locus_correct"] for v in worlds.values()))
    return dict(A1=round(a1, 4), A2=round(a2, 4), A3_fewer=fewer, A3_pooled=pooled, A4=a4)


def removed_detectability(w):
    """{family: (median de, detectable)} for a REMOVED world's absent copies (primary alignments only)."""
    import pysam
    truth = {r["read"]: r for r in csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t")}
    de = collections.defaultdict(list)
    with pysam.AlignmentFile(f"{w}/E0.bam") as bam:
        for rec in bam.fetch(until_eof=True):
            t = truth.get(rec.query_name)
            if t and t["kind"] == "absent" and not (rec.is_unmapped or rec.is_secondary or rec.is_supplementary):
                de[t["family"]].append(rec.get_tag("de"))
    return {f: (statistics.median(v), detectable(v)) for f, v in de.items()}


def cmd_o3(out, worlds):
    subprocess.run([sys.executable, os.path.join(HERE, "score.py"), out] + worlds, check=True, capture_output=True)
    rows = json.load(open(out))["rows"]
    det = {}
    for w in worlds:
        if os.path.basename(w.rstrip("/")) == "REMOVED":
            for f, (_, ok) in removed_detectability(w).items():
                det[(w, f)] = ok
    by = collections.defaultdict(collections.Counter)
    for r in rows:
        if r["arm"] not in ARMS:
            continue
        k = by[(r["scenario"], r["arm"])]
        k["families"] += 1
        k["false"] += r["false_call"]
        if r["scenario"] == "REMOVED":
            d = det.get((r["world"], r["family"]), False)
            k["detectable"] += d
            k["detected"] += bool(d and r["detected"])
        elif r["detected"] is not None:
            k["sources"] += 1
            k["detected"] += bool(r["detected"])
    for sc in sorted({s for s, _ in by}):
        cells = []
        for a in ARMS:
            k = by.get((sc, a))
            if not k:
                continue
            of = k["detectable"] if sc == "REMOVED" else k["sources"]
            cells.append(f"{a}: " + (f"det {k['detected']}/{of} " if of else "") + f"false {k['false']}/{k['families']}")
        print(sc, " | ".join(cells))
    return 0


def cmd_o2bars(prefix, worlds, graph_arm="unionp.m2"):
    res = {}
    for arm in (graph_arm, "locus.m2"):
        p = f"{prefix}_{arm}.json"
        subprocess.run([sys.executable, os.path.join(HERE, "o2graph.py"), "score", p, f"o2.{arm}.tsv"] + worlds, check=True,
                       capture_output=True)
        res[arm] = json.load(open(p))["worlds"]
    per = {}
    for w in worlds:
        u, lo = collections.Counter(), collections.Counter()
        for v in res[graph_arm][w].values():
            u.update(v)
        for v in res["locus.m2"][w].values():
            lo.update(v)
        u["locus_correct"] = lo["graph_correct"]
        per[os.path.basename(w.rstrip("/"))] = u
        print(f"{os.path.basename(w)}: current {u['cur_correct']} padded union {u['graph_correct']}/{u['graph_assigned']} "
              f"locus {lo['graph_correct']} | non-catalogued placed: current {u.get('cur_assigned_noncat', 0)} "
              f"padded union {u.get('graph_assigned_noncat', 0)}")
    print(o2_bars(per))
    return 0


def main(argv):
    if argv[0] == "o3":
        return cmd_o3(argv[1], argv[2:])
    if argv[0] == "o2bars-unionf":                       # round 5: the graph O2 with the locus fallback
        return cmd_o2bars(argv[1], argv[2:], graph_arm="unionf.m2")
    return cmd_o2bars(argv[1], argv[2:])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

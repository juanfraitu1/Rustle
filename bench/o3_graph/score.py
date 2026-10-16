#!/usr/bin/env python3
"""Score graph-O3 worlds (prereg docs/PREREG_o3_graph_sim_2026-10-15.md).

    score.py OUT.json WORLD [WORLD ...]

Per family and arm (CUR = `copy_assign --psv-clusters`: its `unexplained` clusters with est_div > DELTA; R0 / S1 / S12 / S123 /
S3 / S23 / S4n / S4 / S123m / Gn / G / Hn / H = detect.py's calls): the arm's calls are member sets.
  source       the family's copy absent from the genome (kind `absent`); an allele (`allele`), an uncatalogued paralog
               (`outside`) and a twin are no source: every call in such a family is a false call
  detected     some call holds >= RECALL of the source's PLACED reads (those with a copy_assign row of the family) at purity
               (source reads / call reads) >= PURITY  (register 1294's recovery rule)
  false call   some call whose reads are < half the source's
"""
import collections
import csv
import json
import sys

RECALL, PURITY, DELTA = 0.8, 0.9, 0.00958
ARMS = ("CUR", "R0", "S1", "S12", "S123", "S3", "S23", "S4n", "S4", "S123m", "Gn", "G", "Hn", "H")


def source_reads(truth_rows, copy):
    """The absent copy's reads a method can tell apart: for a novel-exon copy (truth column `novel`), those carrying the novel
    exon -- its other reads are the present copy's sequence within the SNP rate; otherwise every read."""
    rows = [r for r in truth_rows if r["copy"] == copy]
    if any(r.get("novel") == "1" for r in rows):          # only a copy that HAS a novel exon is scored on its carriers
        return {r["read"] for r in rows if r["novel"] == "1"}
    return {r["read"] for r in rows}


def family_outcome(calls, source, placed_source):
    """-> (detected: bool, or None without a source; false call: bool)."""
    det = None
    if source:
        det = any(len(c & placed_source) >= RECALL * len(placed_source) and len(c & source) >= PURITY * len(c) for c in calls)
    false = any(len(c & source) < 0.5 * len(c) for c in calls)
    return det, false


def world_calls(w):
    """{arm: {family: [member sets]}}"""
    out = {a: collections.defaultdict(list) for a in ARMS}
    mem = collections.defaultdict(set)
    for r in csv.DictReader(open(f"{w}/cur/A.psv_cluster_reads.tsv"), delimiter="\t"):
        mem[(r["family"], r["host"], r["cluster"])].add(r["read_name"])
    for c in csv.DictReader(open(f"{w}/cur/A.psv_clusters.tsv"), delimiter="\t"):
        if c["status"] == "unexplained" and float(c["est_div"]) > DELTA:
            out["CUR"][c["family"]].append(mem[(c["family"], c["host"], c["cluster"])])
    for r in csv.DictReader(open(f"{w}/graph/calls.tsv"), delimiter="\t"):
        out[r["arm"]][r["family"]].append(set(r["members"].split(",")))
    return out


def score_world(w):
    kinds = json.load(open(f"{w}/kinds.json"))
    scen = json.load(open(f"{w}/world.json"))
    placed = collections.defaultdict(set)
    for r in csv.DictReader(open(f"{w}/cur/A.assignments.tsv"), delimiter="\t"):
        placed[r["family_id"]].add(r["read_name"])
    truth_rows = list(csv.DictReader(open(f"{w}/reads.truth.E0.tsv"), delimiter="\t"))     # the truth table, never read names
    reads_of = collections.defaultdict(set)
    for r in truth_rows:
        reads_of[r["copy"]].add(r["read"])
    import os
    fam_of = json.load(open(f"{w}/copy_family.json")) if os.path.exists(f"{w}/copy_family.json") else {c: c[:-1] for c in kinds}
    calls = world_calls(w)
    rows = []
    for fam in sorted(scen):
        src_copy = [c for c, k in kinds.items() if k == "absent" and fam_of[c] == fam]
        source = source_reads(truth_rows, src_copy[0]) if src_copy else set()
        every = reads_of[src_copy[0]] if src_copy else set()
        for a in ARMS:
            det, false = family_outcome(calls[a].get(fam, []), source, source & placed[fam])
            det_all, _ = family_outcome(calls[a].get(fam, []), every, every & placed[fam])
            rows.append(dict(world=w, family=fam, scenario=scen[fam], arm=a, detected=det, detected_all_reads=det_all,
                             false_call=false, n_calls=len(calls[a].get(fam, [])), n_source=len(source),
                             n_source_placed=len(source & placed[fam]), n_source_all=len(every)))
    return rows


def summary(rows):
    out = collections.defaultdict(dict)
    by = collections.defaultdict(list)
    for r in rows:
        by[(r["scenario"], r["arm"])].append(r)
    for (sc, a), rs in sorted(by.items()):
        d = [r["detected"] for r in rs if r["detected"] is not None]
        out[sc][a] = dict(families=len(rs), detected=sum(d) if d else None, false_call_families=sum(r["false_call"] for r in rs))
    return out


def main(argv):
    out_p, worlds = argv[0], argv[1:]
    rows = [r for w in worlds for r in score_world(w)]
    s = summary(rows)
    json.dump(dict(rows=rows, summary=s), open(out_p, "w"), indent=1)
    for sc, arms in s.items():
        print(sc, " | ".join(f"{a}: det {v['detected']} false {v['false_call_families']}/{v['families']}" for a, v in arms.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

#!/usr/bin/env python3
"""Score the real-data O3 test (prereg docs/PREREG_o3_realdata_2026-10-15.md) from `copy_assign --psv-clusters` products.

    score.py removed   LABELS PREFIX OUT.tsv     per masked copy: recovered, recall, purity, est_div, label, fates of its reads
    score.py addedback LABELS PREFIX OUT.tsv     per masked copy: cleared once the copy is back in the genome
    score.py control   GENES PREFIX OUT.tsv      per single-copy gene: copy / allele_or_copy / - (rd_lib.control_flags)
    score.py verdict   --removed R.tsv --addedback A.tsv [--control C.tsv] [--min-placed 6] [--picked P.tsv]   R1-R3

PREFIX is copy_assign's --out (PREFIX.assignments.tsv, PREFIX.psv_clusters.tsv, PREFIX.psv_cluster_reads.tsv). A read is placed
on a family when it has a row of that family in the assignments (the reads O2 places, every status).
"""
import argparse
import collections
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rd_lib as L  # noqa: E402

R1_BAR, R2_BAR, R3_BAR, MIN_PLACED = 0.70, 0.90, 0.02, 6


def load(prefix):
    placed = collections.defaultdict(set)
    with open(prefix + ".assignments.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            placed[r["family_id"]].add(r["read_name"])
    with open(prefix + ".psv_clusters.tsv") as fh:
        clusters = list(csv.DictReader(fh, delimiter="\t"))
    members = collections.defaultdict(set)
    with open(prefix + ".psv_cluster_reads.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            members[(r["family"], r["host"], r["cluster"])].add(r["read_name"])
    return placed, clusters, members


def absent_reads(labels_path):
    """{family: (masked copy, set of its reads)} from labels.tsv (role D)."""
    out = {}
    with open(labels_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["role"] == "D":
                out.setdefault(r["family"], (r["copy"], set()))[1].add(r["read"])
    return out


def write_rows(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def cmd_removed(a):
    placed, clusters, members = load(a.prefix)
    anywhere = set().union(*placed.values()) if placed else set()
    rows = []
    for fam, (copy, d) in sorted(absent_reads(a.labels).items()):
        s = L.score_absent(fam, d, placed.get(fam, set()), clusters, members)
        elsewhere = {r for r in d - placed.get(fam, set()) if r in anywhere}
        rows.append(dict(family=fam, masked=copy, **s, on_other_family=len(elsewhere),
                         not_placed=len(d) - s["n_placed"] - len(elsewhere),
                         unexplained_clusters=sum(1 for c in clusters if c["family"] == fam and c["status"] == "unexplained")))
    write_rows(a.out, rows)
    print(json.dumps(dict(families=len(rows), recovered=sum(r["recovered"] for r in rows))))


def cmd_addedback(a):
    placed, clusters, members = load(a.prefix)
    rows = []
    for fam, (copy, d) in sorted(absent_reads(a.labels).items()):
        p = placed.get(fam, set())
        rows.append(dict(family=fam, masked=copy, n_placed=len(d & p), cleared=int(L.cleared(fam, d, p, clusters, members)),
                         unexplained_clusters=sum(1 for c in clusters if c["family"] == fam and c["status"] == "unexplained")))
    write_rows(a.out, rows)
    print(json.dumps(dict(families=len(rows), cleared=sum(r["cleared"] for r in rows))))


def cmd_control(a):
    _, clusters, _ = load(a.prefix)
    with open(a.genes) as fh:
        genes = list(csv.DictReader(fh, delimiter="\t"))
    flags = L.control_flags([g["family"] for g in genes], clusters)
    rows = [dict(g, flag=flags[g["family"]]) for g in genes]
    write_rows(a.out, rows)
    print(json.dumps(dict(genes=len(rows), copy=sum(r["flag"] == "copy" for r in rows))))


def verdict(removed, addedback, control, min_placed=MIN_PLACED, copies_of=None):
    elig = [r for r in removed if int(r["n_placed"]) >= min_placed]
    rec = [r for r in elig if int(r["recovered"])]
    ab = {r["family"]: r for r in addedback}
    clr = [r for r in rec if r["family"] in ab and int(ab[r["family"]]["cleared"])]
    out = dict(eligible=len(elig), recovered=len(rec), r1=len(rec) / len(elig) if elig else None,
               cleared=len(clr), r2=len(clr) / len(rec) if rec else None,
               labels=dict(collections.Counter(r["label"] for r in rec)))
    if copies_of:
        strata = collections.defaultdict(lambda: dict(eligible=0, recovered=0, cleared=0))
        for r in elig:
            k = "2" if copies_of[r["family"]] == 2 else ">=3"
            strata[k]["eligible"] += 1
            strata[k]["recovered"] += int(r["recovered"])
            strata[k]["cleared"] += int(r in clr)
        out["by_copies"] = {k: dict(v) for k, v in sorted(strata.items())}
    out["r1_pass"] = out["r1"] is not None and out["r1"] >= R1_BAR
    out["r2_pass"] = out["r2"] is not None and out["r2"] >= R2_BAR
    if control is not None:
        rnd = [r for r in control if r["stratum"] == "random"]
        fc = [r for r in rnd if r["flag"] == "copy"]
        out.update(control_genes=len(rnd), false_copy=len(fc), r3=len(fc) / len(rnd) if rnd else None)
        out["r3_pass"] = out["r3"] is not None and out["r3"] <= R3_BAR
        by = collections.defaultdict(collections.Counter)
        for r in control:
            by[(r["stratum"], r["bin"])][r["flag"]] += 1
        out["control_by_bin"] = {f"{s} {b}": dict(c) for (s, b), c in sorted(by.items())}
    return out


def cmd_verdict(a):
    rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))  # noqa: E731
    copies_of = {r["family"]: 1 + len(r["survivors"].split(",")) for r in rd(a.picked)} if a.picked else None
    res = verdict(rd(a.removed), rd(a.addedback), rd(a.control) if a.control else None, a.min_placed, copies_of)
    print(json.dumps(res, indent=1, sort_keys=True))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("removed", "addedback"):
        p = sub.add_parser(name)
        p.add_argument("labels")
        p.add_argument("prefix")
        p.add_argument("out")
    p = sub.add_parser("control")
    p.add_argument("genes")
    p.add_argument("prefix")
    p.add_argument("out")
    p = sub.add_parser("verdict")
    p.add_argument("--removed", required=True)
    p.add_argument("--addedback", required=True)
    p.add_argument("--control")
    p.add_argument("--min-placed", type=int, default=MIN_PLACED)
    p.add_argument("--picked", help="picked.tsv: results also by copy number (2 vs >= 3)")
    a = ap.parse_args(argv)
    {"removed": cmd_removed, "addedback": cmd_addedback, "control": cmd_control, "verdict": cmd_verdict}[a.cmd](a)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Post hoc, descriptive (not pre-registered): why the real-data O3 test misses (docs/archive/2026-10/merged/PREREG_o3_realdata_2026-10-15.md, Results).

    posthoc.py LABELS REMOVED_PREFIX ADDEDBACK_PREFIX OUT.json [--min-placed 6]

removed arm, per eligible masked copy: its placed reads in PURE unexplained clusters (purity >= 0.9) taken together (the
            pre-registered rule takes the single best cluster), the number of such clusters, and where the rest sit
            (explained / unlinked / unexplained-impure / in no cluster).
added back, everything present: unexplained clusters by label (`copy` = est_div > delta), the families with a `copy`
            cluster, and whether those clusters are mostly the masked copy's reads or a survivor's.
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
import score as S  # noqa: E402


def read_roles(labels):
    role = {}
    with open(labels) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            role[r["read"]] = (r["family"], r["role"])
    return role


def removed_detail(absent, placed, clusters, members, min_placed):
    rows = []
    for fam, (copy, d) in sorted(absent.items()):
        ap = d & placed.get(fam, set())
        if len(ap) < min_placed:
            continue
        where = collections.Counter()
        pure_hit, n_pure = set(), 0
        seen = set()
        for c in clusters:
            if c["family"] != fam:
                continue
            mem = members.get((c["family"], c["host"], c["cluster"]), set())
            got = mem & ap
            if not got:
                continue
            seen |= got
            purity = len(got) / int(c["n_reads"])
            if c["status"] == "unexplained" and purity >= L.PURITY:
                pure_hit |= got
                n_pure += 1
                where["unexplained_pure"] += len(got)
            else:
                where[c["status"] + ("_impure" if c["status"] == "unexplained" else "")] += len(got)
        where["no_cluster"] = len(ap - seen)
        rows.append(dict(family=fam, n_placed=len(ap), combined_recall=round(len(pure_hit) / len(ap), 3),
                         pure_clusters=n_pure, where=dict(where)))
    return rows


def addedback_detail(absent, role, clusters, members):
    out = collections.Counter()
    fams_copy = set()
    for c in clusters:
        if c["status"] != "unexplained":
            continue
        lab = "copy" if float(c["est_div"]) > L.DELTA else "allele_or_copy"
        out[lab] += 1
        if lab == "copy":
            fams_copy.add(c["family"])
            mem = members.get((c["family"], c["host"], c["cluster"]), set())
            rs = collections.Counter(role.get(r, ("?", "?"))[1] for r in mem)
            out["copy_mostly_" + ("D" if rs["D"] > rs["S"] else "S")] += 1
    return dict(unexplained=dict(out), families=len({f for f in absent}), families_with_copy=len(fams_copy))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("labels")
    ap.add_argument("removed")
    ap.add_argument("addedback")
    ap.add_argument("out")
    ap.add_argument("--min-placed", type=int, default=S.MIN_PLACED)
    a = ap.parse_args(argv)
    absent, role = S.absent_reads(a.labels), read_roles(a.labels)
    placed, clusters, members = S.load(a.removed)
    rows = removed_detail(absent, placed, clusters, members, a.min_placed)
    _, cl_ab, mem_ab = S.load(a.addedback)
    res = dict(eligible=len(rows), combined_recall_ge_0_8=sum(r["combined_recall"] >= L.RECALL for r in rows),
               rows=rows, addedback=addedback_detail(absent, role, cl_ab, mem_ab))
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print(json.dumps({k: res[k] for k in ("eligible", "combined_recall_ge_0_8", "addedback")}, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())

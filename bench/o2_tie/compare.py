#!/usr/bin/env python3
"""F2 of docs/PREREG_o2_exact_tie_2026-10-13.md: row-for-row comparison of two copy_assign runs.

    compare.py BEFORE_PREFIX AFTER_PREFIX OUT.json
Every row whose status changes must have had |old margin| < 5e-4 and move from `assigned` to `ambiguous`/`tied`.
"""
import collections, json, sys


def rows(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        out = {}
        for l in fh:
            r = dict(zip(head, l.rstrip("\n").split("\t")))
            out[(r["read_name"], r["family_id"])] = r
        return head, out


def classify(b, a):
    """-> None (identical) or a change kind."""
    if a == b:
        return None
    if a["status"] != b["status"]:
        ok = abs(float(b["margin"])) < 5e-4 and b["status"] == "assigned" and a["status"] in ("ambiguous", "tied")
        return "status_tie_fix" if ok else "status_OTHER"
    if a["assigned_copy"] != b["assigned_copy"]:
        return "copy_only"
    diff = sorted(k for k in a if a[k] != b.get(k))
    return "fields:" + ",".join(diff)


def main(before, after, out):
    hb, B = rows(before + ".assignments.tsv")
    ha, A = rows(after + ".assignments.tsv")
    c = collections.Counter()
    ex = collections.defaultdict(list)
    c["rows_before"], c["rows_after"] = len(B), len(A)
    c["missing_after"] = len(set(B) - set(A))
    c["new_after"] = len(set(A) - set(B))
    for k in set(B) & set(A):
        kind = classify(B[k], A[k])
        if kind is None:
            c["identical"] += 1
            continue
        c[kind] += 1
        if len(ex[kind]) < 10:
            ex[kind].append({"read": k[0], "family": k[1], "before": {f: B[k][f] for f in ("status", "assigned_copy", "margin", "p_value")},
                             "after": {f: A[k][f] for f in ("status", "assigned_copy", "margin", "p_value")}})
    files = {}
    for ext in ("quant", "families"):
        try:
            files[ext] = open(f"{before}.{ext}.tsv").read() == open(f"{after}.{ext}.tsv").read()
        except FileNotFoundError:
            files[ext] = None
    bad = c["status_OTHER"] + c["missing_after"] + c["new_after"] + sum(v for k, v in c.items() if k.startswith("fields:"))
    res = {"counts": dict(c), "identical_files": files, "F2_verdict": "PASS" if bad == 0 else "FAIL", "examples": ex}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k != "examples"}, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:4])

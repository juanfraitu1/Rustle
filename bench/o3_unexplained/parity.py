#!/usr/bin/env python3
"""Parity of `copy_assign --psv-clusters` (Rust) with the frozen prototype psvclust.py, cluster for cluster.

    parity.py PY_DIR RUST_PREFIX     # PY_DIR/psvclust.tsv vs RUST_PREFIX.psv_clusters.tsv + .psv_cluster_reads.tsv

Per (family, host, cluster index): n_reads, the members' source counts (read-name prefix), n_sites, n_alt, n_unexplained_alt, status,
explained_by, and est_div / linkage / editing_frac within rounding. Exit 1 on any difference.
"""
import collections
import csv
import sys


def main():
    py_dir, rp = sys.argv[1], sys.argv[2]
    py = {(r["family"], r["host"], r["cluster"]): r for r in csv.DictReader(open(f"{py_dir}/psvclust.tsv"), delimiter="\t")}
    rs = {(r["family"], r["host"], r["cluster"]): r for r in csv.DictReader(open(f"{rp}.psv_clusters.tsv"), delimiter="\t")}
    src = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(f"{rp}.psv_cluster_reads.tsv"), delimiter="\t"):
        src[(r["family"], r["host"], r["cluster"])][r["read_name"].split("|")[0]] += 1
    diffs = []
    for k in sorted(set(py) | set(rs)):
        a, b = py.get(k), rs.get(k)
        if a is None or b is None:
            diffs.append(f"{k}: only in {'rust' if a is None else 'python'}")
            continue
        pys = dict((x.rsplit(":", 1)[0], int(x.rsplit(":", 1)[1])) for x in a["sources"].split(","))
        for f in ("n_reads", "n_sites", "n_alt", "n_unexplained_alt", "status", "explained_by"):
            if a[f] != b[f]:
                diffs.append(f"{k}: {f} python {a[f]} rust {b[f]}")
        if pys != dict(src[k]):
            diffs.append(f"{k}: sources python {pys} rust {dict(src[k])}")
        for f, tol in (("est_div", 1e-5), ("editing_frac", 1e-3)):
            if abs(float(a[f]) - float(b[f])) > tol:
                diffs.append(f"{k}: {f} python {a[f]} rust {b[f]}")
        if (a["linkage"] == "NA") != (b["linkage"] == "NA") or (a["linkage"] != "NA" and abs(float(a["linkage"]) - float(b["linkage"])) > 1e-3):
            diffs.append(f"{k}: linkage python {a['linkage']} rust {b['linkage']}")
    print(f"{py_dir}\tclusters python {len(py)} rust {len(rs)}\tdifferences {len(diffs)}")
    for d in diffs[:20]:
        print("  " + d)
    return 1 if diffs else 0


if __name__ == "__main__":
    sys.exit(main())

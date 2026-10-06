#!/usr/bin/env python3
"""Own-node flags of the NPIP copies per arm, the input of `bench/copy_support.py --nodes` (the rule of bench/npip_read_pool/pagedata.py).

    nodes.py --copies copies.tsv --truth truth.gtf --family NPIP --arm NAME=loci.gff3,clusters.tsv [--arm ...] --out nodes.json

An arm's NPIP clusters are the clusters (clusters.tsv: one row per member locus, keyed by the span of the locus's gene row in loci.gff3)
that hold a locus overlapping >= 1 exonic bp of a same-strand NPIP copy (the copy's exon union over all its transcripts in --truth). A
copy has an OWN NODE in the arm iff such a locus of an NPIP cluster overlaps it. Two loci on one span, or a cluster row on no locus,
stop the run (the join would otherwise be silent). Writes {"rows": [{"cid", "name", "node": {ARM: bool}}]}.
"""
import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import copy_support as cs  # noqa: E402


def gffx(path):
    """locus name -> dict(c, s1 (1-based gene start), e, st, ex [(s0, e)]) from a loci GFF3 (gene rows + exon rows keyed by `gene=`)."""
    out = {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "gene":
            out[at["Name"]] = dict(c=f[0], s1=int(f[3]), e=int(f[4]), st=f[6], ex=[])
        elif f[2] == "exon":
            out[at["gene"]]["ex"].append((int(f[3]) - 1, int(f[4])))
    return out


def exov(a, b):
    return sum(max(0, min(y, q) - max(x, p)) for x, y in a for p, q in b)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--arm", action="append", required=True, help="NAME=loci.gff3,clusters.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--exons-json", help="npip_read_pool.json: the copy exons as the registered own-node rule took them (copies[].exons, the CAT gene's exon blocks); "
                    "default: the union of the copy's transcripts in --truth (the two differ at PKD1P6-NPIPP1)")
    ap.add_argument("--check", help="a pagedata.json: compare every arm it names (node[arm] per cid) with this run's flags")
    ap.add_argument("--check-out", help="write {compared, mismatches, detail} of --check here")
    a = ap.parse_args()
    copies = [r for r in csv.DictReader(open(a.copies), delimiter="\t") if r["family"] == a.family]
    genes = {r["isoform_gene"] for r in copies} | {r["cid"] for r in copies}
    tx = cs.gtf_transcripts(a.truth, genes)
    cu = {}
    for c in copies:
        texons = tx.get(c["cid"]) or tx.get(c["isoform_gene"], {})
        cu[c["cid"]] = cs.merge([e for ex in texons.values() for e in ex]) or [[int(c["terr_lo0"]), int(c["terr_hi"])]]
    if a.exons_json:
        ej = {c["cid"]: cs.merge([list(e) for e in c["exons"]]) for c in json.load(open(a.exons_json))["copies"]}
        cu = {cid: ej[cid] for cid in cu}
    node = {c["cid"]: {} for c in copies}
    for spec in a.arm:
        name, rest = spec.split("=", 1)
        loci_path, clusters_path = rest.split(",", 1)
        L = gffx(loci_path)
        by = {}
        for n, v in L.items():
            k = (v["c"], v["s1"], v["e"])
            if k in by:
                sys.exit(f"{loci_path}: loci {by[k]} and {n} share the span {k}")
            by[k] = n
        cl = {}
        for r in csv.DictReader(open(clusters_path), delimiter="\t"):
            k = (r["chrom"], int(r["start"]), int(r["end"]))
            if k not in by:
                sys.exit(f"{clusters_path}: the cluster row {k} is on no locus of {loci_path}")
            cl[by[k]] = r["cluster_id"]
        on = {n: [c["cid"] for c in copies if c["strand"] == v["st"] and c["chrom"] == v["c"] and exov(v["ex"], cu[c["cid"]]) > 0]
              for n, v in L.items()}
        fam = {cl[n] for n in L if n in cl and on[n]}
        own = {cid for n in L if cl.get(n) in fam for cid in on[n]}
        for c in copies:
            node[c["cid"]][name] = c["cid"] in own
        print(f"{name}: {len(fam)} NPIP clusters, {len(own)} of {len(copies)} copies with an own node")
    rows = [dict(cid=c["cid"], name=c.get("refseq_name") or c.get("cat_name") or c["name"], node=node[c["cid"]]) for c in copies]
    json.dump(dict(rows=rows), open(a.out, "w"), indent=1)
    if a.check:
        stored = {r["cid"]: r.get("node", {}) for r in json.load(open(a.check))["rows"]}
        compared, detail = 0, []
        for r in rows:
            for arm, v in r["node"].items():
                if arm in stored.get(r["cid"], {}):
                    compared += 1
                    if bool(stored[r["cid"]][arm]) != v:
                        detail.append(f"{r['name']} {arm}: stored {stored[r['cid']][arm]} here {v}")
        rep = dict(compared=compared, mismatches=len(detail), detail=detail)
        print(f"check vs {a.check}: {compared} flags compared, {len(detail)} mismatches")
        if a.check_out:
            json.dump(rep, open(a.check_out, "w"), indent=1)


if __name__ == "__main__":
    main()

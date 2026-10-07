#!/usr/bin/env python3
"""Own-node flags and the composition of the family's clusters, per arm (docs/PREREG_seed_pool_real_reads_2026-10-07.md section 3).

    composition.py --copies copies.tsv --truth truth.gtf --family NPIP --arm NAME=loci.gff3,clusters.tsv [--arm ...]
                   --out nodes.json [--report comp.json] [--exons-json npip_read_pool.json]

The join of bench/default_rescore/nodes.py, with one difference: a cluster file names a locus by its span (chrom, start, end) only, and
nodes.py stops when two loci share a span. Here every locus on a shared span takes the cluster ids of that span (a locus is in a cluster
if any row of its span is), and the number of shared spans is reported. Where nodes.py is defined the flags are equal (test_composition.py).

FAMILY CLUSTERS of an arm = the clusters that hold a locus overlapping >= 1 exonic bp of a same-strand truth copy; their loci are the NODES.
A copy has an OWN NODE iff a node overlaps it. K* = the family cluster holding the most distinct copies; M6 = the copies in K*. A node is ON_COPY (that overlap), else IN_SPAN (overlaps the span of a same-strand copy, off its
exons), else ANTISENSE (overlaps the span of an opposite-strand copy), else ELSEWHERE. NP = on-copy nodes / nodes; nodes that are real but
unlabelled count as false. Writes {"rows": [{"cid", "name", "node": {ARM: bool}}]} (the input of copy_support.py --nodes) and the report.
"""
import argparse
import collections
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import copy_support as cs  # noqa: E402


def read_loci(path):
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


def read_rows(path):
    """clusters.tsv -> [(chrom, start, end, cluster_id)]."""
    return [(r["chrom"], int(r["start"]), int(r["end"]), r["cluster_id"]) for r in csv.DictReader(open(path), delimiter="\t")]


def cluster_ids(loci, rows):
    """({locus: set of cluster ids of its span}, number of spans shared by >= 2 loci). A row on no locus stops the run."""
    by = collections.defaultdict(list)
    for n, v in loci.items():
        by[(v["c"], v["s1"], v["e"])].append(n)
    ids = {n: set() for n in loci}
    for ch, s, e, cid in rows:
        if (ch, s, e) not in by:
            sys.exit(f"the cluster row {(ch, s, e)} is on no locus")
        for n in by[(ch, s, e)]:
            ids[n].add(cid)
    return ids, sum(1 for names in by.values() if len(names) > 1)


def exov(a, b):
    return sum(max(0, min(y, q) - max(x, p)) for x, y in a for p, q in b)


def analyse(loci, rows, copies):
    """copies: cid -> dict(chrom, strand, exons [[s0, e]], span (lo, hi)). Returns the arm's report (lists sorted, so the JSON is stable)."""
    ids, shared = cluster_ids(loci, rows)
    on = {n: [cid for cid, c in copies.items() if c["strand"] == v["st"] and c["chrom"] == v["c"] and exov(v["ex"], c["exons"]) > 0]
          for n, v in loci.items()}
    fam = {i for n in loci if on[n] for i in ids[n]}
    nodes = [n for n in loci if ids[n] & fam]
    own = sorted({cid for n in nodes for cid in on[n]})

    def klass(n):
        v = loci[n]
        if on[n]:
            return "on_copy"
        for want_same, label in ((True, "in_span"), (False, "antisense")):
            for c in copies.values():
                if c["chrom"] == v["c"] and (c["strand"] == v["st"]) == want_same and exov(v["ex"], [list(c["span"])]) > 0:
                    return label
        return "elsewhere"

    classes = collections.Counter(klass(n) for n in nodes)
    classes = {k: classes.get(k, 0) for k in ("on_copy", "in_span", "antisense", "elsewhere")}
    # M6 cohesion: K* = the family cluster holding the most distinct copies (ties to the smaller id); a locus on a shared span is in each of its clusters
    held = collections.defaultdict(set)
    for n in nodes:
        for i in ids[n] & fam:
            held[i].update(on[n])
    kstar = min(held, key=lambda i: (-len(held[i]), i)) if held else None
    return dict(n_loci=len(loci), shared_spans=shared, family_clusters=sorted(fam), nodes=len(nodes), classes=classes,
                np=(classes["on_copy"] / len(nodes)) if nodes else 0.0, own=own, copies_with_node=len(own),
                kstar=kstar, copies_in_kstar=len(held[kstar]) if kstar else 0)


def copy_views(copies_tsv, truth, family, exons_json=None):
    rows = [r for r in csv.DictReader(open(copies_tsv), delimiter="\t") if r["family"] == family]
    genes = {r["isoform_gene"] for r in rows} | {r["cid"] for r in rows}
    tx = cs.gtf_transcripts(truth, genes)
    cu = {}
    for c in rows:
        texons = tx.get(c["cid"]) or tx.get(c["isoform_gene"], {})
        cu[c["cid"]] = cs.merge([e for ex in texons.values() for e in ex]) or [[int(c["terr_lo0"]), int(c["terr_hi"])]]
    if exons_json:
        ej = {c["cid"]: cs.merge([list(e) for e in c["exons"]]) for c in json.load(open(exons_json))["copies"]}
        cu = {cid: ej[cid] for cid in cu}
    return {c["cid"]: dict(chrom=c["chrom"], strand=c["strand"], exons=cu[c["cid"]], span=(int(c["terr_lo0"]), int(c["terr_hi"])),
                           name=c.get("refseq_name") or c.get("cat_name") or c["name"]) for c in rows}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--arm", action="append", required=True, help="NAME=loci.gff3,clusters.tsv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--report")
    ap.add_argument("--exons-json", help="npip_read_pool.json: the copy exons as the registered own-node rule took them (NPIP only)")
    a = ap.parse_args()
    copies = copy_views(a.copies, a.truth, a.family, a.exons_json)
    node = {cid: {} for cid in copies}
    report = {"family": a.family, "copies": len(copies), "arms": {}}
    for spec in a.arm:
        name, rest = spec.split("=", 1)
        loci_path, clusters_path = rest.split(",", 1)
        r = analyse(read_loci(loci_path), read_rows(clusters_path), copies)
        report["arms"][name] = r
        for cid in copies:
            node[cid][name] = cid in r["own"]
        print(f"{name}: {len(r['family_clusters'])} family clusters, {r['nodes']} nodes {r['classes']} NP {r['np']:.3f}, "
              f"{r['copies_with_node']} of {len(copies)} copies with an own node; {r['shared_spans']} shared spans")
    json.dump(dict(rows=[dict(cid=cid, name=c["name"], node=node[cid]) for cid, c in copies.items()]), open(a.out, "w"), indent=1)
    if a.report:
        json.dump(report, open(a.report, "w"), indent=1)


if __name__ == "__main__":
    main()

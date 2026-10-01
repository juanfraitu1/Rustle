#!/usr/bin/env python3
"""Figure data for the NPIP read-pool page (docs/PREREG_npip_read_pool_2026-10-01.md): every locus of each arm inside a few chr16
windows (span, strand, representative exons, class, and what the family step did with it), plus the per-copy locus counts.

    python3 bench/npip_read_pool/figdata.py --dir /mnt/linuxdisk/tmp/readpool_npip --scored npip_read_pool.json --cat-genes genes.tsv \
        --out figdata.json
"""
import argparse
import collections
import csv
import json

ARMS = ("P", "GOOD", "ALL")
WINDOWS = [("NPIPB4", 22_326_000, 22_392_000), ("NPIPA6 and NPIPA7", 16_316_000, 16_412_000), ("NPIPB2", 11_955_000, 12_020_000)]


def gff(path):
    out = {}
    for ln in open(path):
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "gene":
            out[at["Name"]] = dict(c=f[0], s0=int(f[3]) - 1, e=int(f[4]), st=f[6], ex=[])
        elif f[2] == "exon":
            out[at["gene"]]["ex"].append([int(f[3]) - 1, int(f[4])])
    return out


def exov(a, b):
    return sum(max(0, min(y, q) - max(x, p)) for x, y in a for p, q in b)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("dir", "scored", "cat_genes", "out"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    S = json.load(open(a.scored))
    copies = S["copies"]
    catg = [r for r in csv.DictReader(open(a.cat_genes), delimiter="\t") if r["chrom"] == "chr16"]
    out = dict(windows=[], copies=copies, per_copy={arm: S["summary"]["arms"][arm]["loci_per_copy"] for arm in ARMS})
    arms_loci = {}
    for arm in ARMS:
        L = gff(f"{a.dir}/{arm}.gff3")
        by = {(v["c"], v["s0"] + 1, v["e"]): n for n, v in L.items()}
        cl = {}
        for r in csv.DictReader(open(f"{a.dir}/{arm}.fam.clusters.tsv"), delimiter="\t"):
            cl[by[(r["chrom"], int(r["start"]), int(r["end"]))]] = r["cluster_id"]
        fold = {}
        for r in csv.DictReader(open(f"{a.dir}/{arm}.fam.loci.tsv"), delimiter="\t"):
            def key(x):
                c, se = x.rsplit(":", 1)
                s_, e_ = se.split("-")
                return (c, int(s_), int(e_))
            ann, rep = key(r["annotation"]), key(r["representative"])
            if ann in by and rep in by:
                fold[by[ann]] = by[rep]
        # NPIP clusters = every cluster holding a locus whose rep exons overlap an NPIP copy's exons (MCL may split NPIP)
        fam = {cl[n] for n, v in L.items() if n in cl and any(c["strand"] == v["st"] and exov(v["ex"], c["exons"]) > 0 for c in copies)}
        arms_loci[arm] = (L, cl, fold, fam)
    for title, lo, hi in WINDOWS:
        w = dict(title=title, lo=lo, hi=hi, copies=[c for c in copies if c["s0"] < hi and lo < c["e"]],
                 genes=[dict(name=g["gene_name"], st=g["strand"], s0=int(g["start0"]), e=int(g["end"]))
                        for g in catg if int(g["start0"]) < hi and lo < int(g["end"]) and int(g["end"]) - int(g["start0"]) < 400_000],
                 arms={})
        for arm in ARMS:
            L, cl, fold, fam = arms_loci[arm]
            items = []
            for n, v in L.items():
                if not (v["s0"] < hi and lo < v["e"]):
                    continue
                on = [c["name"] for c in copies if c["strand"] == v["st"] and exov(v["ex"], c["exons"]) > 0]
                anti = [c["name"] for c in copies if c["strand"] != v["st"] and c["s0"] < v["e"] and v["s0"] < c["e"]]
                span = [c["name"] for c in copies if c["strand"] == v["st"] and exov([[v["s0"], v["e"]]], c["exons"]) > 0]
                klass = "copy" if on else "antisense" if anti else "offexon" if span else "other"
                if n in cl:
                    fstat = "npip_node" if cl[n] in fam else "other_family"
                elif n in fold:
                    fstat = "folded_npip" if cl.get(fold[n]) in fam else "folded_other"
                else:
                    fstat = "singleton"
                items.append(dict(n=n, s0=v["s0"], e=v["e"], st=v["st"], ex=v["ex"], k=klass, f=fstat, on=on, span=span))
            items.sort(key=lambda x: (x["s0"], x["e"]))
            w["arms"][arm] = items
        out["windows"].append(w)
    json.dump(out, open(a.out, "w"), separators=(",", ":"))
    for w in out["windows"]:
        print(w["title"], {arm: (len(v), dict(collections.Counter(x["f"] for x in v))) for arm, v in w["arms"].items()})


if __name__ == "__main__":
    main()

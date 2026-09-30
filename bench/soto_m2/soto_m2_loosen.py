#!/usr/bin/env python3
"""What loosening our copy-number rule does (2026-09-30): the per-pair MAD cut with OUR copy numbers swept over a grid,
each family re-classified exactly as in soto_m2_families.py, scored on the frozen DEV / HELD-OUT split.

    python3 bench/soto_m2/soto_m2_loosen.py --geneset elig.tsv --full-geneset full.tsv \
        --famcn-ours famcn_ours_allwssd.tsv [--grid 1,1.25,1.5,2,3,4,8]

Soto's rule is MAD < 1 (two genes join when their copy numbers differ by less than 2). Loosening it can only merge
more, so it can recover families we cut (`we_miss`/`mixed` -> `match`) and can break families we had exactly
(`match` -> `soto_smaller`/`mixed`). Both directions are counted per rung; the grid is a display, not a fit.
"""
import argparse
import csv
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "soto"))
sys.path.insert(0, HERE)
import soto_replication as sr  # noqa: E402
from soto_m2_families import classify  # noqa: E402

SOTO = os.path.join(HERE, "..", "soto")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", default=os.path.join(SOTO, "soto_famCN_S1C.tsv"))
    ap.add_argument("--shared", default=os.path.join(SOTO, "shared_exons_5154_exon_mapback.tsv"))
    ap.add_argument("--sedef-edges", nargs="+", default=[os.path.join(SOTO, "shared_exons_2334_finalv1_native.tsv"),
                                                         os.path.join(SOTO, "shared_exons_2334_finalhuman.tsv")])
    ap.add_argument("--split", default=os.path.join(SOTO, "soto_split_2026-09-29.tsv"))
    ap.add_argument("--geneset", required=True)
    ap.add_argument("--full-geneset", required=True)
    ap.add_argument("--famcn-ours", required=True)
    ap.add_argument("--grid", default="1,1.25,1.5,2,3,4,8")
    ap.add_argument("--show", default="ID_41,ID_154,ID_149,ID_356,ID_397,ID_396,ID_400")
    a = ap.parse_args(argv)

    genes, _ = sr.load_geneset(a.geneset)
    full, _ = sr.load_geneset(a.full_geneset)
    edges = sr.read_edges(a.shared)
    sedef = [e for p in a.sedef_edges for e in sr.read_edges(p)]
    clean, _amb = sr.load_truth(a.truth)
    info, fam_all, manual = {}, defaultdict(set), set()
    for r in csv.DictReader(open(a.truth), delimiter="\t"):
        g = r["Gene ID"]
        info.setdefault(g, dict(name=r["Gene Name"], biotype=r["Biotype"], fams=[],
                                backbone=r.get("In Table S1 (SD98 gene set)") == "Yes"))
        if r["Family ID"]:
            info[g]["fams"].append(r["Family ID"])
            if not r["Family ID"].startswith("Unassigned"):
                fam_all[r["Family ID"]].add(g)
        if r.get("Family MAD", "").strip().lower() == "manual merge":
            manual.add(r["Family ID"])
    cn_s1c = {g: v for g, v in sr.load_famcn(a.truth, "Median famCN").items() if g in genes}
    cn_ours = {g: v for g, v in sr.load_famcn(a.famcn_ours, "famCN_sotoiv").items() if g in genes}
    half = sr.load_split(a.split)
    keep = set(info) | genes | full
    show = a.show.split(",")

    base = None
    print("| MAD cut (ours) | ALL ARI / exact | DEV exact | HELD-OUT exact | Soto smaller | we miss · our copy numbers "
          "| recovered vs MAD<1 | broken vs MAD<1 | " + " | ".join(show) + " |")
    print("|---|---|---|---|---|---|---|---|" + "---|" * len(show))
    for t in [float(x) for x in a.grid.split(",")]:
        rows, anc, _ps, _po = classify(keep, genes, full, edges, sedef, clean, fam_all, info, manual, cn_s1c, cn_ours,
                                       mad=t)
        cls = {r["family"]: r for r in rows}
        fam_half = {r["family"]: Counter(half.get(g) for g in r["clean"]).most_common(1)[0][0] for r in rows if r["clean"]}
        ex = Counter(fam_half.get(f) for f, r in cls.items() if r["cls"] == "match")
        k = Counter(r["cls"] for r in rows)
        our_cn = sum(1 for r in rows if r["miss_cause"] == "our_cn")
        if base is None:
            base = cls
        rec = sorted((f for f in cls if cls[f]["cls"] == "match" and base[f]["cls"] != "match"), key=lambda x: int(x[3:]))
        brk = sorted((f for f in cls if cls[f]["cls"] != "match" and base[f]["cls"] == "match"), key=lambda x: int(x[3:]))
        cells = [cls[f]["cls"] + (f" ({cls[f]['miss_cause']})" if cls[f]["miss_cause"] else "") for f in show]
        print(f"| < {t:g} | {anc['ari']:.4f} / {anc['exact']} | {ex['dev']} | {ex['heldout']} | {k['soto_smaller']} | "
              f"{our_cn} | {len(rec)} | {len(brk)} | " + " | ".join(cells) + " |")
        if rec or brk:
            print(f"    recovered: {', '.join(rec) or '-'}", file=sys.stderr)
            print(f"    broken:    {', '.join(brk) or '-'}", file=sys.stderr)


if __name__ == "__main__":
    main()

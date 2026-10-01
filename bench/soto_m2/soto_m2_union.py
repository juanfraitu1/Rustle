#!/usr/bin/env python3
"""Pre-registered KEY=unionedges (docs/PREREG_soto_union_edges_2026-09-30.md): exon map-back edges vs map-back plus every
SEDEF-projected exon link, each family classified exactly as in soto_m2_families.py, with S1C and our copy numbers and the
four biotype filters; exact families on all / dev / held-out, and per family what the union recovers or breaks.

    python3 bench/soto_m2/soto_m2_union.py --geneset elig.tsv --full-geneset full.tsv --famcn-ours famcn_ours_allwssd.tsv
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
from soto_m2_families import COMBOS, classify, is_pseudo  # noqa: E402

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
    a = ap.parse_args(argv)

    genes, _ = sr.load_geneset(a.geneset)
    full, _ = sr.load_geneset(a.full_geneset)
    mapback = sr.read_edges(a.shared)
    sedef = [e for p in a.sedef_edges for e in sr.read_edges(p)]
    seen, union = set(), []
    for x, y in mapback + sedef:
        k = tuple(sorted((x, y)))
        if x != y and k not in seen:
            seen.add(k)
            union.append((x, y))
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
    every = set(info) | genes | full
    nmb = len({tuple(sorted(e)) for e in mapback if e[0] != e[1]})
    print(f"# KEY=unionedges\nedges: map-back {nmb:,}; union {len(union):,} (+{len(union) - nmb:,} SEDEF-only)\n")
    print("| copy numbers | filter | map-back exact (dev / held-out) | union exact (dev / held-out) | ARI map-back -> union "
          "| recovered | broken |")
    print("|---|---|---|---|---|---|---|")
    verdict, details = None, []
    for src, cn_use in (("s1c", cn_s1c), ("ours", cn_ours)):
        for key0, label0, no_p, no_l in COMBOS:
            keep = {g for g in every if not (no_p and is_pseudo(info.get(g, {}).get("biotype", "")))
                    and not (no_l and info.get(g, {}).get("biotype", "") == "lncRNA")}
            res = {}
            for arm, edges in (("mapback", mapback), ("union", union)):
                rows, anc, _ps, _po = classify(keep, genes, full, edges, sedef, clean, fam_all, info, manual, cn_s1c,
                                               cn_use, filtered=bool(key0))
                cls = {r["family"]: r for r in rows}
                fam_half = {r["family"]: Counter(half.get(g) for g in r["clean"]).most_common(1)[0][0]
                            for r in rows if r["clean"]}
                ex = Counter(fam_half.get(f) for f, r in cls.items() if r["cls"] == "match")
                res[arm] = (cls, anc, ex, sum(r["cls"] == "match" for r in rows))
            (b, ba, bex, bn), (u, ua, uex, un) = res["mapback"], res["union"]
            rec = sorted((f for f in u if u[f]["cls"] == "match" and b[f]["cls"] != "match"), key=lambda x: int(x[3:]))
            brk = sorted((f for f in u if u[f]["cls"] != "match" and b[f]["cls"] == "match"), key=lambda x: int(x[3:]))
            print(f"| {src} | {label0} | {bn} ({bex['dev']} / {bex['heldout']}) | {un} ({uex['dev']} / {uex['heldout']}) | "
                  f"{ba['ari']:.4f} -> {ua['ari']:.4f} | {len(rec)} | {len(brk)} |")
            details.append((src, label0, rec, brk, b, u))
            if src == "s1c" and key0 == "":
                if not rec and not brk:
                    verdict = "NO CHANGE"
                elif un > bn and not brk and uex["heldout"] >= bex["heldout"]:
                    verdict = "UNION HELPS"
                elif brk and un <= bn:
                    verdict = "UNION HURTS"
                elif un > bn and brk:
                    verdict = "MIXED"
                else:
                    verdict = "UNION HURTS" if brk else "NO CHANGE (held-out drop)"
    print(f"\nVERDICT (S1C copy numbers, all genes): {verdict}\n")
    for src, label0, rec, brk, b, u in details:
        if rec or brk:
            print(f"- {src}, {label0}:")
            for f in rec:
                print(f"  recovered {f} ({b[f]['cls']}{', ' + b[f]['miss_cause'] if b[f]['miss_cause'] else ''})")
            for f in brk:
                print(f"  broken    {f} -> {u[f]['cls']}{', ' + u[f]['miss_cause'] if u[f]['miss_cause'] else ''}"
                      f" (+{len(u[f]['extra'])} extra genes)")


if __name__ == "__main__":
    main()

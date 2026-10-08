#!/usr/bin/env python3
"""Pre-registered KEY=ourhomology (docs/archive/2026-09/PREREG_our_homology_s1c_cn_2026-09-30.md): our homology method, Soto's copy numbers.

Turns mcl_families' within-family homology edges over the loci of Soto's 2,334 genes into gene-level edges (a locus that folded
several overlapping annotation records links all its genes, and an edge between two loci links every gene of one to every gene of
the other), then scores them exactly like Soto's exon map-back edges with `soto_m2_families.classify`: sequence only, S1C famCN with
Soto's pair rule, our famCN. Reports all / dev / held-out ARI and exact families, nesting and bipartite numbers.

    python3 bench/soto_m2/soto_m2_our_homology.py --geneset elig.tsv --full-geneset full.tsv --famcn-ours famcn_ours_allwssd.tsv \
        --names names.tsv --pairs ours.pairs.tsv --folded ours.loci.tsv
"""
import argparse
import csv
import itertools
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "soto"))
sys.path.insert(0, HERE)
import soto_replication as sr  # noqa: E402
from soto_m2_families import classify  # noqa: E402

SOTO = os.path.join(HERE, "..", "soto")


def gene_edges(names, pairs, folded):
    locus_of = {}
    for line in open(names):
        g, loc = line.rstrip("\n").split("\t")
        locus_of[g] = loc
    rep = {}
    for r in csv.DictReader(open(folded), delimiter="\t"):
        rep[r["annotation"]] = r["representative"]
    genes_at = defaultdict(set)
    for g, loc in locus_of.items():
        genes_at[rep.get(loc, loc)].add(g)
    edges = set()
    for gs in genes_at.values():                      # genes folded into one locus are one node: link them
        for a, b in itertools.combinations(sorted(gs), 2):
            edges.add((a, b))
    for r in csv.DictReader(open(pairs), delimiter="\t"):
        for a in genes_at.get(rep.get(r["a"], r["a"]), ()):
            for b in genes_at.get(rep.get(r["b"], r["b"]), ()):
                if a != b:
                    edges.add(tuple(sorted((a, b))))
    return sorted(edges), sum(len(v) > 1 for v in genes_at.values())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", default=os.path.join(SOTO, "soto_famCN_S1C.tsv"))
    ap.add_argument("--shared", default=os.path.join(SOTO, "shared_exons_5154_exon_mapback.tsv"))
    ap.add_argument("--sedef-edges", nargs="+", default=[os.path.join(SOTO, "shared_exons_2334_finalv1_native.tsv"),
                                                         os.path.join(SOTO, "shared_exons_2334_finalhuman.tsv")])
    ap.add_argument("--split", default=os.path.join(SOTO, "soto_split_2026-09-29.tsv"))
    for k in ("geneset", "full_geneset", "famcn_ours", "names", "pairs", "folded"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)

    genes, _ = sr.load_geneset(a.geneset)
    full, _ = sr.load_geneset(a.full_geneset)
    mapback = sr.read_edges(a.shared)
    sedef = [e for p in a.sedef_edges for e in sr.read_edges(p)]
    ours, nfold = gene_edges(a.names, a.pairs, a.folded)
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
    universe = set(clean)
    uh = {h: {g for g in universe if half.get(g) == h} for h in ("dev", "heldout")}
    print(f"# KEY=ourhomology\nedges: Soto's exon map-back {len(mapback):,}; ours {len(ours):,} gene pairs "
          f"({nfold} loci folding several annotation records)\n")
    print("| homology | copy numbers | ARI all / dev / held-out | exact all (dev / held-out) | nesting | bipartite sens / prec |")
    print("|---|---|---|---|---|---|")
    res = {}
    for hname, edges in (("H0 Soto map-back", mapback), ("H1 ours", ours)):
        for cname, cn in (("none (sequence only)", None), ("S1C (Soto's)", cn_s1c), ("ours", cn_ours)):
            rows, anc, pseq, pcn = classify(keep, genes, full, edges, sedef, clean, fam_all, info, manual, cn_s1c,
                                            cn if cn is not None else cn_s1c)
            pred = pseq if cn is None else pcn
            sc = {k: sr.score(clean, pred, u) for k, u in (("all", universe), ("dev", uh["dev"]), ("heldout", uh["heldout"]))}
            if cn is None:
                ex_txt = f"{sc['all']['n_exact']}"
                bip = "-"
            else:
                fam_half = {r["family"]: Counter(half.get(g) for g in r["clean"]).most_common(1)[0][0] for r in rows if r["clean"]}
                exh = Counter(fam_half.get(r["family"]) for r in rows if r["cls"] == "match")
                ex_txt = f"{anc['exact']} ({exh['dev']} / {exh['heldout']})"
                bip = f"{anc['bip']['sens']:.3f} / {anc['bip']['prec']:.3f}"
            res[(hname, cname)] = sc
            print(f"| {hname} | {cname} | {sc['all']['ari']:.4f} / {sc['dev']['ari']:.4f} / {sc['heldout']['ari']:.4f} | {ex_txt} | "
                  f"{anc['nest']}/{anc['nest_n']} | {bip} |")
    h0 = res[("H0 Soto map-back", "S1C (Soto's)")]["heldout"]["ari"]
    h1 = res[("H1 ours", "S1C (Soto's)")]["heldout"]["ari"]
    verdict = "FINDS THEM" if h1 >= h0 - 0.05 else ("PARTIAL" if h1 >= h0 - 0.15 else "DOES NOT")
    print(f"\nVERDICT (H1 with S1C famCN, held-out ARI {h1:.4f} vs H0 {h0:.4f}): {verdict}")


if __name__ == "__main__":
    main()

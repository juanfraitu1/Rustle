#!/usr/bin/env python3
"""Every Soto S1C family against ours, one row per family, with the cause of each disagreement (2026-09-30).

    python3 bench/soto_m2/soto_m2_families.py --geneset elig.tsv --full-geneset full.tsv \
        --famcn-ours famcn_ours_allwssd.tsv --exons sd98_gene_exons.tsv --cat cat_v4.bed \
        --out-tsv families.tsv --out-json families.json

OURS = the reconciled recipe (exon map-back edges x Soto's per-pair MAD < 1 rule) with OUR copy numbers (268 SGDP
samples, Soto's gene-body ∩ SD98 interval): the ladder rung with ARI 0.9277 / 411 of 491 exact. A Soto family F (its
clean members = genes Soto places in F only) is compared with our families over the ladder's universe (all clean
S1C genes), exactly as `score` counts exact families.

  match          F's clean members are exactly one of our families.
  soto_smaller   all of F is in one of our families, and that family also holds other genes (joined to F by >= 98%
                 identical exons). Cause per extra gene: `cn` (Soto's copy numbers split it off, ours keep it),
                 `soto_rule` (Soto's own rule with Soto's own copy numbers also joins it: their table is narrower than
                 their rule), `unassigned` (Soto left the gene without a family).
  we_miss        F is spread over >= 2 of our families (or a member is left unplaced). Cause, first that applies:
                 `our_cn`      the pieces are joined by >= 98% exon links, our copy numbers cut them (Soto's do not);
                 `soto_rule`   Soto's own rule with Soto's own copy numbers cuts them too (table not reproducible);
                 `our_edges`   no exon map-back link joins the pieces, but a SEDEF-projected link does (our graph gap);
                 `soto_manual` no link in any edge set; S1C marks the family "Manual merge";
                 `soto_no_seq` no link in any edge set, no manual flag: Soto groups genes with no >= 98% exon evidence.
  mixed          both: part of F is cut away, and the rest sits with extra genes.
  gone           (filtered runs only) fewer than 2 of F's members are left once the excluded biotypes are removed.

`narrower_than_homology` flags a family whose members all sit in one sequence-only component (no copy-number gate)
that also holds other genes; it is a property of the copy-number step, independent of the class.

The whole analysis is run four times: all genes, without pseudogenes (any biotype containing "pseudogene"), without
lncRNAs, and without both. Excluded genes are removed everywhere (Soto's families, our graph, the copy-number gate, the
score) before clustering, so a pseudogene can no longer bridge two coding genes.
"""
import argparse
import csv
import json
import os
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "soto"))
import soto_replication as sr  # noqa: E402

SOTO = os.path.join(HERE, "..", "soto")
COMBOS = [("", "all genes", False, False), ("p", "no pseudogenes", True, False), ("l", "no lncRNAs", False, True),
          ("pl", "no pseudogenes or lncRNAs", True, True)]


def components_over(edges, nodes):
    """Union-find over `edges` restricted to `nodes` -> {node: root}."""
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edges:
        if a in parent and b in parent:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
    return {n: find(n) for n in nodes}


def is_pseudo(bt):
    return "pseudogene" in bt


def classify(keep, genes0, full0, edges0, sedef0, clean0, fam_all0, info, manual, cn_s1c0, cn_ours0):
    """One run of the comparison over the genes in `keep`. Returns (rows, anchors, pred_seq, pred_ours)."""
    genes, full = genes0 & keep, full0 & keep
    edges = [(x, y) for x, y in edges0 if x in keep and y in keep]
    sedef = [(x, y) for x, y in sedef0 if x in keep and y in keep]
    clean = {g: f for g, f in clean0.items() if g in keep}
    fam_all = {f: {g for g in m if g in keep} for f, m in fam_all0.items()}
    universe = set(clean)
    cn_s1c = {g: v for g, v in cn_s1c0.items() if g in genes}
    cn_ours = {g: v for g, v in cn_ours0.items() if g in genes}

    def run(famcn, gate):
        cover, kept, leaf_of = sr.pair_families(edges, genes, full, famcn, gate=gate)
        return sr.collapse_cover(cover, leaf_of, genes)
    pred_seq, pred_ours, pred_s1c = run({}, False), run(cn_ours, True), run(cn_s1c, True)

    def fams_over(pred):
        out = defaultdict(set)
        for g in universe:
            if pred.get(g):
                out[pred[g]].add(g)
        return out
    ours_f, seq_f, s1c_f = fams_over(pred_ours), fams_over(pred_seq), fams_over(pred_s1c)
    s1c_sets = {frozenset(v) for v in s1c_f.values()}
    members = defaultdict(set)
    for g, f in clean.items():
        if not f.startswith("Unassigned"):
            members[f].add(g)
    comp = None
    out = []
    for f in sorted(fam_all0, key=lambda x: int(x.split("_")[1])):
        M = members.get(f, set())
        if len(fam_all[f]) < 2 or not M:
            out.append(dict(family=f, cls="gone", miss_cause="", extra_causes={}, n_members=len(fam_all[f]),
                            n_clean=len(M), n_ours_pieces=0, n_extra=0, narrower_than_homology=False, seq_size=0,
                            manual=f in manual, exact_with_s1c_cn=False, clean=sorted(M), extra=[], extra_cause={},
                            seq_members=[], ours_members=[]))
            continue
        labs = {pred_ours.get(g) for g in M}
        placed = [g for g in M if pred_ours.get(g)]
        exact = None not in labs and len(labs) == 1 and ours_f[next(iter(labs))] == M
        split = None in labs or len(labs) > 1
        ours_union = set().union(*(ours_f[l] for l in labs if l)) if placed else set()
        extra = ours_union - M
        extra_cause = {}
        for x in sorted(extra):
            tf = clean.get(x, "")
            if tf.startswith("Unassigned") or not tf:
                extra_cause[x] = "unassigned"
            elif pred_s1c.get(x) and pred_s1c.get(x) in {pred_s1c.get(g) for g in M}:
                extra_cause[x] = "soto_rule"
            else:
                extra_cause[x] = "cn"
        miss_cause = ""
        seq_labs = {pred_seq.get(g) for g in M}
        if split:
            s1c_labs = {pred_s1c.get(g) for g in M}
            if None not in seq_labs and len(seq_labs) == 1:
                miss_cause = "soto_rule" if (None in s1c_labs or len(s1c_labs) > 1) else "our_cn"
            else:
                if comp is None:
                    comp = components_over(edges + sedef, set(full) | set(universe))
                if len({comp[g] for g in M}) == 1:
                    miss_cause = "our_edges"
                elif f in manual:
                    miss_cause = "soto_manual"
                else:
                    miss_cause = "soto_no_seq"
        cls = "match" if exact else ("mixed" if split and extra else ("we_miss" if split else "soto_smaller"))
        narrower = (None not in seq_labs and len(seq_labs) == 1 and len(seq_f[next(iter(seq_labs))]) > len(M))
        seq_union = set().union(*(seq_f[l] for l in seq_labs if l)) if any(seq_labs) else set()
        out.append(dict(
            family=f, cls=cls, miss_cause=miss_cause,
            extra_causes=dict(sorted({c: sum(1 for v in extra_cause.values() if v == c)
                                      for c in set(extra_cause.values())}.items())),
            n_members=len(fam_all[f]), n_clean=len(M),
            n_ours_pieces=len({l for l in labs if l}) + sum(1 for g in M if not pred_ours.get(g)),
            n_extra=len(extra), narrower_than_homology=narrower, seq_size=len(seq_union),
            manual=f in manual, exact_with_s1c_cn=frozenset(M) in s1c_sets,
            clean=sorted(M), extra=sorted(extra), extra_cause=extra_cause, seq_members=sorted(seq_union),
            ours_members=sorted(set().union(*(ours_f[l] for l in labs if l)) if placed else set())))
    # anchors: the ladder's score, plus sequence-only nesting of Soto families with >= 2 clean members
    sc = sr.score(clean, pred_ours, universe)
    sc1 = sr.score(clean, pred_s1c, universe)
    multi = {f: m for f, m in members.items() if len(m) >= 2}
    nest = sum(1 for m in multi.values() if None not in {pred_seq.get(g) for g in m}
               and len({pred_seq.get(g) for g in m}) == 1)
    live = [r for r in out if r["cls"] != "gone"]
    anchors = dict(ari=round(sc["ari"], 4), exact=sc["n_exact"], ari_s1c=round(sc1["ari"], 4), exact_s1c=sc1["n_exact"],
                   families=len(live), gone=len(out) - len(live), nest=nest, nest_n=len(multi),
                   narrower=sum(r["narrower_than_homology"] for r in live), genes=len(keep))
    return out, anchors, pred_seq, pred_ours


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", default=os.path.join(SOTO, "soto_famCN_S1C.tsv"))
    ap.add_argument("--shared", default=os.path.join(SOTO, "shared_exons_5154_exon_mapback.tsv"))
    ap.add_argument("--sedef-edges", nargs="+", default=[os.path.join(SOTO, "shared_exons_2334_finalv1_native.tsv"),
                                                         os.path.join(SOTO, "shared_exons_2334_finalhuman.tsv")])
    ap.add_argument("--geneset", required=True)
    ap.add_argument("--full-geneset", required=True)
    ap.add_argument("--famcn-ours", required=True)
    ap.add_argument("--exons", required=True)
    ap.add_argument("--cat", help="CAT v4 BED12 (gene strand); optional")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-json", required=True)
    a = ap.parse_args(argv)

    genes, _ = sr.load_geneset(a.geneset)
    full, _ = sr.load_geneset(a.full_geneset)
    edges = sr.read_edges(a.shared)
    sedef = [e for p in a.sedef_edges for e in sr.read_edges(p)]
    clean, _ambiguous = sr.load_truth(a.truth)
    info, fam_all, manual = {}, defaultdict(set), set()
    for r in csv.DictReader(open(a.truth), delimiter="\t"):
        g = r["Gene ID"]
        info.setdefault(g, dict(name=r["Gene Name"], biotype=r["Biotype"], fams=[]))
        if r["Family ID"]:
            info[g]["fams"].append(r["Family ID"])
            if not r["Family ID"].startswith("Unassigned"):
                fam_all[r["Family ID"]].add(g)
        if r.get("Family MAD", "").strip().lower() == "manual merge":
            manual.add(r["Family ID"])
    cn_s1c = {g: v for g, v in sr.load_famcn(a.truth, "Median famCN").items() if g in genes}
    cn_ours = {g: v for g, v in sr.load_famcn(a.famcn_ours, "famCN_sotoiv").items() if g in genes}
    every = set(info) | genes | full

    runs = {}
    for key, label, no_p, no_l in COMBOS:
        keep = {g for g in every if not (no_p and is_pseudo(info.get(g, {}).get("biotype", "")))
                and not (no_l and info.get(g, {}).get("biotype", "") == "lncRNA")}
        runs[key] = classify(keep, genes, full, edges, sedef, clean, fam_all, info, manual, cn_s1c, cn_ours)
        rows, anc = runs[key][0], runs[key][1]
        counts = defaultdict(int)
        for r in rows:
            counts[r["cls"]] += 1
        assert counts["match"] == anc["exact"] or key, (counts["match"], anc["exact"])
        print(f"[{label}] genes {anc['genes']}; ARI {anc['ari']:.4f} / {anc['exact']} exact (S1C CN {anc['ari_s1c']:.4f} / "
              f"{anc['exact_s1c']}); nest {anc['nest']}/{anc['nest_n']}; narrower {anc['narrower']}; "
              + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())), file=sys.stderr)
    base = runs[""][0]
    assert sum(r["cls"] == "match" for r in base) == runs[""][1]["exact"]

    with open(a.out_tsv, "w") as fh:
        fh.write("family_id\tclass\tmiss_cause\textra_causes\tn_members\tn_clean\tn_our_pieces\tn_extra\t"
                 "narrower_than_homology\tsequence_family_size\tmanual_merge\texact_with_s1c_cn\texample_genes\t"
                 "class_no_pseudogenes\tclass_no_lncrnas\tclass_no_both\n")
        for i, r in enumerate(base):
            ec = ";".join(f"{k}={v}" for k, v in r["extra_causes"].items())
            top = sorted(r["clean"] or fam_all[r["family"]], key=lambda g: (info[g]["biotype"] != "protein_coding", info[g]["name"]))
            fh.write(f"{r['family']}\t{r['cls']}\t{r['miss_cause']}\t{ec}\t{r['n_members']}\t{r['n_clean']}\t"
                     f"{r['n_ours_pieces']}\t{r['n_extra']}\t{int(r['narrower_than_homology'])}\t{r['seq_size']}\t"
                     f"{int(r['manual'])}\t{int(r['exact_with_s1c_cn'])}\t{','.join(info[g]['name'] for g in top[:3])}\t"
                     f"{runs['p'][0][i]['cls']}\t{runs['l'][0][i]['cls']}\t{runs['pl'][0][i]['cls']}\n")

    # JSON for the viewer: genes (coordinates, exons, copy numbers), links, and per filter: labels and family rows
    exons = defaultdict(list)
    with open(a.exons) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["gene_id"] in info:
                exons[r["gene_id"]].append((r["chrom"], int(r["start"]), int(r["end"])))
    strand = {}
    if a.cat:
        with open(a.cat) as fh:
            for line in fh:
                p = line.split("\t")
                if len(p) > 19 and p[18] in info and p[18] not in strand:
                    strand[p[18]] = p[5]
    gl = sorted(info, key=lambda g: (exons[g][0][0] if exons[g] else "z", min(s for _, s, _ in exons[g]) if exons[g] else 0))
    idx = {g: i for i, g in enumerate(gl)}
    gene_rows = []
    for g in gl:
        merged = []
        for c, s, e in sorted(exons[g], key=lambda t: t[1]):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        bt = info[g]["biotype"]
        gene_rows.append(dict(
            id=g, n=info[g]["name"], bt=bt, k="p" if is_pseudo(bt) else ("l" if bt == "lncRNA" else ""),
            c=exons[g][0][0] if exons[g] else "", s=merged[0][0] if merged else 0,
            e=max(e for _, e in merged) if merged else 0, st=strand.get(g, "."), x=merged, sf=info[g]["fams"],
            cs=round(cn_s1c[g], 1) if g in cn_s1c else None, co=round(cn_ours[g], 1) if g in cn_ours else None))
    link_set = {tuple(sorted((idx[x], idx[y]))) for x, y in edges if x in idx and y in idx and x != y}
    sedef_only = {tuple(sorted((idx[x], idx[y]))) for x, y in sedef if x in idx and y in idx and x != y} - link_set
    combos = {}
    for key, label, _p, _l in COMBOS:
        rows, anc, pred_seq, pred_ours = runs[key]
        fams = {}
        for r in rows:
            show = set(fam_all[r["family"]]) | set(r["extra"]) | set(r["seq_members"]) | set(r["ours_members"])
            fams[r["family"]] = dict(
                k=r["cls"], mc=r["miss_cause"], ec=r["extra_causes"], nm=r["n_members"], nc=r["n_clean"],
                np=r["n_ours_pieces"], nx=r["n_extra"], nh=r["narrower_than_homology"], ns=r["seq_size"],
                x1=r["exact_with_s1c_cn"], g=sorted(idx[g] for g in show if g in idx),
                xc={str(idx[g]): c for g, c in r["extra_cause"].items()})
        combos[key] = dict(label=label, anchor=anc, seq=[pred_seq.get(g, "") for g in gl],
                           our=[pred_ours.get(g, "") for g in gl], fams=fams)
    fam_meta = []
    for r in base:
        top = sorted(r["clean"] or fam_all[r["family"]], key=lambda g: (info[g]["biotype"] != "protein_coding", info[g]["name"]))
        fam_meta.append(dict(f=r["family"], ex=",".join(info[g]["name"] for g in top[:3]), mm=r["manual"]))
    blob = dict(genes=gene_rows, links=sorted(link_set), sedef_links=sorted(sedef_only), families=fam_meta,
                combos=combos)
    with open(a.out_json, "w") as fh:
        json.dump(blob, fh, separators=(",", ":"))


if __name__ == "__main__":
    main()

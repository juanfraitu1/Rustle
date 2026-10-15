#!/usr/bin/env python3
"""Post hoc, descriptive (not pre-registered): why the de novo -copies analogue's bars fail
(docs/PREREG_denovo_copies_2026-10-14.md, Results).

    copies_posthoc.py --species S --base BASE.gtf --run DIR --out J [--liftoff ROOT]

recall  per reachable copy: found by an own template; found by another template off that template's own locus (the
        projection proper) or only by a template at the copy (the copy is itself assembled); whether the best own
        template is gene-like (it covers >= 80% of the source placement's exons and they cover >= 80% of it).
sites   the class `none` sites: template single-exon or not, template at a Liftoff gene or novel, the contigs holding
        the most of them.
"""
import argparse
import collections
import json
import os
import sys
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import copies as K  # noqa: E402

GENE_LIKE = 0.80


def own_at(templates, base_rows):
    """{template: its own locus's exon union}."""
    locus_iv = {(b["contig"], b["locus"]): b["iv"] for b in base_rows}
    return {tid: locus_iv.get((t.contig, t.gene), t.exons) for tid, t in templates.items()}


def is_self(h, templates, own_iv):
    t = templates.get(h["template"])
    return t is not None and h["contig"] == t.contig and K.cov(h["iv"], own_iv[h["template"]]) >= K.LF.MATCH_COV


def recall_detail(loci, base_rows, templates, sc):
    own, _ = K.own_templates(loci, base_rows, templates)
    own_iv = own_at(templates, base_rows)
    placement = {r["source_id"]: r for r in loci if r["cls"] in ("in_place", "moved") and not r["short"]}
    pidx, bidx = K.LF.Index(sc["projected"]), K.LF.Index(base_rows)
    rows = []
    for r in loci:
        if r["cls"] != "extra_copy" or r["short"] or not own.get(r["source_id"]):
            continue
        s, ts = r["source_id"], own[r["source_id"]]
        hits = K._covering(pidx, r)
        p = placement[s]
        like = max(min(K.cov(p["iv"], templates[t].exons), K.cov(templates[t].exons, p["iv"])) for t in ts)
        rows.append(dict(own=any(h["template"] in ts for h in hits),
                         projection=any(not is_self(h, templates, own_iv) for h in hits),
                         assembled=bidx.best_cov(r["contig"], r["iv"])[0] >= K.LF.MATCH_COV,
                         gene_like=like >= GENE_LIKE, any=bool(hits)))
    c = collections.Counter
    out = dict(reachable=len(rows), own=sum(x["own"] for x in rows), any=sum(x["any"] for x in rows),
               any_by_projection=sum(x["projection"] for x in rows),
               any_only_at_assembled_copy=sum(x["any"] and not x["projection"] for x in rows),
               copy_assembled=sum(x["assembled"] for x in rows))
    for key in (True, False):
        sub = [x for x in rows if x["gene_like"] is key]
        out[f"gene_like_{str(key).lower()}"] = [sum(x["own"] for x in sub), len(sub)]
    out["not_own_but_projection"] = sum(x["projection"] and not x["own"] for x in rows)
    out["not_own_any"] = dict(c((x["projection"], x["assembled"]) for x in rows if x["any"] and not x["own"]).most_common())
    out["not_own_any"] = {f"projection={a},assembled={b}": n for (a, b), n in out["not_own_any"].items()}
    return out


def none_sites(loci, base_rows, templates, projected):
    own_iv = own_at(templates, base_rows)
    off = [h for h in projected if not is_self(h, templates, own_iv)]
    extra = K.LF.Index([r for r in loci if r["cls"] == "extra_copy"])
    genes = K.LF.Index([r for r in loci if r["cls"] in ("in_place", "moved")])
    bidx = K.LF.Index(base_rows)
    none = []
    for g in K.site_groups(off):
        best = min(g, key=lambda h: (-h["seq"], -h["cov"], h["template"], h["strand"]))
        if all(idx.best_cov(best["contig"], best["iv"])[0] < K.LF.MATCH_COV for idx in (extra, genes, bidx)):
            none.append(best)
    t_at_gene = {tid: genes.best_cov(t.contig, t.exons)[0] >= K.LF.MATCH_COV for tid, t in templates.items()}
    c = collections.Counter
    return dict(n=len(none),
                template_single_exon=sum(len(templates[h["template"]].exons) == 1 for h in none),
                template_at_liftoff_gene=sum(t_at_gene[h["template"]] for h in none),
                hit_unspliced=sum(len(h["iv"]) == 1 for h in none),
                contigs=c(h["contig"] for h in none).most_common(6),
                templates=len({h["template"] for h in none}))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("--species", "--base", "--run", "--out"):
        ap.add_argument(k, required=True)
    ap.add_argument("--liftoff", default=K.LIFTOFF)
    a = ap.parse_args(argv)
    templates = K.T.read_templates(os.path.join(a.run, "templates.tsv"))
    loci = K.LF.read_loci(Path(a.liftoff) / a.species / "liftoff_loci.tsv")
    base_rows = K.LF.denovo_loci(Path(a.base), Path(os.path.dirname(os.path.abspath(a.out))) / "denovo_cache")
    sc = K.scan(K.batch_pafs(a.run), templates, set())
    res = dict(recall=recall_detail(loci, base_rows, templates, sc),
               none_sites=none_sites(loci, base_rows, templates, sc["projected"]))
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=1, sort_keys=True)
        fh.write("\n")
    print(json.dumps(res, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())

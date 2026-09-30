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

`narrower_than_homology` flags a family whose members all sit in one sequence-only component (no copy-number gate)
that also holds other genes; it is a property of the copy-number step, independent of the class.
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
    clean, ambiguous = sr.load_truth(a.truth)
    universe = set(clean)
    rows = list(csv.DictReader(open(a.truth), delimiter="\t"))
    info, fam_all, manual = {}, defaultdict(set), set()
    for r in rows:
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

    def run(famcn, gate):
        cover, kept, leaf_of = sr.pair_families(edges, genes, full, famcn, gate=gate)
        return sr.collapse_cover(cover, leaf_of, genes), kept
    pred_seq, kept_seq = run({}, False)
    pred_ours, kept_ours = run(cn_ours, True)
    pred_s1c, _ = run(cn_s1c, True)
    anchor = sr.score(clean, pred_ours, universe)
    anchor_s1c = sr.score(clean, pred_s1c, universe)
    print(f"anchor: ours ARI {anchor['ari']:.4f} / {anchor['n_exact']} exact; "
          f"S1C {anchor_s1c['ari']:.4f} / {anchor_s1c['n_exact']} exact", file=sys.stderr)

    def fams_over(pred):
        out = defaultdict(set)
        for g in universe:
            if pred.get(g):
                out[pred[g]].add(g)
        return out
    ours_f, seq_f, s1c_f = fams_over(pred_ours), fams_over(pred_seq), fams_over(pred_s1c)
    sedef = [e for p in a.sedef_edges for e in sr.read_edges(p)]
    members = defaultdict(set)
    for g, f in clean.items():
        if not f.startswith("Unassigned"):
            members[f].add(g)

    out = []
    for f in sorted(members, key=lambda x: int(x.split("_")[1])):
        M = members[f]
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
        if split:
            seq_labs = {pred_seq.get(g) for g in M}
            s1c_labs = {pred_s1c.get(g) for g in M}
            if None not in seq_labs and len(seq_labs) == 1:
                miss_cause = "soto_rule" if (None in s1c_labs or len(s1c_labs) > 1) else "our_cn"
            else:
                comp = components_over(edges + sedef, set(M) | {g for g in full})
                if len({comp[g] for g in M}) == 1:
                    miss_cause = "our_edges"
                elif f in manual:
                    miss_cause = "soto_manual"
                else:
                    miss_cause = "soto_no_seq"
        cls = "match" if exact else ("mixed" if split and extra else ("we_miss" if split else "soto_smaller"))
        seq_labs = {pred_seq.get(g) for g in M}
        narrower = (None not in seq_labs and len(seq_labs) == 1 and len(seq_f[next(iter(seq_labs))]) > len(M))
        seq_union = set().union(*(seq_f[l] for l in seq_labs if l)) if any(seq_labs) else set()
        top = sorted(M, key=lambda g: (info[g]["biotype"] != "protein_coding", info[g]["name"]))
        out.append(dict(
            family=f, cls=cls, miss_cause=miss_cause,
            extra_causes=dict(sorted({c: sum(1 for v in extra_cause.values() if v == c)
                                      for c in set(extra_cause.values())}.items())),
            n_members=len(fam_all[f]), n_clean=len(M), n_ours_pieces=len({l for l in labs if l}) + sum(1 for g in M if not pred_ours.get(g)),
            n_extra=len(extra), narrower_than_homology=narrower, seq_size=len(seq_union),
            manual=f in manual, exact_with_s1c_cn=frozenset(M) in {frozenset(v) for v in s1c_f.values()},
            example=",".join(info[g]["name"] for g in top[:3]),
            clean=sorted(M), extra=sorted(extra), extra_cause=extra_cause, ours_labels=sorted(l or "-" for l in labs),
            seq_members=sorted(seq_union)))
    counts = defaultdict(int)
    for r in out:
        counts[r["cls"]] += 1
    assert counts["match"] == anchor["n_exact"], (counts["match"], anchor["n_exact"])

    with open(a.out_tsv, "w") as fh:
        fh.write("family_id\tclass\tmiss_cause\textra_causes\tn_members\tn_clean\tn_our_pieces\tn_extra\t"
                 "narrower_than_homology\tsequence_family_size\tmanual_merge\texact_with_s1c_cn\texample_genes\n")
        for r in out:
            ec = ";".join(f"{k}={v}" for k, v in r["extra_causes"].items())
            fh.write(f"{r['family']}\t{r['cls']}\t{r['miss_cause']}\t{ec}\t{r['n_members']}\t{r['n_clean']}\t"
                     f"{r['n_ours_pieces']}\t{r['n_extra']}\t{int(r['narrower_than_homology'])}\t{r['seq_size']}\t"
                     f"{int(r['manual'])}\t{int(r['exact_with_s1c_cn'])}\t{r['example']}\n")

    # JSON for the viewer: every S1C gene (coordinates, exons, labels, copy numbers), the links, every family row
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
        ex = sorted(exons[g], key=lambda t: t[1])
        merged = []
        for c, s, e in ex:
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        gene_rows.append(dict(
            id=g, n=info[g]["name"], bt=info[g]["biotype"], c=ex[0][0] if ex else "", s=merged[0][0] if merged else 0,
            e=max(e for _, e in merged) if merged else 0, st=strand.get(g, "."), x=merged,
            sf=info[g]["fams"], seq=pred_seq.get(g, ""), our=pred_ours.get(g, ""), s1=pred_s1c.get(g, ""),
            cs=round(cn_s1c[g], 1) if g in cn_s1c else None, co=round(cn_ours[g], 1) if g in cn_ours else None))
    link_set = set()
    for x, y in edges:
        if x in idx and y in idx and x != y:
            link_set.add(tuple(sorted((idx[x], idx[y]))))
    sedef_only = set()
    for x, y in sedef:
        if x in idx and y in idx and x != y:
            k = tuple(sorted((idx[x], idx[y])))
            if k not in link_set:
                sedef_only.add(k)
    fam_rows = []
    for r in out:
        show = set(fam_all[r["family"]]) | set(r["extra"]) | set(r["seq_members"])
        for g in list(r["clean"]):
            l = pred_ours.get(g)
            if l:
                show |= ours_f[l]
        fam_rows.append(dict(
            f=r["family"], k=r["cls"], mc=r["miss_cause"], ec=r["extra_causes"], nm=r["n_members"], nc=r["n_clean"],
            np=r["n_ours_pieces"], nx=r["n_extra"], nh=r["narrower_than_homology"], ns=r["seq_size"], mm=r["manual"],
            x1=r["exact_with_s1c_cn"], ex=r["example"], g=sorted(idx[g] for g in show),
            xc={str(idx[g]): c for g, c in r["extra_cause"].items()}))
    blob = dict(anchor=dict(ari=round(anchor["ari"], 4), exact=anchor["n_exact"], ari_s1c=round(anchor_s1c["ari"], 4),
                            exact_s1c=anchor_s1c["n_exact"]),
                genes=gene_rows, links=sorted(link_set), sedef_links=sorted(sedef_only), families=fam_rows)
    with open(a.out_json, "w") as fh:
        json.dump(blob, fh, separators=(",", ":"))

    print(f"families {len(out)}: " + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    mc = defaultdict(int)
    for r in out:
        if r["miss_cause"]:
            mc[(r["cls"], r["miss_cause"])] += 1
    print("miss causes: " + ", ".join(f"{k[0]}/{k[1]} {v}" for k, v in sorted(mc.items())))
    ec = defaultdict(int)
    for r in out:
        if r["cls"] in ("soto_smaller", "mixed"):
            for k in r["extra_causes"]:
                ec[(r["cls"], k)] += 1
    print("extra-gene causes (families with >= 1): " + ", ".join(f"{k[0]}/{k[1]} {v}" for k, v in sorted(ec.items())))
    print(f"narrower than homology: {sum(r['narrower_than_homology'] for r in out)}/{len(out)}; "
          f"of the matches: {sum(r['narrower_than_homology'] for r in out if r['cls'] == 'match')}")


if __name__ == "__main__":
    main()

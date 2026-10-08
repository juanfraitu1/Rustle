#!/usr/bin/env python3
"""CAT/Liftoff v2.0 labels for the meeting page's Detection tab (2026-10-01; docs/archive/2026-10/CAT_RERUN_PROTOCOL_2026-10-01.md).

The July detection page (page/july/detection_2026-07-28.html) labels each of Soto's 362 members against NCBI RefSeq and suggests
excluding three classes as segmental-duplication pieces: unannotated, fragment (< 800 bp), piece of another gene. The script that
made those labels (`artifacts/member_annotation.json`, July 2026) was never committed. This script rebuilds the classifier from the
rules the page states and applies the SAME classifier to RefSeq and to CAT, so the two annotations are compared by one scorer:

  fragment       member span < 800 bp
  unannotated    no annotated gene overlaps the member (either strand)
  own gene       a gene of the member's own identity covers >= 75% of it -> that gene's biotype
  piece          otherwise, >= 75% of the member lies inside ONE other gene that is >= 2x the member's length (either strand)
  otherwise      the biotype of the gene overlapping the member most

"Other gene" means a different name for RefSeq (the registered convention) and a different gene id for CAT (R5: Soto's genes are CAT
genes, joined by Gene ID). The 0.75 / 0.75 / 2x constants are the ones that best reproduce the registered RefSeq labels (335 of 362 agree,
and the class counts match: piece 47 = 47, unannotated 9 vs 11, fragment 17 vs 15); they are not tuned on CAT. A member nested in a
longer gene of another identity (the piece test, run even when the own gene is present) is reported beside its CAT label as
information, never as a class. Every Soto member is exactly a CAT v2.0 gene span (checked), so CAT "unannotated" is 0 by construction.

Extra copies (DNA mode, `DEXTRA`): the page flags a copy whose length is > 2x or < 0.5x its family's median member length. That flag
uses no annotation. Each copy is given the CAT genes it overlaps and one more label: a "family-sized" CAT gene is one whose span is
0.5-2x the family median and that overlaps the copy by >= 50% of the shorter of the two. A copy > 2x with such a gene is a block
around a family-sized gene; < 0.5x inside such a gene is a piece of one. Gene content is not homology: the gene is named, not proven.

    python3 bench/soto_m2/soto_m2_cat_labels.py --detection bench/soto_m2/page/july/detection_2026-07-28.html \\
        --cat-genes chm13v2.0_CAT_Liftoff.genes.tsv --s1c bench/soto/soto_famCN_S1C.tsv \\
        --refseq chm13v2.0_RefSeq_full.gff.gz --out cat_labels.json
"""
import argparse
import collections
import csv
import gzip
import json
import re

FRAG_BP, OWN_COV, PIECE_IN, PIECE_X = 800, 0.75, 0.75, 2.0
EXCL = {"unannotated", "fragment", "piece-of-other"}


def page_data(path):
    t = open(path).read()
    dec = json.JSONDecoder()
    D = dec.raw_decode(t[t.index("{", t.index("const D = ")):])[0]
    X = dec.raw_decode(t[t.index("{", t.index("const DEXTRA=")):])[0]
    return D, X


def cls_of(bt):
    if bt == "protein_coding":
        return "coding"
    if "pseudogene" in bt:
        return "transcribed_pseudogene" if bt.startswith("transcribed") else "pseudogene"
    return "lncRNA" if bt == "lncRNA" else "other"


def classify(chrom, s, e, genes, is_self):
    """(class, annotation gene tuple or None, container tuple or None). genes[chrom] = [(start0, end, label, biotype, id), ...]."""
    bp = e - s
    ov = [(min(e, ge) - max(s, gs), gs, ge, lab, bt, gid) for gs, ge, lab, bt, gid in genes.get(chrom, ()) if gs < e and s < ge]
    own = [t for t in ov if is_self(t) and t[0] >= OWN_COV * bp]
    pieces = [t for t in ov if not is_self(t) and t[0] >= PIECE_IN * bp and t[2] - t[1] >= PIECE_X * bp]
    box = max(pieces, key=lambda t: t[2] - t[1]) if pieces else None
    best = max(own or ov, key=lambda t: t[0]) if ov else None
    if bp < FRAG_BP:
        return "fragment", best, box
    if not ov:
        return "unannotated", None, None
    if own:
        return cls_of(best[4]), best, box
    if box:
        return "piece-of-other", box, box
    return cls_of(best[4]), best, box


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("detection", "cat_genes", "s1c", "refseq", "out"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    D, DEXTRA = page_data(a.detection)

    cat, catg = {}, collections.defaultdict(list)
    for r in csv.DictReader(open(a.cat_genes), delimiter="\t"):
        cat[r["gene_id"]] = r
        catg[r["chrom"]].append((int(r["start0"]), int(r["end"]), r["gene_name"], r["gene_biotype"], r["gene_id"]))
    s1c_ids, s1c_fam = collections.defaultdict(list), {}
    for r in csv.DictReader(open(a.s1c), delimiter="\t"):
        s1c_ids[r["Gene Name"]].append(r["Gene ID"])
        s1c_fam[r["Gene ID"]] = r["Family ID"]
    rs = collections.defaultdict(list)
    for line in gzip.open(a.refseq, "rt"):
        if line[0] == "#":
            continue
        f = line.split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene"):
            continue
        at = dict(x.split("=", 1) for x in f[8].rstrip().split(";") if "=" in x)
        rs[f[0]].append((int(f[3]) - 1, int(f[4]), at.get("Name", ""), at.get("gene_biotype", ""), at.get("ID", "")))

    members, agree, cnt = {}, 0, collections.defaultdict(collections.Counter)
    for fam in D["fams"]:
        for m in fam["members"]:
            c, s, e = m["chrom"], m["start"], m["end"]
            own = [cat[i] for i in s1c_ids.get(m["gene"], ()) if i in cat and cat[i]["chrom"] == c
                   and int(cat[i]["start0"]) == s + 1 and int(cat[i]["end"]) == e]
            assert len(own) == 1, (m["gene"], c, s, e, len(own))
            gid = own[0]["gene_id"]
            k_rs, g_rs, _ = classify(c, s, e, rs, lambda t, n=m["gene"]: t[3] == n)
            k_cat, g_cat, box = classify(c, s, e, catg, lambda t, i=gid: t[5] == i)
            agree += k_rs == m["cls"]
            for lab, k in (("registered", m["cls"]), ("refseq_rebuilt", k_rs), ("cat", k_cat)):
                cnt[lab][k] += 1
            if k_cat == "piece-of-other":
                ann = f"{g_cat[3]} ({g_cat[4]}, {(g_cat[2] - g_cat[1]) / 1000:,.1f} kb)"
            elif k_cat == "fragment":
                ann = f"{own[0]['gene_name']} ({own[0]['gene_biotype']}, {e - s} bp)"
            else:
                ann = f"{own[0]['gene_name']} ({own[0]['gene_biotype']})"
            nest = f"{box[3]} ({box[4]}, {(box[2] - box[1]) / 1000:,.1f} kb)" if box and k_cat != "piece-of-other" else ""
            cnt["cat_nested"][bool(nest)] += 1
            members[f"{fam['fam']}:{s}"] = dict(cls=k_cat, ann=ann, excl=k_cat in EXCL, gid=gid, nested=nest, cls_rs=m["cls"],
                                                 ann_rs=m["ann"], cls_rs_rebuilt=k_rs)
    n = len(members)

    fammed = {}
    for fam in D["fams"]:
        b = sorted(m["bp"] for m in fam["members"])
        fammed[fam["fam"]] = b[len(b) // 2] if b else 0
    extras, ecnt = {}, collections.Counter()
    for fam, xs in DEXTRA.items():
        med = fammed.get(fam, 0)
        out = []
        for x in xs:
            c, s, e = x["chrom"], x["start"], x["end"]
            bp = e - s
            r = bp / med if med else 1.0
            raw = "big" if r > 2 else "small" if r < 0.5 else "ok"
            ov = sorted(((min(e, ge) - max(s, gs), gs, ge, lab, bt, gid) for gs, ge, lab, bt, gid in catg.get(c, ()) if gs < e and s < ge),
                        key=lambda t: -t[0])
            fit = [t for t in ov if med and 0.5 <= (t[2] - t[1]) / med <= 2 and t[0] >= 0.5 * min(t[2] - t[1], bp)]
            fit = max(fit, key=lambda t: t[0]) if fit else None
            status = ("nogene" if not ov else "ok" if raw == "ok" else ("block" if raw == "big" else "piece") if fit else "nofit")
            ecnt[(raw, status)] += 1
            gl = lambda t: dict(name=t[3], id=t[5], bt=t[4], bp=t[2] - t[1], ov=t[0], soto=s1c_fam.get(t[5], ""))
            out.append(dict(status=status, raw=raw, genes=[gl(t) for t in ov[:3]], n_genes=len(ov), fit=gl(fit) if fit else None,
                            fit_ratio=round((fit[2] - fit[1]) / med, 2) if fit else None))
        extras[fam] = out

    summ = dict(n_members=n, rebuilt_agreement=agree, counts={k: dict(v) for k, v in cnt.items()},
                excluded={k: sum(v[c] for c in EXCL) for k, v in cnt.items() if k != "cat_nested"},
                extras={f"{a}|{b}": v for (a, b), v in sorted(ecnt.items())}, n_extras=sum(ecnt.values()),
                rules=dict(fragment_bp=FRAG_BP, piece_inside=PIECE_IN, piece_times=PIECE_X))
    json.dump(dict(members=members, extras=extras, summary=summ), open(a.out, "w"), separators=(",", ":"))
    print(f"members {n}; rebuilt RefSeq classifier agrees with the registered label on {agree} ({agree / n:.1%})")
    for k in ("registered", "refseq_rebuilt", "cat"):
        print(f"  {k:15} suggested exclusions {summ['excluded'][k]:3}  {dict(sorted(cnt[k].items()))}")
    print(f"  CAT members nested in a >= 2x longer CAT gene of another id (information only): {cnt['cat_nested'][True]}")
    print(f"extra copies {summ['n_extras']}: {summ['extras']}")


if __name__ == "__main__":
    main()

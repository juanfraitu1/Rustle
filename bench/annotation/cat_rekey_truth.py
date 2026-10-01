#!/usr/bin/env python3
"""Step 1 of docs/CAT_RERUN_PROTOCOL_2026-10-01.md: re-key the NPIP and TBC1D3 truth from RefSeq to CAT/Liftoff v2.0 (rule R1/R2).

NPIP: the 27-row Dishuck-checked table (`docs/lit_subclusters_npip_dishuck_check.tsv`, RefSeq GeneID + coordinates). TBC1D3: the RefSeq
genes whose `description` is "TBC1 domain family member 3..." (the copy-recovery member rule), chr17. Each RefSeq copy goes to the CAT
gene sharing the most exonic bases on the same strand (ties: Jaccard); a RefSeq record without exon lines goes by span overlap; no CAT
gene -> dropped and listed. Other CAT records overlapping a copy are listed, never counted as copies. Output: one TSV per family with
the original columns plus the CAT columns, and a summary.

    python3 bench/annotation/cat_rekey_truth.py --refseq chm13v2.0_RefSeq_full.gff.gz --cat-genes chm13v2.0_CAT_Liftoff.genes.tsv \
        --npip docs/lit_subclusters_npip_dishuck_check.tsv --out-dir docs
"""
import argparse
import collections
import csv
import gzip
import os
import re


def attrs(col):
    return dict(x.split("=", 1) for x in col.rstrip().split(";") if "=" in x)


def merge(ivs):
    out = []
    for a, b in sorted(ivs):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def ov_bp(a, b):
    i = j = t = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if lo < hi:
            t += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return t


def blen(bl):
    return sum(e - s for s, e in bl)


def load_refseq(path, chroms):
    genes, tx_gene, ex = {}, {}, collections.defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in chroms:
                continue
            a = attrs(f[8])
            i = a.get("ID", "")
            if i.startswith("gene-"):
                gid = re.search(r"GeneID:(\d+)", a.get("Dbxref", ""))
                genes[i] = dict(id=i, name=a.get("Name", ""), geneid=f"GeneID:{gid.group(1)}" if gid else "", chrom=f[0],
                                start0=int(f[3]) - 1, end=int(f[4]), strand=f[6], desc=a.get("description", ""),
                                biotype=a.get("gene_biotype", ""))
            elif f[2] == "exon":
                ex[a.get("Parent", "")].append((int(f[3]) - 1, int(f[4])))
            elif a.get("Parent", "").startswith("gene-") and i:
                tx_gene[i] = a["Parent"]
    gex = collections.defaultdict(list)
    for p, v in ex.items():
        g = p if p.startswith("gene-") else tx_gene.get(p)
        if g in genes:
            gex[g].extend(v)
    for g in genes.values():
        g["blocks"] = merge(gex.get(g["id"], []))
    return genes


def load_cat(path, chroms):
    out = []
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["chrom"] in chroms:
            bl = [tuple(map(int, x.split("-"))) for x in r["exon_blocks"].split(",") if x]
            out.append(dict(id=r["gene_id"], name=r["gene_name"], biotype=r["gene_biotype"], source=r["source"], chrom=r["chrom"],
                            start0=int(r["start0"]), end=int(r["end"]), strand=r["strand"], blocks=bl))
    return out


def rekey(rg, cat):
    """(best CAT gene, rule, shared, jaccard, quality, others) for one RefSeq gene."""
    same = [c for c in cat if c["chrom"] == rg["chrom"] and c["strand"] == rg["strand"] and c["start0"] < rg["end"] and rg["start0"] < c["end"]]
    others = [c for c in cat if c["chrom"] == rg["chrom"] and c["start0"] < rg["end"] and rg["start0"] < c["end"]]
    best, rule = None, "none"
    if rg["blocks"]:
        cands = []
        for c in same:
            sh = ov_bp(rg["blocks"], c["blocks"])
            if sh:
                un = blen(rg["blocks"]) + blen(c["blocks"]) - sh
                cands.append((sh, sh / un, c))
        if cands:
            sh, jac, best = max(cands, key=lambda t: (t[0], t[1]))
            rule = "exons"
    else:
        # Amendment 1 (protocol): span Jaccard, so a read-through record that merely contains the copy does not win
        cands = [(min(rg["end"], c["end"]) - max(rg["start0"], c["start0"]),
                  max(rg["end"], c["end"]) - min(rg["start0"], c["start0"]), c) for c in same]
        cands = [t for t in cands if t[0] > 0 and t[2]["blocks"]]
        if cands:
            sh, un, best = max(cands, key=lambda t: t[0] / t[1])
            jac, rule = sh / un, "span"
    if best is None:
        return None, "none", 0, 0.0, "none", others
    if rule == "exons":
        q = ("strong" if sh >= 0.5 * blen(rg["blocks"]) and sh >= 0.5 * blen(best["blocks"]) else
             "partial" if sh >= 0.5 * blen(rg["blocks"]) or sh >= 0.5 * blen(best["blocks"]) else "weak")
    else:
        q = "span"
    return best, rule, sh, jac, q, [c for c in others if c is not best]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("refseq", "cat_genes", "npip", "out_dir"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    chroms = {"chr16", "chr17", "chr18"}
    rs = load_refseq(a.refseq, chroms)
    cat = load_cat(a.cat_genes, chroms)
    by_geneid = {g["geneid"]: g for g in rs.values() if g["geneid"]}

    extra_cols = ["refseq_id", "cat_gene_id", "cat_name", "cat_biotype", "cat_source", "cat_start0", "cat_end", "rule", "shared_bp",
                  "jaccard", "quality", "other_cat_records"]

    def out_row(rg):
        best, rule, sh, jac, q, others = rekey(rg, cat)
        oth = ";".join(f"{c['id']}:{c['name']}" for c in others)
        if best is None:
            return [rg["id"], "", "", "", "", "", "", rule, 0, "0", q, oth], None
        return [rg["id"], best["id"], best["name"], best["biotype"], best["source"], best["start0"], best["end"], rule, sh,
                f"{jac:.3f}", q, oth], best

    summary = []
    # NPIP: the Dishuck-checked table
    rows = list(csv.DictReader(open(a.npip), delimiter="\t"))
    cols = list(rows[0].keys())
    path = os.path.join(a.out_dir, "lit_subclusters_npip_dishuck_check.CAT.tsv")
    q = collections.Counter()
    used = collections.Counter()
    with open(path, "w") as out:
        out.write("\t".join(cols + extra_cols) + "\n")
        for r in rows:
            rg = by_geneid.get(r["gene_id"])
            assert rg is not None and rg["chrom"] == r["chrom"] and rg["start0"] + 1 == int(r["start"]), r["gene_id"]
            ext, best = out_row(rg)
            q[ext[10]] += 1
            if best:
                used[best["id"]] += 1
            out.write("\t".join([r[c] for c in cols] + [str(x) for x in ext]) + "\n")
    shared_cat = [g for g, n in used.items() if n > 1]
    summary.append(f"NPIP (Dishuck table, {len(rows)} rows): " + ", ".join(f"{k} {v}" for k, v in q.most_common())
                   + f"; CAT genes used by more than one copy: {len(shared_cat)} {shared_cat}")

    # TBC1D3: RefSeq description rule (copy-recovery member rule), chr17
    tb = sorted((g for g in rs.values() if g["chrom"] == "chr17" and g["desc"].startswith("TBC1 domain family member 3")
                 and "-" not in g["name"]), key=lambda g: g["start0"])
    path2 = os.path.join(a.out_dir, "tbc1d3_members.CAT.tsv")
    q2 = collections.Counter()
    used2 = collections.Counter()
    with open(path2, "w") as out:
        out.write("\t".join(["refseq_name", "geneid", "chrom", "start0", "end", "strand", "refseq_biotype", "refseq_has_exons"]
                            + extra_cols) + "\n")
        for g in tb:
            ext, best = out_row(g)
            q2[ext[10]] += 1
            if best:
                used2[best["id"]] += 1
            out.write("\t".join([g["name"], g["geneid"], g["chrom"], str(g["start0"]), str(g["end"]), g["strand"], g["biotype"],
                                 str(int(bool(g["blocks"])))] + [str(x) for x in ext]) + "\n")
    cat_tb = sorted(c["name"] for c in cat if c["name"].startswith("TBC1D3") and "-" not in c["name"])
    summary.append(f"TBC1D3 (RefSeq description rule, {len(tb)} genes): " + ", ".join(f"{k} {v}" for k, v in q2.most_common())
                   + f"; CAT genes used by more than one copy: {[g for g, n in used2.items() if n > 1]}; CAT genes named TBC1D3* "
                   f"(reported beside, R2): {len(cat_tb)}")
    print("\n".join(summary))


if __name__ == "__main__":
    main()

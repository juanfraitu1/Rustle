#!/usr/bin/env python3
"""Score the three read-pool arms of docs/PREREG_npip_read_pool_2026-10-01.md on NPIP (human chr16, A119b): locus definition at the
25 CAT NPIP copies, echo loci (zero primary records over the span), annotation class, the NPIP family after the all-vs-all + MCL, and
the advisor's rule f = (ALL-added echo or off-copy loci that end up in the NPIP family) / (such loci).

    python3 bench/npip_read_pool/score.py --dir /mnt/linuxdisk/tmp/readpool_npip --copies copies.hsa.tsv --cat-genes genes.tsv \
        --bam A119b.t2t.bam --out npip_read_pool.json
"""
import argparse
import collections
import csv
import json
import subprocess

ARMS = ("P", "GOOD", "ALL")


def merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def ov(a, b):
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


def load_loci(prefix):
    """gene id -> dict(chrom, s0, e, strand, exons [[s0, e]]) from the arm's locus GFF3 (rep exons, span = transcript hull)."""
    loci = {}
    for ln in open(prefix + ".gff3"):
        if ln[0] == "#":
            continue
        f = ln.rstrip("\n").split("\t")
        at = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "gene":
            loci[at["Name"]] = dict(chrom=f[0], s0=int(f[3]) - 1, e=int(f[4]), strand=f[6], exons=[])
        elif f[2] == "exon":
            loci[at["gene"]]["exons"].append([int(f[3]) - 1, int(f[4])])
    for v in loci.values():
        v["exons"] = merge(v["exons"])
    return loci


def load_clusters(prefix, loci):
    """locus name -> cluster id. mcl_families' clusters.tsv lists multi-member clusters by member coordinates (chrom, GFF start,
    end); fam.loci.tsv lists records folded into a representative (chrom:start-end -> chrom:start-end), which take its cluster."""
    by = {(L["chrom"], L["s0"] + 1, L["e"]): n for n, L in loci.items()}
    out = {}
    for r in csv.DictReader(open(prefix + ".fam.clusters.tsv"), delimiter="\t"):
        key = (r["chrom"], int(r["start"]), int(r["end"]))
        assert key in by, key
        out[by[key]] = r["cluster_id"]
    def parse(x):
        c, se = x.rsplit(":", 1)
        a_, b_ = se.split("-")
        return (c, int(a_), int(b_))
    keycl = {(loci[n]["chrom"], loci[n]["s0"] + 1, loci[n]["e"]): cl for n, cl in out.items()}
    for r in csv.DictReader(open(prefix + ".fam.loci.tsv"), delimiter="\t"):
        a_, rep = parse(r["annotation"]), parse(r["representative"])
        if a_ in by and rep in keycl and by[a_] not in out:
            out[by[a_]] = keycl[rep]
    return out


def primary_counts(bam, loci, chrom):
    """locus -> number of primary records (-F 2308) overlapping its span (one samtools call per locus, region list)."""
    names = sorted(loci, key=lambda n: loci[n]["s0"])
    out = {}
    for n in names:
        L = loci[n]
        r = subprocess.run(["samtools", "view", "-c", "-F", "2308", bam, f"{chrom}:{L['s0'] + 1}-{L['e']}"],
                           capture_output=True, text=True, check=True)
        out[n] = int(r.stdout.strip())
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("dir", "copies", "cat_genes", "bam", "out"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    chrom = "chr16"
    copies = [c for c in csv.DictReader(open(a.copies), delimiter="\t") if c["family"] == "NPIP" and c["chrom"] == chrom]
    cat, catg = {}, []
    for r in csv.DictReader(open(a.cat_genes), delimiter="\t"):
        if r["chrom"] == chrom:
            bl = merge([[int(x) for x in b.split("-")] for b in r["exon_blocks"].split(",") if b])
            cat[r["gene_id"]] = r
            catg.append((r["gene_id"], r["gene_name"], r["strand"], bl))
    cp = [dict(cid=c["cid"], name=c["refseq_name"], cat=c["isoform_gene"], strand=c["strand"],
               exons=merge([[int(x) for x in b.split("-")] for b in cat[c["isoform_gene"]]["exon_blocks"].split(",") if b]))
          for c in copies]
    copy_ids = {c["cat"] for c in cp}
    lo = min(c["exons"][0][0] for c in cp) - 2_000_000
    hi = max(c["exons"][-1][1] for c in cp) + 2_000_000

    res, arms = {}, {}
    for arm in ARMS:
        loci = load_loci(f"{a.dir}/{arm}")
        clus = load_clusters(f"{a.dir}/{arm}", loci)
        on = {}
        for n, L in loci.items():
            on[n] = [c["cid"] for c in cp if c["strand"] == L["strand"] and ov(L["exons"], c["exons"]) > 0]
        # NPIP family = the cluster with the most copy loci
        cc = collections.Counter(clus[n] for n in loci if on[n] and n in clus)
        fam = cc.most_common(1)[0][0] if cc else None
        famset = {n for n in loci if clus.get(n) == fam}
        # loci to characterise: every NPIP-family member and every locus on a copy
        focus = {n for n in loci if on[n]} | famset
        prim = primary_counts(a.bam, {n: loci[n] for n in focus}, chrom)
        rows = {}
        for n in focus:
            L = loci[n]
            genes = [(gid, nm) for gid, nm, st, bl in catg if st == L["strand"] and ov(L["exons"], bl) > 0]
            klass = "npip_copy" if on[n] else ("other_gene" if genes else "no_gene")
            rows[n] = dict(chrom=L["chrom"], s0=L["s0"], e=L["e"], strand=L["strand"], exons=L["exons"], copies=on[n],
                           in_fam=n in famset, cluster=clus.get(n), primaries=prim[n], echo=prim[n] == 0, klass=klass,
                           genes=[nm for _, nm in genes if _ not in copy_ids][:3])
        cov = collections.Counter(c for n in loci for c in on[n])
        arms[arm] = dict(loci_chr16=len(loci), loci_on_copies=sum(1 for n in loci if on[n]),
                         copies_covered=sum(1 for c in cp if cov[c["cid"]]), loci_per_copy=dict(cov),
                         fused=sum(1 for n in loci if len(on[n]) >= 2), echo_on_copies=sum(1 for n in loci if on[n] and rows[n]["echo"]),
                         fam_size=len(famset), fam_by_class=dict(collections.Counter(rows[n]["klass"] for n in famset)),
                         fam_echo=sum(1 for n in famset if rows[n]["echo"]), paf_records=sum(1 for _ in open(f"{a.dir}/{arm}.paf")),
                         clusters=len(set(clus.values())))
        res[arm] = dict(loci=loci, rows=rows)

    # ALL-added loci: rep exons overlap no GOOD locus's rep exons on the same strand (chr16-wide)
    good = res["GOOD"]["loci"]
    gidx = collections.defaultdict(list)
    for n, L in good.items():
        gidx[L["strand"]].append(L)
    added = []
    for n, L in res["ALL"]["loci"].items():
        if not any(g["s0"] < L["e"] and L["s0"] < g["e"] and ov(L["exons"], g["exons"]) > 0 for g in gidx[L["strand"]]):
            added.append(n)
    rowsA = res["ALL"]["rows"]
    # primaries for added loci outside the focus set (needed for the echo status of every added locus)
    extra = {n: res["ALL"]["loci"][n] for n in added if n not in rowsA}
    pe = primary_counts(a.bam, extra, chrom) if extra else {}
    def is_fp(n):
        if n in rowsA:
            return rowsA[n]["echo"] or rowsA[n]["klass"] != "npip_copy"
        return True    # outside the focus set = not on a copy -> off-copy by definition
    fp = [n for n in added if is_fp(n)]
    fp_in_fam = [n for n in fp if n in rowsA and rowsA[n]["in_fam"]]
    f = len(fp_in_fam) / len(fp) if fp else 0.0
    verdict = "claim holds" if f <= 0.10 else "claim fails" if f >= 0.50 else "partial"
    # added loci near NPIP (for the figure), and the family's added loci
    near = [n for n in added if lo <= res["ALL"]["loci"][n]["s0"] <= hi]
    summary = dict(arms=arms, all_added=len(added), all_added_fp=len(fp), all_added_fp_in_npip_family=len(fp_in_fam), f=round(f, 4),
                   verdict=verdict, all_added_in_family=sum(1 for n in added if n in rowsA and rowsA[n]["in_fam"]),
                   all_added_echo=sum(1 for n in added if (rowsA[n]["echo"] if n in rowsA else pe[n] == 0)))
    copies_out = [dict(cid=c["cid"], name=c["name"], cat=c["cat"], strand=c["strand"], s0=c["exons"][0][0], e=c["exons"][-1][1],
                       exons=c["exons"]) for c in cp]
    loci_out = {arm: {n: r for n, r in res[arm]["rows"].items()} for arm in ARMS}
    json.dump(dict(summary=summary, copies=copies_out, loci=loci_out, near_added=near), open(a.out, "w"), separators=(",", ":"))
    for arm in ARMS:
        print(arm, json.dumps({k: v for k, v in arms[arm].items() if k != "loci_per_copy"}))
    print({k: v for k, v in summary.items() if k != "arms"})


if __name__ == "__main__":
    main()

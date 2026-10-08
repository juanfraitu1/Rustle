#!/usr/bin/env python3
"""Chain-level comparison of the reads at annotated copies with their annotated transcripts (2026-10-04, after Amendment B of
docs/archive/2026-10/PREREG_spliced_copy_support_2026-10-04.md): is the annotation what is expressed?

    npip_chains.py --copies copies.tsv --truth truth.gtf [--copies2 copies2.tsv --truth2 truth2.gtf] --bam reads.bam --family NPIP --out PREFIX

Per copy: every same-strand primary read overlapping the exon union (CAT ∪ RefSeq models) gives its in-span junction chain (N >= 50 bp, exact
coordinates). Each distinct chain is classified against each annotation's models:
  FSM  the chain equals a model's intron chain
  ISM  the chain is a contiguous sub-chain of a model's intron chain (the Amendment B support class, with FSM)
  NIC  every donor and acceptor of the chain is an annotated splice site of the copy, but the chain is not a sub-chain of any model
  NNC  at least one splice site is not annotated
  1J   one junction (not enough to be a chain; classified the same way, reported apart); U unspliced.
Reported per copy: reads by class (vs annotation 1, vs annotation 2, vs the union), the top chains by read count with their class and
whether each novel junction recurs (>= 3 reads) or is a singleton, and the dominant chain's class. Writes PREFIX.copies.tsv (one row per
copy), PREFIX.chains.tsv (one row per distinct chain with >= 2 reads) and PREFIX.json.
"""
import argparse
import collections
import csv
import json
import os
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from copy_support import gtf_transcripts, introns_of, merge, inter, read_blocks_junctions  # noqa: E402

MIN_SUPPORT = 3


def classify(chain, models):
    """class of a junction chain against a list of model intron chains (each in order)."""
    if not models:
        return "noann"
    chain = tuple(chain)
    for m in models:
        if chain == tuple(m):
            return "FSM"
    for m in models:
        n, L = len(m), len(chain)
        if L < n and any(tuple(m[i:i + L]) == chain for i in range(n - L + 1)):
            return "ISM"
    sites = {s for m in models for j in m for s in j}
    if all(j[0] in sites and j[1] in sites for j in chain):
        return "NIC"
    return "NNC"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--copies2")
    ap.add_argument("--truth2")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=5)
    a = ap.parse_args(argv)
    copies = [r for r in csv.DictReader(open(a.copies), delimiter="\t") if r["family"] == a.family]
    tx = gtf_transcripts(a.truth, {r["isoform_gene"] for r in copies} | {r["cid"] for r in copies})
    tx2 = {}
    if a.truth2 and a.copies2:
        c2 = list(csv.DictReader(open(a.copies2), delimiter="\t"))
        by_name = {r["name"]: r["cid"] for r in c2}
        t2 = gtf_transcripts(a.truth2, {r["cid"] for r in c2})
        for c in copies:
            nm = c.get("refseq_name") or c["name"]
            if nm in by_name and by_name[nm] in t2:
                tx2[c["cid"]] = t2[by_name[nm]]
    bam = pysam.AlignmentFile(a.bam)
    rows, chain_rows, summary = [], [], {}
    for c in copies:
        label = c.get("refseq_name") or c.get("cat_name") or c["name"]
        ex1 = tx.get(c["cid"]) or tx.get(c["isoform_gene"], {})
        ex2 = tx2.get(c["cid"], {})
        m1 = [introns_of(e) for e in ex1.values() if len(introns_of(e)) >= 1]
        m2 = [introns_of(e) for e in ex2.values() if len(introns_of(e)) >= 1]
        union = merge([e for d in (ex1, ex2) for exs in d.values() for e in exs])
        if not union:
            continue
        lo, hi = union[0][0], union[-1][1]
        chains = collections.Counter()
        n_reads = n_uns = 0
        jcount = collections.Counter()
        per_chain_reads = collections.defaultdict(list)   # chain -> [(mapq, de)] for the tie / divergence split by class
        for rd in bam.fetch(c["chrom"], lo, hi):
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary or ("-" if rd.is_reverse else "+") != c["strand"]:
                continue
            blocks, js = read_blocks_junctions(rd)
            if inter(blocks, union) == 0:
                continue
            js = tuple(j for j in js if j[0] >= lo and j[1] <= hi)
            n_reads += 1
            if not js:
                n_uns += 1
                continue
            chains[js] += 1
            per_chain_reads[js].append((rd.mapping_quality, rd.get_tag("de") if rd.has_tag("de") else float("nan")))
            for j in js:
                jcount[j] += 1
        cls = {}
        by = {"ann1": collections.Counter(), "ann2": collections.Counter(), "union": collections.Counter()}
        for ch, n in chains.items():
            k1, k2, ku = classify(ch, m1), classify(ch, m2), classify(ch, m1 + m2)
            if len(ch) == 1:
                k1, k2, ku = "1J:" + k1, "1J:" + k2, "1J:" + ku
            cls[ch] = (k1, k2, ku)
            by["ann1"][k1] += n
            by["ann2"][k2] += n
            by["union"][ku] += n
        # tie status and divergence by class (vs the union): are the NNC reads placed here by a coin toss (MAPQ 0, higher de)?
        import statistics
        cls_reads = collections.defaultdict(list)
        for ch, lst in per_chain_reads.items():
            k = cls[ch][2].split(":")[-1]
            cls_reads[k].extend(lst)
        tie = {}
        for k in ("FSM", "ISM", "NIC", "NNC"):
            L = cls_reads.get(k, [])
            tie[k] = (round(sum(1 for q, _ in L if q == 0) / len(L), 3) if L else None,
                      round(statistics.median([d for _, d in L if d == d]), 4) if L else None)
        sites_u = {s for m in m1 + m2 for j in m for s in j}
        top = []
        for ch, n in chains.most_common(a.top):
            novel = [j for j in ch if not (j[0] in sites_u and j[1] in sites_u)]
            recurrent = sum(1 for j in novel if jcount[j] >= MIN_SUPPORT)
            top.append(dict(reads=n, n_junctions=len(ch), cls_cat=cls[ch][0], cls_refseq=cls[ch][1], cls_union=cls[ch][2],
                            novel_sites_junctions=len(novel), novel_recurrent=recurrent, chain=";".join(f"{s}-{e}" for s, e in ch)))
            chain_rows.append(dict(copy=label, **top[-1]))
        dom = top[0] if top else None
        n_spliced = n_reads - n_uns
        support_u = by["union"]["FSM"] + by["union"]["ISM"]
        row = dict(cid=c["cid"], copy=label, chrom=c["chrom"], strand=c["strand"], span=f"{lo}-{hi}", models_cat=len(m1), models_refseq=len(m2),
                   cat_chain_lengths=",".join(str(len(m)) for m in m1), refseq_chain_lengths=",".join(str(len(m)) for m in m2),
                   reads=n_reads, unspliced=n_uns, distinct_chains=len(chains), chains_ge2=sum(1 for n in chains.values() if n >= 2),
                   support_union=support_u, support_frac=round(support_u / n_spliced, 3) if n_spliced else 0.0)
        for ann in ("ann1", "ann2", "union"):
            for k in ("FSM", "ISM", "NIC", "NNC"):
                row[f"{ann}_{k}"] = by[ann][k]
            row[f"{ann}_1J"] = sum(v for kk, v in by[ann].items() if kk.startswith("1J:"))
        for k in ("FSM", "ISM", "NIC", "NNC"):
            row[f"{k}_mapq0_frac"] = tie[k][0]
            row[f"{k}_median_de"] = tie[k][1]
        row.update(dominant_reads=dom["reads"] if dom else 0, dominant_junctions=dom["n_junctions"] if dom else 0,
                   dominant_cls_cat=dom["cls_cat"] if dom else "-", dominant_cls_refseq=dom["cls_refseq"] if dom else "-",
                   dominant_cls_union=dom["cls_union"] if dom else "-", dominant_novel_recurrent=dom["novel_recurrent"] if dom else 0)
        rows.append(row)
    with open(a.out + ".copies.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    with open(a.out + ".chains.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(chain_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(chain_rows)
    tot = collections.Counter()
    for r in rows:
        for k in ("FSM", "ISM", "NIC", "NNC", "1J"):
            tot[k] += r[f"union_{k}"]
        tot["unspliced"] += r["unspliced"]
        tot["reads"] += r["reads"]
    summary = dict(copies=len(rows), reads=tot["reads"], by_class_union={k: tot[k] for k in ("FSM", "ISM", "NIC", "NNC", "1J", "unspliced")},
                   copies_dominant_annotated=sum(1 for r in rows if r["dominant_cls_union"] in ("FSM", "ISM")),
                   copies_dominant_NIC=sum(1 for r in rows if r["dominant_cls_union"] == "NIC"),
                   copies_dominant_NNC=sum(1 for r in rows if r["dominant_cls_union"] == "NNC"))
    json.dump(summary, open(a.out + ".json", "w"), indent=1)
    print(json.dumps(summary))
    hdr = ["copy", "reads", "union_FSM", "union_ISM", "union_NIC", "union_NNC", "support_frac", "FSM_mapq0_frac", "ISM_mapq0_frac", "NIC_mapq0_frac", "NNC_mapq0_frac",
           "FSM_median_de", "ISM_median_de", "NNC_median_de", "dominant_reads", "dominant_cls_union", "dominant_novel_recurrent"]
    print("\t".join(hdr))
    for r in rows:
        print("\t".join(str(r[h]) for h in hdr))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Spliced support per annotated copy (docs/PREREG_spliced_copy_support_2026-10-04.md): is a copy FOUND because reads are spliced
transcripts of it, or only because something overlaps it?

    copy_support.py --copies copies.tsv --truth truth.gtf --bam reads.bam --family NPIP --out PREFIX
                    [--loci NAME=loci.gff3 ...] [--nodes pagedata.json]

Per copy (rows of --copies with the given family): the primary reads (-F 2308) overlapping its annotated exon union on its strand;
junctions = N ops >= 50 bp (exact donor/acceptor); supported junction = carried by >= 3 reads at the copy; k = min(2, annotated introns
>= 50 bp); structural-support read = >= k supported junctions (k = 0: aligned blocks cover >= 50% of the exon union); spliced-expressed
= >= 2 such reads. For each --loci set: OLD = a same-strand locus whose rep exons overlap the exon union; STRICT = spliced-expressed and a
same-strand locus whose rep junctions include >= k supported junctions (k = 0: rep exons cover >= 50% of the union). --nodes restricts
the loci of each arm to the page's NPIP-cluster nodes (pagedata.json rows: cid -> node[arm]) and reports found-within-NPIP-clusters too.
Writes PREFIX.copies.tsv and PREFIX.json.
"""
import argparse
import collections
import csv
import json
import re

import pysam

MIN_INTRON = 50
MIN_JUNCTION_READS = 3
FLOOR = 2
COVER = 0.5


def merge(iv):
    iv = sorted(iv)
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def ilen(iv):
    return sum(e - s for s, e in iv)


def inter(a, b):
    i = j = 0
    tot = 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            tot += e - s
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def gtf_transcripts(path, gene_ids):
    """gene_id -> {transcript_id: sorted exon list [[s0, e]]} for the given gene ids (GTF, 1-based closed)."""
    tx = collections.defaultdict(lambda: collections.defaultdict(list))
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        g = re.search(r'gene_id "([^"]+)"', f[8])
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        if not g or g.group(1) not in gene_ids:
            continue
        tx[g.group(1)][t.group(1) if t else g.group(1)].append([int(f[3]) - 1, int(f[4])])
    return {g: {t: sorted(ex) for t, ex in d.items()} for g, d in tx.items()}


def introns_of(exons):
    return [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1) if exons[i + 1][0] - exons[i][1] >= MIN_INTRON]


def read_blocks_junctions(rd):
    """aligned reference blocks [[s0, e]] and junctions [(donor_end, acceptor_start)] from the CIGAR (N >= MIN_INTRON)."""
    blocks, juncs = [], []
    pos = rd.reference_start
    cur_s = pos
    for op, ln in rd.cigartuples:
        if op in (0, 7, 8, 2):          # M, =, X, D consume the reference
            pos += ln
        elif op == 3:                   # N
            if ln >= MIN_INTRON:
                blocks.append([cur_s, pos])
                juncs.append((pos, pos + ln))
                cur_s = pos + ln
            pos += ln
    blocks.append([cur_s, pos])
    return merge([b for b in blocks if b[1] > b[0]]), juncs


def load_loci(path):
    """locus name -> dict(chrom, strand, exons (merged rep exons), juncs set) from a loci GFF3 (gene + exon rows)."""
    loci = {}
    for ln in open(path):
        if ln.startswith("#"):
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] == "gene":
            loci[at["Name"]] = dict(chrom=f[0], strand=f[6], exons=[])
        elif f[2] == "exon":
            loci[at["gene"]]["exons"].append([int(f[3]) - 1, int(f[4])])
    for v in loci.values():
        ex = sorted(v["exons"])
        v["exons"] = merge(ex)
        v["juncs"] = set(introns_of(ex))
    return loci


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--copies", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--family", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--loci", action="append", default=[], help="NAME=loci.gff3 (repeatable)")
    ap.add_argument("--nodes", help="pagedata.json of the read-pool page (its node[arm] per cid)")
    a = ap.parse_args(argv)

    copies = [r for r in csv.DictReader(open(a.copies), delimiter="\t") if r["family"] == a.family]
    genes = {r["isoform_gene"] for r in copies}
    tx = gtf_transcripts(a.truth, genes)
    arms = {}
    for spec in a.loci:
        name, path = spec.split("=", 1)
        arms[name] = load_loci(path)
    nodes = {}
    if a.nodes:
        pd = json.load(open(a.nodes))
        for row in pd.get("rows", pd if isinstance(pd, list) else []):
            nodes[row["cid"]] = row.get("node", {})

    bam = pysam.AlignmentFile(a.bam)
    rows, summary = [], {"copies": len(copies), "spliced_expressed": 0, "arms": {}}
    for arm in arms:
        summary["arms"][arm] = {"old_overlap": 0, "strict_found": 0, "old_overlap_in_npip_nodes": 0, "strict_found_in_npip_nodes": 0}
    for c in copies:
        chrom, strand = c["chrom"], c["strand"]
        texons = tx.get(c["isoform_gene"], {})
        union = merge([e for ex in texons.values() for e in ex]) or [[int(c["terr_lo0"]), int(c["terr_hi"])]]
        ann_introns = {t: set(introns_of(ex)) for t, ex in texons.items()}
        n_intron = max((len(v) for v in ann_introns.values()), default=0)
        k = min(2, n_intron)
        ann_all = set().union(*ann_introns.values()) if ann_introns else set()
        lo, hi = union[0][0], union[-1][1]
        reads = []
        for rd in bam.fetch(chrom, lo, hi):
            if rd.is_unmapped or rd.is_secondary or rd.is_supplementary:
                continue
            if ("-" if rd.is_reverse else "+") != strand:
                continue
            blocks, juncs = read_blocks_junctions(rd)
            if inter(blocks, union) == 0:
                continue
            juncs = [j for j in juncs if j[0] >= lo and j[1] <= hi]
            reads.append((blocks, juncs))
        jcount = collections.Counter(j for _, js in reads for j in js)
        supported = {j for j, n in jcount.items() if n >= MIN_JUNCTION_READS}
        n_unspliced = sum(1 for _, js in reads if not js)
        n_one = sum(1 for _, js in reads if len(js) == 1)
        if k == 0:
            n_support = sum(1 for bl, _ in reads if inter(bl, union) >= COVER * ilen(union))
        else:
            n_support = sum(1 for _, js in reads if sum(1 for j in js if j in supported) >= k)
        n_ann2 = sum(1 for _, js in reads if sum(1 for j in js if j in ann_all) >= min(2, n_intron)) if n_intron else 0
        n_exact = sum(1 for _, js in reads if js and any(set(js) == s for s in ann_introns.values() if s))
        expressed = n_support >= FLOOR
        summary["spliced_expressed"] += expressed
        row = dict(cid=c["cid"], name=c["name"], chrom=chrom, strand=strand, span=f"{lo}-{hi}", n_tx=len(texons), ann_introns=n_intron, k=k,
                   reads=len(reads), unspliced=n_unspliced, one_junction=n_one, supported_junctions=len(supported), support_reads=n_support,
                   ann2_reads=n_ann2, exact_chain_reads=n_exact, spliced_expressed=int(expressed))
        for arm, loci in arms.items():
            same = [L for L in loci.values() if L["chrom"] == chrom and L["strand"] == strand and inter(L["exons"], union) > 0]
            if k == 0:
                strict_loci = [L for L in same if inter(L["exons"], union) >= COVER * ilen(union)]
            else:
                strict_loci = [L for L in same if len(L["juncs"] & supported) >= k]
            old = len(same) > 0
            strict = expressed and len(strict_loci) > 0
            row[f"{arm}_old_overlap_loci"] = len(same)
            row[f"{arm}_strict_loci"] = len(strict_loci)
            row[f"{arm}_strict_found"] = int(strict)
            summary["arms"][arm]["old_overlap"] += old
            summary["arms"][arm]["strict_found"] += strict
            if nodes:
                own = bool(nodes.get(c["cid"], {}).get(arm, False))
                row[f"{arm}_page_own_node"] = int(own)
                summary["arms"][arm]["old_overlap_in_npip_nodes"] += own
                summary["arms"][arm]["strict_found_in_npip_nodes"] += own and strict
        rows.append(row)
    with open(a.out + ".copies.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    json.dump(summary, open(a.out + ".json", "w"), indent=1)
    print(json.dumps(summary))
    hdr = ["name", "reads", "unspliced", "one_junction", "supported_junctions", "support_reads", "ann2_reads", "exact_chain_reads", "spliced_expressed"] + \
          [f"{arm}_{x}" for arm in arms for x in ("strict_found",) + (("page_own_node",) if nodes else ())]
    print("\t".join(hdr))
    for r in rows:
        print("\t".join(str(r.get(h, "")) for h in hdr))


if __name__ == "__main__":
    main()

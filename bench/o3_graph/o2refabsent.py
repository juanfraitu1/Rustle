#!/usr/bin/env python3
"""O3 -> O2 on the real refabsent panel (report-only, docs/REPORT_o2_refabsent_2026-10-24.md): KB3781 fibroblast reads of the 34
families with haplotype-only (B-only) loci, truth from the animal's own haplotypes (refabsent/bonly.tsv, copies_lift.tsv).
  inputs W         W/copies.tsv (refabsent/copies_genes.tsv) and W/net.tsv (refabsent/fam/<family>.fa read names)
  truth W          W/truth.tsv: each read's origin = the lifted copy or B-only locus its best haplotype record overlaps (on the
                   haplotype the family's lift uses; a tie between two loci = unknown)
  score W PREFIX   O3 flags per family vs the expressed B-only loci; O2 on PREFIX.o2.tsv vs truth and vs minimap2 MAPQ >= 1
"""
import collections, csv, os, sys
import pysam

RA = "/mnt/linuxdisk/tmp/rna_allele"
R = f"{RA}/refabsent"


def tsv(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def inputs(W):
    os.makedirs(W, exist_ok=True)
    fams = sorted({r["family"] for r in tsv(f"{R}/bonly.tsv")})
    with open(f"{W}/copies.tsv", "w") as o:
        o.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
        k = collections.Counter()
        for r in tsv(f"{R}/copies_genes.tsv"):
            f = r["gene_id"].split(":")[0]
            if f not in fams:
                continue
            ex = [tuple(map(int, x.split("-"))) for x in r["exons"].split(",")]
            s, e = min(a for a, _ in ex), max(b for _, b in ex)
            o.write(f"{f}\t{k[f]}\t{r['gene_id']}\t{r['chrom']}\t{s}\t{e}\t{len(ex)}\t{r['strand']}\t0\t{r['exons']}\n")
            k[f] += 1
    with open(f"{W}/net.tsv", "w") as o:
        o.write("family\tread\tused\trole\tvia\tclass\tin_reads_tsv\n")
        for f in fams:
            for ln in open(f"{R}/fam/{f}.fa"):
                if ln.startswith(">"):
                    o.write(f"{f}\t{ln[1:].split()[0]}\t1\tS\tA\tNA\t1\n")
    print(f"{len(fams)} families, {sum(k.values())} copies")


def truth(W):
    cm = {}
    for r in tsv(f"{RA}/chrmap.tsv"):
        cm[r["same_name"]] = (r["same_hap"], r["chrom"])
        cm[r["B_name"]] = (r["B_hap"], r["chrom"])
    bams = {h: pysam.AlignmentFile(f"{R}/reads.{h}.bam") for h in ("mat", "pat")}
    name = {h: {} for h in bams}
    for h, b in bams.items():
        for sq in b.references:
            name[h][sq.split("_")[0][3:]] = sq
    loci = collections.defaultdict(list)   # family -> [(label, hap, bam contig, start, end)]
    for r in tsv(f"{R}/copies_lift.tsv"):
        if r["B_chrom"] and r["B_chrom"] in cm:
            h, c = cm[r["B_chrom"]]
            loci[r["gene_id"].split(":")[0]].append((r["gene_id"], h, name[h].get(c), int(r["B_start"]), int(r["B_end"])))
    for r in tsv(f"{R}/bonly.tsv"):
        if r["chrom"] in cm:
            h, c = cm[r["chrom"]]
            loci[r["family"]].append((f"absent:{r['locus']}", h, name[h].get(c), int(r["start"]), int(r["end"])))
    net = collections.defaultdict(set)
    for r in tsv(f"{W}/net.tsv"):
        net[r["family"]].add(r["read"])
    n = collections.Counter()
    with open(f"{W}/truth.tsv", "w") as o:
        o.write("family\tread\ttruth\n")
        for f, L in sorted(loci.items()):
            best = collections.defaultdict(list)    # read -> [(AS, label)]
            for lab, h, ctg, s, e in L:
                if ctg is None:
                    continue
                for rd in bams[h].fetch(ctg, s, e):
                    if rd.query_name in net[f] and not rd.is_unmapped and not rd.is_supplementary and rd.has_tag("AS"):
                        best[rd.query_name].append((rd.get_tag("AS"), lab))
            for rd, xs in best.items():
                xs.sort(reverse=True)
                t = xs[0][1] if len(xs) == 1 or xs[0][0] > xs[1][0] or xs[0][1] == xs[1][1] else "unknown"
                o.write(f"{f}\t{rd}\t{t}\n")
                n["absent" if t.startswith("absent") else t if t == "unknown" else "copy"] += 1
    print(dict(n))


def score(W, prefix):
    tr = {r["read"]: r["truth"] for r in tsv(f"{W}/truth.tsv")}
    cp = {r["tid"]: (r["chrom"], int(r["start"]), int(r["end"])) for r in tsv(f"{W}/copies.tsv")}
    expd = collections.defaultdict(list)
    for r in tsv(f"{R}/bonly_expressed.tsv"):
        if r["expressed"] == "1":
            expd[r["family"]].append(f"{r['locus']}({r['n_reads']})")
    put = collections.Counter(r["family"] for r in tsv(f"{prefix}.putative.tsv"))
    cert = collections.Counter(r["family"] for r in tsv(f"{prefix}.candidates.tsv")) if os.path.exists(f"{prefix}.candidates.tsv") else {}
    print("family\texpressed B-only loci\tcertified\tputative")
    for f in sorted(set(expd) | set(put) | set(cert)):
        print(f"{f}\t{','.join(expd[f]) or '-'}\t{cert.get(f, 0)}\t{put.get(f, 0)}")
    prim, mq = {}, {}
    for rd in pysam.AlignmentFile(f"{R}/R0.bam").fetch(until_eof=True):
        if not (rd.is_secondary or rd.is_supplementary or rd.is_unmapped) and rd.query_name in tr:
            prim[rd.query_name] = next((c for c, (ch, s, e) in cp.items() if ch == rd.reference_name and rd.reference_start < e
                                        and s < rd.reference_end), None)
            mq[rd.query_name] = rd.mapping_quality
    c = collections.Counter()
    for r in tsv(f"{prefix}.o2.tsv"):
        t = tr.get(r["read"])
        if t is None or t == "unknown":
            c["no truth"] += 1
            continue
        k = "absent" if t.startswith("absent") else "copy"
        c[k, "tied"] += 1
        a = r["assigned"]
        if mq.get(r["read"], 0) >= 1:
            c[k, "mq_placed"] += 1
            c[k, "mq_wrong"] += prim.get(r["read"]) != t
        if a in ("ABSTAIN", "HOPELESS"):
            continue
        ok = a.startswith("new:") if k == "absent" else a == t
        c[k, "placed"] += 1
        c[k, "wrong"] += not ok
    for k in ("copy", "absent"):
        print(f"O2 tied reads from {'reference copies' if k == 'copy' else 'B-only loci'}: tied {c[k, 'tied']}; placed "
              f"{c[k, 'placed']} wrong {c[k, 'wrong']}; minimap2 MAPQ>=1 placed {c[k, 'mq_placed']} wrong {c[k, 'mq_wrong']}")
    print(f"tied reads without truth: {c['no truth']}")


if __name__ == "__main__":
    {"inputs": lambda: inputs(sys.argv[2]), "truth": lambda: truth(sys.argv[2]),
     "score": lambda: score(sys.argv[2], sys.argv[3])}[sys.argv[1]]()

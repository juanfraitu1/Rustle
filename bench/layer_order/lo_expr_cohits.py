#!/usr/bin/env python3
"""NPIP/TBC1D3 layer order (integration) — testis read counts for EXPR on the corrected universe.

Rule = heavy/scripts/expr_counts.py as documented in heavy/README.md (EXPR) and re-derived independently by
verify_slim_expr/recount_all.py (0 of 321 genes differ): human_testis.t2t.bam, primary reads only (samtools -F 2308),
read blocks split at N and D (M/=/X consume both), a read counts for gene G if >= 1 block overlaps >= 1 bp of any
exon of G (strand ignored); unique = the read's blocks hit exons of exactly ONE RefSeq gene/pseudogene genome-wide and
that gene is G; exon-less gene/pseudogene records count on a gene-body exon.

This script is the recount_all.py logic with the gene set = heavy/EXPR.counts.tsv (321) + extra ids given on the command
line (universe genes EXPR lacks). It prints any difference against EXPR.counts.tsv (expected: none).

usage: lo_expr_recount.py OUT_TSV EXTRA_ID ...

lo_expr_cohits.py variant (generated from this file): for the genes named on the command line, tabulate which OTHER
genes' exons the same primary reads hit (explains any-vs-unique gaps). usage: lo_expr_cohits.py OUT_TSV GENE ...
"""
import bisect
import csv
import gzip
import re
import subprocess
import sys
from collections import defaultdict

H = "/mnt/linuxdisk/home/juanfraitu/layer_order/npip_tbc1d3/heavy"
BAM = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/human_val/human_testis.t2t.bam"
GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/Reference/chm13v2.0_RefSeq_full.gff.gz"
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def main():
    out_tsv = sys.argv[1]
    extra = set(sys.argv[2:])
    expr = {}
    want = set(extra)
    genes, parent, exraw = {}, {}, []
    with gzip.open(GFF, "rt") as fh:
        for ln in fh:
            if ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            t = f[2]
            if t not in ("gene", "pseudogene", "exon") and "Parent=" not in f[8]:
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            fid, par = a.get("ID"), a.get("Parent")
            if t in ("gene", "pseudogene"):
                s, e = int(f[3]) - 1, int(f[4])
                if fid in genes:
                    c0, s0, e0 = genes[fid]
                    if c0 == f[0]:
                        genes[fid] = (c0, min(s0, s), max(e0, e))
                else:
                    genes[fid] = (f[0], s, e)
                continue
            if t == "exon" and par:
                exraw.append((par.split(",")[0], f[0], int(f[3]) - 1, int(f[4])))
            if fid and par:
                parent.setdefault(fid, par.split(",")[0])

    def up(p):
        k = 0
        while p not in genes:
            p = parent.get(p)
            k += 1
            if p is None or k > 10:
                return None
        return p

    exons = defaultdict(set)
    for p, c, s, e in exraw:
        g = up(p)
        if g:
            exons[g].add((c, s, e))
    for g, (c, s, e) in genes.items():
        if not exons.get(g):
            exons[g] = {(c, s, e)}
    print(f"[gff] gene/pseudogene ids {len(genes)}; distinct exons {sum(len(v) for v in exons.values())}")
    per = defaultdict(list)
    for g, lst in exons.items():
        for c, s, e in lst:
            per[c].append((s, e, g[5:] if g.startswith("gene-") else g))
    idx = {}
    for c, lst in per.items():
        lst.sort()
        idx[c] = ([x[0] for x in lst], lst, max(x[1] - x[0] for x in lst))

    def hits(c, bs, be):
        st, lst, ml = idx[c]
        i, j = bisect.bisect_left(st, bs - ml), bisect.bisect_left(st, be)
        return {lst[k][2] for k in range(i, j) if lst[k][1] > bs and lst[k][0] < be}

    spans = defaultdict(list)
    for g in want:
        c, s, e = genes["gene-" + g]
        spans[c].append((s, e))
    bed = out_tsv + ".windows.bed"
    nwin = 0
    with open(bed, "w") as out:
        for c in sorted(spans):
            iv = sorted(spans[c])
            cs, ce = iv[0]
            for s, e in iv[1:]:
                if s <= ce:
                    ce = max(ce, e)
                else:
                    out.write(f"{c}\t{cs}\t{ce}\n")
                    nwin += 1
                    cs, ce = s, e
            out.write(f"{c}\t{cs}\t{ce}\n")
            nwin += 1
    p = subprocess.Popen(["samtools", "view", "-M", "-F", "2308", "-L", bed, BAM], stdout=subprocess.PIPE, text=True)
    n_any, n_uni = defaultdict(int), defaultdict(int)
    co = defaultdict(lambda: defaultdict(int))
    strand = defaultdict(lambda: defaultdict(int))
    gstrand = {}
    with gzip.open(GFF, "rt") as fh2:
        for ln in fh2:
            f2 = ln.split("\t")
            if len(f2) > 8 and f2[2] in ("gene", "pseudogene"):
                a2 = dict(kv.split("=", 1) for kv in f2[8].strip().split(";") if "=" in kv)
                gstrand[a2.get("ID", "")[5:]] = f2[6]
    seen = set()
    n_rec = n_hit = 0
    for ln in p.stdout:
        f = ln.split("\t", 6)
        key = (f[0], f[1], f[2], f[3])
        if key in seen:
            continue
        seen.add(key)
        n_rec += 1
        c, pos = f[2], int(f[3]) - 1
        hit = set()
        for n, op in CIG.findall(f[5]):
            n = int(n)
            if op in "M=X":
                hit |= hits(c, pos, pos + n)
                pos += n
            elif op in "DN":
                pos += n
        if not hit:
            continue
        n_hit += 1
        for g in hit & want:
            n_any[g] += 1
            for h in hit - {g}:
                co[g][h] += 1
            rs = "-" if int(f[1]) & 16 else "+"
            strand[g]["same" if rs == gstrand.get(g) else "opposite"] += 1
        if len(hit) == 1:
            (g,) = hit
            if g in want:
                n_uni[g] += 1
    assert p.wait() == 0
    print(f"[bam] windows {nwin}; primary records {n_rec}; on >= 1 exon {n_hit}")
    with open(out_tsv, "w") as out:
        out.write("gene\tn_reads_any\tn_reads_unique\tread_strand_same\tread_strand_opposite\tco_hit_genes(n_reads)\n")
        for g in sorted(want):
            cs = ", ".join(f"{h}:{n}" for h, n in sorted(co[g].items(), key=lambda kv: -kv[1]))
            out.write(f"{g}\t{n_any[g]}\t{n_uni[g]}\t{strand[g]['same']}\t{strand[g]['opposite']}\t{cs}\n")
            print(f"{g}: any {n_any[g]} unique {n_uni[g]}; read strand same/opposite {strand[g]['same']}/"
                  f"{strand[g]['opposite']}; co-hit {cs}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""NPIP/TBC1D3 layer order (integration) — testis read counts for EXPR on the corrected universe.

Rule = heavy/scripts/expr_counts.py as documented in heavy/README.md (EXPR) and re-derived independently by
verify_slim_expr/recount_all.py (0 of 321 genes differ): human_testis.t2t.bam, primary reads only (samtools -F 2308),
read blocks split at N and D (M/=/X consume both), a read counts for gene G if >= 1 block overlaps >= 1 bp of any
exon of G (strand ignored); unique = the read's blocks hit exons of exactly ONE RefSeq gene/pseudogene genome-wide and
that gene is G; exon-less gene/pseudogene records count on a gene-body exon.

This script is the recount_all.py logic with the gene set = heavy/EXPR.counts.tsv (321) + extra ids given on the command
line (universe genes EXPR lacks). It prints any difference against EXPR.counts.tsv (expected: none).

usage: lo_expr_recount.py OUT_TSV [--ignore ID,ID,...] EXTRA_ID ...

--ignore (audit 2026-09-16): also count n_reads_unique_mr = reads whose exon hits, after removing the listed records, are
exactly {G} (for a listed record G itself, the other listed records are removed). Used with the RefSeq readthrough
records that overlap a member on the same strand, so 'unique' agrees with the member rule's one-record-per-copy intent.
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
    args = sys.argv[2:]
    ignore = set()
    if args[:1] == ["--ignore"]:
        ignore = set(args[1].split(","))
        args = args[2:]
    extra = set(args)
    expr = {r["gene_id"]: r for r in csv.DictReader(open(f"{H}/EXPR.counts.tsv"), delimiter="\t")}
    want = set(expr) | extra
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
    n_any, n_uni, n_uni_mr = defaultdict(int), defaultdict(int), defaultdict(int)
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
            if ignore and not (hit - ignore - {g}):
                n_uni_mr[g] += 1
        if len(hit) == 1:
            (g,) = hit
            if g in want:
                n_uni[g] += 1
    assert p.wait() == 0
    print(f"[bam] windows {nwin}; primary records {n_rec}; on >= 1 exon {n_hit}")
    diff = 0
    with open(out_tsv, "w") as out:
        out.write("gene_id\tn_reads_any\tn_reads_unique\tin_EXPR_table" + ("\tn_reads_unique_mr" if ignore else "") + "\n")
        for g in sorted(want):
            a, u = n_any[g], n_uni[g]
            if g in expr:
                if (a, u) != (int(expr[g]["n_reads_any"]), int(expr[g]["n_reads_unique"])):
                    diff += 1
                    print("DIFF", g, (a, u), (expr[g]["n_reads_any"], expr[g]["n_reads_unique"]))
            else:
                print(f"EXTRA {g} any {a} unique {u}")
            out.write(f"{g}\t{a}\t{u}\t{'yes' if g in expr else 'no'}" + (f"\t{n_uni_mr[g]}" if ignore else "") + "\n")
    print(f"EXPR.counts.tsv genes recounted {len(expr)}; differing {diff}; extra genes {len(extra - set(expr))}")
    if ignore:
        print(f"unique_mr: ignored records {sorted(ignore)}; genes whose unique count rises: "
              + ", ".join(f"{g} {n_uni[g]}->{n_uni_mr[g]}" for g in sorted(want) if n_uni_mr[g] != n_uni[g]))


if __name__ == "__main__":
    main()

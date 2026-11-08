#!/usr/bin/env python3
"""Register 1583 scorer (docs/PREREG_collapse_copyphasing_2026-11-08.md): copy_phasing hits at r1538's absorbing / null loci.
usage: collapse_cp.py PSV.loci.tsv READ_ORIGIN.tsv NT.bam CP_PREFIX CP.copies.tsv"""
import csv, sys, collections
import pysam
sys.path.insert(0, __import__("os").path.join(__import__("os").path.dirname(__file__), "..", "crossind"))
from copyphasing_lookup import copies, hits, ov

def single_origin(row, origin, bam):
    c = collections.Counter()
    for r in bam.fetch(row["chrom"], int(row["start"]), int(row["end"])):
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            continue
        g = origin.get(r.query_name)
        if g:
            c[g] += 1
    n = sum(c.values())
    if not n or c.most_common(1)[0][0] != row["gene"]:
        return False
    return all(v < 0.10 * n for g, v in c.items() if g != row["gene"])

def main(loci_tsv, origin_tsv, nt, prefix, cp_copies):
    origin = {l.split("\t")[0]: l.split("\t")[1] for l in open(origin_tsv) if not l.startswith("read\t")}
    bam = pysam.AlignmentFile(nt)
    by_key, rows = copies(cp_copies)
    h = hits(prefix, by_key)
    print(f"copies {len(rows)}, certified {sum(k == 'certified' for _, k in h)}, putative {sum(k == 'putative' for _, k in h)}")
    T = collections.defaultdict(lambda: [0, 0])  # (group) -> [testable, hit]
    print("name\tstatus\treason\tkind\tlocus\tstratum\tsingle_origin\ttestable\tcert\tput")
    for r in csv.DictReader(open(loci_tsv), delimiter="\t"):
        if r["kind"] not in ("absorbing", "null"):
            continue
        w = (r["chrom"], int(r["start"]), int(r["end"]))
        xy = "XY" if r["chrom"] in ("chrX", "chrY") else "auto"
        testable = any(ov(c, w) for c in rows)
        cert = sum(k == "certified" and ov(c, w) for c, k in h)
        put = sum(k == "putative" and ov(c, w) for c, k in h)
        so = single_origin(r, origin, bam) if r["kind"] == "null" else None
        groups = [f"{xy} missing", f"{xy} missing {r['reason']}"] if r["kind"] == "absorbing" else \
                 [f"{xy} null all"] + ([f"{xy} null single-origin"] if so else [])
        for g in groups:
            T[g][0] += testable; T[g][1] += testable and (cert + put) > 0
        print(f"{r['name']}\t{r['status']}\t{r['reason']}\t{r['kind']}\t{w[0]}:{w[1]}-{w[2]}\t{xy}\t{so}\t{int(testable)}\t{cert}\t{put}")
    for g in sorted(T):
        t, k = T[g]
        print(f"{g}: hit {k}/{t} testable" + (f" = {k/t:.2f}" if t else ""))
    m, n = T["auto missing"], T["auto null single-origin"]
    if m[0] < 5 or n[0] < 5:
        print("VERDICT UNDERPOWERED"); return
    b1, b2 = m[1] / m[0] >= 0.5, n[1] / n[0] <= 0.25
    print(f"bar1 {'PASS' if b1 else 'FAIL'}  bar2 {'PASS' if b2 else 'FAIL'}  VERDICT {'PASS' if b1 and b2 else 'FAIL'}")

if __name__ == "__main__":
    main(*sys.argv[1:6])

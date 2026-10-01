#!/usr/bin/env python3
"""Truth for docs/PREREG_rna_allele_haplotype_count_2026-10-01.md, step 3: lift each gene's exon union from `_pri` to its B chromosome
(the other haplotype of KB3781) through the primary minimap2 asm5 alignments, and count exonic differences.

Input PAFs: query = `_pri` chromosome in 10 Mb chunks named `<pri>:<offset>`, target = the B chromosome; `--cs` present.
Per gene: lifted = exonic bases covered by aligned columns of primary (tp:A:P) alignments; diffs = mismatching columns + indel events
inside the exons. Classes (lift step only; the "closer B locus elsewhere" refinement of T1 / T? is applied by truth_classes.py):
  lift >= 0.95 and diffs >= 1 -> T2d;  lift >= 0.95 and diffs == 0 -> T2i;  lift < 0.50 -> T1-candidate;  otherwise T?
Genes: a TSV with gene_id, chrom (pri), strand, exon blocks "s0-e,s0-e" (0-based half-open).

    python3 truth_lift.py --chrmap chrmap.tsv --paf-dir out --genes genes.tsv --out lift.tsv
"""
import argparse
import bisect
import collections
import csv
import re

CS = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})")


def merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return out


def walk(rec):
    """Yield (kind, qpos_lo, qpos_hi, t) in global `_pri` coordinates: 'M' aligned run, 'X' mismatch (1 base), 'I'/'D' indel event at a
    query position; t = the target position aligned to the run's first target base. cs follows the target, so for '-' strand
    alignments the query is traversed backwards."""
    q0, q1, strand, cs, t0 = rec
    t = t0
    if strand == "+":
        q = q0
        for op in CS.findall(cs):
            if op[0] == ":":
                n = int(op[1:]); yield ("M", q, q + n, t); q += n; t += n
            elif op[0] == "*":
                yield ("X", q, q + 1, t); q += 1; t += 1
            elif op[0] == "+":          # insertion: bases only in the query
                n = len(op) - 1; yield ("I", q, q + n, t); q += n
            elif op[0] == "-":          # deletion: bases only in the target
                yield ("D", q, q, t); t += len(op) - 1
    else:
        q = q1
        for op in CS.findall(cs):
            if op[0] == ":":
                n = int(op[1:]); yield ("M", q - n, q, t); q -= n; t += n
            elif op[0] == "*":
                yield ("X", q - 1, q, t); q -= 1; t += 1
            elif op[0] == "+":
                n = len(op) - 1; yield ("I", q - n, q, t); q -= n
            elif op[0] == "-":
                yield ("D", q, q, t); t += len(op) - 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for k in ("chrmap", "paf_dir", "genes", "out"):
        ap.add_argument("--" + k.replace("_", "-"), required=True)
    a = ap.parse_args(argv)
    cm = {r["pri"]: r for r in csv.DictReader(open(a.chrmap), delimiter="\t")}
    genes = collections.defaultdict(list)
    for r in csv.DictReader(open(a.genes), delimiter="\t"):
        bl = merge([[int(x) for x in b.split("-")] for b in r["exons"].split(",") if b])
        genes[r["chrom"]].append((bl[0][0], bl[-1][1], r["gene_id"], bl))
    with open(a.out, "w") as out:
        out.write("gene_id\tchrom\texonic_bp\tlifted_bp\tlift_frac\tmismatches\tindels\tclass\tn_alignments\tB_chrom\tB_start\tB_end\n")
        for pri, gl in genes.items():
            if pri not in cm or not cm[pri]["B_name"]:
                for s, e, gid, bl in gl:
                    tot = sum(y - x for x, y in bl)
                    out.write(f"{gid}\t{pri}\t{tot}\t0\t0.000\t0\t0\tT1-candidate\t0\t\t\t\n")   # chrX/chrY in a male: no B
                continue
            chrom = cm[pri]["chrom"]
            runs, mis, ind = [], [], []          # walked ONCE per chromosome: aligned runs, mismatch and indel positions
            tmap = []                            # (q_lo, q_hi, t at the run's first target base, strand) for coordinate lifting
            bname = cm[pri]["B_name"]
            nrec = []
            for ln in open(f"{a.paf_dir}/chr{chrom}.paf"):
                f = ln.rstrip("\n").split("\t")
                if "tp:A:P" not in f[12:]:
                    continue
                off = int(f[0].rsplit(":", 1)[1])
                cs = next(x[5:] for x in f[12:] if x.startswith("cs:Z:"))
                rec = (off + int(f[2]), off + int(f[3]), f[4], cs, int(f[7]))
                nrec.append((rec[0], rec[1]))
                for kind, lo, hi, t in walk(rec):
                    if kind == "M":
                        runs.append([lo, hi]); tmap.append((lo, hi, t, f[4]))
                    elif kind == "X":
                        runs.append([lo, hi]); mis.append(lo)
                    else:
                        ind.append(lo)
            runs = merge(runs)
            rs = [r[0] for r in runs]
            mis.sort(); ind.sort(); nrec.sort(); tmap.sort()
            tq = [r[0] for r in tmap]
            def lift_span(x, y):
                """Target interval covered by the M runs that overlap [x, y) (exact within each run)."""
                lo = hi = None
                i = max(0, bisect.bisect_right(tq, x) - 1)
                while i < len(tmap) and tmap[i][0] < y:
                    q0, q1, t, st = tmap[i]
                    a_, b_ = max(x, q0), min(y, q1)
                    if a_ < b_:
                        if st == "+":
                            u, v = t + (a_ - q0), t + (b_ - q0)
                        else:                      # reverse: the run's first target base pairs with its LAST query base
                            u, v = t + (q1 - b_), t + (q1 - a_)
                        lo = u if lo is None else min(lo, u); hi = v if hi is None else max(hi, v)
                    i += 1
                return lo, hi
            def covered(x, y):
                t, i = 0, max(0, bisect.bisect_right(rs, x) - 1)
                while i < len(runs) and runs[i][0] < y:
                    t += max(0, min(y, runs[i][1]) - max(x, runs[i][0])); i += 1
                return t
            def count(arr, x, y, strict_lo=False):
                return bisect.bisect_left(arr, y) - (bisect.bisect_right(arr, x) if strict_lo else bisect.bisect_left(arr, x))
            for s, e, gid, bl in gl:
                tot = sum(y - x for x, y in bl)
                lb = sum(covered(x, y) for x, y in bl)
                mm = sum(count(mis, x, y) for x, y in bl)
                di = sum(count(ind, x, y) for x, y in bl)
                n = sum(1 for q0, q1 in nrec if q0 < e and s < q1)
                fr = lb / tot if tot else 0.0
                cls = ("T2d" if mm + di else "T2i") if fr >= 0.95 else "T1-candidate" if fr < 0.50 else "T?"
                spans = [lift_span(x, y) for x, y in bl]
                spans = [sp for sp in spans if sp[0] is not None]
                bs = f"{min(u for u, _ in spans)}\t{max(v for _, v in spans)}" if spans else "\t"
                out.write(f"{gid}\t{pri}\t{tot}\t{lb}\t{fr:.3f}\t{mm}\t{di}\t{cls}\t{n}\t{bname if spans else ''}\t{bs}\n")


if __name__ == "__main__":
    main()

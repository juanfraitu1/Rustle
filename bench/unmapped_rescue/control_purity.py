#!/usr/bin/env python3
"""Amendment 36: are the failing control clusters mixtures of reads from different loci? Miniforge python.

    control_purity.py            table over the 12 control runs (needs their classes.a35c.json and the BAMs)"""
import collections
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import discover as D  # noqa: E402
import seeds as SD  # noqa: E402

SPAN = 100_000
FAIL = ("DIVERGED", "UNSUPPORTED", "NOVEL")
GOOD = ("PRESENT", "ALLELE-LIKE")


def is_multi(positions, span=SPAN):
    """positions = [(contig, pos)]: True unless every read lies on one contig within `span` of the others"""
    if len({c for c, _ in positions}) > 1:
        return True
    ps = [p for _, p in positions]
    return bool(ps) and max(ps) - min(ps) > span


def locus_count(positions, gap=SPAN):
    """number of single-linkage groups (a gap over `gap` or another contig starts a new one)"""
    n, prev = 0, None
    for c, p in sorted(positions):
        if prev is None or c != prev[0] or p - prev[1] > gap:
            n += 1
        prev = (c, p)
    return n


def multi_share(rows, classes):
    """share of the reads in clusters of the given classes that lie in multi-locus clusters; None with no such reads"""
    sel = [r for r in rows if r["cls"] in classes]
    tot = sum(r["reads"] for r in sel)
    return sum(r["reads"] for r in sel if r["multi"]) / tot if tot else None


def read_positions(animal, names):
    """name -> (contig, pos) from the primary alignments in the windows `prepare` scans"""
    import pysam
    bam = D.ANIMALS[animal]["bam"]
    h = pysam.AlignmentFile(bam, "rb")
    ctgs = sorted(zip(h.references, h.lengths), key=lambda x: -x[1])[:12]
    regs = [f"{c}:{max(0, L // 2 - 2_500_000) + 1}-{min(L, L // 2 + 2_500_000)}" for c, L in ctgs]
    out = {}
    p = subprocess.Popen(["samtools", "view", "-F", "2308", "-q", "10", bam, *regs], stdout=subprocess.PIPE, text=True)
    for ln in p.stdout:
        f = ln.split("\t", 4)
        if f[0] in names:
            out[f[0]] = (f[2], int(f[3]))
    return out


def run_rows(run, animal, pos):
    d = f"{D.W}/discover_{run}"
    cl = collections.defaultdict(list)
    for r in csv.DictReader(open(f"{d}/clusters.tsv"), delimiter="\t"):
        cl["cl" + r["cluster"]].append(r["read"])
    rows = []
    for r in json.load(open(f"{d}/classes.a35c.json")):
        ps = [pos[n] for n in cl[r["k"]] if n in pos]
        rows.append(dict(run=run, k=r["k"], cls=r["cls"], reads=r["reads"], length=r["length"], R=r["R"], med=r["median_read_divergence"], n_loci=locus_count(ps), multi=is_multi(ps), mapped=len(ps)))
    return rows


def main():
    allrows = []
    for animal in ("ggo_testis", "a119b", "ptr", "ppy"):
        names = set()
        for suf in ("", "_s6", "_s7"):
            fa = f"{D.AD}/{animal}/control{suf}.fa"
            names |= set(SD.read_fa(fa))
        pos = read_positions(animal, names)
        for suf in ("_control", "_control_s6", "_control_s7"):
            allrows += run_rows(animal + suf, animal, pos)
    json.dump(allrows, open(f"{D.W}/control_purity.json", "w"), indent=1)
    for label, classes in (("PRESENT + ALLELE-LIKE", GOOD), ("DIVERGED + UNSUPPORTED + NOVEL", FAIL)):
        sel = [r for r in allrows if r["cls"] in classes]
        sh = multi_share(allrows, classes)
        print(f"{label:32s} clusters {len(sel):4d} reads {sum(r['reads'] for r in sel):5d}  multi-locus clusters {sum(r['multi'] for r in sel):3d}  multi-locus read share {sh:.3f}")
    for c in FAIL + GOOD:
        sel = [r for r in allrows if r["cls"] == c]
        dist = collections.Counter(min(r["n_loci"], 4) for r in sel)
        print(f"  {c:12s} clusters {len(sel):4d} reads {sum(r['reads'] for r in sel):5d} multi share {multi_share(allrows, (c,))}  loci per cluster (4 = 4+): {dict(sorted(dist.items()))}")
    single = [r for r in allrows if r["cls"] in FAIL and not r["multi"]]
    print(f"single-locus failing clusters: {len(single)}; reads {sum(r['reads'] for r in single)}")
    print("  (reads, length, class, R, median read divergence):", [(r["reads"], r["length"], r["cls"][:4], None if r["R"] is None else round(r["R"], 3), None if r["med"] is None else round(r["med"], 4)) for r in sorted(single, key=lambda x: -x["reads"])[:20]])


if __name__ == "__main__":
    main()

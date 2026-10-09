#!/usr/bin/env python3
"""Amendment 41: maternal-assembly truth for the real-read recall test (docs/PREREG_unmapped_rescue_2026-10-08.md). Miniforge python.

Pure helpers (tested in test_mat_truth.py): cs_divergence, mask, assembly_version, group_loci, verdict.
Steps (resumable, each one call): see main()."""
import re

CS = re.compile(r"(:\d+|\*[a-z][a-z]|[+-][a-z]+|~[a-z]{2}\d+[a-z]{2})")
CIG = re.compile(r"(\d+)([MIDNSHP=X])")
DELTA = 0.00958


def cs_divergence(tstart, cs):
    """-> ([(target position, divergent bases)], target end): substitutions (1 base), deletions from the target (their length, at their start), insertions into the
    query (their length, at the target position where they occur); the cs walks the target forward from tstart"""
    pos, ev = tstart, []
    for tok in CS.findall(cs):
        c = tok[0]
        if c == ":":
            pos += int(tok[1:])
        elif c == "*":
            ev.append((pos, 1))
            pos += 1
        elif c == "-":
            ev.append((pos, len(tok) - 1))
            pos += len(tok) - 1
        elif c == "+":
            ev.append((pos, len(tok) - 1))
        else:  # ~ intron (not produced by asm5)
            pos += int(re.findall(r"\d+", tok)[0])
    return ev, pos


def _merge(iv):
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def mask(tlen, covered, events, window=500, max_div=DELTA / 2):
    """target intervals with no counterpart within max_div: not covered by any record, or in a window whose divergent bases exceed max_div x window"""
    cov = _merge([(max(0, a), min(tlen, b)) for a, b in covered if a < b])
    out, prev = [], 0
    for a, b in cov:
        if a > prev:
            out.append((prev, a))
        prev = max(prev, b)
    if prev < tlen:
        out.append((prev, tlen))
    per = {}
    for p, n in events:
        per[p // window] = per.get(p // window, 0) + n
    for w, n in per.items():
        if n > max_div * window:
            out.append((w * window, min(tlen, (w + 1) * window)))
    return _merge(out)


def assembly_version(fetch, pos0, cigar):
    """the reference bases under an alignment's blocks (M = X D), introns (N) skipped, read insertions and clips ignored; fetch(a, b) -> reference[a:b]"""
    out, pos = [], pos0
    for n, op in CIG.findall(cigar):
        n = int(n)
        if op in "M=XD":
            out.append(fetch(pos, pos + n))
            pos += n
        elif op == "N":
            pos += n
    return "".join(out)


def group_loci(reads, gap=1000):
    """reads = [(contig, start, end, name)] -> [dict(contig, start, end, names)], joined when a read starts within `gap` of the current locus end"""
    loci = []
    for c, s, e, n in sorted(reads):
        if loci and loci[-1]["contig"] == c and s <= loci[-1]["end"] + gap:
            loci[-1]["end"] = max(loci[-1]["end"], e)
            loci[-1]["names"].append(n)
        else:
            loci.append(dict(contig=c, start=s, end=e, names=[n]))
    return loci


def verdict(statuses):
    """TRUE when most placed reads (mat / pat / shared) are mat- or pat-specific; WRONG otherwise (including no placed read)"""
    placed = [s for s in statuses if s in ("mat", "pat", "shared")]
    spec = sum(s in ("mat", "pat") for s in placed)
    return "TRUE" if placed and spec > len(placed) / 2 else "WRONG"


# ---------------------------------------------------------------- steps
import json  # noqa: E402
import os  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402

RA = "/mnt/linuxdisk/tmp/rna_allele"
HAPS = "/mnt/linuxdisk/home/juanfraitu/gorilla_haps"
PRI = "/mnt/linuxdisk/home/juanfraitu/_from_wsl/winloci_scratch/GGO.fasta"
FIB = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
OUT = "/mnt/linuxdisk/tmp/o3_rescue/mattruth"


def chrmap():
    rows = [ln.rstrip("\n").split("\t") for ln in open(f"{RA}/chrmap.tsv")][1:]
    return [dict(pri=r[0], chrom=r[1], same=r[2], bhap=r[5] if len(r) > 5 else "", bname=r[6] if len(r) > 6 else "") for r in rows]


def step_mask(hap="mat"):
    """mask.<hap>.bed: <hap> bases with no primary counterpart within delta/2 (chromosomes whose primary is the other haplotype) + every unplaced <hap> contig"""
    os.makedirs(OUT, exist_ok=True)
    lens = {ln.split("\t")[0]: int(ln.split("\t")[1]) for ln in open(f"{HAPS}/{hap}.fa.fai")}
    placed = {r["bname"] for r in chrmap() if r["bhap"] == hap} | {r_ for r in chrmap() for r_ in [r.get("same_name")] if r_}
    rows, tot = [], 0
    for r in chrmap():
        if r["bhap"] != hap:
            continue
        covered, events = [], []
        for ln in open(f"{RA}/out/chr{r['chrom']}.paf"):
            f = ln.rstrip("\n").split("\t")
            if "tp:A:P" not in f[12:]:
                continue
            ts, te = int(f[7]), int(f[8])
            covered.append((ts, te))
            cs = next(t[5:] for t in f[12:] if t.startswith("cs:Z:"))
            ev, _ = cs_divergence(ts, cs)
            events += ev
        for a, b in mask(lens[r["bname"]], covered, events):
            rows.append((r["bname"], a, b, f"chr{r['chrom']}"))
            tot += b - a
    cm = {ln.rstrip("\n").split("\t")[6] for ln in open(f"{RA}/chrmap.tsv") if ln.count("\t") >= 6} | {ln.split("\t")[3] for ln in open(f"{RA}/chrmap.tsv")}
    un = [c for c in lens if c not in cm]
    for c in un:
        rows.append((c, 0, lens[c], "unplaced"))
    with open(f"{OUT}/mask.{hap}.bed", "w") as o:
        for c, a, b, lab in rows:
            o.write(f"{c}\t{a}\t{b}\t{lab}\n")
    print(f"mask.{hap}: {len(rows) - len(un)} intervals, {tot:,} bp on the placed chromosomes; {len(un)} unplaced contigs, {sum(lens[c] for c in un):,} bp")


if __name__ == "__main__":
    {"mask": lambda: step_mask(sys.argv[2] if len(sys.argv) > 2 else "mat")}[sys.argv[1]]()

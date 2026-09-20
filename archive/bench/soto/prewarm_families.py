#!/usr/bin/env python3
"""Precompute and cache everything the mainstay families need, so later analysis is instant.

WHY. Every session re-derives the same objects for the same handful of families: member coordinates,
discovered loci, locus FASTA, the sensitive all-vs-all, E_r, read-supported exons, sec_frac. The BAM
calls dominate — a single sec_frac is ~640 ms cold and ~0.1 ms cached (9,000x); read_exons is ~160 ms
cold. A family with 30 loci is ~25 s of samtools every time somebody asks a question about it.

WHAT IT CACHES. Under $RUSTLE_BENCH_CACHE (default ~/winloci_scratch/benchcache), keyed on each input
file's mtime+size, so a re-aligned BAM or a new reference invalidates itself. Nothing here is a
"result" — it is all recomputable; deleting the cache costs time, never correctness.

MAINSTAY FAMILIES. NPIP and TBC1D3 are the two the thesis narrative rests on; the rest are the curated
Soto families large enough to be scorable. Pass --families to override.

⚠ ONE HEAVY JOB AT A TIME (WSL2). This is serial by construction. Do not run it beside a cargo build or
  another whole-genome pass.
⚠ Uses rustlib's canonical primitives ONLY. The whole point is that nothing here re-implements a rule.

    python3 bench/soto/prewarm_families.py                 # mainstays
    python3 bench/soto/prewarm_families.py --all           # every curated family with >= 3 members
    python3 bench/soto/prewarm_families.py --report        # what is cached, no work
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rustlib as R  # noqa: E402

TRUTH = os.environ.get(
    "RUSTLE_TRUTH", "/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/80_fams.gene_preferred.bed")
BAM = os.environ.get("RUSTLE_BAM", "/home/juanfra/winloci_scratch/chr16_sub.bam")
FULL_BAM = os.environ.get("RUSTLE_FULL_BAM",
                          "/mnt/linuxdisk/home/juanfraitu/winloci_data/A119b.t2t.bam")
WORK = os.environ.get("RUSTLE_PREWARM", "/home/juanfra/winloci_scratch/prewarm")

# NPIP and TBC1D3 carry the narrative; the others are the scorable curated families.
MAINSTAY = ["ID_154", "ID_468", "ID_207", "ID_113", "ID_8", "ID_116", "ID_14", "ID_400", "ID_431"]


def load_families():
    fam = defaultdict(list)
    for line in open(TRUTH):
        f = line.rstrip("\n").split("\t")
        if len(f) >= 4 and "|" in f[3]:
            g, fid = f[3].split("|", 1)
            fam[fid].append((f[0], int(f[1]), int(f[2]), g))
    return fam


def pick_bam(chrom):
    """chr16 has a local whole-chromosome slice; everything else needs the full BAM."""
    if chrom == "chr16" and os.path.exists(BAM):
        return BAM
    return FULL_BAM if os.path.exists(FULL_BAM) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--families", default=",".join(MAINSTAY))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()

    fam = load_families()
    fids = ([f for f, m in fam.items() if len(m) >= 3] if a.all
            else [f for f in a.families.split(",") if f in fam])
    os.makedirs(WORK, exist_ok=True)

    if a.report:
        n = sum(1 for _ in os.scandir(R.CACHE)) if os.path.isdir(R.CACHE) else 0
        sz = sum(e.stat().st_size for e in os.scandir(R.CACHE)) if os.path.isdir(R.CACHE) else 0
        print(f"cache dir {R.CACHE}\n  {n} entries, {sz/1e6:.1f} MB")
        for fid in fids:
            print(f"  {fid:<9} {len(fam[fid]):>3} members  "
                  f"{sorted({c for c, _, _, _ in fam[fid]})}")
        return 0

    print(f"prewarming {len(fids)} families -> {R.CACHE}\n")
    t0 = time.time()
    for fid in fids:
        mem = fam[fid]
        chroms = sorted({c for c, _, _, _ in mem})
        bam = pick_bam(chroms[0])
        if bam is None:
            print(f"  {fid:<9} SKIPPED — no BAM for {chroms} "
                  f"(mount /mnt/linuxdisk for the full A119b)")
            continue
        t = time.time()
        sf, ex, miss = [], 0, 0
        for c, s, e, _g in mem:
            if pick_bam(c) != bam:
                miss += 1
                continue
            sf.append(R.sec_frac(bam, c, s, e)["sec_frac"])
            ex += len(R.read_exons(bam, c, s, e))
        # locus FASTA + the sensitive all-vs-all, both reusable by any later question
        reg = f"{WORK}/{fid}.regions"
        with open(reg, "w") as fh:
            for c, s, e, _g in mem:
                fh.write(f"{c}:{s+1}-{e}\n")
        fa, paf = f"{WORK}/{fid}.fa", f"{WORK}/{fid}_sens.paf"
        if not os.path.exists(fa) or os.path.getsize(fa) == 0:
            with open(fa, "w") as fh:
                subprocess.run(["samtools", "faidx", R.FASTA, "-r", reg],
                               stdout=fh, stderr=subprocess.DEVNULL)
        if not os.path.exists(paf) or os.path.getsize(paf) == 0:
            with open(paf, "w") as fh:
                subprocess.run(["minimap2", "-x", "asm20", "-k11", "-w5", "-c", "--eqx",
                                "-N", "200", "-p", "0.02", "-t", "2", fa, fa],
                               stdout=fh, stderr=subprocess.DEVNULL)
        E = R.er_edges([(paf, R.SHIPPED_TIERS[0][1])])
        names = [f"{c}:{s+1}-{e}" for c, s, e, _g in mem]
        C = R.components(names, E)
        med = sorted(sf)[len(sf) // 2] if sf else float("nan")
        print(f"  {fid:<9} {len(mem):>3} members  {len(E):>5} edges  {len(C):>2} comp  "
              f"largest {len(C[0]):>3}  density {R.density(C[0], E):.3f}  "
              f"med sec_frac {med:.3f}  exon blocks {ex:>4}"
              + (f"  ⚠{miss} members on another chromosome, skipped" if miss else "")
              + f"   [{time.time()-t:.1f}s]")
    print(f"\ndone in {time.time()-t0:.1f}s — rerun is near-instant")
    n = sum(1 for _ in os.scandir(R.CACHE)) if os.path.isdir(R.CACHE) else 0
    sz = sum(e.stat().st_size for e in os.scandir(R.CACHE)) if os.path.isdir(R.CACHE) else 0
    print(f"cache: {n} entries, {sz/1e6:.1f} MB at {R.CACHE}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

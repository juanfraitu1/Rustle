#!/usr/bin/env python3
"""Run a graph-O3 world (prereg docs/PREREG_o3_graph_sim_2026-10-15.md): map its reads, then the current method.

    run.py map WORLD [--threads 4] [--index MMI]   WORLD/E0.bam: reads.E0.fq, the baseline's exact minimap2 command,
                                                   sorted + indexed
    run.py cur WORLD BINARY [--threads 4]    WORLD/cur/A.*: copy_assign --no-as-tied-only --psv-clusters on the catalog
"""
import argparse
import os
import subprocess
import sys

MM2 = "minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes"


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def cmd_map(a):
    w = a.world
    target = a.index or f"{w}/genome.ref.fa"           # a prebuilt .mmi for a whole genome
    sh(f"{MM2} -t {a.threads} {target} {w}/reads.E0.fq 2> {w}/E0.mm2.log | samtools sort -@2 -o {w}/E0.bam - "
       f"&& samtools index {w}/E0.bam")
    print(f"[o3_graph] {w}/E0.bam", flush=True)


def cmd_cur(a):
    w = a.world
    os.makedirs(f"{w}/cur", exist_ok=True)
    if not os.path.exists(f"{w}/regions.txt"):            # a real-genome world writes its families' regions
        with open(f"{w}/regions.txt", "w") as fh:
            for ln in open(f"{w}/genome.ref.fa.fai"):
                c, n = ln.split("\t")[:2]
                fh.write(f"{c}:1-{n}\n")
    sh(f"{a.binary} --bam {w}/E0.bam --fasta {w}/genome.ref.fa --regions {w}/regions.txt --families {w}/P.fam.copies.tsv "
       f"--copies-fa {w}/P.fam.copies.fa --no-as-tied-only --psv-clusters --region-threads 1 --skip-poa-diagnostic "
       f"--threads {a.threads} --out {w}/cur/A > {w}/cur/A.log 2>&1")
    print(f"[o3_graph] {w}/cur/A.psv_clusters.tsv", flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("map")
    p.add_argument("world")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--index", help="a prebuilt minimap2 index of the world's genome")
    p.set_defaults(func=cmd_map)
    p = sub.add_parser("cur")
    p.add_argument("world")
    p.add_argument("binary")
    p.add_argument("--threads", type=int, default=4)
    p.set_defaults(func=cmd_cur)
    a = ap.parse_args(argv)
    a.func(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())

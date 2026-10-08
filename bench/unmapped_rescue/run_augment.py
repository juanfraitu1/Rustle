#!/usr/bin/env python3
"""Augmented-reference realignment on bed H (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 7).

    run_augment.py prep               # W/bedHhalf/{pool.fa,labels.tsv}: building half of the unmapped D reads (seed 1) + the 959 background reads
    (then:  run_bed.py bedHhalf ; consensus.py W/bedHhalf/pool.fa W/bedHhalf/registered/clusters.tsv W/bedHhalf/registered/cons.fa)
    run_augment.py map half|full      # align reads to the consensus sequences alone, combine with the stored genome primaries, write W/bedHhalf/augment_<kind>.json
Heavy step (minimap2 of about 60k reads against a few hundred kb) goes under tools/rlock.sh heavy."""
import csv
import json
import os
import random
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import augment as G  # noqa: E402
import seeds as SD  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
LT = "/mnt/linuxdisk/tmp/rna_allele/linktest"
MM2 = "-ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes"
CIG = re.compile(r"(\d+)([MIDNSHP=X])")


def qcov(cigar, rlen):
    return sum(int(n) for n, op in CIG.findall(cigar) if op in "M=XI") / rlen if rlen else 0.0


def tag(rec, name, default=None):
    for f in rec[11:]:
        if f.startswith(name + ":"):
            return f.split(":", 2)[2]
    return default


def primaries(sam_lines):
    """{read: record dict} of the primary alignment of every mapped read"""
    out = {}
    for ln in sam_lines:
        if ln[0] == "@":
            continue
        f = ln.rstrip("\n").split("\t")
        flag = int(f[1])
        if flag & 2308 or f[2] == "*":
            continue
        rlen = sum(int(n) for n, op in CIG.findall(f[5]) if op in "M=XIS")
        out[f[0]] = dict(score=int(tag(f, "AS", 0)), de=float(tag(f, "de", 1.0)), qcov=qcov(f[5], rlen), ref=f[2], mapq=int(f[4]))
    return out


def prep():
    lab = {r["read"]: r for r in csv.DictReader(open(f"{W}/bedH/labels.tsv"), delimiter="\t")}
    seqs = SD.read_fa(f"{W}/bedH/pool.fa")
    dunm = sorted(n for n, r in lab.items() if r["role"] == "D")
    random.Random(1).shuffle(dunm)
    build = set(dunm[:len(dunm) // 2])
    keep = [n for n in sorted(seqs) if lab[n]["role"] == "bg" or n in build]
    os.makedirs(f"{W}/bedHhalf", exist_ok=True)
    SD.write_fa(f"{W}/bedHhalf/pool.fa", seqs, keep)
    with open(f"{W}/bedHhalf/labels.tsv", "w") as o:
        o.write("read\tfamily\tcopy\trole\n")
        for n in keep:
            r = lab[n]
            o.write(f"{n}\t{r['family']}\t{r['copy']}\t{r['role']}\n")
    json.dump(sorted(build), open(f"{W}/bedHhalf/build_half.json", "w"))
    print(f"building half: {len(build)} D + {len(keep) - len(build)} background; held out {len(dunm) - len(build)} unmapped D")


def cluster_family(clusters_tsv, cons_fa):
    """{consensus name: majority family of its cluster}; consensus headers are '<cluster id>|n=<size>'"""
    maj = {}
    for r in csv.DictReader(open(clusters_tsv), delimiter="\t"):
        maj["cl" + r["cluster"]] = r["majority"] or None
    names = [ln[1:].strip().split()[0] for ln in open(cons_fa) if ln[0] == ">"]
    return {n: maj.get(n.split("|")[0]) for n in names}


def map_reads(kind):
    half = kind == "half"
    cdir = f"{W}/bedHhalf/registered" if half else f"{W}/bedH/registered"
    out = f"{W}/bedHhalf/augment_{kind}"
    os.makedirs(out, exist_ok=True)
    cons_fa = f"{cdir}/cons.fa"
    labs = {r["read"]: r for r in csv.DictReader(open(f"{LT}/labels.tsv"), delimiter="\t")}
    unm = {r["read"]: r for r in csv.DictReader(open(f"{W}/bedH/labels.tsv"), delimiter="\t")}  # unmapped D + bg
    build = set(json.load(open(f"{W}/bedHhalf/build_half.json"))) if half else set()
    genome = primaries(subprocess.run(f"samtools view -F 2308 {LT}/R.bam", shell=True, capture_output=True, text=True, check=True).stdout.splitlines())
    cls, fam = {}, {}
    for n, r in labs.items():
        if r["role"] == "S":
            cls[n] = "S"
        elif n in unm:
            cls[n] = "D_unm"
        else:
            cls[n] = "D_abs"
        fam[n] = r["family"]
    for n, r in unm.items():
        if r["role"] == "bg":
            cls[n], fam[n] = "bg", "bg"
    for n in build:
        cls.pop(n, None)  # reads that built the consensus are not tested in the non-circular run
    reads = SD.read_fa(f"{W}/bedH/pool.fa")
    fq = f"{out}/reads.fa"
    if not os.path.exists(fq):
        subprocess.run(f"samtools fasta -F 2308 {LT}/R.bam > {out}/mapped.fa", shell=True, check=True)
        allseq = SD.read_fa(f"{out}/mapped.fa")
        allseq.update(reads)
        SD.write_fa(fq, allseq, sorted(n for n in cls if n in allseq))
        os.remove(f"{out}/mapped.fa")
    sam = f"{out}/cons.sam"
    if not os.path.exists(sam + ".done"):
        subprocess.run(f"minimap2 {MM2} -t 4 {cons_fa} {fq} > {sam}", shell=True, check=True)
        open(sam + ".done", "w").write("ok")
    cons = primaries(open(sam))
    cf = cluster_family(f"{cdir}/clusters.tsv", cons_fa)
    rows = [(cls[n], fam[n], genome.get(n), cons.get(n), (cons[n]["ref"] if n in cons else None)) for n in sorted(cls)]
    res = dict(kind=kind, cons_sequences=len(cf), cons_with_family=sum(v is not None for v in cf.values()), reads_tested=len(rows),
               all_families=G.move_metrics(rows, cf))
    withcons = {v for v in cf.values() if v}
    res["families_with_consensus"] = G.move_metrics([r for r in rows if r[0] in ("D_unm", "D_abs") and r[1] in withcons], cf)
    json.dump(res, open(f"{W}/bedHhalf/augment_{kind}.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items()}, indent=1))


if __name__ == "__main__":
    prep() if sys.argv[1] == "prep" else map_reads(sys.argv[2])

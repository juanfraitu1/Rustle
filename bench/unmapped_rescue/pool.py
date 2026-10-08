#!/usr/bin/env python3
"""Pools of unmapped reads with truth labels (docs/PREREG_unmapped_rescue_2026-10-08.md section 2).

    pool.py bedA    # W/bedA/{pool.fa,labels.tsv}: unmapped records of o3_excise/masked.bam + the 959 background reads
"""
import csv
import json
import os
import subprocess
import sys

W = "/mnt/linuxdisk/tmp/o3_rescue"
EXC = "/mnt/linuxdisk/home/juanfraitu/o3_excise"
PANEL = "/home/juanfra/winloci_scratch/o3_excise/panel.json"
BG = "/mnt/linuxdisk/tmp/o3_mat/reads/R_unm.fa"


def label_reads(rows, panel):
    """rows: {read: (chrom, start)} the baseline primary placement; panel: [dict(fam, mask_gene, mask=[chrom, s, e], keep_gene, keep=[chrom, s, e])].
    -> {read: (family, copy, 'D'|'S'|'other')}: D = starts inside an erased interval, S = inside a surviving one."""
    out = {}
    for n, (c, s) in rows.items():
        hit = (None, None, "other")
        for p in panel:
            if p["mask"][0] == c and p["mask"][1] <= s < p["mask"][2]:
                hit = (p["fam"], p["mask_gene"], "D")
                break
            if p["keep"][0] == c and p["keep"][1] <= s < p["keep"][2]:
                hit = (p["fam"], p["keep_gene"], "S")
                break
        out[n] = hit
    return out


def bedA():
    os.makedirs(f"{W}/bedA", exist_ok=True)
    panel = json.load(open(PANEL))
    base = {}
    for ln in open(f"{EXC}/baseline.tsv"):
        f = ln.rstrip("\n").split("\t")
        base[f[0]] = (f[1], int(f[2]))
    lab = label_reads(base, panel)
    subprocess.run(f"samtools view -b -f 4 {EXC}/masked.bam | samtools fasta - > {W}/bedA/unmapped.fa", shell=True, check=True)
    n = {"D": 0, "S": 0, "other": 0, "bg": 0}
    with open(f"{W}/bedA/pool.fa", "w") as o, open(f"{W}/bedA/labels.tsv", "w") as t:
        t.write("read\tfamily\tcopy\trole\n")
        for src, bg in ((f"{W}/bedA/unmapped.fa", False), (BG, True)):
            name = None
            for ln in open(src):
                if ln[0] == ">":
                    name = ln[1:].strip().split()[0]
                    fam, cp, role = ("bg", "bg", "bg") if bg else lab.get(name, (None, None, "other"))
                    t.write(f"{name}\t{fam or ''}\t{cp or ''}\t{role}\n")
                    n[role] += 1
                o.write(ln if ln[0] == ">" else ln)
    print("pool:", n, "families with D reads:", len({l[0] for l in lab.values() if l[2] == 'D'}))


COPIES = "/mnt/linuxdisk/tmp/rna_allele/refabsent/copies.fa"


def targets(bed, erased):
    """W/<bed>/targets.fa: the 915 catalog copies minus the erased ones (named FAM:idx)"""
    w, n = False, 0
    with open(f"{W}/{bed}/targets.fa", "w") as o:
        for ln in open(COPIES):
            if ln[0] == ">":
                w = ln[1:].strip().split()[0] not in erased
                n += w
            if w:
                o.write(ln)
    return n


def bedH():
    """A13's 53 held-out multi-copy families (rna_allele/linktest): unmapped records of the masked alignment + the 959 background reads"""
    L = "/mnt/linuxdisk/tmp/rna_allele/linktest"
    os.makedirs(f"{W}/bedH", exist_ok=True)
    lab = {r["read"]: (r["family"], r["copy"], r["role"]) for r in csv.DictReader(open(f"{L}/labels.tsv"), delimiter="\t")}
    subprocess.run(f"samtools view -b -f 4 {L}/R.bam | samtools fasta - > {W}/bedH/unmapped.fa", shell=True, check=True)
    n = {"D": 0, "S": 0, "other": 0, "bg": 0}
    with open(f"{W}/bedH/pool.fa", "w") as o, open(f"{W}/bedH/labels.tsv", "w") as t:
        t.write("read\tfamily\tcopy\trole\n")
        for src, bg in ((f"{W}/bedH/unmapped.fa", False), (BG, True)):
            for ln in open(src):
                if ln[0] == ">":
                    name = ln[1:].strip().split()[0]
                    fam, cp, role = ("bg", "bg", "bg") if bg else lab.get(name, (None, None, "other"))
                    t.write(f"{name}\t{fam or ''}\t{cp or ''}\t{role}\n")
                    n[role] += 1
                o.write(ln)
    erased = {p["mask"][3] for p in json.load(open(f"{L}/panel.json"))}
    print("pool:", n, "targets:", targets("bedH", erased), "erased:", len(erased))


if __name__ == "__main__":
    {"bedA": bedA, "bedH": bedH}[sys.argv[1]]()

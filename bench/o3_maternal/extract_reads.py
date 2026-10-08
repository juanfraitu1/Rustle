#!/usr/bin/env python3
"""R_LRP and R_unm of docs/PREREG_o3_maternal_reference_2026-10-08.md section 4.

    extract_reads.py lrp    # reads with a primary/secondary record on one of the 11 LRPAP1 loci -> W/reads/R_LRP.{fa,names.tsv}
    extract_reads.py unm    # unmapped primaries of the fibroblast BAM -> W/reads/R_unm.fa
"""
import collections
import csv
import os
import random
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

BAM = "/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam"
LRP = "/mnt/linuxdisk/tmp/lrpap1"
R34 = f"{C.TRUTH}/refabsent/labels.tsv"
CAP = 2000


def loci():
    """the 11 LRPAP1 loci (8 full-length + 3 fragments): (cid, name, chrom, lo0, hi), pri coordinates"""
    out = []
    for f in ("lrpap1.copies.tsv", "partial.copies.tsv"):
        for r in csv.DictReader(open(f"{LRP}/{f}"), delimiter="\t"):
            out.append((r["cid"], r["name"], r["chrom"], int(r["terr_lo0"]), int(r["terr_hi"])))
    return out


def net_names(bam_path, loci_, cap=CAP, seed=1):
    """{read: [cid, ...]}: reads with a primary or secondary record on a locus, at most `cap` per locus (seeded shuffle of the sorted names)"""
    rng = random.Random(seed)
    out = collections.defaultdict(list)
    with pysam.AlignmentFile(bam_path) as bam:
        for cid, _name, chrom, lo, hi in loci_:
            names = sorted({rd.query_name for rd in bam.fetch(chrom, lo, hi) if not (rd.is_unmapped or rd.is_supplementary)})
            rng.shuffle(names)
            for n in names[:cap]:
                out[n].append(cid)
    return dict(out)


def merge_labels(r34, lrp):
    """r34: {read: family} of the 34-family set; lrp: {read: [cid]}. -> (rows, overlap): rows = [(read, 'LRPAP1', 'cid,cid')] for reads
    NOT in the 34-family set; overlap = reads in both (they keep their 34-family label and are still LRPAP1-locus reads via names.tsv)"""
    rows, overlap = [], []
    for n, cids in sorted(lrp.items()):
        if n in r34:
            overlap.append(n)
        else:
            rows.append((n, "LRPAP1", ",".join(sorted(cids))))
    return rows, overlap


def sequences(bam_path, names, out_fa):
    """primary-record sequences in the original read orientation (samtools fasta restores the strand); ONE scan of the BAM"""
    nf = out_fa + ".names"
    open(nf, "w").write("\n".join(sorted(names)) + "\n")
    subprocess.run(f"samtools view -b -F 2308 -N {nf} -@ 4 {bam_path} | samtools fasta -@ 2 - > {out_fa}", shell=True, check=True)


def lrp():
    os.makedirs(f"{C.W}/reads", exist_ok=True)
    r34 = {r["read"]: r["family"] for r in csv.DictReader(open(R34), delimiter="\t")}
    names = net_names(BAM, loci())
    rows, overlap = merge_labels(r34, names)
    with open(f"{C.W}/reads/R_LRP.names.tsv", "w") as o:
        o.write("read\tcids\tin_r34\n")
        for n, cids in sorted(names.items()):
            o.write(f"{n}\t{','.join(sorted(cids))}\t{int(n in r34)}\n")
    sequences(BAM, [r[0] for r in rows], f"{C.W}/reads/R_LRP.fa")
    got = sum(1 for ln in open(f"{C.W}/reads/R_LRP.fa") if ln[0] == ">")
    print(f"LRPAP1 net reads {len(names)}; new (not in the 34-family set) {len(rows)}; sequences written {got}; in both {len(overlap)}")


def unm():
    os.makedirs(f"{C.W}/reads", exist_ok=True)
    out = f"{C.W}/reads/R_unm.fa"
    subprocess.run(f"samtools view -b -f 4 {BAM} '*' | samtools fasta - > {out}", shell=True, check=True)
    print("unmapped primaries", sum(1 for ln in open(out) if ln[0] == ">"))


if __name__ == "__main__":
    {"lrp": lrp, "unm": unm}[sys.argv[1]]()

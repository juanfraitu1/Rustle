#!/usr/bin/env python3
"""Gorilla STON1: the mother's chromosome 12 carries it, the father's (and the reference NC_073236.2, which is the father's) does
not -- the data for the gorilla tab of the NPIP Read Pools page (user 2026-10-20: "can we make another tab in this with those
data from gorilla?").

    ston1_fig.py OUT_DIR MAT.fa PAT.fa HSA_TX.fa READS_MAT.bam READS_PAT.bam   ->  OUT_DIR/ston1.json

Windows: mother CM054594.2:95,700,000-96,300,000, father CM054570.2:90,900,000-91,500,000 (as on 2026-10-14). Human mRNAs
(STON1 NM_006873, GTF2A1L NM_006872, LHCGR NM_000233) spliced-aligned to each window place the genes; mother vs father
alignment (minimap2 -x asm5 -c) gives the ribbons; exact 32-mers give the dot plot; the 113 STON1 reads' alignments to each
parent's whole genome (splice:hq) give where they go.
"""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("npip_cnv_fig_s1", os.path.join(HERE, "npip_cnv_fig.py"))
F = importlib.util.module_from_spec(spec)
spec.loader.exec_module(F)

MAT = ("CM054594.2", 95700000, 96300000)
PAT = ("CM054570.2", 90900000, 91500000)


def genes_on(window_seq, tx_fa):
    """Each human mRNA spliced-aligned to the window (minimap2 -x splice -uf, best hit) -> {name: [(s, e), ...] exons}."""
    with tempfile.TemporaryDirectory() as d:
        open(f"{d}/w.fa", "w").write(f">w\n{window_seq}\n")
        sam = subprocess.run(["minimap2", "-a", "-x", "splice", "-uf", "--secondary=no", f"{d}/w.fa", tx_fa],
                             capture_output=True, text=True, check=True).stdout
    out = {}
    for ln in sam.splitlines():
        if ln.startswith("@"):
            continue
        x = ln.split("\t")
        flag = int(x[1])
        if flag & 4 or flag & 2304:
            continue
        name = x[0].split("_")[0]
        pos = int(x[3]) - 1
        import re
        ex, cur = [], pos
        start = cur
        for n, op in re.findall(r"(\d+)([MIDNSHP=X])", x[5]):
            n = int(n)
            if op in "M=XD":
                cur += n
            elif op == "N":
                ex.append((start, cur))
                cur += n
                start = cur
        ex.append((start, cur))
        aligned = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", x[5]) if op in "M=XI")
        qlen = sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", x[5]) if op in "M=XIS")
        out[name] = dict(exons=ex, strand="-" if flag & 16 else "+", qcov=round(aligned / max(1, qlen), 3))
    return out


def main(out_dir, mat_fa, pat_fa, tx_fa, rm_bam, rp_bam):
    os.makedirs(out_dir, exist_ok=True)
    ms = F.fetch(mat_fa, *MAT)
    ps = F.fetch(pat_fa, *PAT)
    rows = F.paf(ps, ms, extra=("-x", "asm5"))                 # mother (query) onto father (target)
    rib = [dict(q0=r["qs"], q1=r["qe"], t0=r["ts"], t1=r["te"], strand=r["strand"], mapq=r["mapq"], primary=r["tp"] == "P",
                identity=round(r["nm"] / r["alen"], 4)) for r in rows if r["qe"] - r["qs"] >= 2000]
    reads = {}
    for side, bam_p, (ch, lo, hi) in (("mother", rm_bam, MAT), ("father", rp_bam, PAT)):
        placed, unmapped, inwin = 0, 0, []
        for rec in pysam.AlignmentFile(bam_p, check_sq=False).fetch(until_eof=True):
            if rec.is_secondary or rec.is_supplementary:
                continue
            if rec.is_unmapped:
                unmapped += 1
                continue
            placed += 1
            if rec.reference_name == ch and lo <= rec.reference_start < hi:
                de = rec.get_tag("de") if rec.has_tag("de") else None
                inwin.append((rec.reference_start - lo, rec.reference_end - lo, de))
        reads[side] = dict(placed=placed, unmapped=unmapped, in_window=sorted(inwin))
    data = dict(mother=dict(contig=MAT[0], start=MAT[1], end=MAT[2], genes=genes_on(ms, tx_fa)),
                father=dict(contig=PAT[0], start=PAT[1], end=PAT[2], genes=genes_on(ps, tx_fa)),
                ribbons=rib, dots=F.dots(ms, ps), reads=reads)
    json.dump(data, open(f"{out_dir}/ston1.json", "w"))
    print("[ston1] genes mother:", {k: (v["exons"][0][0], v["exons"][-1][1], v["qcov"]) for k, v in data["mother"]["genes"].items()})
    print("[ston1] genes father:", {k: (v["exons"][0][0], v["exons"][-1][1], v["qcov"]) for k, v in data["father"]["genes"].items()})
    print("[ston1] ribbons:", [(r["q0"], r["q1"], r["t0"], r["t1"], r["strand"], r["identity"]) for r in rib])
    print("[ston1] reads:", {k: (v["placed"], v["unmapped"], len(v["in_window"])) for k, v in reads.items()})
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:7]))

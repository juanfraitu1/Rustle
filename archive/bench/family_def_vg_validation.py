#!/usr/bin/env python3
"""
family_def_vg_validation.py — VG-style structural validation of read-conflict family
edges, with ISOFORM AGGREGATION done right.

Premise (user): a COPY is a gene; its reads are ISOFORMS. To avoid alt-splicing
masquerading as extra copies, each copy = the EXON-UNION of all its isoform reads
("full gene minus introns"; retained introns appear where reads span them). The VG
backbone is then the shared exonic sequence ACROSS copies (PSVs = inter-copy
variation); alt-splicing is intra-copy (bubbles), not extra copies.

Test: for candidate read-conflict edges (real families + bridges), build each copy's
all-isoform exon-union sequence from the reads, align copies, and measure the backbone:
  - reciprocal backbone coverage (parallel full-length paths = real; short bubble = bridge)
  - identity, and isoform count per copy (the pollution check: many isoforms -> ONE union)
  - compare to the longest-single-isoform cDNA homology (recovers what that truth missed?)

Run: /home/juanfra/miniforge3/bin/python bench/family_def_vg_validation.py
"""
import collections
import os
import subprocess
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

BAM = "/home/juanfra/winloci_scratch/GGO.bam"
GENOME = "/home/juanfra/winloci_scratch/GGO.fasta"
GENES_BED = "/home/juanfra/winloci_scratch/unmapped_poc/genes.bed"
WORK = "/home/juanfra/winloci_scratch/vg_val_seqs.fa"
PAD = 20000          # capture reads extending past the (possibly fragmentary) annotation
MAX_READS = 5000

EDGES = [
    ("RABL2A", "RABL2B", "real_DNA_paralog"),
    ("TBC1D1", "LOC134756953", "real_subbar"),
    ("RABL2B", "LOC134756389", "real_subbar"),
    ("AK6", "LOC115934278", "real_denovo_copy"),
    ("CCDC196", "LOC129526440", "real_denovo_copy"),
    ("OCLN", "SEPTIN7", "bridge_no_homology"),
    ("BCAS4", "CCDC30", "bridge_no_homology"),
    ("GSTM2", "LOC115933235", "bridge_colocated_no_hom"),
]


def load_coords():
    d = {}
    with open(GENES_BED) as f:
        for line in f:
            c, s, e, name = line.rstrip("\n").split("\t")
            if name not in d:
                d[name] = (c, int(s), int(e))
    return d


def merge(iv):
    iv.sort()
    out = []
    for s, e in iv:
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def copy_model(bam, genome, chrom, start, end):
    """exon-union (merged aligned blocks) over reads overlapping [start,end] padded;
    returns (exon_union_intervals, n_reads, n_distinct_isoforms, aggregated_sequence)."""
    blocks = []
    sigs = set()
    n = 0
    s0, e0 = max(0, start - PAD), end + PAD
    for r in bam.fetch(chrom, s0, e0):
        if r.is_unmapped or r.is_supplementary:
            continue
        # keep reads whose footprint overlaps the annotated locus
        if r.reference_end is None or r.reference_end < start or r.reference_start > end:
            continue
        b = r.get_blocks()
        if not b:
            continue
        blocks.extend(b)
        sigs.add(tuple((s, e - s) for s, e in b))   # isoform signature
        n += 1
        if n >= MAX_READS:
            break
    if not blocks:
        return [], 0, 0, ""
    union = merge([list(x) for x in blocks])
    seq = "".join(genome.fetch(chrom, s, e).upper() for s, e in union)
    return union, n, len(sigs), seq


def main():
    coords = load_coords()
    bam = pysam.AlignmentFile(BAM, "rb")
    genome = pysam.FastaFile(GENOME)

    genes = sorted({g for e in EDGES for g in e[:2]})
    models = {}
    print("[build] all-isoform exon-union copy models from reads ...", flush=True)
    with open(WORK, "w") as out:
        for g in genes:
            c, s, e = coords[g]
            union, nr, niso, seq = copy_model(bam, genome, c, s, e)
            ulen = sum(b - a for a, b in union)
            models[g] = dict(reads=nr, isoforms=niso, exon_union_bp=ulen,
                             annot_bp=e - s, seqlen=len(seq))
            print(f"  {g:16} reads={nr:5d} isoforms={niso:4d} exon-union={ulen:6d}bp "
                  f"(annot span {e-s:7d}bp)", flush=True)
            if seq:
                out.write(f">{g}\n{seq}\n")
    bam.close(); genome.close()

    # all-vs-all of the aggregated copy sequences
    print("\n[align] minimap2 all-vs-all of the all-isoform copy models ...", flush=True)
    paf = subprocess.run(["minimap2", "-cx", "asm20", "-N", "20", "-p", "0.05",
                          WORK, WORK], capture_output=True, text=True).stdout
    # best per unordered pair: identity + coverage in each direction
    H = {}
    for line in paf.splitlines():
        f = line.split("\t")
        if len(f) < 12:
            continue
        q, ql, qs, qe, t, m, bl = f[0], int(f[1]), int(f[2]), int(f[3]), f[5], int(f[9]), int(f[10])
        if q == t:
            continue
        ident = m / bl if bl else 0
        cov = (qe - qs) / ql if ql else 0
        a, b = (q, t) if q < t else (t, q)
        r = H.setdefault((a, b), {"id": 0.0, "cov_a": 0.0, "cov_b": 0.0})
        r["id"] = max(r["id"], ident)
        if q == a:
            r["cov_a"] = max(r["cov_a"], cov)
        else:
            r["cov_b"] = max(r["cov_b"], cov)

    print("\n" + "=" * 92)
    print("VG-VALIDATION — all-isoform exon-union backbone per candidate edge")
    print("=" * 92)
    print(f"  {'edge':30} {'class':24} {'id':>5} {'cov_a':>6} {'cov_b':>6} {'min':>5} {'max':>5}  verdict")
    for a, b, cls in EDGES:
        key = (a, b) if a < b else (b, a)
        r = H.get(key, {"id": 0.0, "cov_a": 0.0, "cov_b": 0.0})
        mincov = min(r["cov_a"], r["cov_b"]); maxcov = max(r["cov_a"], r["cov_b"])
        # true bridge: no backbone in EITHER direction. clean family: reciprocal backbone.
        # asymmetric: smaller copy covered but larger copy inflated (coarse annotation vertex).
        if r["id"] < 0.80 or maxcov < 0.20:
            verdict = "NO backbone (BRIDGE)"
        elif mincov >= 0.50:
            verdict = "PARALLEL PATHS (family)"
        else:
            verdict = "ASYMMETRIC (real copy, coarse vertex)"
        print(f"  {a+'~'+b:30} {cls:24} {r['id']:5.2f} {r['cov_a']:6.2f} {r['cov_b']:6.2f} "
              f"{mincov:5.2f} {maxcov:5.2f}  {verdict}")

    print("\n--- isoform-pollution check: copies have many isoforms but ONE exon-union backbone ---")
    for g in genes:
        m = models[g]
        ratio = m["exon_union_bp"] / m["annot_bp"] if m["annot_bp"] else 0
        flag = "  <-- read-union >> annotation (recovered full copy)" if ratio > 2 and m["annot_bp"] < 5000 else ""
        print(f"  {g:16} {m['isoforms']:4d} distinct isoforms -> 1 union ({m['exon_union_bp']}bp){flag}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Prereg Addendum AD-1: is the nucleotide unreachability of old retrogenes synonymous-codon saturation?

For each parent-copy pair: longest-CDS proteins aligned with MAFFT, codon alignment by back-translation, identity at
codon positions 1/2/3, codon difference classes (0/1/2/3 positions) and the synonymous share of multi-difference codons,
dN/dS by PAML yn00, and protein-level reach of the parent protein onto the copy's region with miniprot (coverage and
intron loss).

usage: codon_divergence.py --gff full.gff.gz --genome genome.fa --outdir DIR --mafft mafft --yn00 yn00 --miniprot miniprot
"""
import argparse
import collections
import gzip
import os
import re
import subprocess

import pysam

UNREACHED = [("GK", "GK2"), ("CETN2", "CETN1"), ("NAP1L1", "NAP1L2"), ("NAP1L1", "NAP1L3"), ("CSTF2", "CSTF2T"),
             ("UBL4A", "UBL4B"), ("FAM50A", "FAM50B")]
REACHED = [("MKRN1", "MKRN3"), ("PDHA1", "PDHA2"), ("UTP14A", "UTP14C"), ("PABPC1", "PABPC3"), ("POU5F1", "POU5F1B"),
           ("PGK1", "PGK2"), ("GLUD1", "GLUD2"), ("TAF1", "TAF1L"), ("RPL10", "RPL10L")]
BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k] for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}
COMP = str.maketrans("ACGTN", "TGCAN")


def load(gff, names):
    genes, tx_gene, cds = {}, {}, collections.defaultdict(list)
    with gzip.open(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            if f[2] == "gene" and a.get("Name") in names and a.get("gene_biotype") == "protein_coding":
                genes.setdefault(a["Name"], {"id": a["ID"], "chrom": f[0], "start0": int(f[3]) - 1, "end": int(f[4]),
                                             "strand": f[6]})
            elif f[2] == "mRNA" and "Parent" in a:
                tx_gene[a["ID"]] = a["Parent"]
            elif f[2] == "CDS" and "Parent" in a:
                cds[a["Parent"]].append((int(f[3]) - 1, int(f[4]), int(f[7]) if f[7] != "." else 0))
    by_gene = collections.defaultdict(list)
    for tx, segs in cds.items():
        by_gene[tx_gene.get(tx)].append(segs)
    for g in genes.values():
        txs = by_gene.get(g["id"], [])
        g["cds"] = sorted(max(txs, key=lambda s: sum(e - b for b, e, _ in s)), key=lambda x: x[0]) if txs else []
    return genes


def cds_seq(genome, g):
    segs = g["cds"] if g["strand"] == "+" else g["cds"][::-1]
    parts = [genome.fetch(g["chrom"], s, e).upper() for s, e, _ in g["cds"]]
    seq = "".join(parts)
    if g["strand"] == "-":
        seq = seq.translate(COMP)[::-1]
    phase = segs[0][2]
    seq = seq[phase:]
    seq = seq[:len(seq) - len(seq) % 3]
    junctions, pos = [], -phase
    for s, e, _ in segs[:-1]:
        pos += e - s
        junctions.append(pos / 3)
    prot = "".join(CODE.get(seq[i:i + 3], "X") for i in range(0, len(seq), 3))
    if prot.endswith("*"):
        prot, seq = prot[:-1], seq[:-3]
    return seq, prot, junctions


def mafft_pair(mafft, outdir, tag, p1, p2):
    fa = f"{outdir}/{tag}.prot.fa"
    open(fa, "w").write(f">a\n{p1}\n>b\n{p2}\n")
    out = subprocess.run([mafft, "--auto", "--quiet", fa], capture_output=True, text=True, check=True).stdout
    seqs, name = {}, None
    for line in out.splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return "".join(seqs["a"]), "".join(seqs["b"])


def codon_aln(aln, cds):
    out, k = [], 0
    for c in aln:
        if c == "-":
            out.append("---")
        else:
            out.append(cds[3 * k:3 * k + 3])
            k += 1
    return out


def yn00(yn, outdir, tag, c1, c2):
    keep = [(x, y) for x, y in zip(c1, c2) if "-" not in x + y and len(x) == 3 and len(y) == 3
            and CODE.get(x, "*") != "*" and CODE.get(y, "*") != "*"]
    s1, s2 = "".join(x for x, _ in keep), "".join(y for _, y in keep)
    d = f"{outdir}/{tag}_yn00"
    os.makedirs(d, exist_ok=True)
    open(f"{d}/seq.phy", "w").write(f"  2  {len(s1)}\na          {s1}\nb          {s2}\n")
    open(f"{d}/yn00.ctl", "w").write("seqfile = seq.phy\noutfile = yn.out\nverbose = 0\nicode = 0\nweighting = 0\n"
                                     "commonf3x4 = 0\n")
    subprocess.run([yn, "yn00.ctl"], cwd=d, capture_output=True, text=True)
    txt = open(f"{d}/yn.out").read() if os.path.exists(f"{d}/yn.out") else ""
    m = re.search(r"\(B\) Yang & Nielsen \(2000\).*?\n\s*2\s+1\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+"
                  r"([-\d.]+)\s+\+-\s+[-\d.]+\s+([-\d.]+)", txt, re.S)
    if not m:
        return float("nan"), float("nan")
    return float(m.group(6)), float(m.group(7))


def miniprot_reach(mp, genome, outdir, tag, prot, junctions, g):
    L = g["end"] - g["start0"]
    s, e = max(0, g["start0"] - L), g["end"] + L
    rfa, pfa = f"{outdir}/{tag}.region.fa", f"{outdir}/{tag}.parent.faa"
    open(rfa, "w").write(f">r\n{genome.fetch(g['chrom'], s, e).upper()}\n")
    open(pfa, "w").write(f">p\n{prot}\n")
    out = subprocess.run([mp, "-t", "2", rfa, pfa], capture_output=True, text=True).stdout
    best = None
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 12 or f[0] != "p":
            continue
        AS = next((int(x[5:]) for x in f[12:] if x.startswith("AS:i:")), 0)
        cg = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
        if best is None or AS > best[0]:
            best = (AS, int(f[1]), int(f[2]), int(f[3]), cg, int(f[9]), int(f[10]))
    if best is None:
        return 0.0, None, float("nan")
    AS, qlen, qs, qe, cg, nm, bl = best
    introns = sum(1 for _, op in re.findall(r"(\d+)([MIDFGNUV])", cg) if op in "NUV")
    covers_junction = any(qs < j < qe for j in junctions)
    loss = introns == 0 and covers_junction
    return (qe - qs) / qlen, loss, nm / bl if bl else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--mafft", required=True)
    ap.add_argument("--yn00", required=True)
    ap.add_argument("--miniprot", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    genome = pysam.FastaFile(a.genome)
    names = {x for p in UNREACHED + REACHED for x in p}
    genes = load(a.gff, names)
    missing = sorted(names - set(genes))
    if missing:
        print("genes not found as protein_coding:", missing)
    rows = []
    for group, pairs in (("UNREACHED", UNREACHED), ("REACHED", REACHED)):
        for pa, co in pairs:
            if pa not in genes or co not in genes or not genes[pa]["cds"] or not genes[co]["cds"]:
                print(f"{group} {pa}-{co}: missing CDS")
                continue
            tag = f"{pa}_{co}"
            c1, p1, j1 = cds_seq(genome, genes[pa])
            c2, p2, _ = cds_seq(genome, genes[co])
            a1, a2 = mafft_pair(a.mafft, a.outdir, tag, p1, p2)
            k1, k2 = codon_aln(a1, c1), codon_aln(a2, c2)
            both = [(x, y, u, v) for x, y, u, v in zip(a1, a2, k1, k2) if x != "-" and y != "-"]
            n = len(both)
            prot_id = sum(1 for x, y, _, _ in both if x == y) / n
            pos = [sum(1 for _, _, u, v in both if u[i] == v[i]) / n for i in range(3)]
            nt_id = sum(pos) / 3
            diff = collections.Counter(sum(u[i] != v[i] for i in range(3)) for _, _, u, v in both)
            multi = [(u, v) for _, _, u, v in both if sum(u[i] != v[i] for i in range(3)) >= 2]
            syn_multi = sum(1 for u, v in multi if CODE.get(u) == CODE.get(v) and CODE.get(u) not in (None, "*"))
            dn, ds = yn00(a.yn00, a.outdir, tag, k1, k2)
            cov, loss, mp_id = miniprot_reach(a.miniprot, genome, a.outdir, tag, p1, j1, genes[co])
            rows.append({"group": group, "pair": f"{pa}-{co}", "len_aa": (len(p1), len(p2)), "aligned_codons": n,
                         "internal_stops": (p1.count("*"), p2.count("*")), "prot_id": prot_id, "p1": pos[0], "p2": pos[1],
                         "p3": pos[2], "nt_id": nt_id, "d0": diff[0] / n, "d1": diff[1] / n, "d2": diff[2] / n, "d3": diff[3] / n,
                         "multi_syn": syn_multi / len(multi) if multi else float("nan"), "dN": dn, "dS": ds,
                         "mp_cov": cov, "mp_intron_loss": loss, "mp_id": mp_id})
    hdr = ("group pair aa(parent,copy) codons prot_id p1 p2 p3 nt_id d0 d1 d2 d3 syn_share_of_>=2diff dN dS "
           "miniprot_cov miniprot_intron_loss")
    print(hdr.replace(" ", "\t"))
    for r in rows:
        print("\t".join([r["group"], r["pair"], f"{r['len_aa'][0]},{r['len_aa'][1]}", str(r["aligned_codons"])]
                        + [f"{r[k]:.3f}" for k in ("prot_id", "p1", "p2", "p3", "nt_id", "d0", "d1", "d2", "d3", "multi_syn", "dN", "dS", "mp_cov")]
                        + [str(r["mp_intron_loss"])]))
    U = [r for r in rows if r["group"] == "UNREACHED"]
    sat = sum(1 for r in U if r["p3"] < min(r["p1"], r["p2"]) and r["dS"] >= 1.0 and r["prot_id"] >= 0.70)
    rescue = sum(1 for r in U if r["mp_cov"] >= 0.50 and r["mp_intron_loss"])
    print(f"\nSYNONYMOUS SATURATION: {sat}/{len(U)} unreached pairs (p3 < p1,p2; dS >= 1.0; protein identity >= 0.70) -> "
          f"{'SUPPORTED' if sat >= 5 else 'NOT SUPPORTED'}")
    print(f"PROTEIN LEVEL RESCUE: {rescue}/{len(U)} unreached pairs reached by miniprot (>= 50% of parent protein) with "
          f"intron loss -> {'SUPPORTED' if rescue >= 6 else 'NOT SUPPORTED'}")
    allm = [r["multi_syn"] for r in rows if r["multi_syn"] == r["multi_syn"]]
    print(f"synonymous share of codons differing at >= 2 positions: median {sorted(allm)[len(allm) // 2]:.3f} over {len(allm)} pairs")


if __name__ == "__main__":
    main()

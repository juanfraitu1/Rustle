#!/usr/bin/env python3
"""Prototype: the DNA (segdup) edge oracle via GENOMIC-SPAN projection -- a light, targeted stand-in for
running SEDEF/BISER on the whole assembly (too heavy for this box). Contrast with --dna-family-fallback,
which projects the EXPRESSED transcript: that misses paralogs whose transcript diverged (LRRC37BP1/ANKRD36B).

For each DNA-only floor member we project its GENOMIC span (introns INCLUDED, from soto_members.fa) onto the
genome at segdup identity (asm20, id>=0.90) and count the distinct homologous loci = its DNA-family copies.
This is exactly the homology SEDEF reports (genomic, expression-independent), so it should recover the members
the transcript projection could not, and -- crucially -- surface NOVEL segdup copies not in Soto's annotation
(the extra info a real DNA oracle buys). Copy NUMBER only; per-read resolution still needs parCN.

Light on this box: ONE batched minimap2 index load. Run in background.
Run: /home/juanfra/miniforge3/bin/python bench/soto/soto_dna_oracle_prototype.py
"""
import csv
import os
import subprocess
from collections import defaultdict

D = "/mnt/linuxdisk/home/juanfraitu/winloci_data"
MEMFA = f"{D}/soto_members.fa"
FASTA = f"{D}/chm13v2.0.fa"
BED = "bench/soto/80_fams.chr.bed"
CEIL = "bench/soto/soto_rna_homology_ceiling.tsv"
MM2 = "/home/juanfra/miniforge3/bin/minimap2"
S = "/home/juanfra/winloci_scratch"
ID_FLOOR = 0.90


def load_fa(p):
    seqs = {}; h = None; b = []
    for line in open(p):
        if line.startswith(">"):
            if h:
                seqs[h] = "".join(b)
            h = line[1:].strip(); b = []
        else:
            b.append(line.strip())
    if h:
        seqs[h] = "".join(b)
    return seqs


def merge(iv):
    iv = sorted(iv); out = []
    for s, e in iv:
        if out and s <= out[-1][1] + 1000:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def main():
    seqs = load_fa(MEMFA)
    dna_only = [(r["gene"], r["chrom"], int(r["start"]), int(r["end"]))
                for r in csv.DictReader(open(CEIL), delimiter="\t") if "DNA-family only" in r["class"]]

    qfa = f"{S}/dnaonly_genomic.fa"; n_ex = 0
    with open(qfa, "w") as f:
        for g, c, s, e in dna_only:
            hdr = f"{c}:{s+1}-{e}"
            if hdr in seqs:
                f.write(f">{g}|{c}:{s}-{e}\n{seqs[hdr]}\n"); n_ex += 1
    print(f"extracted {n_ex}/{len(dna_only)} DNA-only genomic spans -> projecting (asm20, id>={ID_FLOOR})...")

    paf = f"{S}/dnaonly_genomic.paf"
    subprocess.run(f"{MM2} -cx asm20 -N50 -p0.1 -t6 {FASTA} {qfa} > {paf} 2>/dev/null", shell=True, check=True)

    mem = defaultdict(list)
    for line in open(BED):
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4 and p[0].startswith("chr"):
            mem[p[0]].append((int(p[1]), int(p[2])))
    def annotated(c, s, e):
        return any(not (a > e or b < s) for a, b in mem.get(c, ()))

    self_span = {g: (c, s, e) for g, c, s, e in dna_only}
    raw = defaultdict(list)
    for line in open(paf):
        f = line.split("\t")
        if len(f) < 12:
            continue
        g = f[0].split("|")[0]
        tgt, ts, te, matches, blk = f[5], int(f[7]), int(f[8]), int(f[9]), int(f[10])
        de = next((float(x[5:]) for x in f[12:] if x.startswith("de:f:")), None)
        idp = (1 - de) if de is not None else matches / blk
        if idp >= ID_FLOOR:
            raw[g].append((tgt, ts, te, idp))

    print("\n=== DNA edge oracle (genomic-span projection) vs transcript projection ===")
    print(f"{'member':13s} {'segdup_copies':>13} {'famCN':>6} {'novel(unannot)':>14}  copies")
    recovered = novel_total = 0
    rows = []
    for g, c, s, e in dna_only:
        hits = raw.get(g, [])
        # exclude the member's own locus; merge the rest into distinct copies
        others = [(t, a, b) for (t, a, b, i) in hits if not (t == c and not (a > e or b < s))]
        bychr = defaultdict(list)
        for t, a, b in others:
            bychr[t].append((a, b))
        copies = [(t, a, b) for t, ivs in bychr.items() for a, b in merge(ivs)]
        nnov = sum(1 for t, a, b in copies if not annotated(t, a, b))
        famcn = len(copies) + 1
        if copies:
            recovered += 1
        novel_total += nnov
        tag = ",".join(f"{t}:{a}{'*' if not annotated(t,a,b) else ''}" for t, a, b in copies[:4])
        rows.append([g, c, s, e, len(copies), famcn, nnov, tag])
        print(f"  {g:11s} {len(copies):>13} {famcn:>6} {nnov:>14}  {tag[:52]}")

    n = len(dna_only)
    print(f"\n  DNA-only members recovered as DNA-families (genomic): {recovered}/{n}  (transcript projection: 10/21)")
    print(f"  NOVEL segdup copies found (not in Soto annotation): {novel_total}  (* marked above -- the extra info a DNA oracle buys)")
    with open("bench/soto/soto_dna_oracle_prototype.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene", "chrom", "start", "end", "segdup_copies", "famCN", "n_novel_unannotated", "copies"])
        w.writerows(rows)
    print("  wrote bench/soto/soto_dna_oracle_prototype.tsv")


if __name__ == "__main__":
    main()

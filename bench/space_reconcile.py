#!/usr/bin/env python3
"""Prereg Addendum AN-2 / AN-3: how protein-, RNA- and DNA-space families relate on one annotation.

usage: space_reconcile.py --nodes REFSEQ_NODES --pfam PROTEIN_FAMILIES --pedges PROTEIN_EDGES --dna DNA_E1_PREFIX
         --paf GENE_SPAN_PAF --expr EXPR_TSV --genome FA --contigs c1,c2 --out DIR [--miniprot BIN]
AN-2: over coding genes of the protein proteome, co-member pairs of the protein families (P) vs the E1 DNA families (D) and
  their RNA restriction (expressed genes, u >= 3). Discordant pairs are described by the best direct protein edge
  (amino-acid identity, coverage) and the best gene-span alignment record (nucleotide identity, aligned bp).
AN-3: pseudogene loci (biotype contains "pseudogene") reached by miniprot from any family protein with an alignment
  covering >= 0.30 of the protein; attachment = the best-scoring protein's family; agreement with the pseudogene's DNA family.
"""
import argparse
import collections
import csv
import itertools
import os
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import protein_families as pf  # noqa: E402


def key1(chrom, start0, end):
    return f"{chrom}:{start0 + 1}-{end}"


def main():
    ap = argparse.ArgumentParser()
    for k in ("--nodes", "--pfam", "--pedges", "--dna", "--paf", "--expr", "--genome", "--contigs", "--out"):
        ap.add_argument(k, required=True)
    ap.add_argument("--miniprot", default="/home/juanfra/miniforge3/envs/prot/bin/miniprot")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    contigs = set(a.contigs.split(","))
    nodes = {r["idx"]: r for r in csv.DictReader(open(a.nodes), delimiter="\t") if r["chrom"] in contigs}
    names = {r["idx"]: (r["name"], r["biotype"]) for r in csv.DictReader(open(a.nodes + ".names.tsv"), delimiter="\t")}
    key_of = {k: key1(r["chrom"], int(r["start"]), int(r["end"])) for k, r in nodes.items()}
    idx_of = {v: k for k, v in key_of.items()}
    P = {r["idx"]: r["family_id"] for r in csv.DictReader(open(a.pfam), delimiter="\t")}
    proteome = {l[1:].strip() for l in open(a.pfam.replace(".families.tsv", ".proteins.faa")) if l.startswith(">")}
    # DNA families: clusters.tsv members (representatives) + loci.tsv folds
    D = {}
    for r in csv.DictReader(open(a.dna + ".clusters.tsv"), delimiter="\t"):
        D[f"{r['chrom']}:{r['start']}-{r['end']}"] = r["cluster_id"]
    if os.path.exists(a.dna + ".loci.tsv"):
        for r in csv.DictReader(open(a.dna + ".loci.tsv"), delimiter="\t"):
            if r["representative"] in D and r["annotation"] not in D:
                D[r["annotation"]] = D[r["representative"]]
    Dg = {k: D[key_of[k]] for k in nodes if key_of[k] in D}
    expr = {}
    for r in csv.DictReader(open(a.expr), delimiter="\t"):
        expr[key1(r["chrom"], int(r["start"]), int(r["end"]))] = int(r["u"])
    expressed = {k for k in nodes if expr.get(key_of[k], 0) >= 3}
    # direct protein edges and best gene-span PAF record per gene pair
    pedge = {}
    for r in csv.DictReader(open(a.pedges), delimiter="\t"):
        pedge[(r["u"], r["v"])] = (float(r["identity"]), float(r["coverage"]))
    coding = sorted((k for k in proteome if k in nodes), key=int)
    cset = set(coding)
    nt = {}
    for line in open(a.paf):
        f = line.split("\t", 12)
        u, v = idx_of.get(f[0]), idx_of.get(f[5])
        if u in cset and v in cset and u != v:
            k = (min(u, v, key=int), max(u, v, key=int))
            ident, bl = int(f[9]) / int(f[10]), int(f[10])
            if bl >= 300 and ident * bl > nt.get(k, (0, 0))[0] * nt.get(k, (0, 0))[1]:
                nt[k] = (ident, bl)

    def pairs(labels, universe):
        g = collections.defaultdict(list)
        for k in universe:
            if k in labels:
                g[labels[k]].append(k)
        return {(min(x, y, key=int), max(x, y, key=int)) for ms in g.values() for x, y in itertools.combinations(ms, 2)}
    print("=== AN-2: protein (P) vs DNA E1 (D) families on coding genes")
    for level, universe in (("DNA", cset), ("RNA (expressed u>=3)", cset & expressed)):
        pp, dp = pairs(P, universe), pairs(Dg, universe)
        inP = len({k for k in universe if k in P})
        inD = len({k for k in universe if k in Dg})
        both = pp & dp
        print(f"[{level}] coding genes {len(universe)}; in P families {inP}, in D families {inD}; "
              f"P pairs {len(pp)}, D pairs {len(dp)}, shared {len(both)}; "
              f"D pairs inside P {len(both) / max(1, len(dp)):.3f}; P pairs inside D {len(both) / max(1, len(pp)):.3f}")
        if level != "DNA":
            continue

        def describe(ps, tag):
            c = collections.Counter()
            ex = collections.defaultdict(list)
            for p in ps:
                pe = pedge.get(p)
                n = nt.get(p)
                pa = "no direct protein edge" if pe is None else ("aa<0.50" if pe[0] < 0.5 else ("aa 0.50-0.70" if pe[0] < 0.7 else "aa>=0.70"))
                na = "no nt record>=300bp id>=0" if n is None else ("nt<0.70" if n[0] < 0.7 else "nt>=0.70")
                c[(pa, na)] += 1
                if len(ex[(pa, na)]) < 4:
                    ex[(pa, na)].append((names[p[0]][0], names[p[1]][0], pe and round(pe[0], 2), n and round(n[0], 2)))
            print(f"  {tag}: {len(ps)} pairs")
            for k, v in c.most_common():
                print(f"    {k[0]:22s} {k[1]:26s} {v:6d}  e.g. {ex[k][:3]}")
        describe(pp - dp, "P-only (protein family, not DNA family)")
        describe(dp - pp, "D-only (DNA family, not protein family)")
        describe(both, "shared")
    # AN-3: pseudogenes reached from family proteins
    print("=== AN-3: pseudogene loci in protein space")
    genome = pysam.FastaFile(a.genome)
    pseudo = [k for k in nodes if "pseudogene" in names[k][1]]
    tfa, pfa, gff = f"{a.out}/pseudo.fa", a.pfam.replace(".families.tsv", ".proteins.faa"), f"{a.out}/pseudo_miniprot.gff"
    with open(tfa, "w") as fh:
        for k in pseudo:
            r = nodes[k]
            fh.write(f">{k}\n{genome.fetch(r['chrom'], int(r['start']), int(r['end'])).upper()}\n")
    fam_prot = [l for l in open(pfa)]
    fam_faa = f"{a.out}/family_proteins.faa"
    with open(fam_faa, "w") as fh:
        keep = False
        for line in fam_prot:
            if line.startswith(">"):
                keep = line[1:].strip() in P
            if keep:
                fh.write(line)
    if not os.path.exists(gff):
        with open(gff + ".tmp", "w") as fh:
            subprocess.run([a.miniprot, "--gff", "-t", "4", "-N", "1000", "--outn=1000", "--outs=0", "-p", "0", tfa, fam_faa],
                           stdout=fh, stderr=subprocess.DEVNULL, check=True)
        os.replace(gff + ".tmp", gff)
    best = {}
    cur = None
    for line in open(gff):
        if line.startswith("##PAF"):
            f = line.rstrip("\n").split("\t")[1:]
            qlen, qs, qe, tname = int(f[1]), int(f[2]), int(f[3]), f[5]
            score = next((int(x[5:]) for x in f[12:] if x.startswith("AS:i:")), 0)
            cur = (f[0], tname, (qe - qs) / qlen, score)
        elif cur and "\tmRNA\t" in line:
            ident = next((float(x.split("=")[1]) for x in line.split("\t")[8].split(";") if x.startswith("Identity=")), 0.0)
            prot, tgt, cov, score = cur
            if cov >= 0.30 and score > best.get(tgt, (0,))[0]:
                best[tgt] = (score, prot, ident, cov)
            cur = None
    pd = [k for k in pseudo if k in Dg]
    reached = [k for k in pseudo if k in best]
    agree = checked = 0
    idents = []
    for k in pd:
        if k not in best:
            continue
        prot = best[k][1]
        idents.append(best[k][2])
        if prot in Dg:
            checked += 1
            agree += Dg[prot] == Dg[k]
    print(f"pseudogene loci {len(pseudo)}; reached by a family protein (cov >= 0.30) {len(reached)} "
          f"({len(reached) / max(1, len(pseudo)):.3f}); in a DNA family {len(pd)}, of which reached "
          f"{sum(1 for k in pd if k in best)}; attachment agrees with the DNA family {agree}/{checked}")
    if idents:
        idents.sort()
        print(f"  identity of attaching alignments: median {idents[len(idents) // 2]:.2f}, "
              f"< 0.50: {sum(1 for x in idents if x < 0.5)}, >= 0.70: {sum(1 for x in idents if x >= 0.7)}")


if __name__ == "__main__":
    main()

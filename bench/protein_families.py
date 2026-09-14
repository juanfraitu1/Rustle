#!/usr/bin/env python3
"""Prereg Addendum AN: protein-space multi-copy gene families, and cross-annotation scoring by CDS overlap.

build: protein_families.py build --nodes NODES_TSV --genome FA --contigs c1,c2 --out PREFIX [--min-ident 0.0] [--threads 4]
  One protein per gene (longest CDS from `<NODES_TSV>.cds.tsv`, >= 10 aa); all-vs-all `blastp -evalue 1e-5` (cached at
  PREFIX.blastp.tsv); per ordered pair, non-overlapping HSPs greedy by bitscore on the longer protein; EDGE iff their union
  covers >= 0.30 of the longer protein (and identity >= --min-ident); weight identity x coverage; MCL I=2.8 prune 1e-9.
  Writes PREFIX[.i<min-ident>].families.tsv (family_id, idx, name, biotype, chrom, strand, cds) and .edges.tsv.
score: protein_families.py score --truth NODES:FAMILIES --test NAME=NODES:FAMILIES ...
  Truth genes in families >= 2; each assigned the family of the test gene with the greatest CDS-base overlap.
"""
import argparse
import bisect
import collections
import csv
import os
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import adjudicated_truth as at  # noqa: E402
import guided_pipeline as gp  # noqa: E402
import mcl_port  # noqa: E402

BLAST = "/home/juanfra/miniforge3/envs/blast/bin"


def excluded(biotype, rule):
    """r1: pseudogenes; r2: r1 + V(D)J recombining antigen-receptor segments (AN revision r2)."""
    if rule >= 1 and "pseudogene" in biotype:
        return True
    if rule >= 2 and (any(x in biotype for x in ("V_segment", "D_segment", "J_segment", "C_region"))
                      or biotype.startswith(("IG_", "TR_"))):
        return True
    return False
MIN_COV = 0.30


def load_genes(nodes, contigs):
    names = {r["idx"]: (r["name"], r["biotype"]) for r in csv.DictReader(open(nodes + ".names.tsv"), delimiter="\t")}
    chrom = {r["idx"]: r["chrom"] for r in csv.DictReader(open(nodes), delimiter="\t")}
    genes = {}
    for r in csv.DictReader(open(nodes + ".cds.tsv"), delimiter="\t"):
        if chrom[r["idx"]] not in contigs:
            continue
        segs = [(int(x.split("-")[0]), int(x.split("-")[1].split(":")[0]), int(x.split(":")[1])) for x in r["cds"].split(",")]
        genes[r["idx"]] = {"idx": r["idx"], "chrom": chrom[r["idx"]], "strand": r["strand"], "cds": segs,
                           "name": names[r["idx"]][0], "biotype": names[r["idx"]][1]}
    return genes


def pair_hsps(path, plen):
    """(q, s) -> list of (bitscore, nident, length, q0, q1, s0, s1) from outfmt 6 lines."""
    hs = collections.defaultdict(list)
    for line in open(path):
        q, s, nid, ln, q0, q1, s0, s1, bits = line.rstrip("\n").split("\t")
        if q != s:
            hs[(q, s)].append((float(bits), int(nid), int(ln), int(q0) - 1, int(q1), int(s0) - 1, int(s1)))
    return hs


def edges_from(hs, plen, min_ident):
    best = {}
    for (q, s), rows in hs.items():
        longer_is_q = plen[q] >= plen[s]
        taken, L, N = [], 0, 0
        for bits, nid, ln, q0, q1, s0, s1 in sorted(rows, reverse=True):
            iv = (q0, q1) if longer_is_q else (s0, s1)
            if any(iv[0] < y and x < iv[1] for x, y in taken):
                continue
            taken.append(iv)
            L += ln
            N += nid
        cov = sum(y - x for x, y in gp.merge(taken)) / max(plen[q], plen[s])
        ident = N / L if L else 0.0
        if cov >= MIN_COV and ident >= min_ident:
            k = (min(q, s, key=int), max(q, s, key=int))
            w = ident * min(cov, 1.0)
            if w > best.get(k, (0.0,))[0]:
                best[k] = (w, ident, cov)
    return best


def cmd_build(a):
    contigs = set(a.contigs.split(","))
    genome = pysam.FastaFile(a.genome)
    genes = load_genes(a.nodes, contigs)
    rule = max(a.rule, 1 if a.no_pseudogenes else 0)
    genes = {k: g for k, g in genes.items() if not excluded(g["biotype"], rule)}
    faa = a.out + ".proteins.faa"
    plen = {}
    with open(faa, "w") as fh:
        for k, g in genes.items():
            p = at.translate(genome, g["chrom"], g["strand"], g["cds"])
            if len(p) >= 10:
                plen[k] = len(p)
                fh.write(f">{k}\n{p}\n")
    bl = a.out + ".blastp.tsv"
    if not os.path.exists(bl):
        subprocess.run([BLAST + "/makeblastdb", "-dbtype", "prot", "-in", faa, "-out", a.out + "_protdb"],
                       stdout=subprocess.DEVNULL, check=True)
        with open(bl + ".tmp", "w") as fh:
            subprocess.run([BLAST + "/blastp", "-query", faa, "-db", a.out + "_protdb", "-evalue", "1e-5",
                            "-max_target_seqs", "100000", "-num_threads", str(a.threads),
                            "-outfmt", "6 qseqid sseqid nident length qstart qend sstart send bitscore"], stdout=fh, check=True)
        os.replace(bl + ".tmp", bl)
    E = edges_from(pair_hsps(bl, plen), plen, a.min_ident)
    fams = [c for c in mcl_port.mcl({k: v[0] for k, v in E.items()}) if len(c) >= 2]
    fams.sort(key=len, reverse=True)
    tag = a.out if a.min_ident == 0 else f"{a.out}.i{a.min_ident:.2f}"
    with open(tag + ".families.tsv", "w") as fh:
        fh.write("family_id\tidx\tname\tbiotype\tchrom\tstrand\tcds\n")
        for i, c in enumerate(fams):
            for k in sorted(c, key=int):
                g = genes[k]
                fh.write(f"P{i}\t{k}\t{g['name']}\t{g['biotype']}\t{g['chrom']}\t{g['strand']}\t"
                         f"{','.join(f'{x}-{y}' for x, y, _ in g['cds'])}\n")
    with open(tag + ".edges.tsv", "w") as fh:
        fh.write("u\tv\tweight\tidentity\tcoverage\n")
        for (u, v), (w, i, c) in sorted(E.items(), key=lambda kv: (int(kv[0][0]), int(kv[0][1]))):
            fh.write(f"{u}\t{v}\t{w:.4f}\t{i:.4f}\t{c:.4f}\n")
    print(f"{tag}: proteins {len(plen)}; edges {len(E)}; families {len(fams)} ({sum(map(len, fams))} genes, "
          f"largest {len(fams[0]) if fams else 0})")


def load_fams(spec):
    nodes, fam = spec.split(":", 1)
    rows = list(csv.DictReader(open(fam), delimiter="\t"))
    return {r["idx"]: r["family_id"] for r in rows}, {
        r["idx"]: (r["chrom"], [tuple(map(int, b.split("-"))) for b in r["cds"].split(",")]) for r in rows}, nodes


def cmd_score(a):
    tfam, tcds, tnodes = load_fams(a.truth)
    cnt = collections.Counter(tfam.values())
    truth = sorted((k for k in tfam if cnt[tfam[k]] >= 2), key=int)
    contigs = {tcds[k][0] for k in truth}
    print(f"truth: {len(truth)} genes in {len({tfam[k] for k in truth})} families, "
          f"{sum(v * (v - 1) // 2 for v in cnt.values() if v >= 2)} pairs")
    print(f"{'catalog':14s} {'pair_sens':>9s} {'pair_prec':>9s} {'bip_R':>6s} {'bip_P':>6s} {'bip_F':>6s} unassigned")
    for spec in a.test:
        name, rest = spec.split("=", 1)
        nodes, fampath = rest.split(":", 1)
        tg = load_genes(nodes, contigs)
        rule = max(a.rule, 1 if a.no_pseudogenes else 0)
        tg = {k: g for k, g in tg.items() if not excluded(g["biotype"], rule)}
        pfam = {r["idx"]: r["family_id"] for r in csv.DictReader(open(fampath), delimiter="\t")}
        by = collections.defaultdict(list)
        for k, g in tg.items():
            segs = sorted((x, y) for x, y, _ in g["cds"])
            by[g["chrom"]].append((segs[0][0], segs[-1][1], k, segs))
        idx = {}
        for c, v in by.items():
            v.sort()
            idx[c] = (v, [x[0] for x in v], max(x[1] - x[0] for x in v))
        pred, un = [], 0
        for i, k in enumerate(truth):
            c, segs = tcds[k]
            s0, e0 = min(x for x, _ in segs), max(y for _, y in segs)
            v, starts, ml = idx.get(c, ([], [], 0))
            lo, hi = bisect.bisect_left(starts, s0 - ml), bisect.bisect_left(starts, e0)
            best = (0, None)
            for a0, a1, kk, tsegs in v[lo:hi]:
                if a1 <= s0:
                    continue
                ov = sum(max(0, min(y, d) - max(x, b)) for x, y in segs for b, d in tsegs)
                if ov > best[0]:
                    best = (ov, kk)
            f = pfam.get(best[1]) if best[1] else None
            if f is None:
                un += 1
            pred.append(f or f"none:{i}")
        true = [tfam[k] for k in truth]
        ps, pp = gp.pairwise(pred, true)
        br, bp = gp.bipartite(pred, true)
        f1 = 2 * br * bp / (br + bp) if br + bp else float("nan")
        print(f"{name:14s} {ps:9.3f} {pp:9.3f} {br:6.3f} {bp:6.3f} {f1:6.3f} {un}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("build")
    for k in ("--nodes", "--genome", "--contigs", "--out"):
        p.add_argument(k, required=True)
    p.add_argument("--min-ident", type=float, default=0.0)
    p.add_argument("--no-pseudogenes", action="store_true", help="AN r1: drop records whose biotype contains 'pseudogene'")
    p.add_argument("--rule", type=int, default=0, help="0 AN, 1 r1 (no pseudogenes), 2 r2 (r1 + V(D)J segments)")
    p.add_argument("--threads", type=int, default=4)
    p = sub.add_parser("score")
    p.add_argument("--truth", required=True)
    p.add_argument("--test", action="append", required=True)
    p.add_argument("--no-pseudogenes", action="store_true")
    p.add_argument("--rule", type=int, default=0)
    a = ap.parse_args()
    {"build": cmd_build, "score": cmd_score}[a.cmd](a)


if __name__ == "__main__":
    main()

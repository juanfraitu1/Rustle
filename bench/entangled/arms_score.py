#!/usr/bin/env python3
"""Window-wide chain recovery, artifact load and gene resolution of transcript sets on the ideal-expression reads (docs/PREREG_entangled_baseline_2026-10-06.md).

    arms_score.py --truth READSPREFIX --arm NAME=GTF [--arm NAME=GTF ...] --out OUTPREFIX

TRUTH = READSPREFIX.transcripts.tsv of bench/ideal_expression/sim_windows.py: every simulated transcript of every gene in the windows (10 full-length reads each), canonical chains.
A truth GENE is entangled (stratum E) iff its exon union shares >= 100 bp with the exon union of another truth gene on the same chromosome and strand (the rule of sim_windows.py); the rest is stratum N.
Only genes with >= 1 multi-exon chain enter the chain and resolution metrics (strata E and N); single-exon genes are counted in OUTPREFIX.json only for the entanglement flag.
An ARM is any GTF with exon records carrying transcript_id and gene_id (StringTie, FLAIR, the assembler's GTF, families.gtf).

Per arm and stratum:
  chains / recovered      distinct truth chains (chrom, strand, ordered intron list) and how many equal the chain of >= 1 arm transcript (exact).
  genes / complete / any  genes with a chain; genes with every chain recovered; genes with >= 1 chain recovered.
  resolved                a truth gene is RESOLVED iff some arm gene_id carries an exact chain of this gene and carries an exact chain of no other truth gene (physical overlap with an entangled partner is allowed).
  artifacts               arm multi-exon transcripts whose chain is no truth chain: fragment (the intron list is a contiguous part of a truth chain), fusion (exons overlap >= 2 truth genes on the strand), other.
  merged_ids              arm gene_ids that carry exact chains of >= 2 truth genes.
Writes OUTPREFIX.json (all numbers) and OUTPREFIX.genes.tsv (one row per arm and gene).
"""
import argparse
import bisect
import collections
import csv
import json
import re

MIN_SHARED = 100
BIN = 2000


def blocks(s):
    return [tuple(map(int, x.split("-"))) for x in s.split(",")] if s else []


def parse_chain(s):
    return tuple(tuple(map(int, x.split("-"))) for x in s.split(",")) if s else ()


def chain_of(exons):
    ex = sorted(exons)
    return tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))


def merge(bl):
    out = []
    for s, e in sorted(bl):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def overlap_bp(a, b):
    i = j = tot = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi > lo:
            tot += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return tot


def load_truth(prefix):
    genes = {}
    with open(prefix + ".transcripts.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["simulated"] != "1":
                continue
            g = genes.setdefault(r["gene"], dict(gene=r["gene"], chrom=r["chrom"], strand=r["strand"], name=r["gene_name"], biotype=r["biotype"],
                                                 copy=r["copy"], blocks=[], chains=set()))
            g["blocks"].extend(blocks(r["canon_blocks"]))
            if r["chain"]:
                g["chains"].add(parse_chain(r["chain"]))
    for g in genes.values():
        g["union"] = merge(g["blocks"])
    # entanglement
    by = collections.defaultdict(list)
    for g in genes.values():
        by[(g["chrom"], g["strand"])].append(g)
    for lst in by.values():
        lst.sort(key=lambda g: g["union"][0][0])
        for i, g in enumerate(lst):
            g["entangled"], g["partner"], g["shared_bp"] = 0, "", 0
            hi = g["union"][-1][1]
            for h in lst:
                if h is g or h["union"][0][0] >= hi or h["union"][-1][1] <= g["union"][0][0]:
                    continue
                ov = overlap_bp(g["union"], h["union"])
                if ov >= MIN_SHARED and ov > g["shared_bp"]:
                    g["entangled"], g["partner"], g["shared_bp"] = 1, h["gene"], ov
    return genes


class Index:
    """Gene exon blocks by chromosome, strand and 2 kb bin."""

    def __init__(self, genes):
        self.bins = collections.defaultdict(list)
        for g in genes.values():
            for s, e in g["union"]:
                for b in range(s // BIN, (e - 1) // BIN + 1):
                    self.bins[(g["chrom"], g["strand"], b)].append((s, e, g["gene"]))

    def genes_hit(self, chrom, strand, exons):
        hit = set()
        for s, e in exons:
            for b in range(s // BIN, (e - 1) // BIN + 1):
                for bs, be, g in self.bins.get((chrom, strand, b), ()):
                    if bs < e and be > s:
                        hit.add(g)
        return hit


def load_arm(path):
    tx = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            t = re.search(r'transcript_id "([^"]+)"', f[8])
            g = re.search(r'gene_id "([^"]+)"', f[8])
            if not t:
                continue
            d = tx.setdefault(t.group(1), dict(chrom=f[0], strand=f[6], gene=g.group(1) if g else t.group(1), exons=[]))
            d["exons"].append((int(f[3]) - 1, int(f[4])))
    return tx


def score_arm(name, tx, genes, idx, truth_chain_set, intron_index):
    arm_chains = set()
    multi = {}
    for t, d in tx.items():
        if len(d["exons"]) >= 2:
            ch = chain_of(d["exons"])
            multi[t] = ch
            arm_chains.add((d["chrom"], d["strand"], ch))
    # per arm gene_id: truth genes overlapped, exact chains carried
    id_genes = collections.defaultdict(set)
    id_exact = collections.defaultdict(set)   # arm gene_id -> truth genes with an exact chain in some transcript
    chain_owner = collections.defaultdict(set)
    for g in genes.values():
        for c in g["chains"]:
            chain_owner[(g["chrom"], g["strand"], c)].add(g["gene"])
    for t, d in tx.items():
        hit = idx.genes_hit(d["chrom"], d["strand"], d["exons"])
        id_genes[d["gene"]] |= hit
        if t in multi:
            for owner in chain_owner.get((d["chrom"], d["strand"], multi[t]), ()):
                id_exact[d["gene"]].add(owner)
    ids_for_gene = collections.defaultdict(set)
    for i, owners in id_exact.items():
        for o in owners:
            ids_for_gene[o].add(i)
    # artifacts
    art = collections.Counter()
    for t, ch in multi.items():
        d = tx[t]
        key = (d["chrom"], d["strand"], ch)
        if key in truth_chain_set:
            continue
        hit = idx.genes_hit(d["chrom"], d["strand"], d["exons"])
        sub = False
        first = ch[0]
        for (cc, ci, pos) in intron_index.get((d["chrom"], d["strand"], first), ()):
            if cc[pos:pos + len(ch)] == ch:
                sub = True
                break
        art["fragment" if sub else "fusion" if len(hit) >= 2 else "other"] += 1
    mono = sum(1 for d in tx.values() if len(d["exons"]) == 1)
    rows, summ = [], {}
    for stratum in ("E", "N"):
        gs = [g for g in genes.values() if g["chains"] and (g["entangled"] == (stratum == "E"))]
        chains = rec = comp = anyr = res = 0
        for g in gs:
            got = [c for c in g["chains"] if (g["chrom"], g["strand"], c) in arm_chains]
            clean = any(id_exact[i] == {g["gene"]} for i in ids_for_gene.get(g["gene"], ()))
            chains += len(g["chains"]); rec += len(got); comp += int(len(got) == len(g["chains"])); anyr += int(len(got) > 0); res += int(clean)
            rows.append(dict(arm=name, gene=g["gene"], name=g["name"], copy=g["copy"], stratum=stratum, partner=g["partner"], shared_bp=g["shared_bp"],
                             chains=len(g["chains"]), recovered=len(got), resolved=int(clean)))
        summ[stratum] = dict(genes=len(gs), chains=chains, recovered=rec, complete=comp, any=anyr, resolved=res)
    merged = sum(1 for i, s in id_exact.items() if len(s) >= 2)
    return dict(transcripts=len(tx), multi_exon=len(multi), mono_exon=mono, arm_gene_ids=len(id_genes), merged_ids=merged, artifacts=dict(art),
                artifact_total=sum(art.values()), strata=summ), rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--arm", action="append", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    genes = load_truth(a.truth)
    idx = Index(genes)
    truth_chain_set = {(g["chrom"], g["strand"], c) for g in genes.values() for c in g["chains"]}
    intron_index = collections.defaultdict(list)
    for g in genes.values():
        for c in g["chains"]:
            for pos, intr in enumerate(c):
                intron_index[(g["chrom"], g["strand"], intr)].append((c, g["gene"], pos))
    out = dict(truth=dict(genes=len(genes), genes_with_chain=sum(1 for g in genes.values() if g["chains"]), entangled=sum(1 for g in genes.values() if g["entangled"] and g["chains"]),
                          chains=len(truth_chain_set)), arms={})
    allrows = []
    for spec in a.arm:
        name, path = spec.split("=", 1)
        res, rows = score_arm(name, load_arm(path), genes, idx, truth_chain_set, intron_index)
        out["arms"][name] = res
        allrows.extend(rows)
        e, n = res["strata"]["E"], res["strata"]["N"]
        print(f"{name:7} tx {res['transcripts']:6} multi {res['multi_exon']:6} | E: chains {e['recovered']}/{e['chains']} genes complete {e['complete']}/{e['genes']} resolved {e['resolved']}/{e['genes']}"
              f" | N: chains {n['recovered']}/{n['chains']} complete {n['complete']}/{n['genes']} resolved {n['resolved']}/{n['genes']} | artifacts {res['artifact_total']} {res['artifacts']} merged ids {res['merged_ids']}")
    json.dump(out, open(a.out + ".json", "w"), indent=1)
    with open(a.out + ".genes.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(allrows)


if __name__ == "__main__":
    main()

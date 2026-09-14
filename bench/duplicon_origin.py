#!/usr/bin/env python3
"""Prereg Addendum W: did a multi-copy gene family arise from duplicons? For each family, every member's region (gene
span +/- one gene-span length) is aligned to a reference member's region; a pair whose alignment chain continues
>= 1 kb past the gene on a matched side in BOTH members co-duplicated its flank (DUPLICON pair), otherwise the homology
stops at the gene (GENE-ONLY pair). A family is DUPLICON-BORNE when most of its aligned pairs co-duplicated flank.

usage: duplicon_origin.py --gff full.gff.gz --genome genome.fa --outdir DIR [--amy-truth truth.tsv] [--threads 4]
"""
import argparse
import collections
import csv
import gzip
import math
import os
import re
import statistics
import subprocess

import pysam

MIN_ID, EXT, CAP, CORE_FRAC = 0.80, 1000, 40, 0.90
FAMILIES = {
    "NPIP": (r"^NPIP[AB]\d+P?$", "core"), "TBC1D3": (r"^TBC1D3[A-Z]?$", "core"),
    "GOLGA8": (r"^GOLGA8[A-Z]*P?\d*$", "core"), "LRRC37A": (r"^LRRC37A\d*P?$", "core"), "RGPD": (r"^RGPD\d+$", "core"),
    "NBPF": (r"^NBPF\d+P?$", "core"), "SPATA31": (r"^SPATA31[A-Z]\d*P?\d*$", "core"), "PMS2": (r"^PMS2(P\d+)?$", "core"),
    "TRIM51": (r"^TRIM51[A-Z]*P?\d*$", "core"), "GUSB": (r"^GUSB(P\d+)?$", "core"),
    "GAPDH": (r"^GAPDH(P\d+)?$", "retro"), "PPIA": (r"^PPIA(P\d+)?$", "retro"), "EEF1A1": (r"^EEF1A1(P\d+)?$", "retro"),
}


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def load_genes(gff):
    genes, tx_parent, exons = {}, {}, collections.defaultdict(int)
    with gzip.open(gff, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] in ("gene", "pseudogene"):
                a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                genes[a["ID"]] = {"name": a.get("Name", "?"), "chrom": f[0], "start0": int(f[3]) - 1, "end": int(f[4]),
                                  "strand": f[6], "biotype": a.get("gene_biotype", f[2])}
            elif f[2] in ("mRNA", "transcript", "lnc_RNA", "ncRNA", "primary_transcript") or f[2] == "exon":
                a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                if f[2] == "exon":
                    exons[a.get("Parent", "")] += 1
                else:
                    tx_parent[a["ID"]] = a.get("Parent", "")
    introns = collections.defaultdict(int)
    for tx, p in tx_parent.items():
        introns[p] = max(introns[p], exons.get(tx, 0) - 1)
    for gid, g in genes.items():
        g["introns"] = max(introns.get(gid, 0), exons.get(gid, 0) - 1, 0)
    return genes


def chains_from_paf(lines):
    by = collections.defaultdict(list)
    for line in lines:
        f = line.split("\t")
        by[(f[0], f[4])].append({"q": f[0], "qlen": int(f[1]), "qs": int(f[2]), "qe": int(f[3]), "strand": f[4],
                                 "ts": int(f[7]), "te": int(f[8]), "nm": int(f[9]), "bl": int(f[10])})
    chains = []
    for (q, strand), rs in by.items():
        rs.sort(key=lambda r: (r["ts"], r["te"]))
        Lq, cur = rs[0]["qlen"], []
        for r in rs + [None]:
            ok = False
            if r is not None and cur:
                gap = r["ts"] - max(x["te"] for x in cur)
                span = r["te"] - min(x["ts"] for x in cur)
                order = r["qs"] >= cur[-1]["qs"] if strand == "+" else r["qs"] <= cur[-1]["qs"]
                ok = gap <= Lq and span <= 2 * Lq and order
            if r is None or (cur and not ok):
                chains.append({"q": q, "strand": strand, "recs": cur, "nm": sum(x["nm"] for x in cur),
                               "bl": sum(x["bl"] for x in cur), "qs": min(x["qs"] for x in cur), "qe": max(x["qe"] for x in cur),
                               "ts": min(x["ts"] for x in cur), "te": max(x["te"] for x in cur),
                               "aligned": sum(e - s for s, e in merge([(x["qs"], x["qe"]) for x in cur]))})
                cur = []
            if r is not None:
                cur.append(r)
    return chains


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--amy-truth")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    genome = pysam.FastaFile(a.genome)
    clen = dict(zip(genome.references, genome.lengths))
    genes = load_genes(a.gff)
    by_name = collections.defaultdict(list)
    for g in genes.values():
        by_name[g["name"]].append(g)

    fams = {}
    for fam, (pat, kind) in FAMILIES.items():
        members = [g for n, gs in by_name.items() if re.match(pat, n) for g in gs if g["chrom"] in clen]
        fams[fam] = (kind, members)
    if a.amy_truth:
        amy = []
        for r in csv.DictReader(open(a.amy_truth), delimiter="\t"):
            gs = [g for g in by_name.get(r["name"], []) if g["chrom"] == r["chrom"]]
            base = gs[0] if gs else {"introns": 0}
            amy.append({"name": r["name"], "chrom": r["chrom"], "start0": int(r["start0"]), "end": int(r["end"]),
                        "strand": r["strand"], "biotype": r["biotype"], "introns": base["introns"]})
        fams["AMY"] = ("reported", amy)

    pair_rows, fam_rows = [], []
    for fam, (kind, members) in fams.items():
        pcs = [g for g in members if g["biotype"] == "protein_coding"] or members
        med = statistics.median(g["end"] - g["start0"] for g in pcs)
        ref = min(pcs, key=lambda g: (abs((g["end"] - g["start0"]) - med), g["name"]))
        others = sorted([g for g in members if g is not ref], key=lambda g: (g["name"], g["chrom"], g["start0"]))
        dropped = max(0, len(others) - CAP)
        others = others[:CAP]

        def region(g):
            L = g["end"] - g["start0"]
            s, e = max(0, g["start0"] - L), min(clen[g["chrom"]], g["end"] + L)
            return s, e

        rs, re_ = region(ref)
        tfa, qfa = f"{a.outdir}/{fam}.ref.fa", f"{a.outdir}/{fam}.members.fa"
        with open(tfa, "w") as fh:
            fh.write(f">ref\n{genome.fetch(ref['chrom'], rs, re_).upper()}\n")
        keys = {}
        with open(qfa, "w") as fh:
            for i, g in enumerate(others):
                s, e = region(g)
                keys[f"m{i}"] = (g, s, e)
                fh.write(f">m{i}\n{genome.fetch(g['chrom'], s, e).upper()}\n")
        paf = subprocess.run(["minimap2", "-c", "-x", "asm20", "-N", "50", "-p", "0.1", "-t", str(a.threads), tfa, qfa],
                             capture_output=True, text=True, check=True).stdout.splitlines()
        chains = chains_from_paf(paf)
        gt0, gt1 = ref["start0"] - rs, ref["end"] - rs
        cov = collections.Counter()
        n_aligned = n_dup = 0
        for key, (g, s, e) in keys.items():
            gq0, gq1 = g["start0"] - s, g["end"] - s
            cand = [c for c in chains if c["q"] == key and c["qs"] < gq1 and gq0 < c["qe"] and c["ts"] < gt1 and gt0 < c["te"]]
            row = {"family": fam, "kind": kind, "reference": ref["name"], "member": g["name"], "chrom": g["chrom"],
                   "member_introns": g["introns"], "ref_introns": ref["introns"], "identity": float("nan"),
                   "ext_q_left": 0, "ext_q_right": 0, "ext_t_left": 0, "ext_t_right": 0, "class": "UNALIGNED"}
            if cand:
                c = max(cand, key=lambda c: c["aligned"])
                ident = c["nm"] / c["bl"]
                row["identity"] = ident
                if ident >= MIN_ID:
                    eql, eqr = gq0 - c["qs"], c["qe"] - gq1
                    etl, etr = gt0 - c["ts"], c["te"] - gt1
                    row.update(ext_q_left=eql, ext_q_right=eqr, ext_t_left=etl, ext_t_right=etr)
                    matched = [(eql, etl), (eqr, etr)] if c["strand"] == "+" else [(eql, etr), (eqr, etl)]
                    dup = any(x >= EXT and y >= EXT for x, y in matched)
                    row["class"] = "DUPLICON" if dup else "GENE-ONLY"
                    n_aligned += 1
                    n_dup += dup
                    for s_, e_ in merge([(r["ts"], r["te"]) for r in c["recs"]]):
                        for p in range(s_, e_):
                            cov[p] += 1
                else:
                    row["class"] = "LOW-IDENTITY"
            pair_rows.append(row)
        need = math.ceil(CORE_FRAC * n_aligned) if n_aligned else 0
        core = sum(1 for v in cov.values() if v >= need) if n_aligned else 0
        call = "DUPLICON-BORNE" if n_aligned >= 2 and n_dup / n_aligned >= 0.5 else "NOT"
        fam_rows.append({"family": fam, "kind": kind, "members": len(members), "dropped_by_cap": dropped,
                         "reference": ref["name"], "ref_span": ref["end"] - ref["start0"], "aligned": n_aligned,
                         "duplicon_pairs": n_dup, "dup_frac": n_dup / n_aligned if n_aligned else float("nan"),
                         "core_bp": core, "core_over_gene": core / (ref["end"] - ref["start0"]), "call": call})

    with open(f"{a.outdir}/pairs.tsv", "w") as fh:
        cols = list(pair_rows[0].keys())
        fh.write("\t".join(cols) + "\n")
        for r in pair_rows:
            fh.write("\t".join(f"{r[k]:.4f}" if isinstance(r[k], float) else str(r[k]) for k in cols) + "\n")
    print("family\tkind\tmembers\tdropped\treference\tref_span\taligned\tduplicon_pairs\tdup_frac\tcore_bp\tcore/gene\tcall")
    for r in fam_rows:
        print(f"{r['family']}\t{r['kind']}\t{r['members']}\t{r['dropped_by_cap']}\t{r['reference']}\t{r['ref_span']}\t"
              f"{r['aligned']}\t{r['duplicon_pairs']}\t{r['dup_frac']:.3f}\t{r['core_bp']}\t{r['core_over_gene']:.2f}\t{r['call']}")
    core_dup = sum(1 for r in fam_rows if r["kind"] == "core" and r["call"] == "DUPLICON-BORNE")
    retro_dup = sum(1 for r in fam_rows if r["kind"] == "retro" and r["call"] == "DUPLICON-BORNE")
    n_core = sum(1 for r in fam_rows if r["kind"] == "core")
    n_retro = sum(1 for r in fam_rows if r["kind"] == "retro")
    verdict = ("SUPPORTED" if core_dup >= 7 and retro_dup <= 1 else
               "UNINFORMATIVE (classifier calls retro families duplicon-borne)" if retro_dup >= 2 else "NOT SUPPORTED")
    print(f"\nREADING: core-duplicon families DUPLICON-BORNE {core_dup}/{n_core}; retrotransposition families DUPLICON-BORNE "
          f"{retro_dup}/{n_retro} -> {verdict}")
    retro_sig = collections.Counter((r["family"], r["member_introns"] == 0 and r["ref_introns"] >= 1) for r in pair_rows)
    print("retro signature (member 0 introns, reference >= 1) per family:",
          {f: f"{retro_sig[(f, True)]}/{retro_sig[(f, True)] + retro_sig[(f, False)]}" for f in fams})


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Ideal-expression read simulation around the copies of one family (docs/PREREG_ideal_expression_2026-10-06.md, substrate S1).

    sim_windows.py --family NPIP|TBC1D3 --out PREFIX [--pad 500000] [--reads-per-tx 10] [--seed 20261006]

Windows: +-pad bp around the territory (terr_lo0..terr_hi) of every CAT/Liftoff copy of the family (copy_recovery_tools_cat/ann/copies.hsa.tsv), merged.
Transcripts: every transcript of every gene (any biotype) whose span overlaps a window. Canonicalization (the 2026-09-18 recipe): an intron is kept only when it is
>= 50 bp and canonical (GT-AG, GC-AG, AT-AC on the transcript's strand, read from the genome); every other gap is merged into its flanking exons. A transcript is
SIMULATED iff its spliced length is >= 120 bp, its molecule is <= 30 kb and no kept intron exceeds 200 kb (minimap2's default maximum intron length); the others are listed.
Reads: `reads_per_tx` per simulated transcript, the model of `bench/sim.py chromosome --arm ideal` (simulate_reads: err .001, indel err/3, no truncation; 0-30 bp end jitter,
skipped when the body would fall under 100 bp), each read from sim.stable_seed(seed, transcript, k); names `transcript|gene|fl|k`.
Writes PREFIX.fq, PREFIX.transcripts.tsv, PREFIX.targets.tsv (the family's copies plus the held-out families' genes, with strata), PREFIX.named_family.tsv, PREFIX.skipped.tsv,
PREFIX.windows.tsv, PREFIX.sim.json. Coordinates 0-based half-open; an intron is (donor_end, acceptor_start).
"""
import argparse
import collections
import csv
import json
import os
import random
import re
import sys

import pysam

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
import lib  # noqa: E402
import sim  # noqa: E402

CAT_GFF = "/mnt/linuxdisk/home/juanfraitu/winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.slim.gff3.gz"
GENOME = "/mnt/linuxdisk/home/juanfraitu/winloci_data/chm13v2.0.fa"
ANN = "/mnt/linuxdisk/tmp/rustle_figures_dev/copy_recovery_tools_cat/ann"
CHROM = {"NPIP": "chr16", "TBC1D3": "chr17"}
HELDOUT = {"NPIP": ("SMG1P", re.compile(r"^SMG1P[0-9]+$")), "TBC1D3": ("KRTAP", re.compile(r"^KRTAP[0-9]+-[0-9]+$"))}
NAMED = re.compile(r"NPIP|TBC1D3")
MIN_INTRON = 50
MAX_MOLECULE = 30_000
MAX_INTRON = 200_000
SHARED_BP = 100
PLUS_OK = {("GT", "AG"), ("GC", "AG"), ("AT", "AC")}
MINUS_OK = {("CT", "AC"), ("CT", "GC"), ("GT", "AT")}   # genomic (first 2, last 2 bases) of the reverse complements


def canon(fa, chrom, strand, blocks):
    """blocks: sorted [(s0, e)] 0-based half-open. Returns (merged blocks, chain [(donor_end, acceptor_start)], n_noncanonical, n_short)."""
    out = [list(blocks[0])]
    chain, n_nc, n_short = [], 0, 0
    for s, e in blocks[1:]:
        pe = out[-1][1]
        keep = False
        if s - pe >= MIN_INTRON:
            first, last = fa.fetch(chrom, pe, pe + 2).upper(), fa.fetch(chrom, s - 2, s).upper()
            keep = (first, last) in (PLUS_OK if strand == "+" else MINUS_OK)
            if not keep:
                n_nc += 1
        else:
            n_short += 1
        if keep:
            chain.append((pe, s))
            out.append([s, e])
        else:
            out[-1][1] = e
    return [tuple(b) for b in out], chain, n_nc, n_short


def load_cat(chrom):
    """gene -> dict(span, strand, name, biotype), transcript -> (gene, strand, sorted 0-based exon blocks) for one chromosome of the CAT slim GFF."""
    genes, tx_gene, tx_strand, ex = {}, {}, {}, collections.defaultdict(list)
    for ln in pysam.TabixFile(CAT_GFF).fetch(chrom):
        f = ln.split("\t")
        at = dict(kv.split("=", 1) for kv in f[8].strip().split(";") if "=" in kv)
        if f[2] == "gene":
            genes[at["ID"]] = dict(span=(int(f[3]) - 1, int(f[4])), strand=f[6], name=at.get("gene_name", ""), biotype=at.get("gene_biotype", ""))
        elif f[2] == "transcript":
            tx_gene[at["ID"]] = at["Parent"]; tx_strand[at["ID"]] = f[6]
        elif f[2] == "exon":
            ex[at["Parent"]].append((int(f[3]) - 1, int(f[4])))
    tx = {t: (tx_gene[t], tx_strand[t], sorted(ex[t])) for t in tx_gene if ex.get(t)}
    return genes, tx


def merge(iv):
    out = []
    for s, e in sorted(iv):
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


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--family", required=True, choices=sorted(CHROM))
    ap.add_argument("--out", required=True)
    ap.add_argument("--pad", type=int, default=500_000)
    ap.add_argument("--reads-per-tx", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20261006)
    ap.add_argument("--err", type=float, default=0.001)
    ap.add_argument("--jitter", type=int, default=30)
    ap.add_argument("--no-reads", action="store_true", help="write the truth tables only (no FASTQ)")
    a = ap.parse_args()
    chrom = CHROM[a.family]
    copies = [r for r in csv.DictReader(open(f"{ANN}/copies.hsa.tsv"), delimiter="\t") if r["family"] == a.family and r["chrom"] == chrom]
    wins = []
    for lo, hi in sorted((int(r["terr_lo0"]), int(r["terr_hi"])) for r in copies):
        w = [max(0, lo - a.pad), hi + a.pad]
        if wins and w[0] <= wins[-1][1]:
            wins[-1][1] = max(wins[-1][1], w[1])
        else:
            wins.append(w)
    genes, tx = load_cat(chrom)
    sel = {g for g, v in genes.items() if any(v["span"][0] < w[1] and v["span"][1] > w[0] for w in wins)}
    fa = pysam.FastaFile(GENOME)
    ho_label, ho_re = HELDOUT[a.family]
    gene_cid = {(c["isoform_gene"] if c["isoform_gene"] in genes else c["cid"]): c["cid"] for c in copies}
    for g in sel:   # held-out family genes (no bars): copy id = family:gene id
        if ho_re.match(genes[g]["name"]) and g not in gene_cid:
            gene_cid[g] = f"{ho_label}:{g}"
    rows, skipped, gene_blocks = [], [], collections.defaultdict(list)
    for t, (g, st, blocks) in sorted(tx.items(), key=lambda kv: (kv[1][2][0][0], kv[0])):
        if g not in sel:
            continue
        cb, chain, n_nc, n_short = canon(fa, chrom, st, blocks)
        spliced = sum(e - s for s, e in cb)
        reason = ""
        if spliced < 120:
            reason = "under 120 bp"
        elif spliced > MAX_MOLECULE:
            reason = f"molecule {spliced} bp > {MAX_MOLECULE}"
        elif any(x - d > MAX_INTRON for d, x in chain):
            reason = "kept intron > 200 kb"
        rows.append(dict(transcript=t, gene=g, chrom=chrom, strand=st, raw_blocks=",".join(f"{s}-{e}" for s, e in blocks),
                         canon_blocks=",".join(f"{s}-{e}" for s, e in cb), chain=",".join(f"{d}-{x}" for d, x in chain),
                         n_raw_introns=len(blocks) - 1, n_noncanonical=n_nc, n_short=n_short, spliced_len=spliced,
                         simulated=int(not reason), truth=int(not reason), skip_reason=reason, copy=gene_cid.get(g, ""), gene_name=genes[g]["name"], biotype=genes[g]["biotype"]))
        if reason:
            skipped.append(dict(transcript=t, gene=g, gene_name=genes[g]["name"], reason=reason, spliced_len=spliced))
        else:
            gene_blocks[g].extend(cb)
    union = {g: merge(b) for g, b in gene_blocks.items()}
    with open(a.out + ".transcripts.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(rows)
    with open(a.out + ".skipped.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["transcript", "gene", "gene_name", "reason", "spliced_len"], delimiter="\t"); w.writeheader(); w.writerows(skipped)
    # named-family genes: simulated non-target genes whose gene_name mentions the family (expected non-target members)
    with open(a.out + ".named_family.tsv", "w") as fh:
        fh.write("gene\tgene_name\tbiotype\tstrand\tchrom\tcanon_blocks\ttarget_cid\n")
        for g in sorted(union):
            if NAMED.search(genes[g]["name"]):
                fh.write(f"{g}\t{genes[g]['name']}\t{genes[g]['biotype']}\t{genes[g]['strand']}\t{chrom}\t" + ",".join(f"{s}-{e}" for s, e in union[g]) + f"\t{gene_cid.get(g, '')}\n")
    # targets (the family's copies, then the held-out family) with strata
    tg = []
    cid_info = {c["cid"]: c for c in copies}
    for g, cid in gene_cid.items():
        c = cid_info.get(cid)
        fam = a.family if c else ho_label
        strand = c["strand"] if c else genes[g]["strand"]
        name = (c.get("cat_name") or c.get("refseq_name") or c["name"]) if c else genes[g]["name"]
        rws = [r for r in rows if r["gene"] == g and r["simulated"] == 1]
        chains = {r["chain"] for r in rws}
        multi = [x for x in chains if x]
        best, best_frac = "", 0.0
        if g in union:
            tot = sum(e - s for s, e in union[g])
            for g2, u2 in union.items():
                if g2 != g and genes[g2]["strand"] == genes[g]["strand"]:
                    ov = overlap_bp(union[g], u2)
                    if ov >= SHARED_BP and ov > best_frac * tot:
                        best, best_frac = g2, ov / tot
        shared_bp = int(round(best_frac * sum(e - s for s, e in union[g]))) if best else 0
        is_E = int(bool(best)); is_C = int(len(multi) == 0); is_X = int(bool(c and c["readthrough"] == "1"))
        tg.append(dict(cid=cid, name=name, gene=g, family=fam, heldout=int(c is None), chrom=chrom, strand=strand, source=genes[g].get("biotype", ""),
                       n_tx=len(rws), n_distinct_chains=len(chains), n_multiexon_chains=len(multi), exon_bp=sum(e - s for s, e in union.get(g, [])),
                       entangled=is_E, entangled_with=best, shared_bp=shared_bp, shared_frac=round(best_frac, 3), mono=is_C, rt_image=is_X,
                       stratum=("E" if is_E else "C" if is_C else "X" if is_X else "R_in")))
    with open(a.out + ".targets.tsv", "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tg[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(sorted(tg, key=lambda r: (r["heldout"], r["cid"])))
    with open(a.out + ".windows.tsv", "w") as fh:
        fh.write("chrom\tstart0\tend\n"); [fh.write(f"{chrom}\t{x}\t{y}\n") for x, y in wins]
    n = 0
    if not a.no_reads:
        with open(a.out + ".fq", "w") as fq:
            for r in rows:
                if not r["simulated"]:
                    continue
                t, g, st = r["transcript"], r["gene"], r["strand"]
                cb = [tuple(map(int, b.split("-"))) for b in r["canon_blocks"].split(",")]
                body0 = lib.spliced1(fa, chrom, [(s + 1, e) for s, e in cb], st)
                for k in range(a.reads_per_tx):
                    rng = random.Random(sim.stable_seed(str(a.seed), t, str(k), "jitter"))
                    lo, hi = rng.randint(0, a.jitter), rng.randint(0, a.jitter)
                    body = body0[lo:len(body0) - hi] if len(body0) - hi > lo + 100 else body0
                    (rd, q), = sim.simulate_reads(body, 1, err=a.err, indel=a.err / 3, seed=sim.stable_seed(str(a.seed), t, str(k), "errors"), trunc_frac=0.0)
                    sim.write_fastq(fq, f"{t}|{g}|fl|{k}", (rd, q)); n += 1
    main_t = [t for t in tg if not t["heldout"]]
    meta = dict(family=a.family, chrom=chrom, windows=wins, pad=a.pad, genes=len(sel), transcripts=len(rows), simulated_transcripts=sum(r["simulated"] for r in rows),
                skipped=collections.Counter(s["reason"].split(" ")[0] + " " + s["reason"].split(" ")[1] if " " in s["reason"] else s["reason"] for s in skipped),
                reads=n, seed=a.seed, err=a.err, jitter=a.jitter, reads_per_tx=a.reads_per_tx, target_copies=len(main_t),
                strata={k: sum(1 for t in main_t if t["stratum"] == k) for k in ("E", "C", "X", "R_in")}, heldout_genes=len(tg) - len(main_t),
                longest_molecule=max((r["spliced_len"] for r in rows if r["simulated"]), default=0))
    json.dump(meta, open(a.out + ".sim.json", "w"), indent=1)
    print(json.dumps(meta))


if __name__ == "__main__":
    main()

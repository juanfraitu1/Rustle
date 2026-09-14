#!/usr/bin/env python3
"""Prereg Addendum W: simulated IsoSeq reads from every amylase locus, aligned to the genome, so de novo O1 can be
tested with expression removed as the limit.

usage: amy_sim.py --truth truth.tsv --gff refseq.gff --genome genome.fa --mmi genome.mmi --outdir DIR --depth N [--threads 4]
Writes DIR/sources.fa (source transcripts), DIR/sim.fq, DIR/sim.bam(.bai), DIR/placement.tsv (read-level placement).
"""
import argparse
import collections
import csv
import os
import re
import subprocess
import sys

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_reads import simulate_reads, write_fastq  # noqa: E402

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--mmi", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--depth", type=int, required=True)
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    genome = pysam.FastaFile(a.genome)
    truth = list(csv.DictReader(open(a.truth), delimiter="\t"))
    for t in truth:
        t["start0"], t["end"] = int(t["start0"]), int(t["end"])
    rec = {t["name"]: t for t in truth}

    gene_id, tx_parent, tx_name, exons = {}, {}, {}, collections.defaultdict(list)
    for line in open(a.gff):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        if f[2] in ("gene", "pseudogene") and at.get("Name") in rec:
            gene_id[at["ID"]] = at["Name"]
        elif f[2] in ("mRNA", "transcript", "lnc_RNA", "ncRNA", "primary_transcript"):
            tx_parent[at["ID"]], tx_name[at["ID"]] = at.get("Parent", ""), at.get("Name", at["ID"])
        elif f[2] == "exon":
            exons[at.get("Parent", "")].append((int(f[3]) - 1, int(f[4])))

    def spliced(t, blocks):
        s = "".join(genome.fetch(t["chrom"], x, y) for x, y in merge(blocks)).upper()
        return s.translate(COMP)[::-1] if t["strand"] == "-" else s

    sources = collections.defaultdict(list)  # locus -> [(tx_name, seq)]
    for gid, name in gene_id.items():
        t = rec[name]
        for tx, p in sorted(tx_parent.items()):
            if p != gid or not exons.get(tx):
                continue
            blocks = [(max(x, t["start0"]), min(y, t["end"])) for x, y in exons[tx] if min(y, t["end"]) > max(x, t["start0"])]
            if name == "LOC124905662":  # v2: exons 2-8 only (the truth span clips the mis-joined exon 1 away)
                pass
            if blocks:
                sources[name].append((tx_name[tx], spliced(t, blocks)))
    if "AMYP1" in rec and not sources.get("AMYP1"):
        t = rec["AMYP1"]
        amy2a = dict(sources.get("AMY2A", []))
        q = amy2a.get("NM_000699.4")
        region = genome.fetch(t["chrom"], t["start0"], t["end"]).upper()
        with open(f"{a.outdir}/amyp1_q.fa", "w") as fh:
            fh.write(f">NM_000699.4\n{q}\n")
        with open(f"{a.outdir}/amyp1_t.fa", "w") as fh:
            fh.write(f">AMYP1\n{region}\n")
        paf = subprocess.run(["minimap2", "-c", "-x", "splice", f"{a.outdir}/amyp1_t.fa", f"{a.outdir}/amyp1_q.fa"],
                             capture_output=True, text=True, check=True).stdout.splitlines()
        best = max((l.split("\t") for l in paf), key=lambda f: int(f[9]), default=None)
        blocks = []
        if best and int(best[9]) / int(best[10]) >= 0.80:
            pos, cur = int(best[7]), int(best[7])
            cg = next(x[5:] for x in best[12:] if x.startswith("cg:Z:"))
            for n, op in re.findall(r"(\d+)([MIDN])", cg):
                n = int(n)
                if op in "MD":
                    pos += n
                elif op == "N":
                    blocks.append((cur, pos))
                    pos += n
                    cur = pos
            blocks.append((cur, pos))
            seq = "".join(region[x:y] for x, y in blocks)
            if best[4] == "-":
                seq = seq.translate(COMP)[::-1]
            model = f"AMY2A-projected ({len(blocks)} exons, identity {int(best[9]) / int(best[10]):.3f})"
        else:
            seq = region.translate(COMP)[::-1] if t["strand"] == "-" else region
            model = "unspliced span"
        sources["AMYP1"].append((model, seq))
        print(f"[sim] AMYP1 model: {model}, {len(seq)} bp", file=sys.stderr)

    with open(f"{a.outdir}/sources.fa", "w") as fs, open(f"{a.outdir}/sim.fq", "w") as fq:
        for li, name in enumerate(sorted(rec)):
            txs = sources.get(name, [])
            if not txs:
                print(f"[sim] {name}: no source transcript", file=sys.stderr)
                continue
            per = [a.depth // len(txs) + (1 if i < a.depth % len(txs) else 0) for i in range(len(txs))]
            for ti, ((tname, seq), n) in enumerate(zip(txs, per)):
                fs.write(f">{name}|{tname}\n{seq}\n")
                for i, rq in enumerate(simulate_reads(seq, n, err=0.003, indel=0.0008, seed=1000 * li + ti + 1, trunc_frac=0.3)):
                    if len(rq[0]) >= 300:
                        write_fastq(fq, f"SIMAMY|{name}|{tname.split()[0]}|{i}", rq)
    bam = f"{a.outdir}/sim.bam"
    subprocess.run(f"minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t {a.threads} {a.mmi} {a.outdir}/sim.fq "
                   f"2>/dev/null | samtools sort -o {bam} - && samtools index {bam}", shell=True, check=True)

    # read-level placement of primaries
    counts = collections.defaultdict(collections.Counter)
    with pysam.AlignmentFile(bam) as B, open(f"{a.outdir}/placement.tsv", "w") as fh:
        fh.write("read\tsource\tprimary_locus\tmapq\n")
        for r in B.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:
                continue
            src = r.query_name.split("|")[1]
            if r.is_unmapped:
                counts[src]["unmapped"] += 1
                fh.write(f"{r.query_name}\t{src}\tunmapped\t0\n")
                continue
            hit = [n for n, t in rec.items() if t["chrom"] == r.reference_name and r.reference_start < t["end"] and t["start0"] < r.reference_end]
            where = max(hit, key=lambda n: min(rec[n]["end"], r.reference_end) - max(rec[n]["start0"], r.reference_start)) if hit else "elsewhere"
            counts[src]["on_source" if where == src else ("other_amylase" if hit else "elsewhere")] += 1
            counts[src]["mapq0"] += r.mapping_quality == 0
            fh.write(f"{r.query_name}\t{src}\t{where}\t{r.mapping_quality}\n")
    print("locus\tlevel2\treads\ton_source\tother_amylase\telsewhere\tunmapped\tmapq0")
    for name in sorted(rec):
        c = counts[name]
        tot = c["on_source"] + c["other_amylase"] + c["elsewhere"] + c["unmapped"]
        print(f"{name}\t{rec[name]['level2']}\t{tot}\t{c['on_source']}\t{c['other_amylase']}\t{c['elsewhere']}\t{c['unmapped']}\t{c['mapq0']}")


if __name__ == "__main__":
    main()

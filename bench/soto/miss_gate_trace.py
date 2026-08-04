#!/usr/bin/env python3
"""For each missed Soto member, replay the pipeline's gates in order and report WHICH one kills it.

The previous triage said "a chain of >=3 reads exists, so a locus is buildable" and concluded the loss must
be at the edge stage. Testing that with the summed-coverage rule recovered 2 of 32, so the inference was
wrong: >=3 reads sharing a chain is a NECESSARY condition for a locus, not a sufficient one. `assemble_gate`
and the transcript filters apply several more tests before an edge is ever considered.

This replays them in the pipeline's own order, per candidate chain, and reports the first gate each locus
fails:

  G1 reads          >= 3 reads share the exact intron chain      (GATE_MIN_READS)
  G2 span           skeleton span <= 3,000,000 bp                (MAX_SPAN)
  G3 canonical      EVERY junction canonical and strand-consistent (build_spliced_seq, strict)
  G4 spliced_len    100 <= exon-sum length <= 300,000            (MIN_SPLICED / MAX_SPLICED)
  G5 readthrough    not a single-exon model engulfing >= 5 junctions
  G6 mis-chain      no intron > 50 kb carried by < 3 reads

A locus survives only if at least one of its candidate chains passes all six. Reporting the FIRST failure
per chain, and whether ANY chain survives, separates "the gate is too strict" from "this never had a model".
"""
import argparse, re, subprocess, sys
from collections import defaultdict

PLUS = {(b"GT", b"AG"), (b"GC", b"AG"), (b"AT", b"AC")}
MINUS = {(b"CT", b"AC"), (b"CT", b"GC"), (b"GT", b"AT")}
GATE_MIN_READS, MAX_SPAN, MIN_SPLICED, MAX_SPLICED = 3, 3_000_000, 100, 300_000
MISCHAIN_BP, MISCHAIN_READS = 50_000, 3


def fetch(fasta, chrom, s, e):
    if e <= s:
        return b""
    r = subprocess.run(["samtools", "faidx", fasta, f"{chrom}:{s+1}-{e}"], capture_output=True).stdout
    return b"".join(r.split(b"\n")[1:]).upper()


def junction_strand(fasta, chrom, d, a):
    if a < d + 4:
        return None
    don, acc = fetch(fasta, chrom, d, d + 2), fetch(fasta, chrom, a - 2, a)
    if (don, acc) in PLUS:
        return "+"
    if (don, acc) in MINUS:
        return "-"
    return None


def chains(bam, chrom, s, e):
    txt = subprocess.run(["samtools", "view", "-F", "2308", bam, f"{chrom}:{s+1}-{e}"],
                         capture_output=True, text=True).stdout
    out = defaultdict(lambda: {"n": 0, "starts": [], "ends": []})
    junc = defaultdict(int)
    for ln in txt.splitlines():
        f = ln.split("\t")
        if len(f) < 6:
            continue
        pos, cig, ex, p = int(f[3]) - 1, f[5], [], int(f[3]) - 1
        for n, o in re.findall(r"(\d+)([MIDNSHP=X])", cig):
            n = int(n)
            if o in "M=X":
                if ex and ex[-1][1] == p:
                    ex[-1][1] = p + n
                else:
                    ex.append([p, p + n])
                p += n
            elif o == "D":
                if ex:
                    ex[-1][1] = p + n
                p += n
            elif o == "N":
                p += n
        if not ex:
            continue
        ch = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        r = out[ch]
        r["n"] += 1
        r["starts"].append(ex[0][0])
        r["ends"].append(ex[-1][1])
        for j in ch:
            junc[j] += 1
    return out, junc


def kth(vals, k, largest):
    v = sorted(vals, reverse=largest)
    return v[min(k - 1, len(v) - 1)]


def verdict(fasta, chrom, ch, rec, junc, all_junc):
    n = rec["n"]
    if n < GATE_MIN_READS:
        return "G1 reads"
    start, end = kth(rec["starts"], 2, False), kth(rec["ends"], 2, True)
    if end - start > MAX_SPAN:
        return "G2 span"
    strand = None
    for d, a in ch:
        st = junction_strand(fasta, chrom, d, a)
        if st is None:
            return "G3 canonical"
        if strand and st != strand:
            return "G3 strand-conflict"
        strand = st
    exonic = (end - start) - sum(a - d for d, a in ch if start <= d and a <= end)
    if exonic < MIN_SPLICED:
        return "G4 too short"
    if exonic > MAX_SPLICED:
        return "G4 too long"
    if not ch:
        engulfed = sum(1 for (d, a) in all_junc if start < d and a < end)
        if engulfed >= 5:
            return "G5 readthrough"
    for d, a in ch:
        if a - d > MISCHAIN_BP and junc.get((d, a), 0) < MISCHAIN_READS:
            return "G6 mis-chain"
    return "PASS"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--loci", nargs="+", required=True, help="name=chrom:start-end")
    a = ap.parse_args()
    for spec in a.loci:
        name, reg = spec.split("=", 1)
        chrom, rng = reg.split(":")
        s, e = (int(x.replace(",", "")) for x in rng.split("-"))
        ch, junc = chains(a.bam, chrom, s, e)
        all_junc = set(junc)
        ranked = sorted(ch.items(), key=lambda kv: -kv[1]["n"])
        counts = defaultdict(int)
        survivors = []
        for c, rec in ranked:
            v = verdict(a.fasta, chrom, c, rec, junc, all_junc)
            counts[v] += 1
            if v == "PASS":
                survivors.append((c, rec["n"]))
        tot = sum(r["n"] for r in ch.values())
        print(f"\n{name}  {chrom}:{s:,}-{e:,}   {tot} reads, {len(ch)} chains")
        for k in sorted(counts, key=lambda x: -counts[x]):
            print(f"    {counts[k]:>4} chains -> {k}")
        if survivors:
            best = max(survivors, key=lambda x: x[1])
            print(f"    => {len(survivors)} chain(s) PASS every gate; best carries {best[1]} reads"
                  f" ({len(best[0])} introns)")
        else:
            print("    => NO chain passes: this locus never produces a transcript model")


if __name__ == "__main__":
    main()

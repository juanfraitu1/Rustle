#!/usr/bin/env python3
"""Final-output parity gate: compare intron-chain sets of two GTFs (Rustle vs StringTie).

This is the HONEST shadow-mode gate (supersedes the pre-collapse transfrag_pre_depl proxy):
it compares the actual emitted transcripts. Reports multi-intron chain overlap and, optionally,
single-exon transcript counts. Drive Rustle-only -> 0 (and understand ST-only) for parity.

Usage: python3 bench/gtf_chain_diff.py <rustle.gtf> <stringtie.gtf>
"""
import re, sys, collections

def load(path):
    tx = collections.defaultdict(list)
    strand = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        tx[m.group(1)].append((int(f[3]), int(f[4])))
        strand[m.group(1)] = f[6]
    multi = set()   # (strand, intron-chain tuple)
    single = 0
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2:
            single += 1
            continue
        introns = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        multi.add((strand[t], introns))
    return multi, single, len(tx)

ru_path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/ru_final.gtf"
st_path = sys.argv[2] if len(sys.argv) > 2 else "/tmp/st_final.gtf"
ru, ru_se, ru_n = load(ru_path)
st, st_se, st_n = load(st_path)
print(f"Rustle: {ru_n} tx ({len(ru)} multi-intron chains, {ru_se} single-exon)")
print(f"ST:     {st_n} tx ({len(st)} multi-intron chains, {st_se} single-exon)")
print(f"  multi-intron in both: {len(ru & st)}  Rustle-only: {len(ru - st)}  ST-only: {len(st - ru)}")
print(f"\nFINAL-PARITY GATE (Rustle-only multi-intron chains): {len(ru - st)} (target 0)")

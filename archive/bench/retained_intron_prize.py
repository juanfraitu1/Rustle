#!/usr/bin/env python3
"""Phase 0a: bound the retained-intron prize. Match ST pred_kill(retained_intron) chains
to RU final chains; net = (RU-only-FP matched) - (RU-TP matched).
Usage: python3 bench/retained_intron_prize.py /tmp/ri_st.jsonl /tmp/ri_ru.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf"""
import json, sys, re
from collections import defaultdict

def gtf_chains(path):
    tx = defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    chains = set()
    for t, ex in tx.items():
        ex.sort()
        ic = tuple((ex[i-1][1] + 1, ex[i][0] - 1) for i in range(1, len(ex)))
        if ic: chains.add((strand[t], ic))
    return chains

def st_killed_chains(path):
    out = set()
    for line in open(path):
        if '"step":"pred_kill"' not in line: continue
        try: r = json.loads(line)
        except ValueError: continue
        p = r.get("payload", {})
        if p.get("reason") != "retained_intron": continue
        ch = p.get("chain", "")
        if not ch: continue
        ic = tuple(tuple(map(int, tok.split("-"))) for tok in ch.split(","))
        out.add((r["strand"], ic))
    return out

def main():
    st_killed = st_killed_chains(sys.argv[1] if len(sys.argv) > 1 else "/tmp/ri_st.jsonl")
    ru = gtf_chains(sys.argv[2] if len(sys.argv) > 2 else "/tmp/ri_ru.gtf")
    ref = gtf_chains(sys.argv[3] if len(sys.argv) > 3 else "/mnt/c/Users/jfris/Desktop/GGO_19.gtf")
    print(f"ST retained_intron-killed chains: {len(st_killed)}")
    print(f"RU final chains: {len(ru)}  (RU-only-FP vs ref: {len(ru - ref)})")
    matched = ru & st_killed
    fp_removed = matched - ref
    tp_lost = matched & ref
    print(f"RU chains matched by an ST retained_intron kill: {len(matched)}")
    print(f"  of which FP (prize): {len(fp_removed)}")
    print(f"  of which TP (cost):  {len(tp_lost)}")
    net = len(fp_removed) - len(tp_lost)
    print(f"NET prize (FP_removed - TP_lost): {net}")
    print(f"ABORT GATE: {'ABORT (shelve Phase 1)' if net <= 0 else ('WEAK (<5)' if net < 5 else 'PROCEED')}")
    for c in list(tp_lost)[:8]:
        print(f"    TP-COST chain {c[0]} {c[1][:3]}...")

if __name__ == "__main__":
    main()

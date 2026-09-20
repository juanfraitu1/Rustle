#!/usr/bin/env python3
"""Classify chimeric_suffix_rescue (CSR) predictions: FP vs TP vs reference.
CSR tx are identified by GTF source="checktrf_rescue" or rescue_class="chimeric_suffix_rescue".
Usage: python3 bench/csr_classify.py /tmp/csr_base.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf"""
import sys, re
from collections import defaultdict

def load(path):
    tx_ex = defaultdict(list); strand = {}; is_csr = {}
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tid = m.group(1)
        if f[2] in ("transcript", "exon"):
            csr = ('checktrf_rescue' in f[1]) or ('chimeric_suffix_rescue' in f[8])
            is_csr[tid] = is_csr.get(tid, False) or csr
        if f[2] == "exon":
            tx_ex[tid].append((int(f[3]), int(f[4]))); strand[tid] = f[6]
    chains = {}
    for t, ex in tx_ex.items():
        ex.sort()
        ic = tuple((ex[i-1][1]+1, ex[i][0]-1) for i in range(1, len(ex)))
        if ic: chains[t] = (strand[t], ic, is_csr.get(t, False))
    return chains

def main():
    ru = load(sys.argv[1] if len(sys.argv) > 1 else "/tmp/csr_base.gtf")
    ref = load(sys.argv[2] if len(sys.argv) > 2 else "/mnt/c/Users/jfris/Desktop/GGO_19.gtf")
    refset = {(s, ic) for (s, ic, _) in ref.values()}
    csr = [(t, s, ic) for t, (s, ic, c) in ru.items() if c]
    fp = [(t, s, ic) for (t, s, ic) in csr if (s, ic) not in refset]
    tp = [(t, s, ic) for (t, s, ic) in csr if (s, ic) in refset]
    print(f"CSR predictions (multi-intron): {len(csr)}")
    print(f"  FP (not in ref) - suppressible prize: {len(fp)}")
    print(f"  TP (in ref) - cost if suppressed:     {len(tp)}")
    net = len(fp) - len(tp)
    print(f"NET (CSR-FP - CSR-TP): {net}")
    print(f"ABORT GATE: {'ABORT (shelve 2B)' if net <= 0 else ('WEAK (<3)' if net < 3 else 'PROCEED')}")
    for t, s, ic in tp[:8]:
        print(f"    CSR-TP {t} {s} {ic[:2]}...")

if __name__ == "__main__":
    main()

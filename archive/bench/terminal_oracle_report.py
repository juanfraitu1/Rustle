#!/usr/bin/env python3
"""Compare oracle-ON vs baseline final GTFs against the reference; report FP/FN/F1 delta
and the causal NET (FP_removed - TP_lost - FP_added).
Usage: python3 bench/terminal_oracle_report.py /tmp/ru_def.gtf /tmp/ru_oracle.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf"""
import sys, re
from collections import defaultdict

def chains(path):
    tx = defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip().split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    out = set()
    for t, ex in tx.items():
        ex.sort()
        ic = tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))
        if ic: out.add((strand[t], ic))
    return out

def main():
    base = chains(sys.argv[1]); orac = chains(sys.argv[2]); ref = chains(sys.argv[3])
    def stats(s):
        tp = len(s & ref); fp = len(s - ref); fn = len(ref - s)
        sn = 100*tp/len(ref) if ref else 0; pr = 100*tp/(tp+fp) if (tp+fp) else 0
        f1 = 2*sn*pr/(sn+pr) if (sn+pr) else 0
        return tp, fp, fn, sn, pr, f1
    bt, bfp, bfn, bsn, bpr, bf1 = stats(base)
    ot, ofp, ofn, osn, opr, of1 = stats(orac)
    print(f"baseline: TP={bt} FP={bfp} FN={bfn}  Sn={bsn:.1f} Pr={bpr:.1f} F1={bf1:.2f}")
    print(f"oracle:   TP={ot} FP={ofp} FN={ofn}  Sn={osn:.1f} Pr={opr:.1f} F1={of1:.2f}")
    print(f"DELTA:    dTP={ot-bt:+d} dFP={ofp-bfp:+d} dFN={ofn-bfn:+d}  dF1={of1-bf1:+.2f}")
    fp_removed = (base - ref) - (orac - ref)
    fp_added = (orac - ref) - (base - ref)
    tp_lost = (base & ref) - (orac & ref)
    print(f"FP removed by oracle: {len(fp_removed)}  FP added: {len(fp_added)}  TP lost: {len(tp_lost)}")
    net = len(fp_removed) - len(tp_lost) - len(fp_added)
    print(f"NET prize (FP_removed - TP_lost - FP_added): {net}")
    print(f"ABORT GATE: {'ABORT (shelve Phase 1)' if net <= 0 else ('WEAK (<5)' if net < 5 else 'PROCEED')}")

if __name__ == "__main__":
    main()

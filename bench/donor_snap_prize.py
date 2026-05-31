#!/usr/bin/env python3
"""2A prize: RU final chains using a weak donor within 25bp of a stronger same-strand donor.
prize = such RU-only-FP chains whose weak intron is NOT a reference intron (snap would remove them);
cost  = chains whose weak intron IS a reference intron (real alt-donor the snap would destroy).
Usage: python3 bench/donor_snap_prize.py /tmp/csr_base.gtf /mnt/c/Users/jfris/Desktop/GGO_19.gtf"""
import sys, re
from collections import defaultdict

def load(path):
    tx_ex = defaultdict(list); strand = {}; cov = {}
    for line in open(path):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tid = m.group(1)
        if f[2] == "transcript":
            cm = re.search(r'cov "([0-9.]+)"', f[8]); cov[tid] = float(cm.group(1)) if cm else 0.0
        if f[2] == "exon":
            tx_ex[tid].append((int(f[3]), int(f[4]))); strand[tid] = f[6]
    out = {}
    for t, ex in tx_ex.items():
        ex.sort()
        introns = [(ex[i-1][1]+1, ex[i][0]-1) for i in range(1, len(ex))]
        out[t] = (strand[t], introns, cov.get(t, 0.0))
    return out

def main():
    ru = load(sys.argv[1] if len(sys.argv) > 1 else "/tmp/csr_base.gtf")
    ref = load(sys.argv[2] if len(sys.argv) > 2 else "/mnt/c/Users/jfris/Desktop/GGO_19.gtf")
    ref_introns = set()
    for (s, introns, _c) in ref.values():
        for (d, a) in introns:
            ref_introns.add((s, d, a))
    donor_cov = defaultdict(float)
    intron_cov = defaultdict(float)
    for (s, introns, c) in ru.values():
        for (d, a) in introns:
            donor_cov[(s, d)] = max(donor_cov[(s, d)], c)
            intron_cov[(s, d, a)] = max(intron_cov[(s, d, a)], c)
    WIN = 25
    prize_chains = set(); cost_chains = set()
    for t, (s, introns, c) in ru.items():
        for (d, a) in introns:
            stronger = any(donor_cov[(s, d2)] > intron_cov[(s, d, a)] * 1.1
                           for d2 in range(d - WIN, d + WIN + 1) if (s, d2) in donor_cov and d2 != d)
            if not stronger:
                continue
            if (s, d, a) in ref_introns:
                cost_chains.add(t)
            else:
                prize_chains.add(t)
    refset = set()
    for (s, ii, c) in ref.values():
        if ii: refset.add((s, tuple(ii)))
    ru_fp = {t for t, (s, ii, c) in ru.items() if ii and (s, tuple(ii)) not in refset}
    prize = prize_chains & ru_fp
    cost = cost_chains
    print(f"chains w/ weak donor near stronger same-strand donor (<=25bp): prize-side {len(prize_chains)}, cost-side {len(cost_chains)}")
    print(f"PRIZE (RU-only-FP w/ weak near-dup donor): {len(prize)}")
    print(f"COST  (chains whose weak intron IS a real reference intron): {len(cost)}")
    net = len(prize) - len(cost)
    print(f"NET (prize - cost): {net}")
    print(f"ABORT GATE: {'ABORT (cost>=prize)' if len(cost) >= len(prize) else ('WEAK (<3)' if net < 3 else 'PROCEED')}")
    for t in list(cost)[:8]:
        print(f"    COST chain {t} strand={ru[t][0]}")

if __name__ == "__main__":
    main()

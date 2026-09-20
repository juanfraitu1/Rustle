#!/usr/bin/env python3
"""For each ST-only multi-intron chain, find rustle's best-matching chain (max shared
introns) and report the exact difference. Reveals near-misses (1-2 intron edits) vs
clean misses. Uses FINAL GTFs.

Usage: python3 nearmiss.py /tmp/m3_ru2.gtf /tmp/m3_st.gtf [lo hi]
"""
import re, sys, collections

def load(path):
    tx = collections.defaultdict(list); strand = {}
    for line in open(path):
        if line.startswith('#'): continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon': continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        tx[m.group(1)].append((int(f[3]), int(f[4]))); strand[m.group(1)] = f[6]
    chains = {}
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2: continue
        introns = tuple((ex[i][1], ex[i+1][0]) for i in range(len(ex)-1))
        chains[(strand[t], introns)] = (t, ex[0][0], ex[-1][1])
    return chains

ru = load(sys.argv[1]); st = load(sys.argv[2])
lo = int(sys.argv[3]) if len(sys.argv) > 3 else 0
hi = int(sys.argv[4]) if len(sys.argv) > 4 else 10**12

st_only = [k for k in st.keys() - ru.keys() if lo <= st[k][1] <= hi]
ru_keys = list(ru.keys())

def best_match(target_introns, strand):
    tset = set(target_introns); best = None; best_shared = -1
    for (s, introns) in ru_keys:
        if s != strand: continue
        shared = len(tset & set(introns))
        if shared > best_shared:
            best_shared = shared; best = (s, introns)
    return best, best_shared

for k in sorted(st_only, key=lambda k: st[k][1]):
    strand, introns = k
    tid, s, e = st[k]
    iset = set(introns)
    bm, shared = best_match(introns, strand)
    print(f"\n=== ST-only {tid} {s}-{e} {strand} nintrons={len(introns)} ===")
    if bm is None:
        print("  NO rustle chain on this strand"); continue
    bm_introns = bm[1]; bset = set(bm_introns)
    rtid, rs, re_ = ru[bm]
    only_st = sorted(iset - bset); only_ru = sorted(bset - iset)
    print(f"  best rustle match: {rtid} {rs}-{re_} nintrons={len(bm_introns)} shared={shared}")
    print(f"  introns only in ST ({len(only_st)}): {only_st}")
    print(f"  introns only in RU ({len(only_ru)}): {only_ru}")

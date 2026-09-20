#!/usr/bin/env python3
"""List the exact divergent multi-intron chains (rustle-vs-StringTie) on mini3,
each tagged with its locus (A/B/C) and span, to drive first-divergence tracing.

Usage: python3 bench/mini3/list_divergent.py /tmp/m3_ru.gtf /tmp/m3_st.gtf
"""
import re, sys, collections

LOCI = [("A", 17160875, 17175957), ("B", 36012544, 36069321), ("C", 68848317, 68871262)]

def locus_of(start, end):
    mid = (start + end) // 2
    for name, lo, hi in LOCI:
        if lo - 5000 <= mid <= hi + 5000:
            return name
    return "?"

def load(path):
    tx = collections.defaultdict(list)
    strand = {}
    tid_of = {}
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
    chains = {}   # (strand, introns) -> (tid, start, end)
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2:
            continue
        introns = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        key = (strand[t], introns)
        chains[key] = (t, ex[0][0], ex[-1][1])
    return chains

ru = load(sys.argv[1] if len(sys.argv) > 1 else "/tmp/m3_ru.gtf")
st = load(sys.argv[2] if len(sys.argv) > 2 else "/tmp/m3_st.gtf")

def dump(label, keys, src):
    print(f"\n### {label} ({len(keys)})")
    rows = []
    for k in keys:
        strand, introns = k
        tid, s, e = src[k]
        rows.append((locus_of(s, e), s, e, strand, len(introns), introns, tid))
    rows.sort()
    for loc, s, e, strand, nin, introns, tid in rows:
        ich = ";".join(f"{a}-{b}" for a, b in introns)
        print(f"  [{loc}] {s}-{e} {strand} nintrons={nin} tid={tid}")
        print(f"       introns: {ich}")

dump("ST-only (rustle MISSES these)", st.keys() - ru.keys(), st)
dump("Rustle-only (rustle EXTRA)", ru.keys() - st.keys(), ru)
print(f"\nSUMMARY: in-both={len(ru.keys() & st.keys())} rustle-only={len(ru.keys()-st.keys())} st-only={len(st.keys()-ru.keys())}")

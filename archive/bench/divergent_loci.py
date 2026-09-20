#!/usr/bin/env python3
"""Cluster rustle-vs-ST divergent multi-intron chains into loci (for the coupled-port work-list).

Usage: python3 bench/divergent_loci.py <rustle.gtf> <stringtie.gtf> [--json]
Emits, per locus (overlapping span cluster): rustle-only and ST-only chains with coords.
"""
import re, sys, collections, json

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
    chains = {}   # (strand, introns) -> (tid, start, end)
    for t, ex in tx.items():
        ex.sort()
        if len(ex) < 2:
            continue
        introns = tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1))
        chains[(strand[t], introns)] = (t, ex[0][0], ex[-1][1])
    return chains

ru = load(sys.argv[1])
st = load(sys.argv[2])
ru_only = set(ru) - set(st)
st_only = set(st) - set(ru)

# Build divergent-chain records with span.
recs = []
for key in ru_only:
    strand, introns = key
    recs.append({'who':'RU','strand':strand,'introns':introns,
                 'start':min(i[0] for i in introns),'end':max(i[1] for i in introns),
                 'nintr':len(introns),'tid':ru[key][0]})
for key in st_only:
    strand, introns = key
    recs.append({'who':'ST','strand':strand,'introns':introns,
                 'start':min(i[0] for i in introns),'end':max(i[1] for i in introns),
                 'nintr':len(introns),'tid':st[key][0]})

# Cluster by overlapping span (greedy interval merge), per strand.
recs.sort(key=lambda r:(r['strand'], r['start'], r['end']))
loci = []
cur = None
for r in recs:
    if cur and r['strand']==cur['strand'] and r['start'] <= cur['end']:
        cur['members'].append(r); cur['end']=max(cur['end'], r['end'])
    else:
        cur = {'strand':r['strand'],'start':r['start'],'end':r['end'],'members':[r]}
        loci.append(cur)

# Summaries
n_ru = sum(1 for r in recs if r['who']=='RU')
n_st = sum(1 for r in recs if r['who']=='ST')
print(f"divergent chains: RU-only={n_ru} ST-only={n_st} total={n_ru+n_st} in {len(loci)} loci")
mixed = [L for L in loci if any(m['who']=='RU' for m in L['members']) and any(m['who']=='ST' for m in L['members'])]
ru_pure = [L for L in loci if all(m['who']=='RU' for m in L['members'])]
st_pure = [L for L in loci if all(m['who']=='ST' for m in L['members'])]
print(f"  MIXED loci (both RU+ST divergent — swap candidates): {len(mixed)}")
print(f"  RU-pure loci (rustle invents, no ST-only here): {len(ru_pure)}  ({sum(len(L['members']) for L in ru_pure)} chains)")
print(f"  ST-pure loci (rustle misses, no RU-only here):  {len(st_pure)}  ({sum(len(L['members']) for L in st_pure)} chains)")

if '--json' in sys.argv:
    out = []
    for L in loci:
        out.append({'strand':L['strand'],'start':L['start'],'end':L['end'],
                    'n_ru':sum(1 for m in L['members'] if m['who']=='RU'),
                    'n_st':sum(1 for m in L['members'] if m['who']=='ST'),
                    'members':[{'who':m['who'],'nintr':m['nintr'],'start':m['start'],'end':m['end'],'tid':m['tid'],
                               'introns':';'.join(f"{a}-{b}" for a,b in m['introns'])} for m in L['members']]})
    json.dump(out, open('/tmp/divergent_loci.json','w'), indent=0)
    print("wrote /tmp/divergent_loci.json")
else:
    # Print mixed loci first (the coupled swap candidates), then sizes
    print("\n=== MIXED loci (sorted by size) ===")
    for L in sorted(mixed, key=lambda L:-(len(L['members']))):
        nr=sum(1 for m in L['members'] if m['who']=='RU'); ns=sum(1 for m in L['members'] if m['who']=='ST')
        print(f"  [{L['strand']}] {L['start']}-{L['end']}  RU={nr} ST={ns}")

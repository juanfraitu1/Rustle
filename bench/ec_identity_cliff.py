#!/usr/bin/env python3
"""E_c support per E_r edge, with the placement-overlap rule swept.

A placement counts as "at" a locus when >= OV of the placement's own reference
span lies inside the locus interval. OV=0 is 'any overlap'. Sweeping it tests
whether the identity cliff is an artefact of the assignment rule.
Tie set = placements scoring >= 0.95 * the read's best AS ANYWHERE in the genome.
"""
import sys, re, bisect, pickle
from collections import defaultdict

DUMP = "/mnt/linuxdisk/home/juanfraitu/npip_cat/arm_f2/dump"
OUT  = "/mnt/linuxdisk/home/juanfraitu/ecband"
TIE  = 0.95

nodes, by_chrom = {}, defaultdict(list)
with open(f"{DUMP}/e.nodes.tsv") as fh:
    h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        k, c = f[ci["node_key"]], f[ci["chrom"]]
        s, e = int(f[ci["start"]]), int(f[ci["end"]])
        nodes[k] = (c, s, e); by_chrom[c].append((s, e, k))
for c in by_chrom: by_chrom[c].sort()
starts = {c: [x[0] for x in v] for c, v in by_chrom.items()}
maxend = {}
for c, v in by_chrom.items():
    acc, m = 0, []
    for s, e, k in v: acc = max(acc, e); m.append(acc)
    maxend[c] = m
keyidx = {k: i for i, k in enumerate(nodes)}
keylist = list(nodes)

def hits(chrom, s, e):
    v = by_chrom.get(chrom)
    if not v: return ()
    out, j = [], bisect.bisect_right(starts[chrom], e) - 1
    span = max(e - s, 1)
    while j >= 0:
        if maxend[chrom][j] <= s: break
        ns, ne, k = v[j]
        ov = min(e, ne) - max(s, ns)
        if ov > 0: out.append((keyidx[k], ov / span))
        j -= 1
    return out

CIG = re.compile(r"(\d+)([MIDNSHP=X])"); RC = set("MDN=X")
best, plac = {}, defaultdict(list)
n = 0
for line in sys.stdin:
    f = line.split("\t", 11)
    flag = int(f[1])
    if flag & 0x4 or flag & 0x800: continue
    tags = f[11] if len(f) > 11 else ""
    p = tags.find("AS:i:")
    if p < 0: continue
    q2 = tags.find("\t", p)
    asv = int(tags[p+5:] if q2 < 0 else tags[p+5:q2])
    q = f[0]
    b = best.get(q)
    if b is None or asv > b: best[q] = asv
    chrom = f[2]
    if chrom in by_chrom:
        pos = int(f[3]) - 1; ref = 0
        for num, op in CIG.findall(f[5]):
            if op in RC: ref += int(num)
        hh = hits(chrom, pos, pos + ref)
        if hh: plac[q].append((asv, hh))
    n += 1
    if n % 500000 == 0: sys.stderr.write(f"  {n//1000}k records\n")
sys.stderr.write(f"records={n} reads={len(best)} reads_on_nodes={len(plac)}\n")

pairs_by_ov = {}
for OV in (0.0, 0.25, 0.50, 0.80):
    pair = defaultdict(int); neg = 0
    for q, pl in plac.items():
        bq = best[q]
        if bq <= 0: neg += 1; continue
        thr = TIE * bq
        ks = set()
        for a, hh in pl:
            if a < thr: continue
            for ki, frac in hh:
                if frac >= OV: ks.add(ki)
        if len(ks) < 2: continue
        ks = sorted(ks)
        for i in range(len(ks)):
            for j in range(i+1, len(ks)): pair[(ks[i], ks[j])] += 1
    pairs_by_ov[OV] = pair
    sys.stderr.write(f"OV={OV}: ec_pairs={len(pair)} (skipped_nonpos={neg})\n")

rows, er = [], set()
with open(f"{DUMP}/e.edges.tsv") as fh:
    h = fh.readline().rstrip("\n").split("\t"); ci = {c: i for i, c in enumerate(h)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        a, b = keyidx[f[ci["node_key_i"]]], keyidx[f[ci["node_key_j"]]]
        key = (a, b) if a < b else (b, a)
        er.add(key)
        rows.append((key, float(f[ci["identity"]]), float(f[ci["coverage"]])))

with open(f"{OUT}/edge_ec2.tsv", "w") as o:
    o.write("node_i\tnode_j\tidentity\tcoverage\t" +
            "\t".join(f"reads_ov{int(v*100):02d}" for v in pairs_by_ov) + "\n")
    for key, idn, cov in rows:
        cts = "\t".join(str(pairs_by_ov[v].get(key, 0)) for v in pairs_by_ov)
        o.write(f"{keylist[key[0]]}\t{keylist[key[1]]}\t{idn:.6f}\t{cov:.6f}\t{cts}\n")

tot = len(rows)
for OV, pair in pairs_by_ov.items():
    k = sum(1 for key, _, _ in rows if pair.get(key, 0) >= 1)
    extra = sum(1 for kk in pair if kk not in er)
    sys.stderr.write(f"OV={OV}: E_r supported {k}/{tot} = {k/tot:.4f} | E_c-only pairs {extra}\n")

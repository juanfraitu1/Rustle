#!/usr/bin/env python3
"""Why are observed (>= 2 molecules) reference transcripts missed? (descriptive; observed.py follow-up)

usage: missed.py REF.gtf STRATUM.tsv COUNTS.tsv TMAP REFMAP OUT.gtf [NOPOLISH.gtf] > missed.tsv
Classes, first that applies: shadowed (a query is `=` to it in the refmap but the tmap gave that query to another
reference with the same chain); chain_out (OUT.gtf has a transcript with its exact intron chain on its contig, yet no
`=`: strand or class mismatch); polish (the chain is in NOPOLISH.gtf only); not_built (in neither).
"""
import re
import sys
from collections import defaultdict

RT = re.compile(r'transcript_id "([^"]+)"')


def chains(gtf):
    ex, key = defaultdict(list), {}
    for l in open(gtf):
        f = l.split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        m = RT.search(f[8])
        if m:
            ex[m.group(1)].append((int(f[3]), int(f[4])))
            key[m.group(1)] = (f[0], f[6])
    out = {}
    for t, e in ex.items():
        if len(e) >= 2:
            e.sort()
            out[t] = (key[t][0], key[t][1], ",".join(f"{e[i][1]}-{e[i + 1][0] - 1}" for i in range(len(e) - 1)), e[0][0], e[-1][1])
    return out


ref, strat, counts, tmap, refmap, outg = sys.argv[1:7]
nop = sys.argv[7] if len(sys.argv) > 7 else None
rc = chains(ref)
stratum = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2] for l in open(strat)}
want = defaultdict(list)
for t, (c, s, ch, a, b) in rc.items():
    want[(c, ch)].append(t)
obs, q0 = {}, {}
for l in open(counts):
    c, ch, n, z = l.rstrip("\n").split("\t")
    for t in want.get((c, ch), ()):
        obs[t], q0[t] = int(n), int(z)
hit = set()
for l in open(tmap):
    f = l.split("\t")
    if f[2] == "=":
        hit.add(f[1])
ref_eq = set()
for l in open(refmap):
    f = l.rstrip("\n").split("\t")
    if len(f) > 3 and f[0] != "ref_gene" and f[2] == "=":
        ref_eq.add(f[1])
oc = defaultdict(list)
for t, (c, s, ch, a, b) in chains(outg).items():
    oc[(c, ch)].append((t, s))
nc = set()
if nop:
    nc = {(c, ch) for c, s, ch, a, b in chains(nop).values()}
print("ref_tx\tstratum\tmolecules\tmapq0\tn_exons\tsame_chain_refs\tclass\tdetail")
for t, (c, s, ch, a, b) in sorted(rc.items()):
    if obs.get(t, 0) < 2 or t in hit:
        continue
    same = len(want[(c, ch)])
    if t in ref_eq:
        cl, d = "shadowed", ""
    elif (c, ch) in oc:
        cl, d = "chain_out", ";".join(f"{q}:{st}" for q, st in oc[(c, ch)][:3]) + f" ref_strand {s}"
    elif (c, ch) in nc:
        cl, d = "polish", ""
    else:
        cl, d = "not_built" if nop else "absent", ""
    print(f"{t}\t{stratum.get(t, '?')}\t{obs[t]}\t{q0[t]}\t{ch.count(',') + 2}\t{same}\t{cl}\t{d}")

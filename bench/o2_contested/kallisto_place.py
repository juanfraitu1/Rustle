#!/usr/bin/env python3
"""r1647 (docs/PREREG_kallisto_long_placement_2026-11-12.md): picks from `kallisto bus --long -n` output. A read is
assigned to a copy iff every sequence of its compatibility class belongs to copies of ONE locus and none is OTHER.
Usage: kallisto_place.py BUS_DIR READS.fa COPIES.tsv > OUT.picks.tsv   (counts on stderr)"""
import collections, csv, struct, sys

def main(bus_dir, reads, copies):
    span = {}
    for r in csv.DictReader(open(copies), delimiter="\t"): span[(r["family_id"], r["copy_idx"])] = (r["chrom"], int(r["start"]), int(r["end"]))
    names = [l[1:].split()[0] for l in open(reads) if l[0] == ">"]
    tx = [tuple(l.split("|")[:2]) for l in open(f"{bus_dir}/transcripts.txt")]
    ec = {}
    for l in open(f"{bus_dir}/matrix.ec"):
        i, ts = l.split("\t"); ec[int(i)] = {tx[int(t)] for t in ts.split(",")}
    b = open(f"{bus_dir}/output.bus", "rb").read()
    assert b[:4] == b"BUS\0"; off = 20 + struct.unpack("<I", b[16:20])[0]
    comp = collections.defaultdict(set)  # read number (flags, from -n) -> copies
    for o in range(off, len(b), 32):
        _, _, e, _, flag, _ = struct.unpack("<QQiIII", b[o:o + 32]); comp[flag] |= ec[e]
    n = collections.Counter(reads=len(names), pseudoaligned=len(comp))
    def same(a, c): return span[a][0] == span[c][0] and span[a][1] < span[c][2] and span[c][1] < span[a][2]
    for i, ks in sorted(comp.items()):
        ks = sorted(ks)
        if any(k[0] == "OTHER" for k in ks): n["abstain_other"] += 1; continue
        if not all(same(ks[0], k) for k in ks[1:]): n["abstain_multi_locus"] += 1; continue
        n["assigned"] += 1; q = names[i]
        print(q[2:] if q.startswith("R|") else q, ks[0][0], ks[0][1], sep="\t")
    print(dict(n), f"pseudoaligned fraction {n['pseudoaligned'] / max(n['reads'], 1):.3f}", file=sys.stderr)

if __name__ == "__main__": main(*sys.argv[1:4])

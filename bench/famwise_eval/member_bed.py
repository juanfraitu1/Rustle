#!/usr/bin/env python3
"""BED of the truth-family members (famtruth recent paralogue components, >= 2 genes) padded by PAD bp, merged: the slice
regions for member_completeness.py on a whole BAM.   member_bed.py OUT.bed [--species human] [--pad 5000]"""
import argparse, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "denovo_eval"))
import famtruth as F
from member_completeness import families
ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("--species", default="human"); ap.add_argument("--pad", type=int, default=5000)
a = ap.parse_args()
G = F.load_genes(a.species); fams, fam_of = families(G, a.species)
iv = sorted((G[g]["chrom"], max(0, min(s for s, _ in G[g]["exons"]) - a.pad), max(e for _, e in G[g]["exons"]) + a.pad) for g in fam_of)
out = []
for c, s, e in iv:
    if out and out[-1][0] == c and s <= out[-1][2]:
        out[-1][2] = max(out[-1][2], e)
    else:
        out.append([c, s, e])
with open(a.out, "w") as fh:
    fh.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in out)
print(f"{len(fam_of)} members, {len(out)} regions, {sum(e - s for _, s, e in out) / 1e6:.1f} Mb")

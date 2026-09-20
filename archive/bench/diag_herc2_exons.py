#!/usr/bin/env python3
import os, sys
if os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)

import pysam
BENCH = os.path.dirname(os.path.abspath(__file__))
SCRATCH = "/home/juanfra/winloci_scratch"
sys.path.insert(0, BENCH)
import recombination_bridge_detector as R

fa = pysam.FastaFile(os.path.join(SCRATCH, "GGO.fasta"))
skel = R.load_skeletons()
strand = R.load_strand()

# HERC2 loci from fam384 and fam385
loci = [
    dict(dn="DN_NC_073240.2_29369420_15", gene="HERC2", chrom="NC_073240.2", start=29369420, end=29393816),
    dict(dn="DN_NC_073240.2_29509629_20", gene="HERC2", chrom="NC_073240.2", start=29509629, end=29577246),
]
recs = R.family_exons(loci, skel, strand, fa)
for r in recs:
    print(r["dn"], r["gene"], r["strand"], "nexons", len(r["exons"]), "exlens", r["exlen"])

best = R.exon_match_tensor(recs)
print("\ncross matches (id, target_j):")
for i in range(len(recs[0]["exons"])):
    idv, j = best[(0, i, 1)]
    if idv >= 0.5:
        print(f"  A exon {i} -> B exon {j} id={idv:.3f}")
print("\nstrict colinear count:", end=" ")
import colinear_multiexon_gate as CM
print(CM.colinear_count(best, 0, 1, len(recs[0]["exons"]), strict=True, thresh=R.ID_THRESH))

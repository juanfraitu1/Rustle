#!/usr/bin/env python3
"""B2-B4 of docs/PREREG_psv_significant_blocks_2026-10-12.md: which labelled (copy, column) pairs survive in a new dump.

    python3 block_eval.py LABELS.tsv NEW_DUMP_PREFIX
A pair (family, copy-0 genome position, copy, allele) is retained iff the new dump has a column at that copy-0 position
in that family where the copy carries the same allele.
"""
import collections, json, os, sys
sys.path.insert(0, os.path.dirname(__file__))
import psv_sun as P

labels, prefix = sys.argv[1], sys.argv[2]
pos = {}
for r in P.tsv(prefix + ".psv_cols.tsv"):
    pos[(r["family_id"], int(r["col_index"]))] = int(r["genome_pos"])
alle = {}
for r in P.tsv(prefix + ".psv_copies.tsv"):
    alle[(r["family_id"], int(r["copy_index"]))] = r["alleles"]
have = set()
for (fam, j), g in pos.items():
    for (f2, i), a in alle.items():
        if f2 == fam and j < len(a) and a[j] != ".":
            have.add((fam, g, i, a[j]))
tot, kept = collections.Counter(), collections.Counter()
for r in P.tsv(labels):
    k = (r["family_id"], int(r["gpos0"]), int(r["copy"]), r["allele"])
    tot[r["cls"]] += 1
    kept[r["cls"]] += k in have
res = {c: {"pairs": tot[c], "retained": kept[c], "retained_frac": round(kept[c] / tot[c], 4)} for c in tot}
nonhom = ["both_unrelated", "terminal_overhang", "no_alignment_unrelated"]
n, k = sum(tot[c] for c in nonhom), sum(kept[c] for c in nonhom)
res["B2_aligned_retained"] = res["aligned"]["retained_frac"]
res["B2_pass"] = res["aligned"]["retained_frac"] >= 0.99
res["B3_nonhomologous_removed"] = round(1 - k / n, 4)
res["B3_pass"] = 1 - k / n >= 0.90
print(json.dumps(res, indent=1))

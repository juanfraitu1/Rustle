#!/usr/bin/env python3
"""Prereg Addendum L1: hybrid (expressed SD atoms) plus RNA-only families.

usage: hybrid_union.py <dna_copies.tsv> <rna_copies.tsv> <expr.tsv> <out_prefix>
<expr.tsv> is `bench/interval_expression.py` output over the atoms (keys `atom:<idx>`). Writes <out>.H.copies.tsv
(expressed atoms, families with >= 2) and <out>.U1/U2/U3.copies.tsv:
  U1 = H + RNA families restricted to copies overlapping NO atom, >= 2 such copies;
  U2 = H + RNA families restricted to copies overlapping no EXPRESSED atom, >= 2 such copies;
  U3 = H + every RNA family, whole, with >= 1 copy overlapping no expressed atom.
RNA family ids are prefixed `RNA:` so they never collide with atom family ids.
"""
import bisect
import collections
import csv
import sys

GATE = 3
dna_path, rna_path, expr_path, out = sys.argv[1:5]
atoms_all, atoms_expr = collections.defaultdict(list), collections.defaultdict(list)
expr = {}
for r in csv.DictReader(open(expr_path), delimiter="\t"):
    if not r["key"].startswith("atom:"):
        continue
    iv = (int(r["start"]), int(r["end"]))
    expr[(r["chrom"],) + iv] = int(r["u"])
    atoms_all[r["chrom"]].append(iv)
    if int(r["u"]) >= GATE:
        atoms_expr[r["chrom"]].append(iv)


def index(d):
    return {c: (sorted(v), [x[0] for x in sorted(v)], max(e - s for s, e in v)) for c, v in d.items()}


IA, IE = index(atoms_all), index(atoms_expr)


def overlaps(I, c, s, e):
    if c not in I:
        return False
    v, st, mx = I[c]
    lo, hi = bisect.bisect_left(st, s - mx), bisect.bisect_left(st, e)
    return any(x[1] > s for x in v[lo:hi])


fields = None
dna = []
for r in csv.DictReader(open(dna_path), delimiter="\t"):
    fields = fields or list(r.keys())
    if expr.get((r["chrom"], int(r["start"]), int(r["end"])), 0) >= GATE:
        dna.append(r)
fc = collections.Counter(r["family_id"] for r in dna)
H = [r for r in dna if fc[r["family_id"]] >= 2]
rna = list(csv.DictReader(open(rna_path), delimiter="\t"))
for r in rna:
    r["_noatom"] = not overlaps(IA, r["chrom"], int(r["start"]), int(r["end"]))
    r["_noexpr"] = not overlaps(IE, r["chrom"], int(r["start"]), int(r["end"]))


def restricted(flag):
    kept = [r for r in rna if r[flag]]
    c = collections.Counter(r["family_id"] for r in kept)
    return [r for r in kept if c[r["family_id"]] >= 2]


whole_fams = {r["family_id"] for r in rna if r["_noexpr"]}
arms = {
    "H": [],
    "U1": restricted("_noatom"),
    "U2": restricted("_noexpr"),
    "U3": [r for r in rna if r["family_id"] in whole_fams],
}
cols = ["family_id", "copy_idx", "tid", "chrom", "start", "end", "n_exon", "strand", "n_reads", "exons"]
for name, extra in arms.items():
    with open(f"{out}.{name}.copies.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in H:
            fh.write("\t".join(str(r.get(k, "NA")) for k in cols) + "\n")
        for r in extra:
            row = dict(r)
            row["family_id"] = "RNA:" + r["family_id"]
            fh.write("\t".join(str(row.get(k, "NA")) for k in cols) + "\n")
    print(f"[union] {name}: {len({r['family_id'] for r in H})} atom families ({len(H)} copies) + "
          f"{len({r['family_id'] for r in extra})} RNA families ({len(extra)} copies)", file=sys.stderr)

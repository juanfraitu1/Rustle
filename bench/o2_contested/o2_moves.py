#!/usr/bin/env python3
"""r1652 (docs/PREREG_o2_moves_same_individual_2026-11-12.md): O2's decidable assigned rows split by whether its copy
keeps the primary record's locus or moves off it. Usage: o2_moves.py COPIES.tsv SPANS.tsv DUMP [DUMP2 ...]
DUMP = strict_truth STRICT_DUMP rows; SPANS = name chrom start end of primary records."""
import collections, csv, sys

def main(copies, spans, *dumps):
    span = {}
    for r in csv.DictReader(open(copies), delimiter="\t"): span[(r["family_id"], r["copy_idx"])] = (r["chrom"], int(r["start"]), int(r["end"]))
    if dumps[0].endswith(".assignments.tsv"):  # r1653 readout: no truth, count assigned rows kept on / moved off the primary record
        rec = {}; a = [r for d in dumps for r in csv.DictReader(open(d), delimiter="\t") if r["status"] == "assigned"]; need = {r["read_name"] for r in a}
        for l in open(spans):
            f = l.split("\t")
            if f[0] in need: rec[f[0]] = (f[1], int(f[2]), int(f[3]))
        c = collections.Counter(); fam = collections.Counter()
        for r in a:
            q = r["read_name"]; k = (r["family_id"], r["catalog_copy_idx"])
            if q not in rec or k not in span: c["no primary span or copy (skipped)"] += 1; continue
            mv = not (rec[q][0] == span[k][0] and rec[q][1] < span[k][2] and span[k][1] < rec[q][2]); c["moved" if mv else "kept"] += 1; fam[r["family_id"]] += mv
        print(dict(c), "moved share %.3f" % (c["moved"] / max(c["moved"] + c["kept"], 1)), "top families:", fam.most_common(5)); return
    rows = []
    for d in dumps:
        for l in open(d):
            f = l.rstrip("\n").split("\t")
            if f[2] == "assigned" and f[6] != "OUTSIDE": rows.append(f)
    need = {f[0] for f in rows}; rec = {}
    for l in open(spans):
        f = l.split("\t")
        if f[0] in need: rec[f[0]] = (f[1], int(f[2]), int(f[3]))
    hit = lambda r, k: r[0] == span[k][0] and r[1] < span[k][2] and span[k][1] < r[2]
    c = collections.Counter(); fam = collections.defaultdict(collections.Counter)
    for q, f, st, v, al, idx, t in rows:
        if q not in rec: c["no primary span (skipped)"] += 1; continue
        T = [tuple(x.split(":")) for x in t.split(";")]; pr = any(hit(rec[q], k) for k in T); mv = "kept" if hit(rec[q], (f, idx)) else "moved"
        c[f"{mv}: O2 {v}"] += 1; c[f"O2 {v}, primary {'right' if pr else 'wrong'}, {mv}"] += 1; fam[f][f"{mv} {v}"] += 1
    for k in sorted(c): print(c[k], k, sep="\t")
    m = c["moved: O2 right"] + c["moved: O2 wrong"]
    print(f"moved rows {m}, right share {c['moved: O2 right'] / m:.3f}" if m else "moved rows 0")
    print("families by moved-wrong:", [(f, dict(v)) for f, v in sorted(fam.items(), key=lambda x: -x[1]["moved wrong"])[:6]])

if __name__ == "__main__": main(*sys.argv[1:])

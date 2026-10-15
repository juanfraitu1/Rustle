#!/usr/bin/env python3
"""Final tally of the `copy` clusters of intact families (copyclusters.py classify + a species screen of the not_in_genome
clusters: the reads' best hit on CHM13 (HSA), mPanTro3 (PTR), mPonPyg2 (PPY) and KB3781's haplotypes; a cluster takes the species
with the highest median identity when it is >= 0.99; below that it is gorilla-like when gorilla is closest by >= 0.02, else
unresolved).

    copyclusters_summary.py X_DIR      (X_DIR = /mnt/linuxdisk/tmp/o3_realdata; writes X_DIR/copyclusters_summary.json)
"""
import collections
import csv
import json
import statistics
import sys

SPECIES_BAR = 0.99
GORILLA_MARGIN = 0.02     # gorilla-like only when gorilla beats every other species by this much


def best(p):
    b = {}
    for ln in open(p):
        f = ln.split("\t")
        if "tp:A:P" not in ln:
            continue
        i, c = int(f[9]) / max(1, int(f[10])), (int(f[3]) - int(f[2])) / max(1, int(f[1]))
        if f[0] not in b or i * c > b[f[0]][0] * b[f[0]][1]:
            b[f[0]] = (i, c)
    return b


def main(x):
    sp = {"GGO": [best(f"{x}/haps/mat.gb.paf"), best(f"{x}/haps/pat.gb.paf")], "HSA": [best(f"{x}/haps/nig.HSA.paf")],
          "PTR": [best(f"{x}/haps/nig.PTR.paf")], "PPY": [best(f"{x}/haps/nig.PPY.paf")]}

    def ident(r, s):
        v = [h[r][0] for h in sp[s] if r in h and h[r][1] >= 0.9]
        return max(v) if v else 0.0
    out = {}
    for s in ("heldout2", "dev"):
        rows = list(csv.DictReader(open(f"{x}/{s}/classified.tsv"), delimiter="\t"))
        mem = collections.defaultdict(list)
        for r in csv.DictReader(open(f"{x}/{s}/ab.psv_cluster_reads.tsv"), delimiter="\t"):
            mem[(r["family"], r["host"], r["cluster"])].append(r["read_name"])
        tally = collections.defaultdict(lambda: [0, 0, set()])
        for r in rows:
            cat = r["category"]
            k = (r["family"], r["host"], r["cluster"])
            if cat == "not_in_genome":
                med = {q: statistics.median(ident(m, q) for m in mem[k]) for q in sp}
                q, v = max(med.items(), key=lambda kv: kv[1])
                margin = med["GGO"] - max(w for q2, w in med.items() if q2 != "GGO")
                cat = (f"other_species_{q}" if q != "GGO" else "gorilla_outside_assemblies") if v >= SPECIES_BAR else \
                      ("gorilla_like_outside_assemblies" if q == "GGO" and v > 0 and margin >= GORILLA_MARGIN
                       else "outside_all_unresolved")
            t = tally[cat]
            t[0] += 1
            t[1] += int(r["n_reads"])
            t[2].add(r["family"])
        out[s] = {c: dict(clusters=v[0], reads=v[1], families=len(v[2])) for c, v in sorted(tally.items())}
        out[s]["_total"] = dict(clusters=len(rows), reads=sum(int(r["n_reads"]) for r in rows),
                                families=len({r["family"] for r in rows}))
    json.dump(out, open(f"{x}/copyclusters_summary.json", "w"), indent=1, sort_keys=True)
    print(json.dumps(out, indent=1, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))

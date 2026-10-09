#!/usr/bin/env python3
"""Three-way O3 call (COPY / ALLELE / PRESENT) of the clusters of an erasure bed: reference = the masked genome, truth = the unmasked genome. Descriptive (Amendment 26).
    o3_bed.py bedH <tag>      (miniforge python; bed H only: the masked-genome index exists for it)"""
import collections
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import apply_trim_real as AT  # noqa: E402
import flagmetric as FM  # noqa: E402
import run_polish as RP  # noqa: E402
import score as S  # noqa: E402
import seeds as SD  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
MASKED_IDX = "/mnt/linuxdisk/tmp/rna_allele/linktest/masked.splice.mmi"


def main(bed, tag):
    d = f"{W}/{bed}"
    out = f"{d}/{tag}"
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{out}/cons.fa").items()}
    fa = f"{out}/cons.o3.fa"
    SD.write_fa(fa, cons, list(cons))
    best = {}
    for which, idx in (("R", MASKED_IDX), ("T", RP.GENOME_IDX)):
        paf = f"{out}/cons.o3.{which}.paf"
        if not os.path.exists(paf + ".done"):
            subprocess.run(f"minimap2 -c -x splice:hq -uf -N 5 -t 4 {idx} {fa} > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
            open(paf + ".done", "w").write("ok")
        b = {}
        for ln in open(paf):
            f = ln.rstrip("\n").split("\t")
            if f[0] not in b or int(f[9]) > b[f[0]][0]:
                b[f[0]] = (int(f[9]), int(f[9]) / max(1, int(f[10])), int(f[2]), int(f[3]), int(f[1]), f[5], int(f[7]), int(f[8]))
        best[which] = b
    gate, _p, _sig = AT.library_gate()
    cl = {}
    for r in csv.DictReader(open(f"{out}/clusters.tsv"), delimiter="\t"):
        cl.setdefault("cl" + r["cluster"], []).append(r["read"])
    lab = {r["read"]: r for r in csv.DictReader(open(f"{d}/labels.tsv"), delimiter="\t")}
    truth = {n: ("bg" if x["role"] == "bg" else (x["family"] or None)) for n, x in lab.items()}
    erased = {p["fam"]: tuple(p["mask"][:3]) for p in json.load(open(RP.PANELS[bed]))}

    def sc(which, k):
        h = best[which].get(k)
        return FM.gate_aware(h[1], h[2], h[3], h[4], cons[k][:h[2]], gate) if h else None
    tab = collections.defaultdict(collections.Counter)
    rows = []
    for k, rs in cl.items():
        fam = S.majority(rs, truth)
        h = best["T"].get(k)
        on = bool(fam in erased and h and h[5] == erased[fam][0] and h[6] < erased[fam][2] and erased[fam][1] < h[7])
        kind = "erased-copy cluster" if on else ("other family cluster" if fam not in (None, "bg") else "background cluster")
        R, T = sc("R", k), sc("T", k)
        c = str(FM.o3_class(R, T))
        tab[kind][c] += 1
        tab[kind]["clusters"] += 1
        tab[kind]["registered_flag"] += FM.registered_flag(R, T)
        rows.append(dict(k=k, reads=len(rs), kind=kind, R=R, T=T, cls=c))
    json.dump(rows, open(f"{out}/o3_bed_classes.json", "w"), indent=1)
    print(f"{bed}/{tag}: reference = masked genome, truth = unmasked genome, gate {'OPEN' if gate else 'closed'}")
    print(f"{'cluster kind':22s} {'clusters':>8s} {'COPY':>5s} {'ALLELE':>6s} {'PRESENT':>7s} {'no call':>7s} {'registered flag':>15s}")
    for kind, c in tab.items():
        print(f"{kind:22s} {c['clusters']:>8d} {c['COPY']:>5d} {c['ALLELE']:>6d} {c['PRESENT']:>7d} {c['None']:>7d} {c['registered_flag']:>15d}")
    cop = [r for r in rows if r["kind"] == "erased-copy cluster" and r["cls"] == "ALLELE"]
    print("erased-copy clusters called ALLELE (the nearest survivor is within the allele cutoff): reads", sorted(r["reads"] for r in cop))


if __name__ == "__main__":
    main(*sys.argv[1:3])

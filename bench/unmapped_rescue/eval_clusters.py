#!/usr/bin/env python3
"""Evaluate a clustering of a real bed (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 22). Miniforge python.
    eval_clusters.py <bedA|bedH> <tag>      # tag = registered (frozen) | chain ; needs W/<bed>/<tag>/{clusters.tsv,metrics.json,cons.fa}
Prints: coverage and purity, clusters on the erased copy, unsupported consensus (Amendment 21), gate-aware identity x coverage >= 0.999, clusters with > 200 edits."""
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import apply_trim_real as AT  # noqa: E402
import consensus_support as CS  # noqa: E402
import flagmetric as FM  # noqa: E402
import run_polish as RP  # noqa: E402
import score as S  # noqa: E402
import seeds as SD  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"


def clusters_of(path):
    cl = {}
    for r in csv.DictReader(open(path), delimiter="\t"):
        cl.setdefault("cl" + r["cluster"], []).append(r["read"])
    return cl


def main(bed, tag):
    d = f"{W}/{bed}"
    out = f"{d}/{tag}"
    cl = clusters_of(f"{out}/clusters.tsv")
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{out}/cons.fa").items()}
    seqs = SD.read_fa(f"{d}/pool.fa")
    lab = {r["read"]: r for r in csv.DictReader(open(f"{d}/labels.tsv"), delimiter="\t")}
    truth = {n: ("bg" if x["role"] == "bg" else (x["family"] or None)) for n, x in lab.items()}
    erased = {p["fam"]: tuple(p["mask"][:3]) for p in json.load(open(RP.PANELS[bed]))}
    m = json.load(open(f"{out}/metrics.json"))
    fa = f"{out}/cons.for_genome.fa"
    SD.write_fa(fa, cons, list(cons))
    paf = f"{out}/cons.genome.paf"
    if not os.path.exists(paf + ".done"):
        subprocess.run(f"minimap2 -c --cs -x splice:hq -uf -N 5 -t 4 {RP.GENOME_IDX} {fa} > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
        open(paf + ".done", "w").write("ok")
    best = RP.best_records(open(paf).read().splitlines(), "")
    gate_ok, _p, _sig = AT.library_gate()
    n_on = ok = big = unsup = unsup_on = 0
    for k, rs in cl.items():
        fam = S.majority(rs, truth)
        names = sorted(rs)[:100]
        med = CS.median_read_divergence(cons[k], {n: seqs[n] for n in names}, f"/tmp/ec_{bed}_{tag}")
        bad = med is None or med > CS.DELTA
        unsup += bad
        if fam is None or fam not in erased or k not in best:
            continue
        h = best[k]
        e = erased[fam]
        if not (h["ref"] == e[0] and h["start"] < e[2] and e[1] < h["end"]):
            continue
        n_on += 1
        unsup_on += bad
        ok += FM.gate_aware(h["ident"], h["qstart"], h["qend"], h["qlen"], cons[k][:h["qstart"]], gate_ok) >= 0.999
        big += h["nm"] > 200
    print(f"{bed}/{tag}: clusters {len(cl)}; D reads clustered {m['d_clustered']} of {m['n_d']} ({m['coverage']:.1%}); purity {m['purity']:.4f}; "
          f"clusters on the erased copy {n_on}; unsupported consensus among them {unsup_on} (all clusters {unsup}); gate-aware >= .999: {ok}; > 200 edits: {big}")


if __name__ == "__main__":
    main(*sys.argv[1:3])

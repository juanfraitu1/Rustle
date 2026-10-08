#!/usr/bin/env python3
"""PART on a real bed, descriptive (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 18, real dev A / H). Run with /home/juanfra/miniforge3/bin/python.

    run_partition_real.py <bedA|bedH> run     # W/<bed>/registered/partition_real.json (resumable: exit 75 = run again)
    run_partition_real.py <bedA|bedH> score   # candidates aligned to the unmasked genome; baseline vs best candidate (oracle on the erased copy), gate-aware metric"""
import collections
import csv
import json
import os
import random
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import apply_trim_real as AT  # noqa: E402
import flagmetric as FM  # noqa: E402
import partition as PT  # noqa: E402
import run_polish as RP  # noqa: E402
import score as S  # noqa: E402
import seeds as SD  # noqa: E402

W = "/mnt/linuxdisk/tmp/o3_rescue"
BATCH = 420


def run(bed):
    d = f"{W}/{bed}"
    reg = f"{d}/registered"
    seqs = SD.read_fa(f"{d}/pool.fa")
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{reg}/cons.fa").items()}
    cl = RP.cluster_reads(d)
    out_path = f"{reg}/partition_real.json"
    done = json.load(open(out_path)) if os.path.exists(out_path) else {}
    align = PT.minimap_align_fn(f"{reg}/part_tmp", preset="splice:hq -uf")
    t0 = time.time()
    for k in sorted(cl, key=lambda x: len(cl[x])):
        if k in done:
            continue
        if time.time() - t0 > BATCH:
            json.dump(done, open(out_path, "w"))
            print(f"{len(done)} of {len(cl)} clusters done; run again")
            sys.exit(75)
        names = sorted(cl[k])
        if len(names) > 400:
            names = sorted(random.Random(1).sample(names, 400))
        leaves = PT.partition({n: seqs[n] for n in names}, align, PT.abpoa_consensus, cons=cons[k], n_as_del=True)
        done[k] = [dict(reads=lf["reads"], cons=lf["cons"]) for lf in sorted(leaves, key=lambda x: -len(x["reads"]))]
    json.dump(done, open(out_path, "w"))
    split = sum(1 for v in done.values() if len(v) > 1)
    print(f"{bed}: {len(done)} clusters, {split} split, {sum(len(v) for v in done.values())} candidates; leaves per split cluster {collections.Counter(len(v) for v in done.values() if len(v) > 1)}")


def score(bed):
    d = f"{W}/{bed}"
    reg = f"{d}/registered"
    part = json.load(open(f"{reg}/partition_real.json"))
    cl = RP.cluster_reads(d)
    cons = {n.split("|")[0]: s for n, s in RP.read_cons(f"{reg}/cons.fa").items()}
    cands = {}
    for k, leaves in part.items():
        for j, lf in enumerate(leaves):
            cands[f"{k}|p{j}"] = lf["cons"]
    fa = f"{reg}/cands.fa"
    SD.write_fa(fa, cands, list(cands))
    paf = f"{reg}/cands.paf"
    if not os.path.exists(paf + ".done"):
        subprocess.run(f"minimap2 -c --cs -x splice:hq -uf -N 5 -t 4 {RP.GENOME_IDX} {fa} > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
        open(paf + ".done", "w").write("ok")
    lines = open(paf).read().splitlines()
    best = {}
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        q = f[0]
        m = int(f[9])
        if q not in best or m > best[q]["m"]:
            best[q] = dict(m=m, ident=m / max(1, int(f[10])), ref=f[5], start=int(f[7]), end=int(f[8]), qs=int(f[2]), qe=int(f[3]), qlen=int(f[1]), nm=int([t for t in f[12:] if t.startswith("NM:i:")][0][5:]))
    lab = {r["read"]: r for r in csv.DictReader(open(f"{d}/labels.tsv"), delimiter="\t")}
    truth = {n: ("bg" if x["role"] == "bg" else (x["family"] or None)) for n, x in lab.items()}
    erased = {p["fam"]: tuple(p["mask"][:3]) for p in json.load(open(RP.PANELS[bed]))}
    gate_ok, _p, _sig = AT.library_gate()
    sc = json.load(open(f"{reg}/polish/score.json"))["clusters"]
    rows = []
    for k, leaves in part.items():
        fam = S.majority(cl[k], truth)
        if fam is None or fam not in erased or k not in sc or not sc[k]["orig"]["on"]:
            continue
        e = erased[fam]

        def idg(name):
            h = best.get(name)
            if not h or not (h["ref"] == e[0] and h["start"] < e[2] and e[1] < h["end"]):
                return None, None
            return FM.gate_aware(h["ident"], h["qs"], h["qe"], h["qlen"], cands[name][:h["qs"]], gate_ok), h["nm"]
        base_id, base_nm = idg(f"{k}|p0") if len(leaves) == 1 else (None, None)
        if len(leaves) > 1:       # the baseline is the frozen consensus, not any leaf: align-free lookup from the Amendment 9 table
            h0 = RP.best_records(open(f"{reg}/polish/both.paf").read().splitlines(), "orig.").get(k)
            base_id = FM.gate_aware(h0["ident"], h0["qstart"], h0["qend"], h0["qlen"], cons[k][:h0["qstart"]], gate_ok)
            base_nm = h0["nm"]
        cand = [idg(f"{k}|p{j}") for j in range(len(leaves))]
        ok = [c for c in cand if c[0] is not None]
        rows.append(dict(k=k, leaves=len(leaves), size=len(cl[k]), base=base_id, base_nm=base_nm, best=max((c[0] for c in ok), default=None), best_nm=min((c[1] for c in ok), default=None)))
    split = [r for r in rows if r["leaves"] > 1]
    ge = lambda x: x is not None and x >= 0.999
    print(f"{bed}: clusters on the erased copy {len(rows)}; split {len(split)} ({collections.Counter(r['leaves'] for r in split)})")
    print(f"  gate-aware identity x coverage >= 0.999: frozen consensus {sum(ge(r['base']) for r in rows)}; best candidate (oracle) {sum(ge(r['best']) or ge(r['base']) for r in rows)}")
    print(f"  among the {len(split)} split clusters: frozen {sum(ge(r['base']) for r in split)}, best candidate {sum(ge(r['best']) for r in split)}; "
          f"edit distance to the erased copy, frozen {sum(r['base_nm'] or 0 for r in split)} vs best candidate {sum((r['best_nm'] if r['best_nm'] is not None else r['base_nm'] or 0) for r in split)}")
    big = [r for r in rows if (r["base_nm"] or 0) > 200]
    print(f"  clusters with > 200 edits for the frozen consensus: {len(big)}, split {sum(r['leaves'] > 1 for r in big)}; best candidate at 0.999: {sum(ge(r['best']) for r in big)}")
    json.dump(rows, open(f"{reg}/partition_real_score.json", "w"), indent=1)


if __name__ == "__main__":
    {"run": run, "score": score}[sys.argv[2]](sys.argv[1])

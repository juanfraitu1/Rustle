#!/usr/bin/env python3
"""Whose reads are the unmapped clusters? (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 34). Miniforge python.

    species_check.py pool              pooled consensus of the ELSEWHERE / NOVEL / UNSUPPORTED clusters of the four real runs
    species_check.py align <assembly>  HSA GGO Tm PTR PPY (resumable, one assembly per call)
    species_check.py report"""
import collections
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import discover as D  # noqa: E402
import run_polish as RP  # noqa: E402
import seeds as SD  # noqa: E402

VARIANT = os.environ.get("SPECIES_VARIANT", "")      # "" = the global-mode consensus of Amendment 34, "l" = the adopted local-mode consensus (Amendment 40)
OUT = f"{D.W}/species" + (f"_{VARIANT}" if VARIANT else "")
ANIMALS = ("ggo_testis", "a119b", "ptr", "ppy")
OWN = {"ggo_testis": "GGO", "a119b": "HSA", "ptr": "PTR", "ppy": "PPY"}
ASSEMBLY = {"HSA": D.CHM13_IDX, "GGO": D.RP.GENOME_IDX, "Tm": D.C.HAP_IDX.format("mat"), "PTR": D.PTR_FA, "PPY": D.PPY_FA}
SPECIES = {"HSA": "HSA", "GGO": "GGO", "Tm": "GGO", "PTR": "PTR", "PPY": "PPY"}
CLASSES = ("ELSEWHERE", "NOVEL", "UNSUPPORTED")
BAR = 0.99


def paf_scores(lines):
    """-> {query: min(1, sum of matching bases over the primary records / query length)}"""
    tot, qlen = collections.defaultdict(int), {}
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        if "tp:A:P" not in f[12:]:
            continue
        tot[f[0]] += int(f[9])
        qlen[f[0]] = int(f[1])
    return {q: min(1.0, tot[q] / qlen[q]) for q in tot}


def best_species(by_species):
    """-> (species, score); a missing score is 0"""
    best = max(by_species, key=lambda s: by_species[s] or 0.0, default=None)
    sc = (by_species.get(best) or 0.0) if best else 0.0
    return (best, sc) if sc > 0 else (None, 0.0)


def same_species(best, own, bar=BAR):
    return best[0] == own and best[1] >= bar


def read_share(rows, pred):
    tot = sum(r["reads"] for r in rows)
    return sum(r["reads"] for r in rows if pred(r)) / tot if tot else None


def pool():
    os.makedirs(OUT, exist_ok=True)
    seqs, names = {}, []
    for a in ANIMALS:
        vp = D.variant_paths(f"{D.W}/discover_{a}", VARIANT)
        cons = {n.split("|")[0]: s for n, s in RP.read_cons(vp["cons"]).items()}
        for r in json.load(open(vp["rescored"] if VARIANT else vp["classes"])):
            if r["cls"] in CLASSES:
                seqs[f"{a}__{r['k']}"] = cons[r["k"]]
                names.append(f"{a}__{r['k']}")
    SD.write_fa(f"{OUT}/pool.fa", seqs, names)
    print(f"pooled {len(seqs)} consensus sequences")


def align(asm):
    paf = f"{OUT}/pool.{asm}.paf"
    if os.path.exists(paf + ".done"):
        print(f"{asm}: done")
        return
    subprocess.run(f"minimap2 -c -x splice -uf -N 5 -t 4 {ASSEMBLY[asm]} {OUT}/pool.fa > {paf}", shell=True, check=True, stderr=subprocess.DEVNULL)
    open(paf + ".done", "w").write("ok")
    print(f"{asm}: aligned")


def report():
    sc = {a: paf_scores(open(f"{OUT}/pool.{a}.paf")) for a in ASSEMBLY}
    rows = []
    for a in ANIMALS:
        vp = D.variant_paths(f"{D.W}/discover_{a}", VARIANT)
        for r in json.load(open(vp["rescored"] if VARIANT else vp["classes"])):
            if r["cls"] not in CLASSES:
                continue
            q = f"{a}__{r['k']}"
            by = {}
            for asm, s in sc.items():
                by[SPECIES[asm]] = max(by.get(SPECIES[asm]) or 0.0, s.get(q, 0.0))
            b = best_species(by)
            rows.append(dict(animal=a, cls=r["cls"], reads=r["reads"], length=r["length"], by=by, best=b, own=same_species(b, OWN[a]), k=r["k"]))
    json.dump(rows, open(f"{OUT}/species.json", "w"), indent=1)
    f2 = lambda v: round(v, 3)
    for a in ANIMALS:
        mine = [r for r in rows if r["animal"] == a]
        for cls in (("ELSEWHERE",), ("NOVEL",), ("UNSUPPORTED",), ("ELSEWHERE", "NOVEL")):
            sel = [r for r in mine if r["cls"] in cls]
            if not sel:
                continue
            comp = collections.Counter()
            for r in sel:
                comp[(r["best"][0] if r["best"][1] >= BAR else ("<0.9" if r["best"][1] < 0.9 else "0.9-0.99"))] += r["reads"]
            print(f"{a:11s} {'+'.join(cls):18s} clusters {len(sel):3d} reads {sum(r['reads'] for r in sel):4d}  same-species share {read_share(sel, lambda r: r['own']):.2f}  "
                  f"best species by reads (>= {BAR}, else band): {dict(comp)}")
    for a in ANIMALS:
        top = sorted((r for r in rows if r["animal"] == a), key=lambda r: -r["reads"])[:6]
        print(a, [(r["reads"], r["length"], r["cls"][:4], {s: f2(v) for s, v in r["by"].items()}) for r in top])
    p1 = read_share([r for r in rows if r["animal"] == "ggo_testis" and r["cls"] == "ELSEWHERE"], lambda r: r["own"])
    hum = [r for r in rows if r["animal"] == "a119b" and r["cls"] in ("ELSEWHERE", "NOVEL")]
    p2 = read_share(hum, lambda r: r["best"][0] not in (None, "HSA") and r["best"][1] >= BAR)
    print(f"P1 (gorilla ELSEWHERE reads in same-species clusters, bar >= 0.80): {p1 if p1 is None else round(p1, 3)}")
    print(f"P2 (human ELSEWHERE+NOVEL reads in clusters best on a NON-human assembly >= {BAR}, bar >= 0.50): {p2 if p2 is None else round(p2, 3)}")


if __name__ == "__main__":
    {"pool": lambda: pool(), "align": lambda: align(sys.argv[2]), "report": lambda: report()}[sys.argv[1]]()

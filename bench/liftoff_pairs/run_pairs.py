#!/usr/bin/env python3
"""Liftoff copy pairs against matched decoys, per sample (prereg docs/PREREG_liftoff_pairs_null_2026-10-14.md).

    run_pairs.py SAMPLE --out DIR [--null 1000] [--seed 7] [--inputs-only]
    SAMPLE = human_A119b | human_testis | gorilla_KB3781

DIR receives inputs.json (claim F1 recomputed with figure 8's own code, the read-supported pair counts and the decoy pools;
all that --inputs-only writes) and test.json (the decoy test per scope and sequence-identity floor; the bar is the
genome-wide row at identity >= 0.95). The run stops when F1 does not reproduce figures/data/fig8_samples.tsv.
"""
import argparse
import csv
import json
import os
import sys
from pathlib import Path

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
import pairs_lib as P  # noqa: E402

LF = P.LF
SAMPLES = {"human_A119b": "human", "human_testis": "human", "gorilla_KB3781": "gorilla"}
LIFTOFF = "/mnt/linuxdisk/tmp/rustle_figures/liftoff"
RUNS = "/mnt/linuxdisk/tmp/rustle_figures/runs"
FIG8 = os.path.join(REPO, "figures", "data", "fig8_samples.tsv")
SCOPES = ("genome", "minus_dev")
PRIMARY = ("genome", "0.95")


def keep_for(species, scope):
    """None (every contig) for the genome; for minus_dev, the contigs the family rule was not developed on."""
    return None if scope == "genome" else (lambda c: c not in LF.DEV_CONTIGS.get(species, set()))


def f1(loci, fam, species):
    """Claim F1 exactly as `figures/fig_loci._pairs` computes it: Liftoff pairs (no read-support condition) whose two loci
    are both covered by family members; k = those sharing a family. -> {scope: {sc: [n, k]}}."""
    out = {}
    for scope in SCOPES:
        keep = keep_for(species, scope) or (lambda c: True)
        out[scope] = {}
        for sc in LF.SC_ROWS:
            res = LF.pair_families(LF.copy_pairs(loci, sc, keep=keep), fam)
            hits = [int(bool(shared)) for _, both, shared in res if both]
            out[scope][f"{sc:.2f}"] = [len(hits), sum(hits)]
    return out


def fig8_f1(path, sample):
    """F1's rows of figures/data/fig8_samples.tsv for one sample. -> {scope: {sc: [n, k]}}."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader((line for line in fh if not line.startswith("#")), delimiter="\t"):
            if r["sample"] == sample and r["rustle_set"] == "families" and r["measure"] == "pair_cofamily":
                out.setdefault(r["scope"], {})[r["sc_min"]] = [int(r["n"]), int(r["k"])]
    return out


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True)
        fh.write("\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample")
    ap.add_argument("--out", required=True)
    ap.add_argument("--null", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--inputs-only", action="store_true", help="write inputs.json and stop (no decoy statistic)")
    ap.add_argument("--liftoff", default=LIFTOFF)
    ap.add_argument("--runs", default=RUNS)
    ap.add_argument("--fig8", default=FIG8)
    a = ap.parse_args(argv)
    if a.sample not in SAMPLES:
        ap.error(f"unknown sample {a.sample}; one of {', '.join(sorted(SAMPLES))}")
    species = SAMPLES[a.sample]
    loci = LF.read_loci(Path(a.liftoff) / species / "liftoff_loci.tsv")
    sup = LF.read_support(Path(a.liftoff) / "support" / f"{a.sample}.tsv")
    fam = LF.catalog_loci(Path(a.runs) / a.sample / f"{a.sample}.fam.copies.tsv")
    got, want = f1(loci, fam, species), fig8_f1(a.fig8, a.sample)
    if got != want:
        raise SystemExit(f"F1 does not reproduce {a.fig8} for {a.sample}: recomputed {got}, table {want}")
    pairs = {scope: {f"{sc:.2f}": LF.copy_pairs(loci, sc, sup=sup, keep=keep_for(species, scope)) for sc in LF.SC_ROWS}
             for scope in SCOPES}
    pools = {scope: P.decoy_pool(loci, sup, keep_for(species, scope)) for scope in SCOPES}
    os.makedirs(a.out, exist_ok=True)
    inputs = dict(species=species, f1=got, pairs={s: {sc: len(v) for sc, v in d.items()} for s, d in pairs.items()},
                  pool={s: len(v) for s, v in pools.items()})
    write_json(os.path.join(a.out, "inputs.json"), dict(args=vars(a), **inputs))
    if a.inputs_only:
        print(json.dumps(inputs, indent=1, sort_keys=True))
        return 0
    rows = {LF.locus_key(r): r for r in pools["genome"]}
    rows.update((LF.locus_key(r), r) for d in pairs.values() for prs in d.values() for pair in prs for r in pair)
    fams = P.cover_families(list(rows.values()), fam)
    out = {}
    for scope in SCOPES:
        pool = P.Pool(pools[scope], fams)
        out[scope] = {sc: P.run_test(prs, pool, fams, a.null, a.seed) for sc, prs in pairs[scope].items()}
    prim = out[PRIMARY[0]][PRIMARY[1]]
    write_json(os.path.join(a.out, "test.json"),
               dict(out, primary={"scope": PRIMARY[0], "sc_min": PRIMARY[1]}, passes=bool(prim["p"] < 0.001)))
    print(json.dumps(dict(sample=a.sample, n=prim["n"], k=prim["k"], rate=prim["rate"], null_mean=prim["null_mean"],
                          null_max=prim["null_max"], p=prim["p"], passes=bool(prim["p"] < 0.001)), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

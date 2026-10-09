#!/usr/bin/env python3
"""Read-level stress harness for PART (docs/PREREG_unmapped_rescue_2026-10-08.md Amendment 25). Run with /home/juanfra/miniforge3/bin/python.

    stress_part.py <outdir> <seed> [n_transcripts=12]      GUARD=0 turns the block-purity guard off (to see the failure it fixes)
Scenarios per transcript (100 reads, end jitter 3, HiFi substitutions 0.001, short indels 0.0003): control (independent errors only); hotspot (75 good reads with p 0.03 and 25
low-quality reads with p at F fixed fragile columns, each a 1-bp deletion); skip (an internal segment of 90-350 bases absent from 33% of the reads); sib (a second
sequence 0.5% apart in 50% of the reads); het (two linked SNPs in 50% of the reads). Prints, for each, how many transcripts were split."""
import collections
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
import partition as PT  # noqa: E402
from sim import simulate_reads  # noqa: E402

MUT = {"A": "C", "C": "G", "G": "T", "T": "A"}


def reads_of(seq, n, seed, prefix, jitter=3, err=0.001, indel=0.0003, fragile=(), p=0.0):
    rng = random.Random(seed)
    out = {}
    for i in range(n):
        drop = {c for c in fragile if rng.random() < p}
        body = "".join(b for j, b in enumerate(seq) if j not in drop)
        lo, hi = rng.randint(0, jitter), len(body) - rng.randint(0, jitter)
        out[f"{prefix}{i}"] = simulate_reads(body[lo:hi], 1, err=err, indel=indel, seed=seed * 100003 + i)[0][0]
    return out


def mutate(seq, cols):
    s = list(seq)
    for c in cols:
        s[c] = MUT[s[c]]
    return "".join(s)


def scenario_reads(kind, t, rng, seed):
    L = len(t)
    if kind == "control":
        return reads_of(t, 100, seed, "c")
    if kind.startswith("hotspot"):
        F, p = (int(x) if i == 0 else float(x) for i, x in enumerate(kind.split("_")[1:3]))
        fragile = sorted(rng.sample(range(30, L - 30), F))
        return {**reads_of(t, 75, seed, "g", fragile=fragile, p=0.03), **reads_of(t, 25, seed + 1, "q", fragile=fragile, p=p)}
    if kind == "skip":
        a = rng.randint(400, L - 600)
        w = rng.randint(90, 350)
        return {**reads_of(t, 67, seed, "a"), **reads_of(t[:a] + t[a + w:], 33, seed + 1, "b")}
    if kind == "sib":
        cols = sorted(rng.sample(range(40, L - 40), max(7, round(0.005 * L))))
        return {**reads_of(t, 50, seed, "a"), **reads_of(mutate(t, cols), 50, seed + 1, "b")}
    if kind == "het":
        cols = sorted(rng.sample(range(40, L - 40), 2))
        return {**reads_of(t, 50, seed, "a"), **reads_of(mutate(t, cols), 50, seed + 1, "b")}
    raise ValueError(kind)


KINDS = ["control", "hotspot_20_0.3", "hotspot_20_0.5", "hotspot_40_0.5", "skip", "sib", "het"]


def main(out, seed, n=12):
    os.makedirs(out, exist_ok=True)
    if os.environ.get("GUARD", "1") == "0":
        PT.pure_blocks = lambda blocks, *a, **k: blocks
    align = PT.minimap_align_fn(f"{out}/tmp", preset="splice:hq -uf")
    res = collections.defaultdict(list)
    for i in range(n):
        rng = random.Random(seed * 1000 + i)
        t = "".join(rng.choice("ACGT") for _ in range(rng.randint(1500, 3000)))
        for kind in KINDS:
            reads = scenario_reads(kind, t, rng, seed * 1000 + i * 10)
            leaves = PT.partition(reads, align, PT.abpoa_consensus, n_as_del=True)
            res[kind].append(len(leaves))
    for kind in KINDS:
        v = res[kind]
        print(f"{kind:16s} split {sum(x > 1 for x in v):2d} of {len(v)}  leaves {sorted(v)}")
    json.dump(res, open(f"{out}/stress.json", "w"))


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 12)

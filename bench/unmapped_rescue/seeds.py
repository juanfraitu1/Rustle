#!/usr/bin/env python3
"""Seed-based clustering of a redundant read pool (docs/PREREG_unmapped_rescue_2026-10-08.md section 1, steps 1-2).

An all-vs-all of a pool in which one transcript has thousands of near-identical reads is quadratic. Instead each round takes `n_seeds` reads spread
over the length order, maps EVERY unclustered read to the seeds only, and keeps the registered edge rule (graph.edges); components of size >= min_size are
done, the rest go to the next round with new seeds. Rounds stop when none is left or one makes no progress.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import graph as G  # noqa: E402


def choose_seeds(lens, k):
    """k read names evenly spread over the descending length order (deterministic): lens = {read: length}"""
    order = sorted(lens, key=lambda r: (-lens[r], r))
    if len(order) <= k:
        return order
    step = len(order) / k
    return [order[int(i * step)] for i in range(k)]


def run_rounds(lens, map_fn, n_seeds=2000, min_size=3, max_rounds=6):
    """lens: {read: length}; map_fn(queries, seed_names) -> [(a, b)] edges that passed the edge rule. Returns {read: component id}."""
    edges = []
    todo = dict(lens)
    comp = G.components([], list(lens))
    for _ in range(max_rounds):
        if not todo:
            break
        seed_names = choose_seeds(todo, n_seeds)
        new = list(map_fn(sorted(todo), seed_names))
        edges += new
        comp = G.components(edges, list(lens))
        cl = G.clusters(comp, min_size)
        done = {r for rs in cl.values() for r in rs}
        left = {r: l for r, l in todo.items() if r not in done}
        if len(left) == len(todo):                      # no progress
            break
        todo = left
    return comp


def write_fa(path, seqs, names):
    with open(path, "w") as o:
        for n in names:
            o.write(f">{n}\n{seqs[n]}\n")


def read_fa(path):
    seqs, cur = {}, None
    for ln in open(path):
        if ln[0] == ">":
            cur = ln[1:].strip().split()[0]
            seqs[cur] = []
        else:
            seqs[cur].append(ln.strip())
    return {k: "".join(v) for k, v in seqs.items()}


class Pause(Exception):
    """a fresh minimap2 round finished and the call's budget is used up: run the program again (finished rounds are read back from disk)"""


def minimap_map_fn(seqs, work, delta, min_frac, threads=4, preset="map-hifi", fresh_budget=1, proper=False):
    """map_fn for run_rounds: minimap2 -x <preset> -c of the queries against the seed set. A round whose PAF and `.done` marker exist is read back;
    after `fresh_budget` fresh rounds in this call the NEXT fresh round raises Pause (each call stays inside the foreground window)."""
    state = {"i": 0, "fresh": 0}

    def f(queries, seed_names):
        i = state["i"] = state["i"] + 1
        sf, qf, pf = f"{work}/seeds{i}.fa", f"{work}/query{i}.fa", f"{work}/round{i}.paf"
        if not os.path.exists(pf + ".done"):
            if state["fresh"] >= fresh_budget:
                raise Pause(i)
            state["fresh"] += 1
            write_fa(sf, seqs, seed_names)
            write_fa(qf, seqs, queries)
            subprocess.run(f"minimap2 -x {preset} -c -N 20 -t {threads} {sf} {qf} > {pf} 2> {pf}.log", shell=True, check=True)
            open(pf + ".done", "w").write("ok")
        with open(pf) as fh:
            return list(G.edges(fh, delta, min_frac, proper))
    return f

#!/usr/bin/env python3
"""Amendment 42: locus-binned clustering of a genome-scale net. Reads are binned by their own primary alignment; inside a bin, the adopted star step on the
all-vs-all of the `cap` longest unassigned reads, the other reads join the first centre they pass the same edge with, leftovers get another pass.

Pure (tested in test_binclust.py): bin_reads, pick_centres, cluster_bin. minimap2-backed helpers: compat_fn, assign_fn."""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chain as CH  # noqa: E402
import graph as G  # noqa: E402

DELTA = 0.00958


def bin_reads(recs):
    """recs = [(name, contig, start, end, strand)] -> [[names]]: reads on one contig and strand whose reference spans overlap (transitively) share a bin"""
    bins, cur, key, end = [], None, None, None
    for n, c, s, e, st in sorted(recs, key=lambda r: (r[1], r[4], r[2], r[3], r[0])):
        if cur is not None and (c, st) == key and s < end:
            cur.append(n)
            end = max(end, e)
        else:
            if cur is not None:
                bins.append(cur)
            cur, key, end = [n], (c, st), e
    if cur is not None:
        bins.append(cur)
    return bins


def pick_centres(compatible, centres):
    """compatible = {read: set of centres it passes the edge with}; -> {read: the first of them in `centres` order}"""
    rank = {c: i for i, c in enumerate(centres)}
    return {r: min(cs, key=lambda c: rank[c]) for r, cs in compatible.items() if cs}


def cluster_bin(names, lens, compat_fn, assign_fn, cap=500, max_pass=20, min_size=3, state=None, stop=lambda: False):
    """-> {centre: [reads]}. compat_fn(head) -> set of compatible sorted pairs among head; assign_fn(rest, centres) -> {read: centre}.
    state (a dict, filled in place) makes a large bin resumable pass by pass: when stop() is true after a pass, the pass is saved in state and None is returned."""
    st = state if state is not None else {}
    if not st:
        st.update(clusters={}, left=sorted(names, key=lambda r: (-lens[r], r)), passes=0, finished=False)
    while not st["finished"]:
        left = st["left"]
        if st["passes"] >= max_pass or len(left) < min_size:
            st["finished"] = True
            break
        head, rest = left[:cap], left[cap:]
        stars = CH.star_clusters({r: lens[r] for r in head}, compat_fn(head), min_size)
        new = {s[0]: list(s) for s in stars}
        st["passes"] += 1
        if not new:
            st["finished"] = True
            break
        if rest:
            for r, c in assign_fn(rest, [s[0] for s in stars]).items():
                new[c].append(r)
        st["clusters"].update(new)
        done = {r for v in new.values() for r in v}
        st["left"] = [r for r in left if r not in done]
        if not rest:
            st["finished"] = True
            break
        if stop():
            return None
    return st["clusters"]


def _write(path, seqs, names):
    with open(path, "w") as o:
        for n in names:
            o.write(f">{n}\n{seqs[n]}\n")


def compat_fn(seqs, work, threads=2):
    """the Amendment 24 edge (proper overlap, NM / block <= delta) over the all-vs-all of a read set"""
    os.makedirs(work, exist_ok=True)
    allvsall = CH.minimap_allvsall(work, threads)

    def f(head):
        return set(G.edges(allvsall({r: seqs[r] for r in head}), DELTA, 0.5, proper=True, edit=True))
    return f


def assign_fn(seqs, work, threads=2):
    """reads mapped to the star centres (map-hifi, -N 500 -p 0.1), the same edge; each read goes to the first passing centre in centre order"""
    os.makedirs(work, exist_ok=True)

    def f(rest, centres):
        _write(f"{work}/c.fa", seqs, centres)
        _write(f"{work}/r.fa", seqs, rest)
        lines = subprocess.run(f"minimap2 -x map-hifi -c -N 500 -p 0.1 -t {threads} {work}/c.fa {work}/r.fa", shell=True, stdout=subprocess.PIPE,
                               stderr=subprocess.DEVNULL, text=True).stdout.splitlines()
        cset, comp = set(centres), {}
        for a, b in G.edges(lines, DELTA, 0.5, proper=True, edit=True):
            r, c = (a, b) if b in cset and a not in cset else (b, a)
            if c in cset and r not in cset:
                comp.setdefault(r, set()).add(c)
        return pick_centres(comp, centres)
    return f

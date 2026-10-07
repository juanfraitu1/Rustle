#!/usr/bin/env python3
"""Pool the arms.json of every run under a work directory (FAMILY/repR/arms.json) per arm and stratum: chains recovered, genes complete / resolved, artifact load.

    pool.py WORKDIR
Prints one block per arm (all four runs summed, then NPIP and TBC1D3 separately) and writes WORKDIR/pooled.json.
"""
import collections
import glob
import json
import os
import sys


def main():
    w = sys.argv[1]
    runs = {}
    for p in sorted(glob.glob(os.path.join(w, "*", "rep*", "arms.json"))):
        fam, rep = p.split(os.sep)[-3], p.split(os.sep)[-2]
        runs[(fam, rep)] = json.load(open(p))
    if not runs:
        raise SystemExit("no arms.json under " + w)
    arms = sorted({a for r in runs.values() for a in r["arms"]})
    out = {}
    for scope, keep in (("ALL", lambda k: True), ("NPIP", lambda k: k[0] == "NPIP"), ("TBC1D3", lambda k: k[0] == "TBC1D3")):
        sel = {k: v for k, v in runs.items() if keep(k)}
        print(f"== {scope} ({len(sel)} runs: {', '.join(f'{f}/{r}' for f, r in sorted(sel))})")
        print(f"{'arm':7} {'stratum':3} {'genes':>6} {'chains':>7} {'recovered':>9} {'%':>6} {'complete':>9} {'resolved':>9}   artifacts(total frag/fusion/other)  merged_ids  tx(multi)")
        out[scope] = {}
        for a in arms:
            tot = collections.defaultdict(lambda: collections.Counter())
            art = collections.Counter(); merged = tx = multi = 0
            for r in sel.values():
                if a not in r["arms"]:
                    continue
                d = r["arms"][a]
                for s in ("E", "N"):
                    for k, v in d["strata"][s].items():
                        tot[s][k] += v
                for k, v in d["artifacts"].items():
                    art[k] += v
                merged += d["merged_ids"]; tx += d["transcripts"]; multi += d["multi_exon"]
            if not tot:
                continue
            out[scope][a] = dict(strata={s: dict(c) for s, c in tot.items()}, artifacts=dict(art), merged_ids=merged, transcripts=tx, multi_exon=multi)
            for s in ("E", "N"):
                c = tot[s]
                pct = 100.0 * c["recovered"] / c["chains"] if c["chains"] else 0.0
                extra = f"   {sum(art.values()):6} {art['fragment']}/{art['fusion']}/{art['other']}   {merged:6}   {tx}({multi})" if s == "E" else ""
                print(f"{a:7} {s:3} {c['genes']:6} {c['chains']:7} {c['recovered']:9} {pct:6.1f} {c['complete']:9} {c['resolved']:9}{extra}")
    json.dump(out, open(os.path.join(w, "pooled.json"), "w"), indent=1)


if __name__ == "__main__":
    main()

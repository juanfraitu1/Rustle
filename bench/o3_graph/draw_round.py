#!/usr/bin/env python3
"""A fresh family draw for a simulation round: N qualifying non-PAR pool families no earlier families*.json set used.
    draw_round.py ROUND N SEED   -> /mnt/linuxdisk/tmp/o3_graph/families<ROUND>.json (set r<ROUND>)"""
import glob
import hashlib
import json
import random
import sys

T = "/mnt/linuxdisk/tmp/o3_graph"


def main(rnd, n, seed):
    q = json.load(open(f"{T}/famident/qualify_nonpar.json"))
    used = set()
    for p in glob.glob(f"{T}/families*.json"):
        try:
            d = json.load(open(p))
        except Exception:
            continue
        for s in d.get("sets", {}).values():
            used |= {f["family"] for f in s}
    pool = json.load(open(f"{T}/r19seq/families19cal_full2_all.json"))
    pf = {f["family"]: f for s in pool["sets"].values() for f in s}
    fresh = sorted(f for f in q if f not in used and f in pf)
    random.seed(int(seed))
    draw = sorted(random.sample(fresh, int(n)))
    out = {k: v for k, v in pool.items() if k != "sets"}
    out["sets"] = {f"r{rnd}": [pf[f] for f in draw]}
    p = f"{T}/families{rnd}.json"
    json.dump(out, open(p, "w"))
    h = hashlib.sha1(open(p, "rb").read()).hexdigest()[:12]
    print(f"{p} sha1 {h}: {len(draw)} of {len(fresh)} fresh families, {sum(len(pf[f]['copies']) for f in draw)} copies")
    print(" ".join(draw))
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))

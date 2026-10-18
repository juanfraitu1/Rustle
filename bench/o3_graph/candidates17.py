#!/usr/bin/env python3
"""Round 17's family pool (2026-10-17), fixed before any build: every CAT protein-coding gene-name STEM (the name less a final
number and optional letter) with EXACTLY 2 genes, neither matched by any earlier group's pattern (families.GROUPS); round 15's
pool (>= 3 genes) is spent, and every such stem is disjoint from it. Built with round 12's rule; a stem forms a family only where
families.py links both copies (identity >= 0.90 over >= 0.8 of a locus).

    candidates17.py OUT_DIR [--build]  -> OUT_DIR/candidates17.json (the stems), OUT_DIR/families17_all.json (the build)
"""
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, f"{name}.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def pairs(all_stems, at_least_3):
    """The stems of exactly 2 genes: those with >= 2 less those with >= 3 (candidates15.stems' two calls)."""
    return {s: p for s, p in all_stems.items() if s not in at_least_3}


def main(argv):
    out = argv[0]
    os.makedirs(out, exist_ok=True)
    F, C = _load("families"), _load("candidates15")
    used = [p for g in F.GROUPS.values() for p in g.values()]
    gff = f"{F.CAT}/chm13v2.0_CAT_Liftoff.slim.gff3.gz"
    cand = pairs(C.stems(gff, used, min_genes=2), C.stems(gff, used, min_genes=3))
    json.dump(cand, open(f"{out}/candidates17.json", "w"), indent=0)
    print(f"[candidates17] {len(cand)} stems -> {out}/candidates17.json", flush=True)
    if "--build" in argv:
        F.GROUPS = {"auto17": cand}
        return F.main(["build", f"{out}/families17_all.json", "--cap", "1000", "--link", "locus", "--min-cover", "0.8"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

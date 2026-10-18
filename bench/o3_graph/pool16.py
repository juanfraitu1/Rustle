"""Round 16's held-out families (2026-10-17): every family of round 15's pool (r15/families15_all.json, set auto15: every
unused CAT protein-coding stem with >= 3 genes, built by round 12's rule) that has >= 3 copies and that round 15 did not draw
(families15.json, set r15). No draw, hence no selection. Writes OUT with the pool's build settings and one set, r16.
Usage: pool16.py POOL ROUND15 OUT"""
import json
import sys


def remaining(pool, used, min_copies=3):
    """The pool's families with >= min_copies copies whose names are not in used, by name."""
    return sorted((f for f in pool if len(f["copies"]) >= min_copies and f["family"] not in used), key=lambda f: f["family"])


def main(pool_p, r15_p, out_p):
    pool = json.load(open(pool_p))
    used = {f["family"] for f in json.load(open(r15_p))["sets"]["r15"]}
    fams = remaining(pool["sets"]["auto15"], used)
    out = {k: v for k, v in pool.items() if k != "sets"}
    out["sets"] = {"r16": fams}
    json.dump(out, open(out_p, "w"), indent=1)
    print(f"{len(fams)} families, {sum(len(f['copies']) for f in fams)} copies: "
          + ", ".join(f"{f['family']} ({len(f['copies'])})" for f in fams))
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))

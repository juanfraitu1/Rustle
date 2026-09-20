#!/usr/bin/env python3
"""o2_emit_materialize.py <fam_id> -- emit a FRESH Python materialize_family for the Rust
e2e parity test, WITHOUT clobbering the shipped cache.

Prints a single JSON object to stdout:
  {"post":        <materialize_family(fam, use_cache=False) return, POST recombinant-abstain>,
   "pre_fresh":   <the pre-abstain VG this fresh run dumped to its private cache>,
   "cached_stale":<the shipped /home/juanfra/winloci_scratch/o2vg/fam<F>.json, or null>}

`post` is what the Rust `materialize_family(...).post_abstain` must equal (both re-run the
SAME current minimap2 -> same alignments -> same columns). `pre_fresh` vs `cached_stale`
discloses any minimap2 drift since the shipped cache was built. The fresh run writes to a
PRIVATE temp cache so the shipped cache is left untouched.

Run: PYTHONHASHSEED=0 /home/juanfra/miniforge3/bin/python bench/o2_emit_materialize.py 39
"""
import json
import os
import sys
import tempfile

os.environ.setdefault("PYTHONHASHSEED", "0")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import o2_vg_visualization as o2   # noqa: E402


def main():
    fam = int(sys.argv[1])
    real_cache = o2.CACHE
    stale_path = os.path.join(real_cache, f"fam{fam}.json")
    stale = json.load(open(stale_path)) if os.path.exists(stale_path) else None
    with tempfile.TemporaryDirectory() as td:
        o2.CACHE = td
        post = o2.materialize_family(fam, use_cache=False)
        pre_fresh = json.load(open(os.path.join(td, f"fam{fam}.json")))
        o2.CACHE = real_cache
    json.dump(dict(post=post, pre_fresh=pre_fresh, cached_stale=stale), sys.stdout)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Diff two parity_decisions JSONL logs (rustle + stringtie).

Usage:
    diff.py rustle.jsonl stringtie.jsonl [--step junction_accept] [--range 22149137-22155355]

Joins events by (step, start, end, strand). Reports:
- Events present in both, payloads identical
- Events present in both, payloads differ → diff fields
- Events only in rustle
- Events only in stringtie
"""
import json, sys, argparse, collections
from typing import Dict, List, Tuple, Any

EventKey = Tuple[str, int, int, str]   # (step, start, end, strand)


def load(path: str, step_filter: str = None, range_filter: Tuple[int, int] = None):
    """Returns dict keyed by EventKey → list of payloads (multi-emit allowed)."""
    by_key: Dict[EventKey, List[dict]] = collections.defaultdict(list)
    n_total = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                ev = json.loads(line)
            except json.JSONDecodeError:
                continue
            step = ev.get("step", "")
            if step.startswith("_"):
                continue
            if step_filter and step != step_filter:
                continue
            s, e = int(ev.get("start", 0)), int(ev.get("end", 0))
            if range_filter:
                rs, re_ = range_filter
                if not (s == 0 and e == 0):
                    if e < rs or s > re_:
                        continue
            strand = ev.get("strand", "?")
            key = (step, s, e, strand)
            by_key[key].append(ev.get("payload", {}))
            n_total += 1
    return by_key, n_total


def diff_payloads(a: dict, b: dict) -> List[str]:
    """Return list of "key: a_val ≠ b_val" strings for differing fields."""
    out = []
    keys = sorted(set(a.keys()) | set(b.keys()))
    for k in keys:
        av, bv = a.get(k), b.get(k)
        if av != bv:
            out.append(f"{k}: rustle={av} stringtie={bv}")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rustle_log")
    ap.add_argument("stringtie_log")
    ap.add_argument("--step", default=None, help="filter to one step name")
    ap.add_argument("--range", default=None, help="filter to coord range LO-HI")
    ap.add_argument("--show-equal", action="store_true", help="also print events where payloads agree")
    ap.add_argument("--max", type=int, default=20, help="max events of each kind to print")
    args = ap.parse_args()

    rng = None
    if args.range:
        lo, hi = args.range.split("-")
        rng = (int(lo), int(hi))

    rustle, n_r = load(args.rustle_log, args.step, rng)
    string, n_s = load(args.stringtie_log, args.step, rng)
    print(f"loaded rustle:    {len(rustle):5d} unique keys ({n_r} events)")
    print(f"loaded stringtie: {len(string):5d} unique keys ({n_s} events)")

    common = set(rustle) & set(string)
    only_r = set(rustle) - set(string)
    only_s = set(string) - set(rustle)
    print(f"\ncommon keys: {len(common)}")
    print(f"only rustle: {len(only_r)}")
    print(f"only stringtie: {len(only_s)}")

    # Within common, count payload divergences
    eq, ne = 0, 0
    diffs = []
    for k in sorted(common):
        rps = rustle[k]
        sps = string[k]
        # Compare first payload of each (events emitted multiple times share key)
        rp, sp = rps[0], sps[0]
        if rp == sp:
            eq += 1
        else:
            ne += 1
            diffs.append((k, rp, sp))

    print(f"common w/ equal payload:    {eq}")
    print(f"common w/ DIFFERING payload: {ne}")

    if args.show_equal and eq > 0:
        print(f"\n--- {min(args.max, eq)} equal-payload events ---")
        n = 0
        for k in sorted(common):
            if rustle[k][0] == string[k][0]:
                print(f"  {k} :: {rustle[k][0]}")
                n += 1
                if n >= args.max:
                    break

    if ne > 0:
        print(f"\n--- {min(args.max, ne)} DIFFERING-payload events ---")
        for (k, rp, sp) in diffs[:args.max]:
            print(f"  {k}")
            for d in diff_payloads(rp, sp):
                print(f"    {d}")

    if only_r:
        print(f"\n--- {min(args.max, len(only_r))} ONLY in rustle ---")
        for k in sorted(only_r)[:args.max]:
            print(f"  {k} :: {rustle[k][0]}")

    if only_s:
        print(f"\n--- {min(args.max, len(only_s))} ONLY in stringtie ---")
        for k in sorted(only_s)[:args.max]:
            print(f"  {k} :: {string[k][0]}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Three-layer junction-set parity: raw / accepted / good_junction.
Answers: do we have exactly ST's junctions? only-RU / only-ST / shared-same / shared-divergent.
Usage: python3 bench/junction_set_diff.py /tmp/ru_parity.jsonl /tmp/st_parity.jsonl"""
import json, sys
from collections import defaultdict

def load(path):
    layers = {"junction_raw": {}, "junction_accept": {}, "good_junction": set()}
    with open(path) as fh:
        for line in fh:
            try:
                r = json.loads(line)
            except ValueError:
                continue
            step = r.get("step")
            if step not in layers:
                continue
            key = (r["start"], r["end"])
            if step == "good_junction":
                layers[step].add(key)
            else:
                layers[step][key] = r.get("payload", {})
    return layers

def report(name, ru_set, st_set):
    only_ru = ru_set - st_set
    only_st = st_set - ru_set
    shared = ru_set & st_set
    print(f"  [{name}] shared={len(shared)}  only-RU={len(only_ru)}  only-ST={len(only_st)}")
    return only_ru, only_st, shared

def main():
    ru = load(sys.argv[1] if len(sys.argv) > 1 else "/tmp/ru_parity.jsonl")
    st = load(sys.argv[2] if len(sys.argv) > 2 else "/tmp/st_parity.jsonl")

    print("=== Layer 1: RAW observed junctions ===")
    report("raw", set(ru["junction_raw"]), set(st["junction_raw"]))

    print("=== Layer 2: ACCEPTED (junction_accept, accepted=true) ===")
    ru_acc = {k for k, p in ru["junction_accept"].items() if p.get("accepted")}
    st_acc = {k for k, p in st["junction_accept"].items() if p.get("accepted")}
    report("accepted", ru_acc, st_acc)

    print("=== Layer 3: good_junction membership ===")
    only_ru_g, only_st_g, shared_g = report("good", ru["good_junction"], st["good_junction"])

    print("=== Shared-divergent treatment (in raw of both, kept by ST, dropped by RU good) ===")
    raw_both = set(ru["junction_raw"]) & set(st["junction_raw"])
    divergent = [k for k in raw_both
                 if k in st["good_junction"] and k not in ru["good_junction"]]
    # bucket the divergent ones by RU's reason (from junction_accept) / mm sign
    buckets = defaultdict(int)
    for k in divergent:
        p = ru["junction_accept"].get(k)
        if p is None:
            buckets["not_in_ru_accept"] += 1
        else:
            reason = p.get("reason", "?")
            mmneg = "mm<0" if float(p.get("mm", 0)) < 0 else "mm>=0"
            buckets[f"{reason}|{mmneg}"] += 1
    print(f"  ST-keeps-but-RU-drops: {len(divergent)}")
    for b, n in sorted(buckets.items(), key=lambda x: -x[1]):
        print(f"    {b}: {n}")

if __name__ == "__main__":
    main()

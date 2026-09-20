#!/usr/bin/env python3
"""Layer-1 gate: compare junction_accept between Rustle and StringTie.

Reports, of the junctions StringTie ACCEPTS, how many Rustle REJECTS and why.
The Layer-1 (mm_negative) gate is satisfied when the 'mm_negative' bucket -> 0.
Inputs: /tmp/ru_ja.jsonl /tmp/st_ja.jsonl (capture with PARITY_FILTER_STEPS=junction_accept).
"""
import json, collections

def accepted_and_rejected(path):
    acc = set()
    rej = {}  # (start,end) -> reason
    for line in open(path):
        try: e = json.loads(line)
        except Exception: continue
        if e.get("step") != "junction_accept": continue
        p = e.get("payload", e)
        k = (e.get("start"), e.get("end"))
        if p.get("accepted") in (True, 1, "true"):
            acc.add(k)
        else:
            rej[k] = str(p.get("reason", "?"))
    return acc, rej

ru_acc, ru_rej = accepted_and_rejected("/tmp/ru_ja.jsonl")
st_acc, _ = accepted_and_rejected("/tmp/st_ja.jsonl")
print(f"Rustle accepted: {len(ru_acc)}  ST accepted: {len(st_acc)}")

# of ST-accepted junctions, which does Rustle reject, by reason
buckets = collections.Counter()
for k in st_acc:
    if k in ru_acc:
        continue
    buckets[ru_rej.get(k, "absent_from_rustle")] += 1
print("=== ST-accepted junctions NOT accepted by Rustle, by Rustle reason ===")
for r, n in buckets.most_common():
    print(f"  {n}\t{r}")
print(f"\nLAYER-1 GATE (mm_negative bucket): {buckets.get('mm_negative', 0)} (target 0)")

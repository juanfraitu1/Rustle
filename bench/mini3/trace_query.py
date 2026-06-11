#!/usr/bin/env python3
"""Tolerant JSONL parity-trace query. Skips malformed lines.

Usage:
  python3 trace_query.py <jsonl> <step> <lo> <hi> [keys...]
  python3 trace_query.py <jsonl> STEPS            # list step counts in range (lo hi optional)
"""
import json, sys

def loadlines(path):
    out = []
    for i, l in enumerate(open(path)):
        l = l.strip()
        if not l:
            continue
        try:
            out.append(json.loads(l))
        except Exception:
            # tolerant: try to salvage step/start/end via regex
            import re
            m_step = re.search(r'"step":"([^"]+)"', l)
            m_s = re.search(r'"start":(-?\d+)', l)
            m_e = re.search(r'"end":(-?\d+)', l)
            out.append({"step": m_step.group(1) if m_step else "?BAD?",
                        "start": int(m_s.group(1)) if m_s else 0,
                        "end": int(m_e.group(1)) if m_e else 0,
                        "_raw": l, "_malformed": True})
    return out

path = sys.argv[1]
step = sys.argv[2]
events = loadlines(path)

if step == "STEPS":
    lo = int(sys.argv[3]) if len(sys.argv) > 3 else None
    hi = int(sys.argv[4]) if len(sys.argv) > 4 else None
    import collections
    c = collections.Counter()
    for o in events:
        s = o.get("start", 0)
        if lo is None or (lo <= s <= hi):
            c[o.get("step", "?")] += 1
    for k, v in c.most_common():
        print(f"  {k}: {v}")
    sys.exit(0)

lo = int(sys.argv[3])
hi = int(sys.argv[4])
keys = sys.argv[5:]
n = 0
for o in events:
    if o.get("step") != step:
        continue
    s = o.get("start", 0)
    if not (lo <= s <= hi):
        continue
    n += 1
    p = o.get("payload", {})
    if keys:
        kv = " ".join(f"{k}={p.get(k)}" for k in keys)
    else:
        kv = json.dumps(p)
    bad = " [MALFORMED]" if o.get("_malformed") else ""
    print(f"  {s}-{o.get('end')} {kv}{bad}")
print(f"  ({n} events)")

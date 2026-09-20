#!/usr/bin/env python3
"""Phase-0b: for extras whose chain ST ALSO extracted, find which ST filter killed it.

Link ST pred_kill (coords+reason) to ST path_extracted (coords+introns) by exact
(start,end,strand), giving chain -> kill_reason. Then for each Rustle extra whose
chain is in ST's extracted set, report the ST kill reason."""
import json, re
from collections import defaultdict, Counter

def parse_introns_field(s):
    if not s: return tuple()
    return tuple((int(a), int(b)) for a, b in (t.split("-") for t in s.split(",")))

# ST path_extracted: (start,end,strand) -> chain ; chain -> set
st_coord_chain = {}
st_chain_set = defaultdict(set)
for line in open("/tmp/st_pe.jsonl"):
    try: e = json.loads(line)
    except: continue
    if e.get("step") != "path_extracted": continue
    p = e.get("payload", e)
    ich = parse_introns_field(p.get("introns", ""))
    if not ich: continue
    st_coord_chain[(e["start"], e["end"], e["strand"])] = ich
    st_chain_set[e["strand"]].add(ich)

# ST pred_kill: (start,end,strand) -> reason
st_kill = {}
for line in open("/tmp/st_pk.jsonl"):
    try: e = json.loads(line)
    except: continue
    if e.get("step") != "pred_kill": continue
    p = e.get("payload", e)
    st_kill[(e["start"], e["end"], e["strand"])] = str(p.get("reason", ""))

# chain -> kill reason (via coord link)
chain_reason = defaultdict(set)
linked = 0
for (s, en, strand), ich in st_coord_chain.items():
    if (s, en, strand) in st_kill:
        chain_reason[(strand, ich)].add(st_kill[(s, en, strand)])
        linked += 1
print(f"ST path_extracted coords linked to a pred_kill: {linked}")

# Rustle GTF extras (reuse phase0 logic, inline)
def introns_from_exons(exons):
    ex = sorted(exons)
    return tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
gtf = {}; exons_by_tid = defaultdict(list)
for line in open("/tmp/full_default.gtf"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9: continue
    m = re.search(r'transcript_id "([^"]+)"', f[8])
    if not m: continue
    tid = m.group(1)
    if f[2] == "transcript":
        src = re.search(r'source "([^"]+)"', f[8])
        gtf[tid] = {"strand": f[6], "source": src.group(1) if src else "?"}
    elif f[2] == "exon":
        exons_by_tid[tid].append((int(f[3]), int(f[4])))
for tid, ex in exons_by_tid.items():
    if tid in gtf:
        gtf[tid]["introns"] = introns_from_exons(ex)
cls = {}
with open("/tmp/cmp_default.full_default.gtf.tmap") as fh:
    next(fh, None)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5: cls[f[4]] = f[2]

# For each non-'=' multi-exon extra whose chain ST extracted: report ST kill reason
report = Counter()
no_kill = 0
for tid, info in gtf.items():
    if cls.get(tid, "?") == "=": continue
    ich = info.get("introns", tuple())
    if not ich: continue
    strand = info["strand"]
    if ich not in st_chain_set[strand]:
        continue   # not co-extracted (the path-enum bucket, handled elsewhere)
    reasons = chain_reason.get((strand, ich))
    if reasons:
        for r in reasons:
            report[r] += 1
    else:
        no_kill += 1
print("\n=== ST kill reason for Rustle extras that ST ALSO extracted ===")
for k, v in report.most_common():
    print(f"  {v:4d}  {k}")
print(f"  {no_kill:4d}  (co-extracted but no coord-linked ST kill found — killed at a stage w/o coord match, or survived in ST as a variant)")

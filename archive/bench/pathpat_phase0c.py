#!/usr/bin/env python3
"""Phase-0c: attribute each Rustle extra (co-extracted or not) to its disposition.

Direct match: Rustle extra (chain) -> is it in ST's extracted set?
  - If yes: find ST pred_kill near its coords (tolerance) -> ST kill reason.
  - If no: path-enum bucket (ST never extracted it); classify vs nearest ST chain.
"""
import json, re
from collections import defaultdict, Counter

TOL = 60  # bp tolerance on start/end for matching extra <-> ST kill

def parse_introns(s):
    if not s: return tuple()
    return tuple((int(a), int(b)) for a, b in (t.split("-") for t in s.split(",")))

# ST extracted chains by strand
st_chain_set = defaultdict(set)
for line in open("/tmp/st_pe.jsonl"):
    try: e = json.loads(line)
    except: continue
    if e.get("step") != "path_extracted": continue
    ich = parse_introns(e.get("payload", e).get("introns", ""))
    if ich: st_chain_set[e["strand"]].add(ich)

# ST pred_kill list by strand: [(start,end,reason)]
st_kill = defaultdict(list)
for line in open("/tmp/st_pk.jsonl"):
    try: e = json.loads(line)
    except: continue
    if e.get("step") != "pred_kill": continue
    st_kill[e["strand"]].append((e["start"], e["end"], str(e.get("payload", e).get("reason", ""))))

# Rustle GTF
def introns_from_exons(ex):
    ex = sorted(ex)
    return tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))
gtf = {}; exons_by = defaultdict(list)
for line in open("/tmp/full_default.gtf"):
    f = line.rstrip("\n").split("\t")
    if len(f) < 9: continue
    m = re.search(r'transcript_id "([^"]+)"', f[8])
    if not m: continue
    tid = m.group(1)
    if f[2] == "transcript":
        src = re.search(r'source "([^"]+)"', f[8])
        gtf[tid] = {"strand": f[6], "start": int(f[3]), "end": int(f[4]),
                    "source": src.group(1) if src else "?"}
    elif f[2] == "exon":
        exons_by[tid].append((int(f[3]), int(f[4])))
for tid, ex in exons_by.items():
    if tid in gtf: gtf[tid]["introns"] = introns_from_exons(ex)
cls = {}
with open("/tmp/cmp_default.full_default.gtf.tmap") as fh:
    next(fh, None)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5: cls[f[4]] = f[2]

def nearest_st(ich, strand):
    eset = set(ich); best = None; bs = -1
    for st in st_chain_set[strand]:
        sh = len(eset & set(st))
        if sh > bs: bs = sh; best = st
    if best is None: return "no_st_overlap"
    eo = eset - set(best); so = set(best) - eset
    if not eo and not so: return "exact"
    if eo and not so: return "superset_extra_junction"
    if so and not eo: return "subset_contained"
    return "alt_splice_diverge"

def find_kill(tid, info):
    st = info["strand"]; s = info["start"]; e = info["end"]
    for (ks, ke, r) in st_kill[st]:
        if abs(ks - s) <= TOL and abs(ke - e) <= TOL:
            return r
    return None

coextracted = Counter()
pathenum = Counter()
n_co = n_pe = 0
for tid, info in gtf.items():
    if cls.get(tid, "?") == "=": continue
    ich = info.get("introns", tuple())
    if not ich: continue
    strand = info["strand"]
    if ich in st_chain_set[strand]:
        n_co += 1
        r = find_kill(tid, info)
        coextracted[r or "no_coord_match"] += 1
    else:
        n_pe += 1
        pathenum[nearest_st(ich, strand)] += 1

print(f"=== 186 multi-exon extras split ===")
print(f"  co-extracted by ST (FILTER divergence): {n_co}")
print(f"  never extracted by ST (PATH-ENUM divergence): {n_pe}")
print(f"\n=== co-extracted: ST kill reason (tol={TOL}bp) ===")
for k, v in coextracted.most_common(): print(f"  {v:4d}  {k}")
print(f"\n=== path-enum: divergence type vs nearest ST chain ===")
for k, v in pathenum.most_common(): print(f"  {v:4d}  {k}")

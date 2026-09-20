#!/usr/bin/env python3
"""Bucket B characterization: alt-splice junction-selection extras.

For each Rustle non-'=' multi-exon extra, find the nearest REFERENCE (GGO_19.gtf =
ST output) intron chain. If they differ by exactly one intron pair sharing one
endpoint (a shifted donor or acceptor), characterize the shift (end, magnitude,
direction) and report whether the corresponding ref is MISSED (→ 2-for-1 fix)."""
import re
from collections import defaultdict, Counter

RU_GTF = "/tmp/full_default.gtf"
TMAP   = "/tmp/cmp_default.full_default.gtf.tmap"
REF    = "../GGO_19.gtf"

def introns(exons):
    ex = sorted(exons)
    return tuple((ex[i][1] + 1, ex[i + 1][0] - 1) for i in range(len(ex) - 1))

def load_gtf(path):
    tx = {}; exb = defaultdict(list)
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9: continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m: continue
        t = m.group(1)
        if f[2] == "transcript":
            tx[t] = {"chrom": f[0], "strand": f[6]}
        elif f[2] == "exon":
            exb[t].append((int(f[3]), int(f[4])))
    for t, ex in exb.items():
        if t in tx:
            tx[t]["ich"] = introns(ex)
    return tx

ru = load_gtf(RU_GTF)
ref = load_gtf(REF)

# ref chains by strand + which refs are matched (have '=' qry) — from refmap
ref_by_strand = defaultdict(list)
for t, info in ref.items():
    if info.get("ich"):
        ref_by_strand[info["strand"]].append((t, info["ich"]))
matched_ref = set()
import glob
rm = glob.glob("/tmp/cmp_default.*.refmap")
if rm:
    with open(rm[0]) as fh:
        next(fh, None)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2] == "=":
                matched_ref.add(f[1])

cls = {}
with open(TMAP) as fh:
    next(fh, None)
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5: cls[f[4]] = f[2]

def nearest_ref(ich, strand):
    es = set(ich); best = None; bid = None; bs = -1
    for rid, rich in ref_by_strand[strand]:
        sh = len(es & set(rich))
        if sh > bs: bs = sh; best = rich; bid = rid
    return bid, best, bs

shifts = []
n_2for1 = 0
shift_dirs = Counter()
shift_mags = []
for t, info in ru.items():
    if cls.get(t, "?") == "=": continue
    ich = info.get("ich")
    if not ich: continue
    strand = info["strand"]
    rid, rich, shared = nearest_ref(ich, strand)
    if rich is None: continue
    eo = sorted(set(ich) - set(rich))   # Rustle-unique introns
    so = sorted(set(rich) - set(ich))   # ref-unique introns
    # alt-splice signature: exactly one unique intron each, sharing ONE endpoint
    if len(eo) == 1 and len(so) == 1:
        (rd, ra) = eo[0]; (sd, sa) = so[0]
        if rd == sd and ra != sa:
            end = "acceptor"; shift = ra - sa
        elif ra == sa and rd != sd:
            end = "donor"; shift = rd - sd
        else:
            continue
        miss = rid not in matched_ref
        if miss: n_2for1 += 1
        shift_dirs["ru>ref (shorter exon/longer intron)" if shift > 0 else "ru<ref"] += 1
        shift_mags.append(abs(shift))
        shifts.append((t, cls.get(t), end, shift, "MISS_ref" if miss else "ref_also_emitted", rid))

print(f"=== single-junction-shift extras (alt-donor/acceptor): {len(shifts)} ===")
print(f"  of which the matching ref is MISSED (true 2-for-1): {n_2for1}")
print(f"\n=== shift end ===")
for k, v in Counter(s[2] for s in shifts).most_common(): print(f"  {v:3d}  {k}")
print(f"\n=== shift direction ===")
for k, v in shift_dirs.most_common(): print(f"  {v:3d}  {k}")
import statistics as st
if shift_mags:
    print(f"\n=== shift magnitude: median={st.median(shift_mags)}bp  <=10bp={sum(1 for m in shift_mags if m<=10)}  <=30bp={sum(1 for m in shift_mags if m<=30)}  >30bp={sum(1 for m in shift_mags if m>30)}")
print("\n=== sample (first 15) ===")
for t, c, end, shift, miss, rid in shifts[:15]:
    print(f"  {t} [{c}] {end} shift={shift:+d}bp {miss} ref={rid}")

#!/usr/bin/env python3
"""Leave-one-out core projection (PREREG_core_projection_2026-09-07, md5 ee7ceb51).
For each member of a family, drop it, recompute the survivors' cores, project them into the genome,
and ask whether the dropped member's locus comes back under the UNCHANGED membership rule.
usage: o1_core_projection.py <units.tsv> <family_id> <truth.bed> <sedef.bed> <genome.mmi> <genome.fa> <workdir>
"""
import sys, csv, os, subprocess, json, collections
import numpy as np

units_p, FAM, truth_p, sedef_p, mmi, gfa, W = sys.argv[1:8]
os.makedirs(W, exist_ok=True)
MIN_ID, MIN_COV, MIN_BP = 0.70, 0.30, 300

truth = []
for l in open(truth_p):
    f = l.rstrip('\n').split('\t')
    if len(f) >= 3 and not f[0].startswith('#'):
        truth.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else f'{f[0]}:{f[1]}'))
rows = [r for r in csv.DictReader(open(units_p), delimiter='\t')
        if r['family_id'] == FAM and r['member_status'] != 'dropped']
mem = [(r['chrom'], int(r['locus_start']), int(r['locus_end']), r['core_hull'], r['copy_idx']) for r in rows]
print(f"family {FAM}: {len(mem)} members", flush=True)

contigs = {m[0] for m in mem}
sed = []
for l in open(sedef_p):
    f = l.split('\t')
    if not f[1].isdigit(): continue
    if f[0] in contigs or f[3] in contigs:
        sed.append((f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])))
print(f"sedef pairs on substrate: {len(sed)}", flush=True)

def ov(a, b, c, d): return max(0, min(b, d) - max(a, c))

def core_profile(c, s, e, others):
    """core segments of interval [s,e) on contig c against a list of (chrom,start,end) partners"""
    half = len(others) / 2.0
    ev = []
    for a in sed:
        for (qc, qs, qe, tc, ts, te) in ((a[0], a[1], a[2], a[3], a[4], a[5]),
                                         (a[3], a[4], a[5], a[0], a[1], a[2])):
            if qc != c or qe <= s or qs >= e: continue
            for oi, o in enumerate(others):
                if o[0] != tc: continue
                if te > o[1] and ts < o[2]:
                    ev.append((max(s, qs), 1, oi)); ev.append((min(e, qe), -1, oi))
    ev.sort()
    cnt = collections.Counter(); depth = 0; md = 0; last = None; segs = []
    for pos, d, oi in ev:
        if last is not None and pos > last and depth >= half:
            if segs and segs[-1][1] == last: segs[-1][1] = pos
            else: segs.append([last, pos])
        if d > 0:
            if cnt[oi] == 0: depth += 1
            cnt[oi] += 1
        else:
            cnt[oi] -= 1
            if cnt[oi] == 0: depth -= 1
        md = max(md, depth); last = pos
    return md, sum(b - a for a, b in segs), segs

def run(cmd, **kw): return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True, **kw)

folds = []
for held in range(len(mem)):
    hm = mem[held]
    S = [m for i, m in enumerate(mem) if i != held]
    Sspan = [(m[0], m[1], m[2]) for m in S]
    # 1. survivors' cores, recomputed over survivors only
    cores = []
    for i, m in enumerate(S):
        others = [x for j, x in enumerate(Sspan) if j != i]
        md, cb, segs = core_profile(m[0], m[1], m[2], others)
        if segs: cores.append((m[0], segs[0][0], segs[-1][1], m[4]))
    med_core = float(np.median([c[2] - c[1] for c in cores])) if cores else 0.0
    fa = f"{W}/f{held}.cores.fa"
    with open(f"{W}/f{held}.regions", 'w') as o:
        for c in cores: o.write(f"{c[0]}:{c[1]+1}-{c[2]}\n")
    run(f"samtools faidx {gfa} -r {W}/f{held}.regions > {fa}")
    # 2. project
    paf = f"{W}/f{held}.paf"
    run(f"minimap2 -x asm20 -c -N 50 -p 0.1 -t 4 {mmi} {fa} > {paf} 2> {W}/f{held}.mm2.log")
    hits = []
    for l in open(paf):
        f = l.rstrip('\n').split('\t')
        qlen, ts, te, tc = int(f[1]), int(f[7]), int(f[8]), f[5]
        blk = int(f[10]); nm = next((int(x[5:]) for x in f[12:] if x.startswith('NM:i:')), blk - int(f[9]))
        ident = 1 - nm / blk if blk else 0
        cov = blk / max(qlen, te - ts)
        if blk >= MIN_BP and ident >= MIN_ID and cov >= MIN_COV: hits.append((tc, ts, te))
    # 3. merge, 4. discard intervals overlapping a survivor
    merged = []
    for c in sorted({h[0] for h in hits}):
        iv = sorted((h[1], h[2]) for h in hits if h[0] == c)
        cur = None
        for a, b in iv:
            if cur and a <= cur[1]: cur[1] = max(cur[1], b)
            else:
                if cur: merged.append((c, cur[0], cur[1]))
                cur = [a, b]
        if cur: merged.append((c, cur[0], cur[1]))
    n_hits, n_merged = len(hits), len(merged)
    cand = [m for m in merged if not any(s[0] == m[0] and ov(m[1], m[2], s[1], s[2]) > 0 for s in Sspan)]
    # 5. admit by the unchanged membership rule
    admitted = []
    for c, s, e in cand:
        if e - s < MIN_BP: continue
        md, cb, segs = core_profile(c, s, e, Sspan)
        if cb == 0: continue
        if cb >= (e - s) / 2 or cb >= med_core / 2:
            admitted.append((c, segs[0][0], segs[-1][1], cb, md, s, e))
    # 6. score
    tm = next((t for t in truth if t[0] == hm[0] and ov(hm[1], hm[2], t[1], t[2]) > 0), None)
    rec = None
    if tm:
        for a in admitted:
            o = ov(a[1], a[2], tm[1], tm[2])
            if o > 0 and (o >= 0.5 * (a[2] - a[1]) or o >= 0.5 * (tm[2] - tm[1])):
                if rec is None or o > rec[0]: rec = (o, a)
    fal = [a for a in admitted if not any(a[0] == t[0] and ov(a[1], a[2], t[1], t[2]) > 0 for t in truth)]
    folds.append(dict(held=hm[4], held_locus=f"{hm[0]}:{hm[1]}-{hm[2]}", truth=tm[3] if tm else None,
                      n_hits=n_hits, n_merged=n_merged, n_cand=len(cand), n_admitted=len(admitted), n_false=len(fal),
                      recovered=rec is not None,
                      rec_interval=f"{rec[1][0]}:{rec[1][1]}-{rec[1][2]}" if rec else None,
                      rec_cov=round(rec[0] / (tm[2] - tm[1]), 3) if rec else None,
                      rec_ratio=round((rec[1][2] - rec[1][1]) / (tm[2] - tm[1]), 3) if rec else None,
                      false_intervals=[f"{a[0]}:{a[1]}-{a[2]}" for a in fal]))
    print(f"  fold {held:2d} copy {hm[4]:>3} {folds[-1]['truth']}: hits {n_hits:4d} merged {n_merged:3d} cand {len(cand):3d} admitted {len(admitted):3d} "
          f"false {len(fal):3d} recovered {folds[-1]['recovered']} cov {folds[-1]['rec_cov']} ratio {folds[-1]['rec_ratio']}", flush=True)
json.dump(folds, open(f"{W}/folds.json", 'w'), indent=1)
R = sum(1 for f in folds if f['recovered'])
covs = [f['rec_cov'] for f in folds if f['recovered']]
rats = [f['rec_ratio'] for f in folds if f['recovered']]
print(f"\nR recovery {R}/{len(folds)}")
print(f"B median truth coverage {np.median(covs):.3f} | in-band 0.5-2x {sum(1 for r in rats if 0.5<=r<=2)}/{len(rats)}")
print(f"F false admissions per fold: median {np.median([f['n_false'] for f in folds]):.1f}, max {max(f['n_false'] for f in folds)}")
print(f"K hits {np.median([f['n_hits'] for f in folds]):.1f} -> merged {np.median([f['n_merged'] for f in folds]):.1f} -> candidates {np.median([f['n_cand'] for f in folds]):.1f} -> admitted {np.median([f['n_admitted'] for f in folds]):.1f} (median)")
fi = collections.Counter(x for f in folds for x in f['false_intervals'])
print("H most frequent false admissions:", fi.most_common(8))

#!/usr/bin/env python3
"""The family as a CLOSURE (PREREG_closure_definition_2026-09-07, md5 8f144fd6).

A family is a set of genomic intervals that all carry copies of the same segment, where the segment is what
at least half of them share. Iterated to a fixed point:
  seed (connected groups of the alignment links) -> core (shared with half the others) -> extend (search the
  genome for the core) -> repeat.
No MCL: the seed is connected components, the annotation only proposes intervals.

usage: o1_closure.py <paf> <sedef.bed> <genome.mmi> <genome.fa> <truth.bed> <workdir> [--seed-mode full|half|neighbours] [--max-iter N]
"""
import sys, os, subprocess, collections, json, statistics
import numpy as np

paf_p, sed_p, mmi, gfa, truth_p, W = sys.argv[1:7]
args = sys.argv[7:]
seed_mode = args[args.index('--seed-mode') + 1] if '--seed-mode' in args else 'full'
seed_bed = args[args.index('--seed-bed') + 1] if '--seed-bed' in args else None
MAX_ITER = int(args[args.index('--max-iter') + 1]) if '--max-iter' in args else 10
os.makedirs(W, exist_ok=True)
MIN_ID, MIN_COV, MIN_BP = 0.70, 0.30, 300
rng = np.random.default_rng(1337)

truth = []
for l in open(truth_p):
    f = l.rstrip('\n').split('\t')
    if len(f) >= 3 and not f[0].startswith('#'):
        truth.append((f[0], int(f[1]), int(f[2]), f[3] if len(f) > 3 else f'{f[0]}:{f[1]}'))

def parse(k):
    c, rest = k.rsplit(':', 1); s, e = rest.split('-'); return (c, int(s), int(e))   # GFF 1-based start
def ov(a, b, c, d): return max(0, min(b, d) - max(a, c))
def tname(m):
    for tc, ts, te, nm in truth:
        if m[0] == tc and ov(m[1] - 1, m[2], ts, te) > 0: return nm
    return None

# ---- step 1: SEED = connected groups of the alignment links
adj = collections.defaultdict(set)
for l in open(paf_p):
    f = l.rstrip('\n').split('\t'); q, t = f[0], f[5]
    if q == t: continue
    blk = int(f[10]); nm = next((int(x[5:]) for x in f[12:] if x.startswith('NM:i:')), None)
    if nm is None: continue
    if blk >= MIN_BP and 1 - nm / blk >= MIN_ID and blk / max(int(f[1]), int(f[6])) >= MIN_COV:
        adj[q].add(t); adj[t].add(q)
seen, comps = set(), []
for n in list(adj):
    if n in seen: continue
    st, c = [n], set()
    while st:
        x = st.pop()
        if x in seen: continue
        seen.add(x); c.add(x); st.extend(adj[x] - seen)
    comps.append(c)
# the group holding the truth
group = max(comps, key=lambda c: sum(1 for t in truth if any(tname(parse(k)) == t[3] for k in c)))
group = sorted(parse(k) for k in group)
if seed_mode == 'half':
    idx = sorted(rng.choice(len(group), size=max(2, len(group) // 2), replace=False)); seed = [group[i] for i in idx]
elif seed_bed:
    seed = sorted((f[0], int(f[1]) + 1, int(f[2])) for f in (l.split('\t') for l in open(seed_bed).read().strip().split('\n')) if len(f) >= 3)
    seed_mode = 'bed'
elif seed_mode == 'neighbours':
    anchor = next(k for k in sorted(adj) if tname(parse(k)))
    seed = sorted({parse(anchor)} | {parse(x) for x in adj[anchor]})
else:
    seed = list(group)
print(f"seed ({seed_mode}): {len(seed)} intervals from a connected group of {len(group)}", flush=True)

# ---- duplication pairs
contigs = {m[0] for m in group}
sed = []
for l in open(sed_p):
    f = l.split('\t')
    if len(f) < 6 or not f[1].isdigit(): continue
    if f[0] in contigs or f[3] in contigs: sed.append((f[0], int(f[1]), int(f[2]), f[3], int(f[4]), int(f[5])))

def core_of(m, others):
    """positions of interval m linked by duplication pairs to >= half of `others`; returns (max_depth, bp, segs)"""
    half = len(others) / 2.0
    c, rs, re = m[0], m[1] - 1, m[2]; ev = []
    for a in sed:
        for (qc, qs, qe, tc, ts, te) in ((a[0], a[1], a[2], a[3], a[4], a[5]), (a[3], a[4], a[5], a[0], a[1], a[2])):
            if qc != c or qe <= rs or qs >= re: continue
            for oi, o in enumerate(others):
                if o[0] != tc: continue
                if te > o[1] - 1 and ts < o[2]: ev.append((max(rs, qs), 1, oi)); ev.append((min(re, qe), -1, oi))
    ev.sort(); cnt = collections.Counter(); depth = md = 0; last = None; segs = []
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

def classify(cands, members):
    """step 2: which candidates are members, judged against `members`"""
    prof = {i: core_of(m, [o for o in members if o != m]) for i, m in enumerate(cands)}
    cores = sorted(p[1] for i, p in prof.items() if cands[i] in members) or [0]
    medc = cores[len(cores) // 2]
    out = {}
    for i, m in enumerate(cands):
        md, cb, segs = prof[i]; span = m[2] - m[1] + 1
        # "shares a core with at least half of the OTHERS" is vacuously true when there are no others: a
        # one-member set must be allowed to extend rather than be pruned to nothing (AMENDMENT 1). Without
        # this, every single-locus seed collapses to the empty set, which is trivially closed and useless.
        if len([o for o in members if o != m]) == 0:
            out[i] = ('member', md, cb, segs if segs else [[m[1] - 1, m[2]]])
            continue
        out[i] = ('member' if (cb >= span / 2 or (cb > 0 and cb >= medc / 2)) else 'candidate', md, cb, segs)
    return out, medc

def run(cmd): return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)

def extend(members, prof, it):
    """step 3: align the members' core hulls to the genome; return intervals not already members"""
    hulls = []
    for i, m in enumerate(members):
        segs = prof[i][3]
        if segs: hulls.append((m[0], segs[0][0], segs[-1][1]))
    if not hulls: return []
    with open(f'{W}/it{it}.regions', 'w') as o:
        for c, s, e in hulls: o.write(f"{c}:{s+1}-{e}\n")
    run(f"samtools faidx {gfa} -r {W}/it{it}.regions > {W}/it{it}.fa")
    run(f"minimap2 -x asm20 -c -N 50 -p 0.1 -t 4 {mmi} {W}/it{it}.fa > {W}/it{it}.paf 2> {W}/it{it}.mm2.log")
    hits = []
    for l in open(f'{W}/it{it}.paf'):
        f = l.rstrip('\n').split('\t')
        qlen, ts, te, tc = int(f[1]), int(f[7]), int(f[8]), f[5]
        blk = int(f[10]); nm = next((int(x[5:]) for x in f[12:] if x.startswith('NM:i:')), blk - int(f[9]))
        if blk >= MIN_BP and (1 - nm / blk) >= MIN_ID and blk / max(qlen, te - ts) >= MIN_COV:
            hits.append((tc, ts, te))
    merged = []
    for c in sorted({h[0] for h in hits}):
        cur = None
        for a, b in sorted((h[1], h[2]) for h in hits if h[0] == c):
            if cur and a <= cur[1]: cur[1] = max(cur[1], b)
            else:
                if cur: merged.append((c, cur[0], cur[1]))
                cur = [a, b]
        if cur: merged.append((c, cur[0], cur[1]))
    return [(c, s + 1, e) for c, s, e in merged
            if not any(m[0] == c and ov(m[1] - 1, m[2], s, e) > 0 for m in members) and e - s >= MIN_BP]

# ---- the closure
members = list(seed); history = []; trace = []
for it in range(MAX_ITER):
    st, medc = classify(members, members)
    members2 = [members[i] for i in range(len(members)) if st[i][0] == 'member']
    left = [m for m in members if m not in members2]
    prof = {i: st[j] for i, (j, m) in enumerate((j, m) for j, m in enumerate(members) if st[j][0] == 'member')}
    prof = {}
    st2, _ = classify(members2, members2)
    for i in range(len(members2)): prof[i] = st2[i]
    new = extend(members2, prof, it)
    admitted = []
    if new:
        cand = members2 + new
        st3, _ = classify(cand, members2)
        admitted = [cand[i] for i in range(len(members2), len(cand)) if st3[i][0] == 'member']
    nxt = sorted(set(members2) | set(admitted))
    trace.append(dict(iteration=it, before=len(members), after_core=len(members2), dropped=len(left),
                      proposed=len(new), admitted=len(admitted), after=len(nxt),
                      truth_members=len({tname(m) for m in nxt if tname(m)})))
    print(f"  it{it}: {len(members)} -> core keeps {len(members2)} (-{len(left)}) -> genome proposes {len(new)}, admits {len(admitted)} -> {len(nxt)}"
          f" | truth loci held {trace[-1]['truth_members']}/{len(truth)}", flush=True)
    history.append(set(nxt))
    if nxt == sorted(members):
        print(f"  FIXED POINT at iteration {it}"); members = nxt; break
    members = nxt
res = dict(seed_mode=seed_mode, seed_n=len(seed), iterations=len(trace), trace=trace,
           members=[f"{m[0]}:{m[1]}-{m[2]}" for m in members],
           n_members=len(members),
           truth_recovered=len({tname(m) for m in members if tname(m)}), truth_total=len(truth),
           precision=round(sum(1 for m in members if tname(m)) / max(len(members), 1), 4),
           monotone=all(history[i] >= history[i+1] for i in range(len(history)-1)))
json.dump(res, open(f'{W}/closure_{seed_mode}.json', 'w'), indent=1)
print(f"\nFIXED POINT: {res['n_members']} members | precision {res['precision']:.3f} | "
      f"truth {res['truth_recovered']}/{res['truth_total']} | iterations {res['iterations']} | monotone {res['monotone']}")

#!/usr/bin/env python3
"""Count primary reads supporting each read-through junction in a BAM (PREREG md5 208648a5).

A read supports a junction when its CIGAR carries an intron whose coordinates match within +-TOL. Junctions
are split into the three classes the cross-species test (§6gc) established, so an aligner setting is scored by
what it does to the CONSERVED ones, not by whether it lowers the read-through count.

usage: junction_support.py <bam> <guard.readthrough.tsv> <conserved_ids.txt> <rejected.tsv> [--tol 5]
"""
import sys, subprocess, re, collections
bam, rt_p, cons_p, rej_p = sys.argv[1:5]
TOL = int(sys.argv[sys.argv.index('--tol') + 1]) if '--tol' in sys.argv else 5
conserved = set(open(cons_p).read().split()) if cons_p != '-' else set()
rejected = set()
if rej_p != '-':
    for l in open(rej_p):
        f = l.rstrip('\n').split('\t')
        if len(f) >= 9 and f[5].isdigit(): rejected.add(f"{f[1]}:{f[5]}-{f[6]}")
junc = []
for l in open(rt_p):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or not f[5].isdigit(): continue
    key = f"{f[1]}:{f[5]}-{f[6]}"
    cls = 'conserved' if key in conserved else ('guard_rejected' if key in rejected else 'not_replicating')
    junc.append((f[1], int(f[5]), int(f[6]), key, cls, int(f[8])))
def introns(pos, cig):
    out = []; p = pos
    for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
        n = int(n)
        if op in 'M=XD': p += n
        elif op == 'N': out.append((p, p + n)); p += n
    return out
sup = collections.Counter()
by_ctg = collections.defaultdict(list)
for c, s, e, k, cls, n in junc: by_ctg[c].append((s, e, k))
for c, v in by_ctg.items():
    lo = min(s for s, _, _ in v) - 200000; hi = max(e for _, e, _ in v) + 200000
    out = subprocess.run(f"samtools view -F 2308 {bam} {c}:{max(lo,1)}-{hi}", shell=True, capture_output=True, text=True).stdout
    for line in out.split('\n'):
        f = line.split('\t')
        if len(f) < 6: continue
        for a, b in introns(int(f[3]) - 1, f[5]):
            for s, e, k in v:
                if abs(a - s) <= TOL and abs(b - e) <= TOL: sup[k] += 1
tot = collections.Counter(); lost = collections.Counter(); reads = collections.Counter()
for c, s, e, k, cls, n in junc:
    tot[cls] += 1; reads[cls] += sup[k]
    if sup[k] == 0: lost[cls] += 1
print(f"{'class':18s} {'junctions':>10} {'with ZERO support':>18} {'supporting reads':>17}")
for cls in ('conserved', 'not_replicating', 'guard_rejected'):
    if tot[cls]: print(f"   {cls:15s} {tot[cls]:10d} {lost[cls]:18d} {reads[cls]:17d}")
print("\nper-junction support:")
for c, s, e, k, cls, n in sorted(junc, key=lambda x: -sup[x[3]]):
    print(f"   {k:34s} {cls:16s} {sup[k]:6d}")

#!/usr/bin/env python3
"""Row 1562: flag satellite families by within-copy k-mer duplication. usage: satseq.py COPIES_FA TRUTH_IDS [K] [THR]"""
import collections, statistics, sys
fa, truth_f = sys.argv[1], sys.argv[2]
k = int(sys.argv[3]) if len(sys.argv) > 3 else 15
thr = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
truth = set(open(truth_f).read().split())
def dup(seq):
    c = collections.Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    n = sum(c.values())
    return sum(v for v in c.values() if v >= 2) / n if n else 0.0
fam = collections.defaultdict(list); gc = collections.defaultdict(list)
name = None; buf = []
def flush():
    if name:
        s = ''.join(buf).upper(); fam[name].append(dup(s)); gc[name].append((s.count('G') + s.count('C')) / max(1, len(s)))
for line in open(fa):
    if line.startswith('>'):
        flush(); name = line[1:].split('|')[0]; buf = []
    else:
        buf.append(line.strip())
flush()
flag = {f for f, v in fam.items() if statistics.median(v) >= thr}
tp = len(flag & truth); fp = len(flag - truth); fn = len(truth - flag)
print(f"families {len(fam)} truth {len(truth)} flagged {len(flag)} tp {tp} fp {fp} fn {fn} precision {tp/max(1,len(flag)):.3f} recall {tp/max(1,len(truth)):.3f}")
for lab, ids in (("truth", truth), ("other", set(fam) - truth)):
    d = [statistics.median(fam[f]) for f in ids if f in fam]; g = [statistics.median(gc[f]) for f in ids if f in fam]
    print(f"{lab}: median dup{k} {statistics.median(d):.3f} q10 {sorted(d)[len(d)//10]:.3f} q90 {sorted(d)[9*len(d)//10]:.3f}  median GC {statistics.median(g):.3f}")
print("false positives:", sorted(flag - truth)[:20])
print("false negatives:", sorted(truth - flag)[:20])

#!/usr/bin/env python3
"""PREREG r1504 scoring: CALLS.fa (only the call ids in IDS) aligned to the T2T panel copies (minimap2 -cx asm20). A call HITS
the copy that is its best target by matching bases, at identity >= 0.99. Recall = expressed missing copies hit / expressed missing;
precision = calls hitting any missing copy / calls; calls on a present copy listed (false alarms on copies the genome has).
  score.py LABEL TRUTH.tsv COPIES.fa CALLS.fa IDS.txt"""
import csv
import subprocess
import sys

label, truth, copies, calls, ids = sys.argv[1:6]
T = {r['gene']: r for r in csv.DictReader(open(truth), delimiter='\t')}
want = {l.strip() for l in open(ids) if l.strip()}
miss = {g for g, r in T.items() if r['status'] == 'missing' and r['expressed'] == '1'}
best = {}
paf = subprocess.run(['minimap2', '-cx', 'asm20', '--secondary=no', copies, calls], capture_output=True, text=True, check=True).stdout
for l in paf.splitlines():
    f = l.split('\t')
    q, t, m, blk = f[0], f[5], int(f[9]), int(f[10])
    if q in want and (q not in best or m > best[q][1]):
        best[q] = (t, m, m / blk)
hit = {q: t for q, (t, m, i) in best.items() if i >= 0.99}
on_miss = {q for q, t in hit.items() if t in miss}
on_present = sorted(f"{q}->{T[t]['name']}" for q, t in hit.items() if T.get(t, {}).get('status') == 'present')
found = {hit[q] for q in on_miss}
rec = len(found) / len(miss) if miss else float('nan')
prec = len(on_miss) / len(want) if want else float('nan')
print(f"{label}\tcalls {len(want)}\ttruth missing+expressed {len(miss)}\tfound {len(found)}\trecall {rec:.3f}\tprecision {prec:.3f}\t"
      f"on present copies {len(on_present)}\tno 99% copy {len(want) - len(hit)}")
print(f"{label}\tfound: {', '.join(sorted(T[g]['name'] for g in found)) or '-'}")
print(f"{label}\tmissed: {', '.join(sorted(T[g]['name'] for g in miss - found)) or '-'}")
print(f"{label}\ton present: {', '.join(on_present[:20]) or '-'}")

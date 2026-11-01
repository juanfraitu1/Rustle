#!/usr/bin/env python3
"""PREREG TIE-MS PORT (register 1483) scorer: arms of `copy_phasing run` on one sample, A (default window) against M10 / M0.
  tie_ms_score.py A_PREFIX M10_PREFIX [M0_PREFIX]
Per arm: tied reads (o2.tsv n_piles >= 2), assigned (a copy), ABSTAIN, HOPELESS, certified copies and putative candidates
(families.tsv / candidates.tsv). Bars (M10 vs A): B1 untied = reads tied under A and not under M10 >= 5% of A's tied;
B2 M10 assigned >= 0.95 x A assigned; B3 conflicting assignments < 1% of the reads both assign."""
import sys, csv

def load(prefix):
    o2 = {(r['family'], r['read']): (r['assigned'], int(r['n_piles'])) for r in csv.DictReader(open(f"{prefix}.o2.tsv"), delimiter='\t')}
    fam = list(csv.DictReader(open(f"{prefix}.families.tsv"), delimiter='\t'))
    cert = sum(int(r['certified']) for r in fam)
    try:
        cand = sum(1 for r in csv.DictReader(open(f"{prefix}.candidates.tsv"), delimiter='\t') if r['flagged'] == '1')
    except FileNotFoundError:
        cand = None
    return o2, cert, cand

def counts(o2):
    tied = {k for k, (a, n) in o2.items() if n >= 2}
    assigned = {k: a for k, (a, n) in o2.items() if a not in ('ABSTAIN', 'HOPELESS')}
    return tied, assigned, sum(1 for a, _ in o2.values() if a == 'ABSTAIN'), sum(1 for a, _ in o2.values() if a == 'HOPELESS')

arms = {n: load(p) for n, p in zip(['A', 'M10', 'M0'], sys.argv[1:4])}
print("arm\trows\ttied\tassigned\tABSTAIN\tHOPELESS\tcertified\tflagged")
C = {}
for n, (o2, cert, cand) in arms.items():
    tied, assigned, ab, ho = counts(o2)
    C[n] = (tied, assigned)
    print(f"{n}\t{len(o2)}\t{len(tied)}\t{len(assigned)}\t{ab}\t{ho}\t{cert}\t{cand if cand is not None else 'NA'}")
if 'M10' in C:
    (ta, aa), (tm, am) = C['A'], C['M10']
    untied = ta - tm
    both = set(aa) & set(am)
    conflict = sum(1 for k in both if aa[k] != am[k])
    b1 = len(untied) >= 0.05 * len(ta)
    b2 = len(am) >= 0.95 * len(aa)
    b3 = conflict < 0.01 * max(len(both), 1)
    print(f"## untied_under_M10 {len(untied)} of A_tied {len(ta)} ({100 * len(untied) / max(len(ta), 1):.1f}%) B1 {'PASS' if b1 else 'FAIL'}")
    print(f"## M10_assigned {len(am)} vs A {len(aa)} ({100 * len(am) / max(len(aa), 1):.1f}%) B2 {'PASS' if b2 else 'FAIL'}")
    print(f"## both_assigned {len(both)} conflicting {conflict} B3 {'PASS' if b3 else 'FAIL'}")
    print(f"## newly_tied_under_M10 {len(tm - ta)}; newly_assigned {len(set(am) - set(aa))}; lost_assignments {len(set(aa) - set(am))}")

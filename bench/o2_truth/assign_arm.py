#!/usr/bin/env python3
"""PREREG r1778 (docs/PREREG_assign_arm_2026-11-17.md): the assign stage's calls as an arm on neutral_truth's read table.

  assign_arm.py READS.tsv PREFIX.assign.assignments.tsv PREFIX.fam.copies.tsv"""
import sys, collections
from allele_graph import table

def call_of(loci, copy):
    """index of the read's placement on the copy's locus, -1 = off the placements"""
    c, s, e, st = copy
    for i, (lc, ls, le, lst, _) in enumerate(loci):
        if lc == c and lst == st and ls < e and s < le: return i
    return -1
assert call_of([('a', 0, 10, '+', 5), ('a', 20, 30, '+', 4)], ('a', 25, 40, '+')) == 1
assert call_of([('a', 0, 10, '+', 5)], ('a', 5, 9, '-')) == -1

def main(tsv, asg, copies):
    T = table(tsv); cp = {}
    it = (l.rstrip('\n').split('\t') for l in open(copies)); h = next(it)
    for f in it:
        d = dict(zip(h, f)); cp[(d['family_id'], d['copy_idx'])] = (d['chrom'], int(d['locus_start']), int(d['locus_end']), d['strand'])
    rows = collections.defaultdict(list)
    it = (l.rstrip('\n').split('\t') for l in open(asg)); h = next(it)
    for f in it:
        if f[0] in T: rows[f[0]].append(dict(zip(h, f)))
    C = collections.Counter()
    for q, d in T.items():
        a = [r for r in rows[q] if r['status'] == 'assigned']
        C[(d['st'], d['status'], 'status', 'assigned' if a else ('no row' if not rows[q] else '/'.join(sorted({r['status'] for r in rows[q]}))))] += 1
        calls = {call_of(d['loci'], cp[(r['family_id'], r['assigned_copy'])]) for r in a}
        call = calls.pop() if len(calls) == 1 else None  # rows that name different placements = no call
        k = (d['st'], d['status'])
        C[k + ('primary', 'right' if d['primary'] in d['origin'] else 'wrong')] += 1
        if call is None: C[k + ('assign', 'abstain')] += 1; fin = d['primary']
        else:
            fin = call; C[k + ('assign', 'right' if call in d['origin'] else 'wrong')] += 1
            if call == -1: C[k + ('assign', 'off')] += 1
            elif call != d['primary']: C[k + ('assign', 'moved', call in d['origin'])] += 1
        C[k + ('assign-or-primary', 'right' if fin in d['origin'] else 'wrong')] += 1
    for st in ('S2', 'S1'):
        n = sum(1 for d in T.values() if d['st'] == st and d['status'] == 'decidable')
        print(f"\n{st}: {n} decidable reads\narm\tright\twrong\tabstain\twrong_share_of_calls\tmoved_right\tmoved_wrong\toff_placements")
        for arm in ('primary', 'assign', 'assign-or-primary'):
            r, x, ab = (C[(st, 'decidable', arm, k)] for k in ('right', 'wrong', 'abstain'))
            print(f"{arm}\t{r}\t{x}\t{ab}\t{x / (r + x) if r + x else 0:.4f}\t{C[(st, 'decidable', arm, 'moved', True)]}\t{C[(st, 'decidable', arm, 'moved', False)]}\t{C[(st, 'decidable', arm, 'off')]}")
        print("class\treads\tassign_calls\tmoves_off_primary\toff_placements\tstatuses")
        for cl in ('decidable', 'outside', 'tie', 'no near-perfect hit'):
            n = sum(1 for d in T.values() if d['st'] == st and d['status'] == cl)
            stt = {k[3]: v for k, v in C.items() if k[:3] == (st, cl, 'status')}
            print(f"{cl}\t{n}\t{n - C[(st, cl, 'assign', 'abstain')]}\t{C[(st, cl, 'assign', 'moved', True)] + C[(st, cl, 'assign', 'moved', False)]}\t{C[(st, cl, 'assign', 'off')]}\t{stt}")

if __name__ == '__main__': main(*sys.argv[1:4])

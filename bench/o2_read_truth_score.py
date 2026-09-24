#!/usr/bin/env python3
"""Per-read scoring of copy_assign --families output (one row per read x family). Three readings of the table:
  OWN     : the row of the read's TRUE family (does the PSV certificate pick the right copy within the family?)
  PRIMARY : the row(s) with primary_local=1 (the family owning the read's primary placement) — the consumer rule
  ANY     : any assigned row; >1 assigned rows pointing to different loci = CONFLICT
usage: CATALOG_TSV=cat.copies.tsv o2_read_truth_score.py PREFIX O2PREFIX"""
import sys, csv, collections, subprocess, os, math
P, O = sys.argv[1], sys.argv[2]; CAT = os.environ['CATALOG_TSV']
cat = {}; div = {}
for r in csv.DictReader(open(CAT), delimiter='\t'):
    cat[(r['family_id'], r['copy_idx'])] = (r['chrom'], int(r['start']), int(r['end']))
    if r.get('max_family_identity') not in (None, '', 'NA'): div[(r['family_id'], r['copy_idx'])] = 1 - float(r['max_family_identity'])
def dbin(d):
    if d is None: return 'NA'
    return '<0.5%' if d < 0.005 else '0.5-1%' if d < 0.01 else '1-2%' if d < 0.02 else '2-5%' if d < 0.05 else '>=5%'
def same_locus(a, b):
    if a is None or b is None or a[0] != b[0]: return False
    o = min(a[2], b[2]) - max(a[1], b[1]); return o >= 0.5 * min(a[2] - a[1], b[2] - b[1])
prim = {}
for ln in subprocess.run(['samtools', 'view', '-F', '2308', P + '.bam'], capture_output=True, text=True).stdout.splitlines():
    f = ln.split('\t'); prim[f[0]] = int(f[4])
by = collections.defaultdict(list)
for r in csv.DictReader(open(O + '.assignments.tsv'), delimiter='\t'): by[r['read_name']].append(r)
def tr(n): return tuple(n.split('|')[:2])
def judge(rows_assigned, t):
    loci = {(r['family_id'], r['catalog_copy_idx']) for r in rows_assigned}
    if not loci: return 'abstain'
    ok = [k == t or same_locus(cat.get(k), cat.get(t)) for k in loci]
    if len(loci) > 1 and not all(ok): return 'conflict' if any(ok) else 'wrong'
    return 'correct' if all(ok) else 'wrong'
S = {v: collections.defaultdict(collections.Counter) for v in ('OWN', 'PRIMARY', 'ANY')}
n_mapq0 = 0; own_status = collections.Counter(); nprim = collections.Counter()
for name, mq in prim.items():
    if mq != 0: continue
    n_mapq0 += 1
    t = tr(name); b = dbin(div.get(t)); rows = by.get(name, [])
    asg = lambda rs: [r for r in rs if r['status'] == 'assigned' and r['origin_rejected'] == '0']
    own = [r for r in rows if r['family_id'] == t[0]]
    own_status[tuple(sorted(r['status'] for r in own)) or ('no_row',)] += 1
    o_own = 'lost' if not rows else ('no_own_row' if not own else judge(asg(own), t))
    pr = [r for r in rows if r['primary_local'] == '1']; nprim[len(pr)] += 1
    o_pr = 'lost' if not rows else ('no_primary_row' if not pr else judge(asg(pr), t))
    o_any = 'lost' if not rows else judge(asg(rows), t)
    for view, o in (('OWN', o_own), ('PRIMARY', o_pr), ('ANY', o_any)):
        for key in ('ALL', b): S[view][key][o] += 1
print(f'MAPQ-0 reads {n_mapq0}; own-family row status combos: {own_status.most_common(6)}; primary_local rows per read: {dict(nprim)}')
order = ['ALL', '<0.5%', '0.5-1%', '1-2%', '2-5%', '>=5%', 'NA']
for view in ('OWN', 'PRIMARY', 'ANY'):
    print(f'== {view}')
    for k in order:
        c = S[view].get(k)
        if not c: continue
        n = sum(c.values()); a = c['correct'] + c['wrong'] + c['conflict']
        acc = c['correct'] / a if a else float('nan')
        print(f"  {k:8s} n={n:5d} correct {c['correct']:4d} wrong {c['wrong']:4d} conflict {c['conflict']:4d} abstain {c['abstain']:4d} other {n - a - c['abstain']:4d} | acc_assigned {acc:.3f} coverage {a / n:.3f}")

#!/usr/bin/env python3
"""O2 LEVERS (register 1481, 2026-11-01): score an assign run against the BASE run by DISTINCT READS, not rows.
  levers_score.py BASE.assign.assignments.tsv ARM.assign.assignments.tsv [--names FILE] [--log ARM.assign.log] [--base-log LOG]
                  [--family-join ARM.assign.family_join.tsv]
Read bucket (one per read, precedence left to right over its rows): assigned > weak (ambiguous, origin_rejected 0, n_candidates >= 1)
> origin_rejected (ambiguous, origin_rejected 1) > orphan (ambiguous, n_candidates 0) > tied. --names: a read list (the < 5%
substitution class of the 2026-10-31 diagnosis) scored as its own block. --family-join: counts assigned rows whose copy tid starts
with `cand_` (a candidate copy). Logs: the `genome-wide unique-mapper agreement` line. Prints one TSV block per input; the
decision rule is in LEDGER_O2 `## PREREG O2 LEVERS`, not here."""
import sys, csv, re, collections
a = sys.argv[1:]
opt = lambda k: a[a.index(k) + 1] if k in a else None
names = set(open(opt('--names')).read().split()) if opt('--names') else None
ORDER = ['assigned', 'weak', 'origin_rejected', 'orphan', 'tied']
def bucket(r):
    s = r['status']
    if s == 'assigned': return 'assigned'
    if s == 'tied': return 'tied'
    if r['n_candidates'] == '0': return 'orphan'
    return 'origin_rejected' if r['origin_rejected'] == '1' else 'weak'
def load(path, fj=None):
    cand = set()
    if fj:
        for r in csv.DictReader(open(fj), delimiter='\t'):
            if r.get('copy_tid', r.get('tid', '')).startswith('cand_'): cand.add((r['family_id'], r['copy_idx']))
    best = {}; rows = 0; on_cand = 0; fams = collections.defaultdict(set)
    for r in csv.DictReader(open(path), delimiter='\t'):
        rows += 1; b = bucket(r); n = r['read_name']
        if n not in best or ORDER.index(b) < ORDER.index(best[n]): best[n] = b
        if b == 'assigned':
            fams[n].add(r['family_id'])
            if (r['family_id'], r['assigned_copy']) in cand: on_cand += 1
    return best, rows, on_cand, sum(1 for n in fams if len(fams[n]) > 1)
def agree(log):
    if not log: return 'NA'
    m = re.findall(r'unique-mapper agreement: (\d+)/(\d+)', open(log).read())
    return ' + '.join(f'{x}/{y}' for x, y in m) if m else 'NA'
pos = [x for i, x in enumerate(a) if not x.startswith('--') and (i == 0 or not a[i - 1].startswith('--'))]
for tag, path, log, fj in [('BASE', pos[0], opt('--base-log'), None), ('ARM', pos[1], opt('--log'), opt('--family-join'))]:
    best, rows, on_cand, dbl = load(path, fj)
    c = collections.Counter(best.values())
    print(f"## {tag}\t{path}\trows\t{rows}\treads\t{len(best)}\tassigned_rows_on_candidate\t{on_cand}\treads_assigned_in_2_families\t{dbl}\tuniq_agree\t{agree(log)}")
    for b in ORDER: print(f"{tag}\t{b}\t{c.get(b, 0)}")
    if names:
        cn = collections.Counter(best[n] for n in names if n in best)
        print(f"{tag}\tlt5_with_row\t{sum(cn.values())}\tof\t{len(names)}")
        for b in ORDER: print(f"{tag}\tlt5_{b}\t{cn.get(b, 0)}")

#!/usr/bin/env python3
"""O2 metric for the ambiguous reads (§6fr): for reads the aligner could not place, to which copy do they most
likely belong? Uses `<out>.posterior.tsv` (copy_assign --posterior). On simulated reads (name SIMNPIP|<true>|<i>):
per class — truth in the tie set, posterior top-1 accuracy, calibration of the max posterior. On real reads: the
distribution of the max posterior and tie-set size.  usage: o2_tied_metric.py <dir> [sim.bam]"""
import sys, csv, collections, statistics as st

def enforce_primary_local(rows, label="rows"):
    """register 734: a per-family read RATE must be taken over molecules with a PRIMARY alignment in a copy.
    `in_copy` is not enough — it fires on any aligned block, so genome-wide multimappers that only visit as
    SECONDARY alignments count as family reads. On DAZ that inflated the denominator 9.4x and turned 34.7%
    assigned into a reported 3.7%. Older tables have no `primary_local` column; those are passed through with
    a loud warning rather than silently rescored."""
    import sys
    if not rows:
        return rows
    if 'primary_local' not in rows[0]:
        print(f"  \u26a0 {label}: no `primary_local` column (table predates register 734) — rates are NOT enforced", file=sys.stderr)
        return rows
    keep = [r for r in rows if r['primary_local'] == '1']
    print(f"  primary-local enforcement: {len(keep)}/{len(rows)} {label} kept, {len(rows)-len(keep)} secondary-only visitors dropped", file=sys.stderr)
    return keep

d = sys.argv[1]; bam = sys.argv[2] if len(sys.argv) > 2 else None
A = {r['read_name']: r for r in enforce_primary_local(list(csv.DictReader(open(f'{d}/A.assignments.tsv'), delimiter='\t')), 'assignment rows')}
J = {r['copy_index']: r['catalog_copy_idx'] for r in csv.DictReader(open(f'{d}/A.family_join.tsv'), delimiter='\t')}
P = {}
for r in csv.DictReader(open(f'{d}/A.posterior.tsv'), delimiter='\t'):
    P[r['read_name']] = {J.get(k, k): float(v) for k, v in (x.split(':') for x in r['posterior'].split(','))}
mapq = {}
if bam:
    import pysam
    for a in pysam.AlignmentFile(bam).fetch(until_eof=True):
        if not (a.is_unmapped or a.is_secondary or a.is_supplementary): mapq[a.query_name] = a.mapping_quality
sim = all(n.startswith('SIMNPIP|') for n in list(A)[:20])
def truth(n): return n.split('|')[1] if sim else None
def report(sel, label):
    if not sel: print(f'{label}: 0 reads'); return
    tops = []; ks = []; in_tie = 0; top1 = 0; cal = collections.defaultdict(lambda: [0, 0])
    for n in sel:
        post = P.get(n)
        if not post: continue
        mx = max(post.values()); tie = [c for c, v in post.items() if abs(v - mx) < 1e-9]; tops.append(mx); ks.append(len(tie))
        if sim:
            t = truth(n); in_tie += t in tie; top1 += (len(tie) == 1 and tie[0] == t)
            b = min(9, int(mx * 10)); cal[b][0] += 1; cal[b][1] += (t in tie) / len(tie)
    msg = f"{label}: {len(tops)} reads | max posterior median {st.median(tops):.3f}, >= 0.9: {sum(1 for t in tops if t >= 0.9)} | tie-set size median {st.median(ks)}"
    if sim: msg += f" | truth in the tie set {in_tie}/{len(tops)} ({100*in_tie/len(tops):.1f} %), top-1 (unique max = truth) {top1} ({100*top1/len(tops):.1f} %)"
    print(msg)
    if sim:
        print('   calibration (max-posterior bin: reads, mean P(truth in tie)/|tie|):', {f'{b/10:.1f}-': (v[0], round(v[1] / v[0], 3)) for b, v in sorted(cal.items())})
inc = [n for n in A if A[n]['in_copy'] == 'true']
eq = [n for n in inc if A[n]['as_margin'] == '0' and A[n]['as_second'] not in ('NA', '')]
con = [n for n in inc if A[n].get('contested') == '1']
report(eq, 'equal best alignment scores (as_margin = 0)')
if mapq: report([n for n in inc if mapq.get(n, 60) == 0], 'MAPQ 0')
report(con, 'contested (MAPQ < 60), all')
report([n for n in con if A[n]['status'] == 'tied'], '  K = 0 ties')
report([n for n in con if A[n]['status'] == 'ambiguous' and A[n]['origin_rejected'] != '1'], '  uncertified pairs (evidence, not enough)')
report([n for n in con if A[n]['origin_rejected'] == '1'], '  certificate-rejected (no candidate explains the read)')
report([n for n in con if A[n]['status'] == 'assigned'], '  assigned')

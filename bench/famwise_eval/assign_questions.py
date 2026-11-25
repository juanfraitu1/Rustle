#!/usr/bin/env python3
"""PREREG r1949 (docs/PREREG_assign_two_questions_2026-11-25.md): statements about a molecule (primary, exact band, within
10, EM top, EM set) and labelled transcript claims against the origin of the reads.

  assign_questions.py --truth READS.tsv --tie-split T.tsv --bam SEL.bam --gtf LABELLED.gtf --shares SHARES.tsv.gz --out OUT.tsv"""
import argparse, collections, gzip, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import records as R
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
from gtf_names import chain

ap = argparse.ArgumentParser(); [ap.add_argument(x, required=True) for x in ('--truth', '--tie-split', '--bam', '--gtf', '--shares', '--out')]
a = ap.parse_args()
truth, recs = R.load_truth(a.truth, a.tie_split), R.load_records(a.bam)
both = sorted(set(truth) & set(recs))
# loci and labelled claims of the GTF
span, ex, at = {}, collections.defaultdict(list), {}
for l in open(a.gtf):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('exon', 'transcript'): continue
    g = re.search(r'gene_id "([^"]+)"', f[8]).group(1); s, e = int(f[3]) - 1, int(f[4])
    c = span.setdefault(g, [f[0], f[6], s, e]); c[2] = min(c[2], s); c[3] = max(c[3], e)
    t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
    if f[2] == 'exon': ex[t].append((int(f[3]), int(f[4])))
    else:
        m = re.search(r'o2_status "([^"]+)"', f[8]); at[t] = (f[0], f[6], m.group(1) if m else 'no_label')
claims = collections.defaultdict(set)   # (contig, strand, chain) -> labels
for t, e in ex.items():
    if len(e) > 1 and t in at: claims[(at[t][0], at[t][1], chain(e))].add(at[t][2])
shares = collections.defaultdict(list)
want = set(both)
with gzip.open(a.shares, 'rt') as fh:
    fh.readline()
    for l in fh:
        q, g, gap, x = l.rstrip('\n').split('\t')
        if q in want: shares[q].append((g, float(x)))

def locus_of(g, t):
    """truth place of a GTF locus (None = never judged, as records.locus_index)"""
    if g not in span: return None
    c, st, s, e = span[g]
    for i, (chrom, a0, b0, strand, sc) in enumerate(t.loci):
        if chrom == c and strand == st and a0 < e and s < b0:
            return i if sc >= R.KEEP * t.loci[0][4] or i in t.origin else None
    return None

def verdict_idx(idx, t): return None if None in idx else (idx == t.origin)
STMTS = ('primary', 'exact band', 'within 10', 'EM top', 'EM set')
n = collections.Counter(); sup = collections.defaultdict(list)
for q in both:
    t, rs = truth[q], recs[q]
    if t.kind in R.JUDGED:
        for r in R.statement(rs, 'ms', 10):
            k, i = (r.chrom, r.strand, r.chain), R.locus_index(r, t)
            if r.chain and k in claims and i is not None: sup[k].append((t.origin, i, t.kind))
    if not (R.near_tied(rs, 'ms') and t.kind in R.JUDGED): continue
    v = {'primary': verdict_idx({R.locus_index(r, t) for r in R.statement(rs, 'ms', None)}, t),
         'exact band': verdict_idx({R.locus_index(r, t) for r in R.statement(rs, 'ms', 0)}, t),
         'within 10': verdict_idx({R.locus_index(r, t) for r in R.statement(rs, 'ms', 10)}, t)}
    sh = shares.get(q)
    if sh:
        top = max(x for _, x in sh); n[t.kind, 'with share rows'] += 1
        v['EM top'] = verdict_idx({locus_of(g, t) for g, x in sh if x >= top - 0.001}, t)
        v['EM set'] = verdict_idx({locus_of(g, t) for g, x in sh if x >= 0.10}, t)
    else: v['EM top'] = v['EM set'] = v['primary']
    n[t.kind, 'molecules'] += 1
    for s in STMTS:
        if v[s] is None: n[t.kind, s, 'no verdict'] += 1
        else: n[t.kind, s, 'n'] += 1; n[t.kind, s, 'equals'] += v[s]
def sh(k, m):
    lo, hi = R.wilson(k, m); return f'{k}\t{m}\t{k / m:.3f}\t{lo:.3f}\t{hi:.3f}' if m else f'{k}\t0\tNA\tNA\tNA'
L = ['kind\tstatement\tequals\tjudged\tshare\tlo\thi\tno_verdict']
for kind in ('one', 'several', 'pooled'):
    for s in STMTS:
        ks = R.JUDGED if kind == 'pooled' else (kind,)
        L.append(f"{kind}\t{s}\t{sh(sum(n[k, s, 'equals'] for k in ks), sum(n[k, s, 'n'] for k in ks))}\t{sum(n[k, s, 'no verdict'] for k in ks)}")
L.append(f"# near-tied judged molecules: one {n['one', 'molecules']} (with share rows {n['one', 'with share rows']}), several {n['several', 'molecules']} (with share rows {n['several', 'with share rows']})")
# claims by label
cl = collections.Counter()
for k, labs in claims.items():
    s = sup.get(k)
    if not s: continue
    vd = 'TRUE' if any(o == {i} for o, i, _ in s) else 'SHARED' if any(i in o for o, i, _ in s) else 'PHANTOM'
    lab = 'expressed' if 'expressed' in labs else sorted(labs)[0]
    cl[lab, vd] += 1; cl[lab, 'several-member isoform'] += any(kd == 'several' for _, _, kd in s)
L.append('label\tjudged_claims\tTRUE\tSHARED\tPHANTOM\tshared_or_phantom\tlo\thi\tcarried_by_a_several_molecule')
for lab in sorted({k[0] for k in cl}):
    tot = sum(cl[lab, x] for x in ('TRUE', 'SHARED', 'PHANTOM')); bad = cl[lab, 'SHARED'] + cl[lab, 'PHANTOM']; lo, hi = R.wilson(bad, tot)
    L.append(f"{lab}\t{tot}\t{cl[lab, 'TRUE']}\t{cl[lab, 'SHARED']}\t{cl[lab, 'PHANTOM']}\t{bad / tot if tot else 0:.3f}\t{lo:.3f}\t{hi:.3f}\t{cl[lab, 'several-member isoform']}")
tt = sum(cl[l, 'TRUE'] for l in {k[0] for k in cl})
L.append(f"# TRUE claims held by expressed: {cl['expressed', 'TRUE']} of {tt}")
open(a.out, 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))

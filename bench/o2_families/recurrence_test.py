#!/usr/bin/env python3
"""IsoCon's decision logic, ported to reference coordinates: which sequences among reference-identical copies are real (register 1926).

    recurrence_test.py ASSEMBLY.gtf BAM AS_TABLE GENOME.fa OUT.tsv [--alpha 0.01] [--min-test-ratio 5] [--threads N]

Groups and tied reads as sample_psv.py (loci joined by reads exactly tied on ms; the tied reads read at one locus).
What is taken from IsoCon (Sahlin et al. 2018; modules hypothesis_test_module.py, functions.py, isocon_statistical_test.py):
 1. candidates = the distinct sequences the reads show after errors are corrected. Here: a difference from the reference
    seen in one read only is an error and is dropped (IsoCon corrects a read towards its nearest neighbours); a read's
    pattern is its alleles at the remaining sites it covers; a pattern contained in exactly one longer pattern is a
    shorter read of it (IsoCon's end-invariance), in several it is ambiguous and takes part in no test.
 2. the test, per candidate c and each nearest candidate t (fewest differing sites), from the least supported up, only
    when c has fewer than --min-test-ratio x t's reads: D = the sites where c and t differ; support = the reads of c or t
    covering all of D that carry c's allele at EVERY site of D; under the null all those reads come from t, and read i
    shows D by error with p_i = product over D of its own error rates (substitution rate / 3, insertion rate / 4,
    deletion rate; an indel in a homopolymer of length u counts u times; rates from the read's own errors, never under
    |D| / length); the p-value is Raghavan's bound on the weighted sum of Bernoullis, multiplied by the number of ways
    D could have been placed ((4 (m + 1)) ^ nI x C(m, nD) x C(3 (m - nD), nS)).
 0. isoforms first: the reads of a group are split by intron chain; a chain that is a contiguous piece of exactly one
    longer chain is a shorter read of it; steps 1-3 run per isoform. A sequence needs >= --min-reads reads (2).
 3. c is kept when every such test gives p < --alpha; otherwise its reads go to its nearest candidate and the tests
    are repeated until nothing changes.
Not taken: IsoCon's all-against-all alignment and consensus (the genome alignment at one locus replaces them).
OUT.tsv: group, loci, read_at, candidate, reads, sites (its alleles), p (largest of its tests, NA untested), ambiguous
reads of the locus, read names, isoform (index of the intron chain at that locus).
It says how many different sequences are expressed among copies the reference cannot separate, never which copy."""
import sys, os, argparse, collections, math
from decimal import Decimal, getcontext
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sample_psv as S
getcontext().prec = 60

def raghavan(prob, support):
    """IsoCon's raghavan_upper_pvalue_bound: prob {read: p_i}, support = reads observed with the variant set"""
    if not support: return 1.0
    lp = {q: -math.log10(p) for q, p in prob.items()}; mx = max(lp.values())
    w = {q: lp[q] / mx for q in lp}
    m = Decimal(sum(w[q] * prob[q] for q in prob)); y = Decimal(sum(w[q] for q in support))
    if m == 0: return 0.0
    d = y / m - 1
    if d <= 0: return 1.0 if d < 0 else 0.5
    k = m * d
    try: return float(k.exp() / (d + 1) ** (k + k / d))
    except Exception: return 0.0

def choose(n, k): return math.comb(max(n, 0), k) if k <= n else 1

def factor(m, D):
    nS = sum(1 for x in D if x[1] == 'S'); nD = sum(1 for x in D if x[1] == 'D'); nI = sum(1 for x in D if x[1] == 'I')
    return (4 * (m + 1)) ** nI * choose(m, nD) * choose(3 * (m - nD), nS)

def p_read(err, length, D, mult):
    """err = (ins, del, sub) of the read; D = variant sites; mult {site: homopolymer multiplicity}"""
    n = float(len(D)); ins, de, sub = err
    pS = max(sub, n) / length / 3.0; pI = max(ins, n) / length / 4.0; pD = max(de, n) / length
    p = 1.0
    for x in D:
        u = mult.get(x, 1)
        p *= pS * u if x[1] == 'S' else min(0.5, pI * u) if x[1] == 'I' else min(0.5, pD * u)
    return min(max(p, 1e-300), 0.99999)

def decide(reads, mult, alpha=0.01, ratio=5):
    """reads {name: (variants, blocks)} -> candidates [{'cov': frozenset, 'pat': {site: bool}, 'reads': [..], 'p': float|None}], n ambiguous"""
    sup = collections.Counter(x for v, _ in reads.values() for x in v)
    sites = sorted(x for x, n in sup.items() if n >= 2)
    def covered(bl, pos): return any(a <= pos < b for a, b in bl)
    pat = {}; err = {}; length = {}
    for q, (v, bl) in reads.items():
        pat[q] = {x: (x in v) for x in sites if covered(bl, x[0])}
        e = [x for x in v if sup[x] < 2]
        err[q] = (sum(1 for x in e if x[1] == 'I'), sum(1 for x in e if x[1] == 'D'), sum(1 for x in e if x[1] == 'S'))
        length[q] = max(1, sum(b - a for a, b in bl))
    by = collections.defaultdict(list)
    for q, p in pat.items(): by[(frozenset(p), tuple(sorted(p.items())))].append(q)
    order = sorted(by, key=lambda k: (-len(k[0]), -len(by[k]), k[1]))
    cands = []; amb = 0
    for k in order:
        cov, items = k; p = dict(items)
        host = [c for c in cands if cov <= c['cov'] and all(c['pat'][x] == a for x, a in p.items())]
        if len(host) == 1: host[0]['reads'] += by[k]
        elif len(host) > 1: amb += len(by[k])
        else: cands.append({'cov': cov, 'pat': p, 'reads': list(by[k]), 'p': None})
    changed = True
    while changed:
        changed = False
        for c in sorted(cands, key=lambda c: (len(c['reads']), sorted(c['pat'].items()))):
            nb = []
            for t in cands:
                if t is c: continue
                D = [x for x in c['cov'] & t['cov'] if c['pat'][x] != t['pat'][x]]
                if D: nb.append((len(D), t, D))
            if not nb: continue
            mn = min(n for n, _, _ in nb); worst = None; tested = False
            for n, t, D in nb:
                if n != mn or len(c['reads']) >= ratio * len(t['reads']): continue
                tested = True
                R = [q for q in c['reads'] + t['reads'] if all(x in pat[q] for x in D)]
                support = [q for q in R if all(pat[q][x] == c['pat'][x] for x in D)]
                if not R: pv = 1.0
                else:
                    prob = {q: p_read(err[q], length[q], D, mult) for q in R}
                    m = sorted(length[q] for q in R)[len(R) // 2]
                    pv = min(1.0, raghavan(prob, support) * factor(m, D))
                if worst is None or pv > worst[0]: worst = (pv, t)
            if not tested: continue
            c['p'] = worst[0]
            if worst[0] >= alpha:
                worst[1]['reads'] += c['reads']; cands.remove(c); changed = True; break
    return cands, amb, sites

def introns(r):
    out = []; pos = r.reference_start
    for op, n in r.cigartuples:
        if op == 3: out.append((pos, pos + n)); pos += n
        elif op in (0, 2, 7, 8): pos += n
    return tuple(out)

def by_isoform(chains):
    """chains {read: intron chain} -> {full chain: [reads]}, ambiguous reads. A chain that is a contiguous piece of exactly
    one longer chain is a shorter read of that isoform (IsoCon's end-invariance); of several, it is ambiguous."""
    cnt = collections.Counter(chains.values()); full = []
    def inside(a, b): return len(a) < len(b) and any(b[i:i + len(a)] == a for i in range(len(b) - len(a) + 1))
    for c in cnt:
        if not any(inside(c, d) for d in cnt): full.append(c)
    out = collections.defaultdict(list); amb = []
    for q, c in chains.items():
        if c in full: out[c].append(q); continue
        host = [d for d in full if inside(c, d)]
        if len(host) == 1: out[host[0]].append(q)
        else: amb.append(q)
    return out, amb

def selftest():
    A = ((10, 20), (30, 40)); B = ((10, 20), (50, 60)); o, amb = by_isoform({'a1': A, 'a2': A, 'b1': B, 'p1': ((30, 40),), 'p2': ((10, 20),), 'm': ()})
    assert sorted(o[A]) == ['a1', 'a2', 'p1'] and o[B] == ['b1'] and sorted(amb) == ['m', 'p2'], (dict(o), amb)
    bl = [(0, 1500)]; R = {}
    for i in range(40): R[f't{i}'] = (set(), bl)
    for i in range(12): R[f'c{i}'] = ({(700, 'S', 'G')}, bl)                 # a real second sequence, 12 reads
    R['e1'] = ({(100, 'S', 'A')}, bl); R['e2'] = ({(100, 'S', 'A')}, bl)       # the same error twice: must not be a copy
    R['e3'] = ({(300, 'S', 'C')}, bl)                                          # a singleton: corrected
    c, amb, s = decide(R, {})
    assert sorted(len(x['reads']) for x in c) == [12, 43], [(len(x['reads']), x['p']) for x in c]
    P = {f'a{i}': ({(10, 'S', 'T')}, [(0, 1000)]) for i in range(20)}; P.update({f'b{i}': (set(), [(0, 1000)]) for i in range(20)})
    P.update({f'p{i}': (set(), [(500, 1000)]) for i in range(6)})             # partial reads covering no site: ambiguous, not a copy
    c, amb, s = decide(P, {})
    assert len(c) == 2 and amb == 6, (len(c), amb)
    assert raghavan({'a': 1e-4, 'b': 1e-4}, []) == 1.0 and raghavan({f'r{i}': 2e-4 for i in range(60)}, [f'r{i}' for i in range(20)]) < 1e-30
    print('recurrence_test selftest ok')

def main():
    if sys.argv[1:] == ['--selftest']: return selftest()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gtf'); ap.add_argument('bam'); ap.add_argument('tab'); ap.add_argument('fasta'); ap.add_argument('out')
    ap.add_argument('--alpha', type=float, default=0.01); ap.add_argument('--min-test-ratio', type=float, default=5); ap.add_argument('--threads', type=int, default=8); ap.add_argument('--min-reads', type=int, default=2)
    a = ap.parse_args()
    import pysam
    meta, sets, groups, gr = S.tie_groups(a.gtf, a.bam, a.tab, a.threads)
    B = pysam.AlignmentFile(a.bam); F = pysam.FastaFile(a.fasta); n = 0
    with open(a.out + '.tmp', 'w') as w:
        w.write('group\tloci\tread_at\tcandidate\treads\talleles\tp\tambiguous_reads\tread_names\tisoform\n')
        for k in sorted(groups, key=lambda k: (meta[k][0], meta[k][1])):
            left = set(gr[k])
            while len(left) >= a.min_reads:      # a chained group is read at as many of its loci as it takes to cover its reads
                reach = collections.Counter(g for q in left for g in sets[q]); rep = max(sorted(reach), key=lambda g: reach[g])
                c, s, e, st = meta[rep][:4]; reads = {}
                for r in B.fetch(c, s - 1, e):
                    q = r.query_name
                    if r.is_unmapped or r.is_supplementary or q not in left or rep not in sets[q]: continue
                    if r.is_secondary and r.get_tag('ms') != S.E.G_TIES.get(q): continue
                    if S.E.tx_strand(r) != st and st in '+-': continue
                    if q not in reads or r.query_alignment_length > reads[q][2]: reads[q] = (*S.variants(r), r.query_alignment_length, introns(r))
                left -= {q for q in left if rep in sets[q]}
                if not reads: continue
                iso, amb_iso = by_isoform({q: v[3] for q, v in reads.items()})
                reads = {q: v[:2] for q, v in reads.items()}
                mult = {}
                for x in {x for v, _ in reads.values() for x in v if x[1] in 'ID'}:
                    seq = F.fetch(c, max(0, x[0] - 20), x[0] + 21).upper(); i = min(20, x[0]); b = seq[i] if i < len(seq) else 'N'; u = 1; j = i + 1
                    while j < len(seq) and seq[j] == b: u += 1; j += 1
                    j = i - 1
                    while j >= 0 and seq[j] == b: u += 1; j -= 1
                    mult[x] = u
                cands = []; amb = len(amb_iso)
                for ci, ch in enumerate(sorted(iso), 1):           # one decision per isoform (intron chain)
                    cs, am, _ = decide({q: reads[q] for q in iso[ch]}, mult, a.alpha, a.min_test_ratio); amb += am
                    for x in cs: x['iso'] = ci
                    cands += [x for x in cs if len(x['reads']) >= a.min_reads]
                for i, cd in enumerate(sorted(cands, key=lambda x: -len(x['reads'])), 1):
                    al = ','.join(f'{x[0] + 1}{x[1]}{x[2]}' for x in sorted(cd['pat']) if cd['pat'][x]) or 'reference'
                    w.write(f"{k}\t{','.join(sorted(groups[k]))}\t{rep}\t{i}\t{len(cd['reads'])}\t{al}\t{'NA' if cd['p'] is None else format(cd['p'], '.3g')}\t{amb}\t{','.join(cd['reads'])}\t{cd['iso']}\n"); n += 1
    os.replace(a.out + '.tmp', a.out)
    print(f'recurrence_test: {len(groups)} groups, {n} sequences kept -> {a.out}')

if __name__ == '__main__':
    main()

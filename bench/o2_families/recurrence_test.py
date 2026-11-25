#!/usr/bin/env python3
"""IsoCon's decision logic, ported to reference coordinates: which sequences among reference-identical copies are real (register 1926).

    recurrence_test.py ASSEMBLY.gtf BAM AS_TABLE GENOME.fa OUT.tsv [--alpha 0.01] [--min-test-ratio 5] [--threads N]
                       [--cache PATH]   keep the reads collected from the BAM; later runs on the same inputs reuse them

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
 1b. --cluster cc (register 1932, opt-in): candidates by correlation clustering of the reads instead of the folding of
    step 1 (pairwise log-likelihood ratio of one sequence against two over shared sites; see correlation_clusters).
 2b. --qual isocon|phred (register 1933, opt-in): IsoCon's quality model in place of the error counts: the probability of
    each variant from the read's base quality there, and no placement factor (see p_read_qual).
 0b. --max-patterns N (default 1000; 0 = no limit): an isoform whose reads show more than N distinct variant patterns is
    not tested and is listed in OUT.skipped.tsv. Such a locus is hypervariable (immunoglobulin and T-cell receptor genes
    are rearranged and mutated in every clone): its many sequences are real and are not copies. HG002's lambda light
    chain locus has 26,776 tied reads with 4,778 patterns in one isoform and took 43 of the tool's 47 minutes.
 3. c is kept when every such test gives p < --alpha; otherwise its reads go to its nearest candidate and the tests
    are repeated until nothing changes.
Not taken: IsoCon's all-against-all alignment and consensus (the genome alignment at one locus replaces them).
OUT.tsv: group, loci, read_at, candidate, reads, sites (its alleles), p (largest of its tests, NA untested), ambiguous
reads of the locus, read names, isoform (index of the intron chain at that locus).
It says how many different sequences are expressed among copies the reference cannot separate, never which copy."""
import sys, os, re, argparse, bisect, collections, math, multiprocessing
from decimal import Decimal, getcontext
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools'))
import eq_classes as E
import tie_table as TT
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

def site_err(err, length, x, mult):
    """probability that a read shows (or hides) the variant x by error; IsoCon's per-variant rates with |D| = 1"""
    ins, de, sub = err; u = mult.get(x, 1)
    return min(0.25, (max(sub, 1) / length / 3.0) * u if x[1] == 'S' else (max(ins, 1) / length / 4.0) * u if x[1] == 'I' else (max(de, 1) / length) * u)

def correlation_clusters(pat, err, length, mult, sites, rounds=10):
    """Correlation clustering of the reads (register 1932). Nodes = reads with the same covered sites and alleles. Two nodes
    are scored by the log-likelihood ratio of 'one sequence' against 'two sequences' over the sites both cover (agree: the
    chance two different sequences agree there, from the allele frequency; disagree: the chance of an error in either).
    Greedy merging of the pair of clusters with the largest positive total, then single-node moves. Reads covering no site
    join the cluster when there is exactly one, else they are ambiguous. -> [[reads]], ambiguous reads"""
    by = collections.defaultdict(list)
    for q, p in pat.items(): by[(frozenset(p), tuple(sorted(p.items())))].append(q)
    blank = [q for k, qs in by.items() if not k[0] for q in qs]
    keys = sorted((k for k in by if k[0]), key=lambda k: (-len(by[k]), k[1]))
    if not keys: return ([blank] if blank else []), []
    cov = collections.Counter(); alt = collections.Counter()
    for q, p in pat.items():
        for x, a in p.items(): cov[x] += 1; alt[x] += a
    f = {x: min(1 - 1e-6, max(1e-6, alt[x] / cov[x])) for x in cov}
    n = len(keys); size = [len(by[k]) for k in keys]; P = [dict(k[1]) for k in keys]
    E = [{x: sum(site_err(err[q], length[q], x, mult) for q in by[k]) / len(by[k]) for x in k[0]} for k in keys]
    w = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            t = 0.0
            for x in P[i].keys() & P[j].keys():
                e = min(0.49, E[i][x] + E[j][x])
                t += math.log((1 - e) / (f[x] ** 2 + (1 - f[x]) ** 2)) if P[i][x] == P[j][x] else math.log(e / (2 * f[x] * (1 - f[x])))
            w[i][j] = w[j][i] = t * size[i] * size[j]
    cl = {i: [i] for i in range(n)}; W = {i: {j: w[i][j] for j in range(n) if j != i and w[i][j] != 0.0} for i in range(n)}   # cluster-to-cluster totals, kept up to date
    while len(cl) > 1:
        best = max(((v, a, b) for a, d in W.items() for b, v in d.items() if a < b), key=lambda z: z[0], default=(0, 0, 0))
        if best[0] <= 0: break
        _, a, b = best; cl[a] += cl.pop(b)
        for c, v in W.pop(b).items():
            if c == a: continue
            W[a][c] = W[a].get(c, 0.0) + v; W[c][a] = W[a][c]; W[c].pop(b, None)
        W[a].pop(b, None)
    cl = list(cl.values())
    for _ in range(rounds):
        moved = False
        for i in range(n):
            cur = next(k for k, c in enumerate(cl) if i in c)
            sc = [sum(w[i][b] for b in c if b != i) for c in cl]
            k = max(range(len(cl)), key=lambda z: sc[z])
            if sc[k] > sc[cur] + 1e-9 and k != cur: cl[cur].remove(i); cl[k].append(i); moved = True
            elif sc[cur] < -1e-9 and max(sc) <= 0 and len(cl[cur]) > 1: cl[cur].remove(i); cl.append([i]); moved = True
            cl = [c for c in cl if c]
        if not moved: break
    out = [[q for i in c for q in by[keys[i]]] for c in cl]
    if len(out) == 1: out[0] += blank; blank = []
    return out, blank

def qual_blocks(r):
    """base qualities of a record for the quality model (register 1933): -> (aligned blocks with their query offset, qualities).
    As IsoCon's fix_quality_values, the qualities inside a homopolymer run of the read are sorted (lowest first): the
    uncertainty of a run belongs to its length, not to one base."""
    q = list(r.query_qualities) if r.query_qualities is not None else None
    if q is None: return None
    seq = r.query_sequence; i = 0
    while i < len(q):
        j = i
        while j + 1 < len(q) and seq[j + 1] == seq[i]: j += 1
        if j > i: q[i:j + 1] = sorted(q[i:j + 1])
        i = j + 1
    bl = []; rp, qp = r.reference_start, 0
    for op, n in r.cigartuples:
        if op in (0, 7, 8): bl.append((rp, rp + n, qp)); rp += n; qp += n
        elif op in (1, 4): qp += n
        elif op in (2, 3): rp += n
    return bl, bytes(min(x, 255) for x in q)

def qual_at(qb, pos):
    """quality of the read base at reference position pos; for a position inside a gap of the read, the lower of the flanking bases"""
    bl, q = qb; lo = None
    for a, b, o in bl:
        if a <= pos < b: return q[o + pos - a]
        if b <= pos: lo = q[o + b - a - 1]
        elif a > pos: return min(lo, q[o]) if lo is not None else q[o]
    return lo if lo is not None else 0

def p_read_qual(qb, D, mult, ratios, mode, max_q=43):
    """IsoCon's get_read_ccs_probabilities: per variant 10^(-Q/10) x the share of that error type among the reads' errors
    (/ 3 for a substitution, / 4 for an insertion); in a homopolymer the full 10^(-Q/10). mode 'isocon': Q rescaled as
    IsoCon does, (Q - 3) x (max_q - 3) / 90 + 3; mode 'phred': Q capped at max_q (for binned qualities that top out near 40)."""
    p = 1.0
    for x in D:
        Q = qual_at(qb, x[0]); Q = (Q - 3) * (max_q - 3.0) / 90.0 + 3 if mode == 'isocon' else min(Q, max_q)
        e = 10 ** (-max(Q, 1) / 10.0)
        p *= e if mult.get(x, 1) > 1 else e * ratios[0] / 3.0 if x[1] == 'S' else e * ratios[1] / 4.0 if x[1] == 'I' else e * ratios[2]
    return min(max(p, 1e-300), 0.99999)

def decide(reads, mult, alpha=0.01, ratio=5, cluster='fold', qinfo=None, qmode=None, max_patterns=0):
    """reads {name: (variants, blocks)} -> candidates [{'cov': frozenset, 'pat': {site: bool}, 'reads': [..], 'p': float|None}], n ambiguous"""
    sup = collections.Counter(x for v, _ in reads.values() for x in v)
    sites = sorted(x for x, n in sup.items() if n >= 2)
    spos = [x[0] for x in sites]
    pat = {}; err = {}; length = {}
    for q, (v, bl) in reads.items():
        p = {}
        for a0, b0 in bl:          # the sites inside each aligned block, by position (blocks are in order and do not overlap)
            for x in sites[bisect.bisect_left(spos, a0):bisect.bisect_left(spos, b0)]: p[x] = x in v
        pat[q] = p
        e = [x for x in v if sup[x] < 2]
        err[q] = (sum(1 for x in e if x[1] == 'I'), sum(1 for x in e if x[1] == 'D'), sum(1 for x in e if x[1] == 'S'))
        length[q] = max(1, sum(b - a for a, b in bl))
    if max_patterns:            # a hypervariable locus (immunoglobulin, T-cell receptor): thousands of real sequences, not copies
        npat = len({(frozenset(p), tuple(sorted(p.items()))) for p in pat.values()})
        if npat > max_patterns: return None, npat, sites
    if cluster == 'cc':
        groups, blank = correlation_clusters(pat, err, length, mult, sites); amb = len(blank)
        def consensus(qs):
            cn = collections.defaultdict(collections.Counter)
            for q in qs:
                for x, v in pat[q].items(): cn[x][v] += 1
            return {x: c[True] > c[False] for x, c in cn.items()}
        groups.sort(key=len, reverse=True); again = True
        while again:            # two clusters whose consensus agrees at every site they share (and share one) are one sequence
            again = False; cons = [consensus(qs) for qs in groups]
            for i in range(len(groups)):
                for j in range(i + 1, len(groups)):
                    sh = cons[i].keys() & cons[j].keys()
                    if sh and all(cons[i][x] == cons[j][x] for x in sh): groups[i] += groups.pop(j); again = True; break
                if again: break
        cands = [{'cov': frozenset(c), 'pat': c, 'reads': list(qs), 'p': None} for qs, c in zip(groups, [consensus(qs) for qs in groups])]
    else:
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
                    if qmode:               # IsoCon with qualities: per-base probabilities, and no placement factor
                        R = [q for q in R if qinfo.get(q) is not None]; Rs = set(R); support = [q for q in support if q in Rs]
                        ti = max(1.0, sum(err[q][0] for q in R)); td = max(1.0, sum(err[q][1] for q in R)); ts = max(1.0, sum(err[q][2] for q in R)); tt = ti + td + ts
                        prob = {q: p_read_qual(qinfo[q], D, mult, (ts / tt, ti / tt, td / tt), qmode) for q in R}
                        pv = min(1.0, raghavan(prob, support)) if R else 1.0
                    else:
                        prob = {q: p_read(err[q], length[q], D, mult) for q in R}
                        m = sorted(length[q] for q in R)[len(R) // 2]
                        rb = raghavan(prob, support)          # the placement factor is an exact integer and can exceed any float: combine in logs
                        pv = 0.0 if rb <= 0 else min(1.0, 10 ** min(0.0, math.log10(rb) + math.log10(factor(m, D))))
                if worst is None or pv > worst[0]: worst = (pv, t)
            if not tested: continue
            c['p'] = worst[0]
            if worst[0] >= alpha:
                worst[1]['reads'] += c['reads']; cands.remove(c); changed = True; break
    return cands, amb, sites

def variants(r):
    """-> set of (pos, kind, alt) on the reference, list of aligned blocks (0-based half open)"""
    out = set(); rp, qp = r.reference_start, 0; seq = r.query_sequence
    for op, n in r.cigartuples:
        if op in (0, 7): rp += n; qp += n
        elif op == 8:
            for i in range(n): out.add((rp + i, 'S', seq[qp + i] if seq else 'N'))
            rp += n; qp += n
        elif op == 1: out.add((rp, 'I', seq[qp:qp + n] if seq else str(n))); qp += n
        elif op == 2: out.add((rp, 'D', n)); rp += n
        elif op == 3: rp += n
        elif op == 4: qp += n
    return out, r.get_blocks()

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
    big = [(i * 7, 'S', 'A') for i in range(400)]                              # 400 separating sites: the factor has > 1,000 digits
    assert factor(3000, big) > 10 ** 400 and math.log10(factor(3000, big)) > 400
    Q = {f'x{i}': (set(big), [(0, 3000)]) for i in range(10)}; Q.update({f'y{i}': (set(), [(0, 3000)]) for i in range(30)})
    c, amb, s_ = decide(Q, {}); assert len(c) == 2, len(c)
    # correlation clustering: two sequences differing at 2 sites, a partial read covering one site joins the right one
    bl = [(0, 1500)]; C = {f'u{i}': ({(100, 'S', 'A'), (900, 'S', 'C')}, bl) for i in range(10)}; C.update({f'v{i}': (set(), bl) for i in range(10)})
    C['part'] = ({(900, 'S', 'C')}, [(800, 1500)]); C['none'] = (set(), [(1200, 1500)])
    c, amb, s_ = decide(C, {}, cluster='cc')
    assert sorted(len(x['reads']) for x in c) == [10, 11] and amb == 1, ([len(x['reads']) for x in c], amb)
    assert any('part' in x['reads'] and 'u0' in x['reads'] for x in c)
    c, amb, s_ = decide({f'r{i}': (set(), bl) for i in range(9)}, {}, cluster='cc'); assert len(c) == 1 and len(c[0]['reads']) == 9 and amb == 0
    # quality model: the same 4 supporting reads pass with high base qualities at the site and fail with low ones
    bl = [(0, 1500)]; G = {f't{i}': (set(), bl) for i in range(40)}; G.update({f'c{i}': ({(700, 'S', 'G')}, bl) for i in range(4)})
    hi = {q: ([(0, 1500, 0)], [90] * 1500) for q in G}; lo = {q: ([(0, 1500, 0)], [8] * 1500) for q in G}
    assert len(decide(G, {}, qinfo=hi, qmode='phred')[0]) == 2 and len(decide(G, {}, qinfo=lo, qmode='phred')[0]) == 1
    assert qual_at(([(0, 10, 0), (12, 20, 10)], list(range(18))), 11) == 9 and qual_at(([(0, 10, 0)], list(range(10))), 3) == 3
    hv = {f'h{i}': ({(10 * i, 'S', 'A'), (10 * i + 5, 'S', 'C')} if i % 2 else {(10 * (i - 1), 'S', 'A'), (10 * (i - 1) + 5, 'S', 'C')}, [(0, 1500)]) for i in range(40)}   # 20 patterns, 2 reads each
    r0 = decide(hv, {}, max_patterns=10); assert r0[0] is None and r0[1] == 20, r0[:2]
    assert decide(hv, {}, max_patterns=0)[0] is not None and decide(hv, {}, max_patterns=1000)[0] is not None
    sig = [['/x.bam', 5, 9], 2, 'v1']
    assert cache_ok({'signature': sig, 'qualities': False}, sig, False) and not cache_ok({'signature': sig, 'qualities': False}, sig, True)
    assert not cache_ok({'signature': [['/x.bam', 6, 9], 2, 'v1'], 'qualities': True}, sig, False) and not cache_ok({}, sig, False)
    print('recurrence_test selftest ok')

def main():
    if sys.argv[1:] == ['--selftest']: return selftest()
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gtf'); ap.add_argument('bam'); ap.add_argument('tab'); ap.add_argument('fasta'); ap.add_argument('out')
    ap.add_argument('--alpha', type=float, default=0.01); ap.add_argument('--min-test-ratio', type=float, default=5); ap.add_argument('--threads', type=int, default=8); ap.add_argument('--min-reads', type=int, default=2); ap.add_argument('--cluster', choices=('fold', 'cc'), default='fold'); ap.add_argument('--qual', choices=('none', 'isocon', 'phred'), default='none'); ap.add_argument('--cache', metavar='PATH'); ap.add_argument('--max-patterns', type=int, default=1000); ap.add_argument('--ties', metavar='FILE', help='the shared tie table (tools/tie_table.py)')
    a = ap.parse_args()
    units, ngroups = collect(a)
    _W['units'] = units; _W['a'] = a
    order = sorted(range(len(units)), key=lambda i: -len(units[i][3]))          # largest groups first, so the pool balances
    if a.threads > 1 and len(units) > 1:
        with multiprocessing.get_context('fork').Pool(a.threads) as pool: res = dict(zip(order, pool.map(_decide_unit, order, chunksize=1)))
    else: res = {i: _decide_unit(i) for i in order}
    n = 0
    with open(a.out + '.tmp', 'w') as w:
        w.write('group\tloci\tread_at\tcandidate\treads\talleles\tp\tambiguous_reads\tread_names\tisoform\n')
        for i in range(len(units)): w.writelines(res[i][0]); n += len(res[i][0])
    os.replace(a.out + '.tmp', a.out)
    sk = [l for i in range(len(units)) for l in res[i][1]]; skp = re.sub(r'\.tsv$', '', a.out) + '.skipped.tsv'
    with open(skp + '.tmp', 'w') as w: w.write('group\tloci\tread_at\tisoform\treads\tsites\tpatterns\treason\n'); w.writelines(sk)
    os.replace(skp + '.tmp', skp)
    print(f'recurrence_test: {ngroups} groups, {n} sequences kept' + (f', {len(sk)} isoforms not tested as hypervariable ({skp})' if sk else '') + f' -> {a.out}')

_W = {}

def _decide_unit(i):
    """the decisions of one locus read (independent of every other: run in parallel) -> its output lines"""
    a = _W['a']; k, loci, rep, rd, mult = _W['units'][i]
    iso, amb_iso = by_isoform({q: v[2] for q, v in rd.items()})
    qinfo = {q: v[3] for q, v in rd.items()} if a.qual != 'none' else None
    reads = {q: v[:2] for q, v in rd.items()}
    cands = []; amb = len(amb_iso); skipped = []
    for ci, ch in enumerate(sorted(iso), 1):           # one decision per isoform (intron chain)
        cs, am, st = decide({q: reads[q] for q in iso[ch]}, mult, a.alpha, a.min_test_ratio, a.cluster, qinfo, None if a.qual == 'none' else a.qual, a.max_patterns)
        if cs is None: skipped.append(f"{k}\t{loci}\t{rep}\t{ci}\t{len(iso[ch])}\t{len(st)}\t{am}\thypervariable: more than {a.max_patterns} distinct variant patterns in one isoform\n"); continue
        amb += am
        for x in cs: x['iso'] = ci
        cands += [x for x in cs if len(x['reads']) >= a.min_reads]
    out = []
    for n, cd in enumerate(sorted(cands, key=lambda x: -len(x['reads'])), 1):      # n, not i: i is this function's unit index
        al = ','.join(f'{x[0] + 1}{x[1]}{x[2]}' for x in sorted(cd['pat']) if cd['pat'][x]) or 'reference'
        out.append(f"{k}\t{loci}\t{rep}\t{n}\t{len(cd['reads'])}\t{al}\t{'NA' if cd['p'] is None else format(cd['p'], '.3g')}\t{amb}\t{','.join(cd['reads'])}\t{cd['iso']}\n")
    return out, skipped

def _collect_group(k):
    """the tied reads of one group, read at as many of its loci as it takes to cover them -> [(group, loci, locus, reads, mult)]"""
    import pysam
    a = _W['a']; meta, sets, groups, gr, wq = _W['meta'], _W['sets'], _W['groups'], _W['gr'], _W['wq']
    B = pysam.AlignmentFile(a.bam); F = pysam.FastaFile(a.fasta); units = []
    left = set(gr[k])
    while len(left) >= a.min_reads:
        reach = collections.Counter(g for q in left for g in sets[q]); rep = max(sorted(reach), key=lambda g: reach[g])
        c, s, e, st = meta[rep][:4]; reads = {}; ln = {}
        for r in B.fetch(c, s - 1, e):
            q = r.query_name
            if r.is_unmapped or r.is_supplementary or q not in left or rep not in sets[q]: continue
            if r.is_secondary and r.get_tag('ms') != E.G_TIES.get(q): continue
            if E.tx_strand(r) != st and st in '+-': continue
            if q not in reads or r.query_alignment_length > ln[q]: reads[q] = (*variants(r), introns(r), qual_blocks(r) if wq else None); ln[q] = r.query_alignment_length
        left -= {q for q in left if rep in sets[q]}
        if not reads: continue
        mult = {}
        for x in {x for v in reads.values() for x in v[0] if x[1] in 'ID'}:
            seq = F.fetch(c, max(0, x[0] - 20), x[0] + 21).upper(); i = min(20, x[0]); b = seq[i] if i < len(seq) else 'N'; u = 1; j = i + 1
            while j < len(seq) and seq[j] == b: u += 1; j += 1
            j = i - 1
            while j >= 0 and seq[j] == b: u += 1; j -= 1
            mult[x] = u
        units.append((k, ','.join(sorted(groups[k])), rep, reads, mult))
    return units

def signature(a):
    """what the collected reads depend on: the four inputs (size and time) and the read floor of the collection"""
    return TT.stamp(a.gtf, a.bam, a.tab, a.fasta) + [a.min_reads, 'v1']

def cache_ok(d, sig, need_qual):
    """a cached collection stands for this run only when it was made from the same inputs and, if the run reads base
    qualities, holds them: anything else is collected again (a stale cache would give another sample's sequences)"""
    return d.get('signature') == sig and (not need_qual or bool(d.get('qualities')))

def collect(a):
    """The costly half (two passes over the BAM): groups of tied loci and, per locus read, each tied read's differences,
    aligned blocks, intron chain and base qualities. With --cache PATH it is written there and reused by any later run on
    the same inputs, so another --qual / --cluster / --alpha takes seconds. -> [(group, loci, locus read at, reads, mult)], groups"""
    import pickle
    if a.cache and os.path.exists(a.cache):
        try:
            with open(a.cache, 'rb') as f: d = pickle.load(f)
            if cache_ok(d, signature(a), a.qual != 'none'):
                print(f'recurrence_test: reads taken from the cache {a.cache}', file=sys.stderr); return d['units'], d['groups']
            print(f'recurrence_test: cache {a.cache} is for other inputs; collecting again', file=sys.stderr)
        except Exception as e: print(f'recurrence_test: cache {a.cache} unreadable ({e}); collecting again', file=sys.stderr)
    meta, sets, groups, gr = TT.tie_groups(a.gtf, a.bam, a.tab, a.threads, a.ties)
    wq = a.qual != 'none' or bool(a.cache)
    _W.update(meta=meta, sets=sets, groups=groups, gr=gr, a=a, wq=wq)
    keys = sorted(groups, key=lambda k: (meta[k][0], meta[k][1])); big = sorted(keys, key=lambda k: -len(gr[k]))
    if a.threads > 1 and len(keys) > 1:
        with multiprocessing.get_context('fork').Pool(a.threads) as pool: got = dict(zip(big, pool.map(_collect_group, big, chunksize=1)))
    else: got = {k: _collect_group(k) for k in big}
    units = [u for k in keys for u in got[k]]
    if a.cache:
        with open(a.cache + '.tmp', 'wb') as f: pickle.dump({'signature': signature(a), 'qualities': wq, 'units': units, 'groups': len(groups)}, f, protocol=4)
        os.replace(a.cache + '.tmp', a.cache)
    return units, len(groups)

if __name__ == '__main__':
    main()

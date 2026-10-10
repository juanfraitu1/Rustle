#!/usr/bin/env python3
"""Allele-specific junctions as an O2 extension (docs/PREREG_asj_o2_2026-10-10.md).

Arm A: reads assigned to the primary's haplotype (S) or the other haplotype (B) at het sites taken from the haplotype assemblies, by O2's certificate
(`copy_pair_significance` / `poisson_binomial_upper_tail`, src/family.rs). Arm R: the June read-anchored split. Both: junction usage compared between the
groups (Fisher, BH, q < 0.05 and |dPSI| >= 0.3). The assemblies classify each arm-R anchor ALLELIC / NON-ALLELIC / UNRESOLVED.
Pure functions (tested in test_asj_o2.py) first; the per-chromosome driver (`run_chrom.py`) and the scorer (`score.py`) use them."""
import bisect
import functools
import math
import re

ERROR_RATE = 0.003            # O2 `error_rate` default (no base qualities)
EPS = ERROR_RATE / 3          # per distinguishing column, as in copy_pair_significance
ALPHA = 1e-3                  # O2 `alpha`; n = 2 haplotypes -> alpha / (n - 1)
INDEL_D = 5
MIN_ALT, MIN_DEPTH, MINOR_LO, MINOR_HI = 5, 12, 0.35, 0.65     # June anchor rule
MIN_JUNCTION_READS, Q_MAX, DPSI_MIN = 3, 0.05, 0.3              # June junction test
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def cs_walk(cs, strand, qs, qe, qoff):
    """PAF cs (query = a primary span starting at genomic `qoff`, target = the other haplotype) -> (snvs {primary pos: (primary base, B base)}
    in primary-strand bases, match runs [(start, end)] and indel positions, all in primary coordinates)"""
    snvs, runs, indels = {}, [], []
    plus = strand == "+"
    q = qs if plus else qe - 1
    for op, val in re.findall(r"([:*+\-~])([0-9]+|[a-z]+)", cs):
        if op == ":":
            n = int(val)
            runs.append((qoff + q, qoff + q + n) if plus else (qoff + q - n + 1, qoff + q + 1))
            q += n if plus else -n
        elif op == "*":
            tb, qb = val[0].upper(), val[1].upper()
            if not plus:
                tb, qb = tb.translate(COMP), qb.translate(COMP)
            snvs[qoff + q] = (qb, tb)
            q += 1 if plus else -1
        elif op == "+":
            indels.append(qoff + q if plus else qoff + q + 1)
            q += len(val) if plus else -len(val)
        elif op == "-":
            indels.append(qoff + q if plus else qoff + q + 1)
    return snvs, sorted(runs), sorted(indels)


def _near(pos, indels, d):
    i = bisect.bisect_left(indels, pos - d)
    return i < len(indels) and indels[i] <= pos + d


def het_sites(snvs, indels, d=INDEL_D):
    return {p: v for p, v in snvs.items() if not _near(p, indels, d)}


def read_obs(cig, rstart, seq, sites):
    """cigartuples (pysam codes) -> ({site: read base}, introns [(start, end)], aligned span, X events [(pos, base)]); sites sorted"""
    obs, introns, xs = {}, [], []
    r, q = rstart, 0
    for op, n in cig:
        if op in (0, 7, 8):
            i = bisect.bisect_left(sites, r)
            while i < len(sites) and sites[i] < r + n:
                obs[sites[i]] = seq[q + sites[i] - r]
                i += 1
            if op == 8:
                xs += [(r + j, seq[q + j]) for j in range(n)]
            r += n
            q += n
        elif op in (1, 4):
            q += n
        elif op == 2:
            r += n
        elif op == 3:
            introns.append((r, r + n))
            r += n
    return obs, introns, (rstart, r), xs


@functools.lru_cache(maxsize=None)
def _dp_tail(k, m, eps):
    """O2's poisson_binomial_upper_tail recursion with equal probabilities"""
    if k == 0:
        return 1.0
    if k > m:
        return 0.0
    dp = [1.0] + [0.0] * m
    for _ in range(m):
        for s in range(m, -1, -1):
            dp[s] = dp[s] * (1.0 - eps) + (dp[s - 1] * eps if s > 0 else 0.0)
    return sum(dp[k:])


@functools.lru_cache(maxsize=None)
def upper_tail(k, m, eps, exact_upto=60):
    """P(X >= k), X ~ Binomial(m, eps): O2's recursion up to `exact_upto` columns (so boundary cases such as k = m = 1 give O2's exact float);
    beyond, the closed-form sum in log space (the recursion is quadratic in m)"""
    if m <= exact_upto or k == 0 or k > m:
        return _dp_tail(k, m, eps)
    le, l1 = math.log(eps), math.log1p(-eps)
    lt = math.lgamma(m + 1) - math.lgamma(k + 1) - math.lgamma(m - k + 1) + k * le + (m - k) * l1
    t, s, x = math.exp(lt), 0.0, k
    while x <= m and t > 0.0:
        s += t
        if t < s * 1e-17:
            break
        t *= (m - x) / (x + 1) * eps / (1.0 - eps)
        x += 1
    return s


def certify(cols):
    """cols = [(read base, S allele, B allele)] at het sites the read covers -> ('S' | 'B' | None, p): strict MLE and p < alpha / (n - 1)"""
    ks = sum(b == s for b, s, _ in cols)
    kb = sum(b == o for b, _, o in cols)
    if ks == kb:
        return None, 1.0
    best, k = ("S", ks) if ks > kb else ("B", kb)
    p = upper_tail(k, len(cols), EPS)
    return (best if p < ALPHA else None), p


def junction_tests(reads, min_reads=MIN_JUNCTION_READS):
    """reads = [(group 0 | 1 | None, introns, (start, end))] -> per tested intron dict(intron, u0, s0, u1, s1, p, dpsi); a read using an intron
    spans it, so `used` counts the group's reads carrying it"""
    import collections
    import numpy as np
    from scipy.stats import fisher_exact
    carried, used = collections.Counter(), {0: collections.Counter(), 1: collections.Counter()}
    for g, ints, _ in reads:
        carried.update(set(ints))
        if g in used:
            used[g].update(set(ints))
    keep = sorted(j for j, n in carried.items() if n >= min_reads)
    if not keep:
        return []
    grp = np.array([-1 if g is None else g for g, _, _ in reads])
    st = np.array([s for _, _, (s, _) in reads])
    en = np.array([e for _, _, (_, e) in reads])
    m0, m1 = grp == 0, grp == 1
    out = []
    for j in keep:
        sp = (st <= j[0]) & (en >= j[1])
        s0, s1 = int(np.count_nonzero(sp & m0)), int(np.count_nonzero(sp & m1))
        u0, u1 = used[0][j], used[1][j]
        if not s0 or not s1 or (u0 == s0 and u1 == s1):
            continue
        p = fisher_exact([[u0, s0 - u0], [u1, s1 - u1]])[1]
        out.append(dict(intron=j, u0=u0, s0=s0, u1=u1, s1=s1, p=p, dpsi=u1 / s1 - u0 / s0))
    return out


def bh(ps):
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i])
    q, run = [0.0] * n, 1.0
    for r in range(n, 0, -1):
        i = order[r - 1]
        run = min(run, ps[i] * n / r)
        q[i] = run
    return q


def x_candidates(xs_per_read, min_alt=MIN_ALT):
    """X events of each read -> {pos: (alt base, reads)} for the top alternative base when it is on >= min_alt reads"""
    cnt = {}
    for xs in xs_per_read:
        for pos, b in set(xs):
            cnt.setdefault(pos, {}).setdefault(b, 0)
            cnt[pos][b] += 1
    out = {}
    for pos, d in cnt.items():
        b, n = max(d.items(), key=lambda kv: (kv[1], kv[0]))
        if n >= min_alt:
            out[pos] = (b, n)
    return out


def pick_anchor(cand):
    """cand = {pos: (ref reads, alt base, alt reads)} -> the anchor (highest depth, then smallest position) or None"""
    ok = []
    for pos, (rn, _, an) in cand.items():
        d = rn + an
        if rn >= MIN_ALT and an >= MIN_ALT and d >= MIN_DEPTH and MINOR_LO <= min(rn, an) / d <= MINOR_HI:
            ok.append((-d, pos))
    return min(ok)[1] if ok else None


def classify(pos, ref, alt, het, snvs, runs, indels, d=INDEL_D):
    if pos in het:
        return "ALLELIC" if {ref, alt} == set(het[pos]) else "UNRESOLVED"
    if pos in snvs or _near(pos, indels, d):
        return "UNRESOLVED"
    i = bisect.bisect_right(runs, (pos, float("inf"))) - 1
    return "NON-ALLELIC" if i >= 0 and runs[i][0] <= pos < runs[i][1] else "UNRESOLVED"


def concordance(pairs):
    """pairs = [(anchor allele 'ref' | 'alt', haplotype 'S' | 'B')] -> best share over the two pairings, None if empty"""
    if not pairs:
        return None
    a = sum((x, h) in (("ref", "S"), ("alt", "B")) for x, h in pairs)
    return max(a, len(pairs) - a) / len(pairs)


def halves(het):
    ps = sorted(het)
    return {p: het[p] for p in ps[0::2]}, {p: het[p] for p in ps[1::2]}


def gene_table(lines):
    """GFF3 lines -> [(chrom, start0, end, gene id, name, strand, most exons of one transcript)] in file order"""
    genes, parent_gene, exons = [], {}, {}
    for ln in lines:
        if ln.startswith("#") or not ln.strip():
            continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9:
            continue
        a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
        if f[2] == "gene":
            genes.append((f[0], int(f[3]) - 1, int(f[4]), a["ID"], a.get("Name", a["ID"]), f[6]))
        elif f[2] == "exon":
            p = a.get("Parent", "")
            exons[p] = exons.get(p, 0) + 1
        elif "Parent" in a and "ID" in a:
            parent_gene[a["ID"]] = a["Parent"]
    best = {}
    for p, n in exons.items():
        g = parent_gene.get(p, p)
        best[g] = max(best.get(g, 0), n)
    return [g + (best.get(g[3], 0),) for g in genes]


def process_gene(reads, het, snvs, runs, indels, pri_base):
    """reads = [(cigartuples, reference start, sequence, MAPQ)] of one gene -> both arms, T1 counts and the arm-R anchor with its referee class"""
    hs = sorted(het)
    h1, h2 = halves(het)
    recs, xs_all = [], []
    for cig, st, seq, mq in reads:
        obs, ints, span, xs = read_obs(cig, st, seq, hs)
        hap = certify([(obs[p],) + het[p] for p in hs if p in obs])[0]
        recs.append((obs, tuple(ints), span, mq, hap))
        xs_all.append(xs)
    both = disc = 0
    for obs, *_ in recs:
        a = certify([(obs[p],) + h1[p] for p in h1 if p in obs])[0]
        b = certify([(obs[p],) + h2[p] for p in h2 if p in obs])[0]
        if a and b:
            both += 1
            disc += a != b
    haps = [r[4] for r in recs]
    out = dict(n_reads=len(reads), n_het=len(het), armA=dict(S=haps.count("S"), B=haps.count("B"), abstain=haps.count(None)),
               testsA=junction_tests([({"S": 0, "B": 1}.get(r[4]), r[1], r[2]) for r in recs]), T1=dict(both=both, discordant=disc),
               anchor=None, testsR=[])
    cand = x_candidates(xs_all)
    if not cand:
        return out
    cs = sorted(cand)
    bases = [read_obs(cig, st, seq, cs)[0] for cig, st, seq, _ in reads]
    full = {}
    for pos in cs:
        ref, alt = pri_base(pos), cand[pos][0]
        full[pos] = (sum(b.get(pos) == ref for b in bases), alt, sum(b.get(pos) == alt for b in bases))
    a = pick_anchor(full)
    if a is None:
        return out
    ref, alt = pri_base(a), full[a][1]
    grp = [0 if b.get(a) == ref else 1 if b.get(a) == alt else None for b in bases]
    out["testsR"] = junction_tests([(g, r[1], r[2]) for g, r in zip(grp, recs)])
    others = [p for p in hs if p != a]
    pairs = []
    for g, r in zip(grp, recs):
        if g is not None:
            hap = certify([(r[0][p],) + het[p] for p in others if p in r[0]])[0]
            if hap:
                pairs.append(("ref" if g == 0 else "alt", hap))

    def mq0(gg):
        m = [r[3] for g, r in zip(grp, recs) if g == gg]
        return sum(x == 0 for x in m) / len(m) if m else None
    out["anchor"] = {"pos": a, "ref": ref, "alt": alt, "ref_n": full[a][0], "alt_n": full[a][2], "class": classify(a, ref, alt, het, snvs, runs, indels),
                     "transition": {ref, alt} in ({"A", "G"}, {"C", "T"}), "conc": concordance(pairs), "n_conc": len(pairs),
                     "mapq0_ref": mq0(0), "mapq0_alt": mq0(1)}
    return out


def allelic_alignments(lines, min_cov=0.9):
    """PAF of primary spans on the other haplotype -> {query: fields} for the primary (tp:A:P) record with the most matching bases, kept when it
    covers >= min_cov of the span"""
    best = {}
    for ln in lines:
        f = ln.rstrip("\n").split("\t")
        if "tp:A:P" not in f[12:]:
            continue
        if f[0] not in best or int(f[9]) > int(best[f[0]][9]):
            best[f[0]] = f
    return {q: f for q, f in best.items() if (int(f[3]) - int(f[2])) / int(f[1]) >= min_cov}

import collections
import importlib.util
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("o3u_psvclust", os.path.join(os.path.dirname(HERE), "o3_unexplained", "psvclust.py"))
P = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(P)


# the pre-2026-10-16 implementations, kept as the reference the faster ones must reproduce exactly
def _old_linkage(sites_u, vec):
    """median over pairs of unexplained alt sites of P(alt at both | alt at either), among reads covering both; None if no pair has
    >= P.MIN_SUPPORT reads with an alt at either"""
    import statistics
    vals = []
    su = sorted(sites_u)
    for i in range(len(su)):
        for j in range(i + 1, len(su)):
            a, b = su[i], su[j]
            either = both = 0
            for v in vec.values():
                if a in v and b in v:
                    xa, xb = v[a] == "alt", v[b] == "alt"
                    if xa or xb:
                        either += 1
                        both += xa and xb
            if either >= P.MIN_SUPPORT:
                vals.append(both / either)
    return statistics.median(vals) if vals else None


def _old_cluster_host(ref, reads_cov):
    """ref: the host locus; reads_cov: {read: {pos: base}} -> (sites, clusters as lists of reads)"""
    alt = collections.defaultdict(collections.Counter)
    for cov in reads_cov.values():
        for p, b in cov.items():
            if b in "ACGT" and b != ref[p]:
                alt[p][b] += 1
    sites = {p: c.most_common(1)[0][0] for p, c in alt.items() if c.most_common(1)[0][1] >= P.MIN_SUPPORT}
    vec = {r: {p: ("alt" if cov[p] == sites[p] else "ref" if cov[p] == ref[p] else "other") for p in sites if p in cov}
           for r, cov in reads_cov.items()}
    order = sorted(vec, key=lambda r: (-len(vec[r]), r))
    clusters = []        # [members, {site: Counter(allele)}]

    def cons(cl):
        return {p: c.most_common(1)[0][0] for p, c in cl[1].items()}

    def fit(v, cl):
        cc = cons(cl)
        co = [p for p in v if p in cc and v[p] != "other" and cc[p] != "other"]
        conf = sum(v[p] != cc[p] for p in co)
        return len(co), conf

    for r in order:
        v = vec[r]
        best = None
        for i, cl in enumerate(clusters):
            co, conf = fit(v, cl)
            if conf <= max(0, (co - 9) // 10) and (best is None or co > best[1]):
                best = (i, co)
        if best is None:
            clusters.append([[r], collections.defaultdict(collections.Counter)])
            i = len(clusters) - 1
        else:
            i = best[0]
            clusters[i][0].append(r)
        for p, a in v.items():
            clusters[i][1][p][a] += 1

    def separated(a, b):
        """a site where each cluster's consensus allele differs, supported by >= P.MIN_SUPPORT reads on each side"""
        ca, cb = cons(a), cons(b)
        return any(ca[p] != cb[p] and a[1][p][ca[p]] >= P.MIN_SUPPORT and b[1][p][cb[p]] >= P.MIN_SUPPORT for p in ca if p in cb)

    changed = True
    while changed:
        changed = False
        big = [c for c in clusters if len(c[0]) >= P.MIN_SUPPORT]
        small = [c for c in clusters if len(c[0]) < P.MIN_SUPPORT]
        for c in small:                                  # dissolve small clusters into their best compatible big cluster
            for r in c[0]:
                v = vec[r]
                fits = [(fit(v, b), k) for k, b in enumerate(big)]
                fits = [(co, k) for (co, conf), k in fits if conf <= max(0, (co - 9) // 10)]
                if fits:
                    k = max(fits)[1]
                    big[k][0].append(r)
                    for p, a in v.items():
                        big[k][1][p][a] += 1
        clusters = big
        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                if not separated(clusters[i], clusters[j]):
                    for p, c in clusters[j][1].items():
                        clusters[i][1][p].update(c)
                    clusters[i][0] += clusters[j][0]
                    del clusters[j]
                    changed = True
                    break
            if changed:
                break
    return sites, [(c[0], cons(c)) for c in clusters], vec


def _random_host(rng, n_reads=40, length=400):
    ref = "".join(rng.choice("ACGT") for _ in range(length))
    copies = []
    for _ in range(rng.randint(1, 4)):
        c = list(ref)
        for p in rng.sample(range(length), rng.randint(2, 12)):
            c[p] = rng.choice([b for b in "ACGT" if b != c[p]])
        copies.append(c)
    reads = {}
    for i in range(n_reads):
        c = rng.choice(copies)
        s = rng.randint(0, length - 50)
        e = rng.randint(s + 30, length)
        cov = {}
        for p in range(s, e):
            if rng.random() < 0.01:
                continue                                          # a deletion: no base
            cov[p] = c[p] if rng.random() > 0.01 else rng.choice("ACGTN")
        reads[f"r{i:03d}"] = cov
    return ref, reads


def test_psv_clustering_reproduces_the_reference_implementation_exactly():
    rng = random.Random(5)
    for k in range(60):
        ref, reads = _random_host(rng, n_reads=rng.randint(3, 60))
        assert P.cluster_host(ref, reads) == _old_cluster_host(ref, reads), k


def test_linkage_reproduces_the_reference_implementation_exactly():
    rng = random.Random(8)
    for k in range(60):
        ref, reads = _random_host(rng, n_reads=rng.randint(3, 60))
        sites, _, vec = _old_cluster_host(ref, reads)
        su = sorted(sites)
        for sub in (su, su[: len(su) // 2], rng.sample(su, min(len(su), 5))):
            assert P.linkage(sub, vec) == _old_linkage(sub, vec), k

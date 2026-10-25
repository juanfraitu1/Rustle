import collections
exec(open(__file__.replace('_robust', '')).read().split("fams = sorted")[0].replace("D = sys.argv[1] if len(sys.argv) > 1 else '.'", "D='.'"))
fams = sorted(n_copies, key=lambda x: int(x[5:]))
def loci2(a, fam, T, C):
    iv = collections.defaultdict(list)
    for fm, q, cov, idn, t, ts, te in hits[a]:
        if fm == fam and cov >= C and idn >= T and not sex(a, t): iv[t].append((ts, te))
    n = 0
    for t, xs in iv.items():
        xs.sort(); end = -1
        for s, e in xs:
            if s > end + 1000: n += 1
            end = max(end, e)
    return n
sets = {}
for C in (0.5, 0.8, 0.95):
    for T in (0.95, 0.98):
        d = {fm: (loci2('ref', fm, T, C), loci2('mat', fm, T, C), loci2('pat', fm, T, C)) for fm in fams}
        S = {fm for fm, (r, m, p) in d.items() if m != p}
        sets[(C, T)] = S
        refside = collections.Counter('ref=mat' if r == m else 'ref=pat' if r == p else 'ref=neither' for fm, (r, m, p) in d.items() if fm in S)
        mag = collections.Counter(min(abs(m - p), 3) for fm, (r, m, p) in d.items() if fm in S)
        tot = sum(m + p for r, m, p in d.values()); dif = sum(abs(m - p) for r, m, p in d.values())
        print(f'qcov>={C} id>={T}: mat!=pat {len(S)}/494 families; {dict(refside)}; |mat-pat| (3=3+) {dict(sorted(mag.items()))}; copy slots differing {dif}/{tot}')
core = set.intersection(*sets.values())
print('stable across all 6 settings:', len(core), sorted(core, key=lambda x: int(x[5:])))
# the non-reference haplotype's version of each autosomal copy
ref_t = {}
for fm, q, cov, idn, t, *_ in hits['ref']:
    if idn > 0.9999 and cov > 0.999: ref_t[q] = t
best = {a: collections.defaultdict(float) for a in ('mat', 'pat')}
for a in best:
    for fm, q, cov, idn, *_ in hits[a]:
        if cov >= 0.8: best[a][q] = max(best[a][q], idn)
c = collections.Counter(); diffs = []
for q, t in ref_t.items():
    if SEX['ref'].get(t): continue
    m, p = best['mat'][q], best['pat'][q]
    other = min(m, p)
    c['identical in both' if other >= 0.99995 else 'other hap >=99.9%' if other >= 0.999 else 'other hap 99-99.9%' if other >= 0.99 else 'other hap 90-99%' if other >= 0.9 else 'other hap <90%/absent'] += 1
print('autosomal reference copies, best match on the OTHER haplotype:', dict(c), 'n', sum(c.values()))

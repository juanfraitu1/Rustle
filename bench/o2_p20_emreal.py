# real fibroblast reads: for exactly tied reads, do the tied loci differ in UNIQUE-read support (the EM's only evidence)?
import collections
BAND = 10
best = {}   # read -> list of (bin, ms, NM)
for g in ("grp00", "grp01", "grp02"):
    for ln in open(f"{g}.tsv"):
        f = ln.rstrip("\n").split("\t"); d = {t[:2]: int(t[5:]) for t in f[5:]}
        best.setdefault(f[0], []).append((f"{f[2]}:{int(f[3]) // 20000}", d["ms"], d["NM"]))
uniq = collections.Counter(); ties = []
for r, rs in best.items():
    m = max(x[1] for x in rs); band = [x for x in rs if x[1] >= m - BAND]
    bins = {x[0] for x in band}
    if len(bins) == 1: uniq[next(iter(bins))] += 1; continue
    top = [x for x in band if x[1] == m]; nm = min(x[2] for x in top)
    tb = sorted({x[0] for x in top if x[2] == nm})
    if len(tb) >= 2: ties.append(tb)
c = collections.Counter(); share = collections.Counter()
for tb in ties:
    u = sorted((uniq[b] for b in tb), reverse=True)
    k = 'no tied locus has unique reads' if u[0] == 0 else 'exactly one tied locus has unique reads' if u[1] == 0 else 'several tied loci have unique reads'
    c[k] += 1
    if u[0] > 0:
        f = u[0] / sum(u); share['top locus >= 99% of unique reads' if f >= 0.99 else 'top locus 90-99%' if f >= 0.9 else 'top locus < 90%'] += 1
print(f'reads {len(best):,}; reads unique to one 20-kb locus {sum(uniq.values()):,}; exactly tied reads (same ms and NM, >= 2 loci) {len(ties):,}')
for k, v in c.items(): print(f'  {k}: {v:,} ({v / len(ties):.1%})')
for k, v in share.items(): print(f'  {k}: {v:,}')

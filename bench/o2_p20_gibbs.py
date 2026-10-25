# P20: collapsed Gibbs over read -> copy assignments, Dirichlet(ALPHA) prior on abundances; posterior = assignment frequency
import sys, random, collections
SPLIT = sys.argv[1]; PFX = sys.argv[2] if len(sys.argv) > 2 else 'p20'; ALPHA = 1.0; BURN, SWEEPS = 200, 2000; POST = 0.99
sys.argv = ['x', SPLIT]
src = open('em19.py').read().replace('p19.', f'{PFX}.').split("units = sorted")[0]
exec(src)
rng = random.Random(f'gibbs:{PFX}:{SPLIT}')
names = sorted(reads); cand = {r: list(reads[r][0].items()) for r in names}
amb = [r for r in names if len(cand[r]) > 1]
n = collections.Counter(); z = {}
for r in names:
    u = max(cand[r], key=lambda x: x[1])[0]; z[r] = u; n[u] += 1
freq = {r: collections.Counter() for r in amb}
for s in range(BURN + SWEEPS):
    for r in amb:
        n[z[r]] -= 1
        w = [(n[u] + ALPHA) * v for u, v in cand[r]]; t = rng.random() * sum(w); acc = 0
        for (u, _), x in zip(cand[r], w):
            acc += x
            if acc >= t: break
        z[r] = u; n[u] += 1
        if s >= BURN: freq[r][u] += 1
c = collections.Counter(); out = []
for r in amb:
    L, pu, tie, tu = reads[r]
    b, k = freq[r].most_common(1)[0]; post = k / SWEEPS
    g = 'LL-tied reads' if tie else 'decidable band reads'
    c[(g, 'reads')] += 1; c[(g, 'primary right')] += pu == tu
    if post >= POST:
        c[(g, 'assigned')] += 1; c[(g, 'assigned wrong')] += b != tu
        c[(g, 'assigned wrong | primary right')] += b != tu and pu == tu; c[(g, 'assigned right | primary wrong')] += b == tu and pu != tu
    else: c[(g, 'abstains')] += 1
    out.append((r, tie, pu == tu, b == tu, round(post, 4)))
for k in sorted(c): print(k, c[k])
open(f'{PFX}.{SPLIT}.gibbs.tsv', 'w').write('read\ttie\tprimary_right\tpick_right\tposterior\n' + '\n'.join('\t'.join(map(str, x)) for x in out) + '\n')
# calibration: posterior bins vs observed right rate (all ambiguous reads)
bins = collections.defaultdict(lambda: [0, 0])
for r, tie, pr, br, p in out:
    k = '>=0.99' if p >= 0.99 else '0.9-0.99' if p >= 0.9 else '0.7-0.9' if p >= 0.7 else '<0.7'
    bins[k][0] += 1; bins[k][1] += br
for k in ('<0.7', '0.7-0.9', '0.9-0.99', '>=0.99'):
    if bins[k][0]: print(f'posterior {k}: {bins[k][0]} reads, right {bins[k][1] / bins[k][0]:.3f}')

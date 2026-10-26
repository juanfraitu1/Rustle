"""WGSCHK scorer (LEDGER_O3 prereg b7894f54). usage: wgschk_score.py WGSCHK_DIR SDCN_DIR"""
import collections, glob, sys
import numpy as np
W, S = sys.argv[1], sys.argv[2]
a = {g: np.fromfile(f"{W}/a_{g}.i32", np.int32) for g in ("ref", "mat", "pat")}
w = sum(np.fromfile(f, np.int32).astype(np.int64) for f in sorted(glob.glob(f"{W}/w_*.i32")))
z = np.load(f"{W}/unit_kmers.npz"); idx, off, names = z["idx"], z["off"], list(z["names"])
one = (a["ref"] == 1) & (a["mat"] == 1) & (a["pat"] == 1)
v = w[one]; h = np.bincount(np.minimum(v, 1000)); h[:2] = 0; mode = int(np.argmax(h))
lam = float(v[v < 6 * mode].mean()) / 2
print(f"lambda_hap {lam:.2f} (single-copy k-mers {one.sum()}, mode {mode})")
# 1416 headline counts (C 0.8, T 0.90), same rules as bench/sdcn_count2.py
SEXREF = {'NC_073247.2', 'NC_073248.2', 'NC_011120.1'}
def sex(g, t): return t in SEXREF if g == 'ref' else t.split('_')[0] in ('chrX', 'chrY', 'chrM')
def count(paf, kind):
    H = collections.defaultdict(lambda: collections.defaultdict(list))
    for l in open(paf):
        f = l.split('\t')
        if f[0].startswith(kind) and (int(f[3]) - int(f[2])) / int(f[1]) >= 0.8 and int(f[9]) / int(f[10]) >= 0.90 and not sex(g, f[5]):
            H[f[0]][f[5]].append((int(f[7]), int(f[8])))
    n = {}
    for u, d in H.items():
        k = 0
        for t, xs in d.items():
            end = -1
            for s0, e0 in sorted(xs):
                if s0 >= end: k += 1
                end = max(end, e0)
        n[u] = k
    return n
N = {}
for kind, files in (("EX", "ex_sr_{}.paf"), ("SD", "{}.paf")):
    for g in ("ref", "mat", "pat"):
        for u, k in count(f"{S}/" + files.format(g), kind).items():
            N.setdefault(u, [0, 0, 0])[("ref", "mat", "pat").index(g)] = k
for kind in ("EX", "SD"):
    res = collections.Counter(); ctrl = []
    for i, u in enumerate(names):
        if not u.startswith(kind): continue
        k = idx[off[i]:off[i + 1]]
        if k.size == 0: continue
        obs = float(np.median(w[k])) / lam; ha = float(np.median(a["mat"][k] + a["pat"][k])); hr = float(np.median(2 * a["ref"][k]))
        n = N.get(u, [0, 0, 0])
        if n[1] != n[2]:
            if ha == hr: res["tested-out (H_asm == H_ref)"] += 1
            elif abs(obs - ha) < abs(obs - hr): res["supports assemblies"] += 1
            elif abs(obs - ha) > abs(obs - hr): res["supports reference"] += 1
            else: res["undecided"] += 1
        else:
            ctrl.append(abs(obs - ha))
    d = res["supports assemblies"] + res["supports reference"]
    print(f"{kind}: {dict(res)}; assemblies supported in {100 * res['supports assemblies'] / max(1, d):.1f}% of decided")
    c = np.array(ctrl); print(f"  control (no mat/pat difference, n={c.size}): |obs - H_asm| median {np.median(c):.2f}, <0.5 in {100 * (c < 0.5).mean():.1f}%, <1 in {100 * (c < 1).mean():.1f}%")

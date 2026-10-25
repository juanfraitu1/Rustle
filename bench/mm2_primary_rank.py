import collections
# which score picks minimap2's primary? per read with >= 1 secondary: primary vs best secondary on ms, AS, NM
prim, best = {}, {}
for g in ("grp00", "grp01", "grp02"):
    for ln, li in zip(open(f"{g}.tsv"), open(f"{g}.intr")):
        f = ln.rstrip("\n").split("\t"); ni = int(li.rsplit("\t", 1)[1])
        d = {t[:2]: int(t[5:]) for t in f[5:]}
        rec = (d["ms"], d["AS"], d["NM"], ni)
        if int(f[1]) & 256:
            b = best.get(f[0])
            best[f[0]] = rec if b is None else (max(b[0], rec[0]), max(b[1], rec[1]), min(b[2], rec[2]), b[3])
            if b is None or rec[0] > b[0]: best[f[0] + "#ms"] = rec
        else: prim[f[0]] = rec
c = collections.Counter(); n = 0
for r, p in prim.items():
    b = best.get(r)
    if b is None: continue
    n += 1
    for i, k in ((0, "ms"), (1, "AS")):
        c[f"{k}: primary {'>' if p[i] > b[i] else '=' if p[i] == b[i] else '<'} best secondary"] += 1
    c[f"NM: primary {'<' if p[2] < b[2] else '=' if p[2] == b[2] else '>'} best secondary"] += 1
    s = best[r + "#ms"]
    if p[0] == s[0]:
        c["ms tie | same NM" if p[2] == s[2] else "ms tie | primary fewer edits" if p[2] < s[2] else "ms tie | primary MORE edits"] += 1
        c["ms tie | same intron count" if p[3] == s[3] else "ms tie | different intron count"] += 1
print(f"reads with >= 1 secondary: {n:,}")
for k in sorted(c): print(f"  {k}: {c[k]:,} ({c[k]/n:.1%})")

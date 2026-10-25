# P15 (LEDGER_O2, register 1404): score sim.bam against truth.tsv for windows A / B / C and the advisor rule R; run in /mnt/linuxdisk/tmp/p15
import collections, re, subprocess, json
PEN = 32
truth = {r[0]: (r[3], int(r[4]), int(r[5])) for r in (l.rstrip("\n").split("\t") for l in list(open("truth.tsv"))[1:])}
sam = subprocess.run("samtools view -F 0x804 sim.bam", shell=True, capture_output=True, text=True, check=True).stdout
prim, secs = {}, collections.defaultdict(list)
for ln in sam.splitlines():
    f = ln.split("\t")
    rl = sum(int(n) for n, o in re.findall(r"(\d+)([MDN=X])", f[5]))
    ni = f[5].count("N")
    d = {t[:2]: int(t[5:]) for t in f[11:] if t[:2] in ("AS", "NM") and t[3] == "i"}
    c, s, e = truth[f[0]]
    right = f[2] == c and int(f[3]) - 1 < e and int(f[3]) - 1 + rl > s
    rec = (d["AS"], d["NM"], ni, int(f[4]), right)
    if int(f[1]) & 256: secs[f[0]].append(rec)
    else: prim[f[0]] = rec
W = {"A": lambda p, s: s[0] >= p[0] - 10,
     "B": lambda p, s: s[0] + PEN * s[2] >= p[0] + PEN * p[2] - 10,
     "C": lambda p, s: s[2] == p[2] and abs(s[0] - p[0]) <= 10}
c = collections.Counter()
c["reads"] = len(prim); c["primary right (all reads)"] = sum(p[4] for p in prim.values())
for n, p in prim.items():
    ss = secs.get(n, [])
    if ss: c["reads with a secondary"] += 1
    inw = {k: [s for s in ss if f(p, s)] for k, f in W.items()}
    for k, v in inw.items():
        if v:
            c[f"{k}: reads in window"] += 1; c[f"{k}: primary right"] += p[4]
            c[f"{k}: recoverable (primary wrong, a right secondary in window)"] += (not p[4]) and any(s[4] for s in v)
    if inw["A"] and not inw["C"]:
        c["A not C: reads"] += 1; c["A not C: primary right"] += p[4]
    if (not p[4]) and any(s[4] for s in inw["A"]):
        c["H3 base: recoverable in A"] += 1
        c["H3: ... C holds a right secondary"] += any(s[4] for s in inw["C"])
        c["   ... B holds a right secondary"] += any(s[4] for s in inw["B"])
    R = [s for s in ss if s[0] >= p[0] + 10]
    if R:
        b = max(R, key=lambda s: (s[0], -s[1]))
        c["R fires: reads"] += 1; c["R: primary right"] += p[4]; c["R: chosen secondary right"] += b[4]
        c["R: chosen secondary has fewer introns"] += b[2] < p[2]
for k, v in c.items(): print(f"{k}: {v:,}")
pct = lambda a, b: c[a] / c[b] if c[b] else float("nan")
h1 = (pct("R: primary right", "R fires: reads"), pct("R: chosen secondary right", "R fires: reads"))
h2 = pct("A not C: primary right", "A not C: reads")
h3 = pct("H3: ... C holds a right secondary", "H3 base: recoverable in A")
print(f"H1 primary right {h1[0]:.3f} (>= .90) / R's secondary right {h1[1]:.3f} (<= .10): {'PASS' if h1[0] >= .9 and h1[1] <= .1 else 'FAIL'}")
print(f"H2 A-not-C primary right {h2:.4f} (>= .98): {'PASS' if h2 >= .98 else 'FAIL'}")
print(f"H3 C keeps {h3:.3f} of A's recoverable (>= .90): {'PASS' if h3 >= .9 else 'FAIL'}")

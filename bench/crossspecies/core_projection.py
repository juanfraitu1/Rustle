#!/usr/bin/env python3
"""A DIFFERENT LEVEL: represent a locus by its PROJECTION ONTO THE FAMILY CORE.

Levels tried so far, and why this one:
  genomic span        -- best representative measured (family recall 0.929)
  read-exon blocks    -- better PROBE, worse REP; max() made it worse, so the denominator
                         explanation is refuted
  k-mer Jaccard       -- already dead
  junction co-thread  -- AUC 0.88 but undefined for 79% of pairs; cannot stand alone

Three of today's results point at the core instead of at pairwise sequence:
  * 27 NPIP loci collapse to ONE 16,118 bp minigraph path -- the family IS a cassette
  * the cassette is size-invariant while annotated spans vary 4.65x
  * true PPIAL4 members cover 98-100% of the seed; processed pseudogenes cover 71-85%

So: build the family core once, project every locus onto it, and relate loci by WHICH PART OF THE
CORE THEY CARRY rather than by aligning them to each other.

  core           = minigraph collapsed path over the family's loci (ref = longest locus)
  node           = the set of CORE INTERVALS a locus covers
  similarity     = Jaccard over core positions
  edge           = Jaccard >= J_FLOOR

WHY THIS MIGHT BEAT PAIRWISE COVERAGE: `min(|u|,|v|)` is length-dependent, which is what let 596 bp
processed pseudogenes reach coverage 0.91 against 760 bp real members. Core positions are a FIXED
frame -- a fragment covers a small share of the core no matter how short it is, so shortness stops
being an advantage.

⚠ Scored at the FAMILY level vs Soto curated membership -- this changes what a node IS.
⚠ J_FLOOR is swept, not chosen. Any single value would be a tuned threshold.
⚠ Precision is a LOWER BOUND (PPIAL4G/H are real members absent from Soto's list).
"""
import os
import subprocess
import sys
from collections import defaultdict

W = "/home/juanfra/winloci_scratch/seedfam/rep2x2"
C = "/home/juanfra/winloci_scratch/seedfam/core"
HS = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0.fa"
TRUTH = "/mnt/c/Users/jfris/Desktop/Rustle/bench/soto/80_fams.gene_preferred.bed"
T = "4"
FAMS = sys.argv[1].split(",") if len(sys.argv) > 1 else [
    "ID_431", "ID_131", "ID_391", "ID_395", "ID_400", "ID_116", "ID_8", "ID_468"]
J_FLOORS = (0.30, 0.50, 0.70)
os.makedirs(C, exist_ok=True)

members = defaultdict(list)
for line in open(TRUTH):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4 and "|" in f[3]:
        g, fid = f[3].split("|", 1)
        if fid in FAMS:
            members[fid].append((f[0], int(f[1]), int(f[2]), g))


def sh(cmd, out=None):
    if out:
        with open(out, "w") as fh:
            return subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL).returncode
    return subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode


def comps(ns, E):
    p = {x: x for x in ns}

    def fi(x):
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x
    for e in E:
        a, b = tuple(e)
        ra, rb = fi(a), fi(b)
        if ra != rb:
            p[ra] = rb
    g = defaultdict(list)
    for x in ns:
        g[fi(x)].append(x)
    return sorted(g.values(), key=len, reverse=True)


def curated_in(fid, comp):
    s = set()
    for reg in comp:
        c, rng = reg.rsplit(":", 1)
        a, b = rng.split("-")
        a, b = int(a) - 1, int(b)
        for cc, gs, ge, g in members[fid]:
            if cc == c and gs < b and a < ge:
                s.add(g)
    return s


rows = {}
for fid in FAMS:
    fa = f"{W}/{fid}_dna.fa"
    if not os.path.exists(fa):
        continue
    seqs, name, buf = [], None, []
    for l in open(fa):
        if l[0] == ">":
            if name:
                seqs.append((name, "".join(buf)))
            name, buf = l[1:].strip(), []
        else:
            buf.append(l.strip())
    if name:
        seqs.append((name, "".join(buf)))
    if len(seqs) < 2:
        continue
    seqs.sort(key=lambda x: -len(x[1]))          # longest locus is the core reference
    d = f"{C}/{fid}"
    os.makedirs(d, exist_ok=True)
    for i, (n, s) in enumerate(seqs):
        open(f"{d}/{i:03d}.fa", "w").write(f">{n}\n{s}\n")
    parts = [f"{d}/{i:03d}.fa" for i in range(len(seqs))]
    sh(["minigraph", "-cxggs", "-t", T] + parts, out=f"{C}/{fid}.gfa")
    core_bp = 0
    with open(f"{C}/{fid}_core.fa", "w") as fh:
        segs = []
        for line in open(f"{C}/{fid}.gfa"):
            if line.startswith("S"):
                p = line.split("\t")
                segs.append(p[2].strip())
        core = "".join(segs)
        core_bp = len(core)
        fh.write(">core\n" + core + "\n")
    if core_bp < 200:
        continue

    # project every locus onto the core
    sh(["minimap2", "-x", "asm20", "-c", "--eqx", "-N", "50", "-p", "0.02", "-t", T,
        f"{C}/{fid}_core.fa", fa], out=f"{C}/{fid}_proj.paf")
    cov = defaultdict(list)
    for line in open(f"{C}/{fid}_proj.paf"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 12:
            continue
        q, ts, te, nm, bl = f[0], int(f[7]), int(f[8]), int(f[9]), int(f[10])
        if bl >= 200 and nm / bl >= 0.80:
            cov[q].append((ts, te))

    def merged(iv):
        iv = sorted(iv)
        out = []
        cs = ce = None
        for s, e in iv:
            if cs is None:
                cs, ce = s, e
            elif s > ce:
                out.append((cs, ce))
                cs, ce = s, e
            else:
                ce = max(ce, e)
        if cs is not None:
            out.append((cs, ce))
        return out

    proj = {n: merged(cov.get(n, [])) for n, _ in seqs}
    ns = [n for n, _ in seqs]

    def ilen(a, b):
        t = 0
        for s1, e1 in a:
            for s2, e2 in b:
                t += max(0, min(e1, e2) - max(s1, s2))
        return t

    def tlen(a):
        return sum(e - s for s, e in a)
    share = {n: tlen(proj[n]) / core_bp for n in ns}
    covered = [n for n in ns if tlen(proj[n]) > 0]
    rows[fid] = {}
    for J in J_FLOORS:
        E = set()
        for i in range(len(covered)):
            for j in range(i + 1, len(covered)):
                a, b = covered[i], covered[j]
                inter = ilen(proj[a], proj[b])
                uni = tlen(proj[a]) + tlen(proj[b]) - inter
                if uni and inter / uni >= J:
                    E.add(frozenset((a, b)))
        cs = comps(ns, E)
        best = max(cs, key=lambda comp: len(curated_in(fid, comp)))
        got = curated_in(fid, best)
        rows[fid][J] = (len(got) / len(members[fid]), len(got) / len(best), len(best))
    sh_true = [share[n] for n in ns if any(
        cc == n.rsplit(":", 1)[0] and gs < int(n.rsplit(":", 1)[1].split("-")[1])
        and int(n.rsplit(":", 1)[1].split("-")[0]) - 1 < ge
        for cc, gs, ge, g in members[fid])]
    sh_oth = [share[n] for n in ns if n not in
              [x for x in ns if any(
                  cc == x.rsplit(":", 1)[0] and gs < int(x.rsplit(":", 1)[1].split("-")[1])
                  and int(x.rsplit(":", 1)[1].split("-")[0]) - 1 < ge
                  for cc, gs, ge, g in members[fid])]]
    sh_true.sort()
    sh_oth.sort()
    print(f"{fid}: core {core_bp:,} bp from {len(seqs)} loci | "
          f"core-share curated median {sh_true[len(sh_true)//2] if sh_true else 0:.2f}, "
          f"other median {sh_oth[len(sh_oth)//2] if sh_oth else 0:.2f}", flush=True)

print("\n" + "=" * 84)
print("CORE-PROJECTION representation, family-level vs Soto curated membership")
print("=" * 84)
print(f"{'family':<9}" + "".join(f"{'J>='+str(J):>22}" for J in J_FLOORS))
agg = defaultdict(list)
for fid in FAMS:
    if fid not in rows:
        continue
    line = f"{fid:<9}"
    for J in J_FLOORS:
        rc, pr, sz = rows[fid][J]
        agg[J].append((rc, pr))
        line += f"{f'{rc:.2f} / {pr:.2f} (n={sz})':>22}"
    print(line)
print()
for J in J_FLOORS:
    v = agg[J]
    if not v:
        continue
    mr = sum(x[0] for x in v) / len(v)
    mp = sum(x[1] for x in v) / len(v)
    print(f"  J>={J}: mean recall {mr:.3f}  mean precision* {mp:.3f}  "
          f"F1 {2*mr*mp/(mr+mp) if mr+mp else 0:.3f}  ({len(v)} families)")
print("\n  compare (same 8 families, pairwise sequence): dna/min recall 0.929 prec* 0.129 F1 0.227")
print("  *precision is a LOWER BOUND.")
print("DONE")

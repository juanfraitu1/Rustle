#!/usr/bin/env python3
"""Where does the 0.54x size deficit come from -- the reads, or the representative we pick from them?

Usage: size_source.py <copies.tsv> <bam> [chr1,chr15]

The bipartite matcher established that predicted copies are about HALF the size of their matched truth
locus (median 0.54, truncated 104 vs over-extended 12). That is a paired per-locus comparison, so it is
not an artifact of families being size-heterogeneous (they are: Soto's own median smallest/largest is 0.33).

This asks whether the deficit is already present in the DATA or is introduced by REPRESENTATIVE SELECTION.
`refine_copy_seq` takes the span from the chosen representative's start/end, and `pick_locus_rep` chooses
ONE intron chain out of many (NPIPB9: 206 reads -> 54 chains, largest 21%). Measured separately, stub reps
average 4,746 bp against spliced reps' 13,065. So the representative may be a fraction of the locus.

Three sizes per detected truth member, all divided by the truth size:
  rep      -- what the pipeline emits today (the chosen representative's span)
  reads    -- union extent of PRIMARY reads overlapping that copy (-F 2308, per the standing rule)
  reads3   -- same, but only reads sharing the locus with >= 3 read support at their own start/end,
              which is the cheap stand-in for a locus GROUP rather than a readthrough pile-up

If `reads` >> `rep`, the size deficit is REPRESENTATIVE SELECTION and is recoverable without new data.
If `reads` ~ `rep`, the reads themselves never covered the locus and no boundary rule can fix it.
"""
import subprocess
import sys
from collections import defaultdict

COP, BAM = sys.argv[1], sys.argv[2]
CHROMS = set(sys.argv[3].split(",")) if len(sys.argv) > 3 else None
BED = "bench/soto/80_fams.gene_preferred.bed"

copies = []
for i, ln in enumerate(open(COP)):
    if i == 0:
        continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 7:
        continue
    if CHROMS and f[3] not in CHROMS:
        continue
    copies.append((f[3], int(f[4]), int(f[5]), int(f[6])))

members = []
for ln in open(BED):
    c, s, e, name, *_ = ln.rstrip("\n").split("\t")
    if CHROMS and c not in CHROMS:
        continue
    members.append((c, int(s), int(e), name.split("|")[0]))

rows = []
for (c, ms, me, gene) in members:
    ov = [(min(me, pe) - max(ms, ps), ps, pe, nx)
          for (pc, ps, pe, nx) in copies if pc == c and ms < pe and me > ps]
    if not ov:
        continue
    _, ps, pe, nx = max(ov)
    out = subprocess.run(["samtools", "view", "-F", "2308", BAM, f"{c}:{ps}-{pe}"],
                         capture_output=True, text=True).stdout
    lo = hi = None
    starts = defaultdict(int)
    ends = defaultdict(int)
    spans = []
    for ln in out.splitlines():
        f = ln.split("\t")
        p = int(f[3])
        q = p
        n = 0
        for ch in f[5]:
            if ch.isdigit():
                n = n * 10 + ord(ch) - 48
            else:
                if ch in "MDN=X":
                    q += n
                n = 0
        spans.append((p, q))
        starts[p // 100] += 1
        ends[q // 100] += 1
        lo = p if lo is None else min(lo, p)
        hi = q if hi is None else max(hi, q)
    if lo is None:
        continue
    # reads3: drop reads whose start AND end bins are both singletons (isolated/readthrough tails)
    keep = [(a, b) for (a, b) in spans if starts[a // 100] >= 3 or ends[b // 100] >= 3]
    lo3 = min((a for a, _ in keep), default=lo)
    hi3 = max((b for _, b in keep), default=hi)
    T = me - ms
    rows.append((gene, nx, (pe - ps) / T, (hi - lo) / T, (hi3 - lo3) / T))


def med(v):
    v = sorted(v)
    return v[len(v) // 2] if v else 0.0


print(f"detected truth members: {len(rows)}\n")
print(f"{'size / truth size':<22}{'median':>9}{'>=0.8x':>9}{'<0.5x':>9}")
for lab, i in (("rep span (shipped)", 2), ("read union extent", 3), ("read extent, >=3 support", 4)):
    v = [r[i] for r in rows]
    print(f"{lab:<22}{med(v):>9.2f}{100*sum(1 for x in v if x>=0.8)/len(v):>8.0f}%"
          f"{100*sum(1 for x in v if x<0.5)/len(v):>8.0f}%")

for lab, pred in (("single-exon rep", lambda r: r[1] == 1), ("spliced rep", lambda r: r[1] > 1)):
    sub = [r for r in rows if pred(r)]
    if sub:
        print(f"\n  {lab} (n={len(sub)}): rep {med([r[2] for r in sub]):.2f}  "
              f"-> reads {med([r[3] for r in sub]):.2f}  -> reads>=3 {med([r[4] for r in sub]):.2f}")

worst = sorted(rows, key=lambda r: r[2])[:8]
print(f"\n{'gene':<14}{'nexon':>6}{'rep':>7}{'reads':>8}{'reads>=3':>10}")
for g, nx, a, b, c in worst:
    print(f"{g:<14}{nx:>6}{a:>7.2f}{b:>8.2f}{c:>10.2f}")

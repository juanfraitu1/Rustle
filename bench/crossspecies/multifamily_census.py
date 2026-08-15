#!/usr/bin/env python3
"""Is the shared-read certificate INFORMATIVE, or does it link everything?

For each gene family: what fraction of WITHIN-family gene pairs is linked by >=k shared reads,
and what fraction of BETWEEN-family pairs is? A certificate that lights up between families as
readily as within them is vacuous, so the between-family number is the load-bearing control --
not a nice-to-have.

⚠ CERTIFICATE ONLY. Read conflict (E_c) is the ambiguity oracle and belongs to O2 (assignment).
  It must never become the O1 edge rule, or O1 and O2 stop being independent.
⚠ Keeps SECONDARY alignments (flag 256) on purpose -- multimapping is the measurement. Drops
  unmapped (4) and supplementary (2048) via -F 2052.
⚠ Shared reads can also arise from a repeat present in BOTH families; that is exactly what the
  between-family control is there to expose.

usage: multifamily_census.py BAM GENE_BED FAM1,FAM2,... [out.tsv]
"""
import subprocess
import sys
from collections import defaultdict
from itertools import combinations

bam, gene_bed, fams_arg = sys.argv[1], sys.argv[2], sys.argv[3]
out_tsv = sys.argv[4] if len(sys.argv) > 4 else None
FAMS = fams_arg.split(",")
THRESH = (1, 3, 10)

genes = []
for line in open(gene_bed):
    f = line.rstrip("\n").split("\t")
    if len(f) >= 4:
        genes.append((f[0], int(f[1]), int(f[2]), f[3]))

fam_of, members = {}, defaultdict(list)
for c, s, e, g in genes:
    for fam in FAMS:
        if g.startswith(fam):
            # longest-prefix wins so NOTCH2NL is not swallowed by a shorter family name
            if g not in fam_of or len(fam) > len(fam_of[g]):
                fam_of[g] = fam
for c, s, e, g in genes:
    if g in fam_of:
        members[fam_of[g]].append((c, s, e, g))

print("family members found:")
for f in FAMS:
    print(f"  {f:<10} {len(members[f]):>3}  " + ", ".join(sorted(x[3] for x in members[f])[:8])
          + (" ..." if len(members[f]) > 8 else ""))

touch = defaultdict(set)
for f in FAMS:
    for c, s, e, g in members[f]:
        p = subprocess.run(["samtools", "view", "-F", "2052", bam, f"{c}:{s+1}-{e}"],
                           capture_output=True, text=True)
        for line in p.stdout.splitlines():
            q = line.split("\t", 1)[0]
            touch[q].add((f, g))

pair = defaultdict(int)
for v in touch.values():
    if len(v) < 2:
        continue
    for a, b in combinations(sorted(v), 2):
        pair[(a, b)] += 1

rows = []
print(f"\n{'family':<10}{'genes':>6}{'pairs':>7}" +
      "".join(f"{'within>='+str(k):>13}" for k in THRESH) +
      "".join(f"{'betw>='+str(k):>12}" for k in THRESH))
for f in FAMS:
    ms = sorted(x[3] for x in members[f])
    if len(ms) < 2:
        continue
    within = list(combinations(ms, 2))
    wc = {k: 0 for k in THRESH}
    for a, b in within:
        n = pair.get(((f, a), (f, b)), 0) or pair.get(((f, b), (f, a)), 0)
        for k in THRESH:
            if n >= k:
                wc[k] += 1
    other = [(g, o) for o in FAMS if o != f for g in sorted(x[3] for x in members[o])]
    bc = {k: 0 for k in THRESH}
    bt = 0
    for a in ms:
        for b, o in other:
            bt += 1
            n = pair.get(((f, a), (o, b)), 0) or pair.get(((o, b), (f, a)), 0)
            for k in THRESH:
                if n >= k:
                    bc[k] += 1
    rows.append((f, len(ms), len(within), wc, bc, bt))
    print(f"{f:<10}{len(ms):>6}{len(within):>7}" +
          "".join(f"{wc[k]}/{len(within)} ({wc[k]/len(within):.0%})".rjust(13) for k in THRESH) +
          "".join(f"{bc[k]/bt:.1%}".rjust(12) for k in THRESH))

tw = {k: sum(r[3][k] for r in rows) for k in THRESH}
tp = sum(r[2] for r in rows)
tb = {k: sum(r[4][k] for r in rows) for k in THRESH}
tbt = sum(r[5] for r in rows)
print(f"\n{'TOTAL':<10}{'':>6}{tp:>7}" +
      "".join(f"{tw[k]}/{tp} ({tw[k]/max(tp,1):.0%})".rjust(13) for k in THRESH) +
      "".join(f"{tb[k]/max(tbt,1):.1%}".rjust(12) for k in THRESH))
for k in THRESH:
    w, b = tw[k] / max(tp, 1), tb[k] / max(tbt, 1)
    print(f"  at >={k:>2} shared reads: within {w:.1%} vs between {b:.1%}"
          + (f"   enrichment {w/b:.0f}x" if b > 0 else "   between = ZERO"))

if out_tsv:
    with open(out_tsv, "w") as fh:
        fh.write("famA\tgeneA\tfamB\tgeneB\tshared_reads\twithin\n")
        for (a, b), n in sorted(pair.items(), key=lambda x: -x[1]):
            fh.write(f"{a[0]}\t{a[1]}\t{b[0]}\t{b[1]}\t{n}\t{int(a[0]==b[0])}\n")
    print(f"\nwrote {out_tsv}")

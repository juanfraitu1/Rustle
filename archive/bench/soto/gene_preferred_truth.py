#!/usr/bin/env python3
"""Build a GENE-PREFERRED truth BED from the Soto member BED.

Soto's member intervals come from the SD98 self-alignment -- they are DUPLICATION BLOCKS, not gene models.
A block can be far larger than the transcribed gene inside it (or, rarely, clip it). Scoring an RNA
prediction's SIZE against a duplicon therefore measures the wrong thing: the RNA can only ever span the gene.

Rule: if a NAMED RefSeq gene whose name matches the member's gene label overlaps the member interval, use the
GENE's span as the truth interval; otherwise keep Soto's block (nothing better is available). Emits the
rewritten BED plus a report of how much each interval moved.
"""
import gzip, re, sys
from collections import defaultdict

GFF = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz"
BED = sys.argv[1] if len(sys.argv) > 1 else "bench/soto/80_fams.chr.bed"
OUT = sys.argv[2] if len(sys.argv) > 2 else "bench/soto/80_fams.gene_preferred.bed"

members = []
for ln in open(BED):
    c, s, e, name, *rest = ln.rstrip("\n").split("\t")
    gene, fam = name.split("|")[0], name.split("|")[1]
    members.append([c, int(s), int(e), name, gene, fam, rest])

# index every named RefSeq gene by (chrom, upper name) -> list of spans
genes = defaultdict(list)
with gzip.open(GFF, "rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("gene", "pseudogene"): continue
        m = re.search(r"(?:^|;)Name=([^;]+)", f[8])
        if not m: continue
        genes[(f[0], m.group(1).upper())].append((int(f[3]) - 1, int(f[4])))

out = open(OUT, "w")
changed = shrunk = grown = moved = kept = 0
rows = []
for c, s, e, name, gene, fam, rest in members:
    cand = [(gs, ge) for (gs, ge) in genes.get((c, gene.upper()), []) if gs < e and s < ge]
    if cand:
        gs, ge = max(cand, key=lambda x: min(e, x[1]) - max(s, x[0]))  # best-overlapping annotation
        if (gs, ge) != (s, e):
            changed += 1
            r = (ge - gs) / max(e - s, 1)
            # r == 1.0 means the interval only MOVED, it did not resize -- counting that as "grown"
            # overstates growth, so it gets its own bucket.
            if   r < 1: shrunk += 1
            elif r > 1: grown  += 1
            else:       moved  += 1
            rows.append((gene, e - s, ge - gs, r))
        s, e = gs, ge
    else:
        kept += 1
    out.write("\t".join([c, str(s), str(e), name] + rest) + "\n")
out.close()

print(f"members            : {len(members)}")
print(f"  matched a named annotation : {len(members)-kept}")
print(f"  no named match (Soto kept) : {kept}")
print(f"  interval CHANGED           : {changed}  (shrunk {shrunk}, grown {grown}, moved-only {moved})")
if rows:
    import statistics
    rr = sorted(r for _, _, _, r in rows)
    print(f"  gene/block size ratio: median {statistics.median(rr):.2f}  "
          f"p10 {rr[len(rr)//10]:.2f}  p90 {rr[9*len(rr)//10]:.2f}")
    # BOTH tails, deliberately. The script's premise is that Soto blocks are LARGER than the gene inside
    # them, so printing only the shrinks would confirm the premise by construction and hide that the
    # majority of corrections go the other way (a block can also CLIP a gene that extends past it).
    print("\n  largest SHRINKS (Soto block >> annotated gene) -- duplicon, not gene model:")
    for g, b, gl, r in sorted(rows, key=lambda x: x[3])[:8]:
        print(f"    {g:14s} block {b:>8}  gene {gl:>8}  {r:.3f}x")
    print("\n  largest GROWTHS (annotated gene extends past the Soto block):")
    for g, b, gl, r in sorted(rows, key=lambda x: -x[3])[:8]:
        print(f"    {g:14s} block {b:>8}  gene {gl:>8}  {r:.1f}x")
    n_gt2 = sum(1 for _, _, _, r in rows if r > 2)
    print(f"\n  direction summary: {sum(1 for _,_,_,r in rows if r<1)} shrink, "
          f"{sum(1 for _,_,_,r in rows if r>1)} grow ({n_gt2} by more than 2x)")
print(f"\nwrote {OUT}")

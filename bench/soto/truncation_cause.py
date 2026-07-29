#!/usr/bin/env python3
"""Is the size truncation OUR fragmentation, or is the RNA genuinely not spanning the locus?

Section 8 found predicted copies are a median 0.54x the truth member's genomic span. Two very different
causes give the same number:
  (a) WE FRAGMENT   - reads cover the whole gene but our copy captures only part of it  -> fixable
  (b) RNA IS SHORT  - no reads outside our span; the copy simply is not transcribed there -> inherent, and
                      the METRIC is at fault for comparing a transcript to a whole gene model

Discriminator: take each truncated matched pair, list the gene's ANNOTATED EXONS that fall OUTSIDE the
predicted span, and count primary reads over them. Exons with real read support that we failed to include
are our fragmentation; exons with no reads are RNA that does not exist.

Note on splicing: a spliced transcript's genomic span still CONTAINS its introns, so being spliced does not
by itself shrink the span. Only a missing 5'/3' portion does.
"""
import subprocess, gzip, re, sys
from collections import defaultdict
import numpy as np
from scipy.optimize import linear_sum_assignment

SAM   = "/home/juanfra/miniforge3/bin/samtools"
BAM   = "/home/juanfra/winloci_scratch/soto_cache/soto_regions.bam"
GFF   = "/mnt/c/Users/jfris/Desktop/Reference/chm13v2.0_RefSeq_full.gff.gz"
CAT   = sys.argv[1] if len(sys.argv) > 1 else "/home/juanfra/winloci_scratch/soto_cache/definitive.copies.tsv"
BED   = "bench/soto/80_fams.gene_preferred.bed"
MIN_CONTAIN = 0.50
MIN_EXON_READS = 3          # an exon counts as EXPRESSED at this many primary reads

truth = []
for ln in open(BED):
    c, s, e, name, *_ = ln.rstrip("\n").split("\t")
    truth.append((name.split("|")[0], c, int(s), int(e)))
pred = []
for i, ln in enumerate(open(CAT)):
    if i == 0: continue
    f = ln.rstrip("\n").split("\t")
    if len(f) >= 6: pred.append((f[3], int(f[4]), int(f[5]), f[2]))

# same 1:1 matching as bipartite_size_match.py
g = np.zeros((len(pred), len(truth)))
for i, (pc, ps, pe, _t) in enumerate(pred):
    for j, (_gn, tc, ts, te) in enumerate(truth):
        if pc != tc: continue
        ov = min(pe, te) - max(ps, ts)
        if ov > 0 and max(ov/max(pe-ps,1), ov/max(te-ts,1)) >= MIN_CONTAIN:
            g[i, j] = ov / max(pe-ps, te-ts, 1)
ri, cj = linear_sum_assignment(-g)
pairs = [(i, j) for i, j in zip(ri, cj) if g[i, j] > 0]

# exons per gene name
exons = defaultdict(list)
with gzip.open(GFF, "rt") as fh:
    for ln in fh:
        if ln.startswith("#"): continue
        f = ln.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon": continue
        m = re.search(r'(?:^|;)gene=([^;]+)', f[8])
        if m: exons[(f[0], m.group(1).upper())].append((int(f[3])-1, int(f[4])))

def reads_over(chrom, s, e):
    out = subprocess.run([SAM, "view", "-c", "-F", "2308", BAM, f"{chrom}:{s+1}-{e}"],
                         capture_output=True, text=True).stdout.strip()
    return int(out) if out.isdigit() else 0

print(f"{'gene':14s} {'pred_bp':>9} {'true_bp':>9} {'ratio':>6} {'ex_in':>6} {'ex_out':>7} "
      f"{'out_expr':>9} {'reads_out':>10}  verdict")
frag = inherent = noann = 0
rows = []
for i, j in pairs:
    plen = pred[i][2] - pred[i][1]; tlen = truth[j][3] - truth[j][2]
    r = plen / max(tlen, 1)
    if r > 0.5: continue                       # only the truncated ones
    gene, tc, ts, te = truth[j]
    ex = [(a, b) for (a, b) in exons.get((tc, gene.upper()), []) if a < te and ts < b]
    if not ex:
        noann += 1; continue
    ps, pe = pred[i][1], pred[i][2]
    inside  = [(a, b) for (a, b) in ex if a < pe and ps < b]
    outside = [(a, b) for (a, b) in ex if not (a < pe and ps < b)]
    # merge outside exons to avoid double-counting overlapping isoform exons
    outside = sorted(set(outside))
    merged, cur = [], None
    for a, b in outside:
        if cur and a <= cur[1]: cur = (cur[0], max(cur[1], b))
        else:
            if cur: merged.append(cur)
            cur = (a, b)
    if cur: merged.append(cur)
    n_expr = sum(1 for a, b in merged if reads_over(tc, a, b) >= MIN_EXON_READS)
    tot = sum(reads_over(tc, a, b) for a, b in merged)
    verdict = "WE FRAGMENT (reads outside)" if n_expr >= 1 else "RNA absent outside (inherent)"
    if n_expr >= 1: frag += 1
    else: inherent += 1
    rows.append((gene, plen, tlen, r, len(inside), len(merged), n_expr, tot, verdict))

for gene, plen, tlen, r, ni, no, ne, tot, v in sorted(rows, key=lambda x: x[3])[:30]:
    print(f"{gene:14s} {plen:>9} {tlen:>9} {r:>6.2f} {ni:>6} {no:>7} {ne:>9} {tot:>10}  {v}")
print(f"\ntruncated pairs analysed: {len(rows)}   (+{noann} with no annotated exon, skipped)")
print(f"  WE FRAGMENT  (>=1 unincluded exon has >={MIN_EXON_READS} reads): {frag}")
print(f"  RNA ABSENT   (no reads on any unincluded exon)                 : {inherent}")

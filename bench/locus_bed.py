#!/usr/bin/env python3
"""Emit predicted loci as BED, and match them one-to-one against an annotation to show size agreement.

Answers "do the loci we call have roughly the extent the annotation gives them, and how far off are we?"
by collapsing each GTF to LOCI (one record per `gene_id`, spanning its transcripts), writing them as BED,
and greedily matching predicted ↔ annotated loci on reciprocal overlap.

⚠ The matching here is EVALUATION ONLY. Loci are never built with bipartite matching — that is a standing
project constraint; this tool only scores loci that were already built without it.

Matching is GREEDY on reciprocal overlap (take the best available pair, remove both, repeat), not optimal
assignment. For the size question that is enough and the tie behaviour is deterministic (ties break on
locus id), but it is a lower bound on the optimal matching, and it is labelled as such in the output.

usage: locus_bed.py PRED.gtf --out PREFIX [--ref REF.gtf] [--min-overlap 0.10]
  PREFIX.loci.bed        predicted loci (BED6; score = summed `reads`, capped at 1000)
  PREFIX.ref_loci.bed    annotated loci, when --ref is given
  PREFIX.locus_match.tsv one row per matched pair, with the size ratio
"""
import argparse, collections, re, sys

def loci(path):
    """gene_id -> (chrom, start0, end, strand, reads, n_tx, spliced_len)"""
    tx = collections.defaultdict(list)
    gene_of, reads_of = {}, {}
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8])
        g = re.search(r'gene_id "([^"]+)"', f[8])
        if not t:
            continue
        t = t.group(1)
        if g:
            gene_of.setdefault(t, g.group(1))
        r = re.search(r'reads "(\d+)"', f[8])
        if r:
            reads_of[t] = max(reads_of.get(t, 0), int(r.group(1)))
        if f[2] == 'exon':
            tx[t].append((f[0], int(f[3]) - 1, int(f[4]), f[6]))
    out = {}
    for t, ex in tx.items():
        g = gene_of.get(t, t)
        ex.sort(key=lambda x: x[1])
        c, s, e, st = ex[0][0], ex[0][1], ex[-1][2], ex[0][3]
        sl = sum(b - a for _, a, b, _ in ex)
        if g in out:
            oc, os_, oe, ost, orr, on, osl = out[g]
            out[g] = (oc, min(os_, s), max(oe, e), ost, orr + reads_of.get(t, 0), on + 1, max(osl, sl))
        else:
            out[g] = (c, s, e, st, reads_of.get(t, 0), 1, sl)
    return out

def write_bed(d, path):
    with open(path, 'w') as fo:
        for g, (c, s, e, st, r, n, _) in sorted(d.items(), key=lambda kv: (kv[1][0], kv[1][1])):
            fo.write(f"{c}\t{s}\t{e}\t{g}\t{min(r,1000)}\t{st if st in '+-' else '.'}\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('pred'); ap.add_argument('--out', required=True)
    ap.add_argument('--ref'); ap.add_argument('--min-overlap', type=float, default=0.10)
    a = ap.parse_args()

    P = loci(a.pred)
    write_bed(P, f"{a.out}.loci.bed")
    print(f"  {len(P)} predicted loci -> {a.out}.loci.bed", file=sys.stderr)
    if not a.ref:
        return
    R = loci(a.ref)
    write_bed(R, f"{a.out}.ref_loci.bed")
    print(f"  {len(R)} annotated loci -> {a.out}.ref_loci.bed", file=sys.stderr)

    # candidate pairs: same contig, any overlap
    bych = collections.defaultdict(list)
    for g, v in R.items():
        bych[v[0]].append((g, v))
    cands = []
    for pg, pv in P.items():
        for rg, rv in bych.get(pv[0], []):
            ov = min(pv[2], rv[2]) - max(pv[1], rv[1])
            if ov <= 0:
                continue
            rec = min(ov / max(1, pv[2] - pv[1]), ov / max(1, rv[2] - rv[1]))
            if rec >= a.min_overlap:
                cands.append((rec, pg, rg, ov))
    cands.sort(key=lambda x: (-x[0], x[1], x[2]))          # deterministic ties
    usedP, usedR, pairs = set(), set(), []
    for rec, pg, rg, ov in cands:
        if pg in usedP or rg in usedR:
            continue
        usedP.add(pg); usedR.add(rg); pairs.append((pg, rg, rec, ov))

    with open(f"{a.out}.locus_match.tsv", 'w') as fo:
        fo.write("pred_locus\tref_locus\tchrom\tpred_start\tpred_end\tref_start\tref_end\t"
                 "pred_span\tref_span\tsize_ratio\trecip_overlap\tpred_reads\tpred_n_tx\n")
        ratios = []
        for pg, rg, rec, ov in sorted(pairs, key=lambda x: x[0]):
            pc, ps, pe, _, pr, pn, _ = P[pg]; _, rs, re_, _, _, _, _ = R[rg]
            psp, rsp = pe - ps, re_ - rs
            ratio = psp / rsp if rsp else float('nan')
            ratios.append(ratio)
            fo.write(f"{pg}\t{rg}\t{pc}\t{ps}\t{pe}\t{rs}\t{re_}\t{psp}\t{rsp}\t{ratio:.4f}\t{rec:.4f}\t{pr}\t{pn}\n")
    ratios.sort()
    n = len(ratios)
    med = ratios[n // 2] if n else float('nan')
    within = lambda lo, hi: sum(1 for r in ratios if lo <= r <= hi)
    print(f"\n  GREEDY one-to-one matching (evaluation only), min reciprocal overlap {a.min_overlap}")
    print(f"    matched pairs            : {n}  of {len(P)} predicted / {len(R)} annotated")
    print(f"    predicted loci unmatched : {len(P)-n}    annotated loci unmatched: {len(R)-n}")
    if n:
        print(f"    size ratio pred/ref      : median {med:.3f}  "
              f"q25 {ratios[n//4]:.3f}  q75 {ratios[3*n//4]:.3f}")
        for lo, hi in ((0.9, 1.1), (0.8, 1.25), (0.5, 2.0)):
            print(f"      within [{lo},{hi}]        : {within(lo,hi)} ({100*within(lo,hi)/n:.1f}%)")
    print(f"    -> {a.out}.locus_match.tsv")
main()

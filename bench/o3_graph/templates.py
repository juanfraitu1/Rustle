#!/usr/bin/env python3
"""Curated templates for the graph-O3 simulation (prereg docs/PREREG_o3_graph_sim_2026-10-15.md): real genes whose curated
isoforms (RefSeq / CAT) are exon-skipping variants of their canonical transcript, so famsim can express every one of them
on a copy (famsim isoforms are skip sets over the model's exons).

    templates.py pick ANNOTATION GENOME OUT.json [--n 24] [--seed 1] [--max-span 40000] [--min-exons 4] [--max-exons 14]

An isoform qualifies when every internal exon is a canonical exon, its first exon ends at a canonical exon's end and its last
exon starts at a canonical exon's start (genomic order); it is then the canonical with the unmapped exons skipped (alternative
ends become skipped end exons). Skip indices are famsim's: 1-based, in TRANSCRIPT orientation. A template needs >= 1
qualifying isoform that skips an INTERNAL exon (a cassette, flanked by kept exons), every canonical intron GT..AG, and a
curated gene (not a transcript-assembly `MSTRG.*` one).
"""
import argparse
import json
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from famsim.template import canonical_transcript, load_genes  # noqa: E402


def cassette_isoforms(rec):
    """-> (canonical transcript id, [(transcript id, skip tuple)]) for the qualifying isoforms, distinct skip sets, sorted."""
    canon_id, canon = canonical_transcript(rec)
    canon = sorted(canon)
    n = len(canon)
    order = list(range(n)) if rec["strand"] == "+" else list(range(n - 1, -1, -1))
    tx_index = {g: i + 1 for i, g in enumerate(order)}           # genomic index -> 1-based transcript-orientation index
    starts = {s: k for k, (s, _) in enumerate(canon)}
    ends = {e: k for k, (_, e) in enumerate(canon)}
    exact = {x: k for k, x in enumerate(canon)}
    out, seen = [], set()
    for tid, ex in sorted(rec["transcripts"].items()):
        if tid == canon_id or len(ex) < 2:
            continue
        ex = sorted(ex)
        first, last, mid = ex[0], ex[-1], ex[1:-1]
        if first[1] not in ends or last[0] not in starts or any(x not in exact for x in mid):
            continue
        mapped = [ends[first[1]]] + [exact[x] for x in mid] + [starts[last[0]]]
        if any(b <= a for a, b in zip(mapped, mapped[1:])):
            continue
        skip = tuple(sorted(tx_index[g] for g in range(n) if g not in set(mapped)))
        if skip and skip not in seen:
            seen.add(skip)
            out.append((tid, skip))
    return canon_id, out


def curated_name(gene):
    """CAT also carries transcript-assembly genes (StringTie `MSTRG.*`); they are not curated."""
    return not gene.startswith("MSTRG.")


def internal_cassettes(skips, n):
    """Exon indices (1-based) some isoform skips while keeping an exon on each side."""
    out = set()
    for sk in skips:
        kept = [k for k in range(1, n + 1) if k not in sk]
        for k in sk:
            if 1 < k < n and any(x < k for x in kept) and any(x > k for x in kept):
                out.add(k)
    return sorted(out)


def gt_ag(fa, chrom, strand, exons):
    from famsim.model import rc
    ex = sorted(exons)
    for (s1, e1), (s2, e2) in zip(ex, ex[1:]):
        intr = fa.fetch(chrom, e1, s2 - 1).upper()            # 1-based inclusive exons -> intron [e1, s2-1)
        if strand == "-":
            intr = rc(intr)
        if not (intr.startswith("GT") and intr.endswith("AG")):
            return False
    return True


def pick(annotation, genome, n, seed, max_span=40000, min_exons=4, max_exons=14):
    import pysam
    fa = pysam.FastaFile(genome)
    genes = load_genes(annotation)
    cands = []
    for g, rec in sorted(genes.items()):
        if rec["strand"] not in "+-" or not rec["transcripts"] or not curated_name(g):
            continue
        cid, isos = cassette_isoforms(rec)
        canon = sorted(rec["transcripts"][cid])
        if not (min_exons <= len(canon) <= max_exons) or canon[-1][1] - canon[0][0] > max_span:
            continue
        cas = internal_cassettes([s for _, s in isos], len(canon))
        if not cas or not gt_ag(fa, rec["chrom"], rec["strand"], canon):
            continue
        cands.append(dict(gene=g, chrom=rec["chrom"], strand=rec["strand"], canonical=cid, n_exons=len(canon),
                          span=canon[-1][1] - canon[0][0] + 1, isoforms=[dict(tid=t, skip=list(s)) for t, s in isos],
                          cassettes=cas))
    rng = random.Random(seed)
    rng.shuffle(cands)
    return cands[:n], len(cands)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pick")
    p.add_argument("annotation")
    p.add_argument("genome")
    p.add_argument("out")
    p.add_argument("--n", type=int, default=24)
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--max-span", type=int, default=40000)
    p.add_argument("--min-exons", type=int, default=4)
    p.add_argument("--max-exons", type=int, default=14)
    a = ap.parse_args(argv)
    picked, n_cand = pick(a.annotation, a.genome, a.n, a.seed, a.max_span, a.min_exons, a.max_exons)
    json.dump(dict(annotation=a.annotation, genome=a.genome, seed=a.seed, candidates=n_cand, templates=picked),
              open(a.out, "w"), indent=1)
    print(f"{len(picked)} of {n_cand} candidate templates -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

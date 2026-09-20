#!/usr/bin/env python3
"""
family_def_denovo_intron.py — does a VERTEX-level intron requirement clean the family
definition?

Hypothesis: the residual "family ⊋ paralogy" contamination (retrocopies / processed
pseudogenes / read-through bridges) shows up as edges touching an INTRONLESS de-novo
locus — because a processed pseudogene is an intronless genomic copy of a mature mRNA,
so reads assemble there as a single-exon locus. Requiring BOTH family members to be
multi-exon (>=1 intron) should therefore strip those bridges while keeping real
multi-exon paralog arrays.

Builds the de-novo read-conflict graph (same de-tie scan as family_def_genomewide, but
over the production de-novo loci, carrying each locus's n_exon), then compares the
baseline component families to the multi-exon-filtered ones:
  - cross-chrom bridge families (>=3 chromosomes = the documented artifact)
  - co-located families (<=2 chrom = real tandem arrays / segdup pairs)
  - exon structure of REMOVED vs KEPT edges (does the filter cut single-exon edges?)
  - cost: real co-located families lost to a single-exon member.

Run: /home/juanfra/miniforge3/bin/python bench/family_def_denovo_intron.py
"""
import collections
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_genomewide import (best_gene, ref_span, de_of, pair_evidence,
                                    components, BAM, DELTA, DE_MAX, MIN_READS)

META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"


def load_loci():
    by_chrom = collections.defaultdict(list)
    nexon = {}; chrom = {}
    with open(META) as f:
        next(f)  # id chrom start end strand n_exon n_reads
        for line in f:
            p = line.rstrip("\n").split("\t")
            vid, c, s, e = p[0], p[1], int(p[2]), int(p[3])
            by_chrom[c].append((s, e, vid))
            nexon[vid] = int(p[5]); chrom[vid] = c
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom, nexon, chrom


def scan(by_chrom):
    mm = collections.defaultdict(dict)
    for flag in ("0x100", None):
        args = ["samtools", "view", "-f", "0x100", BAM] if flag else ["samtools", "view", "-F", "0x900", BAM]
        p = subprocess.Popen(args, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        for line in p.stdout:
            f = line.split("\t")
            if len(f) < 9:
                continue
            de = de_of(f[11:])
            if de is None or de > DE_MAX:
                continue
            g = best_gene(by_chrom, f[2], int(f[3]) - 1, int(f[3]) - 1 + ref_span(f[5]))
            if g is None:
                continue
            d = mm[f[0]]
            if g not in d or de < d[g]:
                d[g] = de
        p.wait()
    return mm


def fam_stats(fams, chrom):
    def nchrom(c): return len({chrom.get(g) for g in c})
    coloc = [c for c in fams if nchrom(c) <= 2]
    bridge = [c for c in fams if nchrom(c) >= 3]
    return len(fams), len(coloc), len(bridge)


def main():
    by_chrom, nexon, chrom = load_loci()
    me = sum(1 for v in nexon if nexon[v] >= 2)
    print(f"[loci] {len(nexon):,} de-novo loci; multi-exon(>=2): {me:,} "
          f"({100*me/len(nexon):.0f}%); single-exon: {len(nexon)-me:,}", flush=True)
    print("[scan] de-tie read-conflict graph over de-novo loci ...", flush=True)
    mm = scan(by_chrom)
    ev = pair_evidence(mm)
    edges, fams = components(ev, DELTA, DE_MAX, MIN_READS)
    edge_set = {(a, b) for a, b, n in edges}

    # multi-exon vertex filter on edges
    def mexon(v): return nexon.get(v, 1) >= 2
    kept_edges = [(a, b, n) for a, b, n in edges if mexon(a) and mexon(b)]
    _, fams_f = components({(a, b): ev[(a, b)] for a, b, n in kept_edges},
                           DELTA, DE_MAX, MIN_READS)

    nB, cB, bB = fam_stats(fams, chrom)
    nF, cF, bF = fam_stats(fams_f, chrom)
    print(f"\n=== BASELINE de-novo families (component) ===")
    print(f"  edges={len(edges)}  families={nB}  co-located(<=2chr)={cB}  cross-chrom-bridge(>=3chr)={bB}")
    print(f"=== MULTI-EXON-FILTERED (both loci >=1 intron) ===")
    print(f"  edges={len(kept_edges)}  families={nF}  co-located={cF}  cross-chrom-bridge={bF}")
    print(f"  -> bridges {bB}->{bF}  ({100*(bB-bF)/max(bB,1):.0f}% removed); "
          f"co-located {cB}->{cF}  ({100*(cB-cF)/max(cB,1):.0f}% lost)")

    # exon structure of removed vs kept edges
    removed = [(a, b, n) for a, b, n in edges if not (mexon(a) and mexon(b))]
    def cls(a, b):
        return "cross_chrom" if chrom.get(a) != chrom.get(b) else "same_chrom"
    print(f"\n--- REMOVED edges ({len(removed)}): how many touch a single-exon locus? ---")
    rc = collections.Counter()
    for a, b, n in removed:
        se = (not mexon(a)) + (not mexon(b))
        rc[(f"{se}_single_exon_member", cls(a, b))] += 1
    for k, v in sorted(rc.items()):
        print(f"  {v:5d}  {k[0]:22} {k[1]}")
    # is the bridge population single-exon-enriched?
    bridge_edges = [(a, b) for a, b, n in edges
                    if chrom.get(a) != chrom.get(b)]
    bridge_se = sum(1 for a, b in bridge_edges if not (mexon(a) and mexon(b)))
    samechr_edges = [(a, b) for a, b, n in edges if chrom.get(a) == chrom.get(b)]
    same_se = sum(1 for a, b in samechr_edges if not (mexon(a) and mexon(b)))
    print(f"\n--- single-exon enrichment: cross-chrom vs same-chrom edges ---")
    print(f"  cross-chrom edges touching single-exon locus: {bridge_se}/{len(bridge_edges)} "
          f"= {bridge_se/max(len(bridge_edges),1):.2f}")
    print(f"  same-chrom  edges touching single-exon locus: {same_se}/{len(samechr_edges)} "
          f"= {same_se/max(len(samechr_edges),1):.2f}")

    out = dict(loci=len(nexon), multiexon=me,
               baseline=dict(edges=len(edges), families=nB, coloc=cB, bridges=bB),
               filtered=dict(edges=len(kept_edges), families=nF, coloc=cF, bridges=bF),
               crosschrom_singleexon_rate=bridge_se / max(len(bridge_edges), 1),
               samechrom_singleexon_rate=same_se / max(len(samechr_edges), 1))
    json.dump(out, open(os.path.join(HERE, "family_def_denovo_intron.json"), "w"), indent=2)
    print(f"\n[+] wrote {os.path.join(HERE,'family_def_denovo_intron.json')}")


if __name__ == "__main__":
    main()

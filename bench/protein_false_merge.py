#!/usr/bin/env python3
"""The FALSE-MERGE measurement that blocks r906's adoption.
Per `docs/PREREG_protein_false_merge_2026-09-22.md` (§6y4, md5 `ca4fa685`).

r906 measured protein edges as ⚠PARTIAL (+13.8 pts held-out pair coverage, 0 cross-family) and recorded
its own gap: "Precision measured only over Soto-labelled genes, so this is NOT a genome-wide false-merge
rate." Every protein edge with an UNLABELLED endpoint is unmeasured, and a sensitive mode's risk lives
exactly there.

§6ko's rule is reused verbatim from `bench/protein_edge_gap.py` -- nothing is re-tuned.
Classification of every PROTEIN-ONLY edge (no nucleotide edge in the shipped graph):
  TRUE     both endpoints Soto-labelled, same family
  FALSE    both endpoints Soto-labelled, different families
  UNKNOWN  at least one endpoint unlabelled  <- r906's blind spot
and each UNKNOWN split by STRUCTURAL corroboration (label-free): does ANY nucleotide PAF record exist for
the pair, even one that failed the gate?
"""
import argparse, collections, csv, os, re, subprocess, sys
sys.path.insert(0, 'bench')
import protein_edge_gap as PEG


def main():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--genome', '--paf', '--chrom', '--soto', '--graph', '--out'):
        ap.add_argument(x, required=True)
    ap.add_argument('--threads', default='4')
    a = ap.parse_args()
    import pysam

    cds = PEG.longest_cds(a.gff, a.chrom)
    fa = pysam.FastaFile(a.genome)
    faa = a.out + '.proteins.faa'
    plen = {}
    with open(faa, 'w') as fh:
        for g, (st, segs) in cds.items():
            p = PEG.translate(fa, a.chrom, st, segs)
            if len(p) >= 10:
                plen[g] = len(p); fh.write(f'>{g}\n{p}\n')
    P = PEG.protein_edges(faa, a.out, plen, a.threads)          # §6ko's rule, unchanged

    # symbol <-> node, from the GFF gene spans (PAF/graph names are chrom:start-end)
    gene_at, span_of = {}, {}
    for line in open(a.gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != a.chrom or f[2] not in ('gene', 'pseudogene', 'ncRNA_gene'):
            continue
        n = re.search(r'Name=([^;]+)', f[8])
        if n:
            s, e = int(f[3]), int(f[4])
            gene_at[f'{a.chrom}:{s}-{e}'] = n.group(1)
            gene_at[f'{a.chrom}:{s-1}-{e}'] = n.group(1)
            span_of[n.group(1)] = (s, e)

    nuc = set()
    for line in open(a.graph):
        f = line.rstrip('\n').split('\t')
        if len(f) == 3 and f[0] != f[1]:
            x, y = gene_at.get(f[0]), gene_at.get(f[1])
            if x and y and x != y:
                nuc.add((x, y) if x < y else (y, x))
    # ANY paf record at all, even one the gate rejected -- label-free corroboration
    anyaln = set()
    for line in open(a.paf):
        f = line.split('\t', 6)
        x, y = gene_at.get(f[0]), gene_at.get(f[5])
        if x and y and x != y:
            anyaln.add((x, y) if x < y else (y, x))

    fam = {}
    for r in csv.DictReader(open(a.soto), delimiter='\t'):
        fi, g = (r.get('Family ID') or '').strip(), (r.get('Gene Name') or '').strip()
        if fi and fi != 'N/A' and g:
            fam.setdefault(g, fi)

    only = {e for e in P if e not in nuc}
    T = F = 0
    unk_corr = unk_bare = 0
    bare_examples = []
    for x, y in sorted(only):
        if x in fam and y in fam:
            if fam[x] == fam[y]:
                T += 1
            else:
                F += 1
        else:
            if (x, y) in anyaln:
                unk_corr += 1
            else:
                unk_bare += 1
                if len(bare_examples) < 6:
                    bare_examples.append((x, y))
    n = len(only)
    opt = F / (T + F) if (T + F) else float('nan')
    pes = (F + unk_bare) / n if n else float('nan')
    print(f"{a.chrom}: protein edges {len(P)} | already nucleotide {len(P)-n} | "
          f"⭐PROTEIN-ONLY (what the mode ADDS) {n}")
    print(f"   TRUE  (both labelled, same family)      {T}")
    print(f"   FALSE (both labelled, diff family)      {F}")
    print(f"   UNKNOWN, sub-threshold alignment exists {unk_corr}")
    print(f"   UNKNOWN, NO alignment at all            {unk_bare}")
    print(f"   ⭐ optimistic false-merge  F/(T+F)      {opt:.4f}"
          if (T + F) else "   ⭐ optimistic false-merge  F/(T+F)      n/a (no labelled pair)")
    print(f"   ⭐ pessimistic (F+bare)/all             {pes:.4f}")
    if bare_examples:
        print(f"   uncorroborated examples: {', '.join(f'{x}~{y}' for x, y in bare_examples)}")


if __name__ == '__main__':
    main()

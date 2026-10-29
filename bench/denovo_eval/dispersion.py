#!/usr/bin/env python3
"""Chromosome-dispersion family filter (PREREG DISPERSION, 2026-10-29). A de novo family whose loci lie on >= T distinct chromosomes is
discarded whole: on HG002 (design set) such families are young transposable-element families (shared segments 99% soft-masked, 93%
inter-chromosomal, identity 0.99) that glue unrelated genes together; true paralog families lie on 1-3 chromosomes.
Sequences shorter than --min-chrom-bp (default 10 Mb: unplaced scaffolds, alts, mito) all count as ONE pseudo-chromosome 'unplaced'.

  dispersion.py FAM_COPIES_TSV GENOME_FAI [--max-chrom 3] [--min-chrom-bp 10000000] [--recluster LOCI_GFF3 GATED_PAF] > filtered.fam.copies.tsv
Keeps families with <= --max-chrom chromosomes (so T = max-chrom + 1 drops). Reports dropped families on stderr.
--recluster: a dispersed family is not discarded but split into the connected components of its INTRA-chromosomal gated-PAF records
(TE-driven inter-chromosomal links removed); components of >= 2 loci are re-emitted as families <id>.r<k>, singletons dropped.
"""
import sys, collections
a = sys.argv[1:]
maxc = int(a[a.index('--max-chrom') + 1]) if '--max-chrom' in a else 3
minbp = int(a[a.index('--min-chrom-bp') + 1]) if '--min-chrom-bp' in a else 10_000_000
copies, fai = a[0], a[1]
big = {l.split('\t')[0] for l in open(fai) if int(l.split('\t')[1]) >= minbp}
rows = [l for l in open(copies)]
head = rows[0]; chrom_i = head.rstrip('\n').split('\t').index('chrom')
chroms = collections.defaultdict(set)
for l in rows[1:]:
    f = l.rstrip('\n').split('\t')
    chroms[f[0]].add(f[chrom_i] if f[chrom_i] in big else 'unplaced')
drop = {k for k, v in chroms.items() if len(v) > maxc}
sys.stdout.write(head)
kept = 0
for l in rows[1:]:
    if l.split('\t')[0] not in drop:
        sys.stdout.write(l); kept += 1
if '--recluster' in a and drop:
    import re
    gff3, paf = a[a.index('--recluster') + 1], a[a.index('--recluster') + 2]
    tid_i = head.rstrip('\n').split('\t').index('tid')
    key_of = {}
    for line in open(gff3):
        if line[0] != '#':
            f = line.rstrip('\n').split('\t')
            if f[2] == 'gene':
                key_of[re.search(r'ID=gene-([^;]+)', f[8]).group(1)] = f"{f[0]}:{f[3]}-{f[4]}"
    fam_rows = collections.defaultdict(list); fam_of = {}
    for l in rows[1:]:
        f = l.rstrip('\n').split('\t')
        if f[0] in drop:
            k = key_of.get(f[tid_i])
            if k: fam_rows[f[0]].append((k, l)); fam_of[k] = f[0]
    parent = {k: k for k in fam_of}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for line in open(paf):
        f = line.split('\t', 6); q, t = f[0], f[5]
        if q in fam_of and fam_of.get(t) == fam_of[q] and q.split(':')[0] == t.split(':')[0]:
            parent[find(q)] = find(t)
    nre = 0
    for fid, lst in fam_rows.items():
        comps = collections.defaultdict(list)
        for k, l in lst: comps[find(k)].append(l)
        for i, (_, ls) in enumerate(sorted(comps.items(), key=lambda kv: -len(kv[1]))):
            if len(ls) >= 2:
                nre += 1
                for l in ls:
                    f = l.split('\t'); f[0] = f"{fid}.r{i}"; sys.stdout.write('\t'.join(f)); kept += 1
    print(f"[dispersion] recluster: {len(drop)} dispersed families -> {nre} intra-chromosomal components of >= 2 loci", file=sys.stderr)
print(f"[dispersion] max-chrom {maxc}: families {len(chroms)} -> kept {len(chroms) - len(drop)}, dropped {len(drop)} "
      f"({sum(1 for l in rows[1:] if l.split(chr(9))[0] in drop)} copies); copies kept {kept}", file=sys.stderr)
for k in sorted(drop, key=lambda k: -len(chroms[k])):
    print(f"[dispersion] drop {k}: {len(chroms[k])} chromosomes", file=sys.stderr)

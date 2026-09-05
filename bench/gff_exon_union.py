#!/usr/bin/env python3
"""Exon-union blocks per gene record, keyed exactly like mcl_families' PAF names (`chrom:start-end`, GFF
1-based verbatim), following Parent chains gene <- rna <- exon.  usage: gff_exon_union.py <gff|gff.gz> [keys_file]
Prints key<TAB>exonic_len<TAB>blocks (0-based half-open, comma-joined)."""
import sys, gzip, collections
gff, keys = sys.argv[1], (set(open(sys.argv[2]).read().split()) if len(sys.argv) > 2 else None)
op = gzip.open if gff.endswith('.gz') else open
gene_key, parent, exons = {}, {}, collections.defaultdict(list)
for l in op(gff, 'rt'):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9: continue
    a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
    if f[2] in ('gene', 'pseudogene'): gene_key[a['ID']] = f"{f[0]}:{f[3]}-{f[4]}"
    elif f[2] == 'exon': exons[a.get('Parent')].append((int(f[3]) - 1, int(f[4])))
    elif 'Parent' in a and 'ID' in a: parent[a['ID']] = a['Parent']
def gene_of(pid):
    seen = 0
    while pid and pid not in gene_key and seen < 5: pid = parent.get(pid); seen += 1
    return pid if pid in gene_key else None
union = collections.defaultdict(list)
for pid, ex in exons.items():
    g = gene_of(pid)
    if g: union[gene_key[g]].extend(ex)
for k, ex in union.items():
    if keys is not None and k not in keys: continue
    ex.sort(); merged = []
    for s, e in ex:
        if merged and s <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], e)
        else: merged.append([s, e])
    print(k, sum(e - s for s, e in merged), ','.join(f'{s}-{e}' for s, e in merged), sep='\t')

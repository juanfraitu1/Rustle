#!/usr/bin/env python3
"""G5 prototype on the gorilla catalog: family units mapped back to the substrate (minimap2 -x splice), hits outside
every annotated gene and every existing unit, admitted when >= min_reads primary reads have a block inside.
usage: rna_admit_gorilla.py <units_prefix> <genome.fa> <bam> <gff> <out.tsv> [min_reads=3]"""
import sys, csv, subprocess, collections, pysam, os
pre, fa, bam_p, gff, out = sys.argv[1:6]; min_reads = int(sys.argv[6]) if len(sys.argv) > 6 else 3
units = list(csv.DictReader(open(pre + '.units.tsv'), delimiter='\t'))
ulen = {f"{u['family_id']}|{u['copy_idx']}": 1 for u in units}
paf = subprocess.run(['minimap2', '-x', 'splice', '-N', '50', '-c', '-t', '4', fa, pre + '.units.fa'], capture_output=True, text=True).stdout
genes = collections.defaultdict(list)
for l in open(gff):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('gene', 'pseudogene'): continue
    a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
    genes[f[0]].append((int(f[3]) - 1, int(f[4]), a.get('gene', '?')))
uni = collections.defaultdict(list)
for u in units: uni[u['chrom']].append((int(u['start']), int(u['end'])))
hits = collections.defaultdict(list)
for l in paf.splitlines():
    f = l.split('\t'); fam = f[0].split('|')[0]; ident = int(f[9]) / max(1, int(f[10])); cov = (int(f[3]) - int(f[2])) / int(f[1])
    if ident < 0.70 or cov < 0.30: continue
    hits[fam].append((f[5], int(f[7]), int(f[8]), ident))
BAM = pysam.AlignmentFile(bam_p); rows = []
for fam, hs in hits.items():
    hs.sort(); merged = []
    for c, s, e, i in hs:
        if merged and merged[-1][0] == c and s <= merged[-1][2]: merged[-1][2] = max(merged[-1][2], e); merged[-1][3] = max(merged[-1][3], i)
        else: merged.append([c, s, e, i])
    for c, s, e, i in merged:
        if any(gs < e and s < ge for gs, ge, _ in genes.get(c, [])): continue
        if any(us < e and s < ue for us, ue in uni.get(c, [])): continue
        n = sum(1 for a in BAM.fetch(c, s, e) if not a.flag & 2308 and any(b1 > s and b0 < e for b0, b1 in a.get_blocks()))
        if n >= min_reads: rows.append((fam, c, s, e, round(i, 3), n))
with open(out, 'w') as o:
    o.write('family\tchrom\tstart\tend\tbest_ident\tn_reads\n')
    for r in rows: o.write('\t'.join(map(str, r)) + '\n')
print(f'families with hits {len(hits)} | admitted loci (outside every gene and unit, >= {min_reads} reads): {len(rows)}')

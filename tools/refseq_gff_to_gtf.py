#!/usr/bin/env python3
"""Extract one chromosome from a RefSeq GFF3 and write a gffread-style GTF (transcript + exon rows).

Stand-in for `gffread -T`, which is not installed on this machine. Validated to reproduce the transcript
count of the bakeoff's `chr20_ref.gtf` exactly (4,574 = 4,574), so the bakeoff chromosomes built with it
are comparable to the 2026-09-15 chr20 reference.

usage: refseq_gff_to_gtf.py REFSEQ.gff[.gz] CHROM OUT.gtf
"""
import sys, gzip, collections, re
src, chrom, out = sys.argv[1], sys.argv[2], sys.argv[3]
op = gzip.open if src.endswith('.gz') else open
exons = collections.defaultdict(list); info = {}
def attrs(s):
    d = {}
    for kv in s.rstrip(';').split(';'):
        if '=' in kv: k, v = kv.split('=', 1); d[k.strip()] = v.strip()
    return d
for line in op(src, 'rt'):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[0] != chrom: continue
    a = attrs(f[8])
    if f[2] == 'exon':
        p = a.get('Parent', '')
        for pid in p.split(','):
            if pid: exons[pid].append((int(f[3]), int(f[4]), f[1], f[6]))
    elif 'ID' in a and f[2] not in ('gene', 'pseudogene', 'CDS', 'region'):
        info[a['ID']] = (a.get('Parent', ''), a.get('gene', a.get('Name', '')))
n = 0
with open(out, 'w') as fo:
    for tid, ex in exons.items():
        ex.sort()
        gid, gname = info.get(tid, ('', ''))
        src_f, strand = ex[0][2], ex[0][3]
        at = f'transcript_id "{tid}"; gene_id "{gid}"; gene_name "{gname}"'
        fo.write(f'{chrom}\t{src_f}\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t{strand}\t.\t{at}\n')
        for i, (s, e, _, _) in enumerate(ex, 1):
            fo.write(f'{chrom}\t{src_f}\texon\t{s}\t{e}\t.\t{strand}\t.\t{at}; exon_number "{i}";\n')
        n += 1
print(f'{chrom}: {n} transcripts', file=sys.stderr)
import sys, gzip, collections, re
src, chrom, out = sys.argv[1], sys.argv[2], sys.argv[3]
op = gzip.open if src.endswith('.gz') else open
exons = collections.defaultdict(list); info = {}
def attrs(s):
    d = {}
    for kv in s.rstrip(';').split(';'):
        if '=' in kv: k, v = kv.split('=', 1); d[k.strip()] = v.strip()
    return d
for line in op(src, 'rt'):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[0] != chrom: continue
    a = attrs(f[8])
    if f[2] == 'exon':
        p = a.get('Parent', '')
        for pid in p.split(','):
            if pid: exons[pid].append((int(f[3]), int(f[4]), f[1], f[6]))
    elif 'ID' in a and f[2] not in ('gene', 'pseudogene', 'CDS', 'region'):
        info[a['ID']] = (a.get('Parent', ''), a.get('gene', a.get('Name', '')))
n = 0
with open(out, 'w') as fo:
    for tid, ex in exons.items():
        ex.sort()
        gid, gname = info.get(tid, ('', ''))
        src_f, strand = ex[0][2], ex[0][3]
        at = f'transcript_id "{tid}"; gene_id "{gid}"; gene_name "{gname}"'
        fo.write(f'{chrom}\t{src_f}\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t{strand}\t.\t{at}\n')
        for i, (s, e, _, _) in enumerate(ex, 1):
            fo.write(f'{chrom}\t{src_f}\texon\t{s}\t{e}\t.\t{strand}\t.\t{at}; exon_number "{i}";\n')
        n += 1
print(f'{chrom}: {n} transcripts', file=sys.stderr)

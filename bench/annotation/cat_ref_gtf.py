# CAT/Liftoff slim GFF3 -> one-chromosome GTF for SQANTI3 (copy of gff2gtf.py's output shape; gene_name from the gene line)
import sys, gzip, collections
src, chrom, out = sys.argv[1], sys.argv[2], sys.argv[3]
def attrs(s):
    return dict(kv.split('=', 1) for kv in s.rstrip(';').split(';') if '=' in kv)
exons = collections.defaultdict(list); tx_gene = {}; gname = {}
for line in gzip.open(src, 'rt'):
    if line.startswith('#'): continue
    f = line.rstrip('\n').split('\t')
    if len(f) < 9 or f[0] != chrom: continue
    a = attrs(f[8])
    if f[2] == 'gene': gname[a['ID']] = a.get('gene_name', a['ID'])
    elif f[2] == 'transcript': tx_gene[a['ID']] = a['Parent']
    elif f[2] == 'exon': exons[a['Parent']].append((int(f[3]), int(f[4]), f[1], f[6]))
assert set(exons) <= set(tx_gene), 'exon without transcript'
with open(out, 'w') as fo:
    for tid, ex in exons.items():
        ex.sort(); gid = tx_gene[tid]
        at = f'transcript_id "{tid}"; gene_id "{gid}"; gene_name "{gname[gid]}"'
        fo.write(f'{chrom}\t{ex[0][2]}\ttranscript\t{ex[0][0]}\t{ex[-1][1]}\t.\t{ex[0][3]}\t.\t{at}\n')
        for i, (s, e, _, _) in enumerate(ex, 1):
            fo.write(f'{chrom}\t{ex[0][2]}\texon\t{s}\t{e}\t.\t{ex[0][3]}\t.\t{at}; exon_number "{i}";\n')
print(chrom, len(exons), 'transcripts', len(set(tx_gene[t] for t in exons)), 'genes')

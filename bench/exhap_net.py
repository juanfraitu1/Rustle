"""EXHAP net (LEDGER_O3): per exon unit, the reads whose primary record overlaps one of its loci; at most CAP per unit, the
CAP smallest sha1(read name). usage: exhap_net.py SPLIT (in the exhap directory)"""
import collections, hashlib, sys
import pysam
BAM = '/mnt/linuxdisk/home/juanfraitu/fibroblasts/GCA_029281585.2_flnc_mm.bam'
CAP = 200
split = sys.argv[1]
units = {l.split('\t')[0] for l in open('units.tsv') if l.split('\t')[2] == split}
loci = collections.defaultdict(list)
for l in open('copies.tsv'):
    f = l.split('\t')
    if f[0] in units: loci[f[0]].append((f[3], int(f[4]), int(f[5])))
bam = pysam.AlignmentFile(BAM)
with open(f'net_{split}.tsv', 'w') as o:
    o.write('family\tread\tused\trole\tvia\tclass\tin_reads_tsv\n')
    for u in sorted(loci):
        rs = {r.query_name for c, s, e in loci[u] for r in bam.fetch(c, s, e)
              if not (r.is_secondary or r.is_supplementary) and r.get_overlap(s, e) >= 0.8 * (e - s)}
        for r in sorted(rs, key=lambda x: hashlib.sha1(x.encode()).hexdigest())[:CAP]:
            o.write(f'{u}\t{r}\t1\tNA\tA\tpass_A\t0\n')

"""EXHAP net (bench/exhap_net.py): per exon unit, the reads whose primary record covers >= 0.8 of one of its loci; at most CAP
per unit, the CAP smallest sha1(read name). usage: exhap_net.py DIR BAM K N  (batch K of N by sha1 of the unit)"""
import collections, hashlib, sys
import pysam
D, BAM, K, N = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
CAP = 200
units = {l.split('\t')[0] for l in open(f'{D}/units.tsv') if not l.startswith('unit\t')
         and int(hashlib.sha1(l.split('\t')[0].encode()).hexdigest(), 16) % N == K}
loci = collections.defaultdict(list)
for l in open(f'{D}/copies.tsv'):
    f = l.split('\t')
    if f[0] in units: loci[f[0]].append((f[3], int(f[4]), int(f[5])))
bam = pysam.AlignmentFile(BAM)
with open(f'{D}/net_{K}.tsv.tmp', 'w') as o:
    for u in sorted(loci):
        rs = {r.query_name for c, s, e in loci[u] for r in bam.fetch(c, s, e)
              if not (r.is_secondary or r.is_supplementary) and r.get_overlap(s, e) >= 0.8 * (e - s)}
        for r in sorted(rs, key=lambda x: hashlib.sha1(x.encode()).hexdigest())[:CAP]:
            o.write(f'{u}\t{r}\t1\tNA\tA\tpass_A\t0\n')

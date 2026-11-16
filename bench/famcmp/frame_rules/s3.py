"""s3.py TAG MARGIN : frozen one-read rules of register 1766 addendum on TAG.agree.tsv (new chains with exactly one full primary read, no longer read)."""
import sys
n, margin = sys.argv[1], float(sys.argv[2]); R = {'M1': [0, 0], 'S3': [0, 0], 'S3b': [0, 0]}
with open(n + '.s3.keys', 'w') as w:
    for l in open(n + '.agree.tsv'):
        f = l.rstrip('\n').split('\t')
        if f[0] == 'chrom' or f[14] != '1' or int(f[16]) > 0: continue
        ag, cf, ok = int(f[15]), int(f[17]), int(f[12]); hit = []
        if ag + cf > 0 and ag / (ag + cf) > 0.5: hit.append('M1')
        if f[19] == '0' and f[5] == '0' and int(f[4]) >= 900:
            hit.append('S3')
            if float(f[18]) > 0.05: hit.append('S3b')
        for h in hit: R[h][0] += 1; R[h][1] += ok; w.write('%s\t%s\t%s\n' % (h, f[0], f[13]))
for h, (t, c) in R.items(): print('%s %s: %d new chains, %d correct, %.3f (margin %.3f)' % (n, h, t, c, c / t, margin))

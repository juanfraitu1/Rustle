"""orang_score.py DIR REF_GTF FASTA : the three bars of the S3b default test on one sample (DIR/o.gtf, o.final.gtf, o.oneread.gtf)."""
import sys, re, collections, pysam
sys.path.insert(0, '/home/hpc/v138da/v138da13/Rustle/tools'); sys.path.insert(0, '/home/hpc/v138da/v138da13/Rustle/bench/famcmp')
import gtf_ptc; from oneread import models, chain
d, ref_p, fa_p = sys.argv[1:4]; fa = pysam.FastaFile(fa_p)
ref = {(c, chain(ex)) for c, s, ex, _ in models(ref_p).values() if len(ex) > 1}
D = models(d + '/o.gtf'); dch = {(c, chain(ex)) for c, s, ex, _ in D.values() if len(ex) > 1}; tier = [0, 0]; seen = set()
for c, s, ex, lines in D.values():
    k = (c, chain(ex))
    if len(ex) < 2 or k in seen or s not in '+-': continue
    seen.add(k); m = re.search(r'reads "(\d+)"', lines[0])
    if m and m.group(1) == '2':
        o, p = gtf_ptc.call(ex, s, lambda a, b: fa.fetch(c, a, b))
        if not p and o >= 900: tier[0] += 1; tier[1] += k in ref
A = {(c, chain(ex)) for c, s, ex, _ in models(d + '/o.oneread.gtf').values()}; ac = len(A & ref); dc = len(dch & ref)
F = models(d + '/o.final.gtf'); fall = {(c, chain(ex)) for c, s, ex, _ in F.values() if len(ex) > 1}
fbase = {(c, chain(ex)) for c, s, ex, l in F.values() if len(ex) > 1 and 'oneread "1"' not in l[0]}
p = lambda S: len(S & ref) / len(S)
b1 = ac / len(A) >= tier[1] / tier[0]; b2 = ac >= .02 * dc; b3 = p(fbase) - p(fall) <= .005
print('%s: default %d / %d / %.3f; admitted %d / %d / %.3f vs tier %.3f (%d chains) %s; yield bar %d %s; final.gtf %d / %d / %.3f -> %d / %d / %.3f %s' % (
    d.split('/')[-1], len(dch), dc, dc / len(dch), len(A), ac, ac / len(A), tier[1] / tier[0], tier[0], 'PASS' if b1 else 'FAIL', round(.02 * dc), 'PASS' if b2 else 'FAIL',
    len(fbase), len(fbase & ref), p(fbase), len(fall), len(fall & ref), p(fall), 'PASS' if b3 else 'FAIL'))

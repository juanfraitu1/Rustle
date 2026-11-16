"""join.py TAG REF_GTF DEFAULT_GTF -> TAG.new.tsv : chains of the 1-read-floor run that the default does not output, with every feature and the RefSeq label."""
import sys, collections
sys.path.insert(0, '/home/hpc/v138da/v138da13/Rustle/bench/famcmp'); from allmiss import chains
n, ref_p, def_p = sys.argv[1:4]; key = lambda c, ch: (c, ','.join('%d-%d' % i for i in ch))
ref = {key(c, ch) for c, s, ch in chains(ref_p)}; dflt = {key(c, ch) for c, s, ch in chains(def_p)}
sis = {}; cand = {}
for l in open(n + '.sis.tsv'):
    f = l.rstrip('\n').split('\t'); sis[(f[0], f[1])] = f[2:4]
for l in open(n + '.cand.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'single': cand[(f[1], f[-1])] = f[3:7]  # min_intron locus longer_reads introns_new
seen = set()
with open(n + '.new.tsv', 'w') as w:
    w.write('chrom\treads\tn_introns\ttx_len\torf_len\tptc\tlifted\tsister_reads\tmin_intron\tlocus\tlonger_reads\tintrons_new\tcorrect\tchain\n')
    for l in open(n + '.orf.tsv'):
        f = l.rstrip('\n').split('\t'); k = (f[0], f[1])
        if f[0] == 'chrom' or k in dflt or k in seen: continue
        seen.add(k)
        w.write('\t'.join([f[0], f[3], str(f[1].count(',') + 1), f[4], f[5], f[6]] + sis.get(k, ['0', '0']) + cand.get(k, ['NA'] * 4) + [str(int(k in ref)), f[1]]) + '\n')

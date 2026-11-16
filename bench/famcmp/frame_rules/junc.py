"""junc.py TAG : per junction of every S3 chain, the primary reads on it by class -> TAG.s3.junctions.tsv
(agree = the read's chain is a contiguous part of the chain, longer = contains it, conflict = neither; intron chains only)."""
import sys, pickle, collections
n = sys.argv[1]; prim = pickle.load(open(n + '.t0.reads.pkl', 'rb')); has = collections.defaultdict(list)
for (c, ch), k in prim.items():
    for i in ch: has[(c, i)].append(ch)
def run(a, b):
    if len(a) > len(b) or a[0] not in b: return False
    k = b.index(a[0]); return b[k:k + len(a)] == a
S3 = {tuple(l.rstrip('\n').split('\t')[1:]) for l in open(n + '.s3.keys') if l.startswith('S3\t')}
with open(n + '.s3.junctions.tsv', 'w') as w:
    w.write('chrom\tchain_start\tn_introns\tcorrect\tjunction\tintron\treads\tfull\tagree\tlonger\tconflict\tcolinear_share\n')
    for l in open(n + '.agree.tsv'):
        f = l.rstrip('\n').split('\t')
        if (f[0], f[13]) not in S3: continue
        c = f[0]; C = tuple(tuple(map(int, x.split('-'))) for x in f[13].split(','))
        for j, i in enumerate(C, 1):
            t = collections.Counter()
            for r in has[(c, i)]: t['full' if r == C else 'agree' if run(r, C) else 'longer' if run(C, r) else 'conflict'] += prim[(c, r)]
            tot = sum(t.values()); w.write('\t'.join(map(str, [c, C[0][0], len(C), f[12], j, '%d-%d' % i, tot, t['full'], t['agree'], t['longer'], t['conflict'], '%.3f' % (1 - t['conflict'] / tot)])) + '\n')

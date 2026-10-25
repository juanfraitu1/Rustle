# P18 simulator: reads from both of Jim's alleles of each catalog copy, HiFi errors from badread's pacbio2021 models
import sys, math, random, hashlib, collections
from badread.error_model import ErrorModel
from badread.qscore_model import QScoreModel
from badread.simulate import sequence_fragment
from badread.identities import Identities
SPLIT = sys.argv[1]                      # dev | test
TAG = 'p20'
N = 10                                   # reads per allele (20 per haploid X/Y copy)
p16 = set(open('p16_families.txt').read().split())
def split(fam): return 'dev' if int(hashlib.sha1(f'p18:{fam}'.encode()).hexdigest(), 16) % 3 == 0 else 'test'
def fasta(p):
    s, cur = {}, None
    for l in open(p):
        if l[0] == '>': cur = l[1:].split()[0]; s[cur] = []
        else: s[cur].append(l.strip())
    return {k: ''.join(v).upper() for k, v in s.items()}
ref = fasta('/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.fa')
alle = {'|'.join(k.split('|')[:5]): v for k, v in fasta('alleles.fa').items()}
sex = {'NC_073247.2', 'NC_073248.2'}
idn = Identities(99.8, 0.15, 100.0, output=open("/dev/null", "w")); em = ErrorModel('pacbio2021', output=open('/dev/null', 'w')); qm = QScoreModel('pacbio2021', output=open('/dev/null', 'w'))
rng = random.Random(f'p20:{SPLIT}'); random.seed(f'p20:{SPLIT}')
import numpy as np; np.random.seed(int(hashlib.sha1(f'p20:{SPLIT}'.encode()).hexdigest(), 16) % 2**32)
k = 0; c = collections.Counter(); ex = open(f'p20.{SPLIT}.expr.tsv', 'w')
with open(f'p20.{SPLIT}.reads.fa', 'w') as fa, open(f'p20.{SPLIT}.truth.tsv', 'w') as tr:
    print('read\tfamily\tcopy\tchrom\tstart\tend\tallele', file=tr)
    for q, s in ref.items():
        fam, idx, loc = q.split('|')[:3]
        if fam in p16 or split(fam) != SPLIT: continue
        chrom, se = loc.split(':'); st, en = se.split('-')
        # uneven expression: 30% of copies silent, else lognormal around 10 reads (sd 1.0 in log), capped at 80
        n = 0 if rng.random() < 0.3 else min(80, max(1, round(math.exp(rng.gauss(math.log(10), 1.0)))))
        if q in alle:
            na = sum(rng.random() < 0.5 for _ in range(n)); srcs = [('ref_hap', s)] * (n - na) + [('other_hap', alle[q])] * na
        elif chrom in sex: srcs = [('ref_hap', s)] * n
        else: c['copies skipped (no orthologous allele)'] += 1; continue
        c['copies'] += 1; c['silent copies'] += n == 0; print(fam, idx, n, sep='\t', file=ex); c['families'] += 0
        for a, t in srcs:
            # 5' truncation, as in Iso-Seq: half the reads full length, half keep a uniform 20-100% from the 3' end
            if rng.random() < 0.5: t = t[int(len(t) * (1 - rng.uniform(0.2, 1.0))):]
            ident = idn.get_identity()
            seq, qual, _, _ = sequence_fragment(t, ident, em, qm)
            name = f'p20{SPLIT[0]}{k:06d}'; k += 1
            fa.write(f'>{name}\n{seq}\n'); print(name, fam, idx, chrom, st, en, a, sep='\t', file=tr)
fams = {l.split('\t')[1] for l in open(f'p20.{SPLIT}.truth.tsv')} - {'family'}
print(SPLIT, dict(c), 'families', len(fams), 'reads', k)

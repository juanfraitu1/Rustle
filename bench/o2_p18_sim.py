# P18 simulator: reads from both of Jim's alleles of each catalog copy, HiFi errors from badread's pacbio2021 models
import sys, random, hashlib, collections
from badread.error_model import ErrorModel
from badread.qscore_model import QScoreModel
from badread.simulate import sequence_fragment
from badread.identities import Identities
SPLIT = sys.argv[1]                      # dev | test
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
rng = random.Random(f'p18:{SPLIT}'); random.seed(f'p18:{SPLIT}')
import numpy as np; np.random.seed(int(hashlib.sha1(f'p18:{SPLIT}'.encode()).hexdigest(), 16) % 2**32)
k = 0; c = collections.Counter()
with open(f'{SPLIT}.reads.fa', 'w') as fa, open(f'{SPLIT}.truth.tsv', 'w') as tr:
    print('read\tfamily\tcopy\tchrom\tstart\tend\tallele', file=tr)
    for q, s in ref.items():
        fam, idx, loc = q.split('|')[:3]
        if fam in p16 or split(fam) != SPLIT: continue
        chrom, se = loc.split(':'); st, en = se.split('-')
        if q in alle: srcs = [('ref_hap', s)] * N + [('other_hap', alle[q])] * N
        elif chrom in sex: srcs = [('ref_hap', s)] * (2 * N)
        else: c['copies skipped (no orthologous allele)'] += 1; continue
        c['copies'] += 1; c['families'] += 0
        for a, t in srcs:
            ident = idn.get_identity()
            seq, qual, _, _ = sequence_fragment(t, ident, em, qm)
            name = f'p18{SPLIT[0]}{k:06d}'; k += 1
            fa.write(f'>{name}\n{seq}\n'); print(name, fam, idx, chrom, st, en, a, sep='\t', file=tr)
fams = {l.split('\t')[1] for l in open(f'{SPLIT}.truth.tsv')} - {'family'}
print(SPLIT, dict(c), 'families', len(fams), 'reads', k)

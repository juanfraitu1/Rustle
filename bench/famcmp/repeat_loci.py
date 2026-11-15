#!/usr/bin/env python3
"""PREREG row 1698: REPEAT = exons >= 0.9 soft-masked or >= 0.5 covered by the rDNA unit, per assembly locus.

  repeat_loci.py GTF FASTA LABELS.tsv|- RDNA.fa OUT_PREFIX      (LABELS = sec_only_loci.py loci.tsv, `-` = none; scaffold = contig JAK*)
Writes OUT_PREFIX.tsv (locus, where, label, mask, rdna_cov, repeat) and OUT_PREFIX.fa; driver stage `repeatloci`."""
import sys, os, collections, subprocess, pysam
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sec_only_loci import read_loci
gtf, fa, tsv, rdna, out = sys.argv[1:6]
loci = read_loci(gtf); f = pysam.FastaFile(fa)
flag = {} if tsv == '-' else {(r[0], int(r[1]), int(r[2])): r[-1] for r in (l.rstrip('\n').split('\t') for l in list(open(tsv))[1:])}
rows = []
with open(out + '.fa', 'w') as w:
    for c, v in loci.items():
        if c not in f.references: continue
        for s, e, st, ts in v:
            t = max(ts, key=lambda t: sum(b - a for a, b in t)); seq = ''.join(f.fetch(c, a, b) for a, b in t)
            name = f'{c}:{s}-{e}'; w.write(f'>{name}\n{seq}\n')
            rows.append((name, 'scaffold' if c.startswith('JAK') else 'chromosome', flag.get((c, s, e), 'not flagged'), sum(x.islower() for x in seq) / max(1, len(seq))))
tool = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'tools', 'sequence_families.py')
cov = dict(l.split('\t') for l in subprocess.run([sys.executable, tool, out + '.fa', rdna, '--per-copy'], check=True, capture_output=True, text=True).stdout.splitlines())
g = collections.defaultdict(lambda: [0, 0, 0, 0])
with open(out + '.tsv', 'w') as w:
    w.write('locus\twhere\tlabel\tmask\trdna_cov\trepeat\n')
    for name, where, lab, m in rows:
        r = float(cov.get(name, 0)); w.write(f'{name}\t{where}\t{lab}\t{m:.3f}\t{r:.3f}\t{int(m >= 0.9 or r >= 0.5)}\n')
        x = g[(where, lab)]; x[0] += 1; x[1] += m >= 0.9; x[2] += r >= 0.5; x[3] += (m >= 0.9 or r >= 0.5)
print('where\tlabel\tn\tmask>=0.9\trdna>=0.5\tREPEAT\tshare')
for k, (n, a, b, c) in sorted(g.items()): print(*k, n, a, b, c, f'{c / n:.3f}', sep='\t')
print(f'{sum(v[3] for v in g.values())} of {sum(v[0] for v in g.values())} loci marked repeat (exons >= 0.9 soft-masked or >= 0.5 covered by {os.path.basename(rdna)})')

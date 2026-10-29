#!/usr/bin/env python3
"""mmseqs families (2026-10-30, user "what about trying mmseqs on the family members"): write every locus's SPLICED exon sequence,
cluster with `mmseqs easy-cluster` (nucleotide, --cov-mode 1 = coverage of the shorter sequence), families = clusters of >= 2 loci,
written in the famtruth/pairscore copies format.
  mmseqs_families.py LOCI_GFF3 GENOME_FA OUT [--min-seq-id 0.8] [--cov 0.5] [--cluster-mode 0] [--threads 16]"""
import sys, os, re, collections, subprocess
a = sys.argv[1:]; gff, fa, out = a[:3]
opt = lambda k, d: a[a.index(k) + 1] if k in a else d
ID = opt('--min-seq-id', '0.8'); COV = opt('--cov', '0.5'); CM = opt('--cluster-mode', '0'); THR = opt('--threads', '16')
loci = {}; exons = collections.defaultdict(list)
for l in open(gff):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t')
    if f[2] == 'gene': loci[re.search(r'ID=gene-([^;]+)', f[8]).group(1)] = (f[0], f[6] if f[6] in '+-' else '+')
    elif f[2] == 'exon': exons[re.search(r'gene=([^;]+)', f[8]).group(1)].append((int(f[3]) - 1, int(f[4])))
by = collections.defaultdict(list)
for t, (c, s) in loci.items(): by[c].append(t)
comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
if not os.path.exists(out + '.spliced.fa'):
    with open(out + '.spliced.fa', 'w') as w:
        for chrom, ts in by.items():
            seq = ''.join(subprocess.run(['samtools', 'faidx', fa, chrom], capture_output=True, text=True).stdout.split('\n')[1:])
            for t in ts:
                s = ''.join(seq[x:y] for x, y in sorted(exons[t]))
                if loci[t][1] == '-': s = s.translate(comp)[::-1]
                if len(s) >= 50: w.write(f">{t}\n{s}\n")
tmp = out + '.mmtmp'
subprocess.run(['mmseqs', 'easy-cluster', out + '.spliced.fa', out + '.mm', tmp, '--min-seq-id', ID, '-c', COV, '--cov-mode', '1', '--cluster-mode', CM, '--threads', THR, '-v', '0'], check=True)
subprocess.run(['rm', '-rf', tmp])
fams = collections.defaultdict(list)
for l in open(out + '.mm_cluster.tsv'):
    r, m = l.split(); fams[r].append(m)
fams = sorted([v for v in fams.values() if len(v) >= 2], key=lambda v: -len(v))
with open(out + '.fam.copies.tsv', 'w') as w:
    w.write("family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\n")
    for i, v in enumerate(fams):
        for j, t in enumerate(sorted(v)):
            c, s = loci[t]; ex = sorted(exons[t])
            w.write(f"MM{i}\t{j}\t{t}\t{c}\t{ex[0][0]+1}\t{ex[-1][1]}\t{len(ex)}\t{s}\t0\t{','.join(f'{p+1}-{q}' for p, q in ex)}\n")
print(f"mmseqs id {ID} cov {COV} mode {CM}: families {len(fams)}, members {sum(len(v) for v in fams)}, largest {[len(v) for v in fams[:5]]}", file=sys.stderr)

#!/usr/bin/env python3
"""O2 read-level truth simulation (PREREG_o2_read_truth_2026-09-23): reads from every copy of every multi-copy
catalog family (the catalog's own spliced sequence), HiFi errors, end jitter, 5'/3' trimming; mapped to the
whole genome with the shipped minimap2 settings. Also writes per-copy closest-sibling identity.
usage: o2_read_truth_sim.py copies.tsv copies.fa INDEX.mmi OUT_PREFIX SEED"""
import sys, random, subprocess, collections, os
sys.path.insert(0, '/mnt/c/Users/jfris/Desktop/Rustle/bench'); from sim_reads import simulate_reads
tsv, fa, idx, out, seed = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
rng = random.Random(seed)
# copies.tsv: family_id copy_idx tid chrom start end n_exon strand n_reads
rows = [l.rstrip('\n').split('\t') for l in open(tsv)][1:]
nfam = collections.Counter(r[0] for r in rows)
seqs = {}
name = None
for l in open(fa):
    if l.startswith('>'):
        h = l[1:].split('|'); name = (h[0], h[1]); seqs[name] = []
    else:
        seqs[name].append(l.strip())
seqs = {k: ''.join(v).upper() for k, v in seqs.items()}
n = 0; ncopies = 0
with open(out + '.fq', 'w') as fq, open(out + '.copies_used.tsv', 'w') as cu:
    cu.write('family_id\tcopy_idx\tchrom\tstart\tend\tn_reads_real\tn_sim\tlen\n')
    for r in rows:
        fam, ci = r[0], r[1]
        if nfam[fam] < 2: continue
        s = seqs.get((fam, ci))
        if not s or len(s) < 300: continue
        k = min(100, max(10, int(r[8])))
        ncopies += 1
        cu.write(f'{fam}\t{ci}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[8]}\t{k}\t{len(s)}\n')
        for i, (rd, q) in enumerate(simulate_reads(s, k, err=0.001, indel=0.0003, seed=seed * 131 + hash(fam + ci) % 1000003, trunc_frac=0.10)):
            j5, j3 = rng.randint(0, 30), rng.randint(0, 30)
            rd = rd[j5:len(rd) - j3]; q = q[j5:len(q) - j3]
            fq.write(f'@{fam}|{ci}|{i}\n{rd}\n+\n{q}\n'); n += 1
print(f'[simO2] {ncopies} copies, {n} reads', flush=True)
subprocess.run(f"minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 4 {idx} {out}.fq 2> {out}.mm2.log | samtools sort -@2 -o {out}.bam - && samtools index {out}.bam", shell=True, check=True)
# closest sibling identity per copy: all-vs-all of the multi-copy families' sequences
with open(out + '.copies.fa', 'w') as f:
    for (fam, ci), s in seqs.items():
        if nfam[fam] >= 2: f.write(f'>{fam}|{ci}\n{s}\n')
subprocess.run(f"minimap2 -x asm20 -c -X -N 200 -p 0 -t 4 {out}.copies.fa {out}.copies.fa > {out}.ava.paf 2>/dev/null", shell=True, check=True)
best = {}
for l in open(out + '.ava.paf'):
    f = l.split('\t'); a, b = f[0], f[5]
    if a == b or a.split('|')[0] != b.split('|')[0]: continue
    idn = int(f[9]) / int(f[10]); cov = (int(f[3]) - int(f[2])) / int(f[1])
    if cov < 0.5: continue
    for x in (a, b):
        if idn > best.get(x, (0, ''))[0]: best[x] = (idn, b if x == a else a)
with open(out + '.sibling.tsv', 'w') as f:
    f.write('family_id\tcopy_idx\tclosest_identity\tclosest_copy\n')
    for (fam, ci) in seqs:
        if nfam[fam] >= 2:
            idn, nb = best.get(f'{fam}|{ci}', (float('nan'), 'NA'))
            f.write(f'{fam}\t{ci}\t{idn:.5f}\t{nb}\n')
print(f'[simO2] wrote {out}.bam, {out}.sibling.tsv ({len(best)} copies with a sibling hit)')

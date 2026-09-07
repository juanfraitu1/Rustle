#!/usr/bin/env python3
"""PREREG docs/PREREG_npip_known_origin_2026-09-06.md — one sensitivity / one specificity for NPIP on reads of
known origin: reads simulated from every NPIP unit's spliced sequence, aligned like the substrate, assigned by
copy_assign with its defaults, scored by the read name. Subcategory = primary MAPQ < 60 (and as_margin = 0).
usage: npip_sim.py <family_dir> <out_dir> [reads_per_unit=200] [trunc_frac=0.3]"""
import sys, os, csv, subprocess, collections, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_reads import simulate_reads, write_fastq
import pysam
fam, out = sys.argv[1], sys.argv[2]; N = int(sys.argv[3]) if len(sys.argv) > 3 else 200; TR = float(sys.argv[4]) if len(sys.argv) > 4 else 0.3
ERR = float(sys.argv[5]) if len(sys.argv) > 5 else 0.003; INDEL = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0008
FA = '/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3_contigs.fa'; BIN = '/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign'
os.makedirs(out, exist_ok=True)
cp = {r['copy_idx']: r for r in csv.DictReader(open(f'{fam}/copies.tsv'), delimiter='\t')}
seqs = {}; name = None
for l in open(f'{fam}/copies.fa'):
    if l.startswith('>'): name = l[1:].split('|')[1]; seqs[name] = []
    else: seqs[name].append(l.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}
fq = f'{out}/sim.fq'
if not os.path.exists(f'{out}/sim.bam'):
    with open(fq, 'w') as fh:
        for ci in sorted(seqs, key=int):
            for i, rq in enumerate(simulate_reads(seqs[ci], N, err=ERR, indel=INDEL, seed=int(ci) + 1, trunc_frac=TR)):
                if len(rq[0]) >= 300: write_fastq(fh, f'SIMNPIP|{ci}|{i}', rq)
    subprocess.run(f'minimap2 -ax splice:hq -uf --eqx -Y -N 50 -p 0.1 --secondary=yes -t 4 {FA} {fq} 2>/dev/null | samtools sort -o {out}/sim.bam - && samtools index {out}/sim.bam', shell=True, check=True)
    for f in ('copies.tsv', 'copies.fa', 'regions', 'forecast.tsv'): subprocess.run(['cp', f'{fam}/{f}', out])
if not os.path.exists(f'{out}/A.assignments.tsv'):
    subprocess.run(f'ulimit -v 12000000; {BIN} --bam {out}/sim.bam --fasta {FA} --families {out}/copies.tsv --copies-fa {out}/copies.fa --regions {out}/regions --out {out}/A > {out}/A.log 2> {out}/A.err', shell=True)
A = {r['read_name']: r for r in csv.DictReader(open(f'{out}/A.assignments.tsv'), delimiter='\t')}
B = pysam.AlignmentFile(f'{out}/sim.bam')
spans = {i: (r['chrom'], int(r['start']), int(r['end'])) for i, r in cp.items()}
prim = {}; mapq = {}; nsec = collections.Counter(); ntot = 0
for a in B.fetch(until_eof=True):
    if a.is_unmapped: continue
    if a.is_secondary: nsec[a.query_name] += 1; continue
    if a.is_supplementary: continue
    best = None
    for i, (c, s, e) in spans.items():
        if a.reference_name != c: continue
        o = sum(min(b1, e) - max(b0, s) for b0, b1 in a.get_blocks() if b1 > s and b0 < e)
        if o > 0 and (best is None or o > best[1]): best = (i, o)
    prim[a.query_name] = best[0] if best else None; mapq[a.query_name] = a.mapping_quality
names = sorted(prim); ntot = len(names)
def tab(sel, label):
    n = len(sel); pr = sum(1 for x in sel if prim[x] == x.split('|')[1]); pw = sum(1 for x in sel if prim[x] is not None and prim[x] != x.split('|')[1])
    st = collections.Counter(A[x]['status'] if x in A else 'absent' for x in sel)
    right = sum(1 for x in sel if x in A and A[x]['status'] == 'assigned' and A[x]['catalog_copy_idx'] == x.split('|')[1])
    wrong = [x for x in sel if x in A and A[x]['status'] == 'assigned' and A[x]['catalog_copy_idx'] != x.split('|')[1]]
    sole = sum(1 for x in sel if x in A and A[x].get('sole_candidate') == '1')
    print(f"{label}: reads {n} | minimap2 primary: right {pr} ({100*pr/max(n,1):.1f} %) wrong {pw} ({100*pw/max(n,1):.1f} %) | O2: right {right} ({100*right/max(n,1):.1f} %) wrong {len(wrong)} ({100*len(wrong)/max(n,1):.2f} %) tied {st['tied']} ambiguous {st['ambiguous']} absent {st['absent']} | O2 precision {right/max(1,right+len(wrong)):.4f} sensitivity {right/max(n,1):.4f} specificity {1-len(wrong)/max(n,1):.4f} | sole {sole}")
    return wrong
print(f'simulated reads with a primary alignment: {ntot}; with >=1 secondary: {sum(1 for x in names if nsec[x])}; in O2 output: {sum(1 for x in names if x in A)}')
tab(names, 'ALL')
low = [x for x in names if mapq[x] < 60]; w = tab(low, 'SUBCATEGORY primary MAPQ<60')
q0 = [x for x in names if mapq[x] == 0]; tab(q0, '  of which MAPQ 0')
ties = [x for x in names if x in A and A[x]['as_second'] not in ('NA', '') and A[x]['as_margin'] == '0']; tab(ties, 'SUBCATEGORY as_margin = 0 (equal best scores)')
tab([x for x in names if mapq[x] >= 60], 'unique (MAPQ 60)')
wc = collections.Counter((x.split('|')[1], A[x]['catalog_copy_idx']) for x in w)
print('wrong calls in the subcategory by (true unit -> called unit):', wc.most_common(10))
ident = {i: r.get('nearest_ident', 'NA') for i, r in csv.DictReader(open(f'{fam}/forecast.tsv'), delimiter='\t') and {} or {}} if False else None
# per true unit: sensitivity and wrong calls
per = collections.defaultdict(lambda: [0, 0, 0, 0])
for x in names:
    t = x.split('|')[1]; per[t][0] += 1
    if x in A and A[x]['status'] == 'assigned': per[t][1 if A[x]['catalog_copy_idx'] == t else 2] += 1
    elif x in A and A[x]['status'] == 'tied': per[t][3] += 1
with open(f'{out}/per_unit.tsv', 'w') as o:
    o.write('unit\tmember_status\tn\tright\twrong\ttied\n')
    for t in sorted(per, key=int): o.write(f"{t}\t{cp[t]['member_status']}\t{per[t][0]}\t{per[t][1]}\t{per[t][2]}\t{per[t][3]}\n")
print('per-unit table:', f'{out}/per_unit.tsv', '| units with any wrong call:', [(t, per[t][2]) for t in sorted(per, key=int) if per[t][2]])

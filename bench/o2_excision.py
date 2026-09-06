#!/usr/bin/env python3
"""PREREG adj/excise: remove copy X from a family, rerun copy_assign (genomic read-star default), follow X's
MAPQ-60 reads, and look for the missing-copy signature: consistent mismatch sites among the origin-rejected reads
that share a best candidate; reconstruct the copy from them.  usage: o2_excision.py <fam_dir> <X> <out_dir>"""
import sys, os, csv, subprocess, collections, shutil, pysam, re
fam_dir, X, out = sys.argv[1], sys.argv[2], sys.argv[3]
BAM_P = '/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3.bam'; FA_P = '/mnt/linuxdisk/home/juanfraitu/npip_cat/npip3_contigs.fa'
B = pysam.AlignmentFile(BAM_P); FA = pysam.FastaFile(FA_P); BIN = '/mnt/linuxdisk/home/juanfraitu/rustle_target/release/copy_assign'
os.makedirs(out, exist_ok=True)
cp = list(csv.DictReader(open(f'{fam_dir}/copies.tsv'), delimiter='\t')); fid = cp[0]['family_id']
keep = [r for r in cp if r['copy_idx'] != X]; remap = {r['copy_idx']: str(i) for i, r in enumerate(keep)}
fa = {}; name = None
for l in open(f'{fam_dir}/copies.fa'):
    if l.startswith('>'): name = l[1:].split('|')[1]; fa[name] = [l]
    else: fa[name].append(l)
with open(f'{out}/copies.tsv', 'w') as o:
    w = csv.DictWriter(o, fieldnames=list(cp[0].keys()), delimiter='\t'); w.writeheader()
    for r in keep: r2 = dict(r); r2['copy_idx'] = remap[r['copy_idx']]; w.writerow(r2)
with open(f'{out}/copies.fa', 'w') as o:
    for r in keep:
        h = fa[r['copy_idx']][0].split('|'); h[1] = remap[r['copy_idx']]; o.write('|'.join(h)); o.write(''.join(fa[r['copy_idx']][1:]))
shutil.copy(f'{fam_dir}/regions', out)
with open(f'{out}/forecast.tsv', 'w') as o:
    for l in open(f'{fam_dir}/forecast.tsv'):
        f = l.split('\t')
        if f[0] == 'copy_idx': o.write(l)
        elif f[0] in remap: o.write('\t'.join([remap[f[0]]] + f[1:]))
if not os.path.exists(f'{out}/A.assignments.tsv'):
    subprocess.run(f'ulimit -v 10000000; {BIN} --bam {BAM_P} --fasta {FA_P} --families {out}/copies.tsv --copies-fa {out}/copies.fa --regions {out}/regions --out {out}/A > {out}/A.log 2> {out}/A.err', shell=True)
A = {r['read_name']: r for r in csv.DictReader(open(f'{out}/A.assignments.tsv'), delimiter='\t')}
A0 = {r['read_name']: r for r in csv.DictReader(open(f'{fam_dir}/A.assignments.tsv'), delimiter='\t')} if os.path.exists(f'{fam_dir}/A.assignments.tsv') else {}
xr = [r for r in cp if r['copy_idx'] == X][0]; xc, xs, xe = xr['chrom'], int(xr['start']), int(xr['end'])
xreads = {}
for a in B.fetch(xc, xs, xe):
    if a.flag & 2308 or a.mapping_quality < 60: continue
    if any(b1 > xs and b0 < xe for b0, b1 in a.get_blocks()): xreads[a.query_name] = a
fate = collections.Counter(); best_of_rej = collections.Counter()
for n in xreads:
    r = A.get(n)
    if r is None: fate['absent'] += 1; continue
    if r['status'] == 'assigned': fate['SILENT misassignment'] += 1
    elif r.get('origin_rejected') == '1': fate['origin-rejected'] += 1; best_of_rej[r['catalog_copy_idx']] += 1
    else: fate[r['status']] += 1
n = len(xreads)
print(f'== {fid} excise copy {X} ({xc}:{xs}-{xe}, nearest_ident {[l.split()[1] for l in open(fam_dir+"/forecast.tsv") if l.split()[0]==X][0]}): {n} MAPQ-60 reads')
print('   fate:', {k: f'{v} ({100*v/n:.0f}%)' for k, v in fate.items()}, '| rejected reads\' best candidate (excised index):', dict(best_of_rej.most_common(3)))
# ---- detector: consistent mismatch sites among origin-rejected reads sharing a best candidate (excised run) vs controls
def consistent_sites(reads, y, tag):
    """align reads to Y's padded locus; return (n_reads, covered_kb, consistent_sites, top allele per site)"""
    yc, ys, ye = y['chrom'], int(y['start']), int(y['end']); pad = max((len(a.query_sequence) for a in reads), default=0)
    ls, le = max(0, ys - pad), ye + pad
    open(f'{out}/{tag}_Y.fa', 'w').write(f'>Y\n{FA.fetch(yc, ls, le)}\n'); open(f'{out}/{tag}_reads.fa', 'w').write(''.join(f'>{a.query_name}\n{a.query_sequence}\n' for a in reads))
    paf = subprocess.run(['minimap2', '-x', 'splice', '-c', '--eqx', '-N', '1', '-t', '4', f'{out}/{tag}_Y.fa', f'{out}/{tag}_reads.fa'], capture_output=True, text=True).stdout
    mism = collections.defaultdict(collections.Counter); cov = collections.Counter(); qseq = {a.query_name: a.query_sequence for a in reads}
    for l in paf.splitlines():
        f = l.split('\t'); cg = [t for t in f[12:] if t.startswith('cg:Z:')][0][5:]; t = int(f[7]); q = int(f[2]) if f[4] == '+' else int(f[1]) - int(f[3]); s = qseq[f[0]]
        if f[4] == '-': s = s.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]
        for num, op in re.findall(r'(\d+)([=XIDNS])', cg):
            num = int(num)
            if op == '=': 
                for k in range(num): cov[t + k] += 1
                t += num; q += num
            elif op == 'X':
                for k in range(num): cov[t + k] += 1; mism[t + k][s[q + k].upper()] += 1
                t += num; q += num
            elif op in 'DN': t += num
            elif op in 'IS': q += num
    cons = {p: c.most_common(1)[0] for p, c in mism.items() if sum(c.values()) >= 3 and sum(c.values()) >= 0.5 * cov[p]}
    covered = sorted(p for p, c in cov.items() if c >= 3); covered_kb = len(covered) / 1000
    return len(reads), covered_kb, cons, (ls, covered)
for yidx, cnt in best_of_rej.most_common(1):
    y = keep[int(yidx)]
    rej = [xreads[n] for n in xreads if A.get(n, {}).get('origin_rejected') == '1' and A[n]['catalog_copy_idx'] == yidx]
    nr, kb, cons, (ls, covered) = consistent_sites(rej, y, 'exc')
    dens = len(cons) / kb if kb else 0
    # control 1: Y's own MAPQ-60 reads vs Y ; control 2: rejected reads of the UN-excised run with best Y (if any)
    yc, ys, ye = y['chrom'], int(y['start']), int(y['end']); yreads = [a for a in B.fetch(yc, ys, ye) if not a.flag & 2308 and a.mapping_quality >= 60 and any(b1 > ys and b0 < ye for b0, b1 in a.get_blocks())][:200]
    nr1, kb1, cons1, _ = consistent_sites(yreads, y, 'ctl1') if yreads else (0, 0, {}, 0)
    orig_idx = [k for k, v in remap.items() if v == yidx][0]
    rej0 = [a for a in B.fetch(xc - 0 if False else yc, max(0, ys - 200000), ye + 200000) if not a.flag & 2308 and A0.get(a.query_name, {}).get('origin_rejected') == '1' and A0[a.query_name]['catalog_copy_idx'] == orig_idx][:200]
    nr2, kb2, cons2, _ = consistent_sites(rej0, y, 'ctl2') if rej0 else (0, 0, {}, 0)
    print(f'   detector at Y=excised idx {yidx} (orig {orig_idx}): rejected reads {nr}, covered {kb:.1f} kb, consistent sites {len(cons)} = {dens:.1f}/kb | control Y-own reads: {nr1} reads, {len(cons1)/kb1 if kb1 else 0:.1f}/kb | control un-excised rejected@Y: {nr2} reads, {len(cons2)/kb2 if kb2 else 0:.1f}/kb')
    # recovery: patch Y with the consistent alleles, compare to X and to Y
    if cons and covered:
        # reconstruct ONLY the covered stretch (the largest run of read-covered positions), Y's bases with the
        # consistent sites patched by the majority read allele; score it against X's locus and Y's locus
        runs = []; 
        for p in covered:
            if runs and p - runs[-1][1] <= 50: runs[-1][1] = p
            else: runs.append([p, p])
        lo, hi = max(runs, key=lambda r: r[1] - r[0]); lo += ls; hi += ls + 1
        seq = list(FA.fetch(yc, lo, hi).upper())
        for p, (allele, _) in cons.items():
            q = p + ls - lo
            if 0 <= q < len(seq): seq[q] = allele
        yorig = FA.fetch(yc, lo, hi).upper()
        open(f'{out}/reconstructed.fa', 'w').write('>reconstructed\n' + ''.join(seq) + '\n>Y_unpatched\n' + yorig + '\n')
        pad = max((len(a.query_sequence) for a in rej), default=0)
        open(f'{out}/targets.fa', 'w').write(f'>X_true\n{FA.fetch(xc, max(0, xs - pad), xe + pad)}\n>Y_locus\n{FA.fetch(yc, max(0, ys - pad), ye + pad)}\n')
        paf = subprocess.run(['minimap2', '-x', 'asm20', '-c', '-N', '5', '-p', '0', '-t', '4', f'{out}/targets.fa', f'{out}/reconstructed.fa'], capture_output=True, text=True).stdout
        best = {}
        for l in paf.splitlines():
            f = l.split('\t'); ident = int(f[9]) / max(1, int(f[10])); k = (f[0], f[5])
            if k not in best or int(f[10]) > best[k][1]: best[k] = (ident, int(f[10]))
        g = lambda q, t: best.get((q, t), (0, 0))
        print(f'   recovery over the covered {hi-lo} bp: reconstructed -> X_true {g("reconstructed","X_true")[0]:.4f} | -> Y {g("reconstructed","Y_locus")[0]:.4f} ; unpatched Y -> X_true {g("Y_unpatched","X_true")[0]:.4f} | -> Y {g("Y_unpatched","Y_locus")[0]:.4f}')

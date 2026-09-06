#!/usr/bin/env python3
"""O3 item 2 (PREREG docs/PREREG_o3_reconstruction_2026-09-06.md): for every `missing_copy` pair of an O3 pass,
patch Y's locus with the consistent alleles over the covered stretch (§6ff recovery, verbatim) and align the
reconstruction to the whole assembly; verdict existing_locus / reference_absent / unresolved.

usage: o3_reconstruct.py <pass_dir> <sweep_dir> --bam B --fasta F --genome-index IDX --units U --gff G
       [--only-reads FILE --x chrom:start-end]   (verification on an excision directory)
Writes <pass_dir>/recon/*.fa, <pass_dir>/recon.tsv.
"""
import sys, os, csv, re, collections, subprocess, tempfile, pysam
args = sys.argv[1:]; pass_dir, sweep = args[0], args[1]
def opt(k, d=None): return args[args.index(k) + 1] if k in args else d
BAM = pysam.AlignmentFile(opt('--bam')); FA = pysam.FastaFile(opt('--fasta')); IDX = opt('--genome-index'); only = set(open(opt('--only-reads')).read().split()) if '--only-reads' in args else None
xspan = None
if '--x' in args:
    c, r = opt('--x').split(':'); a, b = r.split('-'); xspan = (c, int(a), int(b))
genes = collections.defaultdict(list)
for l in open(opt('--gff')):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('gene', 'pseudogene'): continue
    a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv); genes[f[0]].append((int(f[3]) - 1, int(f[4]), a.get('gene', '?')))
units = collections.defaultdict(list)
for r in csv.DictReader(open(opt('--units')), delimiter='\t'): units[r['chrom']].append((int(r['start']), int(r['end']), r['family_id'], r['copy_idx'], r.get('member_status', 'NA')))
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
flags = [p for p in csv.DictReader(open(f'{pass_dir}/flags.tsv'), delimiter='\t') if p['flag'] == 'missing_copy']
rd = f'{pass_dir}/recon'; os.makedirs(rd, exist_ok=True)
rows = []; qmeta = {}
allq = open(f'{rd}/all.fa', 'w')
by_fam = collections.defaultdict(list)
for p in flags: by_fam[p['fam']].append(p)
for fam, ps in by_fam.items():
    d = f'{sweep}/{fam}' if os.path.isdir(f'{sweep}/{fam}') else sweep
    cp = {r['copy_idx']: r for r in csv.DictReader(open(f'{d}/copies.tsv'), delimiter='\t')}
    A = {r['read_name']: r for r in csv.DictReader(open(f'{d}/A.assignments.tsv'), delimiter='\t')}
    truth = {}; rec = {}
    for i, r in cp.items():
        c, s, e = r['chrom'], int(r['start']), int(r['end'])
        for a in BAM.fetch(c, s, e):
            if a.flag & 2308 or a.query_name not in A: continue
            o = sum(min(b1, e) - max(b0, s) for b0, b1 in a.get_blocks() if b1 > s and b0 < e)
            if o > 0 and (a.query_name not in truth or o > truth[a.query_name][1]): truth[a.query_name] = (i, o); rec[a.query_name] = a
    if only is not None:
        for c in {r['chrom'] for r in cp.values()}:
            for a in BAM.fetch(c):
                if not a.flag & 2308 and a.query_name in only and a.query_name in A and a.query_name not in truth: truth[a.query_name] = ('X', 0); rec[a.query_name] = a
    for p in ps:
        y = p['Y']; yr = cp[y]; yc, ys, ye = yr['chrom'], int(yr['start']), int(yr['end'])
        yl = (int(yr['locus_start']), int(yr['locus_end'])) if yr.get('locus_start', 'NA') not in ('NA', '') else None
        names = sorted(n for n in truth if A[n].get('origin_rejected') == '1' and A[n].get('n_candidates', '1') != '0' and A[n]['catalog_copy_idx'] == y and (only is None or n in only))[:500]
        reads = [rec[n] for n in names]
        pad = max((len(a.query_sequence) for a in reads), default=0)
        ls, le = (max(0, ys - pad), ye + pad) if yl is None else (min(yl[0], ys), max(yl[1], ye))
        with tempfile.TemporaryDirectory(dir=rd) as tmp:
            open(f'{tmp}/Y.fa', 'w').write(f'>Y\n{FA.fetch(yc, ls, le)}\n'); open(f'{tmp}/r.fa', 'w').write(''.join(f'>{a.query_name}\n{a.query_sequence}\n' for a in reads))
            paf = subprocess.run(['minimap2', '-x', 'splice', '-c', '--eqx', '-N', '1', '-t', '4', f'{tmp}/Y.fa', f'{tmp}/r.fa'], capture_output=True, text=True).stdout
        mism = collections.defaultdict(collections.Counter); cov = collections.Counter(); qseq = {a.query_name: a.query_sequence for a in reads}
        for l in paf.splitlines():
            f = l.split('\t'); cg = [t for t in f[12:] if t.startswith('cg:Z:')][0][5:]; t = int(f[7]); q = int(f[2]) if f[4] == '+' else int(f[1]) - int(f[3]); s = qseq[f[0]]
            if f[4] == '-': s = s.translate(COMP)[::-1]
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
        cons = {pp: c.most_common(1)[0] for pp, c in mism.items() if sum(c.values()) >= 3 and sum(c.values()) >= 0.5 * cov[pp]}
        covered = sorted(pp for pp, c in cov.items() if c >= 3)
        runs = []
        for pp in covered:
            if runs and pp - runs[-1][1] <= 50: runs[-1][1] = pp
            else: runs.append([pp, pp])
        if not runs: rows.append([fam, y, yr.get('member_status', 'NA'), len(reads), 0, 0, 'nan', 'unresolved', '-', '-', '-', '-']); continue
        lo, hi = max(runs, key=lambda r: r[1] - r[0]); glo, ghi = lo + ls, hi + ls + 1
        seq = list(FA.fetch(yc, glo, ghi).upper()); n_patched = 0
        for pp, (allele, _) in cons.items():
            q = pp + ls - glo
            if 0 <= q < len(seq): seq[q] = allele; n_patched += 1
        ident_y = 1 - n_patched / max(1, len(seq))
        qid = f'{fam}|{y}'; qmeta[qid] = (fam, y, yr.get('member_status', 'NA'), len(reads), len(seq), n_patched, ident_y, yc, min(ls, ys), max(le, ye))
        if n_patched < 3 or len(seq) < 300: rows.append([fam, y, yr.get('member_status', 'NA'), len(reads), len(seq), n_patched, f'{ident_y:.4f}', 'unresolved', '-', '-', '-', '-']); continue
        allq.write(f'>{qid}\n{"".join(seq)}\n'); open(f'{rd}/{fam}_{y}.fa', 'w').write(f'>{qid} {yc}:{glo}-{ghi} patched={n_patched}\n{"".join(seq)}\n')
allq.close()
# one alignment of every reconstruction against the assembly
paf = subprocess.run(['minimap2', '-x', 'asm20', '-c', '-N', '20', '-p', '0', '--secondary=yes', '-t', '4', IDX, f'{rd}/all.fa'], capture_output=True, text=True).stdout
hits = collections.defaultdict(list)
for l in paf.splitlines():
    f = l.split('\t'); qid = f[0]; qlen = int(f[1]); ident = int(f[9]) / max(1, int(f[10])); span = int(f[3]) - int(f[2])
    hits[qid].append((f[5], int(f[7]), int(f[8]), ident, span / qlen))
for qid, (fam, y, st, nr, L, npch, iy, yc, ylo, yhi) in qmeta.items():
    if npch < 3 or L < 300: continue
    outside = [h for h in hits[qid] if h[4] >= 0.5 and not (h[0] == yc and h[1] < yhi and ylo < h[2])]
    inside = [h for h in hits[qid] if h[0] == yc and h[1] < yhi and ylo < h[2]]
    # identity to Y measured the same way as every other hit (aligned, indels included): the best hit inside
    # Y's extent; the patched fraction is only the fallback (§6ff compared patched→X with patched→Y as alignments)
    # §6ff's rule on both sides: only alignments covering ≥ half the reconstruction count; inside Y's locus the
    # LONGEST such alignment gives the identity to Y (a shorter, purer sub-hit would flatter Y)
    inside = [h for h in inside if h[4] >= 0.5]
    if inside: iy = max(inside, key=lambda h: h[4])[3]
    best = max(outside, key=lambda h: h[3]) if outside else None
    if best:
        c, s, e, bi, cov = best
        u = [x for x in units.get(c, []) if x[0] < e and s < x[1]]; g = [x for x in genes.get(c, []) if x[0] < e and s < x[1]]
        where = ';'.join(f'{x[2]}:{x[3]}:{x[4]}' for x in u[:2]) or '-'; gw = ';'.join(x[2] for x in g[:3]) or '-'
        verdict = 'existing_locus' if bi >= iy else 'reference_absent'
        # the true X in verification mode
        xhit = ''
        if xspan: xhit = 'X_hit' if (c == xspan[0] and s < xspan[2] and xspan[1] < e) else 'not_X'
        rows.append([fam, y, st, nr, L, npch, f'{iy:.4f}', verdict, f'{c}:{s}-{e}', f'{bi:.4f}', where, gw + (f' [{xhit}]' if xhit else '')])
    else:
        rows.append([fam, y, st, nr, L, npch, f'{iy:.4f}', 'reference_absent', '-', '-', '-', '-'])
with open(f'{pass_dir}/recon.tsv', 'w') as o:
    o.write('fam\tY\tY_status\tn_rejected\tstretch_bp\tsites_patched\tident_to_Y\tverdict\tbest_hit_outside\thit_ident\tunit_there\tgenes_there\n')
    for r in rows: o.write('\t'.join(map(str, r)) + '\n')
vc = collections.Counter(r[7] for r in rows)
print(f'flags {len(flags)} | {dict(vc)} | reconstructions aligned {len(qmeta)}')

#!/usr/bin/env python3
"""O3 flag pass (PREREG docs/PREREG_o3_flag_pass_2026-09-06.md): every family's origin-rejected unit reads through
the missing-copy detector (§6ff, verbatim) and the admission prototype (§6fc), one flag table per sweep.

usage: o3_flag_pass.py <out_dir> --bam B --fasta F --gff G --units catalog.units.tsv (--sweep DIR | --dirs d1 d2 ...)
       [--budget SEC] [--max-reads 500] [--alpha 0.001]
Per family a cached <out_dir>/<fam>.pairs.tsv / .loci.tsv (skipped when present); the aggregate step (always run
at the end) applies the flag rule over ALL pairs present in <out_dir>: p < alpha / n_pairs.
"""
import sys, os, csv, glob, math, time, collections, subprocess, re, tempfile, pysam
args = sys.argv[1:]; out = args[0]
def opt(k, d=None):
    return args[args.index(k) + 1] if k in args else d
bam_p, fa_p, gff_p, units_p = opt('--bam'), opt('--fasta'), opt('--gff'), opt('--units')
budget = float(opt('--budget', 1e9)); max_reads = int(opt('--max-reads', 500)); alpha = float(opt('--alpha', 0.001))
only = set(open(opt('--only-reads')).read().split()) if '--only-reads' in args else None  # verification: restrict the rejected set
if '--dirs' in args:
    i = args.index('--dirs'); dirs = [a for a in args[i + 1:] if not a.startswith('--')]
else:
    dirs = sorted(glob.glob(f"{opt('--sweep')}/fam_*"))
os.makedirs(out, exist_ok=True)
BAM = pysam.AlignmentFile(bam_p); FA = pysam.FastaFile(fa_p)
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
# RefSeq genes + every catalog unit (for the other-family class)
genes = collections.defaultdict(list)
for l in open(gff_p):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('gene', 'pseudogene'): continue
    a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
    genes[f[0]].append((int(f[3]) - 1, int(f[4]), a.get('gene', '?'), a.get('gene_biotype', '?')))
units = collections.defaultdict(list)
for r in csv.DictReader(open(units_p), delimiter='\t'):
    units[r['chrom']].append((int(r['start']), int(r['end']), r['family_id'], r['copy_idx'], r.get('member_status', 'NA')))
def poisson_tail(k, lam):
    """P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0: return 1.0
    if lam <= 0: return 0.0
    s = 0.0
    for i in range(k):
        s += math.exp(-lam + i * math.log(lam) - math.lgamma(i + 1))
    return max(0.0, min(1.0, 1.0 - s))
def detector(reads, y, tmp):
    """§6ff verbatim: (n_reads, covered_kb, consistent_sites{pos:(allele,n)}, per-read (mismatch, unaligned))."""
    yc, ys, ye, yl = y; pad = max((len(a.query_sequence) for a in reads), default=0)
    # target = O2's own target: the catalog locus extent (L2) when the copies.tsv carries it, else unit ± longest read (§6ff)
    ls, le = (max(0, ys - pad), ye + pad) if yl is None else (min(yl[0], ys), max(yl[1], ye))
    open(f'{tmp}/Y.fa', 'w').write(f'>Y\n{FA.fetch(yc, ls, le)}\n')
    open(f'{tmp}/reads.fa', 'w').write(''.join(f'>{a.query_name}\n{a.query_sequence}\n' for a in reads))
    paf = subprocess.run(['minimap2', '-x', 'splice', '-c', '--eqx', '-N', '1', '-t', '4', f'{tmp}/Y.fa', f'{tmp}/reads.fa'], capture_output=True, text=True).stdout
    mism = collections.defaultdict(collections.Counter); cov = collections.Counter(); qseq = {a.query_name: a.query_sequence for a in reads}
    per_read = {}
    for l in paf.splitlines():
        f = l.split('\t'); cg = [t for t in f[12:] if t.startswith('cg:Z:')][0][5:]; t = int(f[7]); q = int(f[2]) if f[4] == '+' else int(f[1]) - int(f[3]); s = qseq[f[0]]
        if f[4] == '-': s = s.translate(COMP)[::-1]
        nx = 0
        for num, op in re.findall(r'(\d+)([=XIDNS])', cg):
            num = int(num)
            if op == '=':
                for k in range(num): cov[t + k] += 1
                t += num; q += num
            elif op == 'X':
                nx += num
                for k in range(num): cov[t + k] += 1; mism[t + k][s[q + k].upper()] += 1
                t += num; q += num
            elif op in 'DN': t += num
            elif op in 'IS': q += num
        unal = int(f[1]) - (int(f[3]) - int(f[2]))
        if f[0] not in per_read or (int(f[3]) - int(f[2])) > per_read[f[0]][2]: per_read[f[0]] = (nx, unal, int(f[3]) - int(f[2]))
    cons = {p: c.most_common(1)[0] for p, c in mism.items() if sum(c.values()) >= 3 and sum(c.values()) >= 0.5 * cov[p]}
    covered_kb = sum(1 for p, c in cov.items() if c >= 3) / 1000
    return len(reads), covered_kb, cons, per_read
def median(v):
    v = sorted(v); return v[len(v) // 2] if v else float('nan')
t0 = time.time(); done = 0
for d in dirs:
    fam = os.path.basename(d)
    if os.path.exists(f'{out}/{fam}.pairs.tsv'): continue
    if time.time() - t0 > budget: break
    if not os.path.exists(f'{d}/A.assignments.tsv'): continue
    cp = {r['copy_idx']: r for r in csv.DictReader(open(f'{d}/copies.tsv'), delimiter='\t')}
    A = {r['read_name']: r for r in csv.DictReader(open(f'{d}/A.assignments.tsv'), delimiter='\t')}
    truth = {}; mapq = {}; rec = {}
    for i, r in cp.items():
        c, s, e = r['chrom'], int(r['start']), int(r['end'])
        for a in BAM.fetch(c, s, e):
            if a.flag & 2308 or a.query_name not in A: continue
            o = sum(min(b1, e) - max(b0, s) for b0, b1 in a.get_blocks() if b1 > s and b0 < e)
            if o > 0 and (a.query_name not in truth or o > truth[a.query_name][1]):
                truth[a.query_name] = (i, o); mapq[a.query_name] = a.mapping_quality; rec[a.query_name] = a
    if only is not None:
        # verification mode: the named reads are followed wherever their primary lies (an excised copy is no unit)
        for c in {r['chrom'] for r in cp.values()}:
            for a in BAM.fetch(c):
                if not a.flag & 2308 and a.query_name in only and a.query_name in A and a.query_name not in truth:
                    truth[a.query_name] = ('X', 0); mapq[a.query_name] = a.mapping_quality; rec[a.query_name] = a
    rej = collections.defaultdict(list); orphans = []
    for n in truth:
        r = A[n]
        if r.get('n_candidates', '1') == '0': orphans.append(n)  # pre-§6fg outputs have no column
        elif r['origin_rejected'] == '1' and (only is None or n in only): rej[r['catalog_copy_idx']].append(n)
    pairs = []
    with tempfile.TemporaryDirectory(dir=out) as tmp:
        for y, names in sorted(rej.items(), key=lambda kv: -len(kv[1])):
            if len(names) < 3 or y not in cp: continue
            yr = cp[y]; yloc = (yr['chrom'], int(yr['start']), int(yr['end']), (int(yr['locus_start']), int(yr['locus_end'])) if yr.get('locus_start', 'NA') not in ('NA', '', None) else None)
            reads = [rec[n] for n in sorted(names)[:max_reads]]
            nr, kb, cons, pr = detector(reads, yloc, tmp)
            ctl_names = sorted(n for n in truth if truth[n][0] == y and A[n]['origin_rejected'] != '1' and A[n]['catalog_copy_idx'] == y)[:max_reads]  # control = Y's own reads the certificate accepted at Y (any MAPQ, any status: NPIP copies have few MAPQ-60 reads; pre-§6fj outputs label sole candidates `tied`)
            if len(ctl_names) >= 3:
                nc, kbc, consc, _ = detector([rec[n] for n in ctl_names], yloc, tmp)
            else:
                nc, kbc, consc = len(ctl_names), 0.0, {}
            rate = len(cons) / kb if kb else 0.0
            rate_c = (len(consc) / kbc) if kbc else float('nan')
            lam_rate = (max(len(consc), 1) / kbc) if kbc else float('nan')   # floored at 1 site over the control's covered kb
            p = poisson_tail(len(cons), lam_rate * kb) if kb and kbc else float('nan')
            mx = median([v[0] for v in pr.values()]); mu = median([v[1] for v in pr.values()])
            cls = 'divergent' if (not math.isnan(mx) and mx > mu) else 'structural'
            pairs.append([fam, y, yr.get('member_status', 'NA'), len(names), nr, f'{kb:.2f}', len(cons), f'{rate:.2f}', nc, f'{kbc:.2f}', len(consc), f'{rate_c:.2f}', f'{p:.3e}', cls, f'{mx:.0f}', f'{mu:.0f}', len(pr)])
    # admission: rejected reads whose primary has no block in any unit of the family, plus the orphans
    fam_units = [(r['chrom'], int(r['start']), int(r['end'])) for r in cp.values()]
    # admission input = EVERY rejected or orphan row of the family (a unit read's primary is inside a unit by
    # definition; the reads with a primary elsewhere reached O2 through a secondary record at a copy, §6fc)
    targets = {n for n, r in A.items() if r.get('n_candidates', '1') == '0' or r['origin_rejected'] == '1'}
    windows = collections.defaultdict(list)
    for c, s, e in fam_units: windows[c].append((max(0, s - 200000), e + 200000))
    prim = {}
    for c, ws in windows.items():
        ws.sort(); merged = []
        for lo, hi in ws:
            if merged and lo <= merged[-1][1]: merged[-1][1] = max(merged[-1][1], hi)
            else: merged.append([lo, hi])
        for lo, hi in merged:
            for a in BAM.fetch(c, lo, hi):
                if a.flag & 2308 or a.query_name not in targets: continue
                inside = any(a.reference_name == uc and any(b1 > s and b0 < e for b0, b1 in a.get_blocks()) for uc, s, e in fam_units)
                if not inside: prim[a.query_name] = a
    loci = []
    for a in sorted(prim.values(), key=lambda a: (a.reference_name, a.reference_start)):
        if loci and loci[-1][0] == a.reference_name and a.reference_start - loci[-1][2] <= 5000:
            loci[-1][2] = max(loci[-1][2], a.reference_end); loci[-1][3].append(a)
        else: loci.append([a.reference_name, a.reference_start, a.reference_end, [a]])
    locus_rows = []
    for c, s, e, reads in loci:
        if len(reads) < 3: continue
        g = [x for x in genes.get(c, []) if x[0] < e and s < x[1]]
        other = [u for u in units.get(c, []) if u[0] < e and s < u[1] and u[2] != cp[next(iter(cp))]['family_id'].split('_')[0]]
        cls = 'other_family' if other else ('annotated_no_unit' if g else 'unannotated')
        n_orph = sum(1 for a in reads if A[a.query_name].get('n_candidates', '1') == '0')
        locus_rows.append([fam, c, s, e, len(reads), n_orph, cls, ';'.join(f'{x[2]}({x[3]})' for x in g[:3]) or '-', ';'.join(f'{u[2]}:{u[3]}:{u[4]}' for u in other[:3]) or '-'])
    with open(f'{out}/{fam}.pairs.tsv', 'w') as o:
        o.write('fam\tY\tY_status\tn_rejected\tn_aligned\tcovered_kb\tsites\trate_per_kb\tctl_reads\tctl_kb\tctl_sites\tctl_rate\tp_poisson\tclass\tmed_mismatch\tmed_unaligned\tn_paf\n')
        for p in pairs: o.write('\t'.join(map(str, p)) + '\n')
    with open(f'{out}/{fam}.loci.tsv', 'w') as o:
        o.write('fam\tchrom\tstart\tend\tn_reads\tn_orphans\tclass\tgenes\tother_family_units\n')
        for r in locus_rows: o.write('\t'.join(map(str, r)) + '\n')
    with open(f'{out}/{fam}.counts.tsv', 'w') as o:
        o.write(f'fam\tunit_reads\trejected\torphans\trejected_rows_all\tprimaries_outside\n{fam}\t{len(truth)}\t{sum(len(v) for v in rej.values())}\t{len(orphans)}\t{len(targets)}\t{len(prim)}\n')
    done += 1
# ---- aggregate over everything present in out
P = []; L = []; C = []
for f in sorted(glob.glob(f'{out}/*.pairs.tsv')): P += list(csv.DictReader(open(f), delimiter='\t'))
for f in sorted(glob.glob(f'{out}/*.loci.tsv')): L += list(csv.DictReader(open(f), delimiter='\t'))
for f in sorted(glob.glob(f'{out}/*.counts.tsv')): C += list(csv.DictReader(open(f), delimiter='\t'))
n_pairs = sum(1 for p in P if p['p_poisson'] != 'nan'); thr = alpha / max(n_pairs, 1)
for p in P:
    p['flag'] = 'missing_copy' if p['p_poisson'] != 'nan' and float(p['p_poisson']) < thr else ('untestable' if p['p_poisson'] == 'nan' else 'none')
    p['line_5kb'] = '1' if float(p['rate_per_kb']) >= 5 else '0'
with open(f'{out}/flags.tsv', 'w') as o:
    keys = list(P[0].keys()) if P else []
    o.write('\t'.join(keys) + '\n')
    for p in P: o.write('\t'.join(p[k] for k in keys) + '\n')
with open(f'{out}/loci.tsv', 'w') as o:
    keys = list(L[0].keys()) if L else []
    o.write('\t'.join(keys) + '\n')
    for r in L: o.write('\t'.join(r[k] for k in keys) + '\n')
fl = [p for p in P if p['flag'] == 'missing_copy']; ctl_lt1 = sum(1 for p in P if p['ctl_rate'] not in ('nan',) and float(p['ctl_rate']) < 1)
cls = collections.Counter((p['class'], p['flag']) for p in P); lcls = collections.Counter(r['class'] for r in L)
msg = (f"families {len(C)} (this call {done}, {time.time()-t0:.0f} s) | unit reads {sum(int(c['unit_reads']) for c in C)} rejected {sum(int(c['rejected']) for c in C)} orphans {sum(int(c['orphans']) for c in C)} | rejected/orphan rows of any kind {sum(int(c['rejected_rows_all']) for c in C)}, primaries outside every unit {sum(int(c['primaries_outside']) for c in C)}\n"
       f"pairs tested {n_pairs} (threshold {thr:.2e}) | flagged missing_copy {len(fl)} ({100*len(fl)/max(n_pairs,1):.1f} %) | over the 5/kb line {sum(1 for p in P if p['line_5kb']=='1')} | control rate < 1/kb: {ctl_lt1}/{sum(1 for p in P if p['ctl_rate']!='nan')}\n"
       f"class x flag: {dict(cls)}\n"
       f"candidate loci {len(L)}: {dict(lcls)} | reads at them {sum(int(r['n_reads']) for r in L)} (orphans {sum(int(r['n_orphans']) for r in L)})")
open(f'{out}/summary.txt', 'w').write(msg + '\n'); print(msg)

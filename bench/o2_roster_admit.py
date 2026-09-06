#!/usr/bin/env python3
"""§6fc roster admission by read origin: for each sweep family, the origin-rejected unit reads whose primary lies
outside every unit are clustered into loci (>= min_reads within gap bp); each locus becomes an added copy
(`read_admitted`) with an exon chain from those primaries (unit rule: >= min_reads covering AND more covering
than splicing over), sequence from the genome (transcription orientation by the reads' ts/strand majority).
usage: o2_roster_admit.py <sweep_in> <sweep_out> <bam> <genome.fa> [min_reads=3] [gap=5000]"""
import sys, os, glob, csv, collections, shutil, pysam
sin, sout, bam_p, fa_p = sys.argv[1:5]; min_reads = int(sys.argv[5]) if len(sys.argv) > 5 else 3; gap = int(sys.argv[6]) if len(sys.argv) > 6 else 5000
BAM = pysam.AlignmentFile(bam_p); FA = pysam.FastaFile(fa_p)
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
def chain(reads):
    """exon chain from primary records: coverage >= min_reads and > reads splicing over."""
    lo = min(a.reference_start for a in reads); hi = max(a.reference_end for a in reads)
    cov = [0] * (hi - lo); spl = [0] * (hi - lo)
    for a in reads:
        blocks = a.get_blocks()
        for s, e in blocks:
            for x in range(s, e): cov[x - lo] += 1
        for (s1, e1), (s2, e2) in zip(blocks, blocks[1:]):
            for x in range(e1, s2): spl[x - lo] += 1
    ex = []
    for k in range(hi - lo):
        if cov[k] >= min_reads and cov[k] > spl[k]:
            x = lo + k
            if ex and ex[-1][1] == x: ex[-1][1] = x + 1
            else: ex.append([x, x + 1])
    return [(s, e) for s, e in ex]
summary = []
for d in sorted(glob.glob(f'{sin}/fam_*')):
    fam = os.path.basename(d); od = f'{sout}/{fam}'; os.makedirs(od, exist_ok=True)
    for f in ('copies.tsv', 'copies.fa', 'regions', 'forecast.tsv'): shutil.copy(f'{d}/{f}', od)
    if not os.path.exists(f'{d}/A.assignments.tsv'): continue
    cp = list(csv.DictReader(open(f'{d}/copies.tsv'), delimiter='\t')); fid = cp[0]['family_id']
    spans = [(r['chrom'], int(r['start']), int(r['end'])) for r in cp]
    A = {r['read_name']: r for r in csv.DictReader(open(f'{d}/A.assignments.tsv'), delimiter='\t')}
    rej = {n for n, r in A.items() if r.get('origin_rejected') == '1'}
    chrom = spans[0][0]; lo = min(s for _, s, _ in spans) - 200000; hi = max(e for _, _, e in spans) + 200000
    prim = {}
    for a in BAM.fetch(chrom, max(0, lo), hi):
        if a.flag & 2308 or a.query_name not in rej: continue
        inside = any(any(b1 > s and b0 < e for b0, b1 in a.get_blocks()) for _, s, e in spans)
        if not inside: prim[a.query_name] = a
    outside = sorted(prim.values(), key=lambda a: a.reference_start)
    loci = []
    for a in outside:
        if loci and a.reference_start - loci[-1][1] <= gap: loci[-1][1] = max(loci[-1][1], a.reference_end); loci[-1][2].append(a)
        else: loci.append([a.reference_start, a.reference_end, [a]])
    loci = [l for l in loci if len(l[2]) >= min_reads]
    n_added = 0
    if loci:
        rows = cp[:]; fa_lines = open(f'{d}/copies.fa').read()
        idx = len(cp)
        for s, e, reads in loci:
            ex = chain(reads)
            if not ex: continue
            strand = collections.Counter()
            for a in reads:
                ts = a.get_tag('ts') if a.has_tag('ts') else None
                st = ts if ts in ('+', '-') else '+'
                if a.is_reverse and ts in ('+', '-'): st = '-' if ts == '+' else '+'
                strand[st] += 1
            st = strand.most_common(1)[0][0]
            seq = ''.join(FA.fetch(chrom, a_, b_) for a_, b_ in ex).upper()
            if st == '-': seq = seq.translate(COMP)[::-1]
            us, ue = ex[0][0], ex[-1][1]
            rows.append({'family_id': fid, 'copy_idx': str(idx), 'tid': f'ADM_{chrom}_{us}', 'chrom': chrom, 'start': str(us), 'end': str(ue), 'n_exon': str(len(ex)), 'strand': st, 'n_reads': str(len(reads)), 'exons': ','.join(f'{a_}-{b_}' for a_, b_ in ex), 'core_hull': 'NA'})
            fa_lines += f'>{fid}|{idx}|{chrom}:{us}-{ue}|{st}|nexon={len(ex)}\n{seq}\n'; idx += 1; n_added += 1
        with open(f'{od}/copies.tsv', 'w') as o:
            w = csv.DictWriter(o, fieldnames=list(cp[0].keys()), delimiter='\t'); w.writeheader()
            for r in rows: w.writerow({k: r.get(k, 'NA') for k in cp[0].keys()})
        open(f'{od}/copies.fa', 'w').write(fa_lines)
        with open(f'{od}/forecast.tsv', 'a') as o:
            for k in range(len(cp), idx): o.write(f'{k}\tNA\tNA\tNA\tNA\tread_admitted\n')
        rlo = min(min(s for _, s, _ in spans), min(l[0] for l in loci)) - 5000; rhi = max(max(e for _, _, e in spans), max(l[1] for l in loci)) + 5000
        open(f'{od}/regions', 'w').write(f'{chrom}:{max(1, rlo)}-{rhi}\n')
    summary.append((fam, len(cp), len(rej), len(outside), n_added))
n_fam = sum(1 for r in summary if r[4] > 0)
print(f'families {len(summary)} | with admitted loci {n_fam} | admitted copies {sum(r[4] for r in summary)} | families dir: {sout}')
with open(f'{sout}/admitted_summary.tsv', 'w') as o:
    o.write('fam\tn_copies\torigin_rejected\trejected_primary_outside\tadmitted\n')
    for r in summary: o.write('\t'.join(map(str, r)) + '\n')

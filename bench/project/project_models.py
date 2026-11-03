#!/usr/bin/env python3
"""PREREG r1519 (docs/PREREG_project_models_2026-11-03.md): project each family's multi-exon transcript models onto its other
members (structure borrowed, reads not pooled) and tier every projected chain by the reads at the target.

  project_models.py GTF COPIES_TSV FASTA BAM OUTDIR [--targets chr13,chr20] [--threads N] [--min-cov 0.9] [--min-id 0.85] [--as-table TSV]

OUTDIR/templates.fa, targets.fa, aln.sam, projected.gtf (attrs: tier own|tied|none, template, from_copy, cov, ident, motif),
tier_{own,tied,none}.gtf, base.gtf (the arm's own transcripts at the target copies, the SQANTI control), summary.txt.
Tiers at the target copy: own = >= 1 PRIMARY record with the exact projected intron chain; tied = none primary but >= 1 SECONDARY
record with the exact chain at AS >= 0.98 x the molecule's genome-wide best (--as-table; Amendment 1); none = neither.
A chain already present among the target copy's own models is skipped (counted `already`). Single-exon templates are out (v1)."""
import sys, os, re, csv, subprocess, collections
import pysam

def read_gtf(path):
    ex = collections.defaultdict(list); meta = {}
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('exon', 'transcript'): continue
        a = dict(re.findall(r'(\S+) "([^"]*)"', f[8])); t = a['transcript_id']
        if f[2] == 'exon': ex[t].append((int(f[3]) - 1, int(f[4])))
        else: meta[t] = (f[0], f[6], int(a.get('reads', 0)), a.get('gene_id', ''))
    tx = {}
    for t, e in ex.items():
        e.sort(); c, s, r, g = meta[t]
        tx[t] = dict(chrom=c, strand=s, exons=e, reads=r, gene=g, introns=tuple((a[1], b[0]) for a, b in zip(e, e[1:])))
    return tx

def chain_of(rec):
    pos = rec.reference_start; ints = []
    for op, l in rec.cigartuples:
        if op == 3: ints.append((pos, pos + l))
        if op in (0, 2, 3, 7, 8): pos += l
    return tuple(ints)

def tier_target(gn):
    g = _gmeta[gn]; m = g['m']; want = _want[gn]; d = {}
    bam = pysam.AlignmentFile(_bam_path)
    for r in bam.fetch(m['chrom'], g['off'], g['off'] + _glen[gn]):
        if r.is_unmapped or r.is_supplementary or 'N' not in r.cigarstring: continue
        c = chain_of(r)
        if c not in want: continue
        e = d.setdefault(c, [0, 0])
        if not r.is_secondary: e[0] += 1
        else:
            b = _best.get(r.query_name)
            if b and r.has_tag('AS') and r.get_tag('AS') >= 0.98 * b: e[1] += 1
    return gn, {c: tuple(v) for c, v in d.items()}

def main():
    args = sys.argv[1:]; opt = {'--targets': None, '--threads': '16', '--min-cov': '0.9', '--min-id': '0.85', '--as-table': os.path.expanduser('~/runs/human_A119b/a1.molecules.tsv')}
    pos = []
    while args:
        a = args.pop(0)
        if a in opt: opt[a] = args.pop(0)
        else: pos.append(a)
    gtf, copies_tsv, fasta, bam_path, out = pos; os.makedirs(out, exist_ok=True)
    reuse = os.environ.get('PROJECT_REUSE_ALN') and os.path.exists(f'{out}/aln.sam')  # profiling: skip the alignment
    if os.path.exists(f'{out}/aln.sam') and not reuse: os.remove(f'{out}/aln.sam')
    targets_ok = set(opt['--targets'].split(',')) if opt['--targets'] else None
    tx = read_gtf(gtf)
    fams = collections.defaultdict(list)
    for r in csv.DictReader(open(copies_tsv), delimiter='\t'):
        fams[r['family_id']].append(dict(idx=r['copy_idx'], chrom=r['chrom'], start=int(r['locus_start']), end=int(r['locus_end']),
                                        strand=r['strand'], reads=int(r['n_reads']), gene=r['gene_id']))
    fams = {f: m for f, m in fams.items() if len(m) >= 2 and (not targets_ok or any(c['chrom'] in targets_ok for c in m))}
    # transcripts at each copy: same strand, overlapping the locus span
    by_cs = collections.defaultdict(list)
    for t, v in tx.items(): by_cs[(v['chrom'], v['strand'])].append(t)
    for f, ms in fams.items():
        for m in ms:
            m['tx'] = [t for t in by_cs[(m['chrom'], m['strand'])] if tx[t]['exons'][0][0] < m['end'] and m['start'] < tx[t]['exons'][-1][1]]
            m['chains'] = {tx[t]['introns'] for t in m['tx'] if tx[t]['introns']}
    fa = pysam.FastaFile(fasta)
    # templates: unique multi-exon chains per family, spliced sequence in transcript orientation
    tmeta = {}; n_t = 0
    with open(f'{out}/templates.fa', 'w') as w:
        for f, ms in fams.items():
            seen = set()
            for m in ms:
                for t in m['tx']:
                    v = tx[t]
                    if not v['introns'] or v['introns'] in seen: continue
                    seen.add(v['introns'])
                    s = ''.join(fa.fetch(v['chrom'], a, b) for a, b in v['exons']).upper()
                    if v['strand'] == '-': s = s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                    name = f'T{n_t}'; n_t += 1
                    tmeta[name] = dict(fam=f, tid=t, copy=m['idx'], span=v['exons'][-1][1] - v['exons'][0][0], len=len(s), n_int=len(v['introns']))
                    w.write(f'>{name}\n{s}\n')
    # targets: every member on a target contig, padded by the family's longest template span; one minimap2 run PER FAMILY
    # (templates vs that family's targets only, -N = members), so the work is sum over families, not templates x all targets
    gmeta = {}; glen = {}; tfam = collections.defaultdict(list)
    for n, t in tmeta.items(): tfam[t['fam']].append(n)
    tseq = {}; cur = None
    for l in open(f'{out}/templates.fa'):
        if l.startswith('>'): cur = l[1:].strip()
        else: tseq[cur] = l.strip()
    os.makedirs(f'{out}/fam', exist_ok=True)
    jobs = []
    for f, ms in sorted(fams.items(), key=lambda kv: -len(kv[1])):
        if not tfam[f]: continue
        pad = 5000  # a wider window (the longest template span) swallowed neighbouring copies: self-hits and chimeric chains
        tg = [m for m in ms if not targets_ok or m['chrom'] in targets_ok]
        if not tg: continue
        with open(f'{out}/fam/{f}.targets.fa', 'w') as w:
            for m in tg:
                s0 = max(0, m['start'] - pad); e0 = min(fa.get_reference_length(m['chrom']), m['end'] + pad)
                name = f'G{len(gmeta)}'; gmeta[name] = dict(fam=f, m=m, off=s0); glen[name] = e0 - s0
                w.write(f'>{name}\n{fa.fetch(m["chrom"], s0, e0).upper()}\n')
        with open(f'{out}/fam/{f}.templates.fa', 'w') as w:
            for n in tfam[f]: w.write(f'>{n}\n{tseq[n]}\n')
        jobs.append(f"minimap2 -ax splice:hq -uf --eqx -N {max(2, len(tg))} -p 0.1 --secondary=yes -f 2000 -t {max(1, int(opt['--threads']) // 6)} {out}/fam/{f}.targets.fa {out}/fam/{f}.templates.fa 2>> {out}/fam/{f}.mm2.log | grep -v '^@' > {out}/fam/{f}.sam")
    # the chr13 mega-families take ~2 min and 12 GB each: 6 at a time (the node has 256 GB), biggest first
    from multiprocessing.pool import ThreadPool
    if not reuse:
        ThreadPool(6).map(lambda c: subprocess.run(c, shell=True, check=True), jobs)
        subprocess.run(f"cat {out}/fam/*.sam > {out}/aln.sam && cat {out}/fam/*.mm2.log > {out}/mm2.log", shell=True, check=True)
    with open(f'{out}/targets.fa', 'w') as w:
        for f in fams:
            p = f'{out}/fam/{f}.targets.fa'
            if os.path.exists(p): w.write(open(p).read())
    if not reuse: subprocess.run(f"samtools faidx {out}/targets.fa && (samtools dict {out}/targets.fa | grep -v '^@HD'; cat {out}/aln.sam) > {out}/aln.h.sam && mv {out}/aln.h.sam {out}/aln.sam", shell=True, check=True)
    # best alignment per (template, target) within the template's family
    best = {}
    for a in pysam.AlignmentFile(f'{out}/aln.sam', 'r', check_sq=False):
        if a.is_unmapped or a.is_supplementary: continue
        g = gmeta[a.reference_name]; t = tmeta[a.query_name]
        if g['fam'] != t['fam'] or g['m']['idx'] == t['copy']: continue
        AS = a.get_tag('AS')
        if (a.query_name, a.reference_name) not in best or AS > best[(a.query_name, a.reference_name)][0]:
            best[(a.query_name, a.reference_name)] = (AS, a)
    CANON = {('GT', 'AG'), ('GC', 'AG'), ('AT', 'AC')}
    best_as = {}
    with open(opt['--as-table']) as fh:  # molecules with > 1 record: name -> genome-wide best AS (the good pool's tie set)
        for l in fh:
            f = l.rstrip('\n').split('\t')
            if len(f) > 3 and f[0][0] != '#' and f[3].isdigit() and int(f[3]) > 1 and f[1].isdigit(): best_as[f[0]] = int(f[1])
    rows = []; cnt = collections.Counter(); motif_n = motif_ok = 0; gained = collections.defaultdict(set); want = {}
    for (tn, gn), (AS, a) in best.items():
        t = tmeta[tn]; g = gmeta[gn]; m = g['m']
        qlen = a.query_length or t['len']; cov = a.query_alignment_length / t['len']
        ident = 1 - a.get_tag('NM') / max(a.query_alignment_length, 1)
        strand = '-' if a.is_reverse else '+'
        if strand != m['strand']: cnt['wrong_strand'] += 1; continue
        if cov < float(opt['--min-cov']): cnt['low_cov'] += 1; continue
        if ident < float(opt['--min-id']): cnt['low_id'] += 1; continue
        pos = a.reference_start; exons = []; es = pos
        for op, l in a.cigartuples:
            if op == 3: exons.append((es, pos)); pos += l; es = pos
            elif op in (0, 2, 7, 8): pos += l
        exons.append((es, pos))
        exons = [(x + g['off'], y + g['off']) for x, y in exons]
        if len(exons) < 2: cnt['unspliced_at_target'] += 1; continue
        if exons[0][0] < m['start'] - 5000 or exons[-1][1] > m['end'] + 5000: cnt['outside_locus'] += 1; continue
        introns = tuple((x[1], y[0]) for x, y in zip(exons, exons[1:]))
        if introns in m['chains']: cnt['already'] += 1; continue
        ok = 0
        for d, ac in introns:
            dn = fa.fetch(m['chrom'], d, d + 2).upper(); acc = fa.fetch(m['chrom'], ac - 2, ac).upper()
            if strand == '-':
                rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                dn, acc = rc(acc), rc(dn)
            ok += (dn, acc) in CANON
        motif_n += len(introns); motif_ok += ok
        want.setdefault(gn, set()).add(introns)
        rows.append([m['chrom'], strand, exons, None, tn, t['tid'], t['copy'], g['fam'], m['idx'], cov, ident, ok, len(introns), 0, 0, gn, introns])
    # tiers: one BAM pass per target copy in worker processes (the chr13 hotspot copies carry ~1M records each), keeping only
    # the projected chains. tied = a SECONDARY record with the exact chain whose AS >= 0.98 x the molecule's genome-wide best
    # (the good pool's tie: the molecule's best placement is a near-identical copy, i.e. another member).
    from multiprocessing import Pool
    global _bam_path, _want, _gmeta, _glen, _best
    _bam_path, _want, _gmeta, _glen, _best = bam_path, want, gmeta, glen, best_as
    with Pool(16) as pool:
        tiers = dict(pool.imap_unordered(tier_target, sorted(want), chunksize=1))
    for row in rows:
        gn, introns = row[15], row[16]; n_pri, n_tied = tiers[gn].get(introns, (0, 0))
        tier = 'own' if n_pri else 'tied' if n_tied else 'none'
        row[3], row[13], row[14] = tier, n_pri, n_tied
        cnt[tier] += 1
        if gmeta[gn]['m']['reads'] <= 2: gained[(gmeta[gn]['fam'], gmeta[gn]['m']['idx'])].add(tier)
    rows = [tuple(r[:15]) for r in rows]
    ws = {k: open(f'{out}/tier_{k}.gtf', 'w') for k in ('own', 'tied', 'none')}
    with open(f'{out}/projected.gtf', 'w') as w:
        for i, (c, s, ex, tier, tn, tid, fc, f, mi, cov, ident, ok, ni, npri, nsec) in enumerate(rows):
            pid = f'PJ{i}_{f}_{mi}'
            attr = f'gene_id "{f}_{mi}"; transcript_id "{pid}"; tier "{tier}"; template "{tid}"; from_copy "{fc}"; cov "{cov:.3f}"; ident "{ident:.3f}"; motif "{ok}/{ni}"; n_primary "{npri}"; n_tied "{nsec}";'
            lines = [f'{c}\tproject\ttranscript\t{ex[0][0] + 1}\t{ex[-1][1]}\t.\t{s}\t.\t{attr}\n'] + [f'{c}\tproject\texon\t{x + 1}\t{y}\t.\t{s}\t.\tgene_id "{f}_{mi}"; transcript_id "{pid}";\n' for x, y in ex]
            w.writelines(lines); ws[tier].writelines(lines)
    for h in ws.values(): h.close()
    # base control: the arm's own transcripts at the target copies
    with open(f'{out}/base.gtf', 'w') as w, open(f'{out}/base_tids.txt', 'w') as w2:
        tids = {t for g in gmeta.values() for t in g['m']['tx']}
        w2.write('\n'.join(sorted(tids)) + '\n')
        for l in open(gtf):
            f = l.split('\t')
            if len(f) > 8 and re.search(r'transcript_id "([^"]+)"', f[8]).group(1) in tids: w.write(l)
    low = [(f, m['idx']) for f, ms in fams.items() for m in ms if m['reads'] <= 2 and (not targets_ok or m['chrom'] in targets_ok)]
    with open(f'{out}/summary.txt', 'w') as w:
        print(f"families {len(fams)} templates {len(tmeta)} targets {len(gmeta)} pairs_aligned {len(best)}", file=w)
        print(f"filters: {dict(cnt)}", file=w)
        print(f"motif canonical {motif_ok}/{motif_n} = {motif_ok / max(motif_n, 1):.4f} (bar >= 0.95)", file=w)
        print(f"reach: low-read target copies (<= 2 reads) {len(low)}, gaining >= 1 projected chain {len(gained)} = {len(gained) / max(len(low), 1):.3f} (bar >= 0.5); by tier own {sum('own' in v for v in gained.values())} tied {sum('tied' in v for v in gained.values())} none-only {sum(v == {'none'} for v in gained.values())}", file=w)
    print(open(f'{out}/summary.txt').read())

if __name__ == '__main__': main()

#!/usr/bin/env python3
"""PREREG r1537: fate of a T2T copy's reads on a non-T2T genome (docs/PREREG_missing_copy_hallmarks_2026-11-04.md).
  hallmarks.py origin T2T.bam PANEL.genes.tsv TRUTH2.tsv > read_origin.tsv      (primary T2T record over the copy's exons)
  hallmarks.py score read_origin.tsv NT.bam TRUTH2.tsv OUT [thresholds.tsv]      -> OUT.reads.tsv OUT.copies.tsv OUT.thresholds.tsv
Thresholds: 95th percentile of the expressed present copies (development) or the given file (held-out species)."""
import collections
import statistics
import subprocess
import sys

import pysam

FEATS = ['unmapped', 'mapq0', 'tie', 'de_gap', 'clip']


def truth(path):
    return [l.rstrip('\n').split('\t') for l in open(path) if not l.startswith('gene\t')]


def origin(t2t, genes, truth2):
    exp = {r[0]: r[3] for r in truth(truth2) if r[4] == '1'}
    rows = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open(genes) if not l.startswith('gene\t')}
    hits = collections.defaultdict(list)
    for g, st in exp.items():
        d = rows[g]
        reg = ' '.join(f'{d[3]}:{int(s) + 1}-{e}' for s, e in (x.split('-') for x in d[5].split(',')))
        txt = subprocess.run(f'samtools view -F 0x904 {t2t} {reg}', shell=True, capture_output=True, text=True, check=True).stdout
        seen = set()
        for l in txt.splitlines():
            f = l.split('\t')
            if f[0] in seen:
                continue
            seen.add(f[0])
            de = [x[5:] for x in f[11:] if x.startswith('de:f:')]
            hits[f[0]].append((g, d[1], d[2], st, de[0] if de else 'NA', f[4]))
    print('read\tgene\tname\tfamily\tstatus\tde_t2t\tmapq_t2t')
    amb = 0
    for r, hs in hits.items():
        if len(hs) > 1:
            amb += 1
            continue
        print(r, *hs[0], sep='\t')
    print(f'reads {len(hits)} ambiguous (dropped) {amb}', file=sys.stderr)


def score(origin_tsv, nt, truth2, out, thr_file=None):
    org = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1:] for l in open(origin_tsv) if not l.startswith('read')}
    loci = []  # present copy loci on the non-T2T genome
    for r in truth(truth2):
        if r[3] == 'present' and r[7] != '-':
            c, rest = r[7].rsplit(':', 1)
            s, e = rest[:-1].split('-')
            loci.append((c, int(s), int(e), r[0], r[1]))
    own = collections.Counter(g for g, *_ in org.values())  # T2T reads per copy (depth-excess denominator)
    prim, best2 = {}, collections.Counter()
    for a in pysam.AlignmentFile(nt):
        if a.is_supplementary:
            continue
        if a.is_secondary:
            if a.query_name in org:
                best2[a.query_name] = max(best2[a.query_name], a.get_tag('AS'))
            continue
        if a.query_name in org:
            prim[a.query_name] = a
    rd = {}
    with open(out + '.reads.tsv', 'w') as w:
        w.write('read\tgene\tname\tstatus\tunmapped\tmapq0\ttie\tabsorbed_on\tde_gap\tclip\n')
        for r, (g, name, fam, st, de_t2t, _) in org.items():
            a = prim.get(r)
            if a is None or a.is_unmapped:
                f = dict(unmapped=1, mapq0=1, tie=0, absorbed='-', de_gap=float('nan'), clip=1.0)
            else:
                on = '-'
                for c, s, e, gg, nn in loci:
                    if c == a.reference_name and a.reference_start < e and a.reference_end > s:
                        on = gg
                        break
                cl = max([x[1] for x in a.cigartuples if x[0] in (4, 5)] or [0])
                dg = a.get_tag('de') - float(de_t2t) if de_t2t != 'NA' else float('nan')
                f = dict(unmapped=0, mapq0=int(a.mapping_quality == 0), tie=int(best2[r] >= a.get_tag('AS')), absorbed=on,
                         de_gap=dg, clip=cl / a.infer_read_length())
            rd.setdefault(g, []).append(f)
            w.write(f"{r}\t{g}\t{name}\t{st}\t{f['unmapped']}\t{f['mapq0']}\t{f['tie']}\t{f['absorbed']}\t{f['de_gap']:.4f}\t{f['clip']:.3f}\n")
    names = {o[0]: (o[1], o[3]) for o in org.values()}
    cop = {}
    for g, fs in rd.items():
        n = len(fs)
        ab = collections.Counter(f['absorbed'] for f in fs if f['absorbed'] != '-')
        top, k = ab.most_common(1)[0] if ab else ('-', 0)
        cop[g] = dict(n=n, unmapped=sum(f['unmapped'] for f in fs) / n, mapq0=sum(f['mapq0'] for f in fs) / n,
                      tie=sum(f['tie'] for f in fs) / n, absorbed=sum(ab.values()) / n,
                      de_gap=statistics.median([f['de_gap'] for f in fs if f['de_gap'] == f['de_gap']] or [float('nan')]),
                      clip=statistics.median(f['clip'] for f in fs), locus=top, excess=k / own[top] if top in own else float('nan'))
    if thr_file:
        thr = {l.split('\t')[0]: float(l.split('\t')[1]) for l in open(thr_file)}
    else:
        pres = [c for g, c in cop.items() if names[g][1] == 'present']
        thr = {k: sorted(c[k] for c in pres)[int(0.95 * (len(pres) - 1))] for k in FEATS}
        with open(out + '.thresholds.tsv', 'w') as w:
            w.writelines(f'{k}\t{v}\n' for k, v in thr.items())
    with open(out + '.copies.tsv', 'w') as w:
        w.write('gene\tname\tstatus\tn\t' + '\t'.join(FEATS) + '\tabsorbed\tlocus\tdepth_excess\thallmarks\n')
        for g, c in sorted(cop.items(), key=lambda x: (names[x[0]][1], -x[1]['n'])):
            hm = [k for k in FEATS if c[k] > thr[k]]
            w.write(f"{g}\t{names[g][0]}\t{names[g][1]}\t{c['n']}\t" + '\t'.join(f'{c[k]:.4f}' for k in FEATS)
                    + f"\t{c['absorbed']:.3f}\t{c['locus']}\t{c['excess']:.2f}\t{','.join(hm) or '-'}\n")
            c['hm'] = bool(hm)
    tr = {r[0]: r[8] if len(r) > 8 else '' for r in truth(truth2)}
    for st in ('present', 'missing'):
        sc = [(g, c) for g, c in cop.items() if names[g][1] == st and c['n'] >= 10]
        h = sum(c['hm'] for _, c in sc)
        print(f'{st}: scorable {len(sc)} hallmarked {h} ({h / max(1, len(sc)):.2f})')
        if st == 'missing':
            for reason in ('collapsed', 'absent'):
                s2 = [c for g, c in sc if tr.get(g) == reason]
                print(f'  {reason}: {sum(c["hm"] for c in s2)} of {len(s2)}')
    print('thresholds', ' '.join(f'{k}={v:.4f}' for k, v in thr.items()))


if sys.argv[1] == 'origin':
    origin(*sys.argv[2:5])
else:
    score(*sys.argv[2:7])

#!/usr/bin/env python3
"""PREREG r1499 report: SQANTI3 categories / artifact flags / rules-filter removals of a final.gtf, split base (model also in
PREFIX.gtf with the same id and exons) vs layer, and single-copy vs multicopy "true".
  sqf_report.py CLASSIFICATION.txt RULESFILTER_CLASSIFICATION.txt PREFIX.gtf PREFIX.final.gtf"""
import collections
import csv
import re
import sys


def models(path):
    ex, multi = collections.defaultdict(list), {}
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        if f[2] == 'exon':
            ex[t].append((f[0], f[6], int(f[3]), int(f[4])))
        elif f[2] == 'transcript':
            multi[t] = 'multicopy "true"' in f[8]
    return {t: tuple(sorted(v)) for t, v in ex.items()}, multi


def main(cls, filt, base_gtf, final_gtf):
    base, _ = models(base_gtf)
    final, multi = models(final_gtf)
    res = {r['isoform']: r['filter_result'] for r in csv.DictReader(open(filt), delimiter='\t')}
    tab = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(cls), delimiter='\t'):
        t = r['isoform']
        src = 'base' if t in base and base[t] == final.get(t) else 'layer'
        cp = 'multi' if multi.get(t) else 'single'
        art = res.get(t) == 'Artifact'
        for k in (f'{src}', f'{cp}', 'all'):
            c = tab[k]
            c['n'] += 1
            c['cat:' + r['structural_category']] += 1
            c['removed'] += art
            c['RTS'] += r['RTS_stage'] == 'TRUE' or r['RTS_stage'] == 'True'
            c['noncanonical'] += r['all_canonical'] == 'non_canonical'
            try:
                c['intrapriming'] += float(r['perc_A_downstream_TTS']) > 60
            except ValueError:
                pass
    for k in ('all', 'base', 'layer', 'single', 'multi'):
        c = tab[k]
        n = max(c['n'], 1)
        cats = ' '.join(f"{x[4:]}={c[x] / n:.3f}" for x in sorted(c) if x.startswith('cat:'))
        print(f"{k}\tn {c['n']}\tremoved {c['removed']} ({c['removed'] / n:.3f})\tRTS {c['RTS'] / n:.3f}\tnoncanonical "
              f"{c['noncanonical'] / n:.3f}\tintrapriming {c['intrapriming'] / n:.3f}\t{cats}")
    s, m = tab['single'], tab['multi']
    if s['n'] and m['n'] and s['removed']:
        print(f"multi/single removal ratio {(m['removed'] / m['n']) / (s['removed'] / s['n']):.2f} (bar <= 1.5)")


if __name__ == '__main__':
    main(*sys.argv[1:5])

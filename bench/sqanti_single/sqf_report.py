#!/usr/bin/env python3
"""PREREG r1499 report: SQANTI3 categories / artifact flags / rules-filter removals of a final.gtf, split base (model also in
PREFIX.gtf with the same exons, any id) vs layer, and single-copy vs multicopy (exons overlap an SD BED, amendment 1).
  sqf_report.py CLASSIFICATION.txt RULESFILTER_CLASSIFICATION.txt PREFIX.gtf PREFIX.final.gtf SD.bed[.gz]"""
import collections
import bisect
import csv
import gzip
import re
import sys


def models(path):
    ex = collections.defaultdict(list)
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        ex[t].append((f[0], f[6], int(f[3]) - 1, int(f[4])))
    return {t: tuple(sorted(v)) for t, v in ex.items()}


def sd_index(path):
    by = collections.defaultdict(list)
    for line in (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)):
        f = line.split('\t')
        if line[0] != '#' and len(f) >= 3:
            by[f[0]].append((int(f[1]), int(f[2])))
    out = {}
    for c, v in by.items():
        m = []
        for a, b in sorted(v):
            if m and a <= m[-1][1]:
                m[-1][1] = max(m[-1][1], b)
            else:
                m.append([a, b])
        out[c] = ([x[0] for x in m], [x[1] for x in m])
    return out


def in_sd(sd, exons):
    for c, _, a, b in exons:
        if c in sd:
            st, en = sd[c]
            i = bisect.bisect_right(st, b) - 1
            if i >= 0 and en[i] > a:
                return True
    return False


def main(cls, filt, base_gtf, final_gtf, sd_bed):
    base = set(models(base_gtf).values())
    final = models(final_gtf)
    sd = sd_index(sd_bed)
    res = {r['isoform']: r['filter_result'] for r in csv.DictReader(open(filt), delimiter='\t')}
    tab = collections.defaultdict(collections.Counter)
    for r in csv.DictReader(open(cls), delimiter='\t'):
        t = r['isoform']
        src = 'base' if final.get(t) in base else 'layer'
        cp = 'multi' if t in final and in_sd(sd, final[t]) else 'single'
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
    main(*sys.argv[1:6])

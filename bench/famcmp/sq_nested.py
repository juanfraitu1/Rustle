#!/usr/bin/env python3
"""PREREG row 1716 scorer: CAGE support of nested-chain models (SQANTI3 classification + own shifted background).

  sq_nested.py MODELS_GTF CAGE_BED6 SQANTI_classification.txt   |   sq_nested.py --selftest
Model ids: CONTIG_CLASS_GROUP_PEAKED_READS_N (nested_sub.py --models)."""
import sys, csv, bisect, collections

def tss_and_shift(strand, exons, d=100):
    """exons: sorted 0-based half-open. -> (TSS position 0-based, position d bp downstream in the transcript or None)"""
    ex = exons if strand == '+' else [(-e, -s) for s, e in reversed(exons)]   # transcript-order, increasing
    tss = ex[0][0]; left = d; pos = None
    for s, e in ex:
        if e - s > left: pos = s + left; break
        left -= e - s
    f = (lambda p: p) if strand == '+' else (lambda p: -p - 1)
    return f(tss), (None if pos is None else f(pos))

class Peaks:
    def __init__(self, path):
        d = collections.defaultdict(list)
        for l in open(path):
            f = l.split('\t'); d[(f[0], f[5].strip())].append((int(f[1]), int(f[2])))
        self.s = {k: sorted(v) for k, v in d.items()}; self.st = {k: [a for a, _ in v] for k, v in self.s.items()}
    def dist(self, contig, strand, p):
        """0 when p is inside a peak on that strand, else bp to the nearest one (inf when none)."""
        v = self.s.get((contig, strand)); best = float('inf')
        if not v: return best
        i = bisect.bisect_right(self.st[(contig, strand)], p)
        for j in (i - 2, i - 1, i, i + 1):   # shortcut: peaks are narrow (p95 56 bp), neighbours suffice
            if 0 <= j < len(v): a, b = v[j]; best = min(best, 0 if a <= p < b else (a - p if p < a else p - b + 1))
        return best

def selftest():
    assert tss_and_shift('+', [(100, 150), (300, 500)]) == (100, 350) and tss_and_shift('-', [(100, 150), (300, 500)]) == (499, 399)
    assert tss_and_shift('-', [(100, 150), (300, 360)]) == (359, 109) and tss_and_shift('+', [(0, 50)]) == (0, None)
    print('selftest ok')

if __name__ == '__main__':
    if sys.argv[1] == '--selftest': selftest(); sys.exit()
    gtf, bed, cls_p = sys.argv[1:4]; P = Peaks(bed); ex = collections.defaultdict(list); st = {}
    for l in open(gtf):
        f = l.split('\t')
        if f[2] == 'exon': t = f[8].split('transcript_id "')[1].split('"')[0]; ex[t].append((int(f[3]) - 1, int(f[4]))); st[t] = (f[0], f[6])
    sq = {r['isoform']: r for r in csv.DictReader(open(cls_p), delimiter='\t')}
    T = collections.defaultdict(collections.Counter); agree = n_cmp = 0
    for t, e in ex.items():
        contig, strand = st[t]; cls, grp, peaked = t[len(contig) + 1:].split('_')[:3]; tss, sh = tss_and_shift(strand, sorted(e))
        d = P.dist(contig, strand, tss); r = sq.get(t)
        for key in ((cls, grp, 'all'), (cls, grp, 'peaked' if peaked == '1' else 'dispersed')):
            c = T[key]; c['n'] += 1; c['own_in'] += d == 0; c['own_50'] += d <= 50
            if r: c['sq'] += 1; c['sq_in'] += r['within_CAGE_peak'].lower() == 'true'; c['pa'] += r['polyA_motif'] not in ('NA', ''); c['ip'] += r['perc_A_downstream_TTS'] not in ('NA', '') and float(r['perc_A_downstream_TTS']) >= 60
            if sh is not None: ds = P.dist(contig, strand, sh); c['bg_n'] += 1; c['bg_in'] += ds == 0; c['bg_50'] += ds <= 50
        if r: n_cmp += 1; agree += (d == 0) == (r['within_CAGE_peak'].lower() == 'true')
    print(f'own flag vs SQANTI3 within_CAGE_peak: agree {agree} of {n_cmp} ({agree / max(1, n_cmp):.4f})')
    print('class\tgroup\tstratum\tmodels\tin SQANTI3\twithin peak (SQANTI3)\twithin peak (own)\twithin 50 bp\tbackground in peak\tbackground 50 bp\tpolyA motif\tA-rich >= 60%')
    for k in sorted(T):
        c = T[k]; q = lambda a, b: f"{c[a] / max(1, c[b]):.3f}"
        print('\t'.join(k) + f"\t{c['n']}\t{c['sq']}\t{q('sq_in', 'sq')}\t{q('own_in', 'n')}\t{q('own_50', 'n')}\t{q('bg_in', 'bg_n')}\t{q('bg_50', 'bg_n')}\t{q('pa', 'sq')}\t{q('ip', 'sq')}")

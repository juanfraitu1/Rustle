#!/usr/bin/env python3
"""Ends away from the annotated ones: real alternative ends or artifacts? (register 1992, addendum)

    unexpected_ends.py TRUTH.gtf CAGE.bed POLYA.bed NAME=FILE.gtf [...]

Per tool and end: transcripts with an annotated chain by where the end lies against the nearest annotated end of the chain
(expected = within 50 bp, short = > 50 bp inside, past = > 50 bp beyond), and for each class the share within PAD bp of a
CAGE peak (5' end) / a polyA site (3' end). Then the transcripts whose chain is not annotated."""
import sys, os, collections, importlib.util
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread
def _load(n):
    s = importlib.util.spec_from_file_location(n, os.path.join(os.path.dirname(os.path.abspath(__file__)), n + '.py')); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
EA, EC = _load('end_accuracy'), _load('end_choice')
PAD = 25

def main(truth, cage, polya, *arms):
    A = collections.defaultdict(lambda: ([], []))
    for c, s, ex, _ in oneread.models(truth).values():
        if len(ex) > 1 and s in '+-': a = A[(c, s, oneread.chain(ex))]; e = EA.ends(c, s, ex); a[0].append(e[0]); a[1].append(e[1])
    ix = {'5': EC.bed_index(cage, PAD), '3': EC.bed_index(polya, PAD)}
    print("tool\tend\tannotated chains\texpected: n, share, signal\tshort: n, share, signal\tpast: n, share, signal\taway (short + past): signal\tunannotated chains: n, signal")
    for arm in arms:
        name, path = arm.split('=', 1); T = {e: collections.defaultdict(lambda: [0, 0]) for e in '53'}
        for c, s, ex, _ in oneread.models(path).values():
            if len(ex) < 2 or s not in '+-': continue
            k = (c, s, oneread.chain(ex)); e = EA.ends(c, s, ex)
            for i, end in enumerate('53'):
                h = EC.hit(ix[end], c, e[i])
                if k in A:
                    d = EA.signed(e[i], A[k][i], s, end == '5'); cl = 'expected' if abs(d) <= 50 else 'short' if d < 0 else 'past'
                else: cl = 'unannotated'
                t = T[end][cl]; t[0] += 1; t[1] += h
        for end in '53':
            t = T[end]; n = sum(t[k][0] for k in ('expected', 'short', 'past')) or 1; f = lambda k: f"{t[k][0]}, {t[k][0] / n:.3f}, {t[k][1] / max(1, t[k][0]):.3f}"
            aw = (t['short'][1] + t['past'][1]) / max(1, t['short'][0] + t['past'][0])
            print(name, end + "'", n, f('expected'), f('short'), f('past'), f'{aw:.3f}', f"{t['unannotated'][0]}, {t['unannotated'][1] / max(1, t['unannotated'][0]):.3f}", sep='\t')

if __name__ == '__main__': main(*sys.argv[1:])

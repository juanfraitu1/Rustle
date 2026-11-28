#!/usr/bin/env python3
"""Transcript starts and ends against an annotation, per tool (register 1992).

    end_accuracy.py TRUTH.gtf NAME=FILE.gtf [NAME=FILE.gtf ...]

gffcompare matches multi-exon transcripts by intron chain and does not look at their two ends. Here, for every multi-exon
transcript whose chain is annotated: the distance of its 5' and of its 3' end to the NEAREST annotated end of that chain
(signed: negative = the transcript stops short of it, inside; positive = it runs past it). A tool that reports several
transcripts of one chain is scored on each; `chains with > 1 model` counts them.
First lines: what the annotation itself holds (chains with several annotated starts / ends more than 50 bp apart)."""
import sys, os, collections, statistics
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../tools'))
import oneread

def ends(c, s, ex): return (ex[0][0], ex[-1][1]) if s == '+' else (ex[-1][1], ex[0][0])     # (5', 3')

def signed(x, ref, s, five):
    """distance to the nearest of ref; > 0 = past it (longer transcript), < 0 = short of it"""
    d = min(ref, key=lambda r: abs(x - r)); out = (d - x) if (s == '+') == five else (x - d)
    return out

def spread(v): return max(v) - min(v) > 50

def main(truth, *arms):
    A = collections.defaultdict(lambda: ([], []))
    for c, s, ex, _ in oneread.models(truth).values():
        if len(ex) > 1 and s in '+-': a = A[(c, s, oneread.chain(ex))]; e = ends(c, s, ex); a[0].append(e[0]); a[1].append(e[1])
    n = len(A); print(f'annotation: {n} spliced chains; with starts > 50 bp apart {sum(spread(a[0]) for a in A.values()) / n:.3f}; with ends > 50 bp apart {sum(spread(a[1]) for a in A.values()) / n:.3f}')
    print("tool\tmatched transcripts\tchains with > 1 model\t5' within 10 / 50 / 100 bp\t5' short by > 50 / past by > 50\t5' median\t3' within 10 / 50 / 100 bp\t3' short by > 50 / past by > 50\t3' median")
    for arm in arms:
        name, path = arm.split('=', 1); d5 = []; d3 = []; per = collections.Counter()
        for c, s, ex, _ in oneread.models(path).values():
            k = (c, s, oneread.chain(ex)) if len(ex) > 1 else None
            if k in A:
                e = ends(c, s, ex); d5.append(signed(e[0], A[k][0], s, True)); d3.append(signed(e[1], A[k][1], s, False)); per[k] += 1
        m = len(d5)
        if not m: print(name, 0, sep='\t'); continue
        f = lambda d: (f"{sum(abs(x) <= 10 for x in d) / m:.3f} / {sum(abs(x) <= 50 for x in d) / m:.3f} / {sum(abs(x) <= 100 for x in d) / m:.3f}", f"{sum(x < -50 for x in d) / m:.3f} / {sum(x > 50 for x in d) / m:.3f}", f'{statistics.median(d):+.0f}')
        print(name, m, sum(v > 1 for v in per.values()), *f(d5), *f(d3), sep='\t')

def selftest():
    assert signed(100, [90, 300], '+', True) == -10 and signed(80, [90], '+', True) == 10 and signed(100, [90], '+', False) == 10 and signed(100, [90], '-', True) == 10
    print('end_accuracy selftest ok')

if __name__ == '__main__': selftest() if sys.argv[1:] == ['--selftest'] else main(*sys.argv[1:])

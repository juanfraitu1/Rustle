#!/usr/bin/env python3
"""sibling_chains.py emit: the default GTF plus one model per r1628 chain (one primary record, from a no-rival molecule)."""
import os, subprocess, sys, tempfile
d = tempfile.mkdtemp(); here = os.path.dirname(os.path.abspath(__file__))
def tx(t, exons, attr=''):
    a = f'gene_id "g"; transcript_id "{t}"; {attr}'
    return [f'c\tx\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t+\t.\t{a}'] + [f'c\tx\texon\t{s}\t{e}\t.\t+\t.\t{a}' for s, e in exons]
n = lambda p, t, c: f'n_primary "{p}"; n_tied "{t}"; n_clean "{c}";'
open(f'{d}/default.gtf', 'w').write('\n'.join(tx('d1', [(1, 10), (50, 60)])) + '\n')
open(f'{d}/projected.gtf', 'w').write('\n'.join(
    tx('keep', [(1, 10), (20, 30)], n(1, 0, 1)) + tx('dup', [(3, 10), (20, 35)], n(1, 0, 1)) +  # same chain twice: one model
    tx('two', [(1, 10), (70, 80)], n(2, 0, 2)) + tx('rival', [(1, 10), (90, 99)], n(1, 0, 0)) +
    tx('already', [(1, 10), (50, 60)], n(1, 0, 1))) + '\n')
subprocess.run([sys.executable, f'{here}/sibling_chains.py', 'emit', f'{d}/default.gtf', d, f'{d}/out.gtf'], check=True)
out = open(f'{d}/out.gtf').read()
ids = [l.split('transcript_id "')[1].split('"')[0] for l in out.splitlines() if l.split('\t')[2] == 'transcript']
assert ids == ['d1', 'keep'], ids
assert '\tsibling\texon\t20\t30\t' in out and 'sibling_rule "r1628"' in out, out
print('ok')

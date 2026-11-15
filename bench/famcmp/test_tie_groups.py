#!/usr/bin/env python3
"""Self-check of tie_groups.py build: two copies tied by three molecules are one group of two loci, and a long
read whose intron spans both copies does not swallow them (r1612). Run: python3 bench/famcmp/test_tie_groups.py"""
import os, subprocess, sys, tempfile
T = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tie_groups.py')
rec = lambda n, flag, c, pos, cigar: f'{n}\t{flag}\t{c}\t{pos}\t0\t{cigar}\tAS:i:100\tms:i:100\n'
sam = ''.join(rec(m, 0, 'chr1', 1001, '100M100N100M') + rec(m, 256, 'chr1', 50001, '100M100N100M') for m in ('m1', 'm2', 'm3'))
sam += rec('m4', 0, 'chr1', 990, '110M48900N210M') + rec('m4', 256, 'chr2', 1, '320M')
with tempfile.TemporaryDirectory() as d:
    subprocess.run([sys.executable, T, 'build', d + '/t', 'as98'], input=sam, text=True, check=True, stderr=subprocess.DEVNULL)
    rows = [l.rstrip('\n').split('\t') for l in open(d + '/t.groups.tsv')][1:]
got = sorted((int(r[1]), int(r[2]), r[3]) for r in rows)
want = [(1, 1, 'chr1:989-50209+'), (1, 1, 'chr2:0-320+'), (2, 3, 'chr1:1000-1300+,chr1:50000-50300+')]
assert got == want, got
# r1620: two shapes alternating at one place are two loci however many records there are, and the r1612 case holds under cc
alt = ''.join(rec(f'a{i}', 0, 'chr3', 1001, '100M800N100M') + rec(f'a{i}', 256, 'chr4', 1001, '200M') +
              rec(f'b{i}', 0, 'chr3', 1201, '100M100N100M') + rec(f'b{i}', 256, 'chr5', 1001, '200M') for i in range(10))
with tempfile.TemporaryDirectory() as d:
    subprocess.run([sys.executable, T, 'build', d + '/t', 'as98', '0', '1', '', 'stranded', 'cc'], input=sam + alt, text=True, check=True, stderr=subprocess.DEVNULL)
    rows = [l.rstrip('\n').split('\t') for l in open(d + '/t.groups.tsv')][1:]
names = [n for r in rows for n in r[3].split(',')]
assert sum(n.startswith('chr3:') for n in names) == 2, names
# under cc the long read joins the copy its second block covers by half (4 molecules); both models stay copy-sized
got = sorted((int(r[1]), int(r[2]), ','.join(sorted(r[3].split(',')))) for r in rows if not any(x in r[3] for x in ('chr3', 'chr4', 'chr5')))
assert got == [(1, 1, 'chr2:0-320+'), (2, 4, 'chr1:1000-1300+,chr1:50000-50300+')], rows
print('ok')

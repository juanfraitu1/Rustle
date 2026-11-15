#!/usr/bin/env python3
"""RefSeq exact intron-chain precision of GTFs, and new / lost chains of each arm against the first (register 1745).

  chain_prec.py REF_GTF BASE=GTF ARM=GTF [...]   |   chain_prec.py --selftest"""
import sys
from alt_events import load

def chains(path): return {(c, s, ch) for (c, s), d in load(path).items() for ch in d}

def row(name, cs, ref, base=None):
    ok = len(cs & ref); out = f'{name}\t{len(cs)}\t{ok}\t{ok / max(1, len(cs)):.4f}'
    if base is not None:
        new, lost = cs - base, base - cs
        out += f'\t{len(new)}\t{len(new & ref)}\t{len(new & ref) / max(1, len(new)):.4f}\t{len(lost)}\t{len(lost & ref)}'
    return out

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        ref = {1, 2, 3}; assert row('a', {1, 2, 9}, ref) == 'a\t3\t2\t0.6667' and row('b', {1, 3, 8, 7}, ref, {1, 2, 9}).endswith('\t3\t1\t0.3333\t2\t1')
        print('selftest ok'); sys.exit()
    ref = chains(sys.argv[1]); arms = [a.split('=', 1) for a in sys.argv[2:]]; base = None
    print('arm\tmulti-exon chains\tcorrect\tprecision\tnew\tnew correct\tnew precision\tlost\tlost correct')
    for n, p in arms:
        cs = chains(p); print(row(n, cs, ref, base)); base = cs if base is None else base

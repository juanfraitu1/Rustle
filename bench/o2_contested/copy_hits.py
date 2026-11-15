#!/usr/bin/env python3
"""PREREG row 1703: haplotype loci of catalog copies, as haplotype_truth.py writes its `C|family|copy` rows.

  copy_hits.py PREFIX.fam.copies.fa FAMILIES.txt MAT.mmi PAT.mmi OUT_HITS.tsv"""
import sys, os, subprocess, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'o2_truth'))
from haplotype_truth import parse_sam
fa, fams, mat, pat, out = sys.argv[1:6]; fams = set(open(fams).read().split())
with tempfile.TemporaryDirectory() as d:
    q = os.path.join(d, 'q.fa'); on = False
    with open(q, 'w') as w:
        for l in open(fa):
            if l.startswith('>'):
                h = l[1:].split('|'); on = h[0] in fams
                if on: w.write(f'>C|{h[0]}|{h[1]}\n')
            elif on: w.write(l)
    with open(out, 'w') as w:
        w.write('query\thap\tchrom\tstart\tend\tnm\taligned\n')
        for hap, idx in (('mat', mat), ('pat', pat)):
            sam = subprocess.run(['minimap2', '-ax', 'splice:hq', '-uf', '--eqx', '-N', '10', '-t', '16', idx, q], capture_output=True, text=True, check=True).stdout
            for name, hs in parse_sam(sam.splitlines()).items():
                for c, s, e, nm, al in hs: w.write(f'{name}\t{hap}\t{c}\t{s}\t{e}\t{nm}\t{al}\n')

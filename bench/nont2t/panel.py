#!/usr/bin/env python3
"""PREREG r1504 panel: the T2T genes of the MAINSTAYS + YAGS families (gorilla also GGO_MAINSTAYS) -> PREFIX.genes.tsv
(gene, name, family, chrom, strand, exons), PREFIX.exons.bed (per exon, name = gene) and PREFIX.copies.fa (spliced exon union,
transcript orientation).  panel.py human|gorilla T2T.fa PREFIX   (env as famtruth: DENOVO_CAT_*, DENOVO_GGO_GFF)"""
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import famtruth as F  # noqa: E402

species, fa, out = sys.argv[1:4]
G = F.load_cat() if species == 'human' else F.load_ggo(F.load_cat())
pats = {**F.MAINSTAYS, **F.YAGS}
if species == 'gorilla':
    pats.update(F.GGO_MAINSTAYS)
rc = str.maketrans('ACGTacgtN', 'TGCAtgcaN')
rows = []
for g, d in sorted(G.items()):
    for fam, p in pats.items():
        if re.search(p, d['name']) or (species == 'gorilla' and fam in F.GGO_MAINSTAYS and re.search(p, d.get('desc', ''))):
            rows.append((g, d, fam))
            break
with open(out + '.genes.tsv', 'w') as w, open(out + '.exons.bed', 'w') as b, open(out + '.regions', 'w') as r:
    w.write('gene\tname\tfamily\tchrom\tstrand\texons\n')
    for g, d, fam in rows:
        w.write(f"{g}\t{d['name']}\t{fam}\t{d['chrom']}\t{d['strand']}\t{','.join(f'{s}-{e}' for s, e in d['exons'])}\n")
        for s, e in d['exons']:
            b.write(f"{d['chrom']}\t{s}\t{e}\t{g}\n")
            r.write(f"{d['chrom']}:{s + 1}-{e}\n")
seqs = subprocess.run(['samtools', 'faidx', fa, '-r', out + '.regions'], capture_output=True, text=True, check=True).stdout
chunks = [''.join(c.split('\n')[1:]) for c in seqs.split('>')[1:]]
i = 0
with open(out + '.copies.fa', 'w') as w:
    for g, d, fam in rows:
        s = ''.join(chunks[i:i + len(d['exons'])])
        i += len(d['exons'])
        if d['strand'] == '-':
            s = s.translate(rc)[::-1]
        w.write(f'>{g} {d["name"]} {fam}\n{s}\n')
print(f'{species}: {len(rows)} panel genes', file=sys.stderr)

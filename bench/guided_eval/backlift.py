"""GUIDED-S10 reciprocity: which hit genes lift back onto their seed gene (coverage >= 0.5 of the hit gene's CDS, or body when no CDS; identity >= 0.5).

  backlift.py HITS_TSV META_JSON OUT_TSV PAF...   HITS = seed_gene<TAB>hit_gene; META = hit gene -> {cds_len, body_cds}; body PAFs named A*.paf
"""
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cdsseed as C  # noqa: E402
import famtruth as F  # noqa: E402
import guided3 as g3  # noqa: E402

hits, meta, out, pafs = sys.argv[1], json.load(open(sys.argv[2])), sys.argv[3], sys.argv[4:]
SEEDLEN = os.environ.get('BACKLIFT_SEEDLEN')  # GUIDED-S11: JSON {seed gene: CDS length}; coverage over the shorter of hit and seed
G = g3.load_genes('human')
idx = F.exon_index(G)
want = collections.defaultdict(set)
for l in open(hits):
    s, h = l.rstrip('\n').split('\t')[:2]
    want[h].add(s)
body_cds = {h: [tuple(x) for x in m['body_cds']] for h, m in meta.items()}
length = {h: m['cds_len'] for h, m in meta.items()}
for h in want:
    if h not in length:  # no CDS: the whole body is the template
        length[h] = sum(e - s for s, e in G[h]['exons'])  # exon sum, as the body record's aligned bases are exonic+intronic but the gene's own unit is its exons
groups = collections.defaultdict(list)
for q, t, st, bl, qi, ti, cdsb, ident in g3.records(pafs, body_cds):
    if q not in want or not bl or cdsb <= 0:
        continue
    g = C.gene_of(idx, G, t, st, bl)
    if g in want[q]:
        groups[(q, g, st)].append((qi, ti, cdsb, ident))
seedlen = json.load(open(SEEDLEN)) if SEEDLEN else {}
ok = set()
for (h, s, st), rs in groups.items():
    bases, ident = g3.chain(rs, st)
    denom = min(length[h], seedlen.get(s, length[h])) if SEEDLEN else length[h]
    if bases / denom >= g3.COV and ident >= g3.IDENT:
        ok.add((s, h))
with open(out, 'w') as w:
    for s, h in sorted(ok):
        w.write(f'{s}\t{h}\n')
n = sum(len(v) for v in want.values())
print(f'hit pairs {n}, reciprocal {len(ok)} ({len(ok) / max(n, 1):.3f})')

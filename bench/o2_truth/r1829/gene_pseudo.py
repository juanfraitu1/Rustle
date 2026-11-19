#!/usr/bin/env python3
"""r1829 side table (descriptive; user 2026-11-19): the score gap between a protein-coding gene and a pseudogene when one
read aligns to both. Each multimapper's primary and rival place (as_table; rival = best-AS non-primary record at another
place) get the biotype of the CAT gene whose exons contain the record's leftmost position (protein-coding preferred when
several). For pairs coding <-> processed pseudogene, coding <-> unprocessed pseudogene and coding <-> coding, the gap
coding - other on ms and on AS, binned, and how often the non-coding place is the primary.
  gene_pseudo.py SAMPLE AS_TABLE CAT_GENES.tsv > SAMPLE.gene_pseudo.tsv"""
import bisect, collections, sys

GROUP = {'protein_coding': 'coding',
         'processed_pseudogene': 'processed pseudogene', 'transcribed_processed_pseudogene': 'processed pseudogene',
         'translated_processed_pseudogene': 'processed pseudogene',
         'unprocessed_pseudogene': 'unprocessed pseudogene', 'transcribed_unprocessed_pseudogene': 'unprocessed pseudogene',
         'unitary_pseudogene': 'unprocessed pseudogene', 'transcribed_unitary_pseudogene': 'unprocessed pseudogene',
         'polymorphic_pseudogene': 'unprocessed pseudogene'}
RANK = {'coding': 0, 'unprocessed pseudogene': 1, 'processed pseudogene': 2}
BINS = [(-10 ** 9, -51, '< -50'), (-50, -26, '-50..-26'), (-25, -11, '-25..-11'), (-10, -1, '-10..-1'), (0, 0, '0'),
        (1, 5, '1..5'), (6, 10, '6..10'), (11, 25, '11..25'), (26, 50, '26..50'), (51, 100, '51..100'),
        (101, 200, '101..200'), (201, 10 ** 9, '> 200')]

def binof(v):
    return next(l for a, b, l in BINS if a <= v <= b)

s, table, genes = sys.argv[1:4]
iv = collections.defaultdict(list)
for l in open(genes):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'gene_id' or f[2] not in GROUP:
        continue
    for blk in f[9].split(','):
        a, b = blk.split('-'); iv[f[4]].append((int(a), int(b), GROUP[f[2]]))
idx = {}
for c, v in iv.items():
    v.sort(); idx[c] = ([x[0] for x in v], v, max(x[1] - x[0] for x in v))

def group_at(c, p):
    if c not in idx:
        return 'none'
    st, v, mx = idx[c]
    best = None
    for k in range(bisect.bisect_left(st, p - mx), bisect.bisect_right(st, p)):
        a, b, g = v[k]
        if a <= p < b and (best is None or RANK[g] < RANK[best]):
            best = g
    return best or 'none'

C = collections.Counter()
for l in open(table):
    if l[0] == '#':
        continue
    f = l.rstrip('\n').split('\t')
    if f[4] in ('NA', '*') or f[11] == 'NA':
        continue
    pc, rc, pp, rp = f[4], f[14], int(f[15]) - 1, int(f[16]) - 1
    if rc == pc and abs(rp - pp) < 5000:
        continue
    gp, gr = group_at(pc, pp), group_at(rc, rp)
    if 'coding' not in (gp, gr) or gp == 'none' or gr == 'none':
        continue
    other = gr if gp == 'coding' else gp
    pair = 'coding - coding' if other == 'coding' else 'coding - ' + other
    if pair == 'coding - coding':
        cm, om, ca, oa, cs, os_ = int(f[8]), int(f[12]), int(f[5]), int(f[11]), f[10], f[13]   # primary as "coding"
        other_primary = 0
    elif gp == 'coding':
        cm, om, ca, oa, cs, os_ = int(f[8]), int(f[12]), int(f[5]), int(f[11]), f[10], f[13]; other_primary = 0
    else:
        cm, om, ca, oa, cs, os_ = int(f[12]), int(f[8]), int(f[11]), int(f[5]), f[13], f[10]; other_primary = 1
    C[(pair, 'reads')] += 1; C[(pair, 'other_is_primary')] += other_primary
    C[(pair, 'coding_spliced_other_not')] += (cs == '1' and os_ == '0')
    C[(pair, 'ms', binof(cm - om))] += 1; C[(pair, 'AS', binof(ca - oa))] += 1
    # implied divergence: the ms gap read as substitutions (5 each) over the coding place's ms (about its aligned length)
    d = (cm - om) / 5 / max(1, cm)
    C[(pair, 'implied_divergence', '< 0' if d < 0 else '0' if d == 0 else '< 0.5%' if d < 0.005 else '0.5-1%' if d < 0.01 else
       '1-2%' if d < 0.02 else '2-5%' if d < 0.05 else '5-10%' if d < 0.1 else '>= 10%')] += 1
    C[(pair, 'AS_prefers_other_while_ms_prefers_coding')] += (ca < oa and cm > om)
print('sample\tpair\tkey\tbin\tcount')
for k, v in sorted(C.items()):
    print(f"{s}\t{k[0]}\t{k[1]}\t{k[2] if len(k) > 2 else '-'}\t{v}")

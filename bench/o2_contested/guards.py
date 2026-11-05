#!/usr/bin/env python3
"""PREREG r1541 allele-site guard as a post-filter (docs/PREREG_o2_shape_allele_guards_2026-11-05.md).
  guards.py STAR_READS.tsv STRICT_DUMP.tsv BAM COPIES.tsv > report
A PSV column is ALLELIC for a copy when that copy's own uniquely placed reads (primary, MAPQ 60) at its genome position carry
>= 2 alleles with >= 2 reads each and the minor at >= 20%. Decisive columns for the pair (assigned, rival): read base equals one
candidate's base and not the other's; allelic columns (at either candidate's position) do not count. Verdict per decidable row:
KEEP (every rival: non-allelic columns for the assigned copy >= 2 and > columns for the rival), MOVED (some rival wins), else
ABSTAIN. Scored against the strict dump's right / wrong and its truth copies."""
import collections
import csv
import functools
import os
import sys

import pysam

star, dump, bam_path, copies = sys.argv[1:5]
bam = pysam.AlignmentFile(bam_path)
chrom = {}
for r in csv.DictReader(open(copies), delimiter='\t'):
    chrom[(r['family_id'], r['copy_idx'])] = r['chrom']
rows = {}
for l in open(dump):
    f = l.rstrip('\n').split('\t')
    if f[2] == 'assigned':
        rows[(f[0], f[1])] = f  # read, family, status, verdict, aligner, copy, truth


@functools.lru_cache(maxsize=None)
def allelic(c, pos):
    cnt = collections.Counter()
    for col in bam.pileup(c, pos, pos + 1, truncate=True, min_base_quality=10, min_mapping_quality=60, flag_filter=0x904):
        for p in col.pileups:
            if p.is_del or p.is_refskip or p.query_position is None:
                continue
            cnt[p.alignment.query_sequence[p.query_position]] += 1
    top = cnt.most_common(2)
    return len(top) == 2 and top[1][1] >= 2 and top[1][1] >= 0.2 * sum(cnt.values())


res = collections.Counter()
detail = []
for r in csv.DictReader(open(star), delimiter='\t'):
    key = (r['read_name'], r['family_id'])
    if key not in rows or r['status'] != 'assigned':
        continue
    cands = r['candidates'].split(',')
    a = cands.index(r['catalog_copy_idx']) if r['catalog_copy_idx'] in cands else None
    if a is None or r['n_cols'] == '0':
        res[(rows[key][3], 'nostar')] += 1
        continue
    cols = r['columns'].split(',')
    poss = r['positions'].split(',')
    fav = collections.Counter()  # (k) -> columns favouring candidate k over the assigned, and the assigned over k
    n_allelic = 0
    for col, pos in zip(cols, poss):
        _, rb, cb = col.split(':')
        gp = pos.split('/')
        for k in range(len(cands)):
            if k == a or cb[k] == '.' or cb[a] == '.' or cb[k] == cb[a] or rb not in (cb[k], cb[a]):
                continue
            al = False
            for j in (a, k):
                if not os.environ.get('GUARD_NOALLELIC') and gp[j] != '.' and allelic(chrom[(r['family_id'], cands[j])], int(gp[j])):
                    al = True
            if al:
                n_allelic += 1
                continue
            fav[('A', k) if rb == cb[a] else ('K', k)] += 1
    verdict = 'KEEP'
    moved = None
    for k in range(len(cands)):
        if k == a:
            continue
        fa, fk = fav[('A', k)], fav[('K', k)]
        if fk > fa:
            verdict, moved = 'MOVED', cands[k]
            break
        if fa < 2:
            verdict = 'ABSTAIN'
    truth = rows[key][6].split(';')
    new = rows[key][3]
    if verdict == 'ABSTAIN':
        new = 'abstained'
    elif verdict == 'MOVED':
        new = 'right' if f'{r["family_id"]}:{moved}' in truth else 'wrong'
    res[(rows[key][3], verdict, new, 'allelic>0' if n_allelic else 'allelic=0')] += 1
    if rows[key][3] == 'wrong' or verdict != 'KEEP':
        detail.append((r['family_id'], r['read_name'], rows[key][3], verdict, new, n_allelic, dict(fav)))
for k, v in sorted(res.items()):
    print(*k, v, sep='\t')
print('--- rows that are wrong or change')
for d in detail:
    print(*d, sep='\t')

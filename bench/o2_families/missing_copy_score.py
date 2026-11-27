#!/usr/bin/env python3
"""Score locus_sequences.py on a planted missing-copy set (register 1983).

    missing_copy_score.py PLANT_DIR SCENARIO HIDDEN.tsv DETECTOR.tsv CENSUS.tsv [FULL_REFERENCE_DETECTOR.tsv]

DETECTOR.tsv made with --names on the reduced reference; CENSUS.tsv from missing_copy_census.py (class and MAPQ 0 share
of each hidden copy). A group is POSITIVE with >= 10 primaries of hidden copies, NULL with none and >= 20 primaries
(1-9: neither). A hidden copy with >= 20 placed reads is FOUND when the group holding most of its reads is flagged.
FULL_REFERENCE_DETECTOR.tsv = the same reads on the whole genome: every flag there is false.
Prints one line of counts (for pooling) and the tables."""
import sys, csv, collections
csv.field_size_limit(1 << 30)

def score(groups, origin, hidden):
    """groups [(flag, [primary names])] -> positives (flagged, n), null (flagged, n), {hidden copy: (found, reads in its main group)}"""
    pos = [0, 0]; nul = [0, 0]; best = {}
    for f, names in groups:
        c = collections.Counter(origin[q] for q in names if origin.get(q) in hidden); h = sum(c.values())
        if h >= 10: pos[0] += f; pos[1] += 1
        elif h == 0 and len(names) >= 20: nul[0] += f; nul[1] += 1
        for g, k in c.items():
            if k > best.get(g, (0, 0))[1]: best[g] = (f, k)
    return pos, nul, best

def load(p): return [(int(r['flag']), r['primary_names'].split(',') if r['primary_names'] else []) for r in csv.DictReader(open(p), delimiter='\t')]

def main(P, sc, hid, det, census, full=None):
    hidden = {r['copy'] for r in csv.DictReader(open(hid), delimiter='\t')}
    origin = {l.split('\t')[0]: l.split('\t')[1] for l in open(f'{P}/{sc}.truth.tsv')}
    pos, nul, best = score(load(det), origin, hidden)
    C = {r['copy']: r for r in csv.DictReader(open(census), delimiter='\t')}
    T = collections.defaultdict(lambda: [0, 0])
    for g, r in C.items():
        if int(r['reads']) - int(r['unmapped']) < 20: continue
        f = best.get(g, (0, 0))[0]
        for k in (r['class'], 'tied' if float(r['mapq0']) >= 0.2 else 'placed uniquely', 'ALL', '>= 6 edits' if r['class'] in ('6-19', '>= 20') else None):
            if k: T[k][0] += f; T[k][1] += 1
    fl = [0, 0]; one = [0, 0]                      # whole reference: all groups, and those of a single locus
    if full:
        for r in csv.DictReader(open(full), delimiter='\t'):
            if int(r['primaries']) >= 20:
                fl[0] += int(r['flag']); fl[1] += 1
                if r['n_loci'] == '1': one[0] += int(r['flag']); one[1] += 1
    order = ('0 twin', '1-5', '6-19', '>= 20', 'no alignment', '>= 6 edits', 'tied', 'placed uniquely', 'ALL')
    print('COUNTS\t' + '\t'.join(f'{k}={T[k][0]}/{T[k][1]}' for k in order) + f'\tpositive_groups={pos[0]}/{pos[1]}\tnull_groups={nul[0]}/{nul[1]}\tfull_reference_groups={fl[0]}/{fl[1]}\tfull_reference_single_locus={one[0]}/{one[1]}')
    print('hidden copies found (the group holding most of their reads is flagged)')
    for k in order:
        if T[k][1]: print(f'  {k}\t{T[k][0]} of {T[k][1]}\t{T[k][0] / T[k][1]:.3f}')
    print(f'groups with >= 10 hidden-copy reads flagged\t{pos[0]} of {pos[1]}')
    print(f'false flags, reduced reference (groups with no hidden-copy read, >= 20 reads)\t{nul[0]} of {nul[1]}\t{nul[0] / max(1, nul[1]):.3f}')
    if full: print(f'false flags, whole reference (same reads, nothing hidden)\t{fl[0]} of {fl[1]}\t{fl[0] / max(1, fl[1]):.3f}\t(single-locus groups {one[0]} of {one[1]})')

if __name__ == '__main__':
    if sys.argv[1:] == ['--selftest']:
        o = {f'h{i}': 'H' for i in range(12)}; o.update({f'a{i}': 'A' for i in range(30)})
        g = [(1, [f'h{i}' for i in range(10)] + [f'a{i}' for i in range(10)]), (0, ['h10', 'h11'] + [f'a{i}' for i in range(10, 20)]), (1, [f'a{i}' for i in range(20)]), (0, ['a1'] * 5)]
        assert score(g, o, {'H'}) == ([1, 1], [1, 1], {'H': (1, 10)})
        print('missing_copy_score selftest ok')
    else: main(*sys.argv[1:])

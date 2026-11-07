#!/usr/bin/env python3
"""Row 1574: more than one representative per de novo locus (docs/PREREG_locus_reps_2026-11-07.md).

In each gene_id the representatives are chosen greedily: the transcript with the most reads (ties: longer span, then id),
then every transcript with >= MIN_READS reads that shares no exonic base with a representative already chosen, in the
same order. A locus with k >= 2 representatives becomes k groups (gene_id GENE.p1 .. GENE.pk, p1 = the primary): every
transcript joins the representative it shares the most exonic bases with (none or a tie -> the earliest). Families input
only; piece_merge.py / piece_attach.py then attach the groups to the unchanged families.

usage: locus_reps.py IN.families.gtf OUT.gtf
       locus_reps.py --selftest
"""
import collections, re, sys

MIN_READS = 2


def shared(a, b):
    return sum(max(0, min(e, e2) - max(s, s2)) for s, e in a for s2, e2 in b)


def reps(tx, min_reads=MIN_READS):
    """tx: {tid: ([(s, e)], reads)} -> representatives, primary first."""
    span = lambda t: max(e for _, e in tx[t][0]) - min(s for s, _ in tx[t][0])
    order = sorted(tx, key=lambda t: (-tx[t][1], -span(t), t))
    chosen = [order[0]]
    for t in order[1:]:
        if tx[t][1] >= min_reads and all(shared(tx[t][0], tx[c][0]) == 0 for c in chosen):
            chosen.append(t)
    return chosen


def groups(tx, min_reads=MIN_READS):
    r = reps(tx, min_reads)
    out = [[] for _ in r]
    for t in sorted(tx):
        sh = [shared(tx[t][0], tx[c][0]) for c in r]
        out[max(range(len(r)), key=lambda k: (sh[k], -k))].append(t)
    return out


def main(gtf_in, gtf_out):
    loci, lines, order = collections.defaultdict(dict), collections.defaultdict(list), []
    for line in open(gtf_in):
        f = line.split('\t')
        if len(f) < 9:
            continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]).group(1)
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        if gid not in loci:
            order.append(gid)
        lines[(gid, tid)].append(line)
        t = loci[gid].setdefault(tid, [[], 0])
        if f[2] == 'transcript':
            r = re.search(r'\breads "([^"]+)"', f[8])
            t[1] = float(r.group(1)) if r else 0
        elif f[2] == 'exon':
            t[0].append((int(f[3]) - 1, int(f[4])))
    n_multi = n_groups = 0
    with open(gtf_out, 'w') as oh:
        for gid in order:
            tx = loci[gid]
            gs = groups(tx) if len(tx) >= 2 else [sorted(tx)]
            if len(gs) == 1:
                for t in gs[0]:
                    oh.writelines(lines[(gid, t)])
                continue
            n_multi += 1
            n_groups += len(gs)
            for k, g in enumerate(gs, 1):
                for t in g:
                    oh.writelines(l.replace(f'gene_id "{gid}"', f'gene_id "{gid}.p{k}"') for l in lines[(gid, t)])
    print(f"[locus_reps] loci {len(loci)}, {n_multi} with >= 2 representatives -> {n_groups} groups", file=sys.stderr)


def selftest():
    # host gene H (exons 0-100, 900-1000) with a nested gene N in its intron (400-500) and a readthrough R joining them
    tx = {'H1': ([(0, 100), (900, 1000)], 9), 'H2': ([(0, 100), (850, 1000)], 3), 'N1': ([(400, 500)], 4),
          'R': ([(0, 100), (400, 500), (900, 1000)], 20), 'N0': ([(420, 480)], 1)}
    assert reps(tx) == ['R']                                   # the readthrough overlaps both: one representative
    del tx['R']
    assert reps(tx) == ['H1', 'N1']                            # nested gene gets its own; 1-read N0 is never one
    assert groups(tx) == [['H1', 'H2'], ['N0', 'N1']]
    print('selftest ok')


if __name__ == '__main__':
    selftest() if sys.argv[1:] == ['--selftest'] else main(*sys.argv[1:3])

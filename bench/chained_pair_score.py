#!/usr/bin/env python3
"""Arm D of `docs/PREREG_chained_jaccard_2026-09-20.md`: minimap2-style colinear DP chaining as the
pair score, replacing global Jaccard.

Arm B (§6t3) scored a gene pair by `nmatch / (la + lb - nmatch)` on its single best PAF record. But
99.3-99.9% of pairs carry several records, and a global ratio is the wrong shape for a local
phenomenon — r293 (domain-sharers outscore true paralogs) and r359 (short copy vs long parent) are both
that failure. This aggregates every anchor between the pair and rewards COLLINEAR structure.

    f[i] = w_i + max(0, max over colinear j<i of ( f[j] - gamma(gap(j,i)) ))
    colinear : q_start_i >= q_end_j and t_start_i >= t_end_j, same strand
    w_i      = nmatch_i
    gap(j,i) = |(q_start_i - q_end_j) - (t_start_i - t_end_j)|        (diagonal difference)
    gamma(g) = 0 if g == 0 else 0.01 * avg_anchor_len * g + 0.5 * log2(g)

    chain_cov = max_i f[i] / min(len_a, len_b)     -- containment, not a symmetric ratio (r359)

Emits `.clusters.tsv` in `mcl_families`' column order so the same scorer reads every arm.

Usage: chained_pair_score.py --paf P --chrom chrN --out PREFIX --t 0.30
"""
import argparse
import collections
import math


def chain(anchors, avg_len):
    """Best colinear chain score over anchors [(qs, qe, ts, te, w)], already one strand."""
    anchors.sort(key=lambda x: (x[0], x[2]))
    n = len(anchors)
    f = [0.0] * n
    best = 0.0
    for i in range(n):
        qs, qe, ts, te, w = anchors[i]
        fi = float(w)
        # minimap2 caps the predecessor search; 50 back is ample at these anchor counts
        for j in range(max(0, i - 50), i):
            pqs, pqe, pts, pte, pw = anchors[j]
            if qs < pqe or ts < pte:
                continue                      # not colinear / overlapping
            g = abs((qs - pqe) - (ts - pte))
            gamma = 0.0 if g == 0 else 0.01 * avg_len * g + 0.5 * math.log2(g)
            cand = f[j] - gamma + w
            if cand > fi:
                fi = cand
        f[i] = fi
        if fi > best:
            best = fi
    return best


def pair_scores(paf, norm='containment'):
    """(a, b) -> chain_cov, using every PAF record between the pair as an anchor."""
    by_pair = collections.defaultdict(lambda: collections.defaultdict(list))
    length = {}
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        a, b = f[0], f[5]
        if a == b:
            continue
        la, lb, nm = int(f[1]), int(f[6]), int(f[9])
        strand = f[4]
        length[a] = la; length[b] = lb
        # orient the anchor consistently: key on the sorted pair, query = first element
        if a <= b:
            anc = (int(f[2]), int(f[3]), int(f[7]), int(f[8]), nm)
            by_pair[(a, b)][strand].append(anc)
        else:
            anc = (int(f[7]), int(f[8]), int(f[2]), int(f[3]), nm)
            by_pair[(b, a)][strand].append(anc)
    out = {}
    for (a, b), by_strand in by_pair.items():
        la_, lb_ = length.get(a, 1), length.get(b, 1)
        best = 0.0
        for anchors in by_strand.values():
            avg = sum(x[4] for x in anchors) / len(anchors)
            s = chain(list(anchors), avg)
            if s > best:
                best = s
        if norm == 'jaccard':
            # POST-HOC variant (§6t4): the chained score in Jaccard form. Symmetric, so a short gene
            # contained in a long one no longer saturates at 1.0 and cannot become a hub.
            den = la_ + lb_ - best
            out[(a, b)] = (best / den) if den > 0 else 0.0
        else:
            out[(a, b)] = best / (min(la_, lb_) or 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    for x in ('--paf', '--chrom', '--out'):
        ap.add_argument(x, required=True)
    ap.add_argument('--t', type=float, required=True)
    ap.add_argument('--norm', default='containment', choices=('containment', 'jaccard'),
                    help='containment = the pre-registered arm D; jaccard = the post-hoc §6t4 variant')
    a = ap.parse_args()

    sc = pair_scores(a.paf, a.norm)
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    kept = 0
    for (x, y), v in sc.items():
        if v >= a.t:
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[rx] = ry
            kept += 1
    comp = collections.defaultdict(list)
    for n in list(parent):
        comp[find(n)].append(n)
    clusters = {k: v for k, v in comp.items() if len(v) >= 2}

    out = a.out + '.clusters.tsv'
    with open(out, 'w') as fh:
        fh.write('cluster_id\tsize\tdensity\tfrac_in\tcorroborated\tchrom\tstart\tend\n')
        for i, (_, mem) in enumerate(sorted(clusters.items(), key=lambda kv: -len(kv[1]))):
            for name in sorted(mem):
                try:
                    ch, rng = name.rsplit(':', 1); s, e = rng.split('-')
                except ValueError:
                    continue
                fh.write(f'CH{i}\t{len(mem)}\tNA\tNA\tNA\t{ch}\t{int(s)-1}\t{int(e)}\n')
    sizes = sorted((len(v) for v in clusters.values()), reverse=True)
    print(f'{a.chrom} t={a.t}: {kept} pairs kept | {len(clusters)} components | '
          f'{sum(sizes)} members | largest {sizes[0] if sizes else 0}')


if __name__ == '__main__':
    main()

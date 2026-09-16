#!/usr/bin/env python3
"""Choose the `tss` densest-5'-start window W on chr20 (development substrate) and write one simulated GTF per
W for the Rust fidelity check (docs/superpowers/specs/2026-09-16-gtf-refine-and-dedup-fix-design.md, `tss`
addendum).

usage: tss_window_sim.py MODELS_GTF REF_GTF READS_PKL OUT_DIR
       tss_window_sim.py --self-test
MODELS_GTF: fixed-dedup, no-refine `copy_assign --gtf` output (A1). READS_PKL: chr20 primaries (-F 2308, no
placement dedup) as (name, chrom, start, end, is_reverse, ts_strand, chain) tuples, 0-based half-open.
"""
import bisect
import collections
import os
import pickle
import re
import statistics
import sys

GRID = (10, 25, 50, 100)
TOL = 50  # SQANTI3's reference_match tolerance


def densest_five_prime(ends, window, strand):
    """Spec rule. '+': smallest p in ends maximizing #{x: p <= x <= p+W}. '-': largest p maximizing
    #{x: p-W <= x <= p}. None for no ends."""
    if not ends:
        return None
    v = sorted(ends)
    best, best_pos = -1, None
    if strand == '+':
        for lo in v:  # ascending: strict '>' keeps the smallest (most upstream) on ties
            c = bisect.bisect_right(v, lo + window) - bisect.bisect_left(v, lo)
            if c > best:
                best, best_pos = c, lo
    else:
        for hi in reversed(v):  # descending: strict '>' keeps the largest (most upstream on '-') on ties
            c = bisect.bisect_right(v, hi) - bisect.bisect_left(v, hi - window)
            if c > best:
                best, best_pos = c, hi
    return best_pos


def parse_models(path):
    """transcript_id -> dict(chrom, strand, exons, start, end, introns) from exon rows (0-based half-open)."""
    tx = collections.OrderedDict()
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        d = tx.setdefault(tid, dict(chrom=f[0], strand=f[6], exons=[]))
        d['exons'].append((int(f[3]) - 1, int(f[4])))
    for d in tx.values():
        d['exons'].sort()
        d['start'], d['end'] = d['exons'][0][0], d['exons'][-1][1]
        d['introns'] = tuple((a[1], b[0]) for a, b in zip(d['exons'], d['exons'][1:]))
    return tx


def five(m):
    return m['start'] if m['strand'] == '+' else m['end']


def simulate(models, by_chain, window):
    """tid -> new 5' end for every stranded multi-exon model (unchanged when it has no exact reads)."""
    out = {}
    for tid, m in models.items():
        if not m['introns'] or m['strand'] not in ('+', '-'):
            continue
        ex = by_chain.get((m['chrom'], m['introns']), [])
        pos = [s if m['strand'] == '+' else e for s, e in ex]
        out[tid] = densest_five_prime(pos, window, m['strand']) if pos else five(m)
    return out


def metrics(models, new5, ref_by_chain, ref_by_intron):
    within = moved_in = moved_out = guard = changed = n_fsm = 0
    diffs = []
    for tid, p in new5.items():
        m = models[tid]
        cur = five(m)
        changed += p != cur
        refs = ref_by_chain.get((m['chrom'], m['strand'], m['introns']))
        if refs:
            n_fsm += 1
            d_new = min(abs(p - t) for t in refs)
            d_old = min(abs(cur - t) for t in refs)
            diffs.append(d_new)
            within += d_new <= TOL
            moved_in += d_old > TOL and d_new <= TOL
            moved_out += d_old <= TOL and d_new > TOL
        near = set()
        for it in m['introns']:
            near |= ref_by_intron.get((m['chrom'], m['strand'], it), set())
        guard += any(abs(p - t) <= TOL for t in near)
    return dict(n_multi_stranded=len(new5), n_changed=changed, n_fsm_chain=n_fsm, n_within50=within,
                n_moved_in=moved_in, n_moved_out=moved_out,
                median_abs_diff=statistics.median(diffs) if diffs else 'NA', n_guard_within50=guard)


def write_gtf(models_gtf, models, new5, path):
    with open(path, 'w') as out:
        for line in open(models_gtf):
            f = line.rstrip('\n').split('\t')
            if len(f) >= 9 and f[2] in ('transcript', 'exon'):
                tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
                if tid in new5:
                    m = models[tid]
                    if m['strand'] == '+' and int(f[3]) - 1 == m['start']:
                        f[3] = str(new5[tid] + 1)
                    elif m['strand'] == '-' and int(f[4]) == m['end']:
                        f[4] = str(new5[tid])
                    line = '\t'.join(f) + '\n'
            out.write(line)


def main(models_gtf, ref_gtf, reads_pkl, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    models = parse_models(models_gtf)
    ref_by_chain = collections.defaultdict(list)
    ref_by_intron = collections.defaultdict(set)
    for r in parse_models(ref_gtf).values():
        if r['introns']:
            ref_by_chain[(r['chrom'], r['strand'], r['introns'])].append(five(r))
            for it in r['introns']:
                ref_by_intron[(r['chrom'], r['strand'], it)].add(five(r))
    by_chain = collections.defaultdict(list)
    reads = pickle.load(open(reads_pkl, 'rb'))
    for (_name, chrom, s, e, _rev, _ts, chain) in reads:
        if chain:
            by_chain[(chrom, tuple(chain))].append((s, e))
    print(f'models={len(models)} primaries={len(reads)}', file=sys.stderr)
    baseline = {tid: five(m) for tid, m in models.items() if m['introns'] and m['strand'] in ('+', '-')}
    rows = [dict(W='none', **metrics(models, baseline, ref_by_chain, ref_by_intron))]
    for w in GRID:
        new5 = simulate(models, by_chain, w)
        rows.append(dict(W=w, **metrics(models, new5, ref_by_chain, ref_by_intron)))
        write_gtf(models_gtf, models, new5, f'{out_dir}/tss_W{w}.gtf')
    cols = list(rows[0].keys())
    with open(f'{out_dir}/summary.tsv', 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')
    graded = [r for r in rows if r['W'] != 'none']
    chosen = min(graded, key=lambda r: (-r['n_within50'], r['W']))['W']
    open(f'{out_dir}/chosen_W.txt', 'w').write(f'{chosen}\n')
    print(open(f'{out_dir}/summary.tsv').read() + f'chosen_W={chosen}')


def self_test():
    assert densest_five_prime([100, 500, 505, 510], 25, '+') == 500
    assert densest_five_prime([2000, 900, 905, 910], 25, '-') == 910
    assert densest_five_prime([100, 110, 500, 510], 25, '+') == 100
    assert densest_five_prime([100, 110, 500, 510], 25, '-') == 510
    assert densest_five_prime([100], 25, '+') == 100
    assert densest_five_prime([100, 900], 25, '+') == 100
    assert densest_five_prime([100, 900], 25, '-') == 900
    assert densest_five_prime([], 25, '+') is None
    assert densest_five_prime([100, 400, 425], 25, '+') == 400
    assert densest_five_prime([100, 400, 425], 24, '+') == 100
    print('self-test OK')


if __name__ == '__main__':
    if sys.argv[1:] == ['--self-test']:
        self_test()
    else:
        main(*sys.argv[1:5])

#!/usr/bin/env python3
"""
Compare transfrag pre-depletion chain abundances between Rustle and StringTie parity JSONL logs.

Reads 'transfrag_pre_depl' and 'path_extracted' events from each tool's parity log
and produces a comparison table grouped by intron chain fingerprint.

Flags of interest:
  FRAG    - one tool has multiple extractions of the same chain, the other has one or zero
  MISS_R  - chain has pre-depl support in ST (or ST extracted it) but Rustle never extracted it
  MISS_ST - chain has pre-depl support in Rustle (or Rustle extracted it) but ST never extracted it

Usage:
  python3 bench/predepl_chain_compare.py rustle.jsonl st.jsonl [--chrom X] [--range S-E]
"""

import sys
import json
import argparse
from collections import defaultdict


def parse_jsonl(path, steps, chrom_filter=None, range_filter=None):
    events = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            try:
                e = json.loads(line)
            except json.JSONDecodeError:
                continue
            if e.get('step') not in steps:
                continue
            if chrom_filter and e.get('chrom') != chrom_filter:
                continue
            if range_filter:
                lo, hi = range_filter
                ev_s = e.get('start', 0)
                ev_e = e.get('end', 0)
                if ev_e < lo or ev_s > hi:
                    continue
            events.append(e)
    return events


def summarize_jsonl(path, chrom_filter=None, range_filter=None):
    steps = {'transfrag_pre_depl', 'path_extracted'}
    events = parse_jsonl(path, steps, chrom_filter, range_filter)

    # chain fingerprint → pre-depletion totals
    pre_depl: dict = defaultdict(lambda: {
        'abund': 0.0, 'n_tfs': 0, 'strand': '.',
    })
    # chain fingerprint → extraction totals
    extracted: dict = defaultdict(lambda: {
        'longcov': 0.0, 'cov': 0.0, 'n_ext': 0, 'entry_abund': 0.0,
        'sources': [],
    })

    for e in events:
        step = e.get('step')
        p = e.get('payload', {})
        chain = p.get('introns', '')
        strand = e.get('strand', '.')

        if step == 'transfrag_pre_depl':
            pre_depl[chain]['abund'] += float(p.get('abund', 0))
            pre_depl[chain]['n_tfs'] += 1
            if pre_depl[chain]['strand'] == '.':
                pre_depl[chain]['strand'] = strand

        elif step == 'path_extracted':
            extracted[chain]['longcov'] += float(p.get('longcov', 0))
            extracted[chain]['cov'] += float(p.get('cov', 0))
            extracted[chain]['entry_abund'] += float(p.get('entry_abund', 0))
            extracted[chain]['n_ext'] += 1
            src = p.get('source', '')
            extracted[chain]['sources'].append(src)

    return pre_depl, extracted


def chain_summary(introns_str):
    """Return (n_introns, first_intron, last_intron) for display."""
    if not introns_str:
        return 0, '', ''
    parts = introns_str.split(',')
    n = len(parts)
    return n, parts[0], parts[-1]


def main():
    ap = argparse.ArgumentParser(
        description='Compare pre-depletion chain abundances between Rustle and StringTie parity logs'
    )
    ap.add_argument('rustle_jsonl', help='Rustle parity JSONL (RUSTLE_PARITY_LOG output)')
    ap.add_argument('st_jsonl', help='StringTie parity JSONL (STRINGTIE_PARITY_LOG output)')
    ap.add_argument('--chrom', '-c', default=None, help='Filter to this chromosome')
    ap.add_argument('--range', '-r', default=None, help='Genomic range filter, e.g. 38000000-39000000')
    ap.add_argument('--min-abund', '-m', type=float, default=2.0,
                    help='Min pre-depl abund in either tool to include chain (default: 2.0)')
    ap.add_argument('--show-all', '-a', action='store_true',
                    help='Show all chains above min-abund, not just flagged ones')
    ap.add_argument('--sort-by', choices=['abund', 'diff', 'frag'], default='abund',
                    help='Sort order: abund=max pre-depl, diff=|R-ST| pre-depl, frag=fragmentation severity')
    args = ap.parse_args()

    range_filter = None
    if args.range:
        lo, hi = args.range.split('-')
        range_filter = (int(lo), int(hi))

    r_pre, r_ext = summarize_jsonl(args.rustle_jsonl, args.chrom, range_filter)
    s_pre, s_ext = summarize_jsonl(args.st_jsonl, args.chrom, range_filter)

    all_chains = (
        set(r_pre.keys()) | set(s_pre.keys()) |
        set(r_ext.keys()) | set(s_ext.keys())
    )

    rows = []
    for chain in all_chains:
        r_abund = r_pre[chain]['abund']
        s_abund = s_pre[chain]['abund']
        r_tf = r_pre[chain]['n_tfs']
        s_tf = s_pre[chain]['n_tfs']

        r_n_ext = r_ext[chain]['n_ext']
        s_n_ext = s_ext[chain]['n_ext']
        r_lc = r_ext[chain]['longcov']
        s_lc = s_ext[chain]['longcov']
        r_ea = r_ext[chain]['entry_abund']

        max_abund = max(r_abund, s_abund, r_lc, s_lc)
        if max_abund < args.min_abund:
            continue

        # Determine flags
        flags = []
        # Fragmentation: one tool uses many transfrags/extractions, other uses one
        if (r_n_ext > 1 and s_n_ext <= 1 and r_n_ext > 0) or \
           (s_n_ext > 1 and r_n_ext <= 1 and s_n_ext > 0):
            flags.append('FRAG')
        if r_tf > 1 and s_tf <= 1:
            flags.append('FRAG_PREDEPL')  # fragmented even at pre-depl stage

        # Missing: ST extracted but Rustle didn't
        if s_n_ext > 0 and r_n_ext == 0:
            flags.append('MISS_R')
        # Missing: Rustle extracted but ST didn't
        if r_n_ext > 0 and s_n_ext == 0:
            flags.append('MISS_ST')

        # Coverage divergence: both extracted but longcov totals differ significantly
        if r_n_ext > 0 and s_n_ext > 0:
            ratio = max(r_lc, s_lc) / max(min(r_lc, s_lc), 0.01)
            if ratio > 2.0:
                flags.append('COV_DIFF')

        if not args.show_all and not flags:
            continue

        n_int, first_int, last_int = chain_summary(chain)
        rows.append({
            'chain': chain,
            'n_int': n_int,
            'strand': r_pre[chain]['strand'] or s_pre[chain]['strand'],
            'first_int': first_int,
            'last_int': last_int,
            'r_abund': r_abund, 'r_tf': r_tf,
            's_abund': s_abund, 's_tf': s_tf,
            'r_n_ext': r_n_ext, 'r_lc': r_lc,
            's_n_ext': s_n_ext, 's_lc': s_lc,
            'r_ea': r_ea,
            'flags': '|'.join(flags),
            '_max': max_abund,
            '_diff': abs(r_lc - s_lc),
            '_frag': max(r_n_ext, s_n_ext),
        })

    # Sort
    if args.sort_by == 'diff':
        rows.sort(key=lambda x: x['_diff'], reverse=True)
    elif args.sort_by == 'frag':
        rows.sort(key=lambda x: x['_frag'], reverse=True)
    else:
        rows.sort(key=lambda x: x['_max'], reverse=True)

    # Print header
    hdr = (f"{'last_intron':<28} {'ni':>3} {'st':>3} "
           f"{'R_predepl':>10} {'R_tf':>5} {'ST_predepl':>10} {'ST_tf':>5} "
           f"{'R_ext':>6} {'R_lc':>7} {'ST_ext':>6} {'ST_lc':>7}  flags")
    print(hdr)
    print('-' * len(hdr))

    for row in rows:
        li = row['last_int']
        li_short = (li[:26] + '..') if len(li) > 28 else li
        print(
            f"{li_short:<28} {row['n_int']:>3} {row['strand']:>3} "
            f"{row['r_abund']:>10.2f} {row['r_tf']:>5} {row['s_abund']:>10.2f} {row['s_tf']:>5} "
            f"{row['r_n_ext']:>6} {row['r_lc']:>7.2f} {row['s_n_ext']:>6} {row['s_lc']:>7.2f}  "
            f"{row['flags']}"
        )

    print(f"\n{len(rows)} chains shown (min_abund={args.min_abund}, "
          f"show_all={args.show_all})")

    # Summarize flags
    flag_counts: dict = defaultdict(int)
    for row in rows:
        for f in row['flags'].split('|'):
            if f:
                flag_counts[f] += 1
    if flag_counts:
        print("\nFlag summary:")
        for f, cnt in sorted(flag_counts.items(), key=lambda x: -x[1]):
            print(f"  {f:<20} {cnt}")


if __name__ == '__main__':
    main()

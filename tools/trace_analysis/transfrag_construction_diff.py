#!/usr/bin/env python3
"""transfrag_construction_diff.py — cross-tool diff of CONSTRUCTED transfrags.

Both rustle and StringTie emit a `transfrag_define` parity event for every
finalized transfrag, at construction time (immediately after
`process_transfrags`, BEFORE seeding/flow depletion), for ALL transfrags,
with coordinate-aligned intron chains (donor_end+1 .. acceptor convention;
rustle 0-based acceptor == ST 1-based−1, so the emitted strings match —
unlike raw node-coord dumps which have a 0/1-based off-by-one).

This tool joins the two logs by (strand, canonical intron-chain) and reports
where transfrag CONSTRUCTION diverges — the established root cause of the
rustle<StringTie post-flow coverage gap (rustle does not build the
interior-start long-range spanning transfrags StringTie builds; see
project_cov_estimation_divergence_2026_05_15).

Buckets:
  ST_ONLY        chain constructed by StringTie, absent in rustle
                 (the key class — missing capacity-chord transfrags)
  RUSTLE_ONLY    chain constructed by rustle, absent in StringTie
  ABUND_DIVERGE  same chain both sides, abundance differs > --abund-tol
  MATCH          same chain + abundance within tol

Single-exon transfrags (n_introns==0) are bucketed separately (noisy).

Usage:
  transfrag_construction_diff.py --rustle r.jsonl --stringtie st.jsonl \
      [--range LO-HI] [--min-introns N] [--abund-tol 0.5] [--top 40]

Generate the inputs over the same locus, e.g. STRG.15.1:
  RUSTLE_PARITY_LOG=/tmp/r.jsonl RUSTLE_PARITY_FILTER_CHROM=NC_073243.2 \
    RUSTLE_PARITY_FILTER_RANGE=16598646-16622552 RUSTLE_FLOW_RESIDUAL_SE=1 \
    ./Rustle/target/dev-opt/rustle -L GGO_19.bam -o /tmp/r.gtf
  STRINGTIE_PARITY_LOG=/tmp/st.jsonl STRINGTIE_PARITY_FILTER_CHROM=NC_073243.2 \
    STRINGTIE_PARITY_FILTER_RANGE=16598646-16622552 \
    ./stringtie/stringtie -L GGO_19.bam -o /tmp/st.gtf
"""
import argparse, json, sys
from collections import defaultdict, Counter


def load(path):
    """Return dict keyed (strand, introns) -> list of records.
    Each record: dict(span_lo, span_hi, abund, n_introns, n_nodes, longread)."""
    out = defaultdict(list)
    with open(path) as f:
        for line in f:
            try:
                ev = json.loads(line)
            except Exception:
                continue
            if ev.get('step') != 'transfrag_define':
                continue
            p = ev.get('payload', {})
            introns = p.get('introns', '')
            key = (ev.get('strand', '?'), introns)
            out[key].append(dict(
                lo=int(ev.get('start', 0)),
                hi=int(ev.get('end', 0)),
                abund=float(p.get('abund', 0.0)),
                n_introns=int(p.get('n_introns', 0)),
                n_nodes=int(p.get('n_nodes', 0)) if p.get('n_nodes') is not None else -1,
                longread=p.get('longread'),
            ))
    return out


def agg(records):
    """Collapse duplicate (strand,chain) records: sum abund, max span, count."""
    tot_ab = sum(r['abund'] for r in records)
    lo = min(r['lo'] for r in records)
    hi = max(r['hi'] for r in records)
    ni = records[0]['n_introns']
    nn = max(r['n_nodes'] for r in records)
    return tot_ab, lo, hi, ni, nn, len(records)


def in_range(lo, hi, rng):
    if rng is None:
        return True
    return not (hi < rng[0] or lo > rng[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rustle', required=True)
    ap.add_argument('--stringtie', required=True)
    ap.add_argument('--range', help='LO-HI genomic filter (span overlap)')
    ap.add_argument('--min-introns', type=int, default=1,
                    help='ignore transfrags with fewer introns (default 1; '
                         '0 includes single-exon)')
    ap.add_argument('--abund-tol', type=float, default=0.5)
    ap.add_argument('--top', type=int, default=40)
    args = ap.parse_args()

    rng = None
    if args.range:
        a, b = args.range.split('-')
        rng = (int(a), int(b))

    R = load(args.rustle)
    S = load(args.stringtie)
    keys = set(R) | set(S)

    buckets = defaultdict(list)
    for k in keys:
        strand, introns = k
        rr = R.get(k)
        ss = S.get(k)
        recs = rr or ss
        _, lo, hi, ni, _, _ = agg(recs)
        if ni < args.min_introns:
            buckets['_single_or_short'].append(k)
            continue
        if not in_range(lo, hi, rng):
            continue
        if rr and ss:
            ra = agg(rr)
            sa = agg(ss)
            if abs(ra[0] - sa[0]) > args.abund_tol:
                buckets['ABUND_DIVERGE'].append((k, ra, sa))
            else:
                buckets['MATCH'].append((k, ra, sa))
        elif ss and not rr:
            buckets['ST_ONLY'].append((k, agg(ss)))
        else:
            buckets['RUSTLE_ONLY'].append((k, agg(rr)))

    n_st_only = len(buckets['ST_ONLY'])
    n_r_only = len(buckets['RUSTLE_ONLY'])
    n_div = len(buckets['ABUND_DIVERGE'])
    n_match = len(buckets['MATCH'])
    print(f"transfrag_define cross-tool diff "
          f"(range={args.range or 'all'}, min_introns={args.min_introns})")
    print(f"  MATCH={n_match}  ABUND_DIVERGE={n_div}  "
          f"ST_ONLY={n_st_only}  RUSTLE_ONLY={n_r_only}")
    print()

    # ST_ONLY = transfrags StringTie builds that rustle never constructs.
    # Sort by span width (long-range / interior-start chords first) then abund.
    so = sorted(buckets['ST_ONLY'],
                key=lambda x: (-(x[1][2] - x[1][1]), -x[1][0]))
    print(f"== ST_ONLY (rustle never constructs) — top {args.top} by span ==")
    print(f"  {'strand':6} {'span':25} {'nI':>3} {'nN':>3} {'abund':>8}  introns[:60]")
    for (k, a) in so[:args.top]:
        strand, introns = k
        tot_ab, lo, hi, ni, nn, cnt = a
        print(f"  {strand:6} {f'{lo}-{hi}':25} {ni:>3} {nn:>3} "
              f"{tot_ab:>8.3f}  {introns[:60]}")

    if buckets['ABUND_DIVERGE']:
        print()
        print(f"== ABUND_DIVERGE (same chain, |Δabund|>{args.abund_tol}) — "
              f"top {args.top} by |Δ| ==")
        dv = sorted(buckets['ABUND_DIVERGE'],
                    key=lambda x: -abs(x[1][0] - x[2][0]))
        print(f"  {'strand':6} {'span':25} {'nI':>3} "
              f"{'R_ab':>8} {'ST_ab':>8} {'Δ':>8}")
        for (k, ra, sa) in dv[:args.top]:
            strand, introns = k
            print(f"  {strand:6} {f'{ra[1]}-{ra[2]}':25} {ra[3]:>3} "
                  f"{ra[0]:>8.3f} {sa[0]:>8.3f} {sa[0]-ra[0]:>+8.3f}")

    print()
    print(f"(single/short transfrags ignored: "
          f"{len(buckets['_single_or_short'])} chains)")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""seed_only_analysis.py — dissect the `seed_only` missed-ref class.

A `seed_only` ref has its EXACT intron chain registered as a `transfrag_seed`
(flow extraction WAS given the right seed) but no `path_extracted` carries that
chain — flow produced a different path from the same seed. This tool, for each
seed_only ref, aligns:

  ref chain  vs  the seed (confirm registered, abund)
             vs  checktrf_enter/checktrf_result (was it deferred + why)
             vs  the best-overlapping path_extracted chain (what flow built
                 instead) — longest exact consecutive run + first divergence

and classifies the divergence:

  TERMINAL        extracted chain is the ref chain ± contiguous end junctions
                  (ref is a sub/superchain of extracted) — boundary/RI at ends
  INTERNAL_REROUTE shared prefix then flow takes a different junction mid-chain
  READTHR_SKIP    ref chain seen in checktrf with a readthr* reason/outcome
  SHADOWED        a higher-cov path_extracted's chain ⊇ ref chain (sibling won)
  NObest          no overlapping path_extracted (seed consumed, nothing emitted)

Intron-string convention matches no_seed_analysis.py / parity emits:
`donor_exon_end+1 - acceptor_exon_start-1`.

Usage:
  seed_only_analysis.py --parity /tmp/ns.jsonl --ref GGO_19.gtf \
      --refmap /tmp/ns_cmp.ns.gtf.refmap [--tid STRG.93.3]
"""
import argparse, json, re, sys
from collections import defaultdict


def load_ref(ref):
    ex = defaultdict(list)
    info = {}
    for line in open(ref):
        if line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        m = re.search(r'transcript_id "([^"]+)"', f[8])
        if not m:
            continue
        ex[m.group(1)].append((int(f[3]), int(f[4])))
        info[m.group(1)] = (f[0], f[6])
    for v in ex.values():
        v.sort()
    return ex, info


def chain(exons):
    return [f'{exons[i][1]+1}-{exons[i+1][0]-1}' for i in range(len(exons) - 1)]


def matched_set(refmap):
    s = set()
    f = open(refmap)
    next(f, None)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3 and p[2] == '=':
            s.add(p[1])
    return s


def longest_run(a, b):
    """longest exact consecutive run of list a inside list b; returns
    (len, i_in_a, j_in_b)."""
    best = (0, 0, 0)
    for i in range(len(a)):
        for j in range(len(b)):
            k = 0
            while i + k < len(a) and j + k < len(b) and a[i + k] == b[j + k]:
                k += 1
            if k > best[0]:
                best = (k, i, j)
    return best


def index(parity):
    seed_by_chain = {}
    pe = []          # (span_lo, span_hi, strand, chain_list, payload)
    checktrf = defaultdict(list)   # chain_str -> [(kind, payload)]
    for line in open(parity):
        try:
            ev = json.loads(line)
        except Exception:
            continue
        st = ev.get('step', '')
        p = ev.get('payload', {})
        intr = p.get('introns', '')
        if not intr:
            continue
        cl = intr.split(',')
        if st == 'transfrag_seed':
            seed_by_chain[intr] = p
        elif st == 'path_extracted':
            pe.append((ev.get('start', 0), ev.get('end', 0),
                       ev.get('strand', '?'), cl, p))
        elif st in ('checktrf_enter', 'checktrf_result'):
            checktrf[intr].append((st, p))
    return seed_by_chain, pe, checktrf


def classify(rc, pe, strand, lo, hi, checktrf_hits):
    rcs = ','.join(rc)
    for kind, p in checktrf_hits:
        r = (p.get('reason') or p.get('outcome') or '')
        if 'readthr' in r or 'runon' in r:
            return 'READTHR_SKIP', f'{kind}:{r}', None
    # best-overlapping path_extracted
    cands = [(s, e, ch, pp) for (s, e, st, ch, pp) in pe
             if st == strand and not (e < lo or s > hi)]
    if not cands:
        return 'NObest', 'no overlapping path_extracted', None
    scored = []
    for s, e, ch, pp in cands:
        L, i, j = longest_run(rc, ch)
        scored.append((L, s, e, ch, pp, i, j))
    scored.sort(key=lambda x: -x[0])
    L, s, e, ch, pp, i, j = scored[0]
    chs = ','.join(ch)
    note = (f'best run={L}/{len(rc)} src={pp.get("source")} '
            f'cov={pp.get("cov")} nx={pp.get("nexons")} '
            f'extr_span={s}-{e}')
    if rcs in chs:
        # ref chain is a contiguous subchain of an extracted chain
        if pp.get('cov', 0) and pp['cov'] > 1:
            return 'SHADOWED', note + ' (ref ⊂ extracted)', (ch, pp)
        return 'TERMINAL', note + ' (ref ⊂ extracted, lowcov)', (ch, pp)
    if chs in rcs:
        return 'TERMINAL', note + ' (extracted ⊂ ref — truncated)', (ch, pp)
    # Decompose: best run covers ref[i : i+L]. Missing = prefix ref[:i] and
    # suffix ref[i+L:]. Pure prefix-miss or suffix-miss = terminal divergence;
    # both-ends present but inner gap = true internal reroute.
    pre = i               # ref junctions missing at 5' of the run
    suf = len(rc) - (i + L)
    if pre and not suf:
        miss = rc[i - 1] if i - 1 < len(rc) else rc[0]
        return ('TERMINAL_5P',
                note + f' missing {pre} ref 5′ jct(s); first ref jct '
                f'{rc[0]} absent from flow path', (ch, pp))
    if suf and not pre:
        return ('TERMINAL_3P',
                note + f' missing {suf} ref 3′ jct(s); last ref jct '
                f'{rc[-1]} absent from flow path', (ch, pp))
    # both ends missing OR a genuine inner divergence
    div = i + L
    rj = rc[div] if div < len(rc) else '(end)'
    xj = ch[j + L] if (j + L) < len(ch) else '(end)'
    return ('INTERNAL_REROUTE',
            note + f' run=ref[{i}:{i+L}] diverge@ref[{div}]={rj} '
            f'flow_took={xj} (pre_miss={pre} suf_miss={suf})', (ch, pp))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--parity', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--refmap', required=True)
    ap.add_argument('--tid')
    args = ap.parse_args()

    ex, info = load_ref(args.ref)
    matched = matched_set(args.refmap)
    seed_by_chain, pe, checktrf = index(args.parity)

    seed_chains = set(seed_by_chain)
    pe_chains = set(','.join(c) for (_, _, _, c, _) in pe)

    rows = []
    for t in ex:
        if t in matched or len(ex[t]) < 2:
            continue
        rc = chain(ex[t])
        rcs = ','.join(rc)
        if rcs in pe_chains or rcs not in seed_chains:
            continue
        rows.append(t)
    if args.tid:
        rows = [args.tid]

    print(f"seed_only refs: {len(rows)}\n")
    print(f"{'tid':13s} {'st':2s} {'nj':>3s} {'seedAb':>6s} {'class':16s} note")
    print('-' * 110)
    from collections import Counter
    cc = Counter()
    for t in sorted(rows):
        chrom, strand = info[t]
        rc = chain(ex[t])
        rcs = ','.join(rc)
        lo, hi = ex[t][0][0], ex[t][-1][1]
        sp = seed_by_chain.get(rcs, {})
        cls, note, _ = classify(rc, pe, strand, lo, hi,
                                checktrf.get(rcs, []))
        cc[cls] += 1
        print(f"  {t:13s} {strand:2s} {len(rc):>3d} "
              f"{str(sp.get('abund','?')):>6s} {cls:16s} {note}")
    print()
    for k, v in cc.most_common():
        print(f"  {v:3d}  {k}")


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""no_seed_analysis.py — classify missed reference transcripts by the pipeline
stage at which their intron chain disappears, with deep analysis of the
"no_seed" class (chains that never reach transfrag construction).

Inputs:
  --parity   parity_decisions JSONL log (RUSTLE_PARITY_LOG output)
  --ref      reference GTF (the StringTie-de-novo GGO_19.gtf)
  --refmap   gffcompare *.refmap (to determine which refs were matched '=')

Stages a missed multi-exon ref can fall into:
  - matched              : '=' in refmap (not missed; excluded)
  - path_extracted       : exact chain reached path_extracted (killed later;
                           cross-ref with pred_kill stage= to find the killer)
  - seed_only            : chain in transfrag_seed but no path_extracted
  - no_seed              : chain never in transfrag_seed
      sub: in_define_not_seed  — transfrag built, seeding gate dropped it
      sub: not_even_define     — no transfrag with the exact chain

For each no_seed ref it additionally reports:
  - junction acceptance: of the ref's N junctions, how many are accepted
    graph edges (junction_accept accepted=true), rejected (with reason), or
    absent (no junction_accept event — junction never seen in any read CIGAR)
  - longest exact consecutive sub-chain present in any transfrag at the locus,
    and the breakpoint junction (where the combination first diverges)

Intron-string convention (all steps): "donor_exon_end+1 - acceptor_exon_start-1"
(1-based inclusive intron span). junction_accept uses raw GTF
(exon_end, next_exon_start).

Usage:
  no_seed_analysis.py --parity /tmp/p.jsonl --ref GGO_19.gtf \
      --refmap cmp.gtf.refmap [--only-no-seed] [--tid STRG.135.3]
"""
import argparse, json, re, sys
from collections import defaultdict


def load_ref_exons(ref_gtf):
    ex = defaultdict(list)
    info = {}
    with open(ref_gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fs = line.split('\t')
            if len(fs) < 9 or fs[2] != 'exon':
                continue
            m = re.search(r'transcript_id "([^"]+)"', fs[8])
            if not m:
                continue
            tid = m.group(1)
            ex[tid].append((int(fs[3]), int(fs[4])))
            info[tid] = (fs[0], fs[6])
    for v in ex.values():
        v.sort()
    return ex, info


def chain_str(exons):
    return ','.join(f'{exons[i][1]+1}-{exons[i+1][0]-1}'
                     for i in range(len(exons) - 1))


def matched_set(refmap):
    s = set()
    with open(refmap) as f:
        next(f, None)
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 3 and p[2] == '=':
                s.add(p[1])
    return s


def index_parity(parity):
    """Returns: seed_chains, define_chains, pe_chains (sets of intron strings),
    junction_accept dict {(strand,donor,acceptor):(accepted,reason)},
    tf_by_strand {strand:[(lo,hi,[intron_str,...])]}."""
    seed, define, pe = set(), set(), set()
    ja = {}
    tf = defaultdict(list)
    define_meta = {}  # intron_str -> seed-eligibility payload fields
    with open(parity) as f:
        for line in f:
            try:
                ev = json.loads(line)
            except Exception:
                continue
            st = ev.get('step', '')
            p = ev.get('payload', {})
            intr = p.get('introns', '')
            if st == 'junction_accept':
                ja[(ev.get('strand', '?'), ev.get('start', 0), ev.get('end', 0))] = (
                    p.get('accepted'), p.get('reason'))
                continue
            if not intr:
                continue
            if st == 'transfrag_seed':
                seed.add(intr)
            elif st == 'transfrag_define':
                define.add(intr)
                tf[ev.get('strand', '?')].append(
                    (ev.get('start', 0), ev.get('end', 0), intr.split(',')))
                if intr not in define_meta:
                    define_meta[intr] = {
                        k: p.get(k) for k in
                        ('abund', 'longread', 'trflong_seed', 'weak',
                         'usepath', 'n_nodes', 'guide')
                    }
            elif st == 'path_extracted':
                pe.add(intr)
    return seed, define, pe, ja, tf, define_meta


def longest_run(ref_chain, tf_chain):
    best, best_i = 0, 0
    n = len(ref_chain)
    for i in range(n):
        for j in range(len(tf_chain)):
            k = 0
            while (i + k < n and j + k < len(tf_chain)
                   and ref_chain[i + k] == tf_chain[j + k]):
                k += 1
            if k > best:
                best, best_i = k, i
    return best, best_i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--parity', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--refmap', required=True)
    ap.add_argument('--only-no-seed', action='store_true')
    ap.add_argument('--tid', help='deep-dive a single ref id')
    args = ap.parse_args()

    ref_ex, ref_info = load_ref_exons(args.ref)
    matched = matched_set(args.refmap)
    seed_c, define_c, pe_c, ja, tf, define_meta = index_parity(args.parity)

    missed = [t for t in ref_ex
              if t not in matched and len(ref_ex[t]) >= 2]
    rows = []
    for tid in missed:
        ex = ref_ex[tid]
        cs = chain_str(ex)
        if cs in pe_c:
            cls = 'path_extracted'
        elif cs in seed_c:
            cls = 'seed_only'
        elif cs in define_c:
            cls = 'no_seed:in_define_not_seed'
        else:
            cls = 'no_seed:not_even_define'
        rows.append((tid, cls))

    if not args.only_no_seed and not args.tid:
        from collections import Counter
        c = Counter(r[1] for r in rows)
        print(f"Missed multi-exon refs: {len(rows)}")
        for k, v in sorted(c.items(), key=lambda x: -x[1]):
            print(f"  {v:4d}  {k}")
        print()

    ns = [t for t, c in rows if c.startswith('no_seed')]
    if args.tid:
        ns = [args.tid]

    print(f"{'tid':13s} {'chr':3s} {'st':2s} {'nj':>3s} {'acc':>3s} "
          f"{'rej':>3s} {'abs':>3s} {'run':>3s} {'break_idx':>9s} note")
    print('-' * 100)
    for tid in sorted(ns):
        ex = ref_ex[tid]
        chrom, strand = ref_info[tid]
        juncs = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
        rc = chain_str(ex).split(',')
        acc = rej = absent = 0
        rej_r = []
        first_abs = None
        for (d, a) in juncs:
            k = (strand, d, a)
            if k in ja:
                ok, rsn = ja[k]
                if ok:
                    acc += 1
                else:
                    rej += 1
                    rej_r.append(f'{d}-{a}:{rsn}')
            else:
                absent += 1
                if first_abs is None:
                    first_abs = f'{d}-{a}'
        lo, hi = ex[0][0], ex[-1][1]
        best, best_i = 0, 0
        for (s, en, tc) in tf.get(strand, []):
            if en < lo or s > hi:
                continue
            r, ri = longest_run(rc, tc)
            if r > best:
                best, best_i = r, ri
        bp = best_i + best
        if rej:
            note = 'REJECTED_JUNC: ' + ';'.join(rej_r[:2])
        elif absent:
            note = f'ABSENT_JUNC first={first_abs}'
        elif best >= len(rc):
            meta = define_meta.get(','.join(rc), {})
            note = ('EXACT_CHAIN_IN_TF drop: '
                    f"trflong_seed={meta.get('trflong_seed')} "
                    f"weak={meta.get('weak')} usepath={meta.get('usepath')} "
                    f"abund={meta.get('abund')}")
        else:
            bpj = rc[bp] if bp < len(rc) else '(end)'
            note = f'COMBINATORIAL break@{bpj}'
        print(f"  {tid:13s} {chrom[-3:]:3s} {strand:2s} {len(juncs):>3d} "
              f"{acc:>3d} {rej:>3d} {absent:>3d} {best:>3d} {bp:>9d} {note}")


if __name__ == '__main__':
    main()

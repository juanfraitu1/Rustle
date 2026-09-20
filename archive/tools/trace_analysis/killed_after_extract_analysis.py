#!/usr/bin/env python3
"""killed_after_extract_analysis.py — attribute the killing filter for missed
refs whose EXACT intron chain reached `path_extracted` (flow assembled them
fine — no long_max_flow wall) but were removed before GTF emit.

Method: a chain's life is the ordered sequence of `pred_filter_stage` snapshots
(payload.introns == chain) it appears in. The killing stage = the first stage
in pipeline order where it is ABSENT after having been present. Cross-checked
against `pred_kill` events at that stage by (cov, nexons) where available, and
against `path_emit_pre_write` (final survivors).

Inputs: parity JSONL, ref GTF, gffcompare refmap.

Usage:
  killed_after_extract_analysis.py --parity /tmp/mc.jsonl --ref GGO_19.gtf \
      --refmap /tmp/mc_cmp.mc.gtf.refmap [--tid STRG.x.y]
"""
import argparse, json, re
from collections import defaultdict, Counter


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


def chain_str(e):
    return ','.join(f'{e[i][1]+1}-{e[i+1][0]-1}' for i in range(len(e) - 1))


def matched_set(refmap):
    s = set()
    f = open(refmap)
    next(f, None)
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) >= 3 and p[2] == '=':
            s.add(p[1])
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--parity', required=True)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--refmap', required=True)
    ap.add_argument('--tid')
    args = ap.parse_args()

    ex, info = load_ref(args.ref)
    matched = matched_set(args.refmap)

    pe_chains = set()
    emit_chains = set()
    # stage order: first-seen order of pred_filter_stage 'stage' values
    stage_order = []
    stage_seen = set()
    stages_by_chain = defaultdict(list)   # chain -> [stage,...] in file order
    kills_by_stage = defaultdict(list)    # stage -> [(cov,nexons,reason,killer)]
    for line in open(args.parity):
        try:
            ev = json.loads(line)
        except Exception:
            continue
        st = ev.get('step', '')
        p = ev.get('payload', {})
        if st == 'path_extracted':
            if p.get('introns'):
                pe_chains.add(p['introns'])
        elif st == 'path_emit_pre_write':
            if p.get('introns'):
                emit_chains.add(p['introns'])
        elif st == 'pred_filter_stage':
            sg = p.get('stage', '?')
            if sg not in stage_seen:
                stage_seen.add(sg)
                stage_order.append(sg)
            intr = p.get('introns', '')
            if intr:
                stages_by_chain[intr].append(sg)
        elif st == 'pred_kill':
            kills_by_stage[p.get('stage', '?')].append(
                (p.get('cov'), p.get('nexons'), p.get('reason'),
                 (p.get('killer_cov'), p.get('killer_nexons'))))
    stage_rank = {s: i for i, s in enumerate(stage_order)}

    rows = []
    for t in ex:
        if t in matched or len(ex[t]) < 2:
            continue
        cs = chain_str(ex[t])
        if cs in pe_chains:
            rows.append((t, cs))
    if args.tid:
        rows = [(args.tid, chain_str(ex[args.tid]))]

    print(f"missed refs that reached path_extracted: {len(rows)}\n")
    print(f"{'tid':13s} {'st':2s} {'nj':>3s} {'lastAlive':32s} "
          f"{'killedAt':32s} note")
    print('-' * 120)
    bucket = Counter()
    for t, cs in sorted(rows):
        chrom, strand = info[t]
        seq = stages_by_chain.get(cs, [])
        if not seq:
            killed = 'NEVER_IN_pred_filter_stage'
            last_alive = '-'
            note = 'reached path_extracted but no pred_filter_stage row'
        else:
            # last stage (by pipeline rank) where chain present
            present = sorted(set(seq), key=lambda s: stage_rank.get(s, 1e9))
            last_alive = present[-1]
            la_rank = stage_rank.get(last_alive, -1)
            after = [s for s in stage_order
                     if stage_rank[s] > la_rank]
            killed = after[0] if after else (
                'SURVIVED_FILTERS' if cs in emit_chains
                else 'LOST_after_last_stage')
            note = ''
            # try to attribute a pred_kill at killed stage
            kk = kills_by_stage.get(killed, [])
            if kk:
                note = f'{len(kk)} pred_kill@{killed}'
        if cs in emit_chains:
            killed = 'EMITTED(dup/term mismatch — not "=")'
        bucket[killed] += 1
        print(f"  {t:13s} {strand:2s} {len(ex[t])-1:>3d} "
              f"{last_alive:32s} {killed:32s} {note}")
    print()
    for k, v in bucket.most_common():
        print(f"  {v:3d}  {k}")


if __name__ == '__main__':
    main()

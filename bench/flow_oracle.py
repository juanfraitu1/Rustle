#!/usr/bin/env python3
"""Flow-extraction parity oracle (sub-project 2, Task 1).

Confirms the decomposition: are StringTie's final transcripts BUILDABLE in
Rustle's graph (all their introns present + accepted by junction_accept)?
If yes -> extraction parity + selection parity is the complete path to ST's
output. If a large fraction of ST chains need junctions Rustle never accepts,
graph/junction parity is a prerequisite outside this sub-project -> abort.

Coordinate convention (verified): GTF intron (exon_end, next_exon_start) maps
DIRECTLY to junction_accept (start, end). No +/-1 offset.

Inputs:
  st.gtf      StringTie final GTF  (stringtie GGO_19.bam -L -o st.gtf)
  ru.gtf      Rustle final GTF
  ru_ja.jsonl Rustle junction_accept parity log
    (RUSTLE_PARITY_LOG=... RUSTLE_PARITY_FILTER_STEPS=junction_accept rustle ...)

Usage: python3 bench/flow_oracle.py st.gtf ru.gtf ru_ja.jsonl
"""
import json
import sys


def parse_gtf(path):
    tx = {}
    for line in open(path):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        strand = f[6]
        start, end = int(f[3]), int(f[4])
        tid = None
        for part in f[8].split(';'):
            part = part.strip()
            if part.startswith('transcript_id'):
                tid = part.split('"')[1]
        tx.setdefault(tid, {'strand': strand, 'exons': []})
        tx[tid]['exons'].append((start, end))
    for t in tx.values():
        t['exons'].sort()
    return tx


def introns(exons):
    return tuple((exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1))


def load_accepted(path):
    acc = set()
    for line in open(path):
        o = json.loads(line)
        if o.get('step') != 'junction_accept':
            continue
        if o['payload'].get('accepted'):
            acc.add((o['start'], o['end']))  # direct convention
    return acc


def main(st_path, ru_path, ja_path):
    st = parse_gtf(st_path)
    ru = parse_gtf(ru_path)
    acc = load_accepted(ja_path)

    ru_chains = {(t['strand'], introns(t['exons']))
                 for t in ru.values() if len(t['exons']) >= 2}
    st_multi = {tid: t for tid, t in st.items() if len(t['exons']) >= 2}

    def missing(t):
        return [i for i in introns(t['exons']) if i not in acc]

    full = sum(1 for t in st_multi.values() if not missing(t))
    partial = len(st_multi) - full

    fn = [(tid, t) for tid, t in st_multi.items()
          if (t['strand'], introns(t['exons'])) not in ru_chains]
    fn_full = sum(1 for _, t in fn if not missing(t))
    fn_partial = len(fn) - fn_full

    print(f"ST multi-intron chains: {len(st_multi)}")
    print(f"  fully-buildable: {full} ({100 * full / len(st_multi):.1f}%)")
    print(f"  partial (>=1 missing junction): {partial}")
    print(f"Rustle-FN ST chains: {len(fn)}")
    print(f"  fully-buildable: {fn_full} ({100 * fn_full / max(1, len(fn)):.1f}%)")
    print(f"  missing-junction: {fn_partial}")


if __name__ == '__main__':
    main(*sys.argv[1:4]) if len(sys.argv) >= 4 else main(
        '/tmp/st.gtf', '/tmp/ru.gtf', '/tmp/ru_ja.jsonl')

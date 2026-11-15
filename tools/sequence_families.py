#!/usr/bin/env python3
"""Families whose copies are MOSTLY a given sequence, e.g. the rDNA unit: one family id per line (register 1691, opt-in).
usage: sequence_families.py COPIES_FA REFERENCE_FA [FRAC=0.5] [--per-copy] > listed.txt
COPIES_FA = PREFIX.fam.copies.fa (headers FAMILY|copy|...). A copy counts when at least FRAC of its bases are covered by
`minimap2 -x asm20` alignments to REFERENCE_FA; a family is listed when strictly more than half of its copies count.
--per-copy prints `name<TAB>covered fraction` for every sequence instead (any FASTA). Needs minimap2 (RUSTLE_MINIMAP2 or PATH)."""
import collections, os, subprocess, sys

def die(msg):
    print(f"[sequence_families] {msg}", file=sys.stderr); sys.exit(2)

def covered(paf_lines):
    """query name -> fraction of the query covered by the union of its alignments"""
    iv = collections.defaultdict(list); qlen = {}
    for l in paf_lines:
        f = l.split('\t')
        if len(f) < 12: continue
        qlen[f[0]] = int(f[1]); iv[f[0]].append((int(f[2]), int(f[3])))
    out = {}
    for q, v in iv.items():
        v.sort(); tot = 0; end = 0
        for a, b in v:
            if b > end: tot += b - max(a, end); end = b
        out[q] = tot / qlen[q]
    return out

def names(fa):
    return [l[1:].split()[0] for l in open(fa) if l.startswith('>')]

if __name__ == '__main__':
    a = [x for x in sys.argv[1:] if x != '--per-copy']
    if len(a) not in (2, 3): die(__doc__)
    try:
        frac = float(a[2]) if len(a) == 3 else 0.5
    except ValueError:
        frac = -1
    if not 0 < frac <= 1: die(f"FRAC must be a number in (0, 1], got '{a[2]}'")
    for p in a[:2]:
        if not os.path.isfile(p) or not os.path.getsize(p): die(f"{p} is missing or empty")
    mm = os.environ.get('RUSTLE_MINIMAP2', 'minimap2')
    try:  # -N 50: one copy may tile the unit in pieces; the union of its alignments is what counts
        r = subprocess.run([mm, '-x', 'asm20', '-N', '50', '-t', os.environ.get('RUSTLE_SEQFAM_THREADS', '8'), a[1], a[0]], capture_output=True, text=True)
    except OSError as e:
        die(f"cannot run {mm}: {e}")
    if r.returncode: die(f"minimap2 failed: {r.stderr[-300:]}")
    cov = covered(r.stdout.splitlines()); all_names = names(a[0])
    if '--per-copy' in sys.argv:
        print('\n'.join(f"{n}\t{cov.get(n, 0):.3f}" for n in all_names)); sys.exit(0)
    n = collections.Counter(); hit = collections.Counter()
    for q in all_names:
        fam = q.split('|')[0]; n[fam] += 1; hit[fam] += cov.get(q, 0) >= frac
    out = sorted(f for f in n if hit[f] * 2 > n[f])
    print('\n'.join(out))
    print(f"[sequence_families] {len(n)} families, {len(out)} listed (most copies >= {frac:g} covered by {os.path.basename(a[1])})", file=sys.stderr)

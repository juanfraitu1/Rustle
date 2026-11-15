#!/usr/bin/env python3
"""PREREG r1627 (docs/PREREG_sibling_chain_2026-11-11.md): score chains a sibling's model projects onto a family member.

  sibling_chains.py REFSEQ_GTF[.gz] DEFAULT_GTF PROJECT_OUTDIR

Classes of projected chains absent from the default GTF (one row per distinct chain):
  S  >= 1 primary record with the exact chain from a molecule with no rival placement (the rule)
  T  primary records with the chain, all from molecules with a rival placement
  X  no primary record, >= 1 AS-tied secondary record
  N  no record
Baseline: default multi-exon chains carried by exactly 2 reads at the target loci (PROJECT_OUTDIR/base.gtf).
Correct = the intron chain is a RefSeq chain on the same contig and strand."""
import sys, re, gzip, collections

def chains(path, want_attr=False):
    ex = collections.defaultdict(list); meta = {}
    op = gzip.open if path.endswith('.gz') else open
    for l in op(path, 'rt'):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('exon', 'transcript'): continue
        t = re.search(r'transcript_id "([^"]*)"', f[8]).group(1)
        if f[2] == 'exon': ex[(f[0], f[6], t)].append((int(f[3]) - 1, int(f[4])))
        elif want_attr: meta[(f[0], f[6], t)] = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
    out = {}
    for k, e in ex.items():
        e.sort()
        c = tuple((a[1], b[0]) for a, b in zip(e, e[1:]))
        if c: out[k] = (k[0], k[1], c)
    return (out, meta) if want_attr else out

def main():
    ref, default, d = sys.argv[1:4]
    truth = set(chains(ref).values())
    dflt = set(chains(default).values())
    proj, meta = chains(f'{d}/projected.gtf', True)
    cls = {}
    for k, c in proj.items():
        if c in dflt: continue
        a = meta[k]; npri, ntied, nclean = int(a['n_primary']), int(a['n_tied']), int(a['n_clean'])
        v = cls.setdefault(c, [0, 0, 0])
        v[0] = max(v[0], npri); v[1] = max(v[1], ntied); v[2] = max(v[2], nclean)
    rows = collections.defaultdict(lambda: [0, 0])
    for c, (npri, ntied, nclean) in cls.items():
        k = 'S' if nclean else 'T' if npri else 'X' if ntied else 'N'
        rows[k][0] += 1; rows[k][1] += c in truth
        if k == 'S':
            kk = 'S: 1 clean read' if nclean == 1 else 'S: 2+ clean reads'
            rows[kk][0] += 1; rows[kk][1] += c in truth
        if npri == 1 and nclean == 1:  # r1628: one read short of the floor, and that read has no rival
            rows['S1'][0] += 1; rows['S1'][1] += c in truth
    base, bmeta = chains(f'{d}/base.gtf', True)
    seen = set()
    for k, c in base.items():
        if c in seen: continue
        seen.add(c)
        r = bmeta[k].get('reads', '')
        kk = 'base: 2-read default chains' if r == '2' else 'base: other default chains'
        rows[kk][0] += 1; rows[kk][1] += c in truth
    print('class\tchains\tcorrect\tshare')
    for k in sorted(rows):
        n, ok = rows[k]; print(f'{k}\t{n}\t{ok}\t{ok / max(n, 1):.3f}')
    n, b = rows['N'], rows['base: 2-read default chains']
    sh = lambda v: v[1] / max(v[0], 1)
    for name in ('S', 'S1'):
        s = rows[name]
        print(f"{name}: B1 correct {s[1]}; B2 share {sh(s):.3f} >= base 2-read {sh(b):.3f}: {'PASS' if sh(s) >= sh(b) else 'FAIL'}; "
              f"B3 >= 1.5 x N {1.5 * sh(n):.3f}: {'PASS' if sh(s) >= 1.5 * sh(n) else 'FAIL'}")

def emit(default, d, out):
    """r1628 as a product: DEFAULT_GTF + one model per projected chain that is not in it and has exactly one primary record,
    from a molecule with no rival placement (opt-in; yield bar failed, quality bar passed on three samples)."""
    dflt = set(chains(default).values()); proj, meta = chains(f'{d}/projected.gtf', True)
    best = {}
    for k, c in proj.items():
        if c in dflt: continue
        a = meta[k]; v = best.setdefault(c, [0, 0, k])
        v[0] = max(v[0], int(a['n_primary'])); v[1] = max(v[1], int(a['n_clean']))
    keep = {v[2] for v in best.values() if v[0] == 1 and v[1] == 1}
    with open(out, 'w') as w:
        w.write(open(default).read())
        for l in open(f'{d}/projected.gtf'):
            f = l.rstrip('\n').split('\t')
            if len(f) < 9 or (f[0], f[6], re.search(r'transcript_id "([^"]*)"', f[8]).group(1)) not in keep: continue
            f[1] = 'sibling'; f[8] = f[8].rstrip() + ' sibling_rule "r1628";'; w.write('\t'.join(f) + '\n')
    print(f'sibling chains added {len(keep)}', file=sys.stderr)

if __name__ == '__main__':
    if sys.argv[1] == 'emit': emit(*sys.argv[2:5])
    else: main()

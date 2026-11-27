#!/usr/bin/env python3
"""A planted sample that HAS copies the reference lacks: hide copies of a planted set from the genome (register 1982).

    hide_copy.py PLANT_DIR GENOME.fa OUT.fa [--seed 1982] [--share 0.25]

The planted reads stay as they are (every copy expressed as planned); the reads are then mapped to OUT.fa, a copy of
GENOME.fa in which the span of each hidden copy (first to last planted exon) is N. Coordinates do not move, so loci.gtf,
units.tsv and the truth tables stay valid. Per family of >= 2 copies, max(1, round(--share x copies)) copies are hidden,
drawn with --seed among copies whose span overlaps no other planted copy; at least one copy of every family stays.
Needs PLANT_DIR/units.tsv and dist.tsv (units.py) and GENOME.fa.fai. Writes OUT.fa, OUT.fa.fai and OUT.fa.hidden.tsv
(family, copy, name, unit, chrom, start, end, nearest = the remaining copy with the fewest edits, edits; 0 edits with
twin = 1 means a copy of the same unit remains: the reference-identical case).
An N block is not how a real assembly loses a copy (a collapse leaves no gap); for a mapper both mean "no sequence here".
--selftest runs the pure parts."""
import sys, os, random, shutil, argparse, collections

def choose(units, share, rng):
    """units [(family, copy, chrom, start, end)] -> hidden copies (sorted)"""
    fam = collections.defaultdict(list)
    for u in units: fam[u[0]].append(u)
    out = []
    for f in sorted(fam):
        us = sorted(fam[f], key=lambda u: u[1])
        if len(us) < 2: continue
        free = [u[1] for u in us if not any(v[1] != u[1] and v[2] == u[2] and v[3] <= u[4] and u[3] <= v[4] for v in units)]
        k = min(max(1, round(share * len(us))), len(us) - 1, len(free))
        out += rng.sample(free, k)
    return sorted(out)

def byte_runs(off, lb, lw, s, e):
    """0-based half-open [s, e) of a sequence -> [(file offset, length)] one run per FASTA line"""
    o = []
    while s < e:
        n = min(e - s, lb - s % lb); o.append((off + s // lb * lw + s % lb, n)); s += n
    return o

def nearest(copy, hidden, unit, dist):
    """-> (remaining copy, edits, twin) with the fewest edits to `copy`; a remaining copy of its unit counts as 0 edits"""
    tw = sorted(c for c, u in unit.items() if u == unit[copy] and c != copy and c not in hidden)
    if tw: return tw[0], 0, 1
    c = sorted((d, b) for (a, b), d in dist.items() if a == copy and b not in hidden)
    return (c[0][1], c[0][0], 0) if c else ('.', -1, 0)

def main():
    if sys.argv[1:] == ['--selftest']:
        u = [('F', 'a', 'c', 1, 10), ('F', 'b', 'c', 5, 20), ('F', 'c', 'c', 30, 40), ('F', 'd', 'c', 50, 60), ('G', 'x', 'c', 100, 110)]
        h = choose(u, 0.25, random.Random(1)); assert len(h) == 1 and h[0] in ('c', 'd'), h
        assert choose([('F', 'a', 'c', 1, 10), ('F', 'b', 'c', 30, 40)], 0.9, random.Random(1)) in (['a'], ['b'])
        assert byte_runs(100, 60, 61, 58, 63) == [(158, 2), (161, 3)]
        assert nearest('a', {'a'}, {'a': 'u1', 'b': 'u1', 'c': 'u2'}, {}) == ('b', 0, 1)
        assert nearest('a', {'a', 'b'}, {'a': 'u1', 'b': 'u1', 'c': 'u2'}, {('a', 'c'): 7, ('a', 'b'): 0}) == ('c', 7, 0)
        print('hide_copy selftest ok'); return
    ap = argparse.ArgumentParser(); ap.add_argument('plant'); ap.add_argument('genome'); ap.add_argument('out')
    ap.add_argument('--seed', type=int, default=1982); ap.add_argument('--share', type=float, default=0.25); a = ap.parse_args()
    U = [l.rstrip('\n').split('\t') for l in open(f'{a.plant}/units.tsv')][1:]
    hidden = set(choose([(f[0], f[1], f[5], int(f[6]), int(f[7])) for f in U], a.share, random.Random(a.seed)))
    unit = {f[1]: f[3] for f in U}; dist = {}
    for l in list(open(f'{a.plant}/dist.tsv'))[1:]:
        f = l.rstrip('\n').split('\t'); dist[(f[1], f[2])] = dist[(f[2], f[1])] = int(f[3])
    fai = {}
    for l in open(a.genome + '.fai'):
        f = l.rstrip('\n').split('\t'); fai[f[0]] = tuple(map(int, f[1:5]))
    tmp = f'{a.out}.tmp.{os.getpid()}'; shutil.copyfile(a.genome, tmp)
    with open(tmp, 'r+b') as fh, open(a.out + '.hidden.tsv', 'w') as w:
        w.write('family\tcopy\tname\tunit\tchrom\tstart\tend\tnearest\tedits\ttwin\n')
        for f in U:
            if f[1] not in hidden: continue
            n, off, lb, lw = fai[f[5]]
            for o, k in byte_runs(off, lb, lw, int(f[6]) - 1, min(n, int(f[7]))): fh.seek(o); fh.write(b'N' * k)
            w.write('\t'.join([f[0], f[1], f[2], f[3], f[5], f[6], f[7], *map(str, nearest(f[1], hidden, unit, dist))]) + '\n')
    os.replace(tmp, a.out); shutil.copyfile(a.genome + '.fai', a.out + '.fai')
    print(f'hidden {len(hidden)} of {len(U)} copies; reference written to {a.out}')

if __name__ == '__main__': main()

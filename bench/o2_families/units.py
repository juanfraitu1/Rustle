#!/usr/bin/env python3
"""Units of indistinguishable copies and the divergence between copies of a planted set (register 1906).

    units.py PLANT_DIR GENOME.fa [THREADS]

Writes PLANT_DIR/units.tsv (family, copy, name, unit, unit_size, chrom, start, end, strand, isoforms): copies whose
planted isoforms are identical in sequence share a unit; PLANT_DIR/copies.fa (first isoform of each copy, on its strand);
PLANT_DIR/dist.tsv (family, copy_a, copy_b, edits): minimap2 asm20 of the copies of a family against each other,
edits = NM + bases of either sequence outside the alignment, best alignment per pair; a pair with no alignment is absent.
PLANT_DIR/families.tsv: per family copies, units, the smallest and the median pairwise edits, chromosomes, span, strands.
Needs GENOME.fa.fai and minimap2."""
import sys, os, re, subprocess, collections, statistics, hashlib
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plant
d, fa = sys.argv[1:3]; thr = sys.argv[3] if len(sys.argv) > 3 else '8'
FA = plant.Fai(fa); tx = collections.defaultdict(list); meta = {}
for l in open(f'{d}/loci.gtf'):
    f = l.rstrip('\n').split('\t')
    if f[2] != 'exon': continue
    t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); tx[t].append((int(f[3]), int(f[4]))); meta[t] = (f[0], f[6], re.search(r'gene_id "([^"]+)"', f[8]).group(1))
rows = [l.rstrip('\n').split('\t') for l in open(f'{d}/plan.tsv')][1:]
seq = {}
for t, ex in tx.items():
    s = ''.join(FA.fetch(meta[t][0], a - 1, b) for a, b in sorted(ex)).upper(); seq[t] = s[::-1].translate(plant.COMP) if meta[t][1] == '-' else s
copy = collections.OrderedDict()
for r in rows: copy.setdefault(r[1], dict(fam=r[0], name=r[2], tx=[]))['tx'].append(r[3])
key = {g: hashlib.md5('|'.join(sorted(seq[t] for t in c['tx'])).encode()).hexdigest() for g, c in copy.items()}
uid = {}; size = collections.Counter(key.values())
for g in copy: uid.setdefault(key[g], f"{copy[g]['fam']}.u{len({v for v in uid.values() if v.startswith(copy[g]['fam'] + '.u')}) + 1}")
with open(f'{d}/units.tsv', 'w') as w, open(f'{d}/copies.fa', 'w') as q:
    w.write('family\tcopy\tname\tunit\tunit_size\tchrom\tstart\tend\tstrand\tisoforms\n')
    for g, c in copy.items():
        ex = [e for t in c['tx'] for e in tx[t]]; t0 = c['tx'][0]
        w.write(f"{c['fam']}\t{g}\t{c['name']}\t{uid[key[g]]}\t{size[key[g]]}\t{meta[t0][0]}\t{min(a for a, _ in ex)}\t{max(b for _, b in ex)}\t{meta[t0][1]}\t{len(c['tx'])}\n")
        q.write(f">{g}\n{seq[t0]}\n")
L = {g: len(seq[c['tx'][0]]) for g, c in copy.items()}
paf = subprocess.run(['minimap2', '-c', '--eqx', '-x', 'asm20', '-X', '-N', '200', '-p', '0.1', '-t', thr, f'{d}/copies.fa', f'{d}/copies.fa'], capture_output=True, text=True).stdout
best = {}
for l in paf.splitlines():
    f = l.split('\t'); a, b = f[0], f[5]
    if a == b or copy[a]['fam'] != copy[b]['fam'] or f[4] != '+': continue
    nm = int(next(x[5:] for x in f[12:] if x.startswith('NM:i:')))
    e = nm + (int(f[1]) - (int(f[3]) - int(f[2]))) + (int(f[6]) - (int(f[8]) - int(f[7])))
    k = (min(a, b), max(a, b))
    if e < best.get(k, 1 << 60): best[k] = e
for (a, b) in list(best):
    if key[a] == key[b]: best[(a, b)] = 0
with open(f'{d}/dist.tsv', 'w') as w:
    w.write('family\tcopy_a\tcopy_b\tedits\n')
    for (a, b), e in sorted(best.items()): w.write(f"{copy[a]['fam']}\t{a}\t{b}\t{e}\n")
U = [l.rstrip('\n').split('\t') for l in open(f'{d}/units.tsv')][1:]
with open(f'{d}/families.tsv', 'w') as w:
    w.write('family\tcopies\tunits\tlargest_unit\tpairs_aligned\tmin_edits\tmedian_edits\tchroms\tspan_kb\tstrands\tneighbours_within_200kb\n')
    for f in sorted({u[0] for u in U}):
        us = [u for u in U if u[0] == f]; es = [e for (a, b), e in best.items() if copy[a]['fam'] == f]
        ch = collections.Counter(u[5] for u in us); c0 = ch.most_common(1)[0][0]; on = sorted((int(u[6]), int(u[7])) for u in us if u[5] == c0)
        nb = sum(1 for (a, b), (c, e) in zip(on, on[1:]) if c - b <= 200000)
        w.write(f"{f}\t{len(us)}\t{len({u[3] for u in us})}\t{max(int(u[4]) for u in us)}\t{len(es)}\t{min(es) if es else 'NA'}\t{statistics.median(es) if es else 'NA'}\t{','.join(f'{k}:{v}' for k, v in ch.items())}\t{(on[-1][1] - on[0][0]) / 1000:.0f}\t{''.join(sorted({u[8] for u in us}))}\t{nb}\n")
print(open(f'{d}/families.tsv').read(), end='')

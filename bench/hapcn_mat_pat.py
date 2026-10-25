"""Copies per family in the gorilla reference vs Jim's (KB3781) maternal and paternal haplotypes.
Input: <asm>.paf = GGO_gwcat.copies.fa (spliced transcripts) mapped with minimap2 -c -x splice:hq -uf -N 50 -p 0.5.
A locus = merged hits (query coverage >= 0.8, identity = matches/block >= T) of any of the family's copies."""
import collections, sys
D = sys.argv[1] if len(sys.argv) > 1 else '.'
SEX = {'ref': {'NC_073247.2': 'X', 'NC_073248.2': 'Y'}}
def sex(asm, t):
    if asm == 'ref': return SEX['ref'].get(t)
    return 'X' if t.startswith('chrX') else 'Y' if t.startswith('chrY') else None
fam_of = {}; n_copies = collections.Counter()
for l in open('/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.tsv'):
    if l.startswith('family_id'): continue
    f = l.split('\t'); n_copies[f[0]] += 1
hits = {a: [] for a in ('ref', 'mat', 'pat')}
for a in hits:
    for l in open(f'{D}/{a}.paf'):
        f = l.split('\t'); q = f[0]; fam = q.split('|')[0]
        qlen, qs, qe = int(f[1]), int(f[2]), int(f[3]); t, ts, te = f[5], int(f[7]), int(f[8])
        hits[a].append((fam, q, (qe - qs) / qlen, int(f[9]) / int(f[10]), t, ts, te))
def loci(a, fam, T, auto_only):
    iv = collections.defaultdict(list)
    for fm, q, cov, idn, t, ts, te in hits[a]:
        if fm == fam and cov >= 0.8 and idn >= T and not (auto_only and sex(a, t)):
            iv[t].append((ts, te))
    n = 0
    for t, xs in iv.items():
        xs.sort(); end = -1
        for s, e in xs:
            if s > end + 1000: n += 1
            end = max(end, e)
    return n
fams = sorted(n_copies, key=lambda x: int(x[5:]))
sexfam = {fm for a in hits for fm, q, cov, idn, t, *_ in hits[a] if sex(a, t) and cov >= 0.8 and idn >= 0.9}
print(f'families {len(fams)}, copies {sum(n_copies.values())}; families with a hit on X/Y: {len(sexfam)}')
out = open(f'{D}/hapcn.tsv', 'w'); print('family\tT\tcatalog\tref\tmat\tpat\tref_auto\tmat_auto\tpat_auto', file=out)
for T in (0.90, 0.95, 0.98):
    c = collections.Counter(); diffs = []
    for fm in fams:
        r, m, p = (loci(a, fm, T, False) for a in ('ref', 'mat', 'pat'))
        ra, ma, pa = (loci(a, fm, T, True) for a in ('ref', 'mat', 'pat'))
        print(fm, T, n_copies[fm], r, m, p, ra, ma, pa, sep='\t', file=out)
        if ma == pa == ra: c['autosomal: ref = mat = pat'] += 1
        elif ma == pa: c['autosomal: mat = pat, ref differs'] += 1
        else:
            c['autosomal: mat != pat'] += 1; diffs.append((fm, ra, ma, pa))
        if (r, m, p) != (ra, ma, pa): c['has X/Y loci'] += 1
    print(f'\nT={T}:', dict(c))
    print('  mat != pat (ref, mat, pat; autosomes):', ' '.join(f'{f}({a},{b},{d})' for f, a, b, d in diffs[:40]), '...' if len(diffs) > 40 else '')
# per copy: best identity anywhere in each haplotype (query coverage >= 0.8)
best = {a: collections.defaultdict(float) for a in hits}
for a in hits:
    for fm, q, cov, idn, *_ in hits[a]:
        if cov >= 0.8: best[a][q] = max(best[a][q], idn)
qs = [l[1:].split()[0] for l in open('/mnt/linuxdisk/home/juanfraitu/winloci_data/GGO_gwcat.copies.fa') if l.startswith('>')]
def cls(x): return 'identical' if x >= 0.99995 else '>=99.5%' if x >= 0.995 else '99-99.5%' if x >= 0.99 else '90-99%' if x >= 0.9 else '<90% / none'
print('\nper copy, best match anywhere (query coverage >= 0.8):')
for a in ('ref', 'mat', 'pat'):
    print(f'  {a}:', dict(collections.Counter(cls(best[a][q]) for q in qs)))
both = collections.Counter(cls(max(best['mat'][q], best['pat'][q])) for q in qs)
print('  best of mat/pat:', dict(both))
print('  absent from BOTH haplotypes at >=90%:', [q for q in qs if max(best['mat'][q], best['pat'][q]) < 0.9][:20])

#!/usr/bin/env python3
"""PREREG r1538: PSV pile at absorbing loci (docs/PREREG_psv_pile_2026-11-05.md).
  psv_pile.py read_origin.tsv NT.bam GENOME.fa TRUTH2.tsv OUT   -> OUT.loci.tsv OUT.sites.tsv, summary on stdout
Loci by read landing (primary records clustered at <= 10 kb gaps; the copy's locus = the biggest cluster, scorable when it holds
>= 50% of the copy's reads and >= 10 reads). Blind calling: bcftools mpileup over all reads of nt.bam at the locus (SNV records
only); site = DP >= 20, second allele >= 5 carriers at fraction >= 0.10. S = sites per kb at DP >= 20; H = share of reads
covering >= 2 sites that carry the minor allele at >= 80% of them. Purity (labelled) = origin shares of the minor-allele carriers."""
import collections
import os
import subprocess
import sys

import pysam

BCF = os.path.expanduser('~/.pixi/detached-envs/bcftools-13167962409839689241/envs/default/bin/bcftools')
origin_tsv, nt, genome, truth2, out = sys.argv[1:6]
org = {l.split('\t')[0]: l.rstrip('\n').split('\t')[1:] for l in open(origin_tsv) if not l.startswith('read\t')}
reason = {r.split('\t')[0]: r.rstrip('\n').split('\t')[8] for r in open(truth2) if not r.startswith('gene\t')}
bam = pysam.AlignmentFile(nt)
land = collections.defaultdict(list)
for a in bam:
    if a.is_secondary or a.is_supplementary or a.is_unmapped or a.query_name not in org:
        continue
    land[org[a.query_name][0]].append((a.reference_name, a.reference_start, a.reference_end))

loci = {}  # copy -> (chrom, s, e, n_in, n_all)
for g, xs in land.items():
    xs.sort()
    clusters, cur = [], None
    for c, s, e in xs:
        if cur and cur[0] == c and s <= cur[2] + 10000:
            cur[2] = max(cur[2], e)
            cur[3] += 1
        else:
            cur = [c, s, e, 1]
            clusters.append(cur)
    c, s, e, n = max(clusters, key=lambda x: x[3])
    loci[g] = (c, s, e, n, len(xs))
meta = {g: (o[1], o[3]) for o in org.values() for g in [o[0]]}
scor = {g: l for g, l in loci.items() if l[3] >= 10 and l[3] >= 0.5 * l[4]}
absorb = {g: l for g, l in scor.items() if meta[g][1] == 'missing'}
null = {g: l for g, l in scor.items() if meta[g][1] == 'present'
        and not any(l[0] == m[0] and l[1] < m[2] and l[2] > m[1] for m in absorb.values())}


def call(c, s, e):
    """bcftools mpileup parsed directly: SNV sites and the DP >= 20 position count."""
    txt = subprocess.run([BCF, 'mpileup', '-f', genome, '-r', f'{c}:{s + 1}-{e}', '-Q', '10', '-q', '0', '--ff',
                          'UNMAP,SECONDARY,SUPPLEMENTARY', '-a', 'AD,DP', '-d', '100000', nt],
                         capture_output=True, text=True, check=True).stdout
    cov, sites = 0, []
    for l in txt.splitlines():
        if l.startswith('#'):
            continue
        f = l.split('\t')
        if 'INDEL' in f[7]:
            continue
        info = dict(x.split('=') for x in f[7].split(';') if '=' in x)
        dp = int(info['DP'])
        if dp < 20:
            continue
        cov += 1
        alleles = [f[3]] + f[4].split(',')
        ad = list(map(int, f[9].split(':')[-1].split(',')))
        pairs = sorted(zip(ad, alleles), reverse=True)
        pairs = [(n, al) for n, al in pairs if al != '<*>']
        if len(pairs) >= 2 and pairs[1][0] >= 5 and pairs[1][0] / dp >= 0.10:
            sites.append((int(f[1]) - 1, pairs[0][1], pairs[1][1]))
    return cov, sites


rows, siterows = [], []
for g, (c, s, e, n_in, n_all) in sorted(scor.items(), key=lambda x: (meta[x[0]][1], -x[1][3])):
    kind = 'absorbing' if g in absorb else ('null' if g in null else 'present_overlapping')
    cov, sites = call(c, s, e)
    S = len(sites) / cov * 1000 if cov else float('nan')
    per_read = collections.defaultdict(list)  # read -> [carries_minor]
    carriers = collections.Counter()
    for pos, major, minor in sites:
        for col in bam.pileup(c, pos, pos + 1, truncate=True, min_base_quality=10, flag_filter=0x904):
            for p in col.pileups:
                if p.is_del or p.is_refskip or p.query_position is None:
                    continue
                b = p.alignment.query_sequence[p.query_position]
                if b == minor:
                    per_read[p.alignment.query_name].append(1)
                    o = org.get(p.alignment.query_name)
                    carriers['missing_own' if o and o[0] == g else ('ambiguous' if o is None else meta[o[0]][1])] += 1
                elif b == major:
                    per_read[p.alignment.query_name].append(0)
        siterows.append(f'{g}\t{kind}\t{c}\t{pos + 1}\t{major}\t{minor}\n')
    multi = [v for v in per_read.values() if len(v) >= 2]
    H = sum(sum(v) / len(v) >= 0.8 for v in multi) / len(multi) if multi else 0.0
    known = sum(v for k, v in carriers.items() if k != 'ambiguous')
    pur = carriers['missing_own'] / known if known else float('nan')
    rows.append((g, meta[g][0], meta[g][1], reason.get(g, ''), kind, c, s, e, n_in, n_all, cov, len(sites), S, H, len(multi), pur,
                 dict(carriers)))
with open(out + '.loci.tsv', 'w') as w:
    w.write('gene\tname\tstatus\treason\tkind\tchrom\tstart\tend\treads_in\treads_all\tcov20\tsites\tS\tH\treads_ge2\tpurity\tcarriers\n')
    for r in rows:
        w.write('\t'.join(f'{x:.4f}' if isinstance(x, float) else str(x) for x in r) + '\n')
with open(out + '.sites.tsv', 'w') as w:
    w.write('gene\tkind\tchrom\tpos\tmajor\tminor\n')
    w.writelines(siterows)
# thresholds per ploidy stratum from the null loci
for stratum, pick in (('autosomal', lambda c: c not in ('chrX', 'chrY')), ('XY', lambda c: c in ('chrX', 'chrY'))):
    nl = [r for r in rows if r[4] == 'null' and pick(r[5])]
    if not nl:
        print(f'{stratum}: no null loci')
        continue
    q = (lambda xs: sorted(xs)[int(0.95 * (len(xs) - 1))]) if len(nl) >= 5 else max
    tS, tH = q([r[12] for r in nl]), q([r[13] for r in nl])
    ab = [r for r in rows if r[4] == 'absorbing' and pick(r[5])]
    hm = [r for r in ab if r[12] > tS or r[13] > tH]
    nh = [r for r in nl if r[12] > tS or r[13] > tH]
    print(f'{stratum}: null loci {len(nl)} (thresholds S {tS:.2f}/kb H {tH:.3f}; null hallmarked {len(nh)}); '
          f'missing copies {len(ab)} hallmarked {len(hm)}; collapsed {sum(r[3] == "collapsed" for r in hm)}/{sum(r[3] == "collapsed" for r in ab)} '
          f'absent {sum(r[3] == "absent" for r in hm)}/{sum(r[3] == "absent" for r in ab)}; '
          f'unique absorbing loci {len({(r[5], r[6], r[7]) for r in ab})} hallmarked {len({(r[5], r[6], r[7]) for r in hm})}')
    for r in hm:
        print(f'  HM {r[1]} {r[5]}:{r[6]}-{r[7]} S {r[12]:.2f} H {r[13]:.3f} purity {r[15]:.2f} {r[16]}')
print(f'scorable: missing {len(absorb)} of {sum(1 for g in loci if meta[g][1] == "missing")}, null {len(null)}; '
      f'missing copies with reads < 10% of pile (floor): '
      f'{sum(1 for l in absorb.values() if l[3] < 0.1 * sum(1 for a in bam.fetch(l[0], l[1], l[2]) if not a.flag & 0x904))}')

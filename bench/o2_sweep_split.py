#!/usr/bin/env python3
"""Split an `mcl_families --emit-units` catalog into per-(family, contig) O2 input directories.

usage: o2_sweep_split.py <prefix> <out_dir> [family_list]
  <prefix>.units.tsv / .units.fa  ->  <out_dir>/fam_<family>_<contig6>/{copies.tsv,copies.fa,regions,forecast.tsv}
A sweep family = the units of one catalog family on one contig (O2 assumes locality), copies renumbered
0..k-1 in catalog order; `family_list` (one id per line) restricts the output to those sweep families so
two catalogs can be compared PAIRED. Region = contig:min(start)-5000 .. max(end)+5000.
"""
import sys, os, csv, collections
prefix, out = sys.argv[1], sys.argv[2]
only = set(open(sys.argv[3]).read().split()) if len(sys.argv) > 3 else None
rows = list(csv.DictReader(open(prefix + '.units.tsv'), delimiter='\t'))
seqs = {}; name = None
for line in open(prefix + '.units.fa'):
    if line.startswith('>'): name = line[1:].split('|')[0] + '|' + line[1:].split('|')[1]; seqs[name] = []
    else: seqs[name].append(line.strip())
groups = collections.OrderedDict()
for r in rows:
    # contig tag: `NC_073242.2` -> `073242` (the historical form), any other name kept as-is (chrY)
    c = r['chrom']
    fid = f"{r['family_id']}_{c.split('_')[1].split('.')[0] if '_' in c else c}"
    groups.setdefault(fid, []).append(r)
# §6ft partners: for every family unit, the nearest catalog unit of ANOTHER family on each side (the L2 clipping
# neighbours). They enter copies.tsv with member_status = partner: O2 aligns molecules to them so a read-through
# tail is explained, but never assigns to them.
by_chrom = collections.defaultdict(list)
for r in rows: by_chrom[r['chrom']].append(r)
for c in by_chrom: by_chrom[c].sort(key=lambda r: int(r['start']))
def partners_of(rs):
    out = {}
    for r in rs:
        lst = by_chrom[r['chrom']]; s0, e0 = int(r['start']), int(r['end'])
        left = [x for x in lst if x['family_id'] != r['family_id'] and int(x['end']) <= s0]
        right = [x for x in lst if x['family_id'] != r['family_id'] and int(x['start']) >= e0]
        for x in ([max(left, key=lambda x: int(x['end']))] if left else []) + ([min(right, key=lambda x: int(x['start']))] if right else []):
            out[(x['family_id'], x['copy_idx'])] = x
    return [x for k, x in out.items() if not any(x['family_id'] == r['family_id'] for r in rs)]
n = 0
for fid, rs in groups.items():
    if len(rs) < 2 or (only is not None and fid not in only): continue
    partners = partners_of(rs) if '--no-partners' not in sys.argv else []
    d = os.path.join(out, 'fam_' + fid); os.makedirs(d, exist_ok=True)
    with open(f'{d}/copies.tsv', 'w') as ct, open(f'{d}/copies.fa', 'w') as fa, open(f'{d}/forecast.tsv', 'w') as fc:
        extra = [c for c in ('member_status', 'locus_start', 'locus_end') if c in rows[0]]  # L1/L2 columns, when the catalog has them
        ct.write('family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tcore_hull' + ''.join('\t' + c for c in extra) + '\n')
        fc.write('copy_idx\tnearest_ident\tsd_depth\tcore_bp\trep_frac\tsource\n')
        for i, r in enumerate(rs + partners):
            rr = dict(r)
            if i >= len(rs): rr['member_status'] = 'partner'
            ct.write('\t'.join([fid, str(i), r['tid'], r['chrom'], r['start'], r['end'], r['n_exon'], r['strand'], r['n_reads'], r['exons'], r['core_hull']] + [rr.get(c, 'NA') for c in extra]) + '\n')
            fc.write('\t'.join([str(i), r['nearest_ident'], r['sd_depth'], r['core_bp'], r['rep_frac'], r['source'] if i < len(rs) else 'partner']) + '\n')
            fa.write(f">{fid}|{i}|{r['chrom']}:{r['start']}-{r['end']}|{r['strand']}|nexon={r['n_exon']}\n")
            fa.write('\n'.join(seqs[r['family_id'] + '|' + r['copy_idx']]) + '\n')
    lo = min(int(r['start']) for r in rs + partners) - 5000; hi = max(int(r['end']) for r in rs + partners) + 5000  # the region holds the partners too
    open(f'{d}/regions', 'w').write(f"{rs[0]['chrom']}:{lo}-{hi}\n"); n += 1
print(f'wrote {n} sweep families to {out}')

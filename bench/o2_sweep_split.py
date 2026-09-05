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
    fid = f"{r['family_id']}_{r['chrom'].split('_')[1].split('.')[0]}"
    groups.setdefault(fid, []).append(r)
n = 0
for fid, rs in groups.items():
    if len(rs) < 2 or (only is not None and fid not in only): continue
    d = os.path.join(out, 'fam_' + fid); os.makedirs(d, exist_ok=True)
    with open(f'{d}/copies.tsv', 'w') as ct, open(f'{d}/copies.fa', 'w') as fa, open(f'{d}/forecast.tsv', 'w') as fc:
        ct.write('family_id\tcopy_idx\ttid\tchrom\tstart\tend\tn_exon\tstrand\tn_reads\texons\tcore_hull\n')
        fc.write('copy_idx\tnearest_ident\tsd_depth\tcore_bp\trep_frac\tsource\n')
        for i, r in enumerate(rs):
            ct.write('\t'.join([fid, str(i), r['tid'], r['chrom'], r['start'], r['end'], r['n_exon'], r['strand'], r['n_reads'], r['exons'], r['core_hull']]) + '\n')
            fc.write('\t'.join([str(i), r['nearest_ident'], r['sd_depth'], r['core_bp'], r['rep_frac'], r['source']]) + '\n')
            fa.write(f">{fid}|{i}|{r['chrom']}:{r['start']}-{r['end']}|{r['strand']}|nexon={r['n_exon']}\n")
            fa.write('\n'.join(seqs[r['family_id'] + '|' + r['copy_idx']]) + '\n')
    lo = min(int(r['start']) for r in rs) - 5000; hi = max(int(r['end']) for r in rs) + 5000
    open(f'{d}/regions', 'w').write(f"{rs[0]['chrom']}:{lo}-{hi}\n"); n += 1
print(f'wrote {n} sweep families to {out}')

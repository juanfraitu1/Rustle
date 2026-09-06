#!/usr/bin/env python3
"""O3 item 4: the `other_family` admission loci as an MCL-cut certificate. Reads of family F whose primary lies in
family G's unit mean F and G exchange reads; for each (F, G) pair: shared reads, loci, and the best identity between
an F unit and the G unit(s) involved (minimap2 asm20 of units.fa sequences).
usage: o3_cut_certificate.py <pass_dir> <units_prefix> <out.tsv> [min_reads=20]"""
import sys, csv, collections, subprocess, tempfile, os
pass_dir, prefix, out = sys.argv[1:4]; min_reads = int(sys.argv[4]) if len(sys.argv) > 4 else 20
L = [r for r in csv.DictReader(open(f'{pass_dir}/loci.tsv'), delimiter='\t') if r['class'] == 'other_family']
seqs = {}; name = None
for l in open(prefix + '.units.fa'):
    if l.startswith('>'): h = l[1:].split('|'); name = h[0] + '|' + h[1]; seqs[name] = []
    else: seqs[name].append(l.strip())
seqs = {k: ''.join(v) for k, v in seqs.items()}
U = {(r['family_id'], r['copy_idx']): r for r in csv.DictReader(open(prefix + '.units.tsv'), delimiter='\t')}
pairs = collections.defaultdict(lambda: {'reads': 0, 'loci': 0, 'g_units': set(), 'nested_reads': 0, 'nested_loci': 0})
fam_units = collections.defaultdict(list)
for (fid, ci), r in U.items(): fam_units[fid].append((r['chrom'], int(r['locus_start']) if r.get('locus_start', 'NA') not in ('NA', '') else int(r['start']), int(r['locus_end']) if r.get('locus_end', 'NA') not in ('NA', '') else int(r['end'])))
for r in L:
    F = r['fam'].replace('fam_', '').split('_')[0]; gu = r['other_family_units'].split(';')[0].split(':'); G = gu[0]
    if G == F: continue
    k = (F, G); pairs[k]['reads'] += int(r['n_reads']); pairs[k]['loci'] += 1; pairs[k]['g_units'].add(gu[1])
    # nested: the shared locus lies inside (or overlaps) one of F's own loci — units of two families over one place
    c, s, e = r['chrom'], int(r['start']), int(r['end'])
    if any(uc == c and us < e and s < ue for uc, us, ue in fam_units[F]): pairs[k]['nested_reads'] += int(r['n_reads']); pairs[k]['nested_loci'] += 1
rows = []
with tempfile.TemporaryDirectory() as tmp:
    for (F, G), v in sorted(pairs.items(), key=lambda kv: -kv[1]['reads']):
        if v['reads'] < min_reads: continue
        fu = [k for k in seqs if k.startswith(F + '|')]; gu = [f'{G}|{i}' for i in v['g_units'] if f'{G}|{i}' in seqs]
        best = 0.0; blen = 0
        if fu and gu:
            open(f'{tmp}/f.fa', 'w').write(''.join(f'>{k}\n{seqs[k]}\n' for k in fu)); open(f'{tmp}/g.fa', 'w').write(''.join(f'>{k}\n{seqs[k]}\n' for k in gu))
            paf = subprocess.run(['minimap2', '-x', 'asm20', '-c', '-N', '5', '-t', '4', f'{tmp}/g.fa', f'{tmp}/f.fa'], capture_output=True, text=True).stdout
            for l in paf.splitlines():
                f = l.split('\t'); ident = int(f[9]) / max(1, int(f[10]))
                if int(f[10]) >= 300 and (ident > best or (ident == best and int(f[10]) > blen)): best, blen = ident, int(f[10])
        nF = sum(1 for k in U if k[0] == F); nG = sum(1 for k in U if k[0] == G)
        rows.append([F, G, v['reads'], v['loci'], v['nested_reads'], v['nested_loci'], ','.join(sorted(v['g_units'])), nF, nG, f'{best:.4f}' if blen else 'NA', blen, 'nested' if v['nested_reads'] >= 0.5 * v['reads'] else ('cut' if blen and best >= 0.70 else 'disjoint_no_homology')])
with open(out, 'w') as o:
    o.write('F\tG\tshared_reads\tloci\tnested_reads\tnested_loci\tG_units\tF_units\tG_units_total\tbest_unit_identity\taln_bp\tclass\n')
    for r in rows: o.write('\t'.join(map(str, r)) + '\n')
cc = collections.Counter(r[11] for r in rows); rr = collections.Counter(); [rr.__setitem__(r[11], rr[r[11]] + r[2]) for r in rows]
print(f'(F,G) pairs with >= {min_reads} shared reads: {len(rows)} | classes {dict(cc)} | reads by class {dict(rr)} | unit identity >= 0.90: {sum(1 for r in rows if r[9] != "NA" and float(r[9]) >= 0.90)}')

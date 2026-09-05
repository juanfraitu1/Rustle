#!/usr/bin/env python3
"""§6fb: the reads O2's origin certificate rejects (`origin_rejected=1`), per family: how many unit reads, where
their primary placement lies (inside a unit of the family = the unit lacks their content; outside every unit =
a candidate unannotated origin), and the outside placements clustered into candidate loci (>= min_reads
primaries within `gap` bp of each other), annotated with the RefSeq gene there, if any.
usage: o2_unexplained_reads.py <sweep_dir> <bam> <gff> <out_prefix> [min_reads=3] [gap=5000]"""
import sys, os, glob, csv, collections, pysam
sweep, bam_p, gff, out = sys.argv[1:5]; min_reads = int(sys.argv[5]) if len(sys.argv) > 5 else 3; gap = int(sys.argv[6]) if len(sys.argv) > 6 else 5000
BAM = pysam.AlignmentFile(bam_p)
genes = collections.defaultdict(list)
for l in open(gff):
    if l.startswith('#'): continue
    f = l.rstrip('\n').split('\t')
    if len(f) < 9 or f[2] not in ('gene', 'pseudogene'): continue
    a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
    genes[f[0]].append((int(f[3]) - 1, int(f[4]), a.get('gene', '?'), a.get('gene_biotype', '?')))
def gene_at(c, s, e):
    return [g for g in genes.get(c, []) if g[0] < e and s < g[1]]
fam_rows = []; loci_rows = []
for d in sorted(glob.glob(f'{sweep}/fam_*')):
    asg = f'{d}/A.assignments.tsv'
    if not os.path.exists(asg): continue
    fam = os.path.basename(d)
    cp = {r['copy_idx']: (r['chrom'], int(r['start']), int(r['end'])) for r in csv.DictReader(open(f'{d}/copies.tsv'), delimiter='\t')}
    A = {r['read_name']: r for r in csv.DictReader(open(asg), delimiter='\t')}
    rej = {n for n, r in A.items() if r.get('origin_rejected') == '1'}
    # primary placements of the family's unit reads
    unit_reads = {}; prim = {}
    chrom = next(iter(cp.values()))[0]; lo = min(v[1] for v in cp.values()) - 200000; hi = max(v[2] for v in cp.values()) + 200000
    for a in BAM.fetch(chrom, max(0, lo), hi):
        if a.flag & 2308 or a.query_name not in A: continue
        blocks = a.get_blocks(); inside = None
        for i, (c, s, e) in cp.items():
            if any(b1 > s and b0 < e for b0, b1 in blocks): inside = i; break
        prim[a.query_name] = (a.reference_start, a.reference_end, inside, a.mapping_quality)
        if inside is not None: unit_reads[a.query_name] = inside
    n_unit = len(unit_reads); rej_unit = [n for n in unit_reads if n in rej]
    rej_inside = collections.Counter(unit_reads[n] for n in rej_unit)
    # reads the certificate rejected whose primary lies OUTSIDE every unit (the read touched a unit via another record)
    outside = sorted((prim[n][0], prim[n][1], n) for n in rej if n in prim and prim[n][2] is None)
    loci = []
    for s, e, n in outside:
        if loci and s - loci[-1][1] <= gap: loci[-1][1] = max(loci[-1][1], e); loci[-1][2].append(n)
        else: loci.append([s, e, [n]])
    loci = [l for l in loci if len(l[2]) >= min_reads]
    fam_rows.append((fam, len(cp), n_unit, len(rej_unit), sum(rej_inside.values()), len(outside), len(loci)))
    for s, e, ns in loci:
        g = gene_at(chrom, s, e)
        loci_rows.append((fam, chrom, s, e, len(ns), ';'.join(f'{x[2]}[{x[3]}]' for x in g[:3]) or 'NO_GENE', sum(1 for n in ns if prim[n][3] >= 60)))
with open(out + '.families.tsv', 'w') as o:
    o.write('fam\tn_copies\tunit_reads\torigin_rejected_unit_reads\trejected_with_primary_inside_a_unit\trejected_primary_outside_units\tcandidate_loci\n')
    for r in fam_rows: o.write('\t'.join(map(str, r)) + '\n')
with open(out + '.loci.tsv', 'w') as o:
    o.write('fam\tchrom\tstart\tend\tn_reads\tgenes\tn_mapq60\n')
    for r in loci_rows: o.write('\t'.join(map(str, r)) + '\n')
T = [sum(r[i] for r in fam_rows) for i in range(2, 7)]
print(f'families {len(fam_rows)} | unit reads {T[0]} | origin-rejected unit reads {T[1]} ({100*T[1]/max(1,T[0]):.1f}%) — primary inside a unit {T[2]}, outside every unit {T[3]} | candidate loci (>= {min_reads} reads) {T[4]}')

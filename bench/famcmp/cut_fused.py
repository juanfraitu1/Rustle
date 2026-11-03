#!/usr/bin/env python3
"""r1523 step B: cut models that span >= 2 truth-family members (each expressed, T >= 2) at the midpoint between the members.
usage: cut_fused.py IN.gtf GENES.tsv OUT.gtf   (GENES.tsv = family_panel.py's *.genes.tsv: every row is an expressed member)
       cut_fused.py IN.gtf --copies COPIES.tsv OUT.gtf   (r1524: catalog families, x.fam.copies.tsv; a copy counts when n_reads >= 2)
Exons left of the cut form one model, exons right of it another; a side with no exon is dropped. Everything else is copied."""
import sys, os, re, csv, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
try:
    import famtruth as F
except ImportError:  # --copies mode needs no truth
    F = None
members = collections.defaultdict(list)  # (chrom, strand) -> [(lo, hi, exons, family, gene)]
if sys.argv[2] == '--copies':
    gtf, _, copies_tsv, out = sys.argv[1:5]
    for r in csv.DictReader(open(copies_tsv), delimiter='\t'):
        if int(r['n_reads']) < 2: continue
        lo, hi = int(r['locus_start']), int(r['locus_end'])
        members[(r['chrom'], r['strand'])].append((lo, hi, [(lo, hi)], r['family_id'], r['copy_idx']))
else:
    gtf, genes_tsv, out = sys.argv[1:4]
    G = F.load_cat()
    for r in csv.DictReader(open(genes_tsv), delimiter='\t'):
        v = G[r['gene']]; members[(v['chrom'], v['strand'])].append((v['exons'][0][0], v['exons'][-1][1], v['exons'], r['family'], r['gene']))
tx = collections.OrderedDict(); lines = collections.defaultdict(list)
for l in open(gtf):
    f = l.rstrip('\n').split('\t')
    if len(f) < 9: continue
    t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1); lines[t].append(f)
    if f[2] == 'exon': tx.setdefault(t, (f[0], f[6]))
n_cut = 0
with open(out, 'w') as w:
    for t, fs in lines.items():
        ex = sorted((int(f[3]) - 1, int(f[4])) for f in fs if f[2] == 'exon')
        if not ex: continue
        c, s = tx[t]; lo, hi = ex[0][0], ex[-1][1]
        hit = {}
        for mlo, mhi, mex, fam, g in members.get((c, s), []):
            if mlo < hi and lo < mhi and any(a < y and x < b for a, b in ex for x, y in mex): hit.setdefault(fam, []).append((mlo, mhi, g))
        fams = [sorted(v) for v in hit.values() if len(v) >= 2]
        if not fams:
            for f in fs: w.write('\t'.join(f) + '\n')
            continue
        # cut points: midpoints between consecutive members of the same family (members may overlap: then no cut there)
        cuts = sorted({(a[1] + b[0]) // 2 for v in fams for a, b in zip(v, v[1:]) if a[1] < b[0]})
        if not cuts:
            for f in fs: w.write('\t'.join(f) + '\n')
            continue
        n_cut += 1; tf = next(f for f in fs if f[2] == 'transcript'); attr = tf[8]
        parts = [[]]; ci = 0
        for a, b in ex:
            while ci < len(cuts) and a >= cuts[ci]: parts.append([]); ci += 1
            while ci < len(cuts) and b > cuts[ci]:  # an exon spanning the cut is split, not truncated
                parts[-1].append((a, cuts[ci])); a = cuts[ci]; parts.append([]); ci += 1
            parts[-1].append((a, b))
        for i, p in enumerate(p for p in parts if p):
            a2 = attr.replace(f'transcript_id "{t}"', f'transcript_id "{t}.cut{i}"')
            w.write('\t'.join(tf[:3] + [str(p[0][0] + 1), str(p[-1][1])] + tf[5:8] + [a2]) + '\n')
            for a, b in p: w.write('\t'.join(tf[:2] + ['exon', str(a + 1), str(b)] + tf[5:8] + [a2]) + '\n')
print(f'models cut: {n_cut}')

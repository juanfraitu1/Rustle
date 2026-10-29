#!/usr/bin/env python3
"""MINIPROT member recovery (dev, 2026-10-29, gorilla; user "yes lets do that"): take the longest CDS of every protein-coding gene that a
family catalogue places in a family, map it to the genome with miniprot (frameshift-tolerant), and append every placement that lands on
an annotated gene not yet in that family as a new copy of the family. Classification is left to pairscore.py (coding pairs by Compara,
non-coding pairs by the sequence judge).

  miniprot_members.py gorilla GENOME.fa GENOME.mpi ARM.fam.copies.tsv BASE.families.gtf OUTDIR [--min-id 0.85] [--min-cov 0.5]
                      [--mode all|expr|expr_nc] [--threads 32]
  writes OUTDIR/prot.fa (protjudge prep), OUTDIR/mp.gff, OUTDIR/<arm>_mp<id>_<mode>.fam.copies.tsv
  mode: all = every target gene; expr = only targets with a de novo transcript in BASE.families.gtf; expr_nc = expressed and not
  protein-coding (pseudogene / lncRNA only; coding membership stays with the family graph).
"""
import sys, os, collections, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import famtruth as F
a = sys.argv[1:]
species, genome, mpi, arm, gtf, out = a[:6]
minid = float(a[a.index('--min-id') + 1]) if '--min-id' in a else 0.85
mincov = float(a[a.index('--min-cov') + 1]) if '--min-cov' in a else 0.5
mode = a[a.index('--mode') + 1] if '--mode' in a else 'expr_nc'
thr = a[a.index('--threads') + 1] if '--threads' in a else '32'
assert species == 'gorilla'
os.makedirs(out, exist_ok=True)
G = F.load_ggo(F.load_cat()); idx = F.exon_index(G)
fams, _, _ = F.load_catalog(arm, G, idx)
fam_of = {g: f for f, v in fams.items() for g in v if G[g]['biotype'] == 'protein_coding'}
open(f'{out}/genes.txt', 'w').write('\n'.join(sorted(fam_of)) + '\n')
if not os.path.exists(f'{out}/prot.fa'):
    subprocess.run([sys.executable, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'protjudge.py'), 'prep', species, genome, out, f'{out}/genes.txt'], check=True)
if not os.path.exists(f'{out}/mp.gff'):
    with open(f'{out}/mp.gff', 'w') as w:
        subprocess.run(['miniprot', '-t', thr, '--gff', '-N', '50', '--outs', '0.5', mpi, f'{out}/prot.fa'], stdout=w, check=True)
qlen = {}; cur = None
for l in open(f'{out}/prot.fa'):
    if l[0] == '>': cur = l[1:].split()[0]; qlen[cur] = 0
    else: qlen[cur] += len(l.strip())
expressed = F.expressed_from_gtf(gtf, idx)
mr = {}; cds = collections.defaultdict(list)
for l in open(f'{out}/mp.gff'):
    if l[0] == '#': continue
    f = l.rstrip('\n').split('\t'); at = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
    if f[2] == 'mRNA':
        t = at['Target'].split(); mr[at['ID']] = (f[0], f[6], t[0], int(t[1]), int(t[2]), float(at['Identity']))
    elif f[2] == 'CDS': cds[at['Parent']].append((int(f[3]) - 1, int(f[4])))
tag = os.path.basename(arm).replace('.fam.copies.tsv', '') + f'_mp{minid}_{mode}'
seen = set(); n = 0; kinds = collections.Counter()
with open(f'{out}/{tag}.fam.copies.tsv', 'w') as w:
    w.write(open(arm).read())
    for mid, (chrom, strand, q, s, e, ident) in mr.items():
        if (e - s + 1) / qlen[q] < mincov or ident < minid or q not in fam_of: continue
        ex = sorted(cds[mid]); tgt = F.map_copy(idx, chrom, strand, ex)
        if tgt is None or tgt == q or tgt in fams[fam_of[q]] or (fam_of[q], tgt) in seen: continue
        if 'expr' in mode and tgt not in expressed: continue
        if 'nc' in mode and G[tgt]['biotype'] == 'protein_coding': continue
        seen.add((fam_of[q], tgt)); n += 1; kinds[G[tgt]['biotype']] += 1
        w.write(f"{fam_of[q]}\tmp{n}\tMP_{q}_{mid}\t{chrom}\t{ex[0][0]+1}\t{ex[-1][1]}\t{len(ex)}\t{strand}\t0\t{','.join(f'{x+1}-{y}' for x, y in ex)}\n")
print(f"{tag}: placements {len(mr)}, members added {n} {dict(kinds)} -> {out}/{tag}.fam.copies.tsv")

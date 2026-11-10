#!/usr/bin/env python3
"""PREREG r1611 (docs/PREREG_gorilla_protein_truth_2026-11-10.md): human Ensembl gene of each annotated gene by best protein hit.

  protein_truth.py pick GFF PROT_FA > LONGEST_FA      (gffread -y proteins, headers = transcript ids; longest per gene)
  protein_truth.py map HITS_M8 > MAP_TSV              (mmseqs: query,target,fident,alnlen,qlen,bits; gene, ENSG, symbol, identity, cover)
  protein_truth.py check MAP_TSV                      (agreement with symbol-inherited ids on named genes; needs DENOVO_* variables)
"""
import sys, os, collections
MIN_ID, MIN_COVER = 0.80, 0.50
mode = sys.argv[1]
if mode == 'pick':
    gene = {}
    for l in open(sys.argv[2]):
        f = l.rstrip('\n').split('\t')
        if len(f) > 8 and f[2] in ('mRNA', 'transcript'):
            a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
            if 'Parent' in a: gene[a['ID']] = a['Parent']
    seq = collections.defaultdict(str); cur = None
    for l in open(sys.argv[3]):
        if l[0] == '>': cur = l[1:].split()[0]
        else: seq[cur] += l.strip().rstrip('.')
    best = {}
    for t, s in seq.items():
        g = gene.get(t)
        if g and len(s) > len(best.get(g, '')): best[g] = s
    for g, s in best.items(): print(f'>{g}\n{s}')
elif mode == 'map':
    best = {}
    for l in open(sys.argv[2]):
        q, t, fid, aln, qlen, bits = l.split('\t')[:6]
        if q not in best or float(bits) > best[q][0]: best[q] = (float(bits), t, float(fid), int(aln) / int(qlen))
    for q, (_, t, fid, cov) in best.items():
        f = t.split('|')
        if fid >= MIN_ID and cov >= MIN_COVER: print(q, f[2].split('.')[0], f[6], f'{fid:.3f}', f'{cov:.3f}', sep='\t')
elif mode == 'check':
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
    import famtruth as F
    os.environ.pop('DENOVO_GGO_ENS_MAP', None); G = F.load_genes('gorilla'); n = ok = new = 0
    for l in open(sys.argv[2]):
        g, e = l.split('\t')[:2]
        if g not in G: continue
        if G[g]['ens']: n += 1; ok += e in G[g]['ens']
        elif G[g]['biotype'] == 'protein_coding': new += 1
    print(f'named genes with a hit {n}, same Ensembl gene as the symbol {ok} = {ok / max(n, 1):.3f} (bar 0.90); protein-coding genes newly given an id {new}')

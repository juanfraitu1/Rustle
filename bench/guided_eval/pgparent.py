"""Independent pseudogene -> parent table for the guided search truth (PREREG GUIDED-S4, docs/LEDGER_O1.md).

  pgparent.py human   TWOWAY_GTF_GZ V50_GENES_GTF OUT_TSV   GENCODE 2-way consensus pseudogenes (Yale PseudoPipe + UCSC RetroFinder,
                                                            parent_id) -> the GENCODE pseudogene gene they overlap (same strand, >= 0.5 of
                                                            the 2-way record) -> CAT/Liftoff CHM13 gene by ENSG
  pgparent.py gorilla OUT_TSV [GFF [GENE_INFO]]                              RefSeq pseudogene description minus a trailing '-like' -> the one protein-coding
                                                            gene with that description (NCBI Gnomon's homology call)

OUT_TSV: pseudogene gene id (famtruth key) <TAB> parent gene id (famtruth key).
"""
import collections
import gzip
import os
import re
import sys
import urllib.parse

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'denovo_eval'))
import famtruth as F  # noqa: E402

ATTR = re.compile(r'(\S+) "([^"]*)"')


def human(twoway, genes, out):
    G = F.load_cat()
    by_ens = collections.defaultdict(set)
    for g, d in G.items():
        for e in d['ens']:
            by_ens[e.split('.')[0]].add(g)
    pg = collections.defaultdict(list)  # (chrom, strand) -> [(s, e, ensg)]
    for line in open(genes):
        f = line.rstrip('\n').split('\t')
        a = dict(ATTR.findall(f[8]))
        if 'pseudogene' in a.get('gene_type', ''):
            pg[(f[0], f[6])].append((int(f[3]) - 1, int(f[4]), a['gene_id'].split('.')[0]))
    n = linked = 0
    pairs = set()
    for line in gzip.open(twoway, 'rt'):
        if line[0] == '#':
            continue
        f = line.rstrip('\n').split('\t')
        if f[2] != 'transcript':
            continue
        n += 1
        a = dict(ATTR.findall(f[8]))
        s, e = int(f[3]) - 1, int(f[4])
        hit = [x for x in pg[(f[0], f[6])] if min(e, x[1]) - max(s, x[0]) >= 0.5 * (e - s)]
        par = by_ens.get(a['parent_id'].split('.')[0], set())
        if len(hit) == 1 and len(par) == 1:
            for g in by_ens.get(hit[0][2], set()):
                pairs.add((g, next(iter(par))))
            linked += 1
    write(out, pairs)
    print(f'human: 2-way pseudogenes {n}, linked {linked}, CAT pseudogene -> parent pairs {len(pairs)}', file=sys.stderr)


def gorilla(out, gff=None, gene_info=None):
    info = {}  # GUIDED-S6: RefSeq release 102 gene lines carry no description; NCBI Gene's gene_info does (by GeneID)
    if gene_info:
        for line in open(gene_info):
            f = line.split('\t')
            if f[8] != '-':
                info[f[1]] = f[8]
    desc_pc, desc_named = collections.defaultdict(set), collections.defaultdict(set)
    pseudo, loc_pc = {}, {}
    for line in open(gff or F.GGO_GFF):
        if line[0] == '#':
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('gene', 'pseudogene'):
            continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        d = urllib.parse.unquote(a.get('description', ''))
        if not d:
            d = info.get(next((x.split(':')[1] for x in a.get('Dbxref', '').split(',') if x.startswith('GeneID:')), ''), '')
        gid = a.get('ID')
        if not d or not gid:
            continue
        if a.get('gene_biotype') == 'protein_coding':
            desc_pc[d].add(gid)
            if a.get('Name', '').startswith('LOC'):
                loc_pc[gid] = re.sub(r'-like$', '', d)
            else:
                desc_named[d].add(gid)
        elif 'pseudogene' in a.get('gene_biotype', ''):
            pseudo[gid] = re.sub(r'-like$', '', d)
    pairs = {(p, next(iter(desc_pc[d]))) for p, d in pseudo.items() if len(desc_pc.get(d, ())) == 1}
    n_pseudo = len(pairs)
    if gff:  # GUIDED-S6: LOC protein-coding genes -> the single NAMED protein-coding gene with their description
        pairs |= {(p, next(iter(desc_named[d]))) for p, d in loc_pc.items() if len(desc_named.get(d, ())) == 1 and p not in desc_named[d]}
    write(out, pairs)
    print(f'gorilla: RefSeq pseudogenes {len(pseudo)}, pseudogene -> parent pairs {n_pseudo}; LOC protein-coding {len(loc_pc)}, '
          f'LOC -> parent pairs {len(pairs) - n_pseudo}', file=sys.stderr)


def write(out, pairs):
    with open(out, 'w') as w:
        for p, q in sorted(pairs):
            w.write(f'{p}\t{q}\n')


if __name__ == '__main__':
    human(*sys.argv[2:5]) if sys.argv[1] == 'human' else gorilla(*sys.argv[2:5])

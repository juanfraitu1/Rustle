"""DENOVOFAM: score a de novo family catalog (mcl_families copies.tsv) against Ensembl Compara recent paralogues and the named
mainstay / Y-ampliconic (YAG) families. Any species: genes from load_genes (CHM13 CAT by default, or --annotation GTF/GFF3).

  famtruth.py score COPIES_TSV OUT_JSON [--paralogs paralogs.tsv] [--levels recent|TAXON,...] [--species S --annotation A --xref X --ids symbol|own]

Copies map to CAT genes by most shared exonic bp on the same strand. A gene is EXPRESSED when a copy maps to it; de novo can only
see expressed genes, so every truth is restricted to them. Predicted families = the catalog families' gene sets (a gene may sit in
several). Truth families = connected components of Compara within-species paralogue pairs at the chosen duplication levels, over
protein-coding genes with an Ensembl ID that Compara knows (the JUDGEABLE universe).
Metrics (family_score's): pairwise sensitivity / precision / F over gene pairs; bipartite family sensitivity / precision / F
(one-to-one matching maximising overlap); per named family (mainstays, YAGs): members expressed, families holding them, largest
share, foreign judgeable genes in those families. Halves: true families split by the SHA-1 of their smallest Ensembl ID.
"""
import argparse
import sys
import collections
import csv
import hashlib
import json
import re
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

# The truth files come from the environment (no machine defaults; docs/DATA.md: the T2T `chm13.draft_v2.0.gene_annotation.gff3` gene
# lines, the bundle's genes.tsv, a BioMart paralogue export of release 116):
import os as _os
CAT_GFF = _os.environ.get('DENOVO_CAT_GFF')
CAT_TSV = _os.environ.get('DENOVO_CAT_TSV')
PARALOGS = _os.environ.get('DENOVO_PARALOGS')


def need(path, var):
    if not path:
        sys.exit(f'famtruth: set {var}')
    return path
LEVELS = {
    'recent': ['Homo sapiens', 'Homininae', 'Hominidae', 'Hominoidea', 'Catarrhini'],
    'simiiformes': ['Homo sapiens', 'Homininae', 'Hominidae', 'Hominoidea', 'Catarrhini', 'Simiiformes'],
    'primates': ['Homo sapiens', 'Homininae', 'Hominidae', 'Hominoidea', 'Catarrhini', 'Simiiformes', 'Haplorrhini', 'Primates'],
    # family LEVEL of a hierarchical output (2026-10-30, user "could it be that some families are actually subfamilies?"):
    'eutheria': ['Homo sapiens', 'Homininae', 'Hominidae', 'Hominoidea', 'Catarrhini', 'Simiiformes', 'Haplorrhini', 'Primates',
                 'Euarchontoglires', 'Boreoeutheria', 'Eutheria'],
    'mammalia': ['Homo sapiens', 'Homininae', 'Hominidae', 'Hominoidea', 'Catarrhini', 'Simiiformes', 'Haplorrhini', 'Primates',
                 'Euarchontoglires', 'Boreoeutheria', 'Eutheria', 'Theria', 'Mammalia'],
}
# bench/o3_graph/protfam.py KNOWN: the mainstay set and IsoCon's nine Y-ampliconic families
MAINSTAYS = {"NPIP": r"^NPIP[AB]", "TBC1D3": r"^TBC1D3([A-L]|P[0-9]+)?$", "GOLGA8": r"^GOLGA8", "GOLGA6L": r"^GOLGA6L",
             "CT47A": r"^CT47A", "XAGE": r"^XAGE", "SPATA31": r"^SPATA31"}
# Gorilla mainstays (2026-10-30, user: Yoo 2025 ape expansions "MAPKBP1, JMJD7-PLA2G4B, SPTBN5 in chr1 in gorilla ... as well as
# LRPAP1"): the copies are RefSeq LOC genes, so a family is defined by the gene symbol OR the RefSeq description (any biotype: Yoo's
# copies include pseudogene-annotated ones). Matched against `name` and `desc` of famtruth.load_ggo genes.
GGO_MAINSTAYS = {"LRPAP1": r"^LRPAP1$|alpha-2-macroglobulin receptor-associated protein",
                 "MAPKBP1": r"^MAPKBP1$|mitogen-activated protein kinase[- ]binding protein 1",
                 "JMJD7-PLA2G4B": r"^JMJD7$|^PLA2G4B$|lysyl hydroxylase JMJD7|cytosolic phospholipase A2 beta",  # one readthrough unit in Yoo and in our loci
                 "SPTBN5": r"^SPTBN5$|spectrin beta chain, non-erythrocytic 5"}
YAGS = {"BPY2": r"^BPY2[BC]?$", "CDY": r"^CDY[12][AB]?$", "DAZ": r"^DAZ[1-4]$", "HSFY": r"^HSFY[12]$", "PRY": r"^PRY2?$",
        "RBMY": r"^RBMY1[A-J][0-9]?$", "TSPY": r"^TSPY[0-9]+$", "VCY": r"^VCY(1B)?$", "XKRY": r"^XKRY2?$"}


def blocks(s):
    return [tuple(map(int, b.split('-'))) for b in s.split(',') if b]


GENES = {}


def load_cat():
    """gene id -> dict(name, biotype, chrom, strand, exons, ens: set of unversioned Ensembl gene ids)."""
    ens, collapsed = {}, set()
    for line in open(need(CAT_GFF, 'DENOVO_CAT_GFF')):
        f = line.split('\t', 8)
        if len(f) < 9 or f[2] != 'gene':
            continue
        a = dict(x.split('=', 1) for x in f[8].strip().split(';') if '=' in x)
        ids = {a.get('source_gene', 'N/A').split('.')[0]}
        ids |= {x.split('.')[0] for x in a.get('collapsed_gene_ids', 'N/A').split(',')}
        ens[a['gene_id']] = {x for x in ids if x.startswith('ENSG')}
        if a.get('collapsed_gene_ids', 'N/A') != 'N/A':
            collapsed.add(a['gene_id'])
    genes = {}
    for line in open(need(CAT_TSV, 'DENOVO_CAT_TSV')):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'gene_id':
            continue
        genes[f[0]] = dict(name=f[1], biotype=f[2], chrom=f[4], strand=f[7], exons=blocks(f[9]), ens=ens.get(f[0], set()))
    if _os.environ.get('DENOVO_SPLIT_REFSEQ'):
        genes = split_collapsed(genes, collapsed, _os.environ['DENOVO_SPLIT_REFSEQ'])
    GENES.update(genes)
    return genes


def split_collapsed(genes, collapsed, refseq_gtf):
    """register 1494: a protein-coding CAT gene that collapsed several source genes becomes one truth gene per coding RefSeq gene
    (NM_/XM_) on its strand sharing exonic bases with it whose HGNC symbol maps to one of its ENSG ids, when there are >= 2."""
    import gzip
    sym2e = {}
    for r in csv.DictReader(open(_os.environ['DENOVO_HGNC']), delimiter='\t'):
        if r['ensembl_gene_id']:
            sym2e[r['symbol']] = r['ensembl_gene_id']
    ex, meta, coding = collections.defaultdict(list), {}, set()
    for line in (gzip.open(refseq_gtf, 'rt') if refseq_gtf.endswith('.gz') else open(refseq_gtf)):
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        mg = re.search(r'gene_id "([^"]+)"', f[8])
        if not mg:
            continue
        g = mg.group(1)
        ex[g].append((int(f[3]) - 1, int(f[4])))
        meta[g] = (f[0], f[6], (re.search(r'gene_name "([^"]+)"', f[8]) or re.search(r'gene_id "([^"]+)"', f[8])).group(1))
        if re.search(r'transcript_id "rna-[NX]M_', f[8]):
            coding.add(g)
    rs = {}
    for g in coding:
        u = []
        for a, b in sorted(ex[g]):
            if u and a <= u[-1][1]:
                u[-1][1] = max(u[-1][1], b)
            else:
                u.append([a, b])
        rs[g] = dict(chrom=meta[g][0], strand=meta[g][1], name=meta[g][2], exons=[tuple(x) for x in u])
    ridx = exon_index(rs)
    out, n_split, n_sub = {}, 0, 0
    for gid, d in genes.items():
        if gid not in collapsed or d['biotype'] != 'protein_coding':
            out[gid] = d
            continue
        hits = set()
        k = (d['chrom'], d['strand'])
        if k in ridx:
            import bisect
            v, st, m = ridx[k]
            for s0, e0 in d['exons']:
                i = bisect.bisect_left(st, s0 - m)
                while i < len(v) and v[i][0] < e0:
                    if v[i][1] > s0:
                        hits.add(v[i][2])
                    i += 1
        subs = {}
        for r in hits:
            e = sym2e.get(rs[r]['name'])
            if e and e in d['ens']:
                subs.setdefault(rs[r]['name'], (e, r))
        if len(subs) < 2:
            out[gid] = d
            continue
        n_split += 1
        for sym, (e, r) in sorted(subs.items()):
            n_sub += 1
            out[f'{gid}#{sym}'] = dict(name=sym, biotype='protein_coding', chrom=d['chrom'], strand=d['strand'], exons=rs[r]['exons'], ens={e})
    print(f'# split collapsed (register 1494): {n_split} CAT genes -> {n_sub} RefSeq genes', file=sys.stderr)
    return out


GGO_GFF = _os.environ.get('DENOVO_GGO_GFF')


def load_ggo(human, path=None):
    """RefSeq-style GFF3 genes (gorilla first, PREREG DENOVOFAM-G; any RefSeq / NCBI / Liftoff GFF3): exon union per gene; `ens` = the
    Ensembl ids of the human CAT gene(s) of the same symbol, so the human Compara truth applies by symbol. LOC genes get no ids."""
    sym_ens = collections.defaultdict(set)
    for d in human.values():
        sym_ens[d['name']] |= d['ens']
    genes, tx_gene, ex = {}, {}, collections.defaultdict(list)
    for line in _open(need(path or GGO_GFF, 'DENOVO_GGO_GFF or --annotation')):
        if line[0] == '#':
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if f[2] in ('gene', 'pseudogene', 'ncRNA_gene'):
            name = a.get('gene', a.get('gene_name', a.get('Name', a['ID'])))
            genes[a['ID']] = dict(name=name, biotype=a.get('gene_biotype', a.get('biotype', a.get('gene_type', 'other'))), chrom=f[0], strand=f[6], exons=[],
                                  ens=set() if name.startswith('LOC') else sym_ens.get(name, set()),
                                  desc=a.get('description', '').replace('%2C', ','))
        elif f[2] == 'exon':
            ex[a.get('Parent', '')].append((int(f[3]) - 1, int(f[4])))
        elif 'Parent' in a and a['Parent'] in genes:
            tx_gene[a['ID']] = a['Parent']
    by = collections.defaultdict(list)
    for t, v in ex.items():
        g = tx_gene.get(t, t if t in genes else None)
        if g in genes:
            by[g] += v
    for g, v in by.items():
        m = []
        for s0, e0 in sorted(v):
            if m and s0 <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e0)
            else:
                m.append([s0, e0])
        genes[g]['exons'] = [tuple(x) for x in m]
    genes = {g: d for g, d in genes.items() if d['exons']}
    GENES.clear()
    GENES.update(genes)
    return genes


def _open(path):
    import gzip
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def sniff(path):
    """'gtf', 'ensembl' (GFF3 with ID=gene:...) or 'gff3' (RefSeq / NCBI / Liftoff style)."""
    for line in _open(path):
        f = line.split('\t')
        if line[0] == '#' or len(f) < 9:
            continue
        if 'gene_id "' in f[8]:
            return 'gtf'
        if 'ID=gene:' in f[8]:
            return 'ensembl'
        if f[2] in ('gene', 'pseudogene', 'ncRNA_gene'):
            return 'gff3'
    return 'gff3'


def load_gtf(path):
    """GTF genes: exon union per gene_id; name gene_name, biotype gene_biotype / gene_type (else coding by an NM_/XM_ transcript)."""
    genes, ex = {}, collections.defaultdict(list)
    for line in _open(path):
        f = line.rstrip('\n').split('\t')
        if line[0] == '#' or len(f) < 9 or f[2] != 'exon':
            continue
        a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
        g = a['gene_id']
        if g not in genes:
            genes[g] = dict(name=a.get('gene_name', g), biotype=a.get('gene_biotype', a.get('gene_type', 'other')), chrom=f[0],
                            strand=f[6], exons=[], ens=set())
        ex[g].append((int(f[3]) - 1, int(f[4])))
        if genes[g]['biotype'] == 'other' and (re.match(r'(rna-)?[NX]M_', a.get('transcript_id', '')) or
                                               a.get('transcript_biotype', a.get('transcript_type')) == 'protein_coding'):
            genes[g]['biotype'] = 'protein_coding'  # no gene biotype (RefSeq GTF): a coding transcript (NM_/XM_) makes it coding
    for g, v in ex.items():
        m = []
        for s0, e0 in sorted(v):
            if m and s0 <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e0)
            else:
                m.append([s0, e0])
        genes[g]['exons'] = [tuple(x) for x in m]
    GENES.clear()
    GENES.update(genes)
    return genes


def load_ensembl(path):
    """Ensembl GFF3 genes (ID=gene:..., biotype=...; transcripts Parent=gene:...; exons Parent=transcript:...): exon union per gene."""
    import urllib.parse
    genes, tx_gene, ex = {}, {}, collections.defaultdict(list)
    for line in _open(path):
        if line[0] == '#':
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        if a.get('ID', '').startswith('gene:'):
            genes[a['ID']] = dict(name=a['ID'], biotype=a.get('biotype', 'other'), chrom=f[0], strand=f[6], exons=[], ens=set(),
                                  description=urllib.parse.unquote(a.get('description', '')), symbol=a.get('Name'))
        elif a.get('ID', '').startswith('transcript:') and a.get('Parent', '').startswith('gene:'):
            tx_gene[a['ID']] = a['Parent']
        elif f[2] == 'exon':
            ex[a.get('Parent', '')].append((int(f[3]) - 1, int(f[4])))
    for t, v in ex.items():
        g = tx_gene.get(t)
        if g in genes:
            genes[g]['exons'] += v
    for d in genes.values():
        m = []
        for s0, e0 in sorted(d['exons']):
            if m and s0 <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e0)
            else:
                m.append([s0, e0])
        d['exons'] = [tuple(x) for x in m]
    genes = {g: d for g, d in genes.items() if d['exons']}
    GENES.clear()
    GENES.update(genes)
    return genes


def unnamed(d):
    """Outside the by-symbol truth: RefSeq LOC genes, Ensembl genes without a symbol."""
    return d['name'].startswith('LOC') or re.match(r'ENS[A-Z]*G[0-9]{6,}', d['name']) is not None


# Species presets: name -> (annotation env var, xref env var). Any other species: --annotation (+ --xref, --ids).
PRESETS = {'gorilla': ('DENOVO_GGO_GFF', None), 'hg002': ('DENOVO_HG002_GFF', 'DENOVO_HG002_XREF')}


def load_genes(species='human', annotation=None, xref=None, ids='symbol'):
    """The annotated genes a run is scored against (gene id -> name, biotype, chrom, strand, exons, ens).
    species human without --annotation: CHM13 CAT/Liftoff (DENOVO_CAT_GFF/TSV). Otherwise ANNOTATION is any GTF, Ensembl GFF3 or
    RefSeq-style GFF3 of the sample's genome (presets: gorilla = DENOVO_GGO_GFF, hg002 = DENOVO_HG002_GFF + DENOVO_HG002_XREF).
    ids: 'symbol' = a gene inherits the human Ensembl ids of the CAT gene(s) of its symbol (the paralogue table is human Compara);
    'own' = a gene's own unversioned id is its Ensembl id (the paralogue table is that species' Compara export; no CAT needed).
    xref: Ensembl xref TSV (column 1 stable id, 5 symbol, 7 'HGNC Symbol') naming Ensembl genes."""
    if species in PRESETS:
        annotation = annotation or _os.environ.get(PRESETS[species][0])
        if PRESETS[species][1]:
            xref = xref or _os.environ.get(PRESETS[species][1])
    if species == 'human' and not annotation:
        return load_cat()
    path = need(annotation, f'--annotation (or {PRESETS[species][0]})' if species in PRESETS else '--annotation')
    kind = sniff(path)
    if kind == 'gff3' and ids == 'symbol':
        return load_ggo(load_cat(), path)
    G = load_gtf(path) if kind == 'gtf' else load_ensembl(path) if kind == 'ensembl' else load_ggo({}, path)
    if kind == 'ensembl':
        sym = {}
        if xref:
            for line in open(xref):
                f = line.rstrip('\n').split('\t')
                if len(f) > 6 and f[6] == 'HGNC Symbol':
                    sym[f[0]] = f[4]
        for g, d in G.items():
            bare = g.split(':')[-1]
            d['name'] = sym.get(bare, bare) if xref else (d.pop('symbol', None) or bare)
            d.pop('symbol', None)
    if ids == 'own':
        for g, d in G.items():
            d['ens'] = {re.sub(r'^(gene[:-])', '', g).split('.')[0]}
    else:
        sym_ens = collections.defaultdict(set)
        for d in load_cat().values():
            sym_ens[d['name']] |= d['ens']
        for d in G.values():
            d['ens'] = set() if unnamed(d) else sym_ens.get(d['name'], set())
    GENES.clear()
    GENES.update(G)
    return G


SPECIES_FLAGS = ('--species', '--annotation', '--xref', '--ids')


def genes_from_argv(a, default='human'):
    """load_genes from a hand-parsed argv (pairscore.py, seqscore.py)."""
    opt = lambda k, d=None: a[a.index(k) + 1] if k in a else d
    return load_genes(opt('--species', default), opt('--annotation'), opt('--xref'), opt('--ids', 'symbol'))


def levels(name):
    """A LEVELS key, or a comma list of Compara 'last common ancestor' taxa (any species: 'Mus musculus,Murinae,Muridae')."""
    return LEVELS[name] if name in LEVELS else name.split(',')


def exon_index(genes):
    by = collections.defaultdict(list)
    for g, d in genes.items():
        for s, e in d['exons']:
            by[(d['chrom'], d['strand'])].append((s, e, g))
    out = {}
    for k, v in by.items():
        v.sort()
        out[k] = (v, [x[0] for x in v], max(x[1] - x[0] for x in v))
    return out


def map_copy(idx, chrom, strand, exons):
    """CAT gene sharing the most exonic bp with the copy on its strand (None if none)."""
    import bisect
    if (chrom, strand) not in idx:
        return None
    v, st, m = idx[(chrom, strand)]
    share = collections.Counter()
    for s, e in exons:
        i = bisect.bisect_left(st, s - m)
        while i < len(v) and v[i][0] < e:
            a, b, g = v[i]
            if b > s:
                share[g] += min(b, e) - max(a, s)
            i += 1
    if not share:
        return None
    # protein-coding first (PREREG DENOVOFAM Amendment 1): a long non-coding or "unknown_likely_coding" CAT model spanning a
    # coding gene (CHM13_G1000495 over DAZ2, AC091057.3 over GOLGA8Q, PMS2P2 over SPDYE15) otherwise takes the copy away
    # from the gene the truth is about
    pc = [g for g in share if GENES.get(g, {}).get('biotype') == 'protein_coding']
    return max(pc, key=lambda g: (share[g], g)) if pc else max(share, key=lambda g: (share[g], g))


def load_catalog(path, genes, idx):
    fams = collections.defaultdict(set)
    n_copies = n_mapped = 0
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'family_id':
            hdr = {k: i for i, k in enumerate(f)}
            continue
        n_copies += 1
        strand = f[hdr['strand']] if f[hdr['strand']] in '+-' else '+'
        g = map_copy(idx, f[hdr['chrom']], strand, blocks(f[hdr['exons']]))
        if g:
            n_mapped += 1
            fams[f[hdr['family_id']]].add(g)
    return dict(fams), n_copies, n_mapped


def expressed_from_gtf(path, idx):
    """Genes carrying a de novo transcript (exons of each transcript_id mapped like a copy)."""
    tx = collections.defaultdict(list)
    meta = {}
    for line in open(path):
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'exon':
            continue
        t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        tx[t].append((int(f[3]) - 1, int(f[4])))
        meta[t] = (f[0], f[6] if f[6] in '+-' else '+')
    out = set()
    for t, ex in tx.items():
        g = map_copy(idx, meta[t][0], meta[t][1], sorted(ex))
        if g:
            out.add(g)
    return out


def truth_pairs(genes, paralogs, levels):
    by_ens = collections.defaultdict(set)
    for g, d in genes.items():
        if d['biotype'] == 'protein_coding':
            for e in d['ens']:
                by_ens[e].add(g)
    known_ens = set()
    pairs = set()
    lv = set(levels)
    for line in _open(need(paralogs, 'DENOVO_PARALOGS or --paralogs')):
        f = line.rstrip('\n').split('\t')
        if f[0].startswith('Gene'):
            continue
        known_ens.add(f[0])
        if len(f) < 4 or not f[1] or f[3] not in lv:
            continue
        for a in by_ens.get(f[0], ()):
            for b in by_ens.get(f[1], ()):
                if a != b:
                    pairs.add((min(a, b), max(a, b)))
    universe = {g for e, gs in by_ens.items() if e in known_ens for g in gs}
    return pairs, universe


def components(pairs):
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in pairs:
        parent[find(a)] = find(b)
    comp = collections.defaultdict(set)
    for x in parent:
        comp[find(x)].add(x)
    return list(comp.values())


def bipartite(true_fams, pred_fams):
    true_fams = [t for t in true_fams if len(t) >= 2]
    pred_fams = [p for p in pred_fams if len(p) >= 2]
    if not true_fams or not pred_fams:
        return dict(sens=0.0, prec=0.0, f=0.0, overlap=0, true_genes=sum(map(len, true_fams)), matched_pred_genes=0)
    M = np.array([[len(t & p) for p in pred_fams] for t in true_fams])
    r, c = linear_sum_assignment(-M)
    ov = int(M[r, c].sum())
    tg = sum(len(t) for t in true_fams)
    pg = sum(len(pred_fams[j]) for i, j in zip(r, c) if M[i, j] > 0)
    s, p = ov / tg, ov / max(pg, 1)
    return dict(sens=s, prec=p, f=2 * s * p / max(s + p, 1e-12), overlap=ov, true_genes=tg, matched_pred_genes=pg)


def score(true_fams, pred_fams, universe, half=None):
    """Pairwise and bipartite, both restricted to `universe`. With `half` (a gene subset of `universe`), a predicted pair counts
    when it touches >= 1 gene of the half (a pair joining the half to the other half is false and must count against it; audit
    2026-10-27), and the bipartite matching sees predicted families over the whole universe."""
    tp_pairs = {(min(a, b), max(a, b)) for t in true_fams for a in t for b in t if a < b and a in universe and b in universe}
    pr = [p & universe for p in pred_fams]
    pp_pairs = {(min(a, b), max(a, b)) for p in pr for a in p for b in p if a < b}
    if half is not None:
        pp_pairs = {q for q in pp_pairs if q[0] in half or q[1] in half}
        pr = [p for p in pr if p & half]
    hit = len(tp_pairs & pp_pairs)
    s, p = hit / max(len(tp_pairs), 1), hit / max(len(pp_pairs), 1)
    out = dict(pair_true=len(tp_pairs), pair_pred=len(pp_pairs), pair_hit=hit, pair_sens=s, pair_prec=p,
               pair_f=2 * s * p / max(s + p, 1e-12))
    out.update({f'fam_{k}': v for k, v in bipartite([t & universe for t in true_fams], pr).items()})
    return out


def named(genes, expressed, fams, universe, patterns):
    gfam = collections.defaultdict(set)
    for fid, gs in fams.items():
        for g in gs:
            gfam[g].add(fid)
    rep = {}
    for name, pat in patterns.items():
        rx = re.compile(pat)
        members = {g for g, d in genes.items() if rx.match(d['name']) and d['biotype'] in ('protein_coding',) or
                   (rx.match(d['name']) and 'pseudogene' not in d['biotype'] and d['biotype'] != 'lncRNA')}
        members = {g for g in members if rx.match(genes[g]['name'])}
        exp = members & expressed
        held = collections.Counter(f for g in exp for f in gfam[g])
        best = held.most_common(1)[0] if held else (None, 0)
        foreign = {g for f in held for g in fams[f] if g not in members and g in universe}
        rep[name] = dict(members=len(members), expressed=len(exp), families=len(held), largest=best[1],
                         largest_family=best[0], foreign=len(foreign), foreign_names=sorted(genes[g]['name'] for g in foreign)[:12],
                         in_one=len(exp) >= 2 and best[1] == len(exp))
    return rep


def half(fam):
    key = min(e for g in fam for e in G[g]['ens']) if any(G[g]['ens'] for g in fam) else min(fam)
    return int(hashlib.sha1(key.encode()).hexdigest(), 16) % 2


G = {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('cmd')
    ap.add_argument('copies')
    ap.add_argument('out')
    ap.add_argument('--paralogs', default=PARALOGS)
    ap.add_argument('--levels', default='recent')
    ap.add_argument('--show', default='all,dev,heldout', help='halves to PRINT (the JSON always holds all): dev only while tuning')
    ap.add_argument('--species', default='human', help='human (CAT), a preset (gorilla, hg002) or any name with --annotation')
    ap.add_argument('--annotation', '--ggo-gff', dest='annotation', help='GTF / GFF3 of the sample genome (see load_genes)')
    ap.add_argument('--xref', help='Ensembl xref TSV naming Ensembl GFF3 genes (HGNC symbols)')
    ap.add_argument('--ids', default='symbol', choices=['symbol', 'own'], help='symbol: human Compara by symbol; own: the genes\' own ids')
    ap.add_argument('--gtf', help='de novo GTF of the sample: expressed = genes with a transcript (default: genes in a family)')
    a = ap.parse_args()
    global G
    G = load_genes(a.species, a.annotation, a.xref, a.ids)
    idx = exon_index(G)
    fams, n_copies, n_mapped = load_catalog(a.copies, G, idx)
    in_fam = set().union(*fams.values()) if fams else set()
    expressed = expressed_from_gtf(a.gtf, idx) | in_fam if a.gtf else in_fam
    res = dict(copies=n_copies, copies_on_genes=n_mapped, families=len(fams), expressed_genes=len(expressed), genes_in_families=len(in_fam))
    for lvname in [a.levels] + [x for x in LEVELS if x != a.levels and a.levels in LEVELS]:
        pairs, universe = truth_pairs(G, a.paralogs, levels(lvname))
        uexp = universe & expressed
        tf = [c & uexp for c in components({p for p in pairs if p[0] in uexp and p[1] in uexp})]
        tf = [c for c in tf if len(c) >= 2]
        pf = list(fams.values())
        r = dict(universe=len(universe), universe_expressed=len(uexp), true_families=len(tf),
                 all=score(tf, pf, uexp))
        for h in (0, 1):
            th = [t for t in tf if half(t) == h]
            hg = set().union(*th) if th else set()
            # each half: its true families' genes plus singletons hashed to it
            sing = {g for g in uexp - set().union(*tf) if int(hashlib.sha1(min(G[g]['ens']).encode()).hexdigest(), 16) % 2 == h}
            r['dev' if h == 0 else 'heldout'] = score(th, pf, uexp, half=hg | sing)
        res[lvname] = r
    _, universe = truth_pairs(G, a.paralogs, levels('recent'))
    res['mainstays'] = named(G, expressed, fams, universe, MAINSTAYS)
    res['yags'] = named(G, expressed, fams, universe, YAGS)
    json.dump(res, open(a.out, 'w'), indent=1)
    r = res[a.levels]
    print(f"copies {n_copies} on genes {n_mapped}; families {len(fams)}; expressed genes {len(expressed)}; "
          f"judgeable expressed {r['universe_expressed']}; true families ({a.levels}) {r['true_families']}")
    for k in a.show.split(','):
        x = r[k]
        print(f"  {k:8s} pair sens {x['pair_sens']:.3f} prec {x['pair_prec']:.3f} F {x['pair_f']:.3f} ({x['pair_hit']}/{x['pair_true']} | "
              f"{x['pair_pred']})  family sens {x['fam_sens']:.3f} prec {x['fam_prec']:.3f} F {x['fam_f']:.3f}")
    for grp in ('mainstays', 'yags'):
        print(' ', grp, '; '.join(f"{n} {v['largest']}/{v['expressed']}({v['members']}) fam{v['families']} for{v['foreign']}"
                               for n, v in res[grp].items()))


if __name__ == '__main__':
    main()

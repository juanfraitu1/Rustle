"""DENOVOFAM: score a de novo family catalog (mcl_families copies.tsv) against Ensembl Compara recent paralogues and the named
mainstay / Y-ampliconic (YAG) families. Human only (CHM13 v2.0, CAT/Liftoff v2.0).

  famtruth.py score COPIES_TSV OUT_JSON [--paralogs paralogs.tsv] [--levels recent]

Copies map to CAT genes by most shared exonic bp on the same strand. A gene is EXPRESSED when a copy maps to it; de novo can only
see expressed genes, so every truth is restricted to them. Predicted families = the catalog families' gene sets (a gene may sit in
several). Truth families = connected components of Compara within-species paralogue pairs at the chosen duplication levels, over
protein-coding genes with an Ensembl ID that Compara knows (the JUDGEABLE universe).
Metrics (family_score's): pairwise sensitivity / precision / F over gene pairs; bipartite family sensitivity / precision / F
(one-to-one matching maximising overlap); per named family (mainstays, YAGs): members expressed, families holding them, largest
share, foreign judgeable genes in those families. Halves: true families split by the SHA-1 of their smallest Ensembl ID.
"""
import argparse
import collections
import hashlib
import json
import re
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

CAT_GFF = '/mnt/linuxdisk/home/juanfraitu/winloci_data/gencode_chm13/chm13v2.0_gencode.gff3'
CAT_TSV = '/mnt/linuxdisk/home/juanfraitu/winloci_data/gencode_chm13/chm13v2.0_CAT_Liftoff.genes.tsv'
PARALOGS = '/mnt/linuxdisk/tmp/denovo_eval/paralogs.tsv'
LEVELS = {
    'recent': ['Homo sapiens', 'Homininae', 'Hominidae', 'Hominoidea', 'Catarrhini'],
    'simiiformes': ['Homo sapiens', 'Homininae', 'Hominidae', 'Hominoidea', 'Catarrhini', 'Simiiformes'],
    'primates': ['Homo sapiens', 'Homininae', 'Hominidae', 'Hominoidea', 'Catarrhini', 'Simiiformes', 'Haplorrhini', 'Primates'],
}
# bench/o3_graph/protfam.py KNOWN: the mainstay set and IsoCon's nine Y-ampliconic families
MAINSTAYS = {"NPIP": r"^NPIP[AB]", "TBC1D3": r"^TBC1D3([A-L]|P[0-9]+)?$", "GOLGA8": r"^GOLGA8", "GOLGA6L": r"^GOLGA6L",
             "CT47A": r"^CT47A", "XAGE": r"^XAGE", "SPATA31": r"^SPATA31"}
YAGS = {"BPY2": r"^BPY2[BC]?$", "CDY": r"^CDY[12][AB]?$", "DAZ": r"^DAZ[1-4]$", "HSFY": r"^HSFY[12]$", "PRY": r"^PRY2?$",
        "RBMY": r"^RBMY1[A-J][0-9]?$", "TSPY": r"^TSPY[0-9]+$", "VCY": r"^VCY(1B)?$", "XKRY": r"^XKRY2?$"}


def blocks(s):
    return [tuple(map(int, b.split('-'))) for b in s.split(',') if b]


def load_cat():
    """gene id -> dict(name, biotype, chrom, strand, exons, ens: set of unversioned Ensembl gene ids)."""
    ens = {}
    for line in open(CAT_GFF):
        f = line.split('\t', 8)
        if len(f) < 9 or f[2] != 'gene':
            continue
        a = dict(x.split('=', 1) for x in f[8].strip().split(';') if '=' in x)
        ids = {a.get('source_gene', 'N/A').split('.')[0]}
        ids |= {x.split('.')[0] for x in a.get('collapsed_gene_ids', 'N/A').split(',')}
        ens[a['gene_id']] = {x for x in ids if x.startswith('ENSG')}
    genes = {}
    for line in open(CAT_TSV):
        f = line.rstrip('\n').split('\t')
        if f[0] == 'gene_id':
            continue
        genes[f[0]] = dict(name=f[1], biotype=f[2], chrom=f[4], strand=f[7], exons=blocks(f[9]), ens=ens.get(f[0], set()))
    return genes


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
    return share.most_common(1)[0][0] if share else None


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
    for line in open(paralogs):
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


def score(true_fams, pred_fams, universe):
    """Pairwise and bipartite, both restricted to `universe`."""
    tp_pairs = {(min(a, b), max(a, b)) for t in true_fams for a in t for b in t if a < b and a in universe and b in universe}
    pr = [p & universe for p in pred_fams]
    pp_pairs = {(min(a, b), max(a, b)) for p in pr for a in p for b in p if a < b}
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
    ap.add_argument('--gtf', help='de novo GTF of the sample: expressed = genes with a transcript (default: genes in a family)')
    a = ap.parse_args()
    global G
    G = load_cat()
    idx = exon_index(G)
    fams, n_copies, n_mapped = load_catalog(a.copies, G, idx)
    in_fam = set().union(*fams.values()) if fams else set()
    expressed = expressed_from_gtf(a.gtf, idx) | in_fam if a.gtf else in_fam
    res = dict(copies=n_copies, copies_on_genes=n_mapped, families=len(fams), expressed_genes=len(expressed), genes_in_families=len(in_fam))
    for lvname in [a.levels] + [x for x in LEVELS if x != a.levels]:
        pairs, universe = truth_pairs(G, a.paralogs, LEVELS[lvname])
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
            r['dev' if h == 0 else 'heldout'] = score(th, pf, hg | sing)
        res[lvname] = r
    _, universe = truth_pairs(G, a.paralogs, LEVELS['recent'])
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

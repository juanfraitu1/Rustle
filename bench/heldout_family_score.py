#!/usr/bin/env python3
"""Score `mcl_families` clusters against symbol-root truth families, per
`docs/PREREG_heldout_families_2026-09-20.md` (md5 of the pre-registration is quoted in the report).

Implements the PRE-REGISTERED definitions verbatim; nothing here is tunable:

  truth family  = >= 3 RefSeq `gene`/`pseudogene` records on the chromosome whose `Name=` shares a root,
                  root = re.sub(r'(?:P\\d+|\\d+|[A-Z])$', '', symbol) applied ONCE (the amended rule),
                  roots shorter than 3 chars dropped, symbols starting `LOC` excluded from TRUTH but
                  left in the INPUT so they can still cost precision.
  matching      = one-to-one bipartite, maximising total overlap (scipy linear_sum_assignment).
  sensitivity   = matched / truth members;  precision = matched / members of the matched cluster;
  F             = harmonic mean.  Unmatched truth families score 0 and are KEPT in the pooled mean.

Usage:
  heldout_family_score.py --gff chm13.gff --clusters chrN_fam.clusters.tsv --chrom chrN [--json out.json]
"""
import argparse
import collections
import json
import re
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment

SUFFIX = re.compile(r'(?:P\d+|\d+|[A-Z])$')


def symbol_root(sym):
    """The pre-registered root: strip ONE trailing copy-suffix. Not applied repeatedly."""
    return SUFFIX.sub('', sym)


def load_genes(gff, chrom):
    """(start1, end) -> symbol, for gene/pseudogene records on `chrom` carrying a Name."""
    out = {}
    name_re = re.compile(r'Name=([^;]+)')
    with open(gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
                continue
            m = name_re.search(f[8])
            if m:
                out[(int(f[3]), int(f[4]))] = m.group(1)
    return out


def truth_families(genes):
    """root -> [symbols]; >= 3 members, root >= 3 chars, LOC* excluded from truth."""
    by_root = collections.defaultdict(set)
    for sym in genes.values():
        if sym.startswith('LOC'):
            continue
        r = symbol_root(sym)
        if len(r) >= 3:
            by_root[r].add(sym)
    return {r: sorted(v) for r, v in by_root.items() if len(v) >= 3}


def soto_families(s1c, genes, chrom):
    """Soto et al. 2025 published families (S1C `Family ID`), restricted to >= 3 members on `chrom`.

    External, published, SD-derived truth on the same CHM13 assembly. Matched to RefSeq by `Gene Name`.
    A gene carrying more than one distinct Family ID is EXCLUDED (the project's settled rule: a
    partition needs one label per gene; see bench/soto/soto_score_against_truth.py).
    """
    import csv as _csv
    on_chrom = set(genes.values())
    ids = collections.defaultdict(set)
    for r in _csv.DictReader(open(s1c), delimiter='\t'):
        fid = (r.get('Family ID') or '').strip()
        nm = (r.get('Gene Name') or '').strip()
        if fid and fid != 'N/A' and nm:
            ids[nm].add(fid)
    fam = collections.defaultdict(set)
    for nm, fids in ids.items():
        if len(fids) != 1 or nm not in on_chrom:
            continue
        fam[next(iter(fids))].add(nm)
    return {f: sorted(v) for f, v in fam.items() if len(v) >= 3}


def predicted_clusters(clusters_tsv, genes):
    """cluster_id -> [symbols] (members that resolve to a named gene; LOC members are KEPT)."""
    out = collections.defaultdict(list)
    with open(clusters_tsv) as fh:
        for line in fh:
            if line.startswith('cluster_id'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 8:
                continue
            cid, s, e = p[0], int(p[6]), int(p[7])
            sym = genes.get((s + 1, e)) or genes.get((s, e))
            out[cid].append(sym if sym else f'{p[5]}:{s}-{e}')
    return dict(out)


def score(truth, pred):
    """One-to-one bipartite match maximising overlap; unmatched truth families score 0."""
    troots, cids = sorted(truth), sorted(pred)
    if not troots:
        return None
    ov = np.zeros((len(troots), len(cids)), dtype=int)
    for i, r in enumerate(troots):
        tset = set(truth[r])
        for j, c in enumerate(cids):
            ov[i, j] = len(tset & set(pred[c]))
    rows, cols = linear_sum_assignment(-ov) if cids else ([], [])
    matched = {troots[i]: cids[j] for i, j in zip(rows, cols) if ov[i, j] > 0}
    per = {}
    for r in troots:
        t = set(truth[r])
        c = matched.get(r)
        if c is None:
            per[r] = dict(n_truth=len(t), n_pred=0, hit=0, sens=0.0, prec=0.0, f=0.0, cluster=None)
            continue
        p = set(pred[c])
        hit = len(t & p)
        sens = hit / len(t)
        prec = hit / len(p) if p else 0.0
        f = 0.0 if sens + prec == 0 else 2 * sens * prec / (sens + prec)
        per[r] = dict(n_truth=len(t), n_pred=len(p), hit=hit, sens=round(sens, 4),
                      prec=round(prec, 4), f=round(f, 4), cluster=c)
    return per


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--gff', required=True)
    ap.add_argument('--clusters', required=True)
    ap.add_argument('--chrom', required=True)
    ap.add_argument('--json')
    ap.add_argument('--soto', help='score against Soto S1C published families instead of symbol roots')
    a = ap.parse_args()

    genes = load_genes(a.gff, a.chrom)
    truth = soto_families(a.soto, genes, a.chrom) if a.soto else truth_families(genes)
    pred = predicted_clusters(a.clusters, genes)
    per = score(truth, pred)
    if per is None:
        print(f'{a.chrom}: NO TRUTH FAMILIES (>=3 members) — chromosome not scoreable')
        return
    fs = [v['f'] for v in per.values()]
    exact = sum(1 for v in per.values() if v['f'] == 1.0)
    found = sum(1 for v in per.values() if v['cluster'])
    summary = dict(chrom=a.chrom, truth_families=len(truth), pred_clusters=len(pred),
                   truth_families_touched=found, mean_F=round(float(np.mean(fs)), 4),
                   mean_sens=round(float(np.mean([v['sens'] for v in per.values()])), 4),
                   mean_prec=round(float(np.mean([v['prec'] for v in per.values()])), 4),
                   exact_recoveries=exact)
    print(f"{a.chrom}: truth families {len(truth)} | predicted clusters {len(pred)} | "
          f"touched {found} | mean F {summary['mean_F']} "
          f"(sens {summary['mean_sens']} / prec {summary['mean_prec']}) | exact {exact}")
    if a.json:
        with open(a.json, 'w') as fh:
            json.dump(dict(summary=summary, per_family=per), fh, indent=1, sort_keys=True)


if __name__ == '__main__':
    main()

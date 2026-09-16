#!/usr/bin/env python3
"""Pre-registered verdict for the held-out --gtf-refine evaluation (docs/PREREG_gtf_refine_chr17_2026-09-16.md).

usage: gtf_refine_verdict.py BASE_TMAP BASE_STATS BASE_JUNCTIONS BUNDLE_TMAP BUNDLE_STATS BUNDLE_JUNCTIONS \
           BASE_CLASSIFICATION ABL_TSS_CLASSIFICATION
       gtf_refine_verdict.py --self-test
"""
import csv
import json
import re
import sys


def eq_refs(tmap):
    with open(tmap) as fh:
        return {r['ref_id'] for r in csv.DictReader(fh, delimiter='\t') if r['class_code'] == '='}


def precisions(stats):
    txt = open(stats).read()
    tr = float(re.search(r'Transcript level:\s+[\d.]+\s+\|\s+([\d.]+)', txt).group(1))
    ic = float(re.search(r'Intron chain level:\s+[\d.]+\s+\|\s+([\d.]+)', txt).group(1))
    return tr, ic


def novel_canonical_junctions(junctions):
    out = set()
    with open(junctions) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['junction_category'] == 'novel' and r['canonical'] == 'canonical':
                out.add((r['chrom'], r['strand'], r['genomic_start_coord'], r['genomic_end_coord']))
    return out


def verdict(base_eq, base_pr, base_nj, bun_eq, bun_pr, bun_nj):
    e1 = len(bun_eq) >= len(base_eq)
    e2 = bun_pr[0] > base_pr[0] and bun_pr[1] > base_pr[1]
    lost = base_eq - bun_eq
    e3 = len(lost) <= 0.01 * len(base_eq)
    e4 = len(bun_nj) >= len(base_nj)
    if e1 and e2 and e3 and e4:
        v = 'SUPPORTED'
    elif e3 and e4 and (e1 != e2):
        v = 'PARTIAL'
    else:
        v = 'REFUTED'
    return dict(E1_matches=e1, E2_precision=e2, E3_collateral=e3, E4_novel_junctions=e4, verdict=v,
                n_base_eq=len(base_eq), n_bundle_eq=len(bun_eq), n_lost=len(lost), lost_refs=sorted(lost),
                base_tx_ic_pr=base_pr, bundle_tx_ic_pr=bun_pr,
                n_base_novel_junctions=len(base_nj), n_bundle_novel_junctions=len(bun_nj))


def tss_metrics(classification):
    """E5 inputs over SQANTI3 multi-exon full-splice matches: n, p = fraction |diff_to_TSS| <= 50 (unrounded),
    g = count |diff_to_gene_TSS| <= 50. Rows with a non-numeric diff count in n but not in p/g."""
    n = within = gene = 0
    with open(classification) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['structural_category'] != 'full-splice_match' or r['subcategory'] == 'mono-exon':
                continue
            n += 1
            try:
                within += abs(float(r['diff_to_TSS'])) <= 50
            except ValueError:
                pass
            try:
                gene += abs(float(r['diff_to_gene_TSS'])) <= 50
            except ValueError:
                pass
    return dict(n=n, within=within, p=within / n if n else 0.0, g=gene)


def tss_verdict(base, abl):
    e5 = abl['p'] > base['p'] and abl['g'] >= base['g']
    return dict(E5_tss=e5, tss_verdict='SUPPORTED' if e5 else 'REFUTED', baseline_tss=base, abl_tss=abl)


def self_test():
    tb = dict(n=100, within=40, p=0.40, g=60)
    assert tss_verdict(tb, dict(n=100, within=45, p=0.45, g=60))['tss_verdict'] == 'SUPPORTED'
    assert tss_verdict(tb, dict(n=100, within=40, p=0.40, g=70))['tss_verdict'] == 'REFUTED'  # p tie fails
    assert tss_verdict(tb, dict(n=100, within=45, p=0.45, g=59))['tss_verdict'] == 'REFUTED'  # guard fails
    b = set(range(100))
    nj = {('c', '+', '1', '2')}
    assert verdict(b, (35.0, 44.0), nj, b | {100}, (40.0, 48.0), nj)['verdict'] == 'SUPPORTED'
    assert verdict(b, (35.0, 44.0), nj, b, (35.0, 48.0), nj)['verdict'] == 'PARTIAL'  # E2 fails on tx Pr tie
    assert verdict(b, (35.0, 44.0), nj, b, (40.0, 44.0), nj)['verdict'] == 'PARTIAL'  # E2 fails on IC Pr tie
    assert verdict(b, (35.0, 44.0), nj, b - {0}, (40.0, 48.0), nj)['verdict'] == 'PARTIAL'  # E1 fails, 1 lost <= 1%
    assert verdict(b, (35.0, 44.0), nj, b - {0, 1}, (40.0, 48.0), nj)['verdict'] == 'REFUTED'  # E3: 2 lost > 1
    assert verdict(b, (35.0, 44.0), nj, b, (40.0, 48.0), set())['verdict'] == 'REFUTED'  # E4 fails
    assert verdict(b, (35.0, 44.0), nj, b - {0}, (34.0, 44.0), nj)['verdict'] == 'REFUTED'  # neither E1 nor E2
    print('self-test OK')


if __name__ == '__main__':
    if sys.argv[1:] == ['--self-test']:
        self_test()
        sys.exit(0)
    if len(sys.argv) != 9:
        sys.exit(__doc__)
    bt, bs, bj, nt, ns, nj_, bc, tc = sys.argv[1:9]
    res = verdict(eq_refs(bt), precisions(bs), novel_canonical_junctions(bj),
                  eq_refs(nt), precisions(ns), novel_canonical_junctions(nj_))
    res.update(tss_verdict(tss_metrics(bc), tss_metrics(tc)))
    print(json.dumps(res, indent=1))

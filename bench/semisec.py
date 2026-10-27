"""SEMISEC (docs/LEDGER_O1.md): score the candidates of the Rust `semisec` binary (steps 1-4) against the annotation.

  semisec --bam BAM --genome FA --out OUT [--ratio 0.8]     (src/semisec.rs; ~1.5 min testis, 1.6 GB)
  semisec.py report OUT GFF                                  -> OUT/report.md, OUT/silent.json

The annotation is read only here (truth columns, recall, bars). Steps 1-4 were first written in Python in this file (commit
0a12c3d9, held-out run a2e27ff3); the Rust port gives the same candidates on the testis (all 6,630: reads, primaries, class, best).
"""
import bisect
import collections
import csv
import json
import re
import sys

import numpy as np

MIN_READS = 2


class Index:
    """Find intervals overlapping a query (sorted starts + max length)."""

    def __init__(self, ivs):
        self.v = sorted(ivs)
        self.s = [x[0] for x in self.v]
        self.m = max((x[1] - x[0] for x in self.v), default=0)

    def over(self, s, e):
        i = bisect.bisect_left(self.s, s - self.m)
        out = []
        while i < len(self.v) and self.v[i][0] < e:
            if self.v[i][1] > s:
                out.append(self.v[i])
            i += 1
        return out


class PrimCount:
    """Primary records overlapping [s, e), from the binary's prim_iv.bin (u32 start, end pairs sorted per contig)."""

    def __init__(self, out):
        a = np.fromfile(f'{out}/prim_iv.bin', dtype='<u4').reshape(-1, 2)
        self.c = {}
        for line in open(f'{out}/prim_iv.tsv'):
            chrom, off, n = line.split()
            v = a[int(off):int(off) + int(n)]
            self.c[chrom] = (v[:, 0].astype(np.int64), v[:, 1].astype(np.int64), int((v[:, 1] - v[:, 0]).max(initial=0)))

    def count(self, chrom, s, e):
        if chrom not in self.c:
            return 0
        st, en, m = self.c[chrom]
        i, j = np.searchsorted(st, s - m), np.searchsorted(st, e)
        return int((en[i:j] > s).sum())


def load(out):
    cands = []
    for r in csv.DictReader(open(f'{out}/cands.tsv'), delimiter='\t'):
        c = {'chrom': r['chrom'], 'strand': r['strand'], 'start': int(r['start']), 'end': int(r['end']), 'reads': int(r['reads']),
             'primaries': int(r['primaries']), 'silent': r['silent'] == 'true', 'best': r['best'], 'best_seed': r['best_seed']}
        if c['silent']:
            c.update({'class': r['class'], 'orf': r['orf'], 'orf_cov': float(r['orf_cov']), 'processed': r['processed'] == 'true'})
        cands.append(c)
    seeds = {}
    for line in open(f'{out}/seeds.tsv'):
        c, st, s, e, n = line.split()
        seeds[f'{c}:{s}-{e}{st}'] = (c, st, int(s), int(e), int(n))
    return cands, seeds


ISOFORM = re.compile(r' isoform \S+$')


def products(gff):
    """Protein names of coding genes: each protein's CDS span with its `product=` (isoform suffix dropped), by contig."""
    span = {}
    for line in open(gff):
        f = line.split('\t')
        if len(f) < 9 or f[2] != 'CDS':
            continue
        a = dict(x.split('=', 1) for x in f[8].strip().split(';') if '=' in x)
        k = (f[0], a.get('protein_id', a.get('ID', '')))
        s, e = int(f[3]) - 1, int(f[4])
        name = ISOFORM.sub('', a.get('product', '').replace('%2C', ','))
        old = span.get(k)
        span[k] = (min(s, old[0]), max(e, old[1]), name) if old else (s, e, name)
    by = collections.defaultdict(list)
    for (c, _), v in span.items():
        by[c].append(v)
    return {c: Index(v) for c, v in by.items()}


def biotypes(gff):
    genes = collections.defaultdict(list)
    for line in open(gff):
        if line[0] == '#':
            continue
        f = line.split('\t')
        if len(f) < 9 or f[2] not in ('gene', 'pseudogene'):
            continue
        a = dict(x.split('=', 1) for x in f[8].strip().split(';') if '=' in x)
        bt = a.get('gene_biotype', 'other')
        bt = 'pseudogene' if 'pseudogene' in bt else bt if bt in ('protein_coding', 'lncRNA') else 'other'
        genes[f[0]].append((int(f[3]) - 1, int(f[4]), bt, a.get('description', '').replace('%2C', ','), a.get('ID', '')))
    return {c: Index(v) for c, v in genes.items()}


def report(out, gff):
    cands, seeds = load(out)
    pc = PrimCount(out)
    gi = biotypes(gff)

    def top(c, s, e):
        hs = gi[c].over(s, e) if c in gi else []
        if not hs:
            return ('none', '', '')
        h = max(hs, key=lambda h: min(e, h[1]) - max(s, h[0]))
        return h[2], h[3], h[4]

    sil = [c for c in cands if c['silent']]
    for c in sil:
        c['truth'], c['truth_desc'], c['truth_id'] = top(c['chrom'], c['start'], c['end'])
    tab = collections.Counter((c['class'], c['truth']) for c in sil)
    classes = sorted({c['class'] for c in sil})
    bts = ['pseudogene', 'protein_coding', 'lncRNA', 'other', 'none']
    lines = [f"seeds {len(seeds)}; candidates {len(cands)} (silent {len(sil)})", '',
             '| class | n | ' + ' | '.join(bts) + ' |', '|---|---|' + '---|' * len(bts)]
    for k in classes:
        n = sum(tab[(k, b)] for b in bts)
        lines.append(f'| {k} | {n} | ' + ' | '.join(f'{tab[(k, b)]} ({tab[(k, b)] / n:.2f})' for b in bts) + ' |')
    # recall: annotated silent pseudogenes named "X-like" after an expressed seed's protein-coding gene
    seed_desc = set()
    for sid, (c, st, s, e, n) in seeds.items():
        for h in (gi[c].over(s, e) if c in gi else []):
            if h[2] == 'protein_coding' and h[3]:
                seed_desc.add(h[3])
    cidx = collections.defaultdict(list)
    for c in sil:
        cidx[c['chrom']].append((c['start'], c['end']))
    cidx = {k: Index(v) for k, v in cidx.items()}
    found = tot = 0
    for chrom, idx in gi.items():
        for s, e, bt, desc, gid in idx.v:
            if bt != 'pseudogene' or not desc.endswith('-like') or desc[:-5] not in seed_desc:
                continue
            if pc.count(chrom, s, e) >= MIN_READS:
                continue  # expressed copy: not silent, not in the recall denominator
            tot += 1
            found += bool(cidx.get(chrom) and cidx[chrom].over(s, e))
    lines += ['', f'recall: annotated silent pseudogenes "X-like" of an expressed seed gene, overlapped by a silent candidate: '
                  f'{found} / {tot} ({found / max(tot, 1):.3f})']
    # stage-2 bars (LEDGER_O1 "SEMISEC stage-2 bars"): fractions among annotated candidates (truth != none), pooled over PROCESSED/UNPROCESSED
    def frac(bt, orf):
        xs = [c for c in sil if c['truth'] != 'none' and c['class'].endswith('/' + orf)]
        return sum(c['truth'] == bt for c in xs) / max(len(xs), 1), len(xs)
    pd, nd = frac('pseudogene', 'DISRUPTED')
    pi, ni = frac('pseudogene', 'INTACT')
    pdd = [c for c in sil if c['class'] == 'PROCESSED/DISRUPTED']
    ppd = sum(c['truth'] == 'pseudogene' for c in pdd) / max(len(pdd), 1)
    ll, nl = frac('lncRNA', 'LNC-LIKE')
    lc = [c for c in sil if c['truth'] != 'none' and c['class'].split('/')[1] in ('DISRUPTED', 'INTACT')]
    lcf = sum(c['truth'] == 'lncRNA' for c in lc) / max(len(lc), 1)
    rec = found / max(tot, 1)
    bars = {'B1 pseudogene DISRUPTED/INTACT >= 1.5': (pd / max(pi, 1e-9), pd / max(pi, 1e-9) >= 1.5),
            'B2 PROCESSED/DISRUPTED pseudogene >= 0.40': (ppd, ppd >= 0.40),
            'B3 recall >= 0.10': (rec, rec >= 0.10),
            'B4 (separate) lncRNA LNC-LIKE/coding >= 1.25': (ll / max(lcf, 1e-9), ll / max(lcf, 1e-9) >= 1.25)}
    lines += ['', f'annotated-only: pseudogene DISRUPTED {pd:.3f} (n {nd}) INTACT {pi:.3f} (n {ni}); '
                  f'lncRNA LNC-LIKE {ll:.3f} (n {nl}) coding {lcf:.3f} (n {len(lc)})']
    lines += [f'{k}: {v:.3f} {"PASS" if ok else "FAIL"}' for k, (v, ok) in bars.items()]
    lines.append('VERDICT (B1-B3): ' + ('PASS' if all(ok for k, (_, ok) in bars.items() if not k.startswith('B4')) else 'FAIL'))
    # parent agreement (precision on annotated pseudogene hits): a silent candidate whose annotated pseudogene is "X-like" agrees when
    # the best read's seed locus overlaps a protein named X (CDS product; NCBI names "X-like" pseudogenes after the protein, while a
    # gene's description is the gene name: "ribosomal protein S26" makes "small ribosomal subunit protein eS26")
    pi = products(gff)
    agree = n_like = 0
    for c in sil:
        if c['truth'] != 'pseudogene' or not c['truth_desc'].endswith('-like'):
            continue
        ch, iv = c['best_seed'].rsplit(':', 1)
        s0, e0 = map(int, iv.split('-'))
        n_like += 1
        agree += any(h[2] == c['truth_desc'][:-5] for h in (pi[ch].over(s0, e0) if ch in pi else []))
    lines.append(f'parent agreement: {agree} / {n_like} ({agree / max(n_like, 1):.3f}) silent candidates on an "X-like" pseudogene found from a seed on X')
    none = [c for c in sil if c['truth'] == 'none']
    sizes = sorted(c['end'] - c['start'] for c in none)
    lines.append(f"unannotated silent candidates: {len(none)}; median span {sizes[len(sizes) // 2] if sizes else 0} bp")
    open(f'{out}/report.md', 'w').write('\n'.join(lines) + '\n')
    json.dump(sil, open(f'{out}/silent.json', 'w'))
    print('\n'.join(lines))


if __name__ == '__main__':
    report(sys.argv[2], sys.argv[3]) if sys.argv[1] == 'report' else sys.exit(__doc__)

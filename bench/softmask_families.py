#!/usr/bin/env python3
"""PREREG r1690: soft-mask (lowercase) content of every family copy, and the families each way would set aside.

  softmask_families.py COPIES_TSV FASTA OUT_PREFIX [--censat BED] [--rmsk BED] [--assignments TSV] [--filter COPIES_TSV]...

OUT_PREFIX.copies.tsv   per copy: exon / span lowercase fraction, uppercase exonic bp, censat class, RepeatMasker exon class
OUT_PREFIX.families.tsv per family: label, the ways it falls under (majority of copies), assigned rows
OUT_PREFIX.WAY.txt      family ids under the way; --filter T writes OUT_PREFIX.WAY.<basename T> without those families
Ways: W1 exon lowercase >= 0.9; W1b >= 0.5; W2 span lowercase >= 0.9; W3 < 100 uppercase exonic bp."""
import sys, os, csv, collections, numpy as np, pysam
a = sys.argv[1:]; copies_p, fa_p, out = a[:3]
opt = lambda k: a[a.index(k) + 1] if k in a else None
SAT = {'hor', 'hsat1A', 'hsat1B', 'hsat2', 'hsat3', 'mon', 'bsat', 'gsat', 'dhor', 'censat'}
WAYS = {'W1': lambda c: c['exon_low'] >= 0.9, 'W1b': lambda c: c['exon_low'] >= 0.5, 'W2': lambda c: c['span_low'] >= 0.9, 'W3': lambda c: c['upper_exon_bp'] < 100}

def bed(path, col, prefix):
    d = collections.defaultdict(list); names = ['-']
    for l in open(path):
        f = l.rstrip('\n').split('\t')
        if len(f) > col and f[1].isdigit():
            c = f[col].split('_')[0] if prefix else f[col]
            if c not in names: names.append(c)
            d[f[0]].append((int(f[1]), int(f[2]), names.index(c)))
    return d, names
def paint(iv, n):
    x = np.zeros(n, np.uint8)
    for s, e, c in iv: x[s:e] = c
    return x
def top(x, names):  # commonest non-background class of a code slice and its share of the slice
    b = np.bincount(x, minlength=len(names)); b[0] = 0; i = int(b.argmax())
    return (names[i], b[i] / len(x)) if len(x) and b[i] else ('-', 0.0)

cs, csn = bed(opt('--censat'), 3, True) if opt('--censat') else ({}, ['-'])
rm, rmn = bed(opt('--rmsk'), 6, False) if opt('--rmsk') else ({}, ['-'])
rows = list(csv.DictReader(open(copies_p), delimiter='\t')); fa = pysam.FastaFile(fa_p); res = []
for chrom in sorted({r['chrom'] for r in rows}):
    s = np.frombuffer(fa.fetch(chrom).encode(), dtype=np.uint8); low = s >= 97
    C = paint(cs.get(chrom, ()), len(s)); R = paint(rm.get(chrom, ()), len(s))
    for r in rows:
        if r['chrom'] != chrom: continue
        ex = [tuple(map(int, e.split('-'))) for e in r['exons'].split(',')]  # exons as in the table (1-based start, inclusive end)
        idx = np.concatenate([np.arange(x - 1, y) for x, y in ex]); lo, hi = int(r['start']) - 1, int(r['end'])
        cc, cf = top(C[lo:hi], csn); rc, rf = top(R[idx], rmn)
        res.append(dict(family=r['family_id'], copy=r['copy_idx'], chrom=chrom, start=r['start'], end=r['end'], exon_bp=len(idx),
                        exon_low=float(low[idx].mean()), span_low=float(low[lo:hi].mean()), upper_exon_bp=int((~low[idx]).sum()),
                        censat=cc if cf >= 0.5 else '-', rmsk=rc if rf >= 0.5 else '-', rmsk_top=rc, rmsk_top_frac=rf))
with open(out + '.copies.tsv', 'w') as w:
    k = list(res[0]); w.write('\t'.join(k) + '\n')
    for c in res: w.write('\t'.join(f'{c[x]:.3f}' if isinstance(c[x], float) else str(c[x]) for x in k) + '\n')

fam = collections.defaultdict(list)
for c in res: fam[c['family']].append(c)
assigned = collections.Counter()
if opt('--assignments'):
    for r in csv.DictReader(open(opt('--assignments')), delimiter='\t'): assigned[r['family_id']] += r['status'] == 'assigned'
most = lambda cs_, f: sum(map(f, cs_)) * 2 > len(cs_)
def label(v):
    if most(v, lambda c: c['censat'] in SAT or c['rmsk'] == 'Satellite'): return 'satellite'
    if most(v, lambda c: c['censat'] == 'rDNA' or c['rmsk'] == 'rRNA'): return 'rDNA'
    if most(v, lambda c: c['censat'] == 'ct'): return 'ct'
    return 'none'
F = {f: dict(n=len(v), label=label(v), rmsk=collections.Counter(c['rmsk'] for c in v).most_common(1)[0][0], assigned=assigned[f],
             **{w: int(most(v, fn)) for w, fn in WAYS.items()}) for f, v in fam.items()}
with open(out + '.families.tsv', 'w') as w:
    w.write('family\tn_copies\tlabel\trmsk_majority\tassigned\t' + '\t'.join(WAYS) + '\n')
    for f, d in F.items(): w.write(f"{f}\t{d['n']}\t{d['label']}\t{d['rmsk']}\t{d['assigned']}\t" + '\t'.join(str(d[x]) for x in WAYS) + '\n')
labels = ['satellite', 'rDNA', 'ct', 'none']; tot = {l: [d for d in F.values() if d['label'] == l] for l in labels}; A = sum(assigned.values())
print(f"# {copies_p}: {len(F)} families, {len(res)} copies; assigned rows {A}")
print('way\t' + '\t'.join(f'{l} fam (of {len(tot[l])})\t{l} assigned' for l in labels) + '\tnone & no repeat-majority exon class (unexplained)')
for wname in WAYS:
    hit = {f for f, d in F.items() if d[wname]}; open(f'{out}.{wname}.txt', 'w').write(''.join(f + '\n' for f in sorted(hit)))
    print(wname + '\t' + '\t'.join(f"{sum(d[wname] for d in tot[l])} ({sum(d[wname] for d in tot[l]) / max(1, len(tot[l])):.3f})\t{sum(d['assigned'] for d in tot[l] if d[wname])}" for l in labels)
          + f"\t{sum(1 for d in tot['none'] if d[wname] and d['rmsk'] == '-')}")
    for t in [a[i + 1] for i, x in enumerate(a) if x == '--filter']:
        with open(f'{out}.{wname}.{os.path.basename(t)}', 'w') as w:
            for i, l in enumerate(open(t)):
                if i == 0 or l.split('\t', 1)[0] not in hit: w.write(l)

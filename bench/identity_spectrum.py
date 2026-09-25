#!/usr/bin/env python3
"""Identity spectrum (docs/PREREG_identity_spectrum_2026-09-24.md): which alignment tier of the family-edge builder
recovers which Ensembl Compara paralogue pairs, band by band of identity, on one chromosome.

usage: identity_spectrum.py --gtf ASSEMBLED.gtf --ref REF.gtf --fasta GENOME.fa --chrom chr16 --compara compara_chr16.tsv --out PREFIX
       [--mmseqs /path/to/mmseqs] [--threads 4]

Nodes: one spliced representative per expressed locus (gene_id group; most reads, tie longer). Tiers, all-vs-all on
those nodes with the builder's flags: T1 asm20 (k19) identity>=0.80 cov>=0.50; T2 asm20 -k11 -w5 identity>=0.60
cov>=0.50; T3 mmseqs translated (--search-type 2) protein identity>=0.30, e<=1e-5, qcov>=0.50. Loci map to HGNC
symbols by exon overlap with the RefSeq annotation. Truth: Compara pairs (BioMart, columns gene, paralog, perc_id,
perc_id_r1, subtype, paralog_chromosome) with both genes on the chromosome and both expressed."""
import argparse, collections, itertools, os, re, subprocess, sys
import pysam

ap = argparse.ArgumentParser()
ap.add_argument('--gtf', required=True); ap.add_argument('--ref', required=True); ap.add_argument('--fasta', required=True)
ap.add_argument('--chrom', required=True); ap.add_argument('--compara', required=True); ap.add_argument('--out', required=True)
ap.add_argument('--mmseqs', default='mmseqs'); ap.add_argument('--threads', type=int, default=4)
ap.add_argument('--catalog', help='score a gw_family_catalog copies.tsv at the PAIR level against the truth instead of building tiers')
ap.add_argument('--referee', help='truth as a gene->family table (Gene Name, Family ID) instead of Compara; pairs = same family')
ap.add_argument('--gff-genes', help='gene spans from a GFF (gene features, Name=) instead of the RefSeq GTF, for the catalog mapping')
ap.add_argument('--universe', help='catalog mode: a .truth_pairs.tsv from a tier run; recall is ALSO reported over its pairs (both genes expressed), a denominator the catalog cannot move')
a = ap.parse_args()
def attr(s, k):
    m = re.search(k + r' "([^"]+)"', s); return m.group(1) if m else None

if a.catalog:
    # ---- gene spans
    genes = {}
    if a.gff_genes:
        for ln in open(a.gff_genes):
            f = ln.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] != a.chrom or f[2] not in ('gene', 'pseudogene'): continue
            m = re.search(r'(?:^|;)Name=([^;]+)', f[8])
            if m: genes[m.group(1)] = (int(f[3]) - 1, int(f[4]))
    else:
        ex = collections.defaultdict(list)
        for ln in open(a.ref):
            f = ln.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] != a.chrom or f[2] != 'exon': continue
            g = attr(f[8], 'gene_name') or (attr(f[8], 'gene_id') or '').replace('gene-', '', 1)
            if not g: continue  # exon lines without a gene attribute would otherwise form one chromosome-wide span
            ex[g].append((int(f[3]) - 1, int(f[4])))
        genes = {g: (min(s_ for s_, _ in v), max(e for _, e in v)) for g, v in ex.items()}
    glist = sorted((s_, e, g) for g, (s_, e) in genes.items())
    def gene_of(s_, e):
        best = None
        for gs, ge, g in glist:
            if ge <= s_: continue
            if gs >= e: break
            o = min(e, ge) - max(s_, gs)
            if o > 0 and (best is None or o > best[0]): best = (o, g)
        return best[1] if best else None
    # ---- catalog copies -> genes; family pairs
    fam_genes = collections.defaultdict(set); fam_spliced = collections.defaultdict(set); ncop = 0
    for ln in open(a.catalog):
        f = ln.rstrip('\n').split('\t')
        if f[0] == 'family_id' or f[3] != a.chrom: continue
        ncop += 1; g = gene_of(int(f[4]), int(f[5]))
        if g:
            fam_genes[f[0]].add(g)
            if int(f[6]) >= 2: fam_spliced[f[0]].add(g)   # multi-exon copies only (single-exon loci are mostly intronic/unspliced)
    cat_pairs = {frozenset(p) for gs in fam_genes.values() for p in itertools.combinations(sorted(gs), 2)}
    spliced_pairs = {frozenset(p) for gs in fam_spliced.values() for p in itertools.combinations(sorted(gs), 2)}
    cat_genes = set().union(*fam_genes.values()) if fam_genes else set()
    largest = max((len(v) for v in fam_genes.values()), default=0)
    # ---- truth pairs
    if a.referee:
        fam = {}
        for ln in open(a.referee):
            f = ln.rstrip('\n').split('\t')
            if f[0] == 'Gene Name' or len(f) < 2: continue
            fam[f[0]] = f[1]
        truth_all = {frozenset(p) for fid in set(fam.values()) for p in itertools.combinations(sorted(g for g, v in fam.items() if v == fid), 2)}
        judge = set(fam)
        band_of = lambda k: 'all'
    else:
        chrom_num = a.chrom.replace('chr', ''); compara = {}; judge = set()
        for ln in open(a.compara):
            f = ln.rstrip('\n').split('\t')
            if len(f) < 6 or not f[0] or not f[1]: continue
            judge.add(f[0])
            if f[5] != chrom_num or f[0] == f[1]: continue
            try: pid = max(float(f[2] or 0), float(f[3] or 0))
            except ValueError: continue
            k = frozenset((f[0], f[1]))
            if k not in compara or pid > compara[k][0]: compara[k] = (pid, f[4])
        truth_all = set(compara)
        BANDS = [(90, 101, '>=90'), (80, 90, '80-90'), (70, 80, '70-80'), (60, 70, '60-70'), (50, 60, '50-60'), (30, 50, '30-50'), (0, 30, '<30')]
        def band_of(k):
            p = compara[k][0]
            for lo, hi, n in BANDS:
                if lo <= p < hi: return n
            return '<30'
    # recall over truth pairs whose BOTH genes are in the catalog (present as copies) — the catalog cannot join what it did not emit;
    # and, for context, over truth pairs with both genes mapped by any catalog copy or not
    truth_in = {k for k in truth_all if k <= cat_genes}
    print(f'[catalog] {a.catalog}: {ncop} copies on {a.chrom}, {len(fam_genes)} families, largest {largest} genes, {len(cat_pairs)} gene pairs; truth pairs {len(truth_all)}, with both genes in the catalog {len(truth_in)}')
    by_band = collections.defaultdict(lambda: [0, 0])
    for k in truth_in:
        b = band_of(k); by_band[b][1] += 1; by_band[b][0] += (k in cat_pairs)
    for b, (hit, n) in sorted(by_band.items(), key=lambda x: -x[1][1]):
        print(f'   recall {b:6s} {hit}/{n} = {hit/n:.3f}')
    judgeable = {k for k in cat_pairs if k <= judge}
    tp = sum(1 for k in judgeable if k in truth_all)
    print(f'   precision (judgeable catalog pairs): {tp}/{len(judgeable)} = {tp/len(judgeable) if judgeable else float("nan"):.3f}')
    js = {k for k in spliced_pairs if k <= judge}; tps = sum(1 for k in js if k in truth_all)
    print(f'   precision, multi-exon copies only: {tps}/{len(js)} = {tps/len(js) if js else float("nan"):.3f}')
    if a.universe:
        uni = set()
        for ln in open(a.universe):
            f = ln.rstrip('\n').split('\t')
            if f[0] == 'geneA' or len(f) < 2: continue
            uni.add(frozenset((f[0], f[1])))
        uni &= truth_all
        ub = collections.defaultdict(lambda: [0, 0])
        for k in uni:
            b = band_of(k); ub[b][1] += 1; ub[b][0] += (k in cat_pairs)
        print(f'   recall over the fixed universe ({len(uni)} truth pairs, both genes expressed):')
        for b, (hit, n) in sorted(ub.items(), key=lambda x: -x[1][1]):
            print(f'      {b:6s} {hit}/{n} = {hit/n:.3f}')
    sys.exit(0)

fa = pysam.FastaFile(a.fasta)

# ---- nodes: spliced representative per locus
tx = collections.defaultdict(list); reads = {}; gene_of = {}; strand = {}
for ln in open(a.gtf):
    f = ln.rstrip('\n').split('\t')
    if len(f) < 9 or f[0] != a.chrom: continue
    t = attr(f[8], 'transcript_id')
    if f[2] == 'transcript':
        gene_of[t] = attr(f[8], 'gene_id') or t; strand[t] = f[6]; r = attr(f[8], 'reads'); reads[t] = int(r) if r else 0
    elif f[2] == 'exon': tx[t].append((int(f[3]) - 1, int(f[4])))
loci = collections.defaultdict(list)
for t, g in gene_of.items(): loci[g].append(t)
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')
node = {}   # locus -> (span, exons, seq)
with open(a.out + '.nodes.fa', 'w') as fh:
    for g, ts in loci.items():
        ts = [t for t in ts if tx.get(t)]
        if not ts: continue
        rep = max(ts, key=lambda t: (reads.get(t, 0), max(b for _, b in tx[t]) - min(s for s, _ in tx[t])))
        ex = sorted(tx[rep]); seq = ''.join(fa.fetch(a.chrom, s, e) for s, e in ex).upper()
        if strand[rep] == '-': seq = seq.translate(COMP)[::-1]
        if len(seq) < 200: continue
        node[g] = ((ex[0][0], ex[-1][1]), ex, seq); fh.write(f'>{g}\n{seq}\n')
print(f'[spectrum] {len(node)} expressed loci with a representative >= 200 bp', flush=True)

# ---- locus -> gene symbol (RefSeq exon overlap; the symbol is the gene_id with its "gene-" prefix stripped)
ref_ex = collections.defaultdict(list)
for ln in open(a.ref):
    f = ln.rstrip('\n').split('\t')
    if len(f) < 9 or f[0] != a.chrom or f[2] != 'exon': continue
    g = attr(f[8], 'gene_id') or ''
    ref_ex[g.replace('gene-', '', 1)].append((int(f[3]) - 1, int(f[4])))
ref_list = sorted((min(s for s, _ in v), max(e for _, e in v), g, v) for g, v in ref_ex.items())
def symbol_of(exons):
    lo, hi = exons[0][0], exons[-1][1]; best = None
    for s, e, g, v in ref_list:
        if e <= lo: continue
        if s >= hi: break
        o = sum(max(0, min(b, y) - max(a_, x)) for a_, b in exons for x, y in v)
        if o > 0 and (best is None or o > best[0]): best = (o, g)
    return best[1] if best else None
sym = {g: symbol_of(v[1]) for g, v in node.items()}
sym = {g: s for g, s in sym.items() if s}
by_sym = collections.defaultdict(list)
for g, s in sym.items(): by_sym[s].append(g)
print(f'[spectrum] {len(sym)} loci map to {len(by_sym)} RefSeq symbols', flush=True)

# ---- truth: Compara pairs on this chromosome
chrom_num = a.chrom.replace('chr', '')
compara = {}   # frozenset(symbols) -> (max perc_id, subtype)
genes_with_data = set()
for ln in open(a.compara):
    f = ln.rstrip('\n').split('\t')
    if len(f) < 6 or not f[0] or not f[1]: continue
    genes_with_data.add(f[0])
    if f[5] != chrom_num or f[0] == f[1]: continue
    try: pid = max(float(f[2] or 0), float(f[3] or 0))
    except ValueError: continue
    k = frozenset((f[0], f[1]))
    if k not in compara or pid > compara[k][0]: compara[k] = (pid, f[4])
expressed = set(by_sym)
truth = {k: v for k, v in compara.items() if k <= expressed}
unrecoverable = sum(1 for k in compara if not k <= expressed)
print(f'[spectrum] Compara pairs on {a.chrom}: {len(compara)}; both genes expressed: {len(truth)}; with an unexpressed member: {unrecoverable}', flush=True)

# ---- tiers
UNION = {}   # pair -> union-of-records coverage of the shorter sequence (post-hoc variant, addendum 1)
def mm2(flags, out):
    subprocess.run(f"minimap2 {flags} -c -X --no-long-join -N 50 -p 0.1 --secondary=yes -t {a.threads} {a.out}.nodes.fa {a.out}.nodes.fa > {out} 2>/dev/null", shell=True, check=True)
    best = {}; spans = collections.defaultdict(list)
    for l in open(out):
        f = l.split('\t'); q, t = f[0], f[5]
        if q == t: continue
        idn = int(f[9]) / int(f[10]); cov = (int(f[3]) - int(f[2])) / min(int(f[1]), int(f[6]))
        k = frozenset((q, t))
        if k not in best or (idn, cov) > best[k]: best[k] = (idn, cov)
        # union coverage on the SHORTER sequence's coordinates
        if int(f[1]) <= int(f[6]): spans[k].append((int(f[2]), int(f[3]), int(f[1])))
        else: spans[k].append((int(f[7]), int(f[8]), int(f[6])))
    for k, v in spans.items():
        v.sort(); cov = 0; cur = None
        for s_, e_, L in v:
            if cur is None or s_ > cur[1]:
                if cur: cov += cur[1] - cur[0]
                cur = [s_, e_]
            else: cur[1] = max(cur[1], e_)
        cov += cur[1] - cur[0]
        UNION[(out, k)] = cov / v[0][2]
    return best
t1 = mm2('-x asm20', a.out + '.t1.paf'); t2 = mm2('-x asm20 -k11 -w5', a.out + '.t2.paf')
subprocess.run(f"{a.mmseqs} easy-search {a.out}.nodes.fa {a.out}.nodes.fa {a.out}.t3.m8 {a.out}.tmp --search-type 2 --threads {a.threads} -e 1e-5 --format-output query,target,pident,qcov,tcov,evalue > /dev/null 2>&1", shell=True, check=True)
t3 = {}
for l in open(a.out + '.t3.m8'):
    f = l.split('\t')
    if f[0] == f[1]: continue
    idn = float(f[2]) / (100 if float(f[2]) > 1 else 1); cov = max(float(f[3]), float(f[4]))
    k = frozenset((f[0], f[1]))
    if k not in t3 or (idn, cov) > t3[k]: t3[k] = (idn, cov)
def edge(best, k, floor): v = best.get(k); return v is not None and v[0] >= floor and v[1] >= 0.50
def locus_pairs(sa, sb):
    return [frozenset((x, y)) for x in by_sym[sa] for y in by_sym[sb] if x != y]
def tier_hit(sa, sb, tiers):
    return any(edge(b, k, fl) for b, fl in tiers for k in locus_pairs(sa, sb))
TIERS = {'T1': [(t1, 0.80)], 'T1+T2': [(t1, 0.80), (t2, 0.60)], 'T1+T2+T3': [(t1, 0.80), (t2, 0.60), (t3, 0.30)]}
BANDS = [(90, 101, '>=90'), (80, 90, '80-90'), (70, 80, '70-80'), (60, 70, '60-70'), (50, 60, '50-60'), (30, 50, '30-50'), (0, 30, '<30')]
def band(p):
    for lo, hi, n in BANDS:
        if lo <= p < hi: return n
    return '<30'

# ---- recall by Compara band
rows = []; print('\n== RECALL of Compara paralogue pairs (both genes expressed), by Compara protein identity band')
print(f"{'band':8s} {'n_pairs':>7s}  " + '  '.join(f'{t:>9s}' for t in TIERS))
for lo, hi, n in BANDS:
    ks = [k for k, v in truth.items() if lo <= v[0] < hi]
    if not ks: continue
    rec = {t: sum(1 for k in ks if tier_hit(*sorted(k), tiers)) for t, tiers in TIERS.items()}
    print(f"{n:8s} {len(ks):7d}  " + '  '.join(f"{rec[t]/len(ks):9.3f}" for t in TIERS)); rows.append(('recall', n, len(ks), {t: rec[t] / len(ks) for t in TIERS}))
# ---- precision by OUR identity band: aligned symbol pairs (both genes with Compara data) that are Compara paralogues
print('\n== PRECISION of aligned pairs (both genes have Compara data), by the tier\'s own identity band')
def sym_pairs(best, floor):
    out = {}
    for k, (idn, cov) in best.items():
        if idn < floor or cov < 0.50: continue
        x, y = tuple(k); sx, sy = sym.get(x), sym.get(y)
        if not sx or not sy or sx == sy: continue
        kk = frozenset((sx, sy))
        if kk not in out or idn > out[kk]: out[kk] = idn
    return out
for name, best, floor in (('T1 (nt)', t1, 0.80), ('T2 (nt)', t2, 0.60), ('T3 (protein)', t3, 0.30)):
    sp = {k: v for k, v in sym_pairs(best, floor).items() if k <= genes_with_data}
    print(f'-- {name}: {len(sp)} judgeable aligned symbol pairs')
    for lo, hi, n in BANDS:
        ks = [k for k, v in sp.items() if lo <= v * 100 < hi]
        if ks: tp = sum(1 for k in ks if k in compara); print(f"   {n:8s} n={len(ks):5d}  precision {tp/len(ks):.3f}"); rows.append((f'precision {name}', n, len(ks), tp / len(ks)))
missed = [(tuple(sorted(k)), v[0]) for k, v in truth.items() if not tier_hit(*sorted(k), TIERS['T1+T2+T3'])]
# diagnosis of every truth pair: the best record each tier has for ANY locus pair of the two symbols (identity, coverage),
# or none — separates "no seed" (no record) from "coverage clause" (record below 0.50) from "identity floor"
def best_any(best, sa, sb):
    recs = [best[k] for k in locus_pairs(sa, sb) if k in best]
    return max(recs) if recs else None
with open(a.out + '.truth_pairs.tsv', 'w') as fh:
    fh.write('geneA\tgeneB\tcompara_pid\tsubtype\trecovered\tT1_idn\tT1_cov\tT2_idn\tT2_cov\tT3_pid\tT3_cov\tlenA\tlenB\n')
    for k, (pid, sub) in sorted(truth.items(), key=lambda x: -x[1][0]):
        sa, sb = sorted(k); rec = tier_hit(sa, sb, TIERS['T1+T2+T3'])
        cells = []
        for best in (t1, t2, t3):
            b = best_any(best, sa, sb); cells += ([f'{b[0]:.3f}', f'{b[1]:.2f}'] if b else ['-', '-'])
        la = max(len(node[g][2]) for g in by_sym[sa]); lb = max(len(node[g][2]) for g in by_sym[sb])
        fh.write('\t'.join([sa, sb, f'{pid:.1f}', sub, str(int(rec))] + cells + [str(la), str(lb)]) + '\n')
why = collections.defaultdict(collections.Counter)
for k, (pid, sub) in truth.items():
    if pid < 60 or tier_hit(*sorted(k), TIERS['T1+T2']): continue
    sa, sb = sorted(k); b = best_any(t2, sa, sb)
    why[band(pid)]['no nucleotide record (seeding)' if b is None else ('record, coverage < 0.50' if b[1] < 0.50 else 'record, identity < 0.60')] += 1
print('== why truth pairs at >= 60% Compara identity are missed by T1+T2, per band:')
for b_ in ('>=90', '80-90', '70-80', '60-70'):
    if why[b_]: print(f'   {b_:6s} {dict(why[b_])}')
# post-hoc variants of the T2 coverage clause (NOT part of the bar): union-of-records coverage, and a 0.30 floor
def edge_union(k, floor, cov_floor):
    v = t2.get(k); return v is not None and v[0] >= floor and UNION.get((a.out + '.t2.paf', k), 0) >= cov_floor
print('== post-hoc: recall of Compara pairs under T2 coverage variants (identity >= 0.60)')
print(f"{'band':8s} {'n':>5s} {'cov>=0.50':>10s} {'cov>=0.30':>10s} {'union>=0.50':>12s} {'union>=0.30':>12s}")
for lo, hi, n in BANDS:
    ks = [k for k, v in truth.items() if lo <= v[0] < hi]
    if not ks: continue
    def rec(fn): return sum(1 for k in ks if any(fn(kk) for kk in locus_pairs(*sorted(k)))) / len(ks)
    print(f"{n:8s} {len(ks):5d} {rec(lambda kk: edge(t2, kk, 0.60)):10.3f} {rec(lambda kk: t2.get(kk) is not None and t2[kk][0] >= 0.60 and t2[kk][1] >= 0.30):10.3f} {rec(lambda kk: edge_union(kk, 0.60, 0.50)):12.3f} {rec(lambda kk: edge_union(kk, 0.60, 0.30)):12.3f}")
# precision of the union variant, judgeable pairs, by identity band
def sym_pairs_union(floor, cov_floor):
    out = {}
    for k, (idn, cov) in t2.items():
        if idn < floor or UNION.get((a.out + '.t2.paf', k), 0) < cov_floor: continue
        x, y = tuple(k); sx, sy = sym.get(x), sym.get(y)
        if not sx or not sy or sx == sy: continue
        kk = frozenset((sx, sy))
        if kk not in out or idn > out[kk]: out[kk] = idn
    return out
for cf in (0.50, 0.30):
    sp = {k: v for k, v in sym_pairs_union(0.60, cf).items() if k <= genes_with_data}
    tp = sum(1 for k in sp if k in compara)
    print(f"   union>={cf}: {len(sp)} judgeable aligned pairs, precision {tp/len(sp) if sp else float('nan'):.3f}")
print(f'\n== truth pairs recovered by NO tier: {len(missed)} of {len(truth)}; Compara identity median {sorted(p for _, p in missed)[len(missed)//2] if missed else "-"}; examples {missed[:8]}')
with open(a.out + '.spectrum.tsv', 'w') as fh:
    fh.write('metric\tband\tn\tvalue\n')
    for m, b, n, v in rows: fh.write(f'{m}\t{b}\t{n}\t{v}\n')
print(f'wrote {a.out}.spectrum.tsv')

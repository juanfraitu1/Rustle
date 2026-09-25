#!/usr/bin/env python3
"""Scorers: family / pair / copy-assignment scoring against the project's truths (wave 7, 2026-09-24).

One argparse, one subcommand per old script (old files at git tag `notebook-2026-09-24`). Every subcommand keeps its
old stdout format and output files, so `cmp` against the old script proves identity.

Old -> new:
  referee_band_score.py CLUSTERS LOCI.gff3 LABEL  (paths relative to /mnt/linuxdisk/tmp/gw22/sec/)
      -> score.py pairs --members $S/CLUSTERS --genes $S/ref/NC_073244.2.genes.gff --chrom NC_073244.2
           --truth families:$S/ref/NC_073244.2.tsv --expressed $S/ref/NC_073244.2.expressed.tsv
           --bands paf:$S/ref/referee_mrna.paf --label LABEL            (the unused LOCI.gff3 is dropped)
  identity_spectrum.py --gtf x --ref REF.gtf --fasta x --chrom C --compara CMP --out O --catalog COPIES
                       [--universe U] [--referee R] [--gff-genes G]
      -> score.py pairs --members COPIES --genes REF.gtf|G --chrom C --truth compara:CMP|families:R [--universe U]
           (the unused --gtf/--fasta/--out are dropped, B8)
  identity_spectrum.py --gtf A.gtf --ref REF.gtf --fasta FA --chrom C --compara CMP --out O [--mmseqs M] [--threads N]
      -> score.py spectrum (same flags; + --estimator/--coverage, defaults = the old nm/bl and query-over-min)
  heldout_family_score.py ...         -> score.py heldout (same flags; + --exact-only = the B4 fix, off by default)
  soto_vs_us_referee.py ...           -> score.py referee (same flags; + --threads, default 4; B1 NameError fixed)
  protein_edge_gap.py ...             -> score.py edge-gap (same flags)
  rna_truth_from_protein.py ...       -> score.py rna-ceiling (same flags)
  member_completeness.py ARM_GTF REF_GTF UNIVERSE_TSV GENES_GFF CHROM LABEL
      -> score.py members --gtf ARM_GTF --ref REF_GTF --universe UNIVERSE_TSV --chrom CHROM --label LABEL
           (the unused GENES_GFF is dropped)
  adjudicated_truth.py score ...      -> score.py adjudicated (same args)
  protein_families.py score ...       -> score.py protein (same args)
  eichler_compare.py ...              -> score.py eichler (same flags)
  CATALOG_TSV=CAT copy_assign_read_truth.py score PREFIX O2PREFIX
      -> score.py reads --catalog CAT PREFIX O2PREFIX   (CATALOG_TSV is still honoured when --catalog is absent)
  copy_assign_tool_bakeoff.py calls ...   -> score.py bakeoff-calls (same args)
  copy_assign_tool_bakeoff.py compare ... -> score.py bakeoff-compare (same args; samtools output streamed, B6)
  locus_reads.py BAM CHROM START END      -> score.py locus-reads BAM CHROM START END
      (library: `from lib import reads_with_block_in, spanning_genes`)

Determinism: `pairs` breaks ties that the old scripts left to Python's per-process string-hash order (the order of
band lines with equal pair counts; the composition of the largest cluster when counts tie) by a fixed rule, so its
output no longer depends on PYTHONHASHSEED. Everything else is the old code path.

Only the standard library is imported at module top; pysam / numpy / scipy are imported inside the subcommands.
"""
import argparse
import bisect
import collections
import csv
import itertools
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib  # noqa: E402


# ================================================================ pairs (referee_band_score + identity_spectrum --catalog)
BANDS_COMPARA = [(90, 101, '>=90'), (80, 90, '80-90'), (70, 80, '70-80'), (60, 70, '60-70'), (50, 60, '50-60'),
                 (30, 50, '30-50'), (0, 30, '<30')]
BANDS_MRNA = [(0.90, 1.01, '>=90'), (0.80, 0.90, '80-90'), (0.70, 0.80, '70-80'), (0.60, 0.70, '60-70'), (0.0, 0.60, '<60')]


def load_gene_spans(path, chrom, fmt):
    """Gene spans (0-based start, end) on `chrom`.
    gff: gene/pseudogene `Name=` (anchored), last record wins (referee_band_score, identity_spectrum --gff-genes).
    gtf: union of the exon lines keyed by gene_name, else gene_id minus `gene-`; exon lines with neither are skipped
         (identity_spectrum's RefSeq branch; ⚠ D1: same-name copies on one chromosome become ONE span)."""
    genes = {}
    if fmt == 'gff':
        for ln in open(path):
            f = ln.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] != chrom or f[2] not in ('gene', 'pseudogene'):
                continue
            m = re.search(r'(?:^|;)Name=([^;]+)', f[8])
            if m:
                genes[m.group(1)] = (int(f[3]) - 1, int(f[4]))
    else:
        ex = collections.defaultdict(list)
        for ln in open(path):
            f = ln.rstrip('\n').split('\t')
            if len(f) < 9 or f[0] != chrom or f[2] != 'exon':
                continue
            g = lib.gtf_attr(f[8], 'gene_name') or (lib.gtf_attr(f[8], 'gene_id') or '').replace('gene-', '', 1)
            if not g:
                continue  # exon lines without a gene attribute would otherwise form one chromosome-wide span
            ex[g].append((int(f[3]) - 1, int(f[4])))
        genes = {g: (min(s_ for s_, _ in v), max(e for _, e in v)) for g, v in ex.items()}
    return genes


def span_mapper(genes):
    """gene_of(s, e): the gene with the largest SPAN overlap, first maximum in (start, end, name) order (linear scan,
    verbatim from referee_band_score / identity_spectrum)."""
    glist = sorted((s_, e, g) for g, (s_, e) in genes.items())

    def gene_of(s_, e):
        best = None
        for gs, ge, g in glist:
            if ge <= s_:
                continue
            if gs >= e:
                break
            o = min(e, ge) - max(s_, gs)
            if o > 0 and (best is None or o > best[0]):
                best = (o, g)
        return best[1] if best else None
    return gene_of


def load_members(path, chrom, gene_of):
    """Member file -> (group -> genes, group -> genes of multi-exon members or None, members counted).
    A `cluster_id` header (mcl_families clusters.tsv: start/end in columns 7-8) is read WITHOUT a chromosome filter,
    as referee_band_score did; a `family_id` header (gw_family_catalog copies.tsv: chrom/start/end/n_exon in columns
    4-7) is restricted to `chrom`, as identity_spectrum --catalog did."""
    groups = collections.defaultdict(set)
    spliced = collections.defaultdict(set)
    fmt, n = None, 0
    for ln in open(path):
        f = ln.rstrip('\n').split('\t')
        if fmt is None:
            fmt = 'clusters' if f[0] == 'cluster_id' else ('copies' if f[0] == 'family_id' else None)
            if fmt is None:
                sys.exit(f'{path}: header must start with cluster_id (clusters.tsv) or family_id (copies.tsv)')
            continue
        if fmt == 'clusters':
            if f[0] == 'cluster_id':
                continue
            g = gene_of(int(f[6]), int(f[7]))
            if g:
                groups[f[0]].add(g)
        else:
            if f[0] == 'family_id' or f[3] != chrom:
                continue
            n += 1
            g = gene_of(int(f[4]), int(f[5]))
            if g:
                groups[f[0]].add(g)
                if int(f[6]) >= 2:
                    spliced[f[0]].add(g)   # multi-exon copies only (single-exon loci are mostly intronic/unspliced)
    return groups, (spliced if fmt == 'copies' else None), n


def pairs_of(groups):
    return {frozenset(p) for gs in groups.values() for p in itertools.combinations(sorted(gs), 2)}


def family_pairs(fam, keep=None):
    """Same-family gene pairs of a gene -> family table (optionally only genes in `keep`)."""
    by = collections.defaultdict(list)
    for g in sorted(fam):
        if keep is None or g in keep:
            by[fam[g]].append(g)
    return {frozenset(p) for gs in by.values() for p in itertools.combinations(gs, 2)}


def band_lines(counts, order):
    """(band, hit, n) sorted by n descending; ties by the fixed band order (the old scripts left ties to set order)."""
    rank = {b: i for i, b in enumerate(order)}
    return [(b, h, n) for b, (h, n) in sorted(counts.items(), key=lambda x: (-x[1][1], rank.get(x[0], len(rank))))]


def cmd_pairs(a):
    low = a.genes.lower()
    fmt = a.genes_format if a.genes_format != 'auto' else ('gtf' if low.endswith(('.gtf', '.gtf.gz')) else 'gff')
    spans = load_gene_spans(a.genes, a.chrom, fmt)
    if not spans:  # a GTF read as GFF (or the wrong --chrom) finds no genes and would silently score 0 pairs
        sys.exit(f'--genes {a.genes}: no {fmt.upper()} genes on {a.chrom}; check --chrom or pass --genes-format')
    gene_of = span_mapper(spans)
    kind, tpath = a.truth.split(':', 1)
    if kind not in ('families', 'compara'):
        sys.exit('--truth must be families:FILE or compara:FILE')
    if a.bands.startswith('paf:'):
        return pairs_referee_bands(a, gene_of, kind, tpath)
    if a.bands not in ('auto', 'compara') or (a.bands == 'compara' and kind != 'compara'):
        sys.exit("--bands must be auto, compara (with a compara: truth) or paf:FILE")
    return pairs_catalog(a, gene_of, kind, tpath)


def pairs_referee_bands(a, gene_of, kind, tpath):
    """referee_band_score.py: pair recall by ANNOTATED-mRNA identity band over the expressed universe, precision
    against the COMPLETE referee (expression is irrelevant to whether two genes are one family). Truth pairs =
    referee same-family pairs with both genes expressed; band = best minimap2 identity between the two genes'
    annotated mRNAs (nm/bl of any record); 'none' = no alignment record at all."""
    if kind != 'families' or not a.expressed:
        sys.exit('--bands paf: needs --truth families:FILE and --expressed FILE')
    fam = lib.read_referee(tpath)
    expr = set(l.split('\t')[0] for l in open(a.expressed) if not l.startswith('Gene'))
    ident = {}
    for ln in open(a.bands[4:]):
        f = ln.split('\t')
        if f[0] == f[5]:
            continue
        k = frozenset((f[0], f[5])); idn = int(f[9]) / max(1, int(f[10]))
        if k not in ident or idn > ident[k]:
            ident[k] = idn

    def band(k):
        if k not in ident:
            return 'none'
        for lo, hi, n in BANDS_MRNA:
            if lo <= ident[k] < hi:
                return n
    truth = {k: band(k) for k in family_pairs(fam, expr)}
    cg, _, _ = load_members(a.members, a.chrom, gene_of)
    pairs = pairs_of(cg)
    largest = max((len(v) for v in cg.values()), default=0)
    by = collections.defaultdict(lambda: [0, 0])
    for k, b in truth.items():
        by[b][1] += 1; by[b][0] += (k in pairs)
    jd = {k for k in pairs if all(g in fam for g in k)}; tp = sum(1 for k in jd if len({fam[g] for g in k}) == 1)
    big = max(cg.items(), key=lambda x: len(x[1]))[1]
    comp = collections.Counter(fam.get(g, 'not-in-referee') for g in sorted(big))
    order = ['>=90', '80-90', '70-80', '60-70', '<60', 'none']
    rec = ' · '.join(f"{b} {by[b][0]}/{by[b][1]}" for b in order if b in by)
    print(f"{a.label:16s} largest {largest:3d} genes {dict(comp.most_common(3))} | prec {tp}/{len(jd)}={tp/len(jd) if jd else 0:.3f} | recall by annotated-mRNA identity: {rec}")


def pairs_catalog(a, gene_of, kind, tpath):
    """identity_spectrum.py --catalog: a gw_family_catalog copies.tsv at the PAIR level against Compara (recall by
    Compara identity band) or a gene -> family referee (one band, 'all'); precision over judgeable catalog pairs, also
    for multi-exon copies only. ⚠ D11: the first recall block's denominator (truth pairs with BOTH genes in the
    catalog) is conditioned on the prediction; --universe adds recall over a fixed, prediction-independent universe."""
    fam_genes, fam_spliced, ncop = load_members(a.members, a.chrom, gene_of)
    if fam_spliced is None:
        sys.exit('the catalog report needs a gw_family_catalog copies.tsv (family_id header, n_exon column)')
    cat_pairs = pairs_of(fam_genes)
    spliced_pairs = pairs_of(fam_spliced)
    cat_genes = set().union(*fam_genes.values()) if fam_genes else set()
    largest = max((len(v) for v in fam_genes.values()), default=0)
    if kind == 'families':
        fam = lib.read_referee(tpath)
        truth_all = family_pairs(fam)
        judge = set(fam)
        band_of = lambda k: 'all'
        order = ['all']
    else:
        compara, judge = lib.load_compara(tpath, a.chrom)
        truth_all = set(compara)
        order = [n for _, _, n in BANDS_COMPARA]

        def band_of(k):
            p = compara[k][0]
            for lo, hi, n in BANDS_COMPARA:
                if lo <= p < hi:
                    return n
            return '<30'
    # recall over truth pairs whose BOTH genes are in the catalog (present as copies) — the catalog cannot join what
    # it did not emit; and, for context, over truth pairs with both genes mapped by any catalog copy or not
    truth_in = {k for k in truth_all if k <= cat_genes}
    print(f'[catalog] {a.members}: {ncop} copies on {a.chrom}, {len(fam_genes)} families, largest {largest} genes, {len(cat_pairs)} gene pairs; truth pairs {len(truth_all)}, with both genes in the catalog {len(truth_in)}')
    by_band = collections.defaultdict(lambda: [0, 0])
    for k in truth_in:
        b = band_of(k); by_band[b][1] += 1; by_band[b][0] += (k in cat_pairs)
    for b, hit, n in band_lines(by_band, order):
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
            if f[0] == 'geneA' or len(f) < 2:
                continue
            uni.add(frozenset((f[0], f[1])))
        uni &= truth_all
        ub = collections.defaultdict(lambda: [0, 0])
        for k in uni:
            b = band_of(k); ub[b][1] += 1; ub[b][0] += (k in cat_pairs)
        print(f'   recall over the fixed universe ({len(uni)} truth pairs, both genes expressed):')
        for b, hit, n in band_lines(ub, order):
            print(f'      {b:6s} {hit}/{n} = {hit/n:.3f}')


# ================================================================ spectrum (identity_spectrum tier mode)
def cmd_spectrum(a):
    """Identity spectrum (docs/PREREG_identity_spectrum_2026-09-24.md): which alignment tier of the family-edge builder
    recovers which Ensembl Compara paralogue pairs, band by band of identity, on one chromosome.

    Nodes: one spliced representative per expressed locus (gene_id group; most reads, tie longer). Tiers, all-vs-all on
    those nodes: T1 asm20 (k19) identity>=0.80 cov>=0.50; T2 asm20 -k11 -w5 identity>=0.60 cov>=0.50; T3 mmseqs
    translated (--search-type 2) protein identity>=0.30, e<=1e-5, qcov>=0.50. Loci map to HGNC symbols by exon overlap
    with the RefSeq annotation. Truth: Compara pairs (BioMart, columns gene, paralog, perc_id, perc_id_r1, subtype,
    paralog_chromosome) with both genes on the chromosome and both expressed.

    ⚠ B3: the default nucleotide identity/coverage (--estimator nm_bl --coverage query_over_min) are what rows
    1096-1099 were measured with, NOT the builder's (1-de, shorter-axis coverage; query_over_min exceeds 1 when the
    query is the longer sequence). `--estimator de --coverage shorter_axis` is the builder's rule — a separate arm."""
    import pysam
    fa = pysam.FastaFile(a.fasta)

    # ---- nodes: spliced representative per locus
    tx = collections.defaultdict(list); reads = {}; gene_of = {}; strand = {}
    for ln in open(a.gtf):
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != a.chrom:
            continue
        t = lib.gtf_attr(f[8], 'transcript_id')
        if f[2] == 'transcript':
            gene_of[t] = lib.gtf_attr(f[8], 'gene_id') or t; strand[t] = f[6]; r = lib.gtf_attr(f[8], 'reads'); reads[t] = int(r) if r else 0
        elif f[2] == 'exon':
            tx[t].append((int(f[3]) - 1, int(f[4])))
    loci = collections.defaultdict(list)
    for t, g in gene_of.items():
        loci[g].append(t)
    node = {}   # locus -> (span, exons, seq)
    with open(a.out + '.nodes.fa', 'w') as fh:
        for g, ts in loci.items():
            ts = [t for t in ts if tx.get(t)]
            if not ts:
                continue
            rep = max(ts, key=lambda t: (reads.get(t, 0), max(b for _, b in tx[t]) - min(s for s, _ in tx[t])))
            ex = sorted(tx[rep]); seq = ''.join(fa.fetch(a.chrom, s, e) for s, e in ex).upper()
            if strand[rep] == '-':
                seq = lib.rc(seq)
            if len(seq) < 200:
                continue
            node[g] = ((ex[0][0], ex[-1][1]), ex, seq); fh.write(f'>{g}\n{seq}\n')
    print(f'[spectrum] {len(node)} expressed loci with a representative >= 200 bp', flush=True)

    # ---- locus -> gene symbol (RefSeq exon overlap; the symbol is the gene_id with its "gene-" prefix stripped)
    ref_ex = collections.defaultdict(list)
    for ln in open(a.ref):
        f = ln.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != a.chrom or f[2] != 'exon':
            continue
        g = lib.gtf_attr(f[8], 'gene_id') or ''
        ref_ex[g.replace('gene-', '', 1)].append((int(f[3]) - 1, int(f[4])))
    ref_list = sorted((min(s for s, _ in v), max(e for _, e in v), g, v) for g, v in ref_ex.items())

    def symbol_of(exons):
        lo, hi = exons[0][0], exons[-1][1]; best = None
        for s, e, g, v in ref_list:
            if e <= lo:
                continue
            if s >= hi:
                break
            o = sum(max(0, min(b, y) - max(a_, x)) for a_, b in exons for x, y in v)
            if o > 0 and (best is None or o > best[0]):
                best = (o, g)
        return best[1] if best else None
    sym = {g: symbol_of(v[1]) for g, v in node.items()}
    sym = {g: s for g, s in sym.items() if s}
    by_sym = collections.defaultdict(list)
    for g, s in sym.items():
        by_sym[s].append(g)
    print(f'[spectrum] {len(sym)} loci map to {len(by_sym)} RefSeq symbols', flush=True)

    # ---- truth: Compara pairs on this chromosome
    compara, genes_with_data = lib.load_compara(a.compara, a.chrom)   # frozenset(symbols) -> (max perc_id, subtype)
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
            if q == t:
                continue
            idn = lib.paf_identity(f, a.estimator); cov = lib.paf_coverage(f, a.coverage)
            k = frozenset((q, t))
            if k not in best or (idn, cov) > best[k]:
                best[k] = (idn, cov)
            # union coverage on the SHORTER sequence's coordinates
            if int(f[1]) <= int(f[6]):
                spans[k].append((int(f[2]), int(f[3]), int(f[1])))
            else:
                spans[k].append((int(f[7]), int(f[8]), int(f[6])))
        for k, v in spans.items():
            v.sort(); cov = 0; cur = None
            for s_, e_, L in v:
                if cur is None or s_ > cur[1]:
                    if cur:
                        cov += cur[1] - cur[0]
                    cur = [s_, e_]
                else:
                    cur[1] = max(cur[1], e_)
            cov += cur[1] - cur[0]
            UNION[(out, k)] = cov / v[0][2]
        return best
    t1 = mm2('-x asm20', a.out + '.t1.paf'); t2 = mm2('-x asm20 -k11 -w5', a.out + '.t2.paf')
    subprocess.run(f"{a.mmseqs} easy-search {a.out}.nodes.fa {a.out}.nodes.fa {a.out}.t3.m8 {a.out}.tmp --search-type 2 --threads {a.threads} -e 1e-5 --format-output query,target,pident,qcov,tcov,evalue > /dev/null 2>&1", shell=True, check=True)
    t3 = {}
    for l in open(a.out + '.t3.m8'):
        f = l.split('\t')
        if f[0] == f[1]:
            continue
        idn = float(f[2]) / (100 if float(f[2]) > 1 else 1); cov = max(float(f[3]), float(f[4]))
        k = frozenset((f[0], f[1]))
        if k not in t3 or (idn, cov) > t3[k]:
            t3[k] = (idn, cov)

    def edge(best, k, floor):
        v = best.get(k); return v is not None and v[0] >= floor and v[1] >= 0.50

    def locus_pairs(sa, sb):
        return [frozenset((x, y)) for x in by_sym[sa] for y in by_sym[sb] if x != y]

    def tier_hit(sa, sb, tiers):
        return any(edge(b, k, fl) for b, fl in tiers for k in locus_pairs(sa, sb))
    TIERS = {'T1': [(t1, 0.80)], 'T1+T2': [(t1, 0.80), (t2, 0.60)], 'T1+T2+T3': [(t1, 0.80), (t2, 0.60), (t3, 0.30)]}
    BANDS = BANDS_COMPARA

    def band(p):
        for lo, hi, n in BANDS:
            if lo <= p < hi:
                return n
        return '<30'

    # ---- recall by Compara band
    rows = []; print('\n== RECALL of Compara paralogue pairs (both genes expressed), by Compara protein identity band')
    print(f"{'band':8s} {'n_pairs':>7s}  " + '  '.join(f'{t:>9s}' for t in TIERS))
    for lo, hi, n in BANDS:
        ks = [k for k, v in truth.items() if lo <= v[0] < hi]
        if not ks:
            continue
        rec = {t: sum(1 for k in ks if tier_hit(*sorted(k), tiers)) for t, tiers in TIERS.items()}
        print(f"{n:8s} {len(ks):7d}  " + '  '.join(f"{rec[t]/len(ks):9.3f}" for t in TIERS)); rows.append(('recall', n, len(ks), {t: rec[t] / len(ks) for t in TIERS}))
    # ---- precision by OUR identity band: aligned symbol pairs (both genes with Compara data) that are Compara paralogues
    print('\n== PRECISION of aligned pairs (both genes have Compara data), by the tier\'s own identity band')

    def sym_pairs(best, floor):
        out = {}
        for k, (idn, cov) in best.items():
            if idn < floor or cov < 0.50:
                continue
            x, y = tuple(k); sx, sy = sym.get(x), sym.get(y)
            if not sx or not sy or sx == sy:
                continue
            kk = frozenset((sx, sy))
            if kk not in out or idn > out[kk]:
                out[kk] = idn
        return out
    for name, best, floor in (('T1 (nt)', t1, 0.80), ('T2 (nt)', t2, 0.60), ('T3 (protein)', t3, 0.30)):
        sp = {k: v for k, v in sym_pairs(best, floor).items() if k <= genes_with_data}
        print(f'-- {name}: {len(sp)} judgeable aligned symbol pairs')
        for lo, hi, n in BANDS:
            ks = [k for k, v in sp.items() if lo <= v * 100 < hi]
            if ks:
                tp = sum(1 for k in ks if k in compara); print(f"   {n:8s} n={len(ks):5d}  precision {tp/len(ks):.3f}"); rows.append((f'precision {name}', n, len(ks), tp / len(ks)))
    missed = [(tuple(sorted(k)), v[0]) for k, v in truth.items() if not tier_hit(*sorted(k), TIERS['T1+T2+T3'])]

    # diagnosis of every truth pair: the best record each tier has for ANY locus pair of the two symbols (identity,
    # coverage), or none — separates "no seed" (no record) from "coverage clause" (record below 0.50) from "identity floor"
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
        if pid < 60 or tier_hit(*sorted(k), TIERS['T1+T2']):
            continue
        sa, sb = sorted(k); b = best_any(t2, sa, sb)
        why[band(pid)]['no nucleotide record (seeding)' if b is None else ('record, coverage < 0.50' if b[1] < 0.50 else 'record, identity < 0.60')] += 1
    print('== why truth pairs at >= 60% Compara identity are missed by T1+T2, per band:')
    for b_ in ('>=90', '80-90', '70-80', '60-70'):
        if why[b_]:
            print(f'   {b_:6s} {dict(why[b_])}')

    # post-hoc variants of the T2 coverage clause (NOT part of the bar): union-of-records coverage, and a 0.30 floor
    def edge_union(k, floor, cov_floor):
        v = t2.get(k); return v is not None and v[0] >= floor and UNION.get((a.out + '.t2.paf', k), 0) >= cov_floor
    print('== post-hoc: recall of Compara pairs under T2 coverage variants (identity >= 0.60)')
    print(f"{'band':8s} {'n':>5s} {'cov>=0.50':>10s} {'cov>=0.30':>10s} {'union>=0.50':>12s} {'union>=0.30':>12s}")
    for lo, hi, n in BANDS:
        ks = [k for k, v in truth.items() if lo <= v[0] < hi]
        if not ks:
            continue

        def rec(fn):
            return sum(1 for k in ks if any(fn(kk) for kk in locus_pairs(*sorted(k)))) / len(ks)
        print(f"{n:8s} {len(ks):5d} {rec(lambda kk: edge(t2, kk, 0.60)):10.3f} {rec(lambda kk: t2.get(kk) is not None and t2[kk][0] >= 0.60 and t2[kk][1] >= 0.30):10.3f} {rec(lambda kk: edge_union(kk, 0.60, 0.50)):12.3f} {rec(lambda kk: edge_union(kk, 0.60, 0.30)):12.3f}")

    # precision of the union variant, judgeable pairs, by identity band
    def sym_pairs_union(floor, cov_floor):
        out = {}
        for k, (idn, cov) in t2.items():
            if idn < floor or UNION.get((a.out + '.t2.paf', k), 0) < cov_floor:
                continue
            x, y = tuple(k); sx, sy = sym.get(x), sym.get(y)
            if not sx or not sy or sx == sy:
                continue
            kk = frozenset((sx, sy))
            if kk not in out or idn > out[kk]:
                out[kk] = idn
        return out
    for cf in (0.50, 0.30):
        sp = {k: v for k, v in sym_pairs_union(0.60, cf).items() if k <= genes_with_data}
        tp = sum(1 for k in sp if k in compara)
        print(f"   union>={cf}: {len(sp)} judgeable aligned pairs, precision {tp/len(sp) if sp else float('nan'):.3f}")
    print(f'\n== truth pairs recovered by NO tier: {len(missed)} of {len(truth)}; Compara identity median {sorted(p for _, p in missed)[len(missed)//2] if missed else "-"}; examples {missed[:8]}')
    with open(a.out + '.spectrum.tsv', 'w') as fh:
        fh.write('metric\tband\tn\tvalue\n')
        for m, b, n, v in rows:
            fh.write(f'{m}\t{b}\t{n}\t{v}\n')
    print(f'wrote {a.out}.spectrum.tsv')


# ================================================================ heldout (heldout_family_score)
HELDOUT_SUFFIX = re.compile(r'(?:P\d+|\d+|[A-Z])$')


def symbol_root(sym):
    """The pre-registered root: strip ONE trailing copy-suffix. Not applied repeatedly."""
    return HELDOUT_SUFFIX.sub('', sym)


def heldout_load_genes(gff, chrom):
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


def symbol_root_families(genes):
    """root -> [symbols]; >= 3 members, root >= 3 chars, LOC* excluded from truth."""
    by_root = collections.defaultdict(set)
    for sym in genes.values():
        if sym.startswith('LOC'):
            continue
        r = symbol_root(sym)
        if len(r) >= 3:
            by_root[r].add(sym)
    return {r: sorted(v) for r, v in by_root.items() if len(v) >= 3}


def heldout_predicted_clusters(clusters_tsv, genes, exact_only=False):
    """cluster_id -> [symbols] (members that resolve to a named gene; LOC members are KEPT).
    Default: probe (start+1, end) then (start, end), as the old script did; ⚠ B4 (register 966): mcl_families writes
    1-based coordinates verbatim, so --exact-only (probe (start, end) only) is the correct lookup for current files."""
    out = collections.defaultdict(list)
    with open(clusters_tsv) as fh:
        for line in fh:
            if line.startswith('cluster_id'):
                continue
            p = line.rstrip('\n').split('\t')
            if len(p) < 8:
                continue
            cid, s, e = p[0], int(p[6]), int(p[7])
            sym = genes.get((s, e)) if exact_only else (genes.get((s + 1, e)) or genes.get((s, e)))
            out[cid].append(sym if sym else f'{p[5]}:{s}-{e}')
    return dict(out)


def cmd_heldout(a):
    """Score `mcl_families` clusters against symbol-root truth families (docs/PREREG_heldout_families_2026-09-20.md),
    or against Soto S1C families with --soto. Implements the PRE-REGISTERED definitions verbatim:
      truth family  = >= 3 RefSeq `gene`/`pseudogene` records on the chromosome whose `Name=` shares a root,
                      root = re.sub(r'(?:P\\d+|\\d+|[A-Z])$', '', symbol) applied ONCE, roots shorter than 3 chars
                      dropped, `LOC*` symbols excluded from TRUTH but left in the INPUT (they can still cost precision);
      matching      = one-to-one bipartite, maximising total overlap (scipy linear_sum_assignment);
      sensitivity   = matched / truth members; precision = matched / members of the matched cluster; F = harmonic mean;
                      unmatched truth families score 0 and are KEPT in the pooled mean (lib.bipartite_families).
    --soto: Soto et al. 2025 families (S1C `Family ID`), >= 3 members on the chromosome, matched by `Gene Name`; a gene
    with more than one distinct Family ID is EXCLUDED (see bench/soto/soto_replication.py, load_truth)."""
    import json
    import numpy as np
    genes = heldout_load_genes(a.gff, a.chrom)
    if a.soto:
        truth = lib.families_on(lib.soto_gene_family(a.soto), set(genes.values()), 3)
    else:
        truth = symbol_root_families(genes)
    pred = heldout_predicted_clusters(a.clusters, genes, a.exact_only)
    per = lib.bipartite_families(truth, pred)
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


# ================================================================ referee (soto_vs_us_referee)
def cmd_referee(a):
    """Us vs Soto, scored against a NEUTRAL referee.

    Every number this session has quoted used Soto as the truth, so it measures agreement with Soto, not precision.
    To ask whether Soto is more precise ANYWHERE, both have to be scored against a third party.

    Referee: **protein families** (§6ko's rule — longest CDS per gene, translated, all-vs-all blastp e <= 1e-5, edge iff
    non-overlapping HSPs cover >= 0.30 of the longer protein, MCL I = 2.8, r2 exclusions; truth.protein_referee). It
    is independent of BOTH comparators: it never sees our genomic alignment gate, and it never sees Soto's SD/WSSD
    construction. It is amino-acid evidence about the product. The per-chromosome proteome and blastp table are
    cached in --workdir (<chrom>_ref.*).

    ⚠ Register T15: "never consume the comparator's own files and call it replication" — nothing here reads a
    Soto-derived file except the S1C family assignment being SCORED, which is the object under test.

    Reports, pooled and stratified by referee-family size: pairwise precision / recall / F for each comparator against
    the referee, where a "pair" is two genes the comparator places together (lib.pair_scores, fixed universe).

    ⚠ B1: the old soto_vs_us_referee.py crashed with NameError (`re`) on its first call since commit 8db314c7."""
    import pysam
    import truth as truthlib
    os.makedirs(a.workdir, exist_ok=True)
    soto = lib.soto_gene_family(a.soto)
    fa = pysam.FastaFile(a.genome)
    ours, sot, ref = {}, {}, {}
    for chrom in a.chroms.split(','):
        names = lib.gene_key_names(a.gff, chrom)
        fams = truthlib.protein_referee(a.gff, fa, chrom, f'{a.workdir}/{chrom}_ref', a.threads, reuse_faa=True)
        for f, mem in fams.items():
            for g in mem:
                ref[f'{chrom}:{g}'] = f'{chrom}:{f}'
        p = f'{a.clusters}/{chrom}_fam.clusters.tsv'
        for line in open(p):
            if line.startswith('cluster_id'):
                continue
            q = line.rstrip('\n').split('\t')
            g = names.get(f'{q[5]}:{q[6]}-{q[7]}')
            if g:
                ours[f'{chrom}:{g}'] = f'{chrom}:{q[0]}'
        for sp, g in names.items():
            if g in soto:
                sot[f'{chrom}:{g}'] = soto[g]

    print(f"referee: protein families — {len(set(ref.values()))} families over {len(ref)} genes\n")
    print(f"  {'comparator':12s} {'precision':>10} {'recall':>8} {'F':>8} {'pairs called':>13} {'TP':>6}")
    for name, lab in (('OURS (MCL)', ours), ('SOTO', sot)):
        p, r, f, np_, nt, tp = lib.pair_scores(lab, ref)
        print(f"  {name:12s} {p:>10.3f} {r:>8.3f} {f:>8.3f} {np_:>13} {tp:>6}")
    print(f"\n  (referee pairs available: {lib.pair_scores(ours, ref)[4]})")

    print("\nstratified by REFEREE family size — 'in any part':")
    print(f"  {'ref fam size':>13} {'genes':>6} {'OURS prec':>10} {'SOTO prec':>10} {'OURS rec':>9} {'SOTO rec':>9}")
    byr = collections.defaultdict(list)
    for g, f in ref.items():
        byr[f].append(g)
    for lo, hi, lbl in ((2, 2, '2'), (3, 4, '3-4'), (5, 9, '5-9'), (10, 10**9, '>=10')):
        keep = {g for f, v in byr.items() if lo <= len(v) <= hi for g in v}
        sub = {g: f for g, f in ref.items() if g in keep}
        if not sub:
            continue
        po, ro = lib.pair_scores({g: v for g, v in ours.items() if g in keep}, sub)[:2]
        ps, rs = lib.pair_scores({g: v for g, v in sot.items() if g in keep}, sub)[:2]
        print(f"  {lbl:>13} {len(keep):>6} {po:>10.3f} {ps:>10.3f} {ro:>9.3f} {rs:>9.3f}")


# ================================================================ edge-gap (protein_edge_gap)
def nucleotide_edges(paf, gene_at):
    """Shipped gate: identity >= 0.7, cov_longer >= 0.3, >= 300 bp matching (nm / longer, per RECORD)."""
    ed = set()
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        a, b = gene_at.get(f[0]), gene_at.get(f[5])
        if not a or not b or a == b:
            continue
        nmatch, alen = int(f[9]), int(f[10])
        if nmatch < 300 or alen == 0 or nmatch / alen < 0.7:
            continue
        longer = max(int(f[1]), int(f[6]))
        if longer and nmatch / longer >= 0.3:
            ed.add(tuple(sorted((a, b))))
    return ed


def cmd_edge_gap(a):
    """Does a protein-level edge close §6o8's no-edge gap? Per `docs/PREREG_protein_edges_2026-09-20.md`.

    For each chromosome: take Soto's published families as truth (external, unchanged from §6s8), then for every
    within-family PAIR ask whether it carries
      (a) a NUCLEOTIDE edge -- the shipped gene-body gate, read from the same PAF mcl_families consumed
          (identity >= 0.7, cov_longer >= 0.3, >= 300 bp matching), and
      (b) a PROTEIN edge -- §6ko's rule, copied and not re-tuned: one protein per gene (longest CDS, translated),
          all-vs-all blastp -evalue 1e-5, edge iff non-overlapping HSPs cover >= 0.30 of the longer protein.
    Reports the pair rates and the §6o8 statistic: families with NO edge on any member, nucleotide alone vs nucleotide
    union protein."""
    import pysam
    import truth as truthlib
    cds = lib.longest_cds(a.gff, a.chrom)
    truth = lib.families_on(lib.soto_gene_family(a.soto), set(cds), 2)
    if not truth:
        print(f'{a.chrom}: no Soto family with >=2 members carrying CDS'); return

    fa = pysam.FastaFile(a.genome)
    faa = a.out + '.proteins.faa'
    plen = truthlib.write_proteins(fa, a.chrom, cds, faa)

    # PAF sequence names are chrom:start-end; map them to symbols via the GFF gene spans
    gene_at = lib.gene_key_names(a.gff, a.chrom)

    nuc = nucleotide_edges(a.paf, gene_at)
    pro = truthlib.protein_edges(faa, a.out, plen, a.threads)

    pairs = n_nuc = n_pro = n_either = 0
    fam_nuc = fam_either = 0
    for members in truth.values():
        m = [g for g in members if g in plen]
        if len(m) < 2:
            continue
        has_n = has_e = False
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                k = tuple(sorted((m[i], m[j])))
                pairs += 1
                nn, pp = k in nuc, k in pro
                n_nuc += nn; n_pro += pp; n_either += (nn or pp)
                has_n |= nn; has_e |= (nn or pp)
        fam_nuc += not has_n
        fam_either += not has_e
    nf = sum(1 for v in truth.values() if len([g for g in v if g in plen]) >= 2)
    print(f'{a.chrom}: families {nf} | pairs {pairs} | '
          f'nuc {n_nuc} ({100*n_nuc/pairs if pairs else 0:.1f}%) | '
          f'prot {n_pro} ({100*n_pro/pairs if pairs else 0:.1f}%) | '
          f'either {n_either} ({100*n_either/pairs if pairs else 0:.1f}%) || '
          f'NO-EDGE families: nuc {fam_nuc}/{nf} ({100*fam_nuc/nf if nf else 0:.1f}%) -> '
          f'either {fam_either}/{nf} ({100*fam_either/nf if nf else 0:.1f}%)')


# ================================================================ rna-ceiling (rna_truth_from_protein)
def spliced_exons(gff, chrom):
    """gene symbol -> (strand, [(start1,end)]) exon union of the transcript with most exonic bases."""
    by_tx = collections.defaultdict(list)
    tx_gene, tx_strand = {}, {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom or f[2] != 'exon':
            continue
        p = re.search(r'Parent=([^;]+)', f[8]); g = re.search(r'gene=([^;]+)', f[8])
        if not p or not g:
            continue
        by_tx[p.group(1)].append((int(f[3]), int(f[4])))
        tx_gene[p.group(1)] = g.group(1); tx_strand[p.group(1)] = f[6]
    best = {}
    for tx, segs in by_tx.items():
        g = tx_gene[tx]; n = sum(e - s + 1 for s, e in segs)
        if g not in best or n > best[g][0]:
            best[g] = (n, tx_strand[tx], sorted(segs))
    return {g: (st, segs) for g, (n, st, segs) in best.items()}


def cmd_rna_ceiling(a):
    """Build a NON-CIRCULAR RNA-level truth from protein families, and measure the ceiling it implies.

    §6o9 closed the goal question with: *"a non-circular RNA-level truth is the missing ingredient"* — the DNA
    gene-span truth demands pairs that do not exist as RNA (only 6.7% align as spliced RNA, capping pairwise recall at
    0.052), and an alignability-derived truth is circular because it is defined by the same gate that builds the edges.
    Protein families (truth.protein_referee) escape both horns: they are defined on the SPLICED PRODUCT, and built by
    blastp over AMINO ACIDS, a different alphabet, aligner and gate from the nucleotide edges they are used to score.

    Reports, per chromosome: (1) the spliced-RNA ALIGNABLE FRACTION of within-truth-family pairs (§6o9 measured 6.7%,
    5.0% through the shipped gate, for the DNA truth); (2) the implied PAIRWISE RECALL CEILING and the no-edge family
    fraction under this truth. ⚠ the "shipped gate" here is nm/bl >= 0.80 and nmatch / LONGER >= 0.50 (audit D8)."""
    import pysam
    import truth as truthlib
    fa = pysam.FastaFile(a.genome)
    # ---- the TRUTH: protein families (MCL I=2.8 over protein edges), §6ko's rule, not re-tuned
    truth = truthlib.protein_referee(a.gff, fa, a.chrom, a.out, a.threads)
    if not truth:
        print(f'{a.chrom}: no protein family with >= 2 members'); return

    # ---- the spliced RNA of every member, then all-vs-all NUCLEOTIDE alignment (the other alphabet)
    ex = spliced_exons(a.gff, a.chrom)
    members = sorted({g for v in truth.values() for g in v if g in ex})
    rna = a.out + '.rna.fa'
    with open(rna, 'w') as fh:
        for g in members:
            st, segs = ex[g]
            s = lib.spliced1(fa, a.chrom, segs, st)
            if len(s) >= 200:
                fh.write(f'>{g}\n{s}\n')
    paf = a.out + '.rna.paf'
    if not os.path.exists(paf):
        with open(paf, 'w') as fh:
            subprocess.run(['minimap2', '-x', 'asm20', '-c', '-X', '-N', '50', '-p', '0.1',
                            '-t', str(a.threads), rna, rna], stdout=fh,
                           stderr=subprocess.DEVNULL, check=True)

    # §6o9's own two gates: "align at all", and the shipped edge gate id >= 0.80 AND cov >= 0.50
    aligned, gated = set(), set()
    for line in open(paf):
        f = line.rstrip('\n').split('\t')
        if len(f) < 12:
            continue
        q, s = f[0], f[5]
        if q == s:
            continue
        k = tuple(sorted((q, s)))
        aligned.add(k)
        nmatch, alen = int(f[9]), int(f[10])
        if alen and nmatch / alen >= 0.80 and nmatch / max(int(f[1]), int(f[6])) >= 0.50:
            gated.add(k)

    have = {g for g in members}
    pairs = n_al = n_ga = 0
    fam_no = nf = 0
    reach = 0
    for mem in truth.values():
        m = [g for g in mem if g in have]
        if len(m) < 2:
            continue
        nf += 1
        got = False
        # pairwise recall ceiling = pairs joined by ANY path in the gated graph, within the family
        adj = collections.defaultdict(set)
        for i in range(len(m)):
            for j in range(i + 1, len(m)):
                k = tuple(sorted((m[i], m[j])))
                pairs += 1
                if k in aligned:
                    n_al += 1
                if k in gated:
                    n_ga += 1; got = True
                    adj[m[i]].add(m[j]); adj[m[j]].add(m[i])
        fam_no += not got
        seen, comp = set(), []
        for g in m:
            if g in seen:
                continue
            st = [g]; c = []
            while st:
                x = st.pop()
                if x in seen:
                    continue
                seen.add(x); c.append(x); st.extend(adj[x] - seen)
            comp.append(len(c))
        reach += sum(c * (c - 1) // 2 for c in comp)
    print(f'{a.chrom}: protein-family truth {nf} families / {len(have)} members / {pairs} pairs | '
          f'align at all {n_al} ({100*n_al/pairs if pairs else 0:.1f}%) | '
          f'pass shipped gate {n_ga} ({100*n_ga/pairs if pairs else 0:.1f}%) | '
          f'PAIRWISE RECALL CEILING {reach}/{pairs} = {reach/pairs if pairs else 0:.3f} | '
          f'no-edge families {fam_no}/{nf} ({100*fam_no/nf if nf else 0:.1f}%)')


# ================================================================ members (member_completeness)
def cmd_members(a):
    """Member completeness for the family-scoped pool test (docs/PREREG_family_scoped_pool_2026-09-24.md).
    Universe = referee genes (one per line, 'Gene Name' header) fixed before the arms; runs gffcompare (-r REF_GTF) on
    ARM_GTF (outputs next to it, prefix ARM.gc) and reports: complete members (>= 1 transcript of class '='), partial
    ('=', 'c', 'k'), transcripts per member locus, and gffcompare transcript-level sensitivity/precision."""
    gtf, ref, uni, chrom, label = a.gtf, a.ref, a.universe, a.chrom, a.label
    universe = [l.rstrip('\n').split('\t')[0] for l in open(uni) if not l.startswith('Gene Name') and l.strip()]
    # gene name -> ref gene_id in the ref GTF (gffcompare's tmap uses ref_gene_id = gene_id attr, e.g. gene-LOC...)
    name_of_gid = {}
    for l in open(ref):
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom:
            continue
        gid = re.search(r'gene_id "([^"]+)"', f[8]); gn = re.search(r'gene_name "([^"]+)"', f[8])
        if gid and gn:
            name_of_gid[gid.group(1)] = gn.group(1)
    out = os.path.splitext(gtf)[0] + '.gc'
    subprocess.run(['gffcompare', '-r', ref, '-o', out, gtf], capture_output=True, text=True)
    tmap = [p for p in os.listdir(os.path.dirname(gtf) or '.') if p.startswith(os.path.basename(out)) and p.endswith('.tmap')]
    tmap = os.path.join(os.path.dirname(gtf) or '.', tmap[0])
    best = collections.defaultdict(set); ntx_gene = collections.Counter()
    for l in open(tmap):
        f = l.rstrip('\n').split('\t')
        if f[0] == 'ref_gene_id':
            continue
        g = name_of_gid.get(f[0], f[0]); best[g].add(f[2]); ntx_gene[g] += 1
    U = set(universe)
    complete = sum(1 for g in U if '=' in best.get(g, ()))
    partial = sum(1 for g in U if best.get(g, set()) & {'=', 'c', 'k'})
    with_tx = [g for g in U if g in ntx_gene]
    tx_per = sum(ntx_gene[g] for g in with_tx) / max(1, len(with_tx))
    stats = open(out + '.stats').read() if os.path.exists(out + '.stats') else open(out).read()  # a dotted prefix makes gffcompare write the summary to the bare prefix
    prec = re.search(r'Transcript level:\s+([\d.]+)\s+\|\s+([\d.]+)', stats)
    print(f"{label:14s} universe {len(U)} | complete (=) {complete} | partial (=,c,k) {partial} | transcripts per member gene {tx_per:.2f} | gffcompare transcript sens/prec {prec.group(1)}/{prec.group(2)}")


# ================================================================ adjudicated (adjudicated_truth score)
def cmd_adjudicated(a):
    """Score catalogs against the AK adjudicated truth (`truth.py adjudicated`): best-overlap assignment of truth loci
    to method families (ties -> the lexicographically largest family_id, as before); TP/FP on TRUE/FALSE pairs (a pair
    absent from pairs.tsv counts as FALSE), UNSCORED ignored; item-level bipartite (lib.bipartite_items) on the
    connected components of TRUE pairs (with --expr: loci with u >= 3 only, components recomputed)."""
    contigs = set(a.contigs.split(","))
    loci = {r["locus"]: (r["chrom"], int(r["start"]), int(r["end"])) for r in csv.DictReader(open(f"{a.truth}/loci.tsv"), delimiter="\t")
            if r["chrom"] in contigs}
    expressed = None
    if a.expr:
        expressed = {r["locus"] for r in csv.DictReader(open(a.expr), delimiter="\t") if int(r["u"]) >= 3}
    status = {}
    for r in csv.DictReader(open(f"{a.truth}/pairs.tsv"), delimiter="\t"):
        if r["u"] in loci and r["v"] in loci and (expressed is None or (r["u"] in expressed and r["v"] in expressed)):
            status[(r["u"], r["v"])] = r["status"]
    true_pairs = [p for p, s in status.items() if s == "TRUE"]
    uf = lib.UF()
    for u, v in true_pairs:
        uf.union(u, v)
    tl = sorted(uf.p, key=lambda x: int(x[1:]))
    tlabel = [uf.find(x) for x in tl]
    print(f"truth: {len(tl)} loci in {len(set(tlabel))} clusters, {len(true_pairs)} TRUE pairs, "
          f"{sum(1 for s in status.values() if s == 'FALSE')} explicit FALSE, {sum(1 for s in status.values() if s == 'UNSCORED')} UNSCORED")
    # pairs are scored among truth loci only (the AG / `rna_truth.py` convention): span overlap cannot place a nested gene
    # inside a family member's intron, so assigning every joint locus manufactures false pairs (§6kn)
    order = tl
    print(f"{'catalog':16s} {'pair_sens':>9s} {'pair_prec':>9s} {'bip_R':>6s} {'bip_P':>6s} {'bip_F':>6s}  (TP / FP / ignored)")
    for spec in a.catalogs:
        name, path = spec.split("=", 1)
        by = collections.defaultdict(list)
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["chrom"] in contigs:
                by[r["chrom"]].append((int(r["start"]), int(r["end"]), r["family_id"]))
        for c in by:
            by[c].sort()
        pred = {}
        for x in order:
            c, s, e = loci[x]
            h = [(min(e, b) - max(s, a0), f) for a0, b, f in by[c] if a0 < e and s < b]
            if h:
                pred[x] = max(h)[1]
        fam = collections.defaultdict(list)
        for x, f in pred.items():
            fam[f].append(x)
        tp = fp = ign = 0
        for ms in fam.values():
            for u, v in itertools.combinations(sorted(ms, key=lambda x: int(x[1:])), 2):
                s = status.get((u, v), "FALSE")
                if s == "TRUE":
                    tp += 1
                elif s == "FALSE":
                    fp += 1
                else:
                    ign += 1
        sens = tp / max(1, len(true_pairs))
        prec = tp / max(1, tp + fp)
        plabel = [pred.get(x, f"none:{x}") for x in tl]
        br, bp = lib.bipartite_items(plabel, tlabel)
        f1 = 2 * br * bp / (br + bp) if br + bp else float("nan")
        print(f"{name:16s} {sens:9.3f} {prec:9.3f} {br:6.3f} {bp:6.3f} {f1:6.3f}  ({tp} / {fp} / {ign})")


# ================================================================ protein (protein_families score)
def load_fams(spec):
    nodes, fam = spec.split(":", 1)
    rows = list(csv.DictReader(open(fam), delimiter="\t"))
    return {r["idx"]: r["family_id"] for r in rows}, {
        r["idx"]: (r["chrom"], [tuple(map(int, b.split("-"))) for b in r["cds"].split(",")]) for r in rows}, nodes


def cmd_protein(a):
    """Cross-annotation scoring of protein families (Addendum AN): truth genes in families >= 2; each assigned the
    family of the test gene with the greatest CDS-base overlap; pairwise sens/prec (lib.pairwise) and item-level
    bipartite R/P/F (lib.bipartite_items). --truth NODES:FAMILIES, --test NAME=NODES:FAMILIES (repeatable)."""
    import truth as truthlib
    tfam, tcds, tnodes = load_fams(a.truth)
    cnt = collections.Counter(tfam.values())
    truth = sorted((k for k in tfam if cnt[tfam[k]] >= 2), key=int)
    contigs = {tcds[k][0] for k in truth}
    print(f"truth: {len(truth)} genes in {len({tfam[k] for k in truth})} families, "
          f"{sum(v * (v - 1) // 2 for v in cnt.values() if v >= 2)} pairs")
    print(f"{'catalog':14s} {'pair_sens':>9s} {'pair_prec':>9s} {'bip_R':>6s} {'bip_P':>6s} {'bip_F':>6s} unassigned")
    for spec in a.test:
        name, rest = spec.split("=", 1)
        nodes, fampath = rest.split(":", 1)
        tg = truthlib.load_genes(nodes, contigs)
        rule = max(a.rule, 1 if a.no_pseudogenes else 0)
        tg = {k: g for k, g in tg.items() if not truthlib.excluded(g["biotype"], rule)}
        pfam = {r["idx"]: r["family_id"] for r in csv.DictReader(open(fampath), delimiter="\t")}
        by = collections.defaultdict(list)
        for k, g in tg.items():
            segs = sorted((x, y) for x, y, _ in g["cds"])
            by[g["chrom"]].append((segs[0][0], segs[-1][1], k, segs))
        idx = {}
        for c, v in by.items():
            v.sort()
            idx[c] = (v, [x[0] for x in v], max(x[1] - x[0] for x in v))
        pred, un = [], 0
        for i, k in enumerate(truth):
            c, segs = tcds[k]
            s0, e0 = min(x for x, _ in segs), max(y for _, y in segs)
            v, starts, ml = idx.get(c, ([], [], 0))
            lo, hi = bisect.bisect_left(starts, s0 - ml), bisect.bisect_left(starts, e0)
            best = (0, None)
            for a0, a1, kk, tsegs in v[lo:hi]:
                if a1 <= s0:
                    continue
                ovl = sum(max(0, min(y, d) - max(x, b)) for x, y in segs for b, d in tsegs)
                if ovl > best[0]:
                    best = (ovl, kk)
            f = pfam.get(best[1]) if best[1] else None
            if f is None:
                un += 1
            pred.append(f or f"none:{i}")
        true = [tfam[k] for k in truth]
        ps, pp = lib.pairwise(pred, true)
        br, bp = lib.bipartite_items(pred, true)
        f1 = 2 * br * bp / (br + bp) if br + bp else float("nan")
        print(f"{name:14s} {ps:9.3f} {pp:9.3f} {br:6.3f} {bp:6.3f} {f1:6.3f} {un}")


# ================================================================ eichler (eichler_compare)
def cmd_eichler(a):
    """Eichler-style AS-margin assignment, computed alongside ours and compared.

    The method the advisor cites: a multi-mapping read is assigned to its best alignment iff no other alignment scores
    within T units of it (T = 10 by convention); otherwise the read is discarded as ambiguous. It is a MARGIN rule on
    the aligner's own score. `copy_assign` already emits `as_best`, `as_second` and `as_margin` per read alongside our
    `status`, so both calls come from one file (the margin is always recomputed from `as_margin`).

        EICHLER(T):  as_margin >= T          -> assign to the best-AS copy
                     as_margin <  T          -> discard (ambiguous)
        OURS:        status in {assigned, tied, ambiguous}, assign-or-abstain, never 1/k

    ⚠ The two rules do not have the same SUBJECT: our AS-tied gate deliberately selects the reads where the aligner is
    indifferent (margin ~ 0), exactly the population Eichler's rule discards by construction. So "agreement" is not the
    interesting number; the interesting number is what each rule decides on the population the other keeps."""
    rows = list(csv.DictReader(open(a.assignments), delimiter='\t'))
    if not rows:
        raise SystemExit('no rows')

    def num(r, k):
        v = (r.get(k) or '').strip()
        try:
            return float(v)
        except ValueError:
            return None

    joint = collections.Counter()
    margins = collections.Counter()
    n = 0
    for r in rows:
        m = num(r, 'as_margin')
        n += 1
        # ⚠ A read with NO rival placement has margin NA, and nothing is within T of it, so Eichler
        # ASSIGNS it. An earlier version of this tool skipped those rows and undercounted his
        # assignments by 1,522 on the YAG substrate (2,536 instead of 4,058) -- the Rust
        # `--eichler-margin` implementation exposed it.
        eich = 'assign' if (m is None or m >= a.threshold) else 'discard'
        ours = (r.get('status') or '').strip()
        joint[(ours, eich)] += 1
        margins['no rival' if m is None else ('>=T' if m >= a.threshold else ('0' if m == 0 else '0<m<T'))] += 1

    print(f"reads: {n}   (Eichler threshold T = {a.threshold:g})\n")
    print("AS-margin distribution")
    for k in ('0', '0<m<T', '>=T', 'no rival'):
        print(f"  margin {k:6s} {margins[k]:>7}  {100*margins[k]/n if n else 0:>5.1f}%")

    ours_vals = sorted({k[0] for k in joint})
    print(f"\njoint decision table  (rows = OURS, cols = EICHLER T={a.threshold:g})")
    print(f"  {'':12s} {'assign':>9} {'discard':>9} {'total':>8}")
    for o in ours_vals:
        aa, dd = joint[(o, 'assign')], joint[(o, 'discard')]
        print(f"  {o:12s} {aa:>9} {dd:>9} {aa+dd:>8}")
    ta = sum(joint[(o, 'assign')] for o in ours_vals)
    td = sum(joint[(o, 'discard')] for o in ours_vals)
    print(f"  {'TOTAL':12s} {ta:>9} {td:>9} {ta+td:>8}")

    ours_assign = sum(v for k, v in joint.items() if k[0] == 'assigned')
    print(f"\n  Eichler assigns  {ta:>7} / {n} = {100*ta/n if n else 0:.1f}%")
    print(f"  we assign        {ours_assign:>7} / {n} = {100*ours_assign/n if n else 0:.1f}%")
    both = joint[('assigned', 'assign')]
    print(f"  both assign      {both:>7}")
    print(f"  we assign where Eichler discards: {joint[('assigned','discard')]}")
    print(f"  Eichler assigns where we abstain: "
          f"{sum(joint[(o,'assign')] for o in ours_vals if o != 'assigned')}")

    if a.out:
        with open(a.out, 'w') as fh:
            fh.write('ours\teichler\tn\n')
            for (o, e), v in sorted(joint.items()):
                fh.write(f'{o}\t{e}\t{v}\n')
        print(f"\n  wrote {a.out}")


# ================================================================ reads (copy_assign_read_truth score)
def cmd_reads(a):
    """Per-read scoring of `copy_assign --families` output against simulated read truth (`sim.py copies`, read names
    `family|copy|i`): per MAPQ-0 primary read, correct / wrong / conflict / abstain / lost under three readings of the
    per-family table (OWN = the read's true family's row, PRIMARY = rows with primary_local=1, ANY = any assigned row),
    by divergence bin of the source copy (1 - max_family_identity from the catalog).
    docs/PREREG_o2_read_truth_2026-09-23.md."""
    P, O = a.prefix, a.o2prefix
    CAT = a.catalog or os.environ.get('CATALOG_TSV')
    if not CAT:
        sys.exit('score.py reads: give --catalog CAT.copies.tsv (or set CATALOG_TSV)')
    cat = {}; div = {}
    for r in csv.DictReader(open(CAT), delimiter='\t'):
        cat[(r['family_id'], r['copy_idx'])] = (r['chrom'], int(r['start']), int(r['end']))
        if r.get('max_family_identity') not in (None, '', 'NA'):
            div[(r['family_id'], r['copy_idx'])] = 1 - float(r['max_family_identity'])

    def dbin(d):
        if d is None:
            return 'NA'
        return '<0.5%' if d < 0.005 else '0.5-1%' if d < 0.01 else '1-2%' if d < 0.02 else '2-5%' if d < 0.05 else '>=5%'

    def same_locus(x, y):
        if x is None or y is None or x[0] != y[0]:
            return False
        o = min(x[2], y[2]) - max(x[1], y[1]); return o >= 0.5 * min(x[2] - x[1], y[2] - y[1])
    prim = {}
    for ln in lib.sam_lines(['-F', '2308', P + '.bam']):   # streamed (was one captured whole-BAM string)
        f = ln.split('\t'); prim[f[0]] = int(f[4])
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(O + '.assignments.tsv'), delimiter='\t'):
        by[r['read_name']].append(r)

    def tr(n):
        return tuple(n.split('|')[:2])

    def judge(rows_assigned, t):
        loci = {(r['family_id'], r['catalog_copy_idx']) for r in rows_assigned}
        if not loci:
            return 'abstain'
        ok = [k == t or same_locus(cat.get(k), cat.get(t)) for k in loci]
        if len(loci) > 1 and not all(ok):
            return 'conflict' if any(ok) else 'wrong'
        return 'correct' if all(ok) else 'wrong'
    S = {v: collections.defaultdict(collections.Counter) for v in ('OWN', 'PRIMARY', 'ANY')}
    n_mapq0 = 0; own_status = collections.Counter(); nprim = collections.Counter()
    for name, mq in prim.items():
        if mq != 0:
            continue
        n_mapq0 += 1
        t = tr(name); b = dbin(div.get(t)); rows = by.get(name, [])
        asg = lambda rs: [r for r in rs if r['status'] == 'assigned' and r['origin_rejected'] == '0']
        own = [r for r in rows if r['family_id'] == t[0]]
        own_status[tuple(sorted(r['status'] for r in own)) or ('no_row',)] += 1
        o_own = 'lost' if not rows else ('no_own_row' if not own else judge(asg(own), t))
        pr = [r for r in rows if r['primary_local'] == '1']; nprim[len(pr)] += 1
        o_pr = 'lost' if not rows else ('no_primary_row' if not pr else judge(asg(pr), t))
        o_any = 'lost' if not rows else judge(asg(rows), t)
        for view, o in (('OWN', o_own), ('PRIMARY', o_pr), ('ANY', o_any)):
            for key in ('ALL', b):
                S[view][key][o] += 1
    print(f'MAPQ-0 reads {n_mapq0}; own-family row status combos: {own_status.most_common(6)}; primary_local rows per read: {dict(nprim)}')
    order = ['ALL', '<0.5%', '0.5-1%', '1-2%', '2-5%', '>=5%', 'NA']
    for view in ('OWN', 'PRIMARY', 'ANY'):
        print(f'== {view}')
        for k in order:
            c = S[view].get(k)
            if not c:
                continue
            n = sum(c.values()); na = c['correct'] + c['wrong'] + c['conflict']
            acc = c['correct'] / na if na else float('nan')
            print(f"  {k:8s} n={n:5d} correct {c['correct']:4d} wrong {c['wrong']:4d} conflict {c['conflict']:4d} abstain {c['abstain']:4d} other {n - na - c['abstain']:4d} | acc_assigned {acc:.3f} coverage {na / n:.3f}")


# ================================================================ bakeoff-calls / bakeoff-compare (copy_assign_tool_bakeoff)
def cmd_bakeoff_calls(a):
    """Derive a per-molecule COPY CALL from any isoform tool's GTF, by one identical rule (PREREG tool_bakeoff
    2026-09-08, hard_locus_bakeoff 5ca5c7e4).

    Why this exists: scoring tools on "copy attribution accuracy" is vacuous — StringTie, flair and isoseq collapse emit
    no copy attribute, so they score 0 by construction. Every tool DOES emit transcripts with genomic coordinates, and a
    molecule's intron chain either matches a transcript or does not. So the copy call is DERIVABLE for all of them:

      molecule --(exact intron chain)--> transcript(s) --(containment in a copy interval)--> copy

    States, identical for every arm: derived_one (every matching transcript sits in ONE copy) / derived_multi
    (matching transcripts span >= 2 copies: the tool conflated copies) / derived_none (no transcript carries this
    molecule's chain). A tool's own declared state (--own) is reported separately and never enters the derived columns.
    ⚠ B5: a transcript's copy is the copy with the max RAW overlap, strict > over start-sorted copies (ties -> the
    leftmost copy), not the "max reciprocal overlap, ties -> lowest copy_idx" the old docstring claimed."""
    gtf, bam, copies_p = a.gtf, a.bam, a.copies
    label = a.label
    # Junction tolerance. isoseq collapse runs --max-fuzzy-junction 5 by default, so ITS transcript
    # junctions may sit up to 5 bp off the read's, while our GTF is an exact intron-chain collapse and
    # matches by construction. Scoring at fuzz 0 therefore biases `derived_none` in OUR favour. Applied
    # symmetrically to every arm; report both 0 and 5.
    FUZZ = int(a.fuzz)
    own_p = a.own
    out_p = a.out
    # PREREG hard_locus_bakeoff: score only the listed molecules (one read name per line)
    restrict = set(l.strip() for l in open(a.restrict)) if a.restrict else None

    # ---------------------------------------------------------------- copies
    cop = list(csv.DictReader(open(copies_p), delimiter='\t'))
    copies = [(r['chrom'], int(r['start']), int(r['end']), int(r['copy_idx'])) for r in cop]
    by_chrom = collections.defaultdict(list)
    for c, s, e, i in copies:
        by_chrom[c].append((s, e, i))
    for c in by_chrom:
        by_chrom[c].sort()

    def copy_of(chrom, s, e):
        """The copy a transcript belongs to: max raw overlap, ties -> the leftmost copy (B5).
        Returns None when the transcript touches no copy at all."""
        best, best_ov = None, 0
        for cs, ce, ci in by_chrom.get(chrom, ()):
            o = min(e, ce) - max(s, cs)
            if o > best_ov:
                best, best_ov = ci, o
        return best

    # ---------------------------------------------------------------- transcripts
    ex = collections.defaultdict(list)
    for line in open(gtf):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        m = re.search(r'transcript_id[ =]"?([^";]*)"?', f[8])
        if not m:
            continue
        t = m.group(1)
        if f[2] == 'exon':
            ex[t].append((int(f[3]) - 1, int(f[4]), f[0]))

    # chain -> set of copies asserting it; also keep unspliced transcripts by span
    chain_to_copies = collections.defaultdict(set)
    unspliced = []          # (chrom, start, end, copy)
    n_tx = n_tx_in_copy = 0
    tx_multi = []
    copies_hit = set()
    fuzzy_chains = {}
    for t, v in ex.items():
        v.sort()
        chrom = v[0][2]
        s, e = v[0][0], v[-1][1]
        n_tx += 1
        ci = copy_of(chrom, s, e)
        if ci is None:
            continue
        n_tx_in_copy += 1
        # CONFLATION, measured on the transcript rather than on the chain. An intron chain carries
        # genomic coordinates, so two copies can NEVER assert the same chain and a chain-level
        # "spans >= 2 copies" test can never fire (found 2026-09-08, before the flair arm ran;
        # PREREG amendment 1). A transcript that OVERLAPS >= 2 copy intervals is the real thing.
        n_ov = sum(1 for cs, ce, _ in by_chrom.get(chrom, ()) if min(e, ce) - max(s, cs) > 0)
        if n_ov >= 2:
            tx_multi.append((t, chrom, s, e, n_ov))
        copies_hit.add(ci)
        chain = tuple((x[1], y[0]) for x, y in zip(v, v[1:]))
        if chain:
            # index every junction under a rounded key so a fuzzy lookup is O(1) per read
            key = (chrom,) + tuple((x // (2 * FUZZ + 1), y // (2 * FUZZ + 1)) for x, y in chain) if FUZZ else (chrom,) + chain
            chain_to_copies[key].add(ci)
            if FUZZ:
                fuzzy_chains.setdefault(key, []).append((chain, ci))
        else:
            unspliced.append((chrom, s, e, ci))

    # ---------------------------------------------------------------- molecules
    regions = []
    for c in by_chrom:
        cur = None
        for s, e, _ in by_chrom[c]:
            if cur and s <= cur[1]:
                cur[1] = max(cur[1], e)
            else:
                cur = [s, e]
                regions.append((c, cur))

    state = collections.Counter()
    mol_call = {}
    seen = set()
    for c, (lo, hi) in regions:
        # -F 2308: primary, mapped, non-supplementary -- the invariant for per-read statistics
        for ln in lib.sam_lines(['-F', '2308', bam, f'{c}:{lo+1}-{hi}']):
            f = ln.split('\t', 6)
            name, pos, cig = f[0], int(f[3]) - 1, f[5]
            if name in seen or (restrict is not None and name not in restrict):
                continue
            seen.add(name)
            ich = lib.cigar_introns(pos, cig)
            if ich:
                if FUZZ:
                    cps = set()
                    # a read matches a transcript when every junction is within FUZZ bp. The bucket key
                    # can straddle a boundary, so probe the neighbouring bucket on each coordinate too.
                    base = 2 * FUZZ + 1
                    seenk = set()
                    for d1 in (0, -1, 1):
                        for d2 in (0, -1, 1):
                            k = (c,) + tuple(((x // base) + d1, (y // base) + d2) for x, y in ich)
                            if k in seenk:
                                continue
                            seenk.add(k)
                            for cand, ci2 in fuzzy_chains.get(k, ()):
                                if len(cand) == len(ich) and all(
                                        abs(x[0] - y[0]) <= FUZZ and abs(x[1] - y[1]) <= FUZZ
                                        for x, y in zip(cand, ich)):
                                    cps.add(ci2)
                else:
                    cps = chain_to_copies.get((c,) + ich, set())
            else:
                # empty-chain trap (register 757): an unspliced read matches every unspliced
                # transcript's empty chain, so it must be resolved by SPAN CONTAINMENT instead.
                end = pos
                for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cig):
                    if op in 'M=XDN':
                        end += int(n)
                cps = {ci for tc, ts, te, ci in unspliced if tc == c and ts <= pos and end <= te}
            if not cps:
                st = 'derived_none'
            elif len(cps) == 1:
                st = 'derived_one'
            else:
                st = 'derived_multi'
            state[st] += 1
            mol_call[name] = (st, sorted(cps))

    tot = sum(state.values())
    print(f'== {label}')
    print(f'   transcripts {n_tx}  in a copy {n_tx_in_copy}')
    print(f'   transcripts spanning >=2 copies {len(tx_multi)}  ({len(tx_multi)/n_tx_in_copy:.3f} of in-copy)'
          if n_tx_in_copy else '   transcripts spanning >=2 copies 0')
    print(f'   copies with >=1 transcript      {len(copies_hit)} / {len(copies)}')
    print(f'   molecules   {tot}')
    for k in ('derived_one', 'derived_multi', 'derived_none'):
        v = state[k]
        print(f'   {k:14s} {v:7d}  {v/tot:6.3f}' if tot else f'   {k:14s} {v:7d}')

    # ---------------------------------------------------------------- the tool's OWN state, if any
    if own_p:
        own = collections.Counter()
        own_call = {}
        with open(own_p) as fh:
            rd = csv.DictReader(fh, delimiter='\t')
            cols = rd.fieldnames or []
            sc = next((c for c in ('status', 'decision', 'call') if c in cols), None)
            nc = next((c for c in ('read', 'read_name', 'molecule', 'name') if c in cols), None)
            # catalog_copy_idx FIRST: `assigned_copy` is the SWEEP index and does not address
            # copies.tsv (register 756 -- the two were once printed side by side and disagreed).
            cc = next((c for c in ('catalog_copy_idx', 'copy_idx', 'copy') if c in cols), None)
            for r in rd:
                if not sc or not nc:
                    break
                own[r[sc]] += 1
                if cc:
                    own_call[r[nc]] = (r[sc], r[cc])
        if own:
            print(f'   -- own declared states ({sc}):')
            for k, v in own.most_common():
                print(f'      {k:20s} {v:7d}')
        if own_call:
            # Only rows the tool itself CALLED count. Including its abstentions would compare the
            # derived rule against a non-decision and inflate the denominator.
            DECIDED = {'assigned'}
            agree = dis = 0
            for n, (st, cps) in mol_call.items():
                o = own_call.get(n)
                if o and o[0] in DECIDED and st == 'derived_one' and o[1].isdigit():
                    if int(o[1]) == cps[0]:
                        agree += 1
                    else:
                        dis += 1
            if agree + dis:
                print(f'   -- derived vs own copy, where the TOOL decided: {agree}/{agree+dis} = {agree/(agree+dis):.3f}')

    if out_p:
        with open(out_p + '.calls.tsv', 'w') as fh:
            fh.write('molecule\tstate\tcopies\n')
            for n, (st, cps) in sorted(mol_call.items()):
                fh.write(f'{n}\t{st}\t{",".join(map(str, cps))}\n')
        print(f'   wrote {out_p}.calls.tsv')


def load_calls(p):
    out = {}
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["molecule"]] = (r["state"], [int(x) for x in r["copies"].split(",") if x])
    return out


def cmd_bakeoff_compare(a):
    """PREREG hard_locus_bakeoff (5ca5c7e4): compare per-tool derived calls (`score.py bakeoff-calls --out`)."""
    tools = [(t.split("=")[0], load_calls(t.split("=")[1])) for t in a.tools]
    ours_l, ours = tools[0]
    assign = {r["read_name"]: r for r in csv.DictReader(open(a.assign), delimiter="\t")}
    allm = set(ours)
    for _, c in tools:
        allm |= set(c)
    if a.min_mult and a.bam:
        chain = {}
        for ln in lib.sam_lines(["-F", "2308", a.bam]):   # streamed (B6: was one captured whole-BAM string)
            f = ln.split("\t", 6)
            chain.setdefault(f[0], (f[2],) + lib.cigar_introns(int(f[3]) - 1, f[5]))
        mult = collections.Counter(chain.values())
        keep = {m for m in allm if mult.get(chain.get(m), 0) >= a.min_mult}
        print(f"--min-mult {a.min_mult}: keeping {len(keep)} of {len(allm)} molecules whose exact chain has >= {a.min_mult} molecules")
        allm = keep
    hard = {m for m in allm if m in assign}
    easy = allm - hard
    carried = lambda c, m: c.get(m, ("derived_none", []))[0] != "derived_none"
    strata = {
        "hard (all gate rows)": hard,
        "  contested": {m for m in hard if assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "    assigned": {m for m in hard if assign[m]["status"] == "assigned" and assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "    tied": {m for m in hard if assign[m]["status"] == "tied" and assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "    ambiguous": {m for m in hard if assign[m]["status"] == "ambiguous" and assign[m]["origin_rejected"] == "0" and int(assign[m]["n_candidates"]) >= 2},
        "  tie outside catalog": {m for m in hard if assign[m].get("tie_outside_catalog") == "1"},
        "  origin-rejected": {m for m in hard if assign[m]["origin_rejected"] == "1"},
        "easy (not admitted by the gate)": easy,
    }
    print(f"molecules: {len(allm)} total, hard {len(hard)}, easy {len(easy)}")
    print(f"\n{'stratum':34s} {'n':>6} " + " ".join(f"{l:>10}" for l, _ in tools) + "   (fraction carried = derived_one|multi)")
    rates = {}
    for name, S in strata.items():
        row = []
        for l, c in tools:
            k = sum(1 for m in S if carried(c, m))
            rates[(name, l)] = k / len(S) if S else float("nan")
            row.append(f"{k/len(S):10.3f}" if S else f"{'-':>10}")
        print(f"{name:34s} {len(S):>6} " + " ".join(row))
    # P2 discordance on hard and contested
    for name in ("hard (all gate rows)", "  contested"):
        S = strata[name]
        print(f"\n{name.strip()}: discordance vs ours")
        for l, c in tools[1:]:
            on = sum(1 for m in S if carried(ours, m) and not carried(c, m))
            xn = sum(1 for m in S if carried(c, m) and not carried(ours, m))
            both = sum(1 for m in S if carried(ours, m) and carried(c, m))
            nei = len(S) - on - xn - both
            print(f"   {l:10s} ours-not-{l}: {on:5d}   {l}-not-ours: {xn:5d}   ratio {on/max(1,xn):5.2f}   both {both}  neither {nei}")
    # P3 hard vs easy derived_none
    print("\nP3 derived_none hard / easy:")
    for l, c in tools:
        h = 1 - rates[("hard (all gate rows)", l)]; e = 1 - rates[("easy (not admitted by the gate)", l)]
        print(f"   {l:10s} hard {h:.3f}  easy {e:.3f}  ratio {h/e if e else float('nan'):.2f}")
    # P4 per copy: copies where ours carries >=1 hard molecule and X carries none (by the derived copy)
    print("\nP4 per-copy coverage on the hard set (copies where the tool carries >=1 hard molecule, by derived copy):")
    cov = {}
    for l, c in tools:
        cs = set()
        for m in hard:
            st, cps = c.get(m, ("derived_none", []))
            if st != "derived_none":
                cs.update(cps)
        cov[l] = cs
    for l, _ in tools[1:]:
        print(f"   {l:10s} copies {len(cov[l]):2d}; ours-only {sorted(cov[ours_l]-cov[l])}  {l}-only {sorted(cov[l]-cov[ours_l])}")
    print(f"   {ours_l:10s} copies {len(cov[ours_l])}: {sorted(cov[ours_l])}")
    # P6 copy attribution of the O2-assigned
    asg = strata["    assigned"]
    if a.gtf:
        tx_copy = {}   # read but not reported (kept as in the old script)
        for line in open(a.gtf):
            if "\ttranscript\t" not in line:
                continue
            t = re.search(r'transcript_id "([^"]+)"', line); ci = re.search(r'copy_index "([^"]+)"', line)
            if t and ci:
                tx_copy[t.group(1)] = ci.group(1)
    print(f"\nP6 the {len(asg)} O2-assigned molecules (report):")
    for l, c in tools:
        one = [m for m in asg if c.get(m, ("derived_none", []))[0] == "derived_one"]
        agree = sum(1 for m in one if str(c[m][1][0]) == assign[m]["catalog_copy_idx"])
        print(f"   {l:10s} carried {sum(1 for m in asg if carried(c, m)):3d}/{len(asg)}; derived_one {len(one):3d}, of which derived copy == O2 copy: {agree} ({100*agree/max(1,len(one)):.0f}%)")


# ================================================================ locus-reads (locus_reads.py CLI)
def cmd_locus_reads(a):
    """Count primaries (-F 2308) with an ALIGNED BLOCK in [START, END) (lib.reads_with_block_in, the one correct
    count) next to the misleading span-overlap count (lib.reads_overlapping_span). ⛔ ledger §6cm."""
    bam, chrom, start, end = a.bam, a.chrom, a.start, a.end
    tot, over = lib.reads_overlapping_span(bam, chrom, start, end)
    good = lib.reads_with_block_in(bam, chrom, start, end)
    print(f"{chrom}:{start}-{end}")
    print(f"  aligned block inside (USE THIS) : {good}")
    print(f"  overlapping the span            : {tot}")
    print(f"  spliced OVER, no aligned base   : {over} = {over / max(1, tot):.1%}")


# ================================================================ CLI
def _sub(sub, name, func, help_):
    doc = func.__doc__ or help_
    p = sub.add_parser(name, help=help_, description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.set_defaults(func=func)
    return p


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog='old -> new mapping: see the module docstring (python3 -c "import score; help(score)")')
    sub = ap.add_subparsers(dest='cmd', required=True)

    p = _sub(sub, 'pairs', cmd_pairs, 'pair-level family scoring (was referee_band_score.py + identity_spectrum.py --catalog)')
    p.description = (pairs_referee_bands.__doc__ + '\n\n' + pairs_catalog.__doc__ + '\n\nThe report is the referee-band '
                     'one when --bands paf:FILE is given, the catalog one otherwise.')
    p.add_argument('--members', required=True, help='mcl_families clusters.tsv or gw_family_catalog copies.tsv')
    p.add_argument('--genes', required=True, help='gene spans: a GFF (gene/pseudogene Name=) or a GTF (exon lines)')
    p.add_argument('--genes-format', choices=('auto', 'gff', 'gtf'), default='auto', help='auto: .gtf/.gtf.gz suffix (any case) -> gtf, else gff; zero genes found is an error')
    p.add_argument('--chrom', required=True)
    p.add_argument('--truth', required=True, help='families:FILE (Gene Name, Family ID) or compara:FILE (BioMart paralogues)')
    p.add_argument('--bands', default='auto', help="auto (Compara bands, or 'all' for a families truth) | compara | paf:FILE "
                   "(annotated-mRNA identity bands; the referee-band report)")
    p.add_argument('--expressed', help='referee-band report: expressed genes (first column, header starting "Gene")')
    p.add_argument('--universe', help='catalog report: a .truth_pairs.tsv (geneA, geneB); recall is ALSO reported over '
                   'its pairs, a denominator the catalog cannot move')
    p.add_argument('--label', default='', help='row label of the referee-band report')

    p = _sub(sub, 'spectrum', cmd_spectrum, 'edge tiers vs Compara, band by band (was identity_spectrum.py tier mode)')
    for k in ('--gtf', '--ref', '--fasta', '--chrom', '--compara', '--out'):
        p.add_argument(k, required=True)
    p.add_argument('--mmseqs', default='mmseqs'); p.add_argument('--threads', type=int, default=4)
    p.add_argument('--estimator', choices=('nm_bl', 'de'), default='nm_bl', help='nucleotide identity (default: the old nm/bl)')
    p.add_argument('--coverage', choices=('query_over_min', 'shorter_axis'), default='query_over_min',
                   help='nucleotide coverage (default: the old query span / min length, the M1 defect)')

    p = _sub(sub, 'heldout', cmd_heldout, 'mcl_families clusters vs symbol-root or Soto truth (was heldout_family_score.py)')
    p.add_argument('--gff', required=True); p.add_argument('--clusters', required=True); p.add_argument('--chrom', required=True)
    p.add_argument('--json')
    p.add_argument('--soto', help='score against Soto S1C published families instead of symbol roots')
    p.add_argument('--exact-only', action='store_true', help='look up members by (start, end) only (B4 fix; default also '
                   'probes (start+1, end) first, as the old script did)')

    p = _sub(sub, 'referee', cmd_referee, 'ours vs Soto against the protein referee (was soto_vs_us_referee.py)')
    for x in ('--gff', '--genome', '--soto', '--clusters', '--chroms'):
        p.add_argument(x, required=True)
    p.add_argument('--workdir', default='/mnt/linuxdisk/tmp/referee')
    p.add_argument('--threads', default='4', help='blastp threads (the old script hard-coded 4)')

    p = _sub(sub, 'edge-gap', cmd_edge_gap, 'nucleotide vs protein edges on Soto pairs (was protein_edge_gap.py)')
    for x in ('--gff', '--genome', '--paf', '--chrom', '--soto', '--out'):
        p.add_argument(x, required=True)
    p.add_argument('--threads', default='4')

    p = _sub(sub, 'rna-ceiling', cmd_rna_ceiling, 'protein-family truth and its RNA alignability ceiling (was rna_truth_from_protein.py)')
    for x in ('--gff', '--genome', '--chrom', '--out'):
        p.add_argument(x, required=True)
    p.add_argument('--threads', default='4')

    p = _sub(sub, 'members', cmd_members, 'gffcompare "=" per universe gene (was member_completeness.py)')
    p.add_argument('--gtf', required=True, help='the arm GTF (gffcompare writes ARM.gc.* next to it)')
    p.add_argument('--ref', required=True); p.add_argument('--universe', required=True)
    p.add_argument('--chrom', required=True); p.add_argument('--label', required=True)

    p = _sub(sub, 'adjudicated', cmd_adjudicated, 'catalogs vs the AK adjudicated truth (was adjudicated_truth.py score)')
    p.add_argument('--truth', required=True); p.add_argument('--contigs', required=True); p.add_argument('--expr')
    p.add_argument('catalogs', nargs='+', help='name=copies.tsv (chrom, start, end, family_id columns)')

    p = _sub(sub, 'protein', cmd_protein, 'cross-annotation protein-family scoring (was protein_families.py score)')
    p.add_argument('--truth', required=True); p.add_argument('--test', action='append', required=True)
    p.add_argument('--no-pseudogenes', action='store_true'); p.add_argument('--rule', type=int, default=0)

    p = _sub(sub, 'eichler', cmd_eichler, 'Eichler AS-margin rule vs ours (was eichler_compare.py)')
    p.add_argument('--assignments', required=True)
    p.add_argument('--threshold', type=float, default=10.0)
    p.add_argument('--out')

    p = _sub(sub, 'reads', cmd_reads, 'O2 per-read truth scoring (was copy_assign_read_truth.py score)')
    p.add_argument('--catalog', help='the catalog copies.tsv (default: $CATALOG_TSV)')
    p.add_argument('prefix', help='the sim.py copies OUT prefix (reads PREFIX.bam)')
    p.add_argument('o2prefix', help='the copy_assign --out prefix (reads O2PREFIX.assignments.tsv)')

    p = _sub(sub, 'bakeoff-calls', cmd_bakeoff_calls, 'per-molecule copy calls from a tool GTF (was copy_assign_tool_bakeoff.py calls)')
    p.add_argument('gtf'); p.add_argument('bam'); p.add_argument('copies')
    p.add_argument('--label', default='tool'); p.add_argument('--own'); p.add_argument('--out')
    p.add_argument('--restrict', help='score only the listed molecules (one read name per line)')
    p.add_argument('--fuzz', default='0', help='junction tolerance in bp')

    p = _sub(sub, 'bakeoff-compare', cmd_bakeoff_compare, 'compare per-tool calls (was copy_assign_tool_bakeoff.py compare)')
    p.add_argument("--assign", required=True)
    p.add_argument("--gtf", help="our GTF, for copy_index of the transcript carrying each molecule (P6)")
    p.add_argument("--bam", help="with --min-mult: primaries (-F 2308) give each molecule's intron chain")
    p.add_argument("--min-mult", type=int, default=0, help="keep only molecules whose exact chain is carried by >= N molecules (support-policy control)")
    p.add_argument("tools", nargs="+", help="label=calls.tsv; the first is ours")

    p = _sub(sub, 'locus-reads', cmd_locus_reads, 'the correct read count at a locus (was locus_reads.py)')
    p.add_argument('bam'); p.add_argument('chrom'); p.add_argument('start', type=int); p.add_argument('end', type=int)

    a = ap.parse_args(argv)
    a.func(a)


if __name__ == '__main__':
    main()

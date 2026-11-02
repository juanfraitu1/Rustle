#!/usr/bin/env python3
"""units_from_annotation.py ANNOT.gff GTF > list.tsv -- the GUIDED detector of `copy_assign --bridge-regroup f1units --bridge-units-list`.

Which transcripts of an assembled GTF are fusions of two annotated genes, and where to cut them into units
(docs/archive/2026-09/PREREG_container_units_v2_dev_2026-09-30.md Part C section 4.1b; it is the annotation-overlap oracle of
docs/archive/2026-09/PREREG_container_units_mechanism_2026-09-30.md section 3, without its NPIP-specific records).

  * GENES = the `gene` / `pseudogene` / `ncRNA_gene` records of ANNOT.gff (GFF3 with ID / Parent / Name / description, NCBI RefSeq
    style); a gene's exons are the union of every `exon` record
    whose Parent chain (<= 8 levels) reaches it, merged where they overlap or abut; a record without exons is dropped, and so is
    a readthrough gene (the fusions, not their parts): `--readthrough labels` (default) = an attribute of the GENE record contains
    `readthrough` (RefSeq description, GENCODE tag=readthrough_gene); `structural` = a transcript joins exons of two separate,
    non-small-RNA genes (opt-in: GENCODE recall 0.757, RefSeq 0.828, register 1492); `both`; `--readthrough-list FILE` always.
  * For every SPLICED transcript T of GTF (>= 2 exons): G(T) = the genes on T's strand whose exon union shares >= 1 base with an
    exon of T; E_g = the indices of T's exons they overlap; the index intervals [min E_g, max E_g] are merged where they overlap
    (a gene nested in another's interval, or two models over the same exons, are one block and never cut).
  * T has >= 2 blocks: between two consecutive blocks the cut is the WIDEST intron between the last exon of the upstream block and
    the first exon of the downstream one (ties: the first in genomic order). T with one block but >= 2 genes is UNCUTTABLE (shared
    exons or nested genes) and goes to the side report, never to the list.

stdout (and `--out FILE`): `tid chrom strand junctions label gene`, junctions `S-E;S-E` (1-based closed introns, genomic order), the
file `--bridge-units-list` reads. `label` = `--label-prefix` + the gene names, `gene` = T's gene_id in the GTF. `--report FILE`:
one row per uncuttable transcript (`tid chrom strand genes reason`); the counts go to stderr.

ANNOT.gff and GTF may be gzipped (`.gz`). An annotation with no fusion in the sample gives a HEADER-ONLY list, exit 0 and a note on
stderr; `copy_assign --bridge-units-list` refuses a list without a transcript row (`empty list: run without --bridge-units-list`), so a
per-sample loop passes the flag only when the list has a row.

Options: `--contigs a,b` (only these contigs: GFF and GTF), `--replace-genes FILE` (curated gene records that REPLACE a same-named
record or are added: `name chrom strand exons`, exons `S-E,S-E` 1-based closed; how a study supplies records the annotation lacks
or has wrong, e.g. NPIP copies) and `--label-prefix` (default `ann:`; the mechanism test's list was written with `H:`).
Only the standard library. Run `python3 bench/units_from_annotation.py --help`; tests: `python3 bench/test_units_from_annotation.py`.
"""
import argparse
import bisect
import collections
import gzip
import re
import sys

GENE_TYPES = ('gene', 'pseudogene', 'ncRNA_gene')


def open_text(path):
    """a text file, gunzipped when its name ends in `.gz`."""
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


_ATTR = {k: re.compile(r'(?:^|;)' + k + r'=([^;]+)') for k in ('ID', 'Parent', 'Name', 'description')}


def merge(iv):
    """0-based half-open intervals merged where they overlap or abut."""
    out = []
    for a, b in sorted(iv):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [tuple(x) for x in out]


SMALL_RNA = {'mirna', 'snorna', 'snrna', 'scarna', 'trna', 'rrna', 'misc_rna', 'y_rna', 'vault_rna', 'vaultrna'}
_BIOTYPE = re.compile(r'(?:^|;)(?:gene_biotype|gene_type|biotype)=([^;]+)')


def read_genes(gff, contigs=None, mode='labels', rt_list=None, calls=None):
    """[{id, name, chrom, strand, ex (merged 0-based half-open), lo, hi}] of the non-readthrough gene records with exons, in file order.

    Readthrough genes (register 1492, docs/PREREG_readthrough_detect_2026-11-02.md) are dropped by `mode`: `labels` = any attribute
    value of the gene record contains `readthrough` (RefSeq description, GENCODE tag=readthrough_gene, ...);
    `structural` = a transcript of the gene has an intron from an exon of gene A to an exon of gene B (same strand, A and B share no
    exonic base and do not overlap, neither a small RNA: exon union >= 200 bp and biotype not in SMALL_RNA); `both`; `rt_list` = gene
    IDs / names always dropped. `calls` (a dict) receives {gene name|start: reason} of every dropped gene."""
    ftype, parent, exons, genes = {}, {}, collections.defaultdict(list), []
    labelled = set()
    with open_text(gff) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or (contigs and f[0] not in contigs):
                continue
            a = f[8]
            mid, mp = _ATTR['ID'].search(a), _ATTR['Parent'].search(a)
            if mid:
                ftype[mid.group(1)] = f[2]
                if mp:
                    parent[mid.group(1)] = mp.group(1).split(',')[0]
            hit = 'readthrough' in a.lower()
            if f[2] in GENE_TYPES and mid:
                mn, mb = _ATTR['Name'].search(a), _BIOTYPE.search(a)
                genes.append((mid.group(1), mn.group(1) if mn else mid.group(1), f[0], f[6], mb.group(1).lower() if mb else '', f[3]))
                if hit:
                    labelled.add(mid.group(1))
            else:
                if hit and mid:
                    labelled.add(mid.group(1))
                if f[2] == 'exon' and mp:
                    exons[mp.group(1).split(',')[0]].append((int(f[3]) - 1, int(f[4])))
                    if hit:
                        labelled.add(mp.group(1).split(',')[0])

    def gene_of(x):
        for _ in range(8):
            if x is None or ftype.get(x) in GENE_TYPES:
                break
            x = parent.get(x)
        return x if ftype.get(x) in GENE_TYPES else None
    by_gene, tx_of = collections.defaultdict(list), collections.defaultdict(list)
    for p, ex in exons.items():
        g = gene_of(p)
        if g:
            by_gene[g].extend(ex)
            tx_of[g].append(sorted(ex))
    # a gene is labelled by its OWN record only (register 1493: the all-transcripts variant failed its bar as worded)
    lab = {g for g in labelled if ftype.get(g) in GENE_TYPES}
    rt_list = set(rt_list or ())
    all_genes = []
    for gid, name, chrom, strand, bt, start1 in genes:
        ex = merge(by_gene.get(gid, []))
        if ex:
            all_genes.append(dict(gid=gid, id=f'{name}|{start1}', name=name, chrom=chrom, strand=strand, ex=ex, lo=ex[0][0], hi=ex[-1][1], bt=bt))
    struct = set()
    if mode in ('structural', 'both'):
        parts = GeneIndex([g for g in all_genes if g['bt'] not in SMALL_RNA and sum(b - a for a, b in g['ex']) >= 200])
        for g in all_genes:
            for t in tx_of.get(g['gid'], ()):
                if g['gid'] in struct:
                    break
                for d, c in zip(t, t[1:]):
                    left = [h for h in parts.hits(g['chrom'], d[0], d[1], g['strand']) if h is not g]
                    right = [h for h in parts.hits(g['chrom'], c[0], c[1], g['strand']) if h is not g]
                    if any(x is not y and (x['hi'] <= y['lo'] or y['hi'] <= x['lo']) for x in left for y in right):
                        struct.add(g['gid'])
                        break
    out = []
    for g in all_genes:
        why = []
        if g['gid'] in rt_list or g['name'] in rt_list:
            why.append('list')
        if mode in ('labels', 'both') and g['gid'] in lab:
            why.append('label')
        if g['gid'] in struct:
            why.append('structural')
        if why:
            if calls is not None:
                calls[g['id']] = '+'.join(why)
            continue
        out.append(dict(id=g['id'], name=g['name'], chrom=g['chrom'], strand=g['strand'], ex=g['ex'], lo=g['lo'], hi=g['hi']))
    return out


def read_replacements(path):
    out = []
    with open_text(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 4 or f[0] in ('name', '') or line.startswith('#'):
                continue
            ex = merge([(int(a) - 1, int(b)) for a, b in (x.split('-') for x in f[3].split(',') if x)])
            out.append(dict(id=f'copy:{f[0]}', name=f[0], chrom=f[1], strand=f[2], ex=ex, lo=ex[0][0], hi=ex[-1][1]))
    return out


class GeneIndex:
    """genes by contig sorted by start; `hits` = the genes on a strand sharing >= 1 base with an interval."""

    def __init__(self, genes):
        self.by = collections.defaultdict(list)
        for g in genes:
            self.by[g['chrom']].append(g)
        self.lo, self.maxlen = {}, {}
        for c, v in self.by.items():
            v.sort(key=lambda g: g['lo'])
            self.lo[c] = [g['lo'] for g in v]
            self.maxlen[c] = max(g['hi'] - g['lo'] for g in v)

    def hits(self, chrom, a, b, strand):
        v = self.by.get(chrom)
        if not v:
            return []
        out, i = [], bisect.bisect_left(self.lo[chrom], a - self.maxlen[chrom])
        while i < len(v) and self.lo[chrom][i] < b:
            g = v[i]
            i += 1
            if g['strand'] == strand and g['hi'] > a and any(x < b and a < y for x, y in g['ex']):
                out.append(g)
        return out


def read_transcripts(gtf, contigs=None):
    """[(tid, gene_id, chrom, strand, [(s, e)] sorted 1-based closed)] of the spliced transcripts, in file order."""
    txs, order = {}, []
    with open_text(gtf) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] not in ('transcript', 'exon') or (contigs and f[0] not in contigs):
                continue
            tid = re.search(r'transcript_id "([^"]*)"', f[8])
            if not tid:
                continue
            tid = tid.group(1)
            if f[2] == 'transcript':
                g = re.search(r'gene_id "([^"]*)"', f[8])
                txs[tid] = [tid, g.group(1) if g else '', f[0], f[6], []]
                order.append(tid)
            elif tid in txs:
                txs[tid][4].append((int(f[3]), int(f[4])))
    out = []
    for tid in order:
        t = txs[tid]
        t[4].sort()
        if len(t[4]) >= 2:
            out.append(tuple(t))
    return out


def cuts_of(exons, chrom, strand, gix, cut_shared=False):
    """(cut introns, genes, uncuttable reason): see the module docstring. `cut_shared` (register 1496): a transcript whose genes
    share exons is cut at every intron between a PRIVATE exon (one gene) and a SHARED exon (>= 2 genes), so a cluster of
    alternative-promoter genes (PCDHA, UGT1A) leaves each gene's own exons as a unit apart from the shared constant exons."""
    E = collections.defaultdict(set)
    name = {}
    for j, (a, b) in enumerate(exons):
        for g in gix.hits(chrom, a - 1, b, strand):
            E[g['id']].add(j)
            name[g['id']] = g['name']
    if len(E) < 2:
        return [], sorted(name[k] for k in E), None
    iv = sorted((min(s), max(s), k) for k, s in E.items())
    blocks = []
    for lo, hi, k in iv:
        if blocks and lo <= blocks[-1][1]:
            blocks[-1][1] = max(blocks[-1][1], hi)
            blocks[-1][2].append(name[k])
        else:
            blocks.append([lo, hi, [name[k]]])
    genes = [n for _, _, ns in blocks for n in ns]
    if len(blocks) < 2:
        nhit = [sum(1 for s in E.values() if j in s) for j in range(len(exons))]
        shared = any(n >= 2 for n in nhit)
        if shared and cut_shared:
            cuts = [(exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1)
                    if (nhit[i] == 1 and nhit[i + 1] >= 2) or (nhit[i] >= 2 and nhit[i + 1] == 1)]
            if cuts:
                return cuts, genes, None
        return [], genes, 'shared_exon' if shared else 'nested_genes'
    cuts = []
    for (_, hi1, _), (lo2, _, _) in zip(blocks, blocks[1:]):
        best = None
        for i in range(hi1, lo2):                       # the intron between exon i and exon i + 1 (0-based)
            ln = exons[i + 1][0] - exons[i][1] - 1
            if best is None or ln > best[0]:
                best = (ln, i)
        cuts.append((exons[best[1]][1] + 1, exons[best[1] + 1][0] - 1))
    return cuts, genes, None


def run(gff, gtf, contigs=None, replace=None, label_prefix='ann:', mode='labels', rt_list=None, calls=None, cut_shared=False):
    """(list rows, uncuttable rows, counts)."""
    genes = read_genes(gff, contigs, mode, rt_list, calls)
    if replace:
        names = {g['name'] for g in replace}
        genes = [g for g in genes if g['name'] not in names] + [g for g in replace if not contigs or g['chrom'] in contigs]
    gix = GeneIndex(genes)
    rows, unc, n = [], [], collections.Counter()
    for tid, gene, chrom, strand, exons in read_transcripts(gtf, contigs):
        n['spliced'] += 1
        cuts, names, reason = cuts_of(exons, chrom, strand, gix, cut_shared)
        if len(names) >= 2:
            n['ge2_genes'] += 1
        if reason:
            n['uncuttable_' + reason] += 1
            unc.append((tid, chrom, strand, '|'.join(names), reason))
        elif cuts:
            n['transcripts'] += 1
            n['cuts'] += len(cuts)
            rows.append((tid, chrom, strand, ';'.join(f'{s}-{e}' for s, e in cuts), label_prefix + '|'.join(names), gene))
    return rows, unc, n


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('gff')
    ap.add_argument('gtf')
    ap.add_argument('--contigs', default='')
    ap.add_argument('--replace-genes', default='')
    ap.add_argument('--label-prefix', default='ann:')
    ap.add_argument('--report', default='')
    ap.add_argument('--out', default='')
    ap.add_argument('--readthrough', choices=('labels', 'structural', 'both'), default='labels')
    ap.add_argument('--readthrough-list', default='')
    ap.add_argument('--cut-shared-exon', action='store_true', help='cut private-to-shared exon introns (register 1496)')
    ap.add_argument('--readthrough-report', default='', help='dropped readthrough genes: gene reason')
    a = ap.parse_args(argv)
    contigs = set(a.contigs.split(',')) if a.contigs else None
    rt_list = [l.split('\t')[0].strip() for l in open_text(a.readthrough_list) if l.strip() and not l.startswith('#')] if a.readthrough_list else None
    calls = {}
    rows, unc, n = run(a.gff, a.gtf, contigs, read_replacements(a.replace_genes) if a.replace_genes else None, a.label_prefix,
                       a.readthrough, rt_list, calls, a.cut_shared_exon)
    if a.readthrough_report:
        with open(a.readthrough_report, 'w') as fh:
            fh.write('gene\treason\n')
            for k, v in calls.items():
                fh.write(f'{k}\t{v}\n')
    print(f'units_from_annotation: {len(calls)} readthrough genes dropped ({a.readthrough})', file=sys.stderr)
    out = open(a.out, 'w') if a.out else sys.stdout
    out.write('tid\tchrom\tstrand\tjunctions\tlabel\tgene\n')
    for r in rows:
        out.write('\t'.join(r) + '\n')
    if a.out:
        out.close()
    if a.report:
        with open(a.report, 'w') as fh:
            fh.write('tid\tchrom\tstrand\tgenes\treason\n')
            for r in unc:
                fh.write('\t'.join(r) + '\n')
    print(f'units_from_annotation: {n["spliced"]} spliced transcripts, {n["ge2_genes"]} over >= 2 genes: {n["transcripts"]} cut at '
          f'{n["cuts"]} introns, {n["uncuttable_shared_exon"]} uncuttable (shared exon), {n["uncuttable_nested_genes"]} uncuttable (nested genes)',
          file=sys.stderr)
    if not rows:
        print('units_from_annotation: no transcript to cut: the list is HEADER-ONLY (copy_assign --bridge-units-list refuses an empty list: '
              'run without it)', file=sys.stderr)


if __name__ == '__main__':
    main(sys.argv[1:])

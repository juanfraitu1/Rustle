"""Per-exon read support of assembled transcripts from the BAM CIGARs (user question 2026-10-28: exons no read really covers).

  exonsupport.py GTF BAM OUT_TSV [--contigs chr20,chr21]

For every transcript: each exon's cover = primary reads (flag & 2308 == 0) with an aligned block (M/=/X, D-joined) overlapping it, and
skip = primary reads whose N gap spans the whole exon. An exon is UNSUPPORTED when cover < 2 (the assembler's own chain floor) or
skip > cover (the §6er units rule: a base is exonic only when more blocks cover it than reads splice over it). No new numbers.
OUT_TSV: gene_id, transcript_id, reads(attr), n_exon, n_unsupported, min_cover, max_skip_over_cover, unsupported_idx (0-based, 5'->3' genomic).
"""
import argparse
import collections
import re

import pysam


def read_gtf(path, contigs):
    tx = collections.OrderedDict()
    for line in open(path):
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] not in ('transcript', 'exon') or (contigs and f[0] not in contigs):
            continue
        a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
        t = tx.setdefault(a['transcript_id'], dict(gene=a['gene_id'], chrom=f[0], strand=f[6], exons=[], reads=0))
        if f[2] == 'transcript':
            t['reads'] = int(a.get('reads', 0))
        else:
            t['exons'].append((int(f[3]) - 1, int(f[4])))
    for t in tx.values():
        t['exons'] = sorted(set(t['exons']))
    return tx


def blocks_and_gaps(aln):
    """Aligned blocks (M/=/X runs, D joined) and N gaps, in reference coordinates."""
    pos, blocks, gaps, cur = aln.reference_start, [], [], None
    for op, n in aln.cigartuples:
        if op in (0, 7, 8, 2):  # M = X D
            cur = [pos, pos + n] if cur is None else [cur[0], pos + n]
            pos += n
        elif op == 3:  # N
            if cur:
                blocks.append(tuple(cur))
            cur = None
            gaps.append((pos, pos + n))
            pos += n
    if cur:
        blocks.append(tuple(cur))
    return blocks, gaps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('gtf')
    ap.add_argument('bam')
    ap.add_argument('out')
    ap.add_argument('--contigs')
    ap.add_argument('--include-secondary', action='store_true', help='count secondary alignments too (the assembler seed pool); default primaries only')
    a = ap.parse_args()
    contigs = set(a.contigs.split(',')) if a.contigs else None
    tx = read_gtf(a.gtf, contigs)
    by_gene = collections.defaultdict(list)
    for tid, t in tx.items():
        by_gene[t['gene']].append(tid)
    bam = pysam.AlignmentFile(a.bam)
    w = open(a.out, 'w')
    w.write('gene_id\ttranscript_id\treads\tn_exon\tn_unsupported\tmin_cover\tmax_skip_over_cover\tunsupported_idx\n')
    tot = collections.Counter()
    for gene, tids in by_gene.items():
        chrom = tx[tids[0]]['chrom']
        s = min(e[0] for tid in tids for e in tx[tid]['exons'])
        e = max(e[1] for tid in tids for e in tx[tid]['exons'])
        reads = []
        for aln in bam.fetch(chrom, s, e):
            if aln.flag & (2052 if a.include_secondary else 2308) or aln.cigartuples is None:
                continue
            reads.append(blocks_and_gaps(aln))
        for tid in tids:
            t = tx[tid]
            rows = []
            for i, (x, y) in enumerate(t['exons']):
                cover = sum(1 for bl, _ in reads if any(b0 < y and b1 > x for b0, b1 in bl))
                skip = sum(1 for _, gp in reads if any(g0 <= x and g1 >= y for g0, g1 in gp))
                rows.append((cover, skip))
            uns = [i for i, (c, k) in enumerate(rows) if c < 2 or k > c]
            tot['tx'] += 1
            tot['tx_unsupported'] += bool(uns)
            tot['exons'] += len(rows)
            tot['exons_unsupported'] += len(uns)
            w.write(f"{gene}\t{tid}\t{t['reads']}\t{len(rows)}\t{len(uns)}\t{min(c for c, _ in rows)}\t"
                    f"{max((k / max(c, 1) for c, k in rows), default=0):.2f}\t{','.join(map(str, uns))}\n")
    print(f"transcripts {tot['tx']}, with >= 1 unsupported exon {tot['tx_unsupported']} ({tot['tx_unsupported'] / max(tot['tx'], 1):.1%}); "
          f"exons {tot['exons']}, unsupported {tot['exons_unsupported']} ({tot['exons_unsupported'] / max(tot['exons'], 1):.1%})")


if __name__ == '__main__':
    main()

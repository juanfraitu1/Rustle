#!/usr/bin/env python3
"""Does a protein-level edge close §6o8's no-edge gap? Per `docs/PREREG_protein_edges_2026-09-20.md`.

For each chromosome: take Soto's published families as truth (external, unchanged from §6s8), then for
every within-family PAIR ask whether it carries
  (a) a NUCLEOTIDE edge -- the shipped gene-body gate, read from the same PAF mcl_families consumed
      (identity >= 0.7, cov_longer >= 0.3, >= 300 bp matching), and
  (b) a PROTEIN edge -- §6ko's rule, copied and not re-tuned: one protein per gene (longest CDS,
      translated), all-vs-all blastp -evalue 1e-5, edge iff non-overlapping HSPs cover >= 0.30 of the
      longer protein.

Reports the pair rates and the §6o8 statistic: families with NO edge on any member, nucleotide alone
vs nucleotide union protein.

Usage: protein_edge_gap.py --gff G --genome FA --paf P --chrom chrN --soto S1C --out PREFIX
"""
import argparse
import collections
import csv
import os
import re
import subprocess
import sys

BLAST = os.environ.get('BLAST_BIN', '/home/juanfra/miniforge3/envs/blast/bin')
CODON = {}
_B = 'TCAG'
_AA = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for _i, _a in enumerate(_AA):
    CODON[_B[_i // 16] + _B[_i // 4 % 4] + _B[_i % 4]] = _a
COMP = str.maketrans('ACGTacgtN', 'TGCAtgcaN')


def longest_cds(gff, chrom):
    """gene symbol -> (strand, [(start1, end)]) for the transcript with the most CDS bases."""
    by_tx = collections.defaultdict(list)
    tx_gene, tx_strand = {}, {}
    gene_of_id = {}
    for line in open(gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != chrom:
            continue
        if f[2] in ('gene', 'pseudogene'):
            i = re.search(r'ID=([^;]+)', f[8]); n = re.search(r'Name=([^;]+)', f[8])
            if i and n:
                gene_of_id[i.group(1)] = n.group(1)
        elif f[2] == 'CDS':
            p = re.search(r'Parent=([^;]+)', f[8]); g = re.search(r'gene=([^;]+)', f[8])
            if not p:
                continue
            tx = p.group(1)
            by_tx[tx].append((int(f[3]), int(f[4])))
            tx_strand[tx] = f[6]
            if g:
                tx_gene[tx] = g.group(1)
    best = {}
    for tx, segs in by_tx.items():
        g = tx_gene.get(tx)
        if not g:
            continue
        n = sum(e - s + 1 for s, e in segs)
        if g not in best or n > best[g][0]:
            best[g] = (n, tx_strand[tx], sorted(segs))
    return {g: (st, segs) for g, (n, st, segs) in best.items()}


def translate(fa, chrom, strand, segs):
    import pysam
    seq = ''.join(fa.fetch(chrom, s - 1, e) for s, e in segs).upper()
    if strand == '-':
        seq = seq.translate(COMP)[::-1]
    aa = ''.join(CODON.get(seq[i:i + 3], 'X') for i in range(0, len(seq) - len(seq) % 3, 3))
    return aa.split('*')[0] if aa.startswith('M') else aa.replace('*', 'X')


def nucleotide_edges(paf, gene_at):
    """Shipped gate: identity >= 0.7, cov_longer >= 0.3, >= 300 bp matching."""
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


def protein_edges(faa, out, plen, threads):
    bl = out + '.blastp.tsv'
    if not os.path.exists(bl):
        subprocess.run([BLAST + '/makeblastdb', '-dbtype', 'prot', '-in', faa, '-out', out + '_db'],
                       stdout=subprocess.DEVNULL, check=True)
        with open(bl + '.tmp', 'w') as fh:
            subprocess.run([BLAST + '/blastp', '-query', faa, '-db', out + '_db', '-evalue', '1e-5',
                            '-max_target_seqs', '100000', '-num_threads', str(threads), '-outfmt',
                            '6 qseqid sseqid nident length qstart qend sstart send bitscore'],
                           stdout=fh, check=True)
        os.replace(bl + '.tmp', bl)
    hs = collections.defaultdict(list)
    for line in open(bl):
        q, s, nid, ln, q0, q1, s0, s1, bits = line.rstrip('\n').split('\t')
        if q != s:
            hs[tuple(sorted((q, s)))].append((float(bits), int(q0), int(q1), q, s))
    ed = set()
    for (a, b), v in hs.items():
        longer = max(plen.get(a, 0), plen.get(b, 0))
        if not longer:
            continue
        # greedy non-overlapping HSPs by bitscore, projected onto the LONGER protein
        taken = []
        for bits, q0, q1, q, s in sorted(v, reverse=True):
            if plen.get(q, 0) != longer:
                continue
            lo, hi = min(q0, q1), max(q0, q1)
            if all(hi < t0 or lo > t1 for t0, t1 in taken):
                taken.append((lo, hi))
        cov = sum(hi - lo + 1 for lo, hi in taken) / longer
        if cov >= 0.30:
            ed.add((a, b))
    return ed


def soto_truth(s1c, on_chrom):
    ids = collections.defaultdict(set)
    for r in csv.DictReader(open(s1c), delimiter='\t'):
        fid = (r.get('Family ID') or '').strip(); nm = (r.get('Gene Name') or '').strip()
        if fid and fid != 'N/A' and nm:
            ids[nm].add(fid)
    fam = collections.defaultdict(set)
    for nm, f in ids.items():
        if len(f) == 1 and nm in on_chrom:
            fam[next(iter(f))].add(nm)
    return {k: sorted(v) for k, v in fam.items() if len(v) >= 2}


def main():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--genome', '--paf', '--chrom', '--soto', '--out'):
        ap.add_argument(x, required=True)
    ap.add_argument('--threads', default='4')
    a = ap.parse_args()
    import pysam

    cds = longest_cds(a.gff, a.chrom)
    truth = soto_truth(a.soto, set(cds))
    if not truth:
        print(f'{a.chrom}: no Soto family with >=2 members carrying CDS'); return

    fa = pysam.FastaFile(a.genome)
    faa = a.out + '.proteins.faa'
    plen = {}
    with open(faa, 'w') as fh:
        for g, (st, segs) in cds.items():
            p = translate(fa, a.chrom, st, segs)
            if len(p) >= 10:
                plen[g] = len(p); fh.write(f'>{g}\n{p}\n')

    # PAF sequence names are chrom:start-end; map them to symbols via the GFF gene spans
    gene_at = {}
    for line in open(a.gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != a.chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        n = re.search(r'Name=([^;]+)', f[8])
        if n:
            gene_at[f'{f[0]}:{f[3]}-{f[4]}'] = n.group(1)

    nuc = nucleotide_edges(a.paf, gene_at)
    pro = protein_edges(faa, a.out, plen, a.threads)

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


if __name__ == '__main__':
    main()

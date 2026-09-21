#!/usr/bin/env python3
"""Build a NON-CIRCULAR RNA-level truth from protein families, and measure the ceiling it implies.

§6o9 closed the goal question with: *"a non-circular RNA-level truth is the missing ingredient"* — the
DNA gene-span truth demands pairs that do not exist as RNA (only 6.7% align as spliced RNA, capping
pairwise recall at 0.052), and an alignability-derived truth is circular because it is defined by the
same gate that builds the edges.

Protein families escape both horns:
  * they are defined on the SPLICED PRODUCT (longest CDS, translated), not on a genomic span, so their
    pairs are products — the thing an RNA-level definition is about;
  * they are built by blastp over AMINO ACIDS, a different alphabet, aligner and gate from the
    nucleotide spliced-minimap2 edges they are used to score. Nothing in the truth is computed from the
    edge set under test.

What this reports, per chromosome and pooled:
  1. the spliced-RNA ALIGNABLE FRACTION of within-truth-family pairs — §6o9 measured 6.7% (5.0% through
     the shipped gate) for the DNA truth; that number is the ceiling, so this is the headline;
  2. the implied PAIRWISE RECALL CEILING and the no-edge family fraction under this truth.

Usage: rna_truth_from_protein.py --gff G --genome FA --chrom chrN --out PREFIX [--threads 4]
"""
import argparse
import collections
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mcl_port  # noqa: E402
from protein_edge_gap import longest_cds, protein_edges, translate  # noqa: E402
from protein_families import excluded  # noqa: E402  -- §6ko r2: pseudogenes + V(D)J segments


def spliced_exons(gff, chrom):
    """gene symbol -> (strand, [(start1,end)]) exon union of the transcript with most exonic bases."""
    import re
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


COMP = str.maketrans('ACGTacgtN', 'TGCAtgcaN')


def main():
    ap = argparse.ArgumentParser()
    for x in ('--gff', '--genome', '--chrom', '--out'):
        ap.add_argument(x, required=True)
    ap.add_argument('--threads', default='4')
    a = ap.parse_args()
    import pysam

    # ---- the TRUTH: protein families (MCL I=2.8 over protein edges), §6ko's rule, not re-tuned
    cds = longest_cds(a.gff, a.chrom)
    # §6ko rule r2, applied verbatim: exclude pseudogenes AND V(D)J recombining segments. Without it an
    # IGKV cluster merges into one 75-member "family" and dominates every pair-weighted statistic.
    bt = {}
    import re as _re
    for line in open(a.gff):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if len(f) < 9 or f[0] != a.chrom or f[2] not in ('gene', 'pseudogene'):
            continue
        n = _re.search(r'Name=([^;]+)', f[8]); b = _re.search(r'gene_biotype=([^;]+)', f[8])
        if n:
            bt[n.group(1)] = b.group(1) if b else ''
    cds = {g: v for g, v in cds.items() if not excluded(bt.get(g, ''), 2)}
    fa = pysam.FastaFile(a.genome)
    faa = a.out + '.proteins.faa'
    plen = {}
    with open(faa, 'w') as fh:
        for g, (st, segs) in cds.items():
            p = translate(fa, a.chrom, st, segs)
            if len(p) >= 10:
                plen[g] = len(p); fh.write(f'>{g}\n{p}\n')
    pedges = protein_edges(faa, a.out, plen, a.threads)
    fams = mcl_port.mcl({(x, y): 1.0 for x, y in pedges}, inflation=2.8)   # {(a,b): weight}
    truth = {}
    for i, mem in enumerate(fams if isinstance(fams, list) else fams.values()):
        m = sorted(set(mem))
        if len(m) >= 2:
            truth[f'PF{i}'] = m
    if not truth:
        print(f'{a.chrom}: no protein family with >= 2 members'); return

    # ---- the spliced RNA of every member, then all-vs-all NUCLEOTIDE alignment (the other alphabet)
    ex = spliced_exons(a.gff, a.chrom)
    members = sorted({g for v in truth.values() for g in v if g in ex})
    rna = a.out + '.rna.fa'
    with open(rna, 'w') as fh:
        for g in members:
            st, segs = ex[g]
            s = ''.join(fa.fetch(a.chrom, x - 1, y) for x, y in segs).upper()
            if st == '-':
                s = s.translate(COMP)[::-1]
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


if __name__ == '__main__':
    main()

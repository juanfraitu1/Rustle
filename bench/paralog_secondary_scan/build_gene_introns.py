#!/usr/bin/env python3
"""Parse an NCBI RefSeq GFF once into a compact per-gene table:
   gene_id, chrom, start, end, strand, n_tx, intron_chain

intron_chain = union over all the gene's annotated transcripts of
(donor_last_exon_base, acceptor_first_exon_base) junction pairs (1-based),
joined as "d:a,d:a,...". A read "matches the gene's chain" if it shares >=K
of these junctions. Genes with no introns (single-exon) get an empty chain.

Usage: build_gene_introns.py GFF OUT.tsv
"""
import sys, re

def attr(s, key):
    m = re.search(key + r'=([^;]+)', s)
    return m.group(1) if m else None

def main():
    gff, out = sys.argv[1], sys.argv[2]
    gene_meta = {}            # gene_id -> [chrom, start, end, strand]
    rna2gene = {}             # rna_id  -> gene_id
    exons = {}                # gene_id -> list of (start,end)
    txset = {}                # gene_id -> set(rna_id)
    n = 0
    with open(gff) as fh:
        for line in fh:
            if line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            ftype = f[2]
            if ftype == 'gene':
                gid = attr(f[8], 'ID')
                if gid:
                    gene_meta[gid] = [f[0], int(f[3]), int(f[4]), f[6]]
            elif ftype in ('mRNA', 'transcript', 'lnc_RNA', 'ncRNA',
                           'primary_transcript'):
                rid = attr(f[8], 'ID'); par = attr(f[8], 'Parent')
                if rid and par:
                    rna2gene[rid] = par
            elif ftype == 'exon':
                par = attr(f[8], 'Parent')
                if not par:
                    continue
                for p in par.split(','):
                    g = rna2gene.get(p)
                    if g is None:
                        continue
                    exons.setdefault(g, []).append((int(f[3]), int(f[4])))
                    txset.setdefault(g, set()).add(p)
            n += 1
    # derive intron junctions per gene (union across transcripts)
    # group exons per transcript to avoid cross-transcript spurious junctions:
    # re-collect exons keyed by transcript
    # (second light pass not needed; we approximate union via per-transcript)
    # Re-parse exons by transcript for correctness:
    tx_exons = {}             # rna_id -> list of (start,end)
    with open(gff) as fh:
        for line in fh:
            if line[0] == '#':
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon':
                continue
            par = attr(f[8], 'Parent')
            if not par:
                continue
            for p in par.split(','):
                if p in rna2gene:
                    tx_exons.setdefault(p, []).append((int(f[3]), int(f[4])))
    gene_introns = {}         # gene_id -> set((d,a))
    for rid, exs in tx_exons.items():
        g = rna2gene[rid]
        exs = sorted(exs)
        s = gene_introns.setdefault(g, set())
        for i in range(len(exs) - 1):
            d = exs[i][1]            # last base of upstream exon
            a = exs[i + 1][0]        # first base of downstream exon
            if a > d:
                s.add((d, a))
    with open(out, 'w') as o:
        o.write("gene_id\tchrom\tstart\tend\tstrand\tn_tx\tintron_chain\n")
        for gid, (chrom, st, en, strand) in gene_meta.items():
            ints = sorted(gene_introns.get(gid, []))
            chain = ",".join(f"{d}:{a}" for d, a in ints)
            ntx = len(txset.get(gid, set()))
            o.write(f"{gid}\t{chrom}\t{st}\t{en}\t{strand}\t{ntx}\t{chain}\n")
    print(f"parsed {n} feature lines; wrote {len(gene_meta)} genes to {out}",
          file=sys.stderr)
    multi = sum(1 for g in gene_meta if gene_introns.get(g))
    print(f"  {multi} multi-exon genes", file=sys.stderr)

if __name__ == '__main__':
    main()

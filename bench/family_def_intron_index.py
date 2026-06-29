#!/usr/bin/env python3
"""
family_def_intron_index.py — ONE streaming pass over the gorilla RefSeq GFF to build a
per-gene INTRON ARCHITECTURE index (the canonical-transcript intron count + ordered
intron boundary positions), for use as a junction-concordance edge filter on the cDNA
homology graph (angle A2).

Rationale (orthogonal-to-sequence signal): real paralog families arose by WHOLE-GENE
duplication, so the copies share intron-exon ARCHITECTURE (homologous junction positions
=> matching ordered intron-gap structure). Domain-sharers share only a domain (partial /
different architecture). Retrocopies (processed pseudogenes) are INTRONLESS (0 introns)
yet have cDNA id ~1.0. So intron-architecture concordance can cut domain-bridge / retrocopy
edges that pure sequence id cannot.

Canonical transcript per gene = the mRNA/transcript with the MOST exons (ties -> longest
genomic span). We index, per gene Name:
    n_exon, n_intron, strand, exon_lengths (5'->3'), intron_lengths (5'->3')

We DO NOT use absolute genomic coords across genes (paralogs sit at different loci); the
architecture comparison uses the ORDERED exon/intron LENGTH vectors, which are duplication-
preserved and locus-independent.

Output: /home/juanfra/winloci_scratch/gene_intron_index.json  (gene -> dict)

Run: /home/juanfra/miniforge3/bin/python bench/family_def_intron_index.py
"""
import collections
import json
import re

GFF = "/mnt/c/Users/jfris/Desktop/GGO_genomic.gff"
OUT = "/home/juanfra/winloci_scratch/gene_intron_index.json"

GENE_RE = re.compile(r"(?:^|;)gene=([^;]+)")
PARENT_RE = re.compile(r"(?:^|;)Parent=([^;]+)")
ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")

# transcript-like feature types whose children exons define an architecture
TX_TYPES = {"mRNA", "transcript", "lnc_RNA", "ncRNA", "primary_transcript",
            "rRNA", "snRNA", "snoRNA", "tRNA", "miRNA"}


def main():
    # First pass collect: transcript_id -> gene Name ; and exon blocks per transcript_id.
    tx2gene = {}
    tx_exons = collections.defaultdict(list)   # tx_id -> list of (start,end,strand)
    n_lines = 0
    with open(GFF) as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            ftype = p[2]
            attrs = p[8]
            if ftype in TX_TYPES:
                m = ID_RE.search(attrs)
                g = GENE_RE.search(attrs)
                if m and g:
                    tx2gene[m.group(1)] = g.group(1)
            elif ftype == "exon":
                par = PARENT_RE.search(attrs)
                if not par:
                    continue
                start = int(p[3]); end = int(p[4]); strand = p[6]
                tx_exons[par.group(1)].append((start, end, strand))
            n_lines += 1
            if n_lines % 5_000_000 == 0:
                print(f"  ...{n_lines:,} lines", flush=True)
    print(f"[gff] parsed {n_lines:,} feature lines; "
          f"{len(tx2gene):,} transcripts; {len(tx_exons):,} transcripts-with-exons",
          flush=True)

    # Build per-transcript architecture, then pick canonical (most exons) per gene.
    gene_best = {}   # gene -> (n_exon, span, record)
    for tx, exons in tx_exons.items():
        gene = tx2gene.get(tx)
        if gene is None:
            # some exons (pseudogene single-exon) carry gene= directly but no tx parent in TX_TYPES;
            # skip — handled below via direct-gene exons
            continue
        exons.sort()
        strand = exons[0][2]
        # 5'->3' order
        ordered = exons if strand != "-" else list(reversed(exons))
        exon_lengths = [e - s + 1 for (s, e, _) in ordered]
        # introns = gaps between consecutive exons in genomic order, oriented 5'->3'
        gsorted = sorted(exons)
        intron_lengths_g = [gsorted[i + 1][0] - gsorted[i][1] - 1 for i in range(len(gsorted) - 1)]
        intron_lengths = intron_lengths_g if strand != "-" else list(reversed(intron_lengths_g))
        n_exon = len(exons)
        span = gsorted[-1][1] - gsorted[0][0] + 1
        rec = dict(n_exon=n_exon, n_intron=n_exon - 1, strand=strand,
                   exon_lengths=exon_lengths, intron_lengths=intron_lengths, span=span)
        prev = gene_best.get(gene)
        if prev is None or (n_exon, span) > (prev[0], prev[1]):
            gene_best[gene] = (n_exon, span, rec)

    # Also capture genes whose only exon-bearing features have gene= directly (e.g. pseudogene
    # exons with no mRNA parent in TX_TYPES) — gives them n_exon from direct exon grouping.
    # Second light pass: group exons by gene= for genes not yet captured.
    direct = collections.defaultdict(list)
    need = None  # we only do this if there are genes missing; cheaper to just always collect
    with open(GFF) as f:
        for line in f:
            if line[0] == "#":
                continue
            p = line.rstrip("\n").split("\t", 9)
            if len(p) < 9 or p[2] != "exon":
                continue
            g = GENE_RE.search(p[8])
            par = PARENT_RE.search(p[8])
            if not g:
                continue
            gene = g.group(1)
            if gene in gene_best:
                continue
            # only exons whose parent is NOT a known transcript (orphan / pseudogene exons)
            if par and par.group(1) in tx2gene:
                continue
            direct[gene].append((int(p[3]), int(p[4]), p[6]))
    for gene, exons in direct.items():
        exons.sort()
        strand = exons[0][2]
        ordered = exons if strand != "-" else list(reversed(exons))
        exon_lengths = [e - s + 1 for (s, e, _) in ordered]
        gsorted = sorted(exons)
        il_g = [gsorted[i + 1][0] - gsorted[i][1] - 1 for i in range(len(gsorted) - 1)]
        intron_lengths = il_g if strand != "-" else list(reversed(il_g))
        n_exon = len(exons)
        span = gsorted[-1][1] - gsorted[0][0] + 1
        gene_best[gene] = (n_exon, span,
                           dict(n_exon=n_exon, n_intron=n_exon - 1, strand=strand,
                                exon_lengths=exon_lengths, intron_lengths=intron_lengths,
                                span=span))

    index = {g: rec for g, (_, _, rec) in gene_best.items()}
    json.dump(index, open(OUT, "w"))
    nint = collections.Counter()
    for r in index.values():
        nint[min(r["n_intron"], 10)] += 1
    print(f"[+] wrote {OUT}: {len(index):,} genes")
    print("  intron-count distribution (capped at 10): " +
          "  ".join(f"{k}:{nint[k]}" for k in sorted(nint)))
    intronless = sum(1 for r in index.values() if r["n_intron"] == 0)
    print(f"  intronless genes (single-exon, retrocopy/pseudogene candidates): {intronless:,} "
          f"({100*intronless/max(len(index),1):.1f}%)")


if __name__ == "__main__":
    main()

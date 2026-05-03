#!/usr/bin/env python3
"""Look up which genes overlap each dropped family's regions.

For each region in a dropped family, find the gorilla annotation gene
records that overlap. Report the gene NAMES (not IDs) so we can tell
TE-bridge artifacts (different gene families) apart from missed real
paralogs (all overlapping the same gene family description).

Usage: lookup_genes.py <dropped_families.tsv> [n_samples]
"""
import re, sys, os, random
from collections import defaultdict

GFF = "/scratch/jxi21/Assembler/GGO_genomic.gff"

if len(sys.argv) < 2:
    sys.exit("usage: lookup_genes.py <dropped_families.tsv> [n_samples]")
TSV = sys.argv[1]
N_SAMPLES = int(sys.argv[2]) if len(sys.argv) > 2 else 10

# Step 1: index gorilla genes by chrom for fast overlap lookup.
genes_by_chrom = defaultdict(list)
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene": continue
        attrs = p[8]
        m_name = re.search(r"Name=([^;]+)", attrs)
        m_desc = re.search(r"description=([^;]+)", attrs)
        name = m_name.group(1) if m_name else "?"
        desc = m_desc.group(1) if m_desc else ""
        genes_by_chrom[p[0]].append((int(p[3]), int(p[4]), p[6], name, desc))

print(f"# Indexed genes from {len(genes_by_chrom)} chromosomes")
print()

def lookup(chrom, start, end):
    """Return list of (name, desc) for genes overlapping [start, end] on chrom."""
    hits = []
    for (gs, ge, _strand, name, desc) in genes_by_chrom.get(chrom, []):
        if gs <= end and start <= ge:
            hits.append((name, desc[:60]))
    return hits

# Step 2: read dropped families filtered for low_primitive_jaccard, sample.
samples = []
with open(TSV) as fh:
    header = next(fh).rstrip("\n").split("\t")
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 11: continue
        if p[10] != "low_primitive_jaccard": continue
        samples.append(p)

# Random sample with fixed seed for reproducibility
random.seed(42)
chosen = random.sample(samples, min(N_SAMPLES, len(samples)))

# Step 3: for each chosen family, look up gene names per region.
for p in chosen:
    fid = p[0]
    n_copies = p[1]
    n_chroms = p[3]
    n_shared = p[6]
    spc = p[7]
    pj = p[9]
    regions = p[5].split(";")
    print(f"=== family {fid}: {n_copies}c on {n_chroms}chr, {n_shared} shared reads ({spc}/copy), jaccard={pj} ===")
    for r in regions:
        m = re.match(r"^([^:]+):(\d+)-(\d+):(.)$", r)
        if not m: continue
        chrom, s, e, strand = m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)
        hits = lookup(chrom, s, e)
        # Show short name only to keep readable
        names_str = ", ".join(f"{n}({d[:30]})" for n, d in hits[:3])
        more = f" +{len(hits)-3} more" if len(hits) > 3 else ""
        print(f"  {chrom}:{s}-{e}:{strand} → {len(hits)} genes: {names_str}{more}")
    print()

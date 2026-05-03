#!/usr/bin/env python3
"""External-validation proxy for rustle's family detection.

For each rustle KEPT family (from v4_families.tsv), look up which gorilla
genes overlap each region, extract a "paralog family stem" from each
gene's NCBI description, and compute family coherence:
  coherent (≥0.8 share): true paralog family — passes external validation
  partial (0.5-0.8):     plausible — could be true paralogy + 1-2 outliers
  mixed (<0.5):          probable false positive (TE-bridge or noise)
  no-genes:              pseudogene cluster or unannotated region

Family stems come from NCBI descriptions like:
  "golgin subfamily A member 8 N-like"  → "golgin subfamily A"
  "TBC1 domain family member 3K-like"   → "TBC1 domain family"
  "zinc finger protein 226"             → "zinc finger protein"

Usage:
  family_coherence.py <kept_families.tsv> [<dropped_families.tsv>]
"""
import re, sys
from collections import defaultdict, Counter

GFF = "/scratch/jxi21/Assembler/GGO_genomic.gff"

def family_stem(desc):
    """Extract paralog-family stem from an NCBI gene description.
    Strips trailing modifiers (member numbers, "-like", paralog letters).
    """
    if not desc: return None
    # Drop common trailing modifiers
    s = desc.lower()
    # Strip "-like" / "like" suffixes
    s = re.sub(r"[-\s]+like\s*$", "", s)
    # Drop trailing copy designations: "member 8I", "member 3", "member 3K", etc.
    s = re.sub(r"\s+member\s+\S+\s*$", " family", s)
    # Drop trailing numbers: "alpha 1", "subfamily 7"  (but keep "subfamily X" structure)
    s = re.sub(r"\s+\d+\s*$", "", s)
    # Drop trailing single letter or letter+digit: "A", "B1"
    s = re.sub(r"\s+[a-z]\d?\s*$", "", s)
    # Specific cleanups for common families
    s = re.sub(r"^putative\s+", "", s)
    s = re.sub(r"^uncharacterized\s+\w+$", "uncharacterized", s)
    return s.strip() or None

# Step 1: index genes by chrom (start, end, name, family_stem). Only
# protein-coding genes — multi-copy paralog clusters are typically
# protein-coding, and including ncRNAs / pseudogenes adds noise to the
# gene-overlap lookup.
genes_by_chrom = defaultdict(list)
with open(GFF) as fh:
    for line in fh:
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "gene": continue
        attrs = p[8]
        if "gene_biotype=protein_coding" not in attrs: continue
        m_name = re.search(r"Name=([^;]+)", attrs)
        m_desc = re.search(r"description=([^;]+)", attrs)
        name = m_name.group(1) if m_name else "?"
        desc = m_desc.group(1) if m_desc else ""
        stem = family_stem(desc)
        genes_by_chrom[p[0]].append((int(p[3]), int(p[4]), name, desc[:60], stem))

print(f"# Indexed genes from {len(genes_by_chrom)} chromosomes", file=sys.stderr)

def overlapping_genes(chrom, start, end):
    """Return ALL protein-coding genes overlapping [start,end] on chrom.
    Includes (name, stem) for each.
    """
    out = []
    for (gs, ge, name, desc, stem) in genes_by_chrom.get(chrom, []):
        if gs <= end and start <= ge:
            out.append((name, stem))
    return out

def assess_family(regions):
    """Lenient test: does there EXIST a stem S such that every bundle has
    a gene with stem S? This avoids the "best overlap picks wrong gene"
    pitfall.

    real_paralog:    ≥2 distinct genes overall AND a single stem covers
                     ≥80% of bundles
    single_gene:     all bundles overlap the SAME single gene (or its
                     immediate neighbors with same stem) — segmentation
                     artifact
    te_bridge:       no stem covers ≥80% of bundles
    no_protein:      no protein-coding gene overlaps any bundle
    """
    bundle_genes = []  # per-bundle list of (name, stem) tuples
    for (c, s, e, _) in regions:
        bundle_genes.append(overlapping_genes(c, s, e))
    n_bundles = len(bundle_genes)
    if n_bundles == 0:
        return ("no_protein", 0, 0, 0.0, "")
    # Skip bundles with no gene at all
    bundles_with_gene = [bg for bg in bundle_genes if bg]
    if not bundles_with_gene:
        return ("no_protein", 0, 0, 0.0, "")
    # Find best stem coverage
    all_stems = set()
    for bg in bundles_with_gene:
        for (_n, st) in bg:
            if st: all_stems.add(st)
    best_stem = None
    best_coverage = 0
    for st in all_stems:
        cov = sum(1 for bg in bundles_with_gene if any(s == st for (_n, s) in bg))
        if cov > best_coverage:
            best_coverage = cov
            best_stem = st
    coherence = best_coverage / len(bundles_with_gene)
    # Count distinct genes that share the best_stem (among all bundles)
    distinct_genes_in_stem = set()
    for bg in bundles_with_gene:
        for (n, st) in bg:
            if st == best_stem:
                distinct_genes_in_stem.add(n)
    n_distinct = len(distinct_genes_in_stem)
    if coherence >= 0.8 and n_distinct >= 2:
        category = "real_paralog"
    elif coherence >= 0.8 and n_distinct == 1:
        category = "single_gene"
    else:
        category = "te_bridge"
    return (category, n_distinct, len(all_stems), coherence, best_stem or "")

# Step 2: parse kept families and assess.
def parse_tsv_path(path, label):
    fams = []
    with open(path) as fh:
        header = next(fh).rstrip("\n").split("\t")
        col = {h: i for i, h in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5: continue
            try:
                fid = int(p[0]); n_copies = int(p[1])
            except ValueError: continue
            regions_str = p[col.get("regions", 3)]
            regions = []
            for r in regions_str.split(";"):
                m = re.match(r"^([^:]+):(\d+)-(\d+):(.)$", r)
                if m:
                    regions.append((m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)))
            fams.append({"label": label, "fid": fid, "n_copies": n_copies, "regions": regions})
    return fams

if len(sys.argv) < 2:
    sys.exit("usage: family_coherence.py <kept_families.tsv> [<dropped_families.tsv>]")

kept_fams = parse_tsv_path(sys.argv[1], "kept")
dropped_fams = parse_tsv_path(sys.argv[2], "dropped") if len(sys.argv) > 2 else []

def summarize(fams, label):
    cats = Counter()
    examples = defaultdict(list)
    for f in fams:
        cat, n_genes, n_stems, coh, stem = assess_family(f["regions"])
        cats[cat] += 1
        examples[cat].append((f["fid"], f["n_copies"], n_genes, n_stems, coh, stem))
    print(f"## {label} families: {len(fams)}")
    for cat in ["real_paralog", "single_gene", "te_bridge", "no_protein"]:
        n = cats[cat]
        pct = (100.0 * n / max(len(fams), 1))
        print(f"  {cat:14}: {n:4} ({pct:.1f}%)")
        for ex in examples[cat][:3]:
            fid, nc, ng, ns, coh, stem = ex
            print(f"    e.g. fid={fid} n_copies={nc} n_genes={ng} n_stems={ns} coh={coh:.2f} stem='{stem}'")
    print()

summarize(kept_fams, "KEPT")
if dropped_fams:
    # Filter dropped fams to only low_primitive_jaccard
    pj_drops = []
    with open(sys.argv[2]) as fh:
        header = next(fh).rstrip("\n").split("\t")
        col = {h: i for i, h in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 11 or p[col["drop_reason"]] != "low_primitive_jaccard":
                continue
            try:
                fid = int(p[0]); n_copies = int(p[1])
            except ValueError: continue
            regions_str = p[col["regions"]]
            regions = []
            for r in regions_str.split(";"):
                m = re.match(r"^([^:]+):(\d+)-(\d+):(.)$", r)
                if m:
                    regions.append((m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)))
            pj_drops.append({"fid": fid, "n_copies": n_copies, "regions": regions})
    summarize(pj_drops, "DROPPED-by-low_primitive_jaccard")

#!/usr/bin/env python3
"""The two universes of the Soto-subset study: members, our blocks, Soto's labels, and per-member copy number.

guided  Members are Soto's 2,334 genes (S1C). A gene's locus is its gene line in the F_O input GFF
        (`/mnt/linuxdisk/tmp/ourhom/loci.gff3`, `Name=` = CAT v4 gene id), folded through the run's `ours.loci.tsv`.
        Genes on loci of one F_O family share the block `F:<family>`; genes on one unclustered locus share
        `L:<locus>`.
denovo  Members are the catalog loci that a gene of Soto's 491 families reaches by plan 1's exon + strand rule
        (soto_truth.v2.tsv, CHM13 v2.0; a locus folded into a representative counts as the representative). Blocks are
        `F:<family>` from the catalog's hierarchy members table, `L:<locus>` for a reached locus in no family.
A member's Soto labels are its genes' families; `Unassigned_*` genes are Soto singletons (guided only: they are not
in soto_truth.v2.tsv). Copy number per member = the median over its genes that have a value.
"""
import collections
import csv
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "o1_hier"))
from hier_match import SpanIndex, _strand_ok, load_fold, load_loci_gff, load_nodes  # noqa: E402

Universe = collections.namedtuple("Universe", "block_of soto_sets soto_of multi clustered genes_of coverage")


def s1c_labels(path):
    """S1C -> {gene_id: set(Family ID)}, `Unassigned_*` labels included."""
    out = collections.defaultdict(set)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Family ID"]:
                out[r["Gene ID"]].add(r["Family ID"])
    return dict(out)


def soto_views(labels_of, all_labels=()):
    """{member: set(label)} -> (soto_sets {label: frozenset(member)} without `Unassigned_*`, every label of
    `all_labels` present even when empty; soto_of {member: label} for members with exactly one label; multi = the
    members with >= 2 labels)."""
    sets = {lab: set() for lab in all_labels if not lab.startswith("Unassigned_")}
    for m, labs in labels_of.items():
        for lab in labs:
            if not lab.startswith("Unassigned_"):
                sets.setdefault(lab, set()).add(m)
    soto_of = {m: next(iter(labs)) for m, labs in labels_of.items() if len(labs) == 1}
    multi = frozenset(m for m, labs in labels_of.items() if len(labs) > 1)
    return {lab: frozenset(s) for lab, s in sets.items()}, soto_of, multi


def gff_gene_keys(path):
    """Gene lines of a GFF3 -> {Name: 'chrom:start-end'} (GFF 1-based, the hier_match member key)."""
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            at = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            out[at["Name"]] = f"{f[0]}:{f[3]}-{f[4]}"
    return out


def _block(locus, family_of):
    return f"F:{family_of[locus]}" if locus in family_of else f"L:{locus}"


def guided(s1c_path, loci_gff, members_tsv, fold_tsv):
    labels = s1c_labels(s1c_path)
    key = gff_gene_keys(loci_gff)
    missing = sorted(set(labels) - set(key))
    if missing:
        raise ValueError(f"{len(missing)} S1C genes have no gene line in {loci_gff}: {missing[:5]}")
    _, family_of = load_nodes(members_tsv)
    fold = load_fold(fold_tsv)
    locus_of = {g: fold.get(key[g], key[g]) for g in labels}
    block_of = {g: _block(locus_of[g], family_of) for g in labels}
    soto_sets, soto_of, multi = soto_views(labels)
    clustered = frozenset(g for g, b in block_of.items() if b.startswith("F:"))
    per_locus = collections.Counter(locus_of.values())
    coverage = dict(genes=len(labels), genes_on_family_members=len(clustered), loci=len(per_locus),
                    loci_with_several_genes=sum(1 for c in per_locus.values() if c > 1),
                    families=len({b for b in block_of.values() if b.startswith("F:")}))
    return Universe(block_of, soto_sets, soto_of, multi, clustered, {g: frozenset([g]) for g in labels}, coverage)


def read_truth_v2(path):
    """soto_truth.v2.tsv -> ({gene_id: set(label)}, {gene_id: (chrom, strand, [(s, e)] GFF closed)}, n_records)."""
    labels, genes, n = collections.defaultdict(set), {}, 0
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n += 1
            labels[r["gene_id"]].add(r["label"])
            ex = [tuple(int(x) for x in b.split("-")) for b in r["exons"].split(",")]
            genes[r["gene_id"]] = (r["chrom"], r["strand"], ex)
    return dict(labels), genes, n


def denovo(truth_v2, loci_gff, members_tsv, fold_tsv):
    gene_labels, genes, n_records = read_truth_v2(truth_v2)
    loci = load_loci_gff(loci_gff)
    _, family_of = load_nodes(members_tsv)
    fold = load_fold(fold_tsv)
    idx = SpanIndex([(c, s, e, (st, k)) for k, (c, st, ex) in loci.items() for s, e in ex])
    reach = {g: frozenset(fold.get(k, k) for s, e in ex for st, k in idx.query(chrom, s, e)
                          if _strand_ok(st, {strand}))
             for g, (chrom, strand, ex) in sorted(genes.items())}
    labels_of, genes_of = collections.defaultdict(set), collections.defaultdict(set)
    for g, ms in reach.items():
        for m in ms:
            labels_of[m] |= gene_labels[g]
            genes_of[m].add(g)
    block_of = {m: _block(m, family_of) for m in labels_of}
    all_labels = sorted({lab for labs in gene_labels.values() for lab in labs})
    soto_sets, soto_of, multi = soto_views(labels_of, all_labels)
    clustered = frozenset(m for m in block_of if m in family_of)
    rows = lambda gs: sum(len(gene_labels[g]) for g in gs)  # noqa: E731  (one truth row per gene and label)
    on = [g for g, ms in reach.items() if any(m in family_of for m in ms)]
    unclustered = [g for g, ms in reach.items() if ms and not any(m in family_of for m in ms)]
    coverage = dict(genes=len(genes), records=n_records, records_on_members=rows(on),
                    records_on_unclustered_only=rows(unclustered),
                    records_without_locus=rows(g for g, ms in reach.items() if not ms),
                    reached_loci=len(block_of), reached_clustered=len(clustered))
    return Universe(block_of, soto_sets, soto_of, multi, clustered,
                    {m: frozenset(gs) for m, gs in genes_of.items()}, coverage)


def member_cn(gene_cn, genes_of):
    """{gene: value} -> {member: median over its genes that have a value}; members with none are left out, and
    values of genes outside the universe are never used."""
    out = {}
    for m, gs in genes_of.items():
        v = [gene_cn[g] for g in sorted(gs) if g in gene_cn]
        if v:
            out[m] = float(statistics.median(v))
    return out

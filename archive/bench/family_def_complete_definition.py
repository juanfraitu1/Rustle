#!/usr/bin/env python3
"""family_def_complete_definition.py — the COMPLETE, property-based family definition.

A multi-copy gene family is a structural fact: a maximal set of gene copies whose models share a
WHOLE-GENE backbone (~B). Expression is not part of the definition — it is just whether we OBSERVE
the family in this sample. Two properties GRADE each family:
  identifiability (~R, in-silico): would reads be CONFUSABLE (hard, copy-assignment needed) or
      RESOLVABLE (easy)?  proxied by the copies' divergence (id >= 0.90 => read-confusable).
  expression: EXPRESSED (>=1 member overlaps a de-novo locus) or SILENT.

The read-based 196 families are exactly the HARD ∩ EXPRESSED corner. This script enumerates ALL
structural families (so silent / easy ones are included too) and reports the easy/hard x
expressed/silent partition.

Structural ~B edge: id >= 0.80 (asm20-alignable) AND reciprocal coverage min(cov_a,cov_b) >= 0.30
(whole-gene, not a sub-domain), with the retrocopy + antisense precision filters.
Run: /home/juanfra/miniforge3/bin/python bench/family_def_complete_definition.py
"""
import collections
import os
import sys

import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from family_def_read_filters import dna_homology
from family_def_retrocopy_filter import load_intron_index, n_intron
from family_def_strand_probe import edge_strand
from family_def_genomewide import best_gene, GENES_BED

META = "/home/juanfra/winloci_scratch/denovo_transcripts.meta.tsv"
ID_B = 0.80         # asm20 backbone floor
COV_B = 0.30        # reciprocal whole-gene coverage (~B threshold tau)
ID_R = 0.90         # read-confusable (~R potential) threshold


def expressed_genes():
    by = collections.defaultdict(list)
    with open(GENES_BED) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 4:
                by[p[0]].append((int(p[1]), int(p[2]), p[3]))
    for c in by:
        by[c].sort()
    exp = set()
    with open(META) as f:
        next(f)
        for line in f:
            p = line.rstrip("\n").split("\t")
            g = best_gene(by, p[1], int(p[2]), int(p[3]))
            if g:
                exp.add(g)
    return exp


def main():
    H, _ = dna_homology()
    idx = load_intron_index()
    strand = edge_strand()

    # structural ~B graph (id>=0.80, reciprocal whole-gene coverage), retrocopy + antisense filtered
    G = nx.Graph()
    for (a, b), r in H.items():
        if r.get("id", 0) < ID_B:
            continue
        if min(r["cov_a"], r["cov_b"]) < COV_B:
            continue
        na, nb = n_intron(idx, a), n_intron(idx, b)
        if na is not None and nb is not None and ((na == 0) ^ (nb == 0)):
            continue                                   # one-side-intronless retrocopy
        if strand.get((a, b) if a < b else (b, a)) == "-":
            continue                                   # antisense artifact
        G.add_edge(a, b)
    fams = [c for c in nx.connected_components(G) if len(c) >= 2]
    fams.sort(key=len, reverse=True)

    exp = expressed_genes()

    # grade each family
    def is_hard(c):
        for u in c:
            for v in c:
                if u < v:
                    r = H.get((u, v))
                    if r and r.get("id", 0) >= ID_R:
                        return True
        return False
    def is_expressed(c):
        return any(g in exp for g in c)

    cells = collections.Counter()
    sizes = collections.defaultdict(list)
    for c in fams:
        k = ("hard" if is_hard(c) else "easy", "expressed" if is_expressed(c) else "silent")
        cells[k] += 1
        sizes[k].append(len(c))

    tot = len(fams)
    print(f"=== COMPLETE family definition: {tot} structural families (~B, id>=0.80 whole-gene) ===")
    print(f"  edges {G.number_of_edges()}, member loci {sum(len(c) for c in fams)}, "
          f"max family {len(fams[0])}")
    print(f"\n  graded by identifiability (~R) x expression:")
    print(f"  {'':12} {'EXPRESSED':>14} {'SILENT':>10}   (rows: read-identifiability)")
    for r in ("hard", "easy"):
        he, hs = cells[(r, "expressed")], cells[(r, "silent")]
        tag = "read-CONFUSABLE (copy-assignment)" if r == "hard" else "read-RESOLVABLE (trivial)"
        print(f"  {r.upper():12} {he:>14} {hs:>10}   {tag}")
    print(f"\n  HARD ∩ EXPRESSED = {cells[('hard','expressed')]}  "
          f"(the copy-assignment families — what the read-only pipeline observes)")
    print(f"  + EASY ∩ EXPRESSED = {cells[('easy','expressed')]}  (divergent paralogs reads resolve — newly INCLUDED)")
    print(f"  + SILENT (hard {cells[('hard','silent')]} + easy {cells[('easy','silent')]}) "
          f"= {cells[('hard','silent')]+cells[('easy','silent')]}  "
          f"(structurally families; not expressed in this sample — included by the property-based definition)")
    print(f"\n  => the definition now covers ALL {tot} structural families; expression & "
          f"identifiability are LABELS, not gates.")


if __name__ == "__main__":
    main()

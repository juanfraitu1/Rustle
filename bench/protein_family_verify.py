#!/usr/bin/env python
"""Protein-level ORTHOGONAL verifier for the multi-copy gene-family catalog (SEPARATE from the RNA definition).

The thesis defines families at the RNA level (read-conflict + POA core, nucleotide). This is a STANDALONE
protein harness — one of three independent "levels" (RNA / DNA / protein) that should each stand on their own:
  - RNA   = the read-conflict + POA-core catalog (the contribution).
  - DNA   = bench/soto_family_validate.py (Soto-2025 shared-exon + famCN) / SEDEF segdup map.
  - PROTEIN = THIS: 6-frame longest-ORF translation -> mmseqs all-vs-all -> protein homology graph.

It does NOT touch the RNA pipeline. It reads the catalog's emitted SPLICED copy sequences (`<prefix>.copies.fa`,
headers `FAMID|copyIdx|locus|strand|nexon`), translates each copy's longest ORF, runs one mmseqs all-vs-all, and
reports the protein view:

  (1) PRECISION   — what fraction of RNA families are a single connected component in the protein-homology graph
                    (their copies are mutual protein-paralogs, fident>=MIN_FIDENT & min(qcov,tcov)>=MIN_COV) —
                    an orthogonal confirmation that the family is real, independent of the RNA read evidence.
  (2) REACH       — the within-family protein-identity (fident) distribution: how DIVERGENT protein still
                    confirms paralogy (protein reaches synonymous-divergent copies nucleotide cannot, e.g.
                    ~70% genomic but >=87% protein) — the protein analogue of the ~87% nucleotide ceiling.
  (3) EXTENSION   — CROSS-family protein homology: copies in DIFFERENT RNA families (or unattached loci) that
                    are protein-paralogs = loosely-related members the nucleotide RNA definition split apart.

Mirrors the in-pipeline protein tier (`batch_protein_edges`): longest ORF >=40 aa, mmseqs easy-search -s 7.5,
fident>=0.50, min(qcov,tcov)>=0.50. Independent of annotation (de-novo ORFs), so it is a true orthogonal check.

Run: /home/juanfra/miniforge3/bin/python bench/protein_family_verify.py [catalog_prefix]   (default gw_off)
Needs: mmseqs (/home/juanfra/miniforge3/bin/mmseqs), biopython, the genome at GGO_FASTA.
"""
import os
import subprocess
import sys
import tempfile

import pysam
from Bio.Seq import Seq

SCRATCH = "/home/juanfra/winloci_scratch"
GENOME = os.environ.get("GGO_FASTA", f"{SCRATCH}/GGO.fasta")
MMSEQS = os.environ.get("RUSTLE_MMSEQS", "/home/juanfra/miniforge3/bin/mmseqs")
PREFIX = sys.argv[1] if len(sys.argv) > 1 else f"{SCRATCH}/gw_off"
MIN_AA = 40
MIN_FIDENT = 0.50
MIN_COV = 0.50


def longest_orf_aa(nt: str) -> str:
    """Longest stop-free peptide over all 6 frames (>= MIN_AA aa), as the in-pipeline `longest_orf_aa`."""
    best = ""
    s = Seq(nt.upper().replace("N", "A"))
    for seq in (s, s.reverse_complement()):
        for f in range(3):
            aa = str(seq[f:].translate())
            for pep in aa.split("*"):
                if len(pep) > len(best):
                    best = pep
    return best if len(best) >= MIN_AA else ""


def load_families(prefix):
    """famid -> list of (copy_idx, spliced_seq) from the catalog's emitted SPLICED copies.fa."""
    fams = {}
    fa = pysam.FastaFile(f"{prefix}.copies.fa")
    for ref in fa.references:
        parts = ref.split("|")              # FAMID | copyIdx | locus | strand | nexon=N
        fid, ci = parts[0], int(parts[1])
        fams.setdefault(fid, []).append((ci, fa.fetch(ref)))
    return {k: sorted(v) for k, v in fams.items() if len(v) >= 2}


def components(nodes, edges):
    """Connected components (union-find) over `nodes` with undirected `edges`."""
    parent = {n: n for n in nodes}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edges:
        if a in parent and b in parent:
            parent[find(a)] = find(b)
    comp = {}
    for n in nodes:
        comp.setdefault(find(n), []).append(n)
    return list(comp.values())


def main():
    fams = load_families(PREFIX)
    print(f"catalog: {os.path.basename(PREFIX)} — {len(fams)} families (>=2 copies)")

    # translate each copy's longest ORF from its SPLICED sequence; node id = "f{fam}_c{pos}"
    fam_ids = list(fams)
    node_fam = {}            # node id -> family id
    fam_nodes = {}           # family id -> [node ids with an ORF]
    with tempfile.TemporaryDirectory() as td:
        qf = os.path.join(td, "q.fa")
        n_orf = 0
        with open(qf, "w") as out:
            for fi, fid in enumerate(fam_ids):
                fam_nodes[fid] = []
                for pos, (_ci, seq) in enumerate(fams[fid]):
                    aa = longest_orf_aa(seq)
                    if aa:
                        nid = f"f{fi}_c{pos}"
                        out.write(f">{nid}\n{aa}\n")
                        node_fam[nid] = fid
                        fam_nodes[fid].append(nid)
                        n_orf += 1
        print(f"translated {n_orf} ORFs (>= {MIN_AA} aa) from spliced copies; running mmseqs all-vs-all (-s 7.5)...")
        res = os.path.join(td, "res.m8")
        r = subprocess.run(
            [MMSEQS, "easy-search", qf, qf, res, os.path.join(td, "tmp"),
             "-s", "7.5", "--threads", "4", "--format-output", "query,target,fident,qcov,tcov"],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            sys.exit(f"mmseqs failed:\n{r.stderr[-1500:]}")
        within_edges = []          # (node,node) within the same family, passing thresholds
        within_fident = []
        cross = {}                 # (famA,famB) -> max fident, cross-family protein homology
        for line in open(res):
            q, t, fident, qcov, tcov = line.rstrip("\n").split("\t")
            if q == t:
                continue
            fident, qcov, tcov = float(fident), float(qcov), float(tcov)
            if fident < MIN_FIDENT or min(qcov, tcov) < MIN_COV:
                continue
            fq, ft = node_fam.get(q), node_fam.get(t)
            if fq is None or ft is None:
                continue
            if fq == ft:
                within_edges.append((q, t))
                within_fident.append(fident)
            else:
                key = tuple(sorted((fq, ft)))
                cross[key] = max(cross.get(key, 0.0), fident)

    # (1) PRECISION: a family is protein-confirmed if its >=2 ORF-bearing copies form ONE protein component.
    confirmed = testable = 0
    for fid in fam_ids:
        nodes = fam_nodes[fid]
        if len(nodes) < 2:
            continue
        testable += 1
        nset = set(nodes)
        comps = components(nodes, [(a, b) for (a, b) in within_edges if a in nset and b in nset])
        if len(comps) == 1:
            confirmed += 1

    print(f"\n(1) PRECISION — RNA families confirmed as a single protein-paralog component "
          f"(fident>={MIN_FIDENT}, cov>={MIN_COV}):")
    print(f"    {confirmed}/{testable} testable families ({len(fams)-testable} have <2 ORF-bearing copies)")

    if within_fident:
        wf = sorted(within_fident)
        import statistics
        print(f"\n(2) REACH — within-family protein identity (fident) of confirming edges:")
        print(f"    median {100*statistics.median(wf):.0f}%, min {100*min(wf):.0f}%, "
              f"<70%: {sum(1 for x in wf if x<0.70)}, <50.. : down to the {MIN_FIDENT*100:.0f}% gate")
        print(f"    (protein confirms paralogy well BELOW the nucleotide ~81-87% ceiling — the orthogonal reach)")

    print(f"\n(3) EXTENSION — CROSS-family protein homology (loosely-related members the RNA definition split):")
    print(f"    {len(cross)} family-pairs are protein-paralogs (fident>={MIN_FIDENT})")
    for (a, b), fi in sorted(cross.items(), key=lambda kv: -kv[1])[:10]:
        print(f"      {a} ~ {b}  protein fident {100*fi:.0f}%")
    print("\nInterpretation: (1) is the orthogonal CONFIRMATION (protein agrees the RNA family is real), (2) is")
    print("the protein REACH (down to the fident gate, far past nucleotide), (3) flags loosely-related paralogs")
    print("the nucleotide RNA definition could not merge — the protein level's unique contribution. SEPARATE from")
    print("the RNA pipeline by construction; pair with the DNA harness (soto_family_validate.py / SEDEF).")


if __name__ == "__main__":
    main()
